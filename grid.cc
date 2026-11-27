// The model grid and the propagation grid. This file:
// - maps the model cells onto the propagation grid;
// - allocates and holds the state of each cell;
// - assigns the cells to the ranks;
// - reads and writes the grid restart data.
// grid_modelinput.cc reads the input ejecta model, and grid_geometry.cc holds the geometry of the propagation grid.

#include "grid.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <optional>
#include <print>
#include <ranges>
#include <span>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid_internal.h"
#include "kpkt.h"
#include "mpi_logging.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "outputfilestream.h"
#include "radfield.h"
#include "rpkt.h"
#include "sn3d.h"

namespace grid {

namespace {

ptrdiff_t nonempty_npts_model = 0;  // number of allocated non-empty model grid cells

// the model cell of each propagation cell, or -1 for a propagation cell with no matter, e.g. a corner cell outside
// the outermost shell of a 1D or 2D model
std::vector<int> propcell_mgi;
std::vector<int> propcell_nonemptymgi;

std::vector<int> mgi_of_nonemptymgi;

std::vector<double> totmassnuclide{};  // total mass of each nuclide in the ejecta
// the grid-mapped mass of each nuclide in the propagation cells that FORCE_SPHERICAL_ESCAPE_SURFACE
// blanks. The nuclide-mass rescale in init_grid() adds this to the mapped total. The rescale then
// corrects only the discretisation error, and the blanked mass stays removed.
std::vector<double> totmassnuclide_blanked{};

MPI_shared_array<float> initmassfracuntrackedstable_allcells{};
MPI_shared_array<int> elements_uppermost_ion_allcells{};  // Highest ion index that has a significant population
// Lowest ion index that has a significant population. The NLTE ion range reduction can raise it
// above zero, and the paths that replace the NLTE solution with Saha populations reset it to zero.
MPI_shared_array<int> elements_lowermost_ion_allcells{};

// indexed by global rank
std::vector<int> ranks_nstart;
std::vector<int> ranks_nstart_nonempty;
std::vector<int> ranks_ndo;
std::vector<int> ranks_ndo_nonempty;

void set_propcell_modelgridindex(const int cellindex, const int new_modelgridindex) {
  assert_testmodeonly(cellindex >= 0);
  assert_testmodeonly(cellindex < ngrid);
  assert_testmodeonly(new_modelgridindex >= -1);
  assert_testmodeonly(new_modelgridindex < get_npts_model());
  propcell_mgi[cellindex] = new_modelgridindex;
}

// Split each element's mass fraction into its tracked isotopes and an untracked stable remainder, for every
// non-empty cell. Must run after any mapping rescale of the nuclide mass fractions, so that the remainder is
// the one left over from the values everything downstream will use. read_elem_abundances() has already checked
// that the input gives each element at least as much mass as its isotopes, so a negative remainder here is an
// artefact of the rescale rather than bad input, and is clamped away.
void set_untrackedstable_massfracs(std::span<const float> elem_massfracs_of_mgi) {
  if (globals::rank_in_node != 0) {
    return;
  }
  for (ptrdiff_t nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
    const auto tracked_isotope_massfrac_sum_by_atomic_number = get_tracked_isotope_massfrac_sum_by_atomic_number(mgi);
    for (int element = 0; element < get_nelements(); element++) {
      const int atomic_number = get_atomicnumber(element);
      const double elem_massfrac = elem_massfracs_of_mgi[(static_cast<size_t>(mgi) * get_nelements()) + element];
      const double massfrac_allisotopes = tracked_isotope_massfrac_sum_by_atomic_number[atomic_number];

      double massfrac_untrackedstable = elem_massfrac - massfrac_allisotopes;
      if (massfrac_untrackedstable < 0.) {
        //  allow some roundoff error before we complain
        if (std::abs(massfrac_allisotopes - elem_massfrac) > 1e-4) {
          printlnlog("[warning] cell {} Z={} element massfrac is less than the sum of its radioisotope massfracs", mgi,
                     atomic_number);
          printlnlog("  massfrac(Z) {:g} massfrac_radioisotopes(Z) {:g}", elem_massfrac, massfrac_allisotopes);
          printlnlog("  increasing elemental massfrac to {:g} and setting stable isotopic massfrac to zero",
                     massfrac_allisotopes);
        }
        massfrac_untrackedstable = 0.;  // bring up to zero if negative
      }

      initmassfracuntrackedstable_allcells[(nonemptymgi * get_nelements()) + element] =
          static_cast<float>(massfrac_untrackedstable);

      // (massfrac_allisotopes + massfrac_untrackedstable) might not exactly match elem_massfrac if we had to boost
      // it to reach massfrac_allisotopes
      set_elem_massfrac(nonemptymgi, element, static_cast<float>(massfrac_allisotopes + massfrac_untrackedstable));
    }
  }
}

void allocate_nonemptycells_composition_cooling() {
  // allocate the composition, cooling, and NLTE population arrays of the nonempty cells
  const ptrdiff_t nonempty_npts_model_ptrdifft = get_nonempty_npts_model();
  const auto nelements = get_nelements();

  initmassfracuntrackedstable_allcells = MPI_shared_array<float>(nonempty_npts_model_ptrdifft * nelements, 0.);
  elem_meanweight_allcells = MPI_shared_array<float>(nonempty_npts_model_ptrdifft * nelements, 0.);
  elements_uppermost_ion_allcells = MPI_shared_array<int>(nonempty_npts_model_ptrdifft * nelements, -1);
  if constexpr (NLTE_TRACK_SOLUTION_RANGES) {
    elements_lowermost_ion_allcells = MPI_shared_array<int>(nonempty_npts_model_ptrdifft * nelements, 0);
    printlnlog("[info] mem_usage: the lowermost ions of the NLTE solution ranges occupy {:.3f} MB (node shared memory)",
               static_cast<double>(nonempty_npts_model_ptrdifft * nelements * sizeof(int)) / 1024. / 1024.);
  }
  elem_massfracs_allcells = MPI_shared_array<float>(nonempty_npts_model_ptrdifft * nelements, 0.);
  ion_groundlevelpops_allcells = MPI_shared_array<float>(nonempty_npts_model_ptrdifft * get_includedions(), 0.);
  ion_partfuncts_allcells = MPI_shared_array<float>(nonempty_npts_model_ptrdifft * get_includedions(), 0.);
  kpkt::ion_cooling_contribs_allcells = MPI_shared_array<double>(nonempty_npts_model_ptrdifft * get_includedions(), 0.);
  // 1 means that the k-packets keep their full energy (no thermal balance yet)
  kpkt::radiative_energy_factor_allcells = MPI_shared_array<double>(nonempty_npts_model_ptrdifft, 1.);

  // -1 indicates that there is currently no information on the nlte populations
  nltepops_allcells = MPI_shared_array<double>(nonempty_npts_model_ptrdifft * globals::total_nlte_levels, -1.);

  if constexpr (NLTE_TIME_DEPENDENT_FIRST_TIMESTEP.has_value()) {
    nltepop_allocate_time_dependent_arrays();
  }
}

// build the mapping between model cells and propagation cells, then allocate the per-cell
// shared-memory arrays (populations, estimators, temperatures) for the nonempty model cells only
void allocate_nonemptymodelcells() {
  // Determine the number of simulation cells associated with the model cells
  std::ranges::fill(modelgrid_numpropcells, 0);
  if (globals::rank_in_node == 0) {
    for (int mgi = 0; mgi < get_npts_model(); mgi++) {
      modelgrid_input[mgi].initial_radial_pos_sum = 0.;
      modelgrid_input[mgi].initial_radial_pos_squared_sum = 0.;
    }
  }
  MPI_Barrier_node();

  // Record the grid-mapped volume that the spherical escape surface blanks per model cell. This
  // function must record the per-nuclide blanked mass itself. The loop below zeroes the density and
  // the mass fractions of a fully blanked model cell. Only the model-to-grid mapping case needs the
  // record, because only that case rescales the nuclide masses.
  const bool track_blanked_mass = FORCE_SPHERICAL_ESCAPE_SURFACE && (get_modelgridtype() != get_propgridtype());
  std::vector<double> blanked_assocvolume(track_blanked_mass ? get_npts_model() : 0, 0.);
  totmassnuclide_blanked.assign(decay::get_num_nuclides(), 0.);

  if constexpr (FORCE_SPHERICAL_ESCAPE_SURFACE) {
    reserve_resize(propcell_outside_escape_surface, ngrid);
  }

  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const auto radial_pos_mid = get_cellradialposmid(cellindex);

    if constexpr (FORCE_SPHERICAL_ESCAPE_SURFACE) {
      propcell_outside_escape_surface[cellindex] =
          (get_cell_r_inner(cellindex, get_propgridtype()) > globals::rmax) ? 1 : 0;
    }

    if (FORCE_SPHERICAL_ESCAPE_SURFACE && radial_pos_mid > globals::vmax * globals::tmin) {
      // for 1D models, the final shell outer v should already be at vmax
      assert_always(model_type != GridType::SPHERICAL1D || propcell_mgi[cellindex] < 0);
      if (const int mgi_blanked = get_propcell_modelgridindex(cellindex); track_blanked_mass && mgi_blanked >= 0) {
        blanked_assocvolume[mgi_blanked] += get_propcell_volume_tmin(cellindex);
      }
      set_propcell_modelgridindex(cellindex, -1);
    }

    const int mgi = get_propcell_modelgridindex(cellindex);
    if (mgi >= 0) {
      modelgrid_numpropcells[mgi] += 1;
      if (globals::rank_in_node == 0) {
        modelgrid_input[mgi].initial_radial_pos_sum += radial_pos_mid;
        modelgrid_input[mgi].initial_radial_pos_squared_sum += get_cellradialposmeansquared(cellindex);
      }
      assert_always(get_rho_tmin(mgi) > 0);
      // with direct mapping, a non-empty model cell has one propagation cell, or two with a mirror plane
      assert_always(get_modelgridtype() != get_propgridtype() ||
                    modelgrid_numpropcells[mgi] <= (model_mirror_z ? 2 : 1));
    }
  }

  if (track_blanked_mass) {
    // every rank computes the same values from the same node-shared inputs, so no communication is needed
    for (int mgi = 0; mgi < get_npts_model(); mgi++) {
      if (blanked_assocvolume[mgi] > 0.) {
        const double mass_blanked = get_rho_tmin(mgi) * blanked_assocvolume[mgi];
        for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
          totmassnuclide_blanked[nucindex] += mass_blanked * get_modelinitnucmassfrac(mgi, nucindex);
        }
      }
    }
  }

  MPI_Barrier_node();
  // find number of non-empty cells and allocate nonempty list
  nonempty_npts_model = 0;
  for (int mgi = 0; mgi < get_npts_model(); mgi++) {
    if (get_numpropcells(mgi) > 0) {
      nonempty_npts_model++;
    }
  }
  if (nonempty_npts_model == 0) {
    fatal_crash(
        "no model cell with matter maps onto the propagation grid. A model cell keeps its matter only with a positive "
        "density in model.txt and with mass in an element of compositiondata.txt. A propagation cell gets the model "
        "cell that holds the centre of the propagation cell");
  }
  printlnlog("There are {} modelgrid cells with associated propagation cells (nonempty_npts_model)",
             nonempty_npts_model);

  mgi_of_nonemptymgi.assign(nonempty_npts_model, -2);

  propcell_nonemptymgi.assign(ngrid, -1);

  int nonemptymgi = 0;  // index within list of non-empty modelgrid cells

  for (int mgi = 0; mgi < get_npts_model(); mgi++) {
    if (get_numpropcells(mgi) > 0) {
      nonemptymgi_of_mgi[mgi] = nonemptymgi;
      mgi_of_nonemptymgi[nonemptymgi] = mgi;
      nonemptymgi++;
    } else {
      nonemptymgi_of_mgi[mgi] = -1;
      if (globals::rank_in_node == 0) {
        // these arrays are in node-shared memory, so only the node leader writes them
        // (synchronised by the barrier below)
        set_rho_tmin(mgi, 0.);
        for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
          set_modelinitnucmassfrac(mgi, nucindex, 0.);
        }
      }
    }
  }
  MPI_Barrier_allranks();
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    if (const int mgi = get_propcell_modelgridindex(cellindex); mgi >= 0) {
      propcell_nonemptymgi[cellindex] = get_nonemptymgi_of_mgi(mgi);
    }
  }

  assert_always(rho_allcells.empty());
  rho_allcells = MPI_shared_array<float>(nonempty_npts_model, -1.);
  Te_allcells = MPI_shared_array<float>(nonempty_npts_model, -1.);
  TJ_allcells = MPI_shared_array<float>(nonempty_npts_model, -1.);
  TR_allcells = MPI_shared_array<float>(nonempty_npts_model, -1.);
  W_allcells = MPI_shared_array<float>(nonempty_npts_model, -1.);
  nne_allcells = MPI_shared_array<float>(nonempty_npts_model, -1.);
  nnetot_allcells = MPI_shared_array<float>(nonempty_npts_model, -1.);
  kappagrey_allcells = MPI_shared_array<float>(nonempty_npts_model, 0.);
  grey_depth_allcells = MPI_shared_array<float>(nonempty_npts_model, 0.);
  thick_allcells = MPI_shared_array<CellThickness>(nonempty_npts_model, CellThickness::THIN);
  if constexpr (USE_MICROCLUMPING) {
    clumpfactor_allcells = MPI_shared_array<float>(nonempty_npts_model, -1.);
  }
  const auto modelgrid_mem_usage =
      nonempty_npts_model * ((sizeof(float) * (USE_MICROCLUMPING ? 10 : 9)) + sizeof(double) + sizeof(int));
  printlnlog(
      "[info] mem_usage: the modelgrid properties (temperatures and electron densities) occupies {:.3f} MB (node "
      "shared memory)",
      modelgrid_mem_usage / 1024. / 1024.);

  allocate_nonemptycells_composition_cooling();

  if constexpr (expopac_bins_on) {
    allocate_expansionopacities();
  }

  globals::dep_estimator_gamma.assign(nonempty_npts_model, 0.);

  globals::dep_estimator_positron.assign(nonempty_npts_model, 0.);

  globals::dep_estimator_electron.assign(nonempty_npts_model, 0.);

  globals::dep_estimator_alpha.assign(nonempty_npts_model, 0.);

  const auto ionestimcount = nonempty_npts_model * globals::nbfcontinua_ground;
  const auto ionestimsize = ionestimcount * sizeof(double);

  if (ionestimsize > 0) {
    globals::corrphotoionrenorm = MPI_shared_array<double>(ionestimcount, 1.);

    globals::gammaestimator.assign(ionestimcount, 0.);
    globals::col_iongamma_per_gspop_per_clumpednne.assign(ionestimcount, 0.);
  } else {
    globals::corrphotoionrenorm.reset();
    globals::gammaestimator.clear();
    globals::col_iongamma_per_gspop_per_clumpednne.clear();
  }

  if (USE_ION_BFHEATING_ESTIMATORS && ionestimsize > 0) {
    globals::bfheatingestimator.assign(ionestimcount, 0.);
  } else {
    globals::bfheatingestimator.clear();
  }

  globals::ffheatingestimator.assign(nonempty_npts_model, 0.);

  globals::colheatingestimator.assign(COLLISIONAL_HEATING_FROM_LEVELPOPS ? 0 : nonempty_npts_model, 0.);

  MPI_Barrier_allranks();

  printlnlog(
      "[info] mem_usage: NLTE populations for all allocated cells occupy a total of {:.3f} MB (node shared memory)",
      static_cast<ptrdiff_t>(get_nonempty_npts_model()) * globals::total_nlte_levels * sizeof(double) / 1024. / 1024.);
}

// Map 1D spherical model grid onto 3D Cartesian propagation grid
void map_1dmodelto3dgrid() {
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const double cellvmid = get_cellradialposmid(cellindex) / globals::tmin;
    const int mgi = int_index_lowerbound(vout_model, cellvmid);

    if (mgi < get_npts_model() && modelgrid_input[mgi].rhoinit > 0) {
      set_propcell_modelgridindex(cellindex, mgi);
      assert_always(vout_model[mgi] >= cellvmid);
      assert_always((mgi > 0 ? vout_model[mgi - 1] : 0.0) <= cellvmid);
    } else {
      // corner cells outside of the outermost model shell are empty and so are any shells with zero density
      set_propcell_modelgridindex(cellindex, -1);
    }
  }
}

// Map 2D cylindrical model onto 3D Cartesian propagation grid
void map_2dmodelto3dgrid() {
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    // map to 3D Cartesian grid
    const auto pos_mid = Vec3d{
        (get_cellcoordmin(cellindex, 0) + get_cellcoordmax(cellindex, 0)) / 2,
        (get_cellcoordmin(cellindex, 1) + get_cellcoordmax(cellindex, 1)) / 2,
        (get_cellcoordmin(cellindex, 2) + get_cellcoordmax(cellindex, 2)) / 2,
    };

    const double rcylindrical = std::sqrt(pow2(pos_mid[0]) + pow2(pos_mid[1]));

    // 2D grid is uniform so rcyl and z indices can be calculated with no lookup
    const int n_rcyl = static_cast<int>(rcylindrical / globals::tmin / globals::vmax * ncoord_model[0]);
    // the z axis of a mirrored model starts at the plane z = 0, and a cell with z < 0 takes its mirror cell
    const int n_z =
        model_mirror_z
            ? static_cast<int>(std::abs(pos_mid[2]) / globals::tmin / globals::vmax * ncoord_model[1])
            : static_cast<int>(((pos_mid[2] / globals::tmin) + globals::vmax) / (2 * globals::vmax) * ncoord_model[1]);

    if (n_rcyl >= 0 && n_rcyl < ncoord_model[0] && n_z >= 0 && n_z < ncoord_model[1]) {
      const int mgi = (n_z * ncoord_model[0]) + n_rcyl;

      if (modelgrid_input[mgi].rhoinit > 0) {
        set_propcell_modelgridindex(cellindex, mgi);
      } else {
        set_propcell_modelgridindex(cellindex, -1);
      }
    } else {
      set_propcell_modelgridindex(cellindex, -1);
    }
  }
}

// Map a 2D spherical or a 3D spherical model onto the 3D Cartesian propagation grid. Each propagation cell takes the
// model cell that holds its midpoint.
void map_sphericalmodelto3dgrid() {
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const auto pos_mid = Vec3d{
        (get_cellcoordmin(cellindex, 0) + get_cellcoordmax(cellindex, 0)) / 2,
        (get_cellcoordmin(cellindex, 1) + get_cellcoordmax(cellindex, 1)) / 2,
        (get_cellcoordmin(cellindex, 2) + get_cellcoordmax(cellindex, 2)) / 2,
    };
    const int mgi = get_spherical_model_cellindex_tmin(pos_mid);
    // the corner cells outside the outermost shell are empty, and so is a cell of a model cell with zero density
    set_propcell_modelgridindex(cellindex, (mgi >= 0 && modelgrid_input[mgi].rhoinit > 0) ? mgi : -1);
  }
}

// Get the propagation cell with z >= 0 of a model cell of a mirrored model on its own grid type. The mirror axis is
// the axis 1: z of a cylindrical model and cos(theta) of a spherical model.
[[nodiscard]] auto get_upper_propcell_of_mirrored_mgi(const int mgi) -> int {
  const int index0 = mgi % ncoord_model[0];
  const int index1 = (mgi / ncoord_model[0]) % ncoord_model[1];
  const int index2 = mgi / (ncoord_model[0] * ncoord_model[1]);
  // the cells with z >= 0 are the second half of the mirror axis of the propagation grid
  const int propindex1 = ncoord_model[1] + index1;
  return index0 + (ncoordgrid[0] * (propindex1 + (ncoordgrid[1] * index2)));
}

// Map a mirrored model onto a propagation grid of its own type. The propagation grid has the cells of both halves,
// and each propagation cell takes the model cell of the mirror pair.
void map_mirroredmodeltogrid() {
  const int num_half = ncoord_model[1];
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const int index0 = get_cellcoordindex(cellindex, 0);
    const int propindex1 = get_cellcoordindex(cellindex, 1);
    // the model counts the cells of the mirror axis from the plane z = 0, and the propagation grid from the -z end
    const int index1 = (propindex1 >= num_half) ? propindex1 - num_half : num_half - 1 - propindex1;
    const int index2 = (get_modelgridtype() == GridType::SPHERICAL3D) ? get_cellcoordindex(cellindex, 2) : 0;
    const int mgi = index0 + (ncoord_model[0] * (index1 + (ncoord_model[1] * index2)));
    set_propcell_modelgridindex(cellindex, (modelgrid_input[mgi].rhoinit > 0) ? mgi : -1);
  }
}

// here mgi and cellindex are interchangeable (except for empty cells that are associated with mgi == -1)
void map_modeltogrid_direct() {
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const int mgi = (modelgrid_input[cellindex].rhoinit > 0) ? cellindex : -1;
    set_propcell_modelgridindex(cellindex, mgi);
  }
}

auto get_inputcellvolume(const int mgi) -> double {
  switch (get_modelgridtype()) {
    case GridType::SPHERICAL1D: {
      const double v_inner = (mgi == 0) ? 0. : vout_model[mgi - 1];
      return (pow3(vout_model[mgi]) - pow3(v_inner)) * 4 * PI * pow3(globals::tmin) / 3.;
    }

    case GridType::CYLINDRICAL2D: {
      const int n_r = mgi % ncoord_model[0];
      const double delta_rcyl = globals::vmax * t_model / ncoord_model[0];
      // the z axis of a mirrored model is from 0 to vmax, and each cell also fills its mirror cell
      const double delta_z = (model_mirror_z ? 1. : 2.) * globals::vmax * t_model / ncoord_model[1];
      const double mirrorfactor = model_mirror_z ? 2. : 1.;
      return mirrorfactor * pow3(globals::tmin / t_model) * delta_z * PI *
             (pow2((n_r + 1) * delta_rcyl) - pow2(n_r * delta_rcyl));
    }

    case GridType::CARTESIAN3D: {
      // Assumes cells are cubes here - all same volume.
      return pow3(2 * globals::vmax * globals::tmin) / (ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2]);
    }

    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D: {
      // the radial index of a cell changes fastest. A 2D spherical model has one cell of phi from 0 to 2 pi.
      const int index_r = mgi % ncoord_model[0];
      const int index_costheta = (mgi / ncoord_model[0]) % ncoord_model[1];
      const int index_phi = mgi / (ncoord_model[0] * ncoord_model[1]);
      const double radialpart =
          (pow3(model_vel_r_edges[index_r + 1]) - pow3(model_vel_r_edges[index_r])) * pow3(globals::tmin) / 3.;
      const double polarpart = model_costheta_edges[index_costheta + 1] - model_costheta_edges[index_costheta];
      // each cell of a mirrored model also fills its mirror cell
      const double mirrorfactor = model_mirror_z ? 2. : 1.;
      return mirrorfactor * radialpart * polarpart * (model_phi_edges[index_phi + 1] - model_phi_edges[index_phi]);
    }
  }

  fatal_crash("Unknown model grid type {}", static_cast<int>(get_modelgridtype()));
}

void read_grid_restart_data(const int timestep) {
  const auto filename = std::format("gridsave_ts{}.tmp", timestep);

  printlnlog("reading grid restart snapshot from {}", filename);
  FILE* gridsave_file = fopen_required(filename, "rb");

  int ntimesteps_in = -1;
  int nprocs_in = -1;
  double tmin_in = -1.;
  double tmax_in = -1.;
  int timestep_in = -1;
  int nonempty_npts_model_in = -1;
  int nelements_in = -1;
  int includedions_in = -1;
  int nbfcontinua_ground_in = -1;
  read_restart_values(gridsave_file, ntimesteps_in, nprocs_in, tmin_in, tmax_in, timestep_in, nonempty_npts_model_in,
                      nelements_in, includedions_in, nbfcontinua_ground_in);
  if (ntimesteps_in != globals::ntimesteps || nprocs_in != globals::nprocs) {
    fatal_crash("{} has {} timesteps and {} ranks, but this simulation has {} timesteps and {} ranks", filename,
                ntimesteps_in, nprocs_in, globals::ntimesteps, globals::nprocs);
  }
  // the saved per-timestep energies belong to the time grid of the run that wrote the file
  if (tmin_in != globals::tmin || tmax_in != globals::tmax) {
    fatal_crash("{} was written with tmin {:g} tmax {:g} but input.txt gives tmin {:g} tmax {:g}", filename, tmin_in,
                tmax_in, globals::tmin, globals::tmax);
  }
  if (timestep_in != timestep) {
    fatal_crash("{} holds timestep {}, but the expected timestep is {}", filename, timestep_in, timestep);
  }
  if (nonempty_npts_model_in != get_nonempty_npts_model() || nelements_in != get_nelements() ||
      includedions_in != get_includedions() || nbfcontinua_ground_in != globals::nbfcontinua_ground) {
    fatal_crash(
        "{} has {} non-empty cells, {} elements, {} ions and {} ground-level bf continua, but this simulation has {}, "
        "{}, {} and {}",
        filename, nonempty_npts_model_in, nelements_in, includedions_in, nbfcontinua_ground_in,
        get_nonempty_npts_model(), get_nelements(), get_includedions(), globals::nbfcontinua_ground);
  }

  std::vector<int> mgi_of_nonemptymgi_in(mgi_of_nonemptymgi.size());
  read_restart_array(gridsave_file, mgi_of_nonemptymgi_in);
  if (mgi_of_nonemptymgi_in != mgi_of_nonemptymgi) {
    fatal_crash("{} belongs to a model with a different set of non-empty cells", filename);
  }

  for (int element = 0; element < get_nelements(); element++) {
    int atomic_number_in = -1;
    int nions_in = -1;
    read_restart_values(gridsave_file, atomic_number_in, nions_in);
    if (atomic_number_in != get_atomicnumber(element) || nions_in != get_nions(element)) {
      fatal_crash("{} has Z={} with {} ions as element {}, but this simulation has Z={} with {} ions", filename,
                  atomic_number_in, nions_in, element, get_atomicnumber(element), get_nions(element));
    }
  }

  read_restart_array(gridsave_file, globals::timesteps);

  read_restart_array(gridsave_file, TR_allcells);
  read_restart_array(gridsave_file, Te_allcells);
  read_restart_array(gridsave_file, W_allcells);
  read_restart_array(gridsave_file, TJ_allcells);
  read_restart_array(gridsave_file, thick_allcells);
  read_restart_array(gridsave_file, nne_allcells);
  read_restart_array(gridsave_file, nnetot_allcells);
  read_restart_array(gridsave_file, kpkt::radiative_energy_factor_allcells);
  read_restart_array(gridsave_file, globals::dep_estimator_gamma);
  read_restart_array(gridsave_file, globals::dep_estimator_positron);
  read_restart_array(gridsave_file, globals::dep_estimator_electron);
  read_restart_array(gridsave_file, globals::dep_estimator_alpha);
  if constexpr (USE_LUT_PHOTOION) {
    read_restart_array(gridsave_file, globals::corrphotoionrenorm);
    read_restart_array(gridsave_file, globals::gammaestimator);
    read_restart_array(gridsave_file, globals::col_iongamma_per_gspop_per_clumpednne);
  }

  // the order of these calls is very important!
  radfield::read_restart_data(gridsave_file);
  if (globals::rank_in_node == 0) {
    // all data is shared on the node
    nonthermal::read_restart_data(gridsave_file);
    nltepop_read_restart_data(gridsave_file);
    if constexpr (NLTE_TRACK_SOLUTION_RANGES) {
      read_restart_array(gridsave_file, elements_lowermost_ion_allcells);
      read_restart_array(gridsave_file, elements_uppermost_ion_allcells);
    }
    // the file must hold no data after the last section
    assert_always(std::fgetc(gridsave_file) == EOF && std::feof(gridsave_file) != 0);
  }
  MPI_Barrier_node();
  fclose(gridsave_file);
}

// Assign temperatures to the grid cells at the start of the simulation by assuming that all radioactive decay since
// the snapshot time (plus any snapshot initial cell energy) energy is fully trapped
void assign_initial_temperatures() {
  MPI_Barrier_allranks();

  // We assume that for early times the material is so optically thick, that
  // all the radiation is trapped in the cell it originates from. This
  // means furthermore LTE, so that both temperatures can be evaluated
  // according to the local energy density from all radioactive decays and the initial energy q.
  // The dilution factor is W=1 in LTE.

  printlog("Assigning initial temperatures...");

  const double ts0_tmid = globals::timesteps[0].mid;
  int cells_below_mintemp = 0;
  int cells_above_maxtemp = 0;

  // the Bateman factors depend only on the decay path and time, so compute them once and apply them to every
  // cell's initial abundances
  const auto endecay_per_massoftopnuc = decay::calc_energy_per_massoftopnuc_decaypath_withexpansion(ts0_tmid);

  // the decayed energy calculation is expensive and the temperature arrays are in node-shared memory, so stripe
  // the cells across the ranks of each node, with each rank writing its own disjoint subset of the shared arrays
  for (int nonemptymgi = globals::rank_in_node; nonemptymgi < get_nonempty_npts_model();
       nonemptymgi += globals::node_nprocs) {
    const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);

    // Both energies refer to ts0_tmid. The model file gives q at t_model with the adiabatic losses before t_model
    // already applied, and the reader scaled it to tmin. The 1/t law of a comoving mass element takes it from
    // tmin to ts0_tmid. The decay term counts the decays between t_model and ts0_tmid, and a decay at t_decay
    // keeps t_decay / ts0_tmid of its energy under the same law.
    const auto q = INITIAL_PACKETS_ON ? (get_initenergyq(mgi) * globals::tmin / ts0_tmid) : 0.;
    const double decayedenergy_per_mass =
        decay::get_modelcell_endecay_per_mass(nonemptymgi, endecay_per_massoftopnuc) + q;

    auto T_initial = static_cast<float>(std::pow(
        CLIGHT / 4 / STEBO * pow3(globals::tmin / ts0_tmid) * get_rho_tmin(mgi) * decayedenergy_per_mass, 1. / 4.));

    assert_always(std::isfinite(T_initial));
    if (T_initial < MINTEMP) {
      T_initial = MINTEMP;
      cells_below_mintemp++;
    } else if (T_initial > MAXTEMP) {
      T_initial = MAXTEMP;
      cells_above_maxtemp++;
    }

    Te_allcells[nonemptymgi] = T_initial;
    TJ_allcells[nonemptymgi] = T_initial;
    TR_allcells[nonemptymgi] = T_initial;
    W_allcells[nonemptymgi] = 1.;
    thick_allcells[nonemptymgi] = CellThickness::THIN;
  }

  // combine the diagnostic counts over the node communicator (the ranks of each node together cover every cell
  // exactly once)
  MPI_Allreduce_safe(cells_below_mintemp, MPI_SUM, globals::mpi_comm_node);
  MPI_Allreduce_safe(cells_above_maxtemp, MPI_SUM, globals::mpi_comm_node);

  printlnlog("  cells below MINTEMP {:g} [K]: {}. Above MAXTEMP {:g} [K]: {}", MINTEMP, cells_below_mintemp, MAXTEMP,
             cells_above_maxtemp);
  MPI_Barrier_allranks();
}

// find the non-empty model grid index of the first non-empty cell at or after mgi_start, or return -1 if none found
[[nodiscard]] auto get_nonemptymgi_at_or_after(const int mgi_start) -> int {
  for (int mgi = mgi_start; mgi < get_npts_model(); mgi++) {
    if (get_numpropcells(mgi) > 0) {
      return nonemptymgi_of_mgi[mgi];
    }
  }
  return -1;
}

void setup_nstart_ndo() {
  const int nprocesses = globals::nprocs;
  assert_always(nonempty_npts_model > 0);
  const auto min_nonempty_perproc =
      nonempty_npts_model / nprocesses;  // integer division, minimum non-empty cells per process
  const auto n_remainder = nonempty_npts_model % nprocesses;

  // begin with no cell assignments
  ranks_nstart.assign(nprocesses, 0);
  ranks_nstart_nonempty.assign(nprocesses, 0);
  ranks_ndo.assign(nprocesses, 0);
  ranks_ndo_nonempty.assign(nprocesses, 0);

  if (nprocesses >= get_npts_model()) {
    // for convenience, rank == mgi when there is at least one rank per cell
    for (int rank = 0; rank < nprocesses; rank++) {
      if (rank < get_npts_model()) {
        const int mgi = rank;
        ranks_nstart[rank] = mgi;
        ranks_ndo[rank] = 1;
        ranks_nstart_nonempty[rank] = (get_numpropcells(mgi) > 0) ? get_nonemptymgi_of_mgi(mgi) : 0;
        ranks_ndo_nonempty[rank] = (get_numpropcells(mgi) > 0) ? 1 : 0;
      }
    }
  } else {
    // evenly divide up the non-empty cells among the ranks

    int rank = 0;
    for (int mgi = 0; mgi < get_npts_model(); mgi++) {
      const auto target_nonempty_thisrank = (rank < n_remainder) ? min_nonempty_perproc + 1 : min_nonempty_perproc;
      if ((rank < (nprocesses - 1)) && (ranks_ndo_nonempty[rank] >= target_nonempty_thisrank)) {
        // current rank has enough non-empty cells, so start assigning cells to the next rank
        rank++;
        ranks_nstart[rank] = mgi;
        ranks_nstart_nonempty[rank] = get_nonemptymgi_at_or_after(mgi);
        assert_always(ranks_nstart_nonempty[rank] >= 0 || rank >= get_nonempty_npts_model());
      }

      ranks_ndo[rank]++;
      if (get_numpropcells(mgi) > 0) {
        ranks_ndo_nonempty[rank]++;
      }
    }
  }

  int npts_assigned = 0;
  int nonempty_npts_model_assigned = 0;
  for (int r = 0; r < nprocesses; r++) {
    npts_assigned += ranks_ndo[r];
    nonempty_npts_model_assigned += ranks_ndo_nonempty[r];
  }
  assert_always(npts_assigned == get_npts_model());
  assert_always(nonempty_npts_model_assigned == get_nonempty_npts_model());

  if (globals::my_rank == 0) {
    auto fileout = open_output_file("modelgridrankassignments.out");
    fileout << "#rank nstart ndo ndo_nonempty\n";
    for (int r = 0; r < nprocesses; r++) {
      assert_always(ranks_ndo_nonempty[r] <= ranks_ndo[r]);
      fileout << r << ' ' << ranks_nstart[r] << ' ' << ranks_ndo[r] << ' ' << ranks_ndo_nonempty[r] << '\n';
    }
    fileout.close();
    assert_always(!fileout.fail());  // e.g. a full disk
  }
}

// get element mean weight in grams
auto get_element_meanweight(const std::ptrdiff_t nonemptymgi, const int element) -> float {
  if constexpr (USE_CALCULATED_MEANATOMICWEIGHT) {
    const auto mu = elem_meanweight_allcells[(nonemptymgi * get_nelements()) + element];
    assert_always(mu > 0);
    return mu;
  }
  return globals::elements[element].initstablemeannucmass;
}

}  // anonymous namespace

void set_rho_tmin(const int modelgridindex, const float x) { modelgrid_input[modelgridindex].rhoinit = x; }

// sum of the mass fractions of the tracked isotopes of each element in a model cell, indexed by atomic number
[[nodiscard]] auto get_tracked_isotope_massfrac_sum_by_atomic_number(const int mgi)
    -> std::array<double, max_atomic_number + 1> {
  std::array<double, max_atomic_number + 1> massfrac_sum_by_atomic_number{};
  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    const int atomic_number = decay::get_nuc_z(nucindex);
    assert_always(atomic_number <= max_atomic_number);
    massfrac_sum_by_atomic_number[atomic_number] += get_modelinitnucmassfrac(mgi, nucindex);
  }
  return massfrac_sum_by_atomic_number;
}

void calc_modelinit_totmassnuclides() {
  mtot_input = 0.;
  mfegroup = 0.;

  totmassnuclide.assign(decay::get_num_nuclides(), 0.);

  for (int mgi = 0; mgi < get_npts_model(); mgi++) {
    const double mass_in_shell = get_rho_tmin(mgi) * get_inputcellvolume(mgi);
    if (mass_in_shell > 0) {
      mtot_input += mass_in_shell;

      for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
        totmassnuclide[nucindex] += mass_in_shell * get_modelinitnucmassfrac(mgi, nucindex);
      }

      mfegroup += mass_in_shell * get_ffegrp(mgi);
    }
  }
}

// return the model cell volume (when mapped to the propagation cells) at globals::tmin
// for a uniform cubic grid this is constant
[[gnu::pure]] [[nodiscard]] auto get_modelcell_assocvolume_tmin(const int modelgridindex) -> double {
  if (get_propgridtype() == GridType::CARTESIAN3D) {
    // underlying input model could be 1D, 2D, or 3D
    return get_propcell_volume_tmin(0) * get_numpropcells(modelgridindex);
  }

  if (model_mirror_z) {
    // the two propagation cells of a mirror pair have the same volume
    return 2 * get_propcell_volume_tmin(get_upper_propcell_of_mirrored_mgi(modelgridindex));
  }

  // direct mapping, so cellindex = modelgridindex and one propagation cell per model cell
  return get_propcell_volume_tmin(modelgridindex);
}

auto get_rho_tmin(const int modelgridindex) -> float { return modelgrid_input[modelgridindex].rhoinit; }

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_rho(const std::ptrdiff_t nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return rho_allcells[nonemptymgi];
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nne(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return nne_allcells[nonemptymgi];
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_clumpfactor(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());

  if constexpr (USE_MICROCLUMPING) {
    assert_testmodeonly(std::isfinite(clumpfactor_allcells[nonemptymgi]) && clumpfactor_allcells[nonemptymgi] >= 1.F);
    return clumpfactor_allcells[nonemptymgi];
  }

  return 1.;
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nnetot(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return nnetot_allcells[nonemptymgi];
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_ffegrp(const int modelgridindex) -> float {
  return modelgrid_input[modelgridindex].ffegrp;
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_modelcell_mean_radial_pos_tmin(const int modelgridindex) -> double {
  const int assoc_cells = get_numpropcells(modelgridindex);
  return modelgrid_input[modelgridindex].initial_radial_pos_sum / assoc_cells;
}

// volume averaged mean of the squared radius [cm^2] of a model cell at tmin
[[gnu::pure]] [[nodiscard]] auto get_modelcell_mean_radial_pos_squared_tmin(const int modelgridindex) -> double {
  const int assoc_cells = get_numpropcells(modelgridindex);
  return modelgrid_input[modelgridindex].initial_radial_pos_squared_sum / assoc_cells;
}

// mass fraction of an element (all isotopes combined)
[[gnu::pure]] [[nodiscard]] auto get_elem_massfrac(const std::ptrdiff_t nonemptymgi, const int element) -> float {
  const auto massfrac = elem_massfracs_allcells[(nonemptymgi * get_nelements()) + element];
  assert_testmodeonly(massfrac >= 0.0);
  return massfrac;
}

// set the mass fraction of an element (all isotopes combined)
void set_elem_massfrac(const ptrdiff_t nonemptymgi, const int element, const float newmassfrac) {
  elem_massfracs_allcells[(nonemptymgi * get_nelements()) + element] = newmassfrac;
}

// number density [cm^-3] of an element (all isotopes combined)
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_elem_numberdens(const ptrdiff_t nonemptymgi, const int element)
    -> double {
  return get_elem_massfrac(nonemptymgi, element) / static_cast<double>(get_element_meanweight(nonemptymgi, element)) *
         get_rho(nonemptymgi);
}

void set_rho(const int nonemptymgi, const float rho) {
  assert_always(rho >= 0.);
  assert_always(std::isfinite(rho));
  rho_allcells[nonemptymgi] = rho;
}

void set_nne(const int nonemptymgi, const float nne) {
  assert_always(nne >= 0.);
  assert_always(std::isfinite(nne));
  nne_allcells[nonemptymgi] = nne;
}

void set_clumpfactor(const int nonemptymgi, const float clumpfactor) {
  assert_testmodeonly(USE_MICROCLUMPING);
  assert_always(std::isfinite(clumpfactor) && clumpfactor >= 1.F);
  clumpfactor_allcells[nonemptymgi] = clumpfactor;
}

// Calculate and set the total density of electrons (free and bound) in grid cell. These are targets for Compton
// scattering of gamma rays
void set_nnetot(const int nonemptymgi) {
  double nnetot = 0.;
  const auto nelements = get_nelements();
  for (int element = 0; element < nelements; element++) {
    nnetot += get_elem_numberdens(nonemptymgi, element) * get_atomicnumber(element);
  }

  const auto f_nnetot = static_cast<float>(nnetot);
  assert_always(f_nnetot >= 0.);
  nnetot_allcells[nonemptymgi] = f_nnetot;
}

// total number of model grid cells, including empty ones
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_npts_model() -> int {
  assert_testmodeonly(npts_model > 0);
  return static_cast<int>(npts_model);
}

// number of non-empty model grid cells, i.e. the length of the nonemptymgi index space
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nonempty_npts_model() -> int {
  assert_testmodeonly(nonempty_npts_model > 0);
  return static_cast<int>(nonempty_npts_model);
}

// get time at which model input densities are defined
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_t_model() -> double {
  assert_testmodeonly(t_model > 0.);
  return t_model;
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_propcell_modelgridindex(const int cellindex) -> int {
  assert_testmodeonly(cellindex >= 0);
  assert_testmodeonly(cellindex < ngrid);
  const auto mgi = propcell_mgi[cellindex];
  // can return -1 for empty cells, but if not empty then should be a valid modelgridindex
  assert_testmodeonly(mgi >= -1);
  assert_testmodeonly(mgi < get_npts_model());
  return mgi;
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_propcell_nonemptymgi(const int cellindex) -> int {
  const auto nonemptymgi = propcell_nonemptymgi[cellindex];
  assert_testmodeonly(nonemptymgi >= -1);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return nonemptymgi;
}

// number of propagation cells associated with each modelgrid cell
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_numpropcells(const int modelgridindex) -> int {
  assert_testmodeonly(modelgridindex < get_npts_model());
  return modelgrid_numpropcells[modelgridindex];
}

// get the index in the list of non-empty cells for a given model grid cell
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nonemptymgi_of_mgi(const int mgi) -> int {
  assert_always(get_nonempty_npts_model() > 0);
  assert_always(mgi >= 0);
  assert_always(mgi < get_npts_model());

  // the only test for an empty model cell, so an optimised build must also stop here
  const int nonemptymgi = nonemptymgi_of_mgi[mgi];
  assert_always(nonemptymgi >= 0);
  assert_always(nonemptymgi < get_nonempty_npts_model());

  return nonemptymgi;
}

// inverse of get_nonemptymgi_of_mgi(): get the model grid index of an entry in the list of non-empty cells
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_mgi_of_nonemptymgi(const ptrdiff_t nonemptymgi) -> int {
  assert_testmodeonly(get_nonempty_npts_model() > 0);
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());

  const int mgi = mgi_of_nonemptymgi[nonemptymgi];

  assert_always(mgi >= 0);
  return mgi;
}

// the abundances below are initial abundances at t_model

// get the mass fraction of a nuclide in a model grid cell at t=t_model by nuclide index
auto get_modelinitnucmassfrac(const int modelgridindex, const int nucindex) -> float {
  return initnucmassfrac_allcells[(modelgridindex * decay::get_num_nuclides()) + nucindex];
}

// get the t=t_model mass fraction of an element's untracked stable component (i.e. excluding radioactive nuclides)
auto get_elem_untrackedstable_initmassfrac(const std::ptrdiff_t nonemptymgi, const int element) -> float {
  return initmassfracuntrackedstable_allcells[(nonemptymgi * get_nelements()) + element];
}

// set element weight in grams
void set_element_meanweight(const std::ptrdiff_t nonemptymgi, const int element, const float meanweight) {
  assert_always(meanweight > 0.);
  elem_meanweight_allcells[(nonemptymgi * get_nelements()) + element] = meanweight;
}

auto get_electronfrac(const int nonemptymgi) -> double {
  double nucleondens = 0.;
  for (int element = 0; element < get_nelements(); element++) {
    nucleondens += get_elem_numberdens(nonemptymgi, element) * get_element_meanweight(nonemptymgi, element) / AMU;
  }
  return get_nnetot(nonemptymgi) / nucleondens;
}

// q: energy in the model at tmin per gram to use with the INITIAL_PACKETS_ON option [erg/g]
DEVICE_FUNC auto get_initenergyq(const int modelgridindex) -> double {
  return modelgrid_input[modelgridindex].initenergyq;
}

[[nodiscard]] auto get_elements_uppermost_ion(const std::ptrdiff_t nonemptymgi, const int element) -> int {
  const auto uppermost_ion = elements_uppermost_ion_allcells[(nonemptymgi * get_nelements()) + element];
  assert_testmodeonly(uppermost_ion >= -1);  // -1 before the first ion balance of the element in the cell
  assert_testmodeonly(uppermost_ion <= std::max(0, get_nions(element) - 1));
  return uppermost_ion;
}

void set_elements_uppermost_ion(const std::ptrdiff_t nonemptymgi, const int element, const int uppermost_ion) {
  assert_testmodeonly(uppermost_ion <= std::max(0, get_nions(element) - 1));
  elements_uppermost_ion_allcells[(nonemptymgi * get_nelements()) + element] = uppermost_ion;
}

[[nodiscard]] auto get_elements_lowermost_ion(const std::ptrdiff_t nonemptymgi, const int element) -> int {
  if constexpr (!NLTE_TRACK_SOLUTION_RANGES) {
    return 0;  // the array is not allocated when no code reads the solution ranges
  }
  const auto lowermost_ion = elements_lowermost_ion_allcells[(nonemptymgi * get_nelements()) + element];
  assert_testmodeonly(lowermost_ion >= 0);
  assert_testmodeonly(lowermost_ion <= std::max(0, get_nions(element) - 1));
  return lowermost_ion;
}

void set_elements_lowermost_ion(const std::ptrdiff_t nonemptymgi, const int element, const int lowermost_ion) {
  if constexpr (!NLTE_TRACK_SOLUTION_RANGES) {
    return;  // the array is not allocated when no code reads the solution ranges
  }
  assert_testmodeonly(lowermost_ion >= 0);
  assert_testmodeonly(lowermost_ion <= std::max(0, get_nions(element) - 1));
  elements_lowermost_ion_allcells[(nonemptymgi * get_nelements()) + element] = lowermost_ion;
}

void do_MPI_Bcast_nlte_solution_ranges(const ptrdiff_t nstart_nonempty, const ptrdiff_t ndo_nonempty,
                                       const int root_node_id) {
  const auto nelements = get_nelements();
  MPI_Bcast_safe(elements_lowermost_ion_allcells.subspan(nstart_nonempty * nelements, ndo_nonempty * nelements),
                 root_node_id, globals::mpi_comm_internode);
  MPI_Bcast_safe(elements_uppermost_ion_allcells.subspan(nstart_nonempty * nelements, ndo_nonempty * nelements),
                 root_node_id, globals::mpi_comm_internode);
}

[[nodiscard]] auto calculate_cell_kappagrey(const int nonemptymgi) -> float {
  const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
  if (get_rho_tmin(mgi) <= 0.) {
    return 0.;
  }
  double kappa = 0.;
  switch (RPKT_GREY_TYPE) {
    case RpktGreyType::FEGROUP_APPROX: {
      // kappagrey is a simple function of the initial Fe-group mass fraction
      constexpr double GREY_OP = 0.1;
      kappa = ((0.9 * get_ffegrp(mgi)) + 0.1) * GREY_OP / ((0.9 * mfegroup / mtot_input) + 0.1);
      break;
    }

    case RpktGreyType::TANAKA2020_ELECTRONFRAC: {
      // electron-fraction-dependent opacities from table 1 of Tanaka, Kato, Gaigalas & Kawaguchi (2020), MNRAS,
      // 496, 1369-1392, doi:10.1093/mnras/staa1576
      const auto Ye = modelgrid_input[mgi].initelectronfrac;
      if (Ye <= 0.) {
        fatal_crash(
            "model cell {} has no electron fraction. TANAKA2020_ELECTRONFRAC needs a Ye column in model.txt, or "
            "a Ye.txt for a model.txt of format version 1",
            mgi);
      }

      // pairs of (upper Ye limit, kappa [cm^2/g])
      constexpr auto kappa_table = std::to_array<std::pair<double, double>>(
          {{0.1, 19.5}, {0.15, 32.2}, {0.20, 22.3}, {0.25, 5.6}, {0.30, 5.36}, {0.35, 3.3}});
      const auto* const entry =
          std::ranges::find_if(kappa_table, [Ye](const auto& ye_kappa) { return Ye <= ye_kappa.first; });
      kappa = (entry != kappa_table.end()) ? entry->second : 0.96;
      break;
    }

    case RpktGreyType::JUST2022_TEMP_LANTHANIDEFRAC: {
      // grey opacity of Just, Kullmann, Goriely, Bauswein, Janka & Collins (2022), MNRAS, 510, 2820-2840,
      // doi:10.1093/mnras/stab3327. kappa is a simple analytic function of temperature and lanthanide mass fraction,
      // adapted to best fit the light curves of Kasen, Metzger, Barnes, Quataert & Ramirez-Ruiz (2017), Nature, 551,
      // 80-84, doi:10.1038/nature24453, in ALCAR simulations
      const double T_rad = TR_allcells[nonemptymgi];
      double X_lan = 0.;
      for (int element = 0; element < get_nelements(); element++) {
        const int z = get_atomicnumber(element);
        if (z >= 57 && z <= 71) {
          X_lan += get_elem_massfrac(nonemptymgi, element);
        }
      }
      // first step: temperature-independent factor
      if (X_lan < 1e-7) {
        kappa = 0.2;
      } else if (X_lan < 1e-3) {
        kappa = 3 * pow(X_lan / 1e-3, 0.3);
      } else if (X_lan < 1e-1) {
        kappa = 3 * pow(X_lan / 1e-3, 0.5);
      } else {
        kappa = 30 * pow(X_lan / 1e-1, 0.1);
      }
      // second step: multiply temperature-dependent factor
      if (T_rad < 2000.) {
        kappa *= pow5(T_rad / 2000.);
      }
      break;
    }
  }

  const auto kappa_float = static_cast<float>(kappa);
  assert_always(kappa_float >= 0.);
  assert_always(std::isfinite(kappa_float));
  return kappa_float;
}

void write_grid_restart_data(const int timestep) {
  const auto filename = std::format("gridsave_ts{}.tmp", timestep);

  const auto sys_time_start_write_restart = std::chrono::steady_clock::now();
  printlog("Write grid restart data to {}...", filename);

  FILE* gridsave_file = fopen_required(filename, "wb");

  write_restart_values(gridsave_file, globals::ntimesteps, globals::nprocs, globals::tmin, globals::tmax, timestep,
                       get_nonempty_npts_model(), get_nelements(), get_includedions(), globals::nbfcontinua_ground);
  write_restart_array(gridsave_file, mgi_of_nonemptymgi);
  for (int element = 0; element < get_nelements(); element++) {
    write_restart_values(gridsave_file, get_atomicnumber(element), get_nions(element));
  }

  write_restart_array(gridsave_file, globals::timesteps);

  // a negative or NaN deposition estimator stops the run before a restart can keep it
  assert_always(std::ranges::all_of(globals::dep_estimator_gamma, [](const double dep) { return dep >= 0.; }));
  write_restart_array(gridsave_file, TR_allcells);
  write_restart_array(gridsave_file, Te_allcells);
  write_restart_array(gridsave_file, W_allcells);
  write_restart_array(gridsave_file, TJ_allcells);
  write_restart_array(gridsave_file, thick_allcells);
  write_restart_array(gridsave_file, nne_allcells);
  write_restart_array(gridsave_file, nnetot_allcells);
  write_restart_array(gridsave_file, kpkt::radiative_energy_factor_allcells);
  write_restart_array(gridsave_file, globals::dep_estimator_gamma);
  write_restart_array(gridsave_file, globals::dep_estimator_positron);
  write_restart_array(gridsave_file, globals::dep_estimator_electron);
  write_restart_array(gridsave_file, globals::dep_estimator_alpha);
  if constexpr (USE_LUT_PHOTOION) {
    write_restart_array(gridsave_file, globals::corrphotoionrenorm);
    write_restart_array(gridsave_file, globals::gammaestimator);
    write_restart_array(gridsave_file, globals::col_iongamma_per_gspop_per_clumpednne);
  }

  // the order of these calls is very important!
  radfield::write_restart_data(gridsave_file);
  nonthermal::write_restart_data(gridsave_file);
  nltepop_write_restart_data(gridsave_file);
  if constexpr (NLTE_TRACK_SOLUTION_RANGES) {
    // the solved ion range of each element, so that a resumed run starts with the same reactions
    write_restart_array(gridsave_file, elements_lowermost_ion_allcells);
    write_restart_array(gridsave_file, elements_uppermost_ion_allcells);
  }
  // Check earlier writes and the final flush before the caller replaces the previous checkpoint.
  const bool write_failed = (ferror(gridsave_file) != 0);
  const bool close_failed = (fclose(gridsave_file) != 0);
  if (write_failed || close_failed) {
    fatal_crash("Could not write or close {}.", filename);
  }
  const auto write_restart_duration =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_write_restart).count();
  printlnlog("done in {:.1f} seconds.", write_restart_duration);
}

// get lowest modelgridindex assigned to this rank (for update_grid and output files)
auto get_nstart(const int rank) -> int {
  assert_testmodeonly(!ranks_ndo.empty());
  return ranks_nstart[rank];
}

// get lowest nonemptymgi assigned to this rank (for update_grid and output files)
auto get_nstart_nonempty(const int rank) -> int {
  assert_testmodeonly(!ranks_ndo.empty());
  return ranks_nstart_nonempty[rank];
}

// get the count of modelgridindices assigned to this rank (for update_grid and output files)
auto get_ndo(const int rank) -> int {
  assert_testmodeonly(!ranks_ndo.empty());
  return ranks_ndo[rank];
}

// get the count of nonemptymgi assigned to this rank (for update_grid and output files)
auto get_ndo_nonempty(const int rank) -> int {
  assert_testmodeonly(!ranks_ndo.empty());
  return ranks_ndo_nonempty[rank];
}

// Initialise the propagation cells and associate them with modelgrid cells
void init_grid() {
  // The cells will be ordered by x then y, then z. Call a routine that
  // sets up the initial positions and widths of the cells.
  const auto prop_gridtype = get_propgridtype();
  // a model runs on a propagation grid of its own type, or on the 3D Cartesian grid
  if (prop_gridtype != get_modelgridtype() && prop_gridtype != GridType::CARTESIAN3D) {
    fatal_crash(
        "GRID_TYPE_OVERRIDE gives a {}D {} propagation grid, but model.txt holds a {}D {} model. Give "
        "GridType::CARTESIAN3D, or give no value to use the grid type of model.txt",
        get_ndim(prop_gridtype), get_grid_type_name(prop_gridtype), get_ndim(get_modelgridtype()),
        get_grid_type_name(get_modelgridtype()));
  }
  switch (prop_gridtype) {
    case GridType::CARTESIAN3D:
      setup_grid_cartesian_3d();
      break;
    case GridType::CYLINDRICAL2D:
      setup_grid_cylindrical_2d();
      break;
    case GridType::SPHERICAL1D:
      setup_grid_spherical_1d();
      break;
    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D:
      setup_grid_spherical();
      break;
  }
  propcell_mgi.resize(ngrid, -1);

  printlnlog("propagation grid: {}-dimensional {}", get_ndim(prop_gridtype), get_grid_type_name(prop_gridtype));

  for (int d = 0; d < get_ndim(prop_gridtype); d++) {
    printlnlog("    coordinate {} '{}': cells have {} position values", d, get_coordlabel(prop_gridtype, d),
               ncoordgrid[d]);
  }
  printlnlog("    total propagation cells: {}", ngrid);

  // the mapping makes a cell with zero density empty, so it must come after this read
  auto [elem_massfracs_of_mgi, ncells_emptied] = read_elem_abundances();
  if (ncells_emptied > 0) {
    const double mtot_input_with_emptied_cells = mtot_input;
    calc_modelinit_totmassnuclides();
    printlnlog(
        "[warning] {} cells with no mass in the elements of compositiondata.txt became empty. They held {:.3e} "
        "[Msun]. The total input model mass is now {:.3e} [Msun].",
        ncells_emptied, (mtot_input_with_emptied_cells - mtot_input) / MSUN, mtot_input / MSUN);
  }

  if (get_modelgridtype() == prop_gridtype) {
    if (get_modelgridtype() == GridType::CARTESIAN3D) {
      assert_always(ncoord_model[0] == ncoordgrid[0]);
      assert_always(ncoord_model[1] == ncoordgrid[1]);
      assert_always(ncoord_model[2] == ncoordgrid[2]);
    }

    if (model_mirror_z) {
      map_mirroredmodeltogrid();
    } else {
      map_modeltogrid_direct();
    }
  } else if (get_modelgridtype() == GridType::SPHERICAL1D) {
    map_1dmodelto3dgrid();
  } else if (get_modelgridtype() == GridType::CYLINDRICAL2D) {
    map_2dmodelto3dgrid();
  } else if (has_spherical_angle_axes(get_modelgridtype())) {
    assert_always(prop_gridtype == GridType::CARTESIAN3D);
    map_sphericalmodelto3dgrid();
  } else {
    assert_always(false);
  }

  allocate_nonemptymodelcells();

  if (globals::my_rank == 0) {
    auto grid_file = open_output_file("grid.out");
    for (int cellindex = 0; cellindex < ngrid; cellindex++) {
      const int mgi = get_propcell_modelgridindex(cellindex);
      if (mgi >= 0) {
        std::println(grid_file, "{} {}", cellindex, mgi);  // write only non-empty cells to grid file
      }
    }
    grid_file.close();
    assert_always(!grid_file.fail());  // e.g. a full disk
  }

  // the rank ranges need the count of the non-empty cells
  setup_nstart_ndo();

  // when the model grid and the propagation grid differ, rescale the nuclide mass fractions so that the total mass
  // of each nuclide matches the input model again. Every propagation cell takes the model shell that its centre falls
  // in, so the volume associated with a shell is a staircase approximation of the true shell volume and the
  // mapped mass of each nuclide differs from the input.
  //
  // This has to happen before anything derives a quantity from the nuclide mass fractions, which is why it sits
  // ahead of set_untrackedstable_massfracs() (which subtracts them from the elemental mass fractions to get the
  // untracked stable remainder) and assign_initial_temperatures() (which sums their decay energy). It has to
  // come after read_elem_abundances(), which checks the input files against each other while the mass fractions
  // are still the ones they gave.
  //
  // Gated on the two geometries differing rather than on a particular pair of them, because every such mapping
  // loses mass the same way: map_1dmodelto3dgrid() and map_2dmodelto3dgrid() both place a propagation cell by
  // its centre. Only GRID_TYPE_OVERRIDE makes them differ; left alone, get_propgridtype() returns the model's
  // own geometry and map_modeltogrid_direct() gives each model cell exactly its own volume, so there is no
  // discretisation error to correct and the ratio would be one throughout. The "Total grid-mapped mass" line
  // logged below shows the size of the discrepancy for any model.
  //
  // The mapped total includes the mass in the propagation cells that FORCE_SPHERICAL_ESCAPE_SURFACE
  // blanked. allocate_nonemptymodelcells() records that mass in totmassnuclide_blanked. The option
  // removes that mass on purpose, so the ratio must not put it back into the surviving cells. A
  // directly mapped model under the same option also drops the blanked mass, so the two cases agree.
  const bool mapping_loses_nuclide_mass = (get_modelgridtype() != prop_gridtype);
  if (mapping_loses_nuclide_mass && globals::rank_in_node == 0) {
    for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
      if (totmassnuclide[nucindex] <= 0) {
        continue;
      }

      double totmassnuclide_actual = totmassnuclide_blanked[nucindex];
      for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
        const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
        totmassnuclide_actual +=
            get_modelinitnucmassfrac(mgi, nucindex) * get_rho_tmin(mgi) * get_modelcell_assocvolume_tmin(mgi);
      }

      if (totmassnuclide_actual > 0.) {
        const double ratio = totmassnuclide[nucindex] / totmassnuclide_actual;
        for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
          const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
          const double new_massfrac = get_modelinitnucmassfrac(mgi, nucindex) * ratio;
          // a mass fraction above one cannot keep the input mass of the nuclide, so the grid is too coarse. An
          // excess at the rounding level of a float is clamped.
          if (new_massfrac > 1. + 1e-6) {
            fatal_crash(
                "init_grid: the mapping to the propagation grid needs a mass fraction of {:g} for Z={} A={} in cell "
                "{} to keep the input mass of the nuclide. Use a finer grid.",
                new_massfrac, decay::get_nuc_z(nucindex), decay::get_nuc_a(nucindex), mgi);
          }
          set_modelinitnucmassfrac(mgi, nucindex, static_cast<float>(std::min(new_massfrac, 1.)));
        }
      } else {
        fatal_crash(
            "init_grid: the mapping to the propagation grid gives no cell with Z={} A={}, which has mass in the "
            "model. Use a finer grid.",
            decay::get_nuc_z(nucindex), decay::get_nuc_a(nucindex));
      }
    }
  }

  // the untracked stable remainder is whatever the rescaled isotopes leave of each element, so it can only be
  // worked out once the rescale above is done
  set_untrackedstable_massfracs(elem_massfracs_of_mgi);
  elem_massfracs_of_mgi = std::vector<float>{};  // free the memory before the later allocations
  MPI_Barrier_node();

  radfield::init();
  nonthermal::init();

  // and assign a temperature to the cells
  if (globals::simulation_continued_from_saved) {
    // For continuation of an existing simulation we read the temperatures
    // at the end of the simulation and write them to the grid.
    read_grid_restart_data(globals::timestep_initial);
  } else {
    assign_initial_temperatures();
  }

  // this call also fills the cached value before the packet propagation threads read it
  const double mtot_mapped = get_ejecta_mass();
  printlnlog("Total grid-mapped mass: {:9.3e} [Msun] ({:.1f}% of input mass)", mtot_mapped / MSUN,
             mtot_mapped / mtot_input * 100.);

  MPI_Barrier_allranks();

  if (globals::rank_in_node == 0) {
    for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
      kappagrey_allcells[nonemptymgi] = calculate_cell_kappagrey(nonemptymgi);
    }
  }

  MPI_Barrier_node();
}

auto get_totmassnuclide_tmodel(const int z, const int a) -> double { return totmassnuclide[decay::get_nucindex(z, a)]; }

}  // namespace grid
