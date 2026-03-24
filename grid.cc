#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "kpkt.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "radfield.h"
#include "rpkt.h"
#include "sn3d.h"
#include "vectors.h"

namespace grid {

namespace {

std::array<int, 3> ncoordgrid{0, 0, 0};  // propagation grid dimensions

std::optional<GridType> model_type{};

ptrdiff_t npts_model = 0;  // number of model grid cells
ptrdiff_t nonempty_npts_model = 0;  // number of allocated non-empty model grid cells

double t_model = -1.;  // time at which densities in input model are correct.
std::vector<double> vout_model{};
std::array<int, 3> ncoord_model{};  // the model.txt input grid dimensions

double min_den{-1.};  // minimum model density

double mfegroup{0.};  // Total mass of Fe group elements in ejecta

int first_cellindex{-1};  // auto-determine first cell index in model.txt (usually 1 or 0)

// Initial co-ordinates of inner most corner of cell.
std::vector<Vec3d> propcell_pos_min{};

// associate each propagation cell with a model grid cell, or not, if the cell is empty (or doesn't get mapped to
// anything such as 1D/2D to 3D)
std::vector<int> propcell_mgi;
std::vector<int> propcell_nonemptymgi;

std::vector<int> modelgrid_numpropcells;
std::vector<int> nonemptymgi_of_mgi;
std::vector<int> mgi_of_nonemptymgi;

std::vector<double> totmassnuclide{};  // total mass of each nuclide in the ejecta

MPI_Win win_nltepops_allcells = MPI_WIN_NULL;
MPI_Win win_initnucmassfrac_allcells = MPI_WIN_NULL;

std::span<float> initnucmassfrac_allcells{};
std::span<float> initmassfracuntrackedstable_allcells{};
std::span<int> elements_uppermost_ion_allcells{};  // Highest ion index that has a significant population

// indexed by global rank
std::vector<int> ranks_nstart;
std::vector<int> ranks_nstart_nonempty;
std::vector<int> ranks_ndo;
std::vector<int> ranks_ndo_nonempty;

struct ModelGridCellInput {
  float rhoinit = -1.;
  float ffegrp = 0.;
  float initial_radial_pos_sum = 0.;
  float initelectronfrac = 0.4;  // Ye: electrons (or protons) per nucleon
  float initenergyq = 0.;  // q: energy in the model at tmin to use with USE_MODEL_INITIAL_ENERGY [erg/g]
};
std::span<ModelGridCellInput> modelgrid_input{};

enum class BoundaryType : std::uint8_t { INNER, OUTER };

// Get number of dimensions
constexpr auto get_ndim(const GridType gridtype) -> int {
  switch (gridtype) {
    case GridType::SPHERICAL1D:
      return 1;
    case GridType::CYLINDRICAL2D:
      return 2;
    case GridType::CARTESIAN3D:
      return 3;
    default:
      assert_always(false);
      return -1;
  }
}

[[nodiscard]] DEVICE_FUNC constexpr auto get_coordlabel(const GridType gridtype, const int axis) -> char {
  assert_always(axis >= 0 && axis < get_ndim(gridtype));
  switch (gridtype) {
    case GridType::CARTESIAN3D:
      switch (axis) {
        case 0:
          return 'x';
        case 1:
          return 'y';
        case 2:
          return 'z';
        default:
          return '?';
      }
    case GridType::CYLINDRICAL2D:
      switch (axis) {
        case 0:
          return 'r';
        case 1:
          return 'z';
        default:
          return '?';
      }
    case GridType::SPHERICAL1D:
      return 'r';
    default:
      assert_always(false);
      return '?';
  }
}

void set_rho_tmin(const int modelgridindex, const float x) { modelgrid_input[modelgridindex].rhoinit = x; }

void set_initelectronfrac(const int modelgridindex, const float electronfrac) {
  modelgrid_input[modelgridindex].initelectronfrac = electronfrac;
}

void read_possible_yefile() {
  if (!std::filesystem::exists("Ye.txt")) {
    printlnlog("Ye.txt not found");
    return;
  }

  auto filein = fopen_required_uniqueptr("Ye.txt", "r");
  int nlines_in = 0;
  assert_always(fscanf(filein.get(), "%d", &nlines_in) == 1);

  for (int n = 0; n < nlines_in; n++) {
    int mgiplusone = -1;
    float initelecfrac = 0.;
    assert_always(fscanf(filein.get(), "%d %g", &mgiplusone, &initelecfrac) == 2);
    const int mgi = mgiplusone - 1;
    if (mgi >= 0 && mgi < get_npts_model()) {
      set_initelectronfrac(mgi, initelecfrac);
    }
  }
}

// return the inner radius (or equivalent) of a propagation cell at time tmin
auto get_cell_r_inner(const int cellindex, const GridType prop_gridtype) -> double {
  if (prop_gridtype == GridType::SPHERICAL1D) {
    return get_cellcoordmin(cellindex, 0);
  }

  if (prop_gridtype == GridType::CYLINDRICAL2D) {
    const auto rcyl_inner = get_cellcoordmin(cellindex, 0);
    const auto z_inner = std::min(std::abs(get_cellcoordmin(cellindex, 1)), std::abs(get_cellcoordmax(cellindex, 1)));
    return std::sqrt(std::pow(rcyl_inner, 2) + std::pow(z_inner, 2));
  }

  if (prop_gridtype == GridType::CARTESIAN3D) {
    const auto x_inner = std::min(std::abs(get_cellcoordmin(cellindex, 0)), std::abs(get_cellcoordmax(cellindex, 0)));
    const auto y_inner = std::min(std::abs(get_cellcoordmin(cellindex, 1)), std::abs(get_cellcoordmax(cellindex, 1)));
    const auto z_inner = std::min(std::abs(get_cellcoordmin(cellindex, 2)), std::abs(get_cellcoordmax(cellindex, 2)));
    return std::sqrt(std::pow(x_inner, 2) + std::pow(y_inner, 2) + std::pow(z_inner, 2));
  }

  assert_always(false);
  return NAN;
}

void set_ffegrp(const int modelgridindex, float x) {
  if (!(x >= 0.)) {
    printlnlog("WARNING: Fe-group mass fraction {:g} is negative in cell {}", x, modelgridindex);
    assert_always(x > -1e-6);
    x = 0.;
  }

  assert_always(x >= 0);
  assert_always(x <= 1.001);
  modelgrid_input[modelgridindex].ffegrp = x;
}

void set_propcell_modelgridindex(const int cellindex, const int new_modelgridindex) {
  assert_testmodeonly(cellindex >= 0);
  assert_testmodeonly(cellindex < ngrid);
  assert_testmodeonly(new_modelgridindex >= -1);
  assert_testmodeonly(new_modelgridindex < get_npts_model());
  propcell_mgi[cellindex] = new_modelgridindex;
}

void set_modelinitnucmassfrac(const int modelgridindex, const int nucindex, float abund) {
  // set the mass fraction of a nuclide in a model grid cell at t=t_model by nuclide index
  // initnucmassfrac array is in node shared memory
  assert_always(nucindex >= 0);
  if (!(abund >= 0.)) {
    printlnlog("WARNING: nuclear mass fraction for nucindex {} = {:g} is negative in cell {}", nucindex, abund,
               modelgridindex);
    assert_always(abund > -1e-6);
    abund = 0.;
  }

  assert_always(abund >= 0.);
  assert_always(abund <= 1.);
  const ptrdiff_t num_nuclides = decay::get_num_nuclides();

  initnucmassfrac_allcells[(modelgridindex * num_nuclides) + nucindex] = abund;
}

void set_initenergyq(const int modelgridindex, const float initenergyq) {
  modelgrid_input[modelgridindex].initenergyq = initenergyq;
}

void set_elem_untrackedstable_abund_from_total(const int nonemptymgi, const int element, const float elemabundance) {
  // set the stable mass fraction of an element from the total element mass fraction
  // by subtracting the abundances of radioactive isotopes.
  // if the element Z=anumber has no specific stable abundance variable then the function does nothing

  const int atomic_number = get_atomicnumber(element);
  const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);

  double isofracsum = 0.;  // mass fraction sum of radioactive isotopes
  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    if (decay::get_nuc_z(nucindex) == atomic_number) {
      // radioactive isotope of this element
      isofracsum += get_modelinitnucmassfrac(mgi, nucindex);
    }
  }

  double massfrac_untrackedstable = elemabundance - isofracsum;

  if (massfrac_untrackedstable < 0.) {
    //  allow some roundoff error before we complain
    if ((isofracsum - elemabundance - 1.) > 1e-4 && std::abs(isofracsum - elemabundance) > 1e-6) {
      printlnlog("WARNING: cell {} Z={} element abundance is less than the sum of its radioisotope abundances", mgi,
                 atomic_number);
      printlnlog("  massfrac(Z) {:g} massfrac_radioisotopes(Z) {:g}", elemabundance, isofracsum);
      printlnlog("  increasing elemental abundance to {:g} and setting stable isotopic abundance to zero", isofracsum);
    }
    // result is allowed to be slightly negative due to roundoff error
    assert_always(massfrac_untrackedstable >= -1e-2);
    massfrac_untrackedstable = 0.;  // bring up to zero if negative
  }

  // if (globals::rank_in_node == 0)
  {
    initmassfracuntrackedstable_allcells[(nonemptymgi * get_nelements()) + element] =
        static_cast<float>(massfrac_untrackedstable);
  }

  // (isofracsum + massfracstable) might not exactly match elemabundance if we had to boost it to reach isofracsum
  set_elem_abundance(nonemptymgi, element, static_cast<float>(isofracsum + massfrac_untrackedstable));
}

void allocate_nonemptycells_composition_cooling() {
  // Initialise composition dependent cell data for the given cell
  const ptrdiff_t nonempty_npts_model_ptrdifft = get_nonempty_npts_model();
  const auto nelements = get_nelements();

  initmassfracuntrackedstable_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model_ptrdifft * nelements, 0.);
  elem_meanweight_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model_ptrdifft * nelements, 0.);
  elements_uppermost_ion_allcells = MPI_shared_malloc_span<int>(nonempty_npts_model_ptrdifft * nelements, -1);
  elem_massfracs_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model_ptrdifft * nelements, 0.);
  ion_groundlevelpops_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model_ptrdifft * get_includedions(), 0.);
  ion_partfuncts_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model_ptrdifft * get_includedions(), 0.);
  kpkt::ion_cooling_contribs_allcells =
      MPI_shared_malloc_span<double>(nonempty_npts_model_ptrdifft * get_includedions(), 0.);

  // -1 indicates that there is currently no information on the nlte populations
  std::tie(nltepops_allcells, win_nltepops_allcells) =
      MPI_shared_malloc_span_keepwin<double>(nonempty_npts_model_ptrdifft * globals::total_nlte_levels, -1.);
}

void allocate_nonemptymodelcells() {
  // Determine the number of simulation cells associated with the model cells
  std::ranges::fill(modelgrid_numpropcells, 0);
  if (globals::rank_in_node == 0) {
    for (int mgi = 0; mgi < get_npts_model(); mgi++) {
      modelgrid_input[mgi].initial_radial_pos_sum = 0.;
    }
  }
  MPI_Barrier(globals::mpi_comm_node);

  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const auto radial_pos_mid = get_cellradialposmid(cellindex);

    if (FORCE_SPHERICAL_ESCAPE_SURFACE && radial_pos_mid > globals::vmax * globals::tmin) {
      // for 1D models, the final shell outer v should already be at vmax
      assert_always(model_type != GridType::SPHERICAL1D || propcell_mgi[cellindex] < 0);
      set_propcell_modelgridindex(cellindex, -1);
    }

    const int mgi = get_propcell_modelgridindex(cellindex);
    if (mgi >= 0) {
      modelgrid_numpropcells[mgi] += 1;
      if (globals::rank_in_node == 0) {
        modelgrid_input[mgi].initial_radial_pos_sum =
            static_cast<float>(modelgrid_input[mgi].initial_radial_pos_sum + radial_pos_mid);
      }
      assert_always(get_rho_tmin(mgi) > 0);
      // with direct mapping, there must be exactly one propagation cell per non-empty model cell
      assert_always(get_modelgridtype() != get_propgridtype() || modelgrid_numpropcells[mgi] == 1);
    }
  }

  MPI_Barrier(globals::mpi_comm_node);
  // find number of non-empty cells and allocate nonempty list
  nonempty_npts_model = 0;
  for (int mgi = 0; mgi < get_npts_model(); mgi++) {
    if (get_numpropcells(mgi) > 0) {
      nonempty_npts_model++;
    }
  }
  assert_always(nonempty_npts_model > 0);
  printlnlog("There are {} modelgrid cells with associated propagation cells (nonempty_npts_model)",
             nonempty_npts_model);

  resize_exactly(mgi_of_nonemptymgi, nonempty_npts_model);
  std::ranges::fill(mgi_of_nonemptymgi, -2);

  resize_exactly(propcell_nonemptymgi, ngrid);
  std::ranges::fill(propcell_nonemptymgi, -1);

  int nonemptymgi = 0;  // index within list of non-empty modelgrid cells

  for (int mgi = 0; mgi < get_npts_model(); mgi++) {
    if (get_numpropcells(mgi) > 0) {
      assert_always(get_rho_tmin(mgi) >= 0);
      nonemptymgi_of_mgi[mgi] = nonemptymgi;
      mgi_of_nonemptymgi[nonemptymgi] = mgi;
      nonemptymgi++;
    } else {
      nonemptymgi_of_mgi[mgi] = -1;
      set_rho_tmin(mgi, 0.);
      for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
        set_modelinitnucmassfrac(mgi, nucindex, 0.);
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const int mgi = get_propcell_modelgridindex(cellindex);
    if (mgi < 0) {
      propcell_nonemptymgi[cellindex] = -1;
    } else {
      propcell_nonemptymgi[cellindex] = get_nonemptymgi_of_mgi(mgi);
    }
  }

  assert_always(rho_allcells.empty());
  rho_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model, -1.);
  Te_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model, -1.);
  TJ_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model, -1.);
  TR_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model, -1.);
  W_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model, -1.);
  nne_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model, -1.);
  nnetot_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model, -1.);
  kappagrey_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model, 0.);
  grey_depth_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model, 0.);
  thick_allcells = MPI_shared_malloc_span<int>(nonempty_npts_model, 0);
  const auto modelgrid_mem_usage = nonempty_npts_model * ((sizeof(float) * 9) + sizeof(double) + sizeof(int));
  printlnlog(
      "[info] mem_usage: the modelgrid properties (temperatures and electron densities) occupies {:.3f} MB (node "
      "shared memory)",
      modelgrid_mem_usage / 1024. / 1024.);

  allocate_nonemptycells_composition_cooling();

  if constexpr (EXPANSIONOPACITIES_ON || RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value()) {
    allocate_expansionopacities();
  }

  resize_exactly(globals::dep_estimator_gamma, nonempty_npts_model);
  std::ranges::fill(globals::dep_estimator_gamma, 0.);

  resize_exactly(globals::dep_estimator_positron, nonempty_npts_model);
  std::ranges::fill(globals::dep_estimator_positron, 0.);

  resize_exactly(globals::dep_estimator_electron, nonempty_npts_model);
  std::ranges::fill(globals::dep_estimator_electron, 0.);

  resize_exactly(globals::dep_estimator_alpha, nonempty_npts_model);
  std::ranges::fill(globals::dep_estimator_alpha, 0.);

  const auto ionestimcount = nonempty_npts_model * globals::nbfcontinua_ground;
  const auto ionestimsize = ionestimcount * sizeof(double);

  if (ionestimsize > 0) {
    std::tie(globals::corrphotoionrenorm, globals::win_corrphotoionrenorm) =
        MPI_shared_malloc_span_keepwin<double>(ionestimcount, 1.);

    resize_exactly(globals::gammaestimator, ionestimcount);
    std::ranges::fill(globals::gammaestimator, 0.);
#ifdef DO_TITER
    resize_exactly(globals::gammaestimator_save, ionestimcount);
    std::ranges::fill(globals::gammaestimator_save, 0.);
#endif
  } else {
    globals::corrphotoionrenorm = {};
    globals::gammaestimator.clear();
#ifdef DO_TITER
    globals::gammaestimator_save.clear();
#endif
  }

  if (USE_ION_BFHEATING_ESTIMATORS && ionestimsize > 0) {
    resize_exactly(globals::bfheatingestimator, ionestimcount);
    std::ranges::fill(globals::bfheatingestimator, 0.);
#ifdef DO_TITER
    resize_exactly(globals::bfheatingestimator_save, ionestimcount);
    std::ranges::fill(globals::bfheatingestimator_save, 0.);
#endif
  } else {
    globals::bfheatingestimator.clear();
#ifdef DO_TITER
    globals::bfheatingestimator_save.clear();
#endif
  }

  resize_exactly(globals::ffheatingestimator, nonempty_npts_model);
  std::ranges::fill(globals::ffheatingestimator, 0.);

  resize_exactly(globals::colheatingestimator, DIRECT_COL_HEAT ? 0 : nonempty_npts_model);
  std::ranges::fill(globals::colheatingestimator, 0.);

#ifdef DO_TITER
  resize_exactly(globals::ffheatingestimator_save, nonempty_npts_model);
  std::ranges::fill(globals::ffheatingestimator_save, 0.);

  resize_exactly(globals::colheatingestimator_save, DIRECT_COL_HEAT ? 0 : nonempty_npts_model);
  std::ranges::fill(globals::colheatingestimator_save, 0.);
#endif

  // barrier to make sure node master has set abundance values to node shared memory
  MPI_Barrier(MPI_COMM_WORLD);

  printlnlog(
      "[info] mem_usage: NLTE populations for all allocated cells occupy a total of {:.3f} MB (node shared memory)",
      get_nonempty_npts_model() * globals::total_nlte_levels * sizeof(double) / 1024. / 1024.);
}

// Map 1D spherical model grid onto 3D Cartesian propagation grid
void map_1dmodelto3dgrid() {
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const double cellvmid = get_cellradialposmid(cellindex) / globals::tmin;
    const int mgi = static_cast<int>(std::ranges::lower_bound(vout_model, cellvmid) - vout_model.begin());

    if (mgi < get_npts_model() && modelgrid_input[mgi].rhoinit > 0) {
      set_propcell_modelgridindex(cellindex, mgi);
      assert_always(vout_model[mgi] >= cellvmid);
      assert_always((mgi > 0 ? vout_model[mgi - 1] : 0.0) <= cellvmid);
    } else {
      // corner cells outside of the outermost model shell are empty
      // and so are any shells with zero density
      set_propcell_modelgridindex(cellindex, -1);
    }
  }
}

// Map 2D cylindrical model onto 3D Cartesian propagation grid
void map_2dmodelto3dgrid() {
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    // map to 3D Cartesian grid
    const auto pos_mid = Vec3d{get_cellcoordmin(cellindex, 0) + (0.5 * propcell_width_tmin(cellindex, 0)),
                               get_cellcoordmin(cellindex, 1) + (0.5 * propcell_width_tmin(cellindex, 1)),
                               get_cellcoordmin(cellindex, 2) + (0.5 * propcell_width_tmin(cellindex, 2))};

    const double rcylindrical = std::sqrt(std::pow(pos_mid[0], 2) + std::pow(pos_mid[1], 2));

    // 2D grid is uniform so rcyl and z indices can be calculated with no lookup
    const int n_rcyl = static_cast<int>(rcylindrical / globals::tmin / globals::vmax * ncoord_model[0]);
    const int n_z =
        static_cast<int>(((pos_mid[2] / globals::tmin) + globals::vmax) / (2 * globals::vmax) * ncoord_model[1]);

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

// here mgi and cellindex are interchangeable (except for empty cells that are associated with mgi == -1)
void map_modeltogrid_direct() {
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const int mgi = (modelgrid_input[cellindex].rhoinit > 0) ? cellindex : -1;
    set_propcell_modelgridindex(cellindex, mgi);
  }
}

void read_abundances() {
  // barrier to make sure node master has set values in node shared memory
  MPI_Barrier(MPI_COMM_WORLD);
  printlog("reading abundances.txt...");
  const bool threedimensional = (get_modelgridtype() == GridType::CARTESIAN3D);

  // Open the abundances file
  auto abundance_file = fstream_required("abundances.txt", std::ios::in);

  // and process through the grid to read in the abundances per cell
  // The abundance file should only contain information for non-empty
  // cells. Its format must be cellnumber (integer), abundance for
  // element Z=1 (float) up to abundance for element Z=30 (float)
  // i.e. in total one integer and 30 floats.

  static std::string line;
  static std::istringstream ssline;

  // loop over propagation cells for 3D models, or modelgrid cells
  for (int mgi = 0; mgi < get_npts_model(); mgi++) {
    assert_always(get_noncommentline(abundance_file, line));
    ssline.clear();
    ssline.str(line);

    int cellnumberinput = -1;
    assert_always(ssline >> cellnumberinput);
    assert_always(cellnumberinput == mgi + first_cellindex);

    // the abundances.txt file specifies the elemental mass fractions for each model cell
    // (or proportial to mass frac, e.g. element densities because they will be normalised anyway)
    // The abundances begin with hydrogen, helium, etc, going as far up the atomic numbers as required
    double normfactor = 0.;
    std::array<float, 150> elem_abundances_in{};
    std::ranges::fill(elem_abundances_in, 0.);
    double abund_in = 0.;
    for (int elem_z_index = 0; elem_z_index < std::ssize(elem_abundances_in); elem_z_index++) {
      const int atomic_number = elem_z_index + 1;
      if (!(ssline >> abund_in)) {
        // at least one element (hydrogen) should have been specified for nonempty cells
        assert_always(atomic_number > 1 || get_numpropcells(mgi) == 0);
        break;
      }

      if (abund_in < 0. || abund_in < std::numeric_limits<float>::min()) {
        assert_always(abund_in > -1e-6);
        abund_in = 0.;
      }
      elem_abundances_in[elem_z_index] = static_cast<float>(abund_in);
      normfactor += elem_abundances_in[elem_z_index];
    }

    if (get_numpropcells(mgi) > 0) {
      if (threedimensional || normfactor <= 0.) {
        normfactor = 1.;
      }
      const int nonemptymgi = get_nonemptymgi_of_mgi(mgi);

      for (int element = 0; element < get_nelements(); element++) {
        // now set the abundances (by mass) of included elements, i.e.
        // read out the abundances specified in the atomic data file
        const int atomic_number = get_atomicnumber(element);
        const auto elemabundance = static_cast<float>(elem_abundances_in[atomic_number - 1] / normfactor);
        assert_always(elemabundance >= 0.);

        // radioactive nuclide abundances should have already been set by read_??_model
        set_elem_untrackedstable_abund_from_total(nonemptymgi, element, elemabundance);
      }
    }
  }

  // barrier to make sure node master has set values in node shared memory
  MPI_Barrier(MPI_COMM_WORLD);
  printlnlog("done.");
}

void parse_model_headerline(const std::string& line, std::vector<int>& zlist, std::vector<int>& alist,
                            std::vector<std::string>& colnames) {
  // custom header line
  std::istringstream iss(line);
  std::string token;

  int columnindex = -1;

  while (std::getline(iss, token, ' ')) {
    if (std::ranges::all_of(token, isspace)) {  // skip whitespace tokens
      continue;
    }

    columnindex++;

    if (token == "#inputcellid") {
      assert_always(columnindex == 0);
    } else if (token == "velocity_outer") {
      assert_always(columnindex == 1);
    } else if (token == "vel_r_max_kmps") {
      assert_always(columnindex == 1);
    } else if (token.starts_with("pos_")) {
      continue;
    } else if (token == "logrho") {
      // 1D models have log10(rho [g/cm3])
      assert_always(columnindex == 2);
      assert_always(get_modelgridtype() == GridType::SPHERICAL1D);
    } else if (token == "rho") {
      // 2D and 3D models have rho [g/cm3]
      assert_always(get_modelgridtype() != GridType::SPHERICAL1D);
      assert_always((columnindex == 4 && get_modelgridtype() == GridType::CARTESIAN3D) ||
                    (columnindex == 3 && get_modelgridtype() == GridType::CYLINDRICAL2D));
      continue;
    } else if (token.starts_with("X_") && token != "X_Fegroup") {
      colnames.push_back(token);
      const int z = decay::get_nucstring_z(token.substr(2));  // + 2 skips the 'X_'
      const int a = decay::get_nucstring_a(token.substr(2));
      assert_always(z >= 0);
      assert_always(a >= 0);
      zlist.push_back(z);
      alist.push_back(a);
    } else {
      colnames.push_back(token);
      zlist.push_back(-1);
      alist.push_back(-1);
    }
  }
}

auto get_token_count(std::string& line) -> int {
  std::string token;
  int abundcolcount = 0;
  auto ssline = std::istringstream{line};
  while (std::getline(ssline, token, ' ')) {
    if (!std::ranges::all_of(token, isspace)) {  // skip whitespace tokens
      abundcolcount++;
    }
  }
  return abundcolcount;
}

void read_model_radioabundances(std::istream& fmodel, std::istringstream& ssline, const int mgi, const bool keepcell,
                                const std::vector<std::string>& colnames, const std::vector<int>& nucindexlist,
                                const bool one_line_per_cell) {
  if (!one_line_per_cell) {
    static std::string line;
    assert_always(std::getline(fmodel, line));
    ssline.clear();
    ssline.str(line);
  }

  if (!keepcell) {
    return;
  }

  for (auto i = 0Z; i < std::ssize(colnames); i++) {
    double valuein = 0.;
    assert_always(ssline >> valuein);  // usually a mass fraction, but now can be anything

    if (nucindexlist[i] >= 0) {
      assert_testmodeonly(valuein <= 1.);
      set_modelinitnucmassfrac(mgi, nucindexlist[i], static_cast<float>(valuein));
    } else if (colnames[i] == "X_Fegroup") {
      set_ffegrp(mgi, static_cast<float>(valuein));
    } else if (colnames[i] == "cellYe" || colnames[i] == "Ye") {
      set_initelectronfrac(mgi, static_cast<float>(valuein));
    } else if (colnames[i] == "q") {
      // use value for t_model and adjust to tmin with expansion factor
      set_initenergyq(mgi, static_cast<float>(valuein * t_model / globals::tmin));
    } else if (colnames[i] == "tracercount") {
      ;
    } else {
      if (mgi == 0) {
        printlnlog("WARNING: ignoring column '{}' nucindex {} valuein[mgi=0] {:g}", colnames[i], nucindexlist[i],
                   valuein);
      }
    }
  }
  double valuein = 0.;
  assert_always(!(ssline >> valuein));  // should be no tokens left!
}

auto read_model_columns(std::istream& fmodel) -> std::tuple<std::vector<std::string>, std::vector<int>, bool> {
  auto pos_data_start = fmodel.tellg();  // get position in case we need to undo getline

  std::vector<int> zlist;
  std::vector<int> alist;
  std::vector<std::string> colnames;

  std::string line;
  std::getline(fmodel, line);

  std::string headerline;

  const bool header_specified = lineiscommentonly(line);

  if (header_specified) {
    // line is the header
    headerline = line;
    pos_data_start = fmodel.tellg();
    std::getline(fmodel, line);
  } else {
    // line is not a comment, so it must be the first line of data
    // add a default header for unlabelled columns
    switch (get_modelgridtype()) {
      case GridType::SPHERICAL1D:
        headerline = std::string("#inputcellid vel_r_max_kmps logrho");
        break;
      case GridType::CYLINDRICAL2D:
        headerline = std::string("#inputcellid pos_rcyl_mid pos_z_mid rho");
        break;
      case GridType::CARTESIAN3D:
        headerline = std::string("#inputcellid pos_x_min pos_y_min pos_z_min rho");
        break;
    }
    headerline += std::string(" X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48");
  }

  int colcount = get_token_count(line);
  const bool one_line_per_cell = (colcount >= get_token_count(headerline));

  printlnlog("model.txt has {} line per cell format", one_line_per_cell ? "one" : "two");

  if (!one_line_per_cell) {  // add columns from the second line
    std::getline(fmodel, line);
    colcount += get_token_count(line);
  }

  if (!header_specified && colcount > get_token_count(headerline)) {
    headerline += " X_Ni57 X_Co57";
  }

  assert_always(colcount == get_token_count(headerline));

  fmodel.seekg(pos_data_start);  // get back to start of data

  if (header_specified) {
    printlnlog("model.txt has header line: {}", headerline);
  } else {
    printlnlog("model.txt has no header line. Using default: {}", headerline);
  }

  parse_model_headerline(headerline, zlist, alist, colnames);

  decay::init_nuclides(zlist, alist);

  std::vector<int> nucindexlist(zlist.size());
  for (auto i = 0Z; i < std::ssize(zlist); i++) {
    nucindexlist[i] = (zlist[i] > 0) ? decay::get_nucindex(zlist[i], alist[i]) : -1;
  }

  assert_always(npts_model > 0);

  const ptrdiff_t num_nuclides = decay::get_num_nuclides();

  std::tie(initnucmassfrac_allcells, win_initnucmassfrac_allcells) =
      MPI_shared_malloc_span_keepwin<float>((npts_model + 1) * num_nuclides, 0.);
  printlnlog(
      "[info] mem_usage: input abundance data for {} nuclides for {} cells occupies {:.3f} MB (node shared memory)",
      num_nuclides, npts_model, (initnucmassfrac_allcells.size() * sizeof(float)) / 1024. / 1024.);

  return {colnames, nucindexlist, one_line_per_cell};
}

auto get_inputcellvolume(const int mgi) -> double {
  switch (get_modelgridtype()) {
    case GridType::SPHERICAL1D: {
      const double v_inner = (mgi == 0) ? 0. : vout_model[mgi - 1];
      return (pow(vout_model[mgi], 3) - pow(v_inner, 3)) * 4 * PI * pow(globals::tmin, 3) / 3.;
    }

    case GridType::CYLINDRICAL2D: {
      const int n_r = mgi % ncoord_model[0];
      const double delta_rcyl = globals::vmax * t_model / ncoord_model[0];
      const double delta_z = 2. * globals::vmax * t_model / ncoord_model[1];
      return pow(globals::tmin / t_model, 3) * delta_z * PI *
             (pow((n_r + 1) * delta_rcyl, 2) - pow(n_r * delta_rcyl, 2));
    }

    case GridType::CARTESIAN3D: {
      // Assumes cells are cubes here - all same volume.
      return pow((2 * globals::vmax * globals::tmin), 3) / (ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2]);
    }
  }

  assert_always(false);
  return NAN;
}

void calc_modelinit_totmassnuclides() {
  mtot_input = 0.;
  mfegroup = 0.;

  resize_exactly(totmassnuclide, decay::get_num_nuclides());
  std::ranges::fill(totmassnuclide, 0.);

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

void read_grid_restart_data(const int timestep) {
  const auto filename = std::format("gridsave_ts{}.tmp", timestep);

  printlnlog("READIN GRID SNAPSHOT from {}", filename);
  FILE* gridsave_file = fopen_required(filename, "r");

  int ntimesteps_in = -1;
  assert_always(fscanf(gridsave_file, "%d ", &ntimesteps_in) == 1);
  assert_always(ntimesteps_in == globals::ntimesteps);

  int nprocs_in = -1;
  assert_always(fscanf(gridsave_file, "%d ", &nprocs_in) == 1);
  assert_always(nprocs_in == globals::nprocs);

  for (int nts = 0; nts < globals::ntimesteps; nts++) {
    int pellet_decays = 0.;
    assert_always(
        fscanf(gridsave_file, "%la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %d ",
               &globals::timesteps[nts].gamma_dep, &globals::timesteps[nts].gamma_dep_discrete,
               &globals::timesteps[nts].positron_dep, &globals::timesteps[nts].positron_dep_discrete,
               &globals::timesteps[nts].positron_emission, &globals::timesteps[nts].eps_positron_ana_power,
               &globals::timesteps[nts].electron_dep, &globals::timesteps[nts].electron_dep_discrete,
               &globals::timesteps[nts].electron_emission, &globals::timesteps[nts].eps_electron_ana_power,
               &globals::timesteps[nts].alpha_dep, &globals::timesteps[nts].alpha_dep_discrete,
               &globals::timesteps[nts].alpha_emission, &globals::timesteps[nts].eps_alpha_ana_power,
               &globals::timesteps[nts].spfission_dep_discrete, &globals::timesteps[nts].eps_spfission_ana_power,
               &globals::timesteps[nts].qdot_betaminus, &globals::timesteps[nts].qdot_alpha,
               &globals::timesteps[nts].qdot_spfission, &globals::timesteps[nts].qdot_total,
               &globals::timesteps[nts].gamma_emission, &pellet_decays) == 22);
    globals::timesteps[nts].pellet_decays = pellet_decays;
  }

  int timestep_in = 0;
  assert_always(fscanf(gridsave_file, "%d ", &timestep_in) == 1);
  assert_always(timestep_in == timestep);

  for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    int mgi_in = -1;
    float T_R = 0.;
    float T_e = 0.;
    float W = 0.;
    float T_J = 0.;
    int thick = 0;

    assert_always(fscanf(gridsave_file, "%d %a %a %a %a %d %la %la %la %la %a %a", &mgi_in, &T_R, &T_e, &W, &T_J,
                         &thick, &globals::dep_estimator_gamma[nonemptymgi],
                         &globals::dep_estimator_positron[nonemptymgi], &globals::dep_estimator_electron[nonemptymgi],
                         &globals::dep_estimator_alpha[nonemptymgi], &nne_allcells[nonemptymgi],
                         &nnetot_allcells[nonemptymgi]) == 12);

    if (mgi_in != mgi) {
      printlnlog("[fatal] read_grid_restart_data: cell mismatch in reading input gridsave.dat ... abort");
      printlnlog("[fatal] read_grid_restart_data: read cellnumber {}, expected cellnumber {}", mgi_in, mgi);
      assert_always(mgi_in == mgi);
    }

    assert_always(T_R >= 0.);
    assert_always(T_e >= 0.);
    assert_always(W >= 0.);
    assert_always(T_J >= 0.);
    assert_always(globals::dep_estimator_gamma[nonemptymgi] >= 0.);
    assert_always(globals::dep_estimator_positron[nonemptymgi] >= 0.);
    assert_always(globals::dep_estimator_electron[nonemptymgi] >= 0.);
    assert_always(globals::dep_estimator_alpha[nonemptymgi] >= 0.);

    set_TR(nonemptymgi, T_R);
    set_Te(nonemptymgi, T_e);
    set_W(nonemptymgi, W);
    set_TJ(nonemptymgi, T_J);
    thick_allcells[nonemptymgi] = thick;

    if constexpr (USE_LUT_PHOTOION) {
      for (int i = 0; i < globals::nbfcontinua_ground; i++) {
        const int estimindex = (nonemptymgi * globals::nbfcontinua_ground) + i;
        assert_always(fscanf(gridsave_file, " %la %la", &globals::corrphotoionrenorm[estimindex],
                             &globals::gammaestimator[estimindex]) == 2);
      }
    }
  }

  // the order of these calls is very important!
  radfield::read_restart_data(gridsave_file);
  if (globals::rank_in_node == 0) {
    // all data is shared on the node
    nonthermal::read_restart_data(gridsave_file);
    nltepop_read_restart_data(gridsave_file);
  }
  MPI_Barrier(globals::mpi_comm_node);
  fclose(gridsave_file);
}

// Assign temperatures to the grid cells at the start of the simulation by assuming that all radioactive decay since
// the snapshot time (plus any snapshot initial cell energy) energy is fully trapped
void assign_initial_temperatures() {
  MPI_Barrier(MPI_COMM_WORLD);

  // We assume that for early times the material is so optically thick, that
  // all the radiation is trapped in the cell it originates from. This
  // means furthermore LTE, so that both temperatures can be evaluated
  // according to the local energy density resulting from the 56Ni decay.
  // The dilution factor is W=1 in LTE.

  printlnlog("Assigning initial temperatures...");

  const double tstart = globals::timesteps[0].mid;
  int cells_below_mintemp = 0;
  int cells_above_maxtemp = 0;

  for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);

    const auto q = (INITIAL_PACKETS_ON && USE_MODEL_INITIAL_ENERGY) ? get_initenergyq(mgi) : 0.;
    const double decayedenergy_per_mass =
        decay::get_endecay_per_ejectamass_tmodel_to_time_withexpansion(nonemptymgi, tstart) + q;

    auto T_initial = static_cast<float>(
        pow(CLIGHT / 4 / STEBO * pow(globals::tmin / tstart, 3) * get_rho_tmin(mgi) * decayedenergy_per_mass, 1. / 4.));

    if (T_initial < MINTEMP) {
      T_initial = MINTEMP;
      cells_below_mintemp++;
    } else if (T_initial > MAXTEMP) {
      T_initial = MAXTEMP;
      cells_above_maxtemp++;
    } else if (!std::isfinite(T_initial)) {
      printlnlog("mgi {}: T_initial of {:g} is infinite!", mgi, T_initial);
    }

    if (globals::rank_in_node == 0) {
      // set the initial temperatures in the modelgrid
      // this is only done by the node master, so that the values are shared
      // in the node shared memory
      set_Te(nonemptymgi, T_initial);
      set_TJ(nonemptymgi, T_initial);
      set_TR(nonemptymgi, T_initial);
      set_W(nonemptymgi, 1.);
      thick_allcells[nonemptymgi] = 0;
    }
  }
  printlnlog("  cells below MINTEMP {:g}: {}", MINTEMP, cells_below_mintemp);
  printlnlog("  cells above MAXTEMP {:g}: {}", MAXTEMP, cells_above_maxtemp);
  MPI_Barrier(MPI_COMM_WORLD);
}

// start at mgi_start and find the next non-empty cell, or return -1 if none found
[[nodiscard]] auto get_next_nonemptymgi(const int mgi_start) -> int {
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

  ranks_nstart.resize(nprocesses, -1);
  ranks_nstart_nonempty.resize(nprocesses, -1);
  ranks_ndo.resize(nprocesses, 0);
  ranks_ndo_nonempty.resize(nprocesses, 0);

  // begin with no cell assignments
  std::ranges::fill(ranks_nstart, 0);
  std::ranges::fill(ranks_nstart_nonempty, 0);
  std::ranges::fill(ranks_ndo, 0);
  std::ranges::fill(ranks_ndo_nonempty, 0);

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
        ranks_nstart_nonempty[rank] = get_next_nonemptymgi(mgi);
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
    auto fileout = std::ofstream("modelgridrankassignments.out");
    assert_always(fileout.is_open());
    fileout << "#rank nstart ndo ndo_nonempty\n";
    for (int r = 0; r < nprocesses; r++) {
      assert_always(ranks_ndo_nonempty[r] <= ranks_ndo[r]);
      fileout << r << ' ' << ranks_nstart[r] << ' ' << ranks_ndo[r] << ' ' << ranks_ndo_nonempty[r] << '\n';
    }
  }
}

// set up a uniform cuboidal grid.
void setup_grid_cartesian_3d() {
  // vmax is per coordinate, but the simulation volume corners will
  // have a higher expansion velocity than the sides
  const double vmax_corner = sqrt(3 * pow(globals::vmax, 2));
  printlnlog("corner vmax {:g} [cm/s] ({:.2f}c)", vmax_corner, vmax_corner / CLIGHT);
  if (!FORCE_SPHERICAL_ESCAPE_SURFACE) {
    assert_always(vmax_corner < CLIGHT);
  }

  // Set grid size for uniform xyz grid
  if (get_modelgridtype() == GridType::CARTESIAN3D) {
    // if we used in a 3D ejecta model, the propagation grid will match the input grid exactly
    // in case the user specified a grid size, we should ensure that it matches
    assert_always(ncoordgrid[0] > 0);
    assert_always(ncoordgrid[1] > 0);
    assert_always(ncoordgrid[2] > 0);
  } else {
    ncoordgrid = {CUBOID_NCOORDGRID_X, CUBOID_NCOORDGRID_Y, CUBOID_NCOORDGRID_Z};
  }

  // artis assumes in some places that the cells are cubes, not cubioids
  assert_always(ncoordgrid[0] == ncoordgrid[1]);
  assert_always(ncoordgrid[0] == ncoordgrid[2]);

  ngrid = ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2];
  resize_exactly(propcell_pos_min, ngrid);

  std::array<int, 3> nxyz = {0, 0, 0};
  for (int n = 0; n < ngrid; n++) {
    for (int axis = 0; axis < 3; axis++) {
      assert_always(nxyz[axis] == get_cellcoordpointnum(n, axis));
      propcell_pos_min[n][axis] = -globals::rmax + (2 * nxyz[axis] * globals::rmax / ncoordgrid[axis]);
    }

    assert_always(n == ((nxyz[2] * ncoordgrid[1]) * ncoordgrid[2]) + ((nxyz[1] * ncoordgrid[0]) + nxyz[0]));

    nxyz[0]++;  // increment x coordinate
    if (nxyz[0] == ncoordgrid[0]) {
      nxyz[0] = 0;
      nxyz[1]++;  // increment y coordinate
    }
    if (nxyz[1] == ncoordgrid[1]) {
      nxyz[1] = 0;
      nxyz[2]++;  // increment z coordinate
    }
  }
}

void setup_grid_spherical_1d() {
  assert_always(get_modelgridtype() == GridType::SPHERICAL1D);

  ncoordgrid = {get_npts_model(), 1, 1};

  ngrid = ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2];

  resize_exactly(propcell_pos_min, ngrid);

  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const int mgi = cellindex;  // interchangeable in this mode
    const double v_inner = mgi > 0 ? vout_model[mgi - 1] : 0.;
    propcell_pos_min[cellindex] = {v_inner * globals::tmin, 0., 0.};
  }
}

void setup_grid_cylindrical_2d() {
  const double vmax_corner = sqrt(2 * pow(globals::vmax, 2));
  printlnlog("corner vmax {:g} [cm/s] ({:.2f}c)", vmax_corner, vmax_corner / CLIGHT);
  assert_always(vmax_corner < CLIGHT);

  assert_always(get_modelgridtype() == GridType::CYLINDRICAL2D);

  ncoordgrid = ncoord_model;

  ngrid = ncoordgrid[0] * ncoordgrid[1];
  assert_always(ngrid == get_npts_model());

  resize_exactly(propcell_pos_min, ngrid);

  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const int n_rcyl = get_cellcoordpointnum(cellindex, 0);
    const int n_z = get_cellcoordpointnum(cellindex, 1);

    propcell_pos_min[cellindex] = {n_rcyl * globals::rmax / ncoord_model[0],
                                   globals::rmax * (-1 + (n_z * 2. / ncoord_model[1])), 0.};
  }
}

auto get_grid_type_name(const GridType gridtype) -> std::string {
  switch (gridtype) {
    case GridType::SPHERICAL1D:
      return "spherical";
    case GridType::CYLINDRICAL2D:
      return "cylindrical";
    case GridType::CARTESIAN3D:
      return "uniform cuboidal";
    default: {
      return "unknown";
    }
  }
}

// Get the discrete index of the coordinate value (where pos must be position in grid coordinate system, not
// necessarily xyz)
auto get_poscoordpointnum(const double pos, const double time, const int axis) -> int {
  if (get_propgridtype() == GridType::CARTESIAN3D) {
    const auto idx = static_cast<int>(((pos / time) + globals::vmax) / 2 / globals::vmax * ncoordgrid[axis]);
    assert_always(idx >= 0);
    assert_always(idx < ncoordgrid[axis]);
    return idx;
  }

  if (get_propgridtype() == GridType::CYLINDRICAL2D) {
    if (axis == 0) {
      const auto n_rcyl = static_cast<int>(pos / time / globals::vmax * ncoordgrid[axis]);
      assert_always(n_rcyl >= 0);
      assert_always(n_rcyl < ncoordgrid[axis]);
      return n_rcyl;
    }
    if (axis == 1) {
      const auto n_z = static_cast<int>(((pos / time) + globals::vmax) / 2 / globals::vmax * ncoordgrid[axis]);
      assert_always(n_z >= 0);
      assert_always(n_z < ncoordgrid[axis]);
      return n_z;
    }
  }

  if (get_propgridtype() == GridType::SPHERICAL1D) {
    // radial spacing is non-uniform, so we have to do a search
    const auto trat = time / globals::tmin;
    for (int n_r = 0; n_r < ncoordgrid[0]; n_r++) {
      if ((pos < grid::get_cellcoordmax(n_r, 0) * trat) && (pos >= grid::get_cellcoordmin(n_r, 0) * trat)) {
        return n_r;
      }
    }
  }

  assert_always(false);
  return -1;
}

// Convert a position vector from Cartesian xyz to the grid coordinate system
[[nodiscard]] constexpr auto get_gridcoords_from_xyz(const Vec3d& pos_xyz, const GridType gridtype) -> Vec3d {
  switch (gridtype) {
    case GridType::CARTESIAN3D:
      return pos_xyz;
    case GridType::CYLINDRICAL2D:
      return {std::sqrt(pow2(pos_xyz[0]) + pow2(pos_xyz[1])), pos_xyz[2], NAN};
    case GridType::SPHERICAL1D:
      return {vec_len(pos_xyz), NAN, NAN};
  }
  assert_always(false);
  return {NAN, NAN, NAN};
}

// get the velocity in the grid coordinate system from the xyz position and direction
[[nodiscard]] constexpr auto get_gridcoords_vel_from_xyz_pos_dir(const Vec3d& pos_xyz, const Vec3d& dir_xyz,
                                                                 const std::span<const double> pktposgridcoord,
                                                                 const GridType gridtype) -> Vec3d {
  switch (gridtype) {
    case GridType::CARTESIAN3D:
      return {dir_xyz[0] * CLIGHT_PROP, dir_xyz[1] * CLIGHT_PROP, dir_xyz[2] * CLIGHT_PROP};
    case GridType::CYLINDRICAL2D: {
      const double v_rcyl = ((pos_xyz[0] * dir_xyz[0]) + (pos_xyz[1] * dir_xyz[1])) / pktposgridcoord[0] * CLIGHT_PROP;
      const double v_z = dir_xyz[2] * CLIGHT_PROP;
      return {v_rcyl, v_z, NAN};
    }
    case GridType::SPHERICAL1D: {
      const double v_radial = dot(pos_xyz, dir_xyz) / pktposgridcoord[0] * CLIGHT_PROP;
      return {v_radial, NAN, NAN};
    }
  }
  assert_always(false);
  return {NAN, NAN, NAN};
}

// find the closest forward distance to the intersection of a ray with an expanding spherical shell (pos and dir are
// 2-vectors or 3-vectors) or expanding circle (2D vectors)
// returns -1 if there are no forward intersections (or if the intersection
// is tangential to the shell)
template <BoundaryType boundarytype, size_t S1>
[[gnu::pure]] [[nodiscard]] constexpr auto expanding_shell_intersection(const std::array<double, S1>& pos,
                                                                        const std::array<double, S1>& dir,
                                                                        const double speed,
                                                                        const double shellradiuststart,
                                                                        const double tstart) -> double {
  static_assert(S1 == 2 || S1 == 3);
  assert_testmodeonly(shellradiuststart > 0);

  // quadratic equation for intersection of ray with sphere
  // a*d^2 + b*d + c = 0
  const double a = dot(dir, dir) - pow2(shellradiuststart / tstart / speed);
  const double b = 2 * (dot(dir, pos) - (pow2(shellradiuststart) / tstart / speed));
  const double c = dot(pos, pos) - pow2(shellradiuststart);

  const double discriminant = pow2(b) - (4 * a * c);

  if (discriminant < 0) {
    // no intersection
    assert_testmodeonly(boundarytype == BoundaryType::INNER);
    assert_testmodeonly(shellradiuststart < vec_len(pos));
    return -1;
  }

  if (discriminant > 0) {
    // two intersections
    double dist1 = (-b + sqrt(discriminant)) / 2 / a;
    double dist2 = (-b - sqrt(discriminant)) / 2 / a;

    const auto [posfinal1, posfinal2] = [&]() {
      std::array<double, S1> posf1{};
      std::array<double, S1> posf2{};
      for (auto d = 0ZU; d < S1; d++) {
        posf1[d] = pos[d] + (dist1 * dir[d]);
        posf2[d] = pos[d] + (dist2 * dir[d]);
      }
      return std::tuple{posf1, posf2};
    }();

    const double v_rad_shell = shellradiuststart / tstart;
    const double v_rad_final1 = dot(dir, posfinal1) * speed / vec_len(posfinal1);
    const double v_rad_final2 = dot(dir, posfinal2) * speed / vec_len(posfinal2);

    // invalidate any solutions that require entering the boundary from the wrong radial direction
    if constexpr (boundarytype == BoundaryType::INNER) {
      // if the packet's radial velocity at intersection is greater than the inner shell's radial velocity,
      // then it is catching up from below the inner shell and should pass through it
      if (v_rad_final1 > v_rad_shell) {
        dist1 = -1;
      }
      if (v_rad_final2 > v_rad_shell) {
        dist2 = -1;
      }
    } else {
      // if the packet's radial velocity at intersection is less than the outer shell's radial velocity,
      // then it is coming from above the outer shell and should pass through it
      if (v_rad_final1 < v_rad_shell) {
        dist1 = -1;
      }
      if (v_rad_final2 < v_rad_shell) {
        dist2 = -1;
      }
    }

#if (TESTMODE)
    if (dist1 >= 0) {
      const double shellradiusfinal1 = shellradiuststart / tstart * (tstart + (dist1 / speed));
      assert_testmodeonly(fabs((vec_len(posfinal1) / shellradiusfinal1) - 1.) < 1e-3);
    }
    if (dist2 >= 0) {
      const double shellradiusfinal2 = shellradiuststart / tstart * (tstart + (dist2 / speed));
      assert_testmodeonly(fabs((vec_len(posfinal2) / shellradiusfinal2) - 1.) < 1e-3);
    }
#endif

    // negative d means in the reverse direction along the ray
    // ignore negative d values, and if two are positive then return the smaller one
    if (dist1 < 0 && dist2 < 0) {
      return -1;
    }
    if (dist2 < 0) {
      return dist1;
    }
    if (dist1 < 0) {
      return dist2;
    }
    return std::min(dist1, dist2);

  }  // exactly one intersection

  // one intersection
  // ignore this and don't change which cell the packet is in
  assert_testmodeonly(shellradiuststart <= vec_len(pos));
  return -1.;
}

}  // anonymous namespace

// for a uniform grid get the the extent along the x,y,z coordinate (x_2 - x_1, etc.) at time tmin
// for spherical grid get the radial extent (r_outer - r_inner) at time tmin
[[gnu::pure]] [[nodiscard]] auto propcell_width_tmin(const int cellindex, const int axis) -> double {
  switch (get_propgridtype()) {
    case GridType::CARTESIAN3D:
      return 2 * globals::rmax / ncoordgrid[axis];

    case GridType::CYLINDRICAL2D:
      return (axis == 0) ? globals::rmax / ncoordgrid[axis] : 2 * globals::rmax / ncoordgrid[axis];

    case GridType::SPHERICAL1D: {
      const int modelgridindex = cellindex;
      const double v_inner = modelgridindex > 0 ? vout_model[modelgridindex - 1] : 0.;
      return (vout_model[modelgridindex] - v_inner) * globals::tmin;
    }
  }
  assert_always(false);
  return NAN;
}

// return the model cell volume (when mapped to the propagation cells) at globals::tmin
// for a uniform cubic grid this is constant
[[gnu::pure]] [[nodiscard]] auto get_modelcell_assocvolume_tmin(const int modelgridindex) -> double {
  switch (get_propgridtype()) {
    case GridType::CARTESIAN3D:
      return (propcell_width_tmin(modelgridindex, 0) * propcell_width_tmin(modelgridindex, 1) *
              propcell_width_tmin(modelgridindex, 2)) *
             get_numpropcells(modelgridindex);

    case GridType::CYLINDRICAL2D: {
      return propcell_width_tmin(modelgridindex, 1) * PI *
             (pow(get_cellcoordmax(modelgridindex, 0), 2) - pow(get_cellcoordmin(modelgridindex, 0), 2));
    }

    case GridType::SPHERICAL1D: {
      return 4. / 3. * PI * (pow(get_cellcoordmax(modelgridindex, 0), 3) - pow(get_cellcoordmin(modelgridindex, 0), 3));
    }
  }
  assert_always(false);
  return NAN;
}

// return the propagation cell volume at globals::tmin
// for a spherical grid, the cell index is required (and should be equivalent to a modelgridindex)
[[gnu::pure]] [[nodiscard]] auto get_propcell_volume_tmin(const int cellindex) -> double {
  const auto prop_gridtype = get_propgridtype();
  if (prop_gridtype == GridType::CARTESIAN3D) {
    return propcell_width_tmin(cellindex, 0) * propcell_width_tmin(cellindex, 1) * propcell_width_tmin(cellindex, 2);
  }

  // 2D and 1D with direct mapping to propagation cells
  const int mgi = get_propcell_modelgridindex(cellindex);
  return get_modelcell_assocvolume_tmin(mgi);
}

// get the maximum position value of a coordinate axis at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the maximum x position in xyz coords, or the maximum radius in spherical 1D
[[gnu::pure]] [[nodiscard]] auto get_cellcoordmax(const int cellindex, const int axis) -> double {
  return grid::get_cellcoordmin(cellindex, axis) + grid::propcell_width_tmin(cellindex, axis);
}

// get the minimum value of a coordinate at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the minimum x position in xyz coords, or the minimum radius
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordmin(const int cellindex, const int axis) -> double {
  return propcell_pos_min[cellindex][axis];
  // return - coordmax[axis] + (2 * get_cellcoordpointnum(cellindex, axis) * coordmax[axis] / ncoordgrid[axis]);
}

// how much do we change the cellindex to move along a coordinately axis (e.g., the x, y, z directions, or r
// direction)
[[gnu::pure]] [[nodiscard]] auto get_coordcellindexincrement(const int axis) -> int {
  switch (axis) {
    case 0:
      return 1;

    case 1:
      return ncoordgrid[0];

    case 2:
      return ncoordgrid[0] * ncoordgrid[1];

    default:
      if constexpr (TESTMODE) {
        assert_testmodeonly(false);
      } else {
        __builtin_unreachable();
      }
  }
}

// convert a cell index number into an integer (x,y,z or r) coordinate index from 0 to ncoordgrid[axis]
[[gnu::pure]] [[nodiscard]] auto get_cellcoordpointnum(const int cellindex, const int axis) -> int {
  if (get_propgridtype() == GridType::CARTESIAN3D || get_propgridtype() == GridType::CYLINDRICAL2D) {
    switch (axis) {
      // 3D Cartesian: increment x first, then y, then z
      // 2D Cylindrical: increment r first, then z
      case 0:
        return cellindex % ncoordgrid[0];

      case 1:
        return (cellindex / ncoordgrid[0]) % ncoordgrid[1];

      case 2:
        return (cellindex / (ncoordgrid[0] * ncoordgrid[1])) % ncoordgrid[2];

      default:
        if constexpr (TESTMODE) {
          assert_testmodeonly(false);
        } else {
          __builtin_unreachable();
        }
    }
  }

  if (get_propgridtype() == GridType::SPHERICAL1D) {
    return cellindex;
  }

  assert_always(false);
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

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nnetot(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return nnetot_allcells[nonemptymgi];
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_ffegrp(const int modelgridindex) -> float {
  return modelgrid_input[modelgridindex].ffegrp;
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_modelcell_mean_radial_pos(const int modelgridindex,
                                                                           const double tratmid) -> double {
  const int assoc_cells = grid::get_numpropcells(modelgridindex);
  return modelgrid_input[modelgridindex].initial_radial_pos_sum * tratmid / assoc_cells;
}

// mass fraction of an element (all isotopes combined)
[[gnu::pure]] [[nodiscard]] auto get_elem_abundance(const std::ptrdiff_t nonemptymgi, const int element) -> float {
  const auto massfrac = elem_massfracs_allcells[(nonemptymgi * get_nelements()) + element];
  assert_testmodeonly(massfrac >= 0.0);
  return massfrac;
}

// mass fraction of an element (all isotopes combined)
void set_elem_abundance(const ptrdiff_t nonemptymgi, const int element, const float newabundance) {
  elem_massfracs_allcells[(nonemptymgi * get_nelements()) + element] = newabundance;
}

// mass fraction of an element (all isotopes combined)
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_elem_numberdens(const ptrdiff_t nonemptymgi, const int element)
    -> double {
  return get_elem_abundance(nonemptymgi, element) /
         static_cast<double>(grid::get_element_meanweight(nonemptymgi, element)) * grid::get_rho(nonemptymgi);
}

DEVICE_FUNC auto get_kappagrey(const int nonemptymgi) -> float { return kappagrey_allcells[nonemptymgi]; }

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_Te(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return Te_allcells[nonemptymgi];
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_TR(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return TR_allcells[nonemptymgi];
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_TJ(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return TJ_allcells[nonemptymgi];
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_W(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return W_allcells[nonemptymgi];
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

void set_kappagrey(const int nonemptymgi, const float kappagrey) { kappagrey_allcells[nonemptymgi] = kappagrey; }

void set_Te(const int nonemptymgi, const float Te) {
  if (Te > 0.) {
    // ignore the zero initialisation value for this check
    const double nu_peak = 5.879e10 * Te;
    if (nu_peak > NU_MAX_R || nu_peak < NU_MIN_R) {
      const auto modelgridindex = get_mgi_of_nonemptymgi(nonemptymgi);
      printlnlog(
          "[warning] modelgridindex {} B_planck(Te={:g} K) peak at {:g} Hz is outside frequency range NU_MIN_R {:g} "
          "NU_MAX_R {:g}",
          modelgridindex, Te, nu_peak, NU_MIN_R, NU_MAX_R);
    }
  }

  Te_allcells[nonemptymgi] = Te;
}

void set_TR(const int nonemptymgi, const float TR) { TR_allcells[nonemptymgi] = TR; }

void set_TJ(const int nonemptymgi, const float TJ) { TJ_allcells[nonemptymgi] = TJ; }

void set_W(const int nonemptymgi, const float W) { W_allcells[nonemptymgi] = W; }

DEVICE_FUNC auto get_modelgridtype() -> GridType {
  assert_testmodeonly(model_type.has_value());
  return model_type.value();
}

[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_npts_model() -> int
// number of model grid cells
{
  assert_testmodeonly(npts_model > 0);
  return static_cast<int>(npts_model);
}

// number of model grid cells
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
  assert_testmodeonly(modelgridindex <= get_npts_model());
  return modelgrid_numpropcells[modelgridindex];
}

// get the index in the list of non-empty cells for a given model grid cell
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nonemptymgi_of_mgi(const int mgi) -> int {
  assert_testmodeonly(get_nonempty_npts_model() > 0);
  assert_testmodeonly(mgi >= 0);
  assert_testmodeonly(mgi < get_npts_model());

  const int nonemptymgi = nonemptymgi_of_mgi[mgi];
  // assert_testmodeonly(nonemptymgi >= 0 || get_numpropcells(mgi) == 0);
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());

  return nonemptymgi;
}

// get the index in the list of non-empty cells for a given model grid cell
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

auto get_otherstable_initabund(const std::ptrdiff_t nonemptymgi, const int element) -> float {
  return initmassfracuntrackedstable_allcells[(nonemptymgi * get_nelements()) + element];
}

// get element mean weight in grams
auto get_element_meanweight(const std::ptrdiff_t nonemptymgi, const int element) -> float {
  if constexpr (USE_CALCULATED_MEANATOMICWEIGHT) {
    const auto mu = elem_meanweight_allcells[(nonemptymgi * get_nelements()) + element];
    if (mu > 0) {
      return mu;
    }
  }
  return globals::elements[element].initstablemeannucmass;
}

// set element weight in grams
void set_element_meanweight(const std::ptrdiff_t nonemptymgi, const int element, const float meanweight) {
  assert_always(meanweight > 0.);
  elem_meanweight_allcells[(nonemptymgi * get_nelements()) + element] = meanweight;
}

auto get_electronfrac(const int nonemptymgi) -> double {
  double nucleondens = 0.;
  for (int element = 0; element < get_nelements(); element++) {
    nucleondens += get_elem_numberdens(nonemptymgi, element) * get_element_meanweight(nonemptymgi, element) / MH;
  }
  return get_nnetot(nonemptymgi) / nucleondens;
}

// q: energy in the model at tmin per gram to use with USE_MODEL_INITIAL_ENERGY option [erg/g]
DEVICE_FUNC auto get_initenergyq(const int modelgridindex) -> double {
  return modelgrid_input[modelgridindex].initenergyq;
}

// get the radial distance from the origin to the centre of the cell at time tmin
auto get_cellradialposmid(const int cellindex) -> double {
  const auto prop_gridtype = get_propgridtype();
  if (prop_gridtype == GridType::SPHERICAL1D) {
    // volume averaged mean radius is slightly complex for radial shells
    const double r_inner = grid::get_cellcoordmin(cellindex, 0);
    const double r_outer = r_inner + grid::propcell_width_tmin(cellindex, 0);
    return 3. / 4 * (pow(r_outer, 4.) - pow(r_inner, 4.)) / (pow(r_outer, 3) - pow(r_inner, 3.));
  }

  if (prop_gridtype == GridType::CYLINDRICAL2D) {
    const double rcyl_mid = get_cellcoordmin(cellindex, 0) + (0.5 * propcell_width_tmin(cellindex, 0));
    const double z_mid = get_cellcoordmin(cellindex, 1) + (0.5 * propcell_width_tmin(cellindex, 1));
    return std::sqrt(std::pow(rcyl_mid, 2) + std::pow(z_mid, 2));
  }

  // cubic grid requires taking the length of the 3D position vector
  Vec3d dcen{};
  for (int axis = 0; axis < 3; axis++) {
    dcen[axis] = get_cellcoordmin(cellindex, axis) + (0.5 * propcell_width_tmin(cellindex, axis));
  }

  return vec_len(dcen);
}

[[nodiscard]] auto get_elements_uppermost_ion(const int nonemptymgi, const int element) -> int {
  const auto uppermost_ion = elements_uppermost_ion_allcells[(nonemptymgi * get_nelements()) + element];
  assert_testmodeonly(uppermost_ion >= 0);
  assert_testmodeonly(uppermost_ion <= std::max(0, get_nions(element) - 1));
  return uppermost_ion;
}

void set_elements_uppermost_ion(const int nonemptymgi, const int element, const int uppermost_ion) {
  assert_testmodeonly(uppermost_ion <= std::max(0, get_nions(element) - 1));
  elements_uppermost_ion_allcells[(nonemptymgi * get_nelements()) + element] = uppermost_ion;
}

void calculate_kappagrey() {
  double rho_sum = 0.;
  double fe_sum = 0.;
  double opcase3_sum = 0.;
  const int empty_cells = 0;

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const auto mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    rho_sum += get_rho_tmin(mgi);
    fe_sum += get_ffegrp(mgi);

    if (globals::opacity_case == 3) {
      if (get_rho_tmin(mgi) > 0.) {
        const auto kappagrey =
            static_cast<float>(((0.9 * get_ffegrp(mgi)) + 0.1) *
                               ((get_rho_tmin(mgi) > globals::rho_crit) ? globals::rho_crit / get_rho_tmin(mgi) : 1.));

        set_kappagrey(nonemptymgi, kappagrey);
      } else if (get_rho_tmin(mgi) == 0.) {
        set_kappagrey(nonemptymgi, 0.);
      } else if (get_rho_tmin(mgi) < 0.) {
        printlnlog("Error: negative density. Abort.");
        std::abort();
      }
      opcase3_sum += get_kappagrey(nonemptymgi) * get_rho_tmin(mgi);
    }
  }

  // Second pass through allows calculation of normalized chi_grey
  double check1 = 0.;
  double check2 = 0.;
  for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
    if (get_rho_tmin(mgi) > 0) {
      double kappa = 0.;
      if (globals::opacity_case == 0) {
        kappa = globals::GREY_OP;
      } else if (globals::opacity_case == 1 || globals::opacity_case == 4) {
        // kappagrey used for initial grey approximation in case 4
        kappa = ((0.9 * get_ffegrp(mgi)) + 0.1) * globals::GREY_OP / ((0.9 * mfegroup / mtot_input) + 0.1);
      } else if (globals::opacity_case == 2) {
        const double opcase2_normal = globals::GREY_OP * rho_sum / ((0.9 * fe_sum) + (0.1 * (ngrid - empty_cells)));
        kappa = opcase2_normal / get_rho_tmin(mgi) * ((0.9 * get_ffegrp(mgi)) + 0.1);
      } else if (globals::opacity_case == 3) {
        globals::opcase3_normal = globals::GREY_OP * rho_sum / opcase3_sum;
        kappa = get_kappagrey(nonemptymgi) * globals::opcase3_normal;
      } else if (globals::opacity_case == 5) {
        // electron-fraction-dependent opacities
        // values from table 1 of Tanaka et al. (2020).
        const auto Ye = modelgrid_input[mgi].initelectronfrac;
        assert_always(Ye > 0.);

        if (Ye <= 0.1) {
          kappa = 19.5;
        } else if (Ye <= 0.15) {
          kappa = 32.2;
        } else if (Ye <= 0.20) {
          kappa = 22.3;
        } else if (Ye <= 0.25) {
          kappa = 5.6;
        } else if (Ye <= 0.30) {
          kappa = 5.36;
        } else if (Ye <= 0.35) {
          kappa = 3.3;
        } else {
          kappa = 0.96;
        }
      } else if (globals::opacity_case == 6) {
        // grey opacity used in Just+2022, https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.2820J/abstract
        // kappa is a simple analytic function of temperature and lanthanide mass fraction
        // adapted to best fit lightcurves from Kasen+2017 in ALCAR simulations
        const double T_rad = get_TR(nonemptymgi);
        double X_lan = 0.;
        for (int element = 0; element < get_nelements(); element++) {
          const int z = get_atomicnumber(element);
          if (z >= 57 && z <= 71) {
            X_lan += get_elem_abundance(nonemptymgi, element);
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
          kappa *= pow(T_rad / 2000., 5.);
        }
      } else {
        assert_always(false);
      }

      set_kappagrey(nonemptymgi, static_cast<float>(kappa));
    } else {
      set_kappagrey(nonemptymgi, 0.);
    }

    check1 = check1 + (get_kappagrey(nonemptymgi) * get_rho_tmin(mgi));
    check2 = check2 + get_rho_tmin(mgi);
  }

  printlnlog("Grey normalisation check: {:g}", check1 / check2);
}

void read_ejecta_model() {
  auto fmodel = fstream_required("model.txt", std::ios::in);
  std::string line;
  std::optional<GridType> detected_dim{};

  // two integers on the first line of the model file
  int npts_0 = 0;  // total model points for 1D/3D, and number of points in r for 2D
  int npts_1 = 0;  // number of points in z for 2D
  assert_always(get_noncommentline(fmodel, line));
  auto ssline = std::istringstream{line};
  ssline >> npts_0;
  if (ssline >> npts_1) {
    // second number on the line for 2D means the line was n_r n_z
    detected_dim = GridType::CYLINDRICAL2D;
    printlnlog("Detected 2D model");
    ssline >> npts_1;  // r and z (cylindrical polar)
    npts_model = npts_0 * npts_1;
  } else {
    // for 1D and 3D, this was the total number of model cells
    npts_model = npts_0;
  }

  // Now read the time (in days) at which the model is specified.
  double t_model_days{NAN};
  assert_always(get_noncommentline(fmodel, line));
  std::istringstream{line} >> t_model_days;
  t_model = t_model_days * DAY;
  assert_always(globals::tmin >= t_model);

  const auto pos_after_t_model = fmodel.tellg();
  // if the next line is a single float, it is the vmax (so 2D or 3D)
  // otherwise, it is the first line of the model or a header comment (so 1D)
  std::getline(fmodel, line);
  if (!line.starts_with("#")) {
    double num_after_vmax{NAN};
    auto sslinevmax = std::istringstream{line};
    if ((sslinevmax >> globals::vmax) && !(sslinevmax >> num_after_vmax)) {
      // single value on the line is a vmax, so 2D or 3D
      // if it's not already know to be 2D (based on n_r n_z line), then it's 3D
      if (detected_dim != GridType::CYLINDRICAL2D) {
        assert_always(!detected_dim.has_value());
        detected_dim = GridType::CARTESIAN3D;
        printlnlog("Detected 3D model");
      }
    }
  }
  if (!detected_dim.has_value()) {
    assert_always(!detected_dim.has_value());
    detected_dim = GridType::SPHERICAL1D;
    printlnlog("Detected 1D model");
    fmodel.seekg(pos_after_t_model);
  }

  assert_always(detected_dim.has_value());
  model_type = detected_dim;

  assert_always(modelgrid_input.empty());
  modelgrid_input = MPI_shared_malloc_span<ModelGridCellInput>(npts_model + 1, ModelGridCellInput{});
  modelgrid_numpropcells.resize(npts_model + 1, 0);
  nonemptymgi_of_mgi.resize(npts_model + 1, -1);

  if (get_modelgridtype() == GridType::SPHERICAL1D) {
    ncoord_model[0] = npts_0;
    ncoord_model[1] = 0;
    ncoord_model[2] = 0;
    vout_model.resize(get_npts_model(), NAN);

    // Now read in the lines of the model. Each line has 5 entries: the
    // cell number (integer) the velocity at outer boundary of cell (float),
    // the mass density in the cell (float), the abundance of Ni56 by mass
    // in the cell (float) and the total abundance of all Fe-grp elements
    // in the cell (float). For now, the last number is recorded but never
    // used.

    const auto [colnames, nucindexlist, one_line_per_cell] = read_model_columns(fmodel);

    int mgi = 0;
    while (std::getline(fmodel, line)) {
      double vout_kmps{NAN};
      double log_rho{NAN};
      int cellnumberin = 0;
      ssline.clear();
      ssline.str(line);

      if (!(ssline >> cellnumberin >> vout_kmps >> log_rho)) {
        printlnlog("Unexpected number of values in model.txt");
        printlnlog("line: {}", line);
        assert_always(false);
      }

      if (mgi == 0) {
        first_cellindex = cellnumberin;
        printlnlog("first_cellindex {}", first_cellindex);
      }
      assert_always(cellnumberin == mgi + first_cellindex);

      vout_model[mgi] = vout_kmps * 1.e5;

      const auto rho_tmin =
          static_cast<float>(log_rho > -90 ? pow(10., log_rho) * pow(t_model / globals::tmin, 3) : 0.);
      set_rho_tmin(mgi, rho_tmin);
      const bool keepcell = (rho_tmin > 0);
      read_model_radioabundances(fmodel, ssline, mgi, keepcell, colnames, nucindexlist, one_line_per_cell);

      mgi += 1;
      if (mgi == get_npts_model()) {
        break;
      }
    }

    if (mgi != get_npts_model()) {
      printlnlog("ERROR in model.txt. Found only {} cells instead of {} expected.", mgi - 1, get_npts_model());
      std::abort();
    }

    globals::vmax = vout_model[get_npts_model() - 1];
  } else if (get_modelgridtype() == GridType::CYLINDRICAL2D) {
    ncoord_model[0] = npts_0;
    ncoord_model[1] = npts_1;
    ncoord_model[2] = 0;
    const auto [colnames, nucindexlist, one_line_per_cell] = read_model_columns(fmodel);

    // Now read in the model. Each point in the model has two lines of input.
    // First is an index for the cell then its r-mid point then its z-mid point
    // then its total mass density.
    // Second is the total FeG mass, initial 56Ni mass, initial 56Co mass

    int mgi = 0;
    while (std::getline(fmodel, line)) {
      int cellnumberin = 0;
      float cell_r_in{NAN};
      float cell_z_in{NAN};
      double rho_tmodel{NAN};
      ssline.clear();
      ssline.str(line);
      assert_always(ssline >> cellnumberin >> cell_r_in >> cell_z_in >> rho_tmodel);

      if (mgi == 0) {
        first_cellindex = cellnumberin;
      }
      assert_always(cellnumberin == mgi + first_cellindex);

      const int n_rcyl = (mgi % ncoord_model[0]);
      const double pos_r_cyl_mid = (n_rcyl + 0.5) * globals::vmax * t_model / ncoord_model[0];
      assert_always(fabs((cell_r_in / pos_r_cyl_mid) - 1) < 1e-3);
      const int n_z = (mgi / ncoord_model[0]);
      const double pos_z_mid = globals::vmax * t_model * (-1 + (2 * (n_z + 0.5) / ncoord_model[1]));
      assert_always(fabs((cell_z_in / pos_z_mid) - 1) < 1e-3);

      if (rho_tmodel < 0) {
        printlnlog("negative input density {:g} {}", rho_tmodel, mgi);
        std::abort();
      }

      const bool keepcell = (rho_tmodel > 0);
      const auto rho_tmin = static_cast<float>(rho_tmodel * pow(t_model / globals::tmin, 3));
      set_rho_tmin(mgi, rho_tmin);

      read_model_radioabundances(fmodel, ssline, mgi, keepcell, colnames, nucindexlist, one_line_per_cell);

      mgi++;
    }

    if (mgi != get_npts_model()) {
      printlnlog("ERROR in model.txt. Found {} only cells instead of {} expected.", mgi - 1, get_npts_model());
      std::abort();
    }
  } else if (get_modelgridtype() == GridType::CARTESIAN3D) {
    ncoord_model[0] = static_cast<int>(round(std::cbrt(npts_model)));
    ncoord_model[1] = ncoord_model[0];
    ncoord_model[2] = ncoord_model[0];

    assert_always(ncoord_model[0] * ncoord_model[1] * ncoord_model[2] == npts_model);
    // for a 3D input model, the propagation cells will match the input cells exactly
    ncoordgrid = ncoord_model;
    ngrid = npts_model;

    const double xmax_tmodel = globals::vmax * t_model;

    // Now read in the lines of the model.
    min_den = -1.;

    // check if expected positions match in either xyz or zyx column order
    // set false if a problem is detected
    bool posmatch_xyz = true;
    bool posmatch_zyx = true;

    const auto [colnames, nucindexlist, one_line_per_cell] = read_model_columns(fmodel);

    // mgi is the index to the model grid - empty cells are sent to special value get_npts_model(),
    // otherwise each input cell is one modelgrid cell
    int mgi = 0;  // corresponds to model.txt index column, but zero indexed! (model.txt might be 1-indexed)
    while (std::getline(fmodel, line)) {
      int cellnumberin = 0;
      std::array<float, 3> cellpos_in{};
      float rho_model{NAN};
      ssline.clear();
      ssline.str(line);

      assert_always(ssline >> cellnumberin >> cellpos_in[0] >> cellpos_in[1] >> cellpos_in[2] >> rho_model);

      if (mgi == 0) {
        first_cellindex = cellnumberin;
      }
      assert_always(cellnumberin == mgi + first_cellindex);

      if (mgi % (ncoord_model[1] * ncoord_model[2]) == 0) {
        printlnlog("read up to cell mgi {}", mgi);
      }

      // cell coordinates in the 3D model.txt file are sometimes reordered by the scaling script
      // however, the cellindex always should increment X first, then Y, then Z

      for (int axis = 0; axis < 3; axis++) {
        const double cellwidth = 2 * xmax_tmodel / ncoordgrid[axis];
        const double cellpos_expected = -xmax_tmodel + (cellwidth * get_cellcoordpointnum(mgi, axis));
        if (fabs(cellpos_expected - cellpos_in[axis]) > 0.5 * cellwidth) {
          posmatch_xyz = false;
        }
        if (fabs(cellpos_expected - cellpos_in[2 - axis]) > 0.5 * cellwidth) {
          posmatch_zyx = false;
        }
      }

      if (rho_model < 0) {
        printlnlog("negative input density {:g} {}", rho_model, mgi);
        std::abort();
      }

      // in 3D cartesian, cellindex and modelgridindex are interchangeable
      const bool keepcell = (rho_model > 0);
      const auto rho_tmin = static_cast<float>(rho_model * pow(t_model / globals::tmin, 3));
      set_rho_tmin(mgi, rho_tmin);

      if (min_den < 0. || min_den > rho_model) {
        min_den = rho_model;
      }

      read_model_radioabundances(fmodel, ssline, mgi, keepcell, colnames, nucindexlist, one_line_per_cell);

      mgi++;
    }
    if (mgi != npts_model) {
      printlnlog("ERROR in model.txt. Found {} cells instead of {} expected.", mgi, npts_model);
      std::abort();
    }

    //   assert_always(posmatch_zyx ^ posmatch_xyz);  // xor because if both match then probably an infinity occurred
    if (posmatch_xyz) {
      printlnlog("Cell positions in model.txt are consistent with calculated values when x-y-z column order is used.");
    }
    if (posmatch_zyx) {
      printlnlog("Cell positions in model.txt are consistent with calculated values when z-y-x column order is used.");
    }

    if (!posmatch_xyz && !posmatch_zyx) {
      printlnlog(
          "WARNING: Cell positions in model.txt are not consistent with calculated values in either x-y-z or z-y-x "
          "order.");
    }

    printlnlog("min_den {:g} [g/cm3]", min_den);
  }

  assert_always(get_npts_model() ==
                std::max(1, ncoord_model[0]) * std::max(1, ncoord_model[1]) * std::max(1, ncoord_model[2]));
  printlnlog("npts_model: {}", get_npts_model());
  globals::rmax = globals::vmax * globals::tmin;
  printlnlog("vmax {:g} [cm/s] ({:.2f}c)", globals::vmax, globals::vmax / CLIGHT);
  assert_always(globals::vmax < CLIGHT);
  printlnlog("tmin {:g} [s] = {:.2f} [d]", globals::tmin, globals::tmin / 86400.);
  printlnlog("rmax {:g} [cm] (at t=tmin)", globals::rmax);

  calc_modelinit_totmassnuclides();

  printlnlog("Total input model mass: {:9.3e} [Msun]", mtot_input / MSUN);
  printlnlog("Nuclide masses at t=t_model_init [Msun]:");
  printlnlog("  56Ni: {:9.3e}  56Co: {:9.3e}  52Fe: {:9.3e}  48Cr: {:9.3e}", get_totmassnuclide_tmodel(28, 56) / MSUN,
             get_totmassnuclide_tmodel(27, 56) / MSUN, get_totmassnuclide_tmodel(26, 52) / MSUN,
             get_totmassnuclide_tmodel(24, 48) / MSUN);
  printlnlog("  Fe-group: {:9.3e}  57Ni: {:9.3e}  57Co: {:9.3e}", mfegroup / MSUN,
             get_totmassnuclide_tmodel(28, 57) / MSUN, get_totmassnuclide_tmodel(27, 57) / MSUN);

  read_possible_yefile();
}

void write_grid_restart_data(const int timestep) {
  const auto filename = std::format("gridsave_ts{}.tmp", timestep);

  const auto sys_time_start_write_restart = std::time(nullptr);
  printlog("Write grid restart data to {}...", filename);

  FILE* gridsave_file = fopen_required(filename, "w");

  fprintf(gridsave_file, "%d ", globals::ntimesteps);
  fprintf(gridsave_file, "%d ", globals::nprocs);

  for (int nts = 0; nts < globals::ntimesteps; nts++) {
    fprintf(gridsave_file, "%la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %d ",
            globals::timesteps[nts].gamma_dep, globals::timesteps[nts].gamma_dep_discrete,
            globals::timesteps[nts].positron_dep, globals::timesteps[nts].positron_dep_discrete,
            globals::timesteps[nts].positron_emission, globals::timesteps[nts].eps_positron_ana_power,
            globals::timesteps[nts].electron_dep, globals::timesteps[nts].electron_dep_discrete,
            globals::timesteps[nts].electron_emission, globals::timesteps[nts].eps_electron_ana_power,
            globals::timesteps[nts].alpha_dep, globals::timesteps[nts].alpha_dep_discrete,
            globals::timesteps[nts].alpha_emission, globals::timesteps[nts].eps_alpha_ana_power,
            globals::timesteps[nts].spfission_dep_discrete, globals::timesteps[nts].eps_spfission_ana_power,
            globals::timesteps[nts].qdot_betaminus, globals::timesteps[nts].qdot_alpha,
            globals::timesteps[nts].qdot_spfission, globals::timesteps[nts].qdot_total,
            globals::timesteps[nts].gamma_emission, globals::timesteps[nts].pellet_decays);
  }

  fprintf(gridsave_file, "%d ", timestep);

  for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);

    assert_always(globals::dep_estimator_gamma[nonemptymgi] >= 0.);
    fprintf(gridsave_file, "%d %a %a %a %a %d %la %la %la %la %a %a", mgi, get_TR(nonemptymgi), get_Te(nonemptymgi),
            get_W(nonemptymgi), get_TJ(nonemptymgi), thick_allcells[nonemptymgi],
            globals::dep_estimator_gamma[nonemptymgi], globals::dep_estimator_positron[nonemptymgi],
            globals::dep_estimator_electron[nonemptymgi], globals::dep_estimator_alpha[nonemptymgi],
            nne_allcells[nonemptymgi], nnetot_allcells[nonemptymgi]);

    if constexpr (USE_LUT_PHOTOION) {
      for (int i = 0; i < globals::nbfcontinua_ground; i++) {
        const int estimindex = (nonemptymgi * globals::nbfcontinua_ground) + i;
        fprintf(gridsave_file, " %la %la", globals::corrphotoionrenorm[estimindex],
                globals::gammaestimator[estimindex]);
      }
    }
    fprintf(gridsave_file, "\n");
  }

  // the order of these calls is very important!
  radfield::write_restart_data(gridsave_file);
  nonthermal::write_restart_data(gridsave_file);
  nltepop_write_restart_data(gridsave_file);
  fclose(gridsave_file);
  printlnlog("done in {} seconds.", std::time(nullptr) - sys_time_start_write_restart);
}

auto get_nstart(const int rank) -> int {
  if (ranks_ndo.empty()) {
    setup_nstart_ndo();
  }
  return ranks_nstart[rank];
}

auto get_nstart_nonempty(const int rank) -> int {
  if (ranks_ndo.empty()) {
    setup_nstart_ndo();
  }
  return ranks_nstart_nonempty[rank];
}

auto get_ndo(const int rank) -> int {
  if (ranks_ndo.empty()) {
    setup_nstart_ndo();
  }
  return ranks_ndo[rank];
}

auto get_ndo_nonempty(const int rank) -> int {
  if (ranks_ndo.empty()) {
    setup_nstart_ndo();
  }
  return ranks_ndo_nonempty[rank];
}

// Initialise the propagation grid cells and associate them with modelgrid cells
void init_grid(const int my_rank) {
  // The cells will be ordered by x then y, then z. Call a routine that
  // sets up the initial positions and widths of the cells.
  const auto prop_gridtype = get_propgridtype();
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
    default:
      assert_always(false);
  }
  propcell_mgi.resize(ngrid, -1);

  printlnlog("propagation grid: {}-dimensional {}", get_ndim(prop_gridtype), get_grid_type_name(prop_gridtype));

  for (int d = 0; d < get_ndim(prop_gridtype); d++) {
    printlnlog("    coordinate {} '{}': cells have {} position values", d, get_coordlabel(prop_gridtype, d),
               ncoordgrid[d]);
  }
  printlnlog("    total propagation cells: {}", ngrid);

  // Now set up the density in each cell.

  // Calculate the critical opacity at which opacity_case 3 switches from a
  // regime proportional to the density to a regime independent of the density
  // This is done by solving for tau_sobolev == 1
  // tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/(56 * MH) * 3000e-8 * globals::timesteps[m].mid;
  globals::rho_crit = ME * CLIGHT * 56 * MH / (PI * QE * QE * globals::rho_crit_para * 3000e-8 * globals::tmin);
  printlnlog("grid_init: rho_crit = {:g} [g/cm3]", globals::rho_crit);

  if (get_modelgridtype() == prop_gridtype) {
    if (get_modelgridtype() == GridType::CARTESIAN3D) {
      assert_always(ncoord_model[0] == ncoordgrid[0]);
      assert_always(ncoord_model[1] == ncoordgrid[1]);
      assert_always(ncoord_model[2] == ncoordgrid[2]);
    }

    map_modeltogrid_direct();
  } else if (get_modelgridtype() == GridType::SPHERICAL1D) {
    assert_always(prop_gridtype == GridType::CARTESIAN3D);
    map_1dmodelto3dgrid();
  } else if (get_modelgridtype() == GridType::CYLINDRICAL2D) {
    assert_always(prop_gridtype == GridType::CARTESIAN3D);
    map_2dmodelto3dgrid();
  } else {
    assert_always(false);
  }

  if (globals::my_rank == 0) {
    auto grid_file = std::fstream("grid.out", std::ios::out | std::ios::trunc);
    assert_always(grid_file.is_open());
    for (int cellindex = 0; cellindex < ngrid; cellindex++) {
      const int mgi = get_propcell_modelgridindex(cellindex);
      if (mgi >= 0) {
        grid_file << cellindex << ' ' << mgi << '\n';  // write only non-empty cells to grid file
      }
    }
  }

  allocate_nonemptymodelcells();
  calculate_kappagrey();
  read_abundances();

  const int ndo_nonempty = grid::get_ndo_nonempty(my_rank);

  radfield::init(my_rank, ndo_nonempty);
  nonthermal::init();

  // and assign a temperature to the cells
  if (globals::simulation_continued_from_saved) {
    // For continuation of an existing simulation we read the temperatures
    // at the end of the simulation and write them to the grid.
    read_grid_restart_data(globals::timestep_initial);
  } else {
    assign_initial_temperatures();
  }

  // when mapping 1D spherical or 2D cylindrical model onto cubic grid, scale up the
  // radioactive abundances to account for the missing masses in
  // the model cells that are not associated with any propagation cells
  if (prop_gridtype == GridType::CARTESIAN3D && get_modelgridtype() == GridType::SPHERICAL1D &&
      globals::rank_in_node == 0) {
    for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
      if (totmassnuclide[nucindex] <= 0) {
        continue;
      }

      double totmassnuclide_actual = 0.;
      for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
        const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
        totmassnuclide_actual +=
            get_modelinitnucmassfrac(mgi, nucindex) * get_rho_tmin(mgi) * get_modelcell_assocvolume_tmin(mgi);
      }

      if (totmassnuclide_actual > 0.) {
        const double ratio = totmassnuclide[nucindex] / totmassnuclide_actual;
        for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
          const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
          const double prev_abund = get_modelinitnucmassfrac(mgi, nucindex);
          const auto new_abund = static_cast<float>(prev_abund * ratio);
          set_modelinitnucmassfrac(mgi, nucindex, new_abund);
        }
      }
    }
  }

  MPI_Barrier(globals::mpi_comm_node);

  double mtot_mapped = 0.;
  for (int mgi = 0; mgi < get_npts_model(); mgi++) {
    mtot_mapped += get_rho_tmin(mgi) * get_modelcell_assocvolume_tmin(mgi);
  }
  printlnlog("Total grid-mapped mass: {:9.3e} [Msun] ({:.1f}% of input mass)", mtot_mapped / MSUN,
             mtot_mapped / mtot_input * 100.);

  MPI_Barrier(MPI_COMM_WORLD);
}

auto get_totmassnuclide_tmodel(const int z, const int a) -> double { return totmassnuclide[decay::get_nucindex(z, a)]; }

// identify the cell index from an (x,y,z) position and a time.
[[nodiscard]] DEVICE_FUNC auto get_cellindex_from_pos(const Vec3d& pos, const double time) -> int {
  const auto prop_gridtype = get_propgridtype();
  const auto ndim = get_ndim(prop_gridtype);
  auto posgridcoords = get_gridcoords_from_xyz(pos, prop_gridtype);
  int cellindex = 0;
  for (int d = 0; d < ndim; d++) {
    if (std::abs(posgridcoords[d]) > (globals::vmax * time)) {
      // outside grid
      return -99;
    }
    cellindex += get_coordcellindexincrement(d) * get_poscoordpointnum(posgridcoords[d], time, d);
  }
  assert_always(cellindex >= 0);
  assert_always(cellindex < ngrid);

  // do a check that the position is within the cell
  const double trat = time / globals::tmin;
  for (int n = 0; n < ndim; n++) {
    assert_always(posgridcoords[n] >= ((grid::get_cellcoordmin(cellindex, n) * trat) - 10));
    assert_always(posgridcoords[n] <= ((grid::get_cellcoordmax(cellindex, n) * trat) + 10));
  }
  return cellindex;
}

// compute distance to a cell boundary.
[[nodiscard]] DEVICE_FUNC auto boundary_distance(const Vec3d& dir, const Vec3d& pos, const double tstart,
                                                 const int cellindex) -> std::tuple<double, int> {
  const auto prop_gridtype = get_propgridtype();
  if constexpr (FORCE_SPHERICAL_ESCAPE_SURFACE) {
    if (get_cell_r_inner(cellindex, prop_gridtype) > globals::vmax * globals::tmin) {
      return {0., -99};
    }
  }

  // d is used to loop over the coordinate indices 0,1,2 for x,y,z

  // the following vector are in grid coordinates, so either x,y,z (3D) or r (1D), or r_xy, z (2D)

  const auto pktposgridcoord = get_gridcoords_from_xyz(pos, prop_gridtype);

  // dir * CLIGHT_PROP converted from xyz to grid coordinates
  const auto pktvelgridcoord = get_gridcoords_vel_from_xyz_pos_dir(pos, dir, pktposgridcoord, prop_gridtype);

  const auto cellcoordmax = [cellindex, prop_gridtype]() {
    auto _cellcoordmax = std::array<double, 3>{};  // position at time tmin
    for (int d = 0; d < get_ndim(prop_gridtype); d++) {
      _cellcoordmax[d] = grid::get_cellcoordmax(cellindex, d);
    }
    return _cellcoordmax;
  }();

  {
    for (int d = 0; d < get_ndim(prop_gridtype); d++) {
      // pos_component_vel_relative_to_flow is constant along a ray with a given direction in Cartesian coordinates, but
      // for non-Cartesian coordinates, we still need to check at the current position whether the packet is
      // moving in the positive or negative direction in each grid coordinate direction relative to the homologous grid
      // flow, otherwise we might never enter the cell that we're supposed to be in
      const bool pos_component_vel_relative_to_flow = (pktvelgridcoord[d] * tstart) > pktposgridcoord[d];

      bool isoutside_error = false;
      double delta = 0.;
      if (pos_component_vel_relative_to_flow) {
        // check if packet pos is above cell max while moving in the positive direction relative to the grid flow
        const double boundaryposmax = cellcoordmax[d] / globals::tmin * tstart;
        delta = pktposgridcoord[d] - boundaryposmax;
        isoutside_error = pktposgridcoord[d] > (boundaryposmax + 10.);
      } else {
        // check if packet pos is below cell min while moving in the negative direction relative to the grid flow
        const double boundaryposmin = grid::get_cellcoordmin(cellindex, d) / globals::tmin * tstart;
        delta = pktposgridcoord[d] - boundaryposmin;
        isoutside_error = pktposgridcoord[d] < (boundaryposmin - 10.);
      }

      if (isoutside_error) {
        printout(
            "[ERROR] packet outside coord %d %c%c boundary of cell %d. vel %g initpos %g "
            "cellcoordmin %g, cellcoordmax %g\n",
            d, pos_component_vel_relative_to_flow ? '+' : '-', get_coordlabel(prop_gridtype, d), cellindex,
            pktvelgridcoord[d], pktposgridcoord[d], get_cellcoordmin(cellindex, d) / globals::tmin * tstart,
            cellcoordmax[d] / globals::tmin * tstart);
        printout("globals::tmin %g tstart %g tstart/globals::tmin %g\n", globals::tmin, tstart, tstart / globals::tmin);
        printout(" delta %g\n", delta);

        printout("packet dir [%g, %g, %g]\n", dir[0], dir[1], dir[2]);

        const auto snext = grid::get_cellindex_from_pos(pos, tstart);
        if ((grid::get_cellcoordpointnum(cellindex, d) == (grid::ncoordgrid[d] - 1) &&
             pos_component_vel_relative_to_flow) ||
            (grid::get_cellcoordpointnum(cellindex, d) == 0 && !pos_component_vel_relative_to_flow) || (snext < 0)) {
          printout("[warning] escaping packet\n");
          return {0., -99};
        }
        printout("[warning] swapping packet cellindex from %d to %d, which has cellcoordmin %g, cellcoordmax %g\n",
                 cellindex, snext, get_cellcoordmin(snext, d) / globals::tmin * tstart,
                 get_cellcoordmax(snext, d) / globals::tmin * tstart);
        return {0., snext};
      }
    }
  }

  double distance = std::numeric_limits<double>::max();
  int snext{-1};

  if (prop_gridtype == GridType::SPHERICAL1D) {
    // the only coordinate is the radius from the origin

    const double speed = vec_len(dir) * CLIGHT_PROP;  // just in case dir is not normalised

    const double r_outer = cellcoordmax[0] * tstart / globals::tmin;
    const double d_coordmaxboundary =
        expanding_shell_intersection<BoundaryType::OUTER>(pos, dir, speed, r_outer, tstart);

    // upper d coordinate of the current cell
    if ((d_coordmaxboundary >= 0.) && (d_coordmaxboundary < distance)) {
      distance = d_coordmaxboundary;
      snext = (grid::get_cellcoordpointnum(cellindex, 0) == (grid::ncoordgrid[0] - 1))
                  ? -99
                  : cellindex + grid::get_coordcellindexincrement(0);
    }

    const double r_inner = grid::get_cellcoordmin(cellindex, 0) * tstart / globals::tmin;
    if (r_inner > 0.) {
      const double d_coordminboundary =
          expanding_shell_intersection<BoundaryType::INNER>(pos, dir, speed, r_inner, tstart);
      // lower d coordinate of the current cell
      if ((d_coordminboundary >= 0.) && (d_coordminboundary < distance)) {
        distance = d_coordminboundary;
        snext =
            (grid::get_cellcoordpointnum(cellindex, 0) == 0) ? -99 : cellindex - grid::get_coordcellindexincrement(0);
      }
    }
  } else if (prop_gridtype == GridType::CYLINDRICAL2D) {
    // coordinate 0 is cylindrical radius (distance from z=0 in x-y plane), coord 1 is z

    const std::array<double, 2> posnoz = {pos[0], pos[1]};

    // r_cyl component of direction vector
    const double dirxylen = std::sqrt((dir[0] * dir[0]) + (dir[1] * dir[1]));
    // r_cyl component of velocity
    const double xyspeed = dirxylen * CLIGHT_PROP;

    // make a normalised 2D direction vector in the xy plane
    const std::array<double, 2> dirnoz = {dir[0] / dirxylen, dir[1] / dirxylen};

    const double r_outer = cellcoordmax[0] * tstart / globals::tmin;
    const double d_rcyl_coordmaxboundary =
        expanding_shell_intersection<BoundaryType::OUTER>(posnoz, dirnoz, xyspeed, r_outer, tstart);
    if (d_rcyl_coordmaxboundary >= 0.) {
      // how far did the packet travel in the z direction during this time?
      const double d_z_coordmaxboundary = d_rcyl_coordmaxboundary / xyspeed * dir[2] * CLIGHT_PROP;
      // distance from two perpendicular components to the r_cyl upper boundary
      const double d_coordmaxboundary_rcyl = std::sqrt((d_rcyl_coordmaxboundary * d_rcyl_coordmaxboundary) +
                                                       (d_z_coordmaxboundary * d_z_coordmaxboundary));
      if ((d_coordmaxboundary_rcyl > 0) && (d_coordmaxboundary_rcyl < distance)) {
        distance = d_coordmaxboundary_rcyl;
        snext = (grid::get_cellcoordpointnum(cellindex, 0) == (grid::ncoordgrid[0] - 1))
                    ? -99
                    : cellindex + grid::get_coordcellindexincrement(0);
      }
    }

    const double r_inner = grid::get_cellcoordmin(cellindex, 0) * tstart / globals::tmin;
    // don't try to calculate the intersection if the inner radius is zero
    if (r_inner > 0) {
      // calculate the distance in the xy plane to the inner boundary
      const double d_rcyl_coordminboundary =
          expanding_shell_intersection<BoundaryType::INNER>(posnoz, dirnoz, xyspeed, r_inner, tstart);
      if (d_rcyl_coordminboundary >= 0.) {
        const double d_z_coordminboundary = d_rcyl_coordminboundary / xyspeed * dir[2] * CLIGHT_PROP;
        // distance from two perpendicular components to the r_cyl lower boundary
        const double d_coordminboundary_rcyl = std::sqrt((d_rcyl_coordminboundary * d_rcyl_coordminboundary) +
                                                         (d_z_coordminboundary * d_z_coordminboundary));
        if ((d_coordminboundary_rcyl >= 0.) && (d_coordminboundary_rcyl < distance)) {
          distance = d_coordminboundary_rcyl;
          snext =
              (grid::get_cellcoordpointnum(cellindex, 0) == 0) ? -99 : cellindex - grid::get_coordcellindexincrement(0);
        }
      }
    }

    // handle Z boundaries as Cartesian

    if (pktvelgridcoord[1] > (cellcoordmax[1] / globals::tmin)) {
      const double t_zcoordmaxboundary = ((pktposgridcoord[1] - (pktvelgridcoord[1] * tstart)) /
                                          (cellcoordmax[1] - (pktvelgridcoord[1] * globals::tmin)) * globals::tmin) -
                                         tstart;
      const double d_coordmaxboundary_z = CLIGHT_PROP * t_zcoordmaxboundary;

      if ((d_coordmaxboundary_z >= 0.) && (d_coordmaxboundary_z < distance)) {
        distance = d_coordmaxboundary_z;
        snext = (grid::get_cellcoordpointnum(cellindex, 1) == (grid::ncoordgrid[1] - 1))
                    ? -99
                    : cellindex + grid::get_coordcellindexincrement(1);
      }
    } else if (pktvelgridcoord[1] < (grid::get_cellcoordmin(cellindex, 1) / globals::tmin)) {
      const double t_zcoordminboundary =
          ((pktposgridcoord[1] - (pktvelgridcoord[1] * tstart)) /
           ((grid::get_cellcoordmin(cellindex, 1)) - (pktvelgridcoord[1] * globals::tmin)) * globals::tmin) -
          tstart;
      const double d_coordminboundary_z = CLIGHT_PROP * t_zcoordminboundary;

      if ((d_coordminboundary_z >= 0.) && (d_coordminboundary_z < distance)) {
        distance = d_coordminboundary_z;
        snext =
            (grid::get_cellcoordpointnum(cellindex, 1) == 0) ? -99 : cellindex - grid::get_coordcellindexincrement(1);
      }
    }

  } else if (prop_gridtype == GridType::CARTESIAN3D) {
    // There are six possible boundary crossings. Each of the three
    // cartesian coordinates may be taken in turn. For x, the packet
    // trajectory is
    // x = x0 + (dir.x) * c * (t - tstart)
    // the boundaries follow
    // x+/- = x+/-(tmin) * (t/tmin)
    // so the crossing occurs when
    // t = (x0 - (dir.x)*c*tstart)/(x+/-(tmin)/tmin - (dir.x)c)

    // Modified so that it also returns the distance to the closest cell
    // boundary, regardless of direction.

    for (int d = 0; d < 3; d++) {
      if (pktvelgridcoord[d] > (cellcoordmax[d] / globals::tmin)) {
        const double t_coordmaxboundary = ((pktposgridcoord[d] - (pktvelgridcoord[d] * tstart)) /
                                           (cellcoordmax[d] - (pktvelgridcoord[d] * globals::tmin)) * globals::tmin) -
                                          tstart;

        const double d_coordmaxboundary = CLIGHT_PROP * t_coordmaxboundary;

        if ((d_coordmaxboundary >= 0.) && (d_coordmaxboundary < distance)) {
          distance = d_coordmaxboundary;
          snext = (grid::get_cellcoordpointnum(cellindex, d) == (grid::ncoordgrid[d] - 1))
                      ? -99
                      : cellindex + grid::get_coordcellindexincrement(d);
        }
      } else if (pktvelgridcoord[d] < (grid::get_cellcoordmin(cellindex, d) / globals::tmin)) {
        const double t_coordminboundary =
            ((pktposgridcoord[d] - (pktvelgridcoord[d] * tstart)) /
             (grid::get_cellcoordmin(cellindex, d) - (pktvelgridcoord[d] * globals::tmin)) * globals::tmin) -
            tstart;

        const double d_coordminboundary = CLIGHT_PROP * t_coordminboundary;

        // lower d coordinate of the current cell
        if ((d_coordminboundary >= 0.) && (d_coordminboundary < distance)) {
          distance = d_coordminboundary;
          snext =
              (grid::get_cellcoordpointnum(cellindex, d) == 0) ? -99 : cellindex - grid::get_coordcellindexincrement(d);
        }
      }
    }
  } else {
    assert_always(false);
  }

  if constexpr (TESTMODE) {
    if (snext == -1) {
      printout("Something wrong in boundary crossing - didn't find anything.\n");
      printout("packet cell %d\n", cellindex);
      printout("globals::tmin %g tstart %g\n", globals::tmin, tstart);
      for (int d2 = 0; d2 < 3; d2++) {
        printout("coord %d: initpos %g dir %g\n", d2, pos[d2], dir[d2]);
      }
      printout("|initpos| %g |dir| %g |pos.dir| %g\n", vec_len(pos), vec_len(dir), dot(pos, dir));
      for (int d2 = 0; d2 < get_ndim(prop_gridtype); d2++) {
        printout("coord %d: cellcoordmin %g cellcoordmax %g\n", d2,
                 grid::get_cellcoordmin(cellindex, d2) * tstart / globals::tmin,
                 cellcoordmax[d2] * tstart / globals::tmin);
      }
      printout("tstart %g\n", tstart);
    }
  }

  assert_always((snext == -99) || ((snext >= 0) && (snext < grid::ngrid)));

  const double maxsdist =
      (prop_gridtype == GridType::CARTESIAN3D)
          ? std::numbers::sqrt3 * globals::rmax * (tstart + (distance / CLIGHT_PROP)) / globals::tmin
          : 2 * globals::rmax * (tstart + (distance / CLIGHT_PROP)) / globals::tmin;

  assert_always(distance >= 0.);
  assert_always(distance <= maxsdist);

  if (distance > globals::max_path_step) {
    return {globals::max_path_step, cellindex};
  }

  return {distance, snext};
}

}  // namespace grid
