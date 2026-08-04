// The model grid and propagation grid: reads the input ejecta model (model.txt), sets up cell
// densities, abundances, and temperatures, maps model cells onto the propagation grid (1D
// spherical, 2D cylindrical, or 3D Cartesian), and provides the cell geometry and
// boundary-crossing calculations used during packet propagation.

#include "grid.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <print>
#include <set>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "input.h"
#include "kpkt.h"
#include "mpi_logging.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "radfield.h"
#include "random.h"
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

int first_input_cellid{-1};  // auto-determine first cell index in model.txt (usually 1 or 0)

// Initial co-ordinates of inner most corner of cell.
std::array<std::vector<double>, 3> coord_pos_min_tmin{};

// associate each propagation cell with a model grid cell, or not, if the cell is empty (or doesn't get mapped to
// anything such as 1D/2D to 3D)
std::vector<int> propcell_mgi;
std::vector<int> propcell_nonemptymgi;

std::vector<int> modelgrid_numpropcells;
std::vector<int> nonemptymgi_of_mgi;
std::vector<int> mgi_of_nonemptymgi;

std::vector<double> totmassnuclide{};  // total mass of each nuclide in the ejecta

MPI_shared_array<float> initnucmassfrac_allcells{};
MPI_shared_array<float> initmassfracuntrackedstable_allcells{};
MPI_shared_array<int> elements_uppermost_ion_allcells{};  // Highest ion index that has a significant population

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
MPI_shared_array<ModelGridCellInput> modelgrid_input{};

enum class BoundaryType : std::uint8_t { LOWER, UPPER };

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
      return std::array<char, 3>{'x', 'y', 'z'}.at(axis);
    case GridType::CYLINDRICAL2D:
      return std::array<char, 2>{'r', 'z'}.at(axis);
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
    printlnlog("Ye.txt not present, so no initial electron fractions will be applied from it");
    return;
  }

  // the electron fractions are written to node-shared memory, so only the node leaders read the file
  // (synchronised by the barrier below)
  if (globals::rank_in_node == 0) {
    const auto filein = fopen_required_uniqueptr("Ye.txt", "r");
    int nlines_in = 0;
    assert_always(fscanf(filein.get(), "%d", &nlines_in) == 1);

    int cells_set = 0;
    int entries_ignored = 0;
    for (int n = 0; n < nlines_in; n++) {
      int mgiplusone = -1;
      float initelecfrac = 0.;
      assert_always(fscanf(filein.get(), "%d %g", &mgiplusone, &initelecfrac) == 2);
      const int mgi = mgiplusone - 1;
      if (mgi >= 0 && mgi < get_npts_model()) {
        set_initelectronfrac(mgi, initelecfrac);
        cells_set++;
      } else {
        entries_ignored++;
      }
    }
    printlnlog("Ye.txt: set the initial electron fraction for {} of {} model cells ({} out-of-range entries ignored)",
               cells_set, get_npts_model(), entries_ignored);
  }
  MPI_Barrier_allranks();
}

[[gnu::pure]] DEVICE_FUNC auto get_propgridtype() -> GridType {
  if constexpr (GRID_TYPE_OVERRIDE.has_value()) {
    return GRID_TYPE_OVERRIDE.value();
  }
  return get_modelgridtype();
}

// how much do we change the cellindex to move along a coordinately axis (e.g., the x, y, z directions, or r direction)
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_coordcellindexstride(const int axis) -> int {
  int stride = 1;
  for (int a = 0; a < axis; ++a) {
    stride *= ncoordgrid[a];
  }
  return stride;
}

// convert a cell index number into an integer (x,y,z or r) coordinate index from 0 to ncoordgrid[axis]
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordindex(const int cellindex, const int axis) -> int {
  return (cellindex / get_coordcellindexstride(axis)) % ncoordgrid[axis];
}

// get the minimum value of a coordinate at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the minimum x position in xyz coords, or the minimum radius
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordmin(const int cellindex, const int axis) -> double {
  return coord_pos_min_tmin[axis][get_cellcoordindex(cellindex, axis)];
}

// get the maximum position value of a coordinate axis at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the maximum x position in xyz coords, or the maximum radius in spherical 1D
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordmax(const int cellindex, const int axis) -> double {
  const auto idx = get_cellcoordindex(cellindex, axis);
  const auto idxlast = ncoordgrid[axis] - 1;
  return idx < idxlast ? coord_pos_min_tmin[axis][idx + 1] : globals::rmax;
}

// return the inner radius of a propagation cell at time tmin
auto get_cell_r_inner(const int cellindex, const GridType prop_gridtype) -> double {
  if (prop_gridtype == GridType::SPHERICAL1D) {
    return get_cellcoordmin(cellindex, 0);
  }

  // minimum |coordinate| within the cell along one axis: zero if the cell straddles the coordinate plane
  const auto axis_mindist = [](const double coordmin, const double coordmax) {
    return (coordmin <= 0. && coordmax >= 0.) ? 0. : std::min(std::abs(coordmin), std::abs(coordmax));
  };

  if (prop_gridtype == GridType::CYLINDRICAL2D) {
    const auto rcyl_inner = get_cellcoordmin(cellindex, 0);
    const auto z_inner = axis_mindist(get_cellcoordmin(cellindex, 1), get_cellcoordmax(cellindex, 1));
    return std::sqrt(pow2(rcyl_inner) + pow2(z_inner));
  }

  if (prop_gridtype == GridType::CARTESIAN3D) {
    const auto x_inner = axis_mindist(get_cellcoordmin(cellindex, 0), get_cellcoordmax(cellindex, 0));
    const auto y_inner = axis_mindist(get_cellcoordmin(cellindex, 1), get_cellcoordmax(cellindex, 1));
    const auto z_inner = axis_mindist(get_cellcoordmin(cellindex, 2), get_cellcoordmax(cellindex, 2));
    return std::sqrt(pow2(x_inner) + pow2(y_inner) + pow2(z_inner));
  }

  assert_always(false);
  return NAN;
}

// Negative input mass fractions (within roundoff of zero) are counted as they are clamped during
// the model read and summarised once at the end of read_ejecta_model() rather than being logged
// per cell and nuclide, which could produce millions of lines for hydro-derived models.
int n_negative_input_massfrac = 0;
double most_negative_input_massfrac = 0.;
int first_negative_massfrac_mgi = -1;
int first_negative_massfrac_z = -1;  // -1 means the Fe-group mass fraction column
int first_negative_massfrac_a = -1;

// names of model.txt columns that were not recognised, reported after the model read
std::set<std::string> ignored_model_columns;

void count_negative_input_massfrac(const int modelgridindex, const int z, const int a, const double massfrac) {
  n_negative_input_massfrac++;
  most_negative_input_massfrac = std::min(most_negative_input_massfrac, massfrac);
  if (first_negative_massfrac_mgi < 0) {
    first_negative_massfrac_mgi = modelgridindex;
    first_negative_massfrac_z = z;
    first_negative_massfrac_a = a;
  }
}

void report_negative_input_massfracs() {
  if (n_negative_input_massfrac > 0) {
    printlog(
        "[warning] {} negative mass fraction(s) in the input model were clamped to zero (most negative {:g}; first "
        "in cell {} for ",
        n_negative_input_massfrac, most_negative_input_massfrac, first_negative_massfrac_mgi);
    if (first_negative_massfrac_z < 0) {
      printlnlog("X_Fegroup)");
    } else {
      printlnlog("Z={} A={})", first_negative_massfrac_z, first_negative_massfrac_a);
    }
  }
}

void set_ffegrp(const int modelgridindex, float x) {
  if (!(x >= 0.)) {
    if (!(x > -1e-6)) {
      printlnlog("[error] Fe-group mass fraction {:g} is negative in cell {}", x, modelgridindex);
    }
    assert_always(x > -1e-6);
    count_negative_input_massfrac(modelgridindex, -1, -1, x);
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
    if (!(abund > -1e-6)) {
      printlnlog("[error] mass fraction {:g} for Z={} A={} is negative in cell {}", abund, decay::get_nuc_z(nucindex),
                 decay::get_nuc_a(nucindex), modelgridindex);
    }
    assert_always(abund > -1e-6);
    count_negative_input_massfrac(modelgridindex, decay::get_nuc_z(nucindex), decay::get_nuc_a(nucindex), abund);
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

void set_elem_untrackedstable_massfrac(const int nonemptymgi, const int element, const float elem_massfrac) {
  // set the stable mass fraction of an element from the total element mass fraction
  // by subtracting the abundances of radioactive isotopes.
  // if the element Z=anumber has no specific stable abundance variable then the function does nothing

  const int atomic_number = get_atomicnumber(element);
  const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);

  double massfrac_allisotopes = 0.;
  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    if (decay::get_nuc_z(nucindex) == atomic_number) {
      // radioactive isotope of this element
      massfrac_allisotopes += get_modelinitnucmassfrac(mgi, nucindex);
    }
  }

  double massfrac_untrackedstable = elem_massfrac - massfrac_allisotopes;

  if (massfrac_untrackedstable < 0.) {
    //  allow some roundoff error before we complain
    if (std::abs(massfrac_allisotopes - elem_massfrac) > 1e-4) {
      printlnlog("[warning] cell {} Z={} element massfrac is less than the sum of its radioisotope massfracs", mgi,
                 atomic_number);
      printlnlog("  massfrac(Z) {:g} massfrac_radioisotopes(Z) {:g}", elem_massfrac, massfrac_allisotopes);
      printlnlog("  increasing elemental massfrac to {:g} and setting stable isotopic massfrac to zero",
                 massfrac_allisotopes);
    }
    // result is allowed to be slightly negative due to roundoff error
    assert_always(massfrac_untrackedstable >= -1e-2);
    massfrac_untrackedstable = 0.;  // bring up to zero if negative
  }

  initmassfracuntrackedstable_allcells[(nonemptymgi * get_nelements()) + element] =
      static_cast<float>(massfrac_untrackedstable);

  // (massfrac_allisotopes + massfrac_untrackedstable) might not exactly match elem_massfrac if we had to boost it to
  // reach massfrac_allisotopes
  set_elem_massfrac(nonemptymgi, element, static_cast<float>(massfrac_allisotopes + massfrac_untrackedstable));
}

// get the radial distance from the origin to the centre of the cell at time tmin
auto get_cellradialposmid(const int cellindex) -> double {
  const auto prop_gridtype = get_propgridtype();
  if (prop_gridtype == GridType::SPHERICAL1D) {
    // volume averaged mean radius is slightly complex for radial shells
    const double r_inner = get_cellcoordmin(cellindex, 0);
    const double r_outer = get_cellcoordmax(cellindex, 0);
    return 3. / 4 * (pow4(r_outer) - pow4(r_inner)) / (pow3(r_outer) - pow3(r_inner));
  }

  if (prop_gridtype == GridType::CYLINDRICAL2D) {
    const double rcyl_mid = (get_cellcoordmin(cellindex, 0) + get_cellcoordmax(cellindex, 0)) / 2;
    const double z_mid = (get_cellcoordmin(cellindex, 1) + get_cellcoordmax(cellindex, 1)) / 2;
    return std::sqrt(pow2(rcyl_mid) + pow2(z_mid));
  }

  // cubic grid requires taking the length of the 3D position vector
  double lensquared = 0.;
  for (int axis = 0; axis < 3; axis++) {
    lensquared += pow2((get_cellcoordmin(cellindex, axis) + get_cellcoordmax(cellindex, axis)) / 2);
  }

  return std::sqrt(lensquared);
}

void allocate_nonemptycells_composition_cooling() {
  // Initialise composition dependent cell data for the given cell
  const ptrdiff_t nonempty_npts_model_ptrdifft = get_nonempty_npts_model();
  const auto nelements = get_nelements();

  initmassfracuntrackedstable_allcells = MPI_shared_array<float>(nonempty_npts_model_ptrdifft * nelements, 0.);
  elem_meanweight_allcells = MPI_shared_array<float>(nonempty_npts_model_ptrdifft * nelements, 0.);
  elements_uppermost_ion_allcells = MPI_shared_array<int>(nonempty_npts_model_ptrdifft * nelements, -1);
  elem_massfracs_allcells = MPI_shared_array<float>(nonempty_npts_model_ptrdifft * nelements, 0.);
  ion_groundlevelpops_allcells = MPI_shared_array<float>(nonempty_npts_model_ptrdifft * get_includedions(), 0.);
  ion_partfuncts_allcells = MPI_shared_array<float>(nonempty_npts_model_ptrdifft * get_includedions(), 0.);
  kpkt::ion_cooling_contribs_allcells = MPI_shared_array<double>(nonempty_npts_model_ptrdifft * get_includedions(), 0.);

  // -1 indicates that there is currently no information on the nlte populations
  nltepops_allcells = MPI_shared_array<double>(nonempty_npts_model_ptrdifft * globals::total_nlte_levels, -1.);
}

// build the mapping between model cells and propagation cells, then allocate the per-cell
// shared-memory arrays (populations, estimators, temperatures) for the nonempty model cells only
void allocate_nonemptymodelcells() {
  // Determine the number of simulation cells associated with the model cells
  std::ranges::fill(modelgrid_numpropcells, 0);
  if (globals::rank_in_node == 0) {
    for (int mgi = 0; mgi < get_npts_model(); mgi++) {
      modelgrid_input[mgi].initial_radial_pos_sum = 0.;
    }
  }
  MPI_Barrier_node();

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

  MPI_Barrier_node();
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

  reserve_resize(mgi_of_nonemptymgi, nonempty_npts_model);
  std::ranges::fill(mgi_of_nonemptymgi, -2);

  reserve_resize(propcell_nonemptymgi, ngrid);
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
    const int mgi = get_propcell_modelgridindex(cellindex);
    if (mgi < 0) {
      propcell_nonemptymgi[cellindex] = -1;
    } else {
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
  thick_allcells = MPI_shared_array<int>(nonempty_npts_model, 0);
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

  if constexpr (RPKT_USE_EXPANSION_OPACITIES || VPKT_USE_EXPANSION_OPACITIES ||
                RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value()) {
    allocate_expansionopacities();
  }

  reserve_resize(globals::dep_estimator_gamma, nonempty_npts_model);
  std::ranges::fill(globals::dep_estimator_gamma, 0.);

  reserve_resize(globals::dep_estimator_positron, nonempty_npts_model);
  std::ranges::fill(globals::dep_estimator_positron, 0.);

  reserve_resize(globals::dep_estimator_electron, nonempty_npts_model);
  std::ranges::fill(globals::dep_estimator_electron, 0.);

  reserve_resize(globals::dep_estimator_alpha, nonempty_npts_model);
  std::ranges::fill(globals::dep_estimator_alpha, 0.);

  const auto ionestimcount = nonempty_npts_model * globals::nbfcontinua_ground;
  const auto ionestimsize = ionestimcount * sizeof(double);

  if (ionestimsize > 0) {
    globals::corrphotoionrenorm = MPI_shared_array<double>(ionestimcount, 1.);

    reserve_resize(globals::gammaestimator, ionestimcount);
    std::ranges::fill(globals::gammaestimator, 0.);
#ifdef DO_TITER
    reserve_resize(globals::gammaestimator_save, ionestimcount);
    std::ranges::fill(globals::gammaestimator_save, 0.);
#endif
  } else {
    globals::corrphotoionrenorm.reset();
    globals::gammaestimator.clear();
#ifdef DO_TITER
    globals::gammaestimator_save.clear();
#endif
  }

  if (USE_ION_BFHEATING_ESTIMATORS && ionestimsize > 0) {
    reserve_resize(globals::bfheatingestimator, ionestimcount);
    std::ranges::fill(globals::bfheatingestimator, 0.);
#ifdef DO_TITER
    reserve_resize(globals::bfheatingestimator_save, ionestimcount);
    std::ranges::fill(globals::bfheatingestimator_save, 0.);
#endif
  } else {
    globals::bfheatingestimator.clear();
#ifdef DO_TITER
    globals::bfheatingestimator_save.clear();
#endif
  }

  reserve_resize(globals::ffheatingestimator, nonempty_npts_model);
  std::ranges::fill(globals::ffheatingestimator, 0.);

  reserve_resize(globals::colheatingestimator, DIRECT_COL_HEAT ? 0 : nonempty_npts_model);
  std::ranges::fill(globals::colheatingestimator, 0.);

#ifdef DO_TITER
  reserve_resize(globals::ffheatingestimator_save, nonempty_npts_model);
  std::ranges::fill(globals::ffheatingestimator_save, 0.);

  reserve_resize(globals::colheatingestimator_save, DIRECT_COL_HEAT ? 0 : nonempty_npts_model);
  std::ranges::fill(globals::colheatingestimator_save, 0.);
#endif

  // barrier to make sure node master has set abundance values to node shared memory
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

void read_elem_abundances() {
  // barrier to make sure node master has set values in node shared memory
  MPI_Barrier_allranks();
  printlnlog("reading abundances.txt...");
  const bool threedimensional = (get_modelgridtype() == GridType::CARTESIAN3D);

  // Process through the grid to read in the abundances per cell.
  // The abundance file should only contain information for non-empty
  // cells. Its format must be cellnumber (integer), abundance for
  // element Z=1 (float) up to abundance for element Z=30 (float)
  // i.e. in total one integer and 30 floats.

  // the mass fraction arrays are in node-shared memory, so only the node leader of each node parses the file and
  // writes the values (synchronised by the barrier below). The other ranks would just discard everything they read
  if (globals::rank_in_node == 0) {
    auto abundance_file = fstream_required("abundances.txt", std::ios::in);
    std::string line;

    // loop over propagation cells for 3D models, or modelgrid cells
    for (int mgi = 0; mgi < get_npts_model(); mgi++) {
      assert_always(get_noncommentline(abundance_file, line));
      auto remainder = std::string_view{line};

      int cellnumberinput = -1;
      assert_always(parse_next_token(remainder, cellnumberinput));
      assert_always(cellnumberinput == mgi + first_input_cellid);

      // the abundances.txt file specifies the elemental mass fractions for each model cell
      // (or proportional to mass frac, e.g. element densities, because they will be normalised anyway)
      // The abundances begin with hydrogen, helium, etc, going as far up the atomic numbers as required
      double normfactor = 0.;
      std::array<float, 150> elem_massfracs_in{};
      double abund_in = 0.;
      for (int elem_z_index = 0; elem_z_index < std::ssize(elem_massfracs_in); elem_z_index++) {
        const int atomic_number = elem_z_index + 1;
        if (!parse_next_token(remainder, abund_in)) {
          // at least one element (hydrogen) should have been specified for nonempty cells
          assert_always(atomic_number > 1 || get_numpropcells(mgi) == 0);
          break;
        }

        if (abund_in < 0. || abund_in < std::numeric_limits<float>::min()) {
          assert_always(abund_in > -1e-6);
          abund_in = 0.;
        }
        elem_massfracs_in[elem_z_index] = static_cast<float>(abund_in);
        normfactor += elem_massfracs_in[elem_z_index];
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
          const auto elemmassfrac = static_cast<float>(elem_massfracs_in[atomic_number - 1] / normfactor);
          assert_always(elemmassfrac >= 0.);

          // radioactive nuclide abundances should have already been set by read_??_model
          set_elem_untrackedstable_massfrac(nonemptymgi, element, elemmassfrac);
        }
      }
    }
  }

  // barrier to make sure node master has set values in node shared memory
  MPI_Barrier_allranks();
  printlnlog("finished reading abundances.txt");
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

auto get_token_count(std::string const& line) -> int {
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

void read_model_radioabundances(std::istream& fmodel, std::string_view& remainder, const int mgi, const bool keepcell,
                                const std::vector<std::string>& colnames, const std::vector<int>& nucindexlist,
                                const bool one_line_per_cell) {
  if (!one_line_per_cell) {
    static std::string line;
    assert_always(std::getline(fmodel, line));
    remainder = std::string_view{line};
  }

  if (!keepcell) {
    return;
  }

  for (auto i = 0Z; i < std::ssize(colnames); i++) {
    double valuein = 0.;
    assert_always(parse_next_token(remainder, valuein));  // usually a mass fraction, but now can be anything

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
      // reported once after the model read (checking only mgi == 0 could miss it entirely if cell 0 is empty)
      ignored_model_columns.insert(colnames[i]);
    }
  }
  double valuein = 0.;
  assert_always(!parse_next_token(remainder, valuein));  // should be no tokens left!
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
        headerline = "#inputcellid vel_r_max_kmps logrho";
        break;
      case GridType::CYLINDRICAL2D:
        headerline = "#inputcellid pos_rcyl_mid pos_z_mid rho";
        break;
      case GridType::CARTESIAN3D:
        headerline = "#inputcellid pos_x_min pos_y_min pos_z_min rho";
        break;
    }
    headerline += " X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48";
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
    printlnlog("model.txt has a header line.");
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

  initnucmassfrac_allcells = MPI_shared_array<float>((npts_model + 1) * num_nuclides, 0.);
  printlnlog(
      "[info] mem_usage: input abundance data for {} nuclides for {} cells occupies {:.3f} MB (node shared memory)",
      num_nuclides, npts_model, (initnucmassfrac_allcells.size() * sizeof(float)) / 1024. / 1024.);

  return {colnames, nucindexlist, one_line_per_cell};
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
      const double delta_z = 2. * globals::vmax * t_model / ncoord_model[1];
      return pow3(globals::tmin / t_model) * delta_z * PI * (pow2((n_r + 1) * delta_rcyl) - pow2(n_r * delta_rcyl));
    }

    case GridType::CARTESIAN3D: {
      // Assumes cells are cubes here - all same volume.
      return pow3(2 * globals::vmax * globals::tmin) / (ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2]);
    }
  }

  assert_always(false);
  return NAN;
}

void calc_modelinit_totmassnuclides() {
  mtot_input = 0.;
  mfegroup = 0.;

  reserve_resize(totmassnuclide, decay::get_num_nuclides());
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

  printlnlog("reading grid restart snapshot from {}", filename);
  FILE* gridsave_file = fopen_required(filename, "r");

  int ntimesteps_in = -1;
  assert_always(fscanf(gridsave_file, "%d ", &ntimesteps_in) == 1);
  assert_always(ntimesteps_in == globals::ntimesteps);

  int nprocs_in = -1;
  assert_always(fscanf(gridsave_file, "%d ", &nprocs_in) == 1);
  assert_always(nprocs_in == globals::nprocs);

  for (int nts = 0; nts < globals::ntimesteps; nts++) {
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
               &globals::timesteps[nts].gamma_emission, &globals::timesteps[nts].pellet_decays) == 22);
  }

  int timestep_in = 0;
  assert_always(fscanf(gridsave_file, "%d ", &timestep_in) == 1);
  assert_always(timestep_in == timestep);

  for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
    int mgi_in = -1;
    float T_R = 0.;
    float T_e = 0.;
    float W = 0.;
    float T_J = 0.;
    int thick = 0;

    float nne_in = -1.;
    float nnetot_in = -1.;
    assert_always(fscanf(gridsave_file, "%d %a %a %a %a %d %la %la %la %la %a %a", &mgi_in, &T_R, &T_e, &W, &T_J,
                         &thick, &globals::dep_estimator_gamma[nonemptymgi],
                         &globals::dep_estimator_positron[nonemptymgi], &globals::dep_estimator_electron[nonemptymgi],
                         &globals::dep_estimator_alpha[nonemptymgi], &nne_in, &nnetot_in) == 12);

    if (mgi_in != mgi) {
      printlnlog("[error] read_grid_restart_data: cell mismatch in {}: read cellnumber {}, expected {}. aborting",
                 filename, mgi_in, mgi);
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

    if (globals::rank_in_node == 0) {
      // node-shared arrays are written by the node master only (all ranks read identical values from the file)
      TR_allcells[nonemptymgi] = T_R;
      Te_allcells[nonemptymgi] = T_e;
      W_allcells[nonemptymgi] = W;
      TJ_allcells[nonemptymgi] = T_J;
      thick_allcells[nonemptymgi] = thick;
      nne_allcells[nonemptymgi] = nne_in;
      nnetot_allcells[nonemptymgi] = nnetot_in;
    }

    if constexpr (USE_LUT_PHOTOION) {
      for (int i = 0; i < globals::nbfcontinua_ground; i++) {
        const ptrdiff_t estimindex = (static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) + i;
        double corrphotoionrenorm_in = 0.;
        assert_always(fscanf(gridsave_file, " %la %la", &corrphotoionrenorm_in, &globals::gammaestimator[estimindex]) ==
                      2);
        if (globals::rank_in_node == 0) {
          // corrphotoionrenorm is node-shared (gammaestimator is per-rank, so every rank reads into it)
          globals::corrphotoionrenorm[estimindex] = corrphotoionrenorm_in;
        }
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
  // according to the local energy density resulting from the 56Ni decay.
  // The dilution factor is W=1 in LTE.

  printlog("Assigning initial temperatures...");

  const double tstart = globals::timesteps[0].mid;
  int cells_below_mintemp = 0;
  int cells_above_maxtemp = 0;
  int cells_nonfinite_temp = 0;
  int first_nonfinite_mgi = -1;
  double first_nonfinite_rho_tmin = 0.;
  double first_nonfinite_endecay = 0.;

  for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);

    const auto q = (INITIAL_PACKETS_ON && USE_MODEL_INITIAL_ENERGY) ? get_initenergyq(mgi) : 0.;
    const double decayedenergy_per_mass =
        decay::get_endecay_per_ejectamass_tmodel_to_time_withexpansion(nonemptymgi, tstart) + q;

    auto T_initial = static_cast<float>(std::pow(
        CLIGHT / 4 / STEBO * pow3(globals::tmin / tstart) * get_rho_tmin(mgi) * decayedenergy_per_mass, 1. / 4.));

    if (!std::isfinite(T_initial)) {
      // check this first: a NaN would fall through every comparison below and be stored unclamped
      cells_nonfinite_temp++;
      if (first_nonfinite_mgi < 0) {
        first_nonfinite_mgi = mgi;
        first_nonfinite_rho_tmin = get_rho_tmin(mgi);
        first_nonfinite_endecay = decayedenergy_per_mass;
      }
      T_initial = MINTEMP;
    } else if (T_initial < MINTEMP) {
      T_initial = MINTEMP;
      cells_below_mintemp++;
    } else if (T_initial > MAXTEMP) {
      T_initial = MAXTEMP;
      cells_above_maxtemp++;
    }

    if (globals::rank_in_node == 0) {
      // set the initial temperatures in the modelgrid
      // this is only done by the node master, so that the values are shared in the node shared memory
      Te_allcells[nonemptymgi] = T_initial;
      TJ_allcells[nonemptymgi] = T_initial;
      TR_allcells[nonemptymgi] = T_initial;
      W_allcells[nonemptymgi] = 1.;
      thick_allcells[nonemptymgi] = 0;
    }
  }
  printlnlog("  cells below MINTEMP {:g} [K]: {}. Above MAXTEMP {:g} [K]: {}", MINTEMP, cells_below_mintemp, MAXTEMP,
             cells_above_maxtemp);
  if (cells_nonfinite_temp > 0) {
    printlnlog(
        "[warning] {} cells had a non-finite initial temperature and were set to MINTEMP (first was mgi {} with "
        "rho_tmin {:g} [g/cm3] and decayed energy {:g} [erg/g])",
        cells_nonfinite_temp, first_nonfinite_mgi, first_nonfinite_rho_tmin, first_nonfinite_endecay);
  }
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
  // vmax is per coordinate, but the simulation volume corners will have a higher expansion velocity than the sides
  const double vmax_corner = sqrt(3 * pow2(globals::vmax));
  printlnlog("corner vmax {:g} [cm/s] ({:.2f}c)", vmax_corner, vmax_corner / CLIGHT);
  if (!FORCE_SPHERICAL_ESCAPE_SURFACE) {
    assert_always(vmax_corner < CLIGHT);
  }

  // Set grid size for uniform xyz grid
  if (get_modelgridtype() == GridType::CARTESIAN3D) {
    // if we used in a 3D ejecta model, the propagation grid will match the input grid exactly
    assert_always(ncoordgrid[0] > 0);
    assert_always(ncoordgrid[1] > 0);
    assert_always(ncoordgrid[2] > 0);
  } else {
    ncoordgrid = {CUBOID_NCOORDGRID_X, CUBOID_NCOORDGRID_Y, CUBOID_NCOORDGRID_Z};
  }

  // artis assumes in some places that the cells are cubes, not cuboids
  assert_always(ncoordgrid[0] == ncoordgrid[1]);
  assert_always(ncoordgrid[0] == ncoordgrid[2]);

  ngrid = static_cast<ptrdiff_t>(ncoordgrid[0]) * ncoordgrid[1] * ncoordgrid[2];
  reserve_resize(coord_pos_min_tmin[0], ncoordgrid[0]);
  reserve_resize(coord_pos_min_tmin[1], ncoordgrid[1]);
  reserve_resize(coord_pos_min_tmin[2], ncoordgrid[2]);

  for (int axis = 0; axis < 3; axis++) {
    for (int i = 0; i < ncoordgrid[axis]; i++) {
      coord_pos_min_tmin[axis][i] = -globals::rmax + (2 * i * globals::rmax / ncoordgrid[axis]);
    }
  }
}

void setup_grid_spherical_1d() {
  assert_always(get_modelgridtype() == GridType::SPHERICAL1D);

  ncoordgrid = {get_npts_model(), 1, 1};

  ngrid = static_cast<ptrdiff_t>(ncoordgrid[0]) * ncoordgrid[1] * ncoordgrid[2];

  reserve_resize(coord_pos_min_tmin[0], ncoordgrid[0]);

  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const int mgi = cellindex;  // interchangeable in this mode
    const double v_inner = mgi > 0 ? vout_model[mgi - 1] : 0.;
    coord_pos_min_tmin[0][cellindex] = v_inner * globals::tmin;
  }
}

void setup_grid_cylindrical_2d() {
  const double vmax_corner = sqrt(2 * pow2(globals::vmax));
  printlnlog("corner vmax {:g} [cm/s] ({:.2f}c)", vmax_corner, vmax_corner / CLIGHT);
  assert_always(vmax_corner < CLIGHT);

  assert_always(get_modelgridtype() == GridType::CYLINDRICAL2D);

  ncoordgrid = ncoord_model;

  ngrid = ncoordgrid[0] * ncoordgrid[1];
  assert_always(ngrid == get_npts_model());

  reserve_resize(coord_pos_min_tmin[0], ncoordgrid[0]);
  for (int n_rcyl = 0; n_rcyl < ncoordgrid[0]; n_rcyl++) {
    coord_pos_min_tmin[0][n_rcyl] = n_rcyl * globals::rmax / ncoord_model[0];
  }

  reserve_resize(coord_pos_min_tmin[1], ncoordgrid[1]);
  for (int n_z = 0; n_z < ncoordgrid[1]; n_z++) {
    coord_pos_min_tmin[1][n_z] = globals::rmax * (-1 + (n_z * 2. / ncoord_model[1]));
  }
}

constexpr auto get_grid_type_name(const GridType gridtype) -> std::string_view {
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
  switch (get_propgridtype()) {
    case GridType::CARTESIAN3D: {
      const auto idx = static_cast<int>(((pos / time) + globals::vmax) / 2 / globals::vmax * ncoordgrid[axis]);
      assert_always(idx >= 0);
      assert_always(idx < ncoordgrid[axis]);
      return idx;
    }

    case GridType::CYLINDRICAL2D: {
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
      break;
    }

    case GridType::SPHERICAL1D: {
      // radial spacing is non-uniform, so we have to do a search
      const auto trat = time / globals::tmin;
      for (int n_r = 0; n_r < ncoordgrid[0]; n_r++) {
        if ((pos < get_cellcoordmax(n_r, 0) * trat) && (pos >= get_cellcoordmin(n_r, 0) * trat)) {
          return n_r;
        }
      }
      break;
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
    assert_testmodeonly(boundarytype == BoundaryType::LOWER);
    assert_testmodeonly(shellradiuststart < vec_len(pos));
    return -1;
  }

  if (discriminant > 0) {
    // two intersections
    double dist1 = (-b + sqrt(discriminant)) / 2 / a;
    double dist2 = (-b - sqrt(discriminant)) / 2 / a;

    const auto [posfinal1, posfinal2] = [&] {
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
    if constexpr (boundarytype == BoundaryType::LOWER) {
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

// get element mean weight in grams
auto get_element_meanweight(const std::ptrdiff_t nonemptymgi, const int element) -> float {
  if constexpr (USE_CALCULATED_MEANATOMICWEIGHT) {
    const auto mu = elem_meanweight_allcells[(nonemptymgi * get_nelements()) + element];
    assert_always(mu > 0);
    return mu;
  }
  return globals::elements[element].initstablemeannucmass;
}

[[nodiscard]] constexpr auto distance_cartesian_boundary(const double pktposgridcoord, const double pktvelgridcoord,
                                                         const double cellboundarypos, const double tstart) -> double {
  // numerically stable formulation: compute time difference directly as
  // (pos - boundary_at_tstart) / (boundary_velocity - packet_velocity)
  // to avoid catastrophic cancellation when t_crossing ≈ tstart
  return CLIGHT_PROP * (pktposgridcoord - (cellboundarypos / globals::tmin * tstart)) /
         ((cellboundarypos / globals::tmin) - pktvelgridcoord);
}

// maximum out-of-cell position error (in cm) attributable to floating-point rounding of the
// packet position updates. The tolerance scales with the coordinate magnitude because
// rounding errors are proportional to ulp(pos), e.g. about 8 cm at 4e16 cm.
[[nodiscard]] constexpr auto cellbound_tolerance(const double boundarypos) -> double {
  return std::max(10., std::abs(boundarypos) * 1e-12);
}

// Check if the packet is at or (within rounding-error tolerance) past a cell boundary that it
// is moving towards relative to the boundary's homologous motion. In that case the
// intersection calculation finds no forward crossing (a slightly negative distance, or no
// intersection at all for the curved boundaries), and the crossing must instead be treated as
// immediate (zero distance). If it were silently dropped, the crossing could never be
// detected again (the packet is outrunning the boundary), and the out-of-cell position error
// would grow without bound along the rest of the path until the outside-boundary check fails.
template <BoundaryType boundarytype>
[[nodiscard]] constexpr auto is_boundary_overshoot_within_tolerance(const double pktposgridcoord,
                                                                    const double pktvelgridcoord,
                                                                    const double boundarypos_tmin, const double tstart)
    -> bool {
  const double boundaryvel = boundarypos_tmin / globals::tmin;
  const double boundarypos = boundaryvel * tstart;
  const double overshoot =
      boundarytype == BoundaryType::UPPER ? (pktposgridcoord - boundarypos) : (boundarypos - pktposgridcoord);
  const bool movingtowards =
      boundarytype == BoundaryType::UPPER ? (pktvelgridcoord > boundaryvel) : (pktvelgridcoord < boundaryvel);
  return movingtowards && (overshoot >= 0.) && (overshoot <= cellbound_tolerance(boundarypos));
}

}  // anonymous namespace

// for a uniform grid get the extent along the x,y,z coordinate (x_2 - x_1, etc.) at time tmin
// for spherical grid get the radial extent (r_outer - r_inner) at time tmin
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto propcell_width_tmin(const int cellindex, const int axis) -> double {
  return get_cellcoordmax(cellindex, axis) - get_cellcoordmin(cellindex, axis);
}

// return the model cell volume (when mapped to the propagation cells) at globals::tmin
// for a uniform cubic grid this is constant
[[gnu::pure]] [[nodiscard]] auto get_modelcell_assocvolume_tmin(const int modelgridindex) -> double {
  if (get_propgridtype() == GridType::CARTESIAN3D) {
    // underlying input model could be 1D, 2D, or 3D
    return get_propcell_volume_tmin(0) * get_numpropcells(modelgridindex);
  }

  // direct mapping, so cellindex = modelgridindex and one propagation cell per model cell
  return get_propcell_volume_tmin(modelgridindex);
}

// return the propagation cell volume at globals::tmin
// for a spherical grid, the cell index is required (and should be equivalent to a modelgridindex)
[[gnu::pure]] [[nodiscard]] auto get_propcell_volume_tmin(const int cellindex) -> double {
  switch (get_propgridtype()) {
    case GridType::CARTESIAN3D: {
      // volume is constant for uniform grid, so any cell index will do
      const double deltax = get_cellcoordmax(0, 0) - get_cellcoordmin(0, 0);
      // deltax = deltay = deltaz for uniform grid
      return pow3(deltax);
    }

    case GridType::CYLINDRICAL2D: {
      // use cellindex = 0 because spacing is uniform
      const auto delta_z = get_cellcoordmax(0, 1) - get_cellcoordmin(0, 1);
      // delta_rcyl is uniform, but the area of a cylindrical shell is proportional to r^2
      const auto delta_rcylsquared = pow2(get_cellcoordmax(cellindex, 0)) - pow2(get_cellcoordmin(cellindex, 0));
      return delta_z * PI * delta_rcylsquared;
    }

    case GridType::SPHERICAL1D: {
      // the cellindex matters here because the radial spacing is non-uniform
      return 4. / 3. * PI * (pow3(get_cellcoordmax(cellindex, 0)) - pow3(get_cellcoordmin(cellindex, 0)));
    }
  }
  assert_always(false);
  return NAN;
}

[[nodiscard]] auto get_propcell_random_xyz_position_tmin(int cellindex, rngstate_type& rngstate) -> Vec3d {
  switch (get_propgridtype()) {
    case GridType::SPHERICAL1D: {
      const double r_inner = get_cellcoordmin(cellindex, 0);
      const double r_outer = get_cellcoordmax(cellindex, 0);
      // use equal volume probability distribution to select radius
      // (rng_uniform_pos avoids a zero draw placing the packet exactly on the exclusive outer boundary)
      const double radius = std::cbrt(std::lerp(pow3(r_outer), pow3(r_inner), rng_uniform_pos(rngstate)));

      return vec_scale(get_rand_isotropic_unitvec(rngstate), radius);
    }

    case GridType::CYLINDRICAL2D: {
      const double rcyl_inner = get_cellcoordmin(cellindex, 0);
      const double rcyl_outer = get_cellcoordmax(cellindex, 0);
      // use equal area probability distribution to select radius
      const double rcyl_rand = std::sqrt(std::lerp(pow2(rcyl_outer), pow2(rcyl_inner), rng_uniform_pos(rngstate)));
      const double theta_rand = rng_uniform(rngstate) * 2 * PI;
      const double z_rand =
          std::lerp(get_cellcoordmin(cellindex, 1), get_cellcoordmax(cellindex, 1), rng_uniform_pos(rngstate));
      return {std::cos(theta_rand) * rcyl_rand, std::sin(theta_rand) * rcyl_rand, z_rand};
    }

    case GridType::CARTESIAN3D: {
      Vec3d pos;
      for (int axis = 0; axis < 3; axis++) {
        pos[axis] =
            get_cellcoordmin(cellindex, axis) + (rng_uniform_pos(rngstate) * propcell_width_tmin(cellindex, axis));
      }
      return pos;
    }
  }
  assert_always(false);
  return {NAN, NAN, NAN};
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

DEVICE_FUNC auto get_modelgridtype() -> GridType {
  assert_testmodeonly(model_type.has_value());
  return model_type.value();
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
  assert_testmodeonly(modelgridindex <= get_npts_model());
  return modelgrid_numpropcells[modelgridindex];
}

// get the index in the list of non-empty cells for a given model grid cell
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_nonemptymgi_of_mgi(const int mgi) -> int {
  assert_testmodeonly(get_nonempty_npts_model() > 0);
  assert_testmodeonly(mgi >= 0);
  assert_testmodeonly(mgi < get_npts_model());

  const int nonemptymgi = nonemptymgi_of_mgi[mgi];
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());

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
    nucleondens += get_elem_numberdens(nonemptymgi, element) * get_element_meanweight(nonemptymgi, element) / MH;
  }
  return get_nnetot(nonemptymgi) / nucleondens;
}

// q: energy in the model at tmin per gram to use with USE_MODEL_INITIAL_ENERGY option [erg/g]
DEVICE_FUNC auto get_initenergyq(const int modelgridindex) -> double {
  return modelgrid_input[modelgridindex].initenergyq;
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
      // electron-fraction-dependent opacities from Tanaka et al. (2020) table 1.
      const auto Ye = modelgrid_input[mgi].initelectronfrac;
      assert_always(Ye > 0.);

      // pairs of (upper Ye limit, kappa [cm^2/g])
      constexpr auto kappa_table = std::to_array<std::pair<double, double>>(
          {{0.1, 19.5}, {0.15, 32.2}, {0.20, 22.3}, {0.25, 5.6}, {0.30, 5.36}, {0.35, 3.3}});
      const auto* const entry =
          std::ranges::find_if(kappa_table, [Ye](const auto& ye_kappa) { return Ye <= ye_kappa.first; });
      kappa = (entry != kappa_table.end()) ? entry->second : 0.96;
      break;
    }

    case RpktGreyType::JUST2022_TEMP_LANTHANIDEFRAC: {
      // grey opacity used in Just+2022, https://ui.adsabs.harvard.edu/abs/2022MNRAS.510.2820J/abstract
      // kappa is a simple analytic function of temperature and lanthanide mass fraction
      // adapted to best fit lightcurves from Kasen+2017 in ALCAR simulations
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

// read the input ejecta model from model.txt, auto-detecting its dimensionality (1D spherical,
// 2D cylindrical, or 3D Cartesian) and reading the per-cell densities, abundances, and any
// optional extra columns into the model grid
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
    npts_model = static_cast<ptrdiff_t>(npts_0) * npts_1;
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
  if (!line.starts_with('#')) {
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
  modelgrid_input = MPI_shared_array<ModelGridCellInput>(npts_model + 1, ModelGridCellInput{});
  modelgrid_numpropcells.resize(npts_model + 1, 0);
  nonemptymgi_of_mgi.resize(npts_model + 1, -1);

  // geometry-specific setup of the model coordinate counts
  if (get_modelgridtype() == GridType::SPHERICAL1D) {
    ncoord_model[0] = npts_0;
    ncoord_model[1] = 0;
    ncoord_model[2] = 0;
    vout_model.resize(get_npts_model(), NAN);
  } else if (get_modelgridtype() == GridType::CYLINDRICAL2D) {
    ncoord_model[0] = npts_0;
    ncoord_model[1] = npts_1;
    ncoord_model[2] = 0;
  } else {
    ncoord_model[0] = static_cast<int>(round(std::cbrt(npts_model)));
    ncoord_model[1] = ncoord_model[0];
    ncoord_model[2] = ncoord_model[0];
    assert_always(static_cast<ptrdiff_t>(ncoord_model[0]) * ncoord_model[1] * ncoord_model[2] == npts_model);
    // for a 3D input model, the propagation cells will match the input cells exactly
    ncoordgrid = ncoord_model;
    ngrid = npts_model;
    min_den = -1.;
  }

  const auto [colnames, nucindexlist, one_line_per_cell] = read_model_columns(fmodel);

  // 3D only: whether the input cell positions match the expected values in x-y-z or z-y-x
  // column order (set false when a mismatch is detected)
  bool posmatch_xyz = true;
  bool posmatch_zyx = true;

  int mgi = 0;
  // only global rank 0 parses the cell data. The node-shared arrays it fills are broadcast to the leaders of the
  // other nodes below, and the rank-local results are broadcast to all ranks, so that only one rank reads the
  // (potentially very large) file
  if (globals::my_rank == 0) {
    while (mgi < get_npts_model() && std::getline(fmodel, line)) {
      auto remainder = std::string_view{line};
      int cellnumberin = 0;
      double rho_tmodel{NAN};  // the cell density [g/cm3] at the model snapshot time t_model

      // parse the geometry-specific part of the cell line to get the density (and for 1D, the
      // outer velocity), and check that any cell midpoint coordinates match the expected values
      if (get_modelgridtype() == GridType::SPHERICAL1D) {
        // columns: cell number, outer velocity [km/s], log10 of the density
        double vout_kmps{NAN};
        double log_rho{NAN};
        if (!(parse_next_token(remainder, cellnumberin) && parse_next_token(remainder, vout_kmps) &&
              parse_next_token(remainder, log_rho))) {
          printlnlog(
              "[error] model.txt cell {}: expected at least 3 values (inputcellid vel_r_max_kmps log10rho) but could "
              "not parse line: {}",
              mgi, line);
          assert_always(false);
        }
        vout_model[mgi] = vout_kmps * 1.e5;
        // the velocity grid is binary searched by int_index_lowerbound(), so it has to increase
        // outwards. Without this check a non-monotonic model.txt gives silently wrong cell lookups.
        assert_always(mgi == 0 || vout_model[mgi] > vout_model[mgi - 1]);
        rho_tmodel = (log_rho > -90) ? pow(10., log_rho) : 0.;
      } else if (get_modelgridtype() == GridType::CYLINDRICAL2D) {
        // columns: cell number, r midpoint, z midpoint, density
        float cell_r_in{NAN};
        float cell_z_in{NAN};
        assert_always(parse_next_token(remainder, cellnumberin) && parse_next_token(remainder, cell_r_in) &&
                      parse_next_token(remainder, cell_z_in) && parse_next_token(remainder, rho_tmodel));

        const int n_rcyl = (mgi % ncoord_model[0]);
        const double pos_r_cyl_mid = (n_rcyl + 0.5) * globals::vmax * t_model / ncoord_model[0];
        assert_always(fabs((cell_r_in / pos_r_cyl_mid) - 1) < 1e-3);
        const int n_z = (mgi / ncoord_model[0]);
        if (((2 * n_z) + 1) == ncoord_model[1]) {
          // an odd number of z rows puts the middle row at the origin, where a relative comparison
          // has nothing to divide by. Compare against the cell size instead. The row is identified
          // from the integer index rather than by testing the coordinate against zero, because the
          // expression below is not required to evaluate to exactly zero under -ffast-math.
          assert_always(fabs(cell_z_in) < (1e-3 * globals::vmax * t_model / ncoord_model[1]));
        } else {
          const double pos_z_mid = globals::vmax * t_model * (-1 + (2 * (n_z + 0.5) / ncoord_model[1]));
          assert_always(fabs((cell_z_in / pos_z_mid) - 1) < 1e-3);
        }
      } else {
        // columns: cell number, x midpoint, y midpoint, z midpoint, density
        std::array<float, 3> cellpos_in{};
        float rho_model_in{NAN};
        assert_always(parse_next_token(remainder, cellnumberin) && parse_next_token(remainder, cellpos_in[0]) &&
                      parse_next_token(remainder, cellpos_in[1]) && parse_next_token(remainder, cellpos_in[2]) &&
                      parse_next_token(remainder, rho_model_in));
        rho_tmodel = rho_model_in;

        // cell coordinates in the 3D model.txt file are sometimes reordered by the scaling script
        // however, the cellindex always should increment X first, then Y, then Z
        const double xmax_tmodel = globals::vmax * t_model;
        for (int axis = 0; axis < 3; axis++) {
          const double cellwidth = 2 * xmax_tmodel / ncoordgrid[axis];
          const double cellpos_expected = -xmax_tmodel + (cellwidth * get_cellcoordindex(mgi, axis));
          if (fabs(cellpos_expected - cellpos_in[axis]) > 0.5 * cellwidth) {
            posmatch_xyz = false;
          }
          if (fabs(cellpos_expected - cellpos_in[2 - axis]) > 0.5 * cellwidth) {
            posmatch_zyx = false;
          }
        }

        if (min_den < 0. || min_den > rho_tmodel) {
          min_den = rho_tmodel;
        }
      }

      if (mgi == 0) {
        first_input_cellid = cellnumberin;
        printlnlog("first_input_cellid {}", first_input_cellid);
        assert_always(first_input_cellid == 0 || first_input_cellid == 1);
      }
      assert_always(cellnumberin == mgi + first_input_cellid);

      if (rho_tmodel < 0) {
        printlnlog("[error] model.txt cell {} (inputcellid {}) has negative density {:g} [g/cm3] at t_model. aborting",
                   mgi, cellnumberin, rho_tmodel);
        std::abort();
      }

      const bool keepcell = (rho_tmodel > 0);
      set_rho_tmin(mgi, static_cast<float>(rho_tmodel * pow3(t_model / globals::tmin)));

      read_model_radioabundances(fmodel, remainder, mgi, keepcell, colnames, nucindexlist, one_line_per_cell);

      mgi++;
    }

    if (mgi != get_npts_model()) {
      printlnlog("[error] model.txt: found only {} cells instead of {} expected.", mgi, get_npts_model());
      std::abort();
    }

    if (get_modelgridtype() == GridType::SPHERICAL1D) {
      // for 1D models, vmax is the outer velocity of the last cell
      globals::vmax = vout_model[get_npts_model() - 1];
    } else if (get_modelgridtype() == GridType::CARTESIAN3D) {
      // posmatch_xyz and posmatch_zyx should be mutually exclusive: if both column orders matched, an
      // infinity probably occurred in the calculated positions
      if (posmatch_xyz) {
        printlnlog(
            "Cell positions in model.txt are consistent with calculated values when x-y-z column order is used.");
      } else if (posmatch_zyx) {
        printlnlog(
            "Cell positions in model.txt are consistent with calculated values when z-y-x column order is used.");
      } else {
        printlnlog(
            "[warning] Cell positions in model.txt are not consistent with calculated values in either x-y-z or z-y-x "
            "order.");
      }
      printlnlog("minimum model density {:g} [g/cm3]", min_den);
    }
  }

  // share the rank-local results of the cell parsing with the other ranks
  MPI_Bcast_safe(first_input_cellid, 0, MPI_COMM_WORLD);
  MPI_Bcast_safe(globals::vmax, 0, MPI_COMM_WORLD);
  if (get_modelgridtype() == GridType::SPHERICAL1D) {
    MPI_Bcast_safe(vout_model, 0, MPI_COMM_WORLD);
  }

  // send the node-shared model data filled by rank 0 to the node leaders of the other nodes
  if (globals::rank_in_node == 0) {
    MPI_Bcast_safe(modelgrid_input, 0, globals::mpi_comm_internode);
    MPI_Bcast_safe(initnucmassfrac_allcells, 0, globals::mpi_comm_internode);
  }
  MPI_Barrier_allranks();

  assert_always(get_npts_model() ==
                std::max(1, ncoord_model[0]) * std::max(1, ncoord_model[1]) * std::max(1, ncoord_model[2]));
  printlnlog("npts_model: {}", get_npts_model());
  globals::rmax = globals::vmax * globals::tmin;
  printlnlog("vmax {:g} [cm/s] ({:.2f}c)", globals::vmax, globals::vmax / CLIGHT);
  assert_always(globals::vmax < CLIGHT);
  printlnlog("tmin {:g} [s] = {:.2f} [d]", globals::tmin, globals::tmin / DAY);
  printlnlog("rmax {:g} [cm] (at t=tmin)", globals::rmax);

  calc_modelinit_totmassnuclides();

  printlnlog("Total input model mass: {:9.3e} [Msun]", mtot_input / MSUN);
  printlnlog("Nuclide masses at t=t_model_init [Msun]:");
  printlnlog("  56Ni: {:9.3e}  56Co: {:9.3e}  52Fe: {:9.3e}  48Cr: {:9.3e}", get_totmassnuclide_tmodel(28, 56) / MSUN,
             get_totmassnuclide_tmodel(27, 56) / MSUN, get_totmassnuclide_tmodel(26, 52) / MSUN,
             get_totmassnuclide_tmodel(24, 48) / MSUN);
  printlnlog("  Fe-group: {:9.3e}  57Ni: {:9.3e}  57Co: {:9.3e}", mfegroup / MSUN,
             get_totmassnuclide_tmodel(28, 57) / MSUN, get_totmassnuclide_tmodel(27, 57) / MSUN);

  report_negative_input_massfracs();

  if (!ignored_model_columns.empty()) {
    printlog(
        "[warning] ignored model.txt columns not recognised as a known nuclide, X_Fegroup, Ye, q, or tracercount:");
    for (const auto& colname : ignored_model_columns) {
      printlog(" '{}'", colname);
    }
    printlnlog("");
  }

  read_possible_yefile();
}

void write_grid_restart_data(const int timestep) {
  const auto filename = std::format("gridsave_ts{}.tmp", timestep);

  const auto sys_time_start_write_restart = std::chrono::steady_clock::now();
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
    const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);

    assert_always(globals::dep_estimator_gamma[nonemptymgi] >= 0.);
    fprintf(gridsave_file, "%d %a %a %a %a %d %la %la %la %la %a %a", mgi, TR_allcells[nonemptymgi],
            Te_allcells[nonemptymgi], W_allcells[nonemptymgi], TJ_allcells[nonemptymgi],
            static_cast<int>(thick_allcells[nonemptymgi]), globals::dep_estimator_gamma[nonemptymgi],
            globals::dep_estimator_positron[nonemptymgi], globals::dep_estimator_electron[nonemptymgi],
            globals::dep_estimator_alpha[nonemptymgi], nne_allcells[nonemptymgi], nnetot_allcells[nonemptymgi]);

    if constexpr (USE_LUT_PHOTOION) {
      for (int i = 0; i < globals::nbfcontinua_ground; i++) {
        const ptrdiff_t estimindex = (static_cast<ptrdiff_t>(nonemptymgi) * globals::nbfcontinua_ground) + i;
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
  const auto write_restart_duration =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_write_restart).count();
  printlnlog("done in {:.1f} seconds.", write_restart_duration);
}

// get lowest modelgridindex assigned to this rank (for update_grid and output files)
auto get_nstart(const int rank) -> int {
  if (ranks_ndo.empty()) {
    setup_nstart_ndo();
  }
  return ranks_nstart[rank];
}

// get lowest nonemptymgi assigned to this rank (for update_grid and output files)
auto get_nstart_nonempty(const int rank) -> int {
  if (ranks_ndo.empty()) {
    setup_nstart_ndo();
  }
  return ranks_nstart_nonempty[rank];
}

// get the count of modelgridindices assigned to this rank (for update_grid and output files)
auto get_ndo(const int rank) -> int {
  if (ranks_ndo.empty()) {
    setup_nstart_ndo();
  }
  return ranks_ndo[rank];
}

// get the count of nonemptymgi assigned to this rank (for update_grid and output files)
auto get_ndo_nonempty(const int rank) -> int {
  if (ranks_ndo.empty()) {
    setup_nstart_ndo();
  }
  return ranks_ndo_nonempty[rank];
}

// Initialise the propagation grid cells and associate them with modelgrid cells
void init_grid() {
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
    auto grid_file = fstream_required("grid.out", std::ios::out | std::ios::trunc);
    for (int cellindex = 0; cellindex < ngrid; cellindex++) {
      const int mgi = get_propcell_modelgridindex(cellindex);
      if (mgi >= 0) {
        std::println(grid_file, "{} {}", cellindex, mgi);  // write only non-empty cells to grid file
      }
    }
  }

  allocate_nonemptymodelcells();
  read_elem_abundances();

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
        const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
        totmassnuclide_actual +=
            get_modelinitnucmassfrac(mgi, nucindex) * get_rho_tmin(mgi) * get_modelcell_assocvolume_tmin(mgi);
      }

      if (totmassnuclide_actual > 0.) {
        const double ratio = totmassnuclide[nucindex] / totmassnuclide_actual;
        for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
          const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);
          const double prev_massfrac = get_modelinitnucmassfrac(mgi, nucindex);
          const auto new_massfrac = static_cast<float>(prev_massfrac * ratio);
          set_modelinitnucmassfrac(mgi, nucindex, new_massfrac);
        }
      }
    }
  }

  MPI_Barrier_node();

  double mtot_mapped = 0.;
  for (int mgi = 0; mgi < get_npts_model(); mgi++) {
    mtot_mapped += get_rho_tmin(mgi) * get_modelcell_assocvolume_tmin(mgi);
  }
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

// identify the cell index from an (x,y,z) position and a time.
[[nodiscard]] DEVICE_FUNC auto get_cellindex_from_pos(const Vec3d& pos, const double time) -> int {
  const auto prop_gridtype = get_propgridtype();
  const auto ndim = get_ndim(prop_gridtype);
  const auto posgridcoords = get_gridcoords_from_xyz(pos, prop_gridtype);
  int cellindex = 0;
  for (int d = 0; d < ndim; d++) {
    if (std::abs(posgridcoords[d]) > (globals::vmax * time)) {
      // outside grid
      return -99;
    }
    cellindex += get_coordcellindexstride(d) * get_poscoordpointnum(posgridcoords[d], time, d);
  }
  assert_always(cellindex >= 0);
  assert_always(cellindex < ngrid);

  // do a check that the position is within the cell
  const double trat = time / globals::tmin;
  for (int n = 0; n < ndim; n++) {
    const double cellposmin = get_cellcoordmin(cellindex, n) * trat;
    const double cellposmax = get_cellcoordmax(cellindex, n) * trat;
    assert_always(posgridcoords[n] >= (cellposmin - cellbound_tolerance(cellposmin)));
    assert_always(posgridcoords[n] <= (cellposmax + cellbound_tolerance(cellposmax)));
  }
  return cellindex;
}

// Re-establish the invariant that pos lies inside propagation cell cellindex at the given
// time by clamping each grid coordinate into the cell bounds. Called after a boundary
// crossing, this snaps the position exactly onto the crossed cell face(s), so that
// floating-point rounding errors from the position updates (about one ulp per move, i.e.
// up to several cm at the outer grid) cannot accumulate over many crossings.
DEVICE_FUNC void snap_pos_to_cell(Vec3d& pos, const double time, const int cellindex) {
  if (get_propgridtype() != GridType::CARTESIAN3D) {
    return;
  }
  for (int d = 0; d < 3; d++) {
    const int idx = get_cellcoordindex(cellindex, d);
    const double cellposmin = coord_pos_min_tmin[d][idx] / globals::tmin * time;
    // exactly match the boundary used by boundary_distance(): the upper boundary is the
    // lower edge of the neighbouring cell, except at the grid edge
    const double cellposmax = (idx < (ncoordgrid[d] - 1)) ? coord_pos_min_tmin[d][idx + 1] / globals::tmin * time
                                                          : get_cellcoordmax(cellindex, d) / globals::tmin * time;
    const double newpos_d = std::clamp(pos[d], cellposmin, cellposmax);
    // corrections should only ever be at the floating-point rounding error level
    assert_testmodeonly(std::abs(newpos_d - pos[d]) <= cellbound_tolerance(pos[d]));
    pos[d] = newpos_d;
  }
}

// compute distance to a cell boundary.
[[nodiscard]] DEVICE_FUNC auto boundary_distance(const Vec3d& dir, const Vec3d& pos, const double tstart,
                                                 const int cellindex) -> std::tuple<double, int> {
  const auto prop_gridtype = get_propgridtype();
  if constexpr (FORCE_SPHERICAL_ESCAPE_SURFACE) {
    if (get_cell_r_inner(cellindex, prop_gridtype) > globals::rmax) {
      return {0., -99};
    }
  }

  // d is used to loop over the coordinate indices 0,1,2 for x,y,z

  // the following vector are in grid coordinates, so either x,y,z (3D) or r (1D), or r_xy, z (2D)

  const auto pktposgridcoord = get_gridcoords_from_xyz(pos, prop_gridtype);

  // dir * CLIGHT_PROP converted from xyz to grid coordinates
  const auto pktvelgridcoord = get_gridcoords_vel_from_xyz_pos_dir(pos, dir, pktposgridcoord, prop_gridtype);

  const auto [cellcoordidx, cellcoordmin, cellcoordmax] = [cellindex, prop_gridtype] {
    auto _cellcoordidx = std::array<int, 3>{};
    auto _cellcoordmin = std::array<double, 3>{};
    auto _cellcoordmax = std::array<double, 3>{};
    for (int d = 0; d < get_ndim(prop_gridtype); d++) {
      _cellcoordidx[d] = get_cellcoordindex(cellindex, d);

      _cellcoordmin[d] = get_cellcoordmin(cellindex, d);
      _cellcoordmax[d] = get_cellcoordmax(cellindex, d);
    }
    return std::make_tuple(_cellcoordidx, _cellcoordmin, _cellcoordmax);
  }();

  {
    for (int d = 0; d < get_ndim(prop_gridtype); d++) {
      // pos_component_vel_relative_to_flow is constant along a ray with a given direction in Cartesian coordinates, but
      // for non-Cartesian coordinates, we still need to check at the current position whether the packet is
      // moving in the positive or negative direction in each grid coordinate direction relative to the homologous grid
      // flow, otherwise we might never enter the cell that we're supposed to be in
      const bool pos_component_vel_relative_to_flow = (pktvelgridcoord[d] * tstart) > pktposgridcoord[d];

      if constexpr (TESTMODE) {
        bool isoutside_error = false;
        [[maybe_unused]] double delta = 0.;
        if (pos_component_vel_relative_to_flow) {
          // check if packet pos is above cell max while moving in the positive direction relative to the grid flow
          const double boundaryposmax = cellcoordmax[d] / globals::tmin * tstart;
          delta = pktposgridcoord[d] - boundaryposmax;
          isoutside_error = pktposgridcoord[d] > (boundaryposmax + cellbound_tolerance(boundaryposmax));
        } else {
          // check if packet pos is below cell min while moving in the negative direction relative to the grid flow
          const double boundaryposmin = cellcoordmin[d] / globals::tmin * tstart;
          delta = pktposgridcoord[d] - boundaryposmin;
          isoutside_error = pktposgridcoord[d] < (boundaryposmin - cellbound_tolerance(boundaryposmin));
        }

        if (isoutside_error) {
#ifndef GPU_ON
          printlnlog(
              "[error] timestep {}: packet outside coord {} {}{} boundary of cell {} by delta {:g}. vel {:g} initpos "
              "{:g} cellcoordmin {:g} cellcoordmax {:g} dir [{:g}, {:g}, {:g}] tmin {:g} s tstart {:g} s",
              globals::timestep, d, pos_component_vel_relative_to_flow ? '+' : '-', get_coordlabel(prop_gridtype, d),
              cellindex, delta, pktvelgridcoord[d], pktposgridcoord[d], cellcoordmin[d] / globals::tmin * tstart,
              cellcoordmax[d] / globals::tmin * tstart, dir[0], dir[1], dir[2], globals::tmin, tstart);
#endif

          // this should not happen! Leave the check until late 2026 and if it never triggers on any runs, we can remove
          // the check and correction code
          assert_always(!isoutside_error);

          const auto next_cellindex = get_cellindex_from_pos(pos, tstart);
          if ((cellcoordidx[d] == (ncoordgrid[d] - 1) && pos_component_vel_relative_to_flow) ||
              (cellcoordidx[d] == 0 && !pos_component_vel_relative_to_flow) || (next_cellindex < 0)) {
#ifndef GPU_ON
            printlnlog("[warning] treating out-of-boundary packet in cell {} as escaping the grid", cellindex);
#endif
            return {0., -99};
          }
#ifndef GPU_ON
          printlnlog(
              "[warning] swapping packet cellindex from {} to {}, which has cellcoordmin {:g}, cellcoordmax {:g}",
              cellindex, next_cellindex, get_cellcoordmin(next_cellindex, d) / globals::tmin * tstart,
              get_cellcoordmax(next_cellindex, d) / globals::tmin * tstart);
#endif
          return {0., next_cellindex};
        }
      }
    }
  }

  double distance = std::numeric_limits<double>::max();
  int next_cellindex{-1};

  if (prop_gridtype == GridType::SPHERICAL1D) {
    // the only coordinate is the radius from the origin

    const double speed = vec_len(dir) * CLIGHT_PROP;  // just in case dir is not normalised

    const double r_outer = cellcoordmax[0] * tstart / globals::tmin;
    const double d_coordmaxboundary =
        is_boundary_overshoot_within_tolerance<BoundaryType::UPPER>(pktposgridcoord[0], pktvelgridcoord[0],
                                                                    cellcoordmax[0], tstart)
            ? 0.
            : expanding_shell_intersection<BoundaryType::UPPER>(pos, dir, speed, r_outer, tstart);

    // upper d coordinate of the current cell
    if ((d_coordmaxboundary >= 0.) && (d_coordmaxboundary < distance)) {
      distance = d_coordmaxboundary;
      next_cellindex = (cellcoordidx[0] == (ncoordgrid[0] - 1)) ? -99 : cellindex + get_coordcellindexstride(0);
    }

    const double r_inner = cellcoordmin[0] * tstart / globals::tmin;
    if (r_inner > 0.) {
      const double d_coordminboundary =
          is_boundary_overshoot_within_tolerance<BoundaryType::LOWER>(pktposgridcoord[0], pktvelgridcoord[0],
                                                                      cellcoordmin[0], tstart)
              ? 0.
              : expanding_shell_intersection<BoundaryType::LOWER>(pos, dir, speed, r_inner, tstart);
      // lower d coordinate of the current cell
      if ((d_coordminboundary >= 0.) && (d_coordminboundary < distance)) {
        distance = d_coordminboundary;
        next_cellindex = (cellcoordidx[0] == 0) ? -99 : cellindex - get_coordcellindexstride(0);
      }
    }
  } else if (prop_gridtype == GridType::CYLINDRICAL2D) {
    // coordinate 0 is cylindrical radius (distance from z=0 in x-y plane), coord 1 is z

    const std::array<double, 2> posnoz = {pos[0], pos[1]};

    // r_cyl component of direction vector
    const double dirxylen = std::sqrt(pow2(dir[0]) + pow2(dir[1]));
    // r_cyl component of velocity
    const double xyspeed = dirxylen * CLIGHT_PROP;

    if (dirxylen > 0.) {
      // make a normalised 2D direction vector in the xy plane
      const std::array<double, 2> dirnoz = {dir[0] / dirxylen, dir[1] / dirxylen};

      const double r_outer = cellcoordmax[0] * tstart / globals::tmin;
      const double d_rcyl_coordmaxboundary =
          is_boundary_overshoot_within_tolerance<BoundaryType::UPPER>(pktposgridcoord[0], pktvelgridcoord[0],
                                                                      cellcoordmax[0], tstart)
              ? 0.
              : expanding_shell_intersection<BoundaryType::UPPER>(posnoz, dirnoz, xyspeed, r_outer, tstart);
      if (d_rcyl_coordmaxboundary >= 0.) {
        // how far did the packet travel in the z direction during this time?
        const double d_z_coordmaxboundary = d_rcyl_coordmaxboundary / xyspeed * dir[2] * CLIGHT_PROP;
        // distance from two perpendicular components to the r_cyl upper boundary
        const double d_coordmaxboundary_rcyl = std::sqrt(pow2(d_rcyl_coordmaxboundary) + pow2(d_z_coordmaxboundary));
        if ((d_coordmaxboundary_rcyl >= 0.) && (d_coordmaxboundary_rcyl < distance)) {
          distance = d_coordmaxboundary_rcyl;
          next_cellindex = (cellcoordidx[0] == (ncoordgrid[0] - 1)) ? -99 : cellindex + get_coordcellindexstride(0);
        }
      }

      const double r_inner = cellcoordmin[0] * tstart / globals::tmin;
      // don't try to calculate the intersection if the inner radius is zero
      if (r_inner > 0) {
        // calculate the distance in the xy plane to the inner boundary
        const double d_rcyl_coordminboundary =
            is_boundary_overshoot_within_tolerance<BoundaryType::LOWER>(pktposgridcoord[0], pktvelgridcoord[0],
                                                                        cellcoordmin[0], tstart)
                ? 0.
                : expanding_shell_intersection<BoundaryType::LOWER>(posnoz, dirnoz, xyspeed, r_inner, tstart);
        if (d_rcyl_coordminboundary >= 0.) {
          const double d_z_coordminboundary = d_rcyl_coordminboundary / xyspeed * dir[2] * CLIGHT_PROP;
          // distance from two perpendicular components to the r_cyl lower boundary
          const double d_coordminboundary_rcyl = std::sqrt(pow2(d_rcyl_coordminboundary) + pow2(d_z_coordminboundary));
          if ((d_coordminboundary_rcyl >= 0.) && (d_coordminboundary_rcyl < distance)) {
            distance = d_coordminboundary_rcyl;
            next_cellindex = (cellcoordidx[0] == 0) ? -99 : cellindex - get_coordcellindexstride(0);
          }
        }
      }
    } else {
      // The packet is moving exactly along the z axis (reachable when the comoving emission direction
      // is sampled at exactly costheta = +/-1 and the local velocity is parallel to z), so the
      // general-direction code above would divide by zero. The packet's cylindrical radius stays
      // constant while the grid expands: the receding outer r_cyl boundary can never be crossed, and
      // the expanding inner boundary catches up to the packet at t_cross = rcyl / (rcyl_inner_tmin / tmin)
      const double rcyl_inner_tmin = cellcoordmin[0];
      if (rcyl_inner_tmin > 0.) {
        const double d_coordminboundary_rcyl =
            is_boundary_overshoot_within_tolerance<BoundaryType::LOWER>(pktposgridcoord[0], pktvelgridcoord[0],
                                                                        cellcoordmin[0], tstart)
                ? 0.
                : ((pktposgridcoord[0] * globals::tmin / rcyl_inner_tmin) - tstart) * CLIGHT_PROP;
        if ((d_coordminboundary_rcyl >= 0.) && (d_coordminboundary_rcyl < distance)) {
          distance = d_coordminboundary_rcyl;
          next_cellindex = (cellcoordidx[0] == 0) ? -99 : cellindex - get_coordcellindexstride(0);
        }
      }
    }

    // Z boundaries are Cartesian
    constexpr int d = 1;
    if (pktvelgridcoord[d] > (cellcoordmax[d] / globals::tmin)) {
      const double d_coordmaxboundary_z =
          is_boundary_overshoot_within_tolerance<BoundaryType::UPPER>(pktposgridcoord[d], pktvelgridcoord[d],
                                                                      cellcoordmax[d], tstart)
              ? 0.
              : distance_cartesian_boundary(pktposgridcoord[1], pktvelgridcoord[d], cellcoordmax[d], tstart);

      if ((d_coordmaxboundary_z >= 0.) && (d_coordmaxboundary_z < distance)) {
        distance = d_coordmaxboundary_z;
        next_cellindex = (cellcoordidx[d] == (ncoordgrid[d] - 1)) ? -99 : cellindex + get_coordcellindexstride(d);
      }
    } else if (pktvelgridcoord[d] < (cellcoordmin[d] / globals::tmin)) {
      const double d_coordminboundary_z =
          is_boundary_overshoot_within_tolerance<BoundaryType::LOWER>(pktposgridcoord[d], pktvelgridcoord[d],
                                                                      cellcoordmin[d], tstart)
              ? 0.
              : distance_cartesian_boundary(pktposgridcoord[d], pktvelgridcoord[d], cellcoordmin[d], tstart);

      if ((d_coordminboundary_z >= 0.) && (d_coordminboundary_z < distance)) {
        distance = d_coordminboundary_z;
        next_cellindex = (cellcoordidx[d] == 0) ? -99 : cellindex - get_coordcellindexstride(d);
      }
    }

  } else if (prop_gridtype == GridType::CARTESIAN3D) {
    // There are six possible boundary crossings. Each of the three
    // cartesian coordinates may be taken in turn. For x, the packet trajectory is
    // x = x0 + (dir.x) * c * (t - tstart)
    // the boundaries follow
    // x+/- = x+/-(tmin) * (t/tmin)
    // so the crossing occurs when
    // t - tstart = (x0 - x+/-(tmin)/tmin * tstart) / (x+/-(tmin)/tmin - (dir.x)*c)
    // distance = c * (t - tstart)

    // Modified so that it also returns the distance to the closest cell boundary, regardless of direction.

    for (int d = 0; d < 3; d++) {
      if (pktvelgridcoord[d] > (cellcoordmax[d] / globals::tmin)) {
        const double d_coordmaxboundary =
            is_boundary_overshoot_within_tolerance<BoundaryType::UPPER>(pktposgridcoord[d], pktvelgridcoord[d],
                                                                        cellcoordmax[d], tstart)
                ? 0.
                : distance_cartesian_boundary(pktposgridcoord[d], pktvelgridcoord[d], cellcoordmax[d], tstart);

        if ((d_coordmaxboundary >= 0.) && (d_coordmaxboundary < distance)) {
          distance = d_coordmaxboundary;
          next_cellindex = (cellcoordidx[d] == (ncoordgrid[d] - 1)) ? -99 : cellindex + get_coordcellindexstride(d);
        }
      } else if (pktvelgridcoord[d] < (cellcoordmin[d] / globals::tmin)) {
        const double d_coordminboundary =
            is_boundary_overshoot_within_tolerance<BoundaryType::LOWER>(pktposgridcoord[d], pktvelgridcoord[d],
                                                                        cellcoordmin[d], tstart)
                ? 0.
                : distance_cartesian_boundary(pktposgridcoord[d], pktvelgridcoord[d], cellcoordmin[d], tstart);

        // lower d coordinate of the current cell
        if ((d_coordminboundary >= 0.) && (d_coordminboundary < distance)) {
          distance = d_coordminboundary;
          next_cellindex = (cellcoordidx[d] == 0) ? -99 : cellindex - get_coordcellindexstride(d);
        }
      }
    }
  } else {
    assert_always(false);
  }

  assert_always((next_cellindex == -99) || ((next_cellindex >= 0) && (next_cellindex < ngrid)));

  const double maxboundarydist =
      (prop_gridtype == GridType::CARTESIAN3D)
          ? std::numbers::sqrt3 * globals::rmax * (tstart + (distance / CLIGHT_PROP)) / globals::tmin
          : 2 * globals::rmax * (tstart + (distance / CLIGHT_PROP)) / globals::tmin;

  assert_always(distance >= 0.);
  assert_always(distance <= maxboundarydist);

  if (distance > globals::max_path_step) {
    return {globals::max_path_step, cellindex};
  }

  return {distance, next_cellindex};
}

}  // namespace grid
