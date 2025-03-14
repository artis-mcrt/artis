#include "grid.h"

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "input.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "radfield.h"
#include "rpkt.h"
#include "sn3d.h"
#include "vectors.h"

namespace grid {

namespace {

struct ModelGridCellInput {
  float rhoinit = -1.;
  float ffegrp = 0.;
  float initial_radial_pos_sum = 0.;
  float initelectronfrac = -1;  // Ye: electrons (or protons) per nucleon
  float initenergyq = 0.;       // q: energy in the model at tmin to use with USE_MODEL_INITIAL_ENERGY [erg/g]
};

std::array<char, 3> coordlabel{'?', '?', '?'};

std::array<int, 3> ncoordgrid{};  // propagation grid dimensions

GridType model_type = GridType::CARTESIAN3D;
ptrdiff_t npts_model = 0;           // number of model grid cells
ptrdiff_t nonempty_npts_model = 0;  // number of allocated non-empty model grid cells

double t_model = -1.;  // time at which densities in input model are correct.
std::vector<double> vout_model{};
std::array<int, 3> ncoord_model{};  // the model.txt input grid dimensions

double min_den;  // minimum model density

double mfegroup = 0.;  // Total mass of Fe group elements in ejecta

int first_cellindex = -1;  // auto-determine first cell index in model.txt (usually 1 or 0)

// Initial co-ordinates of inner most corner of cell.
std::vector<std::array<double, 3>> propcell_pos_min{};

// associate each propagation cell with a model grid cell, or not, if the cell is empty (or doesn't get mapped to
// anything such as 1D/2D to 3D)
std::vector<int> propcell_mgi;
std::vector<int> propcell_nonemptymgi;

std::vector<int> modelgrid_numpropcells;
std::vector<int> nonemptymgi_of_mgi;
std::vector<int> mgi_of_nonemptymgi;

std::span<double> totmassradionuclide{};  // total mass of each radionuclide in the ejecta

MPI_Win win_nltepops_allcells = MPI_WIN_NULL;
MPI_Win win_initnucmassfrac_allcells = MPI_WIN_NULL;

std::span<float> initnucmassfrac_allcells{};
std::span<float> initmassfracuntrackedstable_allcells{};
std::span<int> elements_uppermost_ion_allcells{};  // Highest ion index that has a significant population

std::vector<int> ranks_nstart;
std::vector<int> ranks_nstart_nonempty;
std::vector<int> ranks_ndo;
std::vector<int> ranks_ndo_nonempty;
inline std::span<ModelGridCellInput> modelgrid_input{};

// Get number of dimensions
consteval auto get_ndim(const GridType gridtype) -> int {
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

void set_rho_tmin(const int modelgridindex, const float x) { modelgrid_input[modelgridindex].rhoinit = x; }

void set_initelectronfrac(const int modelgridindex, const double electronfrac) {
  modelgrid_input[modelgridindex].initelectronfrac = electronfrac;
}

void read_possible_yefile() {
  if (!std::filesystem::exists("Ye.txt")) {
    printout("Ye.txt not found\n");
    return;
  }

  FILE *filein = fopen_required("Ye.txt", "r");
  int nlines_in = 0;
  assert_always(fscanf(filein, "%d", &nlines_in) == 1);

  for (int n = 0; n < nlines_in; n++) {
    int mgiplusone = -1;
    double initelecfrac = 0.;
    assert_always(fscanf(filein, "%d %lg", &mgiplusone, &initelecfrac) == 2);
    const int mgi = mgiplusone - 1;
    if (mgi >= 0 && mgi < get_npts_model()) {
      set_initelectronfrac(mgi, initelecfrac);
    }
  }
  fclose(filein);
}

void allocate_initradiobund() {
  assert_always(npts_model > 0);

  const ptrdiff_t num_nuclides = decay::get_num_nuclides();

  const size_t totalradioabundcount = (npts_model + 1) * num_nuclides;
  std::tie(initnucmassfrac_allcells, win_initnucmassfrac_allcells) =
      MPI_shared_malloc_keepwin_span<float>(totalradioabundcount);
  printout(
      "[info] mem_usage: radioabundance data for %td nuclides for %td cells occupies %.3f MB (node shared memory)\n",
      num_nuclides, npts_model, static_cast<double>(totalradioabundcount * sizeof(float)) / 1024. / 1024.);

  MPI_Barrier(globals::mpi_comm_node);

  assert_always(initnucmassfrac_allcells.data() != nullptr);

  for (ptrdiff_t mgi = 0; mgi < (npts_model + 1); mgi++) {
    if (mgi % static_cast<ptrdiff_t>(globals::node_nprocs) == globals::rank_in_node) {
      std::fill_n(&initnucmassfrac_allcells[mgi * num_nuclides], num_nuclides, 0.);
    }
  }
  MPI_Barrier(globals::mpi_comm_node);
}

auto get_cell_r_inner(const int cellindex) -> double {
  if constexpr (GRID_TYPE == GridType::SPHERICAL1D) {
    return get_cellcoordmin(cellindex, 0);
  }

  if constexpr (GRID_TYPE == GridType::CYLINDRICAL2D) {
    const auto rcyl_inner = get_cellcoordmin(cellindex, 0);
    const auto z_inner = std::min(std::abs(get_cellcoordmin(cellindex, 1)), std::abs(get_cellcoordmax(cellindex, 1)));
    return std::sqrt(std::pow(rcyl_inner, 2) + std::pow(z_inner, 2));
  }

  if constexpr (GRID_TYPE == GridType::CARTESIAN3D) {
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
    printout("WARNING: Fe-group mass fraction %g is negative in cell %d\n", x, modelgridindex);
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
  assert_testmodeonly(new_modelgridindex >= 0);
  assert_testmodeonly(new_modelgridindex <= get_npts_model());
  propcell_mgi[cellindex] = new_modelgridindex;
}

void set_modelinitnucmassfrac(const int modelgridindex, const int nucindex, float abund) {
  // set the mass fraction of a nuclide in a model grid cell at t=t_model by nuclide index
  // initnucmassfrac array is in node shared memory
  assert_always(nucindex >= 0);
  if (!(abund >= 0.)) {
    printout("WARNING: nuclear mass fraction for nucindex %d = %g is negative in cell %d\n", nucindex, abund,
             modelgridindex);
    assert_always(abund > -1e-6);
    abund = 0.;
  }

  assert_always(abund >= 0.);
  assert_always(abund <= 1.);
  const ptrdiff_t num_nuclides = decay::get_num_nuclides();

  initnucmassfrac_allcells[(modelgridindex * num_nuclides) + nucindex] = abund;
}

void set_initenergyq(const int modelgridindex, const double initenergyq) {
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
      printout("WARNING: cell %d Z=%d element abundance is less than the sum of its radioisotope abundances\n", mgi,
               atomic_number);
      printout("  massfrac(Z) %g massfrac_radioisotopes(Z) %g\n", elemabundance, isofracsum);
      printout("  increasing elemental abundance to %g and setting stable isotopic abundance to zero\n", isofracsum);
    }
    // result is allowed to be slightly negative due to roundoff error
    assert_always(massfrac_untrackedstable >= -1e-2);
    massfrac_untrackedstable = 0.;  // bring up to zero if negative
  }

  // if (globals::rank_in_node == 0)
  {
    initmassfracuntrackedstable_allcells[(nonemptymgi * get_nelements()) + element] = massfrac_untrackedstable;
  }

  // (isofracsum + massfracstable) might not exactly match elemabundance if we had to boost it to reach isofracsum
  set_elem_abundance(nonemptymgi, element, isofracsum + massfrac_untrackedstable);
}

void allocate_nonemptycells_composition_cooling()
// Initialise composition dependent cell data for the given cell
{
  const ptrdiff_t nonempty_npts_model_ptrdifft = get_nonempty_npts_model();
  const auto nelements = get_nelements();

  initmassfracuntrackedstable_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model_ptrdifft * nelements);
  elem_meanweight_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model_ptrdifft * nelements);
  elements_uppermost_ion_allcells = MPI_shared_malloc_span<int>(nonempty_npts_model_ptrdifft * nelements);
  elem_massfracs_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model_ptrdifft * nelements);
  ion_groundlevelpops_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model_ptrdifft * get_includedions());
  ion_partfuncts_allcells = MPI_shared_malloc_span<float>(nonempty_npts_model_ptrdifft * get_includedions());
  ion_cooling_contribs_allcells = MPI_shared_malloc_span<double>(nonempty_npts_model_ptrdifft * get_includedions());

  if (globals::total_nlte_levels > 0) {
    std::tie(nltepops_allcells, win_nltepops_allcells) =
        MPI_shared_malloc_keepwin_span<double>(nonempty_npts_model_ptrdifft * globals::total_nlte_levels);

  } else {
    nltepops_allcells = {};
  }

  if (globals::rank_in_node == 0) {
    std::ranges::fill(initmassfracuntrackedstable_allcells, 0.);
    std::ranges::fill(elem_meanweight_allcells, 0.);
    std::ranges::fill(elements_uppermost_ion_allcells, -1);
    std::ranges::fill(elem_massfracs_allcells, -0.);
    std::ranges::fill(ion_groundlevelpops_allcells, 0.);
    std::ranges::fill(ion_partfuncts_allcells, 0.);
    std::ranges::fill(ion_cooling_contribs_allcells, 0.);
    // -1 indicates that there is currently no information on the nlte populations
    std::ranges::fill(grid::nltepops_allcells, -1.);
  }
  MPI_Barrier(globals::mpi_comm_node);
}

void allocate_nonemptymodelcells() {
  // Determine the number of simulation cells associated with the model cells
  std::ranges::fill(modelgrid_numpropcells, 0);
  if (globals::rank_in_node == 0) {
    for (int mgi = 0; mgi < (get_npts_model() + 1); mgi++) {
      modelgrid_input[mgi].initial_radial_pos_sum = 0.;
    }
  }
  MPI_Barrier(globals::mpi_comm_node);

  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const auto radial_pos_mid = get_cellradialposmid(cellindex);

    if (FORCE_SPHERICAL_ESCAPE_SURFACE && radial_pos_mid > globals::vmax * globals::tmin) {
      // for 1D models, the final shell outer v should already be at vmax
      assert_always(model_type != GridType::SPHERICAL1D || propcell_mgi[cellindex] == get_npts_model());
      set_propcell_modelgridindex(cellindex, get_npts_model());
    }

    const int mgi = get_propcell_modelgridindex(cellindex);
    assert_always(!(get_model_type() == GridType::CARTESIAN3D) || (get_rho_tmin(mgi) > 0) || (mgi == get_npts_model()));

    modelgrid_numpropcells[mgi] += 1;
    if (globals::rank_in_node == 0) {
      modelgrid_input[mgi].initial_radial_pos_sum += radial_pos_mid;
    }

    assert_always(!(get_model_type() == GridType::CARTESIAN3D) || (modelgrid_numpropcells[mgi] == 1) ||
                  (mgi == get_npts_model()));
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
    if (mgi >= get_npts_model()) {
      propcell_nonemptymgi[cellindex] = -1;
    } else {
      propcell_nonemptymgi[cellindex] = get_nonemptymgi_of_mgi(mgi);
    }
  }

  assert_always(modelgrid.data() == nullptr);
  modelgrid = MPI_shared_malloc_span<ModelGridCell>(nonempty_npts_model);
  if (globals::rank_in_node == 0) {
    std::ranges::fill(modelgrid, ModelGridCell{});
  }
  MPI_Barrier(globals::mpi_comm_node);
  allocate_nonemptycells_composition_cooling();

  if constexpr (EXPANSIONOPACITIES_ON || RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY > 0.) {
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
        MPI_shared_malloc_keepwin_span<double>(ionestimcount);

    if (globals::rank_in_node == 0) {
      std::ranges::fill(globals::corrphotoionrenorm, 1.);
    }
    MPI_Barrier(globals::mpi_comm_node);

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

  if (USE_LUT_BFHEATING && ionestimsize > 0) {
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

  printout("[info] mem_usage: the modelgrid array occupies %.3f MB\n",
           (get_npts_model() + 1) * sizeof(modelgrid[0]) / 1024. / 1024.);

  printout("There are %td modelgrid cells with associated propagation cells (nonempty_npts_model)\n",
           nonempty_npts_model);

  printout(
      "[info] mem_usage: NLTE populations for all allocated cells occupy a total of %.3f MB (node shared memory)\n",
      get_nonempty_npts_model() * globals::total_nlte_levels * sizeof(double) / 1024. / 1024.);
}

void map_1dmodelto3dgrid()
// Map 1D spherical model grid onto propagation grid
{
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
      set_propcell_modelgridindex(cellindex, get_npts_model());
    }
  }
}

void map_2dmodelto3dgrid()
// Map 2D cylindrical model onto propagation grid
{
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    // map to 3D Cartesian grid
    const auto pos_mid = std::array<double, 3>{get_cellcoordmin(cellindex, 0) + (0.5 * wid_init(cellindex, 0)),
                                               get_cellcoordmin(cellindex, 1) + (0.5 * wid_init(cellindex, 1)),
                                               get_cellcoordmin(cellindex, 2) + (0.5 * wid_init(cellindex, 2))};

    const double rcylindrical = std::sqrt(std::pow(pos_mid[0], 2) + std::pow(pos_mid[1], 2));

    // 2D grid is uniform so rcyl and z indices can be calculated with no lookup
    const int n_rcyl = static_cast<int>(rcylindrical / globals::tmin / globals::vmax * ncoord_model[0]);
    const int n_z =
        static_cast<int>((pos_mid[2] / globals::tmin + globals::vmax) / (2 * globals::vmax) * ncoord_model[1]);

    if (n_rcyl >= 0 && n_rcyl < ncoord_model[0] && n_z >= 0 && n_z < ncoord_model[1]) {
      const int mgi = (n_z * ncoord_model[0]) + n_rcyl;

      if (modelgrid_input[mgi].rhoinit > 0) {
        set_propcell_modelgridindex(cellindex, mgi);
      } else {
        set_propcell_modelgridindex(cellindex, get_npts_model());
      }
    } else {
      set_propcell_modelgridindex(cellindex, get_npts_model());
    }
  }
}

// mgi and cellindex are interchangeable in this mode (except for empty cells that associated with mgi ==
// get_npts_model())
void map_modeltogrid_direct() {
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const int mgi = (modelgrid_input[cellindex].rhoinit > 0) ? cellindex : get_npts_model();
    set_propcell_modelgridindex(cellindex, mgi);
  }
}

void abundances_read() {
  // barrier to make sure node master has set values in node shared memory
  MPI_Barrier(MPI_COMM_WORLD);
  printout("reading abundances.txt...");
  const bool threedimensional = (get_model_type() == GridType::CARTESIAN3D);

  // Open the abundances file
  auto abundance_file = fstream_required("abundances.txt", std::ios::in);

  // and process through the grid to read in the abundances per cell
  // The abundance file should only contain information for non-empty
  // cells. Its format must be cellnumber (integer), abundance for
  // element Z=1 (float) up to abundance for element Z=30 (float)
  // i.e. in total one integer and 30 floats.

  // loop over propagation cells for 3D models, or modelgrid cells
  for (int mgi = 0; mgi < get_npts_model(); mgi++) {
    std::string line;
    assert_always(get_noncommentline(abundance_file, line));
    std::istringstream ssline(line);

    int cellnumberinput = -1;
    assert_always(ssline >> cellnumberinput);
    assert_always(cellnumberinput == mgi + first_cellindex);

    // the abundances.txt file specifies the elemental mass fractions for each model cell
    // (or proportial to mass frac, e.g. element densities because they will be normalised anyway)
    // The abundances begin with hydrogen, helium, etc, going as far up the atomic numbers as required
    double normfactor = 0.;
    float abundances_in[150] = {0.};
    double abund_in = 0.;
    for (int anumber = 1; anumber <= 150; anumber++) {
      abundances_in[anumber - 1] = 0.;
      if (!(ssline >> abund_in)) {
        // at least one element (hydrogen) should have been specified for nonempty cells
        assert_always(anumber > 1 || get_numpropcells(mgi) == 0);
        break;
      }

      if (abund_in < 0. || abund_in < std::numeric_limits<float>::min()) {
        assert_always(abund_in > -1e-6);
        abund_in = 0.;
      }
      abundances_in[anumber - 1] = static_cast<float>(abund_in);
      normfactor += abundances_in[anumber - 1];
    }

    if (get_numpropcells(mgi) > 0) {
      if (threedimensional || normfactor <= 0.) {
        normfactor = 1.;
      }
      const int nonemptymgi = get_nonemptymgi_of_mgi(mgi);

      for (int element = 0; element < get_nelements(); element++) {
        // now set the abundances (by mass) of included elements, i.e.
        // read out the abundances specified in the atomic data file
        const int anumber = get_atomicnumber(element);
        const float elemabundance = abundances_in[anumber - 1] / normfactor;
        assert_always(elemabundance >= 0.);

        // radioactive nuclide abundances should have already been set by read_??_model
        set_elem_untrackedstable_abund_from_total(nonemptymgi, element, elemabundance);
      }
    }
  }

  // barrier to make sure node master has set values in node shared memory
  MPI_Barrier(MPI_COMM_WORLD);
  printout("done.\n");
}

void parse_model_headerline(const std::string &line, std::vector<int> &zlist, std::vector<int> &alist,
                            std::vector<std::string> &colnames) {
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
      assert_always(get_model_type() == GridType::SPHERICAL1D);
    } else if (token == "rho") {
      // 2D and 3D models have rho [g/cm3]
      assert_always(get_model_type() != GridType::SPHERICAL1D);
      assert_always((columnindex == 4 && get_model_type() == GridType::CARTESIAN3D) ||
                    (columnindex == 3 && get_model_type() == GridType::CYLINDRICAL2D));
      continue;
    } else if (token.starts_with("X_") && token != "X_Fegroup") {
      colnames.push_back(token);
      const int z = decay::get_nucstring_z(token.substr(2));  // + 2 skips the 'X_'
      const int a = decay::get_nucstring_a(token.substr(2));
      assert_always(z >= 0);
      assert_always(a >= 0);
      //   printout("Custom column: '%s' Z %d A %d\n", token.c_str(), z, a);
      zlist.push_back(z);
      alist.push_back(a);
    } else {
      //   printout("Custom column: '%s' Z %d A %d\n", token.c_str(), -1, -1);
      colnames.push_back(token);
      zlist.push_back(-1);
      alist.push_back(-1);
    }
  }
}

auto get_token_count(std::string &line) -> int {
  std::string token;
  int abundcolcount = 0;
  auto ssline = std::istringstream(line);
  while (std::getline(ssline, token, ' ')) {
    if (!std::ranges::all_of(token, isspace)) {  // skip whitespace tokens
      abundcolcount++;
    }
  }
  return abundcolcount;
}

void read_model_radioabundances(std::fstream &fmodel, std::istringstream &ssline_in, const int mgi, const bool keepcell,
                                const std::vector<std::string> &colnames, const std::vector<int> &nucindexlist,
                                const bool one_line_per_cell) {
  std::string line;
  if (!one_line_per_cell) {
    assert_always(std::getline(fmodel, line));
  }

  auto ssline = one_line_per_cell ? std::move(ssline_in) : std::istringstream(line);

  if (!keepcell) {
    return;
  }

  for (ptrdiff_t i = 0; i < std::ssize(colnames); i++) {
    double valuein = 0.;
    assert_always(ssline >> valuein);  // usually a mass fraction, but now can be anything

    if (nucindexlist[i] >= 0) {
      assert_testmodeonly(valuein <= 1.);
      set_modelinitnucmassfrac(mgi, nucindexlist[i], valuein);
    } else if (colnames[i] == "X_Fegroup") {
      set_ffegrp(mgi, valuein);
    } else if (colnames[i] == "cellYe") {
      set_initelectronfrac(mgi, valuein);
    } else if (colnames[i] == "q") {
      // use value for t_model and adjust to tmin with expansion factor
      set_initenergyq(mgi, valuein * t_model / globals::tmin);
    } else if (colnames[i] == "tracercount") {
      ;
    } else {
      if (mgi == 0) {
        printout("WARNING: ignoring column '%s' nucindex %d valuein[mgi=0] %lg\n", colnames[i].c_str(), nucindexlist[i],
                 valuein);
      }
    }
  }
  double valuein = 0.;
  assert_always(!(ssline >> valuein));  // should be no tokens left!
}

auto read_model_columns(std::fstream &fmodel) -> std::tuple<std::vector<std::string>, std::vector<int>, bool> {
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
    switch (model_type) {
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

  printout("model.txt has %s line per cell format\n", one_line_per_cell ? "one" : "two");

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
    printout("model.txt has header line: %s\n", headerline.c_str());
  } else {
    printout("model.txt has no header line. Using default: %s\n", headerline.c_str());
  }

  parse_model_headerline(headerline, zlist, alist, colnames);

  decay::init_nuclides(zlist, alist);

  std::vector<int> nucindexlist(zlist.size());
  for (ptrdiff_t i = 0; i < std::ssize(zlist); i++) {
    nucindexlist[i] = (zlist[i] > 0) ? decay::get_nucindex(zlist[i], alist[i]) : -1;
  }

  allocate_initradiobund();

  return {colnames, nucindexlist, one_line_per_cell};
}

auto get_inputcellvolume(const int mgi) -> double {
  if (get_model_type() == GridType::SPHERICAL1D) {
    const double v_inner = (mgi == 0) ? 0. : vout_model[mgi - 1];
    // mass_in_shell = rho_model[mgi] * (pow(vout_model[mgi], 3) - pow(v_inner, 3)) * 4 * PI * pow(t_model, 3) / 3.;
    return (pow(vout_model[mgi], 3) - pow(v_inner, 3)) * 4 * PI * pow(globals::tmin, 3) / 3.;
  }
  if (get_model_type() == GridType::CYLINDRICAL2D) {
    const int n_r = mgi % ncoord_model[0];
    const double dcoord_rcyl = globals::vmax * t_model / ncoord_model[0];    // dr 2D for input model
    const double dcoord_z = 2. * globals::vmax * t_model / ncoord_model[1];  // dz 2D for input model
    return pow(globals::tmin / t_model, 3) * dcoord_z * PI *
           (pow((n_r + 1) * dcoord_rcyl, 2.) - pow(n_r * dcoord_rcyl, 2.));
  }
  if (get_model_type() == GridType::CARTESIAN3D) {
    // Assumes cells are cubes here - all same volume.
    return pow((2 * globals::vmax * globals::tmin), 3.) / (ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2]);
  }
  assert_always(false);
  return NAN;
}

void calc_modelinit_totmassradionuclides() {
  mtot_input = 0.;
  mfegroup = 0.;

  assert_always(totmassradionuclide.data() == nullptr);
  totmassradionuclide =
      std::span(static_cast<double *>(malloc(decay::get_num_nuclides() * sizeof(double))), decay::get_num_nuclides());
  assert_always(totmassradionuclide.data() != nullptr);

  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
    totmassradionuclide[nucindex] = 0.;
  }

  for (int mgi = 0; mgi < get_npts_model(); mgi++) {
    const double mass_in_shell = get_rho_tmin(mgi) * get_inputcellvolume(mgi);
    if (mass_in_shell > 0) {
      mtot_input += mass_in_shell;

      for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
        totmassradionuclide[nucindex] += mass_in_shell * get_modelinitnucmassfrac(mgi, nucindex);
      }

      mfegroup += mass_in_shell * get_ffegrp(mgi);
    }
  }
}

void read_grid_restart_data(const int timestep) {
  char filename[MAXFILENAMELENGTH];
  snprintf(filename, MAXFILENAMELENGTH, "gridsave_ts%d.tmp", timestep);

  printout("READIN GRID SNAPSHOT from %s\n", filename);
  FILE *gridsave_file = fopen_required(filename, "r");

  int ntimesteps_in = -1;
  assert_always(fscanf(gridsave_file, "%d ", &ntimesteps_in) == 1);
  assert_always(ntimesteps_in == globals::ntimesteps);

  int nprocs_in = -1;
  assert_always(fscanf(gridsave_file, "%d ", &nprocs_in) == 1);
  assert_always(nprocs_in == globals::nprocs);

  for (int nts = 0; nts < globals::ntimesteps; nts++) {
    int pellet_decays = 0.;
    assert_always(fscanf(gridsave_file,
                         "%la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %d ",
                         &globals::timesteps[nts].gamma_dep, &globals::timesteps[nts].gamma_dep_discrete,
                         &globals::timesteps[nts].positron_dep, &globals::timesteps[nts].positron_dep_discrete,
                         &globals::timesteps[nts].positron_emission, &globals::timesteps[nts].eps_positron_ana_power,
                         &globals::timesteps[nts].electron_dep, &globals::timesteps[nts].electron_dep_discrete,
                         &globals::timesteps[nts].electron_emission, &globals::timesteps[nts].eps_electron_ana_power,
                         &globals::timesteps[nts].alpha_dep, &globals::timesteps[nts].alpha_dep_discrete,
                         &globals::timesteps[nts].alpha_emission, &globals::timesteps[nts].eps_alpha_ana_power,
                         &globals::timesteps[nts].qdot_betaminus, &globals::timesteps[nts].qdot_alpha,
                         &globals::timesteps[nts].qdot_total, &globals::timesteps[nts].gamma_emission,
                         &globals::timesteps[nts].cmf_lum, &pellet_decays) == 20);
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
                         &globals::dep_estimator_alpha[nonemptymgi], &modelgrid[nonemptymgi].nne,
                         &modelgrid[nonemptymgi].nnetot) == 12);

    if (mgi_in != mgi) {
      printout("[fatal] read_grid_restart_data: cell mismatch in reading input gridsave.dat ... abort\n");
      printout("[fatal] read_grid_restart_data: read cellnumber %d, expected cellnumber %d\n", mgi_in, mgi);
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
    modelgrid[nonemptymgi].thick = thick;

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

// Assign temperatures to the grid cells at the start of the simulation
void assign_initial_temperatures() {
  MPI_Barrier(MPI_COMM_WORLD);  // For a simulation started from scratch we estimate the initial temperatures

  // We assume that for early times the material is so optically thick, that
  // all the radiation is trapped in the cell it originates from. This
  // means furthermore LTE, so that both temperatures can be evaluated
  // according to the local energy density resulting from the 56Ni decay.
  // The dilution factor is W=1 in LTE.

  printout("Assigning initial temperatures...\n");

  const double tstart = globals::timesteps[0].mid;
  int cells_below_mintemp = 0;
  int cells_above_maxtemp = 0;

  for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = get_mgi_of_nonemptymgi(nonemptymgi);

    double decayedenergy_per_mass = decay::get_endecay_per_ejectamass_t0_to_time_withexpansion(nonemptymgi, tstart);
    if constexpr (INITIAL_PACKETS_ON && USE_MODEL_INITIAL_ENERGY) {
      decayedenergy_per_mass += get_initenergyq(mgi);
    }

    double T_initial =
        pow(CLIGHT / 4 / STEBO * pow(globals::tmin / tstart, 3) * get_rho_tmin(mgi) * decayedenergy_per_mass, 1. / 4.);

    if (T_initial < MINTEMP) {
      //   printout("mgi %d: T_initial of %g is below MINTEMP %g K, setting to MINTEMP.\n", mgi, T_initial, MINTEMP);
      T_initial = MINTEMP;
      cells_below_mintemp++;
    } else if (T_initial > MAXTEMP) {
      //   printout("mgi %d: T_initial of %g is above MAXTEMP %g K, setting to MAXTEMP.\n", mgi, T_initial, MAXTEMP);
      T_initial = MAXTEMP;
      cells_above_maxtemp++;
    } else if (!std::isfinite(T_initial)) {
      printout("mgi %d: T_initial of %g is infinite!\n", mgi, T_initial);
    }
    assert_always(std::isfinite(T_initial));

    set_Te(nonemptymgi, T_initial);
    set_TJ(nonemptymgi, T_initial);
    set_TR(nonemptymgi, T_initial);
    set_W(nonemptymgi, 1.);
    modelgrid[nonemptymgi].thick = 0;
  }
  printout("  cells below MINTEMP %g: %d\n", MINTEMP, cells_below_mintemp);
  printout("  cells above MAXTEMP %g: %d\n", MAXTEMP, cells_above_maxtemp);
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
  const int min_nonempty_perproc =
      nonempty_npts_model / nprocesses;  // integer division, minimum non-empty cells per process
  const int n_remainder = nonempty_npts_model % nprocesses;

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
      const int target_nonempty_thisrank = (rank < n_remainder) ? min_nonempty_perproc + 1 : min_nonempty_perproc;
      if ((rank < (nprocesses - 1)) && (ranks_ndo_nonempty[rank] >= target_nonempty_thisrank)) {
        // current rank has enough non-empty cells, so start assigning cells to the next rank
        rank++;
        ranks_nstart[rank] = mgi;
        ranks_nstart_nonempty[rank] = get_next_nonemptymgi(mgi);
        assert_always(ranks_nstart_nonempty[rank] >= 0);
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
      fileout << r << " " << ranks_nstart[r] << " " << ranks_ndo[r] << " " << ranks_ndo_nonempty[r] << "\n";
    }
  }
}

// set up a uniform cuboidal grid.
void setup_grid_cartesian_3d() {
  // vmax is per coordinate, but the simulation volume corners will
  // have a higher expansion velocity than the sides
  const double vmax_corner = sqrt(3 * pow(globals::vmax, 2));
  printout("corner vmax %g [cm/s] (%.2fc)\n", vmax_corner, vmax_corner / CLIGHT);
  if (!FORCE_SPHERICAL_ESCAPE_SURFACE) {
    assert_always(vmax_corner < CLIGHT);
  }

  // Set grid size for uniform xyz grid
  if (get_model_type() == GridType::CARTESIAN3D) {
    // if we used in a 3D ejecta model, the propagation grid must match the input grid exactly
    // in case the user specified a grid size, we should ensure that it matches
    assert_always(ncoordgrid[0] == CUBOID_NCOORDGRID_X || CUBOID_NCOORDGRID_X < 0);
    assert_always(ncoordgrid[1] == CUBOID_NCOORDGRID_Y || CUBOID_NCOORDGRID_Y < 0);
    assert_always(ncoordgrid[2] == CUBOID_NCOORDGRID_Z || CUBOID_NCOORDGRID_Z < 0);
  } else {
    ncoordgrid = {CUBOID_NCOORDGRID_X, CUBOID_NCOORDGRID_Y, CUBOID_NCOORDGRID_Z};
  }

  // artis assumes in some places that the cells are cubes, not cubioids
  assert_always(ncoordgrid[0] == ncoordgrid[1]);
  assert_always(ncoordgrid[0] == ncoordgrid[2]);

  ngrid = ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2];
  resize_exactly(propcell_pos_min, ngrid);

  coordlabel = {'X', 'Y', 'Z'};
  std::array<int, 3> nxyz = {0, 0, 0};
  for (int n = 0; n < ngrid; n++) {
    for (int axis = 0; axis < 3; axis++) {
      assert_always(nxyz[axis] == get_cellcoordpointnum(n, axis));
      propcell_pos_min[n][axis] = -globals::rmax + (2 * nxyz[axis] * globals::rmax / ncoordgrid[axis]);
    }

    assert_always(n == nxyz[2] * ncoordgrid[1] * ncoordgrid[2] + nxyz[1] * ncoordgrid[0] + nxyz[0]);

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
  assert_always(get_model_type() == GridType::SPHERICAL1D);
  coordlabel = {'r', '_', '_'};

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
  printout("corner vmax %g [cm/s] (%.2fc)\n", vmax_corner, vmax_corner / CLIGHT);
  assert_always(vmax_corner < CLIGHT);

  assert_always(get_model_type() == GridType::CYLINDRICAL2D);
  coordlabel = {'r', 'z', '_'};

  ncoordgrid = ncoord_model;

  ngrid = ncoordgrid[0] * ncoordgrid[1];
  assert_always(ngrid == get_npts_model());

  resize_exactly(propcell_pos_min, ngrid);

  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const int n_rcyl = get_cellcoordpointnum(cellindex, 0);
    const int n_z = get_cellcoordpointnum(cellindex, 1);

    propcell_pos_min[cellindex] = {n_rcyl * globals::rmax / ncoord_model[0],
                                   globals::rmax * (-1 + n_z * 2. / ncoord_model[1]), 0.};
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

// Get the discrete index of the coordinate value (where pos must be position in grid coordinate system, not necessarily
// xyz)
auto get_poscoordpointnum(const double pos, const double time, const int axis) -> int {
  if constexpr (GRID_TYPE == GridType::CARTESIAN3D) {
    return static_cast<int>((pos / time + globals::vmax) / 2 / globals::vmax * ncoordgrid[axis]);
  } else if constexpr (GRID_TYPE == GridType::CYLINDRICAL2D) {
    if (axis == 0) {
      return static_cast<int>(pos / time / globals::vmax * ncoordgrid[axis]);
    }
    if (axis == 1) {
      return static_cast<int>((pos / time + globals::vmax) / 2 / globals::vmax * ncoordgrid[axis]);
    }
    assert_always(false);

  } else if constexpr (GRID_TYPE == GridType::SPHERICAL1D) {
    for (int n_r = 0; n_r < ncoordgrid[0]; n_r++) {
      if ((pos >= grid::get_cellcoordmin(n_r, 0)) && (pos < grid::get_cellcoordmax(n_r, 0))) {
        return n_r;
      }
    }
  }

  assert_always(false);
  return -1;
}

// Convert a position in Cartesian xyz to the grid coordinate system (which might the same, or 2D cylindrical or 1D
// spherical)
[[nodiscard]] constexpr auto get_gridcoords_from_xyz(const std::array<double, 3> &pos_xyz) {
  if constexpr (GRID_TYPE == GridType::CARTESIAN3D) {
    return pos_xyz;
  }

  if constexpr (GRID_TYPE == GridType::CYLINDRICAL2D) {
    return std::array<double, 2>{std::sqrt(std::pow(pos_xyz[0], 2) + std::pow(pos_xyz[1], 2)), pos_xyz[2]};
  }

  if constexpr (GRID_TYPE == GridType::SPHERICAL1D) {
    return std::array<double, 1>{vec_len(pos_xyz)};
  }

  assert_always(false);
}

// get the velocity in the grid coordinate system from the xyz position and direction
[[nodiscard]] constexpr auto get_gridcoords_vel_from_xyz_pos_dir(
    const std::array<double, 3> &pos_xyz, const std::array<double, 3> &dir_xyz,
    const std::array<double, get_ndim(GRID_TYPE)> &pktposgridcoord) {
  if constexpr (GRID_TYPE == GridType::CARTESIAN3D) {
    // keep xyz Cartesian coordinates
    return std::array<double, 3>{dir_xyz[0] * CLIGHT_PROP, dir_xyz[1] * CLIGHT_PROP, dir_xyz[2] * CLIGHT_PROP};
  } else if constexpr (GRID_TYPE == GridType::CYLINDRICAL2D) {
    // xy plane radial velocity
    // z velocity
    return std::array<double, 2>{(pos_xyz[0] * dir_xyz[0] + pos_xyz[1] * dir_xyz[1]) / pktposgridcoord[0] * CLIGHT_PROP,
                                 dir_xyz[2] * CLIGHT_PROP};

  } else if constexpr (GRID_TYPE == GridType::SPHERICAL1D) {
    // the only coordinate is radius from the origin
    return std::array<double, 1>{dot(pos_xyz, dir_xyz) / pktposgridcoord[0] * CLIGHT_PROP};
  } else {
    assert_always(false);
  }
}

// find the closest forward distance to the intersection of a ray with an expanding spherical shell (pos and dir are
// 2-vectors or 3-vectors) or expanding circle (2D vectors)
// returns -1 if there are no forward intersections (or if the intersection
// is tangential to the shell)
template <size_t S1>
[[nodiscard]] constexpr auto expanding_shell_intersection(const std::array<double, S1> &pos,
                                                          const std::array<double, S1> &dir, const double speed,
                                                          const double shellradiuststart, const bool isinnerboundary,
                                                          const double tstart) -> double {
  static_assert(S1 == 2 || S1 == 3);
  assert_always(shellradiuststart > 0);

  // quadratic equation for intersection of ray with sphere
  // a*d^2 + b*d + c = 0
  const double a = dot(dir, dir) - pow(shellradiuststart / tstart / speed, 2);
  const double b = 2 * (dot(dir, pos) - pow(shellradiuststart, 2) / tstart / speed);
  const double c = dot(pos, pos) - pow(shellradiuststart, 2);

  const double discriminant = pow(b, 2) - (4 * a * c);

  if (discriminant < 0) {
    // no intersection
    assert_always(isinnerboundary);
    assert_always(shellradiuststart < vec_len(pos));
    return -1;
  }

  if (discriminant > 0) {
    // two intersections
    double dist1 = (-b + sqrt(discriminant)) / 2 / a;
    double dist2 = (-b - sqrt(discriminant)) / 2 / a;

    auto posfinal1 = std::array<double, 3>{0.};
    auto posfinal2 = std::array<double, 3>{0.};

    for (int d = 0; d < std::ssize(pos); d++) {
      posfinal1[d] = pos[d] + dist1 * dir[d];
      posfinal2[d] = pos[d] + dist2 * dir[d];
    }

    const double v_rad_shell = shellradiuststart / tstart;
    const double v_rad_final1 = dot(dir, posfinal1) * speed / vec_len(posfinal1);
    const double v_rad_final2 = dot(dir, posfinal2) * speed / vec_len(posfinal2);

    // invalidate any solutions that require entering the boundary from the wrong radial direction
    if (isinnerboundary) {
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
      const double shellradiusfinal1 = shellradiuststart / tstart * (tstart + dist1 / speed);
      assert_testmodeonly(fabs((vec_len(posfinal1) / shellradiusfinal1) - 1.) < 1e-3);
    }
    if (dist2 >= 0) {
      const double shellradiusfinal2 = shellradiuststart / tstart * (tstart + dist2 / speed);
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
  assert_always(shellradiuststart <= vec_len(pos));
  return -1.;
}

}  // anonymous namespace

// for a uniform grid get the the extent along the x,y,z coordinate (x_2 - x_1, etc.) at time tmin
// for spherical grid get the radial extent (r_outer - r_inner) at time tmin
auto wid_init(const int cellindex, const int axis) -> double {
  if constexpr (GRID_TYPE == GridType::CARTESIAN3D) {
    return 2 * globals::rmax / ncoordgrid[axis];
  }

  if constexpr (GRID_TYPE == GridType::CYLINDRICAL2D) {
    return (axis == 0) ? globals::rmax / ncoordgrid[axis] : 2 * globals::rmax / ncoordgrid[axis];
  }

  if constexpr (GRID_TYPE == GridType::SPHERICAL1D) {
    const int modelgridindex = get_propcell_modelgridindex(cellindex);
    const double v_inner = modelgridindex > 0 ? vout_model[modelgridindex - 1] : 0.;
    return (vout_model[modelgridindex] - v_inner) * globals::tmin;
  }

  assert_always(false);
}

// return the model cell volume (when mapped to the propagation cells) at globals::tmin
// for a uniform cubic grid this is constant
auto get_modelcell_assocvolume_tmin(const int modelgridindex) -> double {
  if constexpr (GRID_TYPE == GridType::CARTESIAN3D) {
    return (wid_init(modelgridindex, 0) * wid_init(modelgridindex, 1) * wid_init(modelgridindex, 2)) *
           get_numpropcells(modelgridindex);
  }

  if constexpr (GRID_TYPE == GridType::CYLINDRICAL2D) {
    return wid_init(modelgridindex, 1) * PI *
           (pow(get_cellcoordmax(modelgridindex, 0), 2) - pow(get_cellcoordmin(modelgridindex, 0), 2));
  }

  if constexpr (GRID_TYPE == GridType::SPHERICAL1D) {
    return 4. / 3. * PI * (pow(get_cellcoordmax(modelgridindex, 0), 3) - pow(get_cellcoordmin(modelgridindex, 0), 3));
  }

  assert_always(false);
}

// return the propagation cell volume at globals::tmin
// for a spherical grid, the cell index is required (and should be equivalent to a modelgridindex)
auto get_propcell_volume_tmin(const int cellindex) -> double {
  if constexpr (GRID_TYPE == GridType::CARTESIAN3D) {
    return (wid_init(cellindex, 0) * wid_init(cellindex, 0) * wid_init(cellindex, 0));
  }

  // 2D and 1D with direct mapping to propagation cells
  const int mgi = get_propcell_modelgridindex(cellindex);
  return get_modelcell_assocvolume_tmin(mgi);
}

// get the minimum value of a coordinate at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the minimum x position in xyz coords, or the minimum radius
auto get_cellcoordmax(const int cellindex, const int axis) -> double {
  if constexpr (GRID_TYPE == GridType::CARTESIAN3D) {
    return grid::get_cellcoordmin(cellindex, axis) + grid::wid_init(0, axis);
  }

  if constexpr (GRID_TYPE == GridType::CYLINDRICAL2D) {
    assert_testmodeonly(axis <= 1);
    return grid::get_cellcoordmin(cellindex, axis) + grid::wid_init(cellindex, axis);
  }

  if constexpr (GRID_TYPE == GridType::SPHERICAL1D) {
    assert_testmodeonly(axis == 0);
    return grid::get_cellcoordmin(cellindex, axis) + grid::wid_init(cellindex, axis);
  }

  assert_always(false);
}

// get the minimum value of a coordinate at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the minimum x position in xyz coords, or the minimum radius
auto get_cellcoordmin(const int cellindex, const int axis) -> double {
  return propcell_pos_min[cellindex][axis];
  // return - coordmax[axis] + (2 * get_cellcoordpointnum(cellindex, axis) * coordmax[axis] / ncoordgrid[axis]);
}

// how much do we change the cellindex to move along a coordinately axis (e.g., the x, y, z directions, or r
// direction)
auto get_coordcellindexincrement(const int axis) -> int {
  switch (axis) {
    case 0:
      return 1;

    case 1:
      return ncoordgrid[0];

    case 2:
      return ncoordgrid[0] * ncoordgrid[1];

    default:
      if constexpr (TESTMODE) {
        printout("invalid coordinate index %d", axis);
        assert_testmodeonly(false);
      } else {
        std::unreachable();
      }
  }
}

// convert a cell index number into an integer (x,y,z or r) coordinate index from 0 to ncoordgrid[axis]
auto get_cellcoordpointnum(const int cellindex, const int axis) -> int {
  if constexpr (GRID_TYPE == GridType::CARTESIAN3D || GRID_TYPE == GridType::CYLINDRICAL2D) {
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
          printout("invalid coordinate index %d", axis);
          assert_testmodeonly(false);
        } else {
          std::unreachable();
        }
    }
  }

  if constexpr (GRID_TYPE == GridType::SPHERICAL1D) {
    return cellindex;
  }

  assert_always(false);
}

auto get_rho_tmin(const int modelgridindex) -> float { return modelgrid_input[modelgridindex].rhoinit; }

__host__ __device__ auto get_rho(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return modelgrid[nonemptymgi].rho;
}

__host__ __device__ auto get_nne(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  const double nne = modelgrid[nonemptymgi].nne;
  assert_testmodeonly(std::isfinite(nne));
  return nne;
}

__host__ __device__ auto get_nnetot(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  const double nnetot = modelgrid[nonemptymgi].nnetot;
  return nnetot;
}

__host__ __device__ auto get_ffegrp(const int modelgridindex) -> float {
  return modelgrid_input[modelgridindex].ffegrp;
}

__host__ __device__ auto get_initial_radial_pos_sum(const int modelgridindex) -> float {
  return modelgrid_input[modelgridindex].initial_radial_pos_sum;
}

auto get_elem_abundance(int nonemptymgi, int element) -> float
// mass fraction of an element (all isotopes combined)
{
  const auto massfrac = elem_massfracs_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_nelements()) + element];
  assert_testmodeonly(massfrac >= 0.0);
  return massfrac;
}

// mass fraction of an element (all isotopes combined)
void set_elem_abundance(const int nonemptymgi, const int element, const float newabundance) {
  elem_massfracs_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_nelements()) + element] = newabundance;
}

// mass fraction of an element (all isotopes combined)
__host__ __device__ auto get_elem_numberdens(const int nonemptymgi, const int element) -> double {
  const double elem_meanweight = grid::get_element_meanweight(nonemptymgi, element);
  return get_elem_abundance(nonemptymgi, element) / elem_meanweight * grid::get_rho(nonemptymgi);
}

__host__ __device__ auto get_kappagrey(const int nonemptymgi) -> float { return modelgrid[nonemptymgi].kappagrey; }

__host__ __device__ auto get_Te(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return modelgrid[nonemptymgi].Te;
}

__host__ __device__ auto get_TR(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return modelgrid[nonemptymgi].TR;
}

__host__ __device__ auto get_TJ(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return modelgrid[nonemptymgi].TJ;
}

__host__ __device__ auto get_W(const int nonemptymgi) -> float {
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return modelgrid[nonemptymgi].W;
}

void set_rho(const int nonemptymgi, const float rho) {
  assert_always(rho >= 0.);
  assert_always(std::isfinite(rho));
  modelgrid[nonemptymgi].rho = rho;
}

void set_nne(const int nonemptymgi, const float nne) {
  assert_always(nne >= 0.);
  assert_always(std::isfinite(nne));
  modelgrid[nonemptymgi].nne = nne;
}

void set_nnetot(const int nonemptymgi, const float nnetot) {
  assert_always(nnetot >= 0.);
  assert_always(std::isfinite(nnetot));
  modelgrid[nonemptymgi].nnetot = nnetot;
}

void set_kappagrey(const int nonemptymgi, const float kappagrey) { modelgrid[nonemptymgi].kappagrey = kappagrey; }

void set_Te(const int nonemptymgi, const float Te) {
  if (Te > 0.) {
    // ignore the zero initialisation value for this check
    const double nu_peak = 5.879e10 * Te;
    if (nu_peak > NU_MAX_R || nu_peak < NU_MIN_R) {
      const auto modelgridindex = get_mgi_of_nonemptymgi(nonemptymgi);
      printout(
          "[warning] modelgridindex %d B_planck(Te=%g K) peak at %g Hz is outside frequency range NU_MIN_R %g NU_MAX_R "
          "%g\n",
          modelgridindex, Te, nu_peak, NU_MIN_R, NU_MAX_R);
    }
  }

  modelgrid[nonemptymgi].Te = Te;
}

void set_TR(const int nonemptymgi, const float TR) { modelgrid[nonemptymgi].TR = TR; }

void set_TJ(const int nonemptymgi, const float TJ) { modelgrid[nonemptymgi].TJ = TJ; }

void set_W(const int nonemptymgi, const float W) { modelgrid[nonemptymgi].W = W; }

auto get_model_type() -> GridType { return model_type; }

void set_model_type(const GridType model_type_value) { model_type = model_type_value; }

__host__ __device__ auto get_npts_model() -> int
// number of model grid cells
{
  assert_testmodeonly(npts_model > 0);
  return npts_model;
}

// number of model grid cells
auto get_nonempty_npts_model() -> int {
  assert_testmodeonly(nonempty_npts_model > 0);
  return nonempty_npts_model;
}

// get time at which model input densities are defined
auto get_t_model() -> double {
  assert_testmodeonly(t_model > 0.);
  return t_model;
}

[[nodiscard]] __host__ __device__ auto get_propcell_modelgridindex(const int cellindex) -> int {
  assert_testmodeonly(cellindex >= 0);
  assert_testmodeonly(cellindex < ngrid);
  const auto mgi = propcell_mgi[cellindex];
  assert_testmodeonly(mgi >= 0);
  assert_testmodeonly(mgi < (get_npts_model() + 1));
  return mgi;
}

[[nodiscard]] __host__ __device__ auto get_propcell_nonemptymgi(const int cellindex) -> int {
  const auto nonemptymgi = propcell_nonemptymgi[cellindex];
  assert_testmodeonly(nonemptymgi >= -1);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());
  return nonemptymgi;
}

// number of propagation cells associated with each modelgrid cell
__host__ __device__ auto get_numpropcells(const int modelgridindex) -> int {
  assert_testmodeonly(modelgridindex <= get_npts_model());
  return modelgrid_numpropcells[modelgridindex];
}

// get the index in the list of non-empty cells for a given model grid cell
__host__ __device__ auto get_nonemptymgi_of_mgi(const int mgi) -> int {
  assert_testmodeonly(get_nonempty_npts_model() > 0);
  assert_testmodeonly(mgi < get_npts_model());

  const int nonemptymgi = nonemptymgi_of_mgi[mgi];
  // assert_testmodeonly(nonemptymgi >= 0 || get_numpropcells(mgi) == 0);
  assert_testmodeonly(nonemptymgi >= 0);
  assert_testmodeonly(nonemptymgi < get_nonempty_npts_model());

  return nonemptymgi;
}

// get the index in the list of non-empty cells for a given model grid cell
__host__ __device__ auto get_mgi_of_nonemptymgi(const int nonemptymgi) -> int {
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
  const ptrdiff_t num_nuclides = decay::get_num_nuclides();

  return initnucmassfrac_allcells[(modelgridindex * num_nuclides) + nucindex];
}

auto get_stable_initabund(const int nonemptymgi, const int element) -> float {
  return initmassfracuntrackedstable_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_nelements()) + element];
}

auto get_element_meanweight(const int nonemptymgi, const int element) -> float
// weight is in grams
{
  if (USE_CALCULATED_MEANATOMICWEIGHT) {
    const double mu = elem_meanweight_allcells[(nonemptymgi * get_nelements()) + element];
    if (mu > 0) {
      return mu;
    }
  }
  return globals::elements[element].initstablemeannucmass;
}

// set element weight in grams
void set_element_meanweight(const int nonemptymgi, const int element, const float meanweight) {
  assert_always(meanweight > 0.);
  elem_meanweight_allcells[(static_cast<ptrdiff_t>(nonemptymgi) * get_nelements()) + element] = meanweight;
}

auto get_electronfrac(const int nonemptymgi) -> double {
  double nucleondens = 0.;
  for (int element = 0; element < get_nelements(); element++) {
    nucleondens += get_elem_numberdens(nonemptymgi, element) * get_element_meanweight(nonemptymgi, element) / MH;
  }
  return get_nnetot(nonemptymgi) / nucleondens;
}

auto get_initelectronfrac(const int modelgridindex) -> double {
  return modelgrid_input[modelgridindex].initelectronfrac;
}

// q: energy in the model at tmin per gram to use with USE_MODEL_INITIAL_ENERGY option [erg/g]
auto get_initenergyq(const int modelgridindex) -> double { return modelgrid_input[modelgridindex].initenergyq; }

// get the radial distance from the origin to the centre of the cell at time tmin
auto get_cellradialposmid(const int cellindex) -> double {
  if (GRID_TYPE == GridType::SPHERICAL1D) {
    // mid point radius
    // return get_cellcoordmin(cellindex, 0) + (0.5 * wid_init(cellindex, 0));
    // volume averaged mean radius is slightly complex for radial shells
    const double r_inner = grid::get_cellcoordmin(cellindex, 0);
    const double r_outer = r_inner + grid::wid_init(cellindex, 0);
    return 3. / 4 * (pow(r_outer, 4.) - pow(r_inner, 4.)) / (pow(r_outer, 3) - pow(r_inner, 3.));
  }

  if (GRID_TYPE == GridType::CYLINDRICAL2D) {
    const double rcyl_mid = get_cellcoordmin(cellindex, 0) + (0.5 * wid_init(cellindex, 0));
    const double z_mid = get_cellcoordmin(cellindex, 1) + (0.5 * wid_init(cellindex, 1));
    return std::sqrt(std::pow(rcyl_mid, 2) + std::pow(z_mid, 2));
  }

  // cubic grid requires taking the length of the 3D position vector
  std::array<double, 3> dcen{};
  for (int axis = 0; axis < 3; axis++) {
    dcen[axis] = get_cellcoordmin(cellindex, axis) + (0.5 * wid_init(cellindex, axis));
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
        double kappagrey = ((0.9 * get_ffegrp(mgi)) + 0.1);

        if (get_rho_tmin(mgi) > globals::rho_crit) {
          kappagrey *= globals::rho_crit / get_rho_tmin(mgi);
        }

        set_kappagrey(nonemptymgi, kappagrey);
      } else if (get_rho_tmin(mgi) == 0.) {
        set_kappagrey(nonemptymgi, 0.);
      } else if (get_rho_tmin(mgi) < 0.) {
        printout("Error: negative density. Abort.\n");
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
        // const double Ye = get_electronfrac(mgi);
        const double Ye = get_initelectronfrac(mgi);
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
        printout("Unknown opacity case. Abort.\n");
        std::abort();
      }

      set_kappagrey(nonemptymgi, kappa);
    } else if (get_rho_tmin(mgi) == 0.) {
      set_kappagrey(nonemptymgi, 0.);
    } else if (get_rho_tmin(mgi) < 0.) {
      printout("Error: negative density. Abort.\n");
      std::abort();
    }

    check1 = check1 + (get_kappagrey(nonemptymgi) * get_rho_tmin(mgi));
    check2 = check2 + get_rho_tmin(mgi);
  }

  printout("Grey normalisation check: %g\n", check1 / check2);
}

void read_ejecta_model() {
  auto fmodel = fstream_required("model.txt", std::ios::in);
  std::string line;
  std::optional<GridType> detected_dim{};

  // two integers on the first line of the model file
  int npts_0 = 0;  // total model points for 1D/3D, and number of points in r for 2D
  int npts_1 = 0;  // number of points in z for 2D
  assert_always(get_noncommentline(fmodel, line));
  auto ssline = std::istringstream(line);
  ssline >> npts_0;
  if (ssline >> npts_1) {
    // second number on the line for 2D means the line was n_r n_z
    detected_dim = GridType::CYLINDRICAL2D;
    printout("Detected 2D model\n");
    ssline >> npts_1;  // r and z (cylindrical polar)
    npts_model = npts_0 * npts_1;
  } else {
    // for 1D and 3D, this was the total number of model cells
    npts_model = npts_0;
  }

  // Now read the time (in days) at which the model is specified.
  double t_model_days{NAN};
  assert_always(get_noncommentline(fmodel, line));
  std::istringstream(line) >> t_model_days;
  t_model = t_model_days * DAY;

  const auto pos_after_t_model = fmodel.tellg();
  // if the next line is a single float, it is the vmax (so 2D or 3D)
  // otherwise, it is the first line of the model or a header comment (so 1D)
  std::getline(fmodel, line);
  if (!line.starts_with("#")) {
    double num_after_vmax{NAN};
    auto sslinevmax = std::istringstream(line);
    if ((sslinevmax >> globals::vmax) && !(sslinevmax >> num_after_vmax)) {
      // single value on the line is a vmax, so 2D or 3D
      // if it's not already know to be 2D (based on n_r n_z line), then it's 3D
      if (detected_dim != GridType::CYLINDRICAL2D) {
        assert_always(!detected_dim.has_value());
        detected_dim = GridType::CARTESIAN3D;
        printout("Detected 3D model\n");
      }
    }
  }
  if (!detected_dim.has_value()) {
    assert_always(!detected_dim.has_value());
    detected_dim = GridType::SPHERICAL1D;
    printout("Detected 1D model\n");
    fmodel.seekg(pos_after_t_model);
  }

  set_model_type(detected_dim.value());

  assert_always(modelgrid_input.data() == nullptr);
  modelgrid_input = MPI_shared_malloc_span<ModelGridCellInput>(npts_model + 1);
  if (globals::rank_in_node == 0) {
    std::ranges::fill(modelgrid_input, ModelGridCellInput{});
  }
  MPI_Barrier(globals::mpi_comm_node);
  modelgrid_numpropcells.resize(npts_model + 1, 0);
  nonemptymgi_of_mgi.resize(npts_model + 1, -1);

  if (get_model_type() == GridType::SPHERICAL1D) {
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
      ssline = std::istringstream(line);

      if (ssline >> cellnumberin >> vout_kmps >> log_rho) {
        if (mgi == 0) {
          first_cellindex = cellnumberin;
          printout("first_cellindex %d\n", first_cellindex);
        }
        assert_always(cellnumberin == mgi + first_cellindex);

        vout_model[mgi] = vout_kmps * 1.e5;

        const double rho_tmin = pow(10., log_rho) * pow(t_model / globals::tmin, 3);
        set_rho_tmin(mgi, rho_tmin);
      } else {
        printout("Unexpected number of values in model.txt\n");
        printout("line: %s\n", line.c_str());
        assert_always(false);
      }
      read_model_radioabundances(fmodel, ssline, mgi, true, colnames, nucindexlist, one_line_per_cell);

      mgi += 1;
      if (mgi == get_npts_model()) {
        break;
      }
    }

    if (mgi != get_npts_model()) {
      printout("ERROR in model.txt. Found only %d cells instead of %d expected.\n", mgi - 1, get_npts_model());
      std::abort();
    }

    globals::vmax = vout_model[get_npts_model() - 1];
  } else if (get_model_type() == GridType::CYLINDRICAL2D) {
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
      ssline = std::istringstream(line);
      assert_always(ssline >> cellnumberin >> cell_r_in >> cell_z_in >> rho_tmodel);

      if (mgi == 0) {
        first_cellindex = cellnumberin;
      }
      assert_always(cellnumberin == mgi + first_cellindex);

      const int n_rcyl = (mgi % ncoord_model[0]);
      const double pos_r_cyl_mid = (n_rcyl + 0.5) * globals::vmax * t_model / ncoord_model[0];
      assert_always(fabs((cell_r_in / pos_r_cyl_mid) - 1) < 1e-3);
      const int n_z = (mgi / ncoord_model[0]);
      const double pos_z_mid = globals::vmax * t_model * (-1 + 2 * (n_z + 0.5) / ncoord_model[1]);
      assert_always(fabs((cell_z_in / pos_z_mid) - 1) < 1e-3);

      if (rho_tmodel < 0) {
        printout("negative input density %g %d\n", rho_tmodel, mgi);
        std::abort();
      }

      const bool keepcell = (rho_tmodel > 0);
      const double rho_tmin = rho_tmodel * pow(t_model / globals::tmin, 3);
      set_rho_tmin(mgi, rho_tmin);

      read_model_radioabundances(fmodel, ssline, mgi, keepcell, colnames, nucindexlist, one_line_per_cell);

      mgi++;
    }

    if (mgi != get_npts_model()) {
      printout("ERROR in model.txt. Found %d only cells instead of %d expected.\n", mgi - 1, get_npts_model());
      std::abort();
    }
  } else if (get_model_type() == GridType::CARTESIAN3D) {
    ncoord_model[0] = static_cast<int>(round(pow(npts_0, 1 / 3.)));
    ncoord_model[1] = ncoord_model[0];
    ncoord_model[2] = ncoord_model[0];
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
      ssline = std::istringstream(line);

      assert_always(ssline >> cellnumberin >> cellpos_in[0] >> cellpos_in[1] >> cellpos_in[2] >> rho_model);

      if (mgi == 0) {
        first_cellindex = cellnumberin;
      }
      assert_always(cellnumberin == mgi + first_cellindex);

      if (mgi % (ncoord_model[1] * ncoord_model[2]) == 0) {
        printout("read up to cell mgi %d\n", mgi);
      }

      // cell coordinates in the 3D model.txt file are sometimes reordered by the scaling script
      // however, the cellindex always should increment X first, then Y, then Z

      for (int axis = 0; axis < 3; axis++) {
        const double cellwidth = 2 * xmax_tmodel / ncoordgrid[axis];
        const double cellpos_expected = -xmax_tmodel + (cellwidth * get_cellcoordpointnum(mgi, axis));
        //   printout("mgi %d coord %d expected %g found %g or %g rmax %g get_cellcoordpointnum(mgi, axis) %d
        //   ncoordgrid %d\n",
        //            mgi, axis, cellpos_expected, cellpos_in[axis], cellpos_in[2 - axis], xmax_tmodel,
        //            get_cellcoordpointnum(mgi, axis), ncoordgrid[axis]);
        if (fabs(cellpos_expected - cellpos_in[axis]) > 0.5 * cellwidth) {
          posmatch_xyz = false;
        }
        if (fabs(cellpos_expected - cellpos_in[2 - axis]) > 0.5 * cellwidth) {
          posmatch_zyx = false;
        }
      }

      if (rho_model < 0) {
        printout("negative input density %g %d\n", rho_model, mgi);
        std::abort();
      }

      // in 3D cartesian, cellindex and modelgridindex are interchangeable
      const bool keepcell = (rho_model > 0);
      const double rho_tmin = rho_model * pow(t_model / globals::tmin, 3);
      set_rho_tmin(mgi, rho_tmin);

      if (min_den < 0. || min_den > rho_model) {
        min_den = rho_model;
      }

      read_model_radioabundances(fmodel, ssline, mgi, keepcell, colnames, nucindexlist, one_line_per_cell);

      mgi++;
    }
    if (mgi != npts_model) {
      printout("ERROR in model.txt. Found %d cells instead of %td expected.\n", mgi, npts_model);
      std::abort();
    }

    //   assert_always(posmatch_zyx ^ posmatch_xyz);  // xor because if both match then probably an infinity occurred
    if (posmatch_xyz) {
      printout("Cell positions in model.txt are consistent with calculated values when x-y-z column order is used.\n");
    }
    if (posmatch_zyx) {
      printout("Cell positions in model.txt are consistent with calculated values when z-y-x column order is used.\n");
    }

    if (!posmatch_xyz && !posmatch_zyx) {
      printout(
          "WARNING: Cell positions in model.txt are not consistent with calculated values in either x-y-z or z-y-x "
          "order.\n");
    }

    printout("min_den %g [g/cm3]\n", min_den);
  }

  assert_always(get_npts_model() ==
                std::max(1, ncoord_model[0]) * std::max(1, ncoord_model[1]) * std::max(1, ncoord_model[2]));
  printout("npts_model: %d\n", get_npts_model());
  globals::rmax = globals::vmax * globals::tmin;
  printout("vmax %g [cm/s] (%.2fc)\n", globals::vmax, globals::vmax / CLIGHT);
  assert_always(globals::vmax < CLIGHT);
  printout("tmin %g [s] = %.2f [d]\n", globals::tmin, globals::tmin / 86400.);
  printout("rmax %g [cm] (at t=tmin)\n", globals::rmax);

  calc_modelinit_totmassradionuclides();

  printout("Total input model mass: %9.3e [Msun]\n", mtot_input / MSUN);
  printout("Nuclide masses at t=t_model_init [Msun]:");
  printout("  56Ni: %9.3e  56Co: %9.3e  52Fe: %9.3e  48Cr: %9.3e\n", get_totmassradionuclide(28, 56) / MSUN,
           get_totmassradionuclide(27, 56) / MSUN, get_totmassradionuclide(26, 52) / MSUN,
           get_totmassradionuclide(24, 48) / MSUN);
  printout("  Fe-group: %9.3e  57Ni: %9.3e  57Co: %9.3e\n", mfegroup / MSUN, get_totmassradionuclide(28, 57) / MSUN,
           get_totmassradionuclide(27, 57) / MSUN);

  read_possible_yefile();
}

void write_grid_restart_data(const int timestep) {
  char filename[MAXFILENAMELENGTH];
  snprintf(filename, MAXFILENAMELENGTH, "gridsave_ts%d.tmp", timestep);

  const auto sys_time_start_write_restart = std::time(nullptr);
  printout("Write grid restart data to %s...", filename);

  FILE *gridsave_file = fopen_required(filename, "w");

  fprintf(gridsave_file, "%d ", globals::ntimesteps);
  fprintf(gridsave_file, "%d ", globals::nprocs);

  for (int nts = 0; nts < globals::ntimesteps; nts++) {
    fprintf(gridsave_file, "%la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %la %d ",
            globals::timesteps[nts].gamma_dep, globals::timesteps[nts].gamma_dep_discrete,
            globals::timesteps[nts].positron_emission, globals::timesteps[nts].positron_dep,
            globals::timesteps[nts].positron_dep_discrete, globals::timesteps[nts].eps_positron_ana_power,
            globals::timesteps[nts].electron_dep, globals::timesteps[nts].electron_dep_discrete,
            globals::timesteps[nts].electron_emission, globals::timesteps[nts].eps_electron_ana_power,
            globals::timesteps[nts].alpha_dep, globals::timesteps[nts].alpha_dep_discrete,
            globals::timesteps[nts].alpha_emission, globals::timesteps[nts].eps_alpha_ana_power,
            globals::timesteps[nts].qdot_betaminus, globals::timesteps[nts].qdot_alpha,
            globals::timesteps[nts].qdot_total, globals::timesteps[nts].gamma_emission, globals::timesteps[nts].cmf_lum,
            globals::timesteps[nts].pellet_decays);
  }

  fprintf(gridsave_file, "%d ", timestep);

  for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);

    assert_always(globals::dep_estimator_gamma[nonemptymgi] >= 0.);
    fprintf(gridsave_file, "%d %a %a %a %a %d %la %la %la %la %a %a", mgi, get_TR(nonemptymgi), get_Te(nonemptymgi),
            get_W(nonemptymgi), get_TJ(nonemptymgi), modelgrid[nonemptymgi].thick,
            globals::dep_estimator_gamma[nonemptymgi], globals::dep_estimator_positron[nonemptymgi],
            globals::dep_estimator_electron[nonemptymgi], globals::dep_estimator_alpha[nonemptymgi],
            modelgrid[nonemptymgi].nne, modelgrid[nonemptymgi].nnetot);

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
  printout("done in %ld seconds.\n", std::time(nullptr) - sys_time_start_write_restart);
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
void grid_init(const int my_rank) {
  // The cells will be ordered by x then y, then z. Call a routine that
  // sets up the initial positions and widths of the cells.

  if (GRID_TYPE == GridType::CARTESIAN3D) {
    setup_grid_cartesian_3d();
  } else if (GRID_TYPE == GridType::CYLINDRICAL2D) {
    setup_grid_cylindrical_2d();
  } else if (GRID_TYPE == GridType::SPHERICAL1D) {
    setup_grid_spherical_1d();
  } else {
    printout("[fatal] grid_init: Error: Unknown grid type. Abort.");
    std::abort();
  }
  propcell_mgi.resize(ngrid, -1);

  printout("propagation grid: %d-dimensional %s\n", get_ndim(GRID_TYPE), get_grid_type_name(GRID_TYPE).c_str());

  for (int d = 0; d < get_ndim(GRID_TYPE); d++) {
    printout("    coordinate %d '%c': cells have %d position values\n", d, coordlabel[d], ncoordgrid[d]);
  }
  printout("    total propagation cells: %d\n", ngrid);

  // Now set up the density in each cell.

  // Calculate the critical opacity at which opacity_case 3 switches from a
  // regime proportional to the density to a regime independent of the density
  // This is done by solving for tau_sobolev == 1
  // tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/nucmass(28, 56) * 3000e-8 * globals::timesteps[m].mid;
  globals::rho_crit =
      ME * CLIGHT * decay::nucmass(28, 56) / (PI * QE * QE * globals::rho_crit_para * 3000e-8 * globals::tmin);
  printout("grid_init: rho_crit = %g [g/cm3]\n", globals::rho_crit);

  if (get_model_type() == GRID_TYPE) {
    if (get_model_type() == GridType::CARTESIAN3D) {
      assert_always(ncoord_model[0] == ncoordgrid[0]);
      assert_always(ncoord_model[1] == ncoordgrid[1]);
      assert_always(ncoord_model[2] == ncoordgrid[2]);
    }

    map_modeltogrid_direct();
  } else if (get_model_type() == GridType::SPHERICAL1D) {
    assert_always(GRID_TYPE == GridType::CARTESIAN3D);
    map_1dmodelto3dgrid();
  } else if (get_model_type() == GridType::CYLINDRICAL2D) {
    assert_always(GRID_TYPE == GridType::CARTESIAN3D);
    map_2dmodelto3dgrid();
  } else {
    printout("[fatal] grid_init: Error: Unknown density type. Abort.");
    std::abort();
  }

  if (globals::my_rank == 0) {
    auto grid_file = std::fstream("grid.out", std::ios::out);
    assert_always(grid_file.is_open());
    for (int n = 0; n < ngrid; n++) {
      const int mgi = get_propcell_modelgridindex(n);
      if (mgi != get_npts_model()) {
        grid_file << n << " " << mgi << "\n";  // write only non-empty cells to grid file
      }
    }
  }

  allocate_nonemptymodelcells();
  calculate_kappagrey();
  abundances_read();

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
  if (GRID_TYPE == GridType::CARTESIAN3D && get_model_type() == GridType::SPHERICAL1D && globals::rank_in_node == 0) {
    for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
      if (totmassradionuclide[nucindex] <= 0) {
        continue;
      }

      double totmassradionuclide_actual = 0.;
      for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
        const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
        totmassradionuclide_actual +=
            get_modelinitnucmassfrac(mgi, nucindex) * get_rho_tmin(mgi) * get_modelcell_assocvolume_tmin(mgi);
      }

      if (totmassradionuclide_actual > 0.) {
        const double ratio = totmassradionuclide[nucindex] / totmassradionuclide_actual;
        for (int nonemptymgi = 0; nonemptymgi < get_nonempty_npts_model(); nonemptymgi++) {
          const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
          const double prev_abund = get_modelinitnucmassfrac(mgi, nucindex);
          const double new_abund = prev_abund * ratio;
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
  printout("Total grid-mapped mass: %9.3e [Msun] (%.1f%% of input mass)\n", mtot_mapped / MSUN,
           mtot_mapped / mtot_input * 100.);

  MPI_Barrier(MPI_COMM_WORLD);
}

auto get_totmassradionuclide(const int z, const int a) -> double {
  return totmassradionuclide[decay::get_nucindex(z, a)];
}

// identify the cell index from an (x,y,z) position and a time.
[[nodiscard]] auto get_cellindex_from_pos(const std::array<double, 3> &pos, const double time) -> int {
  auto posgridcoords = get_gridcoords_from_xyz(pos);
  int cellindex = 0;
  for (int d = 0; d < get_ndim(GRID_TYPE); d++) {
    cellindex += get_coordcellindexincrement(d) * get_poscoordpointnum(posgridcoords[d], time, d);
  }

  // do a check that the position is within the cell
  const double trat = time / globals::tmin;
  for (int n = 0; n < get_ndim(GRID_TYPE); n++) {
    assert_always(posgridcoords[n] >= grid::get_cellcoordmin(cellindex, n) * trat);
    assert_always(posgridcoords[n] <= grid::get_cellcoordmax(cellindex, n) * trat);
  }
  return cellindex;
}

// compute distance to a cell boundary.
[[nodiscard]] __host__ __device__ auto boundary_distance(const std::array<double, 3> &dir,
                                                         const std::array<double, 3> &pos, const double tstart,
                                                         const int cellindex) -> std::tuple<double, int> {
  if constexpr (FORCE_SPHERICAL_ESCAPE_SURFACE) {
    if (get_cell_r_inner(cellindex) > globals::vmax * globals::tmin) {
      return {0., -99};
    }
  }

  // d is used to loop over the coordinate indicies 0,1,2 for x,y,z

  // the following vector are in grid coordinates, so either x,y,z (3D) or r (1D), or r_xy, z (2D)
  static_assert(get_ndim(GRID_TYPE) <= 3);

  const auto pktposgridcoord = get_gridcoords_from_xyz(pos);

  // dir * CLIGHT_PROP converted from xyz to grid coordinates
  const auto pktvelgridcoord = get_gridcoords_vel_from_xyz_pos_dir(pos, dir, pktposgridcoord);

  const auto cellcoordmax = [cellindex] {
    auto _cellcoordmax = std::array<double, get_ndim(GRID_TYPE)>{};  // position at time tmin
    for (int d = 0; d < get_ndim(GRID_TYPE); d++) {
      _cellcoordmax[d] = grid::get_cellcoordmax(cellindex, d);
    }
    return _cellcoordmax;
  }();

  if constexpr (TESTMODE) {
    for (int d = 0; d < get_ndim(GRID_TYPE); d++) {
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
            d, pos_component_vel_relative_to_flow ? '+' : '-', grid::coordlabel[d], cellindex, pktvelgridcoord[d],
            pktposgridcoord[d], grid::get_cellcoordmin(cellindex, d) / globals::tmin * tstart,
            cellcoordmax[d] / globals::tmin * tstart);
        printout("globals::tmin %g tstart %g tstart/globals::tmin %g\n", globals::tmin, tstart, tstart / globals::tmin);
        printout(" delta %g\n", delta);

        printout("packet dir [%g, %g, %g]\n", dir[0], dir[1], dir[2]);
        assert_always(false);
      }
    }
  }

  double distance = std::numeric_limits<double>::max();
  int snext{-1};

  if constexpr (GRID_TYPE == GridType::SPHERICAL1D) {
    // the only coordinate is the radius from the origin

    const double speed = vec_len(dir) * CLIGHT_PROP;  // just in case dir is not normalised

    const double r_outer = cellcoordmax[0] * tstart / globals::tmin;
    const double d_coordmaxboundary = expanding_shell_intersection(pos, dir, speed, r_outer, false, tstart);

    // upper d coordinate of the current cell
    if ((d_coordmaxboundary >= 0.) && (d_coordmaxboundary < distance)) {
      distance = d_coordmaxboundary;
      snext = (grid::get_cellcoordpointnum(cellindex, 0) == (grid::ncoordgrid[0] - 1))
                  ? -99
                  : cellindex + grid::get_coordcellindexincrement(0);
    }

    const double r_inner = grid::get_cellcoordmin(cellindex, 0) * tstart / globals::tmin;
    if (r_inner > 0.) {
      const double d_coordminboundary = expanding_shell_intersection(pos, dir, speed, r_inner, true, tstart);
      // lower d coordinate of the current cell
      if ((d_coordminboundary >= 0.) && (d_coordminboundary < distance)) {
        distance = d_coordminboundary;
        snext =
            (grid::get_cellcoordpointnum(cellindex, 0) == 0) ? -99 : cellindex - grid::get_coordcellindexincrement(0);
      }
    }
  } else if constexpr (GRID_TYPE == GridType::CYLINDRICAL2D) {
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
        expanding_shell_intersection(posnoz, dirnoz, xyspeed, r_outer, false, tstart);
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
          expanding_shell_intersection(posnoz, dirnoz, xyspeed, r_inner, true, tstart);
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

    if ((pktvelgridcoord[1] * tstart) > pktposgridcoord[1]) {
      const double t_zcoordmaxboundary = ((pktposgridcoord[1] - (pktvelgridcoord[1] * tstart)) /
                                          ((cellcoordmax[1]) - (pktvelgridcoord[1] * globals::tmin)) * globals::tmin) -
                                         tstart;
      const double d_coordmaxboundary_z = CLIGHT_PROP * t_zcoordmaxboundary;

      if ((d_coordmaxboundary_z >= 0.) && (d_coordmaxboundary_z < distance)) {
        distance = d_coordmaxboundary_z;
        snext = (grid::get_cellcoordpointnum(cellindex, 1) == (grid::ncoordgrid[1] - 1))
                    ? -99
                    : cellindex + grid::get_coordcellindexincrement(1);
      }
    } else {
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

  } else if constexpr (GRID_TYPE == GridType::CARTESIAN3D) {
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
      if ((pktvelgridcoord[d] * tstart) > pktposgridcoord[d]) {
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
      } else {
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
      for (int d2 = 0; d2 < get_ndim(GRID_TYPE); d2++) {
        printout("coord %d: cellcoordmin %g cellcoordmax %g\n", d2,
                 grid::get_cellcoordmin(cellindex, d2) * tstart / globals::tmin,
                 cellcoordmax[d2] * tstart / globals::tmin);
      }
      printout("tstart %g\n", tstart);
    }
  }

  assert_always((snext == -99) || ((snext >= 0) && (snext < grid::ngrid)));

  const double maxsdist = (GRID_TYPE == GridType::CARTESIAN3D)
                              ? globals::rmax * tstart / globals::tmin
                              : 2 * globals::rmax * (tstart + distance / CLIGHT_PROP) / globals::tmin;

  assert_always(distance >= 0. && distance <= maxsdist);

  if (distance > globals::max_path_step) {
    return {globals::max_path_step, cellindex};
  }

  return {distance, snext};
}

}  // namespace grid
