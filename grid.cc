// The model grid and propagation grid: reads the input ejecta model (model.txt), sets up cell
// densities, abundances, and temperatures, maps model cells onto the propagation grid (1D
// spherical, 2D cylindrical, 3D Cartesian, 2D spherical, or 3D spherical), and provides the cell
// geometry and boundary-crossing calculations used during packet propagation.

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
#include <deque>
#include <format>
#include <iostream>
#include <limits>
#include <numbers>
#include <optional>
#include <print>
#include <ranges>
#include <set>
#include <span>
#include <sstream>
#include <string>
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
#include "input.h"
#include "inputfilestream.h"
#include "kpkt.h"
#include "mpi_logging.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "outputfilestream.h"
#include "radfield.h"
#include "random.h"
#include "rpkt.h"
#include "sn3d.h"
#include "vectors.h"

namespace grid {

namespace {

std::array<int, 3> ncoordgrid{0, 0, 0};  // propagation grid dimensions

std::optional<GridType> model_type{};
// The format version of model.txt. A file of version 2 starts with the line "format_version 2" and has a header of
// keys with the cell edges. A file of version 1 has no version line
int model_format_version{1};
// The columns X_<element>_otherstable of a model.txt of format version 2. Such a column holds the mass fraction of
// the stable isotopes of an element that have no nuclide column. The first list gives the atomic number of each
// column. The node leaders hold the values of each model cell and column until read_elem_abundances() uses them.
std::vector<int> model_otherstable_atomic_numbers;
std::vector<float> model_otherstable_massfracs;
constexpr std::string_view otherstable_suffix = "_otherstable";

// The count of the cell lines that the optional key cell_line_count of a model.txt of format version 2 gives. A
// cell with no line is empty, thus only this count finds a file that lost lines at its end.
std::optional<ptrdiff_t> model_cell_line_count;

// The optional key mirror_plane z=0 of a model.txt of format version 2 gives mirror symmetry in the plane z = 0. The
// file then holds the cells with z >= 0, and each model cell also fills its mirror cell. The propagation grid has the
// cells of both halves, and the two cells of a mirror pair have the same model cell.
bool model_mirror_z{false};

ptrdiff_t npts_model = 0;  // number of model grid cells
ptrdiff_t nonempty_npts_model = 0;  // number of allocated non-empty model grid cells

double t_model = -1.;  // time at which densities in input model are correct.
std::vector<double> vout_model{};
std::array<int, 3> ncoord_model{};  // the model.txt input grid dimensions

// The cell edges of a 2D spherical or a 3D spherical model, which model.txt gives in its header. Each list has one
// value more than the count of the cells on that axis. A 2D spherical model has one phi cell from 0 to 2 pi.
std::vector<double> model_vel_r_edges{};  // radial velocity [cm/s], from 0 to vmax
std::vector<double> model_theta_edges{};  // polar angle from the +z axis [rad], from 0 to pi (pi / 2 with a mirror)
std::vector<double> model_phi_edges{};  // azimuthal angle phi = atan2(y, x) [rad], from 0 to 2 pi

double min_den{-1.};  // minimum model density

double mtot_input{0.};  // Total mass of the input model
double mfegroup{0.};  // Total mass of Fe group elements in ejecta

int first_input_cellid{-1};  // auto-determine first cell index in model.txt (usually 1 or 0)

constexpr int max_atomic_number = 150;  // abundances.txt can give elements up to this atomic number

// Initial co-ordinates of inner most corner of cell. An angular axis of a spherical grid holds the lower angle [rad]
// of each cell, which does not change with time.
std::array<std::vector<double>, 3> coord_pos_min_tmin{};

// The upper bound of the last cell on each axis. A length axis ends at globals::rmax. The polar angle of a spherical
// grid ends at pi, and the azimuthal angle ends at 2 pi.
std::array<double, 3> coord_pos_max_of_grid_tmin{};

// The cosine and the sine of each cell edge of the polar angle and of the azimuthal angle of a spherical propagation
// grid. boundary_distance() runs for every packet step, thus setup_grid_spherical() calculates them one time.
std::vector<double> theta_edge_cos{};
std::vector<double> theta_edge_sin{};
std::vector<double> phi_edge_cos{};
std::vector<double> phi_edge_sin{};

// the model cell of each propagation cell, or -1 for a propagation cell with no matter, e.g. a corner cell outside
// the outermost shell of a 1D or 2D model
std::vector<int> propcell_mgi;
std::vector<int> propcell_nonemptymgi;

// with FORCE_SPHERICAL_ESCAPE_SURFACE: 1 for a propagation cell whose inner corner lies outside the
// spherical escape surface. boundary_distance() runs for every packet step, so
// allocate_nonemptymodelcells() precomputes the per-cell radius comparison.
std::vector<uint8_t> propcell_outside_escape_surface;

std::vector<int> modelgrid_numpropcells;
std::vector<int> nonemptymgi_of_mgi;
std::vector<int> mgi_of_nonemptymgi;

std::vector<double> totmassnuclide{};  // total mass of each nuclide in the ejecta
// the grid-mapped mass of each nuclide in the propagation cells that FORCE_SPHERICAL_ESCAPE_SURFACE
// blanks. The nuclide-mass rescale in init_grid() adds this to the mapped total. The rescale then
// corrects only the discretisation error, and the blanked mass stays removed.
std::vector<double> totmassnuclide_blanked{};

MPI_shared_array<float> initnucmassfrac_allcells{};
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

struct ModelGridCellInput {
  float rhoinit = -1.;
  float ffegrp = 0.;
  // sum of the midpoint radius per propagation cell. A double, because the sum runs over up to millions of
  // propagation cells.
  double initial_radial_pos_sum = 0.;
  // sum of the volume averaged r^2 per propagation cell, for the kinetic energy of the BARNES
  // thermalisation scheme. A double, because the sum runs over up to millions of propagation cells
  // with radii of ~1e15 cm.
  double initial_radial_pos_squared_sum = 0.;
  float initelectronfrac = -1.;  // Ye: electrons (or protons) per nucleon. Negative until the input sets it
  // q: the trapped radiation energy per mass [erg/g], scaled from t_model to tmin (INITIAL_PACKETS_ON)
  float initenergyq = 0.;
};
MPI_shared_array<ModelGridCellInput> modelgrid_input{};

// Get number of dimensions
constexpr auto get_ndim(const GridType gridtype) -> int {
  switch (gridtype) {
    case GridType::SPHERICAL1D:
      return 1;
    case GridType::CYLINDRICAL2D:
    case GridType::SPHERICAL2D:
      return 2;
    case GridType::CARTESIAN3D:
    case GridType::SPHERICAL3D:
      return 3;
  }
  fatal_crash("Unknown grid type {}", static_cast<int>(gridtype));
}

// Get the number of axes that are lengths. The cell boundaries on such an axis move with the homologous expansion.
// These axes come first. Each later axis of a spherical grid is an angle, and its cell boundaries do not move.
constexpr auto get_num_length_axes(const GridType gridtype) -> int {
  switch (gridtype) {
    case GridType::SPHERICAL1D:
    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D:
      return 1;
    case GridType::CYLINDRICAL2D:
      return 2;
    case GridType::CARTESIAN3D:
      return 3;
  }
  fatal_crash("Unknown grid type {}", static_cast<int>(gridtype));
}

// maximum out-of-cell position error (in cm) attributable to floating-point rounding of the
// packet position updates. The tolerance scales with the coordinate magnitude because
// rounding errors are proportional to ulp(pos), e.g. about 8 cm at 4e16 cm.
[[nodiscard]] constexpr auto cellbound_tolerance(const double boundarypos) -> double {
  return std::max(10., std::abs(boundarypos) * 1e-12);
}

constexpr auto has_spherical_angle_axes(const GridType gridtype) -> bool {
  return gridtype == GridType::SPHERICAL2D || gridtype == GridType::SPHERICAL3D;
}

constexpr auto get_grid_type_name(const GridType gridtype) -> std::string_view {
  switch (gridtype) {
    case GridType::SPHERICAL1D:
    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D:
      return "spherical";
    case GridType::CYLINDRICAL2D:
      return "cylindrical";
    case GridType::CARTESIAN3D:
      return "uniform cuboidal";
  }
  return "unknown";
}

[[nodiscard]] DEVICE_FUNC constexpr auto get_coordlabel(const GridType gridtype, const int axis) -> std::string_view {
  assert_always(axis >= 0 && axis < get_ndim(gridtype));
  switch (gridtype) {
    case GridType::CARTESIAN3D:
      return std::array<std::string_view, 3>{"x", "y", "z"}[axis];
    case GridType::CYLINDRICAL2D:
      return std::array<std::string_view, 2>{"r", "z"}[axis];
    case GridType::SPHERICAL1D:
      return "r";
    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D:
      return std::array<std::string_view, 3>{"r", "theta", "phi"}[axis];
  }
  fatal_crash("Unknown grid type {}", static_cast<int>(gridtype));
}

void set_rho_tmin(const int modelgridindex, const float x) { modelgrid_input[modelgridindex].rhoinit = x; }

void set_initelectronfrac(const int modelgridindex, const float electronfrac) {
  if (std::isnan(electronfrac) || electronfrac < 0. || electronfrac > 1.001) {
    fatal_crash("input Ye {:g} for cell {} is outside the physical range [0, 1]", electronfrac, modelgridindex);
  }
  modelgrid_input[modelgridindex].initelectronfrac = electronfrac;
}

void read_possible_yefile() {
  if (model_format_version == 2 && inputfile_exists("Ye.txt")) {
    fatal_crash(
        "Ye.txt is present, but model.txt has format version 2. The column Ye of such a model.txt gives the electron "
        "fraction. Remove Ye.txt");
  }
  if (!inputfile_exists("Ye.txt")) {
    printlnlog("Ye.txt is not present, so the model keeps the electron fractions of model.txt");
    return;
  }

  // the electron fractions are written to node-shared memory, so only the node leaders read the file
  // (synchronised by the barrier below)
  if (globals::rank_in_node == 0) {
    auto filein = istream_required("Ye.txt");
    int nlines_in = 0;
    if (!(filein >> nlines_in) || nlines_in < 0) {
      fatal_crash("Ye.txt must start with the count of its cell lines, but the count is unreadable or negative");
    }

    const int last_mgi = get_npts_model() - 1;
    int cells_set = 0;
    bool has_first_or_last_cell = false;
    for (int n = 0; n < nlines_in; n++) {
      int cellnumberin = -1;
      float initelecfrac = 0.;
      if (!(filein >> cellnumberin >> initelecfrac)) {
        fatal_crash("Ye.txt ends or has an unreadable value after {} of its {} cell lines", n, nlines_in);
      }
      // Ye.txt uses the same cell ids as model.txt. read_ejecta_model() detects the id of the first
      // cell (0 or 1) and stores it in first_input_cellid before this function runs.
      const int mgi = cellnumberin - first_input_cellid;
      if (mgi < 0 || mgi > last_mgi) {
        // An out-of-range id is a sign of a cell id mismatch with model.txt. A silently shifted
        // electron fraction would give wrong grey opacities in every cell.
        fatal_crash("Ye.txt: cell {} is outside the model cell range [0..{}]. The cell ids of model.txt start at {}",
                    mgi, last_mgi, first_input_cellid);
      }
      has_first_or_last_cell = has_first_or_last_cell || mgi == 0 || mgi == last_mgi;
      set_initelectronfrac(mgi, initelecfrac);
      cells_set++;
    }
    if (cells_set > 0 && !has_first_or_last_cell) {
      // A file with either boundary id would hit the out-of-range stop above under the other id
      // convention. A sparse file with neither boundary id gives no such confirmation, e.g. a
      // 1-based file for a 0-based model would read with a one-cell shift.
      printlnlog(
          "[warning] Ye.txt holds neither the first cell 0 nor the last cell {}, so the file cannot confirm "
          "that its ids follow model.txt",
          last_mgi);
    }
    printlnlog("Ye.txt: set the initial electron fraction for {} of {} model cells (ids follow model.txt, first id {})",
               cells_set, get_npts_model(), first_input_cellid);
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

// convert a cell index into the integer coordinate index (x, y, z, or r) in [0, ncoordgrid[axis])
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordindex(const int cellindex, const int axis) -> int {
  return (cellindex / get_coordcellindexstride(axis)) % ncoordgrid[axis];
}

// get the minimum value of a coordinate at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the minimum x position in xyz coords, or the minimum radius
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordmin(const int cellindex, const int axis) -> double {
  return coord_pos_min_tmin[axis][get_cellcoordindex(cellindex, axis)];
}

// get the maximum position value of a coordinate axis at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the maximum x position in xyz coords, or the maximum radius in spherical 1D. An angular axis of a spherical
// grid gives the upper angle [rad] of the cell.
[[gnu::pure]] [[nodiscard]] DEVICE_FUNC auto get_cellcoordmax(const int cellindex, const int axis) -> double {
  const auto idx = get_cellcoordindex(cellindex, axis);
  const auto idxlast = ncoordgrid[axis] - 1;
  return idx < idxlast ? coord_pos_min_tmin[axis][idx + 1] : coord_pos_max_of_grid_tmin[axis];
}

// return the inner radius of a propagation cell at time tmin
auto get_cell_r_inner(const int cellindex, const GridType prop_gridtype) -> double {
  if (prop_gridtype == GridType::SPHERICAL1D || has_spherical_angle_axes(prop_gridtype)) {
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

  fatal_crash("Unknown propagation grid type {}", static_cast<int>(prop_gridtype));
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
      fatal_crash("Fe-group mass fraction {:g} is negative in cell {}", x, modelgridindex);
    }
    count_negative_input_massfrac(modelgridindex, -1, -1, x);
    x = 0.;
  }

  if (x > 1.001) {
    fatal_crash("model.txt gives the Fe-group mass fraction {:g} in cell {}, but a mass fraction must be 1 or less", x,
                modelgridindex);
  }
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
      fatal_crash("mass fraction {:g} for Z={} A={} is negative in cell {}", abund, decay::get_nuc_z(nucindex),
                  decay::get_nuc_a(nucindex), modelgridindex);
    }
    count_negative_input_massfrac(modelgridindex, decay::get_nuc_z(nucindex), decay::get_nuc_a(nucindex), abund);
    abund = 0.;
  }

  if (abund > 1.) {
    fatal_crash("the mass fraction {:.9g} of Z={} A={} in cell {} is above 1", abund, decay::get_nuc_z(nucindex),
                decay::get_nuc_a(nucindex), modelgridindex);
  }
  const ptrdiff_t num_nuclides = decay::get_num_nuclides();

  initnucmassfrac_allcells[(modelgridindex * num_nuclides) + nucindex] = abund;
}

void set_initenergyq(const int modelgridindex, const float initenergyq) {
  modelgrid_input[modelgridindex].initenergyq = initenergyq;
}

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

// radius of the midpoint of a propagation cell at time tmin. The midpoint has the middle value of each grid
// coordinate: r for a spherical grid, (rcyl, z) for a cylindrical grid, and (x, y, z) for a Cartesian grid.
auto get_cellradialposmid(const int cellindex) -> double {
  double radius_squared = 0.;
  for (int axis = 0; axis < get_num_length_axes(get_propgridtype()); axis++) {
    radius_squared += pow2((get_cellcoordmin(cellindex, axis) + get_cellcoordmax(cellindex, axis)) / 2);
  }
  return std::sqrt(radius_squared);
}

// volume averaged mean of the squared radius of a propagation cell. The kinetic energy sum needs the
// mean of r^2, which is larger than the square of the radius that get_cellradialposmid() gives.
auto get_cellradialposmeansquared(const int cellindex) -> double {
  const auto prop_gridtype = get_propgridtype();
  if (prop_gridtype == GridType::SPHERICAL1D || has_spherical_angle_axes(prop_gridtype)) {
    // each cell of a spherical grid is a part of a shell with a constant solid angle, thus the angles cancel
    const double r_inner = get_cellcoordmin(cellindex, 0);
    const double r_outer = get_cellcoordmax(cellindex, 0);
    return 3. / 5. * (pow5(r_outer) - pow5(r_inner)) / (pow3(r_outer) - pow3(r_inner));
  }

  if (prop_gridtype == GridType::CYLINDRICAL2D) {
    const double rcyl_inner = get_cellcoordmin(cellindex, 0);
    const double rcyl_outer = get_cellcoordmax(cellindex, 0);
    const double z_low = get_cellcoordmin(cellindex, 1);
    const double z_high = get_cellcoordmax(cellindex, 1);
    return ((pow2(rcyl_outer) + pow2(rcyl_inner)) / 2.) + ((pow2(z_low) + (z_low * z_high) + pow2(z_high)) / 3.);
  }

  // cubic cell: the mean of x^2 + y^2 + z^2 separates into a sum of the per-axis means
  double meansq = 0.;
  for (int axis = 0; axis < 3; axis++) {
    const double low = get_cellcoordmin(cellindex, axis);
    const double high = get_cellcoordmax(cellindex, axis);
    meansq += (pow2(low) + (low * high) + pow2(high)) / 3.;
  }
  return meansq;
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
  } else {
    globals::corrphotoionrenorm.reset();
    globals::gammaestimator.clear();
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

// Get the index of the cell that holds the value, from the lower edge of each cell of one axis. A cell holds its
// lower edge. A value below the first edge is in the first cell, and the last cell holds each value above its lower
// edge.
[[nodiscard]] DEVICE_FUNC constexpr auto get_cellindex_from_lower_edges(const std::span<const double> lower_edges,
                                                                        const double value) -> int {
  const auto num_cells = static_cast<int>(std::ssize(lower_edges));
  const auto index_above = static_cast<int>(std::ranges::upper_bound(lower_edges, value) - lower_edges.begin());
  return std::clamp(index_above - 1, 0, num_cells - 1);
}

// Get the polar angle theta from the +z axis in [0, pi] and the azimuthal angle phi = atan2(y, x) in [0, 2 pi) of a
// position or a direction. A vector on the z axis has phi = 0, and a vector of zero length has theta = 0.
[[nodiscard]] DEVICE_FUNC inline auto get_theta_phi_from_xyz(const Vec3d& vec_xyz) -> std::array<double, 2> {
  const double theta = std::atan2(std::sqrt(pow2(vec_xyz[0]) + pow2(vec_xyz[1])), vec_xyz[2]);
  double phi = std::atan2(vec_xyz[1], vec_xyz[0]);
  if (phi < 0.) {
    phi += 2 * PI;
  }
  // the sum of 2 pi and a small negative angle can round to 2 pi, which is the angle 0
  if (phi >= 2 * PI) {
    phi = 0.;
  }
  return {theta, phi};
}

// Get the index of the cell of a 2D spherical or a 3D spherical model that holds a position at tmin, or -1 for a
// position outside the outermost shell. The index is i_r + n_r * (i_theta + n_theta * i_phi).
[[nodiscard]] auto get_spherical_model_cellindex_tmin(const Vec3d& pos_tmin) -> int {
  const double vel_r = vec_len(pos_tmin) / globals::tmin;
  if (vel_r > model_vel_r_edges.back()) {
    return -1;
  }
  const auto [theta_pos, phi] = get_theta_phi_from_xyz(pos_tmin);
  // the polar axis of a mirrored model ends at pi / 2, and a position with z < 0 takes its mirror cell
  const double theta = (model_mirror_z && theta_pos > PI / 2) ? PI - theta_pos : theta_pos;
  // the last value of each list is the upper edge of the last cell
  const auto lower_edges = [](const std::vector<double>& edges) {
    return std::span<const double>{edges}.first(edges.size() - 1);
  };
  const int index_r = get_cellindex_from_lower_edges(lower_edges(model_vel_r_edges), vel_r);
  const int index_theta = get_cellindex_from_lower_edges(lower_edges(model_theta_edges), theta);
  const int index_phi = get_cellindex_from_lower_edges(lower_edges(model_phi_edges), phi);
  return index_r + (ncoord_model[0] * (index_theta + (ncoord_model[1] * index_phi)));
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

// The mirror axis of a model with the key mirror_plane is the axis 1: z of a cylindrical model and theta of a
// spherical model.
constexpr int mirror_axis = 1;

// Get the index on the mirror axis of the model grid for an index on that axis of the propagation grid. The model
// grid of a mirrored model holds the half with z >= 0. A cylindrical model counts its z cells from the plane z = 0,
// and a spherical model counts its polar cells from the +z axis.
[[nodiscard]] auto get_model_mirror_axis_index(const int propindex) -> int {
  const int num_half = ncoord_model[mirror_axis];
  if (get_modelgridtype() == GridType::CYLINDRICAL2D) {
    return (propindex >= num_half) ? propindex - num_half : num_half - 1 - propindex;
  }
  return (propindex < num_half) ? propindex : (2 * num_half) - 1 - propindex;
}

// Get the propagation cell with z >= 0 of a model cell of a mirrored model on its own grid type
[[nodiscard]] auto get_upper_propcell_of_mirrored_mgi(const int mgi) -> int {
  const int index0 = mgi % ncoord_model[0];
  const int index1 = (mgi / ncoord_model[0]) % ncoord_model[1];
  const int index2 = mgi / (ncoord_model[0] * ncoord_model[1]);
  const int propindex1 = (get_modelgridtype() == GridType::CYLINDRICAL2D) ? ncoord_model[1] + index1 : index1;
  return index0 + (ncoordgrid[0] * (propindex1 + (ncoordgrid[1] * index2)));
}

// Map a mirrored model onto a propagation grid of its own type. The propagation grid has the cells of both halves,
// and each propagation cell takes the model cell of the mirror pair.
void map_mirroredmodeltogrid() {
  for (int cellindex = 0; cellindex < ngrid; cellindex++) {
    const int index0 = get_cellcoordindex(cellindex, 0);
    const int index1 = get_model_mirror_axis_index(get_cellcoordindex(cellindex, mirror_axis));
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

// Read the element mass fractions of each model cell with a density above zero. A cell with no mass in the
// elements of compositiondata.txt gets zero density, so that the mapping makes it an empty cell. Return the mass
// fractions in model cell order (an empty vector on the ranks other than the node leader) and the number of
// cells that became empty.
auto read_elem_abundances() -> std::tuple<std::vector<float>, int> {
  // a model.txt of format version 2 gives the element mass fractions, and such a model has no abundances.txt
  const bool elements_from_model = (model_format_version == 2);
  const std::string_view massfrac_source = elements_from_model ? "model.txt" : "abundances.txt";
  if (elements_from_model && inputfile_exists("abundances.txt")) {
    fatal_crash(
        "abundances.txt is present, but model.txt has format version 2. Such a model.txt gives the mass fraction of "
        "an element as the sum of its nuclide columns and its column X_<element>_otherstable. Remove abundances.txt");
  }
  if (elements_from_model) {
    printlnlog("calculating the element mass fractions from the columns of model.txt...");
  } else {
    printlnlog("reading abundances.txt...");
  }
  const int nelements = get_nelements();

  // abundances.txt holds one line for each model cell, also for an empty cell. Each line gives the cell number and
  // then the mass fraction of each element from Z=1 upwards. The line can end after the last element that the
  // model needs.

  // the density is in node-shared memory, so only the node leader of each node parses the file and
  // writes the values (synchronised by the barrier below). The other ranks would just discard everything they read
  std::vector<float> elem_massfracs_of_mgi;
  int ncells_emptied = 0;
  if (globals::rank_in_node == 0) {
    elem_massfracs_of_mgi.assign(static_cast<size_t>(get_npts_model()) * nelements, 0.F);
    std::set<int> atomic_numbers_without_atomic_data;  // elements with mass in the cells that became empty
    // Every log line gets a timestamp and a flush, so a large 3D model must not warn per cell.
    constexpr int max_warnings = 10;

    // Check the element mass fractions of one cell with density, and store those of the elements with atomic data.
    // elem_massfracs_in holds the mass fraction of each element from Z=1 upwards.
    const auto store_cell_elem_massfracs = [&](const int mgi,
                                               const std::array<float, max_atomic_number>& elem_massfracs_in,
                                               const double elem_massfrac_sum) {
      // a sum above one is not physical. The limit permits the rounding of the values in the file.
      if (elem_massfrac_sum > 1. + 1e-3) {
        fatal_crash("{}: cell {} has element mass fractions that sum to {:g}, which is more than 1", massfrac_source,
                    mgi, elem_massfrac_sum);
      }

      bool has_mass_in_composition_elements = false;
      const auto tracked_isotope_massfrac_sum_by_atomic_number = get_tracked_isotope_massfrac_sum_by_atomic_number(mgi);

      for (int element = 0; element < nelements; element++) {
        // set the mass fraction of each element that the atomic data includes
        const int atomic_number = get_atomicnumber(element);
        const float elemmassfrac = elem_massfracs_in[atomic_number - 1];
        assert_always(elemmassfrac >= 0.);

        // read_ejecta_model() has already set the nuclide mass fractions. Check here, while
        // the mass fractions are still exactly as the input files gave them, that abundances.txt gives the
        // element at least as much mass as model.txt gives its tracked isotopes. Any later mapping rescale
        // changes the nuclide mass fractions but not the elemental ones, so this cannot be tested afterwards.
        // a small negative remainder is allowed for roundoff error
        const double tracked_isotope_massfrac_sum = tracked_isotope_massfrac_sum_by_atomic_number[atomic_number];
        if ((elemmassfrac - tracked_isotope_massfrac_sum) < -1e-2) {
          fatal_crash(
              "{} cell {} gives the mass fraction {:g} for Z={}, but the nuclide columns of model.txt give this "
              "element the larger mass fraction {:g}",
              massfrac_source, mgi, elemmassfrac, atomic_number, tracked_isotope_massfrac_sum);
        }

        elem_massfracs_of_mgi[(static_cast<size_t>(mgi) * nelements) + element] = elemmassfrac;
        has_mass_in_composition_elements =
            has_mass_in_composition_elements ||
            (get_nions(element) > 0 && (elemmassfrac > 0.F || tracked_isotope_massfrac_sum > 0.));
      }

      if (!has_mass_in_composition_elements) {
        // e.g. a mapped model gives zero abundances to a cell whose mass comes only from particles without
        // nucleosynthesis data. Matter without the elements of compositiondata.txt has no opacity and no cooling
        // processes.
        ncells_emptied++;
        // every element with mass in this cell has no atomic data, so list them for the user
        std::string elements_without_atomic_data;
        for (int elem_z_index = 0; elem_z_index < std::ssize(elem_massfracs_in); elem_z_index++) {
          if (elem_massfracs_in[elem_z_index] > 0.F) {
            atomic_numbers_without_atomic_data.insert(elem_z_index + 1);
            elements_without_atomic_data +=
                std::format(" Z={} ({:.2e})", elem_z_index + 1, elem_massfracs_in[elem_z_index]);
          }
        }
        if (ncells_emptied <= max_warnings) {
          printlnlog(
              "[warning] read_elem_abundances: cell {} has density {:g} [g/cm3] at tmin and no mass in the "
              "elements of compositiondata.txt. The cell becomes empty. Elements with mass and no atomic data:{}",
              mgi, get_rho_tmin(mgi), elements_without_atomic_data.empty() ? " none" : elements_without_atomic_data);
        }
        set_rho_tmin(mgi, 0.);
      }
    };

    if (elements_from_model) {
      // The mass fraction of an element is the sum of its nuclide columns of model.txt and of its column
      // X_<element>_otherstable, which holds the stable isotopes that have no column of their own.
      const auto num_otherstable_columns = std::ssize(model_otherstable_atomic_numbers);
      for (int mgi = 0; mgi < get_npts_model(); mgi++) {
        if (get_rho_tmin(mgi) <= 0.) {
          continue;
        }
        const auto tracked_isotope_massfrac_sum_by_atomic_number =
            get_tracked_isotope_massfrac_sum_by_atomic_number(mgi);
        std::array<double, max_atomic_number> elem_massfracs_double{};
        for (int atomic_number = 1; atomic_number <= max_atomic_number; atomic_number++) {
          elem_massfracs_double[atomic_number - 1] = tracked_isotope_massfrac_sum_by_atomic_number[atomic_number];
        }
        for (auto column = 0Z; column < num_otherstable_columns; column++) {
          elem_massfracs_double[model_otherstable_atomic_numbers[column] - 1] +=
              model_otherstable_massfracs[(static_cast<ptrdiff_t>(mgi) * num_otherstable_columns) + column];
        }
        double elem_massfrac_sum = 0.;
        std::array<float, max_atomic_number> elem_massfracs_in{};
        for (int elem_z_index = 0; elem_z_index < max_atomic_number; elem_z_index++) {
          elem_massfracs_in[elem_z_index] = static_cast<float>(elem_massfracs_double[elem_z_index]);
          elem_massfrac_sum += elem_massfracs_in[elem_z_index];
        }
        store_cell_elem_massfracs(mgi, elem_massfracs_in, elem_massfrac_sum);
      }
      model_otherstable_massfracs.clear();
      model_otherstable_massfracs.shrink_to_fit();
    } else {
      auto abundance_file = istream_required("abundances.txt");
      std::string line;
      for (int mgi = 0; mgi < get_npts_model(); mgi++) {
        if (!get_noncommentline(abundance_file, line)) {
          fatal_crash("abundances.txt ends before cell {} of the {} cells of model.txt", mgi, get_npts_model());
        }
        auto remainder = std::string_view{line};

        int cellnumberinput = -1;
        if (!parse_next_token(remainder, cellnumberinput) || cellnumberinput != mgi + first_input_cellid) {
          fatal_crash("abundances.txt cell {} must start with the cell id {} as in model.txt, but the line is '{}'",
                      mgi, mgi + first_input_cellid, line);
        }

        // the abundances.txt file specifies the elemental mass fractions for each model cell.
        // The values are true mass fractions and get no normalisation. Their sum can be below one.
        // The abundances begin with hydrogen, helium, etc, going as far up the atomic numbers as required
        double elem_massfrac_sum = 0.;
        std::array<float, max_atomic_number> elem_massfracs_in{};
        double abund_in = 0.;
        for (int elem_z_index = 0; elem_z_index < std::ssize(elem_massfracs_in); elem_z_index++) {
          if (!parse_next_token(remainder, abund_in)) {
            // a row with no values is a cell with no elements
            break;
          }

          if (abund_in < std::numeric_limits<float>::min()) {
            if (!(abund_in > -1e-6)) {
              fatal_crash("abundances.txt cell {} gives the negative mass fraction {:g} for Z={}", mgi, abund_in,
                          elem_z_index + 1);
            }
            abund_in = 0.;
          }
          elem_massfracs_in[elem_z_index] = static_cast<float>(abund_in);
          elem_massfrac_sum += elem_massfracs_in[elem_z_index];
        }

        // parse_next_token() gives false at the end of the row and also for a token that is not a
        // number, e.g. "nan". The loop above stops on both, so a bad token would silently zero every
        // later element. Only whitespace may remain here.
        remainder.remove_prefix(std::min(remainder.find_first_not_of(token_whitespace), remainder.size()));
        if (!remainder.empty()) {
          fatal_crash("read_elem_abundances: cell {} has an unreadable token at '{}'", mgi, remainder);
        }

        if (get_rho_tmin(mgi) > 0.) {
          store_cell_elem_massfracs(mgi, elem_massfracs_in, elem_massfrac_sum);
        }
      }
    }

    if (!atomic_numbers_without_atomic_data.empty()) {
      std::string atomic_numbers_text;
      for (const int atomic_number : atomic_numbers_without_atomic_data) {
        atomic_numbers_text += std::format(" {}", atomic_number);
      }
      printlnlog(
          "[warning] read_elem_abundances: add atomic data for these elements to keep the mass of the cells that "
          "became empty. Atomic numbers:{}",
          atomic_numbers_text);
    }
  }

  // barrier to make sure node master has set values in node shared memory
  MPI_Barrier_allranks();
  MPI_Bcast_safe(ncells_emptied, 0, globals::mpi_comm_node);
  printlnlog("finished setting the element mass fractions from {}", massfrac_source);
  return {std::move(elem_massfracs_of_mgi), ncells_emptied};
}

void parse_model_headerline(const std::string& line, std::vector<int>& zlist, std::vector<int>& alist,
                            std::vector<std::string>& colnames) {
  // custom header line
  std::istringstream iss(line);
  std::string token;

  int columnindex = -1;

  while (iss >> token) {
    columnindex++;

    // The header of a file of format version 2 gives the cell edges, thus such a file has no position column, and
    // rho comes directly after the cell id.
    const bool model_has_edge_lists = (model_format_version == 2);

    // -1 for a column that this model type cannot have
    int required_columnindex = columnindex;
    if (token == "#inputcellid") {
      required_columnindex = 0;
    } else if (token == "velocity_outer" || token == "vel_r_max_kmps") {
      required_columnindex = model_has_edge_lists ? -1 : 1;
    } else if (token.starts_with("pos_") && !model_has_edge_lists) {
      continue;
    } else if (token.starts_with("pos_") || token.starts_with("vel_")) {
      required_columnindex = -1;
    } else if (token == "logrho") {
      // 1D models have log10(rho [g/cm3])
      required_columnindex = (get_modelgridtype() == GridType::SPHERICAL1D && !model_has_edge_lists) ? 2 : -1;
    } else if (token == "rho") {
      // 2D and 3D models have rho [g/cm3] after the cell id and the position columns
      if (model_has_edge_lists) {
        required_columnindex = 1;
      } else if (get_modelgridtype() == GridType::SPHERICAL1D) {
        required_columnindex = -1;
      } else {
        required_columnindex = get_ndim(get_modelgridtype()) + 1;
      }
    } else if (token.starts_with("X_") && token.ends_with(otherstable_suffix)) {
      if (!model_has_edge_lists) {
        fatal_crash("the model.txt header line has the column {}, but only format version 2 has such a column", token);
      }
      const int atomic_number =
          decay::get_atomic_number_of_elsymbol(token.substr(2, token.size() - 2 - otherstable_suffix.size()));
      if (atomic_number < 1 || atomic_number > max_atomic_number) {
        fatal_crash("the model.txt header line has the column {}, which does not name an element", token);
      }
      if (std::ranges::find(model_otherstable_atomic_numbers, atomic_number) !=
          model_otherstable_atomic_numbers.end()) {
        fatal_crash("the model.txt header line has two columns {} for the element with Z={}", token, atomic_number);
      }
      model_otherstable_atomic_numbers.push_back(atomic_number);
      colnames.push_back(token);
      zlist.push_back(-1);
      alist.push_back(-1);
    } else if (token.starts_with("X_") && token != "X_Fegroup") {
      const int z = decay::get_nucstring_z(token.substr(2));  // + 2 skips the 'X_'
      const int a = decay::get_nucstring_a(token.substr(2));
      // two names can give one nuclide, e.g. X_Ni56 and X_ni56, and a second column then replaces the mass fraction
      for (auto i = 0Z; i < std::ssize(colnames); i++) {
        if (zlist[i] == z && alist[i] == a) {
          fatal_crash("the model.txt header line has the columns {} and {}, which are the same nuclide (Z={} A={})",
                      colnames[i], token, z, a);
        }
      }
      colnames.push_back(token);
      zlist.push_back(z);
      alist.push_back(a);
    } else {
      if (std::ranges::find(colnames, token) != colnames.end()) {
        fatal_crash("the model.txt header line has the column {} two times", token);
      }
      colnames.push_back(token);
      zlist.push_back(-1);
      alist.push_back(-1);
    }

    if (required_columnindex < 0) {
      fatal_crash(
          "the model.txt header line has the column {}, but a {}D {} model of format version {} cannot have this "
          "column",
          token, get_ndim(get_modelgridtype()), get_grid_type_name(get_modelgridtype()), model_format_version);
    }
    if (columnindex != required_columnindex) {
      fatal_crash("the model.txt header line has the column {} at index {}, but the required index is {}", token,
                  columnindex, required_columnindex);
    }
  }
}

auto get_token_count(std::string const& line) -> int {
  std::string token;
  int abundcolcount = 0;
  auto ssline = std::istringstream{line};
  while (ssline >> token) {
    abundcolcount++;
  }
  return abundcolcount;
}

// The header probes of model.txt read ahead into the data lines. The cell reader takes those lines
// first, so no reader seeks in the file, which a compressed file does not support.
struct ModelFileReader {
  InputFileStream file;
  std::deque<std::string> lines_read_ahead;

  [[nodiscard]] auto getline(std::string& line) -> bool {
    if (!lines_read_ahead.empty()) {
      line = std::move(lines_read_ahead.front());
      lines_read_ahead.pop_front();
      return true;
    }
    return static_cast<bool>(std::getline(file, line));
  }
};

void read_model_radioabundances(ModelFileReader& fmodel, std::string_view& remainder, const int mgi,
                                const bool keepcell, const std::vector<std::string>& colnames,
                                const std::vector<int>& nucindexlist, const bool one_line_per_cell) {
  if (!one_line_per_cell) {
    static std::string line;
    if (!fmodel.getline(line)) {
      fatal_crash("model.txt ends before the second line of cell {}", mgi);
    }
    remainder = std::string_view{line};
  }

  // A line of an empty cell of format version 2 must also have one value for each column. Thus a line that lost its
  // end stops the run, and it does not count as a cell line for the key cell_line_count.
  if (!keepcell && model_format_version != 2) {
    return;
  }

  const auto num_otherstable_columns = std::ssize(model_otherstable_atomic_numbers);
  ptrdiff_t otherstable_column = 0;
  for (auto i = 0Z; i < std::ssize(colnames); i++) {
    double valuein = 0.;  // a mass fraction or another column value, e.g. Ye or q
    if (!parse_next_token(remainder, valuein)) {
      fatal_crash("model.txt cell {} has no readable value for the column {}", mgi, colnames[i]);
    }
    if (!keepcell) {
      continue;
    }

    if (nucindexlist[i] >= 0) {
      set_modelinitnucmassfrac(mgi, nucindexlist[i], static_cast<float>(valuein));
    } else if (colnames[i].starts_with("X_") && colnames[i].ends_with(otherstable_suffix)) {
      // the columns have the sequence of model_otherstable_atomic_numbers
      if (valuein < 0.) {
        if (valuein <= -1e-6) {
          fatal_crash("model.txt cell {} has the negative value {:g} in the column {}", mgi, valuein, colnames[i]);
        }
        valuein = 0.;
      }
      model_otherstable_massfracs[(static_cast<ptrdiff_t>(mgi) * num_otherstable_columns) + otherstable_column] =
          static_cast<float>(valuein);
      otherstable_column++;
    } else if (colnames[i] == "X_Fegroup") {
      set_ffegrp(mgi, static_cast<float>(valuein));
    } else if (colnames[i] == "cellYe" || colnames[i] == "Ye") {
      set_initelectronfrac(mgi, static_cast<float>(valuein));
    } else if (colnames[i] == "q") {
      // The q column holds the trapped radiation energy per mass at t_model. It already includes the adiabatic
      // losses before t_model. The energy of a comoving mass element falls as 1/t, so the value is scaled to tmin.
      if (!(valuein >= 0.)) {
        fatal_crash("model.txt cell {} has the value {:g} in the column q, but the value must be zero or positive", mgi,
                    valuein);
      }
      set_initenergyq(mgi, static_cast<float>(valuein * t_model / globals::tmin));
    } else if (colnames[i] == "tracercount") {
      ;
    } else {
      // reported once after the model read (checking only mgi == 0 could miss it entirely if cell 0 is empty)
      ignored_model_columns.insert(colnames[i]);
    }
  }
  double valuein = 0.;
  if (parse_next_token(remainder, valuein)) {
    fatal_crash("model.txt cell {} has more values than the header has columns", mgi);
  }
  // a cell line of format version 2 has only whitespace after its values. Text there shows a damaged line
  remainder.remove_prefix(std::min(remainder.find_first_not_of(token_whitespace), remainder.size()));
  if (model_format_version == 2 && !remainder.empty()) {
    fatal_crash("model.txt cell {} has the text '{}' after the values of the columns", mgi, remainder);
  }
}

auto read_model_columns(ModelFileReader& fmodel) -> std::tuple<std::vector<std::string>, std::vector<int>, bool> {
  std::vector<int> zlist;
  std::vector<int> alist;
  std::vector<std::string> colnames;

  std::string line;
  if (!fmodel.getline(line)) {
    fatal_crash("model.txt ends before the column header or the first cell");
  }

  std::string headerline;

  const bool header_specified = lineiscommentonly(line);

  if (header_specified) {
    // line is the header
    headerline = line;
    if (!fmodel.getline(line)) {
      fatal_crash("model.txt ends after the column header, before the first cell");
    }
    if (model_format_version == 2) {
      // a comment line or an empty line can come before each cell line of format version 2, also before the first one
      while (lineiscommentonly(line)) {
        if (!fmodel.getline(line)) {
          fatal_crash("model.txt ends after the key columns, before the first cell");
        }
      }
    }
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
      case GridType::SPHERICAL2D:
      case GridType::SPHERICAL3D:
        fatal_crash("model.txt: a 2D spherical or a 3D spherical model needs format version 2");
    }
    headerline += " X_Fegroup X_Ni56 X_Co56 X_Fe52 X_Cr48";
  }

  int colcount = get_token_count(line);
  const bool one_line_per_cell = (colcount >= get_token_count(headerline));

  printlnlog("model.txt has {} line per cell format", one_line_per_cell ? "one" : "two");
  if (!one_line_per_cell && model_format_version == 2) {
    fatal_crash(
        "model.txt: the first cell has {} values, but the key columns has {} names. A cell has one line in format "
        "version 2",
        colcount, get_token_count(headerline));
  }

  std::string secondline;
  if (!one_line_per_cell) {  // add columns from the second line
    if (!fmodel.getline(secondline)) {
      fatal_crash("model.txt ends before the second line of the first cell");
    }
    colcount += get_token_count(secondline);
  }

  // the cell reader takes the lines of the first cell again
  fmodel.lines_read_ahead.push_back(line);
  if (!one_line_per_cell) {
    fmodel.lines_read_ahead.push_back(secondline);
  }

  if (!header_specified && colcount > get_token_count(headerline)) {
    headerline += " X_Ni57 X_Co57";
  }

  if (colcount != get_token_count(headerline)) {
    fatal_crash("model.txt cell 0 has {} values, but the {} '{}' has {} columns", colcount,
                header_specified ? "header line" : "default header for a file without a header line", headerline,
                get_token_count(headerline));
  }

  if (header_specified) {
    printlnlog("model.txt has a header line.");
  } else {
    printlnlog("model.txt has no header line. Using default: {}", headerline);
  }

  parse_model_headerline(headerline, zlist, alist, colnames);

  if (std::ranges::none_of(colnames, [](const std::string& colname) { return colname == "X_Fegroup"; })) {
    printlnlog(
        "model.txt has no column X_Fegroup, thus the Fe-group mass fraction of each cell is 0. This mass fraction "
        "sets the photoelectric opacity and the pair production opacity of the gamma rays, and the grey opacity of "
        "FEGROUP_APPROX");
  }

  // a column that is not a nuclide has Z = -1, and the decay network does not hold the neutron (Z = 0)
  std::vector<int> nuclide_columns_z;
  std::vector<int> nuclide_columns_a;
  for (const auto [z, a] : std::views::zip(zlist, alist)) {
    if (z > 0) {
      nuclide_columns_z.push_back(z);
      nuclide_columns_a.push_back(a);
    }
  }
  decay::init_nuclides(nuclide_columns_z, nuclide_columns_a);

  std::vector<int> nucindexlist(zlist.size());
  for (auto i = 0Z; i < std::ssize(zlist); i++) {
    nucindexlist[i] = (zlist[i] > 0) ? decay::get_nucindex(zlist[i], alist[i]) : -1;
  }

  assert_always(npts_model > 0);

  const ptrdiff_t num_nuclides = decay::get_num_nuclides();

  initnucmassfrac_allcells = MPI_shared_array<float>(npts_model * num_nuclides, 0.);
  if (globals::rank_in_node == 0) {
    model_otherstable_massfracs.assign(npts_model * std::ssize(model_otherstable_atomic_numbers), 0.F);
  }
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
      const int index_theta = (mgi / ncoord_model[0]) % ncoord_model[1];
      const int index_phi = mgi / (ncoord_model[0] * ncoord_model[1]);
      const double radialpart =
          (pow3(model_vel_r_edges[index_r + 1]) - pow3(model_vel_r_edges[index_r])) * pow3(globals::tmin) / 3.;
      const double polarpart = std::cos(model_theta_edges[index_theta]) - std::cos(model_theta_edges[index_theta + 1]);
      // each cell of a mirrored model also fills its mirror cell
      const double mirrorfactor = model_mirror_z ? 2. : 1.;
      return mirrorfactor * radialpart * polarpart * (model_phi_edges[index_phi + 1] - model_phi_edges[index_phi]);
    }
  }

  fatal_crash("Unknown model grid type {}", static_cast<int>(get_modelgridtype()));
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

// Stop the run when the grid corners expand faster than light and nothing removes them. The
// propagation cannot treat a superluminal boundary, so the setup must either cut the corners with
// FORCE_SPHERICAL_ESCAPE_SURFACE or use a smaller vmax.
void require_subluminal_corners(const double vmax_corner) {
  printlnlog("corner vmax {:g} [cm/s] ({:.2f}c)", vmax_corner, vmax_corner / CLIGHT);
  if (!FORCE_SPHERICAL_ESCAPE_SURFACE && vmax_corner >= CLIGHT) {
    fatal_crash(
        "the grid corners expand faster than light. Enable FORCE_SPHERICAL_ESCAPE_SURFACE in "
        "artisoptions.h to remove the superluminal corners, or reduce vmax.");
  }
}

void setup_grid_cartesian_3d() {
  // vmax is per coordinate, but the simulation volume corners will have a higher expansion velocity than the sides
  require_subluminal_corners(sqrt(3 * pow2(globals::vmax)));

  // Set grid size for uniform xyz grid
  // a 3D ejecta model gives the propagation grid its own cell counts
  if (get_modelgridtype() != GridType::CARTESIAN3D) {
    ncoordgrid = {CUBOID_NCOORDGRID_X, CUBOID_NCOORDGRID_Y, CUBOID_NCOORDGRID_Z};
  }

  // artis assumes in some places that the cells are cubes, not cuboids
  if (ncoordgrid[0] <= 0 || ncoordgrid[0] != ncoordgrid[1] || ncoordgrid[0] != ncoordgrid[2]) {
    fatal_crash(
        "the 3D propagation grid has {} x {} x {} cells, but the three counts must be equal and positive. The counts "
        "come from {}",
        ncoordgrid[0], ncoordgrid[1], ncoordgrid[2],
        get_modelgridtype() == GridType::CARTESIAN3D ? "model.txt"
                                                     : "CUBOID_NCOORDGRID_X, CUBOID_NCOORDGRID_Y, and "
                                                       "CUBOID_NCOORDGRID_Z in artisoptions.h");
  }

  ngrid = static_cast<ptrdiff_t>(ncoordgrid[0]) * ncoordgrid[1] * ncoordgrid[2];
  reserve_resize(coord_pos_min_tmin[0], ncoordgrid[0]);
  reserve_resize(coord_pos_min_tmin[1], ncoordgrid[1]);
  reserve_resize(coord_pos_min_tmin[2], ncoordgrid[2]);

  for (int axis = 0; axis < 3; axis++) {
    for (int i = 0; i < ncoordgrid[axis]; i++) {
      coord_pos_min_tmin[axis][i] = -globals::rmax + (2 * i * globals::rmax / ncoordgrid[axis]);
    }
  }
  coord_pos_max_of_grid_tmin = {globals::rmax, globals::rmax, globals::rmax};
}

// The readers of model.txt check the radial edges as velocities. Two edges that differ in the last bits can give
// radii at tmin that are equal or nearly equal. A packet at a boundary is in the position tolerance of that boundary
// (see cellbound_tolerance). In a shell that is not wider than this tolerance, the packet is at once at the next
// boundary, thus the shell gets no path length and no estimators. A grid of radial shells gets this check after the
// multiplication by tmin.
void require_increasing_shell_radii_at_tmin() {
  for (int n_r = 0; n_r < ncoordgrid[0]; n_r++) {
    const double radius_outer_tmin = (n_r < ncoordgrid[0] - 1) ? coord_pos_min_tmin[0][n_r + 1] : globals::rmax;
    const double tolerance = cellbound_tolerance(radius_outer_tmin);
    if (!(radius_outer_tmin - coord_pos_min_tmin[0][n_r] > tolerance)) {
      fatal_crash(
          "model.txt: the radial shell {} of the grid has the inner radius {:.17g} [cm] and the outer radius {:.17g} "
          "[cm] at tmin. Each shell must be wider than the position tolerance of a cell boundary, which is {:g} [cm] "
          "at this radius",
          n_r, coord_pos_min_tmin[0][n_r], radius_outer_tmin, tolerance);
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
  coord_pos_max_of_grid_tmin = {globals::rmax, globals::rmax, globals::rmax};
  require_increasing_shell_radii_at_tmin();
}

void setup_grid_cylindrical_2d() {
  require_subluminal_corners(sqrt(2 * pow2(globals::vmax)));

  assert_always(get_modelgridtype() == GridType::CYLINDRICAL2D);

  ncoordgrid = ncoord_model;
  if (model_mirror_z) {
    // the propagation grid has the z cells of both halves
    ncoordgrid[1] = 2 * ncoord_model[1];
  }

  ngrid = static_cast<ptrdiff_t>(ncoordgrid[0]) * ncoordgrid[1];
  assert_always(ngrid == get_npts_model() * (model_mirror_z ? 2 : 1));

  reserve_resize(coord_pos_min_tmin[0], ncoordgrid[0]);
  for (int n_rcyl = 0; n_rcyl < ncoordgrid[0]; n_rcyl++) {
    coord_pos_min_tmin[0][n_rcyl] = n_rcyl * globals::rmax / ncoord_model[0];
  }

  reserve_resize(coord_pos_min_tmin[1], ncoordgrid[1]);
  for (int n_z = 0; n_z < ncoordgrid[1]; n_z++) {
    coord_pos_min_tmin[1][n_z] = globals::rmax * (-1 + (n_z * 2. / ncoordgrid[1]));
  }
  coord_pos_max_of_grid_tmin = {globals::rmax, globals::rmax, globals::rmax};
}

// Set up a 2D spherical or a 3D spherical propagation grid with the cell edges of the model. A 2D spherical grid is
// a 3D spherical grid with one phi cell from 0 to 2 pi, thus one geometry code serves both.
void setup_grid_spherical() {
  assert_always(has_spherical_angle_axes(get_modelgridtype()));
  assert_always(get_propgridtype() == get_modelgridtype());

  // the propagation grid of a mirrored model has the polar cells of both halves
  std::vector<double> prop_theta_edges = model_theta_edges;
  if (model_mirror_z) {
    for (auto i = std::ssize(model_theta_edges) - 2; i >= 0; i--) {
      prop_theta_edges.push_back(PI - model_theta_edges[i]);
    }
  }
  const std::array<std::span<const double>, 3> model_edges{model_vel_r_edges, prop_theta_edges, model_phi_edges};
  for (int axis = 0; axis < 3; axis++) {
    ncoordgrid[axis] = static_cast<int>(std::ssize(model_edges[axis])) - 1;
    // the last edge of each axis is the upper bound of the grid
    coord_pos_min_tmin[axis].assign(model_edges[axis].begin(), model_edges[axis].end() - 1);
  }
  for (auto& radius_inner_tmin : coord_pos_min_tmin[0]) {
    radius_inner_tmin *= globals::tmin;
  }
  coord_pos_max_of_grid_tmin = {globals::rmax, PI, 2 * PI};
  require_increasing_shell_radii_at_tmin();

  ngrid = static_cast<ptrdiff_t>(ncoordgrid[0]) * ncoordgrid[1] * ncoordgrid[2];
  assert_always(ngrid == get_npts_model() * (model_mirror_z ? 2 : 1));

  // The centre of the grid is a point of every polar cell and every azimuthal cell of the first shell. A packet
  // path within cellbound_tolerance(0.) of the centre changes its cell there in one event, thus that distance must
  // be small against the first shell.
  const double first_shell_outer_radius = (ncoordgrid[0] > 1) ? coord_pos_min_tmin[0][1] : globals::rmax;
  if (!(first_shell_outer_radius > 1e4 * cellbound_tolerance(0.))) {
    fatal_crash(
        "model.txt: the first radial shell of the spherical grid has the outer radius {:g} [cm] at tmin, but it must "
        "be above {:g} [cm]",
        first_shell_outer_radius, 1e4 * cellbound_tolerance(0.));
  }

  // the cosine and the sine at an edge on an axis or in the equatorial plane are exact numbers
  const auto set_cos_sin_of_edges = [](const std::span<const double> angles, std::vector<double>& cos_of_edges,
                                       std::vector<double>& sin_of_edges) {
    cos_of_edges.resize(angles.size());
    sin_of_edges.resize(angles.size());
    for (auto i = 0Z; i < std::ssize(angles); i++) {
      cos_of_edges[i] = std::cos(angles[i]);
      sin_of_edges[i] = std::sin(angles[i]);
    }
    cos_of_edges.front() = 1.;
    sin_of_edges.front() = 0.;
    sin_of_edges.back() = 0.;
  };
  set_cos_sin_of_edges(prop_theta_edges, theta_edge_cos, theta_edge_sin);
  theta_edge_cos.back() = -1.;
  for (auto i = 0Z; i < std::ssize(prop_theta_edges); i++) {
    if (prop_theta_edges[i] == PI / 2) {
      theta_edge_cos[i] = 0.;
      theta_edge_sin[i] = 1.;
    }
  }
  if (model_mirror_z) {
    // the edges of the mirror half have the exact opposite cosines, thus the two cells of a mirror pair have the
    // same volume
    const auto num_half_edges = std::ssize(model_theta_edges);
    for (auto i = num_half_edges; i < std::ssize(prop_theta_edges); i++) {
      const auto mirrorindex = (2 * (num_half_edges - 1)) - i;
      theta_edge_cos[i] = -theta_edge_cos[mirrorindex];
      theta_edge_sin[i] = theta_edge_sin[mirrorindex];
    }
  }
  // A polar edge very near the z axis has the cosine of the axis. The cell between two edges with one cosine has no
  // volume, and its boundary has no cone.
  for (auto i = 1Z; i < std::ssize(theta_edge_cos); i++) {
    if (!(theta_edge_cos[i] < theta_edge_cos[i - 1])) {
      fatal_crash(
          "model.txt: the polar edges {:.17g} and {:.17g} [rad] have the same cosine, thus the polar cell between "
          "them is too thin for the precision of the numbers",
          prop_theta_edges[i - 1], prop_theta_edges[i]);
    }
  }
  set_cos_sin_of_edges(model_phi_edges, phi_edge_cos, phi_edge_sin);
  phi_edge_cos.back() = 1.;
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

    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D:
      // get_cellindex_from_pos() searches the cell edges of these grids itself
      break;
  }

  fatal_crash("Position {:g} on axis {} at time {:g} has no cell index", pos, axis, time);
}

// Convert a position vector from Cartesian xyz to the grid coordinate system
[[nodiscard]] constexpr auto get_gridcoords_from_xyz(const Vec3d& pos_xyz, const GridType gridtype) -> Vec3d {
  switch (gridtype) {
    case GridType::CARTESIAN3D:
      return pos_xyz;
    case GridType::CYLINDRICAL2D:
      return {std::sqrt(pow2(pos_xyz[0]) + pow2(pos_xyz[1])), pos_xyz[2], NAN};
    // only the radius of a spherical grid is a length. boundary_distance() treats the angles of a 2D spherical and a
    // 3D spherical grid separately
    case GridType::SPHERICAL1D:
    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D:
      return {vec_len(pos_xyz), NAN, NAN};
  }
  fatal_crash("Unknown grid type {}", static_cast<int>(gridtype));
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

    case GridType::SPHERICAL1D:
    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D: {
      const double v_radial = dot(pos_xyz, dir_xyz) / pktposgridcoord[0] * CLIGHT_PROP;
      return {v_radial, NAN, NAN};
    }
  }
  fatal_crash("Unknown grid type {}", static_cast<int>(gridtype));
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

// The first crossing of a polar or an azimuthal cell boundary of a spherical grid on a packet path: the distance
// [cm] to the crossing, or -1 if no such crossing lies ahead, and the polar and the azimuthal cell index after it.
struct AngularBoundaryCrossing {
  double distance;
  int index_theta;
  int index_phi;
};

// Find the first crossing of a polar or an azimuthal cell boundary on the path from pos along the unit vector dir,
// for a packet in the cell with the given polar and azimuthal index. The cones of constant theta and the half-planes
// of constant phi do not move with the homologous expansion, thus the result does not depend on the time.
[[nodiscard]] DEVICE_FUNC auto get_first_angular_boundary_crossing(const Vec3d& pos, const Vec3d& dir,
                                                                   const int index_theta, const int index_phi)
    -> AngularBoundaryCrossing {
  const int num_theta = ncoordgrid[1];
  const int num_phi = ncoordgrid[2];
  AngularBoundaryCrossing first_crossing{.distance = -1., .index_theta = index_theta, .index_phi = index_phi};
  const auto take_crossing_if_nearer = [&first_crossing](const double distance, const int index_theta_after,
                                                         const int index_phi_after) {
    if (distance >= 0. && (first_crossing.distance < 0. || distance < first_crossing.distance)) {
      first_crossing = {.distance = distance, .index_theta = index_theta_after, .index_phi = index_phi_after};
    }
  };

  // The angular momentum pos x dir of a unit direction has the length of the smallest distance of the path from the
  // centre. Its z component divided by the length of dir in the xy plane is the smallest distance from the z axis.
  const auto angmom = cross_prod(pos, dir);
  const double dir_xy_squared = pow2(dir[0]) + pow2(dir[1]);
  // every polar cell and every azimuthal cell of a shell meets at the centre, and every azimuthal cell meets at
  // the z axis. A path within this distance [cm] of one of them changes its cell there in one event
  const double path_through_axis_tolerance = cellbound_tolerance(0.);

  if (vec_len(angmom) <= path_through_axis_tolerance) {
    // The path goes through the centre, e.g. the radial path of a new packet. Its theta and phi are constant on each
    // side of the centre, and after the centre they are the angles of dir.
    const auto [theta_of_dir, phi_of_dir] = get_theta_phi_from_xyz(dir);
    const int index_theta_after = get_cellindex_from_lower_edges(coord_pos_min_tmin[1], theta_of_dir);
    const int index_phi_after = get_cellindex_from_lower_edges(coord_pos_min_tmin[2], phi_of_dir);
    if (index_theta_after != index_theta || index_phi_after != index_phi) {
      // A packet after the centre that is not in the cell of its direction changes cell at once. It is within the
      // tolerance of each boundary between the two cells, because the path stays that near to the ray of dir.
      take_crossing_if_nearer(std::max(-dot(pos, dir), 0.), index_theta_after, index_phi_after);
    }
    return first_crossing;
  }

  // the tolerance of a boundary that the packet is at, or just went through (see
  // is_boundary_overshoot_within_tolerance). An angle tolerance has no use near the centre, thus the test uses
  // the distance [cm] from the boundary
  const double tolerance = cellbound_tolerance(vec_len(pos));

  // the distance from the z axis, and its rate of change along the path
  const double rcyl = std::sqrt(pow2(pos[0]) + pow2(pos[1]));
  const double rcyl_rate = (rcyl > 0.) ? ((pos[0] * dir[0]) + (pos[1] * dir[1])) / rcyl : std::sqrt(dir_xy_squared);

  // The polar axis is not a cell boundary, thus the first cell has no lower boundary and the last cell has no upper
  // boundary. theta can have an extremum on the path, thus both boundaries get a test for each direction.
  const auto take_theta_boundary_if_nearer = [&]<BoundaryType boundarytype>(const int index_theta_after) {
    const bool is_upper = (boundarytype == BoundaryType::UPPER);
    const int edgeindex = is_upper ? index_theta + 1 : index_theta;
    // r * sin(theta - theta_edge): the distance of the packet from the cone, positive on the side of larger theta
    const double distance_from_cone = (theta_edge_cos[edgeindex] * rcyl) - (theta_edge_sin[edgeindex] * pos[2]);
    const double distance_from_cone_rate =
        (theta_edge_cos[edgeindex] * rcyl_rate) - (theta_edge_sin[edgeindex] * dir[2]);
    const bool moves_through_boundary = is_upper ? (distance_from_cone_rate > 0.) : (distance_from_cone_rate < 0.);
    const double distance = (moves_through_boundary && std::abs(distance_from_cone) <= tolerance)
                                ? 0.
                                : theta_cone_intersection<boundarytype>(pos, dir, theta_edge_cos[edgeindex]);
    take_crossing_if_nearer(distance, index_theta_after, index_phi);
  };
  if (index_theta < num_theta - 1) {
    take_theta_boundary_if_nearer.template operator()<BoundaryType::UPPER>(index_theta + 1);
  }
  if (index_theta > 0) {
    take_theta_boundary_if_nearer.template operator()<BoundaryType::LOWER>(index_theta - 1);
  }

  // A direction along the z axis keeps phi constant. A 2D spherical grid has one cell of phi.
  if (num_phi == 1 || dir_xy_squared <= 0.) {
    return first_crossing;
  }

  if (std::abs(angmom[2]) <= path_through_axis_tolerance * std::sqrt(dir_xy_squared)) {
    // The path goes through the z axis, where phi changes by pi. After the axis, phi is the angle of dir.
    const int index_phi_after = get_cellindex_from_lower_edges(coord_pos_min_tmin[2], get_theta_phi_from_xyz(dir)[1]);
    if (index_phi_after != index_phi) {
      const double distance_to_axis = -((pos[0] * dir[0]) + (pos[1] * dir[1])) / dir_xy_squared;
      take_crossing_if_nearer(std::max(distance_to_axis, 0.), index_theta, index_phi_after);
    }
    return first_crossing;
  }

  // phi is periodic: the half-plane phi = 2 pi of the last cell is the half-plane phi = 0 of the first cell
  const auto take_phi_boundary_if_nearer = [&]<BoundaryType boundarytype>(const int index_phi_after) {
    const bool is_upper = (boundarytype == BoundaryType::UPPER);
    const int edgeindex = is_upper ? index_phi + 1 : index_phi;
    const double cos_edge = phi_edge_cos[edgeindex];
    const double sin_edge = phi_edge_sin[edgeindex];
    // the distance of the packet from the full plane of the boundary, positive on the side of larger phi
    const double distance_from_plane = (cos_edge * pos[1]) - (sin_edge * pos[0]);
    const double distance_from_plane_rate = (cos_edge * dir[1]) - (sin_edge * dir[0]);
    const bool moves_through_boundary = is_upper ? (angmom[2] > 0. && distance_from_plane_rate > 0.)
                                                 : (angmom[2] < 0. && distance_from_plane_rate < 0.);
    // the other half of the plane has phi = phi_edge + pi
    const bool is_near_boundary_half = ((cos_edge * pos[0]) + (sin_edge * pos[1])) > 0.;
    const double distance =
        (moves_through_boundary && is_near_boundary_half && std::abs(distance_from_plane) <= tolerance)
            ? 0.
            : phi_halfplane_intersection<boundarytype>(pos, dir, cos_edge, sin_edge);
    take_crossing_if_nearer(distance, index_theta, index_phi_after);
  };
  take_phi_boundary_if_nearer.template operator()<BoundaryType::UPPER>((index_phi + 1) % num_phi);
  take_phi_boundary_if_nearer.template operator()<BoundaryType::LOWER>((index_phi + num_phi - 1) % num_phi);

  return first_crossing;
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

  if (model_mirror_z) {
    // the two propagation cells of a mirror pair have the same volume
    return 2 * get_propcell_volume_tmin(get_upper_propcell_of_mirrored_mgi(modelgridindex));
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

    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D: {
      // the integral of r^2 sin(theta) over the cell. A 2D spherical grid has one phi cell from 0 to 2 pi.
      const double radialpart = (pow3(get_cellcoordmax(cellindex, 0)) - pow3(get_cellcoordmin(cellindex, 0))) / 3.;
      const int index_theta = get_cellcoordindex(cellindex, 1);
      const double polarpart = theta_edge_cos[index_theta] - theta_edge_cos[index_theta + 1];
      return radialpart * polarpart * propcell_width_tmin(cellindex, 2);
    }
  }
  fatal_crash("Unknown propagation grid type {}", static_cast<int>(get_propgridtype()));
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

    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D: {
      // The volume element is r^2 dr dcos(theta) dphi, thus r^3, cos(theta), and phi have uniform distributions.
      // rng_uniform_pos() gives no zero and no one, thus the position is not on a cell boundary.
      const double r_inner = get_cellcoordmin(cellindex, 0);
      const double r_outer = get_cellcoordmax(cellindex, 0);
      const double radius = std::cbrt(std::lerp(pow3(r_outer), pow3(r_inner), rng_uniform_pos(rngstate)));
      const int index_theta = get_cellcoordindex(cellindex, 1);
      const double costheta =
          std::lerp(theta_edge_cos[index_theta], theta_edge_cos[index_theta + 1], rng_uniform_pos(rngstate));
      const double sintheta = std::sqrt((1. - costheta) * (1. + costheta));
      const double phi =
          std::lerp(get_cellcoordmin(cellindex, 2), get_cellcoordmax(cellindex, 2), rng_uniform_pos(rngstate));
      return {radius * sintheta * std::cos(phi), radius * sintheta * std::sin(phi), radius * costheta};
    }
  }
  fatal_crash("Unknown propagation grid type {}", static_cast<int>(get_propgridtype()));
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

namespace {

// Set the cell counts of a 2D spherical or a 3D spherical model from its cell edges. The third count of a 2D model
// is zero, as for a 2D cylindrical model.
void set_spherical_model_cell_counts(const GridType spherical_model_type) {
  ncoord_model = {
      static_cast<int>(std::ssize(model_vel_r_edges)) - 1,
      static_cast<int>(std::ssize(model_theta_edges)) - 1,
      (spherical_model_type == GridType::SPHERICAL3D) ? static_cast<int>(std::ssize(model_phi_edges)) - 1 : 0,
  };
  npts_model = static_cast<ptrdiff_t>(ncoord_model[0]) * ncoord_model[1] * std::max(1, ncoord_model[2]);
  if (npts_model > std::numeric_limits<int>::max()) {
    fatal_crash("model.txt: the cell edges give {} cells, which is more than the largest cell index", npts_model);
  }
}

// the largest difference [rad] between the last polar edge and pi, between the last azimuthal edge and 2 pi, and
// between a polar edge and pi / 2, that the reader of model.txt replaces by the exact angle
constexpr double model_angle_edge_tolerance = 1e-8;

// the largest difference of a cylindrical or a Cartesian cell edge from its place on a grid with cells of equal width,
// as a fraction of vmax. The model grid of these two types has cells of equal width only
constexpr double model_uniform_edge_tolerance = 1e-9;

// The smallest angle [rad] between a polar edge and the z axis. A polar edge is a cone with the cosine of the edge,
// and near the axis this cosine moves the cone by approximately 1.1e-16 / theta of the radius. Below this angle, the
// movement exceeds the position tolerance of a cell boundary, which is 1e-12 of the radius.
constexpr double model_min_polar_edge_angle_from_axis = 1e-4;

// Read the header of a model.txt of format version 2, after its first line "format_version 2". Each header line is
// a key and its values, and the keys can come in any sequence. The key columns gives the column names and ends the
// header. The other keys are coordinates, t_model_init_days, and one list of cell edges for each axis of the
// coordinates, e.g. vel_r_edges_kmps. A list has one value more than the count of the cells on its axis. The optional
// key cells_per_axis gives these counts again, in the sequence of the axes, e.g. "5 4 4". A key that
// this version of sn3d does not know gets a warning, thus a later version of sn3d can add a key to format version 2.
// The function gives the column names back to the reader as a header line of format version 1.
auto read_model_keyword_header(ModelFileReader& fmodel) -> GridType {
  using KeyLine = std::pair<std::string, std::string>;  // a key and the text after it, without a comment
  std::vector<KeyLine> keylines;
  const auto find_key = [&keylines](const std::string_view key) {
    return std::ranges::find(keylines, key, &KeyLine::first);
  };
  std::string line;
  while (find_key("columns") == keylines.end()) {
    if (!std::getline(fmodel.file, line)) {
      fatal_crash("model.txt: the header must end with the key columns, which gives the column names");
    }
    if (lineiscommentonly(line)) {
      continue;
    }
    auto text = std::string_view{line};
    text = text.substr(0, text.find('#'));
    text.remove_prefix(std::min(text.find_first_not_of(token_whitespace), text.size()));
    const auto key = text.substr(0, text.find_first_of(token_whitespace));
    if (find_key(key) != keylines.end()) {
      fatal_crash("model.txt: the header has the key {} two times", key);
    }
    keylines.emplace_back(std::string{key}, std::string{text.substr(key.size())});
  }

  const auto get_text_of_key = [&](const std::string_view key) -> std::string_view {
    for (const auto& [linekey, linetext] : keylines) {
      if (linekey == key) {
        auto text = std::string_view{linetext};
        text.remove_prefix(std::min(text.find_first_not_of(token_whitespace), text.size()));
        return text.substr(0, text.find_last_not_of(token_whitespace) + 1);
      }
    }
    fatal_crash("model.txt: the header has no key {}", key);
  };

  const auto get_numbers_of_key = [&](const std::string_view key) {
    auto remainder = get_text_of_key(key);
    std::vector<double> numbers;
    double number{NAN};
    while (parse_next_token(remainder, number)) {
      numbers.push_back(number);
    }
    remainder.remove_prefix(std::min(remainder.find_first_not_of(token_whitespace), remainder.size()));
    if (!remainder.empty()) {
      fatal_crash("model.txt: the key {} has an unreadable number at '{}'", key, remainder);
    }
    return numbers;
  };

  // An edge list increases. A list of r, rcyl, theta, or phi starts at exactly zero. An angular list ends at
  // last_edge_exact (pi or 2 pi)
  const auto require_valid_edges = [](std::vector<double>& edges, const std::string_view key,
                                      const bool first_edge_is_zero, const std::optional<double> last_edge_exact) {
    if (edges.size() < 2) {
      fatal_crash("model.txt: the key {} must have two or more values", key);
    }
    if (first_edge_is_zero && edges.front() != 0.) {
      fatal_crash("model.txt: the first value of {} is {:g}, but it must be 0", key, edges.front());
    }
    if (last_edge_exact.has_value()) {
      if (!(std::abs(edges.back() - last_edge_exact.value()) <= model_angle_edge_tolerance)) {
        fatal_crash("model.txt: the last value of {} is {:.17g}, but it must be {:.17g}", key, edges.back(),
                    last_edge_exact.value());
      }
      edges.back() = last_edge_exact.value();
    }
    for (auto i = 1Z; i < std::ssize(edges); i++) {
      if (!(edges[i] > edges[i - 1])) {
        fatal_crash("model.txt: the values of {} must increase, but {:.17g} follows {:.17g}", key, edges[i],
                    edges[i - 1]);
      }
      // The position tolerance of a cell boundary is 1e-12 of the radius (see cellbound_tolerance). A packet in an
      // angular cell that is thinner than this tolerance goes through the cell with no path length.
      if (last_edge_exact.has_value() && !(edges[i] - edges[i - 1] > model_angle_edge_tolerance)) {
        fatal_crash(
            "model.txt: the values {:.17g} and {:.17g} of {} are the edges of a cell that is too thin. Each angular "
            "cell must be wider than {:g} [rad]",
            edges[i - 1], edges[i], key, model_angle_edge_tolerance);
      }
    }
  };

  // the cell edges of a velocity axis [cm/s]
  const auto get_velocity_edges = [&](const std::string_view key, const bool first_edge_is_zero) {
    auto edges = get_numbers_of_key(key);
    require_valid_edges(edges, key, first_edge_is_zero, std::nullopt);
    for (auto& edge : edges) {
      edge *= 1.e5;  // km/s to cm/s
    }
    // two values in km/s that differ in the last bit can give the same velocity in cm/s
    for (auto i = 1Z; i < std::ssize(edges); i++) {
      if (!(edges[i] > edges[i - 1])) {
        fatal_crash(
            "model.txt: the values of {} at the indices {} and {} give the velocities {:.17g} and {:.17g} [cm/s]. The "
            "values must also increase in cm/s",
            key, i - 1, i, edges[i - 1], edges[i]);
      }
    }
    return edges;
  };

  // sn3d holds a cylindrical grid and a Cartesian grid with cells of equal width from first_edge to vmax [cm/s]
  const auto require_uniform_edges = [](const std::vector<double>& edges, const std::string_view key,
                                        const double first_edge, const double vmax) {
    const auto num_cells = std::ssize(edges) - 1;
    for (auto i = 0Z; i <= num_cells; i++) {
      const double uniform_edge = first_edge + ((vmax - first_edge) * static_cast<double>(i) / num_cells);
      if (!(std::abs(edges[i] - uniform_edge) <= model_uniform_edge_tolerance * vmax)) {
        fatal_crash(
            "model.txt: the values of {} must have equal steps from {:g} to {:g} [km/s], but the value at index {} is "
            "{:.17g}. sn3d reads a cylindrical grid or a Cartesian grid only with cells of equal width and the same "
            "vmax on each axis",
            key, first_edge / 1.e5, vmax / 1.e5, i, edges[i] / 1.e5);
      }
    }
  };

  const auto coordinates = get_text_of_key("coordinates");
  using CoordinatesName = std::pair<std::string_view, GridType>;
  constexpr std::array<CoordinatesName, 5> coordinates_names{
      CoordinatesName{"spherical1d", GridType::SPHERICAL1D}, CoordinatesName{"cylindrical2d", GridType::CYLINDRICAL2D},
      CoordinatesName{"cartesian3d", GridType::CARTESIAN3D}, CoordinatesName{"spherical2d", GridType::SPHERICAL2D},
      CoordinatesName{"spherical3d", GridType::SPHERICAL3D},
  };
  std::optional<GridType> named_model_type;
  for (const auto& [name, gridtype] : coordinates_names) {
    if (name == coordinates) {
      named_model_type = gridtype;
    }
  }
  if (!named_model_type.has_value()) {
    fatal_crash(
        "model.txt: the coordinates '{}' are not known. The choices are spherical1d, cylindrical2d, cartesian3d, "
        "spherical2d, and spherical3d",
        coordinates);
  }
  const auto detected_model_type = named_model_type.value();

  model_mirror_z = false;
  if (find_key("mirror_plane") != keylines.end()) {
    if (get_text_of_key("mirror_plane") != "z=0") {
      fatal_crash("model.txt: the key mirror_plane must have the value z=0, but it has '{}'",
                  get_text_of_key("mirror_plane"));
    }
    if (detected_model_type != GridType::CYLINDRICAL2D && detected_model_type != GridType::SPHERICAL2D &&
        detected_model_type != GridType::SPHERICAL3D) {
      fatal_crash(
          "model.txt: a {} model takes no key mirror_plane. The key applies to cylindrical2d, spherical2d, and "
          "spherical3d models",
          coordinates);
    }
    model_mirror_z = true;
    printlnlog(
        "model.txt has mirror symmetry in the plane z = 0. Each model cell fills its cell with z >= 0 and the mirror "
        "cell with z <= 0");
  }
  printlnlog("Detected {} model of format version 2", coordinates);

  const auto t_model_days = get_numbers_of_key("t_model_init_days");
  if (t_model_days.size() != 1 || !std::isfinite(t_model_days[0]) || t_model_days[0] <= 0.) {
    fatal_crash("model.txt: the key t_model_init_days must have one value above zero, but it has '{}'",
                get_text_of_key("t_model_init_days"));
  }
  t_model = t_model_days[0] * DAY;
  if (t_model > globals::tmin) {
    fatal_crash("The model snapshot time {} d is after tmin {} d", t_model_days[0], globals::tmin / DAY);
  }

  constexpr std::array<std::string_view, 7> edge_keys{
      "vel_r_edges_kmps", "theta_edges_radians", "phi_edges_radians", "vel_rcyl_edges_kmps",
      "vel_x_edges_kmps", "vel_y_edges_kmps",    "vel_z_edges_kmps",
  };
  std::vector<std::string_view> edge_keys_of_model;
  switch (detected_model_type) {
    case GridType::SPHERICAL1D:
      edge_keys_of_model = {"vel_r_edges_kmps"};
      break;
    case GridType::CYLINDRICAL2D:
      edge_keys_of_model = {"vel_rcyl_edges_kmps", "vel_z_edges_kmps"};
      break;
    case GridType::CARTESIAN3D:
      edge_keys_of_model = {"vel_x_edges_kmps", "vel_y_edges_kmps", "vel_z_edges_kmps"};
      break;
    case GridType::SPHERICAL2D:
      edge_keys_of_model = {"vel_r_edges_kmps", "theta_edges_radians"};
      break;
    case GridType::SPHERICAL3D:
      edge_keys_of_model = {"vel_r_edges_kmps", "theta_edges_radians", "phi_edges_radians"};
      break;
  }
  for (const auto& keyline : keylines) {
    const std::string_view key = keyline.first;
    if (key == "coordinates" || key == "cells_per_axis" || key == "mirror_plane" || key == "t_model_init_days" ||
        key == "cell_line_count" || key == "columns") {
      continue;
    }
    if (std::ranges::find(edge_keys, key) == edge_keys.end()) {
      printlnlog("[warning] model.txt: this version of sn3d does not know the key {}, and it ignores the key", key);
    } else if (std::ranges::find(edge_keys_of_model, key) == edge_keys_of_model.end()) {
      fatal_crash("model.txt: a {} model takes no key {}", coordinates, key);
    }
  }

  switch (detected_model_type) {
    case GridType::SPHERICAL1D: {
      model_vel_r_edges = get_velocity_edges("vel_r_edges_kmps", true);
      globals::vmax = model_vel_r_edges.back();
      // a version 1 file gives the outer velocity of each shell in a column
      vout_model.assign(model_vel_r_edges.begin() + 1, model_vel_r_edges.end());
      ncoord_model = {static_cast<int>(std::ssize(vout_model)), 0, 0};
      npts_model = std::ssize(vout_model);
      break;
    }

    case GridType::CYLINDRICAL2D: {
      const auto vel_rcyl_edges = get_velocity_edges("vel_rcyl_edges_kmps", true);
      const auto vel_z_edges = get_velocity_edges("vel_z_edges_kmps", false);
      globals::vmax = vel_rcyl_edges.back();
      require_uniform_edges(vel_rcyl_edges, "vel_rcyl_edges_kmps", 0., globals::vmax);
      // with a mirror plane, the z axis of the file starts at the plane z = 0
      require_uniform_edges(vel_z_edges, "vel_z_edges_kmps", model_mirror_z ? 0. : -globals::vmax, globals::vmax);
      ncoord_model = {
          static_cast<int>(std::ssize(vel_rcyl_edges)) - 1,
          static_cast<int>(std::ssize(vel_z_edges)) - 1,
          0,
      };
      npts_model = static_cast<ptrdiff_t>(ncoord_model[0]) * ncoord_model[1];
      break;
    }

    case GridType::CARTESIAN3D: {
      for (int axis = 0; axis < 3; axis++) {
        const auto vel_edges = get_velocity_edges(edge_keys_of_model[axis], false);
        if (axis == 0) {
          globals::vmax = vel_edges.back();
        }
        require_uniform_edges(vel_edges, edge_keys_of_model[axis], -globals::vmax, globals::vmax);
        ncoord_model[axis] = static_cast<int>(std::ssize(vel_edges)) - 1;
      }
      if (ncoord_model[1] != ncoord_model[0] || ncoord_model[2] != ncoord_model[0]) {
        fatal_crash(
            "model.txt: the lists vel_x_edges_kmps, vel_y_edges_kmps, and vel_z_edges_kmps give {} x {} x {} cells, "
            "but sn3d reads a Cartesian grid with the same count of cells on each axis",
            ncoord_model[0], ncoord_model[1], ncoord_model[2]);
      }
      npts_model = static_cast<ptrdiff_t>(ncoord_model[0]) * ncoord_model[1] * ncoord_model[2];
      break;
    }

    case GridType::SPHERICAL2D:
    case GridType::SPHERICAL3D: {
      model_vel_r_edges = get_velocity_edges("vel_r_edges_kmps", true);
      globals::vmax = model_vel_r_edges.back();

      model_theta_edges = get_numbers_of_key("theta_edges_radians");
      for (auto& theta_edge : model_theta_edges) {
        // the plane z = 0 is a cell boundary with an exact cosine of zero
        if (std::abs(theta_edge - (PI / 2)) <= model_angle_edge_tolerance) {
          theta_edge = PI / 2;
        }
      }
      // this check comes after the change to pi / 2, because two edges that become equal there make a cell with no
      // volume
      // with a mirror plane, the polar axis of the file ends at the plane z = 0
      require_valid_edges(model_theta_edges, "theta_edges_radians", true, model_mirror_z ? PI / 2 : PI);
      for (const double theta_edge : model_theta_edges) {
        if (theta_edge > 0. && theta_edge < PI &&
            std::min(theta_edge, PI - theta_edge) < model_min_polar_edge_angle_from_axis) {
          fatal_crash(
              "model.txt: the polar edge {:.17g} [rad] is less than {:g} rad from the z axis. sn3d takes a polar edge "
              "as "
              "a cone with the cosine of the edge, and the cosine of such an edge does not give its angle to the "
              "position tolerance of a cell boundary",
              theta_edge, model_min_polar_edge_angle_from_axis);
        }
      }

      if (detected_model_type == GridType::SPHERICAL3D) {
        model_phi_edges = get_numbers_of_key("phi_edges_radians");
        require_valid_edges(model_phi_edges, "phi_edges_radians", true, 2 * PI);
      } else {
        model_phi_edges = {0., 2 * PI};
      }
      set_spherical_model_cell_counts(detected_model_type);
      break;
    }
  }
  if (npts_model > std::numeric_limits<int>::max()) {
    fatal_crash("model.txt: the cell edges give {} cells, which is more than the largest cell index", npts_model);
  }

  // An edge list that lost a value or got one more value gives a different number of cells on its axis. The key is
  // optional.
  if (find_key("cells_per_axis") != keylines.end()) {
    const auto cells_per_axis = get_numbers_of_key("cells_per_axis");
    const auto num_axes = std::ssize(edge_keys_of_model);
    std::string edge_cell_counts;
    bool counts_agree = (std::ssize(cells_per_axis) == num_axes);
    for (int axis = 0; axis < num_axes; axis++) {
      edge_cell_counts += std::format("{}{}", (axis > 0) ? " " : "", ncoord_model[axis]);
      counts_agree = counts_agree && (cells_per_axis[axis] == ncoord_model[axis]);
    }
    if (!counts_agree) {
      fatal_crash("model.txt: the key cells_per_axis gives the cell counts '{}', but the cell edges give '{}'",
                  get_text_of_key("cells_per_axis"), edge_cell_counts);
    }
  }

  model_cell_line_count.reset();
  if (find_key("cell_line_count") != keylines.end()) {
    const auto counts = get_numbers_of_key("cell_line_count");
    if (counts.size() != 1 || counts[0] != std::floor(counts[0]) || counts[0] < 1 ||
        counts[0] > static_cast<double>(npts_model)) {
      fatal_crash(
          "model.txt: the key cell_line_count must have one integer from 1 to {}, which is the count of the cells "
          "that the cell edges give, but its text is '{}'",
          npts_model, get_text_of_key("cell_line_count"));
    }
    model_cell_line_count = static_cast<ptrdiff_t>(counts[0]);
  }

  const auto columns = get_text_of_key("columns");
  auto columns_remainder = columns;
  for (const std::string_view required_column : {"inputcellid", "rho"}) {
    columns_remainder.remove_prefix(
        std::min(columns_remainder.find_first_not_of(token_whitespace), columns_remainder.size()));
    const auto column = columns_remainder.substr(0, columns_remainder.find_first_of(token_whitespace));
    if (column != required_column) {
      fatal_crash("model.txt: the first two values of the key columns must be inputcellid and rho, but they are '{}'",
                  columns);
    }
    columns_remainder.remove_prefix(column.size());
  }
  // read_model_columns() takes the column names as a header line of format version 1
  fmodel.lines_read_ahead.push_back(std::format("#{}", columns));

  return detected_model_type;
}

}  // anonymous namespace

// read the input ejecta model from model.txt, auto-detecting its format version and its dimensionality (1D spherical,
// 2D cylindrical, 3D Cartesian, 2D spherical, or 3D spherical) and reading the per-cell densities,
// abundances, and any optional extra columns into the model grid
void read_ejecta_model() {
  auto fmodel = ModelFileReader{.file = istream_required("model.txt"), .lines_read_ahead = {}};
  std::string line;
  std::optional<GridType> detected_model_type{};

  // two integers on the first line of the model file
  int npts_0 = 0;  // total model points for 1D/3D, and number of points in r for 2D
  int npts_1 = 0;  // number of points in z for 2D
  // The first line of a file of format version 2 is "format_version 2". Such a file has a header of keys in place
  // of the header of numbers. A file of version 1 has no version line.
  if (!std::getline(fmodel.file, line)) {
    fatal_crash("model.txt is empty");
  }
  std::string firstword;
  std::string versiontext;
  std::string textafterversion;
  // the version line can end with a comment, and it has no other text after the version
  std::istringstream{line.substr(0, line.find('#'))} >> firstword >> versiontext >> textafterversion;
  model_format_version = (firstword == "format_version") ? 2 : 1;
  if (model_format_version == 2) {
    if (!textafterversion.empty()) {
      // the message gives each value after the key, without the comment
      auto versionvalues = std::string_view{line}.substr(0, line.find('#'));
      versionvalues.remove_prefix(versionvalues.find("format_version") + std::string_view{"format_version"}.size());
      versionvalues.remove_prefix(std::min(versionvalues.find_first_not_of(token_whitespace), versionvalues.size()));
      versiontext = std::string{versionvalues.substr(0, versionvalues.find_last_not_of(token_whitespace) + 1)};
    }
    if (versiontext != "2" || !textafterversion.empty()) {
      fatal_crash(
          "model.txt has the format version '{}'. This version of sn3d reads format version 2, and format version 1, "
          "which has no format_version line",
          versiontext);
    }
    detected_model_type = read_model_keyword_header(fmodel);
  } else {
    if (lineiscommentonly(line)) {
      if (!get_noncommentline(fmodel.file, line)) {
        fatal_crash("model.txt holds only comment lines. Its first line with a value must give the cell count");
      }
      firstword.clear();
      std::istringstream{line} >> firstword;
    }
    if (firstword == "format_version") {
      fatal_crash("model.txt: the line '{}' must be the first line of the file", line);
    }
    if (firstword == "coordinates") {
      fatal_crash("model.txt: the key coordinates needs the line 'format_version 2' as the first line of the file");
    }
    auto ssline = std::istringstream{line};
    ssline >> npts_0;
    if (npts_0 <= 0) {
      fatal_crash("model.txt: could not read a positive cell count from the first line '{}'", line);
    }
    if (ssline >> npts_1) {
      // second number on the line for 2D means the line was n_r n_z
      detected_model_type = GridType::CYLINDRICAL2D;
      printlnlog("Detected 2D model");
      npts_model = static_cast<ptrdiff_t>(npts_0) * npts_1;
    } else {
      npts_model = npts_0;
    }

    // Now read the time (in days) at which the model is specified.
    double t_model_days{NAN};
    if (!get_noncommentline(fmodel.file, line)) {
      fatal_crash("model.txt ends after the cell count, before the line of the snapshot time");
    }
    std::istringstream{line} >> t_model_days;
    // a failed extraction stores zero, which would zero all densities via the (t_model / tmin)^3 scaling
    if (!std::isfinite(t_model_days) || t_model_days <= 0.) {
      fatal_crash("model.txt: could not read a positive snapshot time in days from line '{}'", line);
    }
    t_model = t_model_days * DAY;
    if (t_model > globals::tmin) {
      fatal_crash("The model snapshot time {} d is after tmin {} d", t_model_days, globals::tmin / DAY);
    }

    // if the next line is a single float, it is the vmax (so 2D or 3D)
    // otherwise, it is the first line of the model or a header comment (so 1D)
    if (!fmodel.getline(line)) {
      fatal_crash("model.txt ends after the snapshot time, before the vmax line or the first cell");
    }
    if (!line.starts_with('#')) {
      double num_after_vmax{NAN};
      auto sslinevmax = std::istringstream{line};
      if ((sslinevmax >> globals::vmax) && !(sslinevmax >> num_after_vmax)) {
        // single value on the line is a vmax, so 2D or 3D
        // if it's not already know to be 2D (based on n_r n_z line), then it's 3D
        if (detected_model_type != GridType::CYLINDRICAL2D) {
          detected_model_type = GridType::CARTESIAN3D;
          printlnlog("Detected 3D model");
        }
      }
    }
    if (!detected_model_type.has_value()) {
      detected_model_type = GridType::SPHERICAL1D;
      printlnlog("Detected 1D model");
      fmodel.lines_read_ahead.push_back(line);
    }
  }

  model_type = detected_model_type;

  assert_always(modelgrid_input.empty());
  modelgrid_input = MPI_shared_array<ModelGridCellInput>(npts_model, ModelGridCellInput{});
  modelgrid_numpropcells.resize(npts_model, 0);
  nonemptymgi_of_mgi.resize(npts_model, -1);

  // geometry-specific setup of the model coordinate counts
  if (model_format_version == 2) {
    // read_model_keyword_header() set the cell counts from the cell edges
    printlnlog("model.txt gives {} cells on the first axis{}{}", ncoord_model[0],
               (ncoord_model[1] > 0) ? std::format(", {} on the second axis", ncoord_model[1]) : std::string{},
               (ncoord_model[2] > 0) ? std::format(", and {} on the third axis", ncoord_model[2]) : std::string{});
    if (get_modelgridtype() == GridType::CARTESIAN3D) {
      // for a 3D input model, the propagation cells will match the input cells exactly
      ncoordgrid = ncoord_model;
      ngrid = npts_model;
    }
  } else if (get_modelgridtype() == GridType::SPHERICAL1D) {
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
    if (static_cast<ptrdiff_t>(ncoord_model[0]) * ncoord_model[1] * ncoord_model[2] != npts_model) {
      fatal_crash("model.txt: a 3D model needs a cube number of cells, but it has {} cells", npts_model);
    }
    // for a 3D input model, the propagation cells will match the input cells exactly
    ncoordgrid = ncoord_model;
    ngrid = npts_model;
  }

  const auto [colnames, nucindexlist, one_line_per_cell] = read_model_columns(fmodel);

  // only global rank 0 parses the cell data. The node-shared arrays it fills are broadcast to the leaders of the
  // other nodes below, and the rank-local results are broadcast to all ranks, so that only one rank reads the
  // (potentially very large) file
  if (globals::my_rank == 0) {
    // 3D only: whether the input cell positions match the expected values in x-y-z or z-y-x
    // column order (set false when a mismatch is detected)
    bool posmatch_xyz = true;
    bool posmatch_zyx = true;

    if (model_format_version == 2) {
      first_input_cellid = 0;
    }
    ptrdiff_t cell_lines_read = 0;
    int mgi = 0;
    while (mgi < get_npts_model() && fmodel.getline(line)) {
      auto remainder = std::string_view{line};
      int cellnumberin = 0;
      double rho_tmodel{NAN};  // the cell density [g/cm3] at the model snapshot time t_model

      // parse the geometry-specific part of the cell line to get the density (and for 1D, the
      // outer velocity), and check that any cell midpoint coordinates match the expected values
      if (model_format_version == 2) {
        if (lineiscommentonly(line)) {
          continue;
        }
        // columns: cell number, density. The header gives the cell edges. The cell number is the index of the cell
        // in the grid of the edges, from zero, with the index of the first axis as the fastest index.
        if (!(parse_next_token(remainder, cellnumberin) && parse_next_token(remainder, rho_tmodel))) {
          fatal_crash("model.txt: a cell line must start with 2 numbers (inputcellid rho), but the line is: {}", line);
        }
        if (cellnumberin < mgi || cellnumberin >= get_npts_model()) {
          fatal_crash(
              "model.txt has the cell id {} after {} cells. The cell ids of format version 2 must increase, and the "
              "cell edges of the header give the ids 0 to {}",
              cellnumberin, mgi, get_npts_model() - 1);
        }
        // a cell with no line is empty
        for (; mgi < cellnumberin; mgi++) {
          set_rho_tmin(mgi, 0.F);
        }
        if (get_modelgridtype() == GridType::CARTESIAN3D) {
          // the reader of a version 1 file of a 3D model takes the density as a float, thus the same file in
          // version 2 gives the same density
          const auto rho_model_in = static_cast<float>(rho_tmodel);
          rho_tmodel = rho_model_in;
          if (min_den < 0. || min_den > rho_tmodel) {
            min_den = rho_tmodel;
          }
        }
      } else if (get_modelgridtype() == GridType::SPHERICAL1D) {
        // columns: cell number, outer velocity [km/s], log10 of the density
        double vout_kmps{NAN};
        double log_rho{NAN};
        if (!(parse_next_token(remainder, cellnumberin) && parse_next_token(remainder, vout_kmps) &&
              parse_next_token(remainder, log_rho))) {
          fatal_crash(
              "model.txt cell {}: expected at least 3 values (inputcellid vel_r_max_kmps log10rho) but could "
              "not parse line: {}",
              mgi, line);
        }
        vout_model[mgi] = vout_kmps * 1.e5;
        // the velocity grid is binary searched by int_index_lowerbound(), so it has to increase
        // outwards. Without this check a non-monotonic model.txt gives silently wrong cell lookups.
        if (mgi > 0 && !(vout_model[mgi] > vout_model[mgi - 1])) {
          fatal_crash("model.txt cell {} has vel_r_max_kmps {:g}, but it must be larger than {:g} of the cell before",
                      mgi, vout_kmps, vout_model[mgi - 1] / 1.e5);
        }
        rho_tmodel = (log_rho > -90) ? pow(10., log_rho) : 0.;
      } else if (get_modelgridtype() == GridType::CYLINDRICAL2D) {
        // columns: cell number, r midpoint, z midpoint, density
        float cell_r_in{NAN};
        float cell_z_in{NAN};
        if (!(parse_next_token(remainder, cellnumberin) && parse_next_token(remainder, cell_r_in) &&
              parse_next_token(remainder, cell_z_in) && parse_next_token(remainder, rho_tmodel))) {
          fatal_crash(
              "model.txt cell {} must start with 4 numbers (inputcellid pos_rcyl_mid pos_z_mid rho), but the line is: "
              "{}",
              mgi, line);
        }

        const int n_rcyl = (mgi % ncoord_model[0]);
        const double pos_r_cyl_mid = (n_rcyl + 0.5) * globals::vmax * t_model / ncoord_model[0];
        if (!(fabs((cell_r_in / pos_r_cyl_mid) - 1) < 1e-3)) {
          fatal_crash("model.txt cell {} has pos_rcyl_mid {:g}, but the expected value is {:g}", mgi, cell_r_in,
                      pos_r_cyl_mid);
        }
        const int n_z = (mgi / ncoord_model[0]);
        if (((2 * n_z) + 1) == ncoord_model[1]) {
          // an odd number of z rows puts the middle row at the origin, where a relative comparison
          // has nothing to divide by. Compare against the cell size instead. The row is identified
          // from the integer index rather than by testing the coordinate against zero, because the
          // expression below is not required to evaluate to exactly zero under -ffast-math.
          if (!(fabs(cell_z_in) < (1e-3 * globals::vmax * t_model / ncoord_model[1]))) {
            fatal_crash("model.txt cell {} has pos_z_mid {:g}, but the expected value is 0", mgi, cell_z_in);
          }
        } else {
          const double pos_z_mid = globals::vmax * t_model * (-1 + (2 * (n_z + 0.5) / ncoord_model[1]));
          if (!(fabs((cell_z_in / pos_z_mid) - 1) < 1e-3)) {
            fatal_crash("model.txt cell {} has pos_z_mid {:g}, but the expected value is {:g}", mgi, cell_z_in,
                        pos_z_mid);
          }
        }
      } else {
        // columns: cell number, x minimum, y minimum, z minimum, density. The position columns hold the corner
        // of the cell with the lowest coordinates. The check below permits a difference of half a cell width, so
        // it cannot tell a file of cell midpoints from a file of cell corners.
        std::array<float, 3> cellpos_in{};
        float rho_model_in{NAN};
        if (!(parse_next_token(remainder, cellnumberin) && parse_next_token(remainder, cellpos_in[0]) &&
              parse_next_token(remainder, cellpos_in[1]) && parse_next_token(remainder, cellpos_in[2]) &&
              parse_next_token(remainder, rho_model_in))) {
          fatal_crash(
              "model.txt cell {} must start with 5 numbers (inputcellid pos_x_min pos_y_min pos_z_min rho), but the "
              "line is: {}",
              mgi, line);
        }
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

      if (mgi == 0 && model_format_version == 1) {
        first_input_cellid = cellnumberin;
        printlnlog("model.txt: autodetected cell numbering starting from {}", first_input_cellid);
        if (first_input_cellid != 0 && first_input_cellid != 1) {
          fatal_crash("model.txt: the first cell id is {}, but it must be 0 or 1", first_input_cellid);
        }
      }
      if (cellnumberin != mgi + first_input_cellid) {
        fatal_crash("model.txt cell {} has the cell id {}, but the expected id is {}", mgi, cellnumberin,
                    mgi + first_input_cellid);
      }

      if (rho_tmodel < 0) {
        fatal_crash("model.txt cell {} has negative density {:g} [g/cm3] at t_model", mgi, rho_tmodel);
      }

      const bool keepcell = (rho_tmodel > 0);
      set_rho_tmin(mgi, static_cast<float>(rho_tmodel * pow3(t_model / globals::tmin)));

      read_model_radioabundances(fmodel, remainder, mgi, keepcell, colnames, nucindexlist, one_line_per_cell);

      cell_lines_read++;
      mgi++;
    }

    if (model_format_version == 2) {
      // a cell with no line is empty, also after the last line of the file
      for (; mgi < get_npts_model(); mgi++) {
        set_rho_tmin(mgi, 0.F);
      }
    }
    if (mgi != get_npts_model()) {
      fatal_crash("model.txt: found only {} cells instead of {} expected.", mgi, get_npts_model());
    }

    if (model_format_version == 2) {
      // the header of this format has no cell count, thus a cell after the last one is a sign of a wrong edge list
      while (fmodel.getline(line)) {
        if (!lineiscommentonly(line)) {
          fatal_crash("model.txt has more than {} cells, but the cell edges of its header give {} cells",
                      get_npts_model(), get_npts_model());
        }
      }
      if (model_cell_line_count.has_value() && cell_lines_read != model_cell_line_count.value()) {
        fatal_crash(
            "model.txt: the key cell_line_count gives {} cell lines, but the file has {}. The file is possibly not "
            "complete",
            model_cell_line_count.value(), cell_lines_read);
      }
    }

    if (get_modelgridtype() == GridType::SPHERICAL1D) {
      // for 1D models, vmax is the outer velocity of the last cell
      globals::vmax = vout_model[get_npts_model() - 1];
    } else if (get_modelgridtype() == GridType::CARTESIAN3D && model_format_version == 2) {
      printlnlog("minimum model density {:g} [g/cm3]", min_den);
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
    if (!model_otherstable_massfracs.empty()) {
      MPI_Bcast_safe(model_otherstable_massfracs, 0, globals::mpi_comm_internode);
    }
  }
  MPI_Barrier_allranks();

  assert_always(get_npts_model() ==
                std::max(1, ncoord_model[0]) * std::max(1, ncoord_model[1]) * std::max(1, ncoord_model[2]));
  printlnlog("npts_model: {}", get_npts_model());
  globals::rmax = globals::vmax * globals::tmin;
  printlnlog("vmax {:g} [cm/s] ({:.2f}c)", globals::vmax, globals::vmax / CLIGHT);
  if (!(globals::vmax < CLIGHT)) {
    fatal_crash("model.txt gives vmax = {:g} cm/s, but vmax must be below the speed of light", globals::vmax);
  }
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

void setup_spherical_grid_from_edges(const GridType spherical_grid_type, const std::span<const double> vel_r_edges,
                                     const std::span<const double> theta_edges,
                                     const std::span<const double> phi_edges) {
  assert_always(has_spherical_angle_axes(spherical_grid_type));
  model_type = spherical_grid_type;
  model_vel_r_edges.assign(vel_r_edges.begin(), vel_r_edges.end());
  model_theta_edges.assign(theta_edges.begin(), theta_edges.end());
  if (spherical_grid_type == GridType::SPHERICAL3D) {
    model_phi_edges.assign(phi_edges.begin(), phi_edges.end());
  } else {
    model_phi_edges = {0., 2 * PI};
  }
  set_spherical_model_cell_counts(spherical_grid_type);
  globals::vmax = model_vel_r_edges.back();
  globals::rmax = globals::vmax * globals::tmin;
  setup_grid_spherical();
}

// identify the cell index from an (x,y,z) position and a time.
[[nodiscard]] DEVICE_FUNC auto get_cellindex_from_pos(const Vec3d& pos, const double time) -> int {
  const auto prop_gridtype = get_propgridtype();
  if (has_spherical_angle_axes(prop_gridtype)) {
    const double radius = vec_len(pos);
    if (radius > (globals::vmax * time)) {
      // outside grid
      return -99;
    }
    const double trat = time / globals::tmin;
    const auto [theta, phi] = get_theta_phi_from_xyz(pos);
    const int index_r = get_cellindex_from_lower_edges(coord_pos_min_tmin[0], radius / trat);
    const int index_theta = get_cellindex_from_lower_edges(coord_pos_min_tmin[1], theta);
    const int index_phi = get_cellindex_from_lower_edges(coord_pos_min_tmin[2], phi);
    const int cellindex = index_r + (ncoordgrid[0] * (index_theta + (ncoordgrid[1] * index_phi)));

    // the division by trat can round the radius to the other side of a shell boundary
    const double radius_inner = get_cellcoordmin(cellindex, 0) * trat;
    const double radius_outer = get_cellcoordmax(cellindex, 0) * trat;
    assert_always(radius >= (radius_inner - cellbound_tolerance(radius_inner)));
    assert_always(radius <= (radius_outer + cellbound_tolerance(radius_outer)));
    return cellindex;
  }

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
    const double cellposmin = get_cellcoordmin(cellindex, d) / globals::tmin * time;
    const double cellposmax = get_cellcoordmax(cellindex, d) / globals::tmin * time;
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
    // allocate_nonemptymodelcells() precomputed the radius comparison per cell, because this
    // function runs for every packet step
    if (propcell_outside_escape_surface[cellindex] != 0) {
      return {0., -99};
    }
  }

  // the following vector are in grid coordinates, so either x,y,z (3D) or r (1D), or r_xy, z (2D)

  const auto pktposgridcoord = get_gridcoords_from_xyz(pos, prop_gridtype);

  // dir * CLIGHT_PROP converted from xyz to grid coordinates
  const auto pktvelgridcoord = get_gridcoords_vel_from_xyz_pos_dir(pos, dir, pktposgridcoord, prop_gridtype);

  std::array<int, 3> cellcoordidx{};
  std::array<double, 3> cellcoordmin{};
  std::array<double, 3> cellcoordmax{};
  for (int d = 0; d < get_ndim(prop_gridtype); d++) {
    cellcoordidx[d] = get_cellcoordindex(cellindex, d);
    cellcoordmin[d] = get_cellcoordmin(cellindex, d);
    cellcoordmax[d] = get_cellcoordmax(cellindex, d);
  }

  {
    // the boundaries of an angular axis do not move with the flow, thus the check is for the length axes
    for (int d = 0; d < get_num_length_axes(prop_gridtype); d++) {
      // pos_component_vel_relative_to_flow is constant along a ray with a given direction in Cartesian coordinates, but
      // for non-Cartesian coordinates, we still need to check at the current position whether the packet is
      // moving in the positive or negative direction in each grid coordinate direction relative to the homologous grid
      // flow, otherwise we might never enter the cell that we're supposed to be in
      const bool pos_component_vel_relative_to_flow = (pktvelgridcoord[d] * tstart) > pktposgridcoord[d];

      if constexpr (TESTMODE) {
        bool isoutside_error = false;
        double delta = 0.;
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
          fatal_crash(
              "timestep {}: packet outside coord {} {}{} boundary of cell {} by delta {:g}. vel {:g} initpos {:g} "
              "cellcoordmin {:g} cellcoordmax {:g} dir [{:g}, {:g}, {:g}] tmin {:g} s tstart {:g} s",
              globals::timestep, d, pos_component_vel_relative_to_flow ? '+' : '-', get_coordlabel(prop_gridtype, d),
              cellindex, delta, pktvelgridcoord[d], pktposgridcoord[d], cellcoordmin[d] / globals::tmin * tstart,
              cellcoordmax[d] / globals::tmin * tstart, dir[0], dir[1], dir[2], globals::tmin, tstart);
        }
      }
    }
  }

  double distance = std::numeric_limits<double>::max();
  int next_cellindex{-1};

  // A Cartesian axis d has two possible boundary crossings. For x, the packet trajectory is
  // x = x0 + (dir.x) * c * (t - tstart)
  // the boundaries follow
  // x+/- = x+/-(tmin) * (t/tmin)
  // so the crossing occurs when
  // t - tstart = (x0 - x+/-(tmin)/tmin * tstart) / (x+/-(tmin)/tmin - (dir.x)*c)
  // distance = c * (t - tstart)
  const auto take_cartesian_boundary_if_nearer = [&] [[gnu::always_inline]] (const int d) {
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

      if ((d_coordminboundary >= 0.) && (d_coordminboundary < distance)) {
        distance = d_coordminboundary;
        next_cellindex = (cellcoordidx[d] == 0) ? -99 : cellindex - get_coordcellindexstride(d);
      }
    }
  };

  // The radius from the origin is coordinate 0 of each spherical grid. Its two boundaries are spheres that expand.
  const auto take_radial_shell_boundary_if_nearer = [&] [[gnu::always_inline]] () {
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
  };

  if (prop_gridtype == GridType::SPHERICAL1D) {
    // the only coordinate is the radius from the origin
    take_radial_shell_boundary_if_nearer();
  } else if (has_spherical_angle_axes(prop_gridtype)) {
    // coordinate 0 is the radius, coordinate 1 is the polar angle theta, and coordinate 2 is the azimuthal angle phi
    take_radial_shell_boundary_if_nearer();

    const auto angular_crossing = get_first_angular_boundary_crossing(pos, dir, cellcoordidx[1], cellcoordidx[2]);
    if ((angular_crossing.distance >= 0.) && (angular_crossing.distance < distance)) {
      distance = angular_crossing.distance;
      next_cellindex = cellcoordidx[0] +
                       (ncoordgrid[0] * (angular_crossing.index_theta + (ncoordgrid[1] * angular_crossing.index_phi)));
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
      // is sampled at exactly costheta = +/-1 and the local velocity is parallel to z, or for a virtual
      // packet to an observer at costheta = +/-1), so the
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

    // the z axis is Cartesian
    take_cartesian_boundary_if_nearer(1);
  } else if (prop_gridtype == GridType::CARTESIAN3D) {
    for (int d = 0; d < 3; d++) {
      take_cartesian_boundary_if_nearer(d);
    }
  } else {
    assert_always(false);
  }

  if constexpr (FORCE_SPHERICAL_ESCAPE_SURFACE) {
    if (next_cellindex == -1) {
      // No receding cell boundary lies ahead of the packet, so it escapes through the slower
      // spherical escape surface. The comoving coordinates of a packet converge to a point with
      // speed CLIGHT. Only a cell whose fastest corner is above CLIGHT can hold this state. A
      // slower cell here is a sign of a defective boundary intersection, and the run then stops.
      double cornerspeed_squared = 0.;
      for (int d = 0; d < get_num_length_axes(prop_gridtype); d++) {
        cornerspeed_squared += pow2(std::max(std::abs(cellcoordmin[d]), std::abs(cellcoordmax[d])) / globals::tmin);
      }
      assert_always(cornerspeed_squared > pow2(CLIGHT_PROP));

      // Truncation must act on whole cells to keep the volumes and the packet statistics valid.
      // The trapped packet never leaves this cell. The arrival time is constant along a free ray,
      // so a mid-cell escape is valid only for a cell with no matter. A matter cell can trap
      // packets only when the cell diagonal spans from below vmax to above CLIGHT, so a finer
      // grid removes this stop.
      if (get_propcell_modelgridindex(cellindex) >= 0) {
        fatal_crash(
            "a packet cannot reach any boundary of matter cell {}, because every boundary recedes faster than light. "
            "The cell reaches from inside the escape surface to beyond the light speed, so the grid is too coarse. Use "
            "more propagation cells per axis.",
            cellindex);
      }

      const double r_surface = globals::rmax * tstart / globals::tmin;
      if (dot(pos, pos) >= pow2(r_surface)) {
        // A cell that straddles the sphere keeps positions outside the surface. The intersection
        // then lies behind the packet, so it escapes at its current position.
        distance = 0.;
      } else {
        const double speed = vec_len(dir) * CLIGHT_PROP;  // just in case dir is not normalised
        // A ray that grazes the surface within rounding has no valid forward intersection and
        // gives -1. The packet is then at the surface itself, so the distance clamps to zero.
        distance = std::max(expanding_shell_intersection<BoundaryType::UPPER>(pos, dir, speed, r_surface, tstart), 0.);
      }
      next_cellindex = -99;
    }
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
