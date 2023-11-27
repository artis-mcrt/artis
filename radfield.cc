#include "radfield.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_debye.h>

#include <algorithm>
#include <cmath>
#include <ctime>

#include "atomic.h"
#include "grid.h"
#include "sn3d.h"
#include "vectors.h"

namespace radfield {

static std::vector<double> J_normfactor;

// typedef enum
// {
//   FIT_DILUTE_BLACKBODY = 0,
//   FIT_CONSTANT = 1,
// } enum_bin_fit_type;

struct radfieldbin_solution {
  // these two parameters are used in the current timestep, but were calculated
  // from the values of J and nuJ in the previous timestep
  float W;    // dilution (scaling) factor
  float T_R;  // radiation temperature
  // enum_bin_fit_type fit_type;
};

struct radfieldbin {
  double J_raw;  // value needs to be multipled by J_normfactor to get the true value
  double nuJ_raw;
  int contribcount;
};

static double radfieldbin_nu_upper[RADFIELDBINCOUNT];  // array of upper frequency boundaries of bins
static struct radfieldbin *radfieldbins = nullptr;
static struct radfieldbin_solution *radfieldbin_solutions = nullptr;

#ifdef MPI_ON
MPI_Win win_radfieldbin_solutions = MPI_WIN_NULL;
MPI_Win win_prev_bfrate_normed = MPI_WIN_NULL;
#endif

// ** Detailed lines - Jblue_lu estimators for selected lines

struct Jb_lu_estimator {
  double value = 0.;
  int contribcount = 0;
};

// reallocate the detailed line arrays in units of BLOCKSIZEJBLUE
constexpr int BLOCKSIZEJBLUE = 128;
static int detailed_linecount = 0;

// array of indicies into the linelist[] array for selected lines
static int *detailed_lineindicies;

static struct Jb_lu_estimator **prev_Jb_lu_normed = nullptr;  // value from the previous timestep
static struct Jb_lu_estimator **Jb_lu_raw = nullptr;          // unnormalised estimator for the current timestep

// ** end detailed lines

static float *prev_bfrate_normed = nullptr;  // values from the previous timestep
static double *bfrate_raw = nullptr;         // unnormalised estimators for the current timestep

// expensive debugging mode to track the contributions to each bound-free rate estimator

static std::vector<double> J;  // after normalisation: [ergs/s/sr/cm2/Hz]
#ifdef DO_TITER
static std::vector<double> J_reduced_save;
#endif

// J and nuJ are accumulated and then normalised in-place
// i.e. be sure the normalisation has been applied (exactly once) before using the values here!
static std::vector<double> nuJ;
#ifdef DO_TITER
static std::vector<double> nuJ_reduced_save;
#endif

using enum_prefactor = enum {
  ONE = 0,
  TIMES_NU = 1,
};

using gsl_planck_integral_paras = struct {
  double T_R;
  enum_prefactor prefactor;
};

using gsl_T_R_solver_paras = struct {
  int modelgridindex;
  int binindex;
};

static FILE *radfieldfile = nullptr;

static inline auto get_bin_nu_upper(int binindex) -> double { return radfieldbin_nu_upper[binindex]; }

static void setup_bin_boundaries() {
  // double prev_nu_upper = nu_lower_first_initial;

  // choose between equally spaced in energy/frequency or wavelength (before bf edges shift boundaries around)
  const double delta_nu =
      (nu_upper_last_initial - nu_lower_first_initial) / (RADFIELDBINCOUNT - 1);  // - 1 for the top super bin
  // const double lambda_lower_first_initial = 1e8 * CLIGHT / nu_lower_first_initial;
  // const double lambda_upper_last_initial = 1e8 * CLIGHT / nu_upper_last_initial;
  // const double delta_lambda = (lambda_upper_last_initial - lambda_lower_first_initial) / RADFIELDBINCOUNT;

  for (int binindex = 0; binindex < RADFIELDBINCOUNT - 1; binindex++) {
    // radfieldbin_nu_upper[binindex] = 1e8 * CLIGHT / (lambda_lower_first_initial + (binindex + 1) * delta_lambda);
    radfieldbin_nu_upper[binindex] = nu_lower_first_initial + (binindex + 1) * delta_nu;

    // Align the bin edges with bound-free edges, except for the last one
    // if (binindex < RADFIELDBINCOUNT - 1)
    // {
    //   for (int i = 0; i < globals::nbfcontinua_ground; i++)
    //   {
    //     const double nu_edge = globals::groundcont[i].nu_edge;
    //     const double eV_edge = H * nu_edge / EV;
    //     const double angstrom_edge = 1e8 * CLIGHT / nu_edge;
    //     const int element = globals::groundcont[i].element;
    //     const int ion = globals::groundcont[i].ion;
    //     const int level = globals::groundcont[i].level;
    //     const int phixstargetindex = globals::groundcont[i].phixstargetindex;
    //
    //     const int Z = get_atomicnumber(element);
    //     const int ion_stage = get_ionstage(element, ion);
    //     const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
    //
    //     //printout("bf edge at %g, nu_lower_first %g, nu_upper_last %g\n",nu_edge,nu_lower_first,nu_upper_last);
    //
    //     // this is compares the highest and lowest bins to the bound-free list, only do it once, i.e. binindex == 0
    //     if (binindex == 0 && ((nu_edge < nu_lower_first_initial) || (nu_edge > nu_upper_last_initial)))
    //     {
    //       printout("Missed bf edge at %12.5e Hz (%6.2f eV, %6.1f A), nu_lower_first %11.5e Hz, nu_upper_last %11.5e
    //       Hz, Z=%d ion_stage %d level %d upperionlevel %d\n",
    //                nu_edge, eV_edge, angstrom_edge, nu_lower_first_initial, nu_upper_last_initial, Z, ion_stage,
    //                level, upperionlevel);
    //     }
    //
    //     const double bin_nu_upper = get_bin_nu_upper(binindex);
    //     if ((nu_edge > prev_nu_upper) && (nu_edge < bin_nu_upper))
    //     {
    //       printout("Shifting bin %d nu_upper from %12.5e Hz to bf edge at %12.5e Hz (%6.2f eV, %6.1f A) for Z=%d
    //       ion_stage %d level %d upperionlevel %d\n",
    //                binindex, bin_nu_upper, nu_edge, eV_edge, angstrom_edge, Z, ion_stage, level, upperionlevel);
    //       radfieldbin_nu_upper[binindex] = nu_edge;
    //     }
    //   }
    // }
    // prev_nu_upper = get_bin_nu_upper(binindex);
  }
  radfieldbin_nu_upper[RADFIELDBINCOUNT - 1] = nu_upper_superbin;  // very top end super bin
}

static void realloc_detailed_lines(const int new_size) {
  auto *newptr = static_cast<int *>(realloc(detailed_lineindicies, new_size * sizeof(int)));
  if (newptr == nullptr) {
    printout("ERROR: Not enough memory to reallocate detailed Jblue estimator line list\n");
    abort();
  }
  assert_always(newptr != nullptr);
  detailed_lineindicies = newptr;

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    if (grid::get_numassociatedcells(modelgridindex) > 0) {
      prev_Jb_lu_normed[modelgridindex] = static_cast<struct Jb_lu_estimator *>(
          realloc(prev_Jb_lu_normed[modelgridindex], new_size * sizeof(struct Jb_lu_estimator)));

      Jb_lu_raw[modelgridindex] = static_cast<struct Jb_lu_estimator *>(
          realloc(Jb_lu_raw[modelgridindex], new_size * sizeof(struct Jb_lu_estimator)));

      if (prev_Jb_lu_normed[modelgridindex] == nullptr || Jb_lu_raw[modelgridindex] == nullptr) {
        printout("ERROR: Not enough memory to reallocate detailed Jblue estimator list for cell %d.\n", modelgridindex);
        abort();
      }
    }
  }
}

static void add_detailed_line(const int lineindex)
// associate a Jb_lu estimator with a particular lineindex to be used
// instead of the general radiation field model
{
  if (detailed_linecount % BLOCKSIZEJBLUE == 0) {
    const int new_size = detailed_linecount + BLOCKSIZEJBLUE;
    realloc_detailed_lines(new_size);
  }

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    if (grid::get_numassociatedcells(modelgridindex) > 0) {
      prev_Jb_lu_normed[modelgridindex][detailed_linecount].value = 0;
      prev_Jb_lu_normed[modelgridindex][detailed_linecount].contribcount = 0;

      // zero_estimators should do the next part anyway, but just to be sure:
      Jb_lu_raw[modelgridindex][detailed_linecount].value = 0;
      Jb_lu_raw[modelgridindex][detailed_linecount].contribcount = 0;
    }
  }
  detailed_lineindicies[detailed_linecount] = lineindex;
  detailed_linecount++;
  // printout("Added Jblue estimator for lineindex %d count %d\n", lineindex, detailed_linecount);
}

void init(int my_rank, int ndo_nonempty)
// this should be called only after the atomic data is in memory
{
  const int nonempty_npts_model = grid::get_nonempty_npts_model();

  J_normfactor.resize(nonempty_npts_model + 1);
  J.resize(nonempty_npts_model + 1);

#ifdef DO_TITER
  J_reduced_save.resize(nonempty_npts_model + 1);
#endif

  // J and nuJ are accumulated and then normalised in-place
  // i.e. be sure the normalisation has been applied (exactly once) before using the values here!
  nuJ.resize(nonempty_npts_model + 1);
#ifdef DO_TITER
  nuJ.resize(nonempty_npts_model + 1);
#endif

  prev_Jb_lu_normed =
      static_cast<struct Jb_lu_estimator **>(malloc((grid::get_npts_model() + 1) * sizeof(struct Jb_lu_estimator *)));
  Jb_lu_raw =
      static_cast<struct Jb_lu_estimator **>(malloc((grid::get_npts_model() + 1) * sizeof(struct Jb_lu_estimator *)));

  detailed_linecount = 0;

  detailed_lineindicies = nullptr;
  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    prev_Jb_lu_normed[modelgridindex] = nullptr;
    Jb_lu_raw[modelgridindex] = nullptr;
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    for (int i = 0; i < globals::nlines; i++) {
      const int element = globals::linelist[i].elementindex;
      const int Z = get_atomicnumber(element);
      if (Z == 26) {
        const int lowerlevel = globals::linelist[i].lowerlevelindex;
        // const int upperlevel = linelist[i].upperlevelindex;
        // const int ion = linelist[i].ionindex;
        // const int ionstage = get_ionstage(element, ion);
        const double A_ul = einstein_spontaneous_emission(i);

        bool addline = false;
        // if (ionstage == 1 && lowerlevel == 6 && upperlevel == 55)
        //   addline = true;
        // else if (ionstage == 1 && lowerlevel == 10 && upperlevel == 104)
        //   addline = true;
        // else if (ionstage == 1 && lowerlevel == 10 && upperlevel == 112)
        //   addline = true;
        // else if (ionstage == 2 && lowerlevel == 9 && upperlevel == 64)
        //   addline = true;

        if (lowerlevel <= 15 && A_ul > 0.) {  // ionstage <= 3 && A_ul > 1e3 &&
          addline = true;
        }

        if (addline) {
          // printout("Adding Jblue estimator for lineindex %d Z=%02d ionstage %d lower %d upper %d A_ul %g\n",
          //          i, Z, ionstage, lowerlevel, upperlevel, A_ul);
          add_detailed_line(i);
        }
      }
    }

    // shrink the detailed line list in case detailed_linecount isn't a multiple of BLOCKSIZEJBLUE
    // (important for saving memory if there are a lot of grid cells)
    realloc_detailed_lines(detailed_linecount);

    // these are probably sorted anyway because the previous loop goes in ascending
    // lineindex. But this sorting step is quick and makes sure that the
    // binary searching later will work correctly
    std::sort(detailed_lineindicies, detailed_lineindicies + detailed_linecount);
  }

  printout("There are %d lines with detailed Jblue_lu estimators.\n", detailed_linecount);

  printout("DETAILED_BF_ESTIMATORS %s", DETAILED_BF_ESTIMATORS_ON ? "ON" : "OFF");
  if (DETAILED_BF_ESTIMATORS_ON) {
    printout(" from timestep %d\n", DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP);
  } else {
    printout("\n");
  }

  if (MULTIBIN_RADFIELD_MODEL_ON) {
    printout("The multibin radiation field is being used from timestep %d onwards.\n", FIRST_NLTE_RADFIELD_TIMESTEP);

    printout("Initialising multibin radiation field with %d bins from (%.2f eV, %6.1f A) to (%.2f eV, %6.1f A)\n",
             RADFIELDBINCOUNT, H * nu_lower_first_initial / EV, 1e8 * CLIGHT / nu_lower_first_initial,
             H * nu_upper_last_initial / EV, 1e8 * CLIGHT / nu_upper_last_initial);
    if (ndo_nonempty > 0) {
      char filename[MAXFILENAMELENGTH];
      snprintf(filename, MAXFILENAMELENGTH, "radfield_%.4d.out", my_rank);
      assert_always(radfieldfile == nullptr);
      radfieldfile = fopen_required(filename, "w");
      fprintf(radfieldfile, "%8s %15s %8s %11s %11s %9s %9s %9s %9s %9s %12s\n", "timestep", "modelgridindex",
              "bin_num", "nu_lower", "nu_upper", "nuJ", "J", "J_nu_avg", "ncontrib", "T_R", "W");
      fflush(radfieldfile);
    }

    setup_bin_boundaries();

    const size_t mem_usage_bins = nonempty_npts_model * RADFIELDBINCOUNT * sizeof(struct radfieldbin);
    radfieldbins =
        static_cast<struct radfieldbin *>(malloc(nonempty_npts_model * RADFIELDBINCOUNT * sizeof(struct radfieldbin)));

    const size_t mem_usage_bin_solutions = nonempty_npts_model * RADFIELDBINCOUNT * sizeof(struct radfieldbin_solution);

#ifdef MPI_ON
    {
      int my_rank_cells = nonempty_npts_model / globals::node_nprocs;
      // rank_in_node 0 gets any remainder
      if (globals::rank_in_node == 0) {
        my_rank_cells += nonempty_npts_model - (my_rank_cells * globals::node_nprocs);
      }
      MPI_Aint size = my_rank_cells * RADFIELDBINCOUNT * sizeof(struct radfieldbin_solution);
      int disp_unit = sizeof(struct radfieldbin_solution);
      MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node, &radfieldbin_solutions,
                              &win_radfieldbin_solutions);

      MPI_Win_shared_query(win_radfieldbin_solutions, 0, &size, &disp_unit, &radfieldbin_solutions);
    }
#else
    {
      radfieldbin_solutions = static_cast<struct radfieldbin_solution *>(
          malloc(nonempty_npts_model * RADFIELDBINCOUNT * sizeof(struct radfieldbin_solution)));
    }
#endif

    printout("[info] mem_usage: radiation field bin accumulators for non-empty cells occupy %.3f MB\n",
             mem_usage_bins / 1024. / 1024.);
    printout(
        "[info] mem_usage: radiation field bin solutions for non-empty cells occupy %.3f MB (node shared memory)\n",
        mem_usage_bin_solutions / 1024. / 1024.);
  } else {
    printout("The radiation field model is a full-spectrum fit to a single dilute blackbody TR & W.\n");
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
#ifdef MPI_ON
    {
      int my_rank_cells = nonempty_npts_model / globals::node_nprocs;
      // rank_in_node 0 gets any remainder
      if (globals::rank_in_node == 0) {
        my_rank_cells += nonempty_npts_model - (my_rank_cells * globals::node_nprocs);
      }
      MPI_Aint size = my_rank_cells * globals::nbfcontinua * sizeof(float);
      int disp_unit = sizeof(float);
      MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node, &prev_bfrate_normed,
                              &win_prev_bfrate_normed);
      MPI_Win_shared_query(win_prev_bfrate_normed, 0, &size, &disp_unit, &prev_bfrate_normed);
    }
#else
    {
      prev_bfrate_normed = static_cast<float *>(malloc(nonempty_npts_model * globals::nbfcontinua * sizeof(float)));
    }
#endif
    printout("[info] mem_usage: detailed bf estimators for non-empty cells occupy %.3f MB (node shared memory)\n",
             nonempty_npts_model * globals::nbfcontinua * sizeof(float) / 1024. / 1024.);

    bfrate_raw = static_cast<double *>(malloc(nonempty_npts_model * globals::nbfcontinua * sizeof(double)));

    printout("[info] mem_usage: detailed bf estimator acculumators for non-empty cells occupy %.3f MB\n",
             nonempty_npts_model * globals::nbfcontinua * sizeof(double) / 1024. / 1024.);
  }

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    if (grid::get_numassociatedcells(modelgridindex) > 0) {
      zero_estimators(modelgridindex);

      if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
        if (globals::rank_in_node == 0) {
          const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
          for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
            const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
            radfieldbin_solutions[mgibinindex].W = -1.;
            radfieldbin_solutions[mgibinindex].T_R = -1.;
          }
        }
      }
    }
  }
}

/// Initialise estimator arrays which hold the last time steps values (used to damp out
/// fluctuations over timestep iterations if DO_TITER is defined) to -1.
void initialise_prev_titer_photoionestimators() {
  for (int n = 0; n < grid::get_npts_model(); n++) {
#ifdef DO_TITER
    const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
    J_reduced_save[nonemptymgi] = -1.;
#endif
#ifdef DO_TITER
    nuJ_reduced_save[n] = -1.;
    globals::ffheatingestimator_save[n] = -1.;
    globals::colheatingestimator_save[n] = -1.;
#endif
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++) {
#ifdef DO_TITER
        globals::gammaestimator_save[get_ionestimindex(n, element, ion)] = -1.;
        if constexpr (USE_LUT_BFHEATING) {
          globals::bfheatingestimator_save[get_ionestimindex(n, element, ion)] = -1.;
        }
#endif
      }
    }
  }
}

auto get_Jblueindex(const int lineindex) -> int
// returns -1 if the line does not have a Jblue estimator
{
  // slow linear search
  // for (int i = 0; i < detailed_linecount; i++)
  // {
  //   if (detailed_lineindicies[i] == lineindex)
  //     return i;
  // }

  if constexpr (!DETAILED_LINE_ESTIMATORS_ON) {
    return -1;
  }

  // use a binary search, assuming the list is sorted

  int low = 0;
  int high = detailed_linecount;
  while (low <= high) {
    const int mid = low + ((high - low) / 2);
    if (detailed_lineindicies[mid] < lineindex) {
      low = mid + 1;
    } else if (detailed_lineindicies[mid] > lineindex) {
      high = mid - 1;
    } else {
      return mid;
    }
  }

  // const int element = linelist[lineindex].elementindex;
  // const int ion = linelist[lineindex].ionindex;
  // const int lower = linelist[lineindex].lowerlevelindex;
  // const int upper = linelist[lineindex].upperlevelindex;
  // printout("Could not find lineindex %d among %d items (Z=%02d ionstage %d lower %d upper %d)\n",
  //          lineindex, detailed_linecount, get_atomicnumber(element), get_ionstage(element, ion), lower, upper);

  return -1;
}

auto get_Jb_lu(const int modelgridindex, const int jblueindex) -> double {
  assert_always(jblueindex >= 0);
  assert_always(jblueindex < detailed_linecount);
  return prev_Jb_lu_normed[modelgridindex][jblueindex].value;
}

auto get_Jb_lu_contribcount(const int modelgridindex, const int jblueindex) -> int {
  assert_always(jblueindex >= 0);
  assert_always(jblueindex < detailed_linecount);
  return prev_Jb_lu_normed[modelgridindex][jblueindex].contribcount;
}

static auto get_bin_J(int modelgridindex, int binindex) -> double
// get the normalised J_nu
{
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  assert_testmodeonly(J_normfactor[nonemptymgi] > 0.0);
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
  return radfieldbins[mgibinindex].J_raw * J_normfactor[nonemptymgi];
}

static auto get_bin_nuJ(int modelgridindex, int binindex) -> double {
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  assert_testmodeonly(J_normfactor[nonemptymgi] > 0.0);
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
  return radfieldbins[mgibinindex].nuJ_raw * J_normfactor[nonemptymgi];
}

static inline auto get_bin_nu_bar(int modelgridindex, int binindex) -> double
// importantly, this is average beween the current and previous timestep
{
  const double nuJ_sum = get_bin_nuJ(modelgridindex, binindex);
  const double J_sum = get_bin_J(modelgridindex, binindex);
  return nuJ_sum / J_sum;
}

static inline auto get_bin_nu_lower(int binindex) -> double {
  if (binindex > 0) {
    return radfieldbin_nu_upper[binindex - 1];
  }
  return nu_lower_first_initial;
}

static inline auto get_bin_contribcount(int modelgridindex, int binindex) -> int {
  const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
  return radfieldbins[mgibinindex].contribcount;
}

static inline auto get_bin_W(int modelgridindex, int binindex) -> float {
  const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
  return radfieldbin_solutions[mgibinindex].W;
}

static inline auto get_bin_T_R(int modelgridindex, int binindex) -> float {
  const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
  return radfieldbin_solutions[mgibinindex].T_R;
}

static inline auto select_bin(double nu) -> int {
  if (nu < get_bin_nu_lower(0)) {
    return -2;  // out of range, nu lower than lowest bin's lower boundary
  }

  // find the lowest frequency bin with radfieldbin_nu_upper > nu
  auto *bin = std::upper_bound(&radfieldbin_nu_upper[0], &radfieldbin_nu_upper[RADFIELDBINCOUNT], nu);
  const int binindex = std::distance(&radfieldbin_nu_upper[0], bin);
  if (binindex >= RADFIELDBINCOUNT) {
    // out of range, nu higher than highest bin's upper boundary
    return -1;
  }

  return binindex;
}

void write_to_file(int modelgridindex, int timestep) {
  assert_always(MULTIBIN_RADFIELD_MODEL_ON);
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
#ifdef OPENMP_MT_ON
#pragma omp critical(out_file)
  {
#endif

    int totalcontribs = 0;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      totalcontribs += get_bin_contribcount(modelgridindex, binindex);
    }

    for (int binindex = -1 - detailed_linecount; binindex < RADFIELDBINCOUNT; binindex++) {
      double nu_lower = 0.0;
      double nu_upper = 0.0;
      double nuJ_out = 0.0;
      double J_out = 0.0;
      float T_R = 0.0;
      float W = 0.0;
      double J_nu_bar = 0.0;
      int contribcount = 0;

      const bool skipoutput = false;

      if (binindex >= 0) {
        nu_lower = get_bin_nu_lower(binindex);
        nu_upper = get_bin_nu_upper(binindex);
        nuJ_out = get_bin_nuJ(modelgridindex, binindex);
        J_out = get_bin_J(modelgridindex, binindex);
        T_R = get_bin_T_R(modelgridindex, binindex);
        W = get_bin_W(modelgridindex, binindex);
        J_nu_bar = J_out / (nu_upper - nu_lower);
        contribcount = get_bin_contribcount(modelgridindex, binindex);
      } else if (binindex == -1) {  // bin -1 is the full spectrum fit
        nuJ_out = nuJ[nonemptymgi];
        J_out = J[nonemptymgi];
        T_R = grid::get_TR(modelgridindex);
        W = grid::get_W(modelgridindex);
        contribcount = totalcontribs;
      } else  // use binindex < -1 for detailed line Jb_lu estimators
      {
        const int jblueindex = -2 - binindex;  // -2 is the first detailed line, -3 is the second, etc
        const int lineindex = detailed_lineindicies[jblueindex];
        const double nu_trans = globals::linelist[lineindex].nu;
        nu_lower = nu_trans;
        nu_upper = nu_trans;
        nuJ_out = -1.;
        J_out = -1.;
        T_R = -1.;
        W = -1.;
        J_nu_bar = prev_Jb_lu_normed[modelgridindex][jblueindex].value,
        contribcount = prev_Jb_lu_normed[modelgridindex][jblueindex].contribcount;

        // if (J_nu_bar <= 0.)
        // {
        //   skipoutput = true;
        // }
      }

      if (!skipoutput) {
        fprintf(radfieldfile, "%d %d %d %11.5e %11.5e %9.3e %9.3e %9.3e %d %9.1f %12.5e\n", timestep, modelgridindex,
                binindex, nu_lower, nu_upper, nuJ_out, J_out, J_nu_bar, contribcount, T_R, W);
      }
    }
    fflush(radfieldfile);
#ifdef OPENMP_MT_ON
  }
#endif
}

void close_file() {
  if (radfieldfile != nullptr) {
    fclose(radfieldfile);
    radfieldfile = nullptr;
  }

  if (MULTIBIN_RADFIELD_MODEL_ON) {
    free(radfieldbins);
#ifdef MPI_ON
    if (win_radfieldbin_solutions != MPI_WIN_NULL) {
      MPI_Win_free(&win_radfieldbin_solutions);
    }
#else
    if (radfieldbin_solutions != nullptr) {
      free(radfieldbin_solutions);
    }
#endif
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    free(bfrate_raw);
#ifdef MPI_ON
    if (win_radfieldbin_solutions != MPI_WIN_NULL) {
      MPI_Win_free(&win_prev_bfrate_normed);
    }
#else
    if (prev_bfrate_normed != nullptr) {
      free(prev_bfrate_normed);
    }
#endif
  }
}

void zero_estimators(int modelgridindex)
// set up the new bins and clear the estimators in preparation
// for a timestep
{
  if (grid::get_numassociatedcells(modelgridindex) == 0) {
    return;
  }

  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    assert_always(bfrate_raw != nullptr);
    if (grid::get_numassociatedcells(modelgridindex) > 0) {
      for (int i = 0; i < globals::nbfcontinua; i++) {
        bfrate_raw[nonemptymgi * globals::nbfcontinua + i] = 0.;
      }
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    assert_always(Jb_lu_raw != nullptr);
    assert_always(Jb_lu_raw[modelgridindex] != nullptr);
    for (int i = 0; i < detailed_linecount; i++) {
      Jb_lu_raw[modelgridindex][i].value = 0.;
      Jb_lu_raw[modelgridindex][i].contribcount = 0.;
    }
  }

  J[nonemptymgi] = 0.;  // this is required even if FORCE_LTE is on
  nuJ[nonemptymgi] = 0.;

  if (MULTIBIN_RADFIELD_MODEL_ON) {
    // printout("radfield: zeroing estimators in %d bins in cell %d\n",RADFIELDBINCOUNT,modelgridindex);

    assert_always(radfieldbins != nullptr);
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
      radfieldbins[mgibinindex].J_raw = 0.0;
      radfieldbins[mgibinindex].nuJ_raw = 0.0;
      radfieldbins[mgibinindex].contribcount = 0;
    }
  }
  set_J_normfactor(modelgridindex, -1.0);
}

static void update_bfestimators(const int nonemptymgi, const double distance_e_cmf, const double nu_cmf,
                                const struct packet *const pkt_ptr) {
  assert_testmodeonly(DETAILED_BF_ESTIMATORS_ON);
  assert_always(bfrate_raw != nullptr);

  if (distance_e_cmf == 0) {
    return;
  }

  const int nbfcontinua = globals::nbfcontinua;
  const double dopplerfactor = doppler_packet_nucmf_on_nurf(pkt_ptr->pos, pkt_ptr->dir, pkt_ptr->prop_time);
  // const double dopplerfactor = 1.;

  const int tid = get_thread_num();
  const double distance_e_cmf_over_nu =
      distance_e_cmf / nu_cmf * dopplerfactor;  // TODO: Luke: why did I put a doppler factor here?
  for (int allcontindex = 0; allcontindex < nbfcontinua; allcontindex++) {
    const double nu_edge = globals::allcont_nu_edge[allcontindex];
    const double nu_max_phixs = nu_edge * last_phixs_nuovernuedge;  // nu of the uppermost point in the phixs table

    if (nu_cmf >= nu_edge && nu_cmf <= nu_max_phixs) {
      safeadd(bfrate_raw[nonemptymgi * nbfcontinua + allcontindex],
              globals::phixslist[tid].gamma_contr[allcontindex] * distance_e_cmf_over_nu);

    } else if (nu_cmf < nu_edge) {
      // list is sorted by nu_edge, so all remaining will have nu_cmf < nu_edge
      break;
    }
  }
}

void update_estimators(const int modelgridindex, const double distance_e_cmf, const double nu_cmf,
                       const struct packet *const pkt_ptr) {
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);

  safeadd(J[nonemptymgi], distance_e_cmf);
  safeadd(nuJ[nonemptymgi], distance_e_cmf * nu_cmf);

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    update_bfestimators(nonemptymgi, distance_e_cmf, nu_cmf, pkt_ptr);
  }

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    const int binindex = select_bin(nu_cmf);

    if (binindex >= 0) {
      const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
      safeadd(radfieldbins[mgibinindex].J_raw, distance_e_cmf);
      safeadd(radfieldbins[mgibinindex].nuJ_raw, distance_e_cmf * nu_cmf);
      safeincrement(radfieldbins[mgibinindex].contribcount);
    }
    // else
    // {
    //   printout("WARNING: radfield::update_estimators dropping packet contribution for nu_cmf %g\n",
    //            nu_cmf);
    //   printout("           modelgridindex %d binindex %d nu_lower_first %g nu_upper_last %g \n",
    //            modelgridindex, binindex, nu_lower_first, get_bin_nu_upper(modelgridindex,RADFIELDBINCOUNT - 1));
    // }
  }
}

void update_lineestimator(const int modelgridindex, const int lineindex, const double increment) {
  if constexpr (!DETAILED_LINE_ESTIMATORS_ON) {
    return;
  }

  const int jblueindex = get_Jblueindex(lineindex);
  if (jblueindex >= 0) {
    Jb_lu_raw[modelgridindex][jblueindex].value += increment;
    Jb_lu_raw[modelgridindex][jblueindex].contribcount += 1;
    // const int lineindex = detailed_lineindicies[jblueindex];
    // printout(" increment cell %d lineindex %d Jb_lu_raw %g prev_Jb_lu_normed %g radfield(nu_trans) %g\n",
    //       modelgridindex, lineindex, Jb_lu_raw[modelgridindex][jblueindex],
    //       prev_Jb_lu_normed[modelgridindex][jblueindex].value, radfield(linelist[lineindex].nu, modelgridindex));
  }
}

auto radfield(double nu, int modelgridindex) -> double
// returns mean intensity J_nu [ergs/s/sr/cm2/Hz]
{
  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    if (globals::timestep >= FIRST_NLTE_RADFIELD_TIMESTEP) {
      const int binindex = select_bin(nu);
      if (binindex >= 0) {
        const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
        const struct radfieldbin_solution *const bin = &radfieldbin_solutions[mgibinindex];
        if (bin->W >= 0.) {
          const double J_nu = dbb(nu, bin->T_R, bin->W);
          return J_nu;
        }
      } else {  // binindex < 0
        // if (nu > get_bin_nu_upper(RADFIELDBINCOUNT - 1))
        // {
        //   // undiluted LTE blueward of the bins
        //   const double J_nu_LTE = dbb(nu, grid::get_Te(modelgridindex), 1.0);
        //   return J_nu_LTE;
        // }
        // else
        //   return 0; // no radfield redwards of the bins
        // printout("WARNING: Radfield modelgridindex %d binindex %d nu %g nu_lower_first %g nu_upper_last %g \n",
        //         modelgridindex, binindex, nu, nu_lower_first, nu_upper_last);
      }
      return 0.;
    }
  }

  const float T_R_fullspec = grid::get_TR(modelgridindex);
  const float W_fullspec = grid::get_W(modelgridindex);
  const double J_nu_fullspec = dbb(nu, T_R_fullspec, W_fullspec);
  return J_nu_fullspec;
}

constexpr auto gsl_integrand_planck(const double nu, void *paras) -> double {
  const double T_R = (static_cast<gsl_planck_integral_paras *>(paras))->T_R;
  const enum_prefactor prefactor = (static_cast<gsl_planck_integral_paras *>(paras))->prefactor;

  double integrand = TWOHOVERCLIGHTSQUARED * std::pow(nu, 3) / (std::expm1(HOVERKB * nu / T_R));

  if (prefactor == TIMES_NU) {
    integrand *= nu;
  }

  return integrand;
}

static auto planck_integral(double T_R, double nu_lower, double nu_upper, enum_prefactor prefactor) -> double {
  double integral = 0.;

  double error = 0.;
  const double epsrel = 1e-10;
  const double epsabs = 0.;

  gsl_planck_integral_paras intparas = {.T_R = T_R, .prefactor = prefactor};

  const gsl_function F_planck = {.function = &gsl_integrand_planck, .params = &intparas};

  gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);
  const int status = gsl_integration_qag(&F_planck, nu_lower, nu_upper, epsabs, epsrel, GSLWSIZE, GSL_INTEG_GAUSS61,
                                         gslworkspace, &integral, &error);
  if (status != 0) {
    printout("planck_integral integrator status %d, GSL_FAILURE= %d. Integral value %g, setting to zero.\n", status,
             GSL_FAILURE, integral);
    integral = 0.;
  }
  gsl_set_error_handler(previous_handler);

  return integral;
}

static auto planck_integral_analytic(double T_R, double nu_lower, double nu_upper, enum_prefactor prefactor) -> double {
  double integral = 0.;

  if (prefactor == TIMES_NU) {
    const double debye_upper = gsl_sf_debye_4(HOVERKB * nu_upper / T_R) * pow(nu_upper, 4);
    const double debye_lower = gsl_sf_debye_4(HOVERKB * nu_lower / T_R) * pow(nu_lower, 4);
    integral = TWOHOVERCLIGHTSQUARED * (debye_upper - debye_lower) * T_R / HOVERKB / 4.;
  } else {
    const double debye_upper = gsl_sf_debye_3(HOVERKB * nu_upper / T_R) * pow(nu_upper, 3);
    const double debye_lower = gsl_sf_debye_3(HOVERKB * nu_lower / T_R) * pow(nu_lower, 3);
    integral = TWOHOVERCLIGHTSQUARED * (debye_upper - debye_lower) * T_R / HOVERKB / 3.;

    if (integral == 0.) {
      // double upperexp = exp(HOVERKB * nu_upper / T_R);
      // double upperint = - pow(nu_upper,4) / 4
      //                   + pow(nu_upper,3) * log(1 - upperexp) / HOVERKB
      //                   + 3 * pow(nu_upper,2) * polylog(2,upperexp) / pow(HOVERKB,2)
      //                   - 6 * nu_upper * polylog(3,upperexp) / pow(HOVERKB,3)
      //                   + 6 * polylog(4,upperexp) / pow(HOVERKB,4);
      // double lowerexp = exp(HOVERKB * nu_lower / T_R);
      // double lowerint = - pow(nu_lower,4) / 4
      //                   + pow(nu_lower,3) * log(1 - lowerexp) / HOVERKB
      //                   + 3 * pow(nu_lower,2) * polylog(2,lowerexp) / pow(HOVERKB,2)
      //                   - 6 * nu_lower * polylog(3,lowerexp) / pow(HOVERKB,3)
      //                   + 6 * polylog(4,lowerexp) / pow(HOVERKB,4);
      // double integral2 = TWOHOVERCLIGHTSQUARED * (upperint - lowerint);

      // printout("planck_integral_analytic is zero. debye_upper %g debye_lower %g. Test alternative %g\n",
      //          debye_upper,debye_lower,integral2);
    }
  }

  return integral;
}

static auto delta_nu_bar(double T_R, void *paras) -> double
// difference between the average nu and the average nu of a Planck function
// at temperature T_R, in the frequency range corresponding to a bin
{
  const int modelgridindex = (static_cast<gsl_T_R_solver_paras *>(paras))->modelgridindex;
  const int binindex = (static_cast<gsl_T_R_solver_paras *>(paras))->binindex;

  const double nu_lower = get_bin_nu_lower(binindex);
  const double nu_upper = get_bin_nu_upper(binindex);

  const double nu_bar_estimator = get_bin_nu_bar(modelgridindex, binindex);

  const double nu_times_planck_numerical = planck_integral(T_R, nu_lower, nu_upper, TIMES_NU);
  const double planck_integral_numerical = planck_integral(T_R, nu_lower, nu_upper, ONE);
  const double nu_bar_planck_T_R = nu_times_planck_numerical / planck_integral_numerical;

  /*double nu_times_planck_integral = planck_integral_analytic(T_R, nu_lower, nu_upper, TIMES_NU);
  double planck_integral_result = planck_integral_analytic(T_R, nu_lower, nu_upper, ONE);
  double nu_bar_planck = nu_times_planck_integral / planck_integral_result;

  //printout("nu_bar %g nu_bar_planck(T=%g) %g\n",nu_bar,T_R,nu_bar_planck);

  if (!std::isfinite(nu_bar_planck))
  {
    double nu_times_planck_numerical = planck_integral(T_R, nu_lower, nu_upper, TIMES_NU);
    double planck_integral_numerical = planck_integral(T_R, nu_lower, nu_upper, ONE);
    double nu_bar_planck_numerical = nu_times_planck_numerical / planck_integral_numerical;

    printout("planck_integral_analytic is %g. Replacing with numerical result of
  %g.\n",nu_bar_planck,nu_bar_planck_numerical); nu_bar_planck = nu_bar_planck_numerical;
  }*/

  const double delta_nu_bar = nu_bar_planck_T_R - nu_bar_estimator;

  if (!std::isfinite(delta_nu_bar)) {
    printout(
        "delta_nu_bar is %g. nu_bar_planck_T_R %g nu_times_planck_numerical %g planck_integral_numerical %g "
        "nu_bar_estimator %g\n",
        delta_nu_bar, nu_bar_planck_T_R, nu_times_planck_numerical, planck_integral_numerical, nu_bar_estimator);
  }

  return delta_nu_bar;
}

static auto find_T_R(int modelgridindex, int binindex) -> float {
  double T_R = 0.0;

  gsl_T_R_solver_paras paras;
  paras.modelgridindex = modelgridindex;
  paras.binindex = binindex;

  /// Check whether the equation has a root in [T_min,T_max]
  double delta_nu_bar_min = delta_nu_bar(T_R_min, &paras);
  double delta_nu_bar_max = delta_nu_bar(T_R_max, &paras);

  // printout("find_T_R: bin %4d delta_nu_bar(T_R_min) %g, delta_nu_bar(T_R_max) %g\n",
  //          binindex, delta_nu_bar_min,delta_nu_bar_max);

  if (!std::isfinite(delta_nu_bar_min) || !std::isfinite(delta_nu_bar_max)) {
    delta_nu_bar_max = delta_nu_bar_min = -1;
  }

  if (delta_nu_bar_min * delta_nu_bar_max < 0) {
    /// If there is a root in the interval, solve for T_R

    const double epsrel = 1e-4;
    const double epsabs = 0.;
    const int maxit = 100;

    gsl_function find_T_R_f;
    find_T_R_f.function = &delta_nu_bar;
    find_T_R_f.params = &paras;

    /// one dimensional gsl root solver, bracketing type
    gsl_root_fsolver *T_R_solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(T_R_solver, &find_T_R_f, T_R_min, T_R_max);
    int status = 0;
    for (int iteration_num = 0; iteration_num <= maxit; iteration_num++) {
      gsl_root_fsolver_iterate(T_R_solver);
      T_R = gsl_root_fsolver_root(T_R_solver);

      const double T_R_lower = gsl_root_fsolver_x_lower(T_R_solver);
      const double T_R_upper = gsl_root_fsolver_x_upper(T_R_solver);
      status = gsl_root_test_interval(T_R_lower, T_R_upper, epsabs, epsrel);

      // printout("find_T_R: bin %4d iter %d, T_R is between %7.1f and %7.1f, guess %7.1f, delta_nu_bar %g, status
      // %d\n",
      //          binindex,iteration_num,T_R_lower,T_R_upper,T_R,delta_nu_bar(T_R,&paras),status);
      if (status != GSL_CONTINUE) {
        break;
      }
    }

    if (status == GSL_CONTINUE) {
      printout("[warning] find_T_R: T_R did not converge within %d iterations\n", maxit);
    }

    gsl_root_fsolver_free(T_R_solver);
  } else if (delta_nu_bar_max < 0) {
    /// Thermal balance equation always negative ===> T_R = T_min
    /// Calculate the rates again at this T_e to print them to file
    T_R = T_R_max;
    printout("find_T_R: cell %d bin %4d no solution in interval, clamping to T_R_max=%g\n", modelgridindex, binindex,
             T_R_max);
  } else {
    T_R = T_R_min;
    printout("find_T_R: cell %d bin %4d no solution in interval, clamping to T_R_min=%g\n", modelgridindex, binindex,
             T_R_min);
  }

  return T_R;
}

static void set_params_fullspec(const int modelgridindex, const int timestep) {
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  const double nubar = nuJ[nonemptymgi] / J[nonemptymgi];
  if (!std::isfinite(nubar) || nubar == 0.) {
    printout("[warning] T_R estimator infinite in cell %d, keep T_R, T_J, W of last timestep. J = %g. nuJ = %g\n",
             modelgridindex, J[nonemptymgi], nuJ[nonemptymgi]);
  } else {
    float T_J = pow(J[nonemptymgi] * PI / STEBO, 1 / 4.);
    if (T_J > MAXTEMP) {
      printout("[warning] temperature estimator T_J = %g exceeds T_max %g in cell %d. Setting T_J = T_max!\n", T_J,
               MAXTEMP, modelgridindex);
      T_J = MAXTEMP;
    } else if (T_J < MINTEMP) {
      printout("[warning] temperature estimator T_J = %g below T_min %g in cell %d. Setting T_J = T_min!\n", T_J,
               MINTEMP, modelgridindex);
      T_J = MINTEMP;
    }
    grid::set_TJ(modelgridindex, T_J);

    float T_R = H * nubar / KB / 3.832229494;
    if (T_R > MAXTEMP) {
      printout("[warning] temperature estimator T_R = %g exceeds T_max %g in cell %d. Setting T_R = T_max!\n", T_R,
               MAXTEMP, modelgridindex);
      T_R = MAXTEMP;
    } else if (T_R < MINTEMP) {
      printout("[warning] temperature estimator T_R = %g below T_min %g in cell %d. Setting T_R = T_min!\n", T_R,
               MINTEMP, modelgridindex);
      T_R = MINTEMP;
    }
    grid::set_TR(modelgridindex, T_R);

    const float W = J[nonemptymgi] * PI / STEBO / pow(T_R, 4);
    grid::set_W(modelgridindex, W);

    printout(
        "Full-spectrum fit radfield for cell %d at timestep %d: J %g, nubar %5.1f Angstrom, T_J %g, T_R %g, W %g\n",
        modelgridindex, timestep, J[nonemptymgi], 1e8 * CLIGHT / nubar, T_J, T_R, W);
  }
}

void fit_parameters(int modelgridindex, int timestep)
// finds the best fitting W and temperature parameters in each spectral bin
// using J and nuJ
{
  set_params_fullspec(modelgridindex, timestep);

  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    if (J_normfactor[nonemptymgi] <= 0) {
      printout("radfield: FATAL J_normfactor = %g in cell %d at call to fit_parameters", J_normfactor[nonemptymgi],
               modelgridindex);
      abort();
    }

    double J_bin_sum = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      J_bin_sum += get_bin_J(modelgridindex, binindex);
    }

    printout("radfield bins sum to J of %g (%.1f%% of total J).\n", J_bin_sum, 100. * J_bin_sum / J[nonemptymgi]);
    printout("radfield: Finding parameters for %d bins...\n", RADFIELDBINCOUNT);

    double J_bin_max = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      const double J_bin = get_bin_J(modelgridindex, binindex);
      J_bin_max = std::max(J_bin_max, J_bin);
    }

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      const double nu_lower = get_bin_nu_lower(binindex);
      const double nu_upper = get_bin_nu_upper(binindex);
      const double J_bin = get_bin_J(modelgridindex, binindex);
      float T_R_bin = -1.0;
      double W_bin = -1.0;
      const int contribcount = get_bin_contribcount(modelgridindex, binindex);

      if (contribcount > 0) {
        // // enum_bin_fit_type bin_fit_type = radfieldbin_solutions[modelgridindex][binindex].fit_type;
        // if (bin_fit_type == FIT_DILUTE_BLACKBODY)
        {
          T_R_bin = find_T_R(modelgridindex, binindex);

          if (binindex == RADFIELDBINCOUNT - 1) {
            const auto T_e = grid::get_Te(modelgridindex);
            printout("    replacing bin %d T_R %7.1f with cell T_e = %7.1f\n", binindex,
                     get_bin_T_R(modelgridindex, binindex), T_e);
            T_R_bin = T_e;
          }

          double planck_integral_result = planck_integral(T_R_bin, nu_lower, nu_upper, ONE);
          //          printout("planck_integral(T_R=%g, nu_lower=%g, nu_upper=%g) = %g\n", T_R_bin, nu_lower,
          //          nu_upper, planck_integral_result);

          W_bin = J_bin / planck_integral_result;

          if (W_bin > 1e4) {
            //            printout("T_R_bin %g, nu_lower %g, nu_upper %g\n", T_R_bin, nu_lower, nu_upper);
            printout("W %g too high, trying setting T_R of bin %d to %g. J_bin %g planck_integral %g\n", W_bin,
                     binindex, T_R_max, J_bin, planck_integral_result);
            planck_integral_result = planck_integral(T_R_max, nu_lower, nu_upper, ONE);
            W_bin = J_bin / planck_integral_result;
            if (W_bin > 1e4) {
              printout("W still very high, W=%g. Zeroing bin...\n", W_bin);
              T_R_bin = -99.0;
              W_bin = 0.;
            } else {
              printout("new W is %g. Continuing with this value\n", W_bin);
              T_R_bin = T_R_max;
            }
          }
        }
        // else if (bin_fit_type == FIT_CONSTANT)
        // {
        //   T_R_bin = -1.;
        //   W_bin = J_bin / (nu_upper - nu_lower);
        // }
        // else
        // {
        //   printout("fit_parameters: unknown fit type %d for bin %d\n", bin_fit_type, binindex);
        //   T_R_bin = -1.;
        //   W_bin = -1.;
        // }
      } else {
        T_R_bin = 0.;
        W_bin = 0.;
      }
      // else
      // {
      //   T_R_bin = -1;
      //   W_bin = -1;
      // }
      const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
      radfieldbin_solutions[mgibinindex].T_R = T_R_bin;
      radfieldbin_solutions[mgibinindex].W = W_bin;
    }

    // double prev_nu_upper = nu_lower_first_initial;
    // for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    // {
    //   const double J_bin = get_bin_J(modelgridindex,binindex);
    //   const double T_R_bin = get_bin_T_R(modelgridindex,binindex);
    //   const double W_bin = get_bin_W(modelgridindex,binindex);
    //   const int contribcount = get_bin_contribcount(modelgridindex, binindex);
    //   const double bin_nu_upper = get_bin_nu_upper(binindex);
    //   const double nubar = get_bin_nu_bar(modelgridindex, binindex);
    //
    //   printout("bin %4d (lambda %7.1f Å to %7.1f Å): contribcount %5d J %7.1e T_R %8.1f W %12.5e lambdabar %7.1f
    //   Å\n",
    //          binindex, 1e8 * CLIGHT / prev_nu_upper, 1e8 * CLIGHT / bin_nu_upper, contribcount, J_bin, T_R_bin,
    //          W_bin, 1e8 * CLIGHT / nubar);
    //
    //  prev_nu_upper = get_bin_nu_upper(binindex);
    // }

    write_to_file(modelgridindex, timestep);
  }
}

void set_J_normfactor(int modelgridindex, double normfactor) {
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  J_normfactor[nonemptymgi] = normfactor;
}

void normalise_J(const int modelgridindex, const double estimator_normfactor_over4pi) {
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  assert_always(std::isfinite(J[nonemptymgi]));
  J[nonemptymgi] *= estimator_normfactor_over4pi;
  for (int i = 0; i < detailed_linecount; i++) {
    prev_Jb_lu_normed[modelgridindex][i].value = Jb_lu_raw[modelgridindex][i].value * estimator_normfactor_over4pi;
    prev_Jb_lu_normed[modelgridindex][i].contribcount = Jb_lu_raw[modelgridindex][i].contribcount;
  }
}

void normalise_bf_estimators(const int modelgridindex, const double estimator_normfactor_over_H) {
  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    printout("normalise_bf_estimators for cell %d with factor %g\n", modelgridindex, estimator_normfactor_over_H);
    const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
    assert_always(nonemptymgi >= 0);
    for (int i = 0; i < globals::nbfcontinua; i++) {
      const int mgibfindex = nonemptymgi * globals::nbfcontinua + i;
      prev_bfrate_normed[mgibfindex] = bfrate_raw[mgibfindex] * estimator_normfactor_over_H;
    }
  }
}

static auto get_bfcontindex(const int element, const int lowerion, const int lower, const int phixstargetindex) -> int {
  // simple linear search seems to be faster than the binary search
  // possibly because lower frequency transitions near start of list are more likely to be called?
  const auto &matchbf = std::find_if(globals::allcont, globals::allcont + globals::nbfcontinua, [=](const auto &bf) {
    return (bf.element == element) && (bf.ion == lowerion) && (bf.level == lower) &&
           (bf.phixstargetindex == phixstargetindex);
  });

  const int bfcontindex = std::distance(globals::allcont, matchbf);

  if (bfcontindex < globals::nbfcontinua) {
    return bfcontindex;
  }
  // not found in the continua list
  return -1;
}

auto get_bfrate_estimator(const int element, const int lowerion, const int lower, const int phixstargetindex,
                          const int modelgridindex) -> double {
  if constexpr (!DETAILED_BF_ESTIMATORS_ON) {
    return -1;
  } else {
    const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
    const int allcontindex = get_bfcontindex(element, lowerion, lower, phixstargetindex);
    if (allcontindex >= 0) {
      return prev_bfrate_normed[nonemptymgi * globals::nbfcontinua + allcontindex];
    }

    printout("no bf rate for element Z=%d ion_stage %d lower %d phixstargetindex %d\n", get_atomicnumber(element),
             get_ionstage(element, lowerion), lower, phixstargetindex);
    return -1.;
  }
}

void normalise_nuJ(const int modelgridindex, const double estimator_normfactor_over4pi) {
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  assert_always(std::isfinite(nuJ[nonemptymgi]));
  nuJ[nonemptymgi] *= estimator_normfactor_over4pi;
}

auto get_T_J_from_J(const int modelgridindex) -> double {
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  const double T_J = pow(J[nonemptymgi] * PI / STEBO, 1. / 4.);
  if (!std::isfinite(T_J)) {
    /// keep old value of T_J
    printout("[warning] get_T_J_from_J: T_J estimator infinite in cell %d, use value of last timestep\n",
             modelgridindex);
    return grid::get_TR(modelgridindex);
  }
  /// Make sure that T is in the allowed temperature range.
  if (T_J > MAXTEMP) {
    printout("[warning] get_T_J_from_J: T_J would be %.1f > MAXTEMP. Clamping to MAXTEMP = %.0f K\n", T_J, MAXTEMP);
    return MAXTEMP;
  }
  if (T_J < MINTEMP) {
    printout("[warning] get_T_J_from_J: T_J would be %.1f < MINTEMP. Clamping to MINTEMP = %.0f K\n", T_J, MINTEMP);
    return MINTEMP;
  }
  return T_J;
}

#ifdef DO_TITER
void titer_J(const int modelgridindex) {
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  if (J_reduced_save[nonemptymgi] >= 0) {
    J[nonemptymgi] = (J[nonemptymgi] + J_reduced_save[nonemptymgi]) / 2;
  }
  J_reduced_save[nonemptymgi] = J[nonemptymgi];
}

void titer_nuJ(const int modelgridindex) {
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  if (nuJ_reduced_save[nonemptymgi] >= 0) {
    nuJ[nonemptymgi] = (nuJ[nonemptymgi] + nuJ_reduced_save[nonemptymgi]) / 2;
  }
  nuJ_reduced_save[nonemptymgi] = nuJ[nonemptymgi];
}
#endif

#ifdef MPI_ON
void reduce_estimators()
// reduce and broadcast (allreduce) the estimators for J and nuJ in all bins
{
  const int nonempty_npts_model = grid::get_nonempty_npts_model();

  MPI_Allreduce(MPI_IN_PLACE, J.data(), nonempty_npts_model, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, nuJ.data(), nonempty_npts_model, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    MPI_Allreduce(MPI_IN_PLACE, bfrate_raw, nonempty_npts_model * globals::nbfcontinua, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
  }

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    const time_t sys_time_start_reduction = time(nullptr);
    printout("Reducing binned radiation field estimators");
    assert_always(radfieldbins != nullptr);

    for (int nonemptymgi = 0; nonemptymgi < nonempty_npts_model; nonemptymgi++) {
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
        const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
        MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[mgibinindex].J_raw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[mgibinindex].nuJ_raw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[mgibinindex].contribcount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      }
    }
    const int duration_reduction = time(nullptr) - sys_time_start_reduction;
    printout(" (took %d s)\n", duration_reduction);
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    const time_t sys_time_start_reduction = time(nullptr);
    printout("Reducing detailed line estimators");

    for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
      if (grid::get_numassociatedcells(modelgridindex) > 0) {
        for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
          MPI_Allreduce(MPI_IN_PLACE, &Jb_lu_raw[modelgridindex][jblueindex].value, 1, MPI_DOUBLE, MPI_SUM,
                        MPI_COMM_WORLD);
          MPI_Allreduce(MPI_IN_PLACE, &Jb_lu_raw[modelgridindex][jblueindex].contribcount, 1, MPI_INT, MPI_SUM,
                        MPI_COMM_WORLD);
        }
      }
    }
    const int duration_reduction = time(nullptr) - sys_time_start_reduction;
    printout(" (took %d s)\n", duration_reduction);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void do_MPI_Bcast(const int modelgridindex, const int root, int root_node_id)
// broadcast computed radfield results including parameters
// from the cells belonging to root process to all processes
{
  if (grid::get_numassociatedcells(modelgridindex) == 0) {
    return;
  }

  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  MPI_Bcast(&J_normfactor[nonemptymgi], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
      if (globals::rank_in_node == 0) {
        MPI_Bcast(&radfieldbin_solutions[mgibinindex].W, 1, MPI_FLOAT, root_node_id, globals::mpi_comm_internode);
        MPI_Bcast(&radfieldbin_solutions[mgibinindex].T_R, 1, MPI_FLOAT, root_node_id, globals::mpi_comm_internode);
      }
      MPI_Bcast(&radfieldbins[mgibinindex].J_raw, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&radfieldbins[mgibinindex].nuJ_raw, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&radfieldbins[mgibinindex].contribcount, 1, MPI_INT, root, MPI_COMM_WORLD);
    }
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    if (globals::rank_in_node == 0) {
      MPI_Bcast(&prev_bfrate_normed[nonemptymgi * globals::nbfcontinua], globals::nbfcontinua, MPI_FLOAT, root_node_id,
                globals::mpi_comm_internode);
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
      MPI_Bcast(&prev_Jb_lu_normed[modelgridindex][jblueindex].value, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&prev_Jb_lu_normed[modelgridindex][jblueindex].contribcount, 1, MPI_INT, root, MPI_COMM_WORLD);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
}
#endif

void write_restart_data(FILE *gridsave_file) {
  printout("binned radiation field and detailed lines, ");

  fprintf(gridsave_file, "%d\n", 30490824);  // special number marking the beginning of radfield data

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    fprintf(gridsave_file, "%d %la %la %la %la\n", RADFIELDBINCOUNT, nu_lower_first_initial, nu_upper_last_initial,
            T_R_min, T_R_max);

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      fprintf(gridsave_file, "%d %la\n", binindex, radfieldbin_nu_upper[binindex]);
    }
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    const int nbfcontinua = globals::nbfcontinua;
    fprintf(gridsave_file, "%d\n", nbfcontinua);

    for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
      if (grid::get_numassociatedcells(modelgridindex) > 0) {
        const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
        fprintf(gridsave_file, "%d\n", modelgridindex);
        for (int i = 0; i < nbfcontinua; i++) {
          fprintf(gridsave_file, "%a ", prev_bfrate_normed[nonemptymgi * nbfcontinua + i]);
        }
      }
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    fprintf(gridsave_file, "%d\n", detailed_linecount);

    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
      fprintf(gridsave_file, "%d ", detailed_lineindicies[jblueindex]);
    }
  }

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    if (grid::get_numassociatedcells(modelgridindex) > 0) {
      const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
      assert_testmodeonly(nonemptymgi >= 0);
      fprintf(gridsave_file, "%d %la\n", modelgridindex, J_normfactor[nonemptymgi]);

      if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
        for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
          const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
          fprintf(gridsave_file, "%la %la %a %a %d\n", radfieldbins[mgibinindex].J_raw,
                  radfieldbins[mgibinindex].nuJ_raw, radfieldbin_solutions[mgibinindex].W,
                  radfieldbin_solutions[mgibinindex].T_R, radfieldbins[mgibinindex].contribcount);
          // radfieldbins[mgibinindex].fit_type
        }
      }

      if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
        for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
          fprintf(gridsave_file, "%la %d\n", Jb_lu_raw[modelgridindex][jblueindex].value,
                  Jb_lu_raw[modelgridindex][jblueindex].contribcount);
        }
      }
    }
  }
  fprintf(gridsave_file, "%d\n", 42809403);  // special number marking the end of radfield data
}

void read_restart_data(FILE *gridsave_file) {
  printout("Reading restart data for radiation field\n");

  int code_check = 0;
  assert_always(fscanf(gridsave_file, "%d\n", &code_check) == 1);
  if (code_check != 30490824) {
    printout("ERROR: Beginning of radfield restart data not found! Found %d instead of 30490824\n", code_check);
    abort();
  }

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    double T_R_min_in = NAN;
    double T_R_max_in = NAN;
    double nu_lower_first_initial_in = NAN;
    double nu_upper_last_initial_in = NAN;
    int bincount_in = 0;
    assert_always(fscanf(gridsave_file, "%d %la %la %la %la\n", &bincount_in, &nu_lower_first_initial_in,
                         &nu_upper_last_initial_in, &T_R_min_in, &T_R_max_in) == 5);

    double nu_lower_first_ratio = nu_lower_first_initial_in / nu_lower_first_initial;
    if (nu_lower_first_ratio > 1.0) {
      nu_lower_first_ratio = 1 / nu_lower_first_ratio;
    }

    double nu_upper_last_ratio = nu_upper_last_initial_in / nu_upper_last_initial;
    if (nu_upper_last_ratio > 1.0) {
      nu_upper_last_ratio = 1 / nu_upper_last_ratio;
    }

    if (bincount_in != RADFIELDBINCOUNT || T_R_min_in != T_R_min || T_R_max_in != T_R_max ||
        nu_lower_first_ratio < 0.999 || nu_upper_last_ratio < 0.999) {
      printout(
          "ERROR: gridsave file specifies %d bins, nu_lower_first_initial %lg nu_upper_last_initial %lg T_R_min %lg "
          "T_R_max %lg\n",
          bincount_in, nu_lower_first_initial_in, nu_upper_last_initial_in, T_R_min_in, T_R_max_in);
      printout("require %d bins, nu_lower_first_initial %lg nu_upper_last_initial %lg T_R_min %lg T_R_max %lg\n",
               RADFIELDBINCOUNT, nu_lower_first_initial, nu_upper_last_initial, T_R_min, T_R_max);
      abort();
    }

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      int binindex_in = 0;
      assert_always(fscanf(gridsave_file, "%d %la\n", &binindex_in, &radfieldbin_nu_upper[binindex]) == 2);
      assert_always(binindex_in == binindex);
    }
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    int gridsave_nbf_in = 0;
    assert_always(fscanf(gridsave_file, "%d\n", &gridsave_nbf_in) == 1);
    assert_always(gridsave_nbf_in == globals::nbfcontinua);

    for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
      if (grid::get_numassociatedcells(modelgridindex) > 0) {
        const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
        int mgi_in = 0;
        assert_always(fscanf(gridsave_file, "%d\n", &mgi_in) == 1);
        assert_always(mgi_in == modelgridindex);
        for (int i = 0; i < globals::nbfcontinua; i++) {
          float bfrate_normed = 0;
          assert_always(fscanf(gridsave_file, "%a ", &bfrate_normed) == 1);

          const int mgibfindex = nonemptymgi * globals::nbfcontinua + i;
#ifdef MPI_ON
          if (globals::rank_in_node == 0)
#endif
          {
            prev_bfrate_normed[mgibfindex] = bfrate_normed;
          }
        }
      }
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    int detailed_linecount_in = 0;
    assert_always(fscanf(gridsave_file, "%d\n", &detailed_linecount_in) == 1);

    if (detailed_linecount_in != detailed_linecount) {
      printout("ERROR: gridsave file specifies %d detailed lines but this simulation has %d.\n", detailed_linecount_in,
               detailed_linecount);
      abort();
    }

    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
      assert_always(fscanf(gridsave_file, "%d ", &detailed_lineindicies[jblueindex]) == 1);
    }
  }

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    if (grid::get_numassociatedcells(modelgridindex) > 0) {
      const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
      int mgi_in = 0;
      assert_always(fscanf(gridsave_file, "%d %la\n", &mgi_in, &J_normfactor[nonemptymgi]) == 2);
      if (mgi_in != modelgridindex) {
        printout("ERROR: expected data for cell %d but found cell %d\n", modelgridindex, mgi_in);
        abort();
      }

      if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
        for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
          const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
          float W = 0;
          float T_R = 0;
          assert_always(fscanf(gridsave_file, "%la %la %a %a %d\n", &radfieldbins[mgibinindex].J_raw,
                               &radfieldbins[mgibinindex].nuJ_raw, &W, &T_R,
                               &radfieldbins[mgibinindex].contribcount) == 5);
#ifdef MPI_ON
          if (globals::rank_in_node == 0)
#endif
          {
            radfieldbin_solutions[mgibinindex].W = W;
            radfieldbin_solutions[mgibinindex].T_R = T_R;
          }
        }
      }

      if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
        for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
          assert_always(fscanf(gridsave_file, "%la %d\n", &Jb_lu_raw[modelgridindex][jblueindex].value,
                               &Jb_lu_raw[modelgridindex][jblueindex].contribcount) == 2);
        }
      }
    }
  }
  assert_always(fscanf(gridsave_file, "%d\n", &code_check) == 1);
  if (code_check != 42809403) {
    printout("ERROR: End of radfield restart data not found! Found %d instead of 42809403\n", code_check);
    abort();
  }
}

// not in use, but could potential improve speed and accuracy of integrating
// across the binned radiation field which is discontinuous at the bin boundaries
inline auto integrate(const gsl_function *f, double nu_a, double nu_b, double epsabs, double epsrel, size_t limit,
                      int key, gsl_integration_workspace *workspace, double *result, double *abserr) -> int {
  if (MULTIBIN_RADFIELD_MODEL_ON && (globals::timestep >= FIRST_NLTE_RADFIELD_TIMESTEP)) {
    auto *pts = static_cast<double *>(malloc((RADFIELDBINCOUNT + 3) * sizeof(double)));
    int binindex_a = select_bin(nu_a);
    const int binindex_b = select_bin(nu_b);
    int npts = 0;
    pts[npts++] = nu_a;
    if (binindex_a == binindex_b)  // both higher, both lower, or match the same bin
    {
      // region doesn't contain any bins
      pts[npts++] = nu_b;
    } else {
      if (binindex_a < 0)  // a is below the first bin
      {
        binindex_a = 0;
        pts[npts++] = get_bin_nu_lower(0);
      }

      const int maxbinplusone = (binindex_b < 0) ? RADFIELDBINCOUNT : binindex_b;

      for (int binindex = binindex_a; binindex < maxbinplusone; binindex++) {
        pts[npts++] = get_bin_nu_upper(binindex);
      }

      pts[npts++] = nu_b;
    }
    // for (int e = 0; e < npts; e++)
    // {
    //   printout("radfield::integrate singular point number %d at nu %g, (nu_a %g, nu_b %g), low %g high %g\n",
    //            e, pts[e], nu_a, nu_b, radfield(pts[e] * 0.9999, 0), radfield(pts[e] * 1.0001, 0));
    // }
    const int status = gsl_integration_qagp(f, pts, npts, epsabs, epsrel, limit, workspace, result, abserr);
    free(pts);
    return status;
  }
  return gsl_integration_qag(f, nu_a, nu_b, epsabs, epsrel, limit, key, workspace, result, abserr);
}

}  // namespace radfield