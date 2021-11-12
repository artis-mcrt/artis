#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_debye.h>
#include <gsl/gsl_roots.h>
#include "atomic.h"
#include "grid.h"
#include "ltepop.h"
#include "vectors.h"
#include "radfield.h"
#include "rpkt.h"
#include "sn3d.h"

namespace radfield
{

__managed__ static double *J_normfactor = NULL;

__managed__ static bool initialized = false;

// typedef enum
// {
//   FIT_DILUTE_BLACKBODY = 0,
//   FIT_CONSTANT = 1,
// } enum_bin_fit_type;

struct radfieldbin_solution
{
  // these two parameters are used in the current timestep, but were calculated
  // from the values of J and nuJ in the previous timestep
  float W;                // dilution (scaling) factor
  float T_R;              // radiation temperature
  // enum_bin_fit_type fit_type;
};

struct radfieldbin
{
  double J_raw;           // value needs to be multipled by J_normfactor to get the true value
  double nuJ_raw;
  int contribcount;
};

__managed__ static double radfieldbin_nu_upper[RADFIELDBINCOUNT]; // array of upper frequency boundaries of bins
__managed__ static struct radfieldbin *radfieldbins = NULL;
__managed__ static struct radfieldbin_solution *radfieldbin_solutions = NULL;


#ifdef MPI_ON
  MPI_Win win_radfieldbin_solutions = MPI_WIN_NULL;
  MPI_Win win_prev_bfrate_normed = MPI_WIN_NULL;
#endif

// ** Detailed lines - Jblue_lu estimators for selected lines

struct Jb_lu_estimator
{
  double value;
  int contribcount;
};

// reallocate the detailed line arrays in units of BLOCKSIZEJBLUE
static const int BLOCKSIZEJBLUE = 128;
__managed__ static int detailed_linecount = 0;

// array of indicies into the linelist[] array for selected lines
__managed__ static int *detailed_lineindicies;

__managed__ static struct Jb_lu_estimator **prev_Jb_lu_normed = NULL;  // value from the previous timestep
__managed__ static struct Jb_lu_estimator **Jb_lu_raw = NULL;   // unnormalised estimator for the current timestep

// ** end detailed lines

#if (DETAILED_BF_ESTIMATORS_ON)
__managed__ static float *prev_bfrate_normed = NULL;  // values from the previous timestep
__managed__ static double *bfrate_raw = NULL;   // unnormalised estimators for the current timestep

// expensive debugging mode to track the contributions to each bound-free rate estimator
#if (DETAILED_BF_ESTIMATORS_BYTYPE)
  struct bfratecontrib
  {
    int emissiontype;
    double ratecontrib;
  };

  __managed__ static struct bfratecontrib ***bfrate_raw_bytype;   // unnormalised estimator contributions for stats
  __managed__ static int **bfrate_raw_bytype_size;

  static int compare_bfrate_raw_bytype(const void *p1, const void *p2)
  {
    const struct bfratecontrib *elem1 = (struct bfratecontrib *) p1;
    const struct bfratecontrib *elem2 = (struct bfratecontrib *) p2;

   if (elem1->ratecontrib < elem2->ratecontrib)
      return 1;
   else if (elem1->ratecontrib > elem2->ratecontrib)
      return -1;
   else
      return 0;
  }
  #endif
#endif

__managed__ static double *J = NULL; // after normalisation: [ergs/s/sr/cm2/Hz]
#ifdef DO_TITER
  __managed__ static double *J_reduced_save = NULL;
#endif

// J and nuJ are accumulated and then normalised in-place
// i.e. be sure the normalisation has been applied (exactly once) before using the values here!
#ifndef FORCE_LTE
  __managed__ static double *nuJ = NULL;
  #ifdef DO_TITER
    __managed__ static double *nuJ_reduced_save = NULL;
  #endif
#endif


typedef enum
{
  ONE = 0,
  TIMES_NU = 1,
} enum_prefactor;

typedef struct
{
  double T_R;
  enum_prefactor prefactor;
} gsl_planck_integral_paras;

typedef struct
{
  int modelgridindex;
  int binindex;
} gsl_T_R_solver_paras;

static FILE *radfieldfile = NULL;


__host__ __device__ extern inline double dbb(double nu, float T, float W);


static inline
double get_bin_nu_upper(int binindex)
{
  return radfieldbin_nu_upper[binindex];
}


static void
setup_bin_boundaries(void)
{
  // double prev_nu_upper = nu_lower_first_initial;

  // choose between equally spaced in energy/frequency or wavelength (before bf edges shift boundaries around)
  const double delta_nu = (nu_upper_last_initial - nu_lower_first_initial) / (RADFIELDBINCOUNT - 1); // - 1 for the top super bin
  // const double lambda_lower_first_initial = 1e8 * CLIGHT / nu_lower_first_initial;
  // const double lambda_upper_last_initial = 1e8 * CLIGHT / nu_upper_last_initial;
  // const double delta_lambda = (lambda_upper_last_initial - lambda_lower_first_initial) / RADFIELDBINCOUNT;

  for (int binindex = 0; binindex < RADFIELDBINCOUNT - 1; binindex++)
  {
    // radfieldbin_nu_upper[binindex] = 1e8 * CLIGHT / (lambda_lower_first_initial + (binindex + 1) * delta_lambda);
    radfieldbin_nu_upper[binindex] = nu_lower_first_initial + (binindex + 1) * delta_nu;

    // Align the bin edges with bound-free edges, except for the last one
    // if (binindex < RADFIELDBINCOUNT - 1)
    // {
    //   for (int i = 0; i < globals::nbfcontinua_ground; i++)
    //   {
    //     const double nu_edge = phixslist[tid].groundcont[i].nu_edge;
    //     const double eV_edge = H * nu_edge / EV;
    //     const double angstrom_edge = 1e8 * CLIGHT / nu_edge;
    //     const int element = phixslist[tid].groundcont[i].element;
    //     const int ion = phixslist[tid].groundcont[i].ion;
    //     const int level = phixslist[tid].groundcont[i].level;
    //     const int phixstargetindex = phixslist[tid].groundcont[i].phixstargetindex;
    //
    //     const int Z = get_element(element);
    //     const int ion_stage = get_ionstage(element, ion);
    //     const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
    //
    //     //printout("bf edge at %g, nu_lower_first %g, nu_upper_last %g\n",nu_edge,nu_lower_first,nu_upper_last);
    //
    //     // this is compares the highest and lowest bins to the bound-free list, only do it once, i.e. binindex == 0
    //     if (binindex == 0 && ((nu_edge < nu_lower_first_initial) || (nu_edge > nu_upper_last_initial)))
    //     {
    //       printout("Missed bf edge at %12.5e Hz (%6.2f eV, %6.1f A), nu_lower_first %11.5e Hz, nu_upper_last %11.5e Hz, Z=%d ion_stage %d level %d upperionlevel %d\n",
    //                nu_edge, eV_edge, angstrom_edge, nu_lower_first_initial, nu_upper_last_initial, Z, ion_stage, level, upperionlevel);
    //     }
    //
    //     const double bin_nu_upper = get_bin_nu_upper(binindex);
    //     if ((nu_edge > prev_nu_upper) && (nu_edge < bin_nu_upper))
    //     {
    //       printout("Shifting bin %d nu_upper from %12.5e Hz to bf edge at %12.5e Hz (%6.2f eV, %6.1f A) for Z=%d ion_stage %d level %d upperionlevel %d\n",
    //                binindex, bin_nu_upper, nu_edge, eV_edge, angstrom_edge, Z, ion_stage, level, upperionlevel);
    //       radfieldbin_nu_upper[binindex] = nu_edge;
    //     }
    //   }
    // }
    // prev_nu_upper = get_bin_nu_upper(binindex);
  }
  radfieldbin_nu_upper[RADFIELDBINCOUNT - 1] = nu_upper_superbin;  // very top end super bin
}


static void realloc_detailed_lines(const int new_size)
{
  detailed_lineindicies = (int *) realloc(detailed_lineindicies, new_size * sizeof(int));
  if (detailed_lineindicies == NULL)
  {
    printout("ERROR: Not enough memory to reallocate detailed Jblue estimator line list\n");
    abort();
  }

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++)
  {
    if (grid::get_numassociatedcells(modelgridindex) > 0)
    {
      prev_Jb_lu_normed[modelgridindex] = (struct Jb_lu_estimator *) realloc(
        prev_Jb_lu_normed[modelgridindex], new_size * sizeof(struct Jb_lu_estimator));

      Jb_lu_raw[modelgridindex] = (struct Jb_lu_estimator *) realloc(
        Jb_lu_raw[modelgridindex], new_size * sizeof(struct Jb_lu_estimator));

      if (prev_Jb_lu_normed[modelgridindex] == NULL || Jb_lu_raw[modelgridindex] == NULL)
      {
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
  if (detailed_linecount % BLOCKSIZEJBLUE == 0)
  {
    const int new_size = detailed_linecount + BLOCKSIZEJBLUE;
    realloc_detailed_lines(new_size);
  }

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++)
  {
    if (grid::get_numassociatedcells(modelgridindex) > 0)
    {
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


static int compare_integers(const void* a, const void* b)
{
   int int_a = * ((int*) a);
   int int_b = * ((int*) b);

  return (int_a > int_b) - (int_a < int_b);
}


void init(int my_rank, int ndo)
// this should be called only after the atomic data is in memory
{
  // printout("radfield::init()\n");
  if (initialized)
  {
    printout("ERROR: Tried to initialize radfield twice!\n");
    abort();
  }

#ifdef MPI_ON
  const int rank_in_node = globals::rank_in_node;
#endif

  const int nonempty_npts_model = grid::get_nonempty_npts_model();

  J_normfactor = (double *) malloc((grid::get_npts_model() + 1) * sizeof(double));
  J = (double *) malloc((grid::get_npts_model() + 1) * sizeof(double));
  #ifdef DO_TITER
    J_reduced_save = (double *) malloc((grid::get_npts_model() + 1) * sizeof(double));
  #endif

  // J and nuJ are accumulated and then normalised in-place
  // i.e. be sure the normalisation has been applied (exactly once) before using the values here!
  #ifndef FORCE_LTE
    nuJ = (double *) malloc((grid::get_npts_model() + 1) * sizeof(double));
    #ifdef DO_TITER
    nuJ_reduced_save = (double *) malloc((grid::get_npts_model() + 1) * sizeof(double));
    #endif
  #endif

  prev_Jb_lu_normed = (struct Jb_lu_estimator **) malloc((grid::get_npts_model() + 1) * sizeof(struct Jb_lu_estimator *));
  Jb_lu_raw = (struct Jb_lu_estimator **) malloc((grid::get_npts_model() + 1) * sizeof(struct Jb_lu_estimator *));

  detailed_linecount = 0;

  detailed_lineindicies = NULL;
  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++)
  {
    prev_Jb_lu_normed[modelgridindex] = NULL;
    Jb_lu_raw[modelgridindex] = NULL;
  }

  if (DETAILED_LINE_ESTIMATORS_ON)
  {
    for (int i = 0; i < globals::nlines; i++)
    {
      const int element = globals::linelist[i].elementindex;
      const int Z = get_element(element);
      if (Z == 26)
      {
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

        if (lowerlevel <= 15 && A_ul > 0.) // ionstage <= 3 && A_ul > 1e3 &&
          addline = true;

        if (addline)
        {
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
    qsort(detailed_lineindicies, detailed_linecount, sizeof(int), compare_integers);
  }

  printout("There are %d lines with detailed Jblue_lu estimators.\n", detailed_linecount);

  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    printout("The multibin radiation field estimators are being used instead of the whole-spectrum fit from timestep %d onwards.\n", FIRST_NLTE_RADFIELD_TIMESTEP);

    printout("Initialising multibin radiation field with %d bins from (%.2f eV, %6.1f A) to (%.2f eV, %6.1f A)\n",
             RADFIELDBINCOUNT, H * nu_lower_first_initial / EV, 1e8 * CLIGHT / nu_lower_first_initial,
             H * nu_upper_last_initial / EV, 1e8 * CLIGHT / nu_upper_last_initial);
    if (ndo > 0)
    {
      char filename[100];
      sprintf(filename,"radfield_%.4d.out", my_rank);
      assert_always(radfieldfile == NULL);
      radfieldfile = fopen_required(filename, "w");
      fprintf(radfieldfile,"%8s %15s %8s %11s %11s %9s %9s %9s %9s %9s %12s\n",
              "timestep","modelgridindex","bin_num","nu_lower","nu_upper",
              "nuJ","J","J_nu_avg","ncontrib","T_R","W");
      fflush(radfieldfile);
    }

    setup_bin_boundaries();
  }
  else
  {
    printout("The radiation field model is a whole-spectrum fit to a single diluted blackbody.\n");
  }

  const long mem_usage_bins = nonempty_npts_model * RADFIELDBINCOUNT * sizeof(struct radfieldbin);
  radfieldbins = (struct radfieldbin *) malloc(nonempty_npts_model * RADFIELDBINCOUNT * sizeof(struct radfieldbin));

  const long mem_usage_bin_solutions = nonempty_npts_model * RADFIELDBINCOUNT * sizeof(struct radfieldbin_solution);

  #ifdef MPI_ON
  {
    MPI_Aint size = (rank_in_node == 0) ? nonempty_npts_model * RADFIELDBINCOUNT * sizeof(struct radfieldbin_solution) : 0;
    MPI_Win_allocate_shared(size, sizeof(struct radfieldbin_solution), MPI_INFO_NULL, globals::mpi_comm_node, &radfieldbin_solutions, &win_radfieldbin_solutions);
    if (rank_in_node != 0)
    {
      int disp_unit;
      MPI_Win_shared_query(win_radfieldbin_solutions, MPI_PROC_NULL, &size, &disp_unit, &radfieldbin_solutions);
    }
  }
  #else
  {
    radfieldbin_solutions = (struct radfieldbin_solution *) malloc(nonempty_npts_model * RADFIELDBINCOUNT * sizeof(struct radfieldbin_solution));
  }
  #endif

  #if (DETAILED_BF_ESTIMATORS_ON)
  {
    #ifdef MPI_ON
    {
      MPI_Win win_prev_bfrate_normed;
      MPI_Aint size = (rank_in_node == 0) ? nonempty_npts_model * globals::nbfcontinua * sizeof(float) : 0;
      MPI_Win_allocate_shared(size, sizeof(float), MPI_INFO_NULL, globals::mpi_comm_node, &prev_bfrate_normed, &win_prev_bfrate_normed);
      if (rank_in_node != 0)
      {
        int disp_unit;
        MPI_Win_shared_query(win_prev_bfrate_normed, MPI_PROC_NULL, &size, &disp_unit, &prev_bfrate_normed);
      }
    }
    #else
    {
      prev_bfrate_normed = (float *) malloc(nonempty_npts_model * globals::nbfcontinua * sizeof(float));
    }
    #endif
    printout("[info] mem_usage: detailed bf estimators for non-empty cells occupy %.3f MB (shared node memory)\n",
             nonempty_npts_model * globals::nbfcontinua * sizeof(float) / 1024. / 1024.);

    bfrate_raw = (double *) malloc(nonempty_npts_model * globals::nbfcontinua * sizeof(double));

    printout("[info] mem_usage: detailed bf estimator acculumators for non-empty cells occupy %.3f MB\n",
             nonempty_npts_model * globals::nbfcontinua * sizeof(double) / 1024. / 1024.);

    #if (DETAILED_BF_ESTIMATORS_BYTYPE)
    {
      bfrate_raw_bytype = (struct bfratecontrib ***) malloc((grid::get_npts_model() + 1) * sizeof(struct bfratecontrib **));
      bfrate_raw_bytype_size = (int **) malloc((grid::get_npts_model() + 1) * sizeof(int *));
    }
    #endif
  }
  #endif

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++)
  {
    if (grid::get_numassociatedcells(modelgridindex) > 0)
    {
      const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
      #if (DETAILED_BF_ESTIMATORS_ON && DETAILED_BF_ESTIMATORS_BYTYPE)
      {
        bfrate_raw_bytype[modelgridindex] = (struct bfratecontrib **) malloc(globals::nbfcontinua * sizeof(struct bfratecontrib *));
        bfrate_raw_bytype_size[modelgridindex] = (int *) malloc(globals::nbfcontinua * sizeof(int));
        for (int allcontindex = 0; allcontindex < globals::nbfcontinua; allcontindex++)
        {
          bfrate_raw_bytype[modelgridindex][allcontindex] = NULL;
          bfrate_raw_bytype_size[modelgridindex][allcontindex] = 0.;
        }
      }
      #endif

      zero_estimators(modelgridindex);

      if (MULTIBIN_RADFIELD_MODEL_ON)
      {
        #ifdef MPI_ON
        if (rank_in_node == 0)
        #endif
        {
          for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
          {
            const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
            radfieldbin_solutions[mgibinindex].W = -1.;
            radfieldbin_solutions[mgibinindex].T_R = -1.;
          }
        }
      }
    }
  }
  printout("[info] mem_usage: radiation field bin accumulators for non-empty cells occupy %.3f MB\n", mem_usage_bins / 1024. / 1024.);
  printout("[info] mem_usage: radiation field bin solutions for non-empty cells occupy %.3f MB (shared node memory)\n", mem_usage_bin_solutions / 1024. / 1024.);

  initialized = true;
}


/// Initialise estimator arrays which hold the last time steps values (used to damp out
/// fluctuations over timestep iterations if DO_TITER is defined) to -1.
void initialise_prev_titer_photoionestimators(void)
{
  //for (n = 0; n < ngrid; n++)
  for (int n = 0; n < grid::get_npts_model(); n++)
  {
    //double T_e = grid::get_Te(n);
    #ifdef DO_TITER
      J_reduced_save[n] = -1.;
    #endif
    #ifndef FORCE_LTE
      #ifdef DO_TITER
        nuJ_reduced_save[n] = -1.;
        globals::ffheatingestimator_save[n] = -1.;
        globals::colheatingestimator_save[n] = -1.;
      #endif
      for (int element = 0; element < get_nelements(); element++)
      {
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions - 1; ion++)
        {
          //  double ionpot,Alpha_sp,sw_ratio,Gamma;
          //ionpot = epsilon(element,ion+1,0) - epsilon(element,ion,0);
          //Alpha_sp = interpolate_ions_spontrecombcoeff(element,ion,T_e);
          //sw_ratio = stat_weight(element,ion+1,0)/stat_weight(element,ion,0);
          //Gamma = Alpha_sp * sw_ratio / SAHACONST * pow(T_e,1.5) * exp(-ionpot/KB/T_e);
          ////gamma_lte = interpolate_photoioncoeff_below(element,ion,0,T_e) + interpolate_photoioncoeff_above(element,ion,0,T_e);
          ////zeta = interpolate_zeta(element,ion,T_e);
          //gammaestimator[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = Gamma; //gamma_lte/zeta;
          ////globals::corrphotoionrenorm[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = 1.;
          ////photoionestimator[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = Gamma; //gamma_lte/zeta;

          #ifdef DO_TITER
            globals::gammaestimator_save[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = -1.;
            if (!NO_LUT_BFHEATING)
              globals::bfheatingestimator_save[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = -1.;
            /*
            photoionestimator_save[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = -1.;
            stimrecombestimator_save[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = -1.;
            ionfluxestimator_save[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = -1.;
            */
          #endif
        }
      }
    #endif
  }
}


__host__ __device__
int get_Jblueindex(const int lineindex)
// returns -1 if the line does not have a Jblue estimator
{
  // slow linear search
  // for (int i = 0; i < detailed_linecount; i++)
  // {
  //   if (detailed_lineindicies[i] == lineindex)
  //     return i;
  // }

  if (!DETAILED_LINE_ESTIMATORS_ON)
    return -1;

  // use a binary search, assuming the list is sorted

  int low = 0;
  int high = detailed_linecount;
  while (low <= high)
  {
    int mid = low + ((high - low) / 2);
    if (detailed_lineindicies[mid] < lineindex)
    {
      low = mid + 1;
    }
    else if (detailed_lineindicies[mid] > lineindex)
    {
      high = mid - 1;
    }
    else
    {
      return mid;
    }
   }

   // const int element = linelist[lineindex].elementindex;
   // const int ion = linelist[lineindex].ionindex;
   // const int lower = linelist[lineindex].lowerlevelindex;
   // const int upper = linelist[lineindex].upperlevelindex;
   // printout("Could not find lineindex %d among %d items (Z=%02d ionstage %d lower %d upper %d)\n",
   //          lineindex, detailed_linecount, get_element(element), get_ionstage(element, ion), lower, upper);

  return -1;
}


__host__ __device__
double get_Jb_lu(const int modelgridindex, const int jblueindex)
{
  assert_always(jblueindex >= 0);
  assert_always(jblueindex < detailed_linecount);
  return prev_Jb_lu_normed[modelgridindex][jblueindex].value;
}


__host__ __device__
int get_Jb_lu_contribcount(const int modelgridindex, const int jblueindex)
{
  assert_always(jblueindex >= 0);
  assert_always(jblueindex < detailed_linecount);
  return prev_Jb_lu_normed[modelgridindex][jblueindex].contribcount;
}


static double get_bin_J(int modelgridindex, int binindex)
// get the normalised J_nu
{
  if (J_normfactor[modelgridindex] <= 0.0)
  {
    printout("radfield: Fatal error: get_bin_J called before J_normfactor set for modelgridindex %d, = %g",modelgridindex,J_normfactor[modelgridindex]);
    abort();
  }
  else if (modelgridindex >= grid::get_npts_model())
  {
    printout("radfield: Fatal error: get_bin_J called before on modelgridindex %d >= grid::get_npts_model()",modelgridindex);
    abort();
  }
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
  return radfieldbins[mgibinindex].J_raw * J_normfactor[modelgridindex];
}


__host__ __device__
static double get_bin_nuJ(int modelgridindex, int binindex)
{
  if (J_normfactor[modelgridindex] <= 0.0)
  {
    printout("radfield: Fatal error: get_bin_nuJ called before J_normfactor set for modelgridindex %d",modelgridindex);
    abort();
  }
  else if (modelgridindex >= grid::get_npts_model())
  {
    printout("radfield: Fatal error: get_bin_nuJ called before on modelgridindex %d >= grid::get_npts_model()",modelgridindex);
    abort();
  }
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
  return radfieldbins[mgibinindex].nuJ_raw * J_normfactor[modelgridindex];
}


__host__ __device__
static inline
double get_bin_nu_bar(int modelgridindex, int binindex)
// importantly, this is average beween the current and previous timestep
{
  const double nuJ_sum = get_bin_nuJ(modelgridindex, binindex);
  const double J_sum = get_bin_J(modelgridindex, binindex);
  return nuJ_sum / J_sum;
}


__host__ __device__
static inline
double get_bin_nu_lower(int binindex)
{
  if (binindex > 0)
    return radfieldbin_nu_upper[binindex - 1];
  else
    return nu_lower_first_initial;
}


__host__ __device__
static inline
int get_bin_contribcount(int modelgridindex, int binindex)
{
  const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
  return radfieldbins[mgibinindex].contribcount;
}


static inline
float get_bin_W(int modelgridindex, int binindex)
{
  const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
  return radfieldbin_solutions[mgibinindex].W;
}


static inline
float get_bin_T_R(int modelgridindex, int binindex)
{
  const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
  return radfieldbin_solutions[mgibinindex].T_R;
}


__host__ __device__
static inline
int select_bin(double nu)
{
  // linear search one by one until found
  if (nu >= radfieldbin_nu_upper[RADFIELDBINCOUNT - 1])
    return -1; // out of range, nu higher than highest bin
  else if (nu < get_bin_nu_lower(0))
    return -2; // out of range, nu lower than lowest bin
  else
  {
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      if (radfieldbin_nu_upper[binindex] > nu)
      {
        return binindex;
      }
    }

    // binary search for bin with nu_lower <= nu > nu_upper
    // int low = 0;
    // int high = RADFIELDBINCOUNT - 1;
    // while (low <= high)
    // {
    //   int mid = low + ((high - low) / 2);
    //   if (radfieldbin_nu_upper[mid] <= nu)
    //   {
    //     low = mid + 1;
    //   }
    //   else if (get_bin_nu_lower(mid) > nu)
    //   {
    //     high = mid - 1;
    //   }
    //   else
    //   {
    //     return mid;
    //   }
    //  }
    assert_always(false);
    return -3;
  }
}


#ifndef FORCE_LTE
void write_to_file(int modelgridindex, int timestep)
{
  assert_always(MULTIBIN_RADFIELD_MODEL_ON);

# ifdef _OPENMP
# pragma omp critical (out_file)
  {
# endif
    if (!initialized)
    {
      printout("Call to radfield::write_to_file before radfield::init\n");
      abort();
    }

    int totalcontribs = 0;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      totalcontribs += get_bin_contribcount(modelgridindex, binindex);

    for (int binindex = -1 - detailed_linecount; binindex < RADFIELDBINCOUNT; binindex++)
    {
      double nu_lower = 0.0;
      double nu_upper = 0.0;
      double nuJ_out = 0.0;
      double J_out = 0.0;
      float T_R = 0.0;
      float W = 0.0;
      double J_nu_bar = 0.0;
      int contribcount = 0;

      bool skipoutput = false;

      if (binindex >= 0)
      {
        nu_lower = get_bin_nu_lower(binindex);
        nu_upper = get_bin_nu_upper(binindex);
        nuJ_out = get_bin_nuJ(modelgridindex, binindex);
        J_out = get_bin_J(modelgridindex, binindex);
        T_R = get_bin_T_R(modelgridindex, binindex);
        W = get_bin_W(modelgridindex, binindex);
        J_nu_bar = J_out / (nu_upper - nu_lower);
        contribcount = get_bin_contribcount(modelgridindex, binindex);
      }
      else if (binindex == -1)
      {
        nuJ_out = nuJ[modelgridindex];
        J_out = J[modelgridindex];
        T_R = grid::get_TR(modelgridindex);
        W = grid::get_W(modelgridindex);
        contribcount = totalcontribs;
      }
      else // use binindex < -1 for detailed line Jb_lu estimators
      {
        const int jblueindex = -2 - binindex; // -2 is the first detailed line, -3 is the second, etc
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

      if (!skipoutput)
      {
        fprintf(radfieldfile, "%d %d %d %11.5e %11.5e %9.3e %9.3e %9.3e %d %9.1f %12.5e\n",
                timestep, modelgridindex, binindex, nu_lower, nu_upper, nuJ_out,
                J_out, J_nu_bar, contribcount, T_R, W);
      }
    }
    fflush(radfieldfile);
# ifdef _OPENMP
  }
# endif
}
#endif


void close_file(void)
{
  if (radfieldfile != NULL)
  {
    fclose(radfieldfile);
    radfieldfile = NULL;
  }

  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    free(radfieldbins);
    #ifdef MPI_ON
    if (win_radfieldbin_solutions != MPI_WIN_NULL)
    {
      MPI_Win_free(&win_radfieldbin_solutions);
    }
    #else
    if (radfieldbin_solutions != NULL)
    {
      free(radfieldbin_solutions);
    }
    #endif
  }

  #if (DETAILED_BF_ESTIMATORS_ON)
    free(bfrate_raw);
    #ifdef MPI_ON
    if (win_radfieldbin_solutions != MPI_WIN_NULL)
    {
      MPI_Win_free(&win_prev_bfrate_normed);
    }
    #else
    if (prev_bfrate_normed != NULL)
    {
      free(prev_bfrate_normed);
    }
    #endif
  #endif
}


void zero_estimators(int modelgridindex)
// set up the new bins and clear the estimators in preparation
// for a timestep
{
  #if (DETAILED_BF_ESTIMATORS_ON)
  assert_always(bfrate_raw != NULL);
  if (initialized && (grid::get_numassociatedcells(modelgridindex) > 0))
  {
    const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
    for (int i = 0; i < globals::nbfcontinua; i++)
    {
      bfrate_raw[nonemptymgi * globals::nbfcontinua + i] = 0.;
    }
  }
  #endif

  if (DETAILED_LINE_ESTIMATORS_ON)
  {
    assert_always(Jb_lu_raw != NULL);
    assert_always(Jb_lu_raw[modelgridindex] != NULL);
    for (int i = 0; i < detailed_linecount; i++)
    {
      Jb_lu_raw[modelgridindex][i].value = 0.;
      Jb_lu_raw[modelgridindex][i].contribcount = 0.;
    }
  }

  assert_always(J != NULL);
  J[modelgridindex] = 0.; // this is required even if FORCE_LTE is on
#ifndef FORCE_LTE
  assert_always(nuJ != NULL);
  nuJ[modelgridindex] = 0.;

  if (MULTIBIN_RADFIELD_MODEL_ON && initialized && (grid::get_numassociatedcells(modelgridindex) > 0))
  {
    // printout("radfield: zeroing estimators in %d bins in cell %d\n",RADFIELDBINCOUNT,modelgridindex);

    assert_always(radfieldbins != NULL);
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
      radfieldbins[mgibinindex].J_raw = 0.0;
      radfieldbins[mgibinindex].nuJ_raw = 0.0;
      radfieldbins[mgibinindex].contribcount = 0;
    }
  }
  set_J_normfactor(modelgridindex, -1.0);
#endif
}


#if (DETAILED_BF_ESTIMATORS_ON)
__host__ __device__
static void increment_bfestimators(
  const int modelgridindex, const double distance_e_cmf, const double nu_cmf,
  const PKT *const pkt_ptr, const double t_current)
{
  assert_always(pkt_ptr->prop_time == t_current);
  if (distance_e_cmf == 0)
    return;

  const int nbfcontinua = globals::nbfcontinua;
  const double dopplerfactor = doppler_packetpos(pkt_ptr);
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  // const double dopplerfactor = 1.;

  // const double deltaV = grid::vol_init_modelcell(modelgridindex) * pow(globals::time_step[nts_global].mid / globals::tmin, 3);
  // const double deltat = globals::time_step[nts_global].width;
  // const double estimator_normfactor_over_H = 1 / deltaV / deltat / nprocs / H;

  const int tid = get_thread_num();
  const double distance_e_cmf_over_nu = distance_e_cmf / nu_cmf * dopplerfactor;
  for (int allcontindex = 0; allcontindex < nbfcontinua; allcontindex++)
  {
    const double nu_edge = globals::allcont_nu_edge[allcontindex];
    const double nu_max_phixs = nu_edge * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table

    if (nu_cmf >= nu_edge && nu_cmf <= nu_max_phixs)
    {
      safeadd(bfrate_raw[nonemptymgi * nbfcontinua + allcontindex], globals::phixslist[tid].gamma_contr[allcontindex] * distance_e_cmf_over_nu);

      #if (DETAILED_BF_ESTIMATORS_BYTYPE)
      const int element = allcont[allcontindex].element;
      // const int ion = allcont[allcontindex].ion;
      // const int ionstage = get_ionstage(element, ion);
      const int atomic_number = get_element(element);
      if ((atomic_number == 26)) //  && ionstage == 2
      {
        const int oldlistsize = bfrate_raw_bytype_size[modelgridindex][allcontindex];
        if (oldlistsize % 16 == 0)
        {
          bfrate_raw_bytype[modelgridindex][allcontindex] = (struct bfratecontrib *) realloc(
            bfrate_raw_bytype[modelgridindex][allcontindex], (oldlistsize + 16) * sizeof(struct bfratecontrib));

          assert_always(bfrate_raw_bytype[modelgridindex][allcontindex] != NULL);
        }

        int listindex = oldlistsize;
        for (int i = 0; i < oldlistsize; i++)
        {
          if (bfrate_raw_bytype[modelgridindex][allcontindex][i].emissiontype == pkt_ptr->emissiontype)
          {
            listindex = i;
            break;
          }
        }

        if (listindex == oldlistsize)
        {
          // printout("  bfrate_bytype allcontindex %d expanding list size to %d for emtype %d\n", allcontindex, oldlistsize + 1, pkt_ptr->emissiontype);
          bfrate_raw_bytype_size[modelgridindex][allcontindex] = oldlistsize + 1;
          bfrate_raw_bytype[modelgridindex][allcontindex][listindex].emissiontype = pkt_ptr->emissiontype;
          bfrate_raw_bytype[modelgridindex][allcontindex][listindex].ratecontrib = 0;
        }

        bfrate_raw_bytype[modelgridindex][allcontindex][listindex].ratecontrib += (globals::phixslist[tid].gamma_contr[allcontindex] * distance_e_cmf_over_nu);
      }
      #endif
    }
    else if (nu_cmf < nu_edge)
    {
      // list is sorted by nu_edge, so all remaining will have nu_cmf < nu_edge
      break;
    }
  }
}
#endif


__host__ __device__
void update_estimators(int modelgridindex, double distance_e_cmf, double nu_cmf, const PKT *const pkt_ptr, double t_current)
{
  assert_always(pkt_ptr->prop_time == t_current);
  safeadd(J[modelgridindex], distance_e_cmf);
  #ifdef DEBUG_ON
    if (!std::isfinite(J[modelgridindex]))
    {
      printout("[fatal] update_estimators: estimator becomes non finite: distance_e_cmf %g, nu_cmf %g ... abort\n",distance_e_cmf,nu_cmf);
      abort();
    }
  #endif

#ifndef FORCE_LTE
  safeadd(nuJ[modelgridindex], distance_e_cmf * nu_cmf);
  #ifdef DEBUG_ON
    if (!std::isfinite(nuJ[modelgridindex]))
    {
      printout("[fatal] update_estimators: estimator becomes non finite: distance_e_cmf %g, nu_cmf %g ... abort\n",distance_e_cmf,nu_cmf);
      abort();
    }
  #endif

  #if (DETAILED_BF_ESTIMATORS_ON)
  increment_bfestimators(modelgridindex, distance_e_cmf, nu_cmf, pkt_ptr, t_current);
  #endif

  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    // int binindex = 0;
    // if (nu_cmf <= get_bin_nu_lower(modelgridindex,binindex))
    // {
    //   #ifdef DEBUG_ON
    //   printout("radfield: Extending nu_lower_first from %g down to %g\n",nu_lower_first,nu_cmf);
    //   #endif
    //   nu_lower_first = nu_cmf;
    // }
    // else if (nu_cmf > radfieldbin_nu_upper[modelgridindex][RADFIELDBINCOUNT - 1])
    // {
    //   binindex = RADFIELDBINCOUNT - 1;
    //   #ifdef DEBUG_ON
    //   printout("radfield: Extending nu_upper_last from %g up to %g\n",get_bin_nu_upper(modelgridindex,binindex),nu_cmf);
    //   #endif
    //   get_bin_nu_upper(modelgridindex, binindex) = nu_cmf;
    // }
    // else
    // {
    //   binindex = select_bin(modelgridindex,nu_cmf);
    // }
    const int binindex = select_bin(nu_cmf);

    if (binindex >= 0)
    {
      const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
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
#endif
}


void increment_lineestimator(const int modelgridindex, const int lineindex, const double increment)
{
  if (!DETAILED_LINE_ESTIMATORS_ON) return;

  const int jblueindex = get_Jblueindex(lineindex);
  if (jblueindex >= 0)
  {
    Jb_lu_raw[modelgridindex][jblueindex].value += increment;
    Jb_lu_raw[modelgridindex][jblueindex].contribcount += 1;
    // const int lineindex = detailed_lineindicies[jblueindex];
    // printout(" increment cell %d lineindex %d Jb_lu_raw %g prev_Jb_lu_normed %g radfield(nu_trans) %g\n",
    //       modelgridindex, lineindex, Jb_lu_raw[modelgridindex][jblueindex], prev_Jb_lu_normed[modelgridindex][jblueindex].value, radfield(linelist[lineindex].nu, modelgridindex));
  }
}


__host__ __device__
double dbb_mgi(double nu, int modelgridindex)
{
  const float T_R_fullspec = grid::get_TR(modelgridindex);
  const float W_fullspec   = grid::get_W(modelgridindex);
  return dbb(nu, T_R_fullspec, W_fullspec);
}


__host__ __device__
double radfield(double nu, int modelgridindex)
// returns mean intensity J_nu [ergs/s/sr/cm2/Hz]
{
  if (MULTIBIN_RADFIELD_MODEL_ON && (globals::nts_global >= FIRST_NLTE_RADFIELD_TIMESTEP))
  {
    // const double lambda = 1e8 * CLIGHT / nu;
    // if (lambda < 1085) // Fe II ground state edge
    // {
    //   return dbb(nu, grid::get_TR(modelgridindex), grid::get_W(modelgridindex));
    // }
    const int binindex = select_bin(nu);
    if (binindex >= 0)
    {
      const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
      const struct radfieldbin_solution *const bin = &radfieldbin_solutions[mgibinindex];
      if (bin->W >= 0.)
      {
        // if (bin->fit_type == FIT_DILUTE_BLACKBODY)
        {
          const double J_nu = dbb(nu, bin->T_R, bin->W);
          return J_nu;
        }
        // else
        // {
        //   return bin->W;
        // }
      }
      else //W < 0
      {
        //printout("WARNING: Radfield modelgridindex %d binindex %d has W_bin=%g<0, using W %g T_R %g nu %g\n",
        //         modelgridindex, binindex, W_bin, W_fullspec, T_R_fullspec, nu);
      }
    }
    else //binindex < 0
    {
      // if (nu > get_bin_nu_upper(RADFIELDBINCOUNT - 1))
      // {
      //   // undiluted LTE blueward of the bins
      //   const double J_nu_LTE = dbb(nu, grid::get_Te(modelgridindex), 1.0);
      //   return J_nu_LTE;
      // }
      // else
      //   return 0; // no radfield redwards of the bins
      //printout("WARNING: Radfield modelgridindex %d binindex %d nu %g nu_lower_first %g nu_upper_last %g \n",
      //         modelgridindex, binindex, nu, nu_lower_first, nu_upper_last);
    }
    return 0.;
  }
  /*else
  {
    printout("radfield: WARNING: Radfield called before initialized. Using global T_R %g W %g nu %g modelgridindex %d\n",
             W_fullspec, T_R_fullspec, nu, modelgridindex);
  }*/

  const float T_R_fullspec = grid::get_TR(modelgridindex);
  const float W_fullspec   = grid::get_W(modelgridindex);
  const double J_nu_fullspec = dbb(nu, T_R_fullspec, W_fullspec);
  return J_nu_fullspec;
}


static double gsl_integrand_planck(double nu, void *paras)
{
  const double T_R = ((gsl_planck_integral_paras *) paras)->T_R;
  const enum_prefactor prefactor = ((gsl_planck_integral_paras *) paras)->prefactor;

  double integrand = TWOHOVERCLIGHTSQUARED * pow(nu,3) / (expm1(HOVERKB * nu / T_R));

  if (prefactor == TIMES_NU)
    integrand *= nu;

  return integrand;
}


#ifndef __CUDA_ARCH__
static double planck_integral(double T_R, double nu_lower, double nu_upper, enum_prefactor prefactor)
{
  double integral = 0.;

  double error = 0.;
  const double epsrel = 1e-10;
  const double epsabs = 0.;

  gsl_planck_integral_paras intparas;
  intparas.T_R = T_R;
  intparas.prefactor = prefactor;

  gsl_function F_planck;
  F_planck.function = &gsl_integrand_planck;
  F_planck.params = &intparas;

  gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);
  int status = gsl_integration_qag(&F_planck, nu_lower, nu_upper, epsabs, epsrel, GSLWSIZE, GSL_INTEG_GAUSS61, gslworkspace, &integral, &error);
  if (status != 0)
  {
    printout("planck_integral integrator status %d, GSL_FAILURE= %d. Integral value %g, setting to zero.\n", status,GSL_FAILURE,integral);
    integral = 0.;
  }
  gsl_set_error_handler(previous_handler);

  return integral;
}
#endif


#ifndef __CUDA_ARCH__
static double planck_integral_analytic(double T_R, double nu_lower, double nu_upper, enum_prefactor prefactor)
{
  double integral = 0.;

  if (prefactor == TIMES_NU)
  {
    double debye_upper = gsl_sf_debye_4(HOVERKB * nu_upper / T_R) * pow(nu_upper,4);
    double debye_lower = gsl_sf_debye_4(HOVERKB * nu_lower / T_R) * pow(nu_lower,4);
    integral = TWOHOVERCLIGHTSQUARED * (debye_upper - debye_lower) * T_R / HOVERKB / 4.;
  }
  else
  {
    double debye_upper = gsl_sf_debye_3(HOVERKB * nu_upper / T_R) * pow(nu_upper,3);
    double debye_lower = gsl_sf_debye_3(HOVERKB * nu_lower / T_R) * pow(nu_lower,3);
    integral = TWOHOVERCLIGHTSQUARED * (debye_upper - debye_lower) * T_R / HOVERKB / 3.;

    if (integral == 0.)
    {
      /*double upperexp = exp(HOVERKB * nu_upper / T_R);
      double upperint = - pow(nu_upper,4) / 4
                        + pow(nu_upper,3) * log(1 - upperexp) / HOVERKB
                        + 3 * pow(nu_upper,2) * polylog(2,upperexp) / pow(HOVERKB,2)
                        - 6 * nu_upper * polylog(3,upperexp) / pow(HOVERKB,3)
                        + 6 * polylog(4,upperexp) / pow(HOVERKB,4);
      double lowerexp = exp(HOVERKB * nu_lower / T_R);
      double lowerint = - pow(nu_lower,4) / 4
                        + pow(nu_lower,3) * log(1 - lowerexp) / HOVERKB
                        + 3 * pow(nu_lower,2) * polylog(2,lowerexp) / pow(HOVERKB,2)
                        - 6 * nu_lower * polylog(3,lowerexp) / pow(HOVERKB,3)
                        + 6 * polylog(4,lowerexp) / pow(HOVERKB,4);
      double integral2 = TWOHOVERCLIGHTSQUARED * (upperint - lowerint);

      printout("planck_integral_analytic is zero. debye_upper %g debye_lower %g. Test alternative %g\n",
               debye_upper,debye_lower,integral2);*/
    }
  }

  return integral;
}
#endif


#ifndef __CUDA_ARCH__
static double delta_nu_bar(double T_R, void *paras)
// difference between the average nu and the average nu of a planck function
// at temperature T_R, in the frequency range corresponding to a bin
{
  const int modelgridindex = ((gsl_T_R_solver_paras *) paras)->modelgridindex;
  const int binindex = ((gsl_T_R_solver_paras *) paras)->binindex;

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

    printout("planck_integral_analytic is %g. Replacing with numerical result of %g.\n",nu_bar_planck,nu_bar_planck_numerical);
    nu_bar_planck = nu_bar_planck_numerical;
  }*/

  const double delta_nu_bar = nu_bar_planck_T_R - nu_bar_estimator;

  if (!std::isfinite(delta_nu_bar))
  {
    printout("delta_nu_bar is %g. nu_bar_planck_T_R %g nu_times_planck_numerical %g planck_integral_numerical %g nu_bar_estimator %g\n",
             delta_nu_bar, nu_bar_planck_T_R, nu_times_planck_numerical, planck_integral_numerical, nu_bar_estimator);
  }

  //double delta_nu_bar = nu_bar_planck_T_R / nu_bar_estimator - 1.0;

  //printout("delta_nu_bar %g nu_bar_planck %g\n",delta_nu_bar,nu_bar_planck);

  return delta_nu_bar;
}
#endif


#ifndef __CUDA_ARCH__
static float find_T_R(int modelgridindex, int binindex)
{
  double T_R = 0.0;

  gsl_T_R_solver_paras paras;
  paras.modelgridindex = modelgridindex;
  paras.binindex = binindex;

  /// Check whether the equation has a root in [T_min,T_max]
  double delta_nu_bar_min = delta_nu_bar(T_R_min,&paras);
  double delta_nu_bar_max = delta_nu_bar(T_R_max,&paras);

  // printout("find_T_R: bin %4d delta_nu_bar(T_R_min) %g, delta_nu_bar(T_R_max) %g\n",
  //          binindex, delta_nu_bar_min,delta_nu_bar_max);

  if (!std::isfinite(delta_nu_bar_min) || !std::isfinite(delta_nu_bar_max))
    delta_nu_bar_max = delta_nu_bar_min = -1;

  if (delta_nu_bar_min * delta_nu_bar_max < 0)
  {
    /// If there is a root in the interval, solve for T_R

    const double epsrel = 1e-4;
    const double epsabs = 0.;
    const int maxit = 100;

    gsl_function find_T_R_f;
    find_T_R_f.function = &delta_nu_bar;
    find_T_R_f.params = &paras;

    ///one dimensional gsl root solver, bracketing type
    gsl_root_fsolver *T_R_solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(T_R_solver, &find_T_R_f, T_R_min, T_R_max);
    int iteration_num = 0;
    int status;
    do
    {
      iteration_num++;
      gsl_root_fsolver_iterate(T_R_solver);
      T_R = gsl_root_fsolver_root(T_R_solver);

      const double T_R_lower = gsl_root_fsolver_x_lower(T_R_solver);
      const double T_R_upper = gsl_root_fsolver_x_upper(T_R_solver);
      status = gsl_root_test_interval(T_R_lower,T_R_upper,epsabs,epsrel);

      //printout("find_T_R: bin %4d iter %d, T_R is between %7.1f and %7.1f, guess %7.1f, delta_nu_bar %g, status %d\n",
      //         binindex,iteration_num,T_R_lower,T_R_upper,T_R,delta_nu_bar(T_R,&paras),status);
    }
    while (status == GSL_CONTINUE && iteration_num < maxit);

    if (status == GSL_CONTINUE)
      printout("[warning] find_T_R: T_R did not converge within %d iterations\n", maxit);

    gsl_root_fsolver_free(T_R_solver);
  }
  else if (delta_nu_bar_max < 0)
  {
    /// Thermal balance equation always negative ===> T_R = T_min
    /// Calculate the rates again at this T_e to print them to file
    T_R = T_R_max;
    printout("find_T_R: cell %d bin %4d no solution in interval, clamping to T_R_max=%g\n",
             modelgridindex, binindex, T_R_max);
  }
  else
  {
    T_R = T_R_min;
    printout("find_T_R: cell %d bin %4d no solution in interval, clamping to T_R_min=%g\n",
             modelgridindex, binindex, T_R_min);
  }

  return T_R;
}
#endif


#if (!defined __CUDA_ARCH__ && !defined FORCE_LTE)
static void set_params_fullspec(const int modelgridindex, const int timestep)
{
  const double nubar = nuJ[modelgridindex] / J[modelgridindex];
  if (!std::isfinite(nubar) || nubar == 0.)
  {
    printout("[warning] T_R estimator infinite in cell %d, keep T_R, T_J, W of last timestep. J = %g. nuJ = %g\n",
             modelgridindex, J[modelgridindex], nuJ[modelgridindex]);
  }
  else
  {
    float T_J = pow(J[modelgridindex] * PI / STEBO, 1 / 4.);
    if (T_J > MAXTEMP)
    {
      printout("[warning] temperature estimator T_J = %g exceeds T_max %g in cell %d. Setting T_J = T_max!\n", T_J, MAXTEMP, modelgridindex);
      T_J = MAXTEMP;
    }
    else if (T_J < MINTEMP)
    {
      printout("[warning] temperature estimator T_J = %g below T_min %g in cell %d. Setting T_J = T_min!\n", T_J, MINTEMP, modelgridindex);
      T_J = MINTEMP;
    }
    grid::set_TJ(modelgridindex, T_J);

    float T_R = H * nubar / KB / 3.832229494;
    if (T_R > MAXTEMP)
    {
      printout("[warning] temperature estimator T_R = %g exceeds T_max %g in cell %d. Setting T_R = T_max!\n", T_R, MAXTEMP, modelgridindex);
      T_R = MAXTEMP;
    }
    else if (T_R < MINTEMP)
    {
      printout("[warning] temperature estimator T_R = %g below T_min %g in cell %d. Setting T_R = T_min!\n", T_R, MINTEMP, modelgridindex);
      T_R = MINTEMP;
    }
    grid::set_TR(modelgridindex, T_R);

    const float W = J[modelgridindex] * PI / STEBO / pow(T_R, 4);
    grid::set_W(modelgridindex, W);

    printout("Full-spectrum fit radfield for cell %d at timestep %d: J %g, nubar %5.1f Angstrom, T_J %g, T_R %g, W %g\n",
             modelgridindex, timestep, J[modelgridindex], 1e8 * CLIGHT / nubar,
             T_J, T_R, W);
  }
}
#endif


#if (!defined __CUDA_ARCH__ && !defined FORCE_LTE)
void fit_parameters(int modelgridindex, int timestep)
// finds the best fitting W and temperature parameters in each spectral bin
// using J and nuJ
{
  set_params_fullspec(modelgridindex, timestep);

  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    if (J_normfactor[modelgridindex] <= 0)
    {
      printout("radfield: FATAL J_normfactor = %g in cell %d at call to fit_parameters", J_normfactor[modelgridindex], modelgridindex);
      abort();
    }

    double J_bin_sum = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      J_bin_sum += get_bin_J(modelgridindex, binindex);

    printout("radfield bins sum to J of %g (%.1f%% of total J).\n",
             J_bin_sum, 100. * J_bin_sum / J[modelgridindex]);
    printout("radfield: Finding parameters for %d bins...\n", RADFIELDBINCOUNT);

    double J_bin_max = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      const double J_bin = get_bin_J(modelgridindex, binindex);
      if (J_bin > J_bin_max)
        J_bin_max = J_bin;
    }

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      const double nu_lower = get_bin_nu_lower(binindex);
      const double nu_upper = get_bin_nu_upper(binindex);
      const double J_bin = get_bin_J(modelgridindex, binindex);
      float T_R_bin = -1.0;
      double W_bin = -1.0;
      const int contribcount = get_bin_contribcount(modelgridindex, binindex);

      if (contribcount > 0)
      {
        // // enum_bin_fit_type bin_fit_type = radfieldbin_solutions[modelgridindex][binindex].fit_type;
        // if (bin_fit_type == FIT_DILUTE_BLACKBODY)
        {
          T_R_bin = find_T_R(modelgridindex, binindex);

          if (binindex == RADFIELDBINCOUNT - 1)
          {
            const float T_e = grid::get_Te(modelgridindex);
            printout("    replacing bin %d T_R %7.1f with cell T_e = %7.1f\n",
                     binindex, get_bin_T_R(modelgridindex, binindex), T_e);
            T_R_bin = T_e;
          }

          double planck_integral_result = planck_integral(T_R_bin, nu_lower, nu_upper, ONE);
//          printout("planck_integral(T_R=%g, nu_lower=%g, nu_upper=%g) = %g\n", T_R_bin, nu_lower, nu_upper, planck_integral_result);

          W_bin = J_bin / planck_integral_result;

          if (W_bin > 1e4)
          {
//            printout("T_R_bin %g, nu_lower %g, nu_upper %g\n", T_R_bin, nu_lower, nu_upper);
            printout("W %g too high, trying setting T_R of bin %d to %g. J_bin %g planck_integral %g\n",
                     W_bin, binindex, T_R_max, J_bin, planck_integral_result);
            planck_integral_result = planck_integral(T_R_max, nu_lower, nu_upper, ONE);
            W_bin = J_bin / planck_integral_result;
            if (W_bin > 1e4)
            {
              printout("W still very high, W=%g. Zeroing bin...\n", W_bin);
              T_R_bin = -99.0;
              W_bin = 0.;
            }
            else
            {
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
      }
      else
      {
        T_R_bin = 0.;
        W_bin = 0.;
      }
      // else
      // {
      //   T_R_bin = -1;
      //   W_bin = -1;
      // }
      const int mgibinindex = grid::get_modelcell_nonemptymgi(modelgridindex) * RADFIELDBINCOUNT + binindex;
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
    //   printout("bin %4d (lambda %7.1f Å to %7.1f Å): contribcount %5d J %7.1e T_R %8.1f W %12.5e lambdabar %7.1f Å\n",
    //          binindex, 1e8 * CLIGHT / prev_nu_upper, 1e8 * CLIGHT / bin_nu_upper, contribcount, J_bin, T_R_bin, W_bin, 1e8 * CLIGHT / nubar);
    //
    //  prev_nu_upper = get_bin_nu_upper(binindex);
    // }

    write_to_file(modelgridindex, timestep);
  }
}
#endif


__host__ __device__
void set_J_normfactor(int modelgridindex, double normfactor)
{
  J_normfactor[modelgridindex] = normfactor;
}


__host__ __device__
void normalise_J(const int modelgridindex, const double estimator_normfactor_over4pi)
{
  assert_always(std::isfinite(J[modelgridindex]));
  J[modelgridindex] *= estimator_normfactor_over4pi;
  for (int i = 0; i < detailed_linecount; i++)
  {
    prev_Jb_lu_normed[modelgridindex][i].value = Jb_lu_raw[modelgridindex][i].value * estimator_normfactor_over4pi;
    prev_Jb_lu_normed[modelgridindex][i].contribcount = Jb_lu_raw[modelgridindex][i].contribcount;
  }
}


__host__ __device__
void normalise_bf_estimators(const int modelgridindex, const double estimator_normfactor_over_H)
{
  #if (DETAILED_BF_ESTIMATORS_ON)
  printout("normalise_bf_estimators for cell %d with factor %g\n", modelgridindex, estimator_normfactor_over_H);
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  assert_always(nonemptymgi >= 0);
  for (int i = 0; i < globals::nbfcontinua; i++)
  {
    const int mgibfindex = nonemptymgi * globals::nbfcontinua + i;
    prev_bfrate_normed[mgibfindex] = bfrate_raw[mgibfindex] * estimator_normfactor_over_H;

    #if (DETAILED_BF_ESTIMATORS_BYTYPE)
    const int listsize = bfrate_raw_bytype_size[modelgridindex][i];

    for (int j = 0; j < listsize; j++)
    {
      bfrate_raw_bytype[modelgridindex][i][j].ratecontrib = bfrate_raw_bytype[modelgridindex][i][j].ratecontrib * estimator_normfactor_over_H;
    }

    qsort(bfrate_raw_bytype[modelgridindex][i], listsize, sizeof(struct bfratecontrib), compare_bfrate_raw_bytype);
    #endif
  }
  #endif
}


// TODO: use in get_bfrate_estimator
static int get_bfcontindex(const int element, const int lowerion, const int lower, const int phixstargetindex)
{
  const double nu_edge = get_phixs_threshold(element, lowerion, lower, phixstargetindex);
  for (int i = 0; i < globals::nbfcontinua; i++)
  {
    if ((globals::allcont[i].element == element) && (globals::allcont[i].ion == lowerion) &&
        (globals::allcont[i].level == lower) && (globals::allcont[i].phixstargetindex == phixstargetindex))
    {
      return i;
    }

    if (nu_edge > globals::allcont[i].nu_edge)
      break;
  }
  return -1;
}

#if (DETAILED_BF_ESTIMATORS_BYTYPE)
void reset_bfrate_contributions(const int modelgridindex)
{
  for (int allcontindex = 0; allcontindex < globals::nbfcontinua; allcontindex++)
  {
    free(bfrate_raw_bytype[modelgridindex][allcontindex]);
    bfrate_raw_bytype[modelgridindex][allcontindex] = NULL;
    bfrate_raw_bytype_size[modelgridindex][allcontindex] = 0;
  }
}


void print_bfrate_contributions(const int element, const int lowerion, const int lower, const int phixstargetindex, const int modelgridindex, const double nnlowerlevel, const double nnlowerion)
{
  const int allcontindex = get_bfcontindex(element, lowerion, lower, phixstargetindex);
  if (allcontindex >= 0)
  {
    const int listsize = bfrate_raw_bytype_size[modelgridindex][allcontindex];
    if (listsize > 0)
    {
      printout("  %d contributions found for this bf transition\n", listsize);
    }
    else
    {
      // printout("  no contributions found for this bf transition\n");
    }
    for (int i = 0; i < listsize; i++)
    {
      const int et = bfrate_raw_bytype[modelgridindex][allcontindex][i].emissiontype;
      const double bfcontrib = bfrate_raw_bytype[modelgridindex][allcontindex][i].ratecontrib;
      printout("    Gamma_contrib %7.2e gamma_contrib %7.2e emissiontype ", bfcontrib * nnlowerlevel / nnlowerion, bfcontrib);

      if (et >= 0)
      {
        /// bb-emission
        const int element = globals::linelist[et].elementindex;
        const int ion = globals::linelist[et].ionindex;
        const int upper = globals::linelist[et].upperlevelindex;
        const int lower = globals::linelist[et].lowerlevelindex;
        const double lambda_trans = 1e8 * CLIGHT / globals::linelist[et].nu;
        printout("%7d bound-bound Z=%2d ion_stage %d upper+1 %4d lower+1 %4d lambda %5.1f\n",
                 et, get_element(element), get_ionstage(element, ion), upper + 1, lower + 1, lambda_trans);
      }
      else if (et == -9999999)
      {
        /// ff-emission
        printout("%7d free-free scattering\n", et);
      }
      else
      {
        /// bf-emission
        const int bfindex = -1 - et;
        const int element = globals::bflist[bfindex].elementindex;
        const int ion = globals::bflist[bfindex].ionindex;
        const int lower = globals::bflist[bfindex].levelindex;
        const int phixstargetindex = globals::bflist[bfindex].phixstargetindex;

        const double nuthreshold = get_phixs_threshold(element, ion, lower, phixstargetindex) / H;
        const double lambda_trans = 1e8 * CLIGHT / nuthreshold;
        const int upperionlevel = get_phixsupperlevel(element, ion, lower, phixstargetindex);
        assert_always(get_continuumindex(element, ion, lower, upperionlevel) == et);
        const int lowerionstage = get_ionstage(element, ion);
        printout("%7d bound-free  Z=%2d ion_stage %d->%d upper+1 %4d lower+1 %4d lambda %5.1f\n", et, get_element(element), lowerionstage + 1, lowerionstage, upperionlevel + 1, lower + 1, lambda_trans);
      }
    }

    // reset the list for the next timestep
    free(bfrate_raw_bytype[modelgridindex][allcontindex]);
    bfrate_raw_bytype[modelgridindex][allcontindex] = NULL;
    bfrate_raw_bytype_size[modelgridindex][allcontindex] = 0;
  }
  else
  {
    printout("  no continuum index found for this bf transition\n");
    abort();
  }
}
#endif


double get_bfrate_estimator(const int element, const int lowerion, const int lower, const int phixstargetindex, const int modelgridindex)
{
#if (!DETAILED_BF_ESTIMATORS_ON)
  return -1;
#else
  const int nbfcontinua = globals::nbfcontinua;
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);

  // TODO: speed this up with a binary search
  const double nu_edge = get_phixs_threshold(element, lowerion, lower, phixstargetindex);
  for (int i = 0; i < globals::nbfcontinua; i++)
  {
    if ((globals::allcont[i].element == element) && (globals::allcont[i].ion == lowerion) &&
        (globals::allcont[i].level == lower) && (globals::allcont[i].phixstargetindex == phixstargetindex))
    {
      return prev_bfrate_normed[nonemptymgi * nbfcontinua + i];
    }

    if (nu_edge > globals::allcont[i].nu_edge)
      break;
  }
  printout("no bf rate for element Z=%d ion_stage %d lower %d phixstargetindex %d\n", get_element(element), get_ionstage(element, lowerion), lower, phixstargetindex);
  return 0.;
#endif
}

#ifndef FORCE_LTE
void normalise_nuJ(const int modelgridindex, const double estimator_normfactor_over4pi)
{
  assert_always(std::isfinite(nuJ[modelgridindex]));
  nuJ[modelgridindex] *= estimator_normfactor_over4pi;
}
#endif


double get_T_R_from_J(const int modelgridindex)
{
  const double T_R = pow(PI / STEBO * J[modelgridindex], 1. / 4.);
  if (!std::isfinite(T_R))
  {
    /// keep old value of T_R
    printout("[warning] get_T_R_from_J: T_R estimator infinite in cell %d, use value of last timestep\n", modelgridindex);
    return grid::get_TR(modelgridindex);
  }
  /// Make sure that T is in the allowed temperature range.
  else if (T_R > MAXTEMP)
  {
    printout("[warning] get_T_R_from_J: T_R would be %.1f > MAXTEMP. Clamping to MAXTEMP = %.0f K\n", T_R, MAXTEMP);
    return MAXTEMP;
  }
  else if (T_R < MINTEMP)
  {
    printout("[warning] get_T_R_from_J: T_R would be %.1f < MINTEMP. Clamping to MINTEMP = %.0f K\n", T_R, MINTEMP);
    return MINTEMP;
  }
  else
    return T_R;
}


#ifdef DO_TITER
void titer_J(const int modelgridindex)
{
  if (J_reduced_save[modelgridindex] >= 0)
  {
    J[modelgridindex] = (J[modelgridindex] + J_reduced_save[modelgridindex]) / 2;
  }
  J_reduced_save[modelgridindex] = J[modelgridindex];
}


#ifndef FORCE_LTE
void titer_nuJ(const int modelgridindex)
{
  if (nuJ_reduced_save[modelgridindex] >= 0)
  {
    nuJ[modelgridindex] = (nuJ[modelgridindex] + nuJ_reduced_save[modelgridindex]) / 2;
  }
  nuJ_reduced_save[modelgridindex] = nuJ[modelgridindex];
}
#endif
#endif


#ifdef MPI_ON
void reduce_estimators(void)
// reduce and broadcast (allreduce) the estimators for J and nuJ in all bins
{
  MPI_Allreduce(MPI_IN_PLACE, J, grid::get_npts_model(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #ifndef FORCE_LTE
  MPI_Allreduce(MPI_IN_PLACE, nuJ, grid::get_npts_model(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif

  #if (DETAILED_BF_ESTIMATORS_ON)
  {
    MPI_Allreduce(MPI_IN_PLACE, bfrate_raw, grid::get_nonempty_npts_model() * globals::nbfcontinua, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  }
  #endif

  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    const time_t sys_time_start_reduction = time(NULL);
    printout("Reducing binned radiation field estimators");
    assert_always(radfieldbins != NULL);

    for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++)
    {
      // printout("DEBUGCELLS: cell %d associated_cells %d\n", modelgridindex, grid::get_numassociatedcells(modelgridindex));
      if (grid::get_numassociatedcells(modelgridindex) > 0)
      {
        const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
        for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
        {
          const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
          // printout("MPI: pre-MPI_Allreduce, this process modelgrid %d binindex %d has a individual contribcount of %d\n",modelgridindex,binindex,radfieldbins[mgibinindex].contribcount);
          MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[mgibinindex].J_raw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[mgibinindex].nuJ_raw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[mgibinindex].contribcount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

          // printout("MPI: After MPI_Allreduce: this process modelgrid %d binindex %d has a contribcount of %d\n",modelgridindex,binindex,radfieldbins[mgibinindex].contribcount);
        }
      }
    }
    const int duration_reduction = time(NULL) - sys_time_start_reduction;
    printout(" (took %d s)\n", duration_reduction);
  }

  if (DETAILED_LINE_ESTIMATORS_ON)
  {
    const time_t sys_time_start_reduction = time(NULL);
    printout("Reducing detailed line estimators");

    for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++)
    {
      // printout("DEBUGCELLS: cell %d associated_cells %d\n", modelgridindex, grid::get_numassociatedcells(modelgridindex));
      if (grid::get_numassociatedcells(modelgridindex) > 0)
      {
        for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++)
        {
          MPI_Allreduce(MPI_IN_PLACE, &Jb_lu_raw[modelgridindex][jblueindex].value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          MPI_Allreduce(MPI_IN_PLACE, &Jb_lu_raw[modelgridindex][jblueindex].contribcount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        }
      }
    }
    const int duration_reduction = time(NULL) - sys_time_start_reduction;
    printout(" (took %d s)\n", duration_reduction);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


void do_MPI_Bcast(const int modelgridindex, const int root, int root_node_id)
// broadcast computed radfield results including parameters
// from the cells belonging to root process to all processes
{
  MPI_Bcast(&J_normfactor[modelgridindex], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  if (grid::get_numassociatedcells(modelgridindex) > 0)
  {
    const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
    if (MULTIBIN_RADFIELD_MODEL_ON)
    {
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      {
        const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
        if (globals::rank_in_node == 0)
        {
          MPI_Bcast(&radfieldbin_solutions[mgibinindex].W, 1, MPI_FLOAT, root_node_id, globals::mpi_comm_internode);
          MPI_Bcast(&radfieldbin_solutions[mgibinindex].T_R, 1, MPI_FLOAT, root_node_id, globals::mpi_comm_internode);
        }
        MPI_Bcast(&radfieldbins[mgibinindex].J_raw, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Bcast(&radfieldbins[mgibinindex].nuJ_raw, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Bcast(&radfieldbins[mgibinindex].contribcount, 1, MPI_INT, root, MPI_COMM_WORLD);
      }
    }

    #if (DETAILED_BF_ESTIMATORS_ON)
    {
      if (globals::rank_in_node == 0)
      {
        MPI_Bcast(&prev_bfrate_normed[nonemptymgi * globals::nbfcontinua], globals::nbfcontinua, MPI_FLOAT, root_node_id, globals::mpi_comm_internode);
      }
    }
    #endif

    if (DETAILED_LINE_ESTIMATORS_ON)
    {
      for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++)
      {
        MPI_Bcast(&prev_Jb_lu_normed[modelgridindex][jblueindex].value, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Bcast(&prev_Jb_lu_normed[modelgridindex][jblueindex].contribcount, 1, MPI_INT, root, MPI_COMM_WORLD);
      }
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}
#endif


void write_restart_data(FILE *gridsave_file)
{
  printout("binned radiation field and detailed lines, ");

  fprintf(gridsave_file, "%d\n", 30490824); // special number marking the beginning of radfield data

  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    fprintf(gridsave_file, "%d %lg %lg %lg %lg\n",
            RADFIELDBINCOUNT, nu_lower_first_initial, nu_upper_last_initial,
            T_R_min, T_R_max);

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      fprintf(gridsave_file,"%d %lg\n", binindex, radfieldbin_nu_upper[binindex]);
    }
  }

  #if (DETAILED_BF_ESTIMATORS_ON)
  {
    const int nbfcontinua = globals::nbfcontinua;
    fprintf(gridsave_file, "%d\n", nbfcontinua);

    for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++)
    {
      if (grid::get_numassociatedcells(modelgridindex) > 0)
      {
        const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
        fprintf(gridsave_file, "%d\n", modelgridindex);
        for (int i = 0; i < nbfcontinua; i++)
        {
          fprintf(gridsave_file, "%g ", prev_bfrate_normed[nonemptymgi * nbfcontinua + i]);
        }
      }
    }
  }
  #endif

  if (DETAILED_LINE_ESTIMATORS_ON)
  {
    fprintf(gridsave_file, "%d\n", detailed_linecount);

    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++)
    {
      fprintf(gridsave_file, "%d ", detailed_lineindicies[jblueindex]);
    }
  }

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++)
  {
    if (grid::get_numassociatedcells(modelgridindex) > 0)
    {
      const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
      assert_testmodeonly(nonemptymgi >= 0);
      fprintf(gridsave_file,"%d %lg\n", modelgridindex, J_normfactor[modelgridindex]);

      if (MULTIBIN_RADFIELD_MODEL_ON)
      {
        for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
        {
          const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
          fprintf(gridsave_file, "%lg %lg %g %g %d\n",
                  radfieldbins[mgibinindex].J_raw,
                  radfieldbins[mgibinindex].nuJ_raw,
                  radfieldbin_solutions[mgibinindex].W,
                  radfieldbin_solutions[mgibinindex].T_R,
                  radfieldbins[mgibinindex].contribcount);
                  //radfieldbins[mgibinindex].fit_type
        }
      }

      if (DETAILED_LINE_ESTIMATORS_ON)
      {
        for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++)
        {
          fprintf(gridsave_file, "%lg %d\n",
                  Jb_lu_raw[modelgridindex][jblueindex].value,
                  Jb_lu_raw[modelgridindex][jblueindex].contribcount);
        }
      }
    }
  }
  fprintf(gridsave_file, "%d\n", 42809403); // special number marking the end of radfield data
}


void read_restart_data(FILE *gridsave_file)
{
  printout("Reading restart data for radiation field\n");

  if (!initialized)
  {
    printout("ERROR: Radiation field has not been initialised yet. Can't read saved state.\n");
    abort();
  }

  int code_check;
  fscanf(gridsave_file, "%d\n", &code_check);
  if (code_check != 30490824)
  {
    printout("ERROR: Beginning of radfield restart data not found! Found %d instead of 30490824\n", code_check);
    abort();
  }

  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    int bincount_in;
    double T_R_min_in, T_R_max_in, nu_lower_first_initial_in, nu_upper_last_initial_in;
    fscanf(gridsave_file,"%d %lg %lg %lg %lg\n",
           &bincount_in, &nu_lower_first_initial_in, &nu_upper_last_initial_in,
           &T_R_min_in, &T_R_max_in);

    double nu_lower_first_ratio = nu_lower_first_initial_in / nu_lower_first_initial;
    if (nu_lower_first_ratio > 1.0) nu_lower_first_ratio = 1 / nu_lower_first_ratio;

    double nu_upper_last_ratio = nu_upper_last_initial_in / nu_upper_last_initial;
    if (nu_upper_last_ratio > 1.0) nu_upper_last_ratio = 1 / nu_upper_last_ratio;

    if (bincount_in != RADFIELDBINCOUNT || T_R_min_in != T_R_min || T_R_max_in != T_R_max ||
        nu_lower_first_ratio < 0.999 || nu_upper_last_ratio < 0.999)
    {
      printout("ERROR: gridsave file specifies %d bins, nu_lower_first_initial %lg nu_upper_last_initial %lg T_R_min %lg T_R_max %lg\n",
               bincount_in, nu_lower_first_initial_in, nu_upper_last_initial_in, T_R_min_in, T_R_max_in);
      printout("require %d bins, nu_lower_first_initial %lg nu_upper_last_initial %lg T_R_min %lg T_R_max %lg\n",
              RADFIELDBINCOUNT, nu_lower_first_initial, nu_upper_last_initial, T_R_min, T_R_max);
      abort();
    }

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      int binindex_in;
      fscanf(gridsave_file,"%d %lg\n", &binindex_in, &radfieldbin_nu_upper[binindex]);
      assert_always(binindex_in == binindex);
    }
  }

  #if (DETAILED_BF_ESTIMATORS_ON)
  {
    int gridsave_nbf_in;
    fscanf(gridsave_file, "%d\n", &gridsave_nbf_in);
    assert_always(gridsave_nbf_in == globals::nbfcontinua);

    for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++)
    {
      if (grid::get_numassociatedcells(modelgridindex) > 0)
      {
        const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
        int mgi_in;
        fscanf(gridsave_file, "%d\n", &mgi_in);
        assert_always(mgi_in == modelgridindex);
        for (int i = 0; i < globals::nbfcontinua; i++)
        {
          float bfrate_normed = 0;
          fscanf(gridsave_file, "%g ", &bfrate_normed);

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
  #endif

  if (DETAILED_LINE_ESTIMATORS_ON)
  {
    int detailed_linecount_in;
    fscanf(gridsave_file,"%d\n", &detailed_linecount_in);

    if (detailed_linecount_in != detailed_linecount)
    {
      printout("ERROR: gridsave file specifies %d detailed lines but this simulation has %d.\n",
               detailed_linecount_in, detailed_linecount);
      abort();
    }

    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++)
    {
      fscanf(gridsave_file, "%d ", &detailed_lineindicies[jblueindex]);
    }
  }

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++)
  {
    if (grid::get_numassociatedcells(modelgridindex) > 0)
    {
      const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
      int mgi_in;
      fscanf(gridsave_file,"%d %lg\n", &mgi_in, &J_normfactor[modelgridindex]);
      if (mgi_in != modelgridindex)
      {
        printout("ERROR: expected data for cell %d but found cell %d\n", modelgridindex, mgi_in);
        abort();
      }

      if (MULTIBIN_RADFIELD_MODEL_ON)
      {
        for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
        {
          const int mgibinindex = nonemptymgi * RADFIELDBINCOUNT + binindex;
          float W = 0;
          float T_R = 0;
          fscanf(gridsave_file, "%lg %lg %g %g %d\n",
                 &radfieldbins[mgibinindex].J_raw,
                 &radfieldbins[mgibinindex].nuJ_raw,
                 &W,
                 &T_R,
                 &radfieldbins[mgibinindex].contribcount);
#ifdef MPI_ON
          if (globals::rank_in_node == 0)
#endif
          {
            radfieldbin_solutions[mgibinindex].W = W;
            radfieldbin_solutions[mgibinindex].T_R = T_R;
          }
        }
      }

      if (DETAILED_LINE_ESTIMATORS_ON)
      {
        for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++)
        {
          fscanf(gridsave_file, "%lg %d\n",
                  &Jb_lu_raw[modelgridindex][jblueindex].value,
                  &Jb_lu_raw[modelgridindex][jblueindex].contribcount);
        }
      }
    }
  }
  fscanf(gridsave_file, "%d\n", &code_check);
  if (code_check != 42809403)
  {
    printout("ERROR: End of radfield restart data not found! Found %d instead of 42809403\n", code_check);
    char line[1025];
    fgets(line, 1024, gridsave_file);
    printout("%s\n", line);
    fgets(line, 1024, gridsave_file);
    printout("%s\n", line);
    fgets(line, 1024, gridsave_file);
    printout("%s\n", line);
    abort();
  }
}


// not in use, but could potential improve speed and accuracy of integrating
// across the binned radiation field which is discontinuous at the bin boundaries
inline
int integrate(
  const gsl_function *f, double nu_a, double nu_b, double epsabs,
  double epsrel, size_t limit, int key, gsl_integration_workspace *workspace,
  double *result, double *abserr)
{
  if (MULTIBIN_RADFIELD_MODEL_ON && (globals::nts_global >= FIRST_NLTE_RADFIELD_TIMESTEP))
  {
    double *pts = (double *) malloc((RADFIELDBINCOUNT + 3) * sizeof(double));
    int binindex_a = select_bin(nu_a);
    int binindex_b = select_bin(nu_b);
    int npts = 0;
    pts[npts++] = nu_a;
    if (binindex_a == binindex_b) // both higher, both lower, or match the same bin
    {
      // region doesn't contain any bins
      pts[npts++] = nu_b;
    }
    else
    {
      if (binindex_a < 0) // a is below the first bin
      {
        binindex_a = 0;
        pts[npts++] = get_bin_nu_lower(0);
      }

      const int maxbinplusone = (binindex_b < 0) ? RADFIELDBINCOUNT : binindex_b;

      for (int binindex = binindex_a; binindex < maxbinplusone; binindex++)
        pts[npts++] = get_bin_nu_upper(binindex);

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
  else
    return gsl_integration_qag(f, nu_a, nu_b, epsabs, epsrel, limit, key, workspace, result, abserr);
}

}