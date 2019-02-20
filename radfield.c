#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_debye.h>
#include <gsl/gsl_roots.h>
#include "atomic.h"
#include "grid_init.h"
#include "radfield.h"
#include "sn3d.h"

#define RADFIELDBINCOUNT 168

static const int FIRST_NLTE_RADFIELD_TIMESTEP = 13;

static const double nu_lower_first_initial = (CLIGHT / (10000e-8)); // in Angstroms
static const double nu_upper_last_initial = (CLIGHT /  (50e-8));  // in Angstroms

static const double boost_region_nu_lower = (CLIGHT / (2500e-8)); // in Angstroms
static const double boost_region_nu_upper = (CLIGHT / (2100e-8));  // in Angstroms
static const double boost_region_factor = 1.0;
static const bool boost_region_on = false;

static double J_normfactor[MMODELGRID + 1];

static bool radfield_initialized = false;

static const double T_R_min = 2000;
static const double T_R_max = 280000;

// typedef enum
// {
//   FIT_DILUTED_BLACKBODY = 0,
//   FIT_CONSTANT = 1,
// } enum_bin_fit_type;

struct radfieldbin
{
  double J_raw;           // value needs to be multipled by J_normfactor to get the true value
  double nuJ_raw;
  int contribcount;

  // these two parameters are used in the current timestep, but were calculated
  // from the values of J and nuJ in the previous timestep
  float W;                // dilution (scaling) factor
  float T_R;              // radiation temperature
  // enum_bin_fit_type fit_type;
};

static double radfieldbin_nu_upper[RADFIELDBINCOUNT]; // array of upper frequency boundaries of bins
static struct radfieldbin *radfieldbins[MMODELGRID + 1];

// ** Detailed lines - Jblue_lu estimators for selected lines

struct Jb_lu_estimator
{
  double value;
  int contribcount;
};

// reallocate the detailed line arrays in units of BLOCKSIZEJBLUE
static const int BLOCKSIZEJBLUE = 128;
static int detailed_linecount = 0;

// array of indicies into the linelist[] array for selected lines
static int *detailed_lineindicies;

static struct Jb_lu_estimator *prev_Jb_lu_normed[MMODELGRID + 1];  // value from the previous timestep
static struct Jb_lu_estimator *Jb_lu_raw[MMODELGRID + 1];   // unnormalised estimator for the current timestep

// ** end detailed lines

#if (DETAILED_BF_ESTIMATORS_ON)
static bool normed_bfrates_available = false;
static float *prev_bfrate_normed[MMODELGRID + 1];  // values from the previous timestep
static double *bfrate_raw[MMODELGRID + 1];   // unnormalised estimators for the current timestep
#endif

static double J[MMODELGRID + 1];
#ifdef DO_TITER
  static double J_reduced_save[MMODELGRID + 1];
#endif

// J and nuJ are accumulated and then normalised in-place
// i.e. be sure the normalisation has been applied (exactly once) before using the values here!
#ifndef FORCE_LTE
  static double nuJ[MMODELGRID + 1];
  #ifdef DO_TITER
    static double nuJ_reduced_save[MMODELGRID + 1];
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

static FILE *restrict radfieldfile = NULL;


extern inline double radfield_dbb(double nu, float T, float W);


static inline
double get_bin_nu_upper(int binindex)
{
  return radfieldbin_nu_upper[binindex];
}


static void
setup_bin_boundaries(void)
{
  double prev_nu_upper = nu_lower_first_initial;

  // choose between equally spaced in energy/frequency or wavelength (before bf edges shift boundaries around)
  const double delta_nu = (nu_upper_last_initial - nu_lower_first_initial) / RADFIELDBINCOUNT;
  // const double lambda_lower_first_initial = 1e8 * CLIGHT / nu_lower_first_initial;
  // const double lambda_upper_last_initial = 1e8 * CLIGHT / nu_upper_last_initial;
  // const double delta_lambda = (lambda_upper_last_initial - lambda_lower_first_initial) / RADFIELDBINCOUNT;

  for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
  {
    // radfieldbin_nu_upper[binindex] = 1e8 * CLIGHT / (lambda_lower_first_initial + (binindex + 1) * delta_lambda);
    radfieldbin_nu_upper[binindex] = nu_lower_first_initial + (binindex + 1) * delta_nu;

    // Align the bin edges with bound-free edges, except for the last one
    if (binindex < RADFIELDBINCOUNT - 1)
    {
      for (int i = 0; i < nbfcontinua_ground; i++)
      {
        const double nu_edge = phixslist[tid].groundcont[i].nu_edge;
        const double eV_edge = H * nu_edge / EV;
        const double angstrom_edge = 1e8 * CLIGHT / nu_edge;
        const int element = phixslist[tid].groundcont[i].element;
        const int ion = phixslist[tid].groundcont[i].ion;
        const int level = phixslist[tid].groundcont[i].level;
        const int phixstargetindex = phixslist[tid].groundcont[i].phixstargetindex;

        const int Z = get_element(element);
        const int ion_stage = get_ionstage(element, ion);
        const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);

        //printout("bf edge at %g, nu_lower_first %g, nu_upper_last %g\n",nu_edge,nu_lower_first,nu_upper_last);

        // this is compares the highest and lowest bins to the bound-free list, only do it once, i.e. binindex == 0
        if (binindex == 0 && ((nu_edge < nu_lower_first_initial) || (nu_edge > nu_upper_last_initial)))
        {
          printout("Missed bf edge at %12.5e Hz (%6.2f eV, %6.1f A), nu_lower_first %11.5e Hz, nu_upper_last %11.5e Hz, Z=%d ion_stage %d level %d upperionlevel %d\n",
                   nu_edge, eV_edge, angstrom_edge, nu_lower_first_initial, nu_upper_last_initial, Z, ion_stage, level, upperionlevel);
        }

        const double bin_nu_upper = get_bin_nu_upper(binindex);
        if ((nu_edge > prev_nu_upper) && (nu_edge < bin_nu_upper))
        {
          printout("Shifting bin %d nu_upper from %12.5e Hz to bf edge at %12.5e Hz (%6.2f eV, %6.1f A) for Z=%d ion_stage %d level %d upperionlevel %d\n",
                   binindex, bin_nu_upper, nu_edge, eV_edge, angstrom_edge, Z, ion_stage, level, upperionlevel);
          radfieldbin_nu_upper[binindex] = nu_edge;
        }
      }
    }
    prev_nu_upper = get_bin_nu_upper(binindex);
  }
}


void radfield_jblue_init(void)
// set up the data structures assocated with the Jb_lu estimators
// this is called before/during atomic data input, whereas the
// radfield_init() is called after the atomic data is known
{
  detailed_linecount = 0;

  detailed_lineindicies = NULL;
  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    prev_Jb_lu_normed[modelgridindex] = NULL;
    Jb_lu_raw[modelgridindex] = NULL;
  }
}


static void realloc_detailed_lines(const int new_size)
{
  detailed_lineindicies = realloc(detailed_lineindicies, new_size * sizeof(int));
  if (detailed_lineindicies == NULL)
  {
    printout("ERROR: Not enough memory to reallocate detailed Jblue estimator line list\n");
    abort();
  }

  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    if (get_numassociatedcells(modelgridindex) > 0)
    {
      prev_Jb_lu_normed[modelgridindex] = realloc(
        prev_Jb_lu_normed[modelgridindex], new_size * sizeof(struct Jb_lu_estimator));

      Jb_lu_raw[modelgridindex] = realloc(
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

  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    if (get_numassociatedcells(modelgridindex) > 0)
    {
      prev_Jb_lu_normed[modelgridindex][detailed_linecount].value = 0;
      prev_Jb_lu_normed[modelgridindex][detailed_linecount].contribcount = 0;

      // radfield_zero_estimators should do the next part anyway, but just to be sure:
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


void radfield_init(int my_rank)
// this should be called only after the atomic data is in memory
{
  if (radfield_initialized)
  {
    printout("ERROR: Tried to initialize radfield twice!\n");
    abort();
  }

  radfield_initialized = true;

  if (DETAILED_LINE_ESTIMATORS_ON)
  {
    for (int i = 0; i < nlines; i++)
    {
      const int element = linelist[i].elementindex;
      const int Z = get_element(element);
      if (Z == 26)
      {
        const int lowerlevel = linelist[i].lowerlevelindex;
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

  if (!MULTIBIN_RADFIELD_MODEL_ON)
  {
    printout("The radiation field model is a whole-spectrum fit to a single diluted blackbody.\n");
    return;
  }

  printout("The multibin radiation field estimators are being used instead of the whole-spectrum fit from timestep %d onwards.\n", FIRST_NLTE_RADFIELD_TIMESTEP);

  printout("Initialising multibin radiation field with %d bins from (%.2f eV, %6.1f A) to (%.2f eV, %6.1f A)\n",
           RADFIELDBINCOUNT, H * nu_lower_first_initial / EV, 1e8 * CLIGHT / nu_lower_first_initial,
           H * nu_upper_last_initial / EV, 1e8 * CLIGHT / nu_upper_last_initial);
  char filename[100];
  sprintf(filename,"radfield_%.4d.out", my_rank);
  radfieldfile = fopen_required(filename, "w");
  fprintf(radfieldfile,"%8s %15s %8s %11s %11s %9s %9s %9s %9s %9s %12s\n",
          "timestep","modelgridindex","bin_num","nu_lower","nu_upper",
          "nuJ","J","J_nu_avg","ncontrib","T_R","W");
  fflush(radfieldfile);

  setup_bin_boundaries();

  long radfield_mem_usage = 0;
  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    radfield_set_J_normfactor(modelgridindex, -1.0);
    // printout("DEBUGCELLS: cell %d associated_cells %d\n", modelgridindex, get_numassociatedcells(modelgridindex));
    if (get_numassociatedcells(modelgridindex) > 0)
    {
      #if (DETAILED_BF_ESTIMATORS_ON)
      {
        bfrate_raw[modelgridindex] = calloc(nbfcontinua, sizeof(double));
        prev_bfrate_normed[modelgridindex] = calloc(nbfcontinua, sizeof(float));
      }
      #endif

      radfieldbins[modelgridindex] = (struct radfieldbin *) calloc(RADFIELDBINCOUNT, sizeof(struct radfieldbin));
      radfield_mem_usage += RADFIELDBINCOUNT * sizeof(struct radfieldbin);
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
      {
        radfieldbins[modelgridindex][binindex].J_raw = 0.;
        radfieldbins[modelgridindex][binindex].nuJ_raw = 0.;
        radfieldbins[modelgridindex][binindex].contribcount = 0;
        radfieldbins[modelgridindex][binindex].W = -1.;
        radfieldbins[modelgridindex][binindex].T_R = -1.;
        // radfieldbins[modelgridindex][binindex].fit_type = FIT_DILUTED_BLACKBODY;
      }
    }
  }
  printout("mem_usage: radiation field bins for non-empty cells occupy %.2f MB\n", radfield_mem_usage / 1024. / 1024.);
}


/// Initialise estimator arrays which hold the last time steps values (used to damp out
/// fluctuations over timestep iterations if DO_TITER is defined) to -1.
void initialise_photoionestimators(void)
{
  //for (n = 0; n < ngrid; n++)
  for (int n = 0; n < npts_model; n++)
  {
    //double T_e = get_Te(n);
    #ifdef DO_TITER
      J_reduced_save[n] = -1.;
    #endif
    #ifndef FORCE_LTE
      #ifdef DO_TITER
        nuJ_reduced_save[n] = -1.;
        ffheatingestimator_save[n] = -1.;
        colheatingestimator_save[n] = -1.;
      #endif
      for (int element = 0; element < nelements; element++)
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
          //gammaestimator[n*nelements*maxion+element*maxion+ion] = Gamma; //gamma_lte/zeta;
          ////corrphotoionrenorm[n*nelements*maxion+element*maxion+ion] = 1.;
          ////photoionestimator[n*nelements*maxion+element*maxion+ion] = Gamma; //gamma_lte/zeta;

          #ifdef DO_TITER
            gammaestimator_save[n*nelements*maxion+element*maxion+ion] = -1.;
            if (!NO_LUT_BFHEATING)
              bfheatingestimator_save[n*nelements*maxion+element*maxion+ion] = -1.;
            /*
            photoionestimator_save[n*nelements*maxion+element*maxion+ion] = -1.;
            stimrecombestimator_save[n*nelements*maxion+element*maxion+ion] = -1.;
            ionfluxestimator_save[n*nelements*maxion+element*maxion+ion] = -1.;
            */
          #endif
        }
      }
    #endif
  }
}


inline
int radfield_get_Jblueindex(const int lineindex)
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


inline
double radfield_get_Jb_lu(const int modelgridindex, const int jblueindex)
{
  assert(jblueindex >= 0);
  assert(jblueindex < detailed_linecount);
  return prev_Jb_lu_normed[modelgridindex][jblueindex].value;
}


inline
int radfield_get_Jb_lu_contribcount(const int modelgridindex, const int jblueindex)
{
  assert(jblueindex >= 0);
  assert(jblueindex < detailed_linecount);
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
  else if (modelgridindex >= MMODELGRID)
  {
    printout("radfield: Fatal error: get_bin_J called before on modelgridindex %d >= MMODELGRID",modelgridindex);
    abort();
  }

  return radfieldbins[modelgridindex][binindex].J_raw * J_normfactor[modelgridindex];
}


static void set_bin_J(int modelgridindex, int binindex, double value)
// set the normalised J_nu
{
  if (J_normfactor[modelgridindex] <= 0.0)
  {
    printout("radfield: Fatal error: set_bin_J called before J_normfactor set for modelgridindex %d",modelgridindex);
    abort();
  }
  else if (modelgridindex >= MMODELGRID)
  {
    printout("radfield: Fatal error: set_bin_J called before on modelgridindex %d >= MMODELGRID",modelgridindex);
    abort();
  }
  radfieldbins[modelgridindex][binindex].J_raw = value / J_normfactor[modelgridindex];
}


static double get_bin_nuJ(int modelgridindex, int binindex)
{
  if (J_normfactor[modelgridindex] <= 0.0)
  {
    printout("radfield: Fatal error: get_bin_nuJ called before J_normfactor set for modelgridindex %d",modelgridindex);
    abort();
  }
  else if (modelgridindex >= MMODELGRID)
  {
    printout("radfield: Fatal error: get_bin_nuJ called before on modelgridindex %d >= MMODELGRID",modelgridindex);
    abort();
  }
  return radfieldbins[modelgridindex][binindex].nuJ_raw * J_normfactor[modelgridindex];
}


static inline
double get_bin_nu_bar(int modelgridindex, int binindex)
// importantly, this is average beween the current and previous timestep
{
  const double nuJ_sum = get_bin_nuJ(modelgridindex, binindex);
  const double J_sum = get_bin_J(modelgridindex, binindex);
  return nuJ_sum / J_sum;
}


static inline
double get_bin_nu_lower(int binindex)
{
  if (binindex > 0)
    return radfieldbin_nu_upper[binindex - 1];
  else
    return nu_lower_first_initial;
}


static inline
int get_bin_contribcount(int modelgridindex, int binindex)
{
  return radfieldbins[modelgridindex][binindex].contribcount;
}


static inline
float get_bin_W(int modelgridindex, int binindex)
{
  return radfieldbins[modelgridindex][binindex].W;
}


static inline
float get_bin_T_R(int modelgridindex, int binindex)
{
  return radfieldbins[modelgridindex][binindex].T_R;
}


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
    assert(false);
  }
}


void radfield_write_to_file(int modelgridindex, int timestep)
{
  assert(MULTIBIN_RADFIELD_MODEL_ON);

# ifdef _OPENMP
# pragma omp critical (radfield_out_file)
  {
# endif
    if (!radfield_initialized)
    {
      printout("Call to radfield_write_to_file before radfield_init\n");
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
        T_R = get_TR(modelgridindex);
        W = get_W(modelgridindex);
        contribcount = totalcontribs;
      }
      else // use binindex < -1 for detailed line Jb_lu estimators
      {
        const int jblueindex = -2 - binindex; // -2 is the first detailed line, -3 is the second, etc
        const int lineindex = detailed_lineindicies[jblueindex];
        const double nu_trans = linelist[lineindex].nu;
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


void radfield_close_file(void)
{
  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    fclose(radfieldfile);

    for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
    {
      if (get_numassociatedcells(modelgridindex) > 0)
      {
        free(radfieldbins[modelgridindex]);
        #if (DETAILED_BF_ESTIMATORS_ON)
        free(bfrate_raw[modelgridindex]);
        free(prev_bfrate_normed[modelgridindex]);
        #endif
      }
    }
  }
}


void radfield_zero_estimators(int modelgridindex)
// set up the new bins and clear the estimators in preparation
// for a timestep
{
  #if (DETAILED_BF_ESTIMATORS_ON)
  if (radfield_initialized && (get_numassociatedcells(modelgridindex) > 0))
  {
    for (int i = 0; i < nbfcontinua; i++)
    {
      bfrate_raw[modelgridindex][i] = 0.;
    }
  }
  #endif

  for (int i = 0; i < detailed_linecount; i++)
  {
    Jb_lu_raw[modelgridindex][i].value = 0.;
    Jb_lu_raw[modelgridindex][i].contribcount = 0.;
  }

  J[modelgridindex] = 0.; // this is required even if FORCE_LTE is on
#ifndef FORCE_LTE
  nuJ[modelgridindex] = 0.;

  if (MULTIBIN_RADFIELD_MODEL_ON && radfield_initialized && (get_numassociatedcells(modelgridindex) > 0))
  {
    // printout("radfield: zeroing estimators in %d bins in cell %d\n",RADFIELDBINCOUNT,modelgridindex);

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      radfieldbins[modelgridindex][binindex].J_raw = 0.0;
      radfieldbins[modelgridindex][binindex].nuJ_raw = 0.0;
      radfieldbins[modelgridindex][binindex].contribcount = 0;
    }
    radfield_set_J_normfactor(modelgridindex, -1.0);
  }
#endif
}


#if (DETAILED_BF_ESTIMATORS_ON)
static void radfield_increment_bfestimators(const int modelgridindex, const double distance_e_cmf, const double nu_cmf)
{
  const double distance_e_cmf_over_nu = distance_e_cmf / nu_cmf;
  for (int allcontindex = 0; allcontindex < nbfcontinua; allcontindex++)
  {
    const double nu_edge = phixslist[tid].allcont[allcontindex].nu_edge;
    if (nu_cmf >= nu_edge)
    {
      #ifdef _OPENMP
        #pragma omp atomic
      #endif
      bfrate_raw[modelgridindex][allcontindex] += phixslist[tid].allcont[allcontindex].gamma_contr * distance_e_cmf_over_nu;
    }
    else
    {
      // list is sorted by nu_edge, so all remaining will have nu_cmf < nu_edge
      break;
    }
  }
}
#endif


inline
void radfield_update_estimators(int modelgridindex, double distance_e_cmf, double nu_cmf)
{
  #ifdef _OPENMP
    #pragma omp atomic
  #endif
  J[modelgridindex] += distance_e_cmf;
  #ifdef DEBUG_ON
    if (!isfinite(J[modelgridindex]))
    {
      printout("[fatal] update_estimators: estimator becomes non finite: distance_e_cmf %g, nu_cmf %g ... abort\n",distance_e_cmf,nu_cmf);
      abort();
    }
  #endif

#ifndef FORCE_LTE
  #ifdef _OPENMP
    #pragma omp atomic
  #endif
  nuJ[modelgridindex] += distance_e_cmf * nu_cmf;
  #ifdef DEBUG_ON
    if (!isfinite(nuJ[modelgridindex]))
    {
      printout("[fatal] update_estimators: estimator becomes non finite: distance_e_cmf %g, nu_cmf %g ... abort\n",distance_e_cmf,nu_cmf);
      abort();
    }
  #endif

  #if (DETAILED_BF_ESTIMATORS_ON)
  radfield_increment_bfestimators(modelgridindex, distance_e_cmf, nu_cmf);
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
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      radfieldbins[modelgridindex][binindex].J_raw += distance_e_cmf;
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      radfieldbins[modelgridindex][binindex].nuJ_raw += distance_e_cmf * nu_cmf;
      #ifdef _OPENMP
      #pragma omp atomic
      #endif
      radfieldbins[modelgridindex][binindex].contribcount += 1;
    }
    // else
    // {
    //   printout("WARNING: radfield_update_estimators dropping packet contribution for nu_cmf %g\n",
    //            nu_cmf);
    //   printout("           modelgridindex %d binindex %d nu_lower_first %g nu_upper_last %g \n",
    //            modelgridindex, binindex, nu_lower_first, get_bin_nu_upper(modelgridindex,RADFIELDBINCOUNT - 1));
    // }
  }
#endif
}


void radfield_increment_lineestimator(const int modelgridindex, const int lineindex, const double increment)
{
  if (!DETAILED_LINE_ESTIMATORS_ON) return;

  const int jblueindex = radfield_get_Jblueindex(lineindex);
  if (jblueindex >= 0)
  {
    Jb_lu_raw[modelgridindex][jblueindex].value += increment;
    Jb_lu_raw[modelgridindex][jblueindex].contribcount += 1;
    // const int lineindex = detailed_lineindicies[jblueindex];
    // printout(" increment cell %d lineindex %d Jb_lu_raw %g prev_Jb_lu_normed %g radfield(nu_trans) %g\n",
    //       modelgridindex, lineindex, Jb_lu_raw[modelgridindex][jblueindex], prev_Jb_lu_normed[modelgridindex][jblueindex].value, radfield(linelist[lineindex].nu, modelgridindex));
  }
}

double radfield(double nu, int modelgridindex)
// mean intensity J_nu
{
  // return 0.;
  if (MULTIBIN_RADFIELD_MODEL_ON && (nts_global >= FIRST_NLTE_RADFIELD_TIMESTEP))
  {
    const int binindex = select_bin(nu);
    if (binindex >= 0)
    {
      const struct radfieldbin *restrict const bin = &radfieldbins[modelgridindex][binindex];
      if (bin->W >= 0.)
      {
        // if (bin->fit_type == FIT_DILUTED_BLACKBODY)
        {
          const double J_nu = radfield_dbb(nu, bin->T_R, bin->W);
          /*if (fabs(J_nu / J_nu_fullspec - 1.0) > 0.5)
          {
            printout("WARNING: radfield: significant discrepancy. J_nu_fullspec %g, J_nu %g, nu %g bin->W %g bin->T_R %g\n",
                     J_nu_fullspec, J_nu, nu, bin->W, bin->T_R);
          }*/
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
      //   const double J_nu_LTE = radfield_dbb(nu, get_Te(modelgridindex), 1.0);
      //   return J_nu_LTE;
      // }
      // else
      //   return 0; // no radfield redwards of the bins
      //printout("WARNING: Radfield modelgridindex %d binindex %d nu %g nu_lower_first %g nu_upper_last %g \n",
      //         modelgridindex, binindex, nu, nu_lower_first, nu_upper_last);
    }
  }
  /*else
  {
    printout("radfield: WARNING: Radfield called before initialized. Using global T_R %g W %g nu %g modelgridindex %d\n",
             W_fullspec, T_R_fullspec, nu, modelgridindex);
  }*/

  const float T_R_fullspec = get_TR(modelgridindex);
  const float W_fullspec   = get_W(modelgridindex);
  const double J_nu_fullspec = radfield_dbb(nu, T_R_fullspec, W_fullspec);
  return J_nu_fullspec;
}


static double gsl_integrand_planck(double nu, void *restrict paras)
{
  const double T_R = ((gsl_planck_integral_paras *) paras)->T_R;
  const enum_prefactor prefactor = ((gsl_planck_integral_paras *) paras)->prefactor;

  double integrand = TWOHOVERCLIGHTSQUARED * pow(nu,3) / (expm1(HOVERKB * nu / T_R));

  if (prefactor == TIMES_NU)
    integrand *= nu;

  return integrand;
}


static double planck_integral(double T_R, double nu_lower, double nu_upper, enum_prefactor prefactor)
{
  double integral = 0.;

  double error = 0.;
  const double epsrel = 1e-10;
  const double epsabs = 0.;

  gsl_integration_workspace *restrict w = gsl_integration_workspace_alloc(65536);

  gsl_planck_integral_paras intparas;
  intparas.T_R = T_R;
  intparas.prefactor = prefactor;

  gsl_function F_planck;
  F_planck.function = &gsl_integrand_planck;
  F_planck.params = &intparas;

  gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);
  int status = gsl_integration_qag(&F_planck, nu_lower, nu_upper, epsabs, epsrel, 65536, GSL_INTEG_GAUSS61, w, &integral, &error);
  if (status != 0)
  {
    printout("planck_integral integrator status %d, GSL_FAILURE= %d. Integral value %g, setting to zero.\n", status,GSL_FAILURE,integral);
    integral = 0.;
  }
  gsl_set_error_handler(previous_handler);

  gsl_integration_workspace_free(w);

  return integral;
}


static double planck_integral_analytic(double T_R, double nu_lower, double nu_upper, enum_prefactor prefactor)
{
  double integral = 0.;

  if (prefactor == TIMES_NU)
  {
    double debye_upper = gsl_sf_debye_4(HOVERKB * nu_upper / T_R) * pow(nu_upper,4);
    double debye_lower = gsl_sf_debye_4(HOVERKB * nu_lower / T_R) * pow(nu_lower,4);
    integral = TWOHOVERCLIGHTSQUARED * (debye_upper - debye_lower) * T_R / HOVERKB / 4;
  }
  else
  {
    double debye_upper = gsl_sf_debye_3(HOVERKB * nu_upper / T_R) * pow(nu_upper,3);
    double debye_lower = gsl_sf_debye_3(HOVERKB * nu_lower / T_R) * pow(nu_lower,3);
    integral = TWOHOVERCLIGHTSQUARED * (debye_upper - debye_lower) * T_R / HOVERKB / 3;

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


static double delta_nu_bar(double T_R, void *restrict paras)
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

  if (!isfinite(nu_bar_planck))
  {
    double nu_times_planck_numerical = planck_integral(T_R, nu_lower, nu_upper, TIMES_NU);
    double planck_integral_numerical = planck_integral(T_R, nu_lower, nu_upper, ONE);
    double nu_bar_planck_numerical = nu_times_planck_numerical / planck_integral_numerical;

    printout("planck_integral_analytic is %g. Replacing with numerical result of %g.\n",nu_bar_planck,nu_bar_planck_numerical);
    nu_bar_planck = nu_bar_planck_numerical;
  }*/

  const double delta_nu_bar = nu_bar_planck_T_R - nu_bar_estimator;

  if (!isfinite(delta_nu_bar))
  {
    printout("delta_nu_bar is %g. nu_bar_planck_T_R %g nu_times_planck_numerical %g planck_integral_numerical %g nu_bar_estimator %g\n",
             nu_bar_planck_T_R, nu_times_planck_numerical, planck_integral_numerical, nu_bar_estimator);
  }

  //double delta_nu_bar = nu_bar_planck_T_R / nu_bar_estimator - 1.0;

  //printout("delta_nu_bar %g nu_bar_planck %g\n",delta_nu_bar,nu_bar_planck);

  return delta_nu_bar;
}


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

  if (!isfinite(delta_nu_bar_min) || !isfinite(delta_nu_bar_max))
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


static void set_radfield_params_fullspec(const int modelgridindex, const int timestep)
{
  const double nubar = nuJ[modelgridindex] / J[modelgridindex];
  if (!isfinite(nubar) || nubar == 0.)
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
    set_TJ(modelgridindex, T_J);

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
    set_TR(modelgridindex, T_R);

    const float W = J[modelgridindex] * PI / STEBO / pow(T_R, 4);
    set_W(modelgridindex, W);

    printout("Full-spectrum fit radfield for cell %d at timestep %d: J %g, nubar %5.1f Angstrom, T_J %g, T_R %g, W %g\n",
             modelgridindex, timestep, J[modelgridindex], 1e8 * CLIGHT / nubar,
             T_J, T_R, W);
  }
}


void radfield_fit_parameters(int modelgridindex, int timestep)
// finds the best fitting W and temperature parameters in each spectral bin
// using J and nuJ
{
  set_radfield_params_fullspec(modelgridindex, timestep);

  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    if (J_normfactor[modelgridindex] <= 0)
    {
      printout("radfield: FATAL J_normfactor = %g in cell %d at call to radfield_fit_parameters", J_normfactor[modelgridindex], modelgridindex);
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

    // int contribcount_allbins = 0;
    // for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    //   contribcount_allbins += get_bin_contribcount(modelgridindex, binindex, true);

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      const double nu_lower = get_bin_nu_lower(binindex);
      const double nu_upper = get_bin_nu_upper(binindex);
      if (boost_region_on && boost_region_nu_lower < nu_upper && boost_region_nu_upper > nu_lower)
      {
        printout("Artificially boosting bin %d\n",binindex);
        set_bin_J(modelgridindex, binindex, J_bin_max * boost_region_factor);
      }
      const double J_bin = get_bin_J(modelgridindex, binindex);
      float T_R_bin = -1.0;
      double W_bin = -1.0;
      const int contribcount = get_bin_contribcount(modelgridindex, binindex);

      if (contribcount > 0)
      {
        // // enum_bin_fit_type bin_fit_type = radfieldbins[modelgridindex][binindex].fit_type;
        // if (bin_fit_type == FIT_DILUTED_BLACKBODY)
        {
          T_R_bin = find_T_R(modelgridindex, binindex);

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
        //   printout("radfield_fit_parameters: unknown fit type %d for bin %d\n", bin_fit_type, binindex);
        //   T_R_bin = -1.;
        //   W_bin = -1.;
        // }
      }
      else if (contribcount == 0)
      {
        T_R_bin = 0.;
        W_bin = 0.;
      }
      // else
      // {
      //   T_R_bin = -1;
      //   W_bin = -1;
      // }
      radfieldbins[modelgridindex][binindex].T_R = T_R_bin;
      radfieldbins[modelgridindex][binindex].W = W_bin;
    }

    /*for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      double J_bin = get_bin_J(modelgridindex,binindex);
      double T_R_bin = get_bin_T_R(modelgridindex,binindex);
      double W_bin = get_bin_W(modelgridindex,binindex);
      printout("bin %4d: J %g, T_R %7.1f, W %12.5e\n",
             binindex, J_bin, T_R_bin, W_bin);
    }*/
    // if (timestep % 2 == 0)
    radfield_write_to_file(modelgridindex, timestep);
  }
}


void radfield_set_J_normfactor(int modelgridindex, double normfactor)
{
  J_normfactor[modelgridindex] = normfactor;
}


void radfield_normalise_J(const int modelgridindex, const double estimator_normfactor_over4pi)
{
  assert(isfinite(J[modelgridindex]));
  J[modelgridindex] *= estimator_normfactor_over4pi;
  for (int i = 0; i < detailed_linecount; i++)
  {
    prev_Jb_lu_normed[modelgridindex][i].value = Jb_lu_raw[modelgridindex][i].value * estimator_normfactor_over4pi;
    prev_Jb_lu_normed[modelgridindex][i].contribcount = Jb_lu_raw[modelgridindex][i].contribcount;
  }
}


void radfield_normalise_bf_estimators(const int modelgridindex, const double estimator_normfactor_over_H)
{
  #if (DETAILED_BF_ESTIMATORS_ON)
  // printout("radfield_normalise_bf_estimators for cell %d with factor %g\n", modelgridindex, estimator_normfactor_over_H);
  for (int i = 0; i < nbfcontinua; i++)
  {
    prev_bfrate_normed[modelgridindex][i] = bfrate_raw[modelgridindex][i] * estimator_normfactor_over_H;
  }

  for (int i = 0; i < nbfcontinua; i++)
  {
    if (prev_bfrate_normed[modelgridindex][i] > 0.)
    {
      normed_bfrates_available = true;
      break;
    }
  }
  #endif
}

double get_bfrate_estimator(const int element, const int lowerion, const int lower, const int phixstargetindex, const int modelgridindex)
{
#if (!DETAILED_BF_ESTIMATORS_ON)
  return -1;
#else
  if (!normed_bfrates_available)
    return -1.;

  // TODO: speed this up with a binary search
  const double nu_edge = get_phixs_threshold(element, lowerion, lower, phixstargetindex);
  for (int i = 0; i < nbfcontinua; i++)
  {
    if ((phixslist[tid].allcont[i].element == element) && (phixslist[tid].allcont[i].ion == lowerion) &&
        (phixslist[tid].allcont[i].level == lower) && (phixslist[tid].allcont[i].phixstargetindex == phixstargetindex))
    {
      return prev_bfrate_normed[modelgridindex][i];
    }

    if (nu_edge > phixslist[tid].allcont[i].nu_edge)
      break;
  }
  printout("no bf rate for element Z=%d ion_stage %d lower %d phixstargetindex %d\n", get_element(element), get_ionstage(element, lowerion), lower, phixstargetindex);
  return 0.;
#endif
}

void radfield_normalise_nuJ(const int modelgridindex, const double estimator_normfactor_over4pi)
{
  assert(isfinite(nuJ[modelgridindex]));
  nuJ[modelgridindex] *= estimator_normfactor_over4pi;
}


double get_T_R_from_J(const int modelgridindex)
{
  const double T_R = pow(PI / STEBO * J[modelgridindex], 1. / 4.);
  if (!isfinite(T_R))
  {
    /// keep old value of T_R
    printout("[warning] get_T_R_from_J: T_R estimator infinite in cell %d, use value of last timestep\n", modelgridindex);
    return get_TR(modelgridindex);
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
void radfield_titer_J(const int modelgridindex)
{
  if (J_reduced_save[modelgridindex] >= 0)
  {
    J[modelgridindex] = (J[modelgridindex] + J_reduced_save[modelgridindex]) / 2;
  }
  J_reduced_save[modelgridindex] = J[modelgridindex];
}


#ifndef FORCE_LTE
void radfield_titer_nuJ(const int modelgridindex)
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
void radfield_reduce_estimators(void)
// reduce and broadcast (allreduce) the estimators for J and nuJ in all bins
{
  MPI_Allreduce(MPI_IN_PLACE, J, npts_model, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #ifndef FORCE_LTE
  MPI_Allreduce(MPI_IN_PLACE, nuJ, npts_model, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif

  #if (DETAILED_BF_ESTIMATORS_ON)
  {
    for (int modelgridindex = 0; modelgridindex < npts_model; modelgridindex++)
    {
      if (get_numassociatedcells(modelgridindex) > 0)
      {
        MPI_Allreduce(MPI_IN_PLACE, bfrate_raw[modelgridindex], nbfcontinua, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      }
    }
  }
  #endif

  if (MULTIBIN_RADFIELD_MODEL_ON)
  {
    const time_t sys_time_start_reduction = time(NULL);
    printout("Reducing binned radiation field estimators");

    for (int modelgridindex = 0; modelgridindex < npts_model; modelgridindex++)
    {
      // printout("DEBUGCELLS: cell %d associated_cells %d\n", modelgridindex, get_numassociatedcells(modelgridindex));
      if (get_numassociatedcells(modelgridindex) > 0)
      {
        // MPI_Barrier(MPI_COMM_WORLD);
        for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
        {
          // printout("MPI: pre-MPI_Allreduce, process %d modelgrid %d binindex %d has a individual contribcount of %d\n",my_rank,modelgridindex,binindex,radfieldbins[modelgridindex][binindex].contribcount);
          MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[modelgridindex][binindex].J_raw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[modelgridindex][binindex].nuJ_raw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[modelgridindex][binindex].contribcount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

          // printout("MPI: After MPI_Allreduce: Process %d modelgrid %d binindex %d has a contribcount of %d\n",my_rank,modelgridindex,binindex,radfieldbins[modelgridindex][binindex].contribcount);
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

    for (int modelgridindex = 0; modelgridindex < npts_model; modelgridindex++)
    {
      // printout("DEBUGCELLS: cell %d associated_cells %d\n", modelgridindex, get_numassociatedcells(modelgridindex));
      if (get_numassociatedcells(modelgridindex) > 0)
      {
        MPI_Barrier(MPI_COMM_WORLD);
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
}


void radfield_MPI_Bcast(const int my_rank, const int root, const int root_nstart, const int root_ndo)
// broadcast computed radfield results including parameters
// from the cells belonging to root process to all processes
{
  // double nu_lower_first;
  if (root_ndo > 0)
  {
    // if (root == my_rank)
    // {
    //   printout("radfield_MPI_Bcast root process %d will send data for cells %d to %d\n", my_rank, root_nstart, root_nstart + root_ndo - 1);
    // }
    // else
    // {
    //   printout("radfield_MPI_Bcast process %d will receive data for cells %d to %d\n", my_rank, root_nstart, root_nstart + root_ndo - 1);
    // }
  }

  for (int modelgridindex = root_nstart; modelgridindex < root_nstart + root_ndo; modelgridindex++)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&J_normfactor[modelgridindex], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
    if (get_numassociatedcells(modelgridindex) > 0)
    {
      if (MULTIBIN_RADFIELD_MODEL_ON)
      {
        MPI_Barrier(MPI_COMM_WORLD);
        for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
        {
          // printout("radfield_MPI_Bcast bin %d T_R before: %g\n", binindex, radfieldbins[modelgridindex][binindex].T_R);
          MPI_Bcast(&radfieldbins[modelgridindex][binindex].W, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
          MPI_Bcast(&radfieldbins[modelgridindex][binindex].T_R, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
          MPI_Bcast(&radfieldbins[modelgridindex][binindex].J_raw, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
          MPI_Bcast(&radfieldbins[modelgridindex][binindex].nuJ_raw, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
          MPI_Bcast(&radfieldbins[modelgridindex][binindex].contribcount, 1, MPI_INT, root, MPI_COMM_WORLD);
          // printout("radfield_MPI_Bcast MPI_Bcast radfield bin %d for cell %d from process %d to %d\n", binindex, modelgridindex, root, my_rank);
          // printout("radfield_MPI_Bcast bin %d T_R after: %g\n", binindex, radfieldbins[modelgridindex][binindex].T_R);
        }
      }

      if (DETAILED_LINE_ESTIMATORS_ON)
      {
        for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++)
        {
          MPI_Bcast(&prev_Jb_lu_normed[modelgridindex][jblueindex].value, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
          MPI_Bcast(&prev_Jb_lu_normed[modelgridindex][jblueindex].contribcount, 1, MPI_INT, root, MPI_COMM_WORLD);
        }
      }
    }
  }
}
#endif


void radfield_write_restart_data(FILE *gridsave_file)
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
    fprintf(gridsave_file, "%d\n", nbfcontinua);

    for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
    {
      if (get_numassociatedcells(modelgridindex) > 0)
      {
        fprintf(gridsave_file, "%d\n", modelgridindex);
        for (int i = 0; i < nbfcontinua; i++)
        {
          fprintf(gridsave_file, "%g ", prev_bfrate_normed[modelgridindex][i]);
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

  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    if (get_numassociatedcells(modelgridindex) > 0)
    {
      fprintf(gridsave_file,"%d %lg\n", modelgridindex, J_normfactor[modelgridindex]);

      if (MULTIBIN_RADFIELD_MODEL_ON)
      {
        for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
        {
          fprintf(gridsave_file, "%lg %lg %g %g %d\n",
                  radfieldbins[modelgridindex][binindex].J_raw,
                  radfieldbins[modelgridindex][binindex].nuJ_raw,
                  radfieldbins[modelgridindex][binindex].W,
                  radfieldbins[modelgridindex][binindex].T_R,
                  radfieldbins[modelgridindex][binindex].contribcount);
                  //radfieldbins[modelgridindex][binindex].fit_type
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


void radfield_read_restart_data(FILE *gridsave_file)
{
  printout("Reading restart data for radiation field\n");

  if (!radfield_initialized)
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
      assert(binindex_in == binindex);
    }
  }

  #if (DETAILED_BF_ESTIMATORS_ON)
  {
    int gridave_nbf_in;
    fscanf(gridsave_file, "%d\n", &gridave_nbf_in);
    assert(gridave_nbf_in == nbfcontinua);

    for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
    {
      if (get_numassociatedcells(modelgridindex) > 0)
      {
        int mgi_in;
        fscanf(gridsave_file, "%d\n", &mgi_in);
        assert(mgi_in == modelgridindex);
        for (int i = 0; i < nbfcontinua; i++)
        {
          fscanf(gridsave_file, "%g ", &prev_bfrate_normed[modelgridindex][i]);
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

  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    if (get_numassociatedcells(modelgridindex) > 0)
    {
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
          fscanf(gridsave_file, "%lg %lg %g %g %d\n",
                 &radfieldbins[modelgridindex][binindex].J_raw,
                 &radfieldbins[modelgridindex][binindex].nuJ_raw,
                 &radfieldbins[modelgridindex][binindex].W,
                 &radfieldbins[modelgridindex][binindex].T_R,
                 &radfieldbins[modelgridindex][binindex].contribcount);
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
int radfield_integrate(
  const gsl_function *f, double nu_a, double nu_b, double epsabs,
  double epsrel, size_t limit, int key, gsl_integration_workspace *workspace,
  double *result, double *abserr)
{
  if (MULTIBIN_RADFIELD_MODEL_ON && (nts_global >= FIRST_NLTE_RADFIELD_TIMESTEP))
  {
    double *pts = malloc((RADFIELDBINCOUNT + 3) * sizeof(double));
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
    //   printout("radfield_integrate singular point number %d at nu %g, (nu_a %g, nu_b %g), radfield_low %g radfield_high %g\n",
    //            e, pts[e], nu_a, nu_b, radfield(pts[e] * 0.9999, 0), radfield(pts[e] * 1.0001, 0));
    // }
    const int status = gsl_integration_qagp(f, pts, npts, epsabs, epsrel, limit, workspace, result, abserr);
    free(pts);
    return status;
  }
  else
    return gsl_integration_qag(f, nu_a, nu_b, epsabs, epsrel, limit, key, workspace, result, abserr);
}