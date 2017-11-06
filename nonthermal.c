#include "assert.h"
#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include "atomic.h"
#include "grid_init.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "update_grid.h"
#include "sn3d.h"

#define SFPTS 8192  // number of energy points in the Spencer-Fano solution vector
#define EMAX 32000. // eV
#define EMIN 1.0 // eV

#define BLOCKSIZEEXCITATION 1000    // Realloc the excitation list increasing by this blocksize
#define BLOCKSIZEIONIZATION 100    // Realloc the ionization list increasing by this blocksize

// THESE OPTIONS ARE USED TO TEST THE SF SOLVER
// Compare to Kozma & Fransson (1992) pure-oxygen plasma, nne = 1e8, x_e = 0.01
// #define yscalefactoroverride(mgi) (1e10)
// #define get_tot_nion(x) (1e10)
// #define ionstagepop(modelgridindex, element, ion) ionstagepop_override(modelgridindex, element, ion)
// #define get_nne(x) (1e8)
// #define SFPTS 10000  // number of energy points in the Spencer-Fano solution vector
// #define EMAX 3000. // eV
// #define EMIN 1. // eV
//
// static double ionstagepop_override(const int modelgridindex, const int element, const int ion)
// Fake the composition to test the NT solver
// {
//   const double nntot = get_tot_nion(modelgridindex);
//   if (get_element(element) == 8)
//   {
//     const int ion_stage = get_ionstage(element, ion);
//     if (ion_stage == 1)
//       return 0.99 * nntot;
//     else if (ion_stage == 2)
//       return 0.01 * nntot;
//   }
//   return 0.;
// }

#define STORE_NT_SPECTRUM false // if this is on, the non-thermal energy spectrum will be kept in memory
                                // for every grid cell during packet propagation, which
                                // can take up a lot of memory for large grid sizes
                                // alternatively, just the non-thermal ionization rates can be stored
                                // but we will probably need to re-enable this option to incorporate
                                // non-thermal excitation rates if there are
                                // many more transitions to store than there are spectrum samples

#define NT_EXCITATION_ON true // if this is on, the non-thermal energy spectrum will be kept in memory

static const double DELTA_E = (EMAX - EMIN) / (SFPTS - 1);

static const double minionfraction = 1.e-4;  // minimum number fraction of the total population to include in SF solution
static const double MINDEPRATE = 0.; // minimum deposition rate density (eV/s/cm^3) to solve SF equation

static const double A_naught_squared = 2.800285203e-17; // Bohr radius squared in cm^2

struct collionrow {
  int Z;
  int nelec;
  int n;
  int l;
  double ionpot_ev;
  double A;
  double B;
  double C;
  double D;
};

static struct collionrow *colliondata = NULL;
static int colliondatacount = 0;

static FILE *restrict nonthermalfile = NULL;
static bool nonthermal_initialized = false;

static gsl_vector *envec;            // energy grid on which solution is sampled
static gsl_vector *sourcevec;        // samples of the source function (energy distribution of deposited energy)
static double E_init_ev = 0;         // the energy injection rate density (and mean energy of injected electrons if source integral is one) in eV

struct nt_ionization_struct
{
  double frac_deposition;  // the fraction of the non-thermal deposition energy going to ionizing this ion
  int element;
  int ion;
};

struct nt_excitation_struct
{
  double frac_deposition;  // the fraction of the non-thermal deposition energy going to the excitation transition
  double ratecoeffperdeposition; // the excitation rate coefficient divided by the deposition rate density
  int lineindex;
};

struct nt_solution_struct {
  double E_0;     // the lowest energy ionization or excitation transition in eV
  double *yfunc;  // Samples of the Spencer-Fano solution function. Multiply by energy to get non-thermal electron number flux.
                  // y(E) * dE is the flux of electrons with energy in the range (E, E + dE)

  double deposition_rate_density;

  float frac_heating;              // energy fractions should add up to 1.0 if the solution is good
  float frac_ionization;           // fraction of deposition energy going to ionization
  float frac_excitation;           // fraction of deposition energy going to excitation

  float eff_ionpot[MELEMENTS][MIONS]; // these are used to calculate the non-thermal ionization rate

  int frac_ionizations_list_size;
  struct nt_ionization_struct *frac_ionizations_list;

  int frac_excitations_list_size;
  struct nt_excitation_struct *frac_excitations_list;

  int timestep;                 // the quantities above were calculated for this timestep
};

static struct nt_solution_struct nt_solution[MMODELGRID+1];


static int compare_ionization_fractions(const void *p1, const void *p2)
{
  const struct nt_ionization_struct *elem1 = p1;
  const struct nt_ionization_struct *elem2 = p2;

 if (elem1->frac_deposition < elem2->frac_deposition)
    return 1;
 else if (elem1->frac_deposition > elem2->frac_deposition)
    return -1;
 else
    return 0;
}


static int compare_excitation_fractions(const void *p1, const void *p2)
{
  const struct nt_excitation_struct *elem1 = p1;
  const struct nt_excitation_struct *elem2 = p2;

 if (elem1->frac_deposition < elem2->frac_deposition)
    return 1;
 else if (elem1->frac_deposition > elem2->frac_deposition)
    return -1;
 else
    return 0;
}


#ifndef get_tot_nion
static double get_tot_nion(const int modelgridindex)
{
  double result = 0.;
  for (int element = 0; element < nelements; element++)
  {
    result += modelgrid[modelgridindex].composition[element].abundance / elements[element].mass * get_rho(modelgridindex);

    // alternative method is to add the ion populations
    // const int nions = get_nions(element);
    // for (ion = 0; ion < nions; ion++)
    // {
    //    result += ionstagepop(modelgridindex,element,ion);
    // }
  }

  return result;
}
#endif


static void read_collion_data(void)
{
  printout("Reading collisional ionization data...\n");

  FILE *cifile = fopen("collion.txt", "r");
  if (cifile == NULL)
  {
    printout("Could not open collion.txt\n");
    abort();
  }

  fscanf(cifile, "%d", &colliondatacount);
  printout("Reading %d collisional transition rows\n", colliondatacount);
  colliondata = calloc(colliondatacount, sizeof(struct collionrow));
  for (int n = 0; n < colliondatacount; n++)
  {
    fscanf(cifile, "%2d %2d %1d %1d %lg %lg %lg %lg %lg",
           &colliondata[n].Z, &colliondata[n].nelec, &colliondata[n].n, &colliondata[n].l,
           &colliondata[n].ionpot_ev, &colliondata[n].A, &colliondata[n].B, &colliondata[n].C, &colliondata[n].D);
    // printout("ci row: %2d %2d %1d %1d %lg %lg %lg %lg %lg\n",
    //          colliondata[n].Z, colliondata[n].nelec, colliondata[n].n, colliondata[n].l, colliondata[n].ionpot_ev,
    //          colliondata[n].A, colliondata[n].B, colliondata[n].C, colliondata[n].D);
  }

  fclose(cifile);
}


static void zero_all_effionpot(const int modelgridindex)
{
  for (int element = 0; element < nelements; element++)
  {
    for (int ion = 0; ion < get_nions(element); ion++)
    {
      nt_solution[modelgridindex].eff_ionpot[element][ion] = 0.;
    }
  }
}

void nt_init(const int my_rank)
{
  if (nonthermal_initialized == false)
  {
    printout("Initializing non-thermal solver\n");
    char filename[100];
    sprintf(filename,"nonthermalspec_%.4d.out", my_rank);
    nonthermalfile = fopen(filename, "w");
    if (nonthermalfile == NULL)
    {
      printout("Cannot open %s.\n",filename);
      abort();
    }
    fprintf(nonthermalfile,"%8s %15s %8s %11s %11s %11s\n",
            "timestep","modelgridindex","index","energy_ev","source","y");
    fflush(nonthermalfile);

    for (int modelgridindex = 0; modelgridindex < MMODELGRID + 1; modelgridindex++)
    {
      // should make these negative?
      nt_solution[modelgridindex].frac_heating = 0.98;
      nt_solution[modelgridindex].frac_ionization = 0.02;
      nt_solution[modelgridindex].frac_excitation = 0.0;

      nt_solution[modelgridindex].timestep = -1;
      nt_solution[modelgridindex].E_0 = 0.;
      nt_solution[modelgridindex].deposition_rate_density = -1.;

      if (STORE_NT_SPECTRUM && mg_associated_cells[modelgridindex] > 0)
      {
        nt_solution[modelgridindex].yfunc = calloc(SFPTS, sizeof(double));
      }
      else
      {
        nt_solution[modelgridindex].yfunc = NULL;
      }

      nt_solution[modelgridindex].frac_ionizations_list = NULL;
      nt_solution[modelgridindex].frac_ionizations_list_size = 0;

      nt_solution[modelgridindex].frac_excitations_list = NULL;
      nt_solution[modelgridindex].frac_excitations_list_size = 0;

      zero_all_effionpot(modelgridindex);
    }

    envec = gsl_vector_calloc(SFPTS); // energy grid on which solution is sampled
    sourcevec = gsl_vector_calloc(SFPTS); // energy grid on which solution is sampled

    // const int source_spread_pts = ceil(SFPTS / 20);
    const int source_spread_pts = ceil(SFPTS * 0.03333); // KF92 OXYGEN TEST
    for (int s = 0; s < SFPTS; s++)
    {
      const double energy_ev = EMIN + s * DELTA_E;

      gsl_vector_set(envec, s, energy_ev);

      // spread the source over some energy width
      if (s < SFPTS - source_spread_pts)
        gsl_vector_set(sourcevec, s, 0.);
      else if (s < SFPTS)
        gsl_vector_set(sourcevec, s, 1. / (DELTA_E * source_spread_pts));
    }

    double envec_dot_sourcevec;
    gsl_blas_ddot(envec, sourcevec, &envec_dot_sourcevec);
    E_init_ev = envec_dot_sourcevec * DELTA_E;

    // or put all of the source into one point at EMAX
    // gsl_vector_set_zero(sourcevec);
    // gsl_vector_set(sourcevec, SFPTS-1, 1 / DELTA_E);
    // E_init_ev = EMAX;

    printout("E_init: %14.7e eV\n", E_init_ev);
    double sourceintegral = gsl_blas_dasum(sourcevec) * DELTA_E;
    printout("source vector integral: %14.7e\n", sourceintegral);

    read_collion_data();

    nonthermal_initialized = true;
    printout("Finished initializing non-thermal solver\n");
  }
  else
    printout("Tried to initialize the non-thermal solver more than once!\n");
}


void calculate_deposition_rate_density(const int modelgridindex, const int timestep)
// should be in erg / s / cm^3
{
  const double gamma_deposition = rpkt_emiss[modelgridindex] * 1.e20 * FOURPI;
  // Above is the gamma-ray bit. Below is *supposed* to be the kinetic energy of positrons created by 56Co and 48V. These formulae should be checked, however.

  const double t = time_step[timestep].mid;
  const double rho = get_rho(modelgridindex);

  const double co56_positron_dep = (0.610 * 0.19 * MEV) *
        (exp(-t / TCOBALT) - exp(-t / TNICKEL)) /
        (TCOBALT - TNICKEL) * get_f56ni(modelgridindex) * rho / MNI56;

  const double v48_positron_dep = (0.290 * 0.499 * MEV) *
        (exp(-t / T48V) - exp(-t / T48CR)) /
        (T48V - T48CR) * get_f48cr(modelgridindex) * rho / MCR48;

  //printout("nt_deposition_rate: element: %d, ion %d\n",element,ion);
  //printout("nt_deposition_rate: gammadep: %g, poscobalt %g pos48v %g\n",
  //         gamma_deposition,co56_positron_dep,v48_positron_dep);

  nt_solution[modelgridindex].deposition_rate_density = gamma_deposition + co56_positron_dep + v48_positron_dep;
  nt_solution[modelgridindex].timestep = timestep;
}


double get_deposition_rate_density(const int modelgridindex)
// should be in erg / s / cm^3
{
  // if (nt_solution[modelgridindex].deposition_rate_density <= 0)
  // {
  //   calculate_deposition_rate_density(modelgridindex, nts_global);
  //   printout("No deposition_rate_density for cell %d. Calculated value of %g has been stored.\n",
  //            modelgridindex, nt_solution[modelgridindex].deposition_rate_density);
  // }
  assert(nt_solution[modelgridindex].timestep == nts_global);
  assert(nt_solution[modelgridindex].deposition_rate_density >= 0);
  return nt_solution[modelgridindex].deposition_rate_density;
}


static double get_y_sample(const int modelgridindex, const int index)
{
  if (nt_solution[modelgridindex].yfunc != NULL)
  {
    return nt_solution[modelgridindex].yfunc[index];
  }
  else
  {
    printout("non-thermal: attempted to get y function sample index %d in cell %d, but the y array pointer is null\n",
             index, modelgridindex);
    abort();
  }
}


static void nt_write_to_file(const int modelgridindex, const int timestep)
{
# ifdef _OPENMP
# pragma omp critical (nonthermal_out_file)
  {
# endif
  if (!nonthermal_initialized)
  {
    printout("Call to nonthermal_write_to_file before nonthermal_init");
    abort();
  }

#ifndef yscalefactoroverride // manual override can be defined
  const double yscalefactor = (get_deposition_rate_density(modelgridindex) / (E_init_ev * EV));
#else
  const double yscalefactor = yscalefactoroverride(modelgridindex);
#endif

  for (int s = 0; s < SFPTS; s++)
  {
    fprintf(nonthermalfile,"%8d %15d %8d %11.5e %11.5e %11.5e\n",
            timestep, modelgridindex, s, gsl_vector_get(envec,s),
            gsl_vector_get(sourcevec,s), yscalefactor * get_y_sample(modelgridindex, s));
  }
  fflush(nonthermalfile);
# ifdef _OPENMP
  }
# endif
}


void nt_close_file(void)
{
  fclose(nonthermalfile);
  gsl_vector_free(envec);
  gsl_vector_free(sourcevec);
  if (STORE_NT_SPECTRUM)
  {
    for (int modelgridindex = 0; modelgridindex < MMODELGRID + 1; modelgridindex++)
    {
      if (mg_associated_cells[modelgridindex] > 0)
      {
        free(nt_solution[modelgridindex].yfunc);
        free(nt_solution[modelgridindex].frac_ionizations_list);
        if (nt_solution[modelgridindex].frac_excitations_list_size > 0)
          free(nt_solution[modelgridindex].frac_excitations_list);
      }
    }
  }
  free(colliondata);
  nonthermal_initialized = false;
}


static int get_energyindex_ev(const double energy_ev)
// finds the nearest energy point to energy_ev (may be above or below)
{
  int index = floor((energy_ev - EMIN) / DELTA_E);
  if (index < 0)
    return 0;
  else if (index > SFPTS - 1)
    return SFPTS - 1;
  else
    return index;
}


static int get_y(const int modelgridindex, const double energy_ev)
{
  const int index = (energy_ev - EMIN) / DELTA_E;
  // assert(index > 0);
  if (index < 0)
  {
    return 0.;
  }
  else if (index > SFPTS - 1)
    return 0.;
  else
  {
    const double enbelow = gsl_vector_get(envec, index);
    const double ybelow = get_y_sample(modelgridindex, index);
    const double yabove = get_y_sample(modelgridindex, index + 1);
    const double x = (energy_ev - enbelow) / DELTA_E;
    return (1 - x) * ybelow + x * yabove;
  }
}


static double electron_loss_rate(const double energy, const double nne)
// -dE / dx for fast electrons
// energy is in ergs
// return units are erg / cm
{
  const double omegap = sqrt(4 * PI * nne * pow(QE, 2) / ME);
  const double zetae = H * omegap / 2 / PI;
  const double v = sqrt(2 * energy / ME);
  if (energy > 14 * EV)
  {
    return nne * 2 * PI * pow(QE, 4) / energy * log(2 * energy / zetae);
  }
  else
  {
    const double eulergamma = 0.577215664901532;
    return nne * 2 * PI * pow(QE, 4) / energy * log(ME * pow(v, 3) / (eulergamma * pow(QE, 2) * omegap));
  }
}


static double xs_excitation(const int lineindex, const double epsilon_trans, const double energy)
// collisional excitation cross section in cm^2
// energies are in erg
{
  if (energy < epsilon_trans)
    return 0.;
  const double coll_str = get_coll_str(lineindex);

  if (coll_str >= 0)
  {
    // collision strength is available, so use it
    // Li et al. 2012 equation 11
    return pow(H_ionpot / energy, 2) / statw_lower(lineindex) * coll_str * PI * A_naught_squared;
  }
  else if (!linelist[lineindex].forbidden)
  {
    const double fij = osc_strength(lineindex);
    // permitted E1 electric dipole transitions
    const double U = energy / epsilon_trans;

    // const double g_bar = 0.2;
    const double A = 0.28;
    const double B = 0.15;
    const double g_bar = A * log(U) + B;

    const double prefactor = 45.585750051; // 8 * pi^2/sqrt(3)
    // Eq 4 of Mewe 1972, possibly from Seaton 1962?
    return prefactor * A_naught_squared * pow(H_ionpot / epsilon_trans, 2) * fij * g_bar / U;
  }
  else
    return 0.;
}


static bool get_xs_excitation_vector(gsl_vector *xs_excitation_vec, const int lineindex, const double epsilon_trans)
// vector of collisional excitation cross sections in cm^2
// epsilon_trans is in erg
// returns true if any vector components are (or might be) non-zero
// if it returns false, all vector components are definitely zero
{
  const double coll_str = get_coll_str(lineindex);
  const int en_startindex = 0; // energy point corresponding to epsilon_trans
  // const int en_startindex = get_energyindex_ev(epsilon_trans_ev);

  // make sure every value of xs_excitation_vec is set
  // for (int j = 0; j < en_startindex; j++)
  //   gsl_vector_set(xs_excitation_vec, j, 0.);

  bool hasnonzerovalue = false;
  if (coll_str >= 0)
  {
    // collision strength is available, so use it
    // Li et al. 2012 equation 11
    const double constantfactor = pow(H_ionpot, 2) / statw_lower(lineindex) * coll_str * PI * A_naught_squared;
    for (int j = en_startindex; j < SFPTS; j++)
    {
      const double energy = gsl_vector_get(envec, j) * EV;
      if (energy >= epsilon_trans)
      {
        gsl_vector_set(xs_excitation_vec, j, constantfactor * pow(energy, -2));
        hasnonzerovalue = true;
      }
      else
        gsl_vector_set(xs_excitation_vec, j, 0.);
    }
  }
  else if (!linelist[lineindex].forbidden)
  {
    const double fij = osc_strength(lineindex);
    // permitted E1 electric dipole transitions

    // const double g_bar = 0.2;
    const double A = 0.28;
    const double B = 0.15;

    const double prefactor = 45.585750051; // 8 * pi^2/sqrt(3)
    // Eq 4 of Mewe 1972, possibly from Seaton 1962?
    const double constantfactor = prefactor * A_naught_squared * pow(H_ionpot / epsilon_trans, 2) * fij;
    for (int j = en_startindex; j < SFPTS; j++)
    {
      const double energy = gsl_vector_get(envec, j) * EV;
      if (energy >= epsilon_trans)
      {
        const double U = energy / epsilon_trans;
        const double g_bar = A * log(U) + B;
        gsl_vector_set(xs_excitation_vec, j, constantfactor * g_bar / U);
        hasnonzerovalue = true;
      }
      else
        gsl_vector_set(xs_excitation_vec, j, 0.);
    }
  }
  else
  {
    gsl_vector_set_zero(xs_excitation_vec);
  }

  return hasnonzerovalue;
}


static double xs_impactionization(const double energy_ev, const int collionindex)
// impact ionization cross section in cm^2
// energy and ionization_potential should be in eV
// fitting forumula of Younger 1981
// called Q_i(E) in KF92 equation 7
{
  const double ionpot_ev = colliondata[collionindex].ionpot_ev;
  const double u = energy_ev / ionpot_ev;

  if (u <= 1.)
  {
    return 0;
  }
  else
  {
    const double A = colliondata[collionindex].A;
    const double B = colliondata[collionindex].B;
    const double C = colliondata[collionindex].C;
    const double D = colliondata[collionindex].D;

    return 1e-14 * (A * (1 - 1/u) + B * pow((1 - 1/u), 2) + C * log(u) + D * log(u) / u) / (u * pow(ionpot_ev, 2));
  }
}


static void get_xs_ionization_vector(gsl_vector *xs_vec, const int collionindex)
// xs_vec will be set with impact ionization cross sections for E > ionpot_ev (and zeros below this energy)
{
  const double ionpot_ev = colliondata[collionindex].ionpot_ev;
  const int startindex = get_energyindex_ev(ionpot_ev);

  for (int i = 0; i < startindex; i++)
    gsl_vector_set(xs_vec, i, 0.);

  for (int i = startindex; i < SFPTS; i++)
  {
    const double endash = gsl_vector_get(envec, i);
    gsl_vector_set(xs_vec, i, xs_impactionization(endash, collionindex));
  }
}


static double Psecondary(const double e_p, const double epsilon, const double I, const double J)
// distribution of secondary electron energies for primary electron with energy e_p
// Opal, Peterson, & Beaty (1971)
{
  const double e_s = epsilon - I;
  return 1 / (J * atan((e_p - I) / 2 / J) * (1 + pow(e_s / J, 2)));
}


static double get_J(const int Z, const int ionstage, const double ionpot_ev)
{
  // returns an energy in eV
  // values from Opal et al. 1971 as applied by Kozma & Fransson 1992
  if (ionstage == 1)
  {
    if (Z == 2) // He I
      return 15.8;
    else if (Z == 10) // Ne I
      return 24.2;
    else if (Z == 18) // Ar I
      return 10.0;
  }

  return 0.6 * ionpot_ev;
}


static double N_e(const int modelgridindex, const double energy)
// Kozma & Fransson equation 6.
// Something related to a number of electrons, needed to calculate the heating fraction in equation 3
// possibly not valid for energy > E_0
{
  const double energy_ev = energy / EV;
  const double tot_nion = get_tot_nion(modelgridindex);
  double N_e = 0.;

  for (int element = 0; element < nelements; element++)
  {
    const int Z = get_element(element);
    const int nions = get_nions(element);

    for (int ion = 0; ion < nions; ion++)
    {
      double N_e_ion = 0.;
      const int ionstage = get_ionstage(element, ion);
      const double nnion = ionstagepop(modelgridindex, element, ion);

      if (nnion < minionfraction * tot_nion) // skip negligible ions
        continue;

      // excitation terms
      // for (int level = 0; level < get_nlevels(element,ion); level++)
      const int level = 0; // just consider excitation from the ground level
      {
        const int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
        for (int t = 1; t <= nuptrans; t++)
        {
          const double epsilon_trans = elements[element].ions[ion].levels[level].uptrans[t].epsilon_trans;
          const int lineindex = elements[element].ions[ion].levels[level].uptrans[t].lineindex;
          const double epsilon_trans_ev = epsilon_trans / EV;
          N_e_ion += get_y(modelgridindex, energy_ev + epsilon_trans_ev) * xs_excitation(lineindex, epsilon_trans, energy + epsilon_trans);
        }
      }

      // ionization terms
      for (int n = 0; n < colliondatacount; n++)
      {
        if (colliondata[n].Z == Z && colliondata[n].nelec == Z - ionstage + 1)
        {
          const double ionpot_ev = colliondata[n].ionpot_ev;
          const double J = get_J(Z, ionstage, ionpot_ev);
          const double lambda = fmin(EMAX - energy_ev, energy_ev + ionpot_ev);

          const int integral1startindex = get_energyindex_ev(ionpot_ev);
          const int integral1stopindex = get_energyindex_ev(lambda);
          const int integral2startindex = get_energyindex_ev(2 * energy_ev + ionpot_ev);

          for (int i = 0; i < SFPTS; i++)
          {
            double endash = gsl_vector_get(envec, i);

            if (i >= integral1startindex && i <= integral1stopindex) // integral from ionpot up to lambda
            {
              N_e_ion += get_y(modelgridindex, energy_ev + endash) * xs_impactionization(energy_ev + endash, n) * Psecondary(energy_ev + endash, endash, ionpot_ev, J) * DELTA_E;
            }

            if (i >= integral2startindex) // integral from 2E + I up to E_max
            {
              N_e_ion += get_y_sample(modelgridindex, i) * xs_impactionization(endash, n) * Psecondary(endash, energy_ev + ionpot_ev, ionpot_ev, J) * DELTA_E;
            }
          }

        }
      }

      N_e += nnion * N_e_ion;
    }
  }

  // source term
  N_e += gsl_vector_get(sourcevec, get_energyindex_ev(energy_ev));

  return N_e;
}


static float calculate_frac_heating(const int modelgridindex)
// Kozma & Fransson equation 3
{
  // frac_heating multiplied by E_init, which will be divided out at the end
  double frac_heating_Einit = 0.;

  const float nne = get_nne(modelgridindex);
  const double E_0 = nt_solution[modelgridindex].E_0;

  const int startindex = get_energyindex_ev(E_0);
  for (int i = startindex; i < SFPTS; i++)
  {
    const double endash = gsl_vector_get(envec, i);
    double deltaendash;
    if (i == startindex)
      deltaendash = endash + DELTA_E - E_0;
    else
      deltaendash = DELTA_E;
    // first term
    frac_heating_Einit += get_y_sample(modelgridindex, i) * (electron_loss_rate(endash * EV, nne) / EV) * deltaendash;
  }

  // second term
  frac_heating_Einit += E_0 * get_y(modelgridindex, E_0) * (electron_loss_rate(E_0 * EV, nne) / EV);

  // third term (integral from zero to E_0)
  const int nsteps = 100;
  const double delta_endash = E_0 / nsteps;
  for (int j = 0; j < nsteps; j++)
  {
    const double endash = E_0 * j / nsteps;
    frac_heating_Einit += N_e(modelgridindex, endash * EV) * endash * delta_endash;
  }

  const float frac_heating = frac_heating_Einit / E_init_ev;

  if (!isfinite(frac_heating) || frac_heating < 0 || frac_heating > 1.0)
  {
    printout("WARNING: calculate_frac_heating: invalid result of %g. Setting to 1.0 instead\n");
    return 1.0;
  }

  return frac_heating;
}


float get_nt_frac_heating(const int modelgridindex)
{
  if (!NT_SOLVE_SPENCERFANO)
    return 0.98;

  const float frac_heating = nt_solution[modelgridindex].frac_heating;
  // add any debugging checks here?
  return frac_heating;
}


float get_nt_frac_ionization(const int modelgridindex)
{
  if (!NT_SOLVE_SPENCERFANO)
    return 0.02;

  const float frac_ionization = nt_solution[modelgridindex].frac_ionization;

  if (frac_ionization < 0 || !isfinite(frac_ionization))
  {
    printout("ERROR: get_nt_frac_ionization called with no valid solution stored for cell %d. frac_ionization = %g\n",
             modelgridindex, frac_ionization);
    abort();
  }

  return frac_ionization;
}


float get_nt_frac_excitation(const int modelgridindex)
{
  if (!NT_SOLVE_SPENCERFANO)
    return 0.;

  const float frac_excitation = nt_solution[modelgridindex].frac_excitation;

  if (frac_excitation < 0 || !isfinite(frac_excitation))
  {
    printout("ERROR: get_nt_frac_excitation called with no valid solution stored for cell %d. frac_excitation = %g\n",
             modelgridindex, frac_excitation);
    abort();
  }

  return frac_excitation;
}


static double get_mean_binding_energy(const int element, const int ion)
{
  int q[M_NT_SHELLS];
  double total;

  const int ioncharge = get_ionstage(element,ion) - 1;
  const int nbound = get_element(element) - ioncharge; //number of bound electrons

  if (nbound > 0)
  {
    for (int i = 0; i < M_NT_SHELLS; i++)
    {
      q[i] = 0;
    }

    for (int electron_loop = 0; electron_loop < nbound; electron_loop++)
    {
      if (q[0] < 2) //K 1s
      {
        q[0]++;
      }
      else if (q[1] < 2) //L1 2s
      {
        q[1]++;
      }
      else if (q[2] < 2) //L2 2p[1/2]
      {
        q[2]++;
      }
      else if (q[3] < 4) //L3 2p[3/2]
      {
        q[3]++;
      }
      else if (q[4] < 2) //M1 3s
      {
        q[4]++;
      }
      else if (q[5] < 2) //M2 3p[1/2]
      {
        q[5]++;
      }
      else if (q[6] < 4) //M3 3p[3/2]
      {
        q[6]++;
      }
      else if (ioncharge == 0)
      {
        if (q[9] < 2) //N1 4s
        {
          q[9]++;
        }
        else if (q[7] < 4) //M4 3d[3/2]
        {
          q[7]++;
        }
        else if (q[8] < 6) //M5 3d[5/2]
        {
          q[8]++;
        }
        else
        {
          printout("Going beyond the 4s shell in NT calculation. Abort!\n");
          abort();
        }
      }
      else if (ioncharge == 1)
      {
        if (q[9] < 1) // N1 4s
        {
          q[9]++;
        }
        else if (q[7] < 4) //M4 3d[3/2]
        {
          q[7]++;
        }
        else if (q[8] < 6) //M5 3d[5/2]
        {
          q[8]++;
        }
        else
        {
          printout("Going beyond the 4s shell in NT calculation. Abort!\n");
          abort();
        }
      }
      else if (ioncharge > 1)
      {
        if (q[7] < 4) //M4 3d[3/2]
        {
          q[7]++;
        }
        else if (q[8] < 6) //M5 3d[5/2]
        {
          q[8]++;
        }
        else
        {
          printout("Going beyond the 4s shell in NT calculation. Abort!\n");
          abort();
        }
      }
    }

    //      printout("For element %d ion %d I got q's of: %d %d %d %d %d %d %d %d %d %d\n", element, ion, q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7], q[8], q[9]);
    //printout("%g %g %g %g %g %g %g %g %g %g\n", electron_binding[get_element(element)-1][0], electron_binding[get_element(element)-1][1], electron_binding[get_element(element)-1][2],electron_binding[get_element(element)-1][3],electron_binding[get_element(element)-1][4],electron_binding[get_element(element)-1][5],electron_binding[get_element(element)-1][6],electron_binding[get_element(element)-1][7],electron_binding[get_element(element)-1][8],electron_binding[get_element(element)-1][9]);

    total = 0.0;
    for (int electron_loop = 0; electron_loop < M_NT_SHELLS; electron_loop++)
    {
      const double electronsinshell = q[electron_loop];
      if ((electronsinshell) > 0)
      {
        double use2 = electron_binding[get_element(element) - 1][electron_loop];
        const double use3 = elements[element].ions[ion].ionpot;
        if (use2 <= 0)
        {
          use2 = electron_binding[get_element(element) - 1][electron_loop-1];
          //  to get total += electronsinshell/electron_binding[get_element(element)-1][electron_loop-1];
          //  set use3 = 0.
          if (electron_loop != 8)
          {
            //For some reason in the Lotz data, this is no energy for the M5 shell before Ni. So if the complaint
            //is for 8 (corresponding to that shell) then just use the M4 value
            printout("Huh? I'm trying to use a binding energy when I have no data. element %d ion %d\n",element,ion);
            printout("Z = %d, ion_stage = %d\n", get_element(element), get_ionstage(element, ion));
            abort();
          }
        }
        if (use2 < use3)
        {
          total += electronsinshell / use3;
        }
        else
        {
          total += electronsinshell / use2;
        }
      }
      //printout("total %g\n", total);
    }

  }
  else
  {
    total = 0.0;
  }

  //printout("For element %d ion %d I got mean binding energy of %g (eV)\n", element, ion, 1./total/EV);

  return total;
}


static double get_oneoverw(const int element, const int ion, const int modelgridindex)
{
  // Routine to compute the work per ion pair for doing the NT ionization calculation.
  // Makes use of EXTREMELY SIMPLE approximations - high energy limits only

  // Work in terms of 1/W since this is actually what we want. It is given by sigma/(Latom + Lelec).
  // We are going to start by taking all the high energy limits and ignoring Lelec, so that the
  // denominator is extremely simplified. Need to get the mean Z value.

  double Zbar = 0.0;  // mass-weighted average atomic number
  for (int ielement = 0; ielement < nelements; ielement++)
  {
    Zbar += modelgrid[modelgridindex].composition[ielement].abundance * elements[ielement].anumber;
  }
  //printout("cell %d has Zbar of %g\n", modelgridindex, Zbar);

  const double Aconst = 1.33e-14 * EV * EV;
  const double binding = get_mean_binding_energy(element, ion);
  const double oneoverW = Aconst * binding / Zbar / (2 * 3.14159 * pow(QE, 4));
  //printout("For element %d ion %d I got W of %g (eV)\n", element, ion, 1./oneoverW/EV);

  return oneoverW;
}


static double calculate_nt_frac_excitation_ion(const int modelgridindex, const int element, const int ion)
// Kozma & Fransson equation 4, but summed over all transitions for given ion
// integral in Kozma & Fransson equation 9
{
  gsl_vector *xs_excitation_vec_sum_alltrans = gsl_vector_calloc(SFPTS);
  gsl_vector *xs_excitation_nnlevel_epsilontrans_vec = gsl_vector_calloc(SFPTS);

  bool hasnonzerovalue = false;

  const int maxlevel = 0; // just consider excitation from the ground level
  // const int maxlevel = get_nlevels(element, ion); // excitation from all levels (SLOW)

  for (int level = 0; level <= maxlevel; level++)
  {
    const int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
    const int nnlevel = calculate_exclevelpop(modelgridindex, element, ion, level);

    for (int t = 1; t <= nuptrans; t++)
    {
      const double epsilon_trans = elements[element].ions[ion].levels[level].uptrans[t].epsilon_trans;
      const int lineindex = elements[element].ions[ion].levels[level].uptrans[t].lineindex;

      if (get_xs_excitation_vector(xs_excitation_nnlevel_epsilontrans_vec, lineindex, epsilon_trans))
      {
        hasnonzerovalue = true;
        gsl_blas_daxpy(nnlevel * epsilon_trans / EV, xs_excitation_nnlevel_epsilontrans_vec, xs_excitation_vec_sum_alltrans);
      }
    }
  }

  gsl_vector_free(xs_excitation_nnlevel_epsilontrans_vec);

  if (hasnonzerovalue)
  {
    double y_dot_crosssection = 0.;
    gsl_vector_view yvecview = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);
    gsl_blas_ddot(&yvecview.vector, xs_excitation_vec_sum_alltrans, &y_dot_crosssection);
    gsl_vector_free(xs_excitation_vec_sum_alltrans);

    return y_dot_crosssection * DELTA_E / E_init_ev;
  }
  else
  {
    gsl_vector_free(xs_excitation_vec_sum_alltrans);

    return 0.;
  }
}


static double calculate_nt_frac_ionization_shell(const int modelgridindex, const int element, const int ion, const int collionindex)
// the fraction of deposition energy that goes into ionising electrons in this particular shell
{
  const double nnion = ionstagepop(modelgridindex, element, ion); // hopefully ions per cm^3?
  const double ionpot_ev = colliondata[collionindex].ionpot_ev;

  gsl_vector *cross_section_vec = gsl_vector_alloc(SFPTS);
  get_xs_ionization_vector(cross_section_vec, collionindex);

  gsl_vector_view yvecview = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);

  double y_dot_crosssection = 0.;
  gsl_blas_ddot(&yvecview.vector, cross_section_vec, &y_dot_crosssection);
  gsl_vector_free(cross_section_vec);

  return nnion * ionpot_ev * y_dot_crosssection * DELTA_E / E_init_ev;
}


static double nt_ionization_ratecoeff_wfapprox(const int modelgridindex, const int element, const int ion)
// non-thermal ionization rate coefficient (multiply by population to get rate)
{
  const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
  // to get the non-thermal ionization rate we need to divide the energy deposited
  // per unit volume per unit time in the grid cell (sum of terms above)
  // by the total ion number density and the "work per ion pair"
  return deposition_rate_density / get_tot_nion(modelgridindex) * get_oneoverw(element, ion, modelgridindex);
}


static double calculate_nt_ionization_ratecoeff(const int modelgridindex, const int element, const int ion)
// Integrate the ionization cross section over the electron degradation function to get the ionization rate coefficient
// i.e. multiply this by ion population to get a rate of ionizations per second
// Do not call during packet propagation, as the y vector may not be in memory!
// IMPORTANT: we are dividing by the shell potential, not the valence potential here!
// To change this, include a factor of (ionpot_shell / ionpot_valence)
{
  gsl_vector *cross_section_vec = gsl_vector_alloc(SFPTS);
  gsl_vector *cross_section_vec_allshells = gsl_vector_calloc(SFPTS);

  const int Z = get_element(element);
  const int ionstage = get_ionstage(element, ion);

  for (int collionindex = 0; collionindex < colliondatacount; collionindex++)
  {
    if (colliondata[collionindex].Z == Z && colliondata[collionindex].nelec == Z - ionstage + 1)
    {
      get_xs_ionization_vector(cross_section_vec, collionindex);
      gsl_vector_add(cross_section_vec_allshells, cross_section_vec);
    }
  }

  gsl_vector_free(cross_section_vec);

  double y_dot_crosssection = 0.;
  gsl_vector_view yvecview_thismgi = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);
  gsl_blas_ddot(&yvecview_thismgi.vector, cross_section_vec_allshells, &y_dot_crosssection);

  gsl_vector_free(cross_section_vec_allshells);

  const double deposition_rate_density_ev = get_deposition_rate_density(modelgridindex) / EV;
  const double yscalefactor = deposition_rate_density_ev / E_init_ev;

  return yscalefactor * y_dot_crosssection * DELTA_E;
}


static float calculate_eff_ionpot(const int modelgridindex, const int element, const int ion)
// Kozma & Fransson 1992 equation 12, except modified to be a sum over all shells of an ion
// the result is in ergs
{
  const int Z = get_element(element);
  const int ionstage = get_ionstage(element, ion);
  const double nnion = ionstagepop(modelgridindex, element, ion); // ions/cm^3
  const double tot_nion = get_tot_nion(modelgridindex);
  const double X_ion = nnion / tot_nion; // molar fraction of this ion

  // The ionization rates of all shells of an ion add to make the ion's total ionization rate,
  // i.e., Gamma_ion = Gamma_shell_a + Gamma_shell_b + ...
  // And since the ionization rate is inversely proportional to the effective ion potential,
  // we solve:
  // (eta_ion / ionpot_ion) = (eta_shell_a / ionpot_shell_a) + (eta_shell_b / ionpot_shell_b) + ...
  // where eta is the fraction of the deposition energy going into ionization of the ion or shell

  double eta_over_ionpot_sum = 0.;
  double ionpot_valence = -1;
  for (int collionindex = 0; collionindex < colliondatacount; collionindex++)
  {
    if (colliondata[collionindex].Z == Z && colliondata[collionindex].nelec == Z - ionstage + 1)
    {
      const double frac_ionization_shell = calculate_nt_frac_ionization_shell(modelgridindex, element, ion, collionindex);
      const double ionpot_shell = colliondata[collionindex].ionpot_ev * EV;
      if (ionpot_valence < 0)
        ionpot_valence = ionpot_shell;
      assert(ionpot_shell >= ionpot_valence);

      eta_over_ionpot_sum += frac_ionization_shell / ionpot_valence;
    }
  }

  return X_ion / eta_over_ionpot_sum;
}


static float get_eff_ionpot(const int modelgridindex, const int element, int const ion)
// get the effective ion potential from the stored value
// a value of 0. should be treated as invalid
{
  return nt_solution[modelgridindex].eff_ionpot[element][ion];
  // OR
  // return calculate_eff_ionpot(modelgridindex, element, ion);
}


static double nt_ionization_ratecoeff_sf(const int modelgridindex, const int element, const int ion)
// Kozma & Fransson 1992 equation 13
{
  if (mg_associated_cells[modelgridindex] <= 0)
  {
    printout("ERROR: nt_ionization_ratecoeff_sf called on empty cell %d\n", modelgridindex);
    abort();
  }

  const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
  if (deposition_rate_density > 0.)
  {
    return deposition_rate_density / get_tot_nion(modelgridindex) / get_eff_ionpot(modelgridindex, element, ion);
    // alternatively, if the y vector is still in memory:
    // return calculate_nt_ionization_ratecoeff(modelgridindex, element, ion);
  }
  else
    return 0.;
}


double nt_ionization_ratecoeff(int modelgridindex, int element, int ion)
{
  if (!NT_ON)
  {
    printout("ERROR: NT_ON is false, but nt_ionization_ratecoeff has been called.\n");
    abort();
  }
  if (mg_associated_cells[modelgridindex] <= 0)
  {
    printout("ERROR: nt_ionization_ratecoeff called on empty cell %d\n", modelgridindex);
    abort();
  }

  if (NT_SOLVE_SPENCERFANO)
  {
    double Y_nt = nt_ionization_ratecoeff_sf(modelgridindex, element, ion);
    if (!isfinite(Y_nt))
    {
      // probably because eff_ionpot = 0 because the solver hasn't been run yet, or no impact ionization cross sections exist
      const double Y_nt_wfapprox = nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion);
      // printout("Warning: Spencer-Fano solver gives non-finite ionization rate (%g) for element %d ion_stage %d for cell %d. Using WF approx instead = %g\n",
      //          Y_nt, get_element(element), get_ionstage(element, ion), modelgridindex, Y_nt_wfapprox);
      return Y_nt_wfapprox;
    }
    else if (Y_nt <= 0)
    {
      const double Y_nt_wfapprox = nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion);
      if (Y_nt_wfapprox > 0)
      {
        printout("Warning: Spencer-Fano solver gives negative or zero ionization rate (%g) for element Z=%d ion_stage %d cell %d. Using WF approx instead = %g\n",
                 Y_nt, get_element(element), get_ionstage(element, ion), modelgridindex, Y_nt_wfapprox);
      }
      return Y_nt_wfapprox;
    }
    else
    {
      return Y_nt;
    }
  }
  else
    return nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion);
}


static double calculate_nt_frac_excitation_perlevelpop(const int modelgridindex, const int lineindex, const double epsilon_trans)
// Kozma & Fransson equation 9 divided by level population
{
  if (nt_solution[modelgridindex].yfunc == NULL)
  {
    printout("ERROR: Call to nt_excitation_ratecoeff with no y vector in memory.");
    abort();
  }
  const double epsilon_trans_ev = epsilon_trans / EV;

  gsl_vector *xs_excitation_vec = gsl_vector_alloc(SFPTS);
  if (get_xs_excitation_vector(xs_excitation_vec, lineindex, epsilon_trans))
  {
    double y_dot_crosssection = 0.;
    gsl_vector_view yvecview = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);
    gsl_blas_ddot(xs_excitation_vec, &yvecview.vector, &y_dot_crosssection);
    gsl_vector_free(xs_excitation_vec);

    return epsilon_trans_ev * y_dot_crosssection * DELTA_E / E_init_ev;
  }
  else
  {
    gsl_vector_free(xs_excitation_vec);

    return 0.;
  }
}


double nt_excitation_ratecoeff(const int modelgridindex, const int lineindex)
{
  if (!NT_EXCITATION_ON)
    return 0.;

  const int list_size = nt_solution[modelgridindex].frac_excitations_list_size;
  for (int excitationindex = 0; excitationindex < list_size; excitationindex++)
  {
    if (nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex == lineindex)
    {
      const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
      const double ratecoeffperdeposition = nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition;

      return ratecoeffperdeposition * deposition_rate_density;
    }
  }
  return 0.;
}


void do_ntlepton(PKT *pkt_ptr)
{
  const int modelgridindex = cell[pkt_ptr->where].modelgridindex;
  double zrand = gsl_rng_uniform(rng);
  // const double frac_heating = get_nt_frac_heating(modelgridindex);
  const double frac_excitation = get_nt_frac_excitation(modelgridindex);
  const double frac_ionization = get_nt_frac_ionization(modelgridindex);

  if (zrand < frac_ionization)
  {
    const int frac_ionizations_list_size = nt_solution[modelgridindex].frac_ionizations_list_size;
    for (int allionindex = 0; allionindex < frac_ionizations_list_size; allionindex++)
    {
      const double frac_deposition_ion = nt_solution[modelgridindex].frac_ionizations_list[allionindex].frac_deposition;
      if (zrand < frac_deposition_ion)
      {
        const int element = nt_solution[modelgridindex].frac_ionizations_list[allionindex].element;
        const int lowerion = nt_solution[modelgridindex].frac_ionizations_list[allionindex].ion;

        mastate[tid].element = element;
        mastate[tid].ion = lowerion + 1;
        mastate[tid].level = 0;
        mastate[tid].activatingline = -99;
        pkt_ptr->type = TYPE_MA;
        ma_stat_activation_ntcollion++;
        pkt_ptr->interactions += 1;
        pkt_ptr->last_event = 9;
        pkt_ptr->trueemissiontype = -1; // since this is below zero, macroatom will set it
        pkt_ptr->trueemissionvelocity = -1;

        printout("NTLEPTON packet selected ionization of Z=%d ionstage %d\n", get_element(element), get_ionstage(element, lowerion));

        return;
      }
      zrand -= frac_deposition_ion;
    }
  }
  else if (NT_EXCITATION_ON && zrand < frac_ionization + frac_excitation)
  {
    zrand -= frac_ionization;
    const int frac_excitations_list_size = nt_solution[modelgridindex].frac_excitations_list_size;
    for (int excitationindex = 0; excitationindex < frac_excitations_list_size; excitationindex++)
    {
      const double frac_deposition_exc = nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition;
      if (zrand < frac_deposition_exc)
      {
        const int lineindex = nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex;
        const int element = linelist[lineindex].elementindex;
        const int ion = linelist[lineindex].ionindex;
        const int lower = linelist[lineindex].lowerlevelindex;
        const int upper = linelist[lineindex].upperlevelindex;

        mastate[tid].element = element;
        mastate[tid].ion = ion;
        mastate[tid].level = upper;
        mastate[tid].activatingline = -99;
        pkt_ptr->type = TYPE_MA;
        ma_stat_activation_ntcollexc++;
        pkt_ptr->interactions += 1;
        pkt_ptr->last_event = 8;
        pkt_ptr->trueemissiontype = -1; // since this is below zero, macroatom will set it
        pkt_ptr->trueemissionvelocity = -1;

        printout("NTLEPTON packet selected excitation of Z=%d ionstage %d level %d upperlevel %d\n",
                 get_element(element), get_ionstage(element, ion), lower, upper);

        return;
      }
      zrand -= frac_deposition_exc;
    }
  }


  /*It's an electron - convert to k-packet*/
  //printout("e-minus propagation\n");
  pkt_ptr->type = TYPE_KPKT;
  #ifndef FORCE_LTE
    //kgammadep[pkt_ptr->where] += pkt_ptr->e_cmf;
  #endif
  //pkt_ptr->type = TYPE_PRE_KPKT;
  //pkt_ptr->type = TYPE_GAMMA_KPKT;
  //if (tid == 0) k_stat_from_eminus += 1;
  k_stat_from_eminus += 1;
}


static void analyse_sf_solution(int modelgridindex)
{
  const float nne = get_nne(modelgridindex);
  const double nntot = get_tot_nion(modelgridindex);

  // store the solution properties now while the NT spectrum is in memory (in case we free before packet prop)
  nt_solution[modelgridindex].frac_heating = calculate_frac_heating(modelgridindex);

  double frac_excitation_total = 0.;
  double frac_ionization_total = 0.;

  int allionindex = 0; // unique index for every ion of all elements
#if NT_EXCITATION_ON
  int excitationindex = 0; // unique index for every included excitation transition
#endif
  for (int element = 0; element < nelements; element++)
  {
    const int Z = get_element(element);
    const int nions = get_nions(element);
    for (int ion = 0; ion < get_nions(element); ion++)
    {
      float eff_ionpot = calculate_eff_ionpot(modelgridindex, element, ion);
      if (!isfinite(eff_ionpot))
        eff_ionpot = 0.;
      nt_solution[modelgridindex].eff_ionpot[element][ion] = eff_ionpot;

      const int ionstage = get_ionstage(element, ion);
      const double nnion = ionstagepop(modelgridindex, element, ion);

      // if (nnion < minionfraction * get_tot_nion(modelgridindex)) // skip negligible ions
      if (nnion <= 0.) // skip zero-abundance ions
        continue;

      double frac_ionization_ion = 0.;
      double frac_excitation_ion = 0.;
      printout("  Z=%d ion_stage %d:\n", Z, ionstage);
      // printout("    nnion: %g\n", nnion);
      printout("    nnion/nntot: %g\n", nnion / nntot, get_nne(modelgridindex));
      int matching_nlsubshell_count = 0;
      for (int n = 0; n < colliondatacount; n++)
      {
        if (colliondata[n].Z == Z && colliondata[n].nelec == Z - ionstage + 1)
        {
          const double frac_ionization_ion_shell = calculate_nt_frac_ionization_shell(modelgridindex, element, ion, n);
          frac_ionization_ion += frac_ionization_ion_shell;
          matching_nlsubshell_count++;
          printout("      frac_ionization_shell(n %d l %d): %g (ionpot %.2f eV)\n",
                   colliondata[n].n, colliondata[n].l, frac_ionization_ion_shell, colliondata[n].ionpot_ev);
        }
      }

      if (ion < nions - 1)
      {
        if (allionindex >= nt_solution[modelgridindex].frac_ionizations_list_size)
        {
          nt_solution[modelgridindex].frac_ionizations_list_size += BLOCKSIZEIONIZATION;
          nt_solution[modelgridindex].frac_ionizations_list = realloc(
            nt_solution[modelgridindex].frac_ionizations_list,
            nt_solution[modelgridindex].frac_ionizations_list_size * sizeof(struct nt_ionization_struct));

          if (nt_solution[modelgridindex].frac_ionizations_list == NULL)
          {
            printout("ERROR: Not enough memory to reallocate ionization list.\n");
            abort();
          }
        }

        nt_solution[modelgridindex].frac_ionizations_list[allionindex].frac_deposition = frac_ionization_ion;
        nt_solution[modelgridindex].frac_ionizations_list[allionindex].element = element;
        nt_solution[modelgridindex].frac_ionizations_list[allionindex].ion = ion;
        allionindex++;
        frac_ionization_total += frac_ionization_ion;
      }

      const int maxlevel = 0; // just consider excitation from the ground level
      // const int maxlevel = get_nlevels(element, ion); // excitation from all levels (SLOW)

      for (int level = 0; level <= maxlevel; level++)
      {
        const int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
        const int nnlevel = calculate_exclevelpop(modelgridindex, element, ion, level);

        for (int t = 1; t <= nuptrans; t++)
        {
          const double epsilon_trans = elements[element].ions[ion].levels[level].uptrans[t].epsilon_trans;
          const int lineindex = elements[element].ions[ion].levels[level].uptrans[t].lineindex;

          const double nt_frac_excitation_perlevelpop = calculate_nt_frac_excitation_perlevelpop(modelgridindex, lineindex, epsilon_trans);
          const double frac_excitation_thistrans = nnlevel * nt_frac_excitation_perlevelpop;
          frac_excitation_ion += frac_excitation_thistrans;

#if NT_EXCITATION_ON
          // if (frac_excitation_thistrans > 0.)
          {
            if (excitationindex >= nt_solution[modelgridindex].frac_excitations_list_size)
            {
              nt_solution[modelgridindex].frac_excitations_list_size += BLOCKSIZEEXCITATION;
              nt_solution[modelgridindex].frac_excitations_list = realloc(
                nt_solution[modelgridindex].frac_excitations_list,
                nt_solution[modelgridindex].frac_excitations_list_size * sizeof(struct nt_excitation_struct));

              if (nt_solution[modelgridindex].frac_excitations_list == NULL)
              {
                printout("ERROR: Not enough memory to reallocate excitation list.\n");
                abort();
              }
            }

            const double ratecoeffperdeposition = nt_frac_excitation_perlevelpop / epsilon_trans;
            nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition = frac_excitation_thistrans;
            nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition = ratecoeffperdeposition;
            nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex = lineindex;
            (excitationindex)++;
          }
#endif // NT_EXCITATION_ON
        } // for t
      } // for level

      // alternative way to calculate it
      // const double frac_excitation_ion_2 = calculate_nt_frac_excitation_ion(modelgridindex, element, ion);
      frac_excitation_total += frac_excitation_ion;

      printout("    frac_ionization: %g (%d subshells)\n", frac_ionization_ion, matching_nlsubshell_count);
      printout("    frac_excitation: %g\n", frac_excitation_ion);
      printout("    workfn:       %9.2f eV\n", (1. / get_oneoverw(element, ion, modelgridindex)) / EV);
      printout("    eff_ionpot:   %9.2f eV\n", get_eff_ionpot(modelgridindex, element, ion) / EV);
      printout("    workfn approx Gamma: %9.3e\n",
               nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion));
      printout("    Spencer-Fano Gamma:  %9.3e\n",
               nt_ionization_ratecoeff_sf(modelgridindex, element, ion));
      printout("    SF integral Gamma:   %9.3e\n",
               calculate_nt_ionization_ratecoeff(modelgridindex, element, ion));
    }
  }

  if (allionindex < nt_solution[modelgridindex].frac_ionizations_list_size)
  {
    // shrink the list to match the data
    nt_solution[modelgridindex].frac_ionizations_list = realloc(
      nt_solution[modelgridindex].frac_ionizations_list, (allionindex) * sizeof(struct nt_ionization_struct));
    nt_solution[modelgridindex].frac_ionizations_list_size = allionindex;
  }

  qsort(nt_solution[modelgridindex].frac_ionizations_list,
        nt_solution[modelgridindex].frac_ionizations_list_size, sizeof(struct nt_ionization_struct),
        compare_ionization_fractions);

#if NT_EXCITATION_ON
  if (excitationindex < nt_solution[modelgridindex].frac_excitations_list_size)
  {
    // shrink the list to match the data
    nt_solution[modelgridindex].frac_excitations_list = realloc(
      nt_solution[modelgridindex].frac_excitations_list, (excitationindex) * sizeof(struct nt_excitation_struct));
    nt_solution[modelgridindex].frac_excitations_list_size = excitationindex;
  }

  qsort(nt_solution[modelgridindex].frac_excitations_list,
        nt_solution[modelgridindex].frac_excitations_list_size, sizeof(struct nt_excitation_struct),
        compare_excitation_fractions);

  const float T_e = get_Te(modelgridindex);
  printout("Top non-thermal excitation fractions (total excitations = %d):\n",
           nt_solution[modelgridindex].frac_excitations_list_size);
  int ntransdisplayed = nt_solution[modelgridindex].frac_excitations_list_size;
  ntransdisplayed = (ntransdisplayed > 20) ? 20 : ntransdisplayed;
  for (excitationindex = 0; excitationindex < ntransdisplayed; excitationindex++)
  {
    const double frac_deposition = nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition;
    const int lineindex = nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex;
    const int element = linelist[lineindex].elementindex;
    const int ion = linelist[lineindex].ionindex;
    const int lower = linelist[lineindex].lowerlevelindex;
    const int upper = linelist[lineindex].upperlevelindex;
    const double epsilon_trans = epsilon(element, ion, upper) - epsilon(element, ion, lower);

    const double ntexc_ratecoeff = nt_excitation_ratecoeff(modelgridindex, lineindex);
    const double exc_ratecoeff = col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans);
    if (frac_deposition > 0.)
    {
      printout("frac_deposition %.3e element Z=%d ionstage %d level %d upperlevel %d ratecoeff %.3e ntratecoeff %.3e nt/t %.3e\n",
               frac_deposition, get_element(element), get_ionstage(element, ion), lower, upper,
               exc_ratecoeff, ntexc_ratecoeff, ntexc_ratecoeff / exc_ratecoeff);
    }
  }
#endif // NT_EXCITATION_ON

  const float frac_heating = get_nt_frac_heating(modelgridindex);

  // calculate number density of non-thermal electrons
  const double deposition_rate_density_ev = get_deposition_rate_density(modelgridindex) / EV;
  const double yscalefactor = deposition_rate_density_ev / E_init_ev;

  double nne_nt_max = 0.0;
  for (int i = 0; i < SFPTS; i++)
  {
    const double endash = gsl_vector_get(envec, i);
    const double oneovervelocity = sqrt(9.10938e-31 / 2 / endash / 1.60218e-19) / 100; // in sec/cm
    nne_nt_max += yscalefactor * get_y_sample(modelgridindex, i) * oneovervelocity * DELTA_E;
  }

  nt_solution[modelgridindex].frac_excitation = frac_excitation_total;
  nt_solution[modelgridindex].frac_ionization = frac_ionization_total;

  printout("  E_0:         %9.4f eV\n", nt_solution[modelgridindex].E_0);
  printout("  E_init:      %9.2f eV/s/cm^3\n", E_init_ev);
  printout("  deposition:  %9.2f eV/s/cm^3\n", deposition_rate_density_ev);
  printout("  nne:         %9.3e e-/cm^3\n", nne);
  printout("  nne_nt     < %9.3e e-/cm^3\n", nne_nt_max);
  printout("  nne_nt/nne < %9.3e\n", nne_nt_max / nne);
  printout("  frac_heating_tot:    %g\n", frac_heating);
  printout("  frac_excitation_tot: %g\n", frac_excitation_total);
  printout("  frac_ionization_tot: %g\n", frac_ionization_total);
  const double frac_sum = frac_heating + frac_excitation_total + frac_ionization_total;
  printout("  frac_sum:            %g (should be close to 1.0)\n", frac_sum);

  // const double nnion = ionstagepop(modelgridindex, element, ion);
  // double ntexcit_in_a = 0.;
  // for (int level = 0; level < get_nlevels(0, 1); level++)
  // {
  //   ntexcit_in_a += nnion * nt_excitation_ratecoeff(modelgridindex, 0, 1, level, 75);
  // }
  // printout("  total nt excitation rate into level 75: %g\n", ntexcit_in_a);

  // compensate for lost energy by scaling the solution
  // E_init_ev *= frac_sum;
}


static void sfmatrix_add_excitation(gsl_matrix *sfmatrix, const int modelgridindex, const int element, const int ion, double *E_0)
{
  // excitation terms
  gsl_vector *vec_xs_excitation_nnion_deltae = gsl_vector_alloc(SFPTS);
  const int maxlevel = 0; // just consider excitation from the ground level
  // const int maxlevel = get_nlevels(element, ion); // excitation from all levels (SLOW)

  for (int level = 0; level <= maxlevel; level++)
  {
    const double nnlevel = calculate_exclevelpop(modelgridindex, element, ion, level);
    const int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
    for (int t = 1; t <= nuptrans; t++)
    {
      const double epsilon_trans = elements[element].ions[ion].levels[level].uptrans[t].epsilon_trans;
      const int lineindex = elements[element].ions[ion].levels[level].uptrans[t].lineindex;
      const double epsilon_trans_ev = epsilon_trans / EV;

      if (epsilon_trans / EV < *E_0 || *E_0 <= 0.)
        *E_0 = epsilon_trans / EV;

      if (get_xs_excitation_vector(vec_xs_excitation_nnion_deltae, lineindex, epsilon_trans))
      {
        gsl_blas_dscal(nnlevel * DELTA_E, vec_xs_excitation_nnion_deltae);

        for (int i = 0; i < SFPTS; i++)
        {
          const double en = gsl_vector_get(envec, i);
          const int stopindex = get_energyindex_ev(en + epsilon_trans_ev);
          if (stopindex < SFPTS - 1)
          {
            gsl_vector_view a = gsl_matrix_subrow(sfmatrix, i, i, stopindex - i + 1);
            gsl_vector_const_view b = gsl_vector_const_subvector(vec_xs_excitation_nnion_deltae, i, stopindex - i + 1);
            gsl_vector_add(&a.vector, &b.vector); // add b to a and put the result in a
          }
        }
      }
    }
  }
  gsl_vector_free(vec_xs_excitation_nnion_deltae);
}


static void sfmatrix_add_ionization(gsl_matrix *sfmatrix, const int Z, const int ionstage, const double nnion, double *E_0)
// add the ionization terms to the Spencer-Fano matrix
// also, update the value of E_0, the minimum energy for excitation/ionization
{
  for (int n = 0; n < colliondatacount; n++)
  {
    if (colliondata[n].Z == Z && colliondata[n].nelec == Z - ionstage + 1)
    {
      const double ionpot_ev = colliondata[n].ionpot_ev;
      const double J = get_J(Z, ionstage, ionpot_ev);

      if (ionpot_ev < *E_0 || *E_0 <= 0.)
        *E_0 = ionpot_ev; // this ionization potential is the new minimum energy

      // printout("Z=%2d ion_stage %d n %d l %d ionpot %g eV\n",
      //          Z, ionstage, colliondata[n].n, colliondata[n].l, ionpot_ev);

      for (int i = 0; i < SFPTS; i++)
      {
        // i is the matrix row index, which corresponds to an energy E at which we are solve from y(E)
        const double en = gsl_vector_get(envec, i);

        const int secondintegralstartindex = get_energyindex_ev(2 * en + ionpot_ev);

        for (int j = i; j < SFPTS; j++)
        {
          // j is the matrix column index which corresponds to the piece of the integral at y(E') where E' >= E and E' = envec(j)
          const double endash = gsl_vector_get(envec, j);

          const double prefactor = nnion * xs_impactionization(endash, n) / atan((endash - ionpot_ev) / 2 / J);

          const double epsilon_upper = (endash + ionpot_ev) / 2;
          double epsilon_lower = endash - en;
          // atan bit is the definite integral of 1/[1 + (epsilon - I)/J] in Kozma & Fransson 1992 equation 4
          double ij_contribution = prefactor * (atan((epsilon_upper - ionpot_ev) / J) - atan((epsilon_lower - ionpot_ev) / J)) * DELTA_E;

          if (j >= secondintegralstartindex)
          {
            epsilon_lower = en + ionpot_ev;
            ij_contribution -= prefactor * (atan((epsilon_upper - ionpot_ev) / J) - atan((epsilon_lower - ionpot_ev) / J)) * DELTA_E;
          }
          *gsl_matrix_ptr(sfmatrix, i, j) += ij_contribution;
        }
      }
    }
  }
}


static void sfmatrix_solve(const gsl_matrix *sfmatrix, const gsl_vector *rhsvec, gsl_vector *yvec)
{
  // WARNING: this assumes sfmatrix is in upper triangular form already!
  const gsl_matrix *sfmatrix_LU = sfmatrix;
  gsl_permutation *p = gsl_permutation_calloc(SFPTS); // identity permutation

  // if the matrix is not upper triangular, then do a decomposition
  // printout("Doing LU decomposition of SF matrix\n");
  // make a copy of the matrix for the LU decomp
  // gsl_matrix *sfmatrix_LU = gsl_matrix_alloc(SFPTS, SFPTS);
  // gsl_matrix_memcpy(sfmatrix_LU, sfmatrix);
  // int s; //sign of the transformation
  // gsl_permutation *p = gsl_permutation_alloc(SFPTS);
  // gsl_linalg_LU_decomp(sfmatrix_LU, p, &s);

  // printout("Solving SF matrix equation\n");

  // solve matrix equation: sf_matrix * y_vec = rhsvec for yvec
  gsl_linalg_LU_solve(sfmatrix_LU, p, rhsvec, yvec);
  // printout("Refining solution\n");

  double error_best = -1.;
  gsl_vector *yvec_best = gsl_vector_alloc(SFPTS); // solution vector with lowest error
  gsl_vector *gsl_work_vector = gsl_vector_calloc(SFPTS);
  gsl_vector *residual_vector = gsl_vector_alloc(SFPTS);
  int iteration;
  for (iteration = 0; iteration < 10; iteration++)
  {
    if (iteration > 0)
      gsl_linalg_LU_refine(sfmatrix, sfmatrix_LU, p, rhsvec, yvec, gsl_work_vector); // first argument must be original matrix

    // calculate Ax - b = residual
    gsl_vector_memcpy(residual_vector, rhsvec);
    gsl_blas_dgemv(CblasNoTrans, 1.0, sfmatrix, yvec, -1.0, residual_vector);

    // value of the largest absolute residual
    const double error = fabs(gsl_vector_get(residual_vector, gsl_blas_idamax(residual_vector)));

    if (error < error_best || error_best < 0.)
    {
      gsl_vector_memcpy(yvec_best, yvec);
      error_best = error;
    }
    // printout("Linear algebra solver iteration %d has a maximum residual of %g\n",iteration,error);
  }
  if (error_best >= 0.)
  {
    if (error_best > 1e-10)
      printout("  SF solver LU_refine: After %d iterations, best solution vector has a max residual of %g (WARNING)\n",
               iteration, error_best);
    gsl_vector_memcpy(yvec, yvec_best);
  }
  gsl_vector_free(yvec_best);
  gsl_vector_free(gsl_work_vector);
  gsl_vector_free(residual_vector);

  // gsl_matrix_free(sfmatrix_LU); // if this matrix is different to sfmatrix then free it
  gsl_permutation_free(p);

  if (!gsl_vector_isnonneg(yvec))
  {
    printout("solve_sfmatrix: WARNING: y function goes negative!\n");
  }
}


void nt_solve_spencerfano(const int modelgridindex, const int timestep)
// solve the Spencer-Fano equation to get the non-thermal electron flux energy distribution
// based on Equation (2) of Li et al. (2012)
{
  // timestep now records when the deposition rate density was set
  // nt_solution[modelgridindex].timestep = timestep;
  if (mg_associated_cells[modelgridindex] < 1)
  {
    printout("Associated_cells < 1 in cell %d at timestep %d. Skipping Spencer-Fano solution.\n", modelgridindex, timestep);

    return;
  }
  else
  {
    const double deposition_rate_density_ev = get_deposition_rate_density(modelgridindex) / EV;
    if (deposition_rate_density_ev < MINDEPRATE)
    {
      printout("Non-thermal deposition rate of %g eV/cm/s/cm^3 in cell %d at timestep %d. Skipping Spencer-Fano solution.\n", deposition_rate_density_ev, modelgridindex, timestep);

      // if (!STORE_NT_SPECTRUM)
      // {
      //   nt_solution[modelgridindex].yfunc = calloc(SFPTS, sizeof(double));
      // }

      // nt_write_to_file(modelgridindex, timestep);

      nt_solution[modelgridindex].timestep = timestep;
      nt_solution[modelgridindex].frac_heating = 1.0;
      nt_solution[modelgridindex].E_0 = 0.;

      zero_all_effionpot(modelgridindex);

      return;
    }
  }

  const float nne = get_nne(modelgridindex); // electrons per cm^3

  printout("Setting up Spencer-Fano equation with %d energy points from %g eV to %g eV in cell %d at timestep %d (nne=%g e-/cm^3)\n",
           SFPTS, EMIN, EMAX, modelgridindex, timestep, nne);

  gsl_matrix *const sfmatrix = gsl_matrix_calloc(SFPTS, SFPTS);
  gsl_vector *const rhsvec = gsl_vector_calloc(SFPTS); // constant term (not dependent on y func) in each equation

  // loss terms and source terms
  for (int i = 0; i < SFPTS; i++)
  {
    const double en = gsl_vector_get(envec, i);

    *gsl_matrix_ptr(sfmatrix, i, i) += electron_loss_rate(en * EV, nne) / EV;

    double source_integral_to_emax;
    if (i < SFPTS - 1)
    {
      gsl_vector_const_view source_e_to_emax = gsl_vector_const_subvector(sourcevec, i + 1, SFPTS - i - 1);
      source_integral_to_emax = gsl_blas_dasum(&source_e_to_emax.vector) * DELTA_E;
    }
    else
      source_integral_to_emax = 0;

    gsl_vector_set(rhsvec, i, source_integral_to_emax);
  }
  // gsl_vector_set_all(rhsvec, 1.); // alternative if all electrons are injected at EMAX

  double E_0 = 0.; // reset E_0 so it can be found it again

  for (int element = 0; element < nelements; element++)
  {
    const int Z = get_element(element);
    const int nions = get_nions(element);
    bool first_included_ion_of_element = true;
    for (int ion = 0; ion < nions; ion++)
    {
      const double nnion = ionstagepop(modelgridindex, element, ion); // hopefully ions per cm^3?

      if (nnion < minionfraction * get_tot_nion(modelgridindex)) // skip negligible ions
      {
        continue;
      }

      const int ionstage = get_ionstage(element, ion);
      if (first_included_ion_of_element)
      {
        printout("  including Z=%2d ion_stages: ", Z);
        for (int i = 1; i < get_ionstage(element, ion); i++)
          printout("  ");
        first_included_ion_of_element = false;
      }

      printout("%d ", ionstage);

      sfmatrix_add_excitation(sfmatrix, modelgridindex, element, ion, &E_0);
      if (ion < nions - 1)
        sfmatrix_add_ionization(sfmatrix, Z, ionstage, nnion, &E_0);
    }
    if (!first_included_ion_of_element)
      printout("\n");
  }

  // printout("SF matrix | RHS vector:\n");
  // for (int row = 0; row < 10; row++)
  // {
  //   for (int column = 0; column < 10; column++)
  //   {
  //     char str[15];
  //     sprintf(str, "%+.1e ", gsl_matrix_get(sfmatrix, row, column));
  //     printout(str);
  //   }
  //   printout("| ");
  //   char str[15];
  //   sprintf(str, "%+.1e\n", gsl_vector_get(rhsvec, row));
  //   printout(str);
  // }
  // printout("\n");

  if (!STORE_NT_SPECTRUM)
  {
    nt_solution[modelgridindex].yfunc = calloc(SFPTS, sizeof(double));
  }

  gsl_vector_view yvecview = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);
  gsl_vector *yvec = &yvecview.vector;
  sfmatrix_solve(sfmatrix, rhsvec, yvec);

  // gsl_matrix_free(sfmatrix_LU); // if this matrix is different to sfmatrix

  gsl_matrix_free(sfmatrix);
  gsl_vector_free(rhsvec);

  if (timestep % 10 == 0)
    nt_write_to_file(modelgridindex, timestep);

  nt_solution[modelgridindex].frac_heating = -1.;
  nt_solution[modelgridindex].E_0 = E_0;

  analyse_sf_solution(modelgridindex);

  if (!STORE_NT_SPECTRUM)
  {
    free(nt_solution[modelgridindex].yfunc);
  }
}


void nt_write_restart_data(FILE *gridsave_file)
{
  if (!NT_SOLVE_SPENCERFANO)
    return;

  printout("Writing restart data for non-thermal solver\n");

  fprintf(gridsave_file, "%d\n", 24724518); // special number marking the beginning of NT data
  fprintf(gridsave_file, "%d %lg %lg\n", SFPTS, EMIN, EMAX);

  if (STORE_NT_SPECTRUM)
  {
    printout("nt_write_restart_data not implemented for STORE_NT_SPECTRUM ON");
    abort();
  }
  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    if (mg_associated_cells[modelgridindex] > 0)
    {
      fprintf(gridsave_file, "%d %d %lg %g %lg ",
              modelgridindex,
              nt_solution[modelgridindex].timestep,
              nt_solution[modelgridindex].E_0,
              nt_solution[modelgridindex].frac_heating,
              nt_solution[modelgridindex].deposition_rate_density);

      for (int element = 0; element < nelements; element++)
      {
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions; ion++)
        {
          fprintf(gridsave_file, "%g ", nt_solution[modelgridindex].eff_ionpot[element][ion]);
        }
      }
    }
  }
}


void nt_read_restart_data(FILE *gridsave_file)
{
  if (!NT_SOLVE_SPENCERFANO)
    return;

  printout("Reading restart data for non-thermal solver\n");

  int code_check;
  fscanf(gridsave_file, "%d\n", &code_check);
  if (code_check != 24724518)
  {
    printout("ERROR: Beginning of non-thermal restart data not found!");
    abort();
  }

  int sfpts_in;
  double emin_in;
  double emax_in;
  fscanf(gridsave_file, "%d %lg %lg\n", &sfpts_in, &emin_in, &emax_in);

  if (sfpts_in != SFPTS || emin_in != EMIN || emax_in != EMAX)
  {
    printout("ERROR: gridsave file specifies %d Spencer-Fano samples, emin %lg emax %lg\n",
             sfpts_in, emin_in, emax_in);
    printout("ERROR: This simulation has %d Spencer-Fano samples, emin %lg emax %lg\n",
             SFPTS, EMIN, EMAX);
    abort();
  }

  if (STORE_NT_SPECTRUM)
  {
    printout("nt_write_restart_data not implemented for STORE_NT_SPECTRUM ON");
    abort();
  }
  for (int modelgridindex = 0; modelgridindex < MMODELGRID; modelgridindex++)
  {
    if (mg_associated_cells[modelgridindex] > 0)
    {
      int mgi_in;
      fscanf(gridsave_file, "%d %d %lg %g %lg ",
             &mgi_in,
             &nt_solution[modelgridindex].timestep,
             &nt_solution[modelgridindex].E_0,
             &nt_solution[modelgridindex].frac_heating,
             &nt_solution[modelgridindex].deposition_rate_density);
      if (mgi_in != modelgridindex)
      {
        printout("ERROR: expected data for cell %d but found cell %d\n", modelgridindex, mgi_in);
        abort();
      }

      for (int element = 0; element < nelements; element++)
      {
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions; ion++)
        {
          fscanf(gridsave_file, "%g ", &nt_solution[modelgridindex].eff_ionpot[element][ion]);
        }
      }
    }
  }
}


#ifdef MPI_ON
void nt_MPI_Bcast(const int my_rank, const int root, const int root_nstart, const int root_ndo)
{
  if (!nonthermal_initialized)
    return;

  // const int logged_element_z = 26;
  // const int logged_ion_index = 1;
  // const int logged_element_index = get_elementindex(logged_element_z);
  // const int logged_ion_stage = get_ionstage(logged_element_index, logged_ion_index);

  if (root_ndo > 0)
  {
    if (my_rank == root)
    {
      // printout("nonthermal_MPI_Bcast root process %d will broadcast cells %d to %d\n",
      //          my_rank, root_nstart, root_nstart + root_ndo - 1);
    }
    else
    {
      // printout("nonthermal_MPI_Bcast process %d will receive cells %d to %d from process %d\n",
      //          my_rank, root_nstart, root_nstart + root_ndo - 1, root);
    }
  }

  for (int modelgridindex = root_nstart; modelgridindex < root_nstart + root_ndo; modelgridindex++)
  {
    if (mg_associated_cells[modelgridindex] > 0)
    {
      // printout("nonthermal_MPI_Bcast cell %d before: ratecoeff(Z=%d ion_stage %d): %g, eff_ionpot %g eV\n",
      //          modelgridindex, logged_element_z, logged_ion_stage,
      //          nt_ionization_ratecoeff_sf(modelgridindex, logged_element_index, logged_ion_index),
      //          get_eff_ionpot(modelgridindex, logged_element_index, logged_ion_index) / EV);
      MPI_Barrier(MPI_COMM_WORLD);
      if (STORE_NT_SPECTRUM)
      {
        MPI_Bcast(&nt_solution[modelgridindex].yfunc, SFPTS, MPI_DOUBLE, root, MPI_COMM_WORLD);
        // printout("nonthermal_MPI_Bcast Bcast y vector for cell %d from process %d to %d\n", modelgridindex, root, my_rank);
      }
      MPI_Bcast(&nt_solution[modelgridindex].timestep, 1, MPI_INT, root, MPI_COMM_WORLD);
      MPI_Bcast(&nt_solution[modelgridindex].frac_heating, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
      MPI_Bcast(&nt_solution[modelgridindex].E_0, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&nt_solution[modelgridindex].deposition_rate_density, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      for (int element = 0; element < nelements; element++)
      {
        const int nions = get_nions(element);
        MPI_Bcast(&nt_solution[modelgridindex].eff_ionpot[element], nions, MPI_FLOAT, root, MPI_COMM_WORLD);
      }
      // printout("nonthermal_MPI_Bcast cell %d after: ratecoeff(Z=%d ion_stage %d): %g, eff_ionpot %g eV\n",
      //          modelgridindex, logged_element_z, logged_ion_stage,
      //          nt_ionization_ratecoeff_sf(modelgridindex, logged_element_index, logged_ion_index),
      //          get_eff_ionpot(modelgridindex, logged_element_index, logged_ion_index) / EV);
    }
    else
    {
      // printout("nonthermal_MPI_Bcast Skipping empty grid cell %d.\n", modelgridindex);
    }
  }

}
#endif
