#include <cmath>
#include <cstdbool>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include "atomic.h"
#include "decay.h"
#include "grid.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "stats.h"
#include "update_grid.h"
#include "sn3d.h"


namespace nonthermal
{

// THESE OPTIONS ARE USED TO TEST THE SF SOLVER
// Compare to Kozma & Fransson (1992) pure-oxygen plasma, nne = 1e8, x_e = 0.01
// #define yscalefactoroverride(mgi) (1e10)
// #define get_tot_nion(x) (1e10)
// #define ionstagepop(modelgridindex, element, ion) ionstagepop_override(modelgridindex, element, ion)
// #define get_nne(x) (1e8)
// #define SFPTS 10000  // number of energy points in the Spencer-Fano solution vector
// #define SF_EMAX 3000. // eV
// #define SF_EMIN 1. // eV
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


#define STORE_NT_SPECTRUM false // if this is on, the non-thermal energy spectrum will be kept in memory for
                                // every grid cell during packet propagation, which
                                // can take up a lot of memory for large grid sizes
                                // alternatively, just the non-thermal ionization rates can be stored
                                // but we might want to re-enable this option to incorporate
                                // non-thermal excitation rates if there are
                                // many more transitions to store than there are NT spectrum samples

// minimum number fraction of the total population to include in SF solution
static const double minionfraction = 1.e-8;

// minimum deposition rate density (eV/s/cm^3) to solve SF equation
static const double MINDEPRATE = 0.;

// Bohr radius squared in cm^2
static const double A_naught_squared = 2.800285203e-17;

// specifies max number of shells for which data is known for computing mean binding energies
#define M_NT_SHELLS 10

// maximum number of elements for which binding energy tables are to be used
#define MAX_Z_BINDING 30

static double electron_binding[MAX_Z_BINDING][M_NT_SHELLS];

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
  double auger_g_accumulated; // track the statistical weight represented by the values below, so they can be updated with new g-weighted averaged values
  double prob_num_auger[NT_MAX_AUGER_ELECTRONS + 1];  // probability of 0, 1, ..., NT_MAX_AUGER_ELECTRONS Auger electrons being ejected when the shell is ionised
  float en_auger_ev; // the average kinetic energy released in Auger electrons after making a hole in this shell
  float n_auger_elec_avg;
};

static struct collionrow *colliondata = NULL;
static int colliondatacount = 0;

static FILE *nonthermalfile = NULL;
static bool nonthermal_initialized = false;

static gsl_vector *envec;            // energy grid on which solution is sampled
static gsl_vector *logenvec;         // log of envec
static gsl_vector *sourcevec;        // samples of the source function (energy distribution of deposited energy)
static double E_init_ev = 0;         // the energy injection rate density (and mean energy of injected electrons if source integral is one) in eV

#if (SF_USE_LOG_E_INCREMENT)
  static gsl_vector *delta_envec;
  static double delta_log_e = 0.;
#else
  static const double DELTA_E = (SF_EMAX - SF_EMIN) / (SFPTS - 1);
#endif


// Monte Carlo result - compare to analytical expectation
double nt_energy_deposited;

struct nt_excitation_struct
{
  double frac_deposition;  // the fraction of the non-thermal deposition energy going to the excitation transition
  double ratecoeffperdeposition; // the excitation rate coefficient divided by the deposition rate density
  int lineindex;
};

struct nt_solution_struct {
  double *yfunc;  // Samples of the Spencer-Fano solution function. Multiply by energy to get non-thermal electron number flux.
                  // y(E) * dE is the flux of electrons with energy in the range (E, E + dE)
                  // y has units of particles / cm2 / s / eV

  float frac_heating;              // energy fractions should add up to 1.0 if the solution is good
  float frac_ionization;           // fraction of deposition energy going to ionization
  float frac_excitation;           // fraction of deposition energy going to excitation

  // these points arrays of length includedions
  float *eff_ionpot; // these are used to calculate the non-thermal ionization rate
  double *fracdep_ionization_ion; // the fraction of the non-thermal deposition energy going to ionizing this ion

  // these  point to arrays of length includedions * (NT_MAX_AUGER_ELECTRONS + 1)
  float *prob_num_auger;            // probability that one ionisation of this ion will produce n Auger electrons. elements sum to 1.0 for a given ion
  float *ionenfrac_num_auger;       // like above, but energy weighted. elements sum to 1.0 for an ion

  int frac_excitations_list_size;
  struct nt_excitation_struct *frac_excitations_list;

  int timestep_last_solved;                       // the quantities above were calculated for this timestep
  float nneperion_when_solved;                    // the nne when the solver was last run
};

static struct nt_solution_struct *nt_solution;

static double *deposition_rate_density;
static int *deposition_rate_density_timestep;

// for descending sort
static int compare_excitation_fractions(const void *p1, const void *p2)
{
  const struct nt_excitation_struct *elem1 = (struct nt_excitation_struct *) p1;
  const struct nt_excitation_struct *elem2 = (struct nt_excitation_struct *) p2;

 if (elem1->frac_deposition < elem2->frac_deposition)
    return 1;
 else if (elem1->frac_deposition > elem2->frac_deposition)
    return -1;
 else
    return 0;
}


// for ascending sort
static int compare_excitation_lineindicies(const void *p1, const void *p2)
{
  const struct nt_excitation_struct *elem1 = (struct nt_excitation_struct *) p1;
  const struct nt_excitation_struct *elem2 = (struct nt_excitation_struct *) p2;

 if (elem1->lineindex > elem2->lineindex)
    return 1;
 else if (elem1->lineindex < elem2->lineindex)
    return -1;
 else
    return 0;
}


#ifndef get_tot_nion
static double get_tot_nion(const int modelgridindex)
{
  double result = 0.;
  for (int element = 0; element < get_nelements(); element++)
  {
    result += globals::modelgrid[modelgridindex].composition[element].abundance / globals::elements[element].mass * get_rho(modelgridindex);

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


static void read_binding_energies(void)
{
  FILE *binding = fopen_required("binding_energies.txt", "r");

  int dum1, dum2;
  fscanf(binding, "%d %d", &dum1, &dum2); //dimensions of the table
  if ((dum1 != M_NT_SHELLS) || (dum2 != MAX_Z_BINDING))
  {
    printout("Wrong size for the binding energy tables!\n");
    abort();
  }

  for (int index1 = 0; index1 < dum2; index1++)
  {
    float dum[10];
    fscanf(binding, "%g %g %g %g %g %g %g %g %g %g",
           &dum[0], &dum[1], &dum[2], &dum[3], &dum[4], &dum[5], &dum[6], &dum[7], &dum[8], &dum[9]);

    for (int index2 = 0; index2 < 10; index2++)
    {
      electron_binding[index1][index2] = dum[index2] * EV;
    }
  }

  fclose(binding);
}

static double get_auger_probability(int modelgridindex, int element, int ion, int naugerelec)
{
  assert(naugerelec <= NT_MAX_AUGER_ELECTRONS);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  return nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + naugerelec];
}


static double get_ion_auger_enfrac(int modelgridindex, int element, int ion, int naugerelec)
{
  assert(naugerelec <= NT_MAX_AUGER_ELECTRONS);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  return nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + naugerelec];
}


static void check_auger_probabilities(int modelgridindex)
{
  bool problem_found = false;

  for (int element = 0; element < get_nelements(); element++)
  {
    for (int ion = 0; ion < get_nions(element) - 1; ion++)
    {
      double prob_sum = 0.;
      double ionenfrac_sum = 0.;
      for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
      {
        prob_sum += get_auger_probability(modelgridindex, element, ion, a);
        ionenfrac_sum += get_ion_auger_enfrac(modelgridindex, element, ion, a);
      }

      if (fabs(prob_sum - 1.0) > 0.001)
      {
        printout("Problem with Auger probabilities for cell %d Z=%d ionstage %d prob_sum %g\n", modelgridindex, get_element(element), get_ionstage(element, ion), prob_sum);
        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
        {
          printout("%d: %g\n", a, get_auger_probability(modelgridindex, element, ion, a));
        }
        problem_found = true;
      }

      if (fabs(ionenfrac_sum - 1.0) > 0.001)
      {
        printout("Problem with Auger energy frac sum for cell %d Z=%d ionstage %d ionenfrac_sum %g\n", modelgridindex, get_element(element), get_ionstage(element, ion), ionenfrac_sum);
        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
        {
          printout("%d: %g\n", a, get_ion_auger_enfrac(modelgridindex, element, ion, a));
        }
        problem_found = true;
      }
    }
  }

  if (problem_found)
  {
    abort();
  }
}


static void read_auger_data(void)
{
  printout("Reading Auger effect data...\n");
  FILE *augerfile = fopen_required("auger-km1993-table2.txt", "r");

  char line[151] = "";

  // map x-ray notation shells K L1 L2 L3 M1 M2 M3 to quantum numbers n and l
  const int xrayn[7] = {1, 2, 2, 2, 3, 3, 3};
  const int xrayl[7] = {0, 0, 1, 1, 0, 1, 1};
  const int xrayg[7] = {2, 2, 2, 4, 2, 2, 4}; // g statistical weight = 2j + 1

  while (!feof(augerfile))
  {
    if (line != fgets(line, 151, augerfile))
      break;

    int Z;
    int ionstage;
    int shellnum;

    char *linepos = line;
    int offset = 0;

    assert(sscanf(linepos, "%d %d%n", &Z, &ionstage, &offset) == 2);
    assert(offset == 5);
    linepos += offset;

    const int element = get_elementindex(Z);

    if (element >= 0 && get_ionstage(element, 0) <= ionstage && ionstage < (get_ionstage(element, 0) + get_nions(element)))
    {
      float ionpot_ev;
      float en_auger_ev_total_nocorrection;
      int epsilon_e3;

      assert(sscanf(linepos, "%d %g %g %d%n", &shellnum, &ionpot_ev, &en_auger_ev_total_nocorrection, &epsilon_e3, &offset) == 4);
      assert(offset == 20);
      linepos += offset + 1; // skip one space after so all following columns are exactly 5 characters each

      float n_auger_elec_avg = 0;
      double prob_num_auger[NT_MAX_AUGER_ELECTRONS + 1];
      for (int a = 0; a < 9; a++)
      {
        linepos = line + 26 + a * 5;
        // have to read out exactly 5 characters at a time because the columns are sometimes not separated by a space
        char strprob[6] = "00000";
        assert(sscanf(linepos, "%5c%n", strprob, &offset) == 1);
        assert(offset == 5);
        linepos += offset;
        strprob[5] = '\0';

        int probnaugerelece4;
        assert(sscanf(strprob, "%d", &probnaugerelece4) == 1);

        const double probnaugerelec = probnaugerelece4 / 10000.;

        assert(probnaugerelec <= 1.0);

        n_auger_elec_avg += a * probnaugerelec;

        if (a <= NT_MAX_AUGER_ELECTRONS)
        {
          prob_num_auger[a] = probnaugerelec;
        }
        else
        {
          // add the rates of all higher ionisations to the top one
          prob_num_auger[NT_MAX_AUGER_ELECTRONS] += probnaugerelec;
        }
      }

      // use the epsilon correction factor as in equation 7 of Kaastra & Mewe (1993)
      float en_auger_ev = en_auger_ev_total_nocorrection - (epsilon_e3 / 1000. * ionpot_ev);

      const int n = xrayn[shellnum - 1];
      const int l = xrayl[shellnum - 1];
      const int g = xrayg[shellnum - 1];

      if (!isfinite(en_auger_ev) || en_auger_ev < 0)
      {
        printout("  WARNING: Z=%2d ionstage %2d shellnum %d en_auger_ev is %g. Setting to zero.\n", Z, ionstage, shellnum, en_auger_ev);
        en_auger_ev = 0.;
      }

      // now loop through shells with impact ionisation cross sections and apply Auger data that matches n, l values
      for (int i = 0; i < colliondatacount; i++)
      {
        if (colliondata[i].Z == Z && colliondata[i].nelec == (Z - ionstage + 1) && colliondata[i].n == n && colliondata[i].l == l)
        {
          printout("Z=%2d ionstage %2d shellnum %d n %d l %d ionpot %7.2f E_A %8.1f E_A' %8.1f epsilon %6d <n_Auger> %5.1f P(n_Auger)",
                   Z, ionstage, shellnum, n, l, ionpot_ev, en_auger_ev_total_nocorrection, en_auger_ev, epsilon_e3, n_auger_elec_avg);

          double prob_sum = 0.;
          for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
          {
            prob_sum += prob_num_auger[a];
            printout(" %d: %4.2f", a, prob_num_auger[a]);
          }
          assert(fabs(prob_sum - 1.0) < 0.001);

          printout("\n");
          // printout("ionpot %g %g, g %d\n", colliondata[i].ionpot_ev, ionpot_ev, g);
          bool found_existing_data = (colliondata[i].auger_g_accumulated > 0.);

          // keep existing data but update according to statistical weight represented by existing and new data
          const double oldweight = colliondata[i].auger_g_accumulated / (g + colliondata[i].auger_g_accumulated);
          const double newweight = g / (g + colliondata[i].auger_g_accumulated);
          colliondata[i].auger_g_accumulated += g;

          // update the statistical-weight averaged values
          colliondata[i].en_auger_ev = oldweight * colliondata[i].en_auger_ev + newweight * en_auger_ev;
          colliondata[i].n_auger_elec_avg = oldweight * colliondata[i].n_auger_elec_avg + newweight * n_auger_elec_avg;

          prob_sum = 0.;
          for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
          {
            colliondata[i].prob_num_auger[a] = oldweight * colliondata[i].prob_num_auger[a] + newweight * prob_num_auger[a];
            prob_sum += colliondata[i].prob_num_auger[a];
          }
          assert(fabs(prob_sum - 1.0) < 0.001);

          if (found_existing_data)
          {
            printout("  same NL shell already has data from another X-ray shell. New g-weighted values: P(n_Auger)");

            for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
            {
              printout(" %d: %4.2f", a, colliondata[i].prob_num_auger[a]);
            }
            printout("\n");
          }
        }
      }
    }
  }
  fclose(augerfile);
}


static void read_collion_data(void)
{
  printout("Reading collisional ionization data...\n");

  FILE *cifile = fopen_required("collion.txt", "r");

  fscanf(cifile, "%d", &colliondatacount);
  printout("Reading %d collisional transition rows\n", colliondatacount);
  colliondata = (struct collionrow *) calloc(colliondatacount, sizeof(struct collionrow));
  int n = 0; // the index of kept rows, skipping rows that aren't in the simulation
  for (int i = 0; i < colliondatacount; i++)
  {
    fscanf(cifile, "%2d %2d %1d %1d %lg %lg %lg %lg %lg",
           &colliondata[n].Z, &colliondata[n].nelec, &colliondata[n].n, &colliondata[n].l,
           &colliondata[n].ionpot_ev, &colliondata[n].A, &colliondata[n].B, &colliondata[n].C, &colliondata[n].D);

    colliondata[n].prob_num_auger[0] = 1.;
    for (int a = 1; a <= NT_MAX_AUGER_ELECTRONS; a++)
    {
      colliondata[n].prob_num_auger[a] = 0.;
    }

    colliondata[n].auger_g_accumulated = 0.;
    colliondata[n].en_auger_ev = 0.;
    colliondata[n].n_auger_elec_avg = 0.;

    const int ionstage = colliondata[n].Z - colliondata[n].nelec + 1;
    const int element = get_elementindex(colliondata[n].Z);
    if (element >= 0 && get_ionstage(element, 0) <= ionstage && ionstage < (get_ionstage(element, 0) + get_nions(element)))
    {
      // printout("Including ionisation data for Z=%d ionstage %d\n", colliondata[n].Z, ionstage);
      n++;
    }
    else
    {
      // printout("Excluding ionisation data for Z=%d ionstage %d\n", colliondata[n].Z, ionstage);
    }

    // printout("ci row: %2d %2d %1d %1d %lg %lg %lg %lg %lg\n",
    //        colliondata[n].Z, colliondata[n].nelec, colliondata[n].n, colliondata[n].l, colliondata[n].ionpot_ev,
    //        colliondata[n].A, colliondata[n].B, colliondata[n].C, colliondata[n].D);
  }
  printout("Stored %d of %d input shell cross sections\n", n, colliondatacount);
  colliondatacount = n;
  colliondata = (struct collionrow *) realloc(colliondata, colliondatacount * sizeof(struct collionrow));
  if (colliondata == NULL)
  {
    printout("Could not reallocate colliondata.\n");
    abort();
  }

  fclose(cifile);

  if (NT_MAX_AUGER_ELECTRONS > 0)
  {
    read_auger_data();
  }
}


static void zero_all_effionpot(const int modelgridindex)
{
  assert(nt_solution[modelgridindex].prob_num_auger);
  assert(nt_solution[modelgridindex].ionenfrac_num_auger);

  for (int uniqueionindex = 0; uniqueionindex < globals::includedions; uniqueionindex++)
  {
    nt_solution[modelgridindex].eff_ionpot[uniqueionindex] = 0.;

    nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1)] = 1.;
    nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1)] = 1.;
    for (int a = 1; a <= NT_MAX_AUGER_ELECTRONS; a++)
    {
      nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 0.;
      nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 0.;
    }

    int element = 0;
    int ion = 0;
    get_ionfromuniqueionindex(uniqueionindex, &element, &ion);
    assert(fabs(get_auger_probability(modelgridindex, element, ion, 0) - 1.0) < 1e-3);
    assert(fabs(get_ion_auger_enfrac(modelgridindex, element, ion, 0) - 1.0) < 1e-3);
  }
  check_auger_probabilities(modelgridindex);
}


void init(const int my_rank)
{
  assert(nonthermal_initialized == false);
  nonthermal_initialized = true;

  deposition_rate_density = (double *) calloc(globals::npts_model, sizeof(double));
  deposition_rate_density_timestep = (int *) calloc(globals::npts_model, sizeof(int));

  for (int modelgridindex = 0; modelgridindex < globals::npts_model; modelgridindex++)
  {
    deposition_rate_density[modelgridindex] = -1.;
    deposition_rate_density_timestep[modelgridindex] = -1;
  }

  if (!NT_ON)
  {
    return;
  }

  read_binding_energies();

  if (!NT_SOLVE_SPENCERFANO)
  {
    return;
  }


  printout("Initializing non-thermal solver with:\n");
  printout("  NT_EXCITATION %s\n", NT_EXCITATION_ON ? "on" : "off");
  printout("  MAX_NT_EXCITATIONS_STORED %d\n", MAX_NT_EXCITATIONS_STORED);
  printout("  NTEXCITATION_MAXNLEVELS_LOWER %d\n", NTEXCITATION_MAXNLEVELS_LOWER);
  printout("  NTEXCITATION_MAXNLEVELS_UPPER %d\n", NTEXCITATION_MAXNLEVELS_UPPER);
  printout("  SFPTS %d\n", SFPTS);
  printout("  SF_EMIN %g eV\n", SF_EMIN);
  printout("  SF_EMAX %g eV\n", SF_EMAX);
  printout("  SF_USE_LOG_E_INCREMENT %s\n", SF_USE_LOG_E_INCREMENT ? "on" : "off");
  printout("  STORE_NT_SPECTRUM %s\n", STORE_NT_SPECTRUM ? "on" : "off");
  printout("  NT_USE_VALENCE_IONPOTENTIAL %s\n", NT_USE_VALENCE_IONPOTENTIAL ? "on" : "off");
  printout("  NT_MAX_AUGER_ELECTRONS %d\n", NT_MAX_AUGER_ELECTRONS);
  printout("  SF_AUGER_CONTRIBUTION %s\n", SF_AUGER_CONTRIBUTION_ON ? "on" : "off");
  printout("  SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN %s\n", SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN ? "on" : "off");

  char filename[100];
  sprintf(filename,"nonthermalspec_%.4d.out", my_rank);
  nonthermalfile = fopen_required(filename, "w");
  fprintf(nonthermalfile,"%8s %15s %8s %11s %11s %11s\n",
          "timestep","modelgridindex","index","energy_ev","source","y");
  fflush(nonthermalfile);

  nt_solution = (struct nt_solution_struct *) calloc(globals::npts_model, sizeof(struct nt_solution_struct));

  long mem_usage_yfunc = 0;
  for (int modelgridindex = 0; modelgridindex < globals::npts_model; modelgridindex++)
  {
    // should make these negative?
    nt_solution[modelgridindex].frac_heating = 0.97;
    nt_solution[modelgridindex].frac_ionization = 0.03;
    nt_solution[modelgridindex].frac_excitation = 0.0;

    nt_solution[modelgridindex].nneperion_when_solved = -1.;
    nt_solution[modelgridindex].timestep_last_solved = -1;

    if (get_numassociatedcells(modelgridindex) > 0)
    {
      nt_solution[modelgridindex].eff_ionpot = (float *) calloc(globals::includedions, sizeof(float));
      nt_solution[modelgridindex].fracdep_ionization_ion = (double *) calloc(globals::includedions, sizeof(double));

      nt_solution[modelgridindex].prob_num_auger = (float *) calloc(globals::includedions * (NT_MAX_AUGER_ELECTRONS + 1), sizeof(float));
      nt_solution[modelgridindex].ionenfrac_num_auger = (float *) calloc(globals::includedions * (NT_MAX_AUGER_ELECTRONS + 1), sizeof(float));

      if (STORE_NT_SPECTRUM)
      {
        nt_solution[modelgridindex].yfunc = (double *) calloc(SFPTS, sizeof(double));
        assert(nt_solution[modelgridindex].yfunc != NULL);
        mem_usage_yfunc += SFPTS * sizeof(double);
      }

      zero_all_effionpot(modelgridindex);
    }
    else
    {
      nt_solution[modelgridindex].eff_ionpot = NULL;
      nt_solution[modelgridindex].fracdep_ionization_ion = NULL;

      nt_solution[modelgridindex].prob_num_auger = NULL;
      nt_solution[modelgridindex].ionenfrac_num_auger = NULL;

      nt_solution[modelgridindex].yfunc = NULL;
    }

    nt_solution[modelgridindex].frac_excitations_list = NULL;
    nt_solution[modelgridindex].frac_excitations_list_size = 0;
  }

  if (STORE_NT_SPECTRUM)
    printout("mem_usage: storing non-thermal spectra for all allocated cells occupies %.3f MB\n", mem_usage_yfunc / 1024 / 1024.);;

  envec = gsl_vector_calloc(SFPTS); // energy grid on which solution is sampled
  logenvec = gsl_vector_calloc(SFPTS);
  sourcevec = gsl_vector_calloc(SFPTS); // energy grid on which solution is sampled

  #if (SF_USE_LOG_E_INCREMENT)
  {
    delta_envec = gsl_vector_calloc(SFPTS);
    delta_log_e = log(SF_EMAX / SF_EMIN) / (SFPTS - 1);
  }
  #endif

  #if (SF_USE_LOG_E_INCREMENT)
    const double source_spread_en_target = ceil(SFPTS * 0.03333) * (SF_EMAX - SF_EMIN) / (SFPTS - 1);
    const int sourcelowerindex = round(log((SF_EMAX - source_spread_en_target) / SF_EMIN) / delta_log_e);
    // source_spread_en will be close to source_spread_en_target but lie exactly on an energy grid point
    const double source_spread_en = SF_EMAX - SF_EMIN * exp(sourcelowerindex * delta_log_e);
  #else
    // const int source_spread_pts = ceil(SFPTS / 20);
    const int source_spread_pts = ceil(SFPTS * 0.03333); // KF92 OXYGEN TEST
    const double source_spread_en = source_spread_pts * DELTA_E;
    const int sourcelowerindex = SFPTS - source_spread_pts;
  #endif

  for (int s = 0; s < SFPTS; s++)
  {
    #if (SF_USE_LOG_E_INCREMENT)
      const double energy_ev = SF_EMIN * exp(delta_log_e * s);
      const double energy_ev_next = SF_EMIN * exp(delta_log_e * (s + 1));

      gsl_vector_set(delta_envec, s, energy_ev_next - energy_ev);
    #else
      const double energy_ev = SF_EMIN + s * DELTA_E;
    #endif

    gsl_vector_set(envec, s, energy_ev);
    gsl_vector_set(logenvec, s, log(energy_ev));

    // spread the source over some energy width
    if (s < sourcelowerindex)
      gsl_vector_set(sourcevec, s, 0.);
    else
      gsl_vector_set(sourcevec, s, 1. / source_spread_en);
  }

  gsl_vector *integralvec = gsl_vector_alloc(SFPTS);
  gsl_vector_memcpy(integralvec, sourcevec);
  #if (SF_USE_LOG_E_INCREMENT)
    gsl_vector_mul(integralvec, delta_envec);
  #else
    gsl_vector_scale(integralvec, DELTA_E);
  #endif
  const double sourceintegral = gsl_blas_dasum(integralvec); // integral of S(e) dE

  gsl_vector_mul(integralvec, envec);
  E_init_ev = gsl_blas_dasum(integralvec); // integral of E * S(e) dE
  gsl_vector_free(integralvec);

  // or put all of the source into one point at SF_EMAX
  // gsl_vector_set_zero(sourcevec);
  // gsl_vector_set(sourcevec, SFPTS - 1, 1 / DELTA_E);
  // E_init_ev = SF_EMAX;

  printout("E_init: %14.7e eV\n", E_init_ev);
  printout("source function integral: %14.7e\n", sourceintegral);

  read_collion_data();

  nonthermal_initialized = true;
  printout("Finished initializing non-thermal solver\n");
}


void calculate_deposition_rate_density(const int modelgridindex, const int timestep)
// should be in erg / s / cm^3
{
  const double gamma_deposition = globals::rpkt_emiss[modelgridindex] * 1.e20 * FOURPI;

  const double tmid = globals::time_step[timestep].mid;
  const double positron_deposition = decay::get_positroninjection_rate_density(modelgridindex, tmid);

  deposition_rate_density[modelgridindex] = gamma_deposition + positron_deposition;

  printout("nt_deposition_rate(mgi %d timestep %d): gammadep %9.2f, posdep %9.2f eV/s/cm^3\n", modelgridindex, timestep, gamma_deposition / EV, positron_deposition / EV);

  deposition_rate_density_timestep[modelgridindex] = timestep;
}


double get_deposition_rate_density(const int modelgridindex)
// should be in erg / s / cm^3
{
  // if (deposition_rate_density[modelgridindex] <= 0)
  // {
  //   calculate_deposition_rate_density(modelgridindex, nts_global);
  //   printout("No deposition_rate_density for cell %d. Calculated value of %g has been stored.\n",
  //            modelgridindex, deposition_rate_density[modelgridindex]);
  // }
  assert(deposition_rate_density_timestep[modelgridindex] == globals::nts_global);
  assert(deposition_rate_density[modelgridindex] >= 0);
  return deposition_rate_density[modelgridindex];
}


static double get_y_sample(const int modelgridindex, const int index)
{
  if (nt_solution[modelgridindex].yfunc != NULL)
  {
    if (!isfinite(nt_solution[modelgridindex].yfunc[index]))
    {
      printout("get_y_sample index %d %g\n", index, nt_solution[modelgridindex].yfunc[index]);
    }
    return nt_solution[modelgridindex].yfunc[index];
  }
  else
  {
    printout("non-thermal: attempted to get y function sample index %d in cell %d, but the y array pointer is null\n",
             index, modelgridindex);
    abort();
  }
}


static void nt_write_to_file(const int modelgridindex, const int timestep, const int iteration)
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

  static long nonthermalfile_offset_iteration_zero = 0;
  # ifdef _OPENMP
  # pragma omp threadprivate(nonthermalfile_offset_iteration_zero)
  # endif
  {
    if (iteration == 0)
    {
      nonthermalfile_offset_iteration_zero = ftell(nonthermalfile);
    }
    else
    {
      // overwrite the non-thermal spectrum of a previous iteration of the same timestep and gridcell
      fseek(nonthermalfile, nonthermalfile_offset_iteration_zero, SEEK_SET);
    }
  }

#ifndef yscalefactoroverride // manual override can be defined
  const double yscalefactor = (get_deposition_rate_density(modelgridindex) / (E_init_ev * EV));
#else
  const double yscalefactor = yscalefactoroverride(modelgridindex);
#endif

  for (int s = 0; s < SFPTS; s++)
  {
    fprintf(nonthermalfile, "%8d %15d %8d %11.5e %11.5e %11.5e\n",
            timestep, modelgridindex, s, gsl_vector_get(envec,s),
            gsl_vector_get(sourcevec, s), yscalefactor * get_y_sample(modelgridindex, s));
  }
  fflush(nonthermalfile);
# ifdef _OPENMP
  }
# endif
}


void close_file(void)
{
  nonthermal_initialized = false;

  free(deposition_rate_density);
  free(deposition_rate_density_timestep);

  if (!NT_ON || !NT_SOLVE_SPENCERFANO)
    return;

  fclose(nonthermalfile);
  gsl_vector_free(envec);
  gsl_vector_free(sourcevec);
  for (int modelgridindex = 0; modelgridindex < globals::npts_model; modelgridindex++)
  {
    if (get_numassociatedcells(modelgridindex) > 0)
    {
      if (STORE_NT_SPECTRUM)
      {
        free(nt_solution[modelgridindex].yfunc);
      }
      free(nt_solution[modelgridindex].fracdep_ionization_ion);
      free(nt_solution[modelgridindex].eff_ionpot);
      free(nt_solution[modelgridindex].prob_num_auger);
      free(nt_solution[modelgridindex].ionenfrac_num_auger);
      if (nt_solution[modelgridindex].frac_excitations_list_size > 0)
      {
        free(nt_solution[modelgridindex].frac_excitations_list);
      }
    }
  }
  free(colliondata);
}


inline
static int get_energyindex_ev_lteq(const double energy_ev)
// finds the highest energy point <= energy_ev
{
  #if (SF_USE_LOG_E_INCREMENT)
  const int index = floor(log(energy_ev / SF_EMIN) / delta_log_e);
  #else
  const int index = floor((energy_ev - SF_EMIN) / DELTA_E);
  #endif

  if (index < 0)
    return 0;
  else if (index > SFPTS - 1)
    return SFPTS - 1;
  else
    return index;
}


inline
static int get_energyindex_ev_gteq(const double energy_ev)
// finds the highest energy point <= energy_ev
{
  #if (SF_USE_LOG_E_INCREMENT)
  const int index = ceil(log(energy_ev / SF_EMIN) / delta_log_e);
  #else
  const int index = ceil((energy_ev - SF_EMIN) / DELTA_E);
  #endif

  if (index < 0)
    return 0;
  else if (index > SFPTS - 1)
    return SFPTS - 1;
  else
    return index;
}


static double get_y(const int modelgridindex, const double energy_ev)
{
  if (energy_ev <= 0)
    return 0.;

  #if (SF_USE_LOG_E_INCREMENT)
    const int index = log(energy_ev / SF_EMIN) / delta_log_e;
  #else
    const int index = (energy_ev - SF_EMIN) / DELTA_E;
  #endif

  // assert(index > 0);
  if (index < 0)
  {
    // return 0.;
    assert(isfinite(get_y_sample(modelgridindex, 0)));
    return get_y_sample(modelgridindex, 0);
  }
  else if (index > SFPTS - 1)
    return 0.;
  else
  {
    const double enbelow = gsl_vector_get(envec, index);
    const double enabove = gsl_vector_get(envec, index + 1);
    const double ybelow = get_y_sample(modelgridindex, index);
    const double yabove = get_y_sample(modelgridindex, index + 1);
    const double x = (energy_ev - enbelow) / (enabove - enbelow);
    return (1 - x) * ybelow + x * yabove;

    // or return the nearest neighbour
    // return get_y_sample(modelgridindex, index);
  }
}


static double electron_loss_rate(const double energy, const double nne)
// -dE / dx for fast electrons
// energy is in ergs
// nne is the thermal electron density [cm^-3]
// return units are erg / cm
{
  if (energy <= 0.)
    return 0;
  const double omegap = sqrt(4 * PI * nne * pow(QE, 2) / ME);
  const double zetae = H * omegap / 2 / PI;
  if (energy > 14 * EV)
  {
    return nne * 2 * PI * pow(QE, 4) / energy * log(2 * energy / zetae);
  }
  else
  {
    const double v = sqrt(2 * energy / ME);
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
  else if (!globals::linelist[lineindex].forbidden)
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


static int get_xs_excitation_vector(gsl_vector *const xs_excitation_vec, const int lineindex, const double statweight_lower, const double epsilon_trans)
// vector of collisional excitation cross sections in cm^2
// epsilon_trans is in erg
// returns the index of the first valid cross section point (en >= epsilon_trans)
// all elements below this index are invalid and should not be used
{
  const double coll_str = get_coll_str(lineindex);

  if (coll_str >= 0)
  {
    // collision strength is available, so use it
    // Li et al. 2012 equation 11
    const double constantfactor = pow(H_ionpot, 2) / statweight_lower * coll_str * PI * A_naught_squared;

    const int en_startindex = get_energyindex_ev_gteq(epsilon_trans / EV);

    for (int j = 0; j < en_startindex; j++)
      gsl_vector_set(xs_excitation_vec, j, 0.);

    for (int j = en_startindex; j < SFPTS; j++)
    {
      const double energy = gsl_vector_get(envec, j) * EV;
      gsl_vector_set(xs_excitation_vec, j, constantfactor * pow(energy, -2));
    }
    return en_startindex;
  }
  else if (!globals::linelist[lineindex].forbidden)
  {
    const double fij = osc_strength(lineindex);
    // permitted E1 electric dipole transitions

    // const double g_bar = 0.2;
    const double A = 0.28;
    const double B = 0.15;

    const double prefactor = 45.585750051; // 8 * pi^2/sqrt(3)
    const double epsilon_trans_ev = epsilon_trans / EV;

    // Eq 4 of Mewe 1972, possibly from Seaton 1962?
    const double constantfactor = epsilon_trans_ev * prefactor * A_naught_squared * pow(H_ionpot / epsilon_trans, 2) * fij;

    const int en_startindex = get_energyindex_ev_gteq(epsilon_trans_ev);

    for (int j = 0; j < en_startindex; j++)
      gsl_vector_set(xs_excitation_vec, j, 0.);

    // U = en / epsilon
    // g_bar = A * log(U) + b
    // xs[j] = constantfactor * g_bar / envec[j]

    for (int j = en_startindex; j < SFPTS; j++)
    {
      const double logU = gsl_vector_get(logenvec, j) - log(epsilon_trans_ev);
      const double g_bar = A * logU + B;
      gsl_vector_set(xs_excitation_vec, j, constantfactor * g_bar / gsl_vector_get(envec, j));
    }

    return en_startindex;
  }
  else
  {
    // gsl_vector_set_zero(xs_excitation_vec);
    return -1;
  }
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


static int get_xs_ionization_vector(gsl_vector *const xs_vec, const int collionindex)
// xs_vec will be set with impact ionization cross sections for E > ionpot_ev (and zeros below this energy)
{
  const double ionpot_ev = colliondata[collionindex].ionpot_ev;
  const int startindex = get_energyindex_ev_gteq(ionpot_ev);

  // en points for which en < ionpot
  for (int i = 0; i < startindex; i++)
  {
    gsl_vector_set(xs_vec, i, 0.);
  }

  const double A = colliondata[collionindex].A;
  const double B = colliondata[collionindex].B;
  const double C = colliondata[collionindex].C;
  const double D = colliondata[collionindex].D;

  for (int i = startindex; i < SFPTS; i++)
  {
    const double u = gsl_vector_get(envec, i) / ionpot_ev;
    const double xs_ioniz = 1e-14 * (A * (1 - 1/u) + B * pow((1 - 1/u), 2) + C * log(u) + D * log(u) / u) / (u * pow(ionpot_ev, 2));
    gsl_vector_set(xs_vec, i, xs_ioniz);
  }

  return startindex;
}


static double Psecondary(const double e_p, const double epsilon, const double I, const double J)
// distribution of secondary electron energies for primary electron with energy e_p
// Opal, Peterson, & Beaty (1971)
{
  const double e_s = epsilon - I;

  if (e_p <= I || e_s < 0.)
  {
    return 0.;
  }
  assert(J > 0);
  assert(e_p >= I);
  assert(e_s >= 0);
  assert(isfinite(atan((e_p - I) / 2 / J)));
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
// not valid for energy > SF_EMIN
{
  const double energy_ev = energy / EV;
  const double tot_nion = get_tot_nion(modelgridindex);
  double N_e = 0.;

  for (int element = 0; element < get_nelements(); element++)
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

      const int nlevels_all = get_nlevels(element, ion);
      const int nlevels = (nlevels_all > NTEXCITATION_MAXNLEVELS_LOWER) ? NTEXCITATION_MAXNLEVELS_LOWER : nlevels_all;

      for (int lower = 0; lower < nlevels; lower++)
      {
        const int nuptrans = get_nuptrans(element, ion, lower);
        const double nnlevel = calculate_exclevelpop(modelgridindex, element, ion, lower);
        const double epsilon_lower = epsilon(element, ion, lower);
        for (int t = 0; t < nuptrans; t++)
        {
          const int lineindex = globals::elements[element].ions[ion].levels[lower].uptrans_lineindicies[t];
          const int upper = globals::linelist[lineindex].upperlevelindex;
          if (upper >= NTEXCITATION_MAXNLEVELS_UPPER)
          {
            continue;
          }
          const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
          const double epsilon_trans_ev = epsilon_trans / EV;
          N_e_ion += (nnlevel / nnion) * get_y(modelgridindex, energy_ev + epsilon_trans_ev) * xs_excitation(lineindex, epsilon_trans, energy + epsilon_trans);
        }
      }

      // ionization terms
      for (int n = 0; n < colliondatacount; n++)
      {
        if (colliondata[n].Z == Z && colliondata[n].nelec == Z - ionstage + 1)
        {
          const double ionpot_ev = colliondata[n].ionpot_ev;
          const double J = get_J(Z, ionstage, ionpot_ev);
          const double lambda = fmin(SF_EMAX - energy_ev, energy_ev + ionpot_ev);

          const int integral1startindex = get_energyindex_ev_lteq(ionpot_ev);
          const int integral1stopindex = get_energyindex_ev_lteq(lambda);

          // integral from ionpot up to lambda
          for (int i = integral1startindex; i <= integral1stopindex; i++)
          {
            double endash = gsl_vector_get(envec, i);
            #if (SF_USE_LOG_E_INCREMENT)
            const double delta_endash = gsl_vector_get(delta_envec, i);
            #else
            const double delta_endash = DELTA_E;
            #endif

            N_e_ion += get_y(modelgridindex, energy_ev + endash) * xs_impactionization(energy_ev + endash, n) * Psecondary(energy_ev + endash, endash, ionpot_ev, J) * delta_endash;
          }

          // integral from 2E + I up to E_max
          const int integral2startindex = get_energyindex_ev_lteq(2 * energy_ev + ionpot_ev);
          for (int i = integral2startindex; i < SFPTS; i++)
          {
            double endash = gsl_vector_get(envec, i);
            #if (SF_USE_LOG_E_INCREMENT)
            const double delta_endash = gsl_vector_get(delta_envec, i);
            #else
            const double delta_endash = DELTA_E;
            #endif
            N_e_ion += get_y_sample(modelgridindex, i) * xs_impactionization(endash, n) * Psecondary(endash, energy_ev + ionpot_ev, ionpot_ev, J) * delta_endash;
          }
        }
      }

      N_e += nnion * N_e_ion;
    }
  }

  // source term, should be zero at the low end anyway
  N_e += gsl_vector_get(sourcevec, get_energyindex_ev_lteq(energy_ev));

  assert(isfinite(N_e));
  return N_e;
}


static float calculate_frac_heating(const int modelgridindex)
// Kozma & Fransson equation 3
{
  // frac_heating multiplied by E_init, which will be divided out at the end
  double frac_heating_Einit = 0.;

  const float nne = get_nne(modelgridindex);
  // const float nnetot = get_nnetot(modelgridindex);

  const double E_0 = SF_EMIN;

  const int startindex = get_energyindex_ev_lteq(E_0);
  for (int i = startindex; i < SFPTS; i++)
  {
    const double endash = gsl_vector_get(envec, i);
    #if (SF_USE_LOG_E_INCREMENT)
    double deltaendash = gsl_vector_get(delta_envec, i);
    #else
    double deltaendash = DELTA_E;
    #endif

    if (i == startindex && endash < E_0)
    {
      deltaendash = endash + deltaendash - E_0;
    }
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
    printout("WARNING: calculate_frac_heating: invalid result of %g. Setting to 1.0 instead\n", frac_heating);
    return 1.0;
  }

  return frac_heating;
}


float get_nt_frac_heating(const int modelgridindex)
{
  if (!NT_ON)
    return 1.;
  else if (!NT_SOLVE_SPENCERFANO)
    return 0.97;
  else
  {
    const float frac_heating = nt_solution[modelgridindex].frac_heating;
    // add any debugging checks here?
    return frac_heating;
  }
}


static float get_nt_frac_ionization(const int modelgridindex)
{
  if (!NT_ON)
    return 0.;
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


static float get_nt_frac_excitation(const int modelgridindex)
{
  if (!NT_ON || !NT_SOLVE_SPENCERFANO)
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
        const double use3 = globals::elements[element].ions[ion].ionpot;
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
  for (int ielement = 0; ielement < get_nelements(); ielement++)
  {
    Zbar += globals::modelgrid[modelgridindex].composition[ielement].abundance * globals::elements[ielement].anumber;
  }
  //printout("cell %d has Zbar of %g\n", modelgridindex, Zbar);

  const double Aconst = 1.33e-14 * EV * EV;
  const double binding = get_mean_binding_energy(element, ion);
  const double oneoverW = Aconst * binding / Zbar / (2 * 3.14159 * pow(QE, 4));
  //printout("For element %d ion %d I got W of %g (eV)\n", element, ion, 1./oneoverW/EV);

  return oneoverW;
}


static double calculate_nt_frac_ionization_shell(const int modelgridindex, const int element, const int ion, const int collionindex)
// the fraction of deposition energy that goes into ionising electrons in this particular shell
{
  const double nnion = ionstagepop(modelgridindex, element, ion); // hopefully ions per cm^3?
  const double ionpot_ev = colliondata[collionindex].ionpot_ev;

  gsl_vector *cross_section_vec = gsl_vector_alloc(SFPTS);
  get_xs_ionization_vector(cross_section_vec, collionindex);

  // either multiply by the variable delta_e for LOG_E spacing...
  #if (SF_USE_LOG_E_INCREMENT)
  gsl_vector_mul(cross_section_vec, delta_envec);
  #endif

  gsl_vector_view yvecview = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);

  double y_dot_crosssection_de = 0.;
  gsl_blas_ddot(&yvecview.vector, cross_section_vec, &y_dot_crosssection_de);
  gsl_vector_free(cross_section_vec);

  // or multiply the scalar result by the constant DELTA_E
  #if (!SF_USE_LOG_E_INCREMENT)
  y_dot_crosssection_de *= DELTA_E;
  #endif

  return nnion * ionpot_ev * y_dot_crosssection_de / E_init_ev;
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


static double calculate_nt_ionization_ratecoeff(
  const int modelgridindex, const int element, const int ion, const bool assumeshellpotentialisvalence)
// Integrate the ionization cross section over the electron degradation function to get the ionization rate coefficient
// i.e. multiply this by ion population to get a rate of ionizations per second
// Do not call during packet propagation, as the y vector may not be in memory!
// IMPORTANT: we are dividing by the shell potential, not the valence potential here!
// To change this set assumeshellpotentialisvalence to true
{
  gsl_vector *cross_section_vec = gsl_vector_alloc(SFPTS);
  gsl_vector *cross_section_vec_allshells = gsl_vector_calloc(SFPTS);

  const int Z = get_element(element);
  const int ionstage = get_ionstage(element, ion);
  double ionpot_valence = -1;

  for (int collionindex = 0; collionindex < colliondatacount; collionindex++)
  {
    if (colliondata[collionindex].Z == Z && colliondata[collionindex].nelec == Z - ionstage + 1)
    {
      get_xs_ionization_vector(cross_section_vec, collionindex);

      if (assumeshellpotentialisvalence)
      {
        const double ionpot_shell = colliondata[collionindex].ionpot_ev * EV;
        if (ionpot_valence < 0)
          ionpot_valence = ionpot_shell;

        // ensure that the first shell really was the valence shell (we assumed ascending energy order)
        assert(ionpot_shell >= ionpot_valence);

        // boost the ionization rate by assuming shell vacancy energy is used to eject valence electrons
        gsl_vector_scale(cross_section_vec, ionpot_shell / ionpot_valence);
      }

      gsl_vector_add(cross_section_vec_allshells, cross_section_vec);
    }
  }

  gsl_vector_free(cross_section_vec);

  #if (SF_USE_LOG_E_INCREMENT)
  gsl_vector_mul(cross_section_vec_allshells, delta_envec);
  #endif

  assert(nt_solution[modelgridindex].yfunc != NULL);

  double y_dot_crosssection_de = 0.;
  gsl_vector_view yvecview_thismgi = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);
  gsl_blas_ddot(&yvecview_thismgi.vector, cross_section_vec_allshells, &y_dot_crosssection_de);
  gsl_vector_free(cross_section_vec_allshells);

  #if (!SF_USE_LOG_E_INCREMENT)
  y_dot_crosssection_de *= DELTA_E;
  #endif

  const double deposition_rate_density_ev = get_deposition_rate_density(modelgridindex) / EV;
  const double yscalefactor = deposition_rate_density_ev / E_init_ev;

  return yscalefactor * y_dot_crosssection_de;
}


static void calculate_eff_ionpot_auger_rates(
  const int modelgridindex, const int element, const int ion)
// Kozma & Fransson 1992 equation 12, except modified to be a sum over all shells of an ion
// the result is in ergs
{
  const int Z = get_element(element);
  const int ionstage = get_ionstage(element, ion);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  const double nnion = ionstagepop(modelgridindex, element, ion); // ions/cm^3
  const double tot_nion = get_tot_nion(modelgridindex);
  const double X_ion = nnion / tot_nion; // molar fraction of this ion

  // The ionization rates of all shells of an ion add to make the ion's total ionization rate,
  // i.e., Gamma_ion = Gamma_shell_a + Gamma_shell_b + ...
  // And since the ionization rate is inversely proportional to the effective ion potential,
  // we solve:
  // (eta_ion / ionpot_ion) = (eta_shell_a / ionpot_shell_a) + (eta_shell_b / ionpot_shell_b) + ...
  // where eta is the fraction of the deposition energy going into ionization of the ion or shell

  double eta_nauger_ionize_over_ionpot_sum[NT_MAX_AUGER_ELECTRONS + 1];
  double eta_nauger_ionize_sum[NT_MAX_AUGER_ELECTRONS + 1];

  for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
  {
    eta_nauger_ionize_over_ionpot_sum[a] = 0.;
    nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 0.;

    eta_nauger_ionize_sum[a] = 0.;
    nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 0.;
  }

  double eta_over_ionpot_sum = 0.;
  double eta_sum = 0.;
  double ionpot_valence = -1;
  int matching_nlsubshell_count = 0;
  for (int collionindex = 0; collionindex < colliondatacount; collionindex++)
  {
    if (colliondata[collionindex].Z == Z && colliondata[collionindex].nelec == Z - ionstage + 1)
    {
      matching_nlsubshell_count++;
      const double frac_ionization_shell = calculate_nt_frac_ionization_shell(modelgridindex, element, ion, collionindex);
      eta_sum += frac_ionization_shell;
      const double ionpot_shell = colliondata[collionindex].ionpot_ev * EV;

      if (ionpot_valence < 0)
        ionpot_valence = ionpot_shell;

      // ensure that the first shell really was the valence shell (we assumed ascending energy order)
      assert(ionpot_shell >= ionpot_valence);

      const double ionpot = NT_USE_VALENCE_IONPOTENTIAL ? ionpot_valence : ionpot_shell;
      const double eta_over_ionpot = frac_ionization_shell / ionpot; // this is proportional to rate

      eta_over_ionpot_sum += eta_over_ionpot;

      for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
      {
        eta_nauger_ionize_over_ionpot_sum[a] += eta_over_ionpot * colliondata[collionindex].prob_num_auger[a];
        // printout("test Z=%d ion %d shell %d prob_num_auger(%d) = %g, eta_over_ionpot %g, product %g\n", get_element(element), get_ionstage(element, ion), collionindex, a, colliondata[collionindex].prob_num_auger[a], eta_over_ionpot, eta_over_ionpot * colliondata[collionindex].prob_num_auger[a]);
        eta_nauger_ionize_sum[a] += frac_ionization_shell * colliondata[collionindex].prob_num_auger[a];
      }
    }
  }

  if (NT_MAX_AUGER_ELECTRONS > 0 && matching_nlsubshell_count > 0)
  {
    const int nions = get_nions(element);
    if (ion < nions - 1) // don't try to ionise the top ion
    {
      for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
      {
        // printout("test2 Z=%d ion %d a %d probability %g\n", get_element(element), get_ionstage(element, ion), a, eta_nauger_ionize_over_ionpot_sum[a] / eta_over_ionpot_sum);
        if (ion + 1 + a < nions) // not too many Auger electrons to exceed the top ion of this element
        {
          nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = eta_nauger_ionize_over_ionpot_sum[a] / eta_over_ionpot_sum;
          nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = eta_nauger_ionize_sum[a] / eta_sum;
        }
        else
        {
          // the following ensures that multiple ionisations can't send you to an ion stage that is not in the model
          // could send it to the top one with a = nions - 1 - ion - 1

          nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + nions - 1 - ion - 1] += eta_nauger_ionize_over_ionpot_sum[a] / eta_over_ionpot_sum;
          nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + nions - 1 - ion - 1] += eta_nauger_ionize_sum[a] / eta_sum;

          // printout("test2b going to Z=%d ion %d a %d with new probability %g\n", get_element(element), get_ionstage(element, ion), nions - 1 - ion - 1,  nt_solution[modelgridindex].prob_num_auger[element][ion][nions - 1 - ion - 1]);

          nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 0;
          nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 0.;
        }
      }
    }

    // the following ensures that multiple ionisations can't send you to an ion stage that is not in the model
    // for (int a = NT_MAX_AUGER_ELECTRONS; a > 0; a--)
    // {
    //   if ((ion + a + 1) >= nions)
    //   {
    //     nt_solution[modelgridindex].prob_num_auger[element][ion][a - 1] += nt_solution[modelgridindex].prob_num_auger[element][ion][a];
    //     nt_solution[modelgridindex].prob_num_auger[element][ion][a] = 0.;
    //   }
    // }
  }
  else
  {
    const int a = 0;
    nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 1.;
    nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 1.;
  }

  if (matching_nlsubshell_count > 0)
  {
    double eff_ionpot = X_ion / eta_over_ionpot_sum;
    if (!isfinite(eff_ionpot))
      eff_ionpot = 0.;
    nt_solution[modelgridindex].eff_ionpot[get_uniqueionindex(element, ion)] = eff_ionpot;
  }
  else
  {
    printout("WARNING! No matching subshells in NT impact ionisation cross section data for Z=%d ionstage %d.\n -> Defaulting to work function approximation and ionisation energy is not accounted for in Spencer-Fano solution.\n",
             get_element(element), get_ionstage(element, ion));

    nt_solution[modelgridindex].eff_ionpot[get_uniqueionindex(element, ion)] = 1. / get_oneoverw(element, ion, modelgridindex);
  }
}


static float get_eff_ionpot(const int modelgridindex, const int element, int const ion)
// get the effective ion potential from the stored value
// a value of 0. should be treated as invalid
{
  return nt_solution[modelgridindex].eff_ionpot[get_uniqueionindex(element, ion)];
  // OR
  // return calculate_eff_ionpot(modelgridindex, element, ion);
}


static double nt_ionization_ratecoeff_sf(const int modelgridindex, const int element, const int ion)
// Kozma & Fransson 1992 equation 13
{
  if (get_numassociatedcells(modelgridindex) <= 0)
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


double nt_ionization_upperion_probability(
  const int modelgridindex, const int element, const int lowerion, const int upperion, const bool energyweighted)
{
  assert(upperion > lowerion);
  assert(upperion < get_nions(element));
  assert(upperion <= nt_ionisation_maxupperion(element, lowerion));
  if (NT_SOLVE_SPENCERFANO && NT_MAX_AUGER_ELECTRONS > 0)
  {
    const int numaugerelec = upperion - lowerion - 1; // number of Auger electrons to go from lowerin to upper ion
    const int uniqueionindex = get_uniqueionindex(element, lowerion);

    if (numaugerelec < NT_MAX_AUGER_ELECTRONS)
    {
      if (energyweighted)
      {
        return nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + numaugerelec];
      }
      else
      {
        return nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + numaugerelec];
      }
    }
    else if (numaugerelec == NT_MAX_AUGER_ELECTRONS)
    {
      double prob_remaining = 1.;
      for (int a = 0; a < NT_MAX_AUGER_ELECTRONS; a++)
      {
        if (energyweighted)
        {
          prob_remaining -= nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a];
        }
        else
        {
          prob_remaining -= nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a];
        }
      }
      if (energyweighted)
      {
        assert(fabs(prob_remaining - nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + numaugerelec]) < 0.001);
      }
      else
      {
        if (fabs(prob_remaining - nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + numaugerelec]) >= 0.001)
        {
          printout("Auger probabilities issue for cell %d Z=%02d ionstage %d to %d\n", modelgridindex, get_element(element), get_ionstage(element, lowerion), get_ionstage(element, upperion));
          for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
          {
            printout("  a %d prob %g\n", a, nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a]);
          }
          abort();
        }
      }
      return prob_remaining;
    }
    else
    {
      printout("WARNING: tried to ionise from Z=%02d ionstage %d to %d\n",
        get_element(element), get_ionstage(element, lowerion), get_ionstage(element, upperion));
      return 0.;
    }
  }
  else
  {
    return (upperion == lowerion + 1) ? 1.0 : 0.;
  }
}


int nt_ionisation_maxupperion(const int element, const int lowerion)
{
  const int nions = get_nions(element);
  assert(lowerion < nions - 1);
  int maxupper = lowerion + 1;

  if (NT_SOLVE_SPENCERFANO)
  {
    maxupper = lowerion + 1 + NT_MAX_AUGER_ELECTRONS;
  }

  if (maxupper > nions - 1)
  {
    maxupper = nions - 1;
  }

  return maxupper;
}


int nt_random_upperion(const int modelgridindex, const int element, const int lowerion, const bool energyweighted)
{
  const int nions = get_nions(element);
  assert(lowerion < nions - 1);
  if (NT_SOLVE_SPENCERFANO && NT_MAX_AUGER_ELECTRONS > 0)
  {
    while (true)
    {
      const double zrand = gsl_rng_uniform(rng);

      double prob_sum = 0.;
      for (int upperion = lowerion + 1; upperion <= nt_ionisation_maxupperion(element, lowerion); upperion++)
      {
        prob_sum += nt_ionization_upperion_probability(modelgridindex, element, lowerion, upperion, energyweighted);

        if (zrand <= prob_sum)
        {
          return upperion;
        }
      }

      printout("ERROR: nt_ionization_upperion_probability did not sum to more than zrand = %lg, prob_sum = %lg (Z=%d ionstage %d). Retrying with new random number.\n",
               zrand, prob_sum, get_element(element), get_ionstage(element, lowerion));
      assert(fabs(prob_sum - 1.0) < 1e-3);
    }
  }
  else
  {
    return lowerion + 1;
  }
}



double nt_ionization_ratecoeff(const int modelgridindex, const int element, const int ion)
{
  assert(NT_ON);
  assert(get_numassociatedcells(modelgridindex) > 0);

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


static double calculate_nt_excitation_ratecoeff_perdeposition(
  const int modelgridindex, const int lineindex, const double statweight_lower, const double epsilon_trans)
// Kozma & Fransson equation 9 divided by level population and epsilon_trans
{
  if (nt_solution[modelgridindex].yfunc == NULL)
  {
    printout("ERROR: Call to nt_excitation_ratecoeff with no y vector in memory.");
    abort();
}

  gsl_vector *xs_excitation_vec = gsl_vector_alloc(SFPTS);
  if (get_xs_excitation_vector(xs_excitation_vec, lineindex, statweight_lower, epsilon_trans) >= 0)
  {
    #if (SF_USE_LOG_E_INCREMENT)
    gsl_vector_mul(xs_excitation_vec, delta_envec);
    #endif

    double y_dot_crosssection = 0.;
    gsl_vector_view yvecview = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);
    gsl_blas_ddot(xs_excitation_vec, &yvecview.vector, &y_dot_crosssection);
    gsl_vector_free(xs_excitation_vec);

    #if (!SF_USE_LOG_E_INCREMENT)
    y_dot_crosssection *= DELTA_E;
    #endif

    return y_dot_crosssection / E_init_ev / EV;
  }
  else
  {
    gsl_vector_free(xs_excitation_vec);

    return 0.;
  }
}


double nt_excitation_ratecoeff(const int modelgridindex, const int element, const int ion, const int lower, const int upper, const double epsilon_trans, const int lineindex)
{
#if !NT_EXCITATION_ON
    return 0.;
#endif

  if ((lower >= NTEXCITATION_MAXNLEVELS_LOWER) ||
      (upper >= NTEXCITATION_MAXNLEVELS_UPPER))
    return 0.;

  if (get_numassociatedcells(modelgridindex) <= 0)
  {
    printout("ERROR: nt_excitation_ratecoeff called on empty cell %d\n", modelgridindex);
    abort();
  }

  // if the NT spectrum is stored, we can calculate any non-thermal excitation rate, even if
  // it didn't make the cut to be kept in the stored excitation list
  if (STORE_NT_SPECTRUM)
  {
    const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
    const double statweight_lower = stat_weight(element, ion, lower);

    const double ratecoeffperdeposition = calculate_nt_excitation_ratecoeff_perdeposition(modelgridindex, lineindex, statweight_lower, epsilon_trans);

    return ratecoeffperdeposition * deposition_rate_density;
  }

  const int list_size = nt_solution[modelgridindex].frac_excitations_list_size;

  // linear search for the lineindex
  // for (int excitationindex = 0; excitationindex < list_size; excitationindex++)
  // {
  //   if (nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex == lineindex)
  //   {
  //     const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
  //     const double ratecoeffperdeposition = nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition;
  //
  //     return ratecoeffperdeposition * deposition_rate_density;
  //   }
  // }

  // binary search, assuming the excitation list is sorted by lineindex ascending
  int low = 0;
  int high = list_size - 1;
  while (low <= high)
  {
    const int excitationindex = low + ((high - low) / 2);
    if (nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex < lineindex)
    {
      low = excitationindex + 1;
    }
    else if (nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex > lineindex)
    {
      high = excitationindex - 1;
    }
    else
    {
      const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
      const double ratecoeffperdeposition = nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition;

      return ratecoeffperdeposition * deposition_rate_density;
    }
  }

  return 0.;
}


static void select_nt_ionization(int modelgridindex, int *element, int *lowerion)
// select based on stored frac_deposition for each ion
{
  const double zrand = gsl_rng_uniform(rng);
  double frac_deposition_ion_sum = 0.;
  // zrand is between zero and frac_ionization
  // keep subtracting off deposition fractions of ionizations transitions until we hit the right one
  // e.g. if zrand was less than frac_dep_trans1, then use the first transition
  // e.g. if zrand was between frac_dep_trans1 and frac_dep_trans2 then use the second transition, etc
  for (int allionindex = 0; allionindex < globals::includedions; allionindex++)
  {
    frac_deposition_ion_sum += nt_solution[modelgridindex].fracdep_ionization_ion[allionindex];
    if (frac_deposition_ion_sum >= zrand)
    {
      get_ionfromuniqueionindex(allionindex, element, lowerion);

      return;
    }
  }
  assert(false); // should not reach here
}


static double ion_ntion_energyrate(int modelgridindex, int element, int lowerion)
{
  const double nnlowerion = ionstagepop(modelgridindex, element, lowerion);
  double enrate = 0.;
  for (int upperion = lowerion + 1; upperion <= nt_ionisation_maxupperion(element, lowerion); upperion++)
  {
    const double upperionprobfrac = nt_ionization_upperion_probability(modelgridindex, element, lowerion, upperion, false);
    // for (int lower = 0; lower < get_nlevels(element, lowerion); lower++)
    // {
    //   const double epsilon_trans = epsilon(element, upperion, 0) - epsilon(element, lowerion, lower);
    //   const double nnlower = calculate_exclevelpop(modelgridindex, element, lowerion, lower);
    //   enrate += nnlower * upperionprobfrac * epsilon_trans;
    // }
    const double epsilon_trans = epsilon(element, upperion, 0) - epsilon(element, lowerion, 0);
    enrate += nnlowerion * upperionprobfrac * epsilon_trans;
  }

  const double gamma_nt = nt_ionization_ratecoeff(modelgridindex, element, lowerion);
  return gamma_nt * enrate;
}


static double get_ntion_energyrate(int modelgridindex)
{
  double ratetotal = 0.;
  for (int ielement = 0; ielement < get_nelements(); ielement++)
  {
    const int nions = get_nions(ielement);
    for (int ilowerion = 0; ilowerion < nions - 1; ilowerion++)
    {
      ratetotal += ion_ntion_energyrate(modelgridindex, ielement, ilowerion);
    }
  }
  return ratetotal;
}


static void select_nt_ionization2(int modelgridindex, int *element, int *lowerion)
{
  const double ratetotal = get_ntion_energyrate(modelgridindex);

  const double zrand = gsl_rng_uniform(rng);
  double ratesum = 0.;
  for (int ielement = 0; ielement < get_nelements(); ielement++)
  {
    const int nions = get_nions(ielement);
    for (int ilowerion = 0; ilowerion < nions - 1; ilowerion++)
    {
      ratesum += ion_ntion_energyrate(modelgridindex, ielement, ilowerion);
      if (ratesum >= zrand * ratetotal)
      {
        *element = ielement;
        *lowerion = ilowerion;
        return;
      }
    }
  }
  assert(false);
}


void do_ntlepton(PKT *pkt_ptr)
{
  safeadd(nt_energy_deposited, pkt_ptr->e_cmf);

  const int modelgridindex = globals::cell[pkt_ptr->where].modelgridindex;

  // macroatom should not be activated in thick cells
  if (NT_ON && NT_SOLVE_SPENCERFANO && globals::modelgrid[modelgridindex].thick != 1)
  {
    // here there is some probability to cause ionisation or excitation to a macroatom packet
    // instead of converting directly to k-packet (unless the heating channel is selected)

    double zrand = gsl_rng_uniform(rng);
    // zrand is initially between [0, 1), but we will subtract off each
    // component of the deposition fractions
    // until we end and select transition_ij when zrand < dep_frac_transition_ij

    // const double frac_heating = get_nt_frac_heating(modelgridindex);

    // const double frac_ionization = get_nt_frac_ionization(modelgridindex);
    const double frac_ionization = get_ntion_energyrate(modelgridindex) / get_deposition_rate_density(modelgridindex);
    // printout("frac_ionization compare %g and %g\n", frac_ionization, get_nt_frac_ionization(modelgridindex));
    // const double frac_ionization = 0.;

    // const double frac_excitation = get_nt_frac_excitation(modelgridindex);
    const double frac_excitation = 0.;

    if (zrand < frac_ionization)
    {
      int element = -1;
      int lowerion = -1;
      // select_nt_ionization(modelgridindex, &element, &lowerion);
      select_nt_ionization2(modelgridindex, &element, &lowerion);
      const int upperion = nt_random_upperion(modelgridindex, element, lowerion, true);
      // const int upperion = lowerion + 1;

      pkt_ptr->mastate.element = element;
      pkt_ptr->mastate.ion = upperion;
      pkt_ptr->mastate.level = 0;
      pkt_ptr->mastate.activatingline = -99;
      pkt_ptr->type = TYPE_MA;
      stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_NTCOLLION);
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 20;
      pkt_ptr->trueemissiontype = -1; // since this is below zero, macroatom will set it
      pkt_ptr->trueemissionvelocity = -1;

      stats::increment(stats::COUNTER_NT_STAT_TO_IONIZATION);

      #if (TRACK_ION_STATS)
      assert(upperion < get_nions(element));
      assert(lowerion >= 0);
      const double epsilon_trans = epsilon(element, upperion, 0) - epsilon(element, lowerion, 0);
      stats::increment_ion_stats(modelgridindex, element, lowerion, stats::ION_NTION, pkt_ptr->e_cmf / epsilon_trans);
      stats::increment_ion_stats(modelgridindex, element, upperion, stats::ION_MACROATOM_ENERGYIN_NTCOLLION, pkt_ptr->e_cmf);
      #endif

      // printout("NTLEPTON packet in cell %d selected ionization of Z=%d ionstage %d to %d\n",
      //          modelgridindex, get_element(element), get_ionstage(element, lowerion), get_ionstage(element, upperion));

      return;
    }
    else if (NT_EXCITATION_ON && zrand < frac_ionization + frac_excitation)
    {
      zrand -= frac_ionization;
      // now zrand is between zero and frac_excitation
      // the selection algorithm is the same as for the ionization transitions
      const int frac_excitations_list_size = nt_solution[modelgridindex].frac_excitations_list_size;
      for (int excitationindex = 0; excitationindex < frac_excitations_list_size; excitationindex++)
      {
        const double frac_deposition_exc = nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition;
        if (zrand < frac_deposition_exc)
        {
          const int lineindex = nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex;
          const int element = globals::linelist[lineindex].elementindex;
          const int ion = globals::linelist[lineindex].ionindex;
          // const int lower = linelist[lineindex].lowerlevelindex;
          const int upper = globals::linelist[lineindex].upperlevelindex;

          pkt_ptr->mastate.element = element;
          pkt_ptr->mastate.ion = ion;
          pkt_ptr->mastate.level = upper;
          pkt_ptr->mastate.activatingline = -99;
          pkt_ptr->type = TYPE_MA;
          stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_NTCOLLEXC);
          pkt_ptr->interactions += 1;
          pkt_ptr->last_event = 21;
          pkt_ptr->trueemissiontype = -1; // since this is below zero, macroatom will set it
          pkt_ptr->trueemissionvelocity = -1;

          stats::increment(stats::COUNTER_NT_STAT_TO_EXCITATION);

          // printout("NTLEPTON packet selected in cell %d excitation of Z=%d ionstage %d level %d upperlevel %d\n",
          //          modelgridindex, get_element(element), get_ionstage(element, ion), lower, upper);

          return;
        }
        zrand -= frac_deposition_exc;
      }
      // in case we reached here because the excitation reactions that were stored didn't add up to frac_excitation_ion
      // then just convert it to a kpkt
    }
  }

  pkt_ptr->last_event = 22;
  pkt_ptr->type = TYPE_KPKT;
  stats::increment(stats::COUNTER_NT_STAT_TO_KPKT);
}


static bool realloc_frac_excitations_list(const int modelgridindex, const int newsize)
{
  struct nt_excitation_struct *newptr = (struct nt_excitation_struct *) realloc(
    nt_solution[modelgridindex].frac_excitations_list,
    newsize * sizeof(struct nt_excitation_struct));

  if (newptr == NULL && newsize > 0)
  {
    printout("ERROR: Not enough memory to reallocate NT excitation list for cell %d from size %d to %d.\n",
             modelgridindex, nt_solution[modelgridindex].frac_excitations_list_size, newsize);
    // abort();
    return false;
  }
  else
  {
    nt_solution[modelgridindex].frac_excitations_list = newptr;
    nt_solution[modelgridindex].frac_excitations_list_size = newsize;
    return true;
  }
}


static void analyse_sf_solution(const int modelgridindex, const int timestep)
{
  const float nne = get_nne(modelgridindex);
  const double nntot = get_tot_nion(modelgridindex);
  const double nnetot = get_nnetot(modelgridindex);

  // store the solution properties now while the NT spectrum is in memory (in case we free before packet prop)
  nt_solution[modelgridindex].frac_heating = calculate_frac_heating(modelgridindex);

  double frac_excitation_total = 0.;
  double frac_ionization_total = 0.;

  int excitationindex = 0; // unique index for every included excitation transition
  for (int element = 0; element < get_nelements(); element++)
  {
    const int Z = get_element(element);
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      const int uniqueionindex = get_uniqueionindex(element, ion);

      const int ionstage = get_ionstage(element, ion);
      const double nnion = ionstagepop(modelgridindex, element, ion);

      // if (nnion < minionfraction * get_tot_nion(modelgridindex)) // skip negligible ions
      if (nnion <= 0.) // skip zero-abundance ions
        continue;

      double frac_ionization_ion = 0.;
      double frac_excitation_ion = 0.;
      printout("  Z=%d ion_stage %d:\n", Z, ionstage);
      // printout("    nnion: %g\n", nnion);
      printout("    nnion/nntot: %g\n", nnion / nntot);

      calculate_eff_ionpot_auger_rates(modelgridindex, element, ion);

      int matching_nlsubshell_count = 0;
      for (int n = 0; n < colliondatacount; n++)
      {
        if (colliondata[n].Z == Z && colliondata[n].nelec == Z - ionstage + 1)
        {
          const double frac_ionization_ion_shell = calculate_nt_frac_ionization_shell(modelgridindex, element, ion, n);
          frac_ionization_ion += frac_ionization_ion_shell;
          matching_nlsubshell_count++;
          printout("      shell n %d, l %d, I %5.1f eV: frac_ionization %10.4e",
                   colliondata[n].n, colliondata[n].l, colliondata[n].ionpot_ev, frac_ionization_ion_shell);

          if (NT_MAX_AUGER_ELECTRONS > 0)
          {
            printout("  prob(n Auger elec):");
            for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
            {
              printout(" %d: %.2f", a, colliondata[n].prob_num_auger[a]);
            }
          }
          printout("\n");
        }
      }

      // do not ionize the top ion
      if (ion < nions - 1)
      {
        nt_solution[modelgridindex].fracdep_ionization_ion[uniqueionindex] = frac_ionization_ion;

        frac_ionization_total += frac_ionization_ion;
      }
      else
      {
        nt_solution[modelgridindex].fracdep_ionization_ion[uniqueionindex] = 0.;
      }
      printout("    frac_ionization: %g (%d subshells)\n", frac_ionization_ion, matching_nlsubshell_count);

      // excitation from all levels is very SLOW
      const int nlevels_all = get_nlevels(element, ion);
      // So limit the lower levels to improve performance
      const int nlevels = (nlevels_all > NTEXCITATION_MAXNLEVELS_LOWER) ? NTEXCITATION_MAXNLEVELS_LOWER : nlevels_all;
#if NT_EXCITATION_ON
      const bool above_minionfraction = (nnion >= minionfraction * get_tot_nion(modelgridindex));
#endif

      for (int lower = 0; lower < nlevels; lower++)
      {
        const double statweight_lower = stat_weight(element, ion, lower);
        const int nuptrans = get_nuptrans(element, ion, lower);
        const double nnlevel = calculate_exclevelpop(modelgridindex, element, ion, lower);
        const double epsilon_lower = epsilon(element, ion, lower);

        for (int t = 0; t < nuptrans; t++)
        {
          const int lineindex = globals::elements[element].ions[ion].levels[lower].uptrans_lineindicies[t];
          const int upper = globals::linelist[lineindex].upperlevelindex;
          if (upper >= NTEXCITATION_MAXNLEVELS_UPPER)
          {
            continue;
          }

          const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
          const double nt_frac_excitation_perlevelpop = epsilon_trans * calculate_nt_excitation_ratecoeff_perdeposition(modelgridindex, lineindex, statweight_lower, epsilon_trans);
          const double frac_excitation_thistrans = nnlevel * nt_frac_excitation_perlevelpop;
          frac_excitation_ion += frac_excitation_thistrans;

#if NT_EXCITATION_ON
          // the atomic data set was limited for Fe V, which caused the ground multiplet to be massively
          // depleted, and then almost no recombination happened!
          if (above_minionfraction && nt_frac_excitation_perlevelpop > 0 && !(Z == 26 && ionstage == 5))
          {
            if (excitationindex >= nt_solution[modelgridindex].frac_excitations_list_size)
            {
              const int newsize = nt_solution[modelgridindex].frac_excitations_list_size + NT_BLOCKSIZEEXCITATION;

              realloc_frac_excitations_list(modelgridindex, newsize);
            }

            if (excitationindex < nt_solution[modelgridindex].frac_excitations_list_size)
            {
              double ratecoeffperdeposition = nt_frac_excitation_perlevelpop / epsilon_trans;

              // if (get_coll_str(lineindex) < 0) // if collision strength is not defined, the rate coefficient is unreliable
              //   ratecoeffperdeposition = 0.;

              nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition = frac_excitation_thistrans;
              nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition = ratecoeffperdeposition;
              nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex = lineindex;
              (excitationindex)++;
            }
          }
#endif // NT_EXCITATION_ON
        } // for t
      } // for lower

      frac_excitation_total += frac_excitation_ion;

      printout("    frac_excitation: %g\n", frac_excitation_ion);
      printout("    workfn:       %9.2f eV\n", (1. / get_oneoverw(element, ion, modelgridindex)) / EV);
      printout("    eff_ionpot:   %9.2f eV  (always use valence potential is %s)\n",
               get_eff_ionpot(modelgridindex, element, ion) / EV, (NT_USE_VALENCE_IONPOTENTIAL ? "true" : "false"));

      printout("    workfn approx Gamma:     %9.3e\n",
               nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion));

      printout("    SF integral Gamma:       %9.3e\n",
              calculate_nt_ionization_ratecoeff(modelgridindex, element, ion, false));

      printout("    SF integral(I=Iv) Gamma: %9.3e  (if always use valence potential)\n",
              calculate_nt_ionization_ratecoeff(modelgridindex, element, ion, true));

      printout("    ARTIS using Gamma:       %9.3e\n",
               nt_ionization_ratecoeff(modelgridindex, element, ion));

      // the ion values (unlike shell ones) have been collapsed down to ensure that upperion < nions
      if (ion < nions - 1)
      {
        printout("    probability to ionstage:");
        double prob_sum = 0.;
        for (int upperion = ion + 1; upperion <= nt_ionisation_maxupperion(element, ion); upperion++)
        {
          const double probability = nt_ionization_upperion_probability(modelgridindex, element, ion, upperion, false);
          prob_sum += probability;
          if (probability > 0.)
            printout(" %d: %.3f", get_ionstage(element, upperion), probability);
        }
        printout("\n");
        assert((fabs(prob_sum - 1.0) <= 1e-2) || (nt_ionization_ratecoeff_sf(modelgridindex, element, ion) < 1e-20));

        printout("         enfrac to ionstage:");
        double enfrac_sum = 0.;
        for (int upperion = ion + 1; upperion <= nt_ionisation_maxupperion(element, ion); upperion++)
        {
          const double probability = nt_ionization_upperion_probability(modelgridindex, element, ion, upperion, true);
          enfrac_sum += probability;
          if (probability > 0.)
            printout(" %d: %.3f", get_ionstage(element, upperion), probability);
        }
        printout("\n");
        assert(fabs(enfrac_sum - 1.0) <= 1e-2 || (nt_ionization_ratecoeff_sf(modelgridindex, element, ion) < 1e-20));
      }
    }
  }

  if (NT_EXCITATION_ON && (MAX_NT_EXCITATIONS_STORED > 0))
  {
    if (excitationindex < nt_solution[modelgridindex].frac_excitations_list_size)
    {
      // shrink the list to match the data
      realloc_frac_excitations_list(modelgridindex, excitationindex);
    }

    qsort(nt_solution[modelgridindex].frac_excitations_list,
          nt_solution[modelgridindex].frac_excitations_list_size, sizeof(struct nt_excitation_struct),
          compare_excitation_fractions);

    // the excitation list is now sorted by frac_deposition descending
    const double deposition_rate_density = get_deposition_rate_density(modelgridindex);

    if (nt_solution[modelgridindex].frac_excitations_list_size > MAX_NT_EXCITATIONS_STORED)
    {
      // truncate the sorted list to save memory
      printout("  Truncating non-thermal excitation list from %d to %d transitions.\n",
               nt_solution[modelgridindex].frac_excitations_list_size, MAX_NT_EXCITATIONS_STORED);
      realloc_frac_excitations_list(modelgridindex, MAX_NT_EXCITATIONS_STORED);
    }

    printout("mem_usage: non-thermal excitations for cell %d at this timestep occupy %.1f MB\n",
             modelgridindex, nt_solution[modelgridindex].frac_excitations_list_size *
             (sizeof(nt_solution[modelgridindex].frac_excitations_list) + sizeof(nt_solution[modelgridindex].frac_excitations_list[0])) / 1024. / 1024.);

    const float T_e = get_Te(modelgridindex);
    printout("  Top non-thermal excitation fractions (total excitations = %d):\n",
             nt_solution[modelgridindex].frac_excitations_list_size);
    int ntransdisplayed = nt_solution[modelgridindex].frac_excitations_list_size;
    ntransdisplayed = (ntransdisplayed > 50) ? 50 : ntransdisplayed;
    for (excitationindex = 0; excitationindex < ntransdisplayed; excitationindex++)
    {
      const double frac_deposition = nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition;
      if (frac_deposition > 0.)
      {
        const int lineindex = nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex;
        const int element = globals::linelist[lineindex].elementindex;
        const int ion = globals::linelist[lineindex].ionindex;
        const int lower = globals::linelist[lineindex].lowerlevelindex;
        const int upper = globals::linelist[lineindex].upperlevelindex;
        const double epsilon_trans = epsilon(element, ion, upper) - epsilon(element, ion, lower);

        const double ratecoeffperdeposition = nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition;
        const double ntcollexc_ratecoeff = ratecoeffperdeposition * deposition_rate_density;

        const double t_current = globals::time_step[timestep].start;
        const double radexc_ratecoeff = rad_excitation_ratecoeff(modelgridindex, element, ion, lower, upper, epsilon_trans, lineindex, t_current);

        const double collexc_ratecoeff = col_excitation_ratecoeff(T_e, nne, lineindex, epsilon_trans);

        const double exc_ratecoeff = radexc_ratecoeff + collexc_ratecoeff + ntcollexc_ratecoeff;

        printout("    frac_deposition %.3e Z=%d ionstage %d lower %4d upper %4d rad_exc %.1e coll_exc %.1e nt_exc %.1e nt/tot %.1e collstr %.1e lineindex %d\n",
                 frac_deposition, get_element(element), get_ionstage(element, ion), lower, upper,
                 radexc_ratecoeff, collexc_ratecoeff, ntcollexc_ratecoeff, ntcollexc_ratecoeff / exc_ratecoeff, get_coll_str(lineindex), lineindex);
      }
    }

    // now sort the excitation list by lineindex ascending for fast lookup with a binary search
    qsort(nt_solution[modelgridindex].frac_excitations_list,
          nt_solution[modelgridindex].frac_excitations_list_size, sizeof(struct nt_excitation_struct),
          compare_excitation_lineindicies);

  } // NT_EXCITATION_ON

  const float frac_heating = get_nt_frac_heating(modelgridindex);

  // calculate number density of non-thermal electrons
  const double deposition_rate_density_ev = get_deposition_rate_density(modelgridindex) / EV;
  const double yscalefactor = deposition_rate_density_ev / E_init_ev;

  double nne_nt_max = 0.0;
  for (int i = 0; i < SFPTS; i++)
  {
    const double endash = gsl_vector_get(envec, i);
    #if (SF_USE_LOG_E_INCREMENT)
    const double delta_endash = gsl_vector_get(delta_envec, i);
    #else
    const double delta_endash = DELTA_E;
    #endif
    const double oneovervelocity = sqrt(9.10938e-31 / 2 / endash / 1.60218e-19) / 100; // in sec/cm
    nne_nt_max += yscalefactor * get_y_sample(modelgridindex, i) * oneovervelocity * delta_endash;
  }

  nt_solution[modelgridindex].frac_excitation = frac_excitation_total;
  nt_solution[modelgridindex].frac_ionization = frac_ionization_total;

  printout("  E_init:      %9.2f eV/s/cm^3\n", E_init_ev);
  printout("  deposition:  %9.2f eV/s/cm^3\n", deposition_rate_density_ev);
  printout("  nne:         %9.3e e-/cm^3\n", nne);
  printout("  nnetot:      %9.3e e-/cm^3\n", nnetot);
  printout("  nne_nt     < %9.3e e-/cm^3\n", nne_nt_max);
  printout("  nne_nt/nne < %9.3e\n", nne_nt_max / nne);
  printout("  frac_heating_tot:    %g\n", frac_heating);
  printout("  frac_excitation_tot: %g\n", frac_excitation_total);
  printout("  frac_ionization_tot: %g\n", frac_ionization_total);
  const double frac_sum = frac_heating + frac_excitation_total + frac_ionization_total;
  printout("  frac_sum:            %g (should be close to 1.0)\n", frac_sum);

  nt_solution[modelgridindex].frac_heating = 1. - frac_excitation_total - frac_ionization_total;
  printout("  (replacing calculated frac_heating_tot with %g to make frac_sum = 1.0)\n",
           nt_solution[modelgridindex].frac_heating);

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


static void sfmatrix_add_excitation(gsl_matrix *const sfmatrix, const int modelgridindex, const int element, const int ion)
{
  // excitation terms
  gsl_vector *vec_xs_excitation_deltae = gsl_vector_alloc(SFPTS);

  const int nlevels_all = get_nlevels(element, ion);
  const int nlevels = (nlevels_all > NTEXCITATION_MAXNLEVELS_LOWER) ? NTEXCITATION_MAXNLEVELS_LOWER : nlevels_all;

  for (int lower = 0; lower < nlevels; lower++)
  {
    const double statweight_lower = stat_weight(element, ion, lower);
    const double nnlevel = calculate_exclevelpop(modelgridindex, element, ion, lower);
    const double epsilon_lower = epsilon(element, ion, lower);
    const int nuptrans = get_nuptrans(element, ion, lower);
    for (int t = 0; t < nuptrans; t++)
    {
      const int lineindex = globals::elements[element].ions[ion].levels[lower].uptrans_lineindicies[t];
      const int upper = globals::linelist[lineindex].upperlevelindex;
      if (upper >= NTEXCITATION_MAXNLEVELS_UPPER)
      {
        continue;
      }
      const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
      const double epsilon_trans_ev = epsilon_trans / EV;
      if (epsilon_trans_ev < SF_EMIN)
      {
        continue;
      }

      const int xsstartindex = get_xs_excitation_vector(vec_xs_excitation_deltae, lineindex, statweight_lower, epsilon_trans);
      if (xsstartindex >= 0)
      {
        #if (SF_USE_LOG_E_INCREMENT)
        gsl_vector_mul(vec_xs_excitation_deltae, delta_envec);
        #else
        gsl_blas_dscal(DELTA_E, vec_xs_excitation_deltae);
        #endif

        for (int i = 0; i < SFPTS; i++)
        {
          const double en = gsl_vector_get(envec, i);
          const int stopindex = get_energyindex_ev_lteq(en + epsilon_trans_ev);

          const int startindex = i > xsstartindex ? i : xsstartindex;
          for (int j = startindex; j < stopindex; j++)
          {
            *gsl_matrix_ptr(sfmatrix, i, j) += nnlevel * gsl_vector_get(vec_xs_excitation_deltae, j);
          }

          // do the last bit separately because we're not using the full delta_e interval
          #if (SF_USE_LOG_E_INCREMENT)
          const double delta_en = gsl_vector_get(delta_envec, stopindex);
          #else
          const double delta_en = DELTA_E;
          #endif

          const double delta_en_actual = (en + epsilon_trans_ev - gsl_vector_get(envec, stopindex));

          *gsl_matrix_ptr(sfmatrix, i, stopindex) += nnlevel * gsl_vector_get(vec_xs_excitation_deltae, stopindex) * delta_en_actual / delta_en;
        }
      }
    }
  }
  gsl_vector_free(vec_xs_excitation_deltae);
}


static void sfmatrix_add_ionization(gsl_matrix *const sfmatrix, const int Z, const int ionstage, const double nnion)
// add the ionization terms to the Spencer-Fano matrix
{
  gsl_vector *const vec_xs_ionization = gsl_vector_alloc(SFPTS);
  for (int collionindex = 0; collionindex < colliondatacount; collionindex++)
  {
    if (colliondata[collionindex].Z == Z && colliondata[collionindex].nelec == Z - ionstage + 1)
    {
      const double ionpot_ev = colliondata[collionindex].ionpot_ev;
      const double en_auger_ev = colliondata[collionindex].en_auger_ev;
      // const double n_auger_elec_avg = colliondata[n].n_auger_elec_avg;
      const double J = get_J(Z, ionstage, ionpot_ev);

      assert(ionpot_ev >= SF_EMIN);

      // printout("Z=%2d ion_stage %d n %d l %d ionpot %g eV\n",
      //          Z, ionstage, colliondata[n].n, colliondata[n].l, ionpot_ev);

      const int xsstartindex = get_xs_ionization_vector(vec_xs_ionization, collionindex);

      double atanexp[SFPTS];
      double prefactors[SFPTS];
      for (int j = xsstartindex; j < SFPTS; j++)
      {
        const double endash = gsl_vector_get(envec, j);
        const double epsilon_upper = (endash + ionpot_ev) / 2;
        atanexp[j] = atan((epsilon_upper - ionpot_ev) / J);
        prefactors[j] = gsl_vector_get(vec_xs_ionization, j) * nnion / atan((endash - ionpot_ev) / 2 / J);
      }

      for (int i = 0; i < SFPTS; i++)
      {
        // i is the matrix row index, which corresponds to an energy E at which we are solve from y(E)
        const double en = gsl_vector_get(envec, i);

        // endash ranges from en to SF_EMAX, but skip over the zero-cross section points
        const int jstart = i > xsstartindex ? i : xsstartindex;
        for (int j = jstart; j < SFPTS; j++)
        {
          // j is the matrix column index which corresponds to the piece of the integral at y(E') where E' >= E and E' = envec(j)
          const double endash = gsl_vector_get(envec, j);
          #if (SF_USE_LOG_E_INCREMENT)
          const double deltaendash = gsl_vector_get(delta_envec, j);
          #else
          const double deltaendash = DELTA_E;
          #endif

          // J * atan[(epsilon - ionpot_ev) / J] is the indefinite integral of 1/[1 + (epsilon - ionpot_ev)^2/ J^2]
          // in Kozma & Fransson 1992 equation 4

          const double epsilon_lower = endash - en; // and epsilon_upper = (endash + ionpot_ev) / 2;
          *gsl_matrix_ptr(sfmatrix, i, j) += prefactors[j] * (atanexp[j] - atan((epsilon_lower - ionpot_ev) / J)) * deltaendash;
        }

        // below is atan((epsilon_lower - ionpot_ev) / J) where epsilon_lower = en + ionpot_ev;
        const double atanexp2 = atan(en / J);

        // endash ranges from 2 * en + ionpot_ev to SF_EMAX
        const int secondintegralstartindex = get_energyindex_ev_lteq(2 * en + ionpot_ev);
        for (int j = secondintegralstartindex; j < SFPTS; j++)
        {
          #if (SF_USE_LOG_E_INCREMENT)
          const double deltaendash = gsl_vector_get(delta_envec, j);
          #else
          const double deltaendash = DELTA_E;
          #endif

          // epsilon_lower = en + ionpot_ev;
          // epsilon_upper = (endash + ionpot_ev) / 2;
          *gsl_matrix_ptr(sfmatrix, i, j) -= prefactors[j] * (atanexp[j] - atanexp2) * deltaendash;
        }
      }

      if (SF_AUGER_CONTRIBUTION_ON)
      {
        int augerstopindex = 0;
        if (SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN)
        {
          // en_auger_ev is (if LJS understands it correctly) averaged to include some probability of zero Auger electrons
          // so we need a boost to get the average energy of Auger electrons given that there are one or more
          const double en_boost = 1 / (1. - colliondata[collionindex].prob_num_auger[0]);

          augerstopindex = get_energyindex_ev_gteq(en_auger_ev * en_boost);
        }
        else
        {
          augerstopindex = get_energyindex_ev_gteq(en_auger_ev);
        }

        for (int i = 0; i < augerstopindex; i++)
        {
          const double en = gsl_vector_get(envec, i);
          const int jstart = i > xsstartindex ? i : xsstartindex;
          for (int j = jstart; j < SFPTS; j++)
          {
            const double xs = gsl_vector_get(vec_xs_ionization, j);
            if (SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN)
            {
              const double en_boost = 1 / (1. - colliondata[collionindex].prob_num_auger[0]);
              for (int a = 1; a <= NT_MAX_AUGER_ELECTRONS; a++)
              {
                if (en < (en_auger_ev * en_boost / a))
                {
                  *gsl_matrix_ptr(sfmatrix, i, j) -= nnion * xs * colliondata[collionindex].prob_num_auger[a] * a;
                }
              }
            }
            else
            {
              assert(en < en_auger_ev);
              // printout("SFAuger E %g < en_auger_ev %g so subtracting %g from element with value %g\n", en, en_auger_ev, nnion * xs, ij_contribution);
              *gsl_matrix_ptr(sfmatrix, i, j) -= nnion * xs; // * n_auger_elec_avg; // * en_auger_ev???
            }
          }
        }
      }
    }
  }
  free(vec_xs_ionization);
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


void solve_spencerfano(const int modelgridindex, const int timestep, const int iteration)
// solve the Spencer-Fano equation to get the non-thermal electron flux energy distribution
// based on Equation (2) of Li et al. (2012)
{
  if (get_numassociatedcells(modelgridindex) < 1)
  {
    printout("Associated_cells < 1 in cell %d at timestep %d. Skipping Spencer-Fano solution.\n", modelgridindex, timestep);

    return;
  }
  else if (get_deposition_rate_density(modelgridindex) / EV < MINDEPRATE)
  {
    printout("Non-thermal deposition rate of %g eV/cm/s/cm^3 below  MINDEPRATE %g in cell %d at timestep %d. Skipping Spencer-Fano solution.\n",
    get_deposition_rate_density(modelgridindex) / EV, MINDEPRATE, modelgridindex, timestep);

    nt_solution[modelgridindex].frac_heating = 0.97;
    nt_solution[modelgridindex].frac_ionization = 0.03;
    nt_solution[modelgridindex].frac_excitation = 0.;

    nt_solution[modelgridindex].nneperion_when_solved = -1.;
    nt_solution[modelgridindex].timestep_last_solved = -1;

    free(nt_solution[modelgridindex].frac_excitations_list);
    nt_solution[modelgridindex].frac_excitations_list_size = 0;

    zero_all_effionpot(modelgridindex);
    return;
  }

  const float nne = get_nne(modelgridindex); // electrons per cm^3
  // const double nnetot = get_nnetot(modelgridindex);
  const double nne_per_ion = nne / get_tot_nion(modelgridindex);
  const double nne_per_ion_last = nt_solution[modelgridindex].nneperion_when_solved;
  const double nne_per_ion_fracdiff = fabs((nne_per_ion_last / nne_per_ion) - 1.);
  const int timestep_last_solved = nt_solution[modelgridindex].timestep_last_solved;

  printout("Spencer-Fano solver at timestep %d (last solution was at timestep %d) nne/niontot = %g, at last solution was %g fracdiff %g\n",
           timestep, timestep_last_solved, nne_per_ion, nne_per_ion_last, nne_per_ion_fracdiff);

  if ((nne_per_ion_fracdiff < NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS) && (timestep - timestep_last_solved <= SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS) && timestep_last_solved > globals::n_lte_timesteps)
  {
    printout("Keeping Spencer-Fano solution from timestep %d because x_e fracdiff %g < %g and because timestep %d - %d < %d\n",
             timestep_last_solved, nne_per_ion_fracdiff, NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS, timestep, timestep_last_solved, SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS);

    return;
  }

  printout("Setting up Spencer-Fano equation with %d energy points from %g eV to %g eV in cell %d at timestep %d iteration %d (nne=%g e-/cm^3)\n",
           SFPTS, SF_EMIN, SF_EMAX, modelgridindex, timestep, iteration, nne);

  nt_solution[modelgridindex].nneperion_when_solved = nne_per_ion;
  nt_solution[modelgridindex].timestep_last_solved = timestep;

  bool enable_sfexcitation = true;
  bool enable_sfionization = true;
  if (timestep <= globals::n_lte_timesteps)
  {
    // for the first run of the solver at the first NLTE timestep (which usually requires many iterations),
    // do a fast initial solution but mark it has an invalid nne per ion so it gets replaced at the next timestep
    nt_solution[modelgridindex].nneperion_when_solved = -1.;
    enable_sfexcitation = false;
    enable_sfionization = false;

    printout("Doing a fast initial solution without ionization or excitation in the SF equation for the first NLTE timestep.\n");
  }
  else if (timestep <= globals::n_lte_timesteps + 2)
  {
    // run the solver in a faster mode for the first couple of NLTE timesteps
    // nt_solution[modelgridindex].nneperion_when_solved = -1.;
    enable_sfexcitation = false;
    // enable_sfionization = false;

    printout("Doing a faster solution without excitation in the SF equation for the first couple of NLTE timesteps.\n");
  }

  gsl_matrix *const sfmatrix = gsl_matrix_calloc(SFPTS, SFPTS);
  gsl_vector *const rhsvec = gsl_vector_calloc(SFPTS); // constant term (not dependent on y func) in each equation

  // loss terms and source terms
  for (int i = 0; i < SFPTS; i++)
  {
    const double en = gsl_vector_get(envec, i);

    *gsl_matrix_ptr(sfmatrix, i, i) += electron_loss_rate(en * EV, nne) / EV;

    double source_integral_to_SF_EMAX;
    if (i < SFPTS - 1)
    {
      #if (SF_USE_LOG_E_INCREMENT)
      gsl_vector *sourcevec_de = gsl_vector_alloc(SFPTS);
      gsl_vector_memcpy(sourcevec_de, sourcevec);
      gsl_vector_mul(sourcevec_de, delta_envec);
      gsl_vector_const_view source_de_e_to_SF_EMAX = gsl_vector_const_subvector(sourcevec_de, i + 1, SFPTS - i - 1);
      source_integral_to_SF_EMAX = gsl_blas_dasum(&source_de_e_to_SF_EMAX.vector);
      gsl_vector_free(sourcevec_de);
      #else
      gsl_vector_const_view source_e_to_SF_EMAX = gsl_vector_const_subvector(sourcevec, i + 1, SFPTS - i - 1);
      source_integral_to_SF_EMAX = gsl_blas_dasum(&source_e_to_SF_EMAX.vector) * DELTA_E;
      #endif
    }
    else
      source_integral_to_SF_EMAX = 0;

    gsl_vector_set(rhsvec, i, source_integral_to_SF_EMAX);
  }
  // gsl_vector_set_all(rhsvec, 1.); // alternative if all electrons are injected at SF_EMAX

  if (enable_sfexcitation || enable_sfionization)
  {
    for (int element = 0; element < get_nelements(); element++)
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

        if (enable_sfexcitation)
          sfmatrix_add_excitation(sfmatrix, modelgridindex, element, ion);

        if (enable_sfionization && (ion < nions - 1))
          sfmatrix_add_ionization(sfmatrix, Z, ionstage, nnion);
      }
      if (!first_included_ion_of_element)
        printout("\n");
    }
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
    nt_solution[modelgridindex].yfunc = (double *) calloc(SFPTS, sizeof(double));
  }

  gsl_vector_view yvecview = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);
  gsl_vector *yvec = &yvecview.vector;

  sfmatrix_solve(sfmatrix, rhsvec, yvec);

  // gsl_matrix_free(sfmatrix_LU); // if sfmatrix_LU is different to sfmatrix

  gsl_matrix_free(sfmatrix);
  gsl_vector_free(rhsvec);

  if (timestep % 10 == 0)
    nt_write_to_file(modelgridindex, timestep, iteration);

  analyse_sf_solution(modelgridindex, timestep);

  if (!STORE_NT_SPECTRUM)
  {
    free(nt_solution[modelgridindex].yfunc);
  }
}


void write_restart_data(FILE *gridsave_file)
{
  printout("non-thermal solver, ");

  fprintf(gridsave_file, "%d\n", 24724518); // special number marking the beginning of NT data
  fprintf(gridsave_file, "%d %lg %lg\n", SFPTS, SF_EMIN, SF_EMAX);

  for (int modelgridindex = 0; modelgridindex < globals::npts_model; modelgridindex++)
  {
    if (get_numassociatedcells(modelgridindex) > 0)
    {
      fprintf(gridsave_file, "%d %d %lg ",
              modelgridindex,
              deposition_rate_density_timestep[modelgridindex],
              deposition_rate_density[modelgridindex]);

      if (NT_ON && NT_SOLVE_SPENCERFANO)
      {
        check_auger_probabilities(modelgridindex);

        fprintf(gridsave_file, "%g %g %g %g\n",
                nt_solution[modelgridindex].nneperion_when_solved,
                nt_solution[modelgridindex].frac_heating,
                nt_solution[modelgridindex].frac_ionization,
                nt_solution[modelgridindex].frac_excitation);

        for (int uniqueionindex = 0; uniqueionindex < globals::includedions; uniqueionindex++)
        {
          fprintf(gridsave_file, "%lg ", nt_solution[modelgridindex].fracdep_ionization_ion[uniqueionindex]);
          fprintf(gridsave_file, "%g ", nt_solution[modelgridindex].eff_ionpot[uniqueionindex]);

          for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
          {
            fprintf(gridsave_file, "%g ", nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a]);
            fprintf(gridsave_file, "%g ", nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a]);
          }
        }

        // write NT excitations
        fprintf(gridsave_file, "%d\n", nt_solution[modelgridindex].frac_excitations_list_size);

        const int frac_excitations_list_size = nt_solution[modelgridindex].frac_excitations_list_size;
        for (int excitationindex = 0; excitationindex < frac_excitations_list_size; excitationindex++)
        {
          fprintf(gridsave_file, "%lg %lg %d\n",
                  nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition,
                  nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition,
                  nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex);
        }

        // write non-thermal spectrum
        if (STORE_NT_SPECTRUM)
        {
          for (int s = 0; s < SFPTS; s++)
          {
            fprintf(gridsave_file, "%lg\n", nt_solution[modelgridindex].yfunc[s]);
          }
        }
      }
    }
  }
}


void read_restart_data(FILE *gridsave_file)
{
  printout("Reading restart data for non-thermal solver\n");

  int code_check;
  fscanf(gridsave_file, "%d\n", &code_check);
  if (code_check != 24724518)
  {
    printout("ERROR: Beginning of non-thermal restart data not found! Found %d instead of 24724518\n", code_check);
    abort();
  }

  int sfpts_in;
  double SF_EMIN_in;
  double SF_EMAX_in;
  fscanf(gridsave_file, "%d %lg %lg\n", &sfpts_in, &SF_EMIN_in, &SF_EMAX_in);

  if (sfpts_in != SFPTS || SF_EMIN_in != SF_EMIN || SF_EMAX_in != SF_EMAX)
  {
    printout("ERROR: gridsave file specifies %d Spencer-Fano samples, SF_EMIN %lg SF_EMAX %lg\n",
             sfpts_in, SF_EMIN_in, SF_EMAX_in);
    printout("ERROR: This simulation has %d Spencer-Fano samples, SF_EMIN %lg SF_EMAX %lg\n",
             SFPTS, SF_EMIN, SF_EMAX);
    abort();
  }

  for (int modelgridindex = 0; modelgridindex < globals::npts_model; modelgridindex++)
  {
    if (get_numassociatedcells(modelgridindex) > 0)
    {
      int mgi_in;
      fscanf(gridsave_file, "%d %d %lg ",
             &mgi_in,
             &deposition_rate_density_timestep[modelgridindex],
             &deposition_rate_density[modelgridindex]);

      if (NT_ON && NT_SOLVE_SPENCERFANO)
      {
        fscanf(gridsave_file, "%g %g %g %g\n",
               &nt_solution[modelgridindex].nneperion_when_solved,
               &nt_solution[modelgridindex].frac_heating,
               &nt_solution[modelgridindex].frac_ionization,
               &nt_solution[modelgridindex].frac_excitation);

        if (mgi_in != modelgridindex)
        {
          printout("ERROR: expected data for cell %d but found cell %d\n", modelgridindex, mgi_in);
          abort();
        }

        for (int uniqueionindex = 0; uniqueionindex < globals::includedions; uniqueionindex++)
        {
          fscanf(gridsave_file, "%lg ", &nt_solution[modelgridindex].fracdep_ionization_ion[uniqueionindex]),
          fscanf(gridsave_file, "%g ", &nt_solution[modelgridindex].eff_ionpot[uniqueionindex]);

          for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++)
          {
            fscanf(gridsave_file, "%g ", &nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a]);
            fscanf(gridsave_file, "%g ", &nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a]);
          }
        }

        check_auger_probabilities(modelgridindex);

        // read NT excitations
        int frac_excitations_list_size_in;
        fscanf(gridsave_file, "%d\n", &frac_excitations_list_size_in);

        if (nt_solution[modelgridindex].frac_excitations_list_size != frac_excitations_list_size_in)
        {
          realloc_frac_excitations_list(modelgridindex, frac_excitations_list_size_in);
        }

        for (int excitationindex = 0; excitationindex < frac_excitations_list_size_in; excitationindex++)
        {
          fscanf(gridsave_file, "%lg %lg %d\n",
                  &nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition,
                  &nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition,
                  &nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex);
        }

        // read non-thermal spectrum
        if (STORE_NT_SPECTRUM)
        {
          for (int s = 0; s < SFPTS; s++)
          {
            fscanf(gridsave_file, "%lg\n", &nt_solution[modelgridindex].yfunc[s]);
          }
        }
      }
    }
  }
}


#ifdef MPI_ON
void nt_MPI_Bcast(const int modelgridindex, const int root)
{
  if (get_numassociatedcells(modelgridindex) == 0)
    return;

  // printout("nonthermal_MPI_Bcast cell %d before: ratecoeff(Z=%d ion_stage %d): %g, eff_ionpot %g eV\n",
  //          modelgridindex, logged_element_z, logged_ion_stage,
  //          nt_ionization_ratecoeff_sf(modelgridindex, logged_element_index, logged_ion_index),
  //          get_eff_ionpot(modelgridindex, logged_element_index, logged_ion_index) / EV);

  MPI_Bcast(&deposition_rate_density_timestep[modelgridindex], 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&deposition_rate_density[modelgridindex], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

  if (NT_ON && NT_SOLVE_SPENCERFANO)
  {
    assert(nonthermal_initialized);
    MPI_Bcast(&nt_solution[modelgridindex].nneperion_when_solved, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nt_solution[modelgridindex].timestep_last_solved, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nt_solution[modelgridindex].frac_heating, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nt_solution[modelgridindex].frac_ionization, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nt_solution[modelgridindex].frac_excitation, 1, MPI_FLOAT, root, MPI_COMM_WORLD);

    MPI_Bcast(nt_solution[modelgridindex].fracdep_ionization_ion, globals::includedions, MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(nt_solution[modelgridindex].eff_ionpot, globals::includedions, MPI_FLOAT, root, MPI_COMM_WORLD);

    MPI_Bcast(nt_solution[modelgridindex].prob_num_auger, globals::includedions * (NT_MAX_AUGER_ELECTRONS + 1), MPI_FLOAT, root, MPI_COMM_WORLD);
    MPI_Bcast(nt_solution[modelgridindex].ionenfrac_num_auger, globals::includedions * (NT_MAX_AUGER_ELECTRONS + 1), MPI_FLOAT, root, MPI_COMM_WORLD);

    // communicate NT excitations
    const int frac_excitations_list_size_old = nt_solution[modelgridindex].frac_excitations_list_size;
    MPI_Bcast(&nt_solution[modelgridindex].frac_excitations_list_size, 1, MPI_INT, root, MPI_COMM_WORLD);

    if (nt_solution[modelgridindex].frac_excitations_list_size != frac_excitations_list_size_old)
    {
      assert(realloc_frac_excitations_list(modelgridindex, nt_solution[modelgridindex].frac_excitations_list_size));
    }

    const int frac_excitations_list_size = nt_solution[modelgridindex].frac_excitations_list_size;
    for (int excitationindex = 0; excitationindex < frac_excitations_list_size; excitationindex++)
    {
      MPI_Bcast(&nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex, 1, MPI_INT, root, MPI_COMM_WORLD);
    }

    if (STORE_NT_SPECTRUM)
    {
      assert(nt_solution[modelgridindex].yfunc != NULL);
      MPI_Bcast(nt_solution[modelgridindex].yfunc, SFPTS, MPI_DOUBLE, root, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    check_auger_probabilities(modelgridindex);
  }
}
#endif


void nt_reset_stats(void)
{
  nt_energy_deposited = 0.;
}


void nt_print_stats(const int timestep, const double modelvolume, const double deltat)
{
  const double deposition_rate_density_montecarlo = nt_energy_deposited / EV / modelvolume / deltat;

  // deposition rate density for all cells has not been communicated yet - could change this
  // double total_deposition_rate_density = 0.;
  // for (int mgi = 0; mgi < npts_model; mgi++)
  // {
  //   total_deposition_rate_density += get_deposition_rate_density(mgi) / EV;
  // }
  printout("nt_energy_deposited = %9.2f eV/s/cm^3\n", deposition_rate_density_montecarlo);
}

}  // namespace nonthermal