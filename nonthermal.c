#include <math.h>
#include <stdbool.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "atomic.h"
#include "grid_init.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "sn3d.h"

#define SFPTS 1000  // number of points in Spencer-Fano solution function

#define EMAX 1000. // eV
#define EMIN 1. // eV
const double DELTA_E = (EMAX - EMIN) / SFPTS;

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

struct collionrow *colliondata = NULL;
int colliondatacount = 0;

static FILE *restrict nonthermalfile = NULL;
bool nonthermal_initialized = false;

gsl_vector *yvec = NULL; // Spencer-Fano solution function samples. multiply by energy to get flux
gsl_vector *envec = NULL; // energy grid on which solution is sampled
gsl_vector *sourcevec = NULL; // source distribution function samples (integral normalised to one)
double E_init = 0; // the average energy of the source electrons in eV
double E_0 = 0.; // the lowest ionization or excitation energy in eV

// double cellvolume[MMODELGRID+1];

static void read_collion_data()
{
  printout("Reading collisional ionization data...\n");

  FILE *cifile = fopen("ci.dat", "r");
  if (cifile == NULL)
  {
    printout("Could not open ci.dat\n");
    abort();
  }

  fscanf(cifile, "%d", &colliondatacount);
  printout("Reading %d lines\n", colliondatacount);
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


void nonthermal_init(void)
{
  printout("Initializing non-thermal solver\n");
  if (nonthermal_initialized == false)
  {
    const char filename[100] = "nonthermal.out";
    nonthermalfile = fopen(filename, "w");
    if (nonthermalfile == NULL)
    {
      printout("Cannot open %s.\n",filename);
      exit(0);
    }
    fprintf(nonthermalfile,"%8s %15s %8s %11s %11s %11s\n",
            "timestep","modelgridindex","index","energy_ev","source","y");
    fflush(nonthermalfile);

    yvec = gsl_vector_calloc(SFPTS); // Spencer-Fano solution function samples
                                     // y(E)*  dE is the flux of electrons with energy in the range (E, E + dE)
    envec = gsl_vector_calloc(SFPTS); // energy grid on which solution is sampled
    sourcevec = gsl_vector_calloc(SFPTS); // energy grid on which solution is sampled

    for (int s = 0; s < SFPTS; s++)
    {
      const double energy = EMIN + s * DELTA_E;

      gsl_vector_set(envec, s, energy);

      // spread the source over some energy width
      if (s < SFPTS - 10)
        gsl_vector_set(sourcevec, s, 0.);
      else
        gsl_vector_set(sourcevec, s, 1 / (DELTA_E * 10));
    }
    // put all of the source into one point at EMAX
    // gsl_vector_set(sourcevec, SFPTS-1, 1 / DELTA_E);

    E_init = 0.0;
    for (int i = 0; i < SFPTS; i++)
    {
      E_init += gsl_vector_get(envec, i) * gsl_vector_get(sourcevec, i) * DELTA_E;
    }

    double sourceintegral = 0.0;
    for (int i = 0; i < SFPTS; i++)
    {
      sourceintegral += gsl_vector_get(sourcevec, i) * DELTA_E;
    }
    printout("source vector integral: %14.7e\n", sourceintegral);

    read_collion_data();

    nonthermal_initialized = true;
  }
  printout("Finished initializing non-thermal solver\n");
}


static double get_deposition_rate_density(int modelgridindex)
// should be in erg / s / cm^3
{
  const double gamma_deposition = rpkt_emiss[modelgridindex] * 1.e20 * FOURPI;
  // Above is the gamma-ray bit. Below is *supposed* to be the kinetic energy of positrons created by 56Co and 48V. These formulae should be checked, however.

  const double t = time_step[nts_global].mid;
  const double rho = get_rho(modelgridindex);

  const double co56_positron_dep = (0.610 * 0.19 * MEV) *
        (exp(-t / TCOBALT) - exp(-t / TNICKEL)) /
        (TCOBALT - TNICKEL) * modelgrid[modelgridindex].fni * rho / MNI56;

  const double v48_positron_dep = (0.290 * 0.499 * MEV) *
        (exp(-t / T48V) - exp(-t / T48CR)) /
        (T48V - T48CR) * modelgrid[modelgridindex].f48cr * rho / MCR48;

  //printout("nt_deposition_rate: element: %d, ion %d\n",element,ion);
  //printout("nt_deposition_rate: gammadep: %g, poscobalt %g pos48v %g\n",
  //         gamma_deposition,co56_positron_dep,v48_positron_dep);

  return gamma_deposition + co56_positron_dep + v48_positron_dep;
}


static void nonthermal_write_to_file(int modelgridindex, int timestep)
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

  const double yfactor = (get_deposition_rate_density(modelgridindex) / (E_init * EV));
  for (int s = 0; s < SFPTS; s++)
  {
    fprintf(nonthermalfile,"%8d %15d %8d %11.5e %11.5e %11.5e\n",
            timestep,modelgridindex,s,gsl_vector_get(envec,s),gsl_vector_get(sourcevec,s),gsl_vector_get(yvec,s) * yfactor);
  }
  fflush(nonthermalfile);
# ifdef _OPENMP
  }
# endif
}


void nonthermal_close_file(void)
{
  fclose(nonthermalfile);
  gsl_vector_free(yvec);
  gsl_vector_free(envec);
  gsl_vector_free(sourcevec);
  free(colliondata);
  nonthermal_initialized = false;
}


static int get_energyindex_ev(double energy_ev)
{
  int index = 0;
  for (int i = 1; i < SFPTS; i++) // search through for the grid point matching the energy
  {
    if (gsl_vector_get(yvec, i) <= energy_ev)
      index = i;
  }
  return index;
}


static double xs_excitation(int lineindex, double epsilon_trans, double energy)
// excitation cross section in cm^2
// energies are in erg
{
  double sigma;
  const double coll_str_thisline = coll_str(lineindex);
  const double A_naught_sq = pow(5.2917721067e-9, 2); // Bohr radius squared in cm^2

  if (coll_str_thisline < 0)
  {
    const double fij = osc_strength(lineindex);
    if (coll_str_thisline > -1.5) // to catch -1
    {
      // permitted E1 electric dipole transitions
      const double prefactor = 45.585750051; // 8*pi^2/sqrt(3)
      const double U = energy / epsilon_trans;

      // const double g_bar = 0.2;
      const double A = 0.28;
      const double B = 0.15;
      const double g_bar = A * log(U) + B;

      // Eq 4 of Mewe 1972, possibly from Seaton 1962?
      sigma = prefactor * A_naught_sq * pow(H_ionpot / epsilon_trans, 2) * fij * g_bar / U;
    }
    else if (coll_str_thisline > -3.5) //to catch -2 or -3
    {
      // forbidden transitions: magnetic dipole, electric quadropole...
      sigma = 0.0;
    }
    else
    {
      sigma = 0.0;
    }
    // printout("lineindex %d excitation sigma %g, osc strength %g\n", lineindex, sigma, fij);
  }
  else
  {
    // collision strength is available, so use that
    sigma = pow(H_ionpot / energy, 2) / statw_lower(lineindex) * coll_str_thisline * PI * A_naught_sq;
  }

  return sigma;
}


static double electron_loss_rate(double energy, double nne)
// -dE / dx for fast electrons
// energy is in ergs
// return units are erg / cm
{
  const double omegap = sqrt(4 * PI * nne * pow(QE, 2) / ME);
  const double zetae = H  * omegap / 2 / PI;
  const double v = sqrt(2 * energy / ME);
  if (energy > 14 * EV)
  {
    return nne * 2 * PI * pow(QE,4) / energy * log(2 * energy / zetae);
  }
  else
  {
    const double eulergamma = 0.577215664901532;
    return nne * 2 * PI * pow(QE, 4) / energy * log(ME * pow(v, 3) / (eulergamma * pow(QE, 2) * omegap));
  }
}


static double xs_impactionization(double energy_ev, int collionindex)
// impact ionization cross section in cm^2
// energy and ionization_potential should be in eV
// fitting forumula of Younger 1981
// called Q_i(E) in KF92 equation 7
{
    const double A = colliondata[collionindex].A;
    const double B = colliondata[collionindex].B;
    const double C = colliondata[collionindex].C;
    const double D = colliondata[collionindex].D;
    const double ionpot_ev = colliondata[collionindex].ionpot_ev;
    const double u = energy_ev / ionpot_ev;
    // const double const1 = pow(H * CLIGHT, 2); //convert 1/E^2 to cm^2 where E in in eV
    const double const1 = 1e-14;

    return const1 * (A * (1 - 1/u) + B * pow((1 - 1/u), 2) + C * log(u) + D * log(u) / u) / (u * pow(ionpot_ev, 2));
}


static double Psecondary(double e_p, double epsilon, double I, double J)
// distribution of secondary electron energies for primary electron with energy e_p
// Opal, Peterson, & Beaty (1971)
{
    const double e_s = epsilon - I;
    return 1 / (J * atan((e_p - I) / 2 / J) * (1 + pow(e_s / J, 2)));
}


static double N_e(int modelgridindex, double energy)
// Kozma & Fransson equation 6.
// not sure what quantity means, but it's needed to calculate the heating fraction in equation 3
{
  const double energy_ev = energy / EV;
  double N_e = 0.;

  for (int element = 0; element < nelements; element++)
  {
    const int Z = get_element(element);
    const int nions = get_nions(element);

    for (int ion = 0; ion < nions; ion++)
    {
      double N_e_ion = 0.;
      const int ionstage = get_ionstage(element, ion);
      const double nnion = ionstagepop(modelgridindex,element,ion);

      // excitation terms
      // for (int level = 0; level < get_nlevels(element,ion); level++)
      const int level = 0; // just consider excitation from the ground level
      {
        const double epsilon_current = epsilon(element,ion,level);
        const int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
        for (int t = 1; t <= nuptrans; t++)
        {
          const double epsilon_target = elements[element].ions[ion].levels[level].uptrans[t].epsilon;
          const int lineindex = elements[element].ions[ion].levels[level].uptrans[t].lineindex;
          const double epsilon_trans = epsilon_target - epsilon_current;
          const double epsilon_trans_ev = epsilon_trans / EV;
          N_e_ion += gsl_vector_get(yvec, get_energyindex_ev(energy_ev + epsilon_trans_ev)) * xs_excitation(lineindex, epsilon_trans, energy + epsilon_trans);
        }
      }

      // ionization terms
      for (int n = 0; n < colliondatacount; n++)
      {
        if (colliondata[n].Z == Z && colliondata[n].nelec == Z - ionstage + 1)
        {
          const double ionpot_ev = colliondata[n].ionpot_ev;
          const double J = 0.6 * ionpot_ev;
          const double lambda = fmin(EMAX - energy_ev, energy_ev + ionpot_ev);

          for (int i = 0; i < SFPTS; i++)
          {
            const double endash = gsl_vector_get(envec, i);

            if (endash >= ionpot_ev && endash < lambda)
            {
              N_e_ion += gsl_vector_get(yvec, i) * xs_impactionization(energy_ev + endash, n) * Psecondary(energy_ev + endash, endash, ionpot_ev, J) * DELTA_E;
            }

            if (endash >= 2 * energy_ev + ionpot_ev)
            {
              N_e_ion += gsl_vector_get(yvec, i) * xs_impactionization(endash, n) * Psecondary(endash, energy_ev + ionpot_ev, ionpot_ev, J) * DELTA_E;
            }
          }

          break; // consider only the valence shell
        }
      }

      N_e += nnion * N_e_ion;
    }
  }

  // source term
  N_e += gsl_vector_get(sourcevec, get_energyindex_ev(energy_ev));

  return N_e;
}


static double get_frac_heating(int modelgridindex)
// Kozma & Fransson equation 3
{
  double frac_heating_Einit = 0.;
  const double nne = get_nne(modelgridindex);
  for (int i = 0; i < SFPTS; i++)
  {
    const double endash = gsl_vector_get(envec, i);

    if (endash >= E_0)
    {
      // first term
      frac_heating_Einit += gsl_vector_get(yvec, i) * (electron_loss_rate(endash * EV, nne) / EV) * DELTA_E;
    }
  }

  // second term
  frac_heating_Einit += E_0 * gsl_vector_get(yvec, get_energyindex_ev(E_0)) * (electron_loss_rate(E_0 * EV, nne) / EV);

  // third term (integral from zero to E_0)
  const int nsteps = 100;
  const double delta_endash = E_0 / nsteps;
  for (int j = 0; j < nsteps; j++)
  {
    const double endash = E_0 * j / nsteps;
    frac_heating_Einit += N_e(modelgridindex, endash * EV) * endash * delta_endash;
  }

  return frac_heating_Einit / E_init;
}


static double get_mean_binding_energy(int element, int ion)
{
  int q[M_NT_SHELLS];
  double total;

  const int ioncharge = get_ionstage(element,ion) - 1;
  const int nbound = elements[element].anumber - ioncharge; //number of bound electrons

  if (nbound > 0)
  {
    for (int electron_loop = 0; electron_loop < M_NT_SHELLS; electron_loop++)
    {
      q[electron_loop] = 0;
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
          exit(0);
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
          exit(0);
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
          exit(0);
        }
      }
    }

    //      printout("For element %d ion %d I got q's of: %d %d %d %d %d %d %d %d %d %d\n", element, ion, q[0], q[1], q[2], q[3], q[4], q[5], q[6], q[7], q[8], q[9]);
    //printout("%g %g %g %g %g %g %g %g %g %g\n", electron_binding[elements[element].anumber-1][0], electron_binding[elements[element].anumber-1][1], electron_binding[elements[element].anumber-1][2],electron_binding[elements[element].anumber-1][3],electron_binding[elements[element].anumber-1][4],electron_binding[elements[element].anumber-1][5],electron_binding[elements[element].anumber-1][6],electron_binding[elements[element].anumber-1][7],electron_binding[elements[element].anumber-1][8],electron_binding[elements[element].anumber-1][9]);

    total = 0.0;
    for (int electron_loop = 0; electron_loop < M_NT_SHELLS; electron_loop++)
    {
      const double use1 = q[electron_loop];
      if ((use1) > 0)
      {
        double use2 = electron_binding[elements[element].anumber-1][electron_loop];
        const double use3 = elements[element].ions[ion].ionpot;
        if (use2 > 0)
        {
          if (use2 < use3)
          {
            total += use1 / use3;
          }
          else
          {
            total += use1 / use2;
          }
        }
        else
        {
          use2 = electron_binding[elements[element].anumber-1][electron_loop-1];
          if (use2 < use3)
          {
            total += use1 / use3;
          }
          else
          {
            total += use1 / use2;
          }
          //		  total += use1/electron_binding[elements[element].anumber-1][electron_loop-1];
          if (electron_loop != 8)
          {
            //For some reason in the Lotz data, this is no energy for the M5 shell before Ni. So if the complaint
            //is for 8 (corresponding to that shell) then just use the M4 value
            printout("Huh? I'm trying to use a binding energy when I have no data. element %d ion %d\n",element,ion);
            exit(0);
          }
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


static double get_tot_nion(int modelgridindex)
{
  double result = 0.;
  for (int element = 0; element < nelements; element++)
  {
    result += modelgrid[modelgridindex].composition[element].abundance / elements[element].mass * get_rho(modelgridindex);

    //const int nions = get_nions(element);
    //for (ion = 0; ion < nions; ion++)
    //{
    //  result += ionstagepop(modelgridindex,element,ion);
    //}
  }

  return result;
}


static double get_oneoverw(int element, int ion, int modelgridindex)
{
  // Routine to compute the work per ion pair for doing the NT ionization calculation. Makes use of EXTREMELY SIMPLE approximations - high energy limits only */

  // Work in terms of 1/W since this is actually what we want. It is given by sigma/(Latom + Lelec).
  // We are going to start by taking all the high energy limits and ignoring Lelec, so that the
  // denominator is extremely simplified. Need to get the mean Z value.

  double Zbar = 0.0;
  for (int ielement = 0; ielement < nelements; ielement++)
  {
    Zbar += modelgrid[modelgridindex].composition[ielement].abundance * elements[ielement].anumber;
  }
  //printout("cell %d has Zbar of %g\n", modelgridindex, Zbar);

  double Aconst = 1.33e-14 * EV * EV;
  double binding = get_mean_binding_energy(element, ion);
  double oneoverW = Aconst * binding / Zbar / (2 * 3.14159 * pow(QE, 4));
  //printout("For element %d ion %d I got W of %g (eV)\n", element, ion, 1./oneoverW/EV);

  return oneoverW;
}


static double get_frac_excitation(int modelgridindex, int element, int ion)
// Kozma & Fransson equation 4, but summed over all transitions for given ion
{
  double frac_excitation = 0.;
  const double nnion = ionstagepop(modelgridindex, element, ion);

  // for (int level = 0; level < get_nlevels(element,ion); level++)
  const int level = 0; // just consider excitation from the ground level
  {
    const double epsilon_current = epsilon(element,ion,level);
    const int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;

    for (int t = 1; t <= nuptrans; t++)
    {
      const double epsilon_target = elements[element].ions[ion].levels[level].uptrans[t].epsilon;
      const int lineindex = elements[element].ions[ion].levels[level].uptrans[t].lineindex;
      const double epsilon_trans = epsilon_target - epsilon_current;
      const double epsilon_trans_ev = epsilon_trans / EV;

      double integral = 0.; // integral in Kozma & Fransson equation 9
      const int startindex = get_energyindex_ev(epsilon_trans_ev);
      for (int i = startindex; i < SFPTS; i++)
      {
        const double endash = gsl_vector_get(envec, i);
        integral += gsl_vector_get(yvec, i) * xs_excitation(lineindex, epsilon_trans, endash * EV) * DELTA_E;
      }

      frac_excitation += nnion * epsilon_trans_ev * integral / E_init;
    }
  }

  return frac_excitation;
}


static double get_frac_ionization(int modelgridindex, int element, int ion, int collionindex)
{
  double frac_ionization = 0.;
  const double nnion = ionstagepop(modelgridindex, element, ion); // hopefully ions per cm^3?

  const double ionpot_ev = colliondata[collionindex].ionpot_ev;

  double ionization_integral = 0.; // integral in Kozma & Fransson equation 10
  for (int i = 0; i < SFPTS; i++)
  {
    const double endash = gsl_vector_get(envec, i);
    if (endash >= ionpot_ev)
    {
      ionization_integral += gsl_vector_get(yvec, i) * xs_impactionization(endash, collionindex) * DELTA_E;
    }
  }

  frac_ionization = nnion * ionpot_ev * ionization_integral / E_init;

  return frac_ionization;
}


static double nt_ionization_ratecoeff_wfapprox(int modelgridindex, int element, int ion)
// this is actually a non-thermal ionization rate coefficient (multiply by population to get rate)
{
  const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
  // to get the non-thermal ionization rate we need to divide the energy deposited
  // per unit volume per unit time in the grid cell (sum of terms above)
  // by the total ion number density and the "work per ion pair"
  return deposition_rate_density / get_tot_nion(modelgridindex) * get_oneoverw(element, ion, modelgridindex);
}


static double get_eff_ionpot(int modelgridindex, int element, int ion)
// Kozma & Fransson 1992 equation 12
// result is in erg
{
  double eff_ionpot = 0.;
  const int Z = get_element(element);
  const int ionstage = get_ionstage(element, ion);
  const double nnion = ionstagepop(modelgridindex, element, ion); // hopefully ions per cm^3?
  for (int n = 0; n < colliondatacount; n++)
  {
    if (colliondata[n].Z == Z && colliondata[n].nelec == Z - ionstage + 1)
    {
      const double ionpot = colliondata[n].ionpot_ev * EV;
      const double frac_ionization = get_frac_ionization(modelgridindex, element, ion, n);
      const double tot_nion = get_tot_nion(modelgridindex);

      eff_ionpot = ionpot * (nnion / tot_nion) / frac_ionization;

      break; // consider only the valence shell
    }
  }
  return eff_ionpot;
}


static double nt_ionization_ratecoeff_sf(int modelgridindex, int element, int ion)
// Kozma & Fransson 1992 equation 13
{
  const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
  // const double gamma_deposition = rpkt_emiss[modelgridindex] * 1.e20 * FOURPI;
  return deposition_rate_density / get_tot_nion(modelgridindex) / get_eff_ionpot(modelgridindex, element, ion);
}


double nt_ionization_ratecoeff(int modelgridindex, int element, int ion)
{
  if (NT_SOLVE_SPENCERFANO)
    return nt_ionization_ratecoeff_sf(modelgridindex, element, ion);
  else
    return nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion);
}


void nt_solve_spencerfano(int modelgridindex, int timestep, double deltaV)
// solve the Spencer-Fano equation to get the non-thermal electron flux energy distribution
// based on Equation (2) of Li et al. (2012)
{
  printout("Solving Spencer-Fano equation for cell %d at timestep %d\n", modelgridindex, timestep);

  gsl_matrix *const sfmatrix = gsl_matrix_calloc(SFPTS, SFPTS);
  gsl_vector *const rhsvec = gsl_vector_calloc(SFPTS); // constant term (not dependent on y func) in each equation

  const double nne = get_nne(modelgridindex); // hopefully electrons per cm^3

  // loss terms and source terms
  for (int i = 0; i < SFPTS; i++)
  {
    const double en = gsl_vector_get(envec, i);

    *gsl_matrix_ptr(sfmatrix, i, i) += (electron_loss_rate(en * EV, nne) / EV);

    double source_integral_to_emax = 0.;
    for (int j = i; j < SFPTS; j++) // integral of source function from en to EMAX
    {
      source_integral_to_emax += gsl_vector_get(sourcevec, j) * DELTA_E;
    }
    gsl_vector_set(rhsvec, i, source_integral_to_emax);
  }

  E_0 = -1.; // reset E_0 to find it again

  for (int element = 0; element < nelements; element++)
  {
    const int Z = get_element(element);
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      const int ionstage = get_ionstage(element, ion);
      const double nnion = ionstagepop(modelgridindex,element,ion); // hopefully ions per cm^3?

      // excitation terms
      // for (int level = 0; level < get_nlevels(element,ion); level++)
      const int level = 0; // just consider excitation from the ground level
      {
        const double epsilon_current = epsilon(element,ion,level);
        const int nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
        for (int t = 1; t <= nuptrans; t++)
        {
          const double epsilon_target = elements[element].ions[ion].levels[level].uptrans[t].epsilon;
          const int lineindex = elements[element].ions[ion].levels[level].uptrans[t].lineindex;
          const double epsilon_trans = epsilon_target - epsilon_current;
          const double epsilon_trans_ev = epsilon_trans / EV;

          if (epsilon_trans / EV < E_0 || E_0 < 0.)
            E_0 = epsilon_trans / EV;

          for (int i = 0; i < SFPTS; i++)
          {
            const double en = gsl_vector_get(envec, i);
            const int stopindex = get_energyindex_ev(en + epsilon_trans_ev);

            for (int j = i; j <= stopindex; j++)
            {
              const double endash = gsl_vector_get(envec, j);
              *gsl_matrix_ptr(sfmatrix, i, j) += nnion * xs_excitation(lineindex, epsilon_trans, endash * EV) * DELTA_E;
            }
          }
        }
      }

      // ionization terms
      for (int n = 0; n < colliondatacount; n++)
      {
        if (colliondata[n].Z == Z && colliondata[n].nelec == Z - ionstage + 1)
        {
          const double ionpot_ev = colliondata[n].ionpot_ev;

          if (ionpot_ev < E_0 || E_0 < 0.)
            E_0 = ionpot_ev; // new minimum energy for excitation/ionization

          const double J = 0.6 * ionpot_ev;  // valid for elements other than He, Ne, Ar (Kozma & Fransson 1992)
          printout("Z %d, ion_stage %d, ionpot_ev %g\n", Z, ionstage, ionpot_ev);

          for (int i = 0; i < SFPTS; i++)
          {
            const double en = gsl_vector_get(envec, i);

            const int secondintegralstartindex = get_energyindex_ev(2 * en + ionpot_ev);

            for (int j = i; j < SFPTS; j++)
            {
                const double endash = gsl_vector_get(envec, j);

                const double en_primary = endash - ionpot_ev;
                const double prefactor = xs_impactionization(endash, n) / (J * atan(en_primary / (2 * J)));

                double epsilon_upper = (endash + en) / 2;
                double epsilon_lower = endash - en;
                // atan bit is the definite integral of 1/[1 + (epsilon - I)/J] in Kozma & Fransson 1992 equation 4
                *gsl_matrix_ptr(sfmatrix, i, j) += nnion * prefactor * J * (atan((epsilon_upper - ionpot_ev)/J) - atan((epsilon_lower - ionpot_ev)/J)) * DELTA_E;

                if (j >= secondintegralstartindex)
                {
                  epsilon_upper = (endash + ionpot_ev) / 2;
                  epsilon_lower = en + ionpot_ev;
                  *gsl_matrix_ptr(sfmatrix, i, j) -= nnion * prefactor * J * (atan((epsilon_upper - ionpot_ev)/J) - atan((epsilon_lower - ionpot_ev)/J)) * DELTA_E;
                }
            }
          }
          break; // include only the valence shell
        }
      }
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

  printout("Solving SF matrix equation\n");

  // make a copy of the matrix for the LU decomp
  gsl_matrix *sfmatrix_LU_decomp = gsl_matrix_alloc(SFPTS,SFPTS);
  gsl_matrix_memcpy(sfmatrix_LU_decomp, sfmatrix);

  gsl_permutation *p = gsl_permutation_alloc(SFPTS);

  int s; //sign of the transformation
  gsl_linalg_LU_decomp(sfmatrix_LU_decomp, p, &s);

  // solve matrix equation: sf_matrix * y_vec = constvec for yvec (population vector)
  gsl_linalg_LU_solve(sfmatrix_LU_decomp, p, rhsvec, yvec);

  const double TOLERANCE = 1e-40;

  double error_best = -1.;
  gsl_vector *yvec_best = gsl_vector_alloc(SFPTS); //population solution vector with lowest error
  gsl_vector *gsl_residual_vector = gsl_vector_alloc(SFPTS);
  gsl_vector *residual_vector = gsl_vector_alloc(SFPTS);
  int iteration;
  for (iteration = 0; iteration < 100; iteration++)
  {
    if (iteration > 0)
      gsl_linalg_LU_refine(sfmatrix, sfmatrix_LU_decomp, p, rhsvec, yvec, gsl_residual_vector);

    gsl_vector_memcpy(residual_vector, rhsvec);
    gsl_blas_dgemv(CblasNoTrans, 1.0, sfmatrix, yvec, -1.0, residual_vector); // calculate Ax - b = residual
    const double error = fabs(gsl_vector_get(residual_vector, gsl_blas_idamax(residual_vector))); // value of the largest absolute residual

    if (error < error_best || error_best < 0.)
    {
      gsl_vector_memcpy(yvec_best, yvec);
      error_best = error;
    }
    // printout("Linear algebra solver iteration %d has a maximum residual of %g\n",iteration,error);
    if (error < TOLERANCE)
    {
      break;
    }
  }
  if (error_best >= 0.)
  {
    printout("SF solver LU_refine: After %d iterations, keeping solution vector that had a max residual of %g\n",iteration,error_best);
    gsl_vector_memcpy(yvec,yvec_best);
  }
  gsl_vector_free(yvec_best);
  gsl_vector_free(gsl_residual_vector);
  gsl_vector_free(residual_vector);

  gsl_matrix_free(sfmatrix_LU_decomp);
  gsl_matrix_free(sfmatrix);
  gsl_vector_free(rhsvec);

  nonthermal_write_to_file(modelgridindex, timestep);

  const double deposition_rate_density = get_deposition_rate_density(modelgridindex);

  double frac_excitation_total = 0.;
  double frac_ionization_total = 0.;
  const int element = 0;
  const int Z = get_element(element);
  for (int ion = 0; ion < get_nions(element); ion++)
  {
    const int ionstage = get_ionstage(element, ion);
    int collionindex = 0;
    for (int n = 0; n < colliondatacount; n++)
    {
      if (colliondata[n].Z == Z && colliondata[n].nelec == Z - ionstage + 1)
      {
        collionindex = n;
        break;
      }
    }
    const double frac_ionization_ion = get_frac_ionization(modelgridindex, element, ion, collionindex);
    const double frac_excitation_ion = get_frac_excitation(modelgridindex, element, ion);
    frac_excitation_total += frac_excitation_ion;
    frac_ionization_total += frac_ionization_ion;

    printout("Z=%d ion_stage %d:\n", Z, ionstage);
    printout("  nnion: %g\n", ionstagepop(modelgridindex, element, ion));
    printout("  ionfract: %g\n", ionfract(element, ion, modelgridindex, get_nne(modelgridindex)));
    printout("  frac_ionization: %g\n", frac_ionization_ion);
    printout("  frac_excitation: %g\n", frac_excitation_ion);
    printout("  workfn:     %9.2f eV\n", (1. / get_oneoverw(element, ion, modelgridindex)) / EV);
    printout("  eff_ionpot: %9.2f eV\n", get_eff_ionpot(modelgridindex, element, ion) / EV);
    printout("  workfn approx Gamma: %9.3e\n",
             nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion));
    printout("  Spencer-Fano Gamma:  %9.3e\n",
             nt_ionization_ratecoeff_sf(modelgridindex, element, ion));
  }

  const double frac_heating = get_frac_heating(modelgridindex);

  // calculate number density of non-thermal electrons
  double nne_nt = 0.0;
  for (int i = 0; i < SFPTS; i++)
  {
    nne_nt += (deposition_rate_density / (E_init * EV)) * gsl_vector_get(yvec, i) * 1e2 * sqrt(9.10938e-31 / 2 / gsl_vector_get(envec, i) / 1.60218e-19) * DELTA_E;
  }


  printout("E_init:     %g eV/s/cm^3\n", E_init);
  printout("Deposition: %g eV/s/cm^3\n", deposition_rate_density / EV);
  // printout("Cell volume: %g cm^3\n", deltaV);F
  printout("nne:        %9.3e e-/cm^3\n", nne);
  printout("nne_nt:     %9.3e e-/cm^3\n", nne_nt);
  printout("nne_nt/nne: %9.3e\n", nne_nt / nne);
  printout("frac_heating:    %g\n", frac_heating);
  printout("frac_excitation: %g\n", frac_excitation_total);
  printout("frac_ionization: %g\n", frac_ionization_total);
  printout("Fractions sum to %g\n", frac_heating + frac_excitation_total + frac_ionization_total);
}

