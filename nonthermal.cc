#include "nonthermal.h"

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_double.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <tuple>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "ltepop.h"
#include "macroatom.h"
#include "packet.h"
#include "sn3d.h"
#include "stats.h"

namespace nonthermal {

namespace {

// THESE OPTIONS ARE USED TO TEST THE SF SOLVER
// Compare to Kozma & Fransson (1992) pure-oxygen plasma, nne = 1e8, x_e = 0.01
// #define yscalefactoroverride(mgi) (1e10)
// #define get_nnion_tot(x) (1e10)
// #define get_nnion(modelgridindex, element, ion) get_nnion_override(modelgridindex, element, ion)
// #define grid::get_nne(x) (1e8)
// #define SFPTS 10000  // number of energy points in the Spencer-Fano solution vector
// #define SF_EMAX 3000. // eV
// #define SF_EMIN 1. // eV
//
// double get_nnion_override(const int modelgridindex, const int element, const int ion)
// Fake the composition to test the NT solver
// {
//   const double nntot = get_nnion_tot(modelgridindex);
//   if (get_atomicnumber(element) == 8)
//   {
//     const int ionstage = get_ionstage(element, ion);
//     if (ionstage == 1)
//       return 0.99 * nntot;
//     else if (ionstage == 2)
//       return 0.01 * nntot;
//   }
//   return 0.;
// }

constexpr bool STORE_NT_SPECTRUM = false;  // if this is on, the non-thermal energy spectrum will be kept in memory for
                                           // every grid cell during packet propagation, which
                                           // can take up a lot of memory for large grid sizes
                                           // alternatively, just the non-thermal ionization rates can be stored
                                           // but we might want to re-enable this option to incorporate
                                           // non-thermal excitation rates if there are
                                           // many more transitions to store than there are NT spectrum samples

// minimum number fraction of the total population to include in SF solution
constexpr double minionfraction = 1.e-8;

// minimum deposition rate density (eV/s/cm^3) to solve SF equation
constexpr double MINDEPRATE = 0.;

// Bohr radius squared in cm^2
constexpr double A_naught_squared = 2.800285203e-17;

std::vector<std::vector<double>> electron_binding;
std::vector<std::vector<int>> shells_q;

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
  double auger_g_accumulated;  // track the statistical weight represented by the values below, so they can be updated
                               // with new g-weighted averaged values
  double prob_num_auger[NT_MAX_AUGER_ELECTRONS + 1];  // probability of 0, 1, ..., NT_MAX_AUGER_ELECTRONS Auger
                                                      // electrons being ejected when the shell is ionised
  float en_auger_ev;  // the average kinetic energy released in Auger electrons after making a hole in this shell
  float n_auger_elec_avg;
};

std::vector<collionrow> colliondata;

FILE *nonthermalfile{};
bool nonthermal_initialized = false;

gsl_vector *envec;      // energy grid on which solution is sampled
gsl_vector *logenvec;   // log of envec
gsl_vector *sourcevec;  // samples of the source function (energy distribution of deposited energy)
double E_init_ev =
    0;  // the energy injection rate density (and mean energy of injected electrons if source integral is one) in eV

constexpr double DELTA_E = (SF_EMAX - SF_EMIN) / (SFPTS - 1);

// Monte Carlo result - compare to analytical expectation
double nt_energy_deposited;

struct nt_excitation_struct {
  double frac_deposition;  // the fraction of the non-thermal deposition energy going to the excitation transition
  double ratecoeffperdeposition;  // the excitation rate coefficient divided by the deposition rate density
  int lineindex;
  int loweruptransindex;
};

struct nt_solution_struct {
  double *yfunc{};  // Samples of the Spencer-Fano solution function. Multiply by energy to get non-thermal
                    // electron number flux. y(E) * dE is the flux of electrons with energy in the range (E, E +
                    // dE) y has units of particles / cm2 / s / eV

  float frac_heating = 1.;     // energy fractions should add up to 1.0 if the solution is good
  float frac_ionization = 0.;  // fraction of deposition energy going to ionization
  float frac_excitation = 0.;  // fraction of deposition energy going to excitation

  // these points arrays of length includedions
  float *eff_ionpot{};  // these are used to calculate the non-thermal ionization rate
  double *fracdep_ionization_ion =
      nullptr;  // the fraction of the non-thermal deposition energy going to ionizing this ion

  // these  point to arrays of length includedions * (NT_MAX_AUGER_ELECTRONS + 1)
  float *prob_num_auger{};       // probability that one ionisation of this ion will produce n Auger electrons.
                                 // elements sum to 1.0 for a given ion
  float *ionenfrac_num_auger{};  // like above, but energy weighted. elements sum to 1.0 for an ion

  std::vector<nt_excitation_struct> frac_excitations_list;

  int timestep_last_solved = -1;     // the quantities above were calculated for this timestep
  float nneperion_when_solved{NAN};  // the nne when the solver was last run
};

nt_solution_struct *nt_solution;

double *deposition_rate_density;
int *deposition_rate_density_timestep;

void read_shell_configs() {
  auto shells_file = fstream_required("shells.txt", std::ios::in);

  int nshells = 0;      // number of shell in binding energy file
  int n_z_binding = 0;  // number of elements in file

  std::string line;
  assert_always(get_noncommentline(shells_file, line));
  std::istringstream(line) >> nshells >> n_z_binding;
  printout("Reading shells.txt with %d elements and %d shells\n", n_z_binding, nshells);

  shells_q.resize(n_z_binding, std::vector<int>(nshells, 0.));

  assert_always(shells_q.size() == electron_binding.size());
  assert_always(shells_q[0].size() == electron_binding[0].size());

  int elementcounter = 0;
  while (get_noncommentline(shells_file, line)) {
    std::istringstream ssline(line);

    int z_element = 0;
    assert_always(ssline >> z_element);
    printout("Reading shells Z=%d\n", z_element);

    for (int shell = 0; shell < nshells; shell++) {
      int q = 0;
      assert_always(ssline >> q);
      printout("q of %d in shell %d element number %d Z=%d\n", q, shell, elementcounter, z_element);
      shells_q[elementcounter][shell] = q;
    }
    elementcounter++;
  }
}

void read_binding_energies() {
  bool use_new_format = std::filesystem::exists("bindingenergies_lotz_tab1and2.txt") ||
                        std::filesystem::exists("data/bindingenergies_lotz_tab1and2.txt");

  int nshells = 0;      // number of shell in binding energy file
  int n_z_binding = 0;  // number of elements in binding energy file

  const auto *filename = use_new_format ? "bindingenergies_lotz_tab1and2.txt" : "binding_energies.txt";
  auto binding_energies_file = fstream_required(filename, std::ios::in);

  std::string line;
  assert_always(get_noncommentline(binding_energies_file, line));
  std::istringstream(line) >> nshells >> n_z_binding;
  printout("Reading binding energies file '%s' with %d elements and %d shells\n", filename, n_z_binding, nshells);

  electron_binding.resize(n_z_binding, std::vector<double>(nshells, 0.));

  for (int elemindex = 0; elemindex < n_z_binding; elemindex++) {
    assert_always(get_noncommentline(binding_energies_file, line));
    std::istringstream ssline(line);
    int z_element = elemindex + 1;
    /// new file as an atomic number column
    if (use_new_format) {
      ssline >> z_element;
    }
    for (int shell = 0; shell < nshells; shell++) {
      float bindingenergy = 0.;
      assert_always(ssline >> bindingenergy);
      electron_binding[elemindex][shell] = bindingenergy * EV;
    }
  }

  if (use_new_format) {
    read_shell_configs();
  }
}

auto get_auger_probability(int modelgridindex, int element, int ion, int naugerelec) -> double {
  assert_always(naugerelec <= NT_MAX_AUGER_ELECTRONS);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  return nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + naugerelec];
}

auto get_ion_auger_enfrac(int modelgridindex, int element, int ion, int naugerelec) -> double {
  assert_always(naugerelec <= NT_MAX_AUGER_ELECTRONS);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  return nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + naugerelec];
}

void check_auger_probabilities(int modelgridindex) {
  bool problem_found = false;

  for (int element = 0; element < get_nelements(); element++) {
    for (int ion = 0; ion < get_nions(element) - 1; ion++) {
      double prob_sum = 0.;
      double ionenfrac_sum = 0.;
      for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
        prob_sum += get_auger_probability(modelgridindex, element, ion, a);
        ionenfrac_sum += get_ion_auger_enfrac(modelgridindex, element, ion, a);
      }

      if (fabs(prob_sum - 1.0) > 0.001) {
        printout("Problem with Auger probabilities for cell %d Z=%d ionstage %d prob_sum %g\n", modelgridindex,
                 get_atomicnumber(element), get_ionstage(element, ion), prob_sum);
        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          printout("%d: %g\n", a, get_auger_probability(modelgridindex, element, ion, a));
        }
        problem_found = true;
      }

      if (fabs(ionenfrac_sum - 1.0) > 0.001) {
        printout("Problem with Auger energy frac sum for cell %d Z=%d ionstage %d ionenfrac_sum %g\n", modelgridindex,
                 get_atomicnumber(element), get_ionstage(element, ion), ionenfrac_sum);
        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          printout("%d: %g\n", a, get_ion_auger_enfrac(modelgridindex, element, ion, a));
        }
        problem_found = true;
      }
    }
  }

  if (problem_found) {
    std::abort();
  }
}

void read_auger_data() {
  printout("Reading Auger effect data...\n");
  FILE *augerfile = fopen_required("auger-km1993-table2.txt", "r");

  char line[151] = "";

  // map x-ray notation shells K L1 L2 L3 M1 M2 M3 to quantum numbers n and l
  const int xrayn[7] = {1, 2, 2, 2, 3, 3, 3};
  const int xrayl[7] = {0, 0, 1, 1, 0, 1, 1};
  const int xrayg[7] = {2, 2, 2, 4, 2, 2, 4};  // g statistical weight = 2j + 1

  while (feof(augerfile) == 0) {
    if (line != fgets(line, 151, augerfile)) {
      break;
    }

    int Z = -1;
    int ionstage = -1;
    int shellnum = -1;

    char *linepos = line;
    int offset = 0;

    assert_always(sscanf(linepos, "%d %d%n", &Z, &ionstage, &offset) == 2);
    assert_always(offset == 5);
    linepos += offset;

    const int element = get_elementindex(Z);

    if (element >= 0 && get_ionstage(element, 0) <= ionstage &&
        ionstage < (get_ionstage(element, 0) + get_nions(element))) {
      float ionpot_ev = -1;
      float en_auger_ev_total_nocorrection = -1;
      int epsilon_e3 = -1;

      assert_always(sscanf(linepos, "%d %g %g %d%n", &shellnum, &ionpot_ev, &en_auger_ev_total_nocorrection,
                           &epsilon_e3, &offset) == 4);
      assert_always(offset == 20);

      float n_auger_elec_avg = 0;
      double prob_num_auger[NT_MAX_AUGER_ELECTRONS + 1];
      for (int a = 0; a < 9; a++) {
        linepos = line + static_cast<ptrdiff_t>(26 + a * 5);
        // have to read out exactly 5 characters at a time because the columns are sometimes not separated by a space
        char strprob[6] = "00000";
        assert_always(sscanf(linepos, "%5c%n", strprob, &offset) == 1);
        assert_always(offset == 5);
        strprob[5] = '\0';

        int probnaugerelece4 = -1;
        assert_always(sscanf(strprob, "%d", &probnaugerelece4) == 1);

        const double probnaugerelec = probnaugerelece4 / 10000.;

        assert_always(probnaugerelec <= 1.0);

        n_auger_elec_avg += a * probnaugerelec;

        if (a <= NT_MAX_AUGER_ELECTRONS) {
          prob_num_auger[a] = probnaugerelec;
        } else {
          // add the rates of all higher ionisations to the top one
          prob_num_auger[NT_MAX_AUGER_ELECTRONS] += probnaugerelec;
        }
      }

      // use the epsilon correction factor as in equation 7 of Kaastra & Mewe (1993)
      float en_auger_ev = en_auger_ev_total_nocorrection - (epsilon_e3 / 1000. * ionpot_ev);

      const int n = xrayn[shellnum - 1];
      const int l = xrayl[shellnum - 1];
      const int g = xrayg[shellnum - 1];

      if (!std::isfinite(en_auger_ev) || en_auger_ev < 0) {
        printout("  WARNING: Z=%2d ionstage %2d shellnum %d en_auger_ev is %g. Setting to zero.\n", Z, ionstage,
                 shellnum, en_auger_ev);
        en_auger_ev = 0.;
      }

      // now loop through shells with impact ionisation cross sections and apply Auger data that matches n, l values
      for (auto &collionrow : colliondata) {
        if (collionrow.Z == Z && collionrow.nelec == (Z - ionstage + 1) && collionrow.n == n && collionrow.l == l) {
          printout(
              "Z=%2d ionstage %2d shellnum %d n %d l %d ionpot %7.2f E_A %8.1f E_A' %8.1f epsilon %6d <n_Auger> %5.1f "
              "P(n_Auger)",
              Z, ionstage, shellnum, n, l, ionpot_ev, en_auger_ev_total_nocorrection, en_auger_ev, epsilon_e3,
              n_auger_elec_avg);

          double prob_sum = 0.;
          for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
            prob_sum += prob_num_auger[a];
            printout(" %d: %4.2f", a, prob_num_auger[a]);
          }
          assert_always(fabs(prob_sum - 1.0) < 0.001);

          printout("\n");
          // printout("ionpot %g %g, g %d\n", colliondata[i].ionpot_ev, ionpot_ev, g);
          const bool found_existing_data = (collionrow.auger_g_accumulated > 0.);

          // keep existing data but update according to statistical weight represented by existing and new data
          const double oldweight = collionrow.auger_g_accumulated / (g + collionrow.auger_g_accumulated);
          const double newweight = g / (g + collionrow.auger_g_accumulated);
          collionrow.auger_g_accumulated += g;

          // update the statistical-weight averaged values
          collionrow.en_auger_ev = oldweight * collionrow.en_auger_ev + newweight * en_auger_ev;
          collionrow.n_auger_elec_avg = oldweight * collionrow.n_auger_elec_avg + newweight * n_auger_elec_avg;

          prob_sum = 0.;
          for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
            collionrow.prob_num_auger[a] = oldweight * collionrow.prob_num_auger[a] + newweight * prob_num_auger[a];
            prob_sum += collionrow.prob_num_auger[a];
          }
          assert_always(fabs(prob_sum - 1.0) < 0.001);

          if (found_existing_data) {
            printout("  same NL shell already has data from another X-ray shell. New g-weighted values: P(n_Auger)");

            for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
              printout(" %d: %4.2f", a, collionrow.prob_num_auger[a]);
            }
            printout("\n");
          }
        }
      }
    }
  }
  fclose(augerfile);
}

void read_collion_data() {
  printout("Reading collisional ionization data...\n");

  FILE *cifile = fopen_required("collion.txt", "r");
  int colliondatacount = 0;
  assert_always(fscanf(cifile, "%d", &colliondatacount) == 1);
  printout("Reading %d collisional transition rows\n", colliondatacount);
  assert_always(colliondatacount > 0);

  for (int i = 0; i < colliondatacount; i++) {
    collionrow collionrow{};
    assert_always(fscanf(cifile, "%2d %2d %1d %1d %lg %lg %lg %lg %lg", &collionrow.Z, &collionrow.nelec, &collionrow.n,
                         &collionrow.l, &collionrow.ionpot_ev, &collionrow.A, &collionrow.B, &collionrow.C,
                         &collionrow.D) == 9);

    const int element = get_elementindex(collionrow.Z);
    const int ionstage = collionrow.Z - collionrow.nelec + 1;
    if (element < 0 || ionstage < get_ionstage(element, 0) ||
        ionstage > get_ionstage(element, get_nions(element) - 1)) {
      continue;
    }
    collionrow.prob_num_auger[0] = 1.;
    for (int a = 1; a <= NT_MAX_AUGER_ELECTRONS; a++) {
      collionrow.prob_num_auger[a] = 0.;
    }

    collionrow.auger_g_accumulated = 0.;
    collionrow.en_auger_ev = 0.;
    collionrow.n_auger_elec_avg = 0.;

    colliondata.push_back(collionrow);

    // printout("ci row: %2d %2d %1d %1d %lg %lg %lg %lg %lg\n", collionrow.Z, collionrow.nelec, collionrow.n,
    //          collionrow.l, collionrow.ionpot_ev, collionrow.A, collionrow.B, collionrow.C, collionrow.D);
  }
  printout("Stored %zu of %d input shell cross sections\n", colliondata.size(), colliondatacount);

  fclose(cifile);

  if (NT_MAX_AUGER_ELECTRONS > 0) {
    read_auger_data();
  }
}

void zero_all_effionpot(const int modelgridindex) {
  assert_always(nt_solution[modelgridindex].prob_num_auger);
  assert_always(nt_solution[modelgridindex].ionenfrac_num_auger);

  for (int uniqueionindex = 0; uniqueionindex < get_includedions(); uniqueionindex++) {
    nt_solution[modelgridindex].eff_ionpot[uniqueionindex] = 0.;

    nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1)] = 1.;
    nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1)] = 1.;
    for (int a = 1; a <= NT_MAX_AUGER_ELECTRONS; a++) {
      nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 0.;
      nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 0.;
    }

    const auto [element, ion] = get_ionfromuniqueionindex(uniqueionindex);
    assert_always(fabs(get_auger_probability(modelgridindex, element, ion, 0) - 1.0) < 1e-3);
    assert_always(fabs(get_ion_auger_enfrac(modelgridindex, element, ion, 0) - 1.0) < 1e-3);
  }
  check_auger_probabilities(modelgridindex);
}

auto get_energyindex_ev_lteq(const double energy_ev) -> int
// finds the highest energy point <= energy_ev
{
  const int index = floor((energy_ev - SF_EMIN) / DELTA_E);

  if (index < 0) {
    return 0;
  }
  if (index > SFPTS - 1) {
    return SFPTS - 1;
  }
  return index;
}

auto get_energyindex_ev_gteq(const double energy_ev) -> int
// finds the highest energy point <= energy_ev
{
  const int index = ceil((energy_ev - SF_EMIN) / DELTA_E);

  if (index < 0) {
    return 0;
  }
  if (index > SFPTS - 1) {
    return SFPTS - 1;
  }
  return index;
}

auto get_y_sample(const int modelgridindex, const int index) -> double {
  if (nt_solution[modelgridindex].yfunc != nullptr) {
    if (!std::isfinite(nt_solution[modelgridindex].yfunc[index])) {
      printout("get_y_sample index %d %g\n", index, nt_solution[modelgridindex].yfunc[index]);
    }
    return nt_solution[modelgridindex].yfunc[index];
  }
  printout("non-thermal: attempted to get y function sample index %d in cell %d, but the y array pointer is null\n",
           index, modelgridindex);
  std::abort();
  return -1;
}

auto get_y(const int modelgridindex, const double energy_ev) -> double {
  if (energy_ev <= 0) {
    return 0.;
  }

  const int index = static_cast<int>((energy_ev - SF_EMIN) / DELTA_E);

  // assert_always(index > 0);
  if (index < 0) {
    // return 0.;
    assert_always(std::isfinite(get_y_sample(modelgridindex, 0)));
    return get_y_sample(modelgridindex, 0);
  }
  if (index > SFPTS - 1) {
    return 0.;
  }
  const double enbelow = gsl_vector_get(envec, index);
  const double enabove = gsl_vector_get(envec, index + 1);
  const double ybelow = get_y_sample(modelgridindex, index);
  const double yabove = get_y_sample(modelgridindex, index + 1);
  const double x = (energy_ev - enbelow) / (enabove - enbelow);
  return (1 - x) * ybelow + x * yabove;

  // or return the nearest neighbour
  // return get_y_sample(modelgridindex, index);
}

void nt_write_to_file(const int modelgridindex, const int timestep, const int iteration) {
#ifdef _OPENMP
#pragma omp critical(nonthermal_out_file)
  {
#endif
    if (!nonthermal_initialized || nonthermalfile == nullptr) {
      printout("Call to nonthermal_write_to_file before nonthermal_init");
      std::abort();
    }

    static long nonthermalfile_offset_iteration_zero = 0;
#ifdef _OPENMP
#pragma omp threadprivate(nonthermalfile_offset_iteration_zero)
#endif
    {
      if (iteration == 0) {
        nonthermalfile_offset_iteration_zero = ftell(nonthermalfile);
      } else {
        // overwrite the non-thermal spectrum of a previous iteration of the same timestep and gridcell
        assert_always(fseek(nonthermalfile, nonthermalfile_offset_iteration_zero, SEEK_SET) == 0);
      }
    }

#ifndef yscalefactoroverride  // manual override can be defined
    const double yscalefactor = (get_deposition_rate_density(modelgridindex) / (E_init_ev * EV));
#else
  const double yscalefactor = yscalefactoroverride(modelgridindex);
#endif

    for (int s = 0; s < SFPTS; s++) {
      fprintf(nonthermalfile, "%d %d %d %.5e %.5e %.5e\n", timestep, modelgridindex, s, gsl_vector_get(envec, s),
              gsl_vector_get(sourcevec, s), yscalefactor * get_y_sample(modelgridindex, s));
    }
    fflush(nonthermalfile);
#ifdef _OPENMP
  }
#endif
}

auto get_xs_ionization_vector(gsl_vector *const xs_vec, const collionrow &colliondata) -> int
// xs_vec will be set with impact ionization cross sections for E > ionpot_ev (and zeros below this energy)
{
  const double ionpot_ev = colliondata.ionpot_ev;
  const int startindex = get_energyindex_ev_gteq(ionpot_ev);

  // en points for which en < ionpot
  for (int i = 0; i < startindex; i++) {
    gsl_vector_set(xs_vec, i, 0.);
  }

  const double A = colliondata.A;
  const double B = colliondata.B;
  const double C = colliondata.C;
  const double D = colliondata.D;

  for (int i = startindex; i < SFPTS; i++) {
    const double u = gsl_vector_get(envec, i) / ionpot_ev;
    const double xs_ioniz =
        1e-14 * (A * (1 - 1 / u) + B * pow((1 - 1 / u), 2) + C * log(u) + D * log(u) / u) / (u * pow(ionpot_ev, 2));
    gsl_vector_set(xs_vec, i, xs_ioniz);
  }

  return startindex;
}

auto Psecondary(const double e_p, const double epsilon, const double I, const double J) -> double
// distribution of secondary electron energies for primary electron with energy e_p
// Opal, Peterson, & Beaty (1971)
{
  const double e_s = epsilon - I;

  if (e_p <= I || e_s < 0.) {
    return 0.;
  }
  assert_always(J > 0);
  assert_always(e_p >= I);
  assert_always(e_s >= 0);
  assert_always(std::isfinite(atan((e_p - I) / 2 / J)));
  return 1 / (J * atan((e_p - I) / 2 / J) * (1 + pow(e_s / J, 2)));
}

auto get_J(const int Z, const int ionstage, const double ionpot_ev) -> double {
  // returns an energy in eV
  // values from Opal et al. 1971 as applied by Kozma & Fransson 1992
  if (ionstage == 1) {
    if (Z == 2) {  // He I
      return 15.8;
    }
    if (Z == 10) {  // Ne I
      return 24.2;
    }
    if (Z == 18) {  // Ar I
      return 10.;
    }
  }

  return 0.6 * ionpot_ev;
}

constexpr auto xs_excitation(const int element, const int ion, const int lower, const int uptransindex,
                             const double epsilon_trans, const double lowerstatweight, const double energy) -> double
// collisional excitation cross section in cm^2
// energies are in erg
{
  if (energy < epsilon_trans) {
    return 0.;
  }

  const double coll_strength = globals::elements[element].ions[ion].levels[lower].uptrans[uptransindex].coll_str;
  if (coll_strength >= 0) {
    // collision strength is available, so use it
    // Li et al. 2012 equation 11
    return pow(H_ionpot / energy, 2) / lowerstatweight * coll_strength * PI * A_naught_squared;
  }
  const bool forbidden = globals::elements[element].ions[ion].levels[lower].uptrans[uptransindex].forbidden;
  if (!forbidden) {
    const double trans_osc_strength =
        globals::elements[element].ions[ion].levels[lower].uptrans[uptransindex].osc_strength;
    // permitted E1 electric dipole transitions
    const double U = energy / epsilon_trans;

    // const double g_bar = 0.2;
    const double A = 0.28;
    const double B = 0.15;
    const double g_bar = A * log(U) + B;

    const double prefactor = 45.585750051;  // 8 * pi^2/sqrt(3)
    // Eq 4 of Mewe 1972, possibly from Seaton 1962?
    return prefactor * A_naught_squared * pow(H_ionpot / epsilon_trans, 2) * trans_osc_strength * g_bar / U;
  }
  return 0.;
}

constexpr auto electron_loss_rate(const double energy, const double nne) -> double
// -dE / dx for fast electrons
// energy is in ergs
// nne is the thermal electron density [cm^-3]
// return units are erg / cm
{
  if (energy <= 0.) {
    return 0;
  }

  // normally set to 1.0, but Shingles et al. (2021) boosted this to increase heating
  const double boostfactor = 1.;

  const double omegap = sqrt(4 * PI * nne * pow(QE, 2) / ME);
  const double zetae = H * omegap / 2 / PI;
  if (energy > 14 * EV) {
    return boostfactor * nne * 2 * PI * pow(QE, 4) / energy * log(2 * energy / zetae);
  }
  const double v = sqrt(2 * energy / ME);
  return boostfactor * nne * 2 * PI * pow(QE, 4) / energy * log(ME * pow(v, 3) / (EULERGAMMA * pow(QE, 2) * omegap));
}

constexpr auto xs_impactionization(const double energy_ev, const collionrow &colliondata) -> double
// impact ionization cross section in cm^2
// energy and ionization_potential should be in eV
// fitting forumula of Younger 1981
// called Q_i(E) in KF92 equation 7
{
  const double ionpot_ev = colliondata.ionpot_ev;
  const double u = energy_ev / ionpot_ev;

  if (u <= 1.) {
    return 0;
  }
  const double A = colliondata.A;
  const double B = colliondata.B;
  const double C = colliondata.C;
  const double D = colliondata.D;

  return 1e-14 * (A * (1 - 1 / u) + B * pow((1 - 1 / u), 2) + C * log(u) + D * log(u) / u) / (u * pow(ionpot_ev, 2));
}

auto N_e(const int modelgridindex, const double energy) -> double
// Kozma & Fransson equation 6.
// Something related to a number of electrons, needed to calculate the heating fraction in equation 3
// not valid for energy > SF_EMIN
{
  const double energy_ev = energy / EV;
  const double tot_nion = get_nnion_tot(modelgridindex);
  double N_e = 0.;

  for (int element = 0; element < get_nelements(); element++) {
    const int Z = get_atomicnumber(element);
    const int nions = get_nions(element);

    for (int ion = 0; ion < nions; ion++) {
      double N_e_ion = 0.;
      const int ionstage = get_ionstage(element, ion);
      const double nnion = get_nnion(modelgridindex, element, ion);

      if (nnion < minionfraction * tot_nion) {  // skip negligible ions
        continue;
      }

      // excitation terms

      const int nlevels_all = get_nlevels(element, ion);
      const int nlevels = std::min(NTEXCITATION_MAXNLEVELS_LOWER, nlevels_all);

      for (int lower = 0; lower < nlevels; lower++) {
        const int nuptrans = get_nuptrans(element, ion, lower);
        const double nnlevel = get_levelpop(modelgridindex, element, ion, lower);
        const double epsilon_lower = epsilon(element, ion, lower);
        const auto statweight_lower = stat_weight(element, ion, lower);
        for (int t = 0; t < nuptrans; t++) {
          const int lineindex = globals::elements[element].ions[ion].levels[lower].uptrans[t].lineindex;
          const int upper = globals::linelist[lineindex].upperlevelindex;
          if (upper >= NTEXCITATION_MAXNLEVELS_UPPER) {
            continue;
          }
          const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
          const double epsilon_trans_ev = epsilon_trans / EV;
          N_e_ion += (nnlevel / nnion) * get_y(modelgridindex, energy_ev + epsilon_trans_ev) *
                     xs_excitation(element, ion, lower, t, epsilon_trans, statweight_lower, energy + epsilon_trans);
        }
      }

      // ionization terms
      for (const auto &collionrow : colliondata) {
        if (collionrow.Z == Z && collionrow.nelec == Z - ionstage + 1) {
          const double ionpot_ev = collionrow.ionpot_ev;
          const double J = get_J(Z, ionstage, ionpot_ev);
          const double lambda = std::min(SF_EMAX - energy_ev, energy_ev + ionpot_ev);

          const int integral1startindex = get_energyindex_ev_lteq(ionpot_ev);
          const int integral1stopindex = get_energyindex_ev_lteq(lambda);

          // integral from ionpot up to lambda
          for (int i = integral1startindex; i <= integral1stopindex; i++) {
            const double endash = gsl_vector_get(envec, i);
            const double delta_endash = DELTA_E;

            N_e_ion += get_y(modelgridindex, energy_ev + endash) * xs_impactionization(energy_ev + endash, collionrow) *
                       Psecondary(energy_ev + endash, endash, ionpot_ev, J) * delta_endash;
          }

          // integral from 2E + I up to E_max
          const int integral2startindex = get_energyindex_ev_lteq(2 * energy_ev + ionpot_ev);
          for (int i = integral2startindex; i < SFPTS; i++) {
            const double endash = gsl_vector_get(envec, i);
            const double delta_endash = DELTA_E;
            N_e_ion += get_y_sample(modelgridindex, i) * xs_impactionization(endash, collionrow) *
                       Psecondary(endash, energy_ev + ionpot_ev, ionpot_ev, J) * delta_endash;
          }
        }
      }

      N_e += nnion * N_e_ion;
    }
  }

  // source term, should be zero at the low end anyway
  N_e += gsl_vector_get(sourcevec, get_energyindex_ev_lteq(energy_ev));

  assert_always(std::isfinite(N_e));
  return N_e;
}

auto calculate_frac_heating(const int modelgridindex) -> float
// Kozma & Fransson equation 3
{
  // frac_heating multiplied by E_init, which will be divided out at the end
  double frac_heating_Einit = 0.;

  const float nne = grid::get_nne(modelgridindex);
  // const float nnetot = grid::get_nnetot(modelgridindex);

  for (int i = 0; i < SFPTS; i++) {
    const double endash = gsl_vector_get(envec, i);

    // first term
    frac_heating_Einit += get_y_sample(modelgridindex, i) * (electron_loss_rate(endash * EV, nne) / EV) * DELTA_E;
  }

  // second term
  frac_heating_Einit += SF_EMIN * get_y(modelgridindex, SF_EMIN) * (electron_loss_rate(SF_EMIN * EV, nne) / EV);

  double N_e_contrib = 0.;
  // third term (integral from zero to SF_EMIN)
  const int nsteps = static_cast<int>(ceil(SF_EMIN / DELTA_E) * 10);
  assert_always(nsteps > 0);
  const double delta_endash = SF_EMIN / nsteps;
  for (int j = 0; j < nsteps; j++) {
    const double endash = SF_EMIN * j / nsteps;
    N_e_contrib += N_e(modelgridindex, endash * EV) * endash * delta_endash;
  }
  frac_heating_Einit += N_e_contrib;
  printout(" heating N_e contrib (en < EMIN) %g nsteps %d\n", N_e_contrib / E_init_ev, nsteps);

  const float frac_heating = frac_heating_Einit / E_init_ev;

  if (!std::isfinite(frac_heating) || frac_heating < 0 || frac_heating > 1.0) {
    printout("WARNING: calculate_frac_heating: invalid result of %g. Setting to 1.0 instead\n", frac_heating);
    return 1.;
  }

  return frac_heating;
}

auto get_nt_frac_ionization(const int modelgridindex) -> float {
  if (!NT_ON) {
    return 0.;
  }
  if (!NT_SOLVE_SPENCERFANO) {
    return 0.03;
  }

  const float frac_ionization = nt_solution[modelgridindex].frac_ionization;

  if (frac_ionization < 0 || !std::isfinite(frac_ionization)) {
    printout("ERROR: get_nt_frac_ionization called with no valid solution stored for cell %d. frac_ionization = %g\n",
             modelgridindex, frac_ionization);
    std::abort();
  }

  return frac_ionization;
}

auto get_nt_frac_excitation(const int modelgridindex) -> float {
  if (!NT_ON || !NT_SOLVE_SPENCERFANO) {
    return 0.;
  }

  const float frac_excitation = nt_solution[modelgridindex].frac_excitation;

  if (frac_excitation < 0 || !std::isfinite(frac_excitation)) {
    printout("ERROR: get_nt_frac_excitation called with no valid solution stored for cell %d. frac_excitation = %g\n",
             modelgridindex, frac_excitation);
    std::abort();
  }

  return frac_excitation;
}

auto get_mean_binding_energy(const int element, const int ion) -> double {
  const int ioncharge = get_ionstage(element, ion) - 1;
  const int nbound = get_atomicnumber(element) - ioncharge;  // number of bound electrons

  if (nbound <= 0) {
    return 0.;
  }

  const int num_shells = electron_binding[get_atomicnumber(element) - 1].size();
  auto q = std::vector<int>(num_shells, 0);

  bool use_shells_file = std::filesystem::exists("shells.txt") || std::filesystem::exists("data/shells.txt");
  if (!use_shells_file) {
    for (int electron_loop = 0; electron_loop < nbound; electron_loop++) {
      if (q[0] < 2)  // K 1s
      {
        q[0]++;
      } else if (q[1] < 2)  // L1 2s
      {
        q[1]++;
      } else if (q[2] < 2)  // L2 2p[1/2]
      {
        q[2]++;
      } else if (q[3] < 4)  // L3 2p[3/2]
      {
        q[3]++;
      } else if (q[4] < 2)  // M1 3s
      {
        q[4]++;
      } else if (q[5] < 2)  // M2 3p[1/2]
      {
        q[5]++;
      } else if (q[6] < 4)  // M3 3p[3/2]
      {
        q[6]++;
      } else if (ioncharge == 0) {
        if (q[9] < 2)  // N1 4s
        {
          q[9]++;
        } else if (q[7] < 4)  // M4 3d[3/2]
        {
          q[7]++;
        } else if (q[8] < 6)  // M5 3d[5/2]
        {
          q[8]++;
        } else {
          printout("Going beyond the 4s shell in NT calculation. Abort!\n");
          std::abort();
        }
      } else if (ioncharge == 1) {
        if (q[9] < 1)  // N1 4s
        {
          q[9]++;
        } else if (q[7] < 4)  // M4 3d[3/2]
        {
          q[7]++;
        } else if (q[8] < 6)  // M5 3d[5/2]
        {
          q[8]++;
        } else {
          printout("Going beyond the 4s shell in NT calculation. Abort!\n");
          std::abort();
        }
      } else if (ioncharge > 1) {
        if (q[7] < 4)  // M4 3d[3/2]
        {
          q[7]++;
        } else if (q[8] < 6)  // M5 3d[5/2]
        {
          q[8]++;
        } else {
          printout("Going beyond the 4s shell in NT calculation. Abort!\n");
          std::abort();
        }
      }
    }
  }

  //      printout("For element %d ion %d I got q's of: %d %d %d %d %d %d %d %d %d %d\n", element, ion, q[0], q[1],
  //      q[2], q[3], q[4], q[5], q[6], q[7], q[8], q[9]);
  // printout("%g %g %g %g %g %g %g %g %g %g\n", electron_binding[get_atomicnumber(element)-1][0],
  // electron_binding[get_atomicnumber(element)-1][1],
  // electron_binding[get_atomicnumber(element)-1][2],electron_binding[get_atomicnumber(element)-1][3],electron_binding[get_atomicnumber(element)-1][4],electron_binding[get_atomicnumber(element)-1][5],electron_binding[get_atomicnumber(element)-1][6],electron_binding[get_atomicnumber(element)-1][7],electron_binding[get_atomicnumber(element)-1][8],electron_binding[get_atomicnumber(element)-1][9]);

  double total = 0.;
  for (int shellindex = 0; shellindex < num_shells; shellindex++) {
    double electronsinshell = (use_shells_file ? shells_q[get_atomicnumber(element) - 1] : q)[shellindex];

    if (electronsinshell <= 0) {
      continue;
    }
    double enbinding = electron_binding[get_atomicnumber(element) - 1][shellindex];
    const double ionpot = globals::elements[element].ions[ion].ionpot;
    if (enbinding <= 0) {
      enbinding = electron_binding[get_atomicnumber(element) - 1][shellindex - 1];
      //  to get total += electronsinshell/electron_binding[get_atomicnumber(element)-1][electron_loop-1];
      //  set use3 = 0.
      if (shellindex != 8) {
        // For some reason in the Lotz data, this is no energy for the M5 shell before Ni. So if the complaint
        // is for 8 (corresponding to that shell) then just use the M4 value
        printout("Huh? I'm trying to use a binding energy when I have no data. element %d ion %d\n", element, ion);
        printout("Z = %d, ionstage = %d\n", get_atomicnumber(element), get_ionstage(element, ion));
        std::abort();
      }
    }
    total += electronsinshell / std::max(ionpot, enbinding);

    // printout("total %g\n", total);
  }

  return total;
}

auto get_oneoverw(const int element, const int ion, const int modelgridindex) -> double {
  // Routine to compute the work per ion pair for doing the NT ionization calculation.
  // Makes use of EXTREMELY SIMPLE approximations - high energy limits only

  // Work in terms of 1/W since this is actually what we want. It is given by sigma/(Latom + Lelec).
  // We are going to start by taking all the high energy limits and ignoring Lelec, so that the
  // denominator is extremely simplified. Need to get the mean Z value.

  double Zbar = 0.;  // mass-weighted average atomic number
  for (int ielement = 0; ielement < get_nelements(); ielement++) {
    Zbar += grid::get_elem_abundance(modelgridindex, ielement) * get_atomicnumber(ielement);
  }
  // printout("cell %d has Zbar of %g\n", modelgridindex, Zbar);

  const double Aconst = 1.33e-14 * EV * EV;
  const double binding = get_mean_binding_energy(element, ion);
  const double oneoverW = Aconst * binding / Zbar / (2 * PI * pow(QE, 4));
  // printout("For element %d ion %d I got W of %g (eV)\n", element, ion, 1./oneoverW/EV);

  return oneoverW;
}

auto calculate_nt_frac_ionization_shell(const int modelgridindex, const int element, const int ion,
                                        const collionrow &collionrow) -> double
// the fraction of deposition energy that goes into ionising electrons in this particular shell
{
  const double nnion = get_nnion(modelgridindex, element, ion);  // hopefully ions per cm^3?
  const double ionpot_ev = collionrow.ionpot_ev;

  gsl_vector *cross_section_vec = gsl_vector_alloc(SFPTS);
  get_xs_ionization_vector(cross_section_vec, collionrow);

  // either multiply by the variable delta_e for LOG_E spacing...

  gsl_vector_view const yvecview = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);

  double y_dot_crosssection_de = 0.;
  gsl_blas_ddot(&yvecview.vector, cross_section_vec, &y_dot_crosssection_de);
  gsl_vector_free(cross_section_vec);

  // or multiply the scalar result by the constant DELTA_E
  y_dot_crosssection_de *= DELTA_E;

  return nnion * ionpot_ev * y_dot_crosssection_de / E_init_ev;
}

auto nt_ionization_ratecoeff_wfapprox(const int modelgridindex, const int element, const int ion) -> double
// non-thermal ionization rate coefficient (multiply by population to get rate)
{
  const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
  // to get the non-thermal ionization rate we need to divide the energy deposited
  // per unit volume per unit time in the grid cell (sum of terms above)
  // by the total ion number density and the "work per ion pair"
  return deposition_rate_density / get_nnion_tot(modelgridindex) * get_oneoverw(element, ion, modelgridindex);
}

auto calculate_nt_ionization_ratecoeff(const int modelgridindex, const int element, const int ion,
                                       const bool assumeshellpotentialisvalence) -> double
// Integrate the ionization cross section over the electron degradation function to get the ionization rate coefficient
// i.e. multiply this by ion population to get a rate of ionizations per second
// Do not call during packet propagation, as the y vector may not be in memory!
// IMPORTANT: we are dividing by the shell potential, not the valence potential here!
// To change this set assumeshellpotentialisvalence to true
{
  gsl_vector *cross_section_vec = gsl_vector_alloc(SFPTS);
  gsl_vector *cross_section_vec_allshells = gsl_vector_calloc(SFPTS);

  const int Z = get_atomicnumber(element);
  const int ionstage = get_ionstage(element, ion);
  double ionpot_valence = -1;

  for (auto &collionrow : colliondata) {
    if (collionrow.Z == Z && collionrow.nelec == Z - ionstage + 1) {
      get_xs_ionization_vector(cross_section_vec, collionrow);

      if (assumeshellpotentialisvalence) {
        const double ionpot_shell = collionrow.ionpot_ev * EV;
        if (ionpot_valence < 0) {
          ionpot_valence = ionpot_shell;
        }

        // ensure that the first shell really was the valence shell (we assumed ascending energy order)
        assert_always(ionpot_shell >= ionpot_valence);

        // boost the ionization rate by assuming shell vacancy energy is used to eject valence electrons
        gsl_vector_scale(cross_section_vec, ionpot_shell / ionpot_valence);
      }

      gsl_vector_add(cross_section_vec_allshells, cross_section_vec);
    }
  }

  gsl_vector_free(cross_section_vec);

  assert_always(nt_solution[modelgridindex].yfunc != nullptr);

  double y_dot_crosssection_de = 0.;
  gsl_vector_view const yvecview_thismgi = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);
  gsl_blas_ddot(&yvecview_thismgi.vector, cross_section_vec_allshells, &y_dot_crosssection_de);
  gsl_vector_free(cross_section_vec_allshells);

  y_dot_crosssection_de *= DELTA_E;

  const double deposition_rate_density_ev = get_deposition_rate_density(modelgridindex) / EV;
  const double yscalefactor = deposition_rate_density_ev / E_init_ev;

  return yscalefactor * y_dot_crosssection_de;
}

void calculate_eff_ionpot_auger_rates(const int modelgridindex, const int element, const int ion)
// Kozma & Fransson 1992 equation 12, except modified to be a sum over all shells of an ion
// the result is in ergs
{
  const int Z = get_atomicnumber(element);
  const int ionstage = get_ionstage(element, ion);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  const double nnion = get_nnion(modelgridindex, element, ion);  // ions/cm^3
  const double tot_nion = get_nnion_tot(modelgridindex);
  const double X_ion = nnion / tot_nion;  // molar fraction of this ion

  // The ionization rates of all shells of an ion add to make the ion's total ionization rate,
  // i.e., Gamma_ion = Gamma_shell_a + Gamma_shell_b + ...
  // And since the ionization rate is inversely proportional to the effective ion potential,
  // we solve:
  // (eta_ion / ionpot_ion) = (eta_shell_a / ionpot_shell_a) + (eta_shell_b / ionpot_shell_b) + ...
  // where eta is the fraction of the deposition energy going into ionization of the ion or shell

  std::array<double, NT_MAX_AUGER_ELECTRONS + 1> eta_nauger_ionize_over_ionpot_sum{};
  std::array<double, NT_MAX_AUGER_ELECTRONS + 1> eta_nauger_ionize_sum{};

  std::fill_n(&nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1)],
              NT_MAX_AUGER_ELECTRONS + 1, 0.);

  std::fill_n(&nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1)],
              NT_MAX_AUGER_ELECTRONS + 1, 0.);

  double eta_over_ionpot_sum = 0.;
  double eta_sum = 0.;
  double ionpot_valence = -1;
  int matching_nlsubshell_count = 0;
  for (auto &collionrow : colliondata) {
    if (collionrow.Z == Z && collionrow.nelec == Z - ionstage + 1) {
      matching_nlsubshell_count++;
      const double frac_ionization_shell = calculate_nt_frac_ionization_shell(modelgridindex, element, ion, collionrow);
      eta_sum += frac_ionization_shell;
      const double ionpot_shell = collionrow.ionpot_ev * EV;

      if (ionpot_valence < 0) {
        ionpot_valence = ionpot_shell;
      }

      // ensure that the first shell really was the valence shell (we assumed ascending energy order)
      assert_always(ionpot_shell >= ionpot_valence);

      const double ionpot = NT_USE_VALENCE_IONPOTENTIAL ? ionpot_valence : ionpot_shell;
      const double eta_over_ionpot = frac_ionization_shell / ionpot;  // this is proportional to rate

      eta_over_ionpot_sum += eta_over_ionpot;

      for (size_t a = 0; a < eta_nauger_ionize_sum.size(); a++) {
        eta_nauger_ionize_over_ionpot_sum[a] += eta_over_ionpot * collionrow.prob_num_auger[a];
        eta_nauger_ionize_sum[a] += frac_ionization_shell * collionrow.prob_num_auger[a];
      }
    }
  }

  if (NT_MAX_AUGER_ELECTRONS > 0 && matching_nlsubshell_count > 0) {
    const int nions = get_nions(element);
    const int topion = nions - 1;
    if (ion < topion)  // don't try to ionise the top ion
    {
      for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
        const int upperion = ion + 1 + a;
        if (upperion <= topion)  // not too many Auger electrons to exceed the top ion of this element
        {
          nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] =
              eta_nauger_ionize_over_ionpot_sum[a] / eta_over_ionpot_sum;
          nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] =
              eta_nauger_ionize_sum[a] / eta_sum;
        } else {
          // the following ensures that multiple ionisations can't send you to an ion stage that is not in
          // the model. Send it to the highest ion stage instead
          const int a_replace = topion - ion - 1;

          nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a_replace] +=
              eta_nauger_ionize_over_ionpot_sum[a] / eta_over_ionpot_sum;
          nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a_replace] +=
              eta_nauger_ionize_sum[a] / eta_sum;

          nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 0;
          nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 0.;
        }
      }
    }
  } else {
    const int a = 0;
    nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 1.;
    nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a] = 1.;
  }

  if (matching_nlsubshell_count > 0) {
    double eff_ionpot = X_ion / eta_over_ionpot_sum;
    if (!std::isfinite(eff_ionpot)) {
      eff_ionpot = 0.;
    }
    nt_solution[modelgridindex].eff_ionpot[get_uniqueionindex(element, ion)] = eff_ionpot;
  } else {
    printout("WARNING! No matching subshells in NT impact ionisation cross section data for Z=%d ionstage %d.\n",
             get_atomicnumber(element), get_ionstage(element, ion));
    printout(
        "-> Defaulting to work function approximation and ionisation energy is not accounted for in Spencer-Fano "
        "solution.\n");

    nt_solution[modelgridindex].eff_ionpot[get_uniqueionindex(element, ion)] =
        1. / get_oneoverw(element, ion, modelgridindex);
  }
}

auto get_eff_ionpot(const int modelgridindex, const int element, const int ion) -> float
// get the effective ion potential from the stored value
// a value of 0. should be treated as invalid
{
  return nt_solution[modelgridindex].eff_ionpot[get_uniqueionindex(element, ion)];
  // OR
  // return calculate_eff_ionpot(modelgridindex, element, ion);
}

auto nt_ionization_ratecoeff_sf(const int modelgridindex, const int element, const int ion) -> double
// Kozma & Fransson 1992 equation 13
{
  if (grid::get_numassociatedcells(modelgridindex) <= 0) {
    printout("ERROR: nt_ionization_ratecoeff_sf called on empty cell %d\n", modelgridindex);
    std::abort();
  }

  const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
  if (deposition_rate_density > 0.) {
    return deposition_rate_density / get_nnion_tot(modelgridindex) / get_eff_ionpot(modelgridindex, element, ion);
    // alternatively, if the y vector is still in memory:
    // return calculate_nt_ionization_ratecoeff(modelgridindex, element, ion);
  }
  return 0.;
}

auto get_xs_excitation_vector(gsl_vector *const xs_excitation_vec, const int element, const int ion, const int lower,
                              const int uptransindex, const double statweight_lower, const double epsilon_trans) -> int
// vector of collisional excitation cross sections in cm^2
// epsilon_trans is in erg
// returns the index of the first valid cross section point (en >= epsilon_trans)
// all elements below this index are invalid and should not be used
{
  const double coll_strength = globals::elements[element].ions[ion].levels[lower].uptrans[uptransindex].coll_str;
  if (coll_strength >= 0) {
    // collision strength is available, so use it
    // Li et al. 2012 equation 11
    const double constantfactor = pow(H_ionpot, 2) / statweight_lower * coll_strength * PI * A_naught_squared;

    const int en_startindex = get_energyindex_ev_gteq(epsilon_trans / EV);

    for (int j = 0; j < en_startindex; j++) {
      gsl_vector_set(xs_excitation_vec, j, 0.);
    }

    for (int j = en_startindex; j < SFPTS; j++) {
      const double energy = gsl_vector_get(envec, j) * EV;
      gsl_vector_set(xs_excitation_vec, j, constantfactor * pow(energy, -2));
    }
    return en_startindex;
  }
  const bool forbidden = globals::elements[element].ions[ion].levels[lower].uptrans[uptransindex].forbidden;
  if (!forbidden) {
    const double trans_osc_strength =
        globals::elements[element].ions[ion].levels[lower].uptrans[uptransindex].osc_strength;
    // permitted E1 electric dipole transitions

    // const double g_bar = 0.2;
    const double A = 0.28;
    const double B = 0.15;

    const double prefactor = 45.585750051;  // 8 * pi^2/sqrt(3)
    const double epsilon_trans_ev = epsilon_trans / EV;

    // Eq 4 of Mewe 1972, possibly from Seaton 1962?
    const double constantfactor =
        epsilon_trans_ev * prefactor * A_naught_squared * pow(H_ionpot / epsilon_trans, 2) * trans_osc_strength;

    const int en_startindex = get_energyindex_ev_gteq(epsilon_trans_ev);

    for (int j = 0; j < en_startindex; j++) {
      gsl_vector_set(xs_excitation_vec, j, 0.);
    }

    // U = en / epsilon
    // g_bar = A * log(U) + b
    // xs[j] = constantfactor * g_bar / envec[j]

    for (int j = en_startindex; j < SFPTS; j++) {
      const double logU = gsl_vector_get(logenvec, j) - log(epsilon_trans_ev);
      const double g_bar = A * logU + B;
      gsl_vector_set(xs_excitation_vec, j, constantfactor * g_bar / gsl_vector_get(envec, j));
    }

    return en_startindex;
  }  // gsl_vector_set_zero(xs_excitation_vec);
  return -1;
}

auto calculate_nt_excitation_ratecoeff_perdeposition(const int modelgridindex, const int element, const int ion,
                                                     const int lower, const int uptransindex,
                                                     const double statweight_lower,
                                                     const double epsilon_trans) -> double
// Kozma & Fransson equation 9 divided by level population and epsilon_trans
{
  if (nt_solution[modelgridindex].yfunc == nullptr) {
    printout("ERROR: Call to nt_excitation_ratecoeff with no y vector in memory.");
    std::abort();
  }

  gsl_vector *xs_excitation_vec = gsl_vector_alloc(SFPTS);
  if (get_xs_excitation_vector(xs_excitation_vec, element, ion, lower, uptransindex, statweight_lower, epsilon_trans) >=
      0) {
    double y_dot_crosssection = 0.;
    gsl_vector_view const yvecview = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);
    gsl_blas_ddot(xs_excitation_vec, &yvecview.vector, &y_dot_crosssection);
    gsl_vector_free(xs_excitation_vec);

    y_dot_crosssection *= DELTA_E;

    return y_dot_crosssection / E_init_ev / EV;
  }
  gsl_vector_free(xs_excitation_vec);

  return 0.;
}

auto ion_ntion_energyrate(int modelgridindex, int element, int lowerion) -> double {
  // returns the energy rate [erg/cm3/s] going toward non-thermal ionisation of lowerion
  const double nnlowerion = get_nnion(modelgridindex, element, lowerion);
  double enrate = 0.;
  const auto maxupperion = nt_ionisation_maxupperion(element, lowerion);
  for (int upperion = lowerion + 1; upperion <= maxupperion; upperion++) {
    const double upperionprobfrac =
        nt_ionization_upperion_probability(modelgridindex, element, lowerion, upperion, false);
    // for (int lower = 0; lower < get_nlevels(element, lowerion); lower++)
    // {
    //   const double epsilon_trans = epsilon(element, upperion, 0) - epsilon(element, lowerion, lower);
    //   const double nnlower = get_levelpop(modelgridindex, element, lowerion, lower);
    //   enrate += nnlower * upperionprobfrac * epsilon_trans;
    // }
    const double epsilon_trans = epsilon(element, upperion, 0) - epsilon(element, lowerion, 0);
    enrate += nnlowerion * upperionprobfrac * epsilon_trans;
  }

  const double gamma_nt = nt_ionization_ratecoeff(modelgridindex, element, lowerion);
  return gamma_nt * enrate;
}

auto get_ntion_energyrate(int modelgridindex) -> double {
  // returns the energy rate [erg/s] going toward non-thermal ionisation in a modelgrid cell
  double ratetotal = 0.;
  for (int ielement = 0; ielement < get_nelements(); ielement++) {
    const int nions = get_nions(ielement);
    for (int ilowerion = 0; ilowerion < nions - 1; ilowerion++) {
      ratetotal += ion_ntion_energyrate(modelgridindex, ielement, ilowerion);
    }
  }
  return ratetotal;
}

auto select_nt_ionization(int modelgridindex) -> std::tuple<int, int> {
  const double zrand = rng_uniform();

  // // select based on stored frac_deposition for each ion
  // double frac_deposition_ion_sum = 0.;
  // // zrand is between zero and frac_ionization
  // // keep subtracting off deposition fractions of ionizations transitions until we hit the right one
  // // e.g. if zrand was less than frac_dep_trans1, then use the first transition
  // // e.g. if zrand was between frac_dep_trans1 and frac_dep_trans2 then use the second transition, etc
  // for (int allionindex = 0; allionindex < get_includedions(); allionindex++) {
  //   frac_deposition_ion_sum += nt_solution[modelgridindex].fracdep_ionization_ion[allionindex];
  //   if (frac_deposition_ion_sum >= zrand) {
  //     get_ionfromuniqueionindex(allionindex, element, lowerion);

  //     return;
  //   }
  // }
  // assert_always(false);  // should not reach here

  const double ratetotal = get_ntion_energyrate(modelgridindex);

  // select based on the calcuated energy going to ionisation for each ion
  double ratesum = 0.;
  for (int ielement = 0; ielement < get_nelements(); ielement++) {
    const int nions = get_nions(ielement);
    for (int ilowerion = 0; ilowerion < nions - 1; ilowerion++) {
      ratesum += ion_ntion_energyrate(modelgridindex, ielement, ilowerion);
      if (ratesum >= zrand * ratetotal) {
        return {ielement, ilowerion};
      }
    }
  }
  assert_always(false);
}

void analyse_sf_solution(const int modelgridindex, const int timestep, const bool enable_sfexcitation) {
  const float nne = grid::get_nne(modelgridindex);
  const double nntot = get_nnion_tot(modelgridindex);
  const double nnetot = grid::get_nnetot(modelgridindex);

  double frac_excitation_total = 0.;
  double frac_ionization_total = 0.;

  int excitationindex = 0;  // unique index for every included excitation transition
  nt_solution[modelgridindex].frac_excitations_list.resize(0);
  for (int element = 0; element < get_nelements(); element++) {
    const int Z = get_atomicnumber(element);
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int uniqueionindex = get_uniqueionindex(element, ion);

      const int ionstage = get_ionstage(element, ion);
      const double nnion = get_nnion(modelgridindex, element, ion);

      // if (nnion < minionfraction * get_nnion_tot(modelgridindex)) // skip negligible ions
      if (nnion <= 0.) {  // skip zero-abundance ions
        continue;
      }

      double frac_ionization_ion = 0.;
      double frac_excitation_ion = 0.;
      printout("  Z=%d ionstage %d:\n", Z, ionstage);
      // printout("    nnion: %g\n", nnion);
      printout("    nnion/nntot: %g\n", nnion / nntot);

      calculate_eff_ionpot_auger_rates(modelgridindex, element, ion);

      int matching_nlsubshell_count = 0;
      for (auto &collionrow : colliondata) {
        if (collionrow.Z == Z && collionrow.nelec == Z - ionstage + 1) {
          const double frac_ionization_ion_shell =
              calculate_nt_frac_ionization_shell(modelgridindex, element, ion, collionrow);
          frac_ionization_ion += frac_ionization_ion_shell;
          matching_nlsubshell_count++;
          printout("      shell n %d, l %d, I %5.1f eV: frac_ionization %10.4e", collionrow.n, collionrow.l,
                   collionrow.ionpot_ev, frac_ionization_ion_shell);

          if (NT_MAX_AUGER_ELECTRONS > 0) {
            printout("  prob(n Auger elec):");
            for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
              printout(" %d: %.2f", a, collionrow.prob_num_auger[a]);
            }
          }
          printout("\n");
        }
      }

      // do not ionize the top ion
      if (ion < nions - 1) {
        nt_solution[modelgridindex].fracdep_ionization_ion[uniqueionindex] = frac_ionization_ion;

        frac_ionization_total += frac_ionization_ion;
      } else {
        nt_solution[modelgridindex].fracdep_ionization_ion[uniqueionindex] = 0.;
      }
      printout("    frac_ionization: %g (%d subshells)\n", frac_ionization_ion, matching_nlsubshell_count);

      // excitation from all levels is very SLOW
      const int nlevels_all = get_nlevels(element, ion);
      // So limit the lower levels to improve performance
      int nlevels = (nlevels_all > NTEXCITATION_MAXNLEVELS_LOWER) ? NTEXCITATION_MAXNLEVELS_LOWER : nlevels_all;
      if (!enable_sfexcitation) {
        nlevels = -1;  // disable all excitations
      }
      const bool above_minionfraction = (nnion >= minionfraction * get_nnion_tot(modelgridindex));

      for (int lower = 0; lower < nlevels; lower++) {
        const double statweight_lower = stat_weight(element, ion, lower);
        const int nuptrans = get_nuptrans(element, ion, lower);
        const double nnlevel = get_levelpop(modelgridindex, element, ion, lower);
        const double epsilon_lower = epsilon(element, ion, lower);

        for (int t = 0; t < nuptrans; t++) {
          const int lineindex = globals::elements[element].ions[ion].levels[lower].uptrans[t].lineindex;
          const int upper = globals::linelist[lineindex].upperlevelindex;
          if (upper >= NTEXCITATION_MAXNLEVELS_UPPER) {
            continue;
          }

          const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
          const double nt_frac_excitation_perlevelpop =
              epsilon_trans * calculate_nt_excitation_ratecoeff_perdeposition(modelgridindex, element, ion, lower, t,
                                                                              statweight_lower, epsilon_trans);
          const double frac_excitation_thistrans = nnlevel * nt_frac_excitation_perlevelpop;
          frac_excitation_ion += frac_excitation_thistrans;

          if constexpr (NT_EXCITATION_ON) {
            // the atomic data set was limited for Fe V, which caused the ground multiplet to be massively
            // depleted, and then almost no recombination happened!
            if (above_minionfraction && nt_frac_excitation_perlevelpop > 0 && (Z != 26 || ionstage != 5)) {
              const double ratecoeffperdeposition = nt_frac_excitation_perlevelpop / epsilon_trans;

              assert_always(ratecoeffperdeposition >= 0);
              assert_always(std::isfinite(ratecoeffperdeposition));

              // if (get_coll_str(lineindex) < 0) // if collision strength is not defined, the rate coefficient is
              // unreliable
              //   ratecoeffperdeposition = 0.;

              nt_solution[modelgridindex].frac_excitations_list.push_back({
                  .frac_deposition = frac_excitation_thistrans,
                  .ratecoeffperdeposition = ratecoeffperdeposition,
                  .lineindex = lineindex,
                  .loweruptransindex = t,
              });
              (excitationindex)++;
            }
          }  // NT_EXCITATION_ON
        }  // for t
      }  // for lower

      printout("    frac_excitation: %g\n", frac_excitation_ion);
      if (frac_excitation_ion > 1. || !std::isfinite(frac_excitation_ion)) {
        printout("      WARNING: invalid frac_excitation. Replacing with zero\n");
        frac_excitation_ion = 0.;
      }
      frac_excitation_total += frac_excitation_ion;
      printout("    workfn:       %9.2f eV\n", (1. / get_oneoverw(element, ion, modelgridindex)) / EV);
      printout("    eff_ionpot:   %9.2f eV  (always use valence potential is %s)\n",
               get_eff_ionpot(modelgridindex, element, ion) / EV, (NT_USE_VALENCE_IONPOTENTIAL ? "true" : "false"));

      printout("    workfn approx Gamma:     %9.3e\n", nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion));

      printout("    SF integral Gamma:       %9.3e\n",
               calculate_nt_ionization_ratecoeff(modelgridindex, element, ion, false));

      printout("    SF integral(I=Iv) Gamma: %9.3e  (if always use valence potential)\n",
               calculate_nt_ionization_ratecoeff(modelgridindex, element, ion, true));

      printout("    ARTIS using Gamma:       %9.3e\n", nt_ionization_ratecoeff(modelgridindex, element, ion));

      // the ion values (unlike shell ones) have been collapsed down to ensure that upperion < nions
      if (ion < nions - 1) {
        printout("    probability to ionstage:");
        double prob_sum = 0.;
        for (int upperion = ion + 1; upperion <= nt_ionisation_maxupperion(element, ion); upperion++) {
          const double probability = nt_ionization_upperion_probability(modelgridindex, element, ion, upperion, false);
          prob_sum += probability;
          if (probability > 0.) {
            printout(" %d: %.3f", get_ionstage(element, upperion), probability);
          }
        }
        printout("\n");
        assert_always((fabs(prob_sum - 1.0) <= 1e-2) ||
                      (nt_ionization_ratecoeff_sf(modelgridindex, element, ion) < 1e-20));

        printout("         enfrac to ionstage:");
        double enfrac_sum = 0.;
        for (int upperion = ion + 1; upperion <= nt_ionisation_maxupperion(element, ion); upperion++) {
          const double probability = nt_ionization_upperion_probability(modelgridindex, element, ion, upperion, true);
          enfrac_sum += probability;
          if (probability > 0.) {
            printout(" %d: %.3f", get_ionstage(element, upperion), probability);
          }
        }
        printout("\n");
        assert_always(fabs(enfrac_sum - 1.0) <= 1e-2 ||
                      (nt_ionization_ratecoeff_sf(modelgridindex, element, ion) < 1e-20));
      }
    }
  }

  if constexpr (NT_EXCITATION_ON && (MAX_NT_EXCITATIONS_STORED > 0)) {
    // sort by descending frac_deposition
    std::sort(EXEC_PAR_UNSEQ nt_solution[modelgridindex].frac_excitations_list.begin(),
              nt_solution[modelgridindex].frac_excitations_list.end(),
              [](const auto &a, const auto &b) { return static_cast<bool>(a.frac_deposition > b.frac_deposition); });

    // the excitation list is now sorted by frac_deposition descending
    const double deposition_rate_density = get_deposition_rate_density(modelgridindex);

    if (nt_solution[modelgridindex].frac_excitations_list.size() > MAX_NT_EXCITATIONS_STORED) {
      // truncate the sorted list to save memory
      printout("  Truncating non-thermal excitation list from %zu to %d transitions.\n",
               nt_solution[modelgridindex].frac_excitations_list.size(), MAX_NT_EXCITATIONS_STORED);
      nt_solution[modelgridindex].frac_excitations_list.resize(MAX_NT_EXCITATIONS_STORED);
    }

    printout("[info] mem_usage: non-thermal excitations for cell %d at this timestep occupy %.3f MB\n", modelgridindex,
             nt_solution[modelgridindex].frac_excitations_list.size() *
                 sizeof(nt_solution[modelgridindex].frac_excitations_list[0]) / 1024. / 1024.);

    const auto T_e = grid::get_Te(modelgridindex);
    printout("  Top non-thermal excitation fractions (total excitations = %zu):\n",
             nt_solution[modelgridindex].frac_excitations_list.size());
    const int ntransdisplayed =
        std::min(50, static_cast<int>(nt_solution[modelgridindex].frac_excitations_list.size()));

    for (excitationindex = 0; excitationindex < ntransdisplayed; excitationindex++) {
      const double frac_deposition = nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition;
      if (frac_deposition > 0.) {
        const int lineindex = nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex;
        const TransitionLine *line = &globals::linelist[lineindex];
        const int element = line->elementindex;
        const int ion = line->ionindex;
        const int lower = line->lowerlevelindex;
        const int upper = line->upperlevelindex;
        const int uptransindex = nt_solution[modelgridindex].frac_excitations_list[excitationindex].loweruptransindex;
        const double epsilon_trans = epsilon(element, ion, upper) - epsilon(element, ion, lower);

        const double ratecoeffperdeposition =
            nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition;
        const double ntcollexc_ratecoeff = ratecoeffperdeposition * deposition_rate_density;

        const double t_mid = globals::timesteps[timestep].mid;
        const double radexc_ratecoeff = rad_excitation_ratecoeff(modelgridindex, element, ion, lower, uptransindex,
                                                                 epsilon_trans, lineindex, t_mid);

        const double collexc_ratecoeff = col_excitation_ratecoeff(T_e, nne, element, ion, lower, uptransindex,
                                                                  epsilon_trans, stat_weight(element, ion, lower));

        const double exc_ratecoeff = radexc_ratecoeff + collexc_ratecoeff + ntcollexc_ratecoeff;
        const auto coll_str = globals::elements[element].ions[ion].levels[lower].uptrans[uptransindex].coll_str;

        printout(
            "    frac_deposition %.3e Z=%2d ionstage %d lower %4d upper %4d rad_exc %.1e coll_exc %.1e nt_exc %.1e "
            "nt/tot %.1e collstr %.1e lineindex %d\n",
            frac_deposition, get_atomicnumber(element), get_ionstage(element, ion), lower, upper, radexc_ratecoeff,
            collexc_ratecoeff, ntcollexc_ratecoeff, ntcollexc_ratecoeff / exc_ratecoeff, coll_str, lineindex);
      }
    }

    // sort the excitation list by ascending lineindex for fast lookup with a binary search
    std::sort(EXEC_PAR_UNSEQ nt_solution[modelgridindex].frac_excitations_list.begin(),
              nt_solution[modelgridindex].frac_excitations_list.end(),
              [](const auto &a, const auto &b) { return static_cast<bool>(a.lineindex < b.lineindex); });

  }  // NT_EXCITATION_ON

  // calculate number density of non-thermal electrons
  const double deposition_rate_density_ev = get_deposition_rate_density(modelgridindex) / EV;
  const double yscalefactor = deposition_rate_density_ev / E_init_ev;

  double nne_nt_max = 0.;
  for (int i = 0; i < SFPTS; i++) {
    const double endash = gsl_vector_get(envec, i);
    const double delta_endash = DELTA_E;
    const double oneovervelocity = sqrt(9.10938e-31 / 2 / endash / 1.60218e-19) / 100;  // in sec/cm
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

  // store the solution properties now while the NT spectrum is in memory (in case we free before packet prop)
  nt_solution[modelgridindex].frac_heating = calculate_frac_heating(modelgridindex);

  printout("  frac_heating_tot:    %g\n", nt_solution[modelgridindex].frac_heating);
  printout("  frac_excitation_tot: %g\n", frac_excitation_total);
  printout("  frac_ionization_tot: %g\n", frac_ionization_total);
  const double frac_sum = nt_solution[modelgridindex].frac_heating + frac_excitation_total + frac_ionization_total;
  printout("  frac_sum:            %g (should be close to 1.0)\n", frac_sum);

  nt_solution[modelgridindex].frac_heating = 1. - frac_excitation_total - frac_ionization_total;
  printout("  (replacing calculated frac_heating_tot with %g to make frac_sum = 1.0)\n",
           nt_solution[modelgridindex].frac_heating);

  // const double nnion = get_nnion(modelgridindex, element, ion);
  // double ntexcit_in_a = 0.;
  // for (int level = 0; level < get_nlevels(0, 1); level++)
  // {
  //   ntexcit_in_a += nnion * nt_excitation_ratecoeff(modelgridindex, 0, 1, level, 75);
  // }
  // printout("  total nt excitation rate into level 75: %g\n", ntexcit_in_a);

  // compensate for lost energy by scaling the solution
  // E_init_ev *= frac_sum;
}

void sfmatrix_add_excitation(gsl_matrix *const sfmatrix, const int modelgridindex, const int element, const int ion) {
  // excitation terms
  gsl_vector *vec_xs_excitation_deltae = gsl_vector_alloc(SFPTS);

  const int nlevels_all = get_nlevels(element, ion);
  const int nlevels = (nlevels_all > NTEXCITATION_MAXNLEVELS_LOWER) ? NTEXCITATION_MAXNLEVELS_LOWER : nlevels_all;

  for (int lower = 0; lower < nlevels; lower++) {
    const double statweight_lower = stat_weight(element, ion, lower);
    const double nnlevel = get_levelpop(modelgridindex, element, ion, lower);
    const double epsilon_lower = epsilon(element, ion, lower);
    const int nuptrans = get_nuptrans(element, ion, lower);
    for (int t = 0; t < nuptrans; t++) {
      const int lineindex = globals::elements[element].ions[ion].levels[lower].uptrans[t].lineindex;
      const int upper = globals::linelist[lineindex].upperlevelindex;
      if (upper >= NTEXCITATION_MAXNLEVELS_UPPER) {
        continue;
      }
      const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
      const double epsilon_trans_ev = epsilon_trans / EV;
      if (epsilon_trans_ev < SF_EMIN) {
        continue;
      }

      const int xsstartindex =
          get_xs_excitation_vector(vec_xs_excitation_deltae, element, ion, lower, t, statweight_lower, epsilon_trans);
      if (xsstartindex >= 0) {
        gsl_blas_dscal(DELTA_E, vec_xs_excitation_deltae);

        for (int i = 0; i < SFPTS; i++) {
          const double en = gsl_vector_get(envec, i);
          const int stopindex = get_energyindex_ev_lteq(en + epsilon_trans_ev);

          const int startindex = i > xsstartindex ? i : xsstartindex;
          for (int j = startindex; j < stopindex; j++) {
            *gsl_matrix_ptr(sfmatrix, i, j) += nnlevel * gsl_vector_get(vec_xs_excitation_deltae, j);
          }

          // do the last bit separately because we're not using the full delta_e interval
          const double delta_en = DELTA_E;

          const double delta_en_actual = (en + epsilon_trans_ev - gsl_vector_get(envec, stopindex));

          *gsl_matrix_ptr(sfmatrix, i, stopindex) +=
              nnlevel * gsl_vector_get(vec_xs_excitation_deltae, stopindex) * delta_en_actual / delta_en;
        }
      }
    }
  }
  gsl_vector_free(vec_xs_excitation_deltae);
}

void sfmatrix_add_ionization(gsl_matrix *const sfmatrix, const int Z, const int ionstage, const double nnion)
// add the ionization terms to the Spencer-Fano matrix
{
  gsl_vector *const vec_xs_ionization = gsl_vector_alloc(SFPTS);
  for (auto &collionrow : colliondata) {
    if (collionrow.Z == Z && collionrow.nelec == Z - ionstage + 1) {
      const double ionpot_ev = collionrow.ionpot_ev;
      const double en_auger_ev = collionrow.en_auger_ev;
      // const double n_auger_elec_avg = colliondata[n].n_auger_elec_avg;
      const double J = get_J(Z, ionstage, ionpot_ev);

      assert_always(ionpot_ev >= SF_EMIN);

      // printout("Z=%2d ionstage %d n %d l %d ionpot %g eV\n",
      //          Z, ionstage, colliondata[n].n, colliondata[n].l, ionpot_ev);

      const int xsstartindex = get_xs_ionization_vector(vec_xs_ionization, collionrow);
      // Luke Shingles: the use of min and max on the epsilon limits keeps energies
      // from becoming unphysical. This insight came from reading the
      // CMFGEN Fortran source code (Li, Dessart, Hillier 2012, doi:10.1111/j.1365-2966.2012.21198.x)
      // I had neglected this, so the limits of integration were incorrect. The fix didn't massively affect
      // ionisation rates or spectra, but it was a source of error that led to energy fractions not adding up to 100%
      std::array<double, SFPTS> int_eps_upper = {0};
      std::array<double, SFPTS> prefactors = {0};
      for (int j = xsstartindex; j < SFPTS; j++) {
        const double endash = gsl_vector_get(envec, j);
        const double epsilon_upper = std::min((endash + ionpot_ev) / 2, endash);
        int_eps_upper[j] = atan((epsilon_upper - ionpot_ev) / J);
        prefactors[j] = gsl_vector_get(vec_xs_ionization, j) * nnion / atan((endash - ionpot_ev) / 2 / J);
      }

      for (int i = 0; i < SFPTS; i++) {
        // i is the matrix row index, which corresponds to an energy E at which we are solve from y(E)
        const double en = gsl_vector_get(envec, i);

        // endash ranges from en to SF_EMAX, but skip over the zero-cross section points
        const int jstart = std::max(i, xsstartindex);
        for (int j = jstart; j < SFPTS; j++) {
          // j is the matrix column index which corresponds to the piece of the integral at y(E') where E' >= E and E' =
          // envec(j)
          const double endash = gsl_vector_get(envec, j);

          // J * atan[(epsilon - ionpot_ev) / J] is the indefinite integral of 1/[1 + (epsilon - ionpot_ev)^2/ J^2]
          // in Kozma & Fransson 1992 equation 4

          const double epsilon_lower =
              std::max(endash - en, ionpot_ev);  // and epsilon_upper = (endash + ionpot_ev) / 2;
          const double int_eps_lower = atan((epsilon_lower - ionpot_ev) / J);
          if (int_eps_lower <= int_eps_upper[j]) {
            *gsl_matrix_ptr(sfmatrix, i, j) += prefactors[j] * (int_eps_upper[j] - int_eps_lower) * DELTA_E;
          }
        }

        // below is atan((epsilon_lower - ionpot_ev) / J) where epsilon_lower = en + ionpot_ev;
        const double int_eps_lower2 = atan(en / J);

        // endash ranges from 2 * en + ionpot_ev to SF_EMAX
        if (2 * en + ionpot_ev <= SF_EMAX) {
          const int secondintegralstartindex = std::max(xsstartindex, get_energyindex_ev_lteq(2 * en + ionpot_ev));
          for (int j = secondintegralstartindex; j < SFPTS; j++) {
            // epsilon_lower = en + ionpot_ev;
            // epsilon_upper = (endash + ionpot_ev) / 2;
            if (int_eps_lower2 <= int_eps_upper[j]) {
              *gsl_matrix_ptr(sfmatrix, i, j) -= prefactors[j] * (int_eps_upper[j] - int_eps_lower2) * DELTA_E;
            }
          }
        }
      }

      if constexpr (SF_AUGER_CONTRIBUTION_ON) {
        int augerstopindex = 0;
        if (SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN) {
          // en_auger_ev is (if LJS understands it correctly) averaged to include some probability of zero Auger
          // electrons so we need a boost to get the average energy of Auger electrons given that there are one or more
          const double en_boost = 1 / (1. - collionrow.prob_num_auger[0]);

          augerstopindex = get_energyindex_ev_gteq(en_auger_ev * en_boost);
        } else {
          augerstopindex = get_energyindex_ev_gteq(en_auger_ev);
        }

        for (int i = 0; i < augerstopindex; i++) {
          const double en = gsl_vector_get(envec, i);
          const int jstart = i > xsstartindex ? i : xsstartindex;
          for (int j = jstart; j < SFPTS; j++) {
            const double xs = gsl_vector_get(vec_xs_ionization, j);
            if (SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN) {
              const double en_boost = 1 / (1. - collionrow.prob_num_auger[0]);
              for (int a = 1; a <= NT_MAX_AUGER_ELECTRONS; a++) {
                if (en < (en_auger_ev * en_boost / a)) {
                  *gsl_matrix_ptr(sfmatrix, i, j) -= nnion * xs * collionrow.prob_num_auger[a] * a;
                }
              }
            } else {
              assert_always(en < en_auger_ev);
              // printout("SFAuger E %g < en_auger_ev %g so subtracting %g from element with value %g\n", en,
              // en_auger_ev, nnion * xs, ij_contribution);
              *gsl_matrix_ptr(sfmatrix, i, j) -= nnion * xs;  // * n_auger_elec_avg; // * en_auger_ev???
            }
          }
        }
      }
    }
  }
  gsl_vector_free(vec_xs_ionization);
}

void sfmatrix_solve(const gsl_matrix *sfmatrix, const gsl_vector *rhsvec, gsl_vector *yvec) {
  // WARNING: this assumes sfmatrix is in upper triangular form already!
  const gsl_matrix *sfmatrix_LU = sfmatrix;
  gsl_permutation *p = gsl_permutation_calloc(SFPTS);  // identity permutation

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
  gsl_vector *yvec_best = gsl_vector_alloc(SFPTS);  // solution vector with lowest error
  gsl_vector *gsl_work_vector = gsl_vector_calloc(SFPTS);
  gsl_vector *residual_vector = gsl_vector_alloc(SFPTS);
  int iteration = 0;
  for (iteration = 0; iteration < 10; iteration++) {
    if (iteration > 0) {
      gsl_linalg_LU_refine(sfmatrix, sfmatrix_LU, p, rhsvec, yvec,
                           gsl_work_vector);  // first argument must be original matrix
    }

    // calculate Ax - b = residual
    gsl_vector_memcpy(residual_vector, rhsvec);
    gsl_blas_dgemv(CblasNoTrans, 1.0, sfmatrix, yvec, -1.0, residual_vector);

    // value of the largest absolute residual
    const double error = fabs(gsl_vector_get(residual_vector, gsl_blas_idamax(residual_vector)));

    if (error < error_best || error_best < 0.) {
      gsl_vector_memcpy(yvec_best, yvec);
      error_best = error;
    }
    // printout("Linear algebra solver iteration %d has a maximum residual of %g\n",iteration,error);
  }
  if (error_best >= 0.) {
    if (error_best > 1e-10) {
      printout("  SF solver LU_refine: After %d iterations, best solution vector has a max residual of %g (WARNING)\n",
               iteration, error_best);
    }
    gsl_vector_memcpy(yvec, yvec_best);
  }
  gsl_vector_free(yvec_best);
  gsl_vector_free(gsl_work_vector);
  gsl_vector_free(residual_vector);

  // gsl_matrix_free(sfmatrix_LU); // if this matrix is different to sfmatrix then free it
  gsl_permutation_free(p);

  if (gsl_vector_isnonneg(yvec) == 0) {
    printout("solve_sfmatrix: WARNING: y function goes negative!\n");
  }
}
}  // anonymous namespace

void init(const int my_rank, const int ndo_nonempty) {
  assert_always(nonthermal_initialized == false);
  nonthermal_initialized = true;

  deposition_rate_density = static_cast<double *>(malloc(grid::get_npts_model() * sizeof(double)));
  deposition_rate_density_timestep = static_cast<int *>(malloc(grid::get_npts_model() * sizeof(int)));

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    deposition_rate_density[modelgridindex] = -1.;
    deposition_rate_density_timestep[modelgridindex] = -1;
  }

  if (!NT_ON) {
    return;
  }

  read_binding_energies();

  if (!NT_SOLVE_SPENCERFANO) {
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
  printout("  STORE_NT_SPECTRUM %s\n", STORE_NT_SPECTRUM ? "on" : "off");
  printout("  NT_USE_VALENCE_IONPOTENTIAL %s\n", NT_USE_VALENCE_IONPOTENTIAL ? "on" : "off");
  printout("  NT_MAX_AUGER_ELECTRONS %d\n", NT_MAX_AUGER_ELECTRONS);
  printout("  SF_AUGER_CONTRIBUTION %s\n", SF_AUGER_CONTRIBUTION_ON ? "on" : "off");
  printout("  SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN %s\n", SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN ? "on" : "off");

  if (ndo_nonempty > 0) {
    char filename[MAXFILENAMELENGTH];
    snprintf(filename, MAXFILENAMELENGTH, "nonthermalspec_%.4d.out", my_rank);
    nonthermalfile = fopen_required(filename, "w");
    fprintf(nonthermalfile, "timestep modelgridindex index energy_ev source y\n");
    fflush(nonthermalfile);
  }

  nt_solution = static_cast<nt_solution_struct *>(calloc(grid::get_npts_model(), sizeof(nt_solution_struct)));

  size_t mem_usage_yfunc = 0;
  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    // should make these negative?
    nt_solution[modelgridindex].frac_heating = 0.97;
    nt_solution[modelgridindex].frac_ionization = 0.03;
    nt_solution[modelgridindex].frac_excitation = 0.;

    nt_solution[modelgridindex].nneperion_when_solved = -1.;
    nt_solution[modelgridindex].timestep_last_solved = -1;

    if (grid::get_numassociatedcells(modelgridindex) > 0) {
      nt_solution[modelgridindex].eff_ionpot = static_cast<float *>(calloc(get_includedions(), sizeof(float)));
      nt_solution[modelgridindex].fracdep_ionization_ion =
          static_cast<double *>(calloc(get_includedions(), sizeof(double)));

      nt_solution[modelgridindex].prob_num_auger =
          static_cast<float *>(calloc(get_includedions() * (NT_MAX_AUGER_ELECTRONS + 1), sizeof(float)));
      nt_solution[modelgridindex].ionenfrac_num_auger =
          static_cast<float *>(calloc(get_includedions() * (NT_MAX_AUGER_ELECTRONS + 1), sizeof(float)));

      if (STORE_NT_SPECTRUM) {
        nt_solution[modelgridindex].yfunc = static_cast<double *>(calloc(SFPTS, sizeof(double)));
        assert_always(nt_solution[modelgridindex].yfunc != nullptr);
        mem_usage_yfunc += SFPTS * sizeof(double);
      }

      zero_all_effionpot(modelgridindex);
    } else {
      nt_solution[modelgridindex].eff_ionpot = nullptr;
      nt_solution[modelgridindex].fracdep_ionization_ion = nullptr;

      nt_solution[modelgridindex].prob_num_auger = nullptr;
      nt_solution[modelgridindex].ionenfrac_num_auger = nullptr;

      nt_solution[modelgridindex].yfunc = nullptr;
    }

    nt_solution[modelgridindex].frac_excitations_list.clear();
  }

  if (STORE_NT_SPECTRUM) {
    printout("[info] mem_usage: storing non-thermal spectra for all allocated cells occupies %.3f MB\n",
             mem_usage_yfunc / 1024. / 1024.);
  };

  envec = gsl_vector_calloc(SFPTS);  // energy grid on which solution is sampled
  logenvec = gsl_vector_calloc(SFPTS);
  sourcevec = gsl_vector_calloc(SFPTS);  // energy grid on which solution is sampled

  // const int source_spread_pts = ceil(SFPTS / 20);
  const int source_spread_pts = ceil(SFPTS * 0.03333);  // KF92 OXYGEN TEST
  const double source_spread_en = source_spread_pts * DELTA_E;
  const int sourcelowerindex = SFPTS - source_spread_pts;

  for (int s = 0; s < SFPTS; s++) {
    const double energy_ev = SF_EMIN + s * DELTA_E;

    gsl_vector_set(envec, s, energy_ev);
    gsl_vector_set(logenvec, s, log(energy_ev));

    // spread the source over some energy width
    if (s < sourcelowerindex) {
      gsl_vector_set(sourcevec, s, 0.);
    } else {
      gsl_vector_set(sourcevec, s, 1. / source_spread_en);
    }
  }

  // integrate the source vector to find the assumed injection rate
  gsl_vector *integralvec = gsl_vector_alloc(SFPTS);
  gsl_vector_memcpy(integralvec, sourcevec);
  gsl_vector_scale(integralvec, DELTA_E);
  const double sourceintegral = gsl_blas_dasum(integralvec);  // integral of S(e) dE

  gsl_vector_mul(integralvec, envec);
  E_init_ev = gsl_blas_dasum(integralvec);  // integral of E * S(e) dE
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
  const int nonemptymgi = grid::get_modelcell_nonemptymgi(modelgridindex);
  const double gamma_deposition = globals::dep_estimator_gamma[nonemptymgi] * FOURPI;

  const double tmid = globals::timesteps[timestep].mid;
  const double rho = grid::get_rho(modelgridindex);

  // TODO: calculate thermalisation ratio from the previous timestep either globally (easy) or per cell
  // f = E_dep / E_rad

  // convert from [erg/s/g] to [erg/s/cm3]
  const double positron_deposition =
      rho * decay::get_particle_injection_rate(modelgridindex, tmid, decay::DECAYTYPE_BETAPLUS);

  const double electron_deposition =
      rho * decay::get_particle_injection_rate(modelgridindex, tmid, decay::DECAYTYPE_BETAMINUS);

  const double alpha_deposition =
      rho * decay::get_particle_injection_rate(modelgridindex, tmid, decay::DECAYTYPE_ALPHA);

  deposition_rate_density[modelgridindex] =
      (gamma_deposition + positron_deposition + electron_deposition + alpha_deposition);

  printout(
      "deposition rates [eV/s/cm^3] for timestep %d mgi %d: gamma %8.2e (Monte Carlo), positron %8.2e elec %8.2e alpha "
      "%8.2e (analytic t_mid)\n",
      timestep, modelgridindex, gamma_deposition / EV, positron_deposition / EV, electron_deposition / EV,
      alpha_deposition / EV);

  deposition_rate_density_timestep[modelgridindex] = timestep;
}

auto get_deposition_rate_density(const int modelgridindex) -> double
// should be in erg / s / cm^3
{
  assert_always(deposition_rate_density_timestep[modelgridindex] == globals::timestep);
  assert_always(deposition_rate_density[modelgridindex] >= 0);
  return deposition_rate_density[modelgridindex];
}

void close_file() {
  nonthermal_initialized = false;

  free(deposition_rate_density);
  free(deposition_rate_density_timestep);

  if (!NT_ON || !NT_SOLVE_SPENCERFANO) {
    return;
  }

  if (nonthermalfile != nullptr) {
    fclose(nonthermalfile);
    nonthermalfile = nullptr;
  }
  gsl_vector_free(envec);
  gsl_vector_free(sourcevec);
  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    if (grid::get_numassociatedcells(modelgridindex) > 0) {
      if (STORE_NT_SPECTRUM) {
        free(nt_solution[modelgridindex].yfunc);
      }
      free(nt_solution[modelgridindex].fracdep_ionization_ion);
      free(nt_solution[modelgridindex].eff_ionpot);
      free(nt_solution[modelgridindex].prob_num_auger);
      free(nt_solution[modelgridindex].ionenfrac_num_auger);
    }
  }
  colliondata.clear();
}

auto get_nt_frac_heating(const int modelgridindex) -> float {
  if (!NT_ON) {
    return 1.;
  }
  if (!NT_SOLVE_SPENCERFANO) {
    return 0.97;
  }
  const float frac_heating = nt_solution[modelgridindex].frac_heating;
  // add any debugging checks here?
  return frac_heating;
}

auto nt_ionization_upperion_probability(const int modelgridindex, const int element, const int lowerion,
                                        const int upperion, const bool energyweighted) -> double {
  assert_always(upperion > lowerion);
  assert_always(upperion < get_nions(element));
  assert_always(upperion <= nt_ionisation_maxupperion(element, lowerion));
  if (NT_SOLVE_SPENCERFANO && NT_MAX_AUGER_ELECTRONS > 0) {
    const int numaugerelec = upperion - lowerion - 1;  // number of Auger electrons to go from lowerin to upper ion
    const int uniqueionindex = get_uniqueionindex(element, lowerion);

    if (numaugerelec < NT_MAX_AUGER_ELECTRONS) {
      if (energyweighted) {
        return nt_solution[modelgridindex]
            .ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + numaugerelec];
      }
      return nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + numaugerelec];
    }
    if (numaugerelec == NT_MAX_AUGER_ELECTRONS) {
      double prob_remaining = 1.;
      for (int a = 0; a < NT_MAX_AUGER_ELECTRONS; a++) {
        if (energyweighted) {
          prob_remaining -=
              nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a];
        } else {
          prob_remaining -=
              nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a];
        }
      }
      if (energyweighted) {
        assert_always(fabs(prob_remaining -
                           nt_solution[modelgridindex]
                               .ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + numaugerelec]) <
                      0.001);
      } else {
        if (fabs(prob_remaining - nt_solution[modelgridindex]
                                      .prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + numaugerelec]) >=
            0.001) {
          printout("Auger probabilities issue for cell %d Z=%02d ionstage %d to %d\n", modelgridindex,
                   get_atomicnumber(element), get_ionstage(element, lowerion), get_ionstage(element, upperion));
          for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
            printout("  a %d prob %g\n", a,
                     nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a]);
          }
          std::abort();
        }
      }
      return prob_remaining;
    }
    printout("WARNING: tried to ionise from Z=%02d ionstage %d to %d\n", get_atomicnumber(element),
             get_ionstage(element, lowerion), get_ionstage(element, upperion));
    return 0.;
  }
  return (upperion == lowerion + 1) ? 1.0 : 0.;
}

auto nt_ionisation_maxupperion(const int element, const int lowerion) -> int {
  const int nions = get_nions(element);
  assert_always(lowerion < nions - 1);
  int maxupper = lowerion + 1;

  if (NT_SOLVE_SPENCERFANO) {
    maxupper += NT_MAX_AUGER_ELECTRONS;
  }

  if (maxupper > nions - 1) {
    maxupper = nions - 1;
  }

  return maxupper;
}

auto nt_random_upperion(const int modelgridindex, const int element, const int lowerion,
                        const bool energyweighted) -> int {
  assert_testmodeonly(lowerion < get_nions(element) - 1);
  if (NT_SOLVE_SPENCERFANO && NT_MAX_AUGER_ELECTRONS > 0) {
    while (true) {
      const double zrand = rng_uniform();

      double prob_sum = 0.;
      for (int upperion = lowerion + 1; upperion <= nt_ionisation_maxupperion(element, lowerion); upperion++) {
        prob_sum += nt_ionization_upperion_probability(modelgridindex, element, lowerion, upperion, energyweighted);

        if (zrand <= prob_sum) {
          return upperion;
        }
      }

      printout(
          "ERROR: nt_ionization_upperion_probability did not sum to more than zrand = %lg, prob_sum = %lg (Z=%d "
          "ionstage %d). Retrying with new random number.\n",
          zrand, prob_sum, get_atomicnumber(element), get_ionstage(element, lowerion));
      assert_always(fabs(prob_sum - 1.0) < 1e-3);
    }
  } else {
    return lowerion + 1;
  }
}

auto nt_ionization_ratecoeff(const int modelgridindex, const int element, const int ion) -> double {
  assert_always(NT_ON);
  assert_always(grid::get_numassociatedcells(modelgridindex) > 0);

  if (NT_SOLVE_SPENCERFANO) {
    const double Y_nt = nt_ionization_ratecoeff_sf(modelgridindex, element, ion);
    if (!std::isfinite(Y_nt)) {
      // probably because eff_ionpot = 0 because the solver hasn't been run yet, or no impact ionization cross sections
      // exist
      const double Y_nt_wfapprox = nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion);
      // printout("Warning: Spencer-Fano solver gives non-finite ionization rate (%g) for element %d ionstage %d for
      // cell %d. Using WF approx instead = %g\n",
      //          Y_nt, get_atomicnumber(element), get_ionstage(element, ion), modelgridindex, Y_nt_wfapprox);
      return Y_nt_wfapprox;
    }
    if (Y_nt <= 0) {
      const double Y_nt_wfapprox = nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion);
      if (Y_nt_wfapprox > 0) {
        printout(
            "Warning: Spencer-Fano solver gives negative or zero ionization rate (%g) for element Z=%d ionstage %d "
            "cell %d. Using WF approx instead = %g\n",
            Y_nt, get_atomicnumber(element), get_ionstage(element, ion), modelgridindex, Y_nt_wfapprox);
      }
      return Y_nt_wfapprox;
    }
    return Y_nt;
  }
  return nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion);
}

auto nt_excitation_ratecoeff(const int modelgridindex, const int element, const int ion, const int lowerlevel,
                             const int uptransindex, const double epsilon_trans, const int lineindex) -> double {
  if constexpr (!NT_EXCITATION_ON) {
    return 0.;
  }
  if (lowerlevel >= NTEXCITATION_MAXNLEVELS_LOWER) {
    return 0.;
  }
  const int upperlevel = globals::elements[element].ions[ion].levels[lowerlevel].uptrans[uptransindex].targetlevelindex;
  if (upperlevel >= NTEXCITATION_MAXNLEVELS_UPPER) {
    return 0.;
  }

  if (grid::get_numassociatedcells(modelgridindex) <= 0) {
    printout("ERROR: nt_excitation_ratecoeff called on empty cell %d\n", modelgridindex);
    std::abort();
  }

  // if the NT spectrum is stored, we can calculate any non-thermal excitation rate, even if
  // it didn't make the cut to be kept in the stored excitation list
  if (STORE_NT_SPECTRUM) {
    const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
    const double statweight_lower = stat_weight(element, ion, lowerlevel);

    const double ratecoeffperdeposition = calculate_nt_excitation_ratecoeff_perdeposition(
        modelgridindex, element, ion, lowerlevel, uptransindex, statweight_lower, epsilon_trans);

    return ratecoeffperdeposition * deposition_rate_density;
  }

  // binary search, assuming the excitation list is sorted by lineindex ascending
  auto ntexclist = nt_solution[modelgridindex].frac_excitations_list;
  auto ntexcitation = std::lower_bound(ntexclist.cbegin(), ntexclist.cend(), lineindex,
                                       [](const auto &exc, const int lineindex) { return exc.lineindex < lineindex; });
  if (ntexcitation == ntexclist.end() || ntexcitation->lineindex != lineindex) {
    return 0.;
  }

  const double deposition_rate_density = get_deposition_rate_density(modelgridindex);
  const double ratecoeffperdeposition = ntexcitation->ratecoeffperdeposition;

  return ratecoeffperdeposition * deposition_rate_density;
}

void do_ntlepton(Packet &pkt) {
  atomicadd(nt_energy_deposited, pkt.e_cmf);

  const int modelgridindex = grid::get_cell_modelgridindex(pkt.where);

  // macroatom should not be activated in thick cells
  if (NT_ON && NT_SOLVE_SPENCERFANO && grid::modelgrid[modelgridindex].thick != 1) {
    // here there is some probability to cause ionisation or excitation to a macroatom packet
    // instead of converting directly to k-packet (unless the heating channel is selected)

    double zrand = rng_uniform();
    // zrand is initially between [0, 1), but we will subtract off each
    // component of the deposition fractions
    // until we end and select transition_ij when zrand < dep_frac_transition_ij

    // const double frac_ionization = get_nt_frac_ionization(modelgridindex);
    const double frac_ionization = get_ntion_energyrate(modelgridindex) / get_deposition_rate_density(modelgridindex);
    // printout("frac_ionization compare %g and %g\n", frac_ionization, get_nt_frac_ionization(modelgridindex));
    // const double frac_ionization = 0.;

    if (zrand < frac_ionization) {
      const auto [element, lowerion] = select_nt_ionization(modelgridindex);
      const int upperion = nt_random_upperion(modelgridindex, element, lowerion, true);
      // const int upperion = lowerion + 1;

      pkt.type = TYPE_MA;
      stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_NTCOLLION);
      stats::increment(stats::COUNTER_INTERACTIONS);
      pkt.last_event = 20;
      pkt.trueemissiontype = EMTYPE_NOTSET;
      pkt.trueemissionvelocity = -1;

      stats::increment(stats::COUNTER_NT_STAT_TO_IONIZATION);

      if constexpr (TRACK_ION_STATS) {
        assert_always(upperion < get_nions(element));
        assert_always(lowerion >= 0);
        const double epsilon_trans = epsilon(element, upperion, 0) - epsilon(element, lowerion, 0);
        stats::increment_ion_stats(modelgridindex, element, lowerion, stats::ION_NTION, pkt.e_cmf / epsilon_trans);
        stats::increment_ion_stats(modelgridindex, element, upperion, stats::ION_MACROATOM_ENERGYIN_NTCOLLION,
                                   pkt.e_cmf);
      }

      // printout("NTLEPTON packet in cell %d selected ionization of Z=%d ionstage %d to %d\n",
      //          modelgridindex, get_atomicnumber(element), get_ionstage(element, lowerion), get_ionstage(element,
      //          upperion));
      do_macroatom(pkt, {element, upperion, 0, -99});
      return;
    }

    // const double frac_excitation = NT_EXCITATION_ON ? get_nt_frac_excitation(modelgridindex) : 0;
    const double frac_excitation = 0.;
    if (zrand < (frac_ionization + frac_excitation)) {
      zrand -= frac_ionization;
      // now zrand is between zero and frac_excitation
      // the selection algorithm is the same as for the ionization transitions
      for (const auto &ntexcitation : nt_solution[modelgridindex].frac_excitations_list) {
        const double frac_deposition_exc = ntexcitation.frac_deposition;
        if (zrand < frac_deposition_exc) {
          const int lineindex = ntexcitation.lineindex;
          const int element = globals::linelist[lineindex].elementindex;
          const int ion = globals::linelist[lineindex].ionindex;
          // const int lower = linelist[lineindex].lowerlevelindex;
          const int upper = globals::linelist[lineindex].upperlevelindex;

          pkt.type = TYPE_MA;
          stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_NTCOLLEXC);
          stats::increment(stats::COUNTER_INTERACTIONS);
          pkt.last_event = 21;
          pkt.trueemissiontype = EMTYPE_NOTSET;
          pkt.trueemissionvelocity = -1;

          stats::increment(stats::COUNTER_NT_STAT_TO_EXCITATION);

          // printout("NTLEPTON packet selected in cell %d excitation of Z=%d ionstage %d level %d upperlevel %d\n",
          //          modelgridindex, get_atomicnumber(element), get_ionstage(element, ion), lower, upper);

          do_macroatom(pkt, {element, ion, upper, -99});
          return;
        }
        zrand -= frac_deposition_exc;
      }
      // in case we reached here because the excitation reactions that were stored didn't add up to frac_excitation_ion
      // then just convert it to a kpkt
    }
  }

  pkt.last_event = 22;
  pkt.type = TYPE_KPKT;
  stats::increment(stats::COUNTER_NT_STAT_TO_KPKT);
}

void solve_spencerfano(const int modelgridindex, const int timestep, const int iteration)
// solve the Spencer-Fano equation to get the non-thermal electron flux energy distribution
// based on Equation (2) of Li et al. (2012)
{
  bool skip_solution = false;
  if (grid::get_numassociatedcells(modelgridindex) < 1) {
    printout("Associated_cells < 1 in cell %d at timestep %d. Skipping Spencer-Fano solution.\n", modelgridindex,
             timestep);

    return;
  }
  if (timestep < globals::num_lte_timesteps + 1) {
    printout("Skipping Spencer-Fano solution for first NLTE timestep\n");
    skip_solution = true;
  } else if (get_deposition_rate_density(modelgridindex) / EV < MINDEPRATE) {
    printout(
        "Non-thermal deposition rate of %g eV/cm/s/cm^3 below  MINDEPRATE %g in cell %d at timestep %d. Skipping "
        "Spencer-Fano solution.\n",
        get_deposition_rate_density(modelgridindex) / EV, MINDEPRATE, modelgridindex, timestep);

    skip_solution = true;
  }

  if (skip_solution) {
    nt_solution[modelgridindex].frac_heating = 0.97;
    nt_solution[modelgridindex].frac_ionization = 0.03;
    nt_solution[modelgridindex].frac_excitation = 0.;

    nt_solution[modelgridindex].nneperion_when_solved = -1.;
    nt_solution[modelgridindex].timestep_last_solved = -1;

    nt_solution[modelgridindex].frac_excitations_list.resize(0);

    zero_all_effionpot(modelgridindex);
    return;
  }

  const float nne = grid::get_nne(modelgridindex);  // electrons per cm^3
  // const double nnetot = grid::get_nnetot(modelgridindex);
  const double nne_per_ion = nne / get_nnion_tot(modelgridindex);
  const double nne_per_ion_last = nt_solution[modelgridindex].nneperion_when_solved;
  const double nne_per_ion_fracdiff = fabs((nne_per_ion_last / nne_per_ion) - 1.);
  const int timestep_last_solved = nt_solution[modelgridindex].timestep_last_solved;

  printout(
      "Spencer-Fano solver at timestep %d (last solution was at timestep %d) nne/niontot = %g, at last solution was %g "
      "fracdiff %g\n",
      timestep, timestep_last_solved, nne_per_ion, nne_per_ion_last, nne_per_ion_fracdiff);

  if ((nne_per_ion_fracdiff < NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS) &&
      (timestep - timestep_last_solved <= SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS) &&
      timestep_last_solved > globals::num_lte_timesteps) {
    printout(
        "Keeping Spencer-Fano solution from timestep %d because x_e fracdiff %g < %g and because timestep %d - %d < "
        "%d\n",
        timestep_last_solved, nne_per_ion_fracdiff, NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS, timestep,
        timestep_last_solved, SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS);

    return;
  }
  printout(
      "Setting up Spencer-Fano equation with %d energy points from %g eV to %g eV in cell %d at timestep %d iteration "
      "%d (nne=%g e-/cm^3)\n",
      SFPTS, SF_EMIN, SF_EMAX, modelgridindex, timestep, iteration, nne);

  assert_always(SF_EMIN > 0.);

  nt_solution[modelgridindex].nneperion_when_solved = nne_per_ion;
  nt_solution[modelgridindex].timestep_last_solved = timestep;

  const bool enable_sfexcitation = true;
  const bool enable_sfionization = true;
  // if (timestep <= globals::num_lte_timesteps)
  // {
  //   // for the first run of the solver at the first NLTE timestep (which usually requires many iterations),
  //   // do a fast initial solution but mark it has an invalid nne per ion so it gets replaced at the next timestep
  //   nt_solution[modelgridindex].nneperion_when_solved = -1.;
  //   enable_sfexcitation = false;
  //   enable_sfionization = false;
  //
  //   printout("Doing a fast initial solution without ionization or excitation in the SF equation for the first NLTE
  //   timestep.\n");
  // }
  // if (timestep <= globals::num_lte_timesteps + 2)
  // {
  //   // run the solver in a faster mode for the first couple of NLTE timesteps
  //   // nt_solution[modelgridindex].nneperion_when_solved = -1.;
  //   enable_sfexcitation = false;
  //   // enable_sfionization = false;
  //
  //   printout("Doing a faster solution without excitation in the SF equation for the first couple of NLTE
  //   timesteps.\n");
  // }

  gsl_matrix *const sfmatrix = gsl_matrix_calloc(SFPTS, SFPTS);
  gsl_vector *const rhsvec = gsl_vector_calloc(SFPTS);  // constant term (not dependent on y func) in each equation

  // loss terms and source terms
  for (int i = 0; i < SFPTS; i++) {
    const double en = gsl_vector_get(envec, i);

    *gsl_matrix_ptr(sfmatrix, i, i) += electron_loss_rate(en * EV, nne) / EV;

    double source_integral_to_SF_EMAX{NAN};
    if (i < SFPTS - 1) {
      gsl_vector_const_view source_e_to_SF_EMAX = gsl_vector_const_subvector(sourcevec, i + 1, SFPTS - i - 1);
      source_integral_to_SF_EMAX = gsl_blas_dasum(&source_e_to_SF_EMAX.vector) * DELTA_E;
    } else {
      source_integral_to_SF_EMAX = 0;
    }

    gsl_vector_set(rhsvec, i, source_integral_to_SF_EMAX);
  }
  // gsl_vector_set_all(rhsvec, 1.); // alternative if all electrons are injected at SF_EMAX

  if (enable_sfexcitation || enable_sfionization) {
    for (int element = 0; element < get_nelements(); element++) {
      const int Z = get_atomicnumber(element);
      const int nions = get_nions(element);
      bool first_included_ion_of_element = true;
      for (int ion = 0; ion < nions; ion++) {
        const double nnion = get_nnion(modelgridindex, element, ion);  // hopefully ions per cm^3?

        if (nnion < minionfraction * get_nnion_tot(modelgridindex))  // skip negligible ions
        {
          continue;
        }

        const int ionstage = get_ionstage(element, ion);
        if (first_included_ion_of_element) {
          printout("  including Z=%2d ionstages: ", Z);
          for (int i = 1; i < get_ionstage(element, ion); i++) {
            printout("  ");
          }
          first_included_ion_of_element = false;
        }

        printout("%d ", ionstage);

        if (enable_sfexcitation) {
          sfmatrix_add_excitation(sfmatrix, modelgridindex, element, ion);
        }

        if (enable_sfionization && (ion < nions - 1)) {
          sfmatrix_add_ionization(sfmatrix, Z, ionstage, nnion);
        }
      }
      if (!first_included_ion_of_element) {
        printout("\n");
      }
    }
  }

  // printout("SF matrix | RHS vector:\n");
  // for (int row = 0; row < 10; row++)
  // {
  //   for (int column = 0; column < 10; column++)
  //   {
  //     char str[15];
  //     snprintf(str, 15, "%+.1e ", gsl_matrix_get(sfmatrix, row, column));
  //     printout(str);
  //   }
  //   printout("| ");
  //   char str[15];
  //   snprintf(str, 15, "%+.1e\n", gsl_vector_get(rhsvec, row));
  //   printout(str);
  // }
  // printout("\n");

  if (!STORE_NT_SPECTRUM) {
    nt_solution[modelgridindex].yfunc = static_cast<double *>(calloc(SFPTS, sizeof(double)));
  }

  gsl_vector_view yvecview = gsl_vector_view_array(nt_solution[modelgridindex].yfunc, SFPTS);
  gsl_vector *yvec = &yvecview.vector;

  sfmatrix_solve(sfmatrix, rhsvec, yvec);

  // gsl_matrix_free(sfmatrix_LU); // if sfmatrix_LU is different to sfmatrix

  gsl_matrix_free(sfmatrix);
  gsl_vector_free(rhsvec);

  if (timestep % 10 == 0) {
    nt_write_to_file(modelgridindex, timestep, iteration);
  }

  analyse_sf_solution(modelgridindex, timestep, enable_sfexcitation);

  if (!STORE_NT_SPECTRUM) {
    free(nt_solution[modelgridindex].yfunc);
  }
}

void write_restart_data(FILE *gridsave_file) {
  printout("non-thermal solver, ");

  fprintf(gridsave_file, "%d\n", 24724518);  // special number marking the beginning of NT data
  fprintf(gridsave_file, "%d %la %la\n", SFPTS, SF_EMIN, SF_EMAX);

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    fprintf(gridsave_file, "%d %d %la ", modelgridindex, deposition_rate_density_timestep[modelgridindex],
            deposition_rate_density[modelgridindex]);

    if (NT_ON && NT_SOLVE_SPENCERFANO) {
      check_auger_probabilities(modelgridindex);

      fprintf(gridsave_file, "%a %a %a %a\n", nt_solution[modelgridindex].nneperion_when_solved,
              nt_solution[modelgridindex].frac_heating, nt_solution[modelgridindex].frac_ionization,
              nt_solution[modelgridindex].frac_excitation);

      for (int uniqueionindex = 0; uniqueionindex < get_includedions(); uniqueionindex++) {
        fprintf(gridsave_file, "%la ", nt_solution[modelgridindex].fracdep_ionization_ion[uniqueionindex]);
        fprintf(gridsave_file, "%a ", nt_solution[modelgridindex].eff_ionpot[uniqueionindex]);

        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          fprintf(gridsave_file, "%a %a ",
                  nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a],
                  nt_solution[modelgridindex].ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a]);
        }
      }

      // write NT excitations
      fprintf(gridsave_file, "%d\n", static_cast<int>(nt_solution[modelgridindex].frac_excitations_list.size()));

      for (const auto &excitation : nt_solution[modelgridindex].frac_excitations_list) {
        fprintf(gridsave_file, "%la %la %d\n", excitation.frac_deposition, excitation.ratecoeffperdeposition,
                excitation.lineindex);
      }

      // write non-thermal spectrum
      if (STORE_NT_SPECTRUM) {
        for (int s = 0; s < SFPTS; s++) {
          fprintf(gridsave_file, "%la\n", nt_solution[modelgridindex].yfunc[s]);
        }
      }
    }
  }
}

void read_restart_data(FILE *gridsave_file) {
  printout("Reading restart data for non-thermal solver\n");

  int code_check = 0;
  assert_always(fscanf(gridsave_file, "%d\n", &code_check) == 1);
  if (code_check != 24724518) {
    printout("ERROR: Beginning of non-thermal restart data not found! Found %d instead of 24724518\n", code_check);
    std::abort();
  }

  int sfpts_in = 0;
  double SF_EMIN_in{NAN};
  double SF_EMAX_in{NAN};
  assert_always(fscanf(gridsave_file, "%d %la %la\n", &sfpts_in, &SF_EMIN_in, &SF_EMAX_in) == 3);

  if (sfpts_in != SFPTS || SF_EMIN_in != SF_EMIN || SF_EMAX_in != SF_EMAX) {
    printout("ERROR: gridsave file specifies %d Spencer-Fano samples, SF_EMIN %lg SF_EMAX %lg\n", sfpts_in, SF_EMIN_in,
             SF_EMAX_in);
    printout("ERROR: This simulation has %d Spencer-Fano samples, SF_EMIN %lg SF_EMAX %lg\n", SFPTS, SF_EMIN, SF_EMAX);
    std::abort();
  }

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);

    int mgi_in = 0;
    assert_always(fscanf(gridsave_file, "%d %d %la ", &mgi_in, &deposition_rate_density_timestep[modelgridindex],
                         &deposition_rate_density[modelgridindex]) == 3);

    if (NT_ON && NT_SOLVE_SPENCERFANO) {
      assert_always(fscanf(gridsave_file, "%a %a %a %a\n", &nt_solution[modelgridindex].nneperion_when_solved,
                           &nt_solution[modelgridindex].frac_heating, &nt_solution[modelgridindex].frac_ionization,
                           &nt_solution[modelgridindex].frac_excitation) == 4);

      if (mgi_in != modelgridindex) {
        printout("ERROR: expected data for cell %d but found cell %d\n", modelgridindex, mgi_in);
        std::abort();
      }

      for (int uniqueionindex = 0; uniqueionindex < get_includedions(); uniqueionindex++) {
        assert_always(
            fscanf(gridsave_file, "%la ", &nt_solution[modelgridindex].fracdep_ionization_ion[uniqueionindex]) == 1);
        assert_always(fscanf(gridsave_file, "%a ", &nt_solution[modelgridindex].eff_ionpot[uniqueionindex]) == 1);

        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          assert_always(
              fscanf(gridsave_file, "%a %a ",
                     &nt_solution[modelgridindex].prob_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a],
                     &nt_solution[modelgridindex]
                          .ionenfrac_num_auger[uniqueionindex * (NT_MAX_AUGER_ELECTRONS + 1) + a]) == 2);
        }
      }

      check_auger_probabilities(modelgridindex);

      // read NT excitations
      int frac_excitations_list_size_in = 0;
      assert_always(fscanf(gridsave_file, "%d\n", &frac_excitations_list_size_in) == 1);

      if (static_cast<int>(nt_solution[modelgridindex].frac_excitations_list.size()) != frac_excitations_list_size_in) {
        nt_solution[modelgridindex].frac_excitations_list.resize(frac_excitations_list_size_in);
      }

      for (int excitationindex = 0; excitationindex < frac_excitations_list_size_in; excitationindex++) {
        assert_always(fscanf(gridsave_file, "%la %la %d\n",
                             &nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition,
                             &nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition,
                             &nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex) == 3);
      }

      // read non-thermal spectrum
      if (STORE_NT_SPECTRUM) {
        for (int s = 0; s < SFPTS; s++) {
          assert_always(fscanf(gridsave_file, "%la\n", &nt_solution[modelgridindex].yfunc[s]) == 1);
        }
      }
    }
  }
}

#ifdef MPI_ON
void nt_MPI_Bcast(const int modelgridindex, const int root) {
  if (grid::get_numassociatedcells(modelgridindex) == 0) {
    return;
  }

  // printout("nonthermal_MPI_Bcast cell %d before: ratecoeff(Z=%d ionstage %d): %g, eff_ionpot %g eV\n",
  //          modelgridindex, logged_element_z, logged_ionstage,
  //          nt_ionization_ratecoeff_sf(modelgridindex, logged_element_index, logged_ion_index),
  //          get_eff_ionpot(modelgridindex, logged_element_index, logged_ion_index) / EV);

  MPI_Bcast(&deposition_rate_density_timestep[modelgridindex], 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&deposition_rate_density[modelgridindex], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

  if (NT_ON && NT_SOLVE_SPENCERFANO) {
    assert_always(nonthermal_initialized);
    MPI_Bcast(&nt_solution[modelgridindex].nneperion_when_solved, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nt_solution[modelgridindex].timestep_last_solved, 1, MPI_INT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nt_solution[modelgridindex].frac_heating, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nt_solution[modelgridindex].frac_ionization, 1, MPI_FLOAT, root, MPI_COMM_WORLD);
    MPI_Bcast(&nt_solution[modelgridindex].frac_excitation, 1, MPI_FLOAT, root, MPI_COMM_WORLD);

    MPI_Bcast(nt_solution[modelgridindex].fracdep_ionization_ion, get_includedions(), MPI_DOUBLE, root, MPI_COMM_WORLD);
    MPI_Bcast(nt_solution[modelgridindex].eff_ionpot, get_includedions(), MPI_FLOAT, root, MPI_COMM_WORLD);

    MPI_Bcast(nt_solution[modelgridindex].prob_num_auger, get_includedions() * (NT_MAX_AUGER_ELECTRONS + 1), MPI_FLOAT,
              root, MPI_COMM_WORLD);
    MPI_Bcast(nt_solution[modelgridindex].ionenfrac_num_auger, get_includedions() * (NT_MAX_AUGER_ELECTRONS + 1),
              MPI_FLOAT, root, MPI_COMM_WORLD);

    // communicate NT excitations
    const auto frac_excitations_list_size_old = nt_solution[modelgridindex].frac_excitations_list.size();
    auto frac_excitations_list_size_new = nt_solution[modelgridindex].frac_excitations_list.size();
    MPI_Bcast(&frac_excitations_list_size_new, 1, MPI_INT, root, MPI_COMM_WORLD);

    if (frac_excitations_list_size_new != frac_excitations_list_size_old) {
      nt_solution[modelgridindex].frac_excitations_list.resize(frac_excitations_list_size_new);
    }

    const auto frac_excitations_list_size = nt_solution[modelgridindex].frac_excitations_list.size();
    for (size_t excitationindex = 0; excitationindex < frac_excitations_list_size; excitationindex++) {
      MPI_Bcast(&nt_solution[modelgridindex].frac_excitations_list[excitationindex].frac_deposition, 1, MPI_DOUBLE,
                root, MPI_COMM_WORLD);
      MPI_Bcast(&nt_solution[modelgridindex].frac_excitations_list[excitationindex].ratecoeffperdeposition, 1,
                MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&nt_solution[modelgridindex].frac_excitations_list[excitationindex].lineindex, 1, MPI_INT, root,
                MPI_COMM_WORLD);
    }

    if (STORE_NT_SPECTRUM) {
      assert_always(nt_solution[modelgridindex].yfunc != nullptr);
      MPI_Bcast(nt_solution[modelgridindex].yfunc, SFPTS, MPI_DOUBLE, root, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    check_auger_probabilities(modelgridindex);
  }
}
#endif

void nt_reset_stats() { nt_energy_deposited = 0.; }

void nt_print_stats(const double modelvolume, const double deltat) {
  const double deposition_rate_density_montecarlo = nt_energy_deposited / EV / modelvolume / deltat;

  // deposition rate density for all cells has not been communicated yet - could change this
  // double total_deposition_rate_density = 0.;
  // for (int mgi = 0; mgi < npts_model; mgi++)
  // {
  //   total_deposition_rate_density += get_deposition_rate_density(mgi) / EV;
  // }
  printout("nt_energy_deposited = %g [eV/s/cm^3]\n", deposition_rate_density_montecarlo);
}

}  // namespace nonthermal