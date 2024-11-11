#include "nonthermal.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_vector_double.h>
#include <mpi.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <ios>
#include <numeric>
#include <span>
#include <sstream>
#include <string>
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
#include "thermalbalance.h"

namespace nonthermal {

namespace {

// minimum number fraction of the total population to include in SF solution
constexpr double minionfraction = 1.e-8;

// minimum deposition rate density (eV/s/cm^3) to solve SF equation
constexpr double MINDEPRATE = 0.;

// Bohr radius squared in cm^2
constexpr double A_naught_squared = 2.800285203e-17;

const std::array<std::string, 28> shellnames{"K ", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "N1",
                                             "N2", "N3", "N4", "N5", "N6", "N7", "O1", "O2", "O3", "O4",
                                             "O5", "O6", "O7", "P1", "P2", "P3", "P4", "Q1"};

std::vector<std::vector<double>> elements_electron_binding;
std::vector<std::vector<int>> elements_shells_q;

struct collionrow {
  int Z{-1};
  int ionstage{-1};
  int n = -1;
  int l = -1;
  double ionpot_ev{NAN};
  double A{NAN};
  double B{NAN};
  double C{NAN};
  double D{NAN};
  // track the statistical weight represented by the values below, so they can be updated with new g-weighted averaged
  // values
  double auger_g_accumulated = 0.;

  // probability of 0, 1, ..., NT_MAX_AUGER_ELECTRONS Auger electrons being ejected when the shell is ionised
  std::array<double, NT_MAX_AUGER_ELECTRONS + 1> prob_num_auger{};

  // the average kinetic energy released in Auger electrons after making a hole in this shell
  float en_auger_ev{NAN};
  float n_auger_elec_avg{NAN};

  collionrow() {
    std::ranges::fill(prob_num_auger, 0.);
    prob_num_auger[0] = 1.;
  }
};

std::vector<collionrow> colliondata;

static_assert(SF_EMIN > 0.);
constexpr double DELTA_E = (SF_EMAX - SF_EMIN) / (SFPTS - 1);

// energy grid on which solution is sampled [eV]
constexpr auto engrid(int index) -> double { return SF_EMIN + (index * DELTA_E); }

const auto logengrid = []() {
  std::vector<double> _logengrid(SFPTS);
  for (int i = 0; i < SFPTS; i++) {
    _logengrid[i] = std::log(engrid(i));
  }
  return _logengrid;
}();

// samples of the source function (energy distribution of deposited energy)
constexpr auto sourcevec(const int index) {
  assert_testmodeonly(index >= 0 && index < SFPTS);

  // spread the source over some energy width
  constexpr int source_spread_pts = static_cast<int>(SFPTS * 0.03333) + 1;
  constexpr double source_spread_en = source_spread_pts * DELTA_E;
  constexpr int sourcestartindex = SFPTS - source_spread_pts;

  return (index < sourcestartindex) ? 0. : 1. / source_spread_en;

  // or put all of the source into one point at SF_EMAX
  // return (index < SFPTS - 1) ? 0. : 1. / DELTA_E;
  // so that E_init_ev = SF_EMAX;
};

// the energy injection rate density (integral of E * S(e) dE) in eV/s/cm3 that the Spencer-Fano equation is solved for.
// This is arbitrary and and the solution will be scaled to match the actual energy deposition rate density.
constexpr double E_init_ev = []() {
  double integral = 0.;
  for (int s = 0; s < SFPTS; s++) {
    integral += sourcevec(s) * DELTA_E * engrid(s);
  }
  return integral;
}();

// rhs is the constant term (not dependent on y func) in each equation
constexpr auto rhsvec = []() {
  std::array<double, SFPTS> _rhsvec{};
  double source_integral_to_SF_EMAX = 0.;
  for (int i = SFPTS - 1; i >= 0; i--) {
    _rhsvec[i] = source_integral_to_SF_EMAX * DELTA_E;
    source_integral_to_SF_EMAX += sourcevec(i);
  }
  return _rhsvec;
}();

// Monte Carlo result - compare to analytical expectation
double nt_energy_deposited = 0;

struct NonThermalExcitation {
  double frac_deposition;  // the fraction of the non-thermal deposition energy going to the excitation transition
  double ratecoeffperdeposition;  // the excitation rate coefficient divided by the deposition rate density
  int lineindex;
};

// pointer to either local or node-shared memory excitation list of all cells
std::span<NonThermalExcitation> excitations_list_all_cells{};

// the minimum of MAX_NT_EXCITATIONS_STORED and the number of included excitation transitions in the atomic dataset
int nt_excitations_stored = 0;

struct NonThermalSolutionIon {
  float eff_ionpot{0.};               // these are used to calculate the non-thermal ionization rate
  double fracdep_ionization_ion{0.};  // the fraction of the non-thermal deposition energy going to ionizing each ion

  // probability that one ionisation of this ion will produce n Auger electrons.
  // items sum to 1.0 for a given ion
  std::array<float, NT_MAX_AUGER_ELECTRONS + 1> prob_num_auger{};
  // like prob_num_auger, but energy weighted. items sum to 1.0 for an ion
  std::array<float, NT_MAX_AUGER_ELECTRONS + 1> ionenfrac_num_auger{};
};

std::span<NonThermalSolutionIon> ion_data_all_cells{};

struct NonThermalCellSolution {
  float frac_heating = 1.;     // energy fractions should add up to 1.0 if the solution is good
  float frac_ionization = 0.;  // fraction of deposition energy going to ionization
  float frac_excitation = 0.;  // fraction of deposition energy going to excitation

  int frac_excitations_list_size = 0;

  int timestep_last_solved = -1;     // the quantities above were calculated for this timestep
  float nneperion_when_solved{NAN};  // the nne when the solver was last run
};

std::span<NonThermalCellSolution> nt_solution;

std::span<double> deposition_rate_density_all_cells;

constexpr auto uppertriangular(const int i, const int j) -> int {
  assert_testmodeonly(i >= 0);
  assert_testmodeonly(i < SFPTS);
  // sometimes you might want to get an offset for a row using j = 0 < i, so that j can be added to it.
  // assert_testmodeonly(j >= i);
  assert_testmodeonly(j < SFPTS);
  return (SFPTS * i) - (i * (i + 1) / 2) + j;
}

constexpr void compactify_triangular_matrix(std::vector<double> &matrix) {
  for (int i = 1; i < SFPTS; i++) {
    const int rowoffset = uppertriangular(i, 0);
    for (int j = 0; j < i; j++) {
      assert_always(matrix[(i * SFPTS) + j] == 0.);
    }
    for (int j = i; j < SFPTS; j++) {
      matrix[rowoffset + j] = matrix[(i * SFPTS) + j];
    }
  }
}

constexpr void decompactify_triangular_matrix(std::vector<double> &matrix) {
  for (int i = SFPTS - 1; i > 0; i--) {
    const int rowoffset = uppertriangular(i, 0);
    for (int j = SFPTS - 1; j >= i; j--) {
      matrix[(i * SFPTS) + j] = matrix[rowoffset + j];
    }
    for (int j = i - 1; j >= 0; j--) {
      matrix[(i * SFPTS) + j] = 0.;
    }
  }
}

void read_shell_configs() {
  assert_always(NT_WORKFUNCTION_USE_SHELL_OCCUPANCY_FILE);
  auto shells_file = fstream_required("electron_shell_occupancy.txt", std::ios::in);

  int nshells = 0;      // number of shell in binding energy file
  int n_z_binding = 0;  // number of elements in file

  std::string line;
  assert_always(get_noncommentline(shells_file, line));
  std::istringstream(line) >> nshells >> n_z_binding;
  printout("Reading electron_shell_occupancy.txt with %d elements and %d shells\n", n_z_binding, nshells);

  elements_shells_q.resize(n_z_binding, std::vector<int>(nshells, 0.));

  assert_always(elements_shells_q.size() == elements_electron_binding.size());

  int zminusone = 0;
  while (get_noncommentline(shells_file, line)) {
    std::istringstream ssline(line);

    int z_element = 0;
    assert_always(ssline >> z_element);

    assert_always(elements_shells_q[zminusone].size() == elements_electron_binding[zminusone].size());
    for (int shell = 0; shell < nshells; shell++) {
      int q = 0;
      assert_always(ssline >> q);
      elements_shells_q.at(zminusone).at(shell) = q;
    }
    zminusone++;
  }
}

void read_binding_energies() {
  const bool binding_en_newformat_local = std::filesystem::exists("binding_energies_lotz_tab1and2.txt") ||
                                          std::filesystem::exists("data/binding_energies_lotz_tab1and2.txt");
  bool binding_en_newformat = binding_en_newformat_local;
  // just in case the file system was faulty and the ranks disagree on the existence of the files
  MPI_Allreduce(MPI_IN_PLACE, &binding_en_newformat, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

  int nshells = 0;      // number of shell in binding energy file
  int n_z_binding = 0;  // number of elements in binding energy file

  const auto *filename = binding_en_newformat ? "binding_energies_lotz_tab1and2.txt" : "binding_energies.txt";
  auto binding_energies_file = fstream_required(filename, std::ios::in);

  std::string line;
  assert_always(get_noncommentline(binding_energies_file, line));
  std::istringstream(line) >> nshells >> n_z_binding;
  printout("Reading binding energies file '%s' with %d elements and %d shells\n", filename, n_z_binding, nshells);

  elements_electron_binding.resize(n_z_binding, std::vector<double>(nshells, 0.));

  for (int zminusone = 0; zminusone < n_z_binding; zminusone++) {
    assert_always(get_noncommentline(binding_energies_file, line));
    std::istringstream ssline(line);
    // new file as an atomic number column
    if (binding_en_newformat) {
      int z_element{-1};
      ssline >> z_element;
      assert_always(z_element == (zminusone + 1));
    }
    for (int shell = 0; shell < nshells; shell++) {
      float bindingenergy = 0.;
      assert_always(ssline >> bindingenergy);
      elements_electron_binding.at(zminusone).at(shell) = bindingenergy * EV;
    }
  }

  if constexpr (NT_WORKFUNCTION_USE_SHELL_OCCUPANCY_FILE) {
    if (!binding_en_newformat) {
      printout(
          "NT_WORKFUNCTION_USE_SHELL_OCCUPANCY_FILE is true, but could not find binding_energies_lotz_tab1and2.txt\n");
    }
    assert_always(binding_en_newformat);
    read_shell_configs();
  }
}

[[nodiscard]] auto get_cell_ntexcitations(const int nonemptymgi) {
  return excitations_list_all_cells.subspan(nonemptymgi * nt_excitations_stored,
                                            nt_solution[nonemptymgi].frac_excitations_list_size);
}

[[nodiscard]] auto get_cell_ion_data(const int nonemptymgi) {
  return ion_data_all_cells.subspan(nonemptymgi * get_includedions(), get_includedions());
}

auto get_auger_probability(const int nonemptymgi, const int element, const int ion, const int naugerelec) {
  assert_always(naugerelec <= NT_MAX_AUGER_ELECTRONS);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  return get_cell_ion_data(nonemptymgi)[uniqueionindex].prob_num_auger[naugerelec];
}

auto get_ion_auger_enfrac(const int nonemptymgi, const int element, const int ion, const int naugerelec) {
  assert_always(naugerelec <= NT_MAX_AUGER_ELECTRONS);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  return get_cell_ion_data(nonemptymgi)[uniqueionindex].ionenfrac_num_auger[naugerelec];
}

void check_auger_probabilities(int nonemptymgi) {
  bool problem_found = false;

  for (int element = 0; element < get_nelements(); element++) {
    for (int ion = 0; ion < get_nions(element) - 1; ion++) {
      double prob_sum = 0.;
      double ionenfrac_sum = 0.;
      for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
        prob_sum += get_auger_probability(nonemptymgi, element, ion, a);
        ionenfrac_sum += get_ion_auger_enfrac(nonemptymgi, element, ion, a);
      }

      if (fabs(prob_sum - 1.0) > 0.001) {
        printout("Problem with Auger probabilities for cell %d Z=%d ionstage %d prob_sum %g\n",
                 grid::get_mgi_of_nonemptymgi(nonemptymgi), get_atomicnumber(element), get_ionstage(element, ion),
                 prob_sum);
        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          printout("%d: %g\n", a, get_auger_probability(nonemptymgi, element, ion, a));
        }
        problem_found = true;
      }

      if (fabs(ionenfrac_sum - 1.0) > 0.001) {
        printout("Problem with Auger energy frac sum for cell %d Z=%d ionstage %d ionenfrac_sum %g\n",
                 grid::get_mgi_of_nonemptymgi(nonemptymgi), get_atomicnumber(element), get_ionstage(element, ion),
                 ionenfrac_sum);
        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          printout("%d: %g\n", a, get_ion_auger_enfrac(nonemptymgi, element, ion, a));
        }
        problem_found = true;
      }
    }
  }

  assert_always(!problem_found);
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
        linepos = line + static_cast<ptrdiff_t>(26 + (a * 5));
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
        if (collionrow.Z == Z && collionrow.ionstage == ionstage && collionrow.n == n && collionrow.l == l) {
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

auto get_approx_shell_occupancies(const int nbound, const int ioncharge) {
  assert_always(nbound > 0);
  assert_always(ioncharge >= 0);
  const int Z = nbound + ioncharge;
  std::vector<int> q;
  q.resize(std::max(static_cast<size_t>(10), elements_electron_binding[Z - 1].size()), 0);

  for (int electron_loop = 0; electron_loop < nbound; electron_loop++) {
    if (q[0] < 2) {
      q[0]++;  // K 1s
    } else if (q[1] < 2) {
      q[1]++;  // L1 2s
    } else if (q[2] < 2) {
      q[2]++;  // L2 2p[1/2]
    } else if (q[3] < 4) {
      q[3]++;  // L3 2p[3/2]
    } else if (q[4] < 2) {
      q[4]++;  // M1 3s
    } else if (q[5] < 2) {
      q[5]++;  // M2 3p[1/2]
    } else if (q[6] < 4) {
      q[6]++;  // M3 3p[3/2]
    } else if (ioncharge == 0) {
      if (q[9] < 2) {
        q[9]++;  // N1 4s
      } else if (q[7] < 4) {
        q[7]++;  // M4 3d[3/2]
      } else if (q[8] < 6) {
        q[8]++;  // M5 3d[5/2]
      } else {
        printout("Going beyond the 4s shell in NT calculation. Abort!\n");
        std::abort();
      }
    } else if (ioncharge == 1) {
      if (q[9] < 1) {
        q[9]++;  // N1 4s
      } else if (q[7] < 4) {
        q[7]++;  // M4 3d[3/2]
      } else if (q[8] < 6) {
        q[8]++;  // M5 3d[5/2]
      } else {
        printout("Going beyond the 4s shell in NT calculation. Abort!\n");
        std::abort();
      }
    } else if (ioncharge > 1) {
      if (q[7] < 4) {
        q[7]++;  // M4 3d[3/2]
      } else if (q[8] < 6) {
        q[8]++;  // M5 3d[5/2]
      } else {
        printout("Going beyond the 4s shell in NT calculation. Abort!\n");
        std::abort();
      }
    }
  }
  assert_always(nbound == std::accumulate(q.begin(), q.end(), 0));
  return q;
}

auto get_shell_occupancies(const int nbound, const int ioncharge) {
  assert_always(nbound > 0);
  assert_always(ioncharge >= 0);
  const int Z = nbound + ioncharge;

  if constexpr (!NT_WORKFUNCTION_USE_SHELL_OCCUPANCY_FILE) {
    return get_approx_shell_occupancies(nbound, ioncharge);
  }

  const auto &element_shells_q_neutral = elements_shells_q.at(Z - 1);
  const size_t shellcount = std::min(element_shells_q_neutral.size(), elements_electron_binding[Z - 1].size());
  auto element_shells_q = std::vector<int>(shellcount);

  int electron_count = 0;
  for (size_t shellindex = 0; shellindex < shellcount; shellindex++) {
    const int electronsinshell_neutral = element_shells_q_neutral.at(shellindex);

    int electronsinshell = 0;
    if ((electron_count + electronsinshell_neutral) <= nbound) {
      electronsinshell = electronsinshell_neutral;
    } else {
      electronsinshell = nbound - electron_count;
    }
    assert_always(electronsinshell <= electronsinshell_neutral);
    element_shells_q[shellindex] = electronsinshell;
    electron_count += electronsinshell;
    assert_always(electron_count <= nbound);
  }

  return element_shells_q;
}

auto get_sum_q_over_binding_energy(const int element, const int ion) -> double {
  const int ioncharge = get_ionstage(element, ion) - 1;
  const int nbound = get_atomicnumber(element) - ioncharge;  // number of bound electrons

  if (nbound <= 0) {
    return 0.;
  }

  // get the approximate shell occupancy if we don't have the data file
  const auto shells_q = get_shell_occupancies(nbound, ioncharge);
  const auto &binding_energies = elements_electron_binding.at(get_atomicnumber(element) - 1);

  double total = 0.;
  for (int shellindex = 0; shellindex < std::ssize(shells_q); shellindex++) {
    const int electronsinshell = shells_q[shellindex];

    if (electronsinshell <= 0) {
      continue;
    }
    double enbinding = binding_energies.at(shellindex);
    const double ionpot = globals::elements[element].ions[ion].ionpot;
    if (enbinding <= 0) {
      // if we don't have the shell's binding energy, use the previous one
      enbinding = binding_energies.at(shellindex - 1);
      assert_always(enbinding > 0);
    }
    total += electronsinshell / std::max(ionpot, enbinding);
  }

  return total;
}

void read_collion_data() {
  printout("Reading collisional ionization data from collion.txt...\n");

  FILE *cifile = fopen_required("collion.txt", "r");
  int colliondatacount = 0;
  assert_always(fscanf(cifile, "%d", &colliondatacount) == 1);
  printout("Reading %d collisional transition rows\n", colliondatacount);
  assert_always(colliondatacount >= 0);

  for (int i = 0; i < colliondatacount; i++) {
    collionrow collionrow{};
    int nelec = -1;
    assert_always(fscanf(cifile, "%2d %2d %1d %1d %lg %lg %lg %lg %lg", &collionrow.Z, &nelec, &collionrow.n,
                         &collionrow.l, &collionrow.ionpot_ev, &collionrow.A, &collionrow.B, &collionrow.C,
                         &collionrow.D) == 9);

    assert_always(nelec > 0);
    collionrow.ionstage = collionrow.Z - nelec + 1;

    const int element = get_elementindex(collionrow.Z);
    if (element < 0 || collionrow.ionstage < get_ionstage(element, 0) ||
        collionrow.ionstage > get_ionstage(element, get_nions(element) - 1)) {
      continue;
    }

    std::ranges::fill(collionrow.prob_num_auger, 0.);
    collionrow.prob_num_auger[0] = 1.;

    collionrow.auger_g_accumulated = 0.;
    collionrow.en_auger_ev = 0.;
    collionrow.n_auger_elec_avg = 0.;

    colliondata.push_back(collionrow);
  }
  printout("Stored %zu of %d input shell cross sections\n", colliondata.size(), colliondatacount);
  for (int element = 0; element < get_nelements(); element++) {
    const int Z = get_atomicnumber(element);
    for (int ion = 0; ion < get_nions(element); ion++) {
      const int ionstage = get_ionstage(element, ion);
      const bool any_data_matched = std::ranges::any_of(colliondata, [Z, ionstage](const collionrow &collionrow) {
        return collionrow.Z == Z && collionrow.ionstage == ionstage;
      });
      if (!any_data_matched) {
        const double ionpot_ev = globals::elements[element].ions[ion].ionpot / EV;
        printout("No collisional ionisation data for Z=%d ionstage %d. Using Lotz approximation with ionpot = %g eV\n",
                 Z, ionstage, ionpot_ev);

        const int ioncharge = ionstage - 1;
        const int nbound = Z - ioncharge;  // number of bound electrons
        // get the approximate shell occupancy if we don't have the data file
        auto shells_q = get_shell_occupancies(nbound, ioncharge);
        int electron_count = 0;
        for (int shellindex = 0; shellindex < std::ssize(shells_q); shellindex++) {
          const int electronsinshell = shells_q.at(shellindex);

          electron_count += electronsinshell;

          if (electronsinshell <= 0) {
            continue;
          }
          double enbinding = elements_electron_binding.at(Z - 1).at(shellindex);
          const double ionpot = ionpot_ev * EV;
          if (enbinding <= 0) {
            // if we don't have the shell's binding energy, use the previous one
            enbinding = elements_electron_binding.at(Z - 1).at(shellindex - 1);
            assert_always(enbinding > 0);
          }
          const double p = std::max(ionpot, enbinding);
          collionrow collionrow{};
          collionrow.Z = Z;
          collionrow.ionstage = ionstage;
          collionrow.n = -1;
          collionrow.l = -shellindex;
          collionrow.ionpot_ev = p / EV;
          collionrow.A = -1.;
          collionrow.B = -1.;
          collionrow.C = -1.;
          collionrow.D = -1.;
          std::ranges::fill(collionrow.prob_num_auger, 0.);
          collionrow.prob_num_auger[0] = 1.;
          collionrow.auger_g_accumulated = 0.;
          collionrow.en_auger_ev = 0.;
          collionrow.n_auger_elec_avg = 0.;

          colliondata.push_back(collionrow);
          if (electron_count >= nbound) {
            break;
          }
        }
      }
    }
  }
  colliondata.shrink_to_fit();
  std::ranges::stable_sort(colliondata, [](const collionrow &a, const collionrow &b) {
    return std::tie(a.Z, a.ionstage, a.ionpot_ev, a.n, a.l) < std::tie(b.Z, b.ionstage, b.ionpot_ev, b.n, b.l);
  });

  fclose(cifile);

  if (NT_MAX_AUGER_ELECTRONS > 0) {
    read_auger_data();
  }
}

auto get_possible_nt_excitation_count() -> int {
  // count the number of excitation transitions that pass the MAXNLEVELS_LOWER and MAXNLEVELS_UPPER conditions
  // this count might be higher than the number of stored ratecoeffs due to the MAX_NT_EXCITATIONS_STORED limit
  int ntexcitationcount = 0;
  for (int element = 0; element < get_nelements(); element++) {
    for (int ion = 0; ion < get_nions(element); ion++) {
      const int lower_nlevels = std::min(NTEXCITATION_MAXNLEVELS_LOWER, get_nlevels(element, ion));
      for (int lower = 0; lower < lower_nlevels; lower++) {
        const int nuptrans = get_nuptrans(element, ion, lower);
        for (int t = 0; t < nuptrans; t++) {
          const int upper = get_uptranslist(element, ion, lower)[t].targetlevelindex;
          if (upper < NTEXCITATION_MAXNLEVELS_UPPER) {
            ntexcitationcount++;
          }
        }
      }
    }
  }
  return ntexcitationcount;
}

void zero_all_effionpot(const int nonemptymgi) {
  for (int uniqueionindex = 0; uniqueionindex < get_includedions(); uniqueionindex++) {
    auto &ion_data = get_cell_ion_data(nonemptymgi)[uniqueionindex];
    ion_data.eff_ionpot = 0.;

    std::ranges::fill(get_cell_ion_data(nonemptymgi)[uniqueionindex].prob_num_auger, 0.);
    std::ranges::fill(get_cell_ion_data(nonemptymgi)[uniqueionindex].ionenfrac_num_auger, 0.);
    ion_data.prob_num_auger[0] = 1.;
    ion_data.ionenfrac_num_auger[0] = 1.;

    const auto [element, ion] = get_ionfromuniqueionindex(uniqueionindex);
    assert_always(fabs(get_auger_probability(nonemptymgi, element, ion, 0) - 1.0) < 1e-3);
    assert_always(fabs(get_ion_auger_enfrac(nonemptymgi, element, ion, 0) - 1.0) < 1e-3);
  }
  check_auger_probabilities(nonemptymgi);
}

[[nodiscard]] constexpr auto get_energyindex_ev_lteq(const double energy_ev) -> int
// finds the highest energy point <= energy_ev
{
  const int index = std::floor((energy_ev - SF_EMIN) / DELTA_E);

  return std::clamp(index, 0, SFPTS - 1);
}

[[nodiscard]] constexpr auto get_energyindex_ev_gteq(const double energy_ev) -> int
// finds the highest energy point <= energy_ev
{
  const int index = std::ceil((energy_ev - SF_EMIN) / DELTA_E);

  return std::clamp(index, 0, SFPTS - 1);
}

// interpolate the y flux values to get the value at a given energy
// y has units of particles / cm2 / s / eV
[[nodiscard]] constexpr auto get_y(const std::array<double, SFPTS> &yfunc, const double energy_ev) -> double {
  if (energy_ev <= 0) {
    return 0.;
  }

  const int index = static_cast<int>((energy_ev - SF_EMIN) / DELTA_E);

  // assert_always(index > 0);
  if (index < 0) {
    // return 0.;
    assert_always(std::isfinite(yfunc[0]));
    return yfunc[0];
  }
  if (index >= SFPTS - 1) {
    return 0.;
  }
  const double enbelow = engrid(index);
  const double enabove = engrid(index + 1);
  const double ybelow = yfunc[index];
  const double yabove = yfunc[index + 1];
  const double x = (energy_ev - enbelow) / (enabove - enbelow);
  return ((1 - x) * ybelow) + (x * yabove);

  // or return the nearest neighbour
  // return yfunc[index];
}

auto xs_ionization_lotz(const double en_erg, const collionrow &colliondata_ion) -> double {
  const double ionpot_ev = colliondata_ion.ionpot_ev;
  if (en_erg < (ionpot_ev * EV)) {
    return 0.;
  }
  // const double gamma = (en_erg / (ME * std::pow(CLIGHT, 2))) + 1;
  // const double beta = std::sqrt(1.0 - (1.0 / (std::pow(gamma, 2))));
  const double beta = std::sqrt(2 * en_erg / ME) / CLIGHT;

  const int ioncharge = colliondata_ion.ionstage - 1;
  const int nbound = colliondata_ion.Z - ioncharge;  // number of bound electrons

  if (nbound <= 0) {
    return 0.;
  }

  const int shellindex = -colliondata_ion.l;
  const int electronsinshell = get_shell_occupancies(nbound, ioncharge)[shellindex];

  const double p = colliondata_ion.ionpot_ev * EV;

  if (en_erg > p) {
    const double part_sigma_shell = (electronsinshell / p *
                                     (std::log(std::pow(beta, 2) * ME * std::pow(CLIGHT, 2) / 2.0 / p) -
                                      std::log10(1 - std::pow(beta, 2)) - std::pow(beta, 2)));
    if (part_sigma_shell > 0.) {
      constexpr double Aconst = 1.33e-14 * EV * EV;
      const double sigma = 2 * Aconst / std::pow(beta, 2) / ME / std::pow(CLIGHT, 2) * part_sigma_shell;
      assert_always(sigma >= 0);
      return sigma;
    }
  }

  return 0.;
}

auto get_xs_ionization_vector_lotz(std::array<double, SFPTS> &xs_vec, const collionrow &colliondata_ion) -> int {
  const double ionpot_ev = colliondata_ion.ionpot_ev;
  const int startindex = get_energyindex_ev_gteq(ionpot_ev);

  std::fill_n(xs_vec.begin(), startindex, 0.);

  for (int i = startindex; i < SFPTS; i++) {
    xs_vec[i] = xs_ionization_lotz(engrid(i) * EV, colliondata_ion);
  }

  return startindex;
}

// xs_vec will be set with impact ionization cross sections [cm2] for E > ionpot_ev (and zeros below this energy)
// returns the index of the first energy point >= ionpot_ev
auto get_xs_ionization_vector(std::array<double, SFPTS> &xs_vec, const collionrow &colliondata_ion) -> int {
  const double A = colliondata_ion.A;
  if (A < 0) {
    return get_xs_ionization_vector_lotz(xs_vec, colliondata_ion);
  }

  const double ionpot_ev = colliondata_ion.ionpot_ev;
  const int startindex = get_energyindex_ev_gteq(ionpot_ev);

  std::fill_n(xs_vec.begin(), startindex, 0.);

  const double B = colliondata_ion.B;
  const double C = colliondata_ion.C;
  const double D = colliondata_ion.D;

  for (int i = startindex; i < SFPTS; i++) {
    const double u = engrid(i) / ionpot_ev;
    const double xs_ioniz = 1e-14 *
                            (A * (1 - 1 / u) + B * std::pow((1 - (1 / u)), 2) + C * std::log(u) + D * std::log(u) / u) /
                            (u * std::pow(ionpot_ev, 2));
    xs_vec[i] = xs_ioniz;
  }

  return startindex;
}

// distribution of secondary electron energies for primary electron with energy e_p
// Opal, Peterson, & Beaty (1971)
[[nodiscard]] constexpr auto Psecondary(const double e_p, const double epsilon, const double I, const double J)
    -> double {
  const double e_s = epsilon - I;

  if (e_p <= I || e_s < 0.) {
    return 0.;
  }
  assert_testmodeonly(J > 0);
  assert_testmodeonly(e_p >= I);
  assert_testmodeonly(e_s >= 0);
  assert_testmodeonly(std::isfinite(std::atan((e_p - I) / 2 / J)));
  return 1 / (J * std::atan((e_p - I) / 2 / J) * (1 + std::pow(e_s / J, 2)));
}

[[nodiscard]] constexpr auto get_J(const int Z, const int ionstage, const double ionpot_ev) -> double {
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

// collisional excitation cross section in cm^2
// energies are in erg
constexpr auto xs_excitation(const int element, const int ion, const int lower, const int uptransindex,
                             const double epsilon_trans, const double lowerstatweight, const double energy) -> double {
  if (energy < epsilon_trans) {
    return 0.;
  }

  const auto &uptrans = get_uptranslist(element, ion, lower)[uptransindex];
  if (uptrans.coll_str >= 0) {
    // collision strength is available, so use it
    // Li et al. 2012 equation 11
    return std::pow(H_ionpot / energy, 2) / lowerstatweight * uptrans.coll_str * PI * A_naught_squared;
  }
  if (!uptrans.forbidden) {
    // permitted E1 electric dipole transitions
    const double U = energy / epsilon_trans;

    // constexpr double g_bar = 0.2;
    constexpr double A = 0.28;
    constexpr double B = 0.15;
    const double g_bar = (A * std::log(U)) + B;

    constexpr double prefactor = 45.585750051;  // 8 * pi^2/sqrt(3)
    // Eq 4 of Mewe 1972, possibly from Seaton 1962?
    return prefactor * A_naught_squared * std::pow(H_ionpot / epsilon_trans, 2) * uptrans.osc_strength * g_bar / U;
  }
  return 0.;
}

// -dE / dx for fast electrons
// energy is in ergs
// nne is the thermal electron density [cm^-3]
// return value has units of erg/cm
constexpr auto electron_loss_rate(const double energy, const double nne) -> double {
  if (energy <= 0.) {
    return 0;
  }

  // normally set to 1.0, but Shingles et al. (2021) boosted this to increase heating
  constexpr double boostfactor = 1.;

  const double omegap = std::sqrt(4 * PI * nne * std::pow(QE, 2) / ME);
  const double zetae = H * omegap / 2 / PI;
  if (energy > 14 * EV) {
    return boostfactor * nne * 2 * PI * std::pow(QE, 4) / energy * std::log(2 * energy / zetae);
  }
  const double v = std::sqrt(2 * energy / ME);
  return boostfactor * nne * 2 * PI * std::pow(QE, 4) / energy *
         std::log(ME * std::pow(v, 3) / (EULERGAMMA * std::pow(QE, 2) * omegap));
}

// impact ionization cross section in cm^2
// energy and ionization_potential should be in eV
// fitting formula of Younger 1981
// called Q_i(E) in KF92 equation 7
constexpr auto xs_impactionization(const double energy_ev, const collionrow &colliondata_ion) -> double {
  const double ionpot_ev = colliondata_ion.ionpot_ev;
  const double u = energy_ev / ionpot_ev;

  if (u <= 1.) {
    return 0;
  }
  const double A = colliondata_ion.A;
  if (A < 0) {
    return xs_ionization_lotz(energy_ev / EV, colliondata_ion);
  }

  const double B = colliondata_ion.B;
  const double C = colliondata_ion.C;
  const double D = colliondata_ion.D;

  return 1e-14 * (A * (1 - 1 / u) + B * std::pow((1 - (1 / u)), 2) + C * std::log(u) + D * std::log(u) / u) /
         (u * std::pow(ionpot_ev, 2));
}

// Kozma & Fransson equation 6.
// Something related to a number of electrons, needed to calculate the heating fraction in equation 3
// not valid for energy > SF_EMIN
auto N_e(const int nonemptymgi, const double energy, const std::array<double, SFPTS> &yfunc) -> double {
  const double energy_ev = energy / EV;
  const double tot_nion = get_nnion_tot(nonemptymgi);
  double N_e = 0.;

  for (int element = 0; element < get_nelements(); element++) {
    const int Z = get_atomicnumber(element);
    const int nions = get_nions(element);

    for (int ion = 0; ion < nions; ion++) {
      double N_e_ion = 0.;
      const int ionstage = get_ionstage(element, ion);
      const double nnion = get_nnion(nonemptymgi, element, ion);

      if (nnion < minionfraction * tot_nion) {  // skip negligible ions
        continue;
      }

      // excitation terms

      const int nlevels_all = get_nlevels(element, ion);
      const int nlevels = std::min(NTEXCITATION_MAXNLEVELS_LOWER, nlevels_all);

      for (int lower = 0; lower < nlevels; lower++) {
        const int nuptrans = get_nuptrans(element, ion, lower);
        const auto *const uptranslist = get_uptranslist(element, ion, lower);
        const double nnlevel = get_levelpop(nonemptymgi, element, ion, lower);
        const double epsilon_lower = epsilon(element, ion, lower);
        const auto statweight_lower = stat_weight(element, ion, lower);
        for (int t = 0; t < nuptrans; t++) {
          const int upper = uptranslist[t].targetlevelindex;
          if (upper >= NTEXCITATION_MAXNLEVELS_UPPER) {
            continue;
          }
          const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
          const double epsilon_trans_ev = epsilon_trans / EV;
          N_e_ion += (nnlevel / nnion) * get_y(yfunc, energy_ev + epsilon_trans_ev) *
                     xs_excitation(element, ion, lower, t, epsilon_trans, statweight_lower, energy + epsilon_trans);
        }
      }

      // ionization terms
      for (const auto &collionrow : colliondata) {
        if (collionrow.Z == Z && collionrow.ionstage == ionstage) {
          const double ionpot_ev = collionrow.ionpot_ev;
          const double J = get_J(Z, ionstage, ionpot_ev);
          const double lambda = std::min(SF_EMAX - energy_ev, energy_ev + ionpot_ev);

          const int integral1startindex = get_energyindex_ev_lteq(ionpot_ev);
          const int integral1stopindex = get_energyindex_ev_lteq(lambda);

          // integral from ionpot up to lambda
          for (int i = integral1startindex; i <= integral1stopindex; i++) {
            const double endash = engrid(i);

            N_e_ion += get_y(yfunc, energy_ev + endash) * xs_impactionization(energy_ev + endash, collionrow) *
                       Psecondary(energy_ev + endash, endash, ionpot_ev, J) * DELTA_E;
          }

          // integral from 2E + I up to E_max
          const int integral2startindex = get_energyindex_ev_lteq((2 * energy_ev) + ionpot_ev);
          for (int i = integral2startindex; i < SFPTS; i++) {
            const double endash = engrid(i);
            N_e_ion += yfunc[i] * xs_impactionization(endash, collionrow) *
                       Psecondary(endash, energy_ev + ionpot_ev, ionpot_ev, J) * DELTA_E;
          }
        }
      }

      N_e += nnion * N_e_ion;
    }
  }

  // source term, should be zero at the low end anyway
  N_e += sourcevec(get_energyindex_ev_lteq(energy_ev));

  assert_always(std::isfinite(N_e));
  return N_e;
}

// fraction of deposited energy that goes into heating the thermal electrons
// Kozma & Fransson equation 3
auto calculate_frac_heating(const int nonemptymgi, const std::array<double, SFPTS> &yfunc) -> float {
  // frac_heating multiplied by E_init, which will be divided out at the end
  double frac_heating_Einit = 0.;
  const float nne = grid::get_nne(nonemptymgi);

  for (int i = 0; i < SFPTS; i++) {
    const double endash = engrid(i);

    // first term
    frac_heating_Einit += yfunc[i] * (electron_loss_rate(endash * EV, nne) / EV) * DELTA_E;
  }

  // second term
  frac_heating_Einit += SF_EMIN * get_y(yfunc, SF_EMIN) * (electron_loss_rate(SF_EMIN * EV, nne) / EV);

  double N_e_contrib = 0.;
  // third term (integral from zero to SF_EMIN)
  constexpr int nsteps = (static_cast<int>(SF_EMIN / DELTA_E) + 1) * 10;
  static_assert(nsteps > 0);
  constexpr double delta_endash = SF_EMIN / nsteps;
  for (int j = 1; j < nsteps; j++) {
    const double endash = delta_endash * j;
    N_e_contrib += N_e(nonemptymgi, endash * EV, yfunc) * endash * delta_endash;
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

// fraction of deposited energy that goes into ionization
auto get_nt_frac_ionization(const int modelgridindex) -> float {
  if (!NT_ON) {
    return 0.;
  }
  if (!NT_SOLVE_SPENCERFANO) {
    return 0.03;
  }
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);

  const float frac_ionization = nt_solution[nonemptymgi].frac_ionization;

  if (frac_ionization < 0 || !std::isfinite(frac_ionization)) {
    printout("ERROR: get_nt_frac_ionization called with no valid solution stored for cell %d. frac_ionization = %g\n",
             modelgridindex, frac_ionization);
    std::abort();
  }

  return frac_ionization;
}

// fraction of deposited energy that goes into collisional excitation
auto get_nt_frac_excitation(const int modelgridindex) -> float {
  if (!NT_ON || !NT_SOLVE_SPENCERFANO) {
    return 0.;
  }

  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  const float frac_excitation = nt_solution[nonemptymgi].frac_excitation;

  if (frac_excitation < 0 || !std::isfinite(frac_excitation)) {
    printout("ERROR: get_nt_frac_excitation called with no valid solution stored for cell %d. frac_excitation = %g\n",
             modelgridindex, frac_excitation);
    std::abort();
  }

  return frac_excitation;
}

// compute the work per ion pair for doing the NT ionization calculation.
// Makes use of EXTREMELY SIMPLE approximations - high energy limits only (can be used as an alternative to the
// Spencer-Fano solver)
auto get_oneoverw(const int element, const int ion, const int modelgridindex) -> double {
  // Work in terms of 1/W since this is actually what we want. It is given by sigma/(Latom + Lelec).
  // We are going to start by taking all the high energy limits and ignoring Lelec, so that the
  // denominator is extremely simplified. Need to get the mean Z value.

  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  double Zbar = 0.;  // mass-weighted average atomic number
  for (int ielement = 0; ielement < get_nelements(); ielement++) {
    Zbar += grid::get_elem_abundance(nonemptymgi, ielement) * get_atomicnumber(ielement);
  }

  const double binding = get_sum_q_over_binding_energy(element, ion);
  constexpr double Aconst = 1.33e-14 * EV * EV;
  const double oneoverW = Aconst * binding / Zbar / (2 * PI * std::pow(QE, 4));

  return oneoverW;
}

// the fraction of deposited energy that goes into ionising electrons in a particular shell
auto calculate_nt_frac_ionization_shell(const int nonemptymgi, const int element, const int ion,
                                        const collionrow &collionrow, const std::array<double, SFPTS> &yfunc)
    -> double {
  const double nnion = get_nnion(nonemptymgi, element, ion);
  const double ionpot_ev = collionrow.ionpot_ev;

  std::array<double, SFPTS> cross_section_vec{};
  get_xs_ionization_vector(cross_section_vec, collionrow);

  const double y_dot_crosssection_de = cblas_ddot(SFPTS, yfunc.data(), 1, cross_section_vec.data(), 1) * DELTA_E;

  return nnion * ionpot_ev * y_dot_crosssection_de / E_init_ev;
}

auto nt_ionization_ratecoeff_wfapprox(const int modelgridindex, const int element, const int ion) -> double
// non-thermal ionization rate coefficient (multiply by population to get rate)
{
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  const double deposition_rate_density = get_deposition_rate_density(nonemptymgi);
  // to get the non-thermal ionization rate we need to divide the energy deposited
  // per unit volume per unit time in the grid cell (sum of terms above)
  // by the total ion number density and the "work per ion pair"
  return deposition_rate_density / get_nnion_tot(nonemptymgi) * get_oneoverw(element, ion, modelgridindex);
}

// Integrate the ionization cross section over the electron degradation function to get the ionization rate
// coefficient i.e. multiply this by ion population to get a rate of ionizations per second Do not call during packet
// propagation, as the y vector may not be in memory! IMPORTANT: we are dividing by the shell potential, not the
// valence potential here! To change this set assumeshellpotentialisvalence to true
auto calculate_nt_ionization_ratecoeff(const int modelgridindex, const int element, const int ion,
                                       const bool assumeshellpotentialisvalence, const std::array<double, SFPTS> &yfunc)
    -> double {
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  std::array<double, SFPTS> cross_section_vec{};
  auto gsl_cross_section_vec = gsl_vector_view_array(cross_section_vec.data(), SFPTS).vector;
  std::array<double, SFPTS> cross_section_vec_allshells{};
  auto gsl_cross_section_vec_allshells = gsl_vector_view_array(cross_section_vec_allshells.data(), SFPTS).vector;

  const int Z = get_atomicnumber(element);
  const int ionstage = get_ionstage(element, ion);
  double ionpot_valence = -1;

  for (const auto &collionrow : colliondata) {
    if (collionrow.Z == Z && collionrow.ionstage == ionstage) {
      get_xs_ionization_vector(cross_section_vec, collionrow);

      if (assumeshellpotentialisvalence) {
        const double ionpot_shell = collionrow.ionpot_ev * EV;
        if (ionpot_valence < 0) {
          ionpot_valence = ionpot_shell;
        }

        // ensure that the first shell really was the valence shell (we assumed ascending energy order)
        assert_always(ionpot_shell >= ionpot_valence);

        // boost the ionization rate by assuming shell vacancy energy is used to eject valence electrons
        gsl_vector_scale(&gsl_cross_section_vec, ionpot_shell / ionpot_valence);
      }

      gsl_vector_add(&gsl_cross_section_vec_allshells, &gsl_cross_section_vec);
    }
  }

  const double y_xs_de = cblas_ddot(SFPTS, yfunc.data(), 1, cross_section_vec_allshells.data(), 1) * DELTA_E;

  const double deposition_rate_density_ev = get_deposition_rate_density(nonemptymgi) / EV;
  const double yscalefactor = deposition_rate_density_ev / E_init_ev;

  return yscalefactor * y_xs_de;
}

void calculate_eff_ionpot_auger_rates(const int nonemptymgi, const int element, const int ion,
                                      const std::array<double, SFPTS> &yfunc)
// Kozma & Fransson 1992 equation 12, except modified to be a sum over all shells of an ion
// the result is in ergs
{
  const int Z = get_atomicnumber(element);
  const int ionstage = get_ionstage(element, ion);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  const double nnion = get_nnion(nonemptymgi, element, ion);  // ions/cm^3
  const double tot_nion = get_nnion_tot(nonemptymgi);
  const double X_ion = nnion / tot_nion;  // molar fraction of this ion

  // The ionization rates of all shells of an ion add to make the ion's total ionization rate,
  // i.e., Gamma_ion = Gamma_shell_a + Gamma_shell_b + ...
  // And since the ionization rate is inversely proportional to the effective ion potential,
  // we solve:
  // (eta_ion / ionpot_ion) = (eta_shell_a / ionpot_shell_a) + (eta_shell_b / ionpot_shell_b) + ...
  // where eta is the fraction of the deposition energy going into ionization of the ion or shell

  std::array<double, NT_MAX_AUGER_ELECTRONS + 1> eta_nauger_ionize_over_ionpot_sum{};
  std::array<double, NT_MAX_AUGER_ELECTRONS + 1> eta_nauger_ionize_sum{};
  const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  std::ranges::fill(get_cell_ion_data(nonemptymgi)[uniqueionindex].prob_num_auger, 0.);
  std::ranges::fill(get_cell_ion_data(nonemptymgi)[uniqueionindex].ionenfrac_num_auger, 0.);

  double eta_over_ionpot_sum = 0.;
  double eta_sum = 0.;
  double ionpot_valence = -1;
  int matching_nlsubshell_count = 0;
  for (const auto &collionrow : colliondata) {
    if (collionrow.Z == Z && collionrow.ionstage == ionstage) {
      matching_nlsubshell_count++;
      const double frac_ionization_shell =
          calculate_nt_frac_ionization_shell(nonemptymgi, element, ion, collionrow, yfunc);
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

      for (ptrdiff_t a = 0; a < std::ssize(eta_nauger_ionize_sum); a++) {
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
          get_cell_ion_data(nonemptymgi)[uniqueionindex].prob_num_auger[a] =
              eta_nauger_ionize_over_ionpot_sum[a] / eta_over_ionpot_sum;
          get_cell_ion_data(nonemptymgi)[uniqueionindex].ionenfrac_num_auger[a] = eta_nauger_ionize_sum[a] / eta_sum;
        } else {
          // the following ensures that multiple ionisations can't send you to an ion stage that is not in
          // the model. Send it to the highest ion stage instead
          const int a_replace = topion - ion - 1;

          get_cell_ion_data(nonemptymgi)[uniqueionindex].prob_num_auger.at(a_replace) +=
              eta_nauger_ionize_over_ionpot_sum[a] / eta_over_ionpot_sum;
          get_cell_ion_data(nonemptymgi)[uniqueionindex].ionenfrac_num_auger.at(a_replace) +=
              eta_nauger_ionize_sum[a] / eta_sum;

          get_cell_ion_data(nonemptymgi)[uniqueionindex].prob_num_auger[a] = 0;
          get_cell_ion_data(nonemptymgi)[uniqueionindex].ionenfrac_num_auger[a] = 0.;
        }
      }
    }
  } else {
    const int a = 0;
    get_cell_ion_data(nonemptymgi)[uniqueionindex].prob_num_auger[a] = 1.;
    get_cell_ion_data(nonemptymgi)[uniqueionindex].ionenfrac_num_auger[a] = 1.;
  }

  if (matching_nlsubshell_count > 0) {
    double eff_ionpot = X_ion / eta_over_ionpot_sum;
    if (!std::isfinite(eff_ionpot)) {
      eff_ionpot = 0.;
    }
    get_cell_ion_data(nonemptymgi)[uniqueionindex].eff_ionpot = eff_ionpot;
  } else {
    printout("WARNING! No matching subshells in NT impact ionisation cross section data for Z=%d ionstage %d.\n",
             get_atomicnumber(element), get_ionstage(element, ion));
    printout(
        "-> Defaulting to work function approximation and ionisation energy is not accounted for in Spencer-Fano "
        "solution.\n");

    get_cell_ion_data(nonemptymgi)[uniqueionindex].eff_ionpot = 1. / get_oneoverw(element, ion, modelgridindex);
  }
}

// get the effective ion potential from the stored value
// a value of 0. should be treated as invalid
auto get_eff_ionpot(const int nonemptymgi, const int element, const int ion) {
  return get_cell_ion_data(nonemptymgi)[get_uniqueionindex(element, ion)].eff_ionpot;
  // OR
  // return calculate_eff_ionpot(modelgridindex, element, ion);
}

// Kozma & Fransson 1992 equation 13
// returns the rate coefficient in s^-1
auto nt_ionization_ratecoeff_sf(const int nonemptymgi, const int element, const int ion) -> double {
  const double deposition_rate_density = get_deposition_rate_density(nonemptymgi);
  if (deposition_rate_density > 0.) {
    return deposition_rate_density / get_nnion_tot(nonemptymgi) / get_eff_ionpot(nonemptymgi, element, ion);
    // alternatively, if the y vector is still in memory:
    // return calculate_nt_ionization_ratecoeff(nonemptymgi, element, ion);
  }

  return 0.;
}

// vector of collisional excitation cross sections in cm^2
// epsilon_trans is in erg
// returns the index of the first valid cross section point (en >= epsilon_trans)
// all elements below this index are invalid and should not be used
auto get_xs_excitation_vector(const int element, const int ion, const int lower, const int uptransindex,
                              const double statweight_lower, const double epsilon_trans)
    -> std::tuple<std::array<double, SFPTS>, int> {
  std::array<double, SFPTS> xs_excitation_vec{};
  const auto &uptr = get_uptranslist(element, ion, lower)[uptransindex];
  if (uptr.coll_str >= 0) {
    // collision strength is available, so use it
    // Li et al. 2012 equation 11
    const double constantfactor = std::pow(H_ionpot, 2) / statweight_lower * uptr.coll_str * PI * A_naught_squared;

    const int en_startindex = get_energyindex_ev_gteq(epsilon_trans / EV);

    std::fill_n(xs_excitation_vec.begin(), en_startindex, 0.);

    for (int j = en_startindex; j < SFPTS; j++) {
      const double energy = engrid(j) * EV;
      xs_excitation_vec[j] = constantfactor * std::pow(energy, -2);
    }
    return {xs_excitation_vec, en_startindex};
  }
  if (!uptr.forbidden) {
    const double trans_osc_strength = uptr.osc_strength;
    // permitted E1 electric dipole transitions

    // constexpr double g_bar = 0.2;
    constexpr double A = 0.28;
    constexpr double B = 0.15;

    constexpr double prefactor = 45.585750051;  // 8 * pi^2/sqrt(3)
    const double epsilon_trans_ev = epsilon_trans / EV;

    // Eq 4 of Mewe 1972, possibly from Seaton 1962?
    const double constantfactor =
        epsilon_trans_ev * prefactor * A_naught_squared * std::pow(H_ionpot / epsilon_trans, 2) * trans_osc_strength;

    const int en_startindex = get_energyindex_ev_gteq(epsilon_trans_ev);

    std::fill_n(xs_excitation_vec.begin(), en_startindex, 0.);

    // U = en / epsilon
    // g_bar = A * std::log(U) + b
    // xs[j] = constantfactor * g_bar / engrid(j)
    const double logepsilon = std::log(epsilon_trans_ev);
    for (int j = en_startindex; j < SFPTS; j++) {
      const double logU = logengrid[j] - logepsilon;
      const double g_bar = (A * logU) + B;
      xs_excitation_vec[j] = constantfactor * g_bar / engrid(j);
    }

    return {xs_excitation_vec, en_startindex};
  }

  return {xs_excitation_vec, -1};
}

// Kozma & Fransson equation 9 divided by level population and epsilon_trans
// returns the rate coefficient in s^-1 divided by deposition rate density in erg/cm^3/s
auto calculate_nt_excitation_ratecoeff_perdeposition(const std::array<double, SFPTS> &yvec, const int element,
                                                     const int ion, const int lower, const int uptransindex,
                                                     const double statweight_lower, const double epsilon_trans)
    -> double {
  const auto [xs_excitation_vec, xsstartindex] =
      get_xs_excitation_vector(element, ion, lower, uptransindex, statweight_lower, epsilon_trans);

  if (xsstartindex >= 0) {
    const double y_xs_de =
        cblas_ddot(SFPTS - xsstartindex, xs_excitation_vec.data() + xsstartindex, 1, yvec.data() + xsstartindex, 1) *
        DELTA_E;

    return y_xs_de / E_init_ev / EV;
  }

  return 0.;
}

// returns the energy rate [erg/cm3/s] going toward non-thermal ionisation of lowerion
auto ion_ntion_energyrate(const int modelgridindex, const int element, const int lowerion) -> double {
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  const double nnlowerion = get_nnion(nonemptymgi, element, lowerion);
  double enrate = 0.;
  const auto maxupperion = nt_ionisation_maxupperion(element, lowerion);
  for (int upperion = lowerion + 1; upperion <= maxupperion; upperion++) {
    const double upperionprobfrac = nt_ionization_upperion_probability(nonemptymgi, element, lowerion, upperion, false);
    // for (int lower = 0; lower < get_nlevels(element, lowerion); lower++)
    // {
    //   const double epsilon_trans = epsilon(element, upperion, 0) - epsilon(element, lowerion, lower);
    //   const double nnlower = get_levelpop(nonemptymgi, element, lowerion, lower);
    //   enrate += nnlower * upperionprobfrac * epsilon_trans;
    // }
    const double epsilon_trans = epsilon(element, upperion, 0) - epsilon(element, lowerion, 0);
    enrate += nnlowerion * upperionprobfrac * epsilon_trans;
  }

  const double gamma_nt = nt_ionization_ratecoeff(nonemptymgi, element, lowerion);
  return gamma_nt * enrate;
}

// returns the energy rate [erg/s] going toward non-thermal ionisation in a modelgrid cell
auto get_ntion_energyrate(const int modelgridindex) -> double {
  double ratetotal = 0.;
  for (int ielement = 0; ielement < get_nelements(); ielement++) {
    const int nions = get_nions(ielement);
    for (int ilowerion = 0; ilowerion < nions - 1; ilowerion++) {
      ratetotal += ion_ntion_energyrate(modelgridindex, ielement, ilowerion);
    }
  }
  return ratetotal;
}

auto select_nt_ionization(const int modelgridindex) -> std::tuple<int, int> {
  const double zrand = rng_uniform();

  // // select based on stored frac_deposition for each ion
  // double frac_deposition_ion_sum = 0.;
  // // zrand is between zero and frac_ionization
  // // keep subtracting off deposition fractions of ionizations transitions until we hit the right one
  // // e.g. if zrand was less than frac_dep_trans1, then use the first transition
  // // e.g. if zrand was between frac_dep_trans1 and frac_dep_trans2 then use the second transition, etc
  // for (int allionindex = 0; allionindex < get_includedions(); allionindex++) {
  //   frac_deposition_ion_sum += nt_solution[nonemptymgi].fracdep_ionization_ion[allionindex];
  //   if (frac_deposition_ion_sum >= zrand) {
  //     get_ionfromuniqueionindex(allionindex, element, lowerion);

  //     return;
  //   }
  // }
  // assert_always(false);  // should not reach here

  const double ratetotal = get_ntion_energyrate(modelgridindex);

  // select based on the calculated energy going to ionisation for each ion
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

auto get_uptransindex(const int element, const int ion, const int lower, const int upper) {
  const int nuptrans = get_nuptrans(element, ion, lower);
  const auto *const leveluptrans = get_uptranslist(element, ion, lower);
  for (int t = 0; t < nuptrans; t++) {
    if (upper == leveluptrans[t].targetlevelindex) {
      return t;
    }
  }
  assert_always(false);
  return -1;
}

void analyse_sf_solution(const int nonemptymgi, const int timestep, const bool enable_sfexcitation,
                         const std::array<double, SFPTS> &yfunc) {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const auto nne = grid::get_nne(nonemptymgi);
  const auto nntot = get_nnion_tot(nonemptymgi);
  const auto nnetot = grid::get_nnetot(nonemptymgi);

  double frac_excitation_total = 0.;
  double frac_ionization_total = 0.;

  // temporary storage of the full excitation list for current cell before possible truncation and copying to
  // node-shared memory
  THREADLOCALONHOST std::vector<NonThermalExcitation> tmp_excitation_list;
  tmp_excitation_list.clear();

  for (int element = 0; element < get_nelements(); element++) {
    const int Z = get_atomicnumber(element);
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int uniqueionindex = get_uniqueionindex(element, ion);

      const int ionstage = get_ionstage(element, ion);
      const double nnion = get_nnion(nonemptymgi, element, ion);

      // if (nnion < minionfraction * get_nnion_tot(nonemptymgi)) // skip negligible ions
      if (nnion <= 0.) {  // skip zero-abundance ions
        continue;
      }

      double frac_ionization_ion = 0.;
      double frac_excitation_ion = 0.;
      printout("  Z=%d ionstage %d:\n", Z, ionstage);
      // printout("    nnion: %g\n", nnion);
      printout("    nnion/nntot: %g\n", nnion / nntot);

      calculate_eff_ionpot_auger_rates(nonemptymgi, element, ion, yfunc);

      int matching_subshell_count = 0;
      for (const auto &collionrow : colliondata) {
        if (collionrow.Z == Z && collionrow.ionstage == ionstage) {
          const double frac_ionization_ion_shell =
              calculate_nt_frac_ionization_shell(nonemptymgi, element, ion, collionrow, yfunc);
          frac_ionization_ion += frac_ionization_ion_shell;
          matching_subshell_count++;
          printout("      shell ");
          if (collionrow.n >= 0) {
            printout("n %d, l %d", collionrow.n, collionrow.l);
          } else {
            printout("%s (Lotz)", shellnames.at(-collionrow.l).c_str());
          }
          printout(" I %5.1f eV: frac_ionization %10.4e", collionrow.ionpot_ev, frac_ionization_ion_shell);

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
        get_cell_ion_data(nonemptymgi)[uniqueionindex].fracdep_ionization_ion = frac_ionization_ion;

        frac_ionization_total += frac_ionization_ion;
      } else {
        get_cell_ion_data(nonemptymgi)[uniqueionindex].fracdep_ionization_ion = 0.;
      }
      printout("    frac_ionization: %g (%d subshells)\n", frac_ionization_ion, matching_subshell_count);

      // excitation from all levels is very SLOW
      const int nlevels_all = get_nlevels(element, ion);
      // So limit the lower levels to improve performance
      int nlevels = (nlevels_all > NTEXCITATION_MAXNLEVELS_LOWER) ? NTEXCITATION_MAXNLEVELS_LOWER : nlevels_all;
      if (!enable_sfexcitation) {
        nlevels = -1;  // disable all excitations
      }
      const bool above_minionfraction = (nnion >= minionfraction * get_nnion_tot(nonemptymgi));

      for (int lower = 0; lower < nlevels; lower++) {
        const double statweight_lower = stat_weight(element, ion, lower);
        const int nuptrans = get_nuptrans(element, ion, lower);
        const auto *const uptranslist = get_uptranslist(element, ion, lower);
        const double nnlevel = get_levelpop(nonemptymgi, element, ion, lower);
        const double epsilon_lower = epsilon(element, ion, lower);

        for (int t = 0; t < nuptrans; t++) {
          const int upper = uptranslist[t].targetlevelindex;
          if (upper >= NTEXCITATION_MAXNLEVELS_UPPER) {
            continue;
          }

          const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
          const double ratecoeffperdeposition = calculate_nt_excitation_ratecoeff_perdeposition(
              yfunc, element, ion, lower, t, statweight_lower, epsilon_trans);
          const double frac_excitation_thistrans = nnlevel * epsilon_trans * ratecoeffperdeposition;
          frac_excitation_ion += frac_excitation_thistrans;

          if constexpr (NT_EXCITATION_ON) {
            assert_always(std::isfinite(ratecoeffperdeposition));
            // the atomic data set was limited for Fe V, which caused the ground multiplet to be massively
            // depleted, and then almost no recombination happened!
            if (above_minionfraction && ratecoeffperdeposition > 0 && (Z != 26 || ionstage != 5)) {
              // if (get_coll_str(lineindex) < 0) // if collision strength is not defined, the rate coefficient is
              // unreliable
              //   ratecoeffperdeposition = 0.;

              tmp_excitation_list.push_back({
                  .frac_deposition = frac_excitation_thistrans,
                  .ratecoeffperdeposition = ratecoeffperdeposition,
                  .lineindex = uptranslist[t].lineindex,
              });
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
               get_eff_ionpot(nonemptymgi, element, ion) / EV, (NT_USE_VALENCE_IONPOTENTIAL ? "true" : "false"));

      printout("    workfn approx Gamma:     %9.3e\n", nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion));

      printout("    SF integral Gamma:       %9.3e\n",
               calculate_nt_ionization_ratecoeff(nonemptymgi, element, ion, false, yfunc));

      printout("    SF integral(I=Iv) Gamma: %9.3e  (if always use valence potential)\n",
               calculate_nt_ionization_ratecoeff(nonemptymgi, element, ion, true, yfunc));

      printout("    ARTIS using Gamma:       %9.3e\n", nt_ionization_ratecoeff(nonemptymgi, element, ion));

      // the ion values (unlike shell ones) have been collapsed down to ensure that upperion < nions
      if (ion < nions - 1) {
        printout("    probability to ionstage:");
        double prob_sum = 0.;
        for (int upperion = ion + 1; upperion <= nt_ionisation_maxupperion(element, ion); upperion++) {
          const double probability = nt_ionization_upperion_probability(nonemptymgi, element, ion, upperion, false);
          prob_sum += probability;
          if (probability > 0.) {
            printout(" %d: %.3f", get_ionstage(element, upperion), probability);
          }
        }
        printout("\n");
        assert_always((fabs(prob_sum - 1.0) <= 1e-2) ||
                      (nt_ionization_ratecoeff_sf(nonemptymgi, element, ion) < 1e-20));

        printout("         enfrac to ionstage:");
        double enfrac_sum = 0.;
        for (int upperion = ion + 1; upperion <= nt_ionisation_maxupperion(element, ion); upperion++) {
          const double probability = nt_ionization_upperion_probability(nonemptymgi, element, ion, upperion, true);
          enfrac_sum += probability;
          if (probability > 0.) {
            printout(" %d: %.3f", get_ionstage(element, upperion), probability);
          }
        }
        printout("\n");
        assert_always(fabs(enfrac_sum - 1.0) <= 1e-2 ||
                      (nt_ionization_ratecoeff_sf(nonemptymgi, element, ion) < 1e-20));
      }
    }
  }

  if (nt_excitations_stored > 0) {
    // sort by descending frac_deposition
    std::ranges::SORT_OR_STABLE_SORT(tmp_excitation_list, std::ranges::greater{},
                                     &NonThermalExcitation::frac_deposition);

    // the excitation list is now sorted by frac_deposition descending
    const double deposition_rate_density = get_deposition_rate_density(nonemptymgi);

    if (std::ssize(tmp_excitation_list) > nt_excitations_stored) {
      // truncate the sorted list to save memory
      printout("  Truncating non-thermal excitation list from %zu to %d transitions.\n", tmp_excitation_list.size(),
               nt_excitations_stored);
      tmp_excitation_list.resize(nt_excitations_stored);
    }

    nt_solution[nonemptymgi].frac_excitations_list_size = tmp_excitation_list.size();
    std::ranges::copy(tmp_excitation_list, get_cell_ntexcitations(nonemptymgi).begin());

    printout("[info] mem_usage: non-thermal excitations for cell %d at this timestep occupy %.3f MB\n", modelgridindex,
             nt_solution[nonemptymgi].frac_excitations_list_size * sizeof(NonThermalExcitation) / 1024. / 1024.);

    const auto T_e = grid::get_Te(nonemptymgi);
    printout("  Top non-thermal excitation fractions (total excitations = %d):\n",
             nt_solution[nonemptymgi].frac_excitations_list_size);
    const int ntransdisplayed = std::min(50, nt_solution[nonemptymgi].frac_excitations_list_size);

    for (int excitationindex = 0; excitationindex < ntransdisplayed; excitationindex++) {
      const auto &ntexc = tmp_excitation_list[excitationindex];
      if (ntexc.frac_deposition > 0.) {
        const int lineindex = ntexc.lineindex;
        const auto &line = globals::linelist[lineindex];
        const int element = line.elementindex;
        const int ion = line.ionindex;
        const int lower = line.lowerlevelindex;
        const int upper = line.upperlevelindex;
        const auto nnlevel_lower = get_levelpop(nonemptymgi, element, ion, lower);

        const auto uptransindex = get_uptransindex(element, ion, lower, upper);
        const double epsilon_trans = epsilon(element, ion, upper) - epsilon(element, ion, lower);

        const double ntcollexc_ratecoeff = ntexc.ratecoeffperdeposition * deposition_rate_density;

        const double t_mid = globals::timesteps[timestep].mid;
        const double radexc_ratecoeff = rad_excitation_ratecoeff(nonemptymgi, element, ion, lower, uptransindex,
                                                                 epsilon_trans, nnlevel_lower, lineindex, t_mid);

        const double collexc_ratecoeff = col_excitation_ratecoeff(T_e, nne, element, ion, lower, uptransindex,
                                                                  epsilon_trans, stat_weight(element, ion, lower));

        const double exc_ratecoeff = radexc_ratecoeff + collexc_ratecoeff + ntcollexc_ratecoeff;
        const auto coll_str = get_uptranslist(element, ion, lower)[uptransindex].coll_str;

        printout(
            "    frac_deposition %.3e Z=%2d ionstage %d lower %4d upper %4d rad_exc %.1e coll_exc %.1e nt_exc %.1e "
            "nt/tot %.1e collstr %.1e lineindex %d\n",
            ntexc.frac_deposition, get_atomicnumber(element), get_ionstage(element, ion), lower, upper,
            radexc_ratecoeff, collexc_ratecoeff, ntcollexc_ratecoeff, ntcollexc_ratecoeff / exc_ratecoeff, coll_str,
            lineindex);
      }
    }

    // sort the excitation list by ascending lineindex for fast lookup with a binary search
    std::ranges::SORT_OR_STABLE_SORT(get_cell_ntexcitations(nonemptymgi), std::ranges::less{},
                                     &NonThermalExcitation::lineindex);

  }  // NT_EXCITATION_ON

  // calculate number density of non-thermal electrons
  const double deposition_rate_density_ev = get_deposition_rate_density(nonemptymgi) / EV;
  const double yscalefactor = deposition_rate_density_ev / E_init_ev;

  double nne_nt_max = 0.;
  for (int i = 0; i < SFPTS; i++) {
    const double endash = engrid(i);
    const double delta_endash = DELTA_E;
    const double oneovervelocity = std::sqrt(9.10938e-31 / 2 / endash / 1.60218e-19) / 100;  // in sec/cm
    nne_nt_max += yscalefactor * yfunc[i] * oneovervelocity * delta_endash;
  }

  nt_solution[nonemptymgi].frac_excitation = frac_excitation_total;
  nt_solution[nonemptymgi].frac_ionization = frac_ionization_total;

  printout("  E_init:      %9.2f eV/s/cm^3\n", E_init_ev);
  printout("  deposition:  %9.2f eV/s/cm^3\n", deposition_rate_density_ev);
  printout("  nne:         %9.3e e-/cm^3\n", nne);
  printout("  nnetot:      %9.3e e-/cm^3\n", nnetot);
  printout("  nne_nt     < %9.3e e-/cm^3\n", nne_nt_max);
  printout("  nne_nt/nne < %9.3e\n", nne_nt_max / nne);

  // store the solution properties now while the NT spectrum is in memory (in case we free before packet prop)
  nt_solution[nonemptymgi].frac_heating = calculate_frac_heating(nonemptymgi, yfunc);

  printout("  frac_heating_tot:    %g\n", nt_solution[nonemptymgi].frac_heating);
  printout("  frac_excitation_tot: %g\n", frac_excitation_total);
  printout("  frac_ionization_tot: %g\n", frac_ionization_total);
  const double frac_sum = nt_solution[nonemptymgi].frac_heating + frac_excitation_total + frac_ionization_total;
  printout("  frac_sum:            %g (should be close to 1.0)\n", frac_sum);

  nt_solution[nonemptymgi].frac_heating = 1. - frac_excitation_total - frac_ionization_total;
  printout("  (replacing calculated frac_heating_tot with %g to make frac_sum = 1.0)\n",
           nt_solution[nonemptymgi].frac_heating);
}

void sfmatrix_add_excitation(std::vector<double> &sfmatrixuppertri, const int modelgridindex, const int element,
                             const int ion) {
  // excitation terms
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);

  const int nlevels_all = get_nlevels(element, ion);
  const int nlevels = (nlevels_all > NTEXCITATION_MAXNLEVELS_LOWER) ? NTEXCITATION_MAXNLEVELS_LOWER : nlevels_all;

  for (int lower = 0; lower < nlevels; lower++) {
    const double statweight_lower = stat_weight(element, ion, lower);
    const double nnlevel = get_levelpop(nonemptymgi, element, ion, lower);
    const double epsilon_lower = epsilon(element, ion, lower);
    const int nuptrans = get_nuptrans(element, ion, lower);
    const auto *const uptranslist = get_uptranslist(element, ion, lower);
    for (int t = 0; t < nuptrans; t++) {
      const int upper = uptranslist[t].targetlevelindex;
      if (upper >= NTEXCITATION_MAXNLEVELS_UPPER) {
        continue;
      }
      const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
      const double epsilon_trans_ev = epsilon_trans / EV;
      if (epsilon_trans_ev < SF_EMIN) {
        continue;
      }

      auto [vec_xs_excitation_deltae, xsstartindex] =
          get_xs_excitation_vector(element, ion, lower, t, statweight_lower, epsilon_trans);
      if (xsstartindex >= 0) {
        cblas_dscal(SFPTS - xsstartindex, DELTA_E, vec_xs_excitation_deltae.data() + xsstartindex, 1);

        for (int i = 0; i < SFPTS; i++) {
          const int rowoffset = uppertriangular(i, 0);
          const double en = engrid(i);
          const int stopindex = get_energyindex_ev_lteq(en + epsilon_trans_ev);

          const int startindex = std::max(i, xsstartindex);
          for (int j = startindex; j < stopindex; j++) {
            sfmatrixuppertri[rowoffset + j] += nnlevel * vec_xs_excitation_deltae[j];
          }

          // do the last bit separately because we're not using the full delta_e interval
          const double delta_en_actual = (en + epsilon_trans_ev - engrid(stopindex));
          sfmatrixuppertri[rowoffset + stopindex] +=
              nnlevel * vec_xs_excitation_deltae[stopindex] * delta_en_actual / DELTA_E;
        }
      }
    }
  }
}

void sfmatrix_add_ionization(std::vector<double> &sfmatrixuppertri, const int Z, const int ionstage, const double nnion)
// add the ionization terms to the Spencer-Fano matrix
{
  std::array<double, SFPTS> vec_xs_ionization{};
  for (const auto &collionrow : colliondata) {
    if (collionrow.Z == Z && collionrow.ionstage == ionstage) {
      const double ionpot_ev = collionrow.ionpot_ev;
      const double en_auger_ev = collionrow.en_auger_ev;
      // const double n_auger_elec_avg = colliondata[n].n_auger_elec_avg;
      const double J = get_J(Z, ionstage, ionpot_ev);

      assert_always(ionpot_ev >= SF_EMIN);

      const int xsstartindex = get_xs_ionization_vector(vec_xs_ionization, collionrow);
      assert_always(xsstartindex >= 0);
      // Luke Shingles: the use of min and max on the epsilon limits keeps energies
      // from becoming unphysical. This insight came from reading the
      // CMFGEN Fortran source code (Li, Dessart, Hillier 2012, doi:10.1111/j.1365-2966.2012.21198.x)
      // I had neglected this, so the limits of integration were incorrect. The fix didn't massively affect
      // ionisation rates or spectra, but it was a source of error that led to energy fractions not adding up to
      // 100%
      std::array<double, SFPTS> int_eps_upper = {0};
      std::array<double, SFPTS> prefactors = {0};
      for (int j = xsstartindex; j < SFPTS; j++) {
        const double endash = engrid(j);
        const double epsilon_upper = std::min((endash + ionpot_ev) / 2, endash);
        int_eps_upper[j] = std::atan((epsilon_upper - ionpot_ev) / J);
        prefactors[j] = vec_xs_ionization[j] * nnion / std::atan((endash - ionpot_ev) / 2 / J);
      }

      for (int i = 0; i < SFPTS; i++) {
        // i is the matrix row index, which corresponds to an energy E at which we are solve from y(E)
        const double en = engrid(i);
        const int rowoffset = uppertriangular(i, 0);

        // endash ranges from en to SF_EMAX, but skip over the zero-cross section points
        const int jstart = std::max(i, xsstartindex);
        for (int j = jstart; j < SFPTS; j++) {
          // j is the matrix column index which corresponds to the piece of the integral at y(E') where E' >= E and E'
          // = engrid(j)
          const double endash = engrid(j);

          // J * atan[(epsilon - ionpot_ev) / J] is the indefinite integral of 1/[1 + (epsilon - ionpot_ev)^2/ J^2]
          // in Kozma & Fransson 1992 equation 4

          const double epsilon_lower =
              std::max(endash - en, ionpot_ev);  // and epsilon_upper = (endash + ionpot_ev) / 2;
          const double int_eps_lower = std::atan((epsilon_lower - ionpot_ev) / J);
          if (int_eps_lower <= int_eps_upper[j]) {
            sfmatrixuppertri[rowoffset + j] += prefactors[j] * (int_eps_upper[j] - int_eps_lower) * DELTA_E;
          }
        }

        // below is std::atan((epsilon_lower - ionpot_ev) / J) where epsilon_lower = en + ionpot_ev;
        const double int_eps_lower2 = std::atan(en / J);

        // endash ranges from 2 * en + ionpot_ev to SF_EMAX
        if (2 * en + ionpot_ev <= SF_EMAX) {
          const int secondintegralstartindex = std::max(xsstartindex, get_energyindex_ev_lteq((2 * en) + ionpot_ev));
          for (int j = secondintegralstartindex; j < SFPTS; j++) {
            // epsilon_lower = en + ionpot_ev;
            // epsilon_upper = (endash + ionpot_ev) / 2;
            if (int_eps_lower2 <= int_eps_upper[j]) {
              sfmatrixuppertri[rowoffset + j] -= prefactors[j] * (int_eps_upper[j] - int_eps_lower2) * DELTA_E;
            }
          }
        }
      }

      if constexpr (SF_AUGER_CONTRIBUTION_ON) {
        int augerstopindex = 0;
        if constexpr (SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN) {
          // en_auger_ev is (if LJS understands it correctly) averaged to include some probability of zero Auger
          // electrons so we need a boost to get the average energy of Auger electrons given that there are one or
          // more
          const double en_boost = 1 / (1. - collionrow.prob_num_auger[0]);

          augerstopindex = get_energyindex_ev_gteq(en_auger_ev * en_boost);
        } else {
          augerstopindex = get_energyindex_ev_gteq(en_auger_ev);
        }

        for (int i = 0; i < augerstopindex; i++) {
          const int rowoffset = uppertriangular(i, 0);
          const double en = engrid(i);
          const int jstart = std::max(i, xsstartindex);
          for (int j = jstart; j < SFPTS; j++) {
            const double xs = vec_xs_ionization[j];
            if constexpr (SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN) {
              const double en_boost = 1 / (1. - collionrow.prob_num_auger[0]);
              for (int a = 1; a <= NT_MAX_AUGER_ELECTRONS; a++) {
                if (en < (en_auger_ev * en_boost / a)) {
                  sfmatrixuppertri[rowoffset + j] -= nnion * xs * collionrow.prob_num_auger[a] * a;
                }
              }
            } else {
              assert_always(en < en_auger_ev);
              // printout("SFAuger E %g < en_auger_ev %g so subtracting %g from element with value %g\n", en,
              // en_auger_ev, nnion * xs, ij_contribution);
              sfmatrixuppertri[rowoffset + j] -= nnion * xs;  // * n_auger_elec_avg; // * en_auger_ev???
            }
          }
        }
      }
    }
  }
}

// solve the Spencer-Fano matrix equation and return the y vector (samples of the Spencer-Fano solution function).
// Multiply y by energy interval [eV] to get non-thermal electron number flux. y(E) * dE is the flux of electrons with
// energy in the range (E, E + dE) in units of particles/cm2/s. y has units of particles/cm2/s/eV
auto sfmatrix_solve(const std::vector<double> &sfmatrix) -> std::array<double, SFPTS> {
  std::array<size_t, SFPTS> vec_permutation{};
  gsl_permutation p{.size = SFPTS, .data = vec_permutation.data()};
  gsl_permutation_init(&p);

  const auto gsl_sfmatrix = gsl_matrix_const_view_array(sfmatrix.data(), SFPTS, SFPTS).matrix;

  // sfmatrix must be in upper triangular form
  const auto &gsl_sfmatrix_LU = gsl_sfmatrix;

  // if the matrix is not upper triangular, then do a decomposition
  // make a copy of the matrix for the LU decomp
  // std::array<double, SFPTS> sfmatrix_LU{};
  // auto gsl_sfmatrix_LU = gsl_matrix_view_array(sfmatrix_LU.data(), SFPTS, SFPTS).matrix;
  // gsl_matrix_memcpy(&gsl_sfmatrix_LU, &gsl_sfmatrix);
  // int s{};  // sign of the transformation
  // gsl_linalg_LU_decomp(&gsl_sfmatrix_LU, &p, &s);

  std::array<double, SFPTS> yvec_arr{};
  auto gsl_yvec = gsl_vector_view_array(yvec_arr.data(), SFPTS).vector;

  const auto gsl_rhsvec = gsl_vector_const_view_array(rhsvec.data(), SFPTS).vector;

  // solve matrix equation: sf_matrix * y_vec = rhsvec for yvec
  gsl_linalg_LU_solve(&gsl_sfmatrix_LU, &p, &gsl_rhsvec, &gsl_yvec);

  // refine the solution

  double error_best = -1.;
  std::array<double, SFPTS> yvec_best{};
  auto gsl_yvec_best = gsl_vector_view_array(yvec_best.data(), SFPTS).vector;
  std::array<double, SFPTS> work_vector{};
  auto gsl_work_vector = gsl_vector_view_array(work_vector.data(), SFPTS).vector;
  std::array<double, SFPTS> residual_vector{};
  auto gsl_residual_vector = gsl_vector_view_array(residual_vector.data(), SFPTS).vector;

  int iteration = 0;
  for (iteration = 0; iteration < 10; iteration++) {
    if (iteration > 0) {
      gsl_linalg_LU_refine(&gsl_sfmatrix, &gsl_sfmatrix_LU, &p, &gsl_rhsvec, &gsl_yvec,
                           &gsl_work_vector);  // first argument must be original matrix
    }

    // calculate Ax - b = residual
    gsl_vector_memcpy(&gsl_residual_vector, &gsl_rhsvec);
    gsl_blas_dgemv(CblasNoTrans, 1.0, &gsl_sfmatrix, &gsl_yvec, -1.0, &gsl_residual_vector);

    // value of the largest absolute residual
    const double error = fabs(gsl_vector_get(&gsl_residual_vector, gsl_blas_idamax(&gsl_residual_vector)));

    if (error < error_best || error_best < 0.) {
      gsl_vector_memcpy(&gsl_yvec_best, &gsl_yvec);
      error_best = error;
    }
    // printout("Linear algebra solver iteration %d has a maximum residual of %g\n", iteration, error);
  }
  if (error_best >= 0.) {
    if (error_best > 1e-10) {
      printout("  SF solver LU_refine: After %d iterations, best solution vector has a max residual of %g (WARNING)\n",
               iteration, error_best);
    }
    gsl_vector_memcpy(&gsl_yvec, &gsl_yvec_best);
  }

  if (gsl_vector_isnonneg(&gsl_yvec) == 0) {
    printout("solve_sfmatrix: WARNING: y function goes negative!\n");
  }
  return yvec_arr;
}

}  // anonymous namespace

void init() {
  const ptrdiff_t nonempty_npts_model = grid::get_nonempty_npts_model();

  deposition_rate_density_all_cells = MPI_shared_malloc_span<double>(nonempty_npts_model);

  if (globals::rank_in_node == 0) {
    std::ranges::fill(deposition_rate_density_all_cells, -1.);
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
  printout("  NT_USE_VALENCE_IONPOTENTIAL %s\n", NT_USE_VALENCE_IONPOTENTIAL ? "on" : "off");
  printout("  NT_MAX_AUGER_ELECTRONS %d\n", NT_MAX_AUGER_ELECTRONS);
  printout("  SF_AUGER_CONTRIBUTION %s\n", SF_AUGER_CONTRIBUTION_ON ? "on" : "off");
  printout("  SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN %s\n", SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN ? "on" : "off");

  if (NT_EXCITATION_ON) {
    nt_excitations_stored = std::min(MAX_NT_EXCITATIONS_STORED, get_possible_nt_excitation_count());
    printout("[info] mem_usage: storing %d non-thermal excitations for non-empty cells occupies %.3f MB\n",
             nt_excitations_stored,
             nonempty_npts_model * sizeof(NonThermalExcitation) * nt_excitations_stored / 1024. / 1024.);

    excitations_list_all_cells =
        MPI_shared_malloc_span<NonThermalExcitation>(nonempty_npts_model * nt_excitations_stored);
  }

  ion_data_all_cells = MPI_shared_malloc_span<NonThermalSolutionIon>(nonempty_npts_model * get_includedions());

  nt_solution = MPI_shared_malloc_span<NonThermalCellSolution>(nonempty_npts_model);

  if (globals::rank_in_node == 0) {
    for (ptrdiff_t nonemptymgi = 0; nonemptymgi < nonempty_npts_model; nonemptymgi++) {
      // should make these negative?
      nt_solution[nonemptymgi].frac_heating = 0.97;
      nt_solution[nonemptymgi].frac_ionization = 0.03;
      nt_solution[nonemptymgi].frac_excitation = 0.;

      nt_solution[nonemptymgi].nneperion_when_solved = -1.;
      nt_solution[nonemptymgi].timestep_last_solved = -1;

      zero_all_effionpot(nonemptymgi);

      nt_solution[nonemptymgi].frac_excitations_list_size = 0;
    }
  }
  MPI_Barrier(globals::mpi_comm_node);

  double sourceintegral = 0.;  // integral of S(e) dE
  for (int s = 0; s < SFPTS; s++) {
    sourceintegral += sourcevec(s) * DELTA_E;
  }

  printout("E_init: %14.7e eV/s/cm3\n", E_init_ev);
  printout("source function integral: %14.7e\n", sourceintegral);

  read_collion_data();

  printout("Finished initializing non-thermal solver\n");
}

// set total non-thermal deposition rate from individual gamma/positron/electron/alpha rates. This should be called
// after packet propagation is finished for this timestep and normalise_deposition_estimators() has been done
void calculate_deposition_rate_density(const int nonemptymgi, const int timestep,
                                       HeatingCoolingRates *heatingcoolingrates) {
  heatingcoolingrates->dep_gamma = globals::dep_estimator_gamma[nonemptymgi];

  const double tmid = globals::timesteps[timestep].mid;
  const double rho = grid::get_rho(nonemptymgi);

  // if INSTANT_PARTICLE_DEPOSITION, use the analytic rate at t_mid since it will have no Monte Carlo noise (although
  // strictly, it should be an integral from the timestep start to the end)
  // with time-dependent deposition, we don't have an analytic rate, so we use the Monte Carlo rate
  assert_always(heatingcoolingrates != nullptr);

  heatingcoolingrates->eps_gamma_ana = rho * decay::get_gamma_emission_rate(nonemptymgi, tmid);

  heatingcoolingrates->eps_positron_ana =
      rho * decay::get_particle_injection_rate(nonemptymgi, tmid, decay::DECAYTYPE_BETAPLUS);

  heatingcoolingrates->eps_electron_ana =
      (rho * decay::get_particle_injection_rate(nonemptymgi, tmid, decay::DECAYTYPE_BETAMINUS));

  heatingcoolingrates->eps_alpha_ana =
      rho * decay::get_particle_injection_rate(nonemptymgi, tmid, decay::DECAYTYPE_ALPHA);

  if (PARTICLE_THERMALISATION_SCHEME == ThermalisationScheme::INSTANT) {
    heatingcoolingrates->dep_positron = heatingcoolingrates->eps_positron_ana;
    heatingcoolingrates->dep_electron = heatingcoolingrates->eps_electron_ana;
    heatingcoolingrates->dep_alpha = heatingcoolingrates->eps_alpha_ana;
  } else {
    heatingcoolingrates->dep_positron = globals::dep_estimator_positron[nonemptymgi];
    heatingcoolingrates->dep_electron = globals::dep_estimator_electron[nonemptymgi];
    heatingcoolingrates->dep_alpha = globals::dep_estimator_alpha[nonemptymgi];
  }

  deposition_rate_density_all_cells[nonemptymgi] =
      (heatingcoolingrates->dep_gamma + heatingcoolingrates->dep_positron + heatingcoolingrates->dep_electron);
}

// get non-thermal deposition rate density in erg / s / cm^3 previously stored by calculate_deposition_rate_density()
__host__ __device__ auto get_deposition_rate_density(const int nonemptymgi) -> double {
  assert_always(deposition_rate_density_all_cells[nonemptymgi] >= 0);
  return deposition_rate_density_all_cells[nonemptymgi];
}

void close_file() {
  if (!NT_ON || !NT_SOLVE_SPENCERFANO) {
    return;
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
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  const float frac_heating = nt_solution[nonemptymgi].frac_heating;
  return frac_heating;
}

__host__ __device__ auto nt_ionization_upperion_probability(const int nonemptymgi, const int element,
                                                            const int lowerion, const int upperion,
                                                            const bool energyweighted) -> double {
  assert_always(upperion > lowerion);
  assert_always(upperion < get_nions(element));
  assert_always(upperion <= nt_ionisation_maxupperion(element, lowerion));
  if (NT_SOLVE_SPENCERFANO && NT_MAX_AUGER_ELECTRONS > 0) {
    const int numaugerelec = upperion - lowerion - 1;  // number of Auger electrons to go from lowerin to upper ion
    const int uniqueionindex = get_uniqueionindex(element, lowerion);
    const auto &cell_ion_data = get_cell_ion_data(nonemptymgi)[uniqueionindex];
    if (numaugerelec < NT_MAX_AUGER_ELECTRONS) {
      if (energyweighted) {
        return cell_ion_data.ionenfrac_num_auger[numaugerelec];
      }
      return cell_ion_data.prob_num_auger[numaugerelec];
    }
    if (numaugerelec == NT_MAX_AUGER_ELECTRONS) {
      double prob_remaining = 1.;
      for (int a = 0; a < NT_MAX_AUGER_ELECTRONS; a++) {
        if (energyweighted) {
          prob_remaining -= cell_ion_data.ionenfrac_num_auger[a];
        } else {
          prob_remaining -= cell_ion_data.prob_num_auger[a];
        }
      }
      if (energyweighted) {
        assert_always(fabs(prob_remaining - cell_ion_data.ionenfrac_num_auger[numaugerelec]) < 0.001);
      } else {
        if (fabs(prob_remaining - cell_ion_data.prob_num_auger[numaugerelec]) >= 0.001) {
          printout("Auger probabilities issue for cell %d Z=%02d ionstage %d to %d\n",
                   grid::get_mgi_of_nonemptymgi(nonemptymgi), get_atomicnumber(element),
                   get_ionstage(element, lowerion), get_ionstage(element, upperion));
          for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
            printout("  a %d prob %g\n", a, cell_ion_data.prob_num_auger[a]);
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

__host__ __device__ auto nt_ionisation_maxupperion(const int element, const int lowerion) -> int {
  const int nions = get_nions(element);
  assert_always(lowerion < nions - 1);
  int maxupper = lowerion + 1;

  if (NT_SOLVE_SPENCERFANO) {
    maxupper += NT_MAX_AUGER_ELECTRONS;
  }

  maxupper = std::min(maxupper, nions - 1);

  return maxupper;
}

__host__ __device__ auto nt_random_upperion(const int nonemptymgi, const int element, const int lowerion,
                                            const bool energyweighted) -> int {
  assert_testmodeonly(lowerion < get_nions(element) - 1);
  if (NT_SOLVE_SPENCERFANO && NT_MAX_AUGER_ELECTRONS > 0) {
    while (true) {
      const double zrand = rng_uniform();

      double prob_sum = 0.;
      for (int upperion = lowerion + 1; upperion <= nt_ionisation_maxupperion(element, lowerion); upperion++) {
        prob_sum += nt_ionization_upperion_probability(nonemptymgi, element, lowerion, upperion, energyweighted);

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

__host__ __device__ auto nt_ionization_ratecoeff(const int nonemptymgi, const int element, const int ion) -> double {
  assert_always(NT_ON);
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);

  if (NT_SOLVE_SPENCERFANO) {
    const double Y_nt = nt_ionization_ratecoeff_sf(nonemptymgi, element, ion);
    if (!std::isfinite(Y_nt)) {
      // probably because eff_ionpot = 0 because the solver hasn't been run yet, or no impact ionization cross sections
      // exist
      const double Y_nt_wfapprox = nt_ionization_ratecoeff_wfapprox(modelgridindex, element, ion);
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

#pragma omp declare simd
__host__ __device__ auto nt_excitation_ratecoeff(const int nonemptymgi, const int element, const int ion,
                                                 const int lowerlevel, const int uptransindex, const int lineindex)
    -> double {
  if constexpr (!NT_EXCITATION_ON) {
    return 0.;
  }
  if (lowerlevel >= NTEXCITATION_MAXNLEVELS_LOWER) {
    return 0.;
  }
  const int upperlevel = get_uptranslist(element, ion, lowerlevel)[uptransindex].targetlevelindex;
  if (upperlevel >= NTEXCITATION_MAXNLEVELS_UPPER) {
    return 0.;
  }

  // binary search, assuming the excitation list is sorted by lineindex ascending
  const auto ntexclist = get_cell_ntexcitations(nonemptymgi);
  const auto ntexcitation = std::ranges::lower_bound(ntexclist, lineindex, {}, &NonThermalExcitation::lineindex);
  if (ntexcitation == ntexclist.end() || ntexcitation->lineindex != lineindex) {
    return 0.;
  }

  const double deposition_rate_density = get_deposition_rate_density(nonemptymgi);
  const double ratecoeffperdeposition = ntexcitation->ratecoeffperdeposition;

  return ratecoeffperdeposition * deposition_rate_density;
}

__host__ __device__ void do_ntalpha_deposit(Packet &pkt) {
  // if ionisation by alpha particles is found to be important for the ionisation state, we could do a separate
  // Spencer-Fano solution. For now, just treat alpha deposition as pure heating (even though the alpha deposition rate
  // was calculated from the sum of ionisation and plasma heating)
  atomicadd(nt_energy_deposited, pkt.e_cmf);
  pkt.last_event = 22;
  pkt.type = TYPE_KPKT;
  stats::increment(stats::COUNTER_NT_STAT_TO_KPKT);
}

__host__ __device__ void do_ntlepton_deposit(Packet &pkt) {
  atomicadd(nt_energy_deposited, pkt.e_cmf);

  const int modelgridindex = grid::get_propcell_modelgridindex(pkt.where);
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);

  // macroatom should not be activated in thick cells
  if (NT_ON && NT_SOLVE_SPENCERFANO && grid::modelgrid[nonemptymgi].thick != 1) {
    // here there is some probability to cause ionisation or excitation to a macroatom packet
    // instead of converting directly to k-packet (unless the heating channel is selected)

    double zrand = rng_uniform();
    // zrand is initially between [0, 1), but we will subtract off each
    // component of the deposition fractions
    // until we end and select transition_ij when zrand < dep_frac_transition_ij

    // const double frac_ionization = get_nt_frac_ionization(modelgridindex);
    const double frac_ionization = get_ntion_energyrate(modelgridindex) / get_deposition_rate_density(nonemptymgi);
    // printout("frac_ionization compare %g and %g\n", frac_ionization, get_nt_frac_ionization(modelgridindex));
    // const double frac_ionization = 0.;

    if (zrand < frac_ionization) {
      const auto [element, lowerion] = select_nt_ionization(modelgridindex);
      const int upperion = nt_random_upperion(nonemptymgi, element, lowerion, true);
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
        stats::increment_ion_stats(nonemptymgi, element, lowerion, stats::ION_NTION, pkt.e_cmf / epsilon_trans);
        stats::increment_ion_stats(nonemptymgi, element, upperion, stats::ION_MACROATOM_ENERGYIN_NTCOLLION, pkt.e_cmf);
      }

      do_macroatom(pkt, {.element = element, .ion = upperion, .level = 0, .activatingline = -99});
      return;
    }

    // const double frac_excitation = NT_EXCITATION_ON ? get_nt_frac_excitation(modelgridindex) : 0;
    const double frac_excitation = 0.;
    if (zrand < (frac_ionization + frac_excitation)) {
      zrand -= frac_ionization;
      // now zrand is between zero and frac_excitation
      // the selection algorithm is the same as for the ionization transitions
      for (const auto &ntexcitation : get_cell_ntexcitations(nonemptymgi)) {
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

          do_macroatom(pkt, {.element = element, .ion = ion, .level = upper, .activatingline = -99});
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

// solve the Spencer-Fano equation to get the non-thermal electron flux energy distribution
// based on Equation (2) of Li et al. (2012)
void solve_spencerfano(const int nonemptymgi, const int timestep, const int iteration) {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  bool skip_solution = false;
  if (timestep < globals::num_lte_timesteps + 1) {
    printout("Skipping Spencer-Fano solution for first NLTE timestep\n");
    skip_solution = true;
  } else if (get_deposition_rate_density(nonemptymgi) / EV < MINDEPRATE) {
    printout(
        "Non-thermal deposition rate of %g eV/cm/s/cm^3 below  MINDEPRATE %g in cell %d at timestep %d. Skipping "
        "Spencer-Fano solution.\n",
        get_deposition_rate_density(nonemptymgi) / EV, MINDEPRATE, modelgridindex, timestep);

    skip_solution = true;
  }

  if (skip_solution) {
    // Axelrod values
    nt_solution[nonemptymgi].frac_heating = 0.97;
    nt_solution[nonemptymgi].frac_ionization = 0.03;
    nt_solution[nonemptymgi].frac_excitation = 0.;

    nt_solution[nonemptymgi].nneperion_when_solved = -1.;
    nt_solution[nonemptymgi].timestep_last_solved = -1;

    nt_solution[nonemptymgi].frac_excitations_list_size = 0;

    zero_all_effionpot(nonemptymgi);
    return;
  }

  const auto nne = grid::get_nne(nonemptymgi);  // electrons per cm^3
  const double nne_per_ion = nne / get_nnion_tot(nonemptymgi);
  const double nne_per_ion_last = nt_solution[nonemptymgi].nneperion_when_solved;
  const double nne_per_ion_fracdiff = fabs((nne_per_ion_last / nne_per_ion) - 1.);
  const int timestep_last_solved = nt_solution[nonemptymgi].timestep_last_solved;

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

  nt_solution[nonemptymgi].nneperion_when_solved = nne_per_ion;
  nt_solution[nonemptymgi].timestep_last_solved = timestep;

  const bool enable_sfexcitation = true;
  const bool enable_sfionization = true;
  // if (timestep <= globals::num_lte_timesteps)
  // {
  //   // for the first run of the solver at the first NLTE timestep (which usually requires many iterations),
  //   // do a fast initial solution but mark it has an invalid nne per ion so it gets replaced at the next timestep
  //   nt_solution[nonemptymgi].nneperion_when_solved = -1.;
  //   enable_sfexcitation = false;
  //   enable_sfionization = false;
  //
  //   printout("Doing a fast initial solution without ionization or excitation in the SF equation for the first NLTE
  //   timestep.\n");
  // }
  // if (timestep <= globals::num_lte_timesteps + 2)
  // {
  //   // run the solver in a faster mode for the first couple of NLTE timesteps
  //   // nt_solution[nonemptymgi].nneperion_when_solved = -1.;
  //   enable_sfexcitation = false;
  //   // enable_sfionization = false;
  //
  //   printout("Doing a faster solution without excitation in the SF equation for the first couple of NLTE
  //   timesteps.\n");
  // }

  // sfmatrix will be a compacted upper triangular matrix during construction and then expanded into a full matrix (with
  // lots of zeros) just before the solver is called
  THREADLOCALONHOST std::vector<double> sfmatrix(SFPTS * SFPTS);
  std::fill_n(sfmatrix.begin(), SFPTS * (SFPTS + 1) / 2, 0.);

  // loss terms and source terms
  for (int i = 0; i < SFPTS; i++) {
    sfmatrix[uppertriangular(i, i)] += electron_loss_rate(engrid(i) * EV, nne) / EV;
  }

  if (enable_sfexcitation || enable_sfionization) {
    for (int element = 0; element < get_nelements(); element++) {
      const int Z = get_atomicnumber(element);
      const int nions = get_nions(element);
      bool first_included_ion_of_element = true;
      for (int ion = 0; ion < nions; ion++) {
        const double nnion = get_nnion(nonemptymgi, element, ion);

        // skip negligible ions
        if (nnion < minionfraction * get_nnion_tot(nonemptymgi)) {
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
  // for (int row = 0; row < 10; row++) {
  //   for (int column = 0; column < 10; column++) {
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

  decompactify_triangular_matrix(sfmatrix);
  const auto yfunc = sfmatrix_solve(sfmatrix);

  analyse_sf_solution(nonemptymgi, timestep, enable_sfexcitation, yfunc);
}

void write_restart_data(FILE *gridsave_file) {
  printout("non-thermal solver, ");

  fprintf(gridsave_file, "%d\n", 24724518);  // special number marking the beginning of NT data
  fprintf(gridsave_file, "%d %la %la\n", SFPTS, SF_EMIN, SF_EMAX);

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    fprintf(gridsave_file, "%d %la ", nonemptymgi, deposition_rate_density_all_cells[nonemptymgi]);

    if (NT_ON && NT_SOLVE_SPENCERFANO) {
      check_auger_probabilities(nonemptymgi);

      fprintf(gridsave_file, "%a %a %a %a\n", nt_solution[nonemptymgi].nneperion_when_solved,
              nt_solution[nonemptymgi].frac_heating, nt_solution[nonemptymgi].frac_ionization,
              nt_solution[nonemptymgi].frac_excitation);

      for (int uniqueionindex = 0; uniqueionindex < get_includedions(); uniqueionindex++) {
        fprintf(gridsave_file, "%la ", get_cell_ion_data(nonemptymgi)[uniqueionindex].fracdep_ionization_ion);
        fprintf(gridsave_file, "%a ", get_cell_ion_data(nonemptymgi)[uniqueionindex].eff_ionpot);

        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          fprintf(gridsave_file, "%a %a ", get_cell_ion_data(nonemptymgi)[uniqueionindex].prob_num_auger[a],
                  get_cell_ion_data(nonemptymgi)[uniqueionindex].ionenfrac_num_auger[a]);
        }
      }

      // write NT excitations
      fprintf(gridsave_file, "%d\n", nt_solution[nonemptymgi].frac_excitations_list_size);

      for (const auto &excitation : get_cell_ntexcitations(nonemptymgi)) {
        fprintf(gridsave_file, "%la %la %d\n", excitation.frac_deposition, excitation.ratecoeffperdeposition,
                excitation.lineindex);
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
    int nonemptymgi_in = 0;
    assert_always(fscanf(gridsave_file, "%d %la ", &nonemptymgi_in, &deposition_rate_density_all_cells[nonemptymgi]) ==
                  2);
    assert_always(nonemptymgi_in == nonemptymgi);

    if (NT_ON && NT_SOLVE_SPENCERFANO) {
      assert_always(fscanf(gridsave_file, "%a %a %a %a\n", &nt_solution[nonemptymgi].nneperion_when_solved,
                           &nt_solution[nonemptymgi].frac_heating, &nt_solution[nonemptymgi].frac_ionization,
                           &nt_solution[nonemptymgi].frac_excitation) == 4);

      for (int uniqueionindex = 0; uniqueionindex < get_includedions(); uniqueionindex++) {
        assert_always(
            fscanf(gridsave_file, "%la ", &get_cell_ion_data(nonemptymgi)[uniqueionindex].fracdep_ionization_ion) == 1);
        assert_always(fscanf(gridsave_file, "%a ", &get_cell_ion_data(nonemptymgi)[uniqueionindex].eff_ionpot) == 1);

        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          assert_always(fscanf(gridsave_file, "%a %a ",
                               &get_cell_ion_data(nonemptymgi)[uniqueionindex].prob_num_auger[a],
                               &get_cell_ion_data(nonemptymgi)[uniqueionindex].ionenfrac_num_auger[a]) == 2);
        }
      }

      check_auger_probabilities(nonemptymgi);

      // read NT excitations
      int frac_excitations_list_size_in = 0;
      assert_always(fscanf(gridsave_file, "%d\n", &frac_excitations_list_size_in) == 1);

      nt_solution[nonemptymgi].frac_excitations_list_size = frac_excitations_list_size_in;
      const auto ntexclist = get_cell_ntexcitations(nonemptymgi);

      for (int excitationindex = 0; excitationindex < frac_excitations_list_size_in; excitationindex++) {
        assert_always(fscanf(gridsave_file, "%la %la %d\n", &ntexclist[excitationindex].frac_deposition,
                             &ntexclist[excitationindex].ratecoeffperdeposition,
                             &ntexclist[excitationindex].lineindex) == 3);
      }
    }
  }
}

void nt_MPI_Bcast(const int nonemptymgi, const int root_node_id) {
  MPI_Bcast(&deposition_rate_density_all_cells[nonemptymgi], 1, MPI_DOUBLE, root_node_id, globals::mpi_comm_internode);

  if (NT_ON && NT_SOLVE_SPENCERFANO) {
    if (globals::rank_in_node == 0) {
      MPI_Bcast(&nt_solution[nonemptymgi].nneperion_when_solved, 1, MPI_FLOAT, root_node_id,
                globals::mpi_comm_internode);
      MPI_Bcast(&nt_solution[nonemptymgi].timestep_last_solved, 1, MPI_INT, root_node_id, globals::mpi_comm_internode);
      MPI_Bcast(&nt_solution[nonemptymgi].frac_heating, 1, MPI_FLOAT, root_node_id, globals::mpi_comm_internode);
      MPI_Bcast(&nt_solution[nonemptymgi].frac_ionization, 1, MPI_FLOAT, root_node_id, globals::mpi_comm_internode);
      MPI_Bcast(&nt_solution[nonemptymgi].frac_excitation, 1, MPI_FLOAT, root_node_id, globals::mpi_comm_internode);

      MPI_Bcast(&nt_solution[nonemptymgi].frac_excitations_list_size, 1, MPI_INT, root_node_id,
                globals::mpi_comm_internode);

      MPI_Bcast(get_cell_ntexcitations(nonemptymgi).data(),
                static_cast<size_t>(nt_solution[nonemptymgi].frac_excitations_list_size) * sizeof(NonThermalExcitation),
                MPI_BYTE, root_node_id, globals::mpi_comm_internode);

      const auto ion_data = get_cell_ion_data(nonemptymgi);
      MPI_Bcast(ion_data.data(), ion_data.size() * sizeof(NonThermalSolutionIon), MPI_BYTE, root_node_id,
                globals::mpi_comm_internode);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    check_auger_probabilities(nonemptymgi);
  }
}

void nt_reset_stats() { nt_energy_deposited = 0.; }

void nt_print_stats(const double modelvolume, const double deltat) {
  const double deposition_rate_density_montecarlo = nt_energy_deposited / EV / modelvolume / deltat;

  printout("nt_energy_deposited = %g [eV/s/cm^3]\n", deposition_rate_density_montecarlo);
}

}  // namespace nonthermal
