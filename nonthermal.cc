// Non-thermal energy deposition by the fast leptons produced in radioactive decays: solves the
// Spencer-Fano equation for the electron degradation spectrum (Kozma & Fransson 1992) to obtain
// the fractions of the deposited energy going into heating, ionisation, and excitation, and the
// resulting non-thermal ionisation and excitation rates.

#include "nonthermal.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <ios>
#include <numeric>
#include <ranges>
#include <span>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <Eigen/Core>
#include <Eigen/Dense>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "ltepop.h"
#include "macroatom.h"
#include "mpi_logging.h"
#include "packet.h"
#include "random.h"
#include "sn3d.h"
#include "stats.h"
#include "thermalbalance.h"

namespace nonthermal {

namespace {

// number of energy points in the Spencer-Fano solution vector
constexpr int SFPTS = 4096;

// the minimum and maximum energies for the Spencer-Fano solution vector [eV]
constexpr double SF_EMIN = 0.1;
constexpr double SF_EMAX = 16000;

// number of nodes resolving the integral over the secondary energy epsilon in the first term of Kozma &
// Fransson equation 6. That integral spans at most SF_EMIN, far narrower than one cell of the solution
// energy grid, so it needs a sub-grid of its own rather than being sampled on that grid.
constexpr int NPTS_EPSILON_SUBGRID = 64;

// number of nodes for the integral over E in [0, SF_EMIN] in the third term of Kozma & Fransson equation 3.
// This is a property of that integral alone, not of the solution energy grid. Each node costs a full N_e()
// over every ion and shell, so it dominates the cost of calculate_frac_heating().
constexpr int NPTS_SUB_E0_INTEGRAL = 17;
static_assert(NPTS_SUB_E0_INTEGRAL > 1);

// set true to divide up the mean Auger energy by the number of electrons that come out
constexpr bool SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN = false;

// minimum number fraction of the total population to include in SF solution
constexpr double MIN_ION_OVER_NNTOT = 1.e-8;

// minimum deposition rate density [eV/s/cm3] to solve SF equation
constexpr double MINDEPRATE = 0.;

// Bohr radius squared in cm^2
constexpr double A_naught_squared = 2.800285203e-17;

constexpr std::array shellnames{
    "K ", "L1", "L2", "L3", "M1", "M2", "M3", "M4", "M5", "N1", "N2", "N3", "N4", "N5",
    "N6", "N7", "O1", "O2", "O3", "O4", "O5", "O6", "O7", "P1", "P2", "P3", "P4", "Q1",
};

std::vector<std::vector<double>> elements_electron_binding;
std::vector<std::vector<int>> allions_shell_occupancies;

struct ShellParams {
  int Z{-1};
  int ionstage{-1};
  int n{-1};
  int l{-1};
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

  ShellParams() {
    std::ranges::fill(prob_num_auger, 0.);
    prob_num_auger[0] = 1.;
  }
};

std::vector<ShellParams> colliondata;

static_assert(SF_EMIN > 0.);
static_assert(SF_EMAX > SF_EMIN);
constexpr double DELTA_E = (SF_EMAX - SF_EMIN) / (SFPTS - 1);

// if this is greater than zero, make sure NT_USE_VALENCE_IONPOTENTIAL is false!
static_assert(NT_MAX_AUGER_ELECTRONS == 0 || !NT_USE_VALENCE_IONPOTENTIAL,
              "Overriding the shell potential with the valence potential is not compatible with including Auger "
              "electrons, because the shell potential is used to calculate the energy of Auger electrons.");

// energy grid on which solution is sampled [eV]
constexpr auto engrid(int index) -> double { return SF_EMIN + (index * DELTA_E); }

const auto logengrid = [] {
  std::vector<double> _logengrid(SFPTS);
  for (int i = 0; i < SFPTS; i++) {
    _logengrid[i] = std::log(engrid(i));
  }
  return _logengrid;
}();

// evaluate the source function (distribution of deposited energy) [s^-1 cm^-3 eV^-1] at energy engrid(index)
// This function returns only the spectral shape of the deposited-energy source, spread over a finite energy interval
// and normalized so its integral over energy is 1 (i.e. effective units of eV^-1); the final deposition-rate density
// scaling is applied separately.
constexpr auto sourcevec(const int index) {
  assert_testmodeonly(index >= 0 && index < SFPTS);

  // spread the source over some energy width
  constexpr int source_spread_pts = static_cast<int>(SFPTS * 0.03333) + 1;
  constexpr double source_spread_en = source_spread_pts * DELTA_E;
  constexpr int sourcestartindex = SFPTS - source_spread_pts;

  return (index < sourcestartindex) ? 0. : 1. / source_spread_en;
}

// the energy injection rate density (integral of E * S(e) dE) in eV/s/cm3 that the Spencer-Fano equation is solved for.
// This is arbitrary and and the solution will be scaled to match the actual energy deposition rate density.
constexpr double E_init_ev = [] {
  double integral = 0.;
  for (int s = 0; s < SFPTS; s++) {
    integral += sourcevec(s) * DELTA_E * engrid(s);
  }
  return integral;
}();

// rhs is the constant term (not dependent on y func) in each equation
// rhsvec[i] is the source integral from engrid(i) to SF_EMAX, discretised as a left-endpoint rectangle
// sum that includes point i itself. This matches the convention used for the integrals over y(E') on the
// left-hand side, where sfmatrix_add_excitation() and sfmatrix_add_ionisation() both sum from j = i with
// weight DELTA_E, and the convention used for E_init_ev above.
constexpr auto rhsvec = [] {
  std::array<double, SFPTS> _rhsvec{};
  double source_integral_to_SF_EMAX = 0.;
  for (int i = SFPTS - 1; i >= 0; i--) {
    source_integral_to_SF_EMAX += sourcevec(i);
    _rhsvec[i] = source_integral_to_SF_EMAX * DELTA_E;
  }
  return _rhsvec;
}();

// Monte Carlo result - compare to analytical expectation
// cache-line aligned because it is accumulated atomically by all threads during packet propagation
ALIGNAS_AVOID_FALSE_SHARING double nt_energy_deposited = 0;

struct NonThermalExcitation {
  double frac_deposition;  // the fraction of the non-thermal deposition energy going to the excitation transition
  double ratecoeffperdeposition;  // the excitation rate coefficient divided by the deposition rate density
  int alltransindex;
};

// pointer to either local or node-shared memory excitation list of all cells
MPI_shared_array<NonThermalExcitation> excitations_list_all_cells{};

// the minimum of MAX_NT_EXCITATIONS_STORED and the number of included excitation transitions in the atomic dataset
int nt_excitations_stored = 0;

struct NonThermalSolutionIon {
  float eff_ionpot{0.};  // these are used to calculate the non-thermal ionisation rate
  double fracdep_ionisation_ion{0.};  // the fraction of the non-thermal deposition energy going to ionizing each ion

  // probability that one ionisation of this ion will produce n Auger electrons.
  // items sum to 1.0 for a given ion
  std::array<float, NT_MAX_AUGER_ELECTRONS + 1> prob_num_auger{};
  // like prob_num_auger, but energy weighted. items sum to 1.0 for an ion
  std::array<float, NT_MAX_AUGER_ELECTRONS + 1> ionenfrac_num_auger{};
};

MPI_shared_array<NonThermalSolutionIon> ion_data_all_cells{};

struct NonThermalCellSolution {
  float frac_heating = 1.;  // energy fractions should add up to 1.0 if the solution is good
  float frac_ionisation = 0.;  // fraction of deposition energy going to ionisation
  float frac_excitation = 0.;  // fraction of deposition energy going to excitation

  int frac_excitations_list_size = 0;

  int timestep_last_solved = -1;  // the quantities above were calculated for this timestep
  float nneperion_when_solved{NAN};  // the nne when the solver was last run
};

MPI_shared_array<NonThermalCellSolution> nt_solution;

// Deposition rate density [erg/s/cm3] available to the non-thermal lepton (Spencer-Fano) treatment:
// gamma-ray, positron and electron deposition only. Alpha particles and spontaneous fission fragments
// deposit as pure heating and are accounted for separately in HeatingCoolingRates (see thermalbalance.cc).
MPI_shared_array<double> ntlepton_deposition_rate_density_all_cells;

constexpr auto uppertriangular(const int i, const int j) -> int {
  assert_testmodeonly(i >= 0);
  assert_testmodeonly(i < SFPTS);
  // sometimes you might want to get an offset for a row using j = 0 < i, so that j can be added to it.
  assert_testmodeonly(j >= i || j == 0);
  assert_testmodeonly(j < SFPTS);
  return (SFPTS * i) - (i * (i + 1) / 2) + j;
}

constexpr void decompactify_triangular_matrix(std::span<double> matrix) {
  // the indices work in decreasing order because we're working in-place and
  // don't want to overwrite data that we still need to read
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

auto calculate_ion_shell_occupancies(const int atomic_number, const int nbound,
                                     const std::vector<int>& element_shells_q_neutral) {
  assert_testmodeonly(nbound >= 0);

  const auto shellcount =
      std::min(element_shells_q_neutral.size(), elements_electron_binding[atomic_number - 1].size());
  std::vector<int> element_shells_q;
  reserve_resize(element_shells_q, shellcount);

  int electron_count = 0;
  for (auto shellindex = 0ZU; shellindex < shellcount; shellindex++) {
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

auto read_shell_configs() {
  auto shells_file = fstream_required("electron_shell_occupancy.txt", std::ios::in);

  int nshells = 0;  // number of shell in binding energy file
  int n_z_binding = 0;  // number of elements in file

  std::string line;
  assert_always(get_noncommentline(shells_file, line));
  std::istringstream{line} >> nshells >> n_z_binding;
  printlnlog("Reading electron_shell_occupancy.txt with {} elements and {} shells", n_z_binding, nshells);

  std::vector<std::vector<int>> elements_shells_q;

  elements_shells_q.resize(n_z_binding, std::vector<int>(nshells, 0.));

  assert_always(elements_shells_q.size() == elements_electron_binding.size());

  int zminusone = 0;
  while (get_noncommentline(shells_file, line)) {
    std::istringstream ssline(line);

    int z_element = 0;
    assert_always(ssline >> z_element);
    assert_always(z_element == (zminusone + 1));

    assert_always(elements_shells_q[zminusone].size() == elements_electron_binding[zminusone].size());
    for (int shell = 0; shell < nshells; shell++) {
      int q = 0;
      assert_always(ssline >> q);
      elements_shells_q.at(zminusone).at(shell) = q;
    }
    zminusone++;
  }
  return elements_shells_q;
}

void read_binding_energies() {
  int nshells = 0;  // number of shell in binding energy file
  int n_z_binding = 0;  // number of elements in binding energy file

  constexpr auto filename = "binding_energies_lotz_tab1and2.txt";
  auto binding_energies_file = fstream_required(filename, std::ios::in);

  std::string line;
  assert_always(get_noncommentline(binding_energies_file, line));
  std::istringstream{line} >> nshells >> n_z_binding;
  printlnlog("Reading binding energies file '{}' with {} elements and {} shells", filename, n_z_binding, nshells);

  elements_electron_binding.resize(n_z_binding, std::vector<double>(nshells, 0.));

  for (int zminusone = 0; zminusone < n_z_binding; zminusone++) {
    assert_always(get_noncommentline(binding_energies_file, line));
    std::istringstream ssline(line);
    // new file as an atomic number column
    int z_element{-1};
    ssline >> z_element;
    assert_always(z_element == (zminusone + 1));
    for (int shell = 0; shell < nshells; shell++) {
      float bindingenergy = 0.;
      assert_always(ssline >> bindingenergy);
      elements_electron_binding.at(zminusone).at(shell) = bindingenergy * EV;
    }
  }

  std::vector<std::vector<int>> elements_neutral_shells_q;
  elements_neutral_shells_q = read_shell_configs();

  reserve_resize(allions_shell_occupancies, get_includedions());
  for (int element = 0; element < get_nelements(); element++) {
    for (int ion = 0; ion < get_nions(element); ion++) {
      const int ioncharge = get_ionstage(element, ion) - 1;
      const int atomic_number = get_atomicnumber(element);
      const int nbound = atomic_number - ioncharge;
      if (nbound <= 0) {
        continue;
      }
      allions_shell_occupancies[get_uniqueionindex(element, ion)] =
          calculate_ion_shell_occupancies(atomic_number, nbound, elements_neutral_shells_q.at(atomic_number - 1));
    }
  }
}

[[nodiscard]] auto get_cell_ntexcitations(const ptrdiff_t nonemptymgi) {
  return excitations_list_all_cells.subspan(nonemptymgi * nt_excitations_stored,
                                            nt_solution[nonemptymgi].frac_excitations_list_size);
}

[[nodiscard]] auto get_cell_allions_data(const ptrdiff_t nonemptymgi) {
  return ion_data_all_cells.subspan(nonemptymgi * get_includedions(), get_includedions());
}

auto get_auger_probability(const ptrdiff_t nonemptymgi, const int element, const int ion, const int naugerelec) {
  assert_always(naugerelec <= NT_MAX_AUGER_ELECTRONS);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  return get_cell_allions_data(nonemptymgi)[uniqueionindex].prob_num_auger[naugerelec];
}

auto get_ion_auger_enfrac(const ptrdiff_t nonemptymgi, const int element, const int ion, const int naugerelec) {
  assert_always(naugerelec <= NT_MAX_AUGER_ELECTRONS);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  return get_cell_allions_data(nonemptymgi)[uniqueionindex].ionenfrac_num_auger[naugerelec];
}

void check_auger_probabilities(const ptrdiff_t nonemptymgi) {
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
        printlog(
            "[error] Auger probabilities sum to {:g} (expected 1.0 +/- 0.001) for cell {} Z={} ionstage {}: "
            "P(n_Auger):",
            prob_sum, grid::get_mgi_of_nonemptymgi(nonemptymgi), get_atomicnumber(element), get_ionstage(element, ion));
        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          printlog(" {}:{:g}", a, get_auger_probability(nonemptymgi, element, ion, a));
        }
        printlnlog("");
        problem_found = true;
      }

      if (fabs(ionenfrac_sum - 1.0) > 0.001) {
        printlog(
            "[error] Auger energy fractions sum to {:g} (expected 1.0 +/- 0.001) for cell {} Z={} ionstage {}: "
            "enfrac(n_Auger):",
            ionenfrac_sum, grid::get_mgi_of_nonemptymgi(nonemptymgi), get_atomicnumber(element),
            get_ionstage(element, ion));
        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          printlog(" {}:{:g}", a, get_ion_auger_enfrac(nonemptymgi, element, ion, a));
        }
        printlnlog("");
        problem_found = true;
      }
    }
  }

  assert_always(!problem_found);
}

void read_auger_data() {
  printlnlog("Reading Auger effect data...");
  auto augerfile = fstream_required("auger-km1993-table2.txt", std::ios::in);

  // map x-ray notation shells K L1 L2 L3 M1 M2 M3 to quantum numbers n and l
  constexpr std::array xrayn{1, 2, 2, 2, 3, 3, 3};
  constexpr std::array xrayl{0, 0, 1, 1, 0, 1, 1};
  constexpr std::array xrayg{2, 2, 2, 4, 2, 2, 4};  // g statistical weight = 2j + 1

  std::string strline;
  while (get_noncommentline(augerfile, strline)) {
    int Z = -1;
    int ionstage = -1;
    int shellnum = -1;

    int linepos = 0;
    int offset = 0;

    assert_always(sscanf(strline.c_str(), "%d %d%n", &Z, &ionstage, &offset) == 2);
    assert_always(offset == 5);
    linepos += offset;

    const int element = get_elementindex(Z);

    if (element >= 0 && get_ionstage(element, 0) <= ionstage &&
        ionstage < (get_ionstage(element, 0) + get_nions(element))) {
      float ionpot_ev = -1;
      float en_auger_ev_total_nocorrection = -1;
      int epsilon_e3 = -1;

      assert_always(sscanf(strline.substr(linepos).c_str(), "%d %g %g %d%n", &shellnum, &ionpot_ev,
                           &en_auger_ev_total_nocorrection, &epsilon_e3, &offset) == 4);
      assert_always(offset == 20);

      float n_auger_elec_avg = 0;
      std::array<double, (NT_MAX_AUGER_ELECTRONS + 1)> prob_num_auger{};
      // probability columns are P(1..10 ejected electrons), which includes the primary electron,
      // so column index a corresponds to a Auger electrons
      for (int a = 0; a < 10; a++) {
        linepos = 26 + (a * 5);
        // have to read out exactly 5 characters at a time because the columns are sometimes not separated by a space
        std::array<char, 6> strprob{"00000"};
        assert_always(sscanf(strline.substr(linepos).c_str(), "%5c%n", strprob.data(), &offset) == 1);
        assert_always(offset == 5);
        strprob[5] = '\0';

        int probnaugerelece4 = -1;
        assert_always(sscanf(strprob.data(), "%d", &probnaugerelece4) == 1);

        const double probnaugerelec = probnaugerelece4 / 10000.;

        assert_always(probnaugerelec <= 1.0);

        n_auger_elec_avg += static_cast<float>(a * probnaugerelec);

        if (a <= NT_MAX_AUGER_ELECTRONS) {
          prob_num_auger.at(a) = probnaugerelec;
        } else {
          // add the rates of all higher ionisations to the top one
          prob_num_auger.at(NT_MAX_AUGER_ELECTRONS) += probnaugerelec;
        }
      }

      // use the epsilon correction factor as in equation 7 of Kaastra & Mewe (1993)
      auto en_auger_ev = static_cast<float>(en_auger_ev_total_nocorrection - (epsilon_e3 / 1000. * ionpot_ev));

      assert_always(shellnum > 0);
      assert_always(shellnum <= std::ssize(xrayn));
      const int n = xrayn[shellnum - 1];
      const int l = xrayl[shellnum - 1];
      const int g = xrayg[shellnum - 1];

      if (!std::isfinite(en_auger_ev) || en_auger_ev < 0) {
        printlnlog("  [warning] Z={:2} ionstage {:2} shellnum {} en_auger_ev is {:g}. Setting to zero.", Z, ionstage,
                   shellnum, en_auger_ev);
        en_auger_ev = 0.;
      }

      // now loop through shells with impact ionisation cross sections and apply Auger data that matches n, l values
      for (auto& collionrow : colliondata) {
        if (collionrow.Z == Z && collionrow.ionstage == ionstage && collionrow.n == n && collionrow.l == l) {
          printlog(
              "Z={:2} ionstage {:2} shellnum {} n {} l {} ionpot {:7.2f} E_A {:8.1f} E_A' {:8.1f} epsilon {:6} "
              "<n_Auger> {:5.1f} P(n_Auger)",
              Z, ionstage, shellnum, n, l, ionpot_ev, en_auger_ev_total_nocorrection, en_auger_ev, epsilon_e3,
              n_auger_elec_avg);

          const double prob_sum = std::ranges::fold_left(prob_num_auger, 0.0, std::plus<>());
          for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
            printlog(" {}: {:4.2f}", a, prob_num_auger[a]);
          }
          assert_always(fabs(prob_sum - 1.0) < 0.001);

          printlnlog("");
          const bool found_existing_data = (collionrow.auger_g_accumulated > 0.);

          // keep existing data but update according to statistical weight represented by existing and new data
          const double oldweight = collionrow.auger_g_accumulated / (g + collionrow.auger_g_accumulated);
          const double newweight = g / (g + collionrow.auger_g_accumulated);
          collionrow.auger_g_accumulated += g;

          // update the statistical-weight averaged values
          collionrow.en_auger_ev = static_cast<float>((oldweight * collionrow.en_auger_ev) + (newweight * en_auger_ev));
          collionrow.n_auger_elec_avg =
              static_cast<float>((oldweight * collionrow.n_auger_elec_avg) + (newweight * n_auger_elec_avg));

          for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
            collionrow.prob_num_auger[a] = (oldweight * collionrow.prob_num_auger[a]) + (newweight * prob_num_auger[a]);
          }
          const double prob_sum_final = std::ranges::fold_left(collionrow.prob_num_auger, 0.0, std::plus<>());
          assert_always(fabs(prob_sum_final - 1.0) < 0.001);

          if (found_existing_data) {
            printlog("  same NL shell already has data from another X-ray shell. New g-weighted values: P(n_Auger)");

            for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
              printlog(" {}: {:4.2f}", a, collionrow.prob_num_auger[a]);
            }
            printlnlog("");
          }
        }
      }
    }
  }
}

auto get_sum_q_over_binding_energy(const int element, const int ion) -> double {
  const int ioncharge = get_ionstage(element, ion) - 1;
  const int nbound = get_atomicnumber(element) - ioncharge;  // number of bound electrons

  if (nbound <= 0) {
    return 0.;
  }

  // get the approximate shell occupancy if we don't have the data file
  const auto& shells_q = allions_shell_occupancies[get_uniqueionindex(element, ion)];
  const auto& binding_energies = elements_electron_binding.at(get_atomicnumber(element) - 1);

  double total = 0.;
  const auto num_shells = std::ssize(shells_q);
  for (int shellindex = 0; shellindex < num_shells; shellindex++) {
    const int electronsinshell = shells_q[shellindex];

    if (electronsinshell <= 0) {
      continue;
    }
    double enbinding = binding_energies.at(shellindex);
    if (enbinding <= 0) {
      // if we don't have the shell's binding energy, use the previous one
      assert_always(shellindex > 0);
      enbinding = binding_energies.at(shellindex - 1);
      assert_always(enbinding > 0.);
    }
    total += electronsinshell / std::max(get_ionpot(element, ion), enbinding);
  }

  return total;
}

void read_collion_data() {
  printlnlog("Reading collisional ionisation data from collion.txt...");

  auto cifile = fstream_required("collion.txt", std::ios::in);
  std::string line;
  get_noncommentline(cifile, line);
  std::istringstream ssline(line);
  int colliondatacount = 0;
  assert_always(ssline >> colliondatacount);
  printlnlog("Reading {} collisional transition rows", colliondatacount);
  assert_always(colliondatacount >= 0);

  for (int i = 0; i < colliondatacount; i++) {
    ShellParams collionrow{};
    int nelec = -1;
    get_noncommentline(cifile, line);
    assert_always(std::istringstream{line} >> collionrow.Z >> nelec >> collionrow.n >> collionrow.l >>
                  collionrow.ionpot_ev >> collionrow.A >> collionrow.B >> collionrow.C >> collionrow.D);

    assert_always(nelec > 0);
    collionrow.ionstage = collionrow.Z - nelec + 1;

    const int element = get_elementindex(collionrow.Z);
    if (element < 0 || collionrow.ionstage < get_ionstage(element, 0) ||
        collionrow.ionstage > get_ionstage(element, get_nions(element) - 1)) {
      continue;
    }

    collionrow.en_auger_ev = 0.;
    collionrow.n_auger_elec_avg = 0.;

    colliondata.push_back(collionrow);
  }
  printlnlog("Stored {} of {} input shell cross sections", colliondata.size(), colliondatacount);
  for (int element = 0; element < get_nelements(); element++) {
    const int Z = get_atomicnumber(element);
    for (int ion = 0; ion < get_nions(element); ion++) {
      const int ionstage = get_ionstage(element, ion);
      const int ioncharge = ionstage - 1;
      const int nbound = Z - ioncharge;  // number of bound electrons
      if (nbound <= 0) {
        continue;
      }
      const bool any_data_matched = std::ranges::any_of(colliondata, [Z, ionstage](const ShellParams& collionrow) {
        return collionrow.Z == Z && collionrow.ionstage == ionstage;
      });
      if (!any_data_matched) {
        const double ionpot_ev = get_ionpot(element, ion) / EV;
        printlnlog(
            "No collisional ionisation data for Z={} ionstage {}. Using Lotz approximation with ionpot = {:g} [eV]", Z,
            ionstage, ionpot_ev);

        // get the approximate shell occupancy if we don't have the data file
        const auto& shells_q = allions_shell_occupancies[get_uniqueionindex(element, ion)];
        int electron_count = 0;
        const auto num_shells = std::ssize(shells_q);
        for (int shellindex = 0; shellindex < num_shells; shellindex++) {
          const int electronsinshell = shells_q.at(shellindex);

          electron_count += electronsinshell;

          if (electronsinshell <= 0) {
            continue;
          }
          double enbinding = elements_electron_binding.at(Z - 1).at(shellindex);
          const double ionpot = ionpot_ev * EV;
          if (enbinding <= 0) {
            // if we don't have the shell's binding energy, use the previous one
            assert_always(shellindex > 0);  // first shell binding energy must be positive
            enbinding = elements_electron_binding.at(Z - 1).at(shellindex - 1);
            assert_always(enbinding > 0);
          }
          const double p = std::max(ionpot, enbinding);
          ShellParams collionrow{};
          collionrow.Z = Z;
          collionrow.ionstage = ionstage;
          collionrow.n = -1;
          collionrow.l = -shellindex;
          collionrow.ionpot_ev = p / EV;
          collionrow.A = -1.;
          collionrow.B = -1.;
          collionrow.C = -1.;
          collionrow.D = -1.;
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
  std::ranges::stable_sort(colliondata, [](const ShellParams& a, const ShellParams& b) {
    return std::tie(a.Z, a.ionstage, a.ionpot_ev, a.n, a.l) < std::tie(b.Z, b.ionstage, b.ionpot_ev, b.n, b.l);
  });

  // this condition is checked here once per ion (it does not depend on the cell), so that
  // calculate_eff_ionpot_auger_rates does not have to report it for every cell
  for (int element = 0; element < get_nelements(); element++) {
    const int Z = get_atomicnumber(element);
    for (int ion = 0; ion < get_nions(element); ion++) {
      const int ionstage = get_ionstage(element, ion);
      if ((Z - (ionstage - 1)) <= 0) {
        continue;  // no bound electrons, so no NT impact ionisation
      }
      const bool any_subshells = std::ranges::any_of(colliondata, [Z, ionstage](const ShellParams& collionrow) {
        return collionrow.Z == Z && collionrow.ionstage == ionstage;
      });
      if (!any_subshells) {
        printlnlog(
            "[warning] no NT impact ionisation subshell data for Z={} ionstage {}: the Spencer-Fano solver will "
            "default to the work function approximation and not account for the ionisation energy",
            Z, ionstage);
      }
    }
  }

  if (NT_MAX_AUGER_ELECTRONS > 0) {
    read_auger_data();
  }
}

auto get_possible_nt_excitation_count() -> int {
  // count the number of excitation transitions that pass the MAXNLEVELS_LOWER and MAXNLEVELS_UPPER conditions
  // this count might be higher than the number of stored power_per_source_mass due to the MAX_NT_EXCITATIONS_STORED
  // limit
  int ntexcitationcount = 0;
  for (int element = 0; element < get_nelements(); element++) {
    for (int ion = 0; ion < get_nions(element); ion++) {
      const int lower_nlevels = std::min(NTEXCITATION_MAXNLEVELS_LOWER, get_nlevels(element, ion));
      for (int lower = 0; lower < lower_nlevels; lower++) {
        const int nuptrans = get_nuptrans(element, ion, lower);
        const auto alltrans_startup = get_alltrans_startup(element, ion, lower);
        for (int t = 0; t < nuptrans; t++) {
          const int upper = globals::alltrans.targetlevelindex[alltrans_startup + t];
          if (upper < NTEXCITATION_MAXNLEVELS_UPPER) {
            ntexcitationcount++;
          }
        }
      }
    }
  }
  return ntexcitationcount;
}

void zero_all_effionpot(const ptrdiff_t nonemptymgi) {
  for (int uniqueionindex = 0; uniqueionindex < get_includedions(); uniqueionindex++) {
    auto& celliondata = get_cell_allions_data(nonemptymgi)[uniqueionindex];
    celliondata.eff_ionpot = 0.;

    std::ranges::fill(celliondata.prob_num_auger, 0.);
    celliondata.prob_num_auger[0] = 1.;
    std::ranges::fill(celliondata.ionenfrac_num_auger, 0.);
    celliondata.ionenfrac_num_auger[0] = 1.;

    const auto [element, ion] = get_ionfromuniqueionindex(uniqueionindex);
    assert_always(fabs(get_auger_probability(nonemptymgi, element, ion, 0) - 1.0) < 1e-3);
    assert_always(fabs(get_ion_auger_enfrac(nonemptymgi, element, ion, 0) - 1.0) < 1e-3);
  }
  check_auger_probabilities(nonemptymgi);
}

[[nodiscard]] constexpr auto get_energyindex_ev_lteq(const double energy_ev) -> int
// finds the highest energy point <= energy_ev
{
  const auto index = static_cast<int>(get_linearbinindex(energy_ev, SF_EMIN, DELTA_E));

  return std::clamp(index, 0, SFPTS - 1);
}

[[nodiscard]] constexpr auto get_energyindex_ev_gteq(const double energy_ev) -> int
// finds the lowest energy point >= energy_ev
{
  const int index = std::ceil((energy_ev - SF_EMIN) / DELTA_E);

  return std::clamp(index, 0, SFPTS - 1);
}

// interpolate the y flux values to get the value at a given energy
// y has units of particles / cm2 / s / eV
[[nodiscard]] constexpr auto get_y(const std::array<double, SFPTS>& yfunc, const double energy_ev) -> double {
  // energies below SF_EMIN give a negative index and return zero
  const auto index = static_cast<int>(get_linearbinindex(energy_ev, SF_EMIN, DELTA_E));

  if (index < 0 || index >= (SFPTS - 1)) {
    return 0.;
  }
  const double enbelow = engrid(index);
  const double enabove = engrid(index + 1);
  const double ybelow = yfunc[index];
  const double yabove = yfunc[index + 1];
  return std::lerp(ybelow, yabove, (energy_ev - enbelow) / (enabove - enbelow));
}

constexpr auto xs_ionisation_lotz(const double en_erg, const ShellParams& colliondata_ion, const int electronsinshell)
    -> double {
  const double ionpot = colliondata_ion.ionpot_ev * EV;
  if (en_erg < ionpot) {
    return 0.;
  }
  // relativistic, to match the relativistic correction terms in the Axelrod equation below. The classical
  // betasq = 2 * en_erg / (ME * pow2(CLIGHT)) is 5% high at the 16 keV top of the energy grid, 30% high at 100 keV,
  // and reaches one at 255 keV, where log(1 - betasq) would be undefined and the cross section spuriously zero.
  const double gamma = (en_erg / (ME * pow2(CLIGHT))) + 1.;
  const double betasq = 1. - (1. / pow2(gamma));
  // 0 < betasq < 1 holds for any finite en_erg >= ionpot, so this only fires on a non-finite input. Both log terms
  // below need it, and without the check a NaN would make part_sigma_shell fail its > 0 test and return a silent zero.
  assert_testmodeonly(betasq > 0. && betasq < 1.);

  const int ioncharge = colliondata_ion.ionstage - 1;
  const int nbound = colliondata_ion.Z - ioncharge;  // number of bound electrons

  if (nbound <= 0) {
    return 0.;
  }

  // Equation 3.38 of Axelrod (1980) attributed to Lotz (1967)
  // WARNING: The Axelrod equation uses both ln() and log10(), but the log10() term is likely a typo and has been
  // corrected to ln(). Fortunately, at our typical 16 keV value of EMAX, 511 keV electrons are only mildly
  // relativistic and the log(1 - beta^2) term is small anyway.
  const double part_sigma_shell =
      (electronsinshell / ionpot *
       (std::log(betasq * ME * pow2(CLIGHT) / 2.0 / ionpot) - std::log(1 - betasq) - betasq));
  if (part_sigma_shell > 0.) {
    // See the comment about Aconst in get_oneoverw_approx_axelrod()
    constexpr double Aconst = 1.33e-14 * EV * EV;
    const double sigma = 2 * Aconst / ME / (betasq * pow2(CLIGHT)) * part_sigma_shell;
    assert_always(sigma >= 0);
    return sigma;
  }

  return 0.;
}

auto get_xs_ionisation_vector_lotz(std::array<double, SFPTS>& xs_vec, const ShellParams& colliondata_ion) -> int {
  const int element = get_elementindex(colliondata_ion.Z);
  const int ion = colliondata_ion.ionstage - get_ionstage(element, 0);
  const int shellindex = -colliondata_ion.l;
  const int electronsinshell = allions_shell_occupancies[get_uniqueionindex(element, ion)][shellindex];

  const double ionpot_ev = colliondata_ion.ionpot_ev;
  const int startindex = get_energyindex_ev_gteq(ionpot_ev);
  std::fill_n(xs_vec.begin(), startindex, 0.);

  for (int i = startindex; i < SFPTS; i++) {
    xs_vec[i] = xs_ionisation_lotz(engrid(i) * EV, colliondata_ion, electronsinshell);
  }

  return startindex;
}

// xs_vec will be set with impact ionisation cross sections [cm2] for E > ionpot_ev (and zeros below this energy)
// returns the index of the first energy point >= ionpot_ev, or SFPTS (one past the last index, with xs_vec
// all zeros) when ionpot_ev is above the top of the energy grid
auto get_xs_ionisation_vector(std::array<double, SFPTS>& xs_vec, const ShellParams& colliondata_ion) -> int {
  const double ionpot_ev = colliondata_ion.ionpot_ev;
  if (ionpot_ev > SF_EMAX) {
    // every grid point is below the ionisation threshold, so all cross sections are zero.
    // Without this, the clamped startindex below would evaluate the Younger fitting formula at
    // u = E / I < 1, where it goes negative.
    xs_vec.fill(0.);
    return SFPTS;
  }
  const double A = colliondata_ion.A;
  if (A < 0) {
    return get_xs_ionisation_vector_lotz(xs_vec, colliondata_ion);
  }

  const int startindex = get_energyindex_ev_gteq(ionpot_ev);

  std::fill_n(xs_vec.begin(), startindex, 0.);

  const double B = colliondata_ion.B;
  const double C = colliondata_ion.C;
  const double D = colliondata_ion.D;

  for (int i = startindex; i < SFPTS; i++) {
    const double u = engrid(i) / ionpot_ev;
    const double xs_ioniz =
        1e-14 * ((A * (1 - (1 / u))) + (B * pow2(1 - (1 / u))) + (C * std::log(u)) + (D * std::log(u) / u)) /
        (u * pow2(ionpot_ev));
    xs_vec[i] = xs_ioniz;
  }

  return startindex;
}

// distribution of secondary electron energies for primary electron with energy e_p
// Opal, Peterson, & Beaty (1971)
[[nodiscard]] constexpr auto Psecondary(const double e_p, const double en_epsilon, const double I, const double J)
    -> double {
  const double e_s = en_epsilon - I;

  if (e_p <= I || e_s < 0.) {
    return 0.;
  }
  assert_testmodeonly(J > 0);
  assert_testmodeonly(e_p >= I);
  assert_testmodeonly(e_s >= 0);
  assert_testmodeonly(std::isfinite(std::atan((e_p - I) / 2 / J)));
  return 1 / (J * std::atan((e_p - I) / 2 / J) * (1 + pow2(e_s / J)));
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

  const auto alltrans_startup = get_alltrans_startup(element, ion, lower);
  const auto alltransindex = alltrans_startup + uptransindex;
  if (globals::alltrans.coll_str[alltransindex] >= 0) {
    // collision strength is available, so use it
    // Li et al. 2012 equation 11: sigma = pi * a_0^2 * (I_H / E) * Omega / g_lower,
    // with k_i^2 = E / I_H in units of the inverse Bohr radius squared
    return (H_ionpot / energy) / lowerstatweight * globals::alltrans.coll_str[alltransindex] * PI * A_naught_squared;
  }
  if (!globals::alltrans.forbidden[alltransindex]) {
    // permitted E1 electric dipole transitions
    const double U = energy / epsilon_trans;

    constexpr double A = 0.28;
    constexpr double B = 0.15;
    const double g_bar = (A * std::log(U)) + B;

    constexpr double prefactor = 45.585750051;  // 8 * pi^2/sqrt(3)
    // Eq 4 of Mewe 1972, possibly from Seaton 1962?
    return prefactor * A_naught_squared * pow2(H_ionpot / epsilon_trans) *
           globals::alltrans.osc_strength[alltransindex] * g_bar / U;
  }
  return 0.;
}

// -dE / dx for fast electrons
// energy is in ergs
// nne is the thermal electron density [cm^-3]
// return value has units of erg/cm
constexpr auto electron_loss_rate(const double energy, const double nne) -> double {
  // with no thermal electrons there is no Coulomb energy loss. Without this guard the plasma
  // frequency is zero and the Coulomb logarithm diverges, giving 0 * inf = NaN in a fully neutral
  // cell (nne can reach exactly zero because the MINPOP floor is a float denormal that is flushed
  // to zero by the default -ffast-math build)
  if (energy <= 0. || nne <= 0.) {
    return 0;
  }

  // normally set to 1.0, but Shingles et al. (2021) boosted this to increase heating
  constexpr double boostfactor = 1.;

  const double omegap = std::sqrt(4 * PI * nne * pow2(QE) / ME);
  const double zetae = H * omegap / 2 / PI;
  if (energy > 14 * EV) {
    return boostfactor * nne * 2 * PI * pow4(QE) / energy * std::log(2 * energy / zetae);
  }
  const double v = std::sqrt(2 * energy / ME);
  // Kozma & Fransson (1992) eq. 2 describes the gamma in this Coulomb logarithm as "Euler's constant
  // (Schunk & Hays 1971)", but Schunk & Hays (1971, p. 114) define it by "ln gamma is Euler's
  // constant", i.e. gamma = exp(0.5772) = 1.781 rather than 0.5772 itself.
  return boostfactor * nne * 2 * PI * pow4(QE) / energy * std::log(ME * pow3(v) / (EXP_EULERGAMMA * pow2(QE) * omegap));
}

// impact ionisation cross section in cm^2
// energy and ionisation_potential should be in eV
// fitting formula of Younger 1981
// called Q_i(E) in KF92 equation 7
constexpr auto xs_impactionisation(const double energy_ev, const ShellParams& colliondata_ion) -> double {
  const double ionpot_ev = colliondata_ion.ionpot_ev;
  const double u = energy_ev / ionpot_ev;

  if (u <= 1.) {
    return 0;
  }
  const double A = colliondata_ion.A;
  if (A < 0) {
    const int element = get_elementindex(colliondata_ion.Z);
    const int ion = colliondata_ion.ionstage - get_ionstage(element, 0);
    const int shellindex = -colliondata_ion.l;
    const int electronsinshell = allions_shell_occupancies[get_uniqueionindex(element, ion)][shellindex];
    return xs_ionisation_lotz(energy_ev * EV, colliondata_ion, electronsinshell);
  }

  const double B = colliondata_ion.B;
  const double C = colliondata_ion.C;
  const double D = colliondata_ion.D;

  return 1e-14 * ((A * (1 - (1 / u))) + (B * pow2((1 - (1 / u)))) + (C * std::log(u)) + (D * std::log(u) / u)) /
         (u * pow2(ionpot_ev));
}

// Kozma & Fransson equation 6.
// Something related to a number of electrons, needed to calculate the heating fraction in equation 3
// not valid for energy > SF_EMIN
auto N_e(const int nonemptymgi, const double energy, const std::array<double, SFPTS>& yfunc) -> double {
  const double energy_ev = energy / EV;
  const double nnion_tot = get_nnion_tot(nonemptymgi);
  double N_e_total = 0.;  // named to avoid shadowing the enclosing function N_e()

  for (int element = 0; element < get_nelements(); element++) {
    const int Z = get_atomicnumber(element);
    const int nions = get_nions(element);

    for (int ion = 0; ion < nions; ion++) {
      double N_e_ion = 0.;
      const int ionstage = get_ionstage(element, ion);
      const double nnion = get_nnion(nonemptymgi, element, ion);

      if (nnion < MIN_ION_OVER_NNTOT * nnion_tot) {  // skip negligible ions
        continue;
      }

      // excitation terms

      const int nlevels_all = get_nlevels(element, ion);
      const int nlevels = std::min(NTEXCITATION_MAXNLEVELS_LOWER, nlevels_all);

      for (int lower = 0; lower < nlevels; lower++) {
        const auto uniquelevelindex = get_uniquelevelindex(element, ion, lower);
        const double nnlevel = calculate_levelpop(nonemptymgi, element, ion, lower);
        const double epsilon_lower = epsilon(uniquelevelindex);
        const auto statweight_lower = stat_weight(uniquelevelindex);
        const int nuptrans = get_nuptrans(uniquelevelindex);
        const auto alltrans_startup = get_alltrans_startup(uniquelevelindex);
        for (int t = 0; t < nuptrans; t++) {
          const int upper = globals::alltrans.targetlevelindex[alltrans_startup + t];
          if (upper >= NTEXCITATION_MAXNLEVELS_UPPER) {
            continue;
          }
          const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
          const double epsilon_trans_ev = epsilon_trans / EV;
          if (epsilon_trans_ev < SF_EMIN) {
            // sfmatrix_add_excitation() does not include these transitions in the solution, so skip
            // them here too for a consistent energy accounting
            continue;
          }
          N_e_ion += (nnlevel / nnion) * get_y(yfunc, energy_ev + epsilon_trans_ev) *
                     xs_excitation(element, ion, lower, t, epsilon_trans, statweight_lower, energy + epsilon_trans);
        }
      }

      // ionisation terms
      // sfmatrix_add_ionisation() is only called for ion < nions - 1 (the top ion is not ionised in
      // the solution), so skip the top ion here too for a consistent energy accounting.
      // NB: the Auger-electron source that sfmatrix_add_ionisation() injects has no counterpart here
      if (ion < nions - 1) {
        for (const auto& collionrow : colliondata) {
          if (collionrow.Z == Z && collionrow.ionstage == ionstage) {
            const double ionpot_ev = collionrow.ionpot_ev;
            const double J = get_J(Z, ionstage, ionpot_ev);
            const double lambda = std::min(SF_EMAX - energy_ev, energy_ev + ionpot_ev);

            // integral from ionpot up to lambda, over its own sub-grid. This domain is only
            // (lambda - ionpot_ev) <= energy_ev < SF_EMIN wide, some forty times narrower than one cell of
            // the solution grid, so sampling it on that grid was wrong twice over: it gave a full DELTA_E of
            // weight to the whole domain, and it snapped the limits down onto the grid point at or below
            // ionpot_ev, where Psecondary() sees a negative secondary energy and returns zero. The term was
            // therefore dropped entirely for most shells, and for the few whose threshold sits just below a
            // grid point it instead contributed that oversized weight next to the 1/atan((e_p - I) / 2J)
            // normalisation pole. Midpoint sampling keeps every node strictly inside the domain, away from
            // both of those edges.
            if (lambda > ionpot_ev) {
              const double delta_epsilon = (lambda - ionpot_ev) / NPTS_EPSILON_SUBGRID;
              for (int i = 0; i < NPTS_EPSILON_SUBGRID; i++) {
                const double endash = ionpot_ev + ((i + 0.5) * delta_epsilon);

                N_e_ion += get_y(yfunc, energy_ev + endash) * xs_impactionisation(energy_ev + endash, collionrow) *
                           Psecondary(energy_ev + endash, endash, ionpot_ev, J) * delta_epsilon;
              }
            }

            // integral from 2E + I up to E_max. The lower limit falls between grid points, so the
            // grid-aligned sum starts at the first point at or above it and the leading partial cell is
            // added separately. Starting at the point at or below the limit instead placed a node under the
            // ionisation threshold: usually harmless because Psecondary() returns zero there, but for a
            // shell whose threshold falls within 2E of a grid point that node cleared the e_p > I guard and
            // contributed a spike near the normalisation pole in place of a finite value.
            const double integral2_min = (2 * energy_ev) + ionpot_ev;
            if (integral2_min < SF_EMAX) {
              const int integral2startindex = get_energyindex_ev_gteq(integral2_min);

              const double partialcellwidth = engrid(integral2startindex) - integral2_min;
              if (partialcellwidth > 0.) {
                const double endash = integral2_min + (partialcellwidth / 2.);
                N_e_ion += get_y(yfunc, endash) * xs_impactionisation(endash, collionrow) *
                           Psecondary(endash, energy_ev + ionpot_ev, ionpot_ev, J) * partialcellwidth;
              }

              for (int i = integral2startindex; i < SFPTS; i++) {
                const double endash = engrid(i);
                N_e_ion += yfunc[i] * xs_impactionisation(endash, collionrow) *
                           Psecondary(endash, energy_ev + ionpot_ev, ionpot_ev, J) * DELTA_E;
              }
            }
          }
        }
      }

      N_e_total += nnion * N_e_ion;
    }
  }

  // source term, should be zero at the low end anyway
  N_e_total += sourcevec(get_energyindex_ev_lteq(energy_ev));

  assert_always(std::isfinite(N_e_total));
  return N_e_total;
}

// fraction of deposited energy that goes into heating the thermal electrons
// Kozma & Fransson equation 3
auto calculate_frac_heating(const int nonemptymgi, const std::array<double, SFPTS>& yfunc) -> float {
  // frac_heating multiplied by E_init, which will be divided out at the end
  double frac_heating_Einit = 0.;
  const auto nne = grid::get_nne(nonemptymgi);

  for (int i = 0; i < SFPTS; i++) {
    const double endash = engrid(i);

    // first term
    frac_heating_Einit += yfunc[i] * (electron_loss_rate(endash * EV, nne) / EV) * DELTA_E;
  }

  // second term
  frac_heating_Einit += SF_EMIN * get_y(yfunc, SF_EMIN) * (electron_loss_rate(SF_EMIN * EV, nne) / EV);

  double N_e_contrib_Einit = 0.;
  // third term (integral from zero to SF_EMIN), on a sub-grid of its own. Sizing the node count from
  // DELTA_E tied it to the solution energy grid, which this integral has nothing to do with, and left only
  // nine interior nodes. Integrating by trapezoid rather than by summing interior nodes also recovers the
  // half node that was dropped at the SF_EMIN end; the zero end costs nothing, since the integrand carries
  // a factor of endash and so vanishes there.
  constexpr double delta_endash = SF_EMIN / (NPTS_SUB_E0_INTEGRAL - 1);
  for (int j = 1; j < NPTS_SUB_E0_INTEGRAL; j++) {
    const double endash = delta_endash * j;
    const double weight = (j == (NPTS_SUB_E0_INTEGRAL - 1)) ? 0.5 : 1.;
    N_e_contrib_Einit += weight * N_e(nonemptymgi, endash * EV, yfunc) * endash * delta_endash;
  }
  frac_heating_Einit += N_e_contrib_Einit;

  const auto frac_heating = static_cast<float>(frac_heating_Einit / E_init_ev);

  if (!std::isfinite(frac_heating) || frac_heating < 0 || frac_heating > 1.0) {
    printlnlog("[warning] calculate_frac_heating: cell {} ts {}: invalid result of {:g}. Setting to 1.0 instead",
               grid::get_mgi_of_nonemptymgi(nonemptymgi), globals::timestep, frac_heating);
    return 1.;
  }

  return frac_heating;
}

// fraction of deposited energy that goes into ionisation
auto get_nt_frac_ionisation(const int nonemptymgi) -> float {
  if (!NT_ON) {
    return 0.;
  }
  if (!NT_SOLVE_SPENCERFANO) {
    return 0.03;  // Axelrod 1980 approximation
  }

  assert_always(nt_solution[nonemptymgi].frac_ionisation >= 0.);

  return nt_solution[nonemptymgi].frac_ionisation;
}

// fraction of deposited energy that goes into collisional excitation
auto get_nt_frac_excitation(const int nonemptymgi) -> float {
  if (!NT_ON || !NT_SOLVE_SPENCERFANO) {
    return 0.;
  }

  const float frac_excitation = nt_solution[nonemptymgi].frac_excitation;
  assert_always(frac_excitation >= 0.);
  assert_always(std::isfinite(frac_excitation));

  return frac_excitation;
}

// compute the approximate work per ion pair without using the Spencer-Fano equation
// Makes use of EXTREMELY SIMPLE approximations - high energy limits only (can be used as an alternative to the
// Spencer-Fano solver)
// Luke (2026-02-20): This is so far from the Spencer-Fano eff_ionpot that I think there is a mistake in the calculation
auto get_oneoverw_approx_axelrod(const int element, const int ion, const int nonemptymgi) -> double {
  // Work in terms of 1/W since this is actually what we want. It is given by sigma/(Latom + Lelec).
  // We are going to start by taking all the high energy limits and ignoring Lelec, so that the
  // denominator is extremely simplified. Need to get the mean Z value.

  double nntot = 0.;
  double Zbar = 0.;  // number-weighted average atomic number
  for (int ielement = 0; ielement < get_nelements(); ielement++) {
    const double nnelement = grid::get_elem_numberdens(nonemptymgi, ielement);
    Zbar += nnelement * get_atomicnumber(ielement);
    nntot += nnelement;
  }
  if (nntot > 0) {
    Zbar /= nntot;
  }

  const double binding = get_sum_q_over_binding_energy(element, ion);
  // Axelrod 1980 says the constant A = 1.33e-14 [cm^2 eV^2] has been determined by normalising to the average of the
  // values given by Jacobs et al (1979) and McGuire (1977) at 10 keV. However, this reduces the accuracy of the
  // approximation at lower energies, and since we are mostly interested in the low energy end of the spectrum for
  // calculating the heating and ionisation fractions, it would be better to use Lotz value of A = 4.5e-14 [cm2 eV2].
  constexpr double Aconst = 1.33e-14 * EV * EV;

  return Aconst * binding / Zbar / (2 * PI * pow4(QE));
}

// the fraction of deposited energy that goes into ionising electrons in a particular shell
auto calculate_nt_frac_ionisation_shell(const int nonemptymgi, const int element, const int ion,
                                        const ShellParams& collionrow, const std::array<double, SFPTS>& yfunc)
    -> double {
  const double nnion = get_nnion(nonemptymgi, element, ion);
  const double ionpot_ev = collionrow.ionpot_ev;

  std::array<double, SFPTS> cross_section_vec{};
  get_xs_ionisation_vector(cross_section_vec, collionrow);

  const double y_dot_crosssection_de =
      std::inner_product(yfunc.begin(), yfunc.end(), cross_section_vec.begin(), 0.0) * DELTA_E;

  return nnion * ionpot_ev * y_dot_crosssection_de / E_init_ev;
}

// non-thermal ionisation rate coefficient (multiply by population to get rate)
auto nt_ionisation_ratecoeff_wfapprox(const int nonemptymgi, const int element, const int ion) -> double {
  const double deposition_rate_density = get_ntlepton_deposition_rate_density(nonemptymgi);
  // to get the non-thermal ionisation rate we need to divide the energy deposited
  // per unit volume per unit time in the grid cell (sum of terms above)
  // by the total ion number density and the "work per ion pair"
  return deposition_rate_density / get_nnion_tot(nonemptymgi) * get_oneoverw_approx_axelrod(element, ion, nonemptymgi);
}

// Integrate the ionisation cross section over the electron degradation function to get the ionisation rate
// coefficient i.e. multiply this by ion population to get a rate of ionisations per second Do not call during packet
// propagation, as the y vector may not be in memory! IMPORTANT: we are dividing by the shell potential, not the
// valence potential here! To change this set assumeshellpotentialisvalence to true
auto calculate_nt_ionisation_ratecoeff(const int nonemptymgi, const int element, const int ion,
                                       const bool assumeshellpotentialisvalence, const std::array<double, SFPTS>& yfunc)
    -> double {
  std::array<double, SFPTS> cross_section_vec{};
  std::array<double, SFPTS> cross_section_vec_allshells{};

  const int Z = get_atomicnumber(element);
  const int ionstage = get_ionstage(element, ion);
  double ionpot_valence = -1;

  for (const auto& collionrow : colliondata) {
    if (collionrow.Z == Z && collionrow.ionstage == ionstage) {
      get_xs_ionisation_vector(cross_section_vec, collionrow);

      if (assumeshellpotentialisvalence) {
        const double ionpot_shell = collionrow.ionpot_ev * EV;
        if (ionpot_valence < 0) {
          ionpot_valence = ionpot_shell;
        }

        // ensure that the first shell really was the valence shell (we assumed ascending energy order)
        assert_always(ionpot_shell >= ionpot_valence);

        // boost the ionisation rate by assuming shell vacancy energy is used to eject valence electrons
        for (int i = 0; i < SFPTS; i++) {
          cross_section_vec[i] *= ionpot_shell / ionpot_valence;
        }
      }

      for (int i = 0; i < SFPTS; i++) {
        cross_section_vec_allshells[i] += cross_section_vec[i];
      }
    }
  }

  const double y_xs_de =
      std::inner_product(yfunc.begin(), yfunc.end(), cross_section_vec_allshells.begin(), 0.0) * DELTA_E;

  const double deposition_rate_density_ev = get_ntlepton_deposition_rate_density(nonemptymgi) / EV;
  const double yscalefactor = deposition_rate_density_ev / E_init_ev;

  return yscalefactor * y_xs_de;
}

// Kozma & Fransson 1992 equation 12, except modified to be a sum over all shells of an ion.
// the result is in [erg]
void calculate_eff_ionpot_auger_rates(const int nonemptymgi, const int element, const int ion,
                                      const std::array<double, SFPTS>& yfunc) {
  const int Z = get_atomicnumber(element);
  const int ionstage = get_ionstage(element, ion);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  const double nnion = get_nnion(nonemptymgi, element, ion);  // ions/cm^3
  const double tot_nion = get_nnion_tot(nonemptymgi);
  const double X_ion = nnion / tot_nion;  // number fraction of this ion

  // The ionisation rates of all shells of an ion add to make the ion's total ionisation rate,
  // i.e., Gamma_ion = Gamma_shell_a + Gamma_shell_b + ...
  // And since the ionisation rate is inversely proportional to the effective ion potential, we solve:
  // (eta_ion / ionpot_ion) = (eta_shell_a / ionpot_shell_a) + (eta_shell_b / ionpot_shell_b) + ...
  // where eta is the fraction of the deposition energy going into ionisation of the ion or shell

  std::array<double, NT_MAX_AUGER_ELECTRONS + 1> eta_nauger_ionise_over_ionpot_sum{};
  std::array<double, NT_MAX_AUGER_ELECTRONS + 1> eta_nauger_ionise_sum{};
  auto& celliondata = get_cell_allions_data(nonemptymgi)[uniqueionindex];
  std::ranges::fill(celliondata.prob_num_auger, 0.);
  std::ranges::fill(celliondata.ionenfrac_num_auger, 0.);

  double eta_over_ionpot_sum = 0.;
  double eta_sum = 0.;
  double ionpot_valence = -1;
  int matching_nlsubshell_count = 0;
  for (const auto& collionrow : colliondata) {
    if (collionrow.Z == Z && collionrow.ionstage == ionstage) {
      matching_nlsubshell_count++;
      const double frac_ionisation_shell =
          calculate_nt_frac_ionisation_shell(nonemptymgi, element, ion, collionrow, yfunc);
      eta_sum += frac_ionisation_shell;
      const double ionpot_shell = collionrow.ionpot_ev * EV;

      if (ionpot_valence < 0) {
        ionpot_valence = ionpot_shell;
      }

      // ensure that the first shell really was the valence shell (we assumed ascending energy order)
      assert_always(ionpot_shell >= ionpot_valence);

      const double ionpot = NT_USE_VALENCE_IONPOTENTIAL ? ionpot_valence : ionpot_shell;
      const double eta_over_ionpot = frac_ionisation_shell / ionpot;  // this is proportional to rate

      eta_over_ionpot_sum += eta_over_ionpot;

      for (auto a = 0Z; a < std::ssize(eta_nauger_ionise_sum); a++) {
        eta_nauger_ionise_over_ionpot_sum[a] += eta_over_ionpot * collionrow.prob_num_auger[a];
        eta_nauger_ionise_sum[a] += frac_ionisation_shell * collionrow.prob_num_auger[a];
      }
    }
  }

  if (NT_MAX_AUGER_ELECTRONS > 0 && matching_nlsubshell_count > 0 && eta_over_ionpot_sum > 0.) {
    const int nions = get_nions(element);
    const int topion = nions - 1;
    if (ion < topion)  // don't try to ionise the top ion
    {
      for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
        const int upperion = ion + 1 + a;
        if (upperion <= topion)  // not too many Auger electrons to exceed the top ion of this element
        {
          celliondata.prob_num_auger[a] =
              static_cast<float>(eta_nauger_ionise_over_ionpot_sum[a] / eta_over_ionpot_sum);
          celliondata.ionenfrac_num_auger[a] = static_cast<float>(eta_nauger_ionise_sum[a] / eta_sum);
        } else {
          // the following ensures that multiple ionisations can't send you to an ion stage that is not in
          // the model. Send it to the highest ion stage instead
          const int a_replace = topion - ion - 1;
          celliondata.prob_num_auger.at(a_replace) = static_cast<float>(
              celliondata.prob_num_auger.at(a_replace) + (eta_nauger_ionise_over_ionpot_sum[a] / eta_over_ionpot_sum));
          celliondata.ionenfrac_num_auger.at(a_replace) =
              static_cast<float>(celliondata.ionenfrac_num_auger.at(a_replace) + (eta_nauger_ionise_sum[a] / eta_sum));

          celliondata.prob_num_auger[a] = 0;
          celliondata.ionenfrac_num_auger[a] = 0.;
        }
      }
    } else {
      // the top ion cannot be ionised further; keep the documented invariant that the
      // probabilities sum to one (matching zero_all_effionpot())
      celliondata.prob_num_auger[0] = 1.;
      celliondata.ionenfrac_num_auger[0] = 1.;
    }
  } else {
    const int a = 0;
    celliondata.prob_num_auger[a] = 1.;
    celliondata.ionenfrac_num_auger[a] = 1.;
  }

  if (matching_nlsubshell_count > 0) {
    double eff_ionpot = X_ion / eta_over_ionpot_sum;
    if (!std::isfinite(eff_ionpot)) {
      eff_ionpot = 0.;
    }
    celliondata.eff_ionpot = static_cast<float>(eff_ionpot);
  } else {
    // the absence of matching subshell data is reported once at startup by read_collion_data()
    celliondata.eff_ionpot = static_cast<float>(1. / get_oneoverw_approx_axelrod(element, ion, nonemptymgi));
  }
}

// get the effective ion potential from the stored value
// a value of 0. should be treated as invalid
auto get_eff_ionpot(const int nonemptymgi, const int element, const int ion) {
  return get_cell_allions_data(nonemptymgi)[get_uniqueionindex(element, ion)].eff_ionpot;
}

// Kozma & Fransson 1992 equation 13
// returns the rate coefficient in s^-1
auto nt_ionisation_ratecoeff_sf(const int nonemptymgi, const int element, const int ion) -> double {
  const double deposition_rate_density = get_ntlepton_deposition_rate_density(nonemptymgi);
  if (deposition_rate_density > 0.) {
    return deposition_rate_density / get_nnion_tot(nonemptymgi) / get_eff_ionpot(nonemptymgi, element, ion);
  }

  return 0.;
}

// vector of collisional excitation cross sections in cm^2
// epsilon_trans is in erg
// returns the index of the first valid cross section point (en >= epsilon_trans)
// all elements below this index are invalid and should not be used
auto get_xs_excitation_vector(const int alltransindex, const double statweight_lower, const double epsilon_trans)
    -> std::tuple<std::array<double, SFPTS>, int> {
  std::array<double, SFPTS> xs_excitation_vec{};
  if (epsilon_trans / EV > SF_EMAX) {
    // the excitation threshold is above the top of the energy grid, so the cross section is zero at
    // every grid point. Without this, the clamped startindex below would give a nonzero cross
    // section at the top grid point, which is below the excitation threshold.
    return {xs_excitation_vec, -1};
  }
  if (globals::alltrans.coll_str[alltransindex] >= 0) {
    // collision strength is available, so use it
    // Li et al. 2012 equation 11: sigma = pi * a_0^2 * (I_H / E) * Omega / g_lower,
    // with k_i^2 = E / I_H in units of the inverse Bohr radius squared
    const double constantfactor =
        H_ionpot / statweight_lower * globals::alltrans.coll_str[alltransindex] * PI * A_naught_squared;

    const int en_startindex = get_energyindex_ev_gteq(epsilon_trans / EV);

    std::fill_n(xs_excitation_vec.begin(), en_startindex, 0.);

    for (int j = en_startindex; j < SFPTS; j++) {
      const double energy = engrid(j) * EV;
      xs_excitation_vec[j] = constantfactor / energy;
    }
    return {xs_excitation_vec, en_startindex};
  }
  if (!globals::alltrans.forbidden[alltransindex]) {
    const double trans_osc_strength = globals::alltrans.osc_strength[alltransindex];
    // permitted E1 electric dipole transitions

    constexpr double A = 0.28;
    constexpr double B = 0.15;

    constexpr double prefactor = 45.585750051;  // 8 * pi^2/sqrt(3)
    const double epsilon_trans_ev = epsilon_trans / EV;

    // Eq 4 of Mewe 1972, possibly from Seaton 1962?
    const double constantfactor =
        epsilon_trans_ev * prefactor * A_naught_squared * pow2(H_ionpot / epsilon_trans) * trans_osc_strength;

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
auto calculate_nt_excitation_ratecoeff_perdeposition(const std::array<double, SFPTS>& yvec, const int alltransindex,
                                                     const double statweight_lower, const double epsilon_trans)
    -> double {
  const auto [xs_excitation_vec, xsstartindex] =
      get_xs_excitation_vector(alltransindex, statweight_lower, epsilon_trans);

  if (xsstartindex >= 0) {
    double y_xs_de = 0.;
    for (int i = xsstartindex; i < SFPTS; i++) {
      y_xs_de += yvec[i] * xs_excitation_vec[i];
    }
    // multiply by DELTA_E to get the integral over the energy grid
    y_xs_de *= DELTA_E;

    return y_xs_de / E_init_ev / EV;
  }

  return 0.;
}

// returns the energy rate [erg/cm3/s] going toward non-thermal ionisation of lowerion
auto ion_ntion_energyrate(const int nonemptymgi, const int element, const int lowerion) -> double {
  const double nnlowerion = get_nnion(nonemptymgi, element, lowerion);
  double enrate = 0.;
  const auto maxupperion = nt_ionisation_maxupperion(element, lowerion);
  for (int upperion = lowerion + 1; upperion <= maxupperion; upperion++) {
    const double upperionprobfrac = nt_ionisation_upperion_probability(nonemptymgi, element, lowerion, upperion, false);
    const double epsilon_trans = epsilon(element, upperion, 0) - epsilon(element, lowerion, 0);
    enrate += nnlowerion * upperionprobfrac * epsilon_trans;
  }

  const double gamma_nt = nt_ionisation_ratecoeff(nonemptymgi, element, lowerion);
  return gamma_nt * enrate;
}

// returns the energy rate [erg/s] going toward non-thermal ionisation in a modelgrid cell
auto get_ntion_energyrate(const int nonemptymgi) -> double {
  double ratetotal = 0.;
  for (int ielement = 0; ielement < get_nelements(); ielement++) {
    const int nions = get_nions(ielement);
    for (int ilowerion = 0; ilowerion < nions - 1; ilowerion++) {
      ratetotal += ion_ntion_energyrate(nonemptymgi, ielement, ilowerion);
    }
  }
  return ratetotal;
}

// select an ion to ionise, weighted by each ion's non-thermal ionisation energy rate. Returns {-1, -1}
// if no ion can be selected because the total rate is zero (the caller then deposits the packet as heat).
auto select_nt_ionisation(const int nonemptymgi, rngstate_type& rngstate) -> std::tuple<int, int> {
  const double ratetotal = get_ntion_energyrate(nonemptymgi);

  if (!(ratetotal > 0.)) {
    // No ion has a non-zero non-thermal ionisation energy rate. This happens when the cell's deposition
    // rate density is zero (nt_ionisation_ratecoeff is then zero everywhere), e.g. a cell reached by
    // energy for the first time this timestep, whose stored deposition rate is still from the previous
    // (dark) timestep. Signal that no ion can be selected.
    return {-1, -1};
  }

  const double zrand = rng_uniform(rngstate);

  // select based on the calculated energy going to ionisation for each ion
  double ratesum = 0.;
  for (int ielement = 0; ielement < get_nelements(); ielement++) {
    const int nions = get_nions(ielement);
    for (int ilowerion = 0; ilowerion < nions - 1; ilowerion++) {
      ratesum += ion_ntion_energyrate(nonemptymgi, ielement, ilowerion);
      // Strict comparison ensures that zrand == 0 skips leading ions with zero rate.
      if (ratesum > zrand * ratetotal) {
        return {ielement, ilowerion};
      }
    }
  }
  assert_always(false);
  return {-1, -1};
}

void analyse_sf_solution(const int nonemptymgi, const int timestep, const std::array<double, SFPTS>& yfunc,
                         const bool verbose) {
  const auto nne = grid::get_nne(nonemptymgi);
  const auto nntot = get_nnion_tot(nonemptymgi);
  const auto nnetot = grid::get_nnetot(nonemptymgi);

  double frac_excitation_total = 0.;
  double frac_ionisation_total = 0.;

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
      const int ioncharge = ionstage - 1;
      const int nbound = Z - ioncharge;  // number of bound electrons
      if (nbound <= 0) {
        continue;
      }
      const double nnion = get_nnion(nonemptymgi, element, ion);

      if (nnion <= 0.) {  // skip zero-abundance ions
        continue;
      }

      double frac_ionisation_ion = 0.;
      double frac_excitation_ion = 0.;
      if (verbose) {
        printlnlog("  Z={} ionstage {}:", Z, ionstage);
        printlnlog("    nnion/nntot: {:g}", nnion / nntot);
      }

      calculate_eff_ionpot_auger_rates(nonemptymgi, element, ion, yfunc);

      int matching_subshell_count = 0;
      for (const auto& collionrow : colliondata) {
        if (collionrow.Z != Z || collionrow.ionstage != ionstage) {
          continue;
        }
        const double frac_ionisation_ion_shell =
            calculate_nt_frac_ionisation_shell(nonemptymgi, element, ion, collionrow, yfunc);
        // with SF_AUGER_CONTRIBUTION_ON, the mean Auger energy per ionisation is re-injected into the
        // electron pool and gets counted in the heating/excitation fractions, so only the net energy
        // removed per ionisation (shell potential minus Auger energy) counts as ionisation here
        const double frac_auger_recycled =
            SF_AUGER_CONTRIBUTION_ON ? frac_ionisation_ion_shell * collionrow.en_auger_ev / collionrow.ionpot_ev : 0.;
        frac_ionisation_ion += frac_ionisation_ion_shell - frac_auger_recycled;
        matching_subshell_count++;

        if (verbose) {
          printlog("      shell ");
          if (collionrow.n >= 0) {
            printlog("n {}, l {}", collionrow.n, collionrow.l);
          } else {
            printlog("{} (Lotz)", shellnames.at(-collionrow.l));
          }
          printlog(" I {:5.1f} [eV]: frac_ionisation {:10.4e}", collionrow.ionpot_ev, frac_ionisation_ion_shell);

          if (NT_MAX_AUGER_ELECTRONS > 0) {
            printlog("  prob(n Auger elec):");
            for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
              printlog(" {}: {:.2f}", a, collionrow.prob_num_auger[a]);
            }
          }
          printlnlog("");
        }
      }

      // do not ionise the top ion
      if (ion < nions - 1) {
        get_cell_allions_data(nonemptymgi)[uniqueionindex].fracdep_ionisation_ion = frac_ionisation_ion;

        frac_ionisation_total += frac_ionisation_ion;
      } else {
        get_cell_allions_data(nonemptymgi)[uniqueionindex].fracdep_ionisation_ion = 0.;
      }
      if (verbose) {
        printlnlog("    frac_ionisation: {:g} ({} subshells)", frac_ionisation_ion, matching_subshell_count);
      }

      // excitation from all levels is expensive, so we limit it to a maximum number of levels
      const int nlevels = std::min(get_nlevels(element, ion), NTEXCITATION_MAXNLEVELS_LOWER);
      const bool above_minionfraction = (nnion >= MIN_ION_OVER_NNTOT * get_nnion_tot(nonemptymgi));

      // mark the current size so that this ion's entries can be removed again if its
      // frac_excitation_ion total turns out to be invalid below
      const auto tmp_excitation_list_size_before_ion = std::ssize(tmp_excitation_list);

      for (int lower = 0; lower < nlevels; lower++) {
        const auto uniquelevelindex = get_uniquelevelindex(element, ion, lower);
        const double statweight_lower = stat_weight(uniquelevelindex);
        const auto alltrans_startup = get_alltrans_startup(uniquelevelindex);
        const int nuptrans = get_nuptrans(uniquelevelindex);
        const double nnlevel = calculate_levelpop(nonemptymgi, element, ion, lower);
        const double epsilon_lower = epsilon(uniquelevelindex);

        for (int t = 0; t < nuptrans; t++) {
          const int alltransindex = alltrans_startup + t;
          const int upper = globals::alltrans.targetlevelindex[alltransindex];
          if (upper >= NTEXCITATION_MAXNLEVELS_UPPER) {
            continue;
          }

          const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
          if (epsilon_trans / EV < SF_EMIN) {
            // sfmatrix_add_excitation() does not include these transitions in the solution (their
            // energy sink is below the grid), so exclude them from the excitation fractions and the
            // stored excitation list too. Otherwise they would contribute to frac_excitation (making
            // frac_sum differ from 1.0) and receive packet excitation events and NLTE rates that the
            // Spencer-Fano solution never accounted for.
            continue;
          }
          const double ratecoeffperdeposition =
              calculate_nt_excitation_ratecoeff_perdeposition(yfunc, alltransindex, statweight_lower, epsilon_trans);
          const double frac_excitation_thistrans = nnlevel * epsilon_trans * ratecoeffperdeposition;
          frac_excitation_ion += frac_excitation_thistrans;

          if constexpr (NT_EXCITATION_ON) {
            assert_always(std::isfinite(ratecoeffperdeposition));
            // the atomic data set was limited for Fe V, which caused the ground multiplet to be massively
            // depleted, and then almost no recombination happened!
            if (above_minionfraction && ratecoeffperdeposition > 0 && (Z != 26 || ionstage != 5)) {
              tmp_excitation_list.push_back({
                  .frac_deposition = frac_excitation_thistrans,
                  .ratecoeffperdeposition = ratecoeffperdeposition,
                  .alltransindex = alltransindex,
              });
            }
          }  // NT_EXCITATION_ON
        }  // for t
      }  // for lower

      if (verbose) {
        printlnlog("    frac_excitation: {:g}", frac_excitation_ion);
      }
      if (frac_excitation_ion > 1. || !std::isfinite(frac_excitation_ion)) {
        printlnlog(
            "      [warning] cell {} ts {}: Z={} ionstage {}: invalid frac_excitation {:g}. Replacing with zero and "
            "dropping {} stored excitations for this ion",
            grid::get_mgi_of_nonemptymgi(nonemptymgi), timestep, Z, ionstage, frac_excitation_ion,
            std::ssize(tmp_excitation_list) - tmp_excitation_list_size_before_ion);
        frac_excitation_ion = 0.;
        // remove this ion's transitions from the excitation list. Otherwise their (invalid, possibly
        // huge) frac_deposition and ratecoeffperdeposition values would remain in the stored list and
        // distort the packet excitation-channel selection and the NLTE non-thermal excitation rates,
        // even though the ion's contribution has been excluded from frac_excitation_total.
        tmp_excitation_list.resize(tmp_excitation_list_size_before_ion);
      }
      frac_excitation_total += frac_excitation_ion;

      if (verbose) {
        printlnlog("    approxworkfn: {:9.2f} [eV]  (without using the Spencer-Fano solution)",
                   (1. / get_oneoverw_approx_axelrod(element, ion, nonemptymgi)) / EV);
        printlnlog("    eff_ionpot:   {:9.2f} [eV]  (always use valence potential is {})",
                   get_eff_ionpot(nonemptymgi, element, ion) / EV, (NT_USE_VALENCE_IONPOTENTIAL ? "true" : "false"));

        printlnlog("    approxworkfn Gamma:      {:9.3e}", nt_ionisation_ratecoeff_wfapprox(nonemptymgi, element, ion));

        if constexpr (NT_USE_VALENCE_IONPOTENTIAL) {
          printlnlog("    SF integral Gamma:       {:9.3e} (alternative if NT_USE_VALENCE_IONPOTENTIAL was disabled)",
                     calculate_nt_ionisation_ratecoeff(nonemptymgi, element, ion, false, yfunc));
        } else {
          printlnlog("    SF integral(I=Iv) Gamma: {:9.3e}  (alternative if NT_USE_VALENCE_IONPOTENTIAL was enabled)",
                     calculate_nt_ionisation_ratecoeff(nonemptymgi, element, ion, true, yfunc));
        }

        printlnlog("    ARTIS using Gamma:       {:9.3e}", nt_ionisation_ratecoeff(nonemptymgi, element, ion));
      }

      // the ion values (unlike shell ones) have been collapsed down to ensure that upperion < nions
      if (ion < nions - 1) {
        if (verbose) {
          printlog("    probability to ionstage:");
        }
        double prob_sum = 0.;
        for (int upperion = ion + 1; upperion <= nt_ionisation_maxupperion(element, ion); upperion++) {
          const double probability = nt_ionisation_upperion_probability(nonemptymgi, element, ion, upperion, false);
          prob_sum += probability;
          if (probability > 0. && verbose) {
            printlog(" {}: {:.3f}", get_ionstage(element, upperion), probability);
          }
        }
        assert_always((fabs(prob_sum - 1.0) <= 1e-2) ||
                      (nt_ionisation_ratecoeff_sf(nonemptymgi, element, ion) < 1e-20));

        if (verbose) {
          printlnlog("");
          printlog("         enfrac to ionstage:");
        }
        double enfrac_sum = 0.;
        for (int upperion = ion + 1; upperion <= nt_ionisation_maxupperion(element, ion); upperion++) {
          const double probability = nt_ionisation_upperion_probability(nonemptymgi, element, ion, upperion, true);
          enfrac_sum += probability;
          if (probability > 0. && verbose) {
            printlog(" {}: {:.3f}", get_ionstage(element, upperion), probability);
          }
        }
        assert_always(fabs(enfrac_sum - 1.0) <= 1e-2 ||
                      (nt_ionisation_ratecoeff_sf(nonemptymgi, element, ion) < 1e-20));
        if (verbose) {
          printlnlog("");
        }
      }
    }
  }

  if (nt_excitations_stored > 0) {
    // sort by descending frac_deposition
    std::ranges::SORT_OR_STABLE_SORT(tmp_excitation_list, std::ranges::greater{},
                                     &NonThermalExcitation::frac_deposition);

    // the excitation list is now sorted by frac_deposition descending
    const double deposition_rate_density = get_ntlepton_deposition_rate_density(nonemptymgi);

    if (std::ssize(tmp_excitation_list) > nt_excitations_stored) {
      // truncate the sorted list to save memory
      printlnlog("  Truncating non-thermal excitation list from {} to {} transitions.", tmp_excitation_list.size(),
                 nt_excitations_stored);
      tmp_excitation_list.resize(nt_excitations_stored);
    }

    nt_solution[nonemptymgi].frac_excitations_list_size = static_cast<int>(std::ssize(tmp_excitation_list));
    std::ranges::copy(tmp_excitation_list, get_cell_ntexcitations(nonemptymgi).begin());

    const auto T_e = grid::Te_allcells[nonemptymgi];
    if (verbose) {
      printlnlog("  Top non-thermal excitation fractions (total excitations = {}):",
                 nt_solution[nonemptymgi].frac_excitations_list_size);
      const int ntransdisplayed = std::min(10, nt_solution[nonemptymgi].frac_excitations_list_size);

      for (int excitationindex = 0; excitationindex < ntransdisplayed; excitationindex++) {
        const auto& ntexc = tmp_excitation_list[excitationindex];
        if (ntexc.frac_deposition > 0.) {
          const auto alltransindex = ntexc.alltransindex;
          const auto lineindex = globals::alltrans.lineindex[alltransindex];
          const int element = globals::linelist.elementindex[lineindex];
          const int ion = globals::linelist.ionindex[lineindex];
          const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
          const auto lower_uniquelevelindex = globals::linelist.uniquelevelindex_lower[lineindex];
          const auto upper_uniquelevelindex = globals::linelist.uniquelevelindex_upper[lineindex];
          const int lower = lower_uniquelevelindex - ionuniquelevelindexstart;
          const int upper = upper_uniquelevelindex - ionuniquelevelindexstart;
          const auto nnlevel_lower = calculate_levelpop(nonemptymgi, element, ion, lower);
          const auto nnlevel_upper = calculate_levelpop(nonemptymgi, element, ion, upper);
          const auto statweight_lower = stat_weight(lower_uniquelevelindex);

          const double epsilon_trans = epsilon(upper_uniquelevelindex) - epsilon(lower_uniquelevelindex);

          const double ntcollexc_ratecoeff = ntexc.ratecoeffperdeposition * deposition_rate_density;
          const auto statweight_upper = stat_weight(upper_uniquelevelindex);

          const double t_mid = globals::timesteps[timestep].mid;
          const double radexc_ratecoeff = rad_excitation_ratecoeff(
              nonemptymgi, statweight_upper, globals::alltrans.einstein_A[alltransindex], epsilon_trans, nnlevel_lower,
              nnlevel_upper, statweight_lower, alltransindex, t_mid);

          const double collexc_ratecoeff =
              col_excitation_ratecoeff(T_e, grid::get_clumpfactor(nonemptymgi) * nne, statweight_upper, alltransindex,
                                       epsilon_trans, statweight_lower);

          const double exc_ratecoeff = radexc_ratecoeff + collexc_ratecoeff + ntcollexc_ratecoeff;
          const auto coll_str = globals::alltrans.coll_str[alltransindex];

          printlnlog(
              "    frac_deposition {:.3e} Z={:2} ionstage {} lower {:4} upper {:4} rad_exc {:.1e} coll_exc {:.1e} "
              "nt_exc {:.1e} nt/tot {:.1e} collstr {:.1e} lineindex {}",
              ntexc.frac_deposition, get_atomicnumber(element), get_ionstage(element, ion), lower, upper,
              radexc_ratecoeff, collexc_ratecoeff, ntcollexc_ratecoeff, ntcollexc_ratecoeff / exc_ratecoeff, coll_str,
              lineindex);
        }
      }
    }

    // sort the excitation list by ascending alltransindex for fast lookup with a binary search
    std::ranges::SORT_OR_STABLE_SORT(get_cell_ntexcitations(nonemptymgi), std::ranges::less{},
                                     &NonThermalExcitation::alltransindex);

  }  // NT_EXCITATION_ON

  // calculate number density of non-thermal electrons
  const double deposition_rate_density_ev = get_ntlepton_deposition_rate_density(nonemptymgi) / EV;
  const double yscalefactor = deposition_rate_density_ev / E_init_ev;

  double nne_nt_max = 0.;
  for (int i = 0; i < SFPTS; i++) {
    const double endash = engrid(i);
    const double delta_endash = DELTA_E;
    // 1 / v for a non-relativistic electron of kinetic energy endash [eV], in sec/cm
    const double oneovervelocity = std::sqrt(ME / (2 * endash * EV));
    nne_nt_max += yscalefactor * yfunc[i] * oneovervelocity * delta_endash;
  }

  nt_solution[nonemptymgi].frac_excitation = static_cast<float>(frac_excitation_total);
  nt_solution[nonemptymgi].frac_ionisation = static_cast<float>(frac_ionisation_total);
  const auto frac_heating_calculated = calculate_frac_heating(nonemptymgi, yfunc);
  const double frac_sum = frac_heating_calculated + frac_excitation_total + frac_ionisation_total;

  if (verbose) {
    printlnlog("  deposition:  {:9.2f} [eV/s/cm^3]", deposition_rate_density_ev);
    printlnlog("  nne:         {:9.3e} [e-/cm^3]", nne);
    printlnlog("  nnetot:      {:9.3e} [e-/cm^3]", nnetot);
    printlnlog("  nne_nt     < {:9.3e} [e-/cm^3]", nne_nt_max);
    printlnlog("  nne_nt/nne < {:9.3e}", nne_nt_max / nne);

    printlnlog("  frac_heating_tot:    {:g}", frac_heating_calculated);
    printlnlog("  frac_excitation_tot: {:g}", frac_excitation_total);
    printlnlog("  frac_ionisation_tot: {:g}", frac_ionisation_total);
    printlnlog("  frac_sum:            {:g} (should be close to 1.0)", frac_sum);
  }

  // force frac_sum to be 1.0 by adjusting frac_heating. Clamp to [0, 1] so a bad solution (where
  // excitation + ionisation exceed 1) cannot feed a negative heating fraction into the T_e solver.
  // When non-thermal excitation is not being modelled, do_ntlepton_deposit() sends the excitation
  // share to k-packets instead, so it has to be counted as heating here for the energy that the
  // packets carry to match the heating rate that the T_e solver is given.
  const double frac_excitation_nonheating = NT_EXCITATION_ON ? frac_excitation_total : 0.;
  nt_solution[nonemptymgi].frac_heating =
      static_cast<float>(std::clamp(1. - frac_excitation_nonheating - frac_ionisation_total, 0., 1.));

  if (!ftol<0.02>(frac_sum, 1.0)) {
    printlnlog("[warning] frac_sum is {:g}, but should be 1.0", frac_sum);
    printlnlog("  (replacing calculated frac_heating_tot {:g} with {:g} to make frac_sum = 1.0)",
               frac_heating_calculated, nt_solution[nonemptymgi].frac_heating);
  }
}

void sfmatrix_add_excitation(std::span<double> sfmatrixuppertri, const int nonemptymgi, const int element,
                             const int ion) {
  // excitation terms

  const int nlevels_all = get_nlevels(element, ion);
  const int nlevels = (nlevels_all > NTEXCITATION_MAXNLEVELS_LOWER) ? NTEXCITATION_MAXNLEVELS_LOWER : nlevels_all;

  const auto lowers = std::ranges::iota_view{0, nlevels};
  std::for_each(lowers.begin(), lowers.end(), [&](const int lower) {
    const auto uniquelevelindex = get_uniquelevelindex(element, ion, lower);
    const double statweight_lower = stat_weight(uniquelevelindex);
    const double nnlevel = calculate_levelpop(nonemptymgi, element, ion, lower);
    const double epsilon_lower = epsilon(uniquelevelindex);
    const auto alltrans_startup = get_alltrans_startup(uniquelevelindex);
    const int nuptrans = get_nuptrans(uniquelevelindex);
    for (int t = 0; t < nuptrans; t++) {
      const int alltransindex = alltrans_startup + t;
      const int upper = globals::alltrans.targetlevelindex[alltransindex];
      if (upper >= NTEXCITATION_MAXNLEVELS_UPPER) {
        continue;
      }
      const double epsilon_trans = epsilon(element, ion, upper) - epsilon_lower;
      const double epsilon_trans_ev = epsilon_trans / EV;
      if (epsilon_trans_ev < SF_EMIN) {
        continue;
      }

      auto [vec_xs_excitation_deltae, xsstartindex] =
          get_xs_excitation_vector(alltransindex, statweight_lower, epsilon_trans);
      if (xsstartindex >= 0) {
        for (int j = xsstartindex; j < SFPTS; j++) {
          vec_xs_excitation_deltae[j] *= DELTA_E;
        }

        for (int i = 0; i < SFPTS; i++) {
          const int rowoffset = uppertriangular(i, 0);
          const double en = engrid(i);
          const int stopindex = get_energyindex_ev_lteq(en + epsilon_trans_ev);

          const int startindex = std::max(i, xsstartindex);
          for (int j = startindex; j < stopindex; j++) {
            atomicadd(sfmatrixuppertri[rowoffset + j], nnlevel * vec_xs_excitation_deltae[j]);
          }

          // do the last bit separately because we're not using the full delta_e interval.
          // clamp to DELTA_E: when en + epsilon_trans_ev exceeds SF_EMAX, get_energyindex_ev_lteq()
          // clamps stopindex to SFPTS-1, so the raw remainder can exceed DELTA_E and would over-weight
          // the top energy point. Capping it gives that bin a full (not inflated) contribution.
          const double delta_en_actual = std::min(en + epsilon_trans_ev - engrid(stopindex), DELTA_E);
          atomicadd(sfmatrixuppertri[rowoffset + stopindex],
                    nnlevel * vec_xs_excitation_deltae[stopindex] * delta_en_actual / DELTA_E);
        }
      }
    }
  });
}

void sfmatrix_add_ionisation(std::span<double> sfmatrixuppertri, const int Z, const int ionstage, const double nnion)
// add the ionisation terms to the Spencer-Fano matrix
{
  std::array<double, SFPTS> vec_xs_ionisation{};
  for (const auto& collionrow : colliondata) {
    if (collionrow.Z == Z && collionrow.ionstage == ionstage) {
      const double ionpot_ev = collionrow.ionpot_ev;
      const double en_auger_ev = collionrow.en_auger_ev;
      const double J = get_J(Z, ionstage, ionpot_ev);

      assert_always(ionpot_ev >= SF_EMIN);

      const int xsstartindex = get_xs_ionisation_vector(vec_xs_ionisation, collionrow);
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
        prefactors[j] = vec_xs_ionisation[j] * nnion / std::atan((endash - ionpot_ev) / 2 / J);
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
        if ((2 * en) + ionpot_ev <= SF_EMAX) {
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

      // shells with no Auger data have en_auger_ev == 0 (and prob_num_auger[0] == 1, which would make the
      // energy boost factor below infinite) and inject no Auger electrons
      if (SF_AUGER_CONTRIBUTION_ON && en_auger_ev > 0.) {
        int augerstopindex = 0;
        if constexpr (SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN) {
          // en_auger_ev is (if LJS understands it correctly) averaged to include some probability of zero Auger
          // electrons so we need a boost to get the average energy of Auger electrons given that there are one or more
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
            // the Auger electron source is an integral of the ionisation rate over E', so each column
            // carries a DELTA_E integration weight like the other ionisation integrals above
            const double xs = vec_xs_ionisation[j];
            if constexpr (SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN) {
              const double en_boost = 1 / (1. - collionrow.prob_num_auger[0]);
              for (int a = 1; a <= NT_MAX_AUGER_ELECTRONS; a++) {
                if (en < (en_auger_ev * en_boost / a)) {
                  sfmatrixuppertri[rowoffset + j] -= nnion * xs * collionrow.prob_num_auger[a] * a * DELTA_E;
                }
              }
            } else {
              // inject a single electron per ionisation carrying the mean total Auger energy
              assert_always(en < en_auger_ev);
              sfmatrixuppertri[rowoffset + j] -= nnion * xs * DELTA_E;
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
auto sfmatrix_solve(const std::span<const double> sfmatrix) -> std::array<double, SFPTS> {
  // solve the matrix-vector equation sfmatrix * yvec = rhsvec for yvec

  THREADLOCALONHOST std::array<double, SFPTS> yvec_arr{};

  // The solve below reads only the upper triangle (via triangularView<Upper>()) while the residual is formed from the
  // full matrix, so anything left below the diagonal would make the refinement iterate on a different system than the
  // one that was solved.
  assert_testmodeonly([sfmatrix] {
    for (int i = 1; i < SFPTS; i++) {
      for (int j = 0; j < i; j++) {
        if (sfmatrix[(i * SFPTS) + j] != 0.) {
          return false;
        }
      }
    }
    return true;
  }());

  const Eigen::Map<const Eigen::Vector<double, SFPTS>> eigen_rhsvec{rhsvec.data()};
  const Eigen::Map<const Eigen::Matrix<double, SFPTS, SFPTS, Eigen::RowMajor>> eigen_sfmatrix{sfmatrix.data()};

  const auto eigen_sfmatrix_upper = eigen_sfmatrix.triangularView<Eigen::Upper>();
  Eigen::Map<Eigen::Vector<double, SFPTS>> eigen_yvec{yvec_arr.data()};
  eigen_yvec = eigen_sfmatrix_upper.solve(eigen_rhsvec);

  // refine the solution

  THREADLOCALONHOST std::array<double, SFPTS> yvec_best{};
  THREADLOCALONHOST Eigen::Vector<double, SFPTS> eigen_residual_vec{};

  double error_best = -1.;
  for (int iteration = 0; iteration < 10; iteration++) {
    if (iteration > 0) {
      eigen_yvec += eigen_sfmatrix_upper.solve(eigen_residual_vec);
    }
    eigen_residual_vec = eigen_rhsvec - eigen_sfmatrix * eigen_yvec;

    // PropagateNaN is required for the non-finite test below: Eigen documents the default (PropagateFast) as
    // undefined in the presence of a NaN, and its scalar and SIMD paths genuinely differ there, so a NaN residual
    // could otherwise be scored as a finite value and accepted. For an all-finite residual the two agree exactly.
    const double error = eigen_residual_vec.cwiseAbs().maxCoeff<Eigen::PropagateNaN>();

    // only ever store a finite error, so that a non-finite iteration cannot latch error_best and block a later
    // iteration from being accepted. If every iteration is non-finite then error_best stays negative and the assertion
    // below fails.
    if (std::isfinite(error) && (error_best < 0. || error < error_best)) {
      std::ranges::copy(yvec_arr, yvec_best.begin());
      error_best = error;
    }
  }

  assert_always(error_best >= 0.);
  std::ranges::copy(yvec_best, yvec_arr.begin());

  if (error_best > 1e-10) {
    printlnlog("  SF solver iterative refinement: best solution vector has a max residual of {:g}", error_best);
  }

  assert_always(std::ranges::all_of(yvec_arr, [](double y) { return y >= 0.; }));

  return yvec_arr;
}

}  // anonymous namespace

void init() {
  const ptrdiff_t nonempty_npts_model = grid::get_nonempty_npts_model();

  ntlepton_deposition_rate_density_all_cells = MPI_shared_array<double>(nonempty_npts_model);

  if (globals::rank_in_node == 0) {
    std::ranges::fill(ntlepton_deposition_rate_density_all_cells, -1.);
  }

  if (!NT_ON) {
    return;
  }

  read_binding_energies();

  if (!NT_SOLVE_SPENCERFANO) {
    return;
  }

  printlnlog("Initializing non-thermal solver with:");
  printlnlog("  NT_EXCITATION {}", NT_EXCITATION_ON ? "on" : "off");
  printlnlog("  MAX_NT_EXCITATIONS_STORED {}", MAX_NT_EXCITATIONS_STORED);
  printlnlog("  NTEXCITATION_MAXNLEVELS_LOWER {}", NTEXCITATION_MAXNLEVELS_LOWER);
  printlnlog("  NTEXCITATION_MAXNLEVELS_UPPER {}", NTEXCITATION_MAXNLEVELS_UPPER);
  printlnlog("  SFPTS {}", SFPTS);
  printlnlog("  SF_EMIN {:g} [eV]", SF_EMIN);
  printlnlog("  SF_EMAX {:g} [eV]", SF_EMAX);
  printlnlog("  NT_USE_VALENCE_IONPOTENTIAL {}", NT_USE_VALENCE_IONPOTENTIAL ? "on" : "off");
  printlnlog("  NT_MAX_AUGER_ELECTRONS {}", NT_MAX_AUGER_ELECTRONS);
  printlnlog("  SF_AUGER_CONTRIBUTION {}", SF_AUGER_CONTRIBUTION_ON ? "on" : "off");
  printlnlog("  SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN {}", SF_AUGER_CONTRIBUTION_DISTRIBUTE_EN ? "on" : "off");

  if (NT_EXCITATION_ON) {
    nt_excitations_stored = std::min(MAX_NT_EXCITATIONS_STORED, get_possible_nt_excitation_count());
    printlnlog("[info] mem_usage: storing {} non-thermal excitations for non-empty cells occupies {:.3f} MB",
               nt_excitations_stored,
               nonempty_npts_model * sizeof(NonThermalExcitation) * nt_excitations_stored / 1024. / 1024.);

    excitations_list_all_cells = MPI_shared_array<NonThermalExcitation>(nonempty_npts_model * nt_excitations_stored);
  }

  ion_data_all_cells = MPI_shared_array<NonThermalSolutionIon>(nonempty_npts_model * get_includedions());

  nt_solution = MPI_shared_array<NonThermalCellSolution>(nonempty_npts_model);

  if (globals::rank_in_node == 0) {
    for (auto nonemptymgi = 0Z; nonemptymgi < nonempty_npts_model; nonemptymgi++) {
      // should make these negative?
      nt_solution[nonemptymgi].frac_heating = 0.97;
      nt_solution[nonemptymgi].frac_ionisation = 0.03;
      nt_solution[nonemptymgi].frac_excitation = 0.;

      nt_solution[nonemptymgi].nneperion_when_solved = -1.;
      nt_solution[nonemptymgi].timestep_last_solved = -1;

      zero_all_effionpot(nonemptymgi);

      nt_solution[nonemptymgi].frac_excitations_list_size = 0;
    }
  }
  MPI_Barrier_node();

  const double sourceintegral = std::ranges::fold_left(
      std::views::iota(0, SFPTS), 0.0, [](double sum, int s) { return sum + (sourcevec(s) * DELTA_E); });

  printlnlog("E_init: {:14.7e} [eV/s/cm^3]", E_init_ev);
  printlnlog("source function integral: {:14.7e}", sourceintegral);

  read_collion_data();

  printlnlog("Finished initializing non-thermal solver");
}

// set total non-thermal deposition rate from individual gamma/positron/electron/alpha rates. This should be called
// after packet propagation is finished for this timestep and normalise_deposition_estimators() has been done
void calculate_deposition_rate_density(const int nonemptymgi, HeatingCoolingRates& heatingcoolingrates,
                                       const decay::AnaEmissionPowerPerMass& emission_power_per_mass) {
  heatingcoolingrates.dep_gamma = globals::dep_estimator_gamma[nonemptymgi];

  const double rho = grid::get_rho(nonemptymgi);

  heatingcoolingrates.eps_gamma_ana =
      rho * decay::get_modelcell_decaypower_per_mass(nonemptymgi, emission_power_per_mass.gamma);

  heatingcoolingrates.eps_positron_ana =
      rho * decay::get_modelcell_decaypower_per_mass(nonemptymgi, emission_power_per_mass.positron);

  heatingcoolingrates.eps_electron_ana =
      rho * decay::get_modelcell_decaypower_per_mass(nonemptymgi, emission_power_per_mass.electron);

  heatingcoolingrates.eps_alpha_ana =
      rho * decay::get_modelcell_decaypower_per_mass(nonemptymgi, emission_power_per_mass.alpha);

  heatingcoolingrates.eps_spfission_ana =
      rho * decay::get_modelcell_decaypower_per_mass(nonemptymgi, emission_power_per_mass.spfission);

  if (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::INSTANTFULLDEPOSITION) {
    // for instant full deposition, the deposition rate is the same as the emission rate, which we know analytically
    // without Monte Carlo noise (although strictly, it should be an integral from the timestep start to the end divided
    // by timestep duration instead of the instantaneous rate at tmid)
    heatingcoolingrates.dep_positron = heatingcoolingrates.eps_positron_ana;
    heatingcoolingrates.dep_electron = heatingcoolingrates.eps_electron_ana;
    heatingcoolingrates.dep_alpha = heatingcoolingrates.eps_alpha_ana;
  } else {
    // with time-dependent deposition, we don't have an analytic rate, so we use the Monte Carlo rate
    heatingcoolingrates.dep_positron = globals::dep_estimator_positron[nonemptymgi];
    heatingcoolingrates.dep_electron = globals::dep_estimator_electron[nonemptymgi];
    heatingcoolingrates.dep_alpha = globals::dep_estimator_alpha[nonemptymgi];
  }

  // spontaneous fission contribution is always treated as instant deposition
  heatingcoolingrates.dep_spfission = heatingcoolingrates.eps_spfission_ana;

  ntlepton_deposition_rate_density_all_cells[nonemptymgi] =
      (heatingcoolingrates.dep_gamma + heatingcoolingrates.dep_positron + heatingcoolingrates.dep_electron);
}

// get the non-thermal lepton deposition rate density (gamma + positron + electron, excluding alpha and
// spontaneous fission) in erg / s / cm^3, previously stored by calculate_deposition_rate_density()
DEVICE_FUNC auto get_ntlepton_deposition_rate_density(const int nonemptymgi) -> double {
  assert_always(ntlepton_deposition_rate_density_all_cells[nonemptymgi] >= 0);
  return ntlepton_deposition_rate_density_all_cells[nonemptymgi];
}

auto get_nt_frac_heating(const int nonemptymgi) -> float {
  if (!NT_ON) {
    return 1.;
  }
  if (!NT_SOLVE_SPENCERFANO) {
    return 0.97;
  }
  const float frac_heating = nt_solution[nonemptymgi].frac_heating;
  return frac_heating;
}

DEVICE_FUNC auto nt_ionisation_upperion_probability(const int nonemptymgi, const int element, const int lowerion,
                                                    const int upperion, const bool energyweighted) -> double {
  assert_always(upperion > lowerion);
  assert_always(upperion < get_nions(element));
  assert_always(upperion <= nt_ionisation_maxupperion(element, lowerion));
  if (NT_SOLVE_SPENCERFANO && NT_MAX_AUGER_ELECTRONS > 0) {
    const int numaugerelec = upperion - lowerion - 1;  // number of Auger electrons to go from lowerin to upper ion
    const int uniqueionindex = get_uniqueionindex(element, lowerion);
    const auto& cell_ion_data = get_cell_allions_data(nonemptymgi)[uniqueionindex];
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
        assert_always(fabs(prob_remaining - cell_ion_data.prob_num_auger[numaugerelec]) < 0.001);
      }
      return prob_remaining;
    }

    return 0.;
  }
  return (upperion == lowerion + 1) ? 1.0 : 0.;
}

DEVICE_FUNC auto nt_ionisation_maxupperion(const int element, const int lowerion) -> int {
  const int nions = get_nions(element);
  assert_always(lowerion < nions - 1);
  int maxupper = lowerion + 1;

  if (NT_SOLVE_SPENCERFANO) {
    maxupper += NT_MAX_AUGER_ELECTRONS;
  }

  maxupper = std::min(maxupper, nions - 1);

  return maxupper;
}

DEVICE_FUNC auto nt_random_upperion(const int nonemptymgi, const int element, const int lowerion,
                                    const bool energyweighted, rngstate_type& rngstate) -> int {
  assert_testmodeonly(lowerion < get_nions(element) - 1);
  if (NT_SOLVE_SPENCERFANO && NT_MAX_AUGER_ELECTRONS > 0) {
    const double zrand = rng_uniform(rngstate);

    double prob_sum = 0.;
    for (int upperion = lowerion + 1; upperion <= nt_ionisation_maxupperion(element, lowerion); upperion++) {
      prob_sum += nt_ionisation_upperion_probability(nonemptymgi, element, lowerion, upperion, energyweighted);

      // Strict comparison ensures that zrand == 0 skips leading zero-probability outcomes.
      if (zrand < prob_sum) {
        return upperion;
      }
    }

    // the channel probabilities are stored as floats and the tabulated Auger yields are only
    // normalised to about one part in a thousand, so a draw very close to one can fall past the end
    // of the cumulative sum. Absorb that small deficit into the highest reachable ion rather than
    // failing the whole run, but still reject a distribution that is badly wrong rather than
    // silently returning the top ion for every draw.
    assert_always(prob_sum > 0.99);
    return nt_ionisation_maxupperion(element, lowerion);
  }
  return lowerion + 1;
}

DEVICE_FUNC auto nt_ionisation_ratecoeff(const int nonemptymgi, const int element, const int ion) -> double {
  assert_always(NT_ON);

  if (NT_SOLVE_SPENCERFANO) {
    const double Y_nt = nt_ionisation_ratecoeff_sf(nonemptymgi, element, ion);
    if (!std::isfinite(Y_nt)) {
      // probably because eff_ionpot = 0 because the solver hasn't been run yet, or no impact ionisation cross sections
      // exist
      const double Y_nt_wfapprox = nt_ionisation_ratecoeff_wfapprox(nonemptymgi, element, ion);
      return Y_nt_wfapprox;
    }
    assert_always(Y_nt >= 0.);

    return Y_nt;
  }
  return nt_ionisation_ratecoeff_wfapprox(nonemptymgi, element, ion);
}

DEVICE_FUNC auto nt_excitation_ratecoeff(const int nonemptymgi, const int lowerlevel, const int upperlevel,
                                         const int alltransindex) -> double {
  if constexpr (!NT_EXCITATION_ON) {
    return 0.;
  }
  if (lowerlevel >= NTEXCITATION_MAXNLEVELS_LOWER) {
    return 0.;
  }
  if (upperlevel >= NTEXCITATION_MAXNLEVELS_UPPER) {
    return 0.;
  }

  // binary search, assuming the excitation list is sorted by alltransindex ascending
  const auto ntexclist = get_cell_ntexcitations(nonemptymgi);
  const auto ntexcitation =
      std::ranges::lower_bound(ntexclist, alltransindex, {}, &NonThermalExcitation::alltransindex);
  if (ntexcitation == ntexclist.end() || ntexcitation->alltransindex != alltransindex) {
    return 0.;
  }

  const double deposition_rate_density = get_ntlepton_deposition_rate_density(nonemptymgi);
  const double ratecoeffperdeposition = ntexcitation->ratecoeffperdeposition;

  return ratecoeffperdeposition * deposition_rate_density;
}

DEVICE_FUNC void do_ntalpha_fisprod_deposit(Packet& pkt) {
  // if ionisation by alpha particles is found to be important for the ionisation state, we could do a separate
  // Spencer-Fano solution. For now, just treat alpha deposition as pure heating (even though the alpha deposition rate
  // was calculated from the sum of ionisation and plasma heating)
  atomicadd(nt_energy_deposited, pkt.e_cmf);
  pkt.type = TYPE_KPKT;
  stats::increment(stats::Counter::NT_STAT_TO_KPKT);
}

DEVICE_FUNC void do_ntlepton_deposit(Packet& pkt) {
  atomicadd(nt_energy_deposited, pkt.e_cmf);

  const int modelgridindex = grid::get_propcell_modelgridindex(pkt.cellindex);
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);

  // macroatom should not be activated in thick cells
  if (NT_ON && NT_SOLVE_SPENCERFANO && grid::thick_allcells[nonemptymgi] != 1) {
    // here there is some probability to cause ionisation or excitation to a macroatom packet
    // instead of converting directly to k-packet (unless the heating channel is selected)

    double zrand = rng_uniform(get_rngstate(pkt));
    // zrand is initially between [0, 1), but we will subtract off each component of the deposition fractions
    // until we end and select transition_ij when zrand < dep_frac_transition_ij

    // Gate on the stored ionisation fraction so that the k-packet probability below is exactly
    // 1 - frac_ionisation - frac_excitation = frac_heating, the deposition partition that
    // analyse_sf_solution() stores and the T_e solver applies. (The live get_ntion_energyrate() /
    // deposition estimate used previously is a different estimator and does not sum with frac_excitation
    // to give frac_heating.) The live per-ion rates still choose which ion is ionised in
    // select_nt_ionisation(), via a separate random draw that is renormalised internally, so the gate
    // and the ion selection do not need to share a normalisation.
    const double frac_ionisation = get_nt_frac_ionisation(nonemptymgi);

    if (zrand < frac_ionisation) {
      const auto [element, lowerion] = select_nt_ionisation(nonemptymgi, get_rngstate(pkt));

      if (lowerion >= 0) {
        const int upperion = nt_random_upperion(nonemptymgi, element, lowerion, true, get_rngstate(pkt));

        stats::increment(stats::Counter::MA_STAT_ACTIVATION_NTCOLLION);
        stats::increment(stats::Counter::INTERACTIONS);
        pkt.trueemissiontype = EMTYPE_NOTSET;
        pkt.trueem_pos = {NAN, NAN, NAN};

        stats::increment(stats::Counter::NT_STAT_TO_IONISATION);

        do_macroatom(pkt, {.element = element, .ion = upperion, .level = 0, .activatingline = -99});
        return;
      }

      // No ion could be selected because the deposition rate density is zero. Deposit the packet as heat
      // rather than falling through to the excitation block below, where the outer zrand (< frac_ionisation)
      // would be an invalid position in the excitation list.
      pkt.type = TYPE_KPKT;
      stats::increment(stats::Counter::NT_STAT_TO_KPKT);
      return;
    }

    // Route the excitation share of the deposition to macroatoms. Whatever is left over after the
    // ionisation and excitation channels becomes a k-packet (heating) below, so the k-packet
    // probability is 1 - frac_ionisation - frac_excitation, matching the frac_heating that
    // analyse_sf_solution() stores and the T_e solver applies to the deposition rate.
    const double frac_excitation = NT_EXCITATION_ON ? get_nt_frac_excitation(nonemptymgi) : 0.;
    if (zrand < (frac_ionisation + frac_excitation)) {
      zrand -= frac_ionisation;
      // now zrand is between zero and frac_excitation
      // the selection algorithm is the same as for the ionisation transitions
      for (const auto& ntexcitation : get_cell_ntexcitations(nonemptymgi)) {
        const double frac_deposition_exc = ntexcitation.frac_deposition;
        if (zrand < frac_deposition_exc) {
          const auto lineindex = globals::alltrans.lineindex[ntexcitation.alltransindex];
          const auto [element, ion, upper] =
              get_levelfromuniquelevelindex(globals::linelist.uniquelevelindex_upper[lineindex]);

          stats::increment(stats::Counter::MA_STAT_ACTIVATION_NTCOLLEXC);
          stats::increment(stats::Counter::INTERACTIONS);
          pkt.trueemissiontype = EMTYPE_NOTSET;
          pkt.trueem_pos = {NAN, NAN, NAN};

          stats::increment(stats::Counter::NT_STAT_TO_EXCITATION);

          do_macroatom(pkt, {.element = element, .ion = ion, .level = upper, .activatingline = -99});
          return;
        }
        zrand -= frac_deposition_exc;
      }
      // Reaching here means zrand landed in the part of frac_excitation that the stored list does not
      // cover: the list is truncated to MAX_NT_EXCITATIONS_STORED and excludes some transitions (e.g.
      // Fe V, and ions below MIN_ION_OVER_NNTOT), while frac_excitation counts every transition. That
      // unresolved remainder falls through to a k-packet, i.e. it is treated as heating.
    }
  }

  pkt.type = TYPE_KPKT;
  stats::increment(stats::Counter::NT_STAT_TO_KPKT);
}

// The matrix discretisation follows Equation (2) of Li et al. (2012); the Auger-electron extension is described by
// Shingles et al. (2020), Section 2.5, doi:10.1093/mnras/stz3412.
void solve_spencerfano(const int nonemptymgi, const int timestep, const int iteration) {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  bool skip_solution = false;
  if (timestep < globals::num_lte_timesteps + 1) {
    // this global condition is reported once per timestep by update_grid(), not once per cell
    skip_solution = true;
  } else if (get_ntlepton_deposition_rate_density(nonemptymgi) / EV < MINDEPRATE) {
    printlnlog(
        "Non-thermal deposition rate of {:g} [eV/s/cm^3] below MINDEPRATE {:g} [eV/s/cm^3] in cell {} at timestep {}. "
        "Skipping Spencer-Fano solution.",
        get_ntlepton_deposition_rate_density(nonemptymgi) / EV, MINDEPRATE, modelgridindex, timestep);

    skip_solution = true;
  }

  if (skip_solution) {
    // Axelrod values
    nt_solution[nonemptymgi].frac_heating = 0.97;
    nt_solution[nonemptymgi].frac_ionisation = 0.03;
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

  printlnlog(
      "Spencer-Fano solver at timestep {} (last solution was at timestep {}) nne/niontot = {:g}, at last solution was "
      "{:g} fracdiff {:g}",
      timestep, timestep_last_solved, nne_per_ion, nne_per_ion_last, nne_per_ion_fracdiff);

  if ((nne_per_ion_fracdiff < NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS) &&
      (timestep - timestep_last_solved <= SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS) &&
      timestep_last_solved > globals::num_lte_timesteps) {
    printlnlog(
        "Keeping Spencer-Fano solution from timestep {} because x_e fracdiff {:g} < {:g} and because timestep {} - {} "
        "<= {}",
        timestep_last_solved, nne_per_ion_fracdiff, NT_MAX_FRACDIFF_NNEPERION_BETWEEN_SOLUTIONS, timestep,
        timestep_last_solved, SF_MAX_TIMESTEPS_BETWEEN_SOLUTIONS);

    return;
  }
  printlnlog(
      "Setting up Spencer-Fano equation with {} energy points from {:g} [eV] to {:g} [eV] in cell {} at timestep {} "
      "iteration {} (nne={:g} [e-/cm^3])",
      SFPTS, SF_EMIN, SF_EMAX, modelgridindex, timestep, iteration, nne);

  nt_solution[nonemptymgi].nneperion_when_solved = static_cast<float>(nne_per_ion);
  nt_solution[nonemptymgi].timestep_last_solved = timestep;

  // sfmatrix will be a compacted upper triangular matrix during construction and then expanded into a full matrix (with
  // lots of zeros) just before the solver is called
  THREADLOCALONHOST std::vector<double> sfmatrix(SFPTS * SFPTS);
  // while we are filling the matrix, we only need to fill the upper triangular part in a compacted form (to reduce
  // cache misses). Later when we go to solve it, we will expand it back to a full square matrix with zeros in the lower
  // triangle
  std::fill_n(sfmatrix.begin(), SFPTS * (SFPTS + 1) / 2, 0.);
  // loss terms and source terms
  for (int i = 0; i < SFPTS; i++) {
    sfmatrix[uppertriangular(i, i)] += electron_loss_rate(engrid(i) * EV, nne) / EV;
  }

  for (int element = 0; element < get_nelements(); element++) {
    const int Z = get_atomicnumber(element);
    const int nions = get_nions(element);
    bool first_included_ion_of_element = true;
    for (int ion = 0; ion < nions; ion++) {
      const double nnion = get_nnion(nonemptymgi, element, ion);

      // skip negligible ions
      if (nnion < MIN_ION_OVER_NNTOT * get_nnion_tot(nonemptymgi)) {
        continue;
      }

      const int ionstage = get_ionstage(element, ion);
      if (first_included_ion_of_element) {
        printlog("  including Z={:2} ionstages: ", Z);
        for (int i = 1; i < get_ionstage(element, ion); i++) {
          printlog("  ");
        }
        first_included_ion_of_element = false;
      }

      printlog("{} ", ionstage);

      sfmatrix_add_excitation(sfmatrix, nonemptymgi, element, ion);

      if ((ion < nions - 1)) {
        sfmatrix_add_ionisation(sfmatrix, Z, ionstage, nnion);
      }
    }
    if (!first_included_ion_of_element) {
      printlnlog("");
    }
  }

  decompactify_triangular_matrix(sfmatrix);
  const auto yfunc = sfmatrix_solve(sfmatrix);
  constexpr bool verbose = false;
  analyse_sf_solution(nonemptymgi, timestep, yfunc, verbose);
}

void write_restart_data(FILE* gridsave_file) {
  printlog("non-thermal solver, ");

  fprintf(gridsave_file, "%d\n", 24724518);  // special number marking the beginning of NT data
  fprintf(gridsave_file, "%d %la %la\n", SFPTS, SF_EMIN, SF_EMAX);

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    fprintf(gridsave_file, "%d %la ", nonemptymgi, ntlepton_deposition_rate_density_all_cells[nonemptymgi]);

    if (NT_ON && NT_SOLVE_SPENCERFANO) {
      check_auger_probabilities(nonemptymgi);

      fprintf(gridsave_file, "%a %a %a %a\n", nt_solution[nonemptymgi].nneperion_when_solved,
              nt_solution[nonemptymgi].frac_heating, nt_solution[nonemptymgi].frac_ionisation,
              nt_solution[nonemptymgi].frac_excitation);

      for (int uniqueionindex = 0; uniqueionindex < get_includedions(); uniqueionindex++) {
        const auto& celliondata = get_cell_allions_data(nonemptymgi)[uniqueionindex];
        fprintf(gridsave_file, "%la ", celliondata.fracdep_ionisation_ion);
        fprintf(gridsave_file, "%a ", celliondata.eff_ionpot);

        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          fprintf(gridsave_file, "%a %a ", celliondata.prob_num_auger[a], celliondata.ionenfrac_num_auger[a]);
        }
      }

      // write NT excitations
      fprintf(gridsave_file, "%d\n", nt_solution[nonemptymgi].frac_excitations_list_size);

      for (const auto& excitation : get_cell_ntexcitations(nonemptymgi)) {
        fprintf(gridsave_file, "%la %la %d\n", excitation.frac_deposition, excitation.ratecoeffperdeposition,
                excitation.alltransindex);
      }
    }
  }
}

void read_restart_data(FILE* gridsave_file) {
  printlnlog("Reading restart data for non-thermal solver");

  int code_check = 0;
  assert_always(fscanf(gridsave_file, "%d\n", &code_check) == 1);
  assert_always(code_check == 24724518);  // special number marking the beginning of NT data

  int sfpts_in = 0;
  double SF_EMIN_in{NAN};
  double SF_EMAX_in{NAN};
  assert_always(fscanf(gridsave_file, "%d %la %la\n", &sfpts_in, &SF_EMIN_in, &SF_EMAX_in) == 3);

  if (sfpts_in != SFPTS || SF_EMIN_in != SF_EMIN || SF_EMAX_in != SF_EMAX) {
    printlnlog("[error] gridsave file specifies {} Spencer-Fano samples, SF_EMIN {:g} SF_EMAX {:g}", sfpts_in,
               SF_EMIN_in, SF_EMAX_in);
    printlnlog("[error] This simulation has {} Spencer-Fano samples, SF_EMIN {:g} SF_EMAX {:g}", SFPTS, SF_EMIN,
               SF_EMAX);
    std::abort();
  }

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    int nonemptymgi_in = 0;
    assert_always(fscanf(gridsave_file, "%d %la ", &nonemptymgi_in,
                         &ntlepton_deposition_rate_density_all_cells[nonemptymgi]) == 2);
    assert_always(nonemptymgi_in == nonemptymgi);

    if (NT_ON && NT_SOLVE_SPENCERFANO) {
      assert_always(fscanf(gridsave_file, "%a %a %a %a\n", &nt_solution[nonemptymgi].nneperion_when_solved,
                           &nt_solution[nonemptymgi].frac_heating, &nt_solution[nonemptymgi].frac_ionisation,
                           &nt_solution[nonemptymgi].frac_excitation) == 4);

      for (int uniqueionindex = 0; uniqueionindex < get_includedions(); uniqueionindex++) {
        auto& celliondata = get_cell_allions_data(nonemptymgi)[uniqueionindex];
        assert_always(fscanf(gridsave_file, "%la ", &celliondata.fracdep_ionisation_ion) == 1);
        assert_always(fscanf(gridsave_file, "%a ", &celliondata.eff_ionpot) == 1);

        for (int a = 0; a <= NT_MAX_AUGER_ELECTRONS; a++) {
          assert_always(fscanf(gridsave_file, "%a %a ", &celliondata.prob_num_auger[a],
                               &celliondata.ionenfrac_num_auger[a]) == 2);
        }
      }

      check_auger_probabilities(nonemptymgi);

      // read NT excitations
      int frac_excitations_list_size_in = 0;
      assert_always(fscanf(gridsave_file, "%d\n", &frac_excitations_list_size_in) == 1);

      // gridsave file must not have been written with a larger per-cell excitation list capacity
      assert_always(frac_excitations_list_size_in >= 0);
      assert_always(frac_excitations_list_size_in <= nt_excitations_stored);

      nt_solution[nonemptymgi].frac_excitations_list_size = frac_excitations_list_size_in;

      for (auto& excitation : get_cell_ntexcitations(nonemptymgi)) {
        assert_always(fscanf(gridsave_file, "%la %la %d\n", &excitation.frac_deposition,
                             &excitation.ratecoeffperdeposition, &excitation.alltransindex) == 3);
      }
    }
  }
}

void nt_MPI_Bcast(const ptrdiff_t nonemptymgi, const int root_node_id) {
  if (globals::rank_in_node == 0) {
    // node-shared memory, so only node leaders participate in the internode broadcast
    // (root_node_id is only a valid root rank within the rank_in_node == 0 communicator)
    MPI_Bcast_safe(ntlepton_deposition_rate_density_all_cells[nonemptymgi], root_node_id, globals::mpi_comm_internode);
  }

  if (NT_ON && NT_SOLVE_SPENCERFANO) {
    if (globals::rank_in_node == 0) {
      MPI_Bcast_safe(nt_solution[nonemptymgi].nneperion_when_solved, root_node_id, globals::mpi_comm_internode);
      MPI_Bcast_safe(nt_solution[nonemptymgi].timestep_last_solved, root_node_id, globals::mpi_comm_internode);
      MPI_Bcast_safe(nt_solution[nonemptymgi].frac_heating, root_node_id, globals::mpi_comm_internode);
      MPI_Bcast_safe(nt_solution[nonemptymgi].frac_ionisation, root_node_id, globals::mpi_comm_internode);
      MPI_Bcast_safe(nt_solution[nonemptymgi].frac_excitation, root_node_id, globals::mpi_comm_internode);

      MPI_Bcast_safe(nt_solution[nonemptymgi].frac_excitations_list_size, root_node_id, globals::mpi_comm_internode);

      MPI_Bcast_safe(get_cell_ntexcitations(nonemptymgi), root_node_id, globals::mpi_comm_internode);

      MPI_Bcast_safe(get_cell_allions_data(nonemptymgi), root_node_id, globals::mpi_comm_internode);
    }

    MPI_Barrier_allranks();

    check_auger_probabilities(nonemptymgi);
  }
}

void reset_stats() { nt_energy_deposited = 0.; }

void print_stats(const double modelvolume, const double deltat) {
  const double deposition_rate_density_montecarlo = nt_energy_deposited / EV / modelvolume / deltat;

  printlnlog("nt_energy_deposited = {:g} [eV/s/cm^3]", deposition_rate_density_montecarlo);
}

}  // namespace nonthermal
