// Charge transfer reactions between the ions of the included elements. A reaction moves one electron
// from a donor ion to an acceptor ion, e.g. A+2 + B0 -> A+1 + B+1. The reactions change the ionisation
// balance of both elements and keep the free electron density constant.
//
// The rates come from three sources:
// - Published analytic fits for reactions with hydrogen and helium, in the fit form of
//   Kingdon & Ferland (1996, ApJS, 106, 205), hereafter KF96. The file data/chargetransfer.txt holds
//   these fits, and its header names the sources. KF96 give no rate for the exothermic ionisation
//   of Ca+, Sc+, and Ti+ by protons, so those reactions are absent from the file as well.
// - A flat estimate for the other electron captures from a neutral donor by a singly charged ion.
//   A singly charged ion has no Coulomb curve crossing, so the energy release does not predict the
//   rate. The tabulated rates of such reactions with an energy release up to 4 eV spread from
//   1e-15 to 1e-9 cm3/s, with a median near 1e-12 cm3/s at 1e4 K. These reactions get that median,
//   and a larger energy release gets the radiative floor.
// - Landau-Zener estimates for the other electron captures from a neutral donor by an ion with a
//   charge of two or more. The method is the Landau-Zener approach of Butler & Dalgarno (1980,
//   ApJ, 241, 838), hereafter BD80, which KF96 also used for their reactions without quantal data.
//   Each level of the lower acceptor ion is one capture channel. The trajectory crosses the
//   channels in sequence from the outside inwards and back, with the classical two-state mixing of
//   Landau and Zener at each crossing, so the total transfer probability stays below one. The
//   approach probability is one. These estimates carry an accuracy of a factor of a few at best.
// init() generates the estimates at startup from the ionisation energies and the level energies of
// the loaded atomic dataset, so any composition (e.g. r-process ejecta) gets rates. The
// Landau-Zener rate coefficients go into tables on the temperature grid of the rate coefficients.
//
// The reverse (ionisation) direction of each generated reaction comes from detailed balance with the
// partition functions of the four ions. The solver applies a forward rate to the whole reactant ion
// populations, and the partition function ratio makes the forward and the reverse fluxes balance in
// LTE. A Boltzmann level distribution inside the lower ion turns the sum over the capture channels
// into exactly this per-ion factor. The smaller energy release of an excited channel cancels against
// the Boltzmann factor of that level.
//
// The NLTE population solver applies the rates as per-ion coefficients between neighbouring ion
// stages (see nltepop.cc). A reaction is active only when both elements have NLTE levels and a free
// ionisation balance. With only one side in the NLTE matrix, the partner element would keep its
// populations, and the reaction would break the charge conservation. The solver does not add the
// heat of a reaction (the energy release, typically ~1 eV) to the thermal balance.

#include "chargetransfer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <ios>
#include <span>
#include <sstream>
#include <string>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "ltepop.h"
#include "mpi_logging.h"
#include "nltepop.h"

namespace {

// velocity grid of the Landau-Zener cross sections [cm/s]
constexpr int LZ_VGRIDSIZE = 72;
constexpr double LZ_VGRID_MIN_CMS = 1e4;
constexpr double LZ_VGRID_MAX_CMS = 1e9;
const double LZ_VGRID_DLOG = std::log(LZ_VGRID_MAX_CMS / LZ_VGRID_MIN_CMS) / (LZ_VGRIDSIZE - 1);

// the temperature grid of the rate coefficient tables, the same grid as in ratecoeff.cc
const double T_STEP_LOG = (std::log(MAXTEMP) - std::log(MINTEMP)) / (TABLESIZE - 1.);

// nodes of the impact parameter integral and of the thermal average integral
constexpr int LZ_BNODES = 48;
constexpr int LZ_THERMALNODES = 48;

// keep only the crossing channels inside this radius window [Bohr radii]. Outside the window the
// coupling makes the channel contribution negligible (Sterling & Stancil 2011 name the range
// 5 to 15-20 Bohr radii as the favourable window).
constexpr double LZ_RX_MIN_BOHR = 0.5;
constexpr double LZ_RX_MAX_BOHR = 40.;

// rate floor [cm3/s] from radiative charge transfer, which operates for every exoergic reaction
// (Butler, Guberman & Dalgarno 1977; adopted as a floor by Sterling & Stancil 2011)
constexpr double RADIATIVE_CT_FLOOR = 1e-14;

// the flat estimate [cm3/s] for a capture by a singly charged ion with an energy release up to
// SINGLY_CHARGED_MAX_DEFECT [erg]: the median of the tabulated rates of such reactions at 1e4 K
constexpr double SINGLY_CHARGED_RATE = 1e-12;
constexpr double SINGLY_CHARGED_MAX_DEFECT = 4. * EV;

struct CTReaction {
  int partner_element{-1};  // the ion that reacts with the owner ion of the list entry
  int partner_ion{-1};  // the partner ion before the reaction. The list role gives its product ion.

  // KF96 fit parameters (a in cm3/s). Used when ratetable_index is negative.
  double a{0.};
  double b{0.};
  double c{0.};
  double d{0.};
  double eexp{0.};
  double tmin{0.};
  double tmax{0.};

  int ratetable_index{-1};  // >= 0 selects a tabulated Landau-Zener rate coefficient

  // a reverse (endothermic) entry carries the energy release of its exothermic forward reaction [erg]
  bool is_reverse{false};
  double delta_e{0.};
};

// reaction lists per unique ion index of the owner ion:
// - recombination list: the owner ion is the electron acceptor, transition owner ion -> (ion - 1)
// - ionisation list: the owner ion is the electron donor, transition owner ion -> (ion + 1)
// The list role gives the step of the owner ion, -1 or +1, and the partner ion takes the opposite step.
std::vector<std::vector<CTReaction>> reactions_rec_perion;
std::vector<std::vector<CTReaction>> reactions_ion_perion;
constexpr int STEP_RECOMBINATION = -1;
constexpr int STEP_IONISATION = 1;

// tabulated Landau-Zener rate coefficients [cm3/s] on the temperature grid
std::vector<std::array<double, TABLESIZE>> lz_ratecoeff_tables;

// one capture channel at its avoided crossing, in atomic units. The diabatic passage probability at
// the radial velocity v_r is exp(-w / v_r).
struct LZChannel {
  double rx_au{0.};
  double w_au{0.};
};

auto lz_vgrid_v_cms(const int i) -> double { return LZ_VGRID_MIN_CMS * std::exp(i * LZ_VGRID_DLOG); }

auto ratetable_temperature(const int i) -> double { return MINTEMP * std::exp(i * T_STEP_LOG); }

// a reaction is active only when both elements solve their ionisation balance in the NLTE matrix.
// With only one side in the matrix, the partner element would keep its populations, and the
// reaction would change the total ionic charge, which charge transfer must conserve.
[[nodiscard]] auto element_has_ct_balance(const int element) -> bool {
  return elem_has_nlte_levels(element) && !FORCE_SAHA_ION_BALANCE(get_atomicnumber(element));
}

// reduced mass of the collision pair [g], from the mean nuclear masses of compositiondata.txt
auto reduced_mass_g(const int element_a, const int element_b) -> double {
  const double ma = globals::elements[element_a].initstablemeannucmass;
  const double mb = globals::elements[element_b].initstablemeannucmass;
  assert_always(ma > 0. && mb > 0.);
  return ma * mb / (ma + mb);
}

// the crossing of one capture channel (BD80): the entrance channel (ion + neutral) is flat, and the
// exit channel has the Coulomb repulsion of the two product ions minus the energy release. The
// crossing sits where the repulsion equals the energy release, so rx = zeta / deltae with zeta the
// product of the product ion charges. This is the zeroth-order crossing radius of BD80, which they
// found accurate for most cases of interest. The adiabatic splitting at the crossing is the BD80 fit
// dU = rx^2 exp(-alpha rx) hartree, with alpha = sqrt(2 IP) from the ionisation potential of the
// neutral donor (BD80 give alpha = 1.0 for H and 1.34 for He, and state the sqrt(2 IP) rule).
auto make_lz_channel(const int ioncharge, const double deltae_erg, const double ip_donor_erg) -> LZChannel {
  // the crossing radius scales with (charge - 1) and needs an energy release
  assert_always(ioncharge >= 2);
  assert_always(deltae_erg > 0.);
  const double deltae_au = deltae_erg / E_HARTREE;
  const int zeta = ioncharge - 1;
  const double rx_au = zeta / deltae_au;

  // slope difference of the two diabatic curves at the crossing (BD80 with their mu = 0)
  const double dfdr_au = zeta / (rx_au * rx_au);

  const double alpha = std::sqrt(2. * ip_donor_erg / E_HARTREE);
  const double du_au = rx_au * rx_au * std::exp(-alpha * rx_au);

  // Landau-Zener single passage survival p = exp(-w / v_radial), with w = pi^2 dU^2 / (h v F)
  // and h = 2 pi in atomic units
  return {.rx_au = rx_au, .w_au = PI / 2. * du_au * du_au / dfdr_au};
}

// cross section [cm2] of a set of capture channels, sorted by the crossing radius from the outside
// inwards, at the relative velocity v. The trajectory at impact parameter b crosses every channel
// with a radius above b on the way in and again on the way out. At each crossing the classical
// two-state mixing moves probability between the entrance channel and that capture channel. The
// transfer probability is the probability that the trajectory leaves the entrance channel.
auto sigma_lz_sorted_channels(const std::span<const LZChannel> channels, const double v_cms) -> double {
  if (channels.empty() || v_cms <= 0.) {
    return 0.;
  }
  const double v_au = v_cms / V_ATOMIC_CMS;
  const double bmax_au = channels.front().rx_au;
  const double db_au = bmax_au / LZ_BNODES;

  std::vector<double> p_pass(channels.size());
  std::vector<double> captured(channels.size());

  double integral_au2 = 0.;
  for (int i = 0; i < LZ_BNODES; i++) {
    const double b_au = (i + 0.5) * db_au;

    // inbound: a crossing takes the fraction (1 - p) out of the entrance channel
    double p_entrance = 1.;
    std::size_t ncrossed = 0;
    for (; ncrossed < channels.size() && channels[ncrossed].rx_au > b_au; ncrossed++) {
      const auto& channel = channels[ncrossed];
      const double t = std::sqrt(1. - ((b_au * b_au) / (channel.rx_au * channel.rx_au)));
      const double p = std::exp(-channel.w_au / (v_au * t));
      p_pass[ncrossed] = p;
      captured[ncrossed] = p_entrance * (1. - p);
      p_entrance *= p;
    }

    // outbound: the same crossings in reverse order mix the entrance channel with each capture channel
    for (std::size_t k = 0; k < ncrossed; k++) {
      const std::size_t j = ncrossed - 1 - k;
      const double p = p_pass[j];
      const double p_entrance_new = (p_entrance * p) + (captured[j] * (1. - p));
      captured[j] = (captured[j] * p) + (p_entrance * (1. - p));
      p_entrance = p_entrance_new;
    }

    integral_au2 += (1. - p_entrance) * b_au * db_au;
  }

  return 2. * PI * integral_au2 * A_BOHR_CM * A_BOHR_CM;
}

// thermal average k(T) = <sigma * v> [cm3/s] over a Maxwell distribution of the relative velocity,
// with sigma(v) interpolated from the tabulated cross section on the velocity grid
auto lz_thermal_average(const std::array<double, LZ_VGRIDSIZE>& sigmatable, const double T, const double mu_grams)
    -> double {
  const double kT = KB * T;

  // integrate sigma(v(x)) * x^2 * exp(-x) over ln x, with x = E / kT
  constexpr double XMIN = 1e-3;
  constexpr double XMAX = 40.;
  const double dlnx = std::log(XMAX / XMIN) / (LZ_THERMALNODES - 1);

  double integral = 0.;
  for (int i = 0; i < LZ_THERMALNODES; i++) {
    const double x = XMIN * std::exp(i * dlnx);
    const double v = std::sqrt(2. * kT * x / mu_grams);

    double sigma = 0.;
    if (v <= LZ_VGRID_MIN_CMS) {
      sigma = sigmatable[0];
    } else if (v >= LZ_VGRID_MAX_CMS) {
      sigma = sigmatable[LZ_VGRIDSIZE - 1];
    } else {
      const double findex = std::log(v / LZ_VGRID_MIN_CMS) / LZ_VGRID_DLOG;
      // findex can round up to the last grid point, so keep ilow + 1 inside the table
      const int ilow = std::min(static_cast<int>(findex), LZ_VGRIDSIZE - 2);
      const double frac = findex - ilow;
      sigma = ((1. - frac) * sigmatable[ilow]) + (frac * sigmatable[ilow + 1]);
    }

    const double weight = (i == 0 || i == LZ_THERMALNODES - 1) ? 0.5 : 1.;
    integral += weight * sigma * x * x * std::exp(-x) * dlnx;
  }

  return std::sqrt(8. * kT / (PI * mu_grams)) * integral;
}

// tabulate the rate coefficient of a reaction on the temperature grid, with the radiative floor
auto make_lz_ratecoeff_table(const std::span<const LZChannel> channels, const double mu_grams)
    -> std::array<double, TABLESIZE> {
  std::array<double, LZ_VGRIDSIZE> sigmatable{};
  for (int i = 0; i < LZ_VGRIDSIZE; i++) {
    sigmatable[i] = sigma_lz_sorted_channels(channels, lz_vgrid_v_cms(i));
  }

  std::array<double, TABLESIZE> ratetable{};
  for (int i = 0; i < TABLESIZE; i++) {
    ratetable[i] = std::max(lz_thermal_average(sigmatable, ratetable_temperature(i), mu_grams), RADIATIVE_CT_FLOOR);
  }
  return ratetable;
}

// interpolate a tabulated rate coefficient at the temperature T, with clamping at the grid ends
auto lz_ratecoeff(const int ratetable_index, const double T) -> double {
  const auto& ratetable = lz_ratecoeff_tables[ratetable_index];
  if (T <= MINTEMP) {
    return ratetable[0];
  }
  if (T >= MAXTEMP) {
    return ratetable[TABLESIZE - 1];
  }
  const double findex = std::log(T / MINTEMP) / T_STEP_LOG;
  const int ilow = std::min(static_cast<int>(findex), TABLESIZE - 2);
  const double frac = findex - ilow;
  return ((1. - frac) * ratetable[ilow]) + (frac * ratetable[ilow + 1]);
}

// the rate coefficient of a list entry of the owner ion (element, ion) with the list step
auto eval_reaction_ratecoeff(const CTReaction& r, const int nonemptymgi, const double T_e, const int element,
                             const int ion, const int step) -> double {
  double k = (r.ratetable_index >= 0) ? lz_ratecoeff(r.ratetable_index, T_e)
                                      : chargetransfer::evaluate_ctfit(r.a, r.b, r.c, r.d, r.eexp, r.tmin, r.tmax, T_e);
  if (r.is_reverse) {
    // detailed balance for the per-ion rates: the products of this entry are the reactants of the
    // exothermic forward reaction, so the partition functions of the products go over those of the
    // reactants, with the Boltzmann factor of the energy release
    const double partfunc_ratio = (get_ion_partfunct(nonemptymgi, element, ion + step) *
                                   get_ion_partfunct(nonemptymgi, r.partner_element, r.partner_ion - step)) /
                                  (get_ion_partfunct(nonemptymgi, element, ion) *
                                   get_ion_partfunct(nonemptymgi, r.partner_element, r.partner_ion));
    k *= partfunc_ratio * std::exp(-r.delta_e / (KB * T_e));
  }
  return k;
}

// sum the rates of the list entries of the owner ion (element, ion) with the list step, while the
// NLTE matrix of the element covers the ion range [first_ion_used, first_ion_used + nions_used)
auto sum_reaction_list(const std::vector<CTReaction>& list, const int nonemptymgi, const int element, const int ion,
                       const int step, const int first_ion_used, const int nions_used) -> double {
  if (list.empty()) {
    return 0.;
  }
  const auto T_e = grid::Te_allcells[nonemptymgi];
  double total = 0.;
  for (const auto& r : list) {
    // the partner ion and its product ion must both be inside the NLTE matrix solution of the
    // partner element. An edge ion that the ion range reduction removed, or an element that fell
    // back to LTE, gets no reciprocal transition, and the reaction would then change the total
    // ionic charge. A self-reaction has the owner element as its partner. The matrix under
    // assembly then gives the range, because the stored range belongs to the previous solve.
    const bool partner_is_owner = (r.partner_element == element);
    const auto ion_in_partner_solution = [&](const int partner_ion) {
      return partner_is_owner ? (partner_ion >= first_ion_used && partner_ion < (first_ion_used + nions_used))
                              : ion_in_nlte_solution(nonemptymgi, r.partner_element, partner_ion);
    };
    if (!ion_in_partner_solution(r.partner_ion) || !ion_in_partner_solution(r.partner_ion - step)) {
      continue;
    }
    const double nnpartner = get_nnion(nonemptymgi, r.partner_element, r.partner_ion);
    if (nnpartner > 0.) {
      total += eval_reaction_ratecoeff(r, nonemptymgi, T_e, element, ion, step) * nnpartner;
    }
  }
  // with microclumping, the reactions happen at the in-clump partner density
  return grid::get_clumpfactor(nonemptymgi) * total;
}

// register the two forward list entries of a reaction, and the two reverse entries when autoreverse
// is set. The acceptor ion captures the electron (ion_acc -> ion_acc - 1) and the donor ion loses it
// (ion_don -> ion_don + 1). Return true when the reverse entries exist.
auto add_reaction(const int element_acc, const int ion_acc, const int element_don, const int ion_don,
                  const CTReaction& params, const bool autoreverse) -> bool {
  CTReaction fwd = params;
  fwd.is_reverse = false;

  fwd.partner_element = element_don;
  fwd.partner_ion = ion_don;
  reactions_rec_perion[get_uniqueionindex(element_acc, ion_acc)].push_back(fwd);

  fwd.partner_element = element_acc;
  fwd.partner_ion = ion_acc;
  reactions_ion_perion[get_uniqueionindex(element_don, ion_don)].push_back(fwd);

  if (!autoreverse) {
    return false;
  }

  // energy release of the forward reaction: the captured electron binds more strongly on the
  // acceptor than it did on the donor. An endothermic file entry gets no detailed balance
  // reverse, because its listed rate is a floor value and the balance would amplify it.
  const double delta_e = get_ionpot(element_acc, ion_acc - 1) - get_ionpot(element_don, ion_don);
  if (delta_e <= 0.) {
    return false;
  }

  CTReaction rev = params;
  rev.is_reverse = true;
  rev.delta_e = delta_e;

  rev.partner_element = element_acc;
  rev.partner_ion = ion_acc - 1;
  reactions_rec_perion[get_uniqueionindex(element_don, ion_don + 1)].push_back(rev);

  rev.partner_element = element_don;
  rev.partner_ion = ion_don + 1;
  reactions_ion_perion[get_uniqueionindex(element_acc, ion_acc - 1)].push_back(rev);

  return true;
}

// read the published fits from data/chargetransfer.txt. Return the list of (Z, ionstage) pairs of
// (acceptor, donor) that the file covered, with the reverse pair of each reaction that got a
// detailed balance reverse, so that the estimate generator can skip them.
auto read_reaction_file() -> std::vector<std::array<int, 4>> {
  auto covered = std::vector<std::array<int, 4>>();

  printlnlog("Reading charge transfer reactions from chargetransfer.txt...");
  auto ctfile = fstream_required("chargetransfer.txt", std::ios::in);
  std::string line;
  assert_always(get_noncommentline(ctfile, line));
  int filereactioncount = 0;
  assert_always(std::istringstream{line} >> filereactioncount);
  assert_always(filereactioncount >= 0);

  int used = 0;
  for (int n = 0; n < filereactioncount; n++) {
    assert_always(get_noncommentline(ctfile, line));
    int z_acc{-1};
    int ionstage_acc{-1};
    int z_don{-1};
    int ionstage_don{-1};
    int autoreverse{0};
    CTReaction params{};
    double a_e9 = 0.;
    assert_always(std::istringstream{line} >> z_acc >> ionstage_acc >> z_don >> ionstage_don >> a_e9 >> params.b >>
                  params.c >> params.d >> params.eexp >> params.tmin >> params.tmax >> autoreverse);
    assert_always(z_acc > 0 && ionstage_acc > 1 && z_don > 0 && ionstage_don > 0);
    assert_always(a_e9 >= 0.);
    assert_always(params.tmin > 0. && params.tmin <= params.tmax);
    params.a = a_e9 * 1e-9;

    covered.push_back({z_acc, ionstage_acc, z_don, ionstage_don});

    const int element_acc = get_elementindex(z_acc);
    const int element_don = get_elementindex(z_don);
    if (element_acc < 0 || element_don < 0) {
      continue;
    }
    if (!element_has_ct_balance(element_acc) || !element_has_ct_balance(element_don)) {
      continue;
    }
    // the reaction needs both ion stages of both species in the loaded dataset
    const int ion_acc = ionstage_acc - get_ionstage(element_acc, 0);
    const int ion_don = ionstage_don - get_ionstage(element_don, 0);
    if (ion_acc < 1 || ion_acc >= get_nions(element_acc) || ion_don < 0 || ion_don >= (get_nions(element_don) - 1)) {
      continue;
    }

    if (add_reaction(element_acc, ion_acc, element_don, ion_don, params, autoreverse != 0)) {
      covered.push_back({z_don, ionstage_don + 1, z_acc, ionstage_acc - 1});
    }
    used++;
  }
  printlnlog("Stored {} of {} charge transfer reactions from the file", used, filereactioncount);
  if (used == 0 && filereactioncount > 0) {
    printlnlog(
        "[warning] No charge transfer reaction of the file applies to this model. The file holds reactions with "
        "hydrogen and helium, which need NLTE levels and a free ionisation balance.");
  }

  return covered;
}

// generate the estimated reactions: electron capture from every neutral donor by every ion, for the
// pairs that the data file does not cover. A singly charged ion gets the flat estimate, and an ion
// with a charge of two or more gets a Landau-Zener rate.
void generate_estimated_reactions(const std::vector<std::array<int, 4>>& covered) {
  int singly_charged_count = 0;
  int lz_reaction_count = 0;
  int lz_channel_count = 0;

  for (int element_acc = 0; element_acc < get_nelements(); element_acc++) {
    if (!element_has_ct_balance(element_acc)) {
      continue;
    }
    for (int ion_acc = 1; ion_acc < get_nions(element_acc); ion_acc++) {
      const int ioncharge = get_ionstage(element_acc, ion_acc) - 1;
      if (ioncharge < 1) {
        continue;
      }
      for (int element_don = 0; element_don < get_nelements(); element_don++) {
        if (get_nions(element_don) < 2 || get_ionstage(element_don, 0) != 1 || !element_has_ct_balance(element_don)) {
          continue;
        }
        const int ion_don = 0;
        const auto key = std::array<int, 4>{
            get_atomicnumber(element_acc),
            get_ionstage(element_acc, ion_acc),
            get_atomicnumber(element_don),
            1,
        };
        if (std::ranges::find(covered, key) != covered.end()) {
          continue;
        }
        // the file can hold the reverse reaction as its own fit, and the estimate then gets no
        // detailed balance reverse
        const auto reverse_key = std::array<int, 4>{
            get_atomicnumber(element_don),
            2,
            get_atomicnumber(element_acc),
            get_ionstage(element_acc, ion_acc) - 1,
        };
        const bool autoreverse = (std::ranges::find(covered, reverse_key) == covered.end());

        const double ip_donor = get_ionpot(element_don, ion_don);
        const double ip_final_ground = get_ionpot(element_acc, ion_acc - 1);

        // the reaction must release energy for a curve crossing and for the radiative floor
        const double deltae_ground = ip_final_ground - ip_donor;
        if (deltae_ground <= 0.) {
          continue;
        }

        if (ioncharge == 1) {
          CTReaction params{};
          params.a = chargetransfer::singly_charged_ratecoeff(deltae_ground);
          add_reaction(element_acc, ion_acc, element_don, ion_don, params, autoreverse);
          singly_charged_count++;
          continue;
        }

        // one capture channel per final level of the lower acceptor ion, sorted by the crossing
        // radius from the outside inwards
        std::vector<LZChannel> channels;
        const int nlevels_lower = get_nlevels(element_acc, ion_acc - 1);
        for (int level = 0; level < nlevels_lower; level++) {
          const double excitation = epsilon(element_acc, ion_acc - 1, level) - epsilon(element_acc, ion_acc - 1, 0);
          const double ip_final = ip_final_ground - excitation;
          const double deltae = ip_final - ip_donor;
          if (ip_final <= 0. || deltae <= 0.) {
            continue;
          }
          const double rx_bohr = (ioncharge - 1) * E_HARTREE / deltae;
          if (rx_bohr < LZ_RX_MIN_BOHR || rx_bohr > LZ_RX_MAX_BOHR) {
            continue;
          }
          channels.push_back(make_lz_channel(ioncharge, deltae, ip_donor));
        }
        // equal radii give equal channels, so the order of the ties does not change the result
        std::ranges::stable_sort(channels, std::ranges::greater{}, &LZChannel::rx_au);
        lz_channel_count += static_cast<int>(channels.size());

        CTReaction params{};
        params.ratetable_index = static_cast<int>(lz_ratecoeff_tables.size());
        lz_ratecoeff_tables.push_back(make_lz_ratecoeff_table(channels, reduced_mass_g(element_acc, element_don)));

        add_reaction(element_acc, ion_acc, element_don, ion_don, params, autoreverse);
        lz_reaction_count++;
      }
    }
  }
  printlnlog("Generated {} charge transfer estimates for singly charged ions (flat rate {} cm3/s up to {} eV)",
             singly_charged_count, SINGLY_CHARGED_RATE, SINGLY_CHARGED_MAX_DEFECT / EV);
  printlnlog("Generated {} Landau-Zener charge transfer reactions with {} capture channels", lz_reaction_count,
             lz_channel_count);
}

}  // anonymous namespace

namespace chargetransfer {

auto singly_charged_ratecoeff(const double deltae_erg) -> double {
  if (deltae_erg <= 0.) {
    return 0.;
  }
  return (deltae_erg <= SINGLY_CHARGED_MAX_DEFECT) ? SINGLY_CHARGED_RATE : RADIATIVE_CT_FLOOR;
}

auto evaluate_ctfit(const double a, const double b, const double c, const double d, const double eexp,
                    const double tmin, const double tmax, const double T) -> double {
  const double t_clamped = std::clamp(T, tmin, tmax);
  const double t4 = t_clamped / 1e4;
  double k = a * std::pow(t4, b) * (1. + (c * std::exp(d * t4)));
  if (eexp > 0.) {
    // the Boltzmann suppression of an endothermic fit uses the true temperature
    k *= std::exp(-eexp / T);
  }
  return std::max(k, 0.);
}

auto sigma_lz_channels(const int ioncharge, const std::span<const double> deltae_erg_list, const double ip_donor_erg,
                       const double v_cms) -> double {
  std::vector<LZChannel> channels;
  channels.reserve(deltae_erg_list.size());
  for (const double deltae_erg : deltae_erg_list) {
    channels.push_back(make_lz_channel(ioncharge, deltae_erg, ip_donor_erg));
  }
  std::ranges::stable_sort(channels, std::ranges::greater{}, &LZChannel::rx_au);
  return sigma_lz_sorted_channels(channels, v_cms);
}

auto sigma_lz_channel(const int ioncharge, const double deltae_erg, const double ip_donor_erg, const double v_cms)
    -> double {
  const auto deltae_list = std::array<double, 1>{deltae_erg};
  return sigma_lz_channels(ioncharge, deltae_list, ip_donor_erg, v_cms);
}

void init() {
  if (!ENABLE_CHARGE_TRANSFER_REACTIONS) {
    return;
  }
  reactions_rec_perion.resize(get_includedions());
  reactions_ion_perion.resize(get_includedions());

  const auto covered = read_reaction_file();
  generate_estimated_reactions(covered);
}

auto ct_recombination_rate(const int nonemptymgi, const int element, const int upperion, const int first_ion_used,
                           const int nions_used) -> double {
  if (!ENABLE_CHARGE_TRANSFER_REACTIONS) {
    return 0.;
  }
  return sum_reaction_list(reactions_rec_perion[get_uniqueionindex(element, upperion)], nonemptymgi, element, upperion,
                           STEP_RECOMBINATION, first_ion_used, nions_used);
}

auto ct_ionisation_rate(const int nonemptymgi, const int element, const int ion, const int first_ion_used,
                        const int nions_used) -> double {
  if (!ENABLE_CHARGE_TRANSFER_REACTIONS) {
    return 0.;
  }
  return sum_reaction_list(reactions_ion_perion[get_uniqueionindex(element, ion)], nonemptymgi, element, ion,
                           STEP_IONISATION, first_ion_used, nions_used);
}

}  // namespace chargetransfer
