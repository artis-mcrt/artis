// Charge transfer reactions between the ions of the included elements. A reaction moves one electron
// from a donor ion to an acceptor ion, e.g. A+2 + B0 -> A+1 + B+1. The reactions change the ionisation
// balance of both elements and keep the free electron density constant.
//
// The rates come from two sources:
// - Published analytic fits for reactions with hydrogen and helium, in the fit form of
//   Kingdon & Ferland (1996, ApJS, 106, 205), hereafter KF96. The file data/chargetransfer.txt holds
//   these fits (see data/chargetransfer-reference.txt for the sources).
// - Landau-Zener estimates for electron capture from a neutral donor by an ion with a charge of two
//   or more. init() generates these reactions at startup from the ionisation energies and the level
//   energies of the loaded atomic dataset, so any composition (e.g. r-process ejecta) gets rates.
//   The method is the Landau-Zener approach of Butler & Dalgarno (1980, ApJ, 241, 838), which KF96
//   also used for their reactions without quantal data. The channels get an independent-crossing
//   sum with an approach probability of one, so the estimates carry an accuracy of a factor of a
//   few at best.
//
// The reverse (ionisation) direction of each generated reaction comes from detailed balance with the
// partition functions of the four ions. The solver applies a forward rate to the whole reactant
// ion populations, and the partition function ratio makes the forward and the reverse fluxes
// balance in LTE. A Boltzmann level distribution inside the lower ion turns the sum over the
// capture channels into exactly this per-ion factor, because the smaller energy defect of an
// excited channel cancels against the Boltzmann factor of that level. The NLTE population solver
// applies the rates as per-ion
// coefficients between neighbouring ion stages (see nltepop.cc). A reaction is active only when
// both elements have NLTE levels and a free ionisation balance. With only one side in the NLTE
// matrix, the partner element would keep its populations, and the reaction would break the charge
// conservation. The heat released by a reaction (the energy defect, typically ~1 eV) is not added
// to the thermal balance.

#include "chargetransfer.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <ios>
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

// atomic units
constexpr double A_BOHR_CM = 5.29177211e-9;  // Bohr radius [cm]
constexpr double E_HARTREE = 4.35974472e-11;  // Hartree energy [erg]
constexpr double V_ATOMIC_CMS = 2.18769126e8;  // atomic velocity unit [cm/s]

// velocity grid for the tabulated Landau-Zener cross sections
constexpr int LZ_VGRIDSIZE = 72;
constexpr double LZ_VGRID_MIN_CMS = 1e4;
constexpr double LZ_VGRID_MAX_CMS = 1e9;

// nodes for the impact parameter integral and for the thermal average integral
constexpr int LZ_BNODES = 48;
constexpr int LZ_THERMALNODES = 48;

// keep only the crossing channels inside this radius window [Bohr radii]. Outside the window the
// coupling makes the channel contribution negligible (Sterling & Stancil 2011 name the range
// 5 to 15-20 Bohr radii as the favourable window).
constexpr double LZ_RX_MIN_BOHR = 0.5;
constexpr double LZ_RX_MAX_BOHR = 40.;

// rate floor [cm3/s] from radiative charge transfer, which operates for every exoergic reaction
// (Butler, Guberman & Dalgarno 1977; adopted as a floor by Sterling & Stancil 2011)
constexpr double LZ_RADIATIVE_CT_FLOOR = 1e-14;

struct CTReaction {
  int partner_element{-1};  // the ion that reacts with the owner ion of the list entry
  int partner_ion{-1};
  int partner_product_ion{-1};  // the partner ion after the reaction
  double mu_grams{0.};  // reduced mass of the collision pair

  // KF96 fit parameters (a in cm3/s). Used when lz_tableindex is negative.
  double a{0.};
  double b{0.};
  double c{0.};
  double d{0.};
  double eexp{0.};
  double tmin{0.};
  double tmax{0.};

  int lz_tableindex{-1};  // >= 0 selects a tabulated Landau-Zener cross section

  // detailed balance data of a reverse (endothermic) entry: the energy release of the forward
  // reaction and the identity of its reactant ions (the acceptor captures from the donor)
  bool is_reverse{false};
  double delta_e{0.};  // [erg]
  int element_acc{-1};
  int ion_acc{-1};
  int element_don{-1};
  int ion_don{-1};
};

// reaction lists per unique ion index of the owner ion:
// - recombination list: the owner ion is the electron acceptor, transition owner ion -> (ion - 1)
// - ionisation list: the owner ion is the electron donor, transition owner ion -> (ion + 1)
std::vector<std::vector<CTReaction>> reactions_rec_perion;
std::vector<std::vector<CTReaction>> reactions_ion_perion;

// tabulated cross sections sigma(v) [cm2] on the log-spaced velocity grid
std::vector<std::array<double, LZ_VGRIDSIZE>> lz_sigma_tables;

int reaction_count = 0;

auto lz_vgrid_v_cms(const int i) -> double {
  const double logstep = std::log(LZ_VGRID_MAX_CMS / LZ_VGRID_MIN_CMS) / (LZ_VGRIDSIZE - 1);
  return LZ_VGRID_MIN_CMS * std::exp(i * logstep);
}

// a reaction is active only when both elements solve their ionisation balance in the NLTE matrix.
// With only one side in the matrix, the partner element would keep its populations, and the
// reaction would change the total ionic charge, which charge transfer must conserve.
[[nodiscard]] auto element_has_ct_balance(const int element) -> bool {
  return elem_has_nlte_levels(element) && !FORCE_SAHA_ION_BALANCE(get_atomicnumber(element));
}

// approximate mass of one atom of the element [g]
auto element_mass_g(const int element) -> double {
  const double m = globals::elements[element].initstablemeannucmass;
  if (m > 0.) {
    return m;
  }
  // the model holds no stable isotope of the element, so estimate the mass from the atomic number
  return 2.3 * get_atomicnumber(element) * MH;
}

auto reduced_mass_g(const int element_a, const int element_b) -> double {
  const double ma = element_mass_g(element_a);
  const double mb = element_mass_g(element_b);
  return ma * mb / (ma + mb);
}

// thermal average k(T) = <sigma * v> [cm3/s] over a Maxwell distribution of the relative velocity,
// with sigma(v) interpolated from a tabulated cross section
auto lz_thermal_ratecoeff(const int lz_tableindex, const double T, const double mu_grams) -> double {
  const auto& sigmatable = lz_sigma_tables[lz_tableindex];
  const double kT = KB * T;

  // integrate sigma(v(x)) * x^2 * exp(-x) over ln x, with x = E / kT
  constexpr double XMIN = 1e-3;
  constexpr double XMAX = 40.;
  const double dlnx = std::log(XMAX / XMIN) / (LZ_THERMALNODES - 1);
  const double logstep_v = std::log(LZ_VGRID_MAX_CMS / LZ_VGRID_MIN_CMS) / (LZ_VGRIDSIZE - 1);

  double integral = 0.;
  for (int i = 0; i < LZ_THERMALNODES; i++) {
    const double x = XMIN * std::exp(i * dlnx);
    const double v = std::sqrt(2. * kT * x / mu_grams);

    // interpolate log sigma on the log velocity grid, with clamping at the grid ends
    double sigma = 0.;
    if (v <= LZ_VGRID_MIN_CMS) {
      sigma = sigmatable[0];
    } else if (v >= LZ_VGRID_MAX_CMS) {
      sigma = sigmatable[LZ_VGRIDSIZE - 1];
    } else {
      const double findex = std::log(v / LZ_VGRID_MIN_CMS) / logstep_v;
      const int ilow = static_cast<int>(findex);
      const double frac = findex - ilow;
      sigma = ((1. - frac) * sigmatable[ilow]) + (frac * sigmatable[ilow + 1]);
    }

    const double weight = (i == 0 || i == LZ_THERMALNODES - 1) ? 0.5 : 1.;
    integral += weight * sigma * x * x * std::exp(-x) * dlnx;
  }

  const double k_lz = std::sqrt(8. * kT / (PI * mu_grams)) * integral;

  // every exoergic reaction keeps at least the radiative charge transfer rate
  return std::max(k_lz, LZ_RADIATIVE_CT_FLOOR);
}

auto eval_reaction_ratecoeff(const CTReaction& r, const int nonemptymgi, const double T_e) -> double {
  double k = (r.lz_tableindex >= 0) ? lz_thermal_ratecoeff(r.lz_tableindex, T_e, r.mu_grams)
                                    : chargetransfer::evaluate_ctfit(r.a, r.b, r.c, r.d, r.eexp, r.tmin, r.tmax, T_e);
  if (r.is_reverse) {
    // detailed balance for the per-ion rates: the partition functions of the forward reactants
    // over those of the forward products, with the Boltzmann factor of the energy release
    const double partfunc_ratio = (get_ion_partfunct(nonemptymgi, r.element_acc, r.ion_acc) *
                                   get_ion_partfunct(nonemptymgi, r.element_don, r.ion_don)) /
                                  (get_ion_partfunct(nonemptymgi, r.element_acc, r.ion_acc - 1) *
                                   get_ion_partfunct(nonemptymgi, r.element_don, r.ion_don + 1));
    k *= partfunc_ratio * std::exp(-r.delta_e / (KB * T_e));
  }
  return k;
}

auto sum_reaction_list(const std::vector<CTReaction>& list, const int nonemptymgi) -> double {
  if (list.empty()) {
    return 0.;
  }
  const auto T_e = grid::Te_allcells[nonemptymgi];
  double total = 0.;
  for (const auto& r : list) {
    // after a failed matrix solve, the partner element holds LTE populations in this cell and gets
    // no reciprocal transition, so skip the reaction to keep the total ionic charge constant
    // the partner ion and its product ion must both be inside the NLTE matrix solution of the
    // partner element. An edge ion that the ion range reduction removed, or an element that fell
    // back to LTE, gets no reciprocal transition, and the reaction would then change the total
    // ionic charge.
    if (!ion_in_nlte_solution(nonemptymgi, r.partner_element, r.partner_ion) ||
        !ion_in_nlte_solution(nonemptymgi, r.partner_element, r.partner_product_ion)) {
      continue;
    }
    const double nnpartner = get_nnion(nonemptymgi, r.partner_element, r.partner_ion);
    if (nnpartner > 0.) {
      total += eval_reaction_ratecoeff(r, nonemptymgi, T_e) * nnpartner;
    }
  }
  // with microclumping, the reactions happen at the in-clump partner density
  return grid::get_clumpfactor(nonemptymgi) * total;
}

// register the two forward list entries of a reaction, and the two reverse entries when autoreverse
// is set. The acceptor ion captures the electron (ion_acc -> ion_acc - 1) and the donor ion loses it
// (ion_don -> ion_don + 1).
void add_reaction(const int element_acc, const int ion_acc, const int element_don, const int ion_don,
                  const CTReaction& params, const bool autoreverse) {
  CTReaction fwd = params;
  fwd.mu_grams = reduced_mass_g(element_acc, element_don);
  fwd.is_reverse = false;

  fwd.partner_element = element_don;
  fwd.partner_ion = ion_don;
  fwd.partner_product_ion = ion_don + 1;
  reactions_rec_perion[get_uniqueionindex(element_acc, ion_acc)].push_back(fwd);

  fwd.partner_element = element_acc;
  fwd.partner_ion = ion_acc;
  fwd.partner_product_ion = ion_acc - 1;
  reactions_ion_perion[get_uniqueionindex(element_don, ion_don)].push_back(fwd);

  reaction_count++;

  if (!autoreverse) {
    return;
  }

  // energy release of the forward reaction: the captured electron binds more strongly on the
  // acceptor than it did on the donor. An endothermic file entry gets no detailed balance
  // reverse, because its listed rate is a floor value and the balance would amplify it.
  const double delta_e = get_ionpot(element_acc, ion_acc - 1) - get_ionpot(element_don, ion_don);
  if (delta_e <= 0.) {
    return;
  }

  CTReaction rev = params;
  rev.mu_grams = fwd.mu_grams;
  rev.is_reverse = true;
  rev.delta_e = delta_e;
  rev.element_acc = element_acc;
  rev.ion_acc = ion_acc;
  rev.element_don = element_don;
  rev.ion_don = ion_don;

  rev.partner_element = element_acc;
  rev.partner_ion = ion_acc - 1;
  rev.partner_product_ion = ion_acc;
  reactions_rec_perion[get_uniqueionindex(element_don, ion_don + 1)].push_back(rev);

  rev.partner_element = element_don;
  rev.partner_ion = ion_don + 1;
  rev.partner_product_ion = ion_don;
  reactions_ion_perion[get_uniqueionindex(element_acc, ion_acc - 1)].push_back(rev);

  reaction_count++;
}

// read the published fits from data/chargetransfer.txt. Return the list of (Z, ionstage) pairs of
// (acceptor, donor) that the file covered, so that the Landau-Zener generator can skip them.
auto read_reaction_file() -> std::vector<std::array<int, 4>> {
  auto covered = std::vector<std::array<int, 4>>();

  printlnlog("Reading charge transfer reactions from chargetransfer.txt...");
  auto ctfile = fstream_required("chargetransfer.txt", std::ios::in);
  std::string line;
  get_noncommentline(ctfile, line);
  int filereactioncount = 0;
  assert_always(std::istringstream{line} >> filereactioncount);
  assert_always(filereactioncount >= 0);

  int used = 0;
  // store one reaction when the loaded dataset holds both ion stages of both species
  auto store = [&](const int z_acc, const int ionstage_acc, const int z_don, const int ionstage_don,
                   const CTReaction& params, const int autoreverse) {
    covered.push_back({z_acc, ionstage_acc, z_don, ionstage_don});

    const int element_acc = get_elementindex(z_acc);
    const int element_don = get_elementindex(z_don);
    if (element_acc < 0 || element_don < 0) {
      return;
    }
    if (!element_has_ct_balance(element_acc) || !element_has_ct_balance(element_don)) {
      return;
    }
    // the reaction needs both ion stages of both species in the loaded dataset
    const int ion_acc = ionstage_acc - get_ionstage(element_acc, 0);
    const int ion_don = ionstage_don - get_ionstage(element_don, 0);
    if (ion_acc < 1 || ion_acc >= get_nions(element_acc) || ion_don < 0 || ion_don >= (get_nions(element_don) - 1)) {
      return;
    }

    add_reaction(element_acc, ion_acc, element_don, ion_don, params, autoreverse != 0);
    used++;
  };

  // the first table: each line holds a reaction with its own fit
  for (int n = 0; n < filereactioncount; n++) {
    get_noncommentline(ctfile, line);
    int z_acc{-1};
    int ionstage_acc{-1};
    int z_don{-1};
    int ionstage_don{-1};
    int autoreverse{0};
    CTReaction params{};
    double a_e9 = 0.;
    assert_always(std::istringstream{line} >> z_acc >> ionstage_acc >> z_don >> ionstage_don >> a_e9 >> params.b >>
                  params.c >> params.d >> params.eexp >> params.tmin >> params.tmax >> autoreverse);
    params.a = a_e9 * 1e-9;
    store(z_acc, ionstage_acc, z_don, ionstage_don, params, autoreverse);
  }

  // the second table: the estimates, which share one line of coefficients
  get_noncommentline(ctfile, line);
  int fileestimatecount = 0;
  assert_always(std::istringstream{line} >> fileestimatecount);
  assert_always(fileestimatecount >= 0);
  CTReaction estimateparams{};
  {
    get_noncommentline(ctfile, line);
    double a_e9 = 0.;
    assert_always(std::istringstream{line} >> a_e9 >> estimateparams.b >> estimateparams.c >> estimateparams.d >>
                  estimateparams.eexp >> estimateparams.tmin >> estimateparams.tmax);
    estimateparams.a = a_e9 * 1e-9;
  }
  for (int n = 0; n < fileestimatecount; n++) {
    get_noncommentline(ctfile, line);
    int z_acc{-1};
    int ionstage_acc{-1};
    int z_don{-1};
    int ionstage_don{-1};
    int autoreverse{0};
    assert_always(std::istringstream{line} >> z_acc >> ionstage_acc >> z_don >> ionstage_don >> autoreverse);
    store(z_acc, ionstage_acc, z_don, ionstage_don, estimateparams, autoreverse);
  }
  printlnlog("Stored {} of {} charge transfer reactions from the file", used, filereactioncount + fileestimatecount);

  return covered;
}

// generate the Landau-Zener reactions: electron capture from every neutral donor by every ion with a
// charge of two or more, for the pairs that the data file does not cover
void generate_lz_reactions(const std::vector<std::array<int, 4>>& covered) {
  int lz_reaction_count = 0;
  int lz_channel_count = 0;

  for (int element_acc = 0; element_acc < get_nelements(); element_acc++) {
    if (!element_has_ct_balance(element_acc)) {
      continue;
    }
    for (int ion_acc = 1; ion_acc < get_nions(element_acc); ion_acc++) {
      const int ioncharge = get_ionstage(element_acc, ion_acc) - 1;
      if (ioncharge < 2) {
        continue;
      }
      for (int element_don = 0; element_don < get_nelements(); element_don++) {
        if (get_ionstage(element_don, 0) != 1 || get_nions(element_don) < 2 || !element_has_ct_balance(element_don)) {
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

        const double ip_donor = get_ionpot(element_don, ion_don);
        const double ip_final_ground = get_ionpot(element_acc, ion_acc - 1);

        // the reaction must release energy for a curve crossing and for the radiative floor
        if ((ip_final_ground - ip_donor) <= 0.) {
          continue;
        }

        // one crossing channel per final level of the lower acceptor ion. The channel cross
        // sections add, which follows the independent-crossing treatment of BD80 and KF96.
        auto sigmatable = std::array<double, LZ_VGRIDSIZE>{};
        int nchannels = 0;
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
          nchannels++;
          for (int i = 0; i < LZ_VGRIDSIZE; i++) {
            sigmatable[i] += chargetransfer::sigma_lz_channel(ioncharge, deltae, ip_donor, lz_vgrid_v_cms(i));
          }
        }
        lz_channel_count += nchannels;

        CTReaction params{};
        params.lz_tableindex = static_cast<int>(lz_sigma_tables.size());
        lz_sigma_tables.push_back(sigmatable);

        add_reaction(element_acc, ion_acc, element_don, ion_don, params, true);
        lz_reaction_count++;
      }
    }
  }
  printlnlog("Generated {} Landau-Zener charge transfer reactions with {} capture channels", lz_reaction_count,
             lz_channel_count);
}

}  // anonymous namespace

namespace chargetransfer {

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

auto sigma_lz_channel(const int ioncharge, const double deltae_erg, const double ip_donor_erg, const double v_cms)
    -> double {
  // diabatic curves in atomic units: the entrance channel (ion + neutral) is taken as flat, and the
  // exit channel has the Coulomb repulsion of the two product ions plus the energy release. The
  // crossing sits where the repulsion equals the energy release, so rx = zeta / deltae with
  // zeta the product of the product ion charges (Butler & Dalgarno 1980, hereafter BD80). This is
  // the zeroth-order crossing radius of BD80, which they found accurate for most cases of interest.
  const double deltae_au = deltae_erg / E_HARTREE;
  const int zeta = ioncharge - 1;
  const double rx_au = zeta / deltae_au;

  // slope difference of the two diabatic curves at the crossing (BD80 with their mu = 0)
  const double dfdr_au = zeta / (rx_au * rx_au);

  // adiabatic splitting at the crossing: the BD80 fit dU = rx^2 exp(-alpha rx) hartree, with
  // alpha = sqrt(2 IP) from the ionisation potential of the neutral donor (BD80 give alpha = 1.0
  // for H and 1.34 for He, and state the sqrt(2 IP) rule for the exponent)
  const double alpha = std::sqrt(2. * ip_donor_erg / E_HARTREE);
  const double du_au = rx_au * rx_au * std::exp(-alpha * rx_au);

  const double v_au = v_cms / V_ATOMIC_CMS;
  if (v_au <= 0.) {
    return 0.;
  }

  // Landau-Zener single passage survival p = exp(-w / t) at radial velocity fraction
  // t = v_radial / v, with w = pi^2 dU^2 / (h v F) and h = 2 pi in atomic units.
  // The net transfer probability over the inbound and the outbound crossing is 2 p (1 - p).
  const double w = PI / 2. * du_au * du_au / (dfdr_au * v_au);

  // sigma = 4 pi rx^2 * integral of p (1 - p) t dt over t from 0 to 1, which is the impact
  // parameter average of BD80 equation (3) with their lambda = 0 and approach probability 1
  const double dt = 1. / LZ_BNODES;
  double integral = 0.;
  for (int i = 0; i < LZ_BNODES; i++) {
    const double t = (i + 0.5) * dt;
    const double p = std::exp(-w / t);
    integral += p * (1. - p) * t * dt;
  }

  const double rx_cm = rx_au * A_BOHR_CM;
  return 4. * PI * rx_cm * rx_cm * integral;
}

void init() {
  if (!ENABLE_CHARGE_TRANSFER_REACTIONS) {
    return;
  }
  reactions_rec_perion.resize(get_includedions());
  reactions_ion_perion.resize(get_includedions());

  const auto covered = read_reaction_file();
  generate_lz_reactions(covered);

  printlnlog("Charge transfer setup stored {} reactions (each with a reverse where detailed balance applies)",
             reaction_count);
}

auto get_reaction_count() -> int { return reaction_count; }

auto ct_recombination_rate(const int nonemptymgi, const int element, const int upperion) -> double {
  if (!ENABLE_CHARGE_TRANSFER_REACTIONS) {
    return 0.;
  }
  return sum_reaction_list(reactions_rec_perion[get_uniqueionindex(element, upperion)], nonemptymgi);
}

auto ct_ionisation_rate(const int nonemptymgi, const int element, const int ion) -> double {
  if (!ENABLE_CHARGE_TRANSFER_REACTIONS) {
    return 0.;
  }
  return sum_reaction_list(reactions_ion_perion[get_uniqueionindex(element, ion)], nonemptymgi);
}

}  // namespace chargetransfer
