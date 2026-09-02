// The per-timestep packet update: sorts the packet population and advances every packet to the
// end of the timestep, handling pellet decays and non-thermal deposition and dispatching gamma,
// r-, and k-packets to their propagation and interaction routines.

#include "update_packets.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <optional>
#ifdef GPU_ON
#include <ranges>
#endif
#include <span>
#include <tuple>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "gammapkt.h"
#include "globals.h"
#include "grid.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "mpi_logging.h"
#include "nonthermal.h"
#include "packet.h"
#include "random.h"
#include "rpkt.h"
#include "stats.h"
#include "vectors.h"

namespace {

void do_nonthermal_predeposit(Packet& pkt, const int nts, const double ts_end) {
  // handle deposition by non-thermal alpha and beta particles that are emitted from pellets and then
  // deposit some or all of their energy locally in the ejecta (possibly after some time delay).
  double e_cmf_deposited = pkt.e_cmf;
  const auto nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);
  const auto priortype = pkt.type;
  const double ts = pkt.prop_time;
  const auto deposit_type =
      (pkt.type == TYPE_NONTHERMAL_PREDEPOSIT_ALPHA) ? TYPE_NTALPHA_FISPROD_DEPOSITED : TYPE_NTLEPTON_DEPOSITED;

  // deposit the packet with probability f_p, otherwise discard its deposited energy and let it escape
  [[maybe_unused]] const auto deposit_or_escape = [&](const double f_p) {
    assert_always(f_p >= 0.);
    assert_always(f_p <= 1.);
    if (rng_uniform(get_rngstate(pkt)) < f_p) {
      pkt.type = deposit_type;
    } else {
      e_cmf_deposited = 0.;
      // change_cell_or_escape() records pkt.escape_type = pkt.type before setting the type to
      // TYPE_ESCAPE, so leave pkt.type as the predeposit particle type here to preserve it.
      grid::change_cell_or_escape(pkt, -99);
    }
  };

  if constexpr (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::INSTANTFULLDEPOSITION) {
    // absorption happens
    pkt.type = deposit_type;
  } else if constexpr (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::BARNES) {
    const double E_kin = grid::get_ejecta_kinetic_energy();
    const double m_ej = grid::get_ejecta_mass();
    const double v_ej = std::sqrt(E_kin * 2 / m_ej);

    const double prefactor = (pkt.type == TYPE_NONTHERMAL_PREDEPOSIT_ALPHA) ? 7.74 : 7.4;
    const double tau_ineff =
        prefactor * DAY * std::sqrt(m_ej / (5.e-3 * MSUN)) * std::pow((0.2 * CLIGHT) / v_ej, 3. / 2.);
    const double f_p = std::log1p(2. * ts * ts / tau_ineff / tau_ineff) / (2. * ts * ts / tau_ineff / tau_ineff);
    deposit_or_escape(f_p);
  } else if constexpr (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::WOLLAEGER) {
    // particle thermalisation from Wollaeger+2018, similar to Barnes but using a slightly different expression
    const double A = (pkt.type == TYPE_NONTHERMAL_PREDEPOSIT_ALPHA) ? 1.2 * 1.e-11 : 1.3 * 1.e-11;
    const double aux_term = 2 * A / (ts * grid::get_rho(nonemptymgi));
    // In Bulla 2023 (arXiv:2211.14348), the following line contains (<-> eq. 7) contains a typo. The way implemented
    // here is the original from Wollaeger paper without the typo
    const double f_p = std::log1p(aux_term) / aux_term;
    deposit_or_escape(f_p);
  } else if constexpr (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::TIMEDEPENDENT ||
                       PARTICLE_THERMALISATION_SCHEME ==
                           ParticleThermalisationScheme::TIMEDEPENDENT_WITH_ADIABATIC_LOSS ||
                       PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::TIMEDEPENDENTWITHGAMMAPRODUCTS) {
    // local time-dependent absorption described by Shingles et al. (2023)
    const double rho = grid::get_rho(nonemptymgi);

    const double particle_en = H * pkt.nu_cmf;  // energy of the particles in the packet

    // the positive energy loss rate per particle [erg/s] from Barnes et al. (2016). see their figure 6.
    const double endot_collisional =
        (pkt.type == TYPE_NONTHERMAL_PREDEPOSIT_ALPHA) ? 5.e11 * MEV * rho : 4.e10 * MEV * rho;
    // positive energy loss rate from adiabatic expansion in [erg/s], assuming homologous expansion
    const double endot_adiabatic =
        (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::TIMEDEPENDENT_WITH_ADIABATIC_LOSS)
            ? particle_en / ts
            : 0.;
    const double endot = endot_collisional + endot_adiabatic;

    // time of deposition is the smaller out of (a) the time until which the particle loses all its energy according to
    // the loss rates above and (b) the time remaining until the end of the current time step
    // Only collisional losses count as energy deposited into the gas (not adiabatic losses)
    e_cmf_deposited = pkt.e_cmf * endot_collisional * std::min(ts_end - ts, particle_en / endot) / particle_en;

    // A discrete absorption event should occur somewhere along the continuous track from initial kinetic energy to zero
    // KE with equal probability of happening at any energy in between. So we can just randomly
    // select an energy at which the absorption happens between particle_en and 0 and then calculate the corresponding
    // time.
    const double rnd_en_absorb = rng_uniform(get_rngstate(pkt)) * particle_en;
    const double t_absorb = ts + (rnd_en_absorb / endot);

    // if absorption happens beyond the end of the current timestep,
    // just reduce the particle energy up to the end of this timestep
    const auto t_new = std::min(t_absorb, ts_end);

    const bool absorbed = (t_absorb <= ts_end);
    if (absorbed) {
      pkt.type = deposit_type;
    } else {
      pkt.nu_cmf -= (endot * (ts_end - ts)) / H;
    }

    pkt.pos = vec_scale(pkt.pos, t_new / ts);
    pkt.prop_time = t_new;
    if constexpr (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::TIMEDEPENDENT_WITH_ADIABATIC_LOSS) {
      if (absorbed) {
        // Only the collisional part of the energy loss heats the gas, so the packet gives up that
        // share of its energy and the remainder goes to the expansion.
        //
        // The packet energy is deliberately not degraded adiabatically while the particle is still
        // travelling. The absorption energy is drawn uniformly, so a packet survives to time t with
        // probability E(t) / E_start, and the expected total of the trajectory estimator above is
        // (endot_collisional / E_start) times the integral of e_cmf over time. The energy that the
        // particles really surrender to the gas is endot_collisional * (t_end - t_start) each, so
        // the two agree only while e_cmf is held constant. The adiabatic loss is already carried by
        // that survival weighting, because it shortens the particle's life; degrading e_cmf as well
        // counted it a second time and made the deposition low by a factor
        // t_start * ln(t_end / t_start) / (t_end - t_start), which falls to 0.42 once the adiabatic
        // losses dominate.
        pkt.e_cmf *= endot_collisional / endot;
      }
    }
    assert_testmodeonly(grid::get_cellindex_from_pos(pkt.pos, pkt.prop_time) == pkt.cellindex);
  } else {
    assert_always(false);  // unhandled thermalisation scheme
  }

  // contribute to the trajectory integrated deposition estimator
  // and if a deposition event occurred, also the discrete Monte Carlo count deposition rate
  // for TIMEDEPENDENTWITHGAMMAPRODUCTS, gamma-ray deposition will lead to predeposit beta particles, but they will
  // count toward "gamma deposition" not particle deposition
  if (pkt.originated_from_particlenotgamma) {
    if (priortype == TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS) {
      atomicadd(globals::dep_estimator_electron[nonemptymgi], e_cmf_deposited);
      if (pkt.type == deposit_type) {
        atomicadd(globals::timesteps[nts].electron_dep_discrete, pkt.e_cmf);
      }
    } else if (priortype == TYPE_NONTHERMAL_PREDEPOSIT_BETAPLUS) {
      atomicadd(globals::dep_estimator_positron[nonemptymgi], e_cmf_deposited);
      if (pkt.type == deposit_type) {
        atomicadd(globals::timesteps[nts].positron_dep_discrete, pkt.e_cmf);
      }
    } else if (priortype == TYPE_NONTHERMAL_PREDEPOSIT_ALPHA) {
      atomicadd(globals::dep_estimator_alpha[nonemptymgi], e_cmf_deposited);
      if (pkt.type == deposit_type) {
        atomicadd(globals::timesteps[nts].alpha_dep_discrete, pkt.e_cmf);
      }
    }
  } else if constexpr (PARTICLE_THERMALISATION_SCHEME == ParticleThermalisationScheme::TIMEDEPENDENTWITHGAMMAPRODUCTS) {
    atomicadd(globals::dep_estimator_gamma[nonemptymgi], e_cmf_deposited);
    if (pkt.type == TYPE_NTLEPTON_DEPOSITED) {
      atomicadd(globals::timesteps[nts].gamma_dep_discrete, pkt.e_cmf);
    }
  }
}

// Handle inactive pellets. Need to do two things (a) check if it
// decays in this time step and if it does handle that. (b) if it doesn't decay in
// this time step then just move the packet along with the matter for the start of the next time step.
void update_pellet(Packet& pkt, const int nts, const double t2) {
  assert_always(pkt.prop_time < t2);
  const double ts = pkt.prop_time;

  const double tdecay = pkt.tdecay;  // after packet_init(), this value never changes
  if (tdecay > t2) {
    // It won't decay in this timestep, so just need to move it on with the flow.
    pkt.pos = vec_scale(pkt.pos, t2 / ts);
    pkt.prop_time = t2;
    assert_testmodeonly(grid::get_cellindex_from_pos(pkt.pos, pkt.prop_time) == pkt.cellindex);

    // That's all that needs to be done for the inactive pellet.
  } else if (tdecay > ts) {
    // The packet decays in the current timestep.
    atomicadd(globals::timesteps[nts].pellet_decays, 1);

    pkt.prop_time = tdecay;
    pkt.pos = vec_scale(pkt.pos, tdecay / ts);
    assert_testmodeonly(grid::get_cellindex_from_pos(pkt.pos, pkt.prop_time) == pkt.cellindex);

    if (pkt.originated_from_particlenotgamma) {
      // decay to non-thermal particle
      if (pkt.pellet_decaytype == decay::DECAYTYPE_BETAPLUS) {
        pkt.type = TYPE_NONTHERMAL_PREDEPOSIT_BETAPLUS;
        atomicadd(globals::timesteps[nts].positron_emission, pkt.e_cmf);
      } else if (pkt.pellet_decaytype == decay::DECAYTYPE_BETAMINUS) {
        pkt.type = TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS;
        atomicadd(globals::timesteps[nts].electron_emission, pkt.e_cmf);
      } else if (pkt.pellet_decaytype == decay::DECAYTYPE_ALPHA) {
        atomicadd(globals::timesteps[nts].alpha_emission, pkt.e_cmf);
        pkt.type = TYPE_NONTHERMAL_PREDEPOSIT_ALPHA;
      } else if (pkt.pellet_decaytype == decay::DECAYTYPE_SPONTFISSION) {
        atomicadd(globals::timesteps[nts].spfission_dep_discrete, pkt.e_cmf);
        pkt.type = TYPE_NTALPHA_FISPROD_DEPOSITED;
      } else if constexpr (TESTMODE) {
        printlnlog(
            "[error] pellet marked as particle emission is for decaytype {} != any of (alpha, beta+, beta-, spfission)",
            pkt.pellet_decaytype);
        assert_always(false);
      } else {
        __builtin_unreachable();
      }
      pkt.em_time = static_cast<float>(pkt.prop_time);
      pkt.absorptiontype = ABSTYPE_PELLET_PARTICLEDECAY;
    } else {
      // decay to gamma-ray packet
      atomicadd(globals::timesteps[nts].gamma_emission, pkt.e_cmf);
      gammapkt::pellet_gamma_decay(pkt);
    }
  } else if ((tdecay > 0) && (nts == 0)) {
    // These are pellets whose decay times were before the first time step. They become pre-k-packets, which
    // do_packet() immediately re-emits as r-packets with a blackbody frequency. The energy is reduced by
    // tdecay / tmin to account for the work done on the ejecta by the trapped radiation between the decay and
    // the start of the simulation (equation 18 of Lucy 2005, as also applied in decay.cc).
    // The position is already set at globals::tmin so don't need to move it. Assume
    // that it is fixed in place from decay to globals::tmin - i.e. short mfp.

    pkt.e_cmf *= tdecay / globals::tmin;
    pkt.type = TYPE_PRE_KPKT;
    pkt.absorptiontype = ABSTYPE_PELLET_BEFORESIMSTART;
    stats::increment(stats::Counter::K_STAT_FROM_EARLIERDECAY);

    pkt.prop_time = globals::tmin;
  } else if constexpr (TESTMODE) {
    printlnlog("[error] Something wrong with decaying pellets. tdecay {:g} ts {:g} (ts + tw) {:g}", tdecay, ts, t2);
    assert_testmodeonly(false);
  } else {
    __builtin_unreachable();
  }
}

// update a packet no further than time t2
void do_packet(Packet& pkt, const double t2, const int nts, ContinuumOpacity& chi_rpkt_cont) {
  switch (pkt.type) {
    case TYPE_RADIOACTIVE_PELLET: {
      update_pellet(pkt, nts, t2);
      break;
    }

    case TYPE_GAMMA: {
      gammapkt::do_gamma(pkt, nts, t2);
      break;
    }

    case TYPE_RPKT: {
      do_rpkt(pkt, t2, chi_rpkt_cont);
      break;
    }

    case TYPE_NONTHERMAL_PREDEPOSIT_ALPHA:
    case TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS:
    case TYPE_NONTHERMAL_PREDEPOSIT_BETAPLUS: {
      do_nonthermal_predeposit(pkt, nts, t2);
      break;
    }

    case TYPE_NTLEPTON_DEPOSITED: {
      nonthermal::do_ntlepton_deposit(pkt);
      break;
    }

    case TYPE_NTALPHA_FISPROD_DEPOSITED: {
      nonthermal::do_ntalpha_fisprod_deposit(pkt);
      break;
    }

    case TYPE_PRE_KPKT: {
      kpkt::do_kpkt_blackbody(pkt);
      break;
    }

    case TYPE_KPKT: {
      const int nonemptymgi = grid::get_propcell_nonemptymgi(pkt.cellindex);
      assert_testmodeonly(nonemptymgi >= 0);  // k-packets can only exist in non-empty cells
      if (grid::thick_allcells[nonemptymgi] == grid::CellThickness::THICK ||
          RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value()) {
        kpkt::do_kpkt_blackbody(pkt);
      } else {
        kpkt::do_kpkt(pkt, t2, nts);
      }
      break;
    }

    default: {
      if constexpr (TESTMODE) {
        printlnlog("[error] Unknown packet type {}", static_cast<int>(pkt.type));
        assert_testmodeonly(false);
      } else {
        __builtin_unreachable();
      }
    }
  }
}

// Return true if the packet still needs propagating this timestep, i.e. it has not escaped and its propagation
// time has not yet reached the end of the timestep.
constexpr auto packetprop_update_required(const Packet& pkt, const double ts_end) -> bool {
  if (pkt.type == TYPE_ESCAPE) {
    return false;
  }
  return pkt.prop_time < ts_end;
}

// Return the id of the cell cache group that this packet belongs to, or an empty std::optional if the
// packet does not use the cell cache at all (pellet/gamma/predeposit types, empty cells, thick cells).
// In single-slot mode the group id is the nonemptymgi, because packets of a given cell must be processed
// together while that cell occupies the rank's one cache slot. In multi-slot mode every cell has its own
// persistent slot, so no partitioning is needed and all cache-using packets share group 0.
auto get_packet_cellcachegroupid(const Packet& pkt) -> std::optional<int> {
  constexpr auto nocache_packettypes = std::array{
      TYPE_RADIOACTIVE_PELLET,
      TYPE_GAMMA,
      TYPE_PRE_KPKT,
      TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS,
      TYPE_NONTHERMAL_PREDEPOSIT_BETAPLUS,
      TYPE_NONTHERMAL_PREDEPOSIT_ALPHA,
      TYPE_NTALPHA_FISPROD_DEPOSITED,
  };
  if (std::ranges::find(nocache_packettypes, pkt.type) != nocache_packettypes.end()) {
    return std::nullopt;  // these types do not use the cell cache
  }
  const auto mgi = grid::get_propcell_modelgridindex(pkt.cellindex);
  if (mgi < 0) {
    return std::nullopt;  // for empty cell, no cell cache required
  }
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
  if (grid::thick_allcells[nonemptymgi] == grid::CellThickness::THICK) {
    return std::nullopt;  // for thick cell, no cell cache required
  }
  if (!cellcache_singleslot) {
    return 0;  // all cell caches are available, so no partitioning is required
  }
  return nonemptymgi;
}

// Strict-weak ordering used to sort packets before propagation: packets still needing updates come first, then
// they are grouped by cell-cache group and ordered by descending cell density and descending nu_cmf so that
// packets sharing a cell (and cache state) are processed consecutively for cache efficiency.
auto compare_packet_order(const Packet& p1, const Packet& p2, const double ts_end) -> bool {
  // return true if packet p1 goes before p2

  // first order by whether the packet has reached the end of the timestep or escaped (both of which mean it won't be
  // updated anymore by update_packets in this timestep)

  const auto pktactive1 = packetprop_update_required(p1, ts_end);
  const auto pktactive2 = packetprop_update_required(p2, ts_end);
  // if one packet is active and the other is not, then the active one goes first
  if (pktactive1 && !pktactive2) {
    return true;
  }
  if (!pktactive1 && pktactive2) {
    return false;
  }
  if (!pktactive1 && !pktactive2) {
    // both are inactive - order doesn't matter
    return false;
  }

  // all packets in empty or thick cells can be grouped together since they don't use the cell cache
  const int cellcachenonemptymgi1 = get_packet_cellcachegroupid(p1).value_or(-1);
  const int cellcachenonemptymgi2 = get_packet_cellcachegroupid(p2).value_or(-1);
  const auto rho1 = cellcachenonemptymgi1 >= 0 ? grid::get_rho(cellcachenonemptymgi1) : 0.0;
  const auto rho2 = cellcachenonemptymgi2 >= 0 ? grid::get_rho(cellcachenonemptymgi2) : 0.0;

  // rho1 and rho2 are swapped here so that higher density cells are first (descending order of density)
  // nu_cmf 1 and 2 are also swapped so that higher energy packets are first (descending order of nu_cmf)
  // other fields are ascending order
  return std::tie(rho2, cellcachenonemptymgi1, p1.type, p2.nu_cmf) <
         std::tie(rho1, cellcachenonemptymgi2, p2.type, p1.nu_cmf);
}

// fill the cellcache with values for the current cell
void cellcacheslot_populate(globals::CellCache& cacheslot, const int nonemptymgi) {
  assert_always(nonemptymgi >= 0);
  stats::increment(stats::Counter::UPDATECELL);

  cacheslot.nonemptymgi = nonemptymgi;
  cacheslot.chi_ff_nnionpart[0] = calculate_chi_ffheat_nnionpart(nonemptymgi);

  const int nelements = get_nelements();
  for (int element = 0; element < nelements; element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      cacheslot.cooling_contrib[kpkt::get_coolinglistoffset(element, ion)] = -99.;
    }

    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels(element, ion);
      const auto uniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int level = 0; level < nlevels; level++) {
        cacheslot.alllevels_pops[uniquelevelindexstart + level] = calculate_levelpop(nonemptymgi, element, ion, level);
      }
    }
  }

  std::ranges::fill(cacheslot.allphixstargets_corrphotoioncoeff, -99.);

  std::ranges::fill(cacheslot.allcont_modified_departureratios, -1.);
  std::ranges::fill(cacheslot.allcont_stimfactor_edgepart, -1.);

  std::ranges::fill(cacheslot.allcont_keepbits, 0);

  const auto nnetot = grid::get_nnetot(nonemptymgi);
  for (int i = 0; i < globals::nbfcontinua; i++) {
    const auto nnlevel = cacheslot.alllevels_pops[globals::allcont.uniquelevelindex[i]];
    // an unpopulated level is skipped either way, so skip keep_this_cont's division for those
    const bool keep = nnlevel > 0 && keep_this_cont(globals::allcont.element[i], globals::allcont.ion[i],
                                                    globals::allcont.level[i], nonemptymgi, nnetot);
    cacheslot.allcont_nnlevel[i] = nnlevel;
    if (keep) {
      set_allcont_keepbit(cacheslot.allcont_keepbits, i);
    }
  }

  if constexpr (!cellcache_singleslot) {
    // The cache slot is shared between node ranks and persists for the whole timestep, so the macroatom
    // transition rates and k-packet cooling rates cannot be filled in lazily under a per-cell mutex during
    // propagation (that would race across ranks, and on the GPU there is no mutex at all). Instead
    // pre-calculate every cell's rates here, up front. The populate calls are distributed across the node's
    // ranks (see update_packets), so each rank only does this for a subset of cells.
    const double t_mid = globals::timesteps[globals::timestep].mid;
    const auto alllevelindices = std::ranges::iota_view{0, get_includedlevels()};
    std::for_each(EXEC_PAR_UNSEQ alllevelindices.begin(), alllevelindices.end(),
                  [nonemptymgi, t_mid](const int uniquelevelindex) {
                    calculate_cellcache_macroatom_transitionrates(nonemptymgi, uniquelevelindex, t_mid);
                  });

    const auto allionindices = std::ranges::iota_view{0, get_includedions()};
    std::for_each(EXEC_PAR_UNSEQ allionindices.begin(), allionindices.end(), [nonemptymgi](const int uniqueionindex) {
      kpkt::calculate_cellcache_cooling_rates_ion(nonemptymgi, uniqueionindex);
    });
  } else {
    for (int uniquelevelindex = 0; uniquelevelindex < get_includedlevels(); uniquelevelindex++) {
      cacheslot.alllevels_maprocessrates[uniquelevelindex * MA_ACTION_COUNT] = -99.;
    }
  }
}

// Propagate a group of packets that share a cell-cache group until the end of the timestep, (re)populating the
// single shared cache slot for that group's cell first when cellcache_singleslot is enabled.
void update_packet_cellcache_group(const int cellcache_groupid, std::span<Packet> packetgroup, const int nts,
                                   const double ts_end) {
  if (cellcache_singleslot) {
    // in this case, a positive groupid is a nonemptymgi, and -1 is the no-cache-required group
    auto& cacheslot = globals::cellcache[globals::rank_in_node];
    if (cellcache_groupid >= 0 && cacheslot.nonemptymgi != cellcache_groupid) {
      cellcacheslot_populate(cacheslot, cellcache_groupid);
    }
  }

#ifdef GPU_ON
  // we don't know how many GPU threads will exist, and we can't use a thread_local variables on device.
  // So, we assume the worst case that each packet is handled simultaneously by different GPU threads.
  static std::vector<ContinuumOpacity> chi_rpkt_cont_vec;
  if (cellcache_groupid >= 0) {
    if (chi_rpkt_cont_vec.size() < packetgroup.size()) {
      const auto size_mb = (static_cast<ptrdiff_t>(globals::nbfcontinua_ground) +
                            static_cast<ptrdiff_t>(globals::nbfcontinua) + std::ssize(globals::bfestim_nu_edge)) *
                           std::ssize(packetgroup) * sizeof(double) / 1024. / 1024.;
      printlog("Resizing chi_rpkt_cont_vec from {} to {} ({:g} MB)...", chi_rpkt_cont_vec.size(), packetgroup.size(),
               size_mb);
      chi_rpkt_cont_vec.resize(packetgroup.size());
      printlnlog("done.");
    }
  } else if (chi_rpkt_cont_vec.empty()) {
    // we're not going to use this, but we need to pass a reference to something
    chi_rpkt_cont_vec.resize(1);
  }
#endif

  const auto update_packet = [cellcache_groupid, ts_end, nts, &packetgroup](const ptrdiff_t index_in_group) {
#ifdef GPU_ON
    auto& chi_rpkt_cont = chi_rpkt_cont_vec[(cellcache_groupid >= 0) ? index_in_group : 0];
#else
    // thread_local lets us reuse this allocation on each CPU thread
    THREADLOCALONHOST auto chi_rpkt_cont = ContinuumOpacity{};
#endif
    auto& pkt = packetgroup[index_in_group];
    while (packetprop_update_required(pkt, ts_end) &&
           (get_packet_cellcachegroupid(pkt).value_or(cellcache_groupid) == cellcache_groupid)) {
      do_packet(pkt, ts_end, nts, chi_rpkt_cont);
    }
  };

#if defined(STDPAR_ON) || !defined(_OPENMP)
  const auto pktgroupindices = std::ranges::iota_view{0Z, std::ssize(packetgroup)};
  std::for_each(EXEC_PAR pktgroupindices.begin(), pktgroupindices.end(), update_packet);
#else
#ifdef GPU_ON
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(nonmonotonic : dynamic)
#endif
  for (auto i = 0Z; i < std::ssize(packetgroup); i++) {
    update_packet(i);
  }
#endif
}

}  // anonymous namespace

// Move and update packets during the current timestep (nts)
void update_packets(const int nts, std::span<Packet> packets) {
  // At the start, the packets have all either just been initialised or have already been
  // processed for one or more timesteps. Those that are pellets will just be sitting in the
  // matter. Those that are photons (or one sort or another) will already have a position and a direction.
  const double ts = globals::timesteps[nts].start;
  const double tw = globals::timesteps[nts].width;
  const double ts_end = ts + tw;

#ifdef GPU_ON
  // to avoid failed allocation of chi_rpkt_cont_vec on GPU
  // ~1e6 allocation was 155GB for classic model, so limit to 1e5 packets per group for now, which is ~15GB on GPU
  constexpr auto MAX_PACKETGROUP_SIZE = 100'000Z;
#else
  constexpr auto MAX_PACKETGROUP_SIZE = static_cast<ptrdiff_t>(std::numeric_limits<int>::max());
#endif

  const auto time_update_packets_start = std::chrono::steady_clock::now();
  printlnlog("timestep {}: start update_packets", nts);
  if constexpr (cellcache_singleslot) {
    // invalidate the cell cache, since it might match the nonemptymgi but with old values from the previous timestep
    globals::cellcache[globals::rank_in_node].nonemptymgi = -2;
  } else {
    // The cell cache is shared between the node's MPI ranks, so distribute the populate (and up-front rate
    // pre-calculation) across those ranks. Each rank fills the cells in its chunk of the shared memory and
    // then we wait at a node barrier so that every rank sees the fully-populated cache before propagation.
    const auto nonempty_npts_model = grid::get_nonempty_npts_model();
    const auto [nonemptymgi_start, nonemptymgi_count] =
        get_range_chunk(nonempty_npts_model, globals::node_nprocs, globals::rank_in_node);
    for (auto nonemptymgi = nonemptymgi_start; nonemptymgi < (nonemptymgi_start + nonemptymgi_count); nonemptymgi++) {
      cellcacheslot_populate(globals::cellcache.at(nonemptymgi), static_cast<int>(nonemptymgi));
    }
    MPI_Barrier_node();
    printlnlog("timestep {}: all {} cellcaches set (took {:.1f} s)", nts, nonempty_npts_model,
               std::chrono::duration<double>(std::chrono::steady_clock::now() - time_update_packets_start).count());
  }
  int prev_cellcache_groupid = -2;
  int passnumber = 0;
  while (true) {
    const auto sys_time_start_pass = std::chrono::steady_clock::now();

    std::SORT_OR_STABLE_SORT(
        EXEC_PAR_UNSEQ packets.begin(), packets.end(),
        [ts_end](const Packet& p1, const Packet& p2) { return compare_packet_order(p1, p2, ts_end); });

    static std::vector<std::tuple<int, std::span<Packet>>> packet_groups;
    packet_groups.clear();
    auto pass_packets_updated{0Z};
    const auto updatecellcounter_beforepass = stats::get_counter(stats::Counter::UPDATECELL);
    auto packetgroupstart{0Z};
    // because of the sort, we don't need to check the entire packet list
    auto pktindex = 0Z;
    for (; pktindex < std::ssize(packets); pktindex++) {
      const auto& pkt = packets[pktindex];
      if (!packetprop_update_required(pkt, ts_end)) {
        // due to the sorting, all following packets will also not require updating, so can break out of the loop
        break;
      }
      const auto cellcache_groupid = get_packet_cellcachegroupid(pkt).value_or(-1);
      const auto packetgroupsize = pktindex - packetgroupstart;  // not counting the current packet, which would be in
                                                                 // the next group if a new group is started
      if (packetgroupsize > 0 &&
          (cellcache_groupid != prev_cellcache_groupid || packetgroupsize >= MAX_PACKETGROUP_SIZE)) {
        packet_groups.emplace_back(prev_cellcache_groupid, packets.subspan(packetgroupstart, packetgroupsize));
        packetgroupstart = pktindex;
      }
      prev_cellcache_groupid = cellcache_groupid;
    }
    if (packetgroupstart != pktindex) {
      // finish the last group of packets that needed updating
      packet_groups.emplace_back(prev_cellcache_groupid,
                                 packets.subspan(packetgroupstart, pktindex - packetgroupstart));
    }
    if (packet_groups.empty()) {
      // if no packets needed updating, then this timestep is complete
      break;
    }

    // process the packets grouped by their required cell cache, which should minimise the number of times we need to
    // change the cell cache during the packet updates
    for (auto [cellcache_groupid, grouppackets] : packet_groups) {
      update_packet_cellcache_group(cellcache_groupid, grouppackets, nts, ts_end);
      pass_packets_updated += std::ssize(grouppackets);
    }

    const auto cellcacheresets = stats::get_counter(stats::Counter::UPDATECELL) - updatecellcounter_beforepass;
    const auto pass_duration_s =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_pass).count();
    printlnlog("  update_packets timestep {} pass {:3d}: packetsupdated {:7d} cellcacheresets {:7d} (took {:.1f}s)",
               nts, passnumber, pass_packets_updated, cellcacheresets, pass_duration_s);

    passnumber++;
  }

  stats::pkt_action_counters_printout(nts);

  const auto time_update_packets_end_thisrank = std::chrono::steady_clock::now();
  const auto rank_process_time =
      std::chrono::duration<double>(time_update_packets_end_thisrank - time_update_packets_start).count();
  printlnlog("timestep {}: finished update_packets for rank {} (took {:.1f} seconds)", nts, globals::my_rank,
             rank_process_time);

  MPI_Barrier_allranks();  // hold all processes once the packets are updated
  const auto time_update_packets_end_allranks = std::chrono::steady_clock::now();
  const auto rank_wait_time =
      std::chrono::duration<double>(time_update_packets_end_allranks - time_update_packets_end_thisrank).count();
  const auto rank_total_time =
      std::chrono::duration<double>(time_update_packets_end_allranks - time_update_packets_start).count();
  printlnlog(
      "timestep {}: time after update packets for all processes (rank {} took {:.1f}s, waited {:.1f}s, total {:.1f}s)",
      nts, globals::my_rank, rank_process_time, rank_wait_time, rank_total_time);
}
