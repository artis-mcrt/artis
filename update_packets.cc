#include "update_packets.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <optional>
#include <span>
#include <tuple>
#include <vector>

#include "atomic.h"
#include "ltepop.h"

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "constants.h"
#include "decay.h"
#include "gammapkt.h"
#include "globals.h"
#include "grid.h"
#include "kpkt.h"
#include "nonthermal.h"
#include "packet.h"
#include "random.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "vectors.h"

namespace {

void do_nonthermal_predeposit(Packet& pkt, const int nts, const double t2) {
  double en_deposited = pkt.e_cmf;
  const auto mgi = grid::get_propcell_modelgridindex(pkt.where);
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
  const auto priortype = pkt.type;
  const double ts = pkt.prop_time;
  const auto deposit_type =
      (pkt.type == TYPE_NONTHERMAL_PREDEPOSIT_ALPHA) ? TYPE_NTALPHA_FISPROD_DEPOSITED : TYPE_NTLEPTON_DEPOSITED;

  if constexpr (PARTICLE_THERMALISATION_SCHEME == ThermalisationScheme::INSTANT) {
    // absorption happens
    pkt.type = deposit_type;
  } else if constexpr (PARTICLE_THERMALISATION_SCHEME == ThermalisationScheme::BARNES) {
    const double E_kin = grid::get_ejecta_kinetic_energy();
    const double v_ej = std::sqrt(E_kin * 2 / grid::mtot_input);

    const double prefactor = (pkt.type == TYPE_NONTHERMAL_PREDEPOSIT_ALPHA) ? 7.74 : 7.4;
    const double tau_ineff = prefactor * 86400 * std::sqrt(grid::mtot_input / (5.e-3 * 1.989 * 1.e33)) *
                             std::pow((0.2 * 29979200000) / v_ej, 3. / 2.);
    const double f_p = std::log1p(2. * ts * ts / tau_ineff / tau_ineff) / (2. * ts * ts / tau_ineff / tau_ineff);
    assert_always(f_p >= 0.);
    assert_always(f_p <= 1.);
    if (rng_uniform() < f_p) {
      pkt.type = deposit_type;
    } else {
      en_deposited = 0.;
      pkt.type = TYPE_ESCAPE;
      grid::change_cell(pkt, -99);
    }
  } else if constexpr (PARTICLE_THERMALISATION_SCHEME == ThermalisationScheme::WOLLAEGER) {
    // particle thermalisation from Wollaeger+2018, similar to Barnes but using a slightly different expression
    const double A = (pkt.type == TYPE_NONTHERMAL_PREDEPOSIT_ALPHA) ? 1.2 * 1.e-11 : 1.3 * 1.e-11;
    const double aux_term = 2 * A / (ts * grid::get_rho(nonemptymgi));
    // In Bulla 2023 (arXiv:2211.14348), the following line contains (<-> eq. 7) contains a typo. The way implemented
    // here is the original from Wollaeger paper without the typo
    const double f_p = std::log1p(aux_term) / aux_term;
    assert_always(f_p >= 0.);
    assert_always(f_p <= 1.);
    if (rng_uniform() < f_p) {
      pkt.type = deposit_type;
    } else {
      en_deposited = 0.;
      pkt.type = TYPE_ESCAPE;
      grid::change_cell(pkt, -99);
    }
  } else {
    // ThermalisationScheme::DETAILED or ThermalisationScheme::DETAILEDWITHGAMMAPRODUCTS
    // local, detailed absorption following Shingles+2023
    const double rho = grid::get_rho(nonemptymgi);

    // endot is energy loss rate (positive) in [erg/s]
    // endot [erg/s] from Barnes et al. (2016). see their figure 6.
    const double endot = (pkt.type == TYPE_NONTHERMAL_PREDEPOSIT_ALPHA) ? 5.e11 * MEV * rho : 4.e10 * MEV * rho;

    const double particle_en = H * pkt.nu_cmf;  // energy of the particles in the packet

    // for endot independent of energy, the next line is trivial (for E dependent endot, an integral would be needed)

    const double t_enzero = ts + (particle_en / endot);  // time at which zero energy is reached
    en_deposited = pkt.e_cmf * (std::min(t2, t_enzero) - ts) / (particle_en / endot);

    // A discrete absorption event should occur somewhere along the
    // continuous track from initial kinetic energy to zero KE.
    // The probability of being absorbed in energy range [E, E+delta_E] is proportional to
    // endot(E) * delta_t = endot(E) * delta_E / endot(E) = delta_E (delta_t is the time spent in the bin range)
    // so all final energies are equally likely.
    // Choose random en_absorb [0, particle_en]

    const double rnd_en_absorb = rng_uniform() * particle_en;
    const double t_absorb = ts + (rnd_en_absorb / endot);

    // if absorption happens beyond the end of the current timestep,
    // just reduce the particle energy up to the end of this timestep
    const auto t_new = std::min(t_absorb, t2);

    if (t_absorb <= t2) {
      pkt.type = deposit_type;
    } else {
      pkt.nu_cmf = (particle_en - (endot * (t_new - ts))) / H;
    }

    pkt.pos = vec_scale(pkt.pos, t_new / ts);
    pkt.prop_time = t_new;
    // pkt.e_cmf *= ts / t_new;
    assert_testmodeonly(grid::get_cellindex_from_pos(pkt.pos, pkt.prop_time) == pkt.where);
  }

  // contribute to the trajectory integrated deposition estimator
  // and if a deposition event occurred, also the discrete Monte Carlo count deposition rate
  // for DETAILEDWITHGAMMAPRODUCTS, gamma-ray deposition will lead to predeposit beta particles, but they will count
  // toward "gamma deposition" not particle deposition
  if (pkt.originated_from_particlenotgamma) {
    if (priortype == TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS) {
      atomicadd(globals::dep_estimator_electron[nonemptymgi], en_deposited);
      if (pkt.type == deposit_type) {
        atomicadd(globals::timesteps[nts].electron_dep_discrete, pkt.e_cmf);
      }
    } else if (priortype == TYPE_NONTHERMAL_PREDEPOSIT_BETAPLUS) {
      atomicadd(globals::dep_estimator_positron[nonemptymgi], en_deposited);
      if (pkt.type == deposit_type) {
        atomicadd(globals::timesteps[nts].positron_dep_discrete, pkt.e_cmf);
      }
    } else if (priortype == TYPE_NONTHERMAL_PREDEPOSIT_ALPHA) {
      atomicadd(globals::dep_estimator_alpha[nonemptymgi], en_deposited);
      if (pkt.type == deposit_type) {
        atomicadd(globals::timesteps[nts].alpha_dep_discrete, pkt.e_cmf);
      }

    } else if constexpr (PARTICLE_THERMALISATION_SCHEME == ThermalisationScheme::DETAILEDWITHGAMMAPRODUCTS) {
      atomicadd(globals::dep_estimator_gamma[nonemptymgi], en_deposited);
      if (pkt.type == TYPE_NTLEPTON_DEPOSITED) {
        atomicadd(globals::timesteps[nts].gamma_dep_discrete, pkt.e_cmf);
      }
    }
  }
}

// Handle inactive pellets. Need to do two things (a) check if it
// decays in this time step and if it does handle that. (b) if it doesn't decay in
// this time step then just move the packet along with the matter for the
// start of the next time step.
void update_pellet(Packet& pkt, const int nts, const double t2) {
  assert_always(pkt.prop_time < t2);
  const double ts = pkt.prop_time;

  const double tdecay = pkt.tdecay;  // after packet_init(), this value never changes
  if (tdecay > t2) {
    // It won't decay in this timestep, so just need to move it on with the flow.
    pkt.pos = vec_scale(pkt.pos, t2 / ts);
    pkt.prop_time = t2;
    assert_testmodeonly(grid::get_cellindex_from_pos(pkt.pos, pkt.prop_time) == pkt.where);

    // That's all that needs to be done for the inactive pellet.
  } else if (tdecay > ts) {
    // The packet decays in the current timestep.
    atomicadd(globals::timesteps[nts].pellet_decays, 1);

    pkt.prop_time = tdecay;
    pkt.pos = vec_scale(pkt.pos, tdecay / ts);
    assert_testmodeonly(grid::get_cellindex_from_pos(pkt.pos, pkt.prop_time) == pkt.where);

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
        assert_testmodeonly(DECAY_SPONTFISSION_ON);
        atomicadd(globals::timesteps[nts].spfission_dep_discrete, pkt.e_cmf);
        pkt.type = TYPE_NTALPHA_FISPROD_DEPOSITED;
      } else if constexpr (TESTMODE) {
        printlnlog(
            "ERROR: pellet marked as particle emission is for decaytype {} != any of (alpha, beta+, beta-, spfission)",
            pkt.pellet_decaytype);
        std::abort();
      } else {
        __builtin_unreachable();
      }
      pkt.em_time = static_cast<float>(pkt.prop_time);
      pkt.absorptiontype = -10;
    } else {
      // decay to gamma-ray packet
      atomicadd(globals::timesteps[nts].gamma_emission, pkt.e_cmf);
      gammapkt::pellet_gamma_decay(pkt);
    }
  } else if ((tdecay > 0) && (nts == 0)) {
    // These are pellets whose decay times were before the first time step
    // They will be made into r-packets with energy reduced for doing work on the
    // ejecta following Lucy 2004.
    // The position is already set at globals::tmin so don't need to move it. Assume
    // that it is fixed in place from decay to globals::tmin - i.e. short mfp.

    pkt.e_cmf *= tdecay / globals::tmin;
    pkt.type = TYPE_PRE_KPKT;
    pkt.absorptiontype = -7;
    stats::increment(stats::Counter::K_STAT_FROM_EARLIERDECAY);

    pkt.prop_time = globals::tmin;
  } else if constexpr (TESTMODE) {
    printlnlog("ERROR: Something wrong with decaying pellets. tdecay {:g} ts {:g} (ts + tw) {:g}", tdecay, ts, t2);
    assert_testmodeonly(false);
  } else {
    __builtin_unreachable();
  }
}

// update a packet no further than time t2
void do_packet(Packet& pkt, const double t2, const int nts) {
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
      do_rpkt(pkt, t2);
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
      const int mgi = grid::get_propcell_modelgridindex(pkt.where);
      const int nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
      if (grid::thick_allcells[nonemptymgi] == 1 ||
          (EXPANSIONOPACITIES_ON && RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY.has_value())) {
        kpkt::do_kpkt_blackbody(pkt);
      } else {
        kpkt::do_kpkt(pkt, t2, nts);
      }
      break;
    }

    default: {
      if constexpr (TESTMODE) {
        printlnlog("ERROR: Unknown packet type {}", static_cast<int>(pkt.type));
        assert_testmodeonly(false);
      } else {
        __builtin_unreachable();
      }
    }
  }
}

constexpr auto packetprop_update_required(const Packet& pkt, const double ts_end) -> bool {
  if (pkt.type == TYPE_ESCAPE) {
    return false;
  }
  return pkt.prop_time < ts_end;
}

// Return the nonemptymgi for the cell cache if required (non-empty, non-thick cell),
// otherwise return an empty std::optional to indicate that no cell cache is used
auto get_packet_cellcachenonemptymgi(const Packet& pkt) -> std::optional<int> {
  constexpr auto nocache_packettypes = std::array<packet_type, 7>{TYPE_RADIOACTIVE_PELLET,
                                                                  TYPE_GAMMA,
                                                                  TYPE_PRE_KPKT,
                                                                  TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS,
                                                                  TYPE_NONTHERMAL_PREDEPOSIT_BETAPLUS,
                                                                  TYPE_NONTHERMAL_PREDEPOSIT_ALPHA,
                                                                  TYPE_NTALPHA_FISPROD_DEPOSITED};
  if (std::ranges::find(nocache_packettypes, pkt.type) != nocache_packettypes.end()) {
    return {};  // these types do not use the cell cache
  }
  const auto mgi = grid::get_propcell_modelgridindex(pkt.where);
  if (mgi < 0) {
    return {};  // for empty cell, no cell cache required
  }
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
  if (grid::thick_allcells[nonemptymgi] == 1) {
    return {};  // for thick cell, no cell cache required
  }
  return {nonemptymgi};
}

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
  const int cellcachenonemptymgi1 = get_packet_cellcachenonemptymgi(p1).value_or(-1);
  const int cellcachenonemptymgi2 = get_packet_cellcachenonemptymgi(p2).value_or(-1);
  const auto rho1 = cellcachenonemptymgi1 >= 0 ? grid::get_rho(cellcachenonemptymgi1) : 0.0;
  const auto rho2 = cellcachenonemptymgi2 >= 0 ? grid::get_rho(cellcachenonemptymgi2) : 0.0;

  // rho 1 and 2 are swapped here since we want higher density cells to come first
  return std::tie(rho2, cellcachenonemptymgi1, p1.type, p2.nu_cmf) <
         std::tie(rho1, cellcachenonemptymgi2, p2.type, p1.nu_cmf);
}

// fill the cellcache with values for the current cell
void cellcache_change_cell(globals::CellCache& cacheslot, const int nonemptymgi) {
  assert_always(nonemptymgi >= 0);
  stats::increment(stats::Counter::UPDATECELL);

  cacheslot.nonemptymgi = nonemptymgi;
  cacheslot.chi_ff_nnionpart = -1.;

  const int nelements = get_nelements();
  for (int element = 0; element < nelements; element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      cacheslot.cooling_contrib[kpkt::get_coolinglistoffset(element, ion)] = COOLING_UNDEFINED;
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

  for (int uniquelevelindex = 0; uniquelevelindex < std::ssize(cacheslot.alllevels_maprocessrates);
       uniquelevelindex++) {
    cacheslot.alllevels_maprocessrates[uniquelevelindex][0] = -99.;
  }

  std::ranges::fill(cacheslot.allcont_departureratios, -1.);

  const auto nnetot = grid::get_nnetot(nonemptymgi);
  for (int i = 0; i < globals::nbfcontinua; i++) {
    const auto nnlevel = cacheslot.alllevels_pops[globals::allcont.uniquelevelindex[i]];
    cacheslot.allcont_nnlevel[i] = nnlevel;
    cacheslot.allcont_keep[i] = nnlevel > 0 && keep_this_cont(globals::allcont.element[i], globals::allcont.ion[i],
                                                              globals::allcont.level[i], nonemptymgi, nnetot);
  }
}

void update_packet_cellcache_group(const int cellcache_nonemptymgi, std::span<Packet> packets, const int nts,
                                   const double ts_end) {
  if (cellcache_nonemptymgi >= 0 && globals::cellcache[cellcacheslotid].nonemptymgi != cellcache_nonemptymgi) {
    cellcache_change_cell(globals::cellcache[cellcacheslotid], cellcache_nonemptymgi);
  }

  auto update_packet = [cellcache_nonemptymgi, ts_end, nts](auto& pkt) {
    while (packetprop_update_required(pkt, ts_end) &&
           (get_packet_cellcachenonemptymgi(pkt).value_or(cellcache_nonemptymgi) == cellcache_nonemptymgi)) {
      do_packet(pkt, ts_end, nts);
    }
  };

#if defined(STDPAR_ON) || !defined(_OPENMP)
  std::for_each(EXEC_PAR packets.begin(), packets.end(), update_packet);
#else
#ifdef GPU_ON
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(nonmonotonic : dynamic)
#endif
  for (auto i = 0Z; i < std::ssize(packets); i++) {
    update_packet(packets[i]);
  }
#endif
}

}  // anonymous namespace

// Move and update packets during the current timestep (nts)
void update_packets(const int nts, std::span<Packet> packets) {
  // At the start, the packets have all either just been initialised or have already been
  // processed for one or more timesteps. Those that are pellets will just be sitting in the
  // matter. Those that are photons (or one sort or another) will already have a position and
  // a direction.
  const double ts = globals::timesteps[nts].start;
  const double tw = globals::timesteps[nts].width;
  const double ts_end = ts + tw;

  const auto time_update_packets_start = std::time(nullptr);
  printlnlog("timestep {}: start update_packets at time {}", nts, time_update_packets_start);
  // first group will probably be the -1 no-cache required group, so -2 triggers the first update
  globals::cellcache[cellcacheslotid].nonemptymgi = -2;
  int prevpkt_cellcache_nonemptymgi = -2;
  int passnumber = 0;
  while (true) {
    const auto sys_time_start_pass = std::time(nullptr);

    std::ranges::SORT_OR_STABLE_SORT(
        packets, [ts_end](const Packet& p1, const Packet& p2) { return compare_packet_order(p1, p2, ts_end); });

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
      const auto cellcache_nonemptymgi = get_packet_cellcachenonemptymgi(pkt).value_or(-1);
      if (cellcache_nonemptymgi != prevpkt_cellcache_nonemptymgi && packetgroupstart != pktindex) {
        packet_groups.emplace_back(prevpkt_cellcache_nonemptymgi,
                                   packets.subspan(packetgroupstart, pktindex - packetgroupstart));
        packetgroupstart = pktindex;
      }
      prevpkt_cellcache_nonemptymgi = cellcache_nonemptymgi;
    }
    if (packetgroupstart != pktindex) {
      // finish the last group of packets that needed updating
      packet_groups.emplace_back(prevpkt_cellcache_nonemptymgi,
                                 packets.subspan(packetgroupstart, pktindex - packetgroupstart));
    }
    if (packet_groups.empty()) {
      // if no packets needed updating, then this timestep is complete
      break;
    }

    // process the packets grouped by their required cell cache, which should minimise the number of times we need to
    // change the cell cache during the packet updates
    for (auto [cellcache_nonemptymgi, grouppackets] : packet_groups) {
      update_packet_cellcache_group(cellcache_nonemptymgi, grouppackets, nts, ts_end);
      pass_packets_updated += std::ssize(grouppackets);
    }

    const auto cellcacheresets = stats::get_counter(stats::Counter::UPDATECELL) - updatecellcounter_beforepass;
    printlnlog("  update_packets timestep {} pass {:3d}: packetsupdated {:7d} cellcacheresets {:7d} (took {}s)", nts,
               passnumber, pass_packets_updated, cellcacheresets, std::time(nullptr) - sys_time_start_pass);

    passnumber++;
  }

  stats::pkt_action_counters_printout(nts);

  const auto time_update_packets_end_thisrank = std::time(nullptr);
  printlnlog("timestep {}: finished update_packets for rank {} (took {} seconds)", nts, globals::my_rank,
             time_update_packets_end_thisrank - time_update_packets_start);

  MPI_Barrier(MPI_COMM_WORLD);  // hold all processes once the packets are updated
  const auto time_update_packets_end_allranks = std::time(nullptr);
  printlnlog("timestep {}: time after update packets for all processes (rank {} took {}s, waited {}s, total {}s)", nts,
             globals::my_rank, time_update_packets_end_thisrank - time_update_packets_start,
             time_update_packets_end_allranks - time_update_packets_end_thisrank,
             time_update_packets_end_allranks - time_update_packets_start);
}
