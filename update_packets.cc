#include "update_packets.h"

#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iterator>
#include <span>

#include "artisoptions.h"
#include "constants.h"
#include "decay.h"
#include "gammapkt.h"
#include "globals.h"
#include "grid.h"
#include "kpkt.h"
#ifdef MPI_ON
#include "mpi.h"
#endif
#include "nonthermal.h"
#include "packet.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "update_grid.h"
#include "vectors.h"

static void do_nonthermal_predeposit(Packet &pkt, const int nts, const double t2) {
  double en_deposited = pkt.e_cmf;

  if constexpr (INSTANT_PARTICLE_DEPOSITION) {
    // absorption happens
    pkt.type = TYPE_NTLEPTON;
  } else {
    const double rho = grid::get_rho(grid::get_cell_modelgridindex(pkt.where));

    // endot is energy loss rate (positive) in [erg/s]
    // endot [erg/s] from Barnes et al. (2016). see their figure 6.
    const double endot = (pkt.pellet_decaytype == decay::DECAYTYPE_ALPHA) ? 5.e11 * MEV * rho : 4.e10 * MEV * rho;

    const double ts = pkt.prop_time;
    const double particle_en = H * pkt.nu_cmf;  // energy of the particles in the packet

    // for endot independent of energy, the next line is trival (for E dependent endot, an integral would be needed)

    const double t_enzero = ts + particle_en / endot;  // time at which zero energy is reached
    if (t_enzero > t2) {
      en_deposited = pkt.e_cmf * (t2 - ts) / (particle_en / endot);
    } else {
      en_deposited = pkt.e_cmf * (t_enzero - ts) / (particle_en / endot);
    }

    // A discrete absorption event should occur somewhere along the
    // continuous track from initial kinetic energy to zero KE.
    // The probability of being absorbed in energy range [E, E+delta_E] is proportional to
    // endot(E) * delta_t = endot(E) * delta_E / endot(E) = delta_E (delta_t is the time spent in the bin range)
    // so all final energies are equally likely.
    // Choose random en_absorb [0, particle_en]

    const double rnd_en_absorb = rng_uniform() * particle_en;
    const double t_absorb = ts + rnd_en_absorb / endot;

    // if absorption happens beyond the end of the current timestep,
    // just reduce the particle energy up to the end of this timestep
    const auto t_new = std::min(t_absorb, t2);

    if (t_absorb <= t2) {
      pkt.type = TYPE_NTLEPTON;
    } else {
      pkt.nu_cmf = (particle_en - endot * (t_new - ts)) / H;
    }

    vec_scale(pkt.pos, t_new / ts);
    pkt.prop_time = t_new;
  }

  if (pkt.pellet_decaytype == decay::DECAYTYPE_ALPHA) {
    safeadd(globals::timesteps[nts].alpha_dep, en_deposited);
  } else if (pkt.pellet_decaytype == decay::DECAYTYPE_BETAMINUS) {
    safeadd(globals::timesteps[nts].electron_dep, en_deposited);
  } else if (pkt.pellet_decaytype == decay::DECAYTYPE_BETAPLUS) {
    safeadd(globals::timesteps[nts].positron_dep, en_deposited);
  }
}

static void update_pellet(Packet &pkt, const int nts, const double t2) {
  // Handle inactive pellets. Need to do two things (a) check if it
  // decays in this time step and if it does handle that. (b) if it doesn't decay in
  // this time step then just move the packet along with the matter for the
  // start of the next time step.
  assert_always(pkt.prop_time < t2);
  const double ts = pkt.prop_time;

  const double tdecay = pkt.tdecay;  // after packet_init(), this value never changes
  if (tdecay > t2) {
    // It won't decay in this timestep, so just need to move it on with the flow.
    vec_scale(pkt.pos, t2 / ts);
    pkt.prop_time = t2;

    // That's all that needs to be done for the inactive pellet.
  } else if (tdecay > ts) {
    // The packet decays in the current timestep.
    globals::timesteps[nts].pellet_decays++;

    pkt.prop_time = tdecay;
    vec_scale(pkt.pos, tdecay / ts);

    if (pkt.originated_from_particlenotgamma)  // will decay to non-thermal particle
    {
      if (pkt.pellet_decaytype == decay::DECAYTYPE_BETAPLUS) {
        safeadd(globals::timesteps[nts].positron_dep, pkt.e_cmf);
        pkt.type = TYPE_NTLEPTON;
        pkt.absorptiontype = -10;
      } else if (pkt.pellet_decaytype == decay::DECAYTYPE_BETAMINUS) {
        safeadd(globals::timesteps[nts].electron_emission, pkt.e_cmf);
        pkt.em_time = pkt.prop_time;
        pkt.type = TYPE_NONTHERMAL_PREDEPOSIT;
        pkt.absorptiontype = -10;
      } else if (pkt.pellet_decaytype == decay::DECAYTYPE_ALPHA) {
        safeadd(globals::timesteps[nts].alpha_emission, pkt.e_cmf);
        pkt.em_time = pkt.prop_time;
        pkt.type = TYPE_NONTHERMAL_PREDEPOSIT;
        pkt.absorptiontype = -10;
      }
    } else {
      safeadd(globals::timesteps[nts].gamma_emission, pkt.e_cmf);
      // decay to gamma-ray, kpkt, or ntlepton
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
    stats::increment(stats::COUNTER_K_STAT_FROM_EARLIERDECAY);

    // printout("already decayed packets and propagation by packet_prop\n");
    pkt.prop_time = globals::tmin;
  } else {
    printout("ERROR: Something gone wrong with decaying pellets. tdecay %g ts %g (ts + tw) %g\n", tdecay, ts, t2);
    assert_testmodeonly(false);
    std::abort();
  }
}

static void do_packet(Packet &pkt, const double t2, const int nts)
// update a packet no further than time t2
{
  switch (pkt.type) {
    case TYPE_RADIOACTIVE_PELLET: {
      update_pellet(pkt, nts, t2);
      break;
    }

    case TYPE_GAMMA: {
      gammapkt::do_gamma(pkt, t2);

      if (pkt.type != TYPE_GAMMA && pkt.type != TYPE_ESCAPE) {
        safeadd(globals::timesteps[nts].gamma_dep, pkt.e_cmf);
      }
      break;
    }

    case TYPE_RPKT: {
      do_rpkt(pkt, t2);

      if (pkt.type == TYPE_ESCAPE) {
        safeadd(globals::timesteps[nts].cmf_lum, pkt.e_cmf);
      }
      break;
    }

    case TYPE_NONTHERMAL_PREDEPOSIT: {
      do_nonthermal_predeposit(pkt, nts, t2);
      break;
    }

    case TYPE_NTLEPTON: {
      nonthermal::do_ntlepton(pkt);
      break;
    }

    case TYPE_PRE_KPKT: {
      kpkt::do_kpkt_blackbody(pkt);
      break;
    }

    case TYPE_KPKT: {
      if (grid::modelgrid[grid::get_cell_modelgridindex(pkt.where)].thick == 1 || EXPANSIONOPACITIES_ON) {
        kpkt::do_kpkt_blackbody(pkt);
      } else {
        kpkt::do_kpkt(pkt, t2, nts);
      }
      break;
    }

    default: {
      printout("packet_prop: Unknown packet type %d. Abort.\n", pkt.type);
      assert_testmodeonly(false);
      std::abort();
    }
  }
}

static auto std_compare_packets_bymodelgriddensity(const Packet &p1, const Packet &p2) -> bool {
  // return true if packet p1 goes before p2

  // move escaped packets to the end of the list for better performance
  const bool esc1 = (p1.type == TYPE_ESCAPE);
  const bool esc2 = (p2.type == TYPE_ESCAPE);

  if (!esc1 && esc2) {
    return true;
  }
  if (esc1) {
    return false;
  }

  // const auto ts_end = globals::timesteps[globals::timestep].start + globals::timesteps[globals::timestep].width;

  // const bool pktdone1 = (p1.prop_time >= ts_end);
  // const bool pktdone2 = (p2.prop_time >= ts_end);

  // if (!pktdone1 && pktdone2) {
  //   return true;
  // }
  // if (pktdone1) {
  //   return false;
  // }

  // for both non-escaped packets, order by descending cell density
  const int mgi1 = grid::get_cell_modelgridindex(p1.where);
  const int mgi2 = grid::get_cell_modelgridindex(p2.where);
  const auto rho1 = mgi1 < grid::get_npts_model() ? grid::get_rho(mgi1) : 0.0;
  const auto rho2 = mgi2 < grid::get_npts_model() ? grid::get_rho(mgi2) : 0.0;

  if (rho1 > rho2) {
    return true;
  }

  if (rho1 == rho2 && (mgi1 < mgi2)) {
    return true;
  }

  // same cell, order by type
  if ((mgi1 == mgi2) && (p1.type < p2.type)) {
    return true;
  }

  // same cell and type, order by decreasing frequency
  if ((mgi1 == mgi2) && (p1.type == p2.type) && (p1.nu_cmf > p2.nu_cmf)) {
    return true;
  }

  return false;
}

static void do_cell_packet_updates(std::span<Packet> packets, const int nts, const double ts_end) {
  auto update_packet = [ts_end, nts](auto &pkt) {
    const int mgi = grid::get_cell_modelgridindex(pkt.where);
    int newmgi = mgi;
    while (pkt.prop_time < ts_end && pkt.type != TYPE_ESCAPE && (newmgi == mgi || newmgi == grid::get_npts_model())) {
      do_packet(pkt, ts_end, nts);
      newmgi = grid::get_cell_modelgridindex(pkt.where);
    }
  };

#if defined(STDPAR_ON) || !defined(_OPENMP)
  std::for_each(EXEC_PAR packets.begin(), packets.end(), update_packet);
#else
#pragma omp parallel for schedule(nonmonotonic : dynamic)
  for (ptrdiff_t i = 0; i < std::ssize(packets); i++) {
    update_packet(packets[i]);
  }
#endif
}

void update_packets(const int my_rank, const int nts, std::span<Packet> packets)
// Subroutine to move and update packets during the current timestep (nts)
{
  // At the start, the packets have all either just been initialised or have already been
  // processed for one or more timesteps. Those that are pellets will just be sitting in the
  // matter. Those that are photons (or one sort or another) will already have a position and
  // a direction.
  const double ts = globals::timesteps[nts].start;
  const double tw = globals::timesteps[nts].width;
  const double ts_end = ts + tw;

  const auto time_update_packets_start = std::time(nullptr);
  printout("timestep %d: start update_packets at time %ld\n", nts, time_update_packets_start);
  bool timestepcomplete = false;
  int passnumber = 0;
  while (!timestepcomplete) {
    const auto sys_time_start_pass = std::time(nullptr);

    // printout("sorting packets...");

    std::sort(EXEC_PAR_UNSEQ std::begin(packets), std::end(packets), std_compare_packets_bymodelgriddensity);

    // printout("took %lds\n", std::time(nullptr) - sys_time_start_pass);

    printout("  update_packets timestep %d pass %3d: started at %ld\n", nts, passnumber, sys_time_start_pass);

    const int count_pktupdates = static_cast<int>(std::ranges::count_if(
        packets, [ts_end](const auto &pkt) { return pkt.prop_time < ts_end && pkt.type != TYPE_ESCAPE; }));
    const int updatecellcounter_beforepass = stats::get_counter(stats::COUNTER_UPDATECELL);
    auto *packetgroupstart = packets.data();

    for (auto &pkt : packets) {
      if ((pkt.type != TYPE_ESCAPE && pkt.prop_time < ts_end)) {
        const int mgi = grid::get_cell_modelgridindex(pkt.where);
        const bool cellcache_change_cell_required =
            (mgi != grid::get_npts_model() && globals::cellcache[cellcacheslotid].cellnumber != mgi &&
             grid::modelgrid[mgi].thick != 1);

        if (cellcache_change_cell_required) {
          if (packetgroupstart != &pkt) {
            do_cell_packet_updates(std::span(packetgroupstart, std::distance(packetgroupstart, &pkt)), nts, ts_end);
          }

#ifdef _OPENMP
#pragma omp critical(cellchange)
#endif
          {
            stats::increment(stats::COUNTER_UPDATECELL);
            cellcache_change_cell(mgi);
          }
          packetgroupstart = &pkt;
        }
      }
    }
    const auto packets_remaining = std::distance(packetgroupstart, packets.data() + packets.size());
    if (packets_remaining > 0) {
      do_cell_packet_updates(std::span(packetgroupstart, packets_remaining), nts, ts_end);
    }

    timestepcomplete = std::ranges::all_of(
        packets, [ts_end](const auto &pkt) { return pkt.prop_time >= ts_end || pkt.type == TYPE_ESCAPE; });

    const int cellcacheresets = stats::get_counter(stats::COUNTER_UPDATECELL) - updatecellcounter_beforepass;
    printout(
        "  update_packets timestep %d pass %3d: finished at %ld packetsupdated %7d cellcacheresets %7d (took %lds)\n",
        nts, passnumber, std::time(nullptr), count_pktupdates, cellcacheresets,
        std::time(nullptr) - sys_time_start_pass);

    passnumber++;
  }

  stats::pkt_action_counters_printout(nts);

  const auto time_update_packets_end_thisrank = std::time(nullptr);
  printout("timestep %d: end of update_packets for this rank at time %ld\n", nts, time_update_packets_end_thisrank);

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);  // hold all processes once the packets are updated
#endif
  printout(
      "timestep %d: time after update packets for all processes %ld (rank %d took %lds, waited %lds, total %lds)\n",
      nts, std::time(nullptr), my_rank, time_update_packets_end_thisrank - time_update_packets_start,
      std::time(nullptr) - time_update_packets_end_thisrank, std::time(nullptr) - time_update_packets_start);
}
