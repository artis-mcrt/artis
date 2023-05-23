#include "update_packets.h"

#include <algorithm>

#include "decay.h"
#include "gammapkt.h"
#include "grid.h"
#include "kpkt.h"
#include "nonthermal.h"
#include "packet.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "update_grid.h"

static void do_nonthermal_predeposit(struct packet *pkt_ptr, const int nts, const double t2) {
  if constexpr (!INSTANT_PARTICLE_DEPOSITION) {
    const double ts = pkt_ptr->prop_time;
    const int mgi = grid::get_cell_modelgridindex(pkt_ptr->where);
    const double rho = grid::get_rho(mgi);

    // endot is energy loss rate (positive) in [erg/s]
    // endot [erg/s] from Barnes et al. (2016). see their figure 6.
    const double endot = (pkt_ptr->pellet_decaytype == decay::DECAYTYPE_ALPHA) ? 5.e11 * MEV * rho : 4.e10 * MEV * rho;

    // A discrete absorption event should occur somewhere along the
    // continuous track from initial kinetic energy to zero KE.
    // The probability of being absorbed in energy range [E, E+delta_E] is proportional to
    // endot(E) * delta_t = endot(E) * delta_E / endot(E) = delta_E (delta_t is the time spent in the bin range)
    // so all final energies are equally likely.
    // Choose random en_absorb [0, particle_en]
    const double zrand = rng_uniform();

    const double particle_en = H * pkt_ptr->nu_cmf;
    const double en_absorb = zrand * particle_en;

    // for endot independent of energy, the next line is trival (for E dependent endot, an integral would be needed)
    const double t_absorb = ts + en_absorb / endot;

    // const double deltat_zeroen = particle_en / endot;
    // const double t_sim_zeroen = ts + deltat_zeroen;
    // printout("%s packet: nts %d energy_mev %g ts %g deltat_zeroen %g t_sim_zeroen %g t2 %g t_absorb %g\n",
    //          pkt_ptr->pellet_decaytype == decay::DECAYTYPE_ALPHA ? "alpha" : "beta", nts,
    //          particle_en / MEV,
    //          ts / 86400,
    //          deltat_zeroen / 86400, t_sim_zeroen / 86400, t2 / 86400, t_absorb / 86400);

    if (t_absorb > t2) {
      // absorption happens beyond the end of the current timestep,
      // so reduce the particle energy for the end of this timestep
      pkt_ptr->nu_cmf = (particle_en - endot * (t2 - ts)) / H;
      vec_scale(pkt_ptr->pos, t2 / ts);
      pkt_ptr->prop_time = t2;
      return;
    }

    // absorption happen part way through this timestep
    vec_scale(pkt_ptr->pos, t_absorb / ts);
    pkt_ptr->prop_time = t_absorb;
  }

  // absorption happens

  if (pkt_ptr->pellet_decaytype == decay::DECAYTYPE_ALPHA) {
    safeadd(globals::time_step[nts].alpha_dep, pkt_ptr->e_cmf);
  } else if (pkt_ptr->pellet_decaytype == decay::DECAYTYPE_BETAMINUS) {
    safeadd(globals::time_step[nts].electron_dep, pkt_ptr->e_cmf);
  } else if (pkt_ptr->pellet_decaytype == decay::DECAYTYPE_BETAPLUS) {
    safeadd(globals::time_step[nts].positron_dep, pkt_ptr->e_cmf);
  }

  pkt_ptr->type = TYPE_NTLEPTON;
}

static void update_pellet(struct packet *pkt_ptr, const int nts, const double t2) {
  // Handle inactive pellets. Need to do two things (a) check if it
  // decays in this time step and if it does handle that. (b) if it doesn't decay in
  // this time step then just move the packet along with the matter for the
  // start of the next time step.
  assert_always(pkt_ptr->prop_time < t2);
  const double ts = pkt_ptr->prop_time;

  const double tdecay = pkt_ptr->tdecay;  // after packet_init(), this value never changes
  if (tdecay > t2) {
    // It won't decay in this timestep, so just need to move it on with the flow.
    vec_scale(pkt_ptr->pos, t2 / ts);
    pkt_ptr->prop_time = t2;

    // That's all that needs to be done for the inactive pellet.
  } else if (tdecay > ts) {
    // The packet decays in the current timestep.
    safeincrement(globals::time_step[nts].pellet_decays);

    pkt_ptr->prop_time = tdecay;
    vec_scale(pkt_ptr->pos, tdecay / ts);

    if (pkt_ptr->originated_from_particlenotgamma)  // will decay to non-thermal particle
    {
      if (pkt_ptr->pellet_decaytype == decay::DECAYTYPE_BETAPLUS) {
        safeadd(globals::time_step[nts].positron_dep, pkt_ptr->e_cmf);
        pkt_ptr->type = TYPE_NTLEPTON;
        pkt_ptr->absorptiontype = -10;
      } else if (pkt_ptr->pellet_decaytype == decay::DECAYTYPE_BETAMINUS) {
        safeadd(globals::time_step[nts].electron_emission, pkt_ptr->e_cmf);
        pkt_ptr->em_time = pkt_ptr->prop_time;
        pkt_ptr->type = TYPE_NONTHERMAL_PREDEPOSIT;
        pkt_ptr->absorptiontype = -10;
      } else if (pkt_ptr->pellet_decaytype == decay::DECAYTYPE_ALPHA) {
        safeadd(globals::time_step[nts].alpha_emission, pkt_ptr->e_cmf);
        pkt_ptr->em_time = pkt_ptr->prop_time;
        pkt_ptr->type = TYPE_NONTHERMAL_PREDEPOSIT;
        pkt_ptr->absorptiontype = -10;
      }
    } else {
      safeadd(globals::time_step[nts].gamma_emission, pkt_ptr->e_cmf);
      // decay to gamma-ray, kpkt, or ntlepton
      gammapkt::pellet_gamma_decay(pkt_ptr);
    }
  } else if ((tdecay > 0) && (nts == 0)) {
    // These are pellets whose decay times were before the first time step
    // They will be made into r-packets with energy reduced for doing work on the
    // ejecta following Lucy 2004.
    // The position is already set at globals::tmin so don't need to move it. Assume
    // that it is fixed in place from decay to globals::tmin - i.e. short mfp.

    pkt_ptr->e_cmf *= tdecay / globals::tmin;
    pkt_ptr->type = TYPE_PRE_KPKT;
    pkt_ptr->absorptiontype = -7;
    stats::increment(stats::COUNTER_K_STAT_FROM_EARLIERDECAY);

    // printout("already decayed packets and propagation by packet_prop\n");
    pkt_ptr->prop_time = globals::tmin;
  } else {
    printout("ERROR: Something gone wrong with decaying pellets. tdecay %g ts %g (ts + tw) %g\n", tdecay, ts, t2);
    abort();
  }
}

static void do_packet(struct packet *const pkt_ptr, const double t2, const int nts)
// update a packet no further than time t2
{
  const int pkt_type = pkt_ptr->type;

  switch (pkt_type) {
    case TYPE_RADIOACTIVE_PELLET: {
      update_pellet(pkt_ptr, nts, t2);
      break;
    }

    case TYPE_GAMMA: {
      gammapkt::do_gamma(pkt_ptr, t2);

      if (pkt_ptr->type != TYPE_GAMMA && pkt_ptr->type != TYPE_ESCAPE) {
        safeadd(globals::time_step[nts].gamma_dep, pkt_ptr->e_cmf);
      }
      break;
    }

    case TYPE_RPKT: {
      do_rpkt(pkt_ptr, t2);

      if (pkt_ptr->type == TYPE_ESCAPE) {
        safeadd(globals::time_step[nts].cmf_lum, pkt_ptr->e_cmf);
      }
      break;
    }

    case TYPE_NONTHERMAL_PREDEPOSIT: {
      do_nonthermal_predeposit(pkt_ptr, nts, t2);
      break;
    }

    case TYPE_NTLEPTON: {
      nonthermal::do_ntlepton(pkt_ptr);
      break;
    }

    case TYPE_PRE_KPKT: {
      kpkt::do_kpkt_bb(pkt_ptr);
      break;
    }

    case TYPE_KPKT: {
      if (grid::modelgrid[grid::get_cell_modelgridindex(pkt_ptr->where)].thick == 1) {
        kpkt::do_kpkt_bb(pkt_ptr);
      } else {
        kpkt::do_kpkt(pkt_ptr, t2, nts);
      }
      break;
    }

    case TYPE_MA: {
      do_macroatom(pkt_ptr, nts);
      break;
    }

    default:
      printout("packet_prop: Unknown packet type %d. Abort.\n", pkt_ptr->type);
      abort();
  }
}

static auto std_compare_packets_bymodelgriddensity(const struct packet &p1, const struct packet &p2) -> bool {
  // return true if packet p1 goes before p2

  // move escaped packets to the end of the list for better performance
  const bool esc1 = (p1.type == TYPE_ESCAPE);
  const bool esc2 = (p2.type == TYPE_ESCAPE);

  if (!esc1 && esc2) {
    return true;
  }
  if (esc1 && !esc2) {
    return false;
  }
  if (esc1 && esc2) {
    return false;
  }

  // for both non-escaped packets, order by descending cell density
  const int mgi1 = grid::get_cell_modelgridindex(p1.where);
  const int mgi2 = grid::get_cell_modelgridindex(p2.where);
  if (grid::get_rho(mgi1) > grid::get_rho(mgi2)) {
    return true;
  }

  if (grid::get_rho(mgi1) == grid::get_rho(mgi2) && (mgi1 < mgi2)) {
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

void update_packets(const int my_rank, const int nts, struct packet *packets)
// Subroutine to move and update packets during the current timestep (nts)
{
  // At the start, the packets have all either just been initialised or have already been
  // processed for one or more timesteps. Those that are pellets will just be sitting in the
  // matter. Those that are photons (or one sort or another) will already have a position and
  // a direction.

  const double ts = globals::time_step[nts].start;
  const double tw = globals::time_step[nts].width;

  const time_t time_update_packets_start = time(nullptr);
  printout("timestep %d: start update_packets at time %ld\n", nts, time_update_packets_start);
  bool timestepcomplete = false;
  int passnumber = 0;
  while (!timestepcomplete) {
    timestepcomplete = true;  // will be set false if any packets did not finish propagating in this pass

    const time_t sys_time_start_pass = time(nullptr);

    // printout("sorting packets...");

    std::sort(packets, packets + globals::npkts, std_compare_packets_bymodelgriddensity);

    // printout("took %lds\n", time(nullptr) - sys_time_start_pass);

    printout("  update_packets timestep %d pass %3d: started at %ld\n", nts, passnumber, sys_time_start_pass);

    int count_pktupdates = 0;
    const int updatecellcounter_beforepass = stats::get_counter(stats::COUNTER_UPDATECELL);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int n = 0; n < globals::npkts; n++) {
      struct packet *pkt_ptr = &packets[n];

      // if (pkt_ptr->type == TYPE_ESCAPE)
      // {
      //   printout("packet index %d already escaped. Skipping rest of packets (which are all escaped).\n", n);
      //   // for (int n2 = n; n2 < globals::npkts; n2++)
      //   // {
      //   //   assert_always(packets[n2].type == TYPE_ESCAPE);
      //   // }
      //   break;
      // }
      // pkt_ptr->timestep = nts;

      if (passnumber == 0) {
        pkt_ptr->interactions = 0;
      }

      if (pkt_ptr->type != TYPE_ESCAPE && pkt_ptr->prop_time < (ts + tw)) {
        const int cellindex = pkt_ptr->where;
        const int mgi = grid::get_cell_modelgridindex(cellindex);
        /// for non empty cells update the global available level populations and cooling terms
        /// Reset cellhistory if packet starts up in another than the last active cell
        if (mgi != grid::get_npts_model() && globals::cellhistory[tid].cellnumber != mgi) {
          stats::increment(stats::COUNTER_UPDATECELL);
          cellhistory_reset(mgi, false);
        }

        // enum packet_type oldtype = pkt_ptr->type;
        int newmgi = mgi;
        bool workedonpacket = false;
        while ((newmgi == mgi || newmgi == grid::get_npts_model()) && pkt_ptr->prop_time < (ts + tw) &&
               pkt_ptr->type != TYPE_ESCAPE) {
          workedonpacket = true;
          do_packet(pkt_ptr, ts + tw, nts);
          const int newcellnum = pkt_ptr->where;
          newmgi = grid::get_cell_modelgridindex(newcellnum);
        }
        count_pktupdates += workedonpacket ? 1 : 0;

        if (pkt_ptr->type != TYPE_ESCAPE && pkt_ptr->prop_time < (ts + tw)) {
          timestepcomplete = false;
        }
      }
    }
    const int cellhistresets = stats::get_counter(stats::COUNTER_UPDATECELL) - updatecellcounter_beforepass;
    printout(
        "  update_packets timestep %d pass %3d: finished at %ld packetsupdated %7d cellhistoryresets %7d (took %lds)\n",
        nts, passnumber, time(nullptr), count_pktupdates, cellhistresets, time(nullptr) - sys_time_start_pass);

    passnumber++;
  }

  stats::pkt_action_counters_printout(packets, nts);

  const time_t time_update_packets_end_thisrank = time(nullptr);
  printout("timestep %d: end of update_packets for this rank at time %ld\n", nts, time_update_packets_end_thisrank);

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);  // hold all processes once the packets are updated
#endif
  printout(
      "timestep %d: time after update packets for all processes %ld (rank %d took %lds, waited %lds, total %lds)\n",
      nts, time(nullptr), my_rank, time_update_packets_end_thisrank - time_update_packets_start,
      time(nullptr) - time_update_packets_end_thisrank, time(nullptr) - time_update_packets_start);
}
