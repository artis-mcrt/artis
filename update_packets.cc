#include "update_packets.h"

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <span>

#include "artisoptions.h"
#include "constants.h"
#include "decay.h"
#include "gammapkt.h"
#include "globals.h"
#include "grid.h"
#include "kpkt.h"
#include "nonthermal.h"
#include "packet.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "update_grid.h"
#include "vectors.h"

namespace {

void do_nonthermal_predeposit(Packet &pkt, const int nts, const double t2) {
  double en_deposited = pkt.e_cmf;
  const auto mgi = grid::get_cell_modelgridindex(pkt.where);
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
  const auto priortype = pkt.type;
  const double ts = pkt.prop_time;
  const auto deposit_type =
      (pkt.type == TYPE_NONTHERMAL_PREDEPOSIT_ALPHA) ? TYPE_NTALPHA_DEPOSITED : TYPE_NTLEPTON_DEPOSITED;

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
      pkt.nu_cmf = (particle_en - endot * (t_new - ts)) / H;
    }

    pkt.pos = vec_scale(pkt.pos, t_new / ts);
    pkt.prop_time = t_new;
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
void update_pellet(Packet &pkt, const int nts, const double t2) {
  assert_always(pkt.prop_time < t2);
  const double ts = pkt.prop_time;

  const double tdecay = pkt.tdecay;  // after packet_init(), this value never changes
  if (tdecay > t2) {
    // It won't decay in this timestep, so just need to move it on with the flow.
    pkt.pos = vec_scale(pkt.pos, t2 / ts);
    pkt.prop_time = t2;

    // That's all that needs to be done for the inactive pellet.
  } else if (tdecay > ts) {
    // The packet decays in the current timestep.
    atomicadd(globals::timesteps[nts].pellet_decays, 1);

    pkt.prop_time = tdecay;
    pkt.pos = vec_scale(pkt.pos, tdecay / ts);

    if (pkt.originated_from_particlenotgamma)  // will decay to non-thermal particle
    {
      if (pkt.pellet_decaytype == decay::DECAYTYPE_BETAPLUS) {
        pkt.type = TYPE_NONTHERMAL_PREDEPOSIT_BETAPLUS;
        atomicadd(globals::timesteps[nts].positron_emission, pkt.e_cmf);
      } else if (pkt.pellet_decaytype == decay::DECAYTYPE_BETAMINUS) {
        pkt.type = TYPE_NONTHERMAL_PREDEPOSIT_BETAMINUS;
        atomicadd(globals::timesteps[nts].electron_emission, pkt.e_cmf);
      } else if (pkt.pellet_decaytype == decay::DECAYTYPE_ALPHA) {
        atomicadd(globals::timesteps[nts].alpha_emission, pkt.e_cmf);
        pkt.type = TYPE_NONTHERMAL_PREDEPOSIT_ALPHA;
      } else if constexpr (TESTMODE) {
        printout("ERROR: pellet marked as particle emission is for decaytype %d != any of (alpha, beta+, beta-)\n",
                 pkt.pellet_decaytype);
        std::abort();
      } else {
        __builtin_unreachable();
      }
      pkt.em_time = pkt.prop_time;
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
    stats::increment(stats::COUNTER_K_STAT_FROM_EARLIERDECAY);

    pkt.prop_time = globals::tmin;
  } else if constexpr (TESTMODE) {
    printout("ERROR: Something wrong with decaying pellets. tdecay %g ts %g (ts + tw) %g\n", tdecay, ts, t2);
    assert_testmodeonly(false);
  } else {
    __builtin_unreachable();
  }
}

void do_packet(Packet &pkt, const double t2, const int nts)
// update a packet no further than time t2
{
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

      if (pkt.type == TYPE_ESCAPE) {
        atomicadd(globals::timesteps[nts].cmf_lum, pkt.e_cmf);
      }
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

    case TYPE_NTALPHA_DEPOSITED: {
      nonthermal::do_ntalpha_deposit(pkt);
      break;
    }

    case TYPE_PRE_KPKT: {
      kpkt::do_kpkt_blackbody(pkt);
      break;
    }

    case TYPE_KPKT: {
      const int mgi = grid::get_cell_modelgridindex(pkt.where);
      const int nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
      if (grid::modelgrid[nonemptymgi].thick == 1 ||
          (EXPANSIONOPACITIES_ON && RPKT_BOUNDBOUND_THERMALISATION_PROBABILITY >= 0.)) {
        kpkt::do_kpkt_blackbody(pkt);
      } else {
        kpkt::do_kpkt(pkt, t2, nts);
      }
      break;
    }

    default: {
      if constexpr (TESTMODE) {
        printout("ERROR: Unknown packet type %d\n", pkt.type);
        assert_testmodeonly(false);
      } else {
        __builtin_unreachable();
      }
    }
  }
}

auto std_compare_packets_bymodelgriddensity(const Packet &p1, const Packet &p2) -> bool {
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
  const auto rho1 = mgi1 < grid::get_npts_model() ? grid::get_rho(grid::get_nonemptymgi_of_mgi(mgi1)) : 0.0;
  const auto rho2 = mgi2 < grid::get_npts_model() ? grid::get_rho(grid::get_nonemptymgi_of_mgi(mgi2)) : 0.0;

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

void do_cell_packet_updates(std::span<Packet> packets, const int nts, const double ts_end) {
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
#ifdef GPU_ON
#pragma omp target teams distribute parallel for
#else
#pragma omp parallel for schedule(nonmonotonic : dynamic)
#endif
  for (ptrdiff_t i = 0; i < std::ssize(packets); i++) {
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
  printout("timestep %d: start update_packets at time %ld\n", nts, time_update_packets_start);
  bool timestepcomplete = false;
  int passnumber = 0;
  while (!timestepcomplete) {
    const auto sys_time_start_pass = std::time(nullptr);

    std::ranges::SORT_OR_STABLE_SORT(packets, std_compare_packets_bymodelgriddensity);

    printout("  update_packets timestep %d pass %3d: started at %ld\n", nts, passnumber, sys_time_start_pass);

    const int count_pktupdates = static_cast<int>(std::ranges::count_if(
        packets, [ts_end](const auto &pkt) { return pkt.prop_time < ts_end && pkt.type != TYPE_ESCAPE; }));
    const int updatecellcounter_beforepass = stats::get_counter(stats::COUNTER_UPDATECELL);
    auto *packetgroupstart = packets.data();

    for (auto &pkt : packets) {
      if ((pkt.type != TYPE_ESCAPE && pkt.prop_time < ts_end)) {
        const int mgi = grid::get_cell_modelgridindex(pkt.where);
        const int nonemptymgi = (mgi < grid::get_npts_model()) ? grid::get_nonemptymgi_of_mgi(mgi) : -1;
        const bool cellcache_change_cell_required =
            (nonemptymgi >= 0 && globals::cellcache[cellcacheslotid].nonemptymgi != nonemptymgi &&
             grid::modelgrid[nonemptymgi].thick != 1);

        if (cellcache_change_cell_required) {
          if (packetgroupstart != &pkt) {
            do_cell_packet_updates(std::span(packetgroupstart, &pkt - packetgroupstart), nts, ts_end);
          }

#ifdef _OPENMP
#pragma omp critical(cellchange)
#endif
          {
            stats::increment(stats::COUNTER_UPDATECELL);
            cellcache_change_cell(nonemptymgi);
          }
          packetgroupstart = &pkt;
        }
      }
    }
    const auto packets_remaining = packets.data() + packets.size() - packetgroupstart;
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
      nts, std::time(nullptr), globals::my_rank, time_update_packets_end_thisrank - time_update_packets_start,
      std::time(nullptr) - time_update_packets_end_thisrank, std::time(nullptr) - time_update_packets_start);
}
