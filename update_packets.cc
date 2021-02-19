#include "sn3d.h"
#include "gamma.h"
#include "grid.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "update_packets.h"
#include "rpkt.h"
#include "stats.h"
#include "vectors.h"

#include <algorithm>

static void update_pellet(
  PKT *pkt_ptr, const bool decay_to_kpkt, const bool decay_to_ntlepton, const int nts, const double t2)
{
  // Handle inactive pellets. Need to do two things (a) check if it
  // decays in this time step and if it does handle that. (b) if it doesn't decay in
  // this time step then just move the packet along with the matter for the
  // start of the next time step.
  assert_always(pkt_ptr->prop_time < t2);
  assert_always(!decay_to_kpkt || !decay_to_ntlepton); // can't decay to both!
  const double ts = pkt_ptr->prop_time;

  const double tdecay = pkt_ptr->tdecay; // after packet_init(), this value never changes
  if (tdecay > t2)
  {
    // It won't decay in this timestep, so just need to move it on with the flow.
    vec_scale(pkt_ptr->pos, t2 / ts);
    pkt_ptr->prop_time = t2;

    // That's all that needs to be done for the inactive pellet.
  }
  else if (tdecay > ts)
  {
    // The packet decays in the current timestep.
    safeincrement(globals::time_step[nts].pellet_decays);

    pkt_ptr->prop_time = tdecay;
    vec_scale(pkt_ptr->pos, tdecay / ts);

    if (decay_to_kpkt)
    {
      pkt_ptr->type = TYPE_KPKT;
      pkt_ptr->absorptiontype = -6;
    }
    else if (decay_to_ntlepton)
    {
      safeadd(globals::time_step[nts].positron_dep, pkt_ptr->e_cmf);

      pkt_ptr->type = TYPE_NTLEPTON;
      pkt_ptr->absorptiontype = -10;
    }
    else
    {
      // decay to gamma ray
      pellet_decay(nts, pkt_ptr);
      //printout("pellet to photon packet and propagation by packet_prop\n");
    }
  }
  else if ((tdecay > 0) && (nts == 0))
  {
    // These are pellets whose decay times were before the first time step
    // They will be made into r-packets with energy reduced for doing work on the
    // ejecta following Lucy 2004.
    // The position is already set at globals::tmin so don't need to move it. Assume
    // that it is fixed in place from decay to globals::tmin - i.e. short mfp.

    pkt_ptr->e_cmf *= tdecay / globals::tmin;
    //pkt_ptr->type = TYPE_KPKT;
    pkt_ptr->type = TYPE_PRE_KPKT;
    pkt_ptr->absorptiontype = -7;
    stats::increment(stats::COUNTER_K_STAT_FROM_EARLIERDECAY);

    //printout("already decayed packets and propagation by packet_prop\n");
    pkt_ptr->prop_time = globals::tmin;
  }
  else
  {
    printout("ERROR: Something gone wrong with decaying pellets. tdecay %g ts %g (ts + tw) %g\n", tdecay, ts, t2);
    abort();
  }
}


static void do_packet(PKT *const pkt_ptr, const double t2, const int nts)
// update a packet no further than time t2
{
  const int pkt_type = pkt_ptr->type; // avoid dereferencing multiple times

  switch (pkt_type)
  {
    case TYPE_56NI_PELLET:
    case TYPE_56CO_PELLET:
    case TYPE_57NI_PELLET:
    case TYPE_57CO_PELLET:
    case TYPE_48CR_PELLET:
    case TYPE_48V_PELLET:
      // decay to gamma ray
      update_pellet(pkt_ptr, false, false, nts, t2);
      break;

    case TYPE_52FE_PELLET:
    case TYPE_52MN_PELLET:
      // convert to kpkts
      update_pellet(pkt_ptr, true, false, nts, t2);
      break;

    case TYPE_57NI_POSITRON_PELLET:
    case TYPE_56CO_POSITRON_PELLET:
      // convert to to non-thermal leptons
      update_pellet(pkt_ptr, false, true, nts, t2);
      break;

    case TYPE_GAMMA:
      do_gamma(pkt_ptr, t2);
	    /* This returns a flag if the packet gets to t2 without
      changing to something else. If the packet does change it
      returns the time of change and sets everything for the
      new packet.*/
      if (pkt_ptr->type != TYPE_GAMMA && pkt_ptr->type != TYPE_ESCAPE)
      {
        safeadd(globals::time_step[nts].gamma_dep, pkt_ptr->e_cmf);
      }
      break;

    case TYPE_RPKT:
      do_rpkt(pkt_ptr, t2);

      if (pkt_ptr->type == TYPE_ESCAPE)
      {
        safeadd(globals::time_step[nts].cmf_lum, pkt_ptr->e_cmf);
      }
      break;

    case TYPE_NTLEPTON:
      nonthermal::do_ntlepton(pkt_ptr);
      break;

    case TYPE_KPKT:
    case TYPE_PRE_KPKT:
    case TYPE_GAMMA_KPKT:
      /*It's a k-packet - convert to r-packet (low freq).*/
      //printout("k-packet propagation\n");

      //t_change_type = do_kpkt(pkt_ptr, t_current, t2);
      if (pkt_type == TYPE_PRE_KPKT || globals::modelgrid[get_cell_modelgridindex(pkt_ptr->where)].thick == 1)
      {
        do_kpkt_bb(pkt_ptr);
      }
      else if (pkt_type == TYPE_KPKT)
      {
        do_kpkt(pkt_ptr, t2, nts);
      }
      else
      {
        printout("kpkt not of type TYPE_KPKT or TYPE_PRE_KPKT\n");
        abort();
        //t_change_type = do_kpkt_ffonly(pkt_ptr, t_current, t2);
      }
      break;

    case TYPE_MA:
      do_macroatom(pkt_ptr, nts);
      break;

    default:
      printout("packet_prop: Unknown packet type %d. Abort.\n", pkt_ptr->type);
      abort();
  }
}


static bool std_compare_packets_bymodelgriddensity(const PKT &p1, const PKT &p2)
{
  // return true if packet p1 goes before p2

  // move escaped packets to the end of the list for better performance
  const bool esc1 = (p1.type == TYPE_ESCAPE);
  const bool esc2 = (p2.type == TYPE_ESCAPE);

  if (!esc1 && esc2)
    return true;
  else if (esc1 && !esc2)
    return false;
  else if (esc1 && esc2)
    return false;

  // for both non-escaped packets, order by descending cell density
  const int a1_where = p1.where;
  const int a2_where = p2.where;

  const int mgi1 = get_cell_modelgridindex(a1_where);
  const int mgi2 = get_cell_modelgridindex(a2_where);
  if (get_rho(mgi1) > get_rho(mgi2))
    return true;

  if (get_rho(mgi1) == get_rho(mgi2) && (mgi1 < mgi2))
    return true;

  return false;
  // return (p1.type > p2.type);
}


void update_packets(const int my_rank, const int nts, PKT *pkt)
// Subroutine to move and update packets during the current timestep (nts)
{
  /** At the start, the packets have all either just been initialised or have already been
  processed for one or more timesteps. Those that are pellets will just be sitting in the
  matter. Those that are photons (or one sort or another) will already have a position and
  a direction.*/

  const double ts = globals::time_step[nts].start;
  const double tw = globals::time_step[nts].width;

  printout("start of parallel update_packets loop %ld\n", time(NULL));
  bool timestepcomplete = false;
  int passnumber = 0;
  while (!timestepcomplete)
  {
    timestepcomplete = true;

    const time_t sys_time_start_sort = time(NULL);

    std::sort(pkt, pkt + globals::npkts, std_compare_packets_bymodelgriddensity);

    const int duration_sortpackets = time(NULL) - sys_time_start_sort;
    if (duration_sortpackets > 1)
    {
      printout("sorting packets took %ds\n", duration_sortpackets);
    }

    int count_pktupdates = 0;
    const int updatecellcounter_beforepass = stats::get_counter(stats::COUNTER_UPDATECELL);

    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
    #endif
    for (int n = 0; n < globals::npkts; n++)
    {
      PKT *pkt_ptr = &pkt[n];

      // if (pkt_ptr->type == TYPE_ESCAPE)
      // {
      //   printout("packet index %d already escaped. Skipping rest of packets (which are all escaped).\n", n);
      //   // for (int n2 = n; n2 < globals::npkts; n2++)
      //   // {
      //   //   assert_always(pkt[n2].type == TYPE_ESCAPE);
      //   // }
      //   break;
      // }
      //pkt_ptr->timestep = nts;

      if (passnumber == 0)
      {
        pkt_ptr->interactions = 0;
        pkt_ptr->scat_count = 0;
      }

      if (pkt_ptr->type != TYPE_ESCAPE && pkt_ptr->prop_time < (ts + tw))
      {
        const int cellindex = pkt_ptr->where;
        const int mgi = get_cell_modelgridindex(cellindex);
        /// for non empty cells update the global available level populations and cooling terms
        /// Reset cellhistory if packet starts up in another than the last active cell
        if (mgi != MMODELGRID && globals::cellhistory[tid].cellnumber != mgi)
        {
          stats::increment(stats::COUNTER_UPDATECELL);
          cellhistory_reset(mgi, false);
        }

        // enum packet_type oldtype = pkt_ptr->type;
        int newmgi = mgi;
        bool workedonpacket = false;
        while ((newmgi == mgi || newmgi == MMODELGRID) && pkt_ptr->prop_time < (ts + tw) && pkt_ptr->type != TYPE_ESCAPE)
        {
          workedonpacket = true;
          do_packet(pkt_ptr, ts + tw, nts);
          const int newcellnum = pkt_ptr->where;
          newmgi = get_cell_modelgridindex(newcellnum);
        }
        count_pktupdates += workedonpacket ? 1 : 0;

        if (pkt_ptr->type != TYPE_ESCAPE && pkt_ptr->prop_time < (ts + tw))
        {
          timestepcomplete = false;
        }
      }
    }
    const int cellhistresets = stats::get_counter(stats::COUNTER_UPDATECELL) - updatecellcounter_beforepass;
    printout("  update_packets timestep %d pass %3d: updated packets %7d cellhistoryresets %7d at %ld\n",
             nts, passnumber, count_pktupdates, cellhistresets, time(NULL));

    passnumber++;
  }

  printout("end of update_packets parallel for loop %ld\n", time(NULL));
}


/*static int compare_packets_byposition(const void *p1, const void *p2)
/// Helper function to sort the phixslist by ascending threshold frequency.
{
  const PKT *a1 = (PKT *)(p1);
  const PKT *a2 = (PKT *)(p2);

  int cell_diff = a1->where - a2->where;
  if (cell_diff < 0)
    return -1;
  else if (cell_diff > 0)
    return 1;
  else
    return 0;
}


static int compare_packets_bymodelgridposition(const void *p1, const void *p2)
/// Helper function to sort the phixslist by ascending threshold frequency.
{
  const PKT *a1 = (PKT *)(p1);
  const PKT *a2 = (PKT *)(p2);

  int mgi_diff = get_cell_modelgridindex(a1->where) - get_cell_modelgridindex(a2->where);
  if (mgi_diff < 0)
    return -1;
  else if (mgi_diff > 0)
    return 1;
  else
    return 0;
}*/
