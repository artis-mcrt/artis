#include "sn3d.h"
#include "gamma.h"
#include "grid_init.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nonthermal.h"
#include "update_packets.h"
#include "rpkt.h"
#include "vectors.h"

#include <algorithm>

static void update_pellet(
  PKT *pkt_ptr, const bool decay_to_kpkt, const bool decay_to_ntlepton, const int nts, const double t2)
{
  // Handle inactive pellets. Need to do two things (a) check if it
  // decays in this time step and if it does handle that. (b) if it doesn't decay in
  // this time step then just move the packet along with the matter for the
  // start of the next time step.
  assert(pkt_ptr->prop_time < t2);
  assert(!decay_to_kpkt || !decay_to_ntlepton); // can't decay to both!
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
    #ifdef _OPENMP
      #pragma omp atomic
    #endif
    time_step[nts].pellet_decays++;

    pkt_ptr->prop_time = tdecay;
    vec_scale(pkt_ptr->pos, tdecay / ts);

    if (decay_to_kpkt)
    {
      pkt_ptr->type = TYPE_KPKT;
      pkt_ptr->absorptiontype = -6;
    }
    else if (decay_to_ntlepton)
    {
      #ifdef _OPENMP
        #pragma omp atomic
      #endif
      time_step[nts].positron_dep += pkt_ptr->e_cmf;

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
    // The position is already set at tmin so don't need to move it. Assume
    // that it is fixed in place from decay to tmin - i.e. short mfp.

    pkt_ptr->e_cmf *= tdecay / tmin;
    //pkt_ptr->type = TYPE_KPKT;
    pkt_ptr->type = TYPE_PRE_KPKT;
    pkt_ptr->absorptiontype = -7;
    //if (tid == 0) k_stat_from_earlierdecay++;
    k_stat_from_earlierdecay++;

    //printout("already decayed packets and propagation by packet_prop\n");
    pkt_ptr->prop_time = tmin;
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
        #ifdef _OPENMP
          #pragma omp atomic
        #endif
        time_step[nts].gamma_dep += pkt_ptr->e_cmf;
      }
      break;

    case TYPE_RPKT:
      while (do_rpkt(pkt_ptr, t2))
      {
        ;
      }

      if (pkt_ptr->type == TYPE_ESCAPE)
      {
        #ifdef _OPENMP
          #pragma omp atomic
        #endif
        time_step[nts].cmf_lum += pkt_ptr->e_cmf;
      }
      break;

    case TYPE_NTLEPTON:
      do_ntlepton(pkt_ptr);
      break;

    case TYPE_KPKT:
    case TYPE_PRE_KPKT:
    case TYPE_GAMMA_KPKT:
      /*It's a k-packet - convert to r-packet (low freq).*/
      //printout("k-packet propagation\n");

      //t_change_type = do_kpkt(pkt_ptr, t_current, t2);
      if (pkt_type == TYPE_PRE_KPKT || modelgrid[cell[pkt_ptr->where].modelgridindex].thick == 1)
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


static int compare_packets_bymodelgriddensity(const void *p1, const void *p2)
{
  // <0 The element pointed by p1 goes before the element pointed by p2
  // 0  The element pointed by p1 is equivalent to the element pointed by p2
  // >0 The element pointed by p1 goes after the element pointed by p2

  // move escaped packets to the end of the list for better performance
  const bool esc1 = (((PKT *) p1)->type == TYPE_ESCAPE);
  const bool esc2 = (((PKT *) p2)->type == TYPE_ESCAPE);
  if (esc1 && !esc2)
    return 1;
  else if (!esc1 && esc2)
    return -1;
  else if (esc1 && esc2)
    return 0;

  // for both non-escaped packets, order by descending cell density
  const int a1_where = ((PKT *) p1)->where;
  const int a2_where = ((PKT *) p2)->where;

  const int mgi1 = cell[a1_where].modelgridindex;
  const int mgi2 = cell[a2_where].modelgridindex;
  const double rho_diff = get_rho(mgi1) - get_rho(mgi2);
  if (rho_diff < 0)
    return 1;
  else if (rho_diff > 0)
    return -1;
  else
    return (mgi1 - mgi2);
}

static bool std_compare_packets_bymodelgriddensity(const PKT &p1, const PKT &p2)
{
  // true if p1 goes before the element pointed by p2

  // move escaped packets to the end of the list for better performance
  const bool esc1 = (p1.type == TYPE_ESCAPE);
  const bool esc2 = (p2.type == TYPE_ESCAPE);
  if (!esc1 && esc2)
    return true;
  else if (esc1 && !esc2)
    return false;

  // for both non-escaped packets, order by descending cell density
  const int a1_where = p1.where;
  const int a2_where = p2.where;

  const int mgi1 = cell[a1_where].modelgridindex;
  const int mgi2 = cell[a2_where].modelgridindex;
  if (get_rho(mgi1) > get_rho(mgi2))
    return true;
  else if (get_rho(mgi1) < get_rho(mgi2))
    return false;
  else
    return (mgi1 < mgi2);
}


void update_packets(const int nts, PKT *pkt)
// Subroutine to move and update packets during the current timestep (nts)
{
  /** At the start, the packets have all either just been initialised or have already been
  processed for one or more timesteps. Those that are pellets will just be sitting in the
  matter. Those that are photons (or one sort or another) will already have a position and
  a direction.*/

  const double ts = time_step[nts].start;
  const double tw = time_step[nts].width;

  printout("start of parallel update_packets loop %ld\n", time(NULL));
  /// Initialise the OpenMP reduction target to zero
  bool timestepcomplete = false;
  int passnumber = 0;
  #ifdef _OPENMP
    #pragma omp parallel
    //copyin(debuglevel,nuJ,J)
  #endif
  while (!timestepcomplete)
  {
    timestepcomplete = true;

    const bool photonpkt_pass = (passnumber % 2 == 1);

    // after a !photonpkt_pass, packets have not propagated and changed cells
    if (passnumber == 0 || !photonpkt_pass)
    {
      const time_t sys_time_start_sort = time(NULL);

      qsort(pkt, npkts, sizeof(PKT), compare_packets_bymodelgriddensity);
      // std::sort(pkt, pkt + npkts, std_compare_packets_bymodelgriddensity);

      const int duration_sortpackets = time(NULL) - sys_time_start_sort;
      if (duration_sortpackets > 1)
      {
        printout("sorting packets took %ds\n", duration_sortpackets);
      }
    }

    int count_photpktupdates = 0;
    int count_otherupdates = 0;
    const int updatecellcounter_beforepass = updatecellcounter;

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic) reduction(+:escounter,resonancescatterings,cellcrossings,nesc,updatecellcounter,coolingratecalccounter,upscatter,downscatter,ma_stat_activation_collexc,ma_stat_activation_collion,ma_stat_activation_ntcollexc,ma_stat_activation_ntcollion,ma_stat_activation_bb,ma_stat_activation_bf,ma_stat_activation_fb,ma_stat_deactivation_colldeexc,ma_stat_deactivation_collrecomb,ma_stat_deactivation_bb,ma_stat_deactivation_fb,k_stat_to_ma_collexc,k_stat_to_ma_collion,k_stat_to_r_ff,k_stat_to_r_fb,k_stat_from_ff,k_stat_from_bf,nt_stat_from_gamma,k_stat_from_earlierdecay)
    #endif
    for (int n = 0; n < npkts; n++)
    {
      PKT *pkt_ptr = &pkt[n];

      // if (pkt_ptr->type == TYPE_ESCAPE)
      // {
      //   printout("packet index %d already escaped. Skipping rest of packets (which are all escaped).\n", n);
      //   // for (int n2 = n; n2 < npkts; n2++)
      //   // {
      //   //   assert(pkt[n2].type == TYPE_ESCAPE);
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
        const int mgi = cell[cellindex].modelgridindex;
        /// for non empty cells update the global available level populations and cooling terms
        /// Reset cellhistory if packet starts up in another than the last active cell
        if (mgi != MMODELGRID && cellhistory[tid].cellnumber != mgi)
        {
          updatecellcounter++;
          cellhistory_reset(mgi, false);
        }

        timestepcomplete = false;
        if (!photonpkt_pass && pkt_ptr->type != TYPE_RPKT && pkt_ptr->type != TYPE_GAMMA)
        {
          bool workedonpacket = false;
          while (pkt_ptr->type != TYPE_RPKT && pkt_ptr->type != TYPE_GAMMA && pkt_ptr->prop_time < (ts + tw) && pkt_ptr->type != TYPE_ESCAPE)
          {
            workedonpacket = true;
            do_packet(pkt_ptr, ts + tw, nts);
          }

          count_otherupdates += workedonpacket ? 1 : 0;
        }
        else if (photonpkt_pass && (pkt_ptr->type == TYPE_RPKT || pkt_ptr->type == TYPE_GAMMA))
        {
          count_photpktupdates++;
          if (pkt_ptr->type == TYPE_RPKT && modelgrid[mgi].thick != 1 && mgi != MMODELGRID)
          {
            // printout("calculate_kappa_rpkt_cont(mgi %d)...", mgi);
            calculate_kappa_rpkt_cont(pkt_ptr, mgi);
            // printout("done\n");
          }
          enum packet_type oldtype = pkt_ptr->type;
          int newmgi = mgi;
          while ((newmgi == mgi || newmgi == MMODELGRID) && pkt_ptr->prop_time < (ts + tw) && pkt_ptr->type == oldtype)
          {
            do_packet(pkt_ptr, ts + tw, nts);
            const int newcellnum = pkt_ptr->where;
            newmgi = cell[newcellnum].modelgridindex;
          }
        }
      }
    }
    printout("[debug] update_packets pass %d: updated %d photon packets and %d other packets cellhistoryresets %d at %ld\n",
             passnumber, count_photpktupdates, count_otherupdates, updatecellcounter - updatecellcounter_beforepass, time(NULL));

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

  int mgi_diff = cell[a1->where].modelgridindex - cell[a2->where].modelgridindex;
  if (mgi_diff < 0)
    return -1;
  else if (mgi_diff > 0)
    return 1;
  else
    return 0;
}*/
