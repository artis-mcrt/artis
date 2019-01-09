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


extern void update_cell(const int mgi);


static void packet_prop(PKT *restrict const pkt_ptr, const double t1, const double t2, const int nts)
// Master routine for moving packets around. When it called,
//   it is given the time at start of inverval and at end - when it finishes,
//   everything the packet does during this time should be sorted out.
{
  double t_current = t1;

  /* 0 the scatter counter for the packet. */
  pkt_ptr->scat_count = 0;

  // t_current == PACKET_SAME which is < 0 after t2 has been reached;
  while (t_current >= 0)
  {
    /* Start by sorting out what sort of packet it is.*/
    //printout("start of packet_prop loop %d\n", pkt_ptr->type );
    const int pkt_type = pkt_ptr->type; // avoid dereferencing multiple times

    switch (pkt_type)
    {
      case TYPE_GAMMA:
        /*It's a gamma-ray packet.*/
        /* Call do_gamma. */
        //printout("gamma propagation\n");
        t_current = do_gamma(pkt_ptr, t_current, t2);
  	    /* This returns a flag if the packet gets to t2 without
        changing to something else. If the packet does change it
        returns the time of change and sets everything for the
        new packet.*/
        if (t_current >= 0)
        {
          #ifdef _OPENMP
            #pragma omp atomic
          #endif
          time_step[nts].gamma_dep += pkt_ptr->e_cmf; ///***??
        }
        break;

      case TYPE_RPKT:
        /*It's an r-packet. */
        //printout("r-pkt propagation\n");
        t_current = do_rpkt(pkt_ptr, t_current, t2);
  //       if (modelgrid[cell[pkt_ptr->where].modelgridindex].thick == 1)
  //         t_change_type = do_rpkt_thickcell( pkt_ptr, t_current, t2);
  //       else
  //         t_change_type = do_rpkt( pkt_ptr, t_current, t2);

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
          t_current = do_kpkt_bb(pkt_ptr, t_current);
        else if (pkt_type == TYPE_KPKT)
          t_current = do_kpkt(pkt_ptr, t_current, t2, nts);
        else
        {
          printout("kpkt not of type TYPE_KPKT or TYPE_PRE_KPKT\n");
          abort();
          //t_change_type = do_kpkt_ffonly(pkt_ptr, t_current, t2);
        }
        break;

      case TYPE_MA:
        // It's an active macroatom - apply transition probabilities
        //printout("MA-packet handling\n");

        t_current = do_macroatom(pkt_ptr, t_current, t2, nts);
        break;

      default:
        printout("Unknown packet type - abort\n");
        abort();
    }
  }
  assert(t_current == PACKET_SAME);
}


static void update_pellet(
  PKT *restrict pkt_ptr, const bool decay_to_kpkt, const bool decay_to_ntlepton, const int nts, const double ts, const double tw)
{
  // Handle inactive pellets. Need to do two things (a) check if it
  // decays in this time step and if it does handle that. (b) if it doesn't decay in
  // this time step then just move the packet along with the matter for the
  // start of the next time step.

  assert(!decay_to_kpkt || !decay_to_ntlepton); // can't decay to both!

  const double tdecay = pkt_ptr->tdecay; // after packet_init(), this value never changes
  if (tdecay > (ts + tw))
  {
    // It won't decay in this timestep, so just need to move it on with the flow.
    vec_scale(pkt_ptr->pos, (ts + tw) / ts);

    // That's all that needs to be done for the inactive pellet.
  }
  else if (tdecay > ts)
  {
    // The packet decays in the current timestep.
    if (decay_to_kpkt)
    {
      vec_scale(pkt_ptr->pos, tdecay / ts);

      pkt_ptr->type = TYPE_KPKT;
      pkt_ptr->absorptiontype = -6;
      packet_prop(pkt_ptr, tdecay, ts + tw, nts);
    }
    else if (decay_to_ntlepton)
    {
      vec_scale(pkt_ptr->pos, tdecay / ts);

      pkt_ptr->type = TYPE_NTLEPTON;
      pkt_ptr->absorptiontype = -10;
      place_ntlepton(pkt_ptr, tdecay);
      packet_prop(pkt_ptr, tdecay, ts + tw, nts);
    }
    else
    {
      // decay to gamma ray
      pellet_decay(nts, pkt_ptr);
      //printout("pellet to photon packet and propagation by packet_prop\n");
      packet_prop(pkt_ptr, tdecay, ts + tw, nts);
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
    packet_prop(pkt_ptr, tmin, ts + tw, nts);
  }
  else
  {
    printout("Something gone wrong with decaying pellets. Abort.\n");
    abort();
  }
}


static int compare_packets_bymodelgriddensity(const void *restrict p1, const void *restrict p2)
/// Helper function to sort the phixslist by descending cell density.
{
  const int a1_where = ((PKT *)p1)->where;
  const int a2_where = ((PKT *)p2)->where;

  const double rho_diff = get_rho(cell[a1_where].modelgridindex) - get_rho(cell[a2_where].modelgridindex);
  if (rho_diff < 0)
    return 1;
  else if (rho_diff > 0)
    return -1;
  else
    return 0;
}


void update_packets(const int nts, PKT *pkt)
/** Subroutine to move the packets and update them during the currect timestep. */
///nts the time step we're doing
{
  //double n_1;

  /** At the start, the packets have all either just been initialised or have already been
  processed for one or more timesteps. Those that are pellets will just be sitting in the
  matter. Those that are photons (or one sort or another) will already have a position and
  a direction.*/

  const double ts = time_step[nts].start;
  const double tw = time_step[nts].width;

  /*for (n = 0; n < npkts; n++)
  {
    printout("pkt[%d].where = %d\n",n,pkt[n].where);
  }*/

  /// If we want to do that with this version, sorting should be by modelgridcell
  //qsort(pkt,npkts,sizeof(PKT),compare_packets_bymodelgridposition);
  /// For 2D and 3D models sorting by the modelgrid cell's density should be most efficient
  qsort(pkt, npkts, sizeof(PKT), compare_packets_bymodelgriddensity);
  /*for (n = 0; n < npkts; n++)
  {
    printout("pkt[%d].where = %d, mgi %d\n",n,pkt[n].where,cell[pkt[n].where].modelgridindex);
  }*/

  //printout("before update packets\n");

  printout("start of parallel update_packets loop %d\n",time(NULL));
  /// Initialise the OpenMP reduction target to zero
  #ifdef _OPENMP
    #pragma omp parallel
    //copyin(debuglevel,nuJ,J)
  #endif
  {
    // time_t time_of_last_packet_printout = 0;
    #ifdef _OPENMP
    #pragma omp for schedule(dynamic) reduction(+:escounter,resonancescatterings,cellcrossings,nesc,updatecellcounter,coolingratecalccounter,upscatter,downscatter,ma_stat_activation_collexc,ma_stat_activation_collion,ma_stat_activation_ntcollexc,ma_stat_activation_ntcollion,ma_stat_activation_bb,ma_stat_activation_bf,ma_stat_deactivation_colldeexc,ma_stat_deactivation_collrecomb,ma_stat_deactivation_bb,ma_stat_deactivation_fb,k_stat_to_ma_collexc,k_stat_to_ma_collion,k_stat_to_r_ff,k_stat_to_r_fb,k_stat_from_ff,k_stat_from_bf,k_stat_from_gamma,k_stat_from_eminus,k_stat_from_earlierdecay)
    #endif
    for (int n = 0; n < npkts; n++)
    {
      // if ((time(NULL) - time_of_last_packet_printout > 1) || n == npkts - 1)
      if ((n % 10000 == 0) || n == npkts - 1)
      {
        // time_of_last_packet_printout = time(NULL);
        printout("[debug] update_packets: updating packet %d for timestep %d at time %d...\n",n,nts, time(NULL));
      }

      PKT *restrict pkt_ptr = &pkt[n];
      pkt_ptr->interactions = 0;
      //pkt_ptr->timestep = nts;

//        if (pkt_ptr->number == debug_packet) debuglevel = 2;
//        else debuglevel = 4;
      //if (n % 10000 == 0) debuglevel = 2;

      if (debuglevel == 2) printout("[debug] update_packets: updating packet %d for timestep %d __________________________\n",n,nts);

      const int cellindex = pkt_ptr->where;
      const int mgi = cell[cellindex].modelgridindex;
      /// for non empty cells update the global available level populations and cooling terms
      if (mgi != MMODELGRID)
      {
        //printout("thread%d _ pkt %d in cell %d with density %g\n",tid,n,pkt_ptr->where,cell[pkt_ptr->where].rho);
        /// Reset cellhistory if packet starts up in another than the last active cell
        if (cellhistory[tid].cellnumber != mgi)
          update_cell(mgi);

        /*
        histindex = 0;
        /// in the case of non isothermal or non homogenous grids this must be done for each packet
        /// unless we do some special treatment
        if (cellhistory[histindex].cellnumber != mgi)
        {
          histindex = search_cellhistory(mgi);
          if (histindex < 0)
          {
            //histindex = find_farthestcell_initial(mgi);
            histindex = 0;
            //printout("thread%d _ update_packets histindex %d\n",tid,histindex);
            update_cell(mgi);
          }
        }
        */

        //copy_populations_to_phixslist();
        /// rpkt's continuum opacity depends on nu, therefore it must be calculated by packet
        if (pkt_ptr->type == TYPE_RPKT && modelgrid[mgi].thick != 1)
        {
          calculate_kappa_rpkt_cont(pkt_ptr,ts);
        }
      }

      //printout("[debug] update_packets: current position of packet %d (%g, %g, %g)\n",n,pkt_ptr->pos[0],pkt_ptr->pos[1],pkt_ptr->pos[2]);
      //printout("[debug] update_packets: target position of homologous flow (%g, %g, %g)\n",pkt_ptr->pos[0]*(ts + tw)/ts,pkt_ptr->pos[1]*(ts + tw)/ts,pkt_ptr->pos[2]*(ts + tw)/ts);
      //printout("[debug] update_packets: current direction of packet %d (%g, %g, %g)\n",n,pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);

      switch (pkt_ptr->type)
      {
        case TYPE_56NI_PELLET:
        case TYPE_56CO_PELLET:
        case TYPE_57NI_PELLET:
        case TYPE_57CO_PELLET:
        case TYPE_48CR_PELLET:
        case TYPE_48V_PELLET:
          update_pellet(pkt_ptr, false, false, nts, ts, tw);
          break;

        case TYPE_52FE_PELLET:
        case TYPE_52MN_PELLET:
          // convert to kpts
          update_pellet(pkt_ptr, true, false, nts, ts, tw);
          break;

        case TYPE_57NI_POSITRON_PELLET:
        case TYPE_56CO_POSITRON_PELLET:
          // covert to to non-thermal leptons
          update_pellet(pkt_ptr, false, true, nts, ts, tw);
          break;

        case TYPE_GAMMA:
        case TYPE_RPKT:
        case TYPE_KPKT:
          /**Stuff for processing photons. */
          //printout("further propagate a photon packet via packet_prop\n");

          packet_prop(pkt_ptr, ts, ts + tw, nts);
          break;

        case TYPE_ESCAPE:
          break;

        default:
          printout("Unknown packet type %d %d. Abort.\n", pkt_ptr->type, n);
          abort();
      }

      if (debuglevel == 10 || debuglevel == 2)
        printout("[debug] update_packets: packet %d had %d interactions during timestep %d\n",
                 n, pkt_ptr->interactions, nts);

      if (n == npkts - 1)
        printout("last packet updated at %d\n",time(NULL));
    }
  }

  printout("end of update_packets parallel for loop %d\n",time(NULL));
  //printout("[debug] update_packets: packet %d updated for timestep %d\n",n,nts);
}





///***************************************************************************/
/*int search_cellhistory(int cellnumber)
/// Returns the historyindex of cellnumber if available in cellhistory.
/// Otherwise a negative value is returned as not found flag.
{
  int i;

  for (i = 0; i < CELLHISTORYSIZE; i++)
  {
    //printout("search_cellhistory: cellnumber %d in stack entry %d\n",cellhistory[i].cellnumber,i);
    if (cellnumber == cellhistory[i].cellnumber) return i;
  }

  return -99;
}*/


///***************************************************************************/
/*int find_farthestcell_initial(int cellnumber)
/// Searches the cellhistory for the cell with largest distance to cellnumber
/// and returns this cells historyindex.
/// This version may replace any element and should be used only for calls
/// from update_packets.
{
  int i,n,overwrite_histindex;
  double previous_distance,current_distance;

  n = cellhistory[0].cellnumber;
  overwrite_histindex = 0;
  if (n >= 0)
  {
    previous_distance = pow(cell[n].pos_init[0] - cell[cellnumber].pos_init[0], 2) + pow(cell[n].pos_init[1] - cell[cellnumber].pos_init[1], 2) + pow(cell[n].pos_init[2] - cell[cellnumber].pos_init[2], 2);

    for (i = 1; i < CELLHISTORYSIZE; i++)
    {
      n = cellhistory[i].cellnumber;
      if (n >= 0)
      {
        current_distance = pow(cell[n].pos_init[0] - cell[cellnumber].pos_init[0], 2) + pow(cell[n].pos_init[1] - cell[cellnumber].pos_init[1], 2) + pow(cell[n].pos_init[2] - cell[cellnumber].pos_init[2], 2);
        if (current_distance > previous_distance)
        {
          previous_distance = current_distance;
          overwrite_histindex = i;
        }
      }
      else
      {
        overwrite_histindex = i;
        break;
      }
    }
  }

  return overwrite_histindex;
}*/


///***************************************************************************/
/*int find_farthestcell(int cellnumber)
/// Searches the cellhistory for the cell with largest distance to cellnumber
/// and returns this cells historyindex.
/// This version which keeps the zeroth element of cellhistory should be used
/// for calls from change_cell.
{
  int i,n,overwrite_histindex;
  double previous_distance,current_distance;

  n = cellhistory[1].cellnumber;
  overwrite_histindex = 1;
  if (n >= 0)
  {
    previous_distance = pow(cell[n].pos_init[0] - cell[cellnumber].pos_init[0], 2) + pow(cell[n].pos_init[1] - cell[cellnumber].pos_init[1], 2) + pow(cell[n].pos_init[2] - cell[cellnumber].pos_init[2], 2);

    for (i = 2; i < CELLHISTORYSIZE; i++)
    {
      n = cellhistory[i].cellnumber;
      if (n >= 0)
      {
        current_distance = pow(cell[n].pos_init[0] - cell[cellnumber].pos_init[0], 2) + pow(cell[n].pos_init[1] - cell[cellnumber].pos_init[1], 2) + pow(cell[n].pos_init[2] - cell[cellnumber].pos_init[2], 2);
        if (current_distance > previous_distance)
        {
          previous_distance = current_distance;
          overwrite_histindex = i;
        }
      }
      else
      {
        overwrite_histindex = i;
        break;
      }
    }
  }

  return overwrite_histindex;
}*/


///***************************************************************************/
/*void copy_populations_to_phixslist()
{
  //double nnlevel;
  int i;
  for (i=0; i < importantbfcontinua; i++)
  {
    phixslist[i].nnlevel = get_levelpop(phixslist[i].element, phixslist[i].ion, phixslist[i].level);
    //nnlevel = get_levelpop(phixslist[i].element, phixslist[i].ion, phixslist[i].level);
    //printout("phixlistnnlevel i%d, element %d, ion %d, level %d, nnlevel %g\n",i,phixslist[i].element,phixslist[i].ion,phixslist[i].level,nnlevel);
    //phixslist[i].nnlevel = nnlevel;
  }
}*/


/*static int compare_packets_byposition(const void *restrict p1, const void *restrict p2)
/// Helper function to sort the phixslist by ascending threshold frequency.
{
  const PKT *restrict a1 = (PKT *)(p1);
  const PKT *restrict a2 = (PKT *)(p2);

  int cell_diff = a1->where - a2->where;
  if (cell_diff < 0)
    return -1;
  else if (cell_diff > 0)
    return 1;
  else
    return 0;
}


static int compare_packets_bymodelgridposition(const void *restrict p1, const void *restrict p2)
/// Helper function to sort the phixslist by ascending threshold frequency.
{
  const PKT *restrict a1 = (PKT *)(p1);
  const PKT *restrict a2 = (PKT *)(p2);

  int mgi_diff = cell[a1->where].modelgridindex - cell[a2->where].modelgridindex;
  if (mgi_diff < 0)
    return -1;
  else if (mgi_diff > 0)
    return 1;
  else
    return 0;
}*/
