#include "sn3d.h"
#include "atomic.h"
#include "gamma.h"
#include "kpkt.h"
#include "ltepop.h"
#include "packet_prop.h"
#include "update_packets.h"
#include "rpkt.h"


//private functions
int compare_packets_byposition(const void *p1, const void *p2);
int compare_packets_bymodelgridposition(const void *p1, const void *p2);
int compare_packets_bymodelgriddensity(const void *p1, const void *p2);


/** Subroutine to move the packets and update them during the currect timestep. */

int update_packets(int nts)
///nts the time step we're doing
{
  //double n_1;

  /** At the start, the packets have all either just been initialised or have already been
  processed for one or more timesteps. Those that are pellets will just be sitting in the
  matter. Those that are photons (or one sort or another) will already have a position and
  a direction.*/

  double ts = time_step[nts].start;
  double tw = time_step[nts].width;

  /*for (n = 0; n < npkts; n++)
  {
    printout("pkt[%d].where = %d\n",n,pkt[n].where);
  }*/

  /// If we want to do that with this version, sorting should be by modelgridcell
  //qsort(pkt,npkts,sizeof(PKT),compare_packets_bymodelgridposition);
  /// For 2D and 3D models sorting by the modelgrid cell's density should be most efficient
  qsort(pkt,npkts,sizeof(PKT),compare_packets_bymodelgriddensity);
  /*for (n = 0; n < npkts; n++)
  {
    printout("pkt[%d].where = %d, mgi %d\n",n,pkt[n].where,cell[pkt[n].where].modelgridindex);
  }*/

  //printout("before update packets\n");

  printout("start of parallel update_packets loop %d\n",time(NULL));
  int time_of_last_packet_printout = 0;
  /// Initialise the OpenMP reduction target to zero
  #ifdef _OPENMP
    #pragma omp parallel
    //copyin(debuglevel,nuJ,J)
    {
      #pragma omp for schedule(dynamic) reduction(+:escounter,resonancescatterings,cellcrossings,nesc,updatecellcounter,coolingratecalccounter,upscatter,downscatter,ma_stat_activation_collexc,ma_stat_activation_collion,ma_stat_activation_bb,ma_stat_activation_bf,ma_stat_deactivation_colldeexc,ma_stat_deactivation_collrecomb,ma_stat_deactivation_bb,ma_stat_deactivation_fb,k_stat_to_ma_collexc,k_stat_to_ma_collion,k_stat_to_r_ff,k_stat_to_r_fb,k_stat_from_ff,k_stat_from_bf,k_stat_from_gamma,k_stat_from_eminus,k_stat_from_earlierdecay)
  #endif
      for (int n = 0; n < npkts; n++)
      {
        //printout("[debug] update_packets: updating packet %d for timestep %d...\n",n,nts);
        if (time(NULL) - time_of_last_packet_printout > 2)
        {
          time_of_last_packet_printout = time(NULL);
          printout("[debug] update_packets: updating packet %d for timestep %d...\n",n,nts);
        }
        //if (n == 5000) exit(0);

        PKT *restrict pkt_ptr = &pkt[n];
        pkt_ptr->interactions = 0;
        //pkt_ptr->timestep = nts;

//        if (pkt_ptr->number == debug_packet) debuglevel = 2;
//        else debuglevel = 4;
        //if (n % 10000 == 0) debuglevel = 2;

        if (debuglevel == 2) printout("[debug] update_packets: updating packet %d for timestep %d __________________________\n",n,nts);

        int mgi = cell[pkt_ptr->where].modelgridindex;
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

        if ((pkt_ptr->type == TYPE_NICKEL_PELLET) || (pkt_ptr->type == TYPE_COBALT_PELLET) || (pkt_ptr->type == TYPE_48CR_PELLET) || (pkt_ptr->type == TYPE_48V_PELLET)|| (pkt_ptr->type == TYPE_52FE_PELLET) || (pkt_ptr->type == TYPE_52MN_PELLET) || (pkt_ptr->type == TYPE_COBALT_POSITRON_PELLET))
        {
          //printout("inactive pellet\n");
          /**It's still an inactive pellet. Need to do two things (a) check if it
          decays in this time step and if it does handle that. (b) if it doesn't decay in
          this time step then just move the packet along with the matter for the
          start of the next time step. */

          if (pkt_ptr->tdecay > (ts + tw))
          {
            /**It won't decay in this timestep just need to move it on.*/

            pkt_ptr->pos[0] = pkt_ptr->pos[0] * (ts + tw) / ts;
            pkt_ptr->pos[1] = pkt_ptr->pos[1] * (ts + tw) / ts;
            pkt_ptr->pos[2] = pkt_ptr->pos[2] * (ts + tw) / ts;

            /**That's all that needs to be done for the inactive pellet. */
          }
          else if (pkt_ptr->tdecay > ts)
          {
            /**These are the packets decaying in this timestep.*/
            if ((pkt_ptr->type == TYPE_52FE_PELLET) || (pkt_ptr->type == TYPE_52MN_PELLET) || (pkt_ptr->type == TYPE_COBALT_POSITRON_PELLET))
            {
              pkt_ptr->pos[0] = pkt_ptr->pos[0] * pkt_ptr->tdecay / ts;
              pkt_ptr->pos[1] = pkt_ptr->pos[1] * pkt_ptr->tdecay / ts;
              pkt_ptr->pos[2] = pkt_ptr->pos[2] * pkt_ptr->tdecay / ts;

              pkt_ptr->type = TYPE_KPKT;
              pkt_ptr->absorptiontype = -6;
              packet_prop(pkt_ptr,pkt_ptr->tdecay,ts+tw,nts);
            }
            else
            {
              pellet_decay(nts,pkt_ptr);
              //printout("pellet to photon packet and propagation by packet_prop\n");
              packet_prop(pkt_ptr,pkt_ptr->tdecay,ts+tw,nts);
            }
          }
          else if ((pkt_ptr->tdecay > 0) && (nts == 0))
          {
            /** These are pellets whose decay times were before the first time step
            They will be made into r-packets with energy reduced for doing work on the
            ejects following Lucy 2004. */
            /** The position is already set at tmin so don't need to move it. Assume
            that it is fixed in place from decay to tmin - i.e. short mfp. */

            pkt_ptr->e_cmf = pkt_ptr->e_cmf *pkt_ptr->tdecay/tmin;
            //pkt_ptr->type = TYPE_KPKT;
            pkt_ptr->type = TYPE_PRE_KPKT;
            pkt_ptr->absorptiontype = -7;
            //if (tid == 0) k_stat_from_earlierdecay += 1;
            k_stat_from_earlierdecay += 1;

            //printout("already decayed packets and propagation by packet_prop\n");
            packet_prop(pkt_ptr, tmin, ts+tw, nts);
          }
          else
          {
            printout("Something gone wrong with decaying pellets. Abort.\n");
            exit(0);
          }
        }
        else if (pkt_ptr->type == TYPE_GAMMA || pkt_ptr->type == TYPE_RPKT || pkt_ptr->type == TYPE_KPKT)
        {
          /**Stuff for processing photons. */
          //printout("further propagate a photon packet via packet_prop\n");

          packet_prop(pkt_ptr, ts, ts+tw, nts);

        }
        else if (pkt_ptr->type != TYPE_ESCAPE)
        {
          printout("Unknown packet type %d %d. Abort.\n", pkt_ptr->type, n);
          exit(0);
        }

        if (debuglevel == 10 || debuglevel == 2) printout("[debug] update_packets: packet %d had %d interactions during timestep %d\n",n,pkt_ptr->interactions,nts);

        if (n == npkts-1) printout("last packet updated at %d\n",time(NULL));

      }
  #ifdef _OPENMP
    }
  #endif

  printout("end of update_packets parallel for loop %d\n",time(NULL));
  //printout("[debug] update_packets: packet %d updated for timestep %d\n",n,nts);
  return 0;
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
void update_cell(int cellnumber)
///=calculate_levelpops for non isothermal homogeneous grids
///
{
  updatecellcounter += 1;

  /// Make known that cellhistory[tid] contains information about the
  /// cell given by cellnumber.
  cellhistory[tid].cellnumber = cellnumber;

  /// Calculate the level populations for this cell, and flag the other entries
  /// as empty.
  //printout("update cell %d at histindex %d\n",cellnumber,histindex);
  for (int element = 0; element < nelements; element++)
  {
    int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      cellhistory[tid].coolinglist[get_coolinglistoffset(element,ion)].contribution = COOLING_UNDEFINED;
      int nlevels = get_nlevels(element,ion);
      for (int level = 0; level < nlevels; level++)
      {
        double population = calculate_exclevelpop(cellnumber,element,ion,level);
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].population = population;

        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
        {
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].sahafact = -99.;
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].spontaneousrecombrate = -99.;
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].bfcooling = -99.;
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].bfheating = -99.;
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].corrphotoioncoeff = -99.;
        }

        /// This is the only flag needed for all of the following MA stuff!
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_deexc = -99.;
        /*
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_deexc = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_recomb = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_recomb = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_same = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_same = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_lower = -99.;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_higher = -99.;

        ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
        nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
        for (i = 0; i < ndowntrans; i++)
        {
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc[i] = -99.;
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[i] = -99.;
        }
        for (i = 0; i < nuptrans; i++)
        {
          cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same[i] = -99.;
        }
        */
      }
    }
  }
  //cellhistory[tid].totalcooling = COOLING_UNDEFINED;
  //cellhistory[tid].phixsflag = PHIXS_UNDEFINED;
}

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


///****************************************************************************
int compare_packets_byposition(const void *p1, const void *p2)
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


///****************************************************************************
int compare_packets_bymodelgridposition(const void *p1, const void *p2)
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
}


///****************************************************************************
int compare_packets_bymodelgriddensity(const void *p1, const void *p2)
/// Helper function to sort the phixslist by descending cell density.
{
  const PKT *a1 = (PKT *)(p1);
  const PKT *a2 = (PKT *)(p2);

  double rho_diff = modelgrid[cell[a1->where].modelgridindex].rho - modelgrid[cell[a2->where].modelgridindex].rho;
  if (rho_diff < 0)
    return 1;
  else if (rho_diff > 0)
    return -1;
  else
    return 0;
}
