#include "sn3d.h"

/* Master routine for moving packets around. When it called,
   it is given the time at start of inverval and at end - when it finishes,
   everything the packet does during this time should be sorted out. */

int packet_prop(pkt_ptr, t1, t2, nts)
     PKT *pkt_ptr;
     double t1, t2;
     int nts; //the time step we are doing
{
  double t_current; 
  double do_gamma();
  double do_rpkt();
  double do_rpkt_thickcell();
  double do_kpkt();
  double do_kpkt_bb();
  double do_kpkt_ffonly();
  double do_ma();
  int end_packet;
  double t_change_type;

  end_packet = 0;  //means "keep working"

  t_current = t1;

  /* 0 the scatter counter for the packet. */
  pkt_ptr->scat_count = 0;


  while (end_packet == 0)
  {
    /* Start by sorting out what sort of packet it is.*/
    //printout("start of packet_prop loop %d\n", pkt_ptr->type );
    if (pkt_ptr->type == TYPE_GAMMA)
    {
      /*It's a gamma-ray packet.*/
      /* Call do_gamma. */
      //printout("gamma propagation\n");
      t_change_type = do_gamma( pkt_ptr, t_current, t2);
	  /* This returns a flag if the packet gets to t2 without
      changing to something else. If the packet does change it
      returns the time of change and sets everything for the
      new packet.*/
      if (t_change_type < 0)
      {
        end_packet = 1;
      }
      else
      {
        #ifdef _OPENMP 
          #pragma omp atomic
        #endif
        time_step[nts].gamma_dep += pkt_ptr->e_cmf; ///***??
        t_current = t_change_type;
      }
    }
    else if (pkt_ptr->type == TYPE_RPKT)
    {
      /*It's an r-packet. */
      //printout("r-pkt propagation\n");
      t_change_type = do_rpkt( pkt_ptr, t_current, t2);
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
	  
      if (t_change_type < 0)
      {
        end_packet = 1;
      }
      else
      {
        t_current = t_change_type;
      }
    }
    else if (pkt_ptr->type == TYPE_EMINUS)
    {
      /*It's an electron - convert to k-packet*/
      //printout("e-minus propagation\n");
      pkt_ptr->type = TYPE_KPKT;
      #ifndef FORCE_LTE
        //kgammadep[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif
      //pkt_ptr->type = TYPE_PRE_KPKT;
      //pkt_ptr->type = TYPE_GAMMA_KPKT;
      //if (tid == 0) k_stat_from_eminus += 1;
      k_stat_from_eminus += 1;
    }
    else if (pkt_ptr->type == TYPE_KPKT || pkt_ptr->type == TYPE_PRE_KPKT || pkt_ptr->type == TYPE_GAMMA_KPKT)
    {
      /*It's a k-packet - convert to r-packet (low freq).*/
      //printout("k-packet propagation\n");
	  
      //t_change_type = do_kpkt(pkt_ptr, t_current, t2);
      
      if (pkt_ptr->type == TYPE_PRE_KPKT || modelgrid[cell[pkt_ptr->where].modelgridindex].thick == 1)
        t_change_type = do_kpkt_bb(pkt_ptr, t_current, t2);
      else if (pkt_ptr->type == TYPE_KPKT)
        t_change_type = do_kpkt(pkt_ptr, t_current, t2, nts);
      else
      {
        printout("kpkt not of type TYPE_KPKT or TYPE_PRE_KPKT\n");
        //t_change_type = do_kpkt_ffonly(pkt_ptr, t_current, t2);
      }
      
      if (t_change_type < 0)
      {
        end_packet = 1;
      }
      else
      {
        t_current = t_change_type;
      }
    }
    else if (pkt_ptr->type == TYPE_MA)
    {
      /*It's an active macroatom - apply transition probabilities*/
      //printout("MA-packet handling\n");
        
      t_change_type = do_ma(pkt_ptr, t_current, t2, nts);

      if (t_change_type < 0)
      {
        /// shouldn't be possible in this case
        end_packet = 1;
      }
      else
      {
        /// should be the same as original t_current, as MA processes 
        /// happen instantaneously so far
        t_current = t_change_type;
      }
    }
    else
    {
      printout("Unknown packet type - abort\n");
      exit(0);
    }

  }

  return(0);
}
