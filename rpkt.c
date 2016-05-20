#include "sn3d.h"
#include "atomic.h"
#include "boundary.h"
#include "grey_emissivities.h"
#include "grid_init.h"
#include "ltepop.h"
#include "move.h"
#include "polarization.h"
#include "rpkt.h"
#include "update_grid.h"
#include "vectors.h"


// private functions
double get_event(int modelgridindex, PKT *pkt_ptr, int *rpkt_eventtype, double t_current, double tau_rnd, double abort_dist);
int rpkt_event(PKT *pkt_ptr, int rpkt_eventtype, double t_current);
void rpkt_event_thickcell(PKT *pkt_ptr, double t_current);

static inline
double min(double a, double b)
// returns minimum of a and b
{
  if (a >= b)
  {
    return(b);
  }
  else
  {
    return(a);
  }
}


// Material for handing r-packet propagation.
///****************************************************************************
double do_rpkt(PKT *pkt_ptr, double t1, double t2)
/** Routine for moving an r-packet. Similar to do_gamma in objective.*/
{
  double vel_vec[3];

  double zrand, tau_next, tau_current, tdist;
  double sdist, edist;
  int snext;
  int find_nextline;
  int rpkt_eventtype;

  double kappa, kappa_cont, sigma, kappa_ff, kappa_bf;
  double *kappacont_ptr, *sigma_ptr, *kappaff_ptr, *kappabf_ptr;
  kappacont_ptr = &kappa_cont;
  sigma_ptr = &sigma;
  kappaff_ptr = &kappa_ff;
  kappabf_ptr = &kappa_bf;

  int mgi = cell[pkt_ptr->where].modelgridindex;

  double t_current = t1; ///this will keep track of time in the calculation
  //printout("[debug] r-pkt propagation init\n");
  //int it = 1;

  int end_packet = 0; ///means "keep working"
  while (end_packet == 0)
  {
    find_nextline = 0;

    /*
    int i;
    for (i=0; i < CELLHISTORYSIZE; i++)
      printout("thread%d _ do_rpkt: cellhistory[%d].cellnumber = %d\n",tid,i,cellhistory[i].cellnumber);
    */

    #ifdef DEBUG_ON
      if (pkt_ptr->next_trans > 0)
      {
        if (debuglevel == 2) printout("[debug] do_rpkt: init: pkt_ptr->nu_cmf %g, nu(pkt_ptr->next_trans=%d) %g, nu(pkt_ptr->next_trans-1=%d) %g, pkt_ptr->where %d\n", pkt_ptr->nu_cmf, pkt_ptr->next_trans, linelist[pkt_ptr->next_trans].nu, pkt_ptr->next_trans-1, linelist[pkt_ptr->next_trans-1].nu, pkt_ptr->where );
        if (debuglevel == 2) printout("[debug] do_rpkt: init: (pkt_ptr->nu_cmf - nu(pkt_ptr->next_trans-1))/pkt_ptr->nu_cmf %g\n", (pkt_ptr->nu_cmf-linelist[pkt_ptr->next_trans-1].nu)/pkt_ptr->nu_cmf);
      }
    #endif

    //printout("[debug] r-pkt propagation iteration %d\n",it);
    //it += 1;
    /** Assign optical depth to next physical event. And start counter of
    optical depth for this path.*/
    zrand = gsl_rng_uniform(rng);
    tau_next = -1. * log(zrand);
    tau_current = 0.0;

    /** Start by finding the distance to the crossing of the grid cell
    boundaries. sdist is the boundary distance and snext is the
    grid cell into which we pass.*/
    sdist = boundary_cross(pkt_ptr, t_current, &snext);

    if (sdist == 0)
    {
      change_cell(pkt_ptr, snext, &end_packet, t_current);
      mgi = cell[pkt_ptr->where].modelgridindex;
    }
    else
    {
      if (sdist > (rmax * t_current/tmin))
      {
        printout("[fatal] do_rpkt: Unreasonably large sdist. Rpkt. Abort. %g %g %g\n", rmax, t_current/tmin, sdist);
        abort();
      }

      if (sdist < 1)
      {
        printout("[warning] r_pkt: Negative distance (sdist = %g). Abort.\n", sdist);
        printout("[warning] r_pkt: cell %d snext %d\n", pkt_ptr->where, snext);
        printout("[warning] r_pkt: pos %g %g %g\n", pkt_ptr->pos[0], pkt_ptr->pos[1], pkt_ptr->pos[2]);
        printout("[warning] r_pkt: dir %g %g %g\n", pkt_ptr->dir[0], pkt_ptr->dir[1], pkt_ptr->dir[2]);
        printout("[warning] r_pkt: cell corner %g %g %g\n",cell[pkt_ptr->where].pos_init[0]*t_current/tmin, cell[pkt_ptr->where].pos_init[1]*t_current/tmin,  cell[pkt_ptr->where].pos_init[2]*t_current/tmin);
        printout("[warning] r_pkt: cell width %g %g %g\n",wid_init*t_current/tmin, wid_init*t_current/tmin,  wid_init*t_current/tmin);
        //abort();
      }
      if (((snext != -99) && (snext < 0)) || (snext >= ngrid))
      {
        printout("[fatal] r_pkt: Heading for inappropriate grid cell. Abort.\n");
        printout("[fatal] r_pkt: Current cell %d, target cell %d.\n", pkt_ptr->where, snext);
        abort();
      }
      if (sdist > max_path_step)
      {
        sdist = max_path_step;
        snext = pkt_ptr->where;
      }


      /** At present there is no scattering/destruction process so all that needs to
      happen is that we determine whether the packet reaches the boundary during the timestep. */

      /** Find how far it can travel during the time inverval. */

      tdist = (t2 - t_current) * CLIGHT_PROP;

      if (tdist < 0)
      {
        printout("[fatal] do_rpkt: Negative distance (tdist). Abort. \n");
        abort();
      }

      //if (cell[pkt_ptr->where].nne < 1e-40)
      //if (get_nne(cell[pkt_ptr->where].modelgridindex) < 1e-40)
      if (mgi == MMODELGRID)
      {
        /// for empty cells no physical event occurs. The packets just propagate.
        edist = 1e99;
        find_nextline = 1;
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_rpkt: propagating through empty cell, set edist=1e99\n");
        #endif
      }
      else if (modelgrid[mgi].thick == 1)
      {
        /// In the case ot optically thick cells, we treat the packets in grey approximation to speed up the calculation
        /// Get distance to the next physical event in this case only electron scattering
        //kappa = SIGMA_T*get_nne(mgi);
        kappa = get_kappagrey(mgi)*get_rho(mgi);
        get_velocity(pkt_ptr->pos, vel_vec, t_current);
        kappa = kappa * doppler(pkt_ptr->dir, vel_vec);
        edist = (tau_next - tau_current) / kappa;
        find_nextline = 1;
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_rpkt: propagating through grey cell, edist  %g\n",edist);
        #endif
      }
      else
      {
        /// get distance to the next physical event (continuum or bound-bound)
        edist = get_event(mgi,pkt_ptr, &rpkt_eventtype, t_current, tau_next, min(tdist,sdist)); //, kappacont_ptr, sigma_ptr, kappaff_ptr, kappabf_ptr);
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_rpkt: after edist: pkt_ptr->nu_cmf %g, nu(pkt_ptr->next_trans=%d) %g\n", pkt_ptr->nu_cmf, pkt_ptr->next_trans, linelist[pkt_ptr->next_trans].nu);
        #endif
      }
      if (edist < 0)
      {
        printout("[fatal] do_rpkt: Negative distance (edist). Abort. \n");
        printout("[fatal] do_rpkt: Trouble was due to packet number %d.\n", pkt_ptr->number);
        abort();
      }
      //printout("[debug] do_rpkt: sdist, tdist, edist %g, %g, %g\n",sdist,tdist,edist);

      if ((sdist < tdist) && (sdist < edist))
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_rpkt: sdist < tdist && sdist < edist\n");
        #endif
        /** Move it into the new cell. */
        sdist = sdist / 2.;
        t_current += sdist / CLIGHT_PROP;
        move_pkt(pkt_ptr,sdist,t_current);
        update_estimators(pkt_ptr,sdist*2);
        if (do_rlc_est != 0 && do_rlc_est != 3)
        {
          sdist = sdist * 2.;
          rlc_emiss_rpkt(pkt_ptr, sdist, t_current);
          sdist = sdist / 2.;
        }
        t_current += sdist / CLIGHT_PROP;
        move_pkt(pkt_ptr,sdist,t_current);
        sdist = sdist * 2.;

        if (snext != pkt_ptr->where)
        {
          change_cell(pkt_ptr, snext, &end_packet, t_current);
          mgi = cell[pkt_ptr->where].modelgridindex;
        }
        #ifdef DEBUG_ON
          /** New cell so reset the scat_counter */
          pkt_ptr->scat_count = 0;
          //if (debuglevel == 2) printout("[debug] do_rpkt:   pkt_ptr->last_event %d\n",pkt_ptr->last_event);
          pkt_ptr->last_event = pkt_ptr->last_event + 100;
        #endif
        /// For empty or grey cells a photon can travel over several bb-lines. Thus we need to
        /// find the next possible line interaction.
        if (find_nextline == 1)
        {
          /// However, this is only required if the new cell is non-empty or non-grey
          if (mgi != MMODELGRID && modelgrid[mgi].thick != 1)
          {
            closest_transition_empty(pkt_ptr);
          }
        }
      }
      else if ((tdist < sdist) && (tdist < edist))
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_rpkt: tdist < sdist && tdist < edist\n");
        #endif
        /** Doesn't reach boundary. */
        tdist = tdist / 2.;
        t_current += tdist / CLIGHT_PROP;
        move_pkt(pkt_ptr,tdist,t_current);
        update_estimators(pkt_ptr,tdist*2);
        if (do_rlc_est != 0 && do_rlc_est != 3)
        {
          tdist = tdist * 2.;
          rlc_emiss_rpkt(pkt_ptr, tdist, t_current);
          tdist = tdist / 2.;
        }
        t_current = t2;
        move_pkt(pkt_ptr,tdist,t_current);
        tdist = tdist * 2.;
        #ifdef DEBUG_ON
          pkt_ptr->last_event = pkt_ptr->last_event + 1000;
        #endif
        /// For empty or grey cells a photon can travel over several bb-lines. Thus we need to
        /// find the next possible line interaction.
        if (find_nextline == 1) closest_transition_empty(pkt_ptr);
        end_packet = 1;
      }
      else if ((edist < sdist) && (edist < tdist))
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_rpkt: edist < sdist && edist < tdist\n");
        #endif
        edist = edist / 2.;
        t_current += edist / CLIGHT_PROP;
        move_pkt(pkt_ptr,edist,t_current);
        update_estimators(pkt_ptr,edist*2);
        if (do_rlc_est != 0 && do_rlc_est != 3)
        {
          edist = edist * 2.;
          rlc_emiss_rpkt(pkt_ptr, edist, t_current);
          edist = edist / 2.;
        }
        t_current += edist / CLIGHT_PROP;
        move_pkt(pkt_ptr,edist,t_current);
        edist = edist * 2.;

        /** The previously selected and in pkt_ptr stored event occurs. Handling is done by rpkt_event*/
        if (modelgrid[mgi].thick == 1)
          rpkt_event_thickcell(pkt_ptr,t_current);
        else
          rpkt_event(pkt_ptr,rpkt_eventtype,t_current);
        if (pkt_ptr->type != TYPE_RPKT)
        {
          /** It's not an r-packet any more - return.*/
          return(t_current);
        }
      }
      else
      {
        printout("[fatal] do_rpkt: Failed to identify event . Rpkt. edist %g, sdist %g, tdist %g Abort.\n", edist, sdist, tdist);
        printout("[fatal] do_rpkt: Trouble was due to packet number %d.\n", pkt_ptr->number);
        abort();
      }
    }
  }

  return PACKET_SAME;
}



///****************************************************************************
double get_event(int modelgridindex, PKT *pkt_ptr, int *rpkt_eventtype, double t_current, double tau_rnd, double abort_dist)
///     PKT *pkt_ptr;      //pointer to packet object
///     double t_current;  //current time
///     double tau_rnd;    //random optical depth until which the packet travels
///     double abort_dist;      //maximal travel distance before packet leaves cell or time step ends
/// BE AWARE THAT THIS PROCEDURE SHOULD BE ONLY CALLED FOR NON EMPTY CELLS!!!
{
  /// initialize loop variables
  double tau = 0.;        ///initial optical depth along path
  double dist = 0.;       ///initial position on path
  double edist = 0.;

  PKT dummypkt;
  dummypkt = *pkt_ptr;
  PKT *dummypkt_ptr;
  dummypkt_ptr = &dummypkt;
  //propagationcounter = 0;
  int endloop = 0;
  while (endloop == 0)
  {
    /// calculate distance to next line encounter ldist
    /// first select the closest transition in frequency
    /// we need its frequency nu_trans, the element/ion and the corresponding levels
    /// create therefore new variables in packet, which contain next_lowerlevel, ...
    double nu_trans = closest_transition(dummypkt_ptr);  ///returns negative value if nu_cmf > nu_trans
    int element = mastate[tid].element;
    int ion = mastate[tid].ion;
    int upper = mastate[tid].level;
    int lower = mastate[tid].activatedfromlevel;
    #ifdef DEBUG_ON
      //if (debuglevel == 2) printout("[debug] get_event: propagationcounter %d\n",propagationcounter);
      if (debuglevel == 2) printout("[debug] get_event:   dummypkt_ptr->nu_cmf %g, dummypkt_ptr->nu_rf %g, dummypkt_ptr->where %d\n", dummypkt_ptr->nu_cmf,dummypkt_ptr->nu_rf,dummypkt_ptr->where);
      if (debuglevel == 2) printout("[debug] get_event:   closest_transition %g, pkt_ptr->next_trans %d, nu_next_trans %g\n",nu_trans,dummypkt_ptr->next_trans,linelist[dummypkt_ptr->next_trans].nu);
    #endif

    if (nu_trans > 0)
    {
      /// line interaction in principle possible (nu_cmf > nu_trans)
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] get_event:   line interaction possible\n");
      #endif

      double ldist;
      if (dummypkt_ptr->nu_cmf < nu_trans)
      {
        //printout("dummypkt_ptr->nu_cmf %g < nu_trans %g, next_trans %d, element %d, ion %d, lower%d, upper %d\n",dummypkt_ptr->nu_cmf,nu_trans,dummypkt_ptr->next_trans,element,ion,lower,upper);
        ldist = 0;  /// photon was propagated too far, make sure that we don't miss a line
      }
      else
        ldist = CLIGHT * t_current * (dummypkt_ptr->nu_cmf/nu_trans - 1);
      //fprintf(ldist_file,"%25.16e %25.16e\n",dummypkt_ptr->nu_cmf,ldist);
      if (ldist < 0.) printout("[warning] get_event: ldist < 0 %g\n",ldist);

      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] get_event:     ldist %g\n",ldist);
      #endif

      //A_ul = einstein_spontaneous_emission(element,ion,upper,lower);
      double A_ul = einstein_spontaneous_emission(dummypkt_ptr->next_trans-1);
      double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans,3) * A_ul;
      double B_lu = stat_weight(element,ion,upper)/stat_weight(element,ion,lower) * B_ul;

      double n_u = get_levelpop(modelgridindex,element,ion,upper);
      double n_l = get_levelpop(modelgridindex,element,ion,lower);

      double tau_line = (B_lu*n_l - B_ul*n_u) * HCLIGHTOVERFOURPI * t_current;
      //if (element == 7) fprintf(tau_file,"%g %g %d\n",nu_trans,tau_line,ion);
      if (tau_line < 0)
      {
        if (SILENT == 0) printout("[warning] get_event: tau_line %g < 0, n_l %g, n_u %g, B_lu %g, B_ul %g, W %g, T_R %g, element %d, ion %d, upper %d, lower %d ... abort\n",tau_line, n_l,n_u,B_lu,B_ul,get_W(cell[pkt_ptr->where].modelgridindex),get_TR(cell[pkt_ptr->where].modelgridindex),element,ion,upper,lower);
        if (SILENT == 0) printout("[warning] get_event: set tau_line = 0\n");
        tau_line = 0.;
        //printout("[fatal] get_event: tau_line < 0 ... abort\n");
        //abort();
      }

      //calculate_kappa_rpkt_cont(dummypkt_ptr, t_current);
      ///restore values which were changed by calculate_kappa_rpkt_cont to those set by closest_transition
      //mastate[tid].element = element;
      //mastate[tid].ion = ion;
      //mastate[tid].level = upper;
      double kap_cont = kappa_rpkt_cont[tid].total;
      double tau_cont = kap_cont*ldist;

      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] get_event:     tau_rnd %g, tau %g, tau_cont %g, tau_line %g\n",tau_rnd,tau,tau_cont,tau_line);
      #endif
      if (tau_rnd - tau > tau_cont)
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] get_event:       tau_rnd - tau > tau_cont\n");
        #endif
        if (tau_rnd - tau > tau_cont + tau_line)
        {
          /// total optical depth still below tau_rnd: continue
          #ifdef DEBUG_ON
            if (debuglevel == 2) printout("[debug] get_event:         tau_rnd - tau > tau_cont + tau_line ... proceed this packets propagation\n");
          #endif

          dist = dist + ldist;
          #ifdef DEBUG_ON
            if (debuglevel == 2) printout("[debug] get_event:         dist %g, abort_dist %g, dist-abort_dist %g\n", dist, abort_dist, dist-abort_dist);
          #endif
          if (dist > abort_dist)
          {
            dummypkt_ptr->next_trans -= 1;
            pkt_ptr->next_trans = dummypkt_ptr->next_trans;
            #ifdef DEBUG_ON
              if (debuglevel == 2) printout("[debug] get_event:         leave propagation loop (dist %g > abort_dist %g) ... dummypkt_ptr->next_trans %d\n", dist,abort_dist,dummypkt_ptr->next_trans);
            #endif
            return abort_dist+1e20;
          }
          tau = tau + tau_cont + tau_line;
          //dummypkt_ptr->next_trans += 1;
          t_current += ldist / CLIGHT_PROP;
          move_pkt(dummypkt_ptr,ldist,t_current);

          #ifdef DEBUG_ON
            if (debuglevel == 2) printout("[debug] get_event:         dummypkt_ptr->nu_cmf %g, nu(dummypkt_ptr->next_trans=%d) %g, nu(dummypkt_ptr->next_trans-1=%d) %g\n", dummypkt_ptr->nu_cmf, dummypkt_ptr->next_trans, linelist[dummypkt_ptr->next_trans].nu, dummypkt_ptr->next_trans-1, linelist[dummypkt_ptr->next_trans-1].nu);
            if (debuglevel == 2) printout("[debug] get_event:         (dummypkt_ptr->nu_cmf - nu(dummypkt_ptr->next_trans-1))/dummypkt_ptr->nu_cmf %g\n", (dummypkt_ptr->nu_cmf-linelist[dummypkt_ptr->next_trans-1].nu)/dummypkt_ptr->nu_cmf);

            if (debuglevel == 2)
            {
              if (dummypkt_ptr->nu_cmf >= linelist[dummypkt_ptr->next_trans].nu && dummypkt_ptr->nu_cmf < linelist[dummypkt_ptr->next_trans-1].nu)
                printout("[debug] get_event:           nu(next_trans-1) > nu_cmf >= nu(next_trans)\n");
              else if (dummypkt_ptr->nu_cmf < linelist[dummypkt_ptr->next_trans].nu)
                printout("[debug] get_event:           nu_cmf < nu(next_trans)\n");
              else
                printout("[debug] get_event:           nu_cmf >= nu(next_trans-1)\n");
            }
          #endif
        }
        else
        {
          /// bound-bound process occurs
          #ifdef DEBUG_ON
            if (debuglevel == 2) printout("[debug] get_event:         tau_rnd - tau <= tau_cont + tau_line: bb-process occurs\n");
          #endif
          edist = dist+ldist;
          if (edist > abort_dist) dummypkt_ptr->next_trans -= 1;
          *rpkt_eventtype = RPKT_EVENTTYPE_BB;
          /// the line and its parameters were already selected by closest_transition!
          endloop = 1;
          #ifdef DEBUG_ON
            if (debuglevel == 2) printout("[debug] get_event:         edist %g, abort_dist %g, edist-abort_dist %g, endloop   %d\n",edist,abort_dist,edist-abort_dist,endloop);
          #endif
        }
      }
      else
      {
        /// continuum process occurs
        edist = dist + (tau_rnd-tau)/kap_cont;
        dummypkt_ptr->next_trans -= 1;
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] get_event:        distance to the occuring continuum event %g, abort_dist %g\n",edist, abort_dist);
        #endif
        *rpkt_eventtype = RPKT_EVENTTYPE_CONT;
        endloop = 1;
      }
    }
    else
    {
      /// no line interaction possible - check whether continuum process occurs in cell
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] get_event:     line interaction impossible\n");
      #endif
      //calculate_kappa_rpkt_cont(dummypkt_ptr, t_current);
      ///no need to restore values set by closest_transition, as nothing was set in this case
      double kap_cont = kappa_rpkt_cont[tid].total;
      double tau_cont = kap_cont*(abort_dist-dist);
      //printout("nu_cmf %g, opticaldepths in ff %g, es %g\n",pkt_ptr->nu_cmf,kappa_rpkt_cont[tid].ff*(abort_dist-dist),kappa_rpkt_cont[tid].es*(abort_dist-dist));
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] get_event:     tau_rnd %g, tau %g, tau_cont %g\n",tau_rnd,tau,tau_cont);
      #endif

      if (tau_rnd - tau > tau_cont)
      {
        /// travel out of cell or time step
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] get_event:       travel out of cell or time step\n");
        #endif
        edist = abort_dist+1e20;
        endloop = 1;
      }
      else
      {
        /// continuum process occurs at edist
        edist = dist + (tau_rnd-tau)/kap_cont;
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] get_event:       continuum process occurs at edist %g\n",edist);
        #endif
        *rpkt_eventtype = RPKT_EVENTTYPE_CONT;
        endloop = 1;
      }
    }
    //propagationcounter += 1;
  }

  pkt_ptr->next_trans = dummypkt_ptr->next_trans;
  #ifdef DEBUG_ON
    if (edist >= 0 || edist < 0)
      return edist;
    else
    {
      printout("edist NaN %g... abort\n",edist);
      abort();
    }
  #else
    return edist;
  #endif
}





///****************************************************************************
int rpkt_event(PKT *pkt_ptr, int rpkt_eventtype, double t_current) //, double kappa_cont, double sigma, double kappa_ff, double kappa_bf)
{
  double zrand;
  //double nnionlevel,nnlevel,nne;
  double nu_edge;
  //double ma_prob,p_maactivate,p_bf,prob;
  //double departure_ratio,corr_photoion;
  //double sigma_bf;
  double kappa_bf_sum;

  //calculate_kappa_rpkt_cont(pkt_ptr, t_current);

  int modelgridindex = cell[pkt_ptr->where].modelgridindex;

  //double nne = get_nne(modelgridindex);
  //double T_e = get_Te(modelgridindex);
  double nu = pkt_ptr->nu_cmf;

  double kappa_cont = kappa_rpkt_cont[tid].total;;
  double sigma = kappa_rpkt_cont[tid].es;
  double kappa_ff = kappa_rpkt_cont[tid].ff;
  double kappa_bf = kappa_rpkt_cont[tid].bf;

  int element,ion,upper,level;
  int i;

  if (rpkt_eventtype == RPKT_EVENTTYPE_BB)
  {
    /// bound-bound transition occured
    /// activate macro-atom in corresponding upper-level. Actually all the information
    /// about the macro atoms state has already been set by closest_transition, so
    /// we need here just the activation!
    #ifdef DEBUG_ON
      if (debuglevel == 2) printout("[debug] rpkt_event: bound-bound activation of macroatom\n");
      //if (tid == 0) ma_stat_activation_bb += 1;
      ma_stat_activation_bb += 1;
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 1;
    #endif
    pkt_ptr->absorptiontype = mastate[tid].activatingline;
    pkt_ptr->absorptionfreq = pkt_ptr->nu_rf;//pkt_ptr->nu_cmf;
    pkt_ptr->absorptiondir[0] = pkt_ptr->dir[0];//pkt_ptr->nu_cmf;
    pkt_ptr->absorptiondir[1] = pkt_ptr->dir[1];//pkt_ptr->nu_cmf;
    pkt_ptr->absorptiondir[2] = pkt_ptr->dir[2];//pkt_ptr->nu_cmf;
    pkt_ptr->type = TYPE_MA;
    #ifndef FORCE_LTE
      //maabs[pkt_ptr->where] += pkt_ptr->e_cmf;
    #endif
    #ifdef RECORD_LINESTAT
      if (tid == 0) acounter[pkt_ptr->next_trans-1] += 1;  /// This way we will only record line statistics from OMP-thread 0
                                                           /// With an atomic pragma or a thread-private structure with subsequent
                                                           /// reduction this could be extended to all threads. However, I'm not
                                                           /// sure if this is worth the additional computational expenses.
    #endif
    //mastate[tid].element = pkt_ptr->nextrans_element;   //store all these nextrans data to MA to save memory!!!!
    //mastate[tid].ion     = pkt_ptr->nextrans_ion;       //MA info becomes important just after activating!
    //mastate[tid].level   = pkt_ptr->nextrans_uppper;
  }
  else
  {
    /// else: continuum process happens. select due to its probabilities sigma/kappa_cont, kappa_ff/kappa_cont, kappa_bf/kappa_cont
    zrand = gsl_rng_uniform(rng);
    #ifdef DEBUG_ON
      if (debuglevel == 2) printout("[debug] rpkt_event: r-pkt undergoes a continuum transition\n");
      if (debuglevel == 2) printout("[debug] rpkt_event:   zrand*kappa_cont %g, sigma %g, kappa_ff %g, kappa_bf %g\n", zrand*kappa_cont,sigma,kappa_ff,kappa_bf);
    #endif
    if (zrand*kappa_cont < sigma)
    {
      /// electron scattering occurs
      /// in this case the packet stays a R_PKT of same nu_cmf than before (coherent scattering)
      /// but with different direction
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] rpkt_event:   electron scattering\n");
        pkt_ptr->interactions += 1;
        pkt_ptr->nscatterings += 1;
        pkt_ptr->last_event = 12;
        escounter += 1;
      #endif

      //pkt_ptr->nu_cmf = 3.7474058e+14;
      escat_rpkt(pkt_ptr,t_current);
      /// Electron scattering does not modify the last emission flag
      //pkt_ptr->emissiontype = get_continuumindex(element,ion-1,lower);
      /// but it updates the last emission position
      pkt_ptr->em_pos[0] = pkt_ptr->pos[0];
      pkt_ptr->em_pos[1] = pkt_ptr->pos[1];
      pkt_ptr->em_pos[2] = pkt_ptr->pos[2];
      pkt_ptr->em_time = t_current;

      /// Set some flags
      //pkt_ptr->next_trans = 0;   ///packet's comoving frame frequency is conserved during electron scattering
                                   ///don't touch the value of next_trans to save transition history
    }
    else if (zrand*kappa_cont < sigma+kappa_ff)
    {
      /// ff: transform to k-pkt
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] rpkt_event:   free-free transition\n");
        //if (tid == 0) k_stat_from_ff += 1;
        k_stat_from_ff += 1;
        pkt_ptr->interactions += 1;
        pkt_ptr->last_event = 5;
      #endif
      pkt_ptr->type = TYPE_KPKT;
      pkt_ptr->absorptiontype = -1;
      #ifndef FORCE_LTE
        //kffabs[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif
    }
    else
    {
      /// bf: transform to k-pkt or activate macroatom correspondig to probabilities
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] rpkt_event:   bound-free transition\n");
      #endif
      pkt_ptr->absorptiontype = -2;

      /// Update the bf-opacity for the packets current frequency
      //calculate_kappa_rpkt_cont(pkt_ptr, t_current);
      double kappa_bf_inrest = kappa_rpkt_cont[tid].bf_inrest;

      /// Determine in which continuum the bf-absorption occurs
      zrand = gsl_rng_uniform(rng);
      kappa_bf_sum = 0.;
      for (i = 0; i < nbfcontinua; i++)
      {
        kappa_bf_sum += phixslist[tid].allcont[i].kappa_bf_contr;
        if (kappa_bf_sum > zrand*kappa_bf_inrest)
        {
          nu_edge = phixslist[tid].allcont[i].nu_edge;
          //if (nu < nu_edge) printout("does this ever happen?\n");
          element = phixslist[tid].allcont[i].element;
          ion = phixslist[tid].allcont[i].ion;
          level = phixslist[tid].allcont[i].level;
          upper = 0;

          #ifdef DEBUG_ON
            if (debuglevel == 2) printout("[debug] rpkt_event:   bound-free: element %d, ion+1 %d, upper %d, ion %d, lower %d\n",element,ion+1,0,ion,level);
            if (debuglevel == 2) printout("[debug] rpkt_event:   bound-free: nu_edge %g, nu %g\n",nu_edge,nu);
          #endif

          /// and decide whether we go to ionisation energy
          zrand = gsl_rng_uniform(rng);
          if (zrand < nu_edge/nu)
          {
            #ifdef DEBUG_ON
              //if (tid == 0) ma_stat_activation_bf += 1;
              ma_stat_activation_bf += 1;
              pkt_ptr->interactions += 1;
              pkt_ptr->last_event = 3;
            #endif
            pkt_ptr->type = TYPE_MA;
            #ifndef FORCE_LTE
              //maabs[pkt_ptr->where] += pkt_ptr->e_cmf;
            #endif
            mastate[tid].element = element;
            mastate[tid].ion     = ion+1;
            mastate[tid].level   = upper;
            mastate[tid].nnlevel = get_levelpop(modelgridindex,element,ion+1,upper);
            mastate[tid].activatingline = -99;
            //if (element == 6) cell[pkt_ptr->where].photoion[ion] += pkt_ptr->e_cmf/pkt_ptr->nu_cmf/H;
          }
          /// or to the thermal pool
          else
          {
            /// transform to k-pkt
            #ifdef DEBUG_ON
              if (debuglevel == 2) printout("[debug] rpkt_event:   bound-free: transform to k-pkt\n");
              //if (tid == 0) k_stat_from_bf += 1;
              k_stat_from_bf += 1;
              pkt_ptr->interactions += 1;
              pkt_ptr->last_event = 4;
            #endif
            pkt_ptr->type = TYPE_KPKT;
            #ifndef FORCE_LTE
              //kbfabs[pkt_ptr->where] += pkt_ptr->e_cmf;
            #endif
            //if (element == 6) cell[pkt_ptr->where].bfabs[ion] += pkt_ptr->e_cmf/pkt_ptr->nu_cmf/H;
          }
          break;
        }
      }

      #ifdef DEBUG_ON
        if (i >= nbfcontinua) printout("[warning] rpkt_event: problem in selecting bf-continuum zrand %g, kappa_bf_sum %g, kappa_bf_inrest %g\n",zrand,kappa_bf_sum,kappa_bf_inrest);
      #endif
    }
  }

  return 0;
}



///****************************************************************************
double closest_transition(PKT *pkt_ptr)
/// for the propagation through non empty cells
{
  double nu_trans;
  int match;

  //int left = 0;
  int left = pkt_ptr->next_trans;
  //printout("[debug] closest_transition: initial left %d\n",left);
  int right = nlines-1;
  int middle = 1;

  //printout("[debug] ___closest_transition___: initial left %d, right %d, nu_cmf %g\n",left,right,pkt_ptr->nu_cmf);
  //printout("[debug] ___closest_transition___: nu_left %g, nu_right%g\n",linelist[left].nu,linelist[right].nu);
  /// if nu_cmf is smaller than the lowest frequency in the linelist,
  /// no line interaction is possible: return negative value as a flag
  if (pkt_ptr->nu_cmf < linelist[right].nu)
  {
    pkt_ptr->next_trans   = nlines+1;  ///helper variable to overcome numerical problems after line scattering
    return -1;
  }
  if (left > right)
  {
    //printout("[debug] pp should have no line interaction anymore\n");
    pkt_ptr->next_trans   = nlines+1;  ///helper variable to overcome numerical problems after line scattering
    return -1;
  }

  if (left > 0)
  {
    /// if left = pkt_ptr->next_trans > 0 we know the next line we should interact with, independent of the packets
    /// current nu_cmf which might be smaller than linelist[left].nu due to propagation errors
    match = left;
  }
  else if (pkt_ptr->nu_cmf >= linelist[0].nu)
  {
    /// if nu_cmf is larger than the highest frequency in the the linelist,
    /// interaction with the first line occurs - no search
    match = 0;
  }
  else
  {
    /// otherwise go through the list until nu_cmf is located between two
    /// entries in the line list and get the index of the closest line
    /// to lower frequencies
    /*
    if (left > 0)
    {
      if ((linelist[left].nu-pkt_ptr->nu_cmf)/pkt_ptr->nu_cmf > 1e-8)
      {
        printout("[warning] closest_transition: search was executed for initial left > 0\n");
        printout("[warning] closest_transition:   left %d, right %d\n",left,right);
        //printf("[warning] closest_transition:   pkt_ptr->nu_cmf %25.16e, linelist[left].nu %25.16e,\n",pkt_ptr->nu_cmf, linelist[left].nu);
        printout("[warning] closest_transition:   pkt_ptr->nu_cmf %g, linelist[left].nu %g,\n",pkt_ptr->nu_cmf, linelist[left].nu);
        printout("[warning] closest_transition:   (linelist[left].nu-pkt_ptr->nu_cmf)/pkt_ptr->nu_cmf %g,\n",(linelist[left].nu-pkt_ptr->nu_cmf)/pkt_ptr->nu_cmf);
        printout("[warning] closest_transition:   linelist[left-1].nu %g, linelist[left+1].nu %g,\n",linelist[left-1].nu, linelist[left+1].nu);
        printout("[warning] closest_transition:   linelist[right].nu %g,\n",linelist[right].nu);
        //printout("[warning] closest_transition:   propagationcounter %d,\n",propagationcounter);
        printout("[warning] closest_transition:   pkt_ptr->last_event %d,\n",pkt_ptr->last_event);
        printout("[warning] closest_transition:   pkt_ptr->number %d,\n",pkt_ptr->number);
        if (pkt_ptr->nu_cmf >= linelist[left+1].nu) printout("[warning] closest_transition:   nu_cmf larger than next line\n");
        else printout("[warning] closest_transition:   jump over several lines\n");
      }
    }
    */
    while (left <= right)  ///must be a "<=" to obtain proper search results!!!
                          ///access to negative array indices is prevented by the upper check
    {
      middle = left + ((right-left) / 2);

      //printout("[debug] middle %d, left %d, right %d, nlines %d\n",middle,left,right,nlines);
      //printout("[debug] linelist[middle].nu %g, linelist[middle-1].nu %g\n",linelist[middle].nu,linelist[middle-1].nu);
      if (pkt_ptr->nu_cmf >= linelist[middle].nu && pkt_ptr->nu_cmf < linelist[middle-1].nu) break;

      if (pkt_ptr->nu_cmf >= linelist[middle].nu)
        right = middle - 1;
      else
        left = middle + 1;
    }
    match = middle;
  }

  /// read transition data out of the linelist and store it as the
  /// next transition for this packet. To save memory it is stored
  /// to the macro atoms state variables. This has no influence until
  /// the macro atom becomes activated by rpkt_event.
  nu_trans = linelist[match].nu;
  mastate[tid].element = linelist[match].elementindex;
  mastate[tid].ion     = linelist[match].ionindex;
  mastate[tid].level   = linelist[match].upperlevelindex;  ///if the MA will be activated it must be in the transitions upper level
  mastate[tid].activatedfromlevel = linelist[match].lowerlevelindex;
  mastate[tid].activatingline = match;
  pkt_ptr->next_trans   = match+1;  ///helper variable to overcome numerical problems after line scattering
                                     ///further scattering events should be located at lower frequencies to prevent
                                     ///multiple scattering events of one pp in a single line
  //printout("[debug] closest_transition: next_trans %d\n",pkt_ptr->next_trans);
  //printout("[debug] closest_transition: linelistindex %d corresponds to transition from level %d to level %d of ion %d of element %d\n",match,mastate[tid].level,pkt_ptr->nextrans_lower,mastate[tid].ion,mastate[tid].element);
  //printout("[debug] closest_transition:   frequency of this transiton %g\n",nu_trans);

  /// return the transitions frequency
  return nu_trans;
}


///****************************************************************************
void emitt_rpkt(PKT *pkt_ptr, double t_current)
{
  /// now make the packet a r-pkt and set further flags
  pkt_ptr->type = TYPE_RPKT;
  pkt_ptr->last_cross = NONE;  /// allow all further cell crossings

  /// Need to assign a new direction. Assume isotropic emission in the cmf
  double zrand = gsl_rng_uniform(rng);
  double zrand2 = gsl_rng_uniform(rng);

  double mu = -1 + (2.*zrand);
  double phi = zrand2 * 2 * PI;
  double sintheta = sqrt(1. - (mu * mu));

  pkt_ptr->dir[0] = sintheta * cos(phi);
  pkt_ptr->dir[1] = sintheta * sin(phi);
  pkt_ptr->dir[2] = mu;
  //printout("[debug] pkt_ptr->dir in CMF: %g %g %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);

  double dummy_dir[3], vel_vec[3];
  /// This direction is in the cmf - we want to convert it to the rest
  /// frame - use aberation of angles. We want to convert from cmf to
  /// rest so need -ve velocity.
  get_velocity(pkt_ptr->pos, vel_vec, (-1*(t_current)));
  ///negative time since we want the backwards transformation here

  angle_ab(pkt_ptr->dir, vel_vec, dummy_dir);
  pkt_ptr->dir[0] = dummy_dir[0];
  pkt_ptr->dir[1] = dummy_dir[1];
  pkt_ptr->dir[2] = dummy_dir[2];
  //printout("[debug] pkt_ptr->dir in RF: %g %g %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);

  /// Check unit vector.
  #ifdef DEBUG_ON
    if (fabs(vec_len(pkt_ptr->dir) - 1) > 1.e-8)
    {
      printout("[fatal] do_ma: Not a unit vector. Abort.\n");
      abort();
    }
  #endif

  /// Finally we want to put in the rest frame energy and frequency. And record
  /// that it's now a r-pkt.
  get_velocity(pkt_ptr->pos, vel_vec, t_current);
  pkt_ptr->nu_rf = pkt_ptr->nu_cmf / doppler(pkt_ptr->dir, vel_vec);

  #ifdef DEBUG_ON
    if (pkt_ptr->e_cmf >1e50)
    {
      printout("[fatal] emitt_rpkt: here %g\n", pkt_ptr->e_cmf);
      abort();
    }
  #endif
  pkt_ptr->e_rf = pkt_ptr->e_cmf * pkt_ptr->nu_rf /pkt_ptr->nu_cmf;

  /**Rest polarization information */
  pkt_ptr->stokes_qu[0]=pkt_ptr->stokes_qu[1]=0.0;
  dummy_dir[0]=dummy_dir[1]=0.0;
  dummy_dir[2]=1.0;
  cross_prod(pkt_ptr->dir,dummy_dir,pkt_ptr->pol_dir);
  if ((dot(pkt_ptr->pol_dir,pkt_ptr->pol_dir)) < 1.e-8)
  {
    dummy_dir[0]=dummy_dir[2]=0.0;
    dummy_dir[1]=1.0;
    cross_prod(pkt_ptr->dir,dummy_dir,pkt_ptr->pol_dir);

  }

  vec_norm(pkt_ptr->pol_dir, pkt_ptr->pol_dir);
  //printout("initialise pol state of packet %g, %g, %g, %g, %g\n",pkt_ptr->stokes_qu[0],pkt_ptr->stokes_qu[1],pkt_ptr->pol_dir[0],pkt_ptr->pol_dir[1],pkt_ptr->pol_dir[2]);
  //printout("pkt direction %g, %g, %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
}


///****************************************************************************
double closest_transition_empty(PKT *pkt_ptr)
/// for the propagation through empty cells
/// here its possible that the packet jumps over several lines
{
  int match;
  //int left = 0;
  int left = pkt_ptr->next_trans;
  //printout("[debug] closest_transition: initial left %d\n",left);
  int right = nlines-1;
  int middle = 1;

  //printout("[debug] ___closest_transition___: initial left %d, right %d, nu_cmf %g\n",left,right,pkt_ptr->nu_cmf);
  //printout("[debug] ___closest_transition___: nu_left %g, nu_right%g\n",linelist[left].nu,linelist[right].nu);
  /// if nu_cmf is smaller than the lowest frequency in the linelist,
  /// no line interaction is possible: return negative value as a flag
  if (pkt_ptr->nu_cmf < linelist[right].nu)
  {
    pkt_ptr->next_trans   = nlines+1;  ///helper variable to overcome numerical problems after line scattering
    return -1;
  }
  if (left > right)
  {
    //printout("[debug] pp should have no line interaction anymore\n");
    pkt_ptr->next_trans   = nlines+1;  ///helper variable to overcome numerical problems after line scattering
    return -1;
  }

  /// no check for left > 0 in the empty case as it is possible that the packet is moved over
  /// several lines through the empty cell
  if (pkt_ptr->nu_cmf >= linelist[left].nu)
  {
    /// if nu_cmf is larger than the highest frequency in the allowed part of the linelist,
    /// interaction with the first line of this part of the list occurs
    match = left;
  }
  else
  {
    /// otherwise go through the list until nu_cmf is located between two
    /// entries in the line list and get the index of the closest line
    /// to lower frequencies

    while (left <= right)  ///must be a "<=" to obtain proper search results!!!
                          ///access to negative array indices is prevented by the upper check
    {
      middle = left + ((right-left) / 2);

      //printout("[debug] middle %d, left %d, right %d, nlines %d\n",middle,left,right,nlines);
      //printout("[debug] linelist[middle].nu %g, linelist[middle-1].nu %g\n",linelist[middle].nu,linelist[middle-1].nu);
      if (pkt_ptr->nu_cmf >= linelist[middle].nu && pkt_ptr->nu_cmf < linelist[middle-1].nu) break;

      if (pkt_ptr->nu_cmf >= linelist[middle].nu)
        right = middle - 1;
      else
        left = middle + 1;
    }
    match = middle;
  }

  /// read transition data out of the linelist and store it as the
  /// next transition for this packet. To save memory it is stored
  /// to the macro atoms state variables. This has no influence until
  /// the macro atom becomes activated by rpkt_event.
  double nu_trans = linelist[match].nu;
  //mastate[tid].element = linelist[match].elementindex;
  //mastate[tid].ion     = linelist[match].ionindex;
  //mastate[tid].level   = linelist[match].upperlevelindex;  ///if the MA will be activated it must be in the transitions upper level
  //mastate[tid].activatedfromlevel   = linelist[match].lowerlevelindex;  ///helper variable for the transitions lower level

  /// For the empty case it's match not match+1: a line interaction is only possible in the next iteration
  /// of the propagation loop. We just have to make sure that the next "normal" line search knows about the
  /// current position of the photon in the frequency list.
  pkt_ptr->next_trans   = match;  /// helper variable to overcome numerical problems after line scattering
                                     ///further scattering events should be located at lower frequencies to prevent
                                     ///multiple scattering events of one pp in a single line

  /// return the transitions frequency
  return nu_trans;
}


///****************************************************************************
void calculate_kappa_rpkt_cont(const PKT *pkt_ptr, double t_current)
{
  double kappa_ffheating = 0.;
  double kappa_ff = 0.;
  double kappa_bf = 0.;

  int modelgridindex = cell[pkt_ptr->where].modelgridindex;

  double sigma = 0.0;
  if (do_r_lc == 1)
  {
    if (opacity_case == 4)
    {
      double nne = get_nne(modelgridindex);
      double T_e = get_Te(modelgridindex);
      //double T_R = get_TR(modelgridindex);
      double nu = pkt_ptr->nu_cmf;
      double g_ff = 1;
      //double g_bf = 1;

      /// First contribution: Thomson scattering on free electrons
      sigma = SIGMA_T * nne;
      //reduced e/s for debugging
      //sigma = 1e-30*sigma;
      //switched off e/s for debugging
      //sigma_cmf = 0. * nne;
      //sigma *= 0.1;

      /// Second contribution: free-free absorption
      kappa_ff = 0.;
      /// Estimator for bound-free heating
      kappa_ffheating = 0.;
      //kappa_bfheating = 0.;
      double ionpops_local[MELEMENTS][MIONS];
      for (int element = 0; element < nelements; element++)
      {
        int nions = get_nions(element);
        for (int ion = 0; ion < nions; ion++)
        {
          ///calculate population of ionstage ...
          double nnion = ionstagepop(modelgridindex,element,ion); ///partfunct needs to be adjusted
          ionpops_local[element][ion] = nnion / get_nnetot(modelgridindex);
          //Z = get_element(element);  ///atomic number
          //if (get_ionstage(element,ion) > 1)
          /// Z is ionic charge in the following formula
          int Z = get_ionstage(element,ion)-1;
          if (Z > 0)
          {
            //kappa_ff += 3.69255e8 * pow(Z,2) / sqrt(T_e) * pow(nu,-3) * g_ff * nne * nnion * (1-exp(-HOVERKB*nu/T_e));
            kappa_ff += pow(Z,2) * g_ff * nnion;
            //kappa_ffheating += pow(Z,2) * g_ff * nnion;
            /// heating with level dependence
            //kappa_ffheating += 3.69255e8 * pow(Z,2) / sqrt(T_e) * pow(nu,-3) * g_ff * nne * nnion * (1-exp(-HOVERKB*nu/T_e));
            /// heating without level dependence
            //kappa_ffheating += 3.69255e8 * pow(Z,2) * pow(nu,-3) * g_ff * (1-exp(-HOVERKB*nu/T_e));
          }
        }
      }
      kappa_ff *= 3.69255e8 / sqrt(T_e) * pow(nu,-3) * nne * (1 - exp(-HOVERKB*nu/T_e));
      //kappa_ffheating *= 3.69255e8 / sqrt(T_e) * pow(nu,-3) * nne * (1-exp(-HOVERKB*nu/T_e));
      kappa_ffheating = kappa_ff;
      //kappa_ff *= 1e5;

      /// Third contribution: bound-free absorption
      kappa_bf = 0.;
      for (int i = 0; i < nbfcontinua; i++)
      {
        int element = phixslist[tid].allcont[i].element;
        int ion = phixslist[tid].allcont[i].ion;
        int level = phixslist[tid].allcont[i].level;
        /// The bf process happens only if the current cell contains
        /// the involved atomic species
        if ((get_abundance(modelgridindex,element) > 0) && ((ionpops_local[element][ion] > 1.e-6) || (level == 0)))
        {
          double nu_edge = phixslist[tid].allcont[i].nu_edge;
          //printout("i %d, nu_edge %g\n",i,nu_edge);
          if (nu >= nu_edge)
          {
            //ion = phixslist[tid].allcont[i].ion;
            //level = phixslist[tid].allcont[i].level;
            //nnlevel = samplegrid[samplecell].phixslist[i].nnlevel;
            //printout("element %d, ion %d, level %d, nnlevel %g\n",element,ion,level,nnlevel);

            double nnlevel = get_levelpop(modelgridindex,element,ion,level);
            //if (fabs(nnlevel - phixslist[tid].allcont[i].nnlevel) > 0)
            //{
            //  printout("history value %g, phixslist value %g\n",nnlevel,phixslist[tid].allcont[i].nnlevel);
            //  printout("pkt_ptr->number %d\n",pkt_ptr->number);
            //}
            mastate[tid].element = element;
            mastate[tid].ion = ion;
            mastate[tid].level = level;
            double sigma_bf = photoionization_crosssection(nu_edge,nu);
            double check = 0.0; //corrfactor
            for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
            {
              int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
              double nnionlevel = calculate_exclevelpop(modelgridindex,element,ion+1,upper);
              double sf = get_sahafact(element,ion,level,phixstargetindex,T_e,nu_edge*H);
              double helper = nnlevel * sigma_bf * get_phixsprobability(element,ion,level,phixstargetindex);
              double departure_ratio = nnionlevel / nnlevel * nne * sf; // put that to phixslist
              if (nnlevel == 0.)
                check += 0.;
              else
                check += helper * (1 - departure_ratio * exp(-HOVERKB*nu/T_e));
              if ((level == 0) && (phixstargetindex == 0))
              {
                int gphixsindex = phixslist[tid].allcont[i].index_in_groundphixslist;
                double corrfactor = 1 - departure_ratio * exp(-HOVERKB*nu/T_e);
                if (corrfactor < 0)
                  corrfactor = 1;
                phixslist[tid].groundcont[gphixsindex].gamma_contr = sigma_bf * corrfactor * get_phixsprobability(element,ion,level,phixstargetindex);
              }
            }

            if (check <= 0)
            {
              #ifdef DEBUG_ON
                //printout("[warning] calculate_kappa_rpkt_cont: kappa_bf has negative contribution T_e %g, T_R %g, W %g, E_threshold %g, check %g\n", T_e,cell[pkt_ptr->where].T_R,cell[pkt_ptr->where].W,nu_edge*H,nnlevel*sigma_bf);
                //printout("[warning] calculate_kappa_rpkt_cont: set this contribution to zero\n");
                if (fabs(check) > 1e-10) printout("[warning] calculate_kappa_rpkt_cont: kappa_bf has negative contribution %g for element %d ion %d level %d (nnlevel %g, nne %g, sigma_bf %g)\n", check,element,ion,level,nnlevel,nne,sigma_bf);
              #endif
              check = 0.;
              //phixslist[tid].allcont[i].kappa_bf_contr = check;
              //phixslist[tid].allcont[i].photoion_contr = 0.;
              //phixslist[tid].allcont[i].stimrecomb_contr = 0.;
            }
/*
            else
            {
              phixslist[tid].allcont[i].kappa_bf_contr = check;
              phixslist[tid].allcont[i].photoion_contr = nnlevel * sigma_bf;
              phixslist[tid].allcont[i].stimrecomb_contr = sf * sigma_bf;
            }*/
            //check *= 2;
            phixslist[tid].allcont[i].kappa_bf_contr = check;
            kappa_bf += check;
            #ifdef DEBUG_ON
              if (!isfinite(check))
              {
                printout("[fatal] calculate_kappa_rpkt_cont: non-finite contribution to kappa_bf %g ... abort\n",check);
                printout("[fatal] phixslist index %d, element %d, ion %d, level %d\n",i,element,ion,level);
                printout("[fatal] cell[%d].composition[%d].abundance = %g\n",modelgridindex,element,get_abundance(modelgridindex,element));
                printout("[fatal] nne %g, nnlevel %g, \n",nne,nnlevel);
                printout("[fatal] sigma_bf %g, T_e %g, nu %g, nu_edge %g\n",sigma_bf,T_e,nu,nu_edge);
                abort();
              }
            #endif
          }
          /// The important part of phixslist is sorted by nu_edge in ascending order
          /// If nu < phixslist[tid].allcont[i].nu_edge no absorption in any of the following continua
          /// is possible, therefore leave the loop.
          else
            break;
          //  {
          //  /// Set photoion_contr to zero for continua with nu < nu_edge
          //  /// to get the correct estimators for the photoionisation rate coefficients
          //  for (ii = i; ii < importantbfcontinua; ii++) phixslist[tid].allcont[ii].photoion_contr = 0;
          //  break;
          //}
        }
        else
        {
          phixslist[tid].allcont[i].kappa_bf_contr = 0.;
          if (phixslist[tid].allcont[i].level == 0)
          {
            int gphixsindex = phixslist[tid].allcont[i].index_in_groundphixslist;
            //phixslist[tid].groundcont[gphixsindex].photoion_contr = 0.;
            phixslist[tid].groundcont[gphixsindex].gamma_contr = 0.;
            //phixslist[tid].groundcont[gphixsindex].stimrecomb_contr = 0.;
            //phixslist[tid].groundcont[gphixsindex].bfheating_contr = 0.;
          }
        }
      }
      kappa_rpkt_cont[tid].bf_inrest = kappa_bf;
    }
    else
    {
      double nne = get_nne(modelgridindex);
      double T_e = get_Te(modelgridindex);
      //double T_R = get_TR(modelgridindex);
      double g_ff = 1;
      double nu = pkt_ptr->nu_cmf;
      /// in the other cases kappa_grey is an mass absorption coefficient
      /// therefore use the mass density
      //sigma = cell[pkt_ptr->where].kappa_grey * cell[pkt_ptr->where].rho;
      sigma = SIGMA_T * nne;

      sigma = 0.;
      /*
      kappa_ff = 0.9*sigma;
      sigma *= 0.1;
      kappa_bf = 0.;
      */

      /// Second contribution: free-free absorption
      kappa_ff = 0.;
      for (int element = 0; element < nelements; element++)
      {
        int nions = get_nions(element);
        for (int ion = 0; ion < nions; ion++)
        {
          ///calculate population of ionstage ...
          double nnion = ionstagepop(modelgridindex,element,ion); ///partfunct needs to be adjusted
          /// Z is ionic charge in the following formula
          int Z = get_ionstage(element,ion)-1;
          if (Z > 0)
          {
            kappa_ff += pow(Z,2) * g_ff * nnion;
          }
        }
      }
      kappa_ff *= 1e5 * 3.69255e8 / sqrt(T_e) * pow(nu,-3) * nne * (1-exp(-HOVERKB*nu/T_e));
      kappa_bf = 0.;
    }

    /// Now need to convert between frames.
    double vel_vec[3];
    get_velocity(pkt_ptr->pos, vel_vec, t_current);
    double dopplerfac = doppler(pkt_ptr->dir, vel_vec);
    sigma = sigma * dopplerfac;
    kappa_ff = kappa_ff * dopplerfac;
    kappa_bf = kappa_bf * dopplerfac;
  }
  else
  {
    sigma = 0.;
    kappa_ff = 0.;
    kappa_bf = 0.;
    //kappa_ffheating = 0.;
    //kappa_bfheating = 0.;
  }

  kappa_rpkt_cont[tid].total = sigma + kappa_bf + kappa_ff;
  #ifdef DEBUG_ON
    if (debuglevel == 2) printout("[debug]  ____kappa_rpkt____: kappa_cont %g, sigma %g, kappa_ff %g, kappa_bf %g\n",kappa_rpkt_cont[tid].total,sigma,kappa_ff,kappa_bf);
  #endif
  kappa_rpkt_cont[tid].es = sigma;
  kappa_rpkt_cont[tid].ff = kappa_ff;
  kappa_rpkt_cont[tid].bf = kappa_bf;
  kappa_rpkt_cont[tid].ffheating = kappa_ffheating;
  //kappa_rpkt_cont[tid].bfheating = kappa_bfheating;


  #ifdef DEBUG_ON
    if (!isfinite(kappa_rpkt_cont[tid].total))
    {
      printout("[fatal] calculate_kappa_rpkt_cont: resulted in non-finite kappa_rpkt_cont.total ... abort\n");
      printout("[fatal] es %g, ff %g, bf %g\n",kappa_rpkt_cont[tid].es,kappa_rpkt_cont[tid].ff,kappa_rpkt_cont[tid].bf);
      printout("[fatal] nbfcontinua %d\n",nbfcontinua);
      printout("[fatal] in cell %d with density %g\n",modelgridindex,get_rho(modelgridindex));
      printout("[fatal] pkt_ptr->nu_cmf %g\n",pkt_ptr->nu_cmf);
      if (isfinite(kappa_rpkt_cont[tid].es))
      {
        kappa_rpkt_cont[tid].ff = 0.;
        kappa_rpkt_cont[tid].bf = 0.;
        kappa_rpkt_cont[tid].total = kappa_rpkt_cont[tid].es;
      }
      else
      {
        abort();
      }
    }
  #endif

}



///****************************************************************************
void calculate_kappa_vpkt_cont(const PKT *pkt_ptr, double t_current)
{
    double vel_vec[3];
    double sigma;
    double kappa_ffheating = 0.;//,kappa_bfheating;
    double nne,T_e,T_R,nu;
    double g_ff,g_bf;
    double nnion,nnionlevel,nnlevel,departure_ratio;
    int element,ion,level;//,samplecell;
    int Z;
    double nu_edge;
    int i,nions;
    double sf,check;
    double helper;
    int gphixsindex;
    double corrfactor;

    double bef;

    double kappa_ff = 0.;
    double kappa_bf = 0.;

    int modelgridindex = cell[pkt_ptr->where].modelgridindex;

    if (do_r_lc == 1)
    {
        if (opacity_case == 4)
        {
            nne = get_nne(modelgridindex);
            T_e = get_Te(modelgridindex);
            T_R = get_TR(modelgridindex);
            nu = pkt_ptr->nu_cmf;
            g_ff = 1;
            g_bf = 1;

            /// First contribution: Thomson scattering on free electrons
            sigma = SIGMA_T * nne;
            //reduced e/s for debugging
            //sigma = 1e-30*sigma;
            //switched off e/s for debugging
            //sigma_cmf = 0. * nne;
            //sigma *= 0.1;

            /// Second contribution: free-free absorption
            kappa_ff = 0.;
            /// Estimator for bound-free heating
            kappa_ffheating = 0.;
            //kappa_bfheating = 0.;
            double ionpops_local[MELEMENTS][MIONS];
            for (int element = 0; element < nelements; element++)
            {
                int nions = get_nions(element);
                for (int ion = 0; ion < nions; ion++)
                {
                    ///calculate population of ionstage ...
                    nnion = ionstagepop(modelgridindex,element,ion); ///partfunct needs to be adjusted
                    ionpops_local[element][ion] = nnion/get_nnetot(modelgridindex);
                    //Z = get_element(element);  ///atomic number
                    //if (get_ionstage(element,ion) > 1)
                    /// Z is ionic charge in the following formula
                    Z = get_ionstage(element,ion)-1;
                    if (Z > 0)
                    {
                        //kappa_ff += 3.69255e8 * pow(Z,2) / sqrt(T_e) * pow(nu,-3) * g_ff * nne * nnion * (1-exp(-HOVERKB*nu/T_e));
                        kappa_ff += pow(Z,2) * g_ff * nnion;
                        //kappa_ffheating += pow(Z,2) * g_ff * nnion;
                        /// heating with level dependence
                        //kappa_ffheating += 3.69255e8 * pow(Z,2) / sqrt(T_e) * pow(nu,-3) * g_ff * nne * nnion * (1-exp(-HOVERKB*nu/T_e));
                        /// heating without level dependence
                        //kappa_ffheating += 3.69255e8 * pow(Z,2) * pow(nu,-3) * g_ff * (1-exp(-HOVERKB*nu/T_e));
                    }
                }
            }
            kappa_ff *= 3.69255e8 / sqrt(T_e) * pow(nu,-3) * nne * (1-exp(-HOVERKB*nu/T_e));
            //kappa_ffheating *= 3.69255e8 / sqrt(T_e) * pow(nu,-3) * nne * (1-exp(-HOVERKB*nu/T_e));
            kappa_ffheating = kappa_ff;
            //kappa_ff *= 1e5;

            /// Third contribution: bound-free absorption
            kappa_bf = 0.;
            for (i = 0; i < nbfcontinua; i++)
            {
                element = phixslist[tid].allcont[i].element;
                ion = phixslist[tid].allcont[i].ion;
                level = phixslist[tid].allcont[i].level;

                /// The bf process happens only if the current cell contains
                /// the involved atomic species
                if ((ionpops_local[element][ion] > 1.e-6) || (level == 0))
                    ///if (get_abundance(modelgridindex,element) > 0)
                {
                    nu_edge = phixslist[tid].allcont[i].nu_edge;
                    //printout("i %d, nu_edge %g\n",i,nu_edge);
                    if (nu >= nu_edge)
                    {
                        //ion = phixslist[tid].allcont[i].ion;
                        //level = phixslist[tid].allcont[i].level;
                        //nnlevel = samplegrid[samplecell].phixslist[i].nnlevel;
                        //printout("element %d, ion %d, level %d, nnlevel %g\n",element,ion,level,nnlevel);

                        nnlevel = calculate_exclevelpop(modelgridindex,element,ion,level);
                        //if (fabs(nnlevel - phixslist[tid].allcont[i].nnlevel) > 0)
                        //{
                        //  printout("history value %g, phixslist value %g\n",nnlevel,phixslist[tid].allcont[i].nnlevel);
                        //  printout("pkt_ptr->number %d\n",pkt_ptr->number);
                        //}
                        nnionlevel = get_groundlevelpop(modelgridindex,element,ion+1);

                        mastate[tid].element = element;
                        mastate[tid].ion = ion;
                        mastate[tid].level = level;
                        double sigma_bf = photoionization_crosssection(nu_edge,nu);

                        bef = cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[0].sahafact;

                        sf = calculate_sahafact(element,ion,level,0,T_e,nu_edge*H);
                        helper = nnlevel * sigma_bf;
                        departure_ratio = nnionlevel/nnlevel * nne*sf; ///put that to phixslist
                        if (nnlevel == 0.) check = 0.;
                        else check = helper * (1 - departure_ratio * exp(-HOVERKB*nu/T_e));

                        //printout("[warning] calculate_kappa_rpkt_cont: kappa_bf has negative contribution %g for element %d ion %d level %d (nnionlevel %g, nnlevel %g, nne %g, sf %g, sigma_bf %g) before %g after %g \n", check,element,ion,level,nnionlevel,nnlevel,nne,sf,sigma_bf,bef,cellhistory[tid].chelements[element].chions[ion].chlevels[level].sahafact);

                        if (check <= 0)
                        {
                            #ifdef DEBUG_ON
                            //printout("[warning] calculate_kappa_rpkt_cont: kappa_bf has negative contribution T_e %g, T_R %g, W %g, E_threshold %g, check %g\n", T_e,cell[pkt_ptr->where].T_R,cell[pkt_ptr->where].W,nu_edge*H,nnlevel*sigma_bf);
                            //printout("[warning] calculate_kappa_rpkt_cont: set this contribution to zero\n");
                            if (fabs(check) > 1e-10) printout("[warning] calculate_kappa_rpkt_cont: kappa_bf has negative contribution %g for element %d ion %d level %d (nnionlevel %g, nnlevel %g, nne %g, sf %g, sigma_bf %g) before %g after %g \n", check,element,ion,level,nnionlevel,nnlevel,nne,sf,sigma_bf,bef,cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[0].sahafact);
                            #endif
                            check = 0.;
                            //phixslist[tid].allcont[i].kappa_bf_contr = check;
                            //phixslist[tid].allcont[i].photoion_contr = 0.;
                            //phixslist[tid].allcont[i].stimrecomb_contr = 0.;
                        }
                        /*            else
                         {
                         phixslist[tid].allcont[i].kappa_bf_contr = check;
                         phixslist[tid].allcont[i].photoion_contr = nnlevel * sigma_bf;
                         phixslist[tid].allcont[i].stimrecomb_contr = sf * sigma_bf;
                         }*/
                        //check *= 2;
                        phixslist[tid].allcont[i].kappa_bf_contr = check;
                        kappa_bf += check;
                        if (level == 0)
                        {
                            gphixsindex = phixslist[tid].allcont[i].index_in_groundphixslist;
                            //groundphixslist[gphixsindex].photoion_contr = helper;
                            corrfactor = 1 - departure_ratio * exp(-HOVERKB*nu/T_e);
                            if (corrfactor < 0) corrfactor = 1;
                            phixslist[tid].groundcont[gphixsindex].gamma_contr = sigma_bf * corrfactor;
                            //phixslist[tid].groundcont[gphixsindex].stimrecomb_contr = sf * sigma_bf;
                            //phixslist[tid].groundcont[gphixsindex].bfheating_contr = helper * nu_edge;
                        }
                        #ifdef DEBUG_ON
                        if (!isfinite(check))
                        {
                            printout("[fatal] calculate_kappa_rpkt_cont: non-finite contribution to kappa_bf %g ... abort\n",check);
                            printout("[fatal] phixslist index %d, element %d, ion %d, level %d\n",i,element,ion,level);
                            printout("[fatal] cell[%d].composition[%d].abundance = %g\n",modelgridindex,element,get_abundance(modelgridindex,element));
                            printout("[fatal] nne %g, nnlevel %g, nnionlevel %g, departure_ratio %g\n",nne,nnlevel,nnionlevel,departure_ratio);
                            printout("[fatal] sigma_bf %g, T_e %g, nu %g, nu_edge %g\n",sigma_bf,T_e,nu,nu_edge);
                            abort();
                        }
                        #endif
                    }
                    /// The important part of phixslist is sorted by nu_edge in ascending order
                    /// If nu < phixslist[tid].allcont[i].nu_edge no absorption in any of the following continua
                    /// is possible, therefore leave the loop.
                    else break;
                    //  {
                    //  /// Set photoion_contr to zero for continua with nu < nu_edge
                    //  /// to get the correct estimators for the photoionisation rate coefficients
                    //  for (int ii = i; ii < importantbfcontinua; ii++) phixslist[tid].allcont[ii].photoion_contr = 0;
                    //  break;
                    //}
                }
                else
                {
                    phixslist[tid].allcont[i].kappa_bf_contr = 0.;
                    if (phixslist[tid].allcont[i].level == 0)
                    {
                        gphixsindex = phixslist[tid].allcont[i].index_in_groundphixslist;
                        //phixslist[tid].groundcont[gphixsindex].photoion_contr = 0.;
                        phixslist[tid].groundcont[gphixsindex].gamma_contr = 0.;
                        //phixslist[tid].groundcont[gphixsindex].stimrecomb_contr = 0.;
                        //phixslist[tid].groundcont[gphixsindex].bfheating_contr = 0.;
                    }
                }
            }
            kappa_rpkt_cont[tid].bf_inrest = kappa_bf;
        }
        else
        {
            nne = get_nne(modelgridindex);
            T_e = get_Te(modelgridindex);
            T_R = get_TR(modelgridindex);
            g_ff = 1;
            nu = pkt_ptr->nu_cmf;
            /// in the other cases kappa_grey is an mass absorption coefficient
            /// therefore use the mass density
            //sigma = cell[pkt_ptr->where].kappa_grey * cell[pkt_ptr->where].rho;
            sigma = SIGMA_T*nne;

            sigma = 0.;
            /*
             kappa_ff = 0.9*sigma;
             sigma *= 0.1;
             kappa_bf = 0.;
             */

            /// Second contribution: free-free absorption
            kappa_ff = 0.;
            for (element = 0; element < nelements; element++)
            {
                nions = get_nions(element);
                for (ion = 0; ion < nions; ion++)
                {
                    ///calculate population of ionstage ...
                    nnion = ionstagepop(modelgridindex,element,ion); ///partfunct needs to be adjusted
                    /// Z is ionic charge in the following formula
                    Z = get_ionstage(element,ion)-1;
                    if (Z > 0)
                    {
                        kappa_ff += pow(Z,2) * g_ff * nnion;
                    }
                }
            }
            kappa_ff *= 1e5*3.69255e8 / sqrt(T_e) * pow(nu,-3) * nne * (1-exp(-HOVERKB*nu/T_e));
            kappa_bf = 0.;
        }

        /// Now need to convert between frames.
        get_velocity(pkt_ptr->pos, vel_vec, t_current);
        sigma = sigma * doppler(pkt_ptr->dir, vel_vec);
        kappa_ff = kappa_ff * doppler(pkt_ptr->dir, vel_vec);
        kappa_bf = kappa_bf * doppler(pkt_ptr->dir, vel_vec);
    }
    else
    {
        sigma = 0.0;
        kappa_ff = 0.;
        kappa_bf = 0.;
        //kappa_ffheating = 0.;
        //kappa_bfheating = 0.;
    }

    kappa_rpkt_cont[tid].total = sigma + kappa_bf + kappa_ff;
    #ifdef DEBUG_ON
    if (debuglevel == 2) printout("[debug]  ____kappa_rpkt____: kappa_cont %g, sigma %g, kappa_ff %g, kappa_bf %g\n",kappa_rpkt_cont[tid].total,sigma,kappa_ff,kappa_bf);
    #endif
    kappa_rpkt_cont[tid].es = sigma;
    kappa_rpkt_cont[tid].ff = kappa_ff;
    kappa_rpkt_cont[tid].bf = kappa_bf;
    kappa_rpkt_cont[tid].ffheating = kappa_ffheating;
    //kappa_rpkt_cont[tid].bfheating = kappa_bfheating;


    #ifdef DEBUG_ON
    if (!isfinite(kappa_rpkt_cont[tid].total))
    {
        printout("[fatal] calculate_kappa_rpkt_cont: resulted in non-finite kappa_rpkt_cont.total ... abort\n");
        printout("[fatal] es %g, ff %g, bf %g\n",kappa_rpkt_cont[tid].es,kappa_rpkt_cont[tid].ff,kappa_rpkt_cont[tid].bf);
        printout("[fatal] nbfcontinua %d\n",nbfcontinua);
        printout("[fatal] in cell %d with density %g\n",modelgridindex,get_rho(modelgridindex));
        printout("[fatal] pkt_ptr->nu_cmf %g, T_e %g, nne %g\n",pkt_ptr->nu_cmf,T_e,nne);
        if (isfinite(kappa_rpkt_cont[tid].es))
        {
            kappa_rpkt_cont[tid].ff = 0.;
            kappa_rpkt_cont[tid].bf = 0.;
            kappa_rpkt_cont[tid].total = kappa_rpkt_cont[tid].es;
        }
        else
        {
            abort();
        }
    }
    #endif

    return;
}



///****************************************************************************
int compare_phixslistentry_bynuedge(const void *p1, const void *p2)
/// Helper function to sort the phixslist by ascending threshold frequency.
{
  const fullphixslist_t *a1 = (fullphixslist_t *)(p1);
  const fullphixslist_t *a2 = (fullphixslist_t *)(p2);

  double edge_diff = a1->nu_edge - a2->nu_edge;
  if (edge_diff < 0)
    return -1;
  else if (edge_diff > 0)
    return 1;
  else
    return 0;
}



///****************************************************************************
int compare_groundphixslistentry_bynuedge(const void *p1, const void *p2)
/// Helper function to sort the groundphixslist by ascending threshold frequency.
{
  const groundphixslist_t *a1 = (groundphixslist_t *)(p1);
  const groundphixslist_t *a2 = (groundphixslist_t *)(p2);

  double edge_diff = a1->nu_edge - a2->nu_edge;
  if (edge_diff < 0)
    return -1;
  else if (edge_diff > 0)
    return 1;
  else
    return 0;
}


///****************************************************************************
double do_rpkt_thickcell(PKT *pkt_ptr, double t1, double t2)
/** Routine for moving an r-packet. Similar to do_gamma in objective.*/
{
  double vel_vec[3];

  double tdist;
  double edist;
  int snext;

  double kappa, kappa_cont, sigma, kappa_ff, kappa_bf;
  double *kappacont_ptr, *sigma_ptr, *kappaff_ptr, *kappabf_ptr;
  kappacont_ptr = &kappa_cont;
  sigma_ptr = &sigma;
  kappaff_ptr = &kappa_ff;
  kappabf_ptr = &kappa_bf;

  double t_current = t1; ///this will keep track of time in the calculation
  //printout("[debug] r-pkt propagation init\n");
  //int it = 1;

  int end_packet = 0; ///means "keep working"
  while (end_packet == 0)
  {

    /*
    int i;
    for (i=0; i < CELLHISTORYSIZE; i++)
      printout("thread%d _ do_rpkt: cellhistory[%d].cellnumber = %d\n",tid,i,cellhistory[i].cellnumber);
    */

    #ifdef DEBUG_ON
      if (pkt_ptr->next_trans > 0)
      {
        if (debuglevel == 2) printout("[debug] do_rpkt: init: pkt_ptr->nu_cmf %g, nu(pkt_ptr->next_trans=%d) %g, nu(pkt_ptr->next_trans-1=%d) %g, pkt_ptr->where %d\n", pkt_ptr->nu_cmf, pkt_ptr->next_trans, linelist[pkt_ptr->next_trans].nu, pkt_ptr->next_trans-1, linelist[pkt_ptr->next_trans-1].nu, pkt_ptr->where );
        if (debuglevel == 2) printout("[debug] do_rpkt: init: (pkt_ptr->nu_cmf - nu(pkt_ptr->next_trans-1))/pkt_ptr->nu_cmf %g\n", (pkt_ptr->nu_cmf-linelist[pkt_ptr->next_trans-1].nu)/pkt_ptr->nu_cmf);
      }
    #endif

    //printout("[debug] r-pkt propagation iteration %d\n",it);
    //it += 1;
    /** Assign optical depth to next physical event. And start couter of
    optical depth for this path.*/
    double zrand = gsl_rng_uniform(rng);
    double tau_next = -1. * log(zrand);
    double tau_current = 0.0;

    /** Start by finding the distance to the crossing of the grid cell
    boundaries. sdist is the boundary distance and snext is the
    grid cell into which we pass.*/
    double sdist = boundary_cross(pkt_ptr, t_current, &snext);

    if (sdist == 0)
      change_cell(pkt_ptr, snext, &end_packet, t_current);
    else
    {
      if (sdist > (rmax * t_current/tmin))
      {
        printout("[fatal] do_rpkt: Unreasonably large sdist. Rpkt. Abort. %g %g %g\n", rmax, t_current/tmin, sdist);
        abort();
      }

      if (sdist < 1)
      {
        printout("[warning] r_pkt: Negative distance (sdist = %g). Abort.\n", sdist);
        printout("[warning] r_pkt: cell %d snext %d\n", pkt_ptr->where, snext);
        printout("[warning] r_pkt: pos %g %g %g\n", pkt_ptr->pos[0], pkt_ptr->pos[1], pkt_ptr->pos[2]);
        printout("[warning] r_pkt: dir %g %g %g\n", pkt_ptr->dir[0], pkt_ptr->dir[1], pkt_ptr->dir[2]);
        printout("[warning] r_pkt: cell corner %g %g %g\n",cell[pkt_ptr->where].pos_init[0]*t_current/tmin, cell[pkt_ptr->where].pos_init[1]*t_current/tmin,  cell[pkt_ptr->where].pos_init[2]*t_current/tmin);
        printout("[warning] r_pkt: cell width %g %g %g\n",wid_init*t_current/tmin, wid_init*t_current/tmin,  wid_init*t_current/tmin);
        //abort();
      }
      if (((snext != -99) && (snext < 0)) || (snext >= ngrid))
      {
        printout("[fatal] r_pkt: Heading for inappropriate grid cell. Abort.\n");
        printout("[fatal] r_pkt: Current cell %d, target cell %d.\n", pkt_ptr->where, snext);
        abort();
      }
      if (sdist > max_path_step)
      {
        sdist = max_path_step;
        snext = pkt_ptr->where;
      }


      /// Find how far it can travel during the time inverval.
      tdist = (t2 - t_current) * CLIGHT_PROP;
      if (tdist < 0)
      {
        printout("[fatal] do_rpkt: Negative distance (tdist). Abort. \n");
        abort();
      }

      /// Get distance to the next physical event in this case only electron scattering
      //kappa = SIGMA_T*get_nne(cell[pkt_ptr->where].modelgridindex);
      kappa = get_kappagrey(cell[pkt_ptr->where].modelgridindex) * get_rho(cell[pkt_ptr->where].modelgridindex);
      get_velocity(pkt_ptr->pos, vel_vec, t_current);
      kappa = kappa * doppler(pkt_ptr->dir, vel_vec);
      edist = (tau_next - tau_current) / kappa;
      if (edist < 0)
      {
        printout("[fatal] do_rpkt: Negative distance (edist). Abort. \n");
        printout("[fatal] do_rpkt: Trouble was due to packet number %d.\n", pkt_ptr->number);
        abort();
      }

      if ((sdist < tdist) && (sdist < edist))
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_rpkt: sdist < tdist && sdist < edist\n");
        #endif
        /** Move it into the new cell. */
        sdist = sdist / 2.;
        t_current += sdist / CLIGHT_PROP;
        move_pkt(pkt_ptr,sdist,t_current);
        update_estimators(pkt_ptr,sdist*2);
        if (do_rlc_est != 0 && do_rlc_est != 3)
        {
          sdist = sdist * 2.;
          rlc_emiss_rpkt(pkt_ptr, sdist, t_current);
          sdist = sdist / 2.;
        }
        t_current += sdist / CLIGHT_PROP;
        move_pkt(pkt_ptr,sdist,t_current);
        sdist = sdist * 2.;

        if (snext != pkt_ptr->where)
        {
          change_cell(pkt_ptr, snext, &end_packet, t_current);
        }
        #ifdef DEBUG_ON
          /** New cell so reset the scat_counter */
          pkt_ptr->scat_count = 0;
          //if (debuglevel == 2) printout("[debug] do_rpkt:   pkt_ptr->last_event %d\n",pkt_ptr->last_event);
          pkt_ptr->last_event = pkt_ptr->last_event + 100;
        #endif
      }
      else if ((tdist < sdist) && (tdist < edist))
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_rpkt: tdist < sdist && tdist < edist\n");
        #endif
        /** Doesn't reach boundary. */
        tdist = tdist / 2.;
        t_current += tdist / CLIGHT_PROP;
        move_pkt(pkt_ptr,tdist,t_current);
        update_estimators(pkt_ptr,tdist*2);
        if (do_rlc_est != 0 && do_rlc_est != 3)
        {
          tdist = tdist * 2.;
          rlc_emiss_rpkt(pkt_ptr, tdist, t_current);
          tdist = tdist / 2.;
        }
        t_current = t2;
        move_pkt(pkt_ptr,tdist,t_current);
        tdist = tdist * 2.;
        #ifdef DEBUG_ON
          pkt_ptr->last_event = pkt_ptr->last_event + 1000;
        #endif
        end_packet = 1;
      }
      else if ((edist < sdist) && (edist < tdist))
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] do_rpkt: edist < sdist && edist < tdist\n");
        #endif
        edist = edist / 2.;
        t_current += edist / CLIGHT_PROP;
        move_pkt(pkt_ptr,edist,t_current);
        update_estimators(pkt_ptr,edist*2);
        if (do_rlc_est != 0 && do_rlc_est != 3)
        {
          edist = edist * 2.;
          rlc_emiss_rpkt(pkt_ptr, edist, t_current);
          edist = edist / 2.;
        }
        t_current += edist / CLIGHT_PROP;
        move_pkt(pkt_ptr,edist,t_current);
        edist = edist * 2.;

        /// electron scattering occurs
        /// in this case the packet stays a R_PKT of same nu_cmf than before (coherent scattering)
        /// but with different direction
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] rpkt_event:   electron scattering\n");
          pkt_ptr->interactions += 1;
          pkt_ptr->nscatterings += 1;
          pkt_ptr->last_event = 12;
          escounter += 1;
        #endif

        //pkt_ptr->nu_cmf = 3.7474058e+14;
        emitt_rpkt(pkt_ptr,t_current);
        /// Electron scattering does not modify the last emission flag
        //pkt_ptr->emissiontype = get_continuumindex(element,ion-1,lower);
        /// but it updates the last emission position
        pkt_ptr->em_pos[0] = pkt_ptr->pos[0];
        pkt_ptr->em_pos[1] = pkt_ptr->pos[1];
        pkt_ptr->em_pos[2] = pkt_ptr->pos[2];
        pkt_ptr->em_time = t_current;

        /*/// The previously selected and in pkt_ptr stored event occurs. Handling is done by rpkt_event
        rpkt_event(pkt_ptr,rpkt_eventtype,t_current);
        if (pkt_ptr->type != TYPE_RPKT)
        {
          /// It's not an r-packet any more - return.
          return(t_current);
        }*/
      }
      else
      {
        printout("[fatal] do_rpkt: Failed to identify event . Rpkt. edist %g, sdist %g, tdist %g Abort.\n", edist, sdist, tdist);
        printout("[fatal] do_rpkt: Trouble was due to packet number %d.\n", pkt_ptr->number);
        abort();
      }
    }
  }

  return PACKET_SAME;
}


void rpkt_event_thickcell(PKT *pkt_ptr, double t_current)
/// Event handling for optically thick cells. Those cells are treated in a grey
/// approximation with electron scattering only.
/// The packet stays an R_PKT of same nu_cmf than before (coherent scattering)
/// but with different direction.
{
  #ifdef DEBUG_ON
    if (debuglevel == 2) printout("[debug] rpkt_event_thickcell:   electron scattering\n");
    pkt_ptr->interactions += 1;
    pkt_ptr->nscatterings += 1;
    pkt_ptr->last_event = 12;
    escounter += 1;
  #endif

  //pkt_ptr->nu_cmf = 3.7474058e+14;
  emitt_rpkt(pkt_ptr,t_current);
  /// Electron scattering does not modify the last emission flag
  //pkt_ptr->emissiontype = get_continuumindex(element,ion-1,lower);
  /// but it updates the last emission position
  pkt_ptr->em_pos[0] = pkt_ptr->pos[0];
  pkt_ptr->em_pos[1] = pkt_ptr->pos[1];
  pkt_ptr->em_pos[2] = pkt_ptr->pos[2];
  pkt_ptr->em_time = t_current;
}
