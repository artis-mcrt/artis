#include "sn3d.h"
#include <gsl/gsl_integration.h>


double do_ma(PKT *pkt_ptr, double t1, double t2, int timestep)
/// Material for handling activated macro atoms.
{
  double t_current;
  int end_packet;
  double zrand;
  double boundary_cross();
  int locate();
  int change_cell();
  int rpkt_event();
  int rlc_emiss_rpkt();

  double get_individ_rad_deexc(int i);
  double get_individ_internal_down_same(int i);
  double get_individ_internal_up_same(int i);

  void calculate_kappa_rpkt_cont(PKT *pkt_ptr, double t_current);
  double rad_deexcitation(PKT *pkt_ptr, int lower, double epsilon_trans, double statweight_target, int lineindex, double t_current);
  double rad_recombination(int modelgridindex, int lower, double epsilon_trans);
  double rad_excitation(PKT *pkt_ptr, int upper, double epsilon_trans, double statweight_target, int lineindex, double t_current);//, double T_R, double W);
  double photoionization(int modelgridindex, int phixstargetindex, double epsilon_trans);
  double col_excitation(int modelgridindex, int upper, int lineindex, double epsilon_trans);
  double col_ionization(int modelgridindex, int phixstargetindex, double epsilon_trans);
  double col_deexcitation(int modelgridindex, int lower, double epsilon_trans, double statweight_target, int lineindex);
  double col_recombination(int modelgridindex, int lower, double epsilon_trans);
  double get_levelpop(int element, int ion, int level);
  double stat_weight(int element, int ion, int level);
  double epsilon(int element, int ion, int level);
  int get_element(int element);
  int get_ionstage(int element, int ion);
  void emitt_rpkt(PKT *pkt_ptr, double t_current);
  int get_continuumindex(int element, int ion, int level);
  int get_bfcontinua(int element, int ion);
  int get_nphixstargets(int element, int ion, int level);
  int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);

  double rad_deexc,rad_recomb,col_deexc,col_recomb;
  double internal_down_same,internal_down_lower,internal_up_same,internal_up_higher;
  double individ_rad_deexc,individ_col_deexc,individ_internal_down_same,individ_internal_up_same;
  double R,C;
  double total_transitions;
  double rate;
  double t_mid;
  double oldnucmf;

  double nu_threshold, nu_max_phixs;
  double epsilon_current;
  double epsilon_target;
  double epsilon_trans;
  double statweight_target;

  int i,lineindex;
  int linelistindex = -99;
  int element,ion,level,upper,lower,phixstargetindex;
  int ndowntrans,nuptrans;
  int nlevels,ionisinglevels;

  end_packet = 0; ///means "keep working"
  t_current = t1; ///this will keep track of time in the calculation
  t_mid = time_step[timestep].mid;

  gsl_integration_workspace *wsp;
  gslintegration_paras intparas;
  double alpha_sp_integrand_gsl(double nu, void *paras);
  double alpha_sp_E_integrand_gsl(double nu, void *paras);
  gsl_function F_alpha_sp;
  //F_alpha_sp.function = &alpha_sp_integrand_gsl;
  F_alpha_sp.function = &alpha_sp_E_integrand_gsl;
  double intaccuracy = 1e-2;        /// Fractional accuracy of the integrator
  double error;
  double nu_lower,deltanu,alpha_sp,total_alpha_sp,alpha_sp_old,nuoffset;
  double calculate_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold);
  double get_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
  wsp = gsl_integration_workspace_alloc(1000);

  //printout("[debug] do MA\n");

  int modelgridindex = cell[pkt_ptr->where].modelgridindex;

  /// calculate occupation number for active MA level ////////////////////////////////////
  /// general QUESTION: is it better to calculate the n_1 (later the n_ionstage and
  /// U_ionstage) here where we need them or once in update_grid for each grid cell
  /// not sure whether this reduces the number of calculations, as number of grid cells
  /// is much larger than number of pellets (next question: connection to number of
  /// photons)
  double T_e = get_Te(modelgridindex);
  //double T_R = cell[pkt_ptr->where].T_R;
  //double W = cell[pkt_ptr->where].W;
  element = mastate[tid].element;

  /// dummy-initialize these to nonsense values, if something goes wrong with the real
  /// initialization we should see errors
  epsilon_trans = -100.;
  upper = -100;
  lower = -100;

  //debuglevel = 2;
  #ifdef DEBUG_ON
    if (debuglevel == 2) printout("[debug] =============entering do_ma\n");
    int jumps = 0;
    int jump = -99;
  #endif


  //debuglevel = 2;
  while (end_packet == 0)
  {
    ion = mastate[tid].ion;
    level = mastate[tid].level;
    mastate[tid].statweight = stat_weight(element,ion,level);
    ionisinglevels = get_bfcontinua(element,ion);
    //ionisinglevels = get_ionisinglevels(element,ion);

    /// Set this here to 1 to overcome problems in cells which have zero population
    /// in some ionisation stage. This is possible because the dependence on the
    /// originating levels population cancels out in the macroatom transition probabilities
    /// which are based on detailed balance.
    mastate[tid].nnlevel = 1.;
    //mastate[tid].nnlevel = get_levelpop(element,ion,level);

    #ifdef DEBUG_ON
      //if (element == 1 && ion == 1 && level == 149) debuglevel = 2;
      //if (element == 1 && ion == 1 && level == 150) debuglevel = 2;
      //if (element == 1 && ion == 1 && level == 151) debuglevel = 2;
      //if (element == 1 && ion == 1 && level == 61) debuglevel = 2;
      if (get_ionstage(element,ion) == 0 && level == 0)
      {
        printout("element %d, ion %d, level %d, ionstage %d\n",element,ion,level,get_ionstage(element,ion));
        debuglevel = 2;
      }
      if (ion >= get_nions(element))
      {
        printout("[fatal] do_ma: ion does not exist ... abort\n");
        abort();
      }
      if (debuglevel == 2) printout("[debug] do_ma: element %d, ion %d, level %d:\n",element,ion,level);
      if (debuglevel == 2) printout("[debug] do_ma: jumps = %d\n",jumps);
    #endif

    epsilon_current = epsilon(element,ion,level);
    ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
    nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
    //nlevels_nextion  ///not needed as long we only ionise to the ground state
    //nlevels_lowerion ///only needed if level = 0, this won't happen too often

    #ifdef DEBUG_ON
      if (debuglevel == 2) printout("[debug] do_ma: ndowntrans %d, nuptrans %d\n",ndowntrans,nuptrans);
      if (debuglevel == 2) printout("[debug] do_ma: col_deexc_stored %g\n",cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_deexc);
    #endif
    if (cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_deexc >= 0 && level != 0)
    {
      /// Take MA event rates from memory
      rad_deexc = cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_deexc;
      col_deexc = cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_deexc;
      rad_recomb = cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_recomb;
      col_recomb = cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_recomb;
      internal_down_same = cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_same;
      internal_up_same = cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_same;
      internal_down_lower = cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_lower;
      internal_up_higher = cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_higher;
    }
    else
    {
      /// If there are no precalculated rates available we must calculate them

      /// Downward transitions within the current ionisation stage:
      /// radiative/collisional deexcitation and internal downward jumps
      rad_deexc = 0.;
      col_deexc = 0.;
      internal_down_same = 0.;
      for (i = 1; i <= ndowntrans; i++)
      {
        lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
        epsilon_target = elements[element].ions[ion].levels[level].downtrans[i].epsilon;
        statweight_target = elements[element].ions[ion].levels[level].downtrans[i].stat_weight;
        lineindex = elements[element].ions[ion].levels[level].downtrans[i].lineindex;
        epsilon_trans = epsilon_current - epsilon_target;
        R = rad_deexcitation(pkt_ptr,lower,epsilon_trans,statweight_target,lineindex,t_mid);
        C = col_deexcitation(modelgridindex,lower,epsilon_trans,statweight_target,lineindex);

        individ_rad_deexc = R * epsilon_trans;
        individ_col_deexc = C * epsilon_trans;
        individ_internal_down_same = (R + C) * epsilon_target;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc[i] = individ_rad_deexc;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[i] = individ_internal_down_same;

        rad_deexc += individ_rad_deexc;
        col_deexc += individ_col_deexc;
        internal_down_same += individ_internal_down_same;

        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("checking downtrans %d to level %d: R %g, C %g, epsilon_trans %g\n",i,lower,R,C,epsilon_trans);
        #endif
      }

      /// Downward transitions to lower ionisation stages:
      /// radiative/collisional recombination and internal downward jumps
      rad_recomb = 0.;
      col_recomb = 0.;
      internal_down_lower = 0.;
      if (ion > 0) ///checks only if there is a lower ion, doesn't make sure that Z(ion)=Z(ion-1)+1
      {
        //nlevels = get_nlevels(element,ion-1);
        nlevels = get_bfcontinua(element,ion-1);
        //nlevels = get_ionisinglevels(element,ion-1);
        for (lower = 0; lower < nlevels; lower++)
        {
          epsilon_target = epsilon(element,ion-1,lower);
          epsilon_trans = epsilon_current - epsilon_target;
          R = rad_recombination(modelgridindex,lower,epsilon_trans);
          //printout("rad recombination of element %d, ion %d, level %d, to lower level %d has rate %g\n",element,ion,level,lower,R);
          C = col_recombination(modelgridindex,lower,epsilon_trans);
          rad_recomb += R * epsilon_trans;
          col_recomb += C * epsilon_trans;
          internal_down_lower += (R + C) * epsilon_target;
        }
      }

      /// Calculate sum for upward internal transitions
      /// transitions within the current ionisation stage
      internal_up_same = 0.;
      for (i = 1; i <= nuptrans; i++)
      {
        upper = elements[element].ions[ion].levels[level].uptrans[i].targetlevel;
        epsilon_target = elements[element].ions[ion].levels[level].uptrans[i].epsilon;
        statweight_target = elements[element].ions[ion].levels[level].uptrans[i].stat_weight;
        lineindex = elements[element].ions[ion].levels[level].uptrans[i].lineindex;
        epsilon_trans = epsilon_target - epsilon_current;
        R = rad_excitation(pkt_ptr,upper,epsilon_trans,statweight_target,lineindex,t_mid);//,T_R,W);
        C = col_excitation(modelgridindex,upper,lineindex,epsilon_trans);

        //individ_internal_up_same = (C) * epsilon_current;
        individ_internal_up_same = (R + C) * epsilon_current;
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same[i] = individ_internal_up_same;

        internal_up_same += individ_internal_up_same;

        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("checking uptrans %d to level %d: R %g, C %g, epsilon_trans %g\n",i,upper,R,C,epsilon_trans);
          if (!isfinite(internal_up_same)) {printout("fatal: internal_up_same has nan contribution\n");}
        #endif
      }

      /// Transitions to higher ionisation stages
      internal_up_higher = 0.;
      if (ion < get_nions(element)-1 && level < ionisinglevels)  //&& get_ionstage(element,ion) < get_element(element)+1)
      {
        for (phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
        {
          upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
          epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
          R = photoionization(modelgridindex,phixstargetindex,epsilon_trans);
          C = col_ionization(modelgridindex,phixstargetindex,epsilon_trans);
          internal_up_higher += (R + C) * epsilon_current;
        }
      }

      /// and store them to memory
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_deexc = rad_deexc;
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_deexc = col_deexc;
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].rad_recomb = rad_recomb;
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_recomb = col_recomb;
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_same = internal_down_same;
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_same = internal_up_same;
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_down_lower = internal_down_lower;
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_higher = internal_up_higher;
    }

    /// select transition according to probabilities /////////////////////////////////////
    zrand = gsl_rng_uniform(rng);
    //printout("zrand %g\n",zrand);

    //internal_down_same = internal_down_lower = internal_up_same = internal_up_higher = 0.; ///DEBUG ONLY
    total_transitions = rad_deexc + col_deexc + internal_down_same + rad_recomb + col_recomb + internal_down_lower + internal_up_same + internal_up_higher;
    #ifdef DEBUG_ON
      if (debuglevel == 2) printout("[debug] do_ma: element %d, ion %d, level %d\n",element,ion,level);
      if (debuglevel == 2) printout("[debug] do_ma:   rad_deexc %g, col_deexc %g, internal_down_same %g, rad_recomb %g, col_recomb %g, internal_down_lower %g, internal_up_same %g, internal_up_higher %g\n",rad_deexc,col_deexc,internal_down_same,rad_recomb,col_recomb, internal_down_lower, internal_up_same,internal_up_higher);
      if (total_transitions <= 0.)
      {
        printout("[debug] do_ma: element %d, ion %d, level %d, total_transitions = %g\n",element,ion,level,total_transitions);
        printout("[debug]    col_deexc %g, col_recomb %g, rad_deexc %g, rad_recomb %g\n",col_deexc,col_recomb,rad_deexc,rad_recomb);
        printout("[debug]    internal_down_same %g, internal_down_lower %g\n",internal_down_same,internal_down_lower);
        printout("[debug]    internal_up_same %g, internal_up_higher %g\n",internal_up_same,internal_up_higher);
        printout("[debug]    zrand %g\n",zrand);
        printout("[debug]    jumps %d, jump %d\n",jumps,jump);
        printout("[debug]    pkt_ptr->number %d, pkt_ptr->where %d\n",pkt_ptr->number,pkt_ptr->where);
        printout("[debug]    groundlevelpop of current ion in current cell %g\n",modelgrid[cell[pkt_ptr->where].modelgridindex].composition[element].groundlevelpop[ion]);
        printout("[debug]    levelpop %g\n",mastate[tid].nnlevel);

        R = 0.0;
        C = 0.0;
        for (phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
        {
          if (get_phixsupperlevel(element,ion,level,phixstargetindex) == upper)
          {
            epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            R = photoionization(modelgridindex,phixstargetindex,epsilon_trans);
            C = col_ionization(modelgridindex,phixstargetindex,epsilon_trans);
            printout("epsilon_current %g, epsilon_trans %g, photion %g, colion %g, internal_up_higher %g, saved_internal_up_higher %g\n",epsilon_current,epsilon_trans,R,C,(R + C) * epsilon_current,cellhistory[tid].chelements[element].chions[ion].chlevels[level].internal_up_higher);
            break;
          }
        }

        double T_R = get_TR(modelgridindex);
        double W = get_W(modelgridindex);
        double interpolate_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, double T);
        printout("modelgridindex %d, T_R %g, T_e %g, W %g, T_J %g\n",modelgridindex,T_R,T_e,W,get_TJ(modelgridindex));
        double gammacorr = W*interpolate_corrphotoioncoeff(element,ion,level,0,T_R);
        int index_in_groundlevelcontestimor = elements[element].ions[ion].levels[level].closestgroundlevelcont;
        double renorm  = corrphotoionrenorm[modelgridindex*nelements*maxion+index_in_groundlevelcontestimor];
        printout("gammacorr %g, index %d, renorm %g, total %g\n",gammacorr,index_in_groundlevelcontestimor,renorm,gammacorr*renorm);


        //abort();
      }
    #endif
    if (zrand*total_transitions < rad_deexc)
    {
      ///radiative deexcitation of MA: emitt rpkt
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   radiative deexcitation\n");
        if (debuglevel == 2) printout("[debug] do_ma:   jumps = %d\n",jumps);
      #endif

      ///randomly select which line transitions occurs
      zrand = gsl_rng_uniform(rng);
      //zrand = 1. - 1e-14; /// ONLY FOR DEBUG!!!
      rate = 0.;
      for (i = 1; i <= ndowntrans; i++)
      {
        rate += get_individ_rad_deexc(i);
        if (zrand*rad_deexc < rate)
        {
          #ifdef DEBUG_ON
            if (debuglevel == 2) printout("[debug] do_ma:   jump to level %d\n",lower);
          #endif
          lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
          epsilon_target = elements[element].ions[ion].levels[level].downtrans[i].epsilon;
          linelistindex = elements[element].ions[ion].levels[level].downtrans[i].lineindex;
          epsilon_trans = epsilon_current - epsilon_target;
          //linelistindex = elements[element].ions[ion].levels[level].transitions[level-lower-1].linelistindex;
          #ifdef RECORD_LINESTAT
            if (tid == 0) ecounter[linelistindex] += 1;    /// This way we will only record line statistics from OMP-thread 0
                                                           /// With an atomic pragma or a thread-private structure with subsequent
                                                           /// reduction this could be extended to all threads. However, I'm not
                                                           /// sure if this is worth the additional computational expenses.
          #endif
          break;
        }
      }
      if (pkt_ptr->last_event == 1) oldnucmf = pkt_ptr->nu_cmf;
      pkt_ptr->nu_cmf = epsilon_trans/H;
      //pkt_ptr->nu_cmf = 3.7474058e+14;
      //if (tid == 0)
      //{
        if (pkt_ptr->last_event == 1)
        {
          if (oldnucmf < pkt_ptr->nu_cmf) upscatter += 1;
          else downscatter += 1;
        }
      //}

      #ifdef DEBUG_ON
        if (!isfinite(pkt_ptr->nu_cmf))
        {
          printout("[fatal] rad deexcitation of MA: selected frequency not finite ... abort\n");
          abort();
        }
        if (i > ndowntrans)
        {
          printout("[fatal] problem in selecting radiative downward transition of MA zrand %g, rate %g, rad_deexc %g, i %d, ndowntrans %d\n",zrand,rate,rad_deexc,i,ndowntrans);
          printout("[fatal] total_transitions %g, element %d, ion %d, level %d\n",total_transitions,element,ion,level);
          abort();
        }
        if (debuglevel == 2)printout("[debug] do_ma: calculate_kappa_rpkt_cont after MA deactivation\n");
        //if (tid == 0) ma_stat_deactivation_bb += 1;
        ma_stat_deactivation_bb += 1;
        mastate[tid].lastaction = MA_ACTION_RADDEEXC;
        pkt_ptr->interactions += 1;
        pkt_ptr->last_event = 0;
      #endif

      /// Emitt the rpkt in a random direction
      emitt_rpkt(pkt_ptr,t_current);
      if (linelistindex == mastate[tid].activatingline)
      {
        resonancescatterings += 1;
      }
      else calculate_kappa_rpkt_cont(pkt_ptr,t_current);
      /// NB: the r-pkt can only interact with lines redder than the current one
      pkt_ptr->next_trans = linelistindex + 1;
      pkt_ptr->emissiontype = linelistindex;
      pkt_ptr->em_pos[0] = pkt_ptr->pos[0];
      pkt_ptr->em_pos[1] = pkt_ptr->pos[1];
      pkt_ptr->em_pos[2] = pkt_ptr->pos[2];
      pkt_ptr->em_time = t_current;
      pkt_ptr->nscatterings = 0;
      //printout("next possible line encounter %d\n",pkt_ptr->next_trans);
      end_packet = 1;
      #ifndef FORCE_LTE
        //matotem[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif
    }
    else if (zrand*total_transitions < rad_deexc + col_deexc)
    {
      ///collisional deexcitation of macro atom => convert the packet into a k-packet
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   collisonal deexcitation\n");
        if (debuglevel == 2) printout("[debug] do_ma: jumps = %d\n",jumps);
        //if (tid == 0) ma_stat_deactivation_colldeexc += 1;
        ma_stat_deactivation_colldeexc += 1;
        mastate[tid].lastaction = MA_ACTION_COLDEEXC;
        pkt_ptr->interactions += 1;
        pkt_ptr->last_event = 10;
      #endif
      pkt_ptr->type = TYPE_KPKT;
      end_packet = 1;
      #ifndef FORCE_LTE
        //matotem[pkt_ptr->where] += pkt_ptr->e_cmf;
        #ifdef _OPENMP
          #pragma omp atomic
        #endif
        colheatingestimator[modelgridindex] += pkt_ptr->e_cmf;
      #endif
    }
    else if (zrand*total_transitions < rad_deexc + col_deexc + internal_down_same)
    {
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   internal downward jump within current ionstage\n");
        mastate[tid].lastaction = MA_ACTION_INTERNALDOWNSAME;
        pkt_ptr->interactions += 1;
        jumps += 1;
        jump = 0;
      #endif

      /// Randomly select the occuring transition
      zrand = gsl_rng_uniform(rng);
      lower = -99;
      rate = 0.;
      for (i = 1; i <= ndowntrans; i++)
      {
        rate += get_individ_internal_down_same(i);
        if (zrand*internal_down_same < rate)
        {
          lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
          break;
        }
      }
      /// and set the macroatom's new state
      mastate[tid].ion = ion;
      mastate[tid].level = lower;

      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   to level %d\n",lower);
        if (get_ionstage(element,ion) == 0 && lower == 0)
        {
          printout("internal downward transition to ground level occured ... abort\n");
          printout("element %d, ion %d, level %d, lower %d\n",element,ion,level,lower);
          printout("Z %d, ionstage %d, energy %g\n",elements[element].anumber,get_ionstage(element,ion),elements[element].ions[ion].levels[lower].epsilon);
          printout("[debug] do_ma:   internal downward jump within current ionstage\n");
          abort();
        }
      #endif
    }
    else if (zrand*total_transitions < rad_deexc + col_deexc + internal_down_same + rad_recomb)
    {
      /// Radiative recombination of MA: emitt a continuum-rpkt
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   radiative recombination\n");
        if (debuglevel == 2) printout("[debug] do_ma:   jumps = %d\n",jumps);
        if (debuglevel == 2) printout("[debug] do_ma:   element %d, ion %d, level %d\n",element,ion,level);
      #endif

      /// Randomly select a continuum
      zrand = gsl_rng_uniform(rng);
      //zrand = 1. - 1e-14;
      rate = 0;
      //nlevels = get_nlevels(element,ion-1);
      nlevels = get_bfcontinua(element,ion-1);
      //nlevels = get_ionisinglevels(element,ion-1);
      for (lower=0; lower < nlevels; lower++)
      {
        epsilon_trans = epsilon_current - epsilon(element,ion-1,lower);
        R = rad_recombination(modelgridindex,lower,epsilon_trans);
        rate += R * epsilon_trans;
        #ifdef DEBUG_ON
          if (debuglevel == 2)
          {
            printout("[debug] do_ma:   R %g, deltae %g\n",R,(epsilon(element,ion,level)-epsilon(element,ion-1,lower)));
            printout("[debug] do_ma:   rate to level %d of ion %d = %g\n",lower,ion-1,rate);
            printout("[debug] do_ma:   zrand*rad_recomb = %g\n",zrand*rad_recomb);
          }
        #endif
        if (zrand*rad_recomb < rate) break;
      }
      /// and set its threshold frequency
      nu_threshold = epsilon_trans/H;
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   going to level %d of ion %d of element %d\n",lower,ion-1,element);
      #endif

      /// Then randomly sample the packets frequency according to the continuums
      /// energy distribution and set some flags
      //zrand = gsl_rng_uniform(rng);
      //zrand = 1. - zrand;  /// Make sure that 0 < zrand <= 1
      //pkt_ptr->nu_cmf = nu_threshold * (1 - KB*T_e/H/nu_threshold*log(zrand));
      //pkt_ptr->nu_cmf = nu_threshold;


      zrand = gsl_rng_uniform(rng);
      zrand = 1. - zrand;  /// Make sure that 0 < zrand <= 1
      mastate[tid].element = element;
      mastate[tid].ion = ion-1;
      mastate[tid].level = lower;
      intparas.T = T_e;
      intparas.nu_edge = nu_threshold;   /// Global variable which passes the threshold to the integrator
      F_alpha_sp.params = &intparas;
      deltanu = nu_threshold * NPHIXSNUINCREMENT;
      nu_max_phixs = nu_threshold * (1.0 + NPHIXSNUINCREMENT * (NPHIXSPOINTS - 1)); //nu of the uppermost point in the phixs table
      gsl_integration_qag(&F_alpha_sp, nu_threshold, nu_max_phixs, 0, intaccuracy, 1000, 6, wsp, &total_alpha_sp, &error);
      alpha_sp = total_alpha_sp;
      for (i = 0; i < NPHIXSPOINTS; i++)
      // LJS: this loop could probably be made a bit faster
      // use the overlap with the previous integral and add on a piece each time instead of recalculating the
      // integral over the entire region
      {
        // the reason the lower limit of integration is incremented is that most of the probability is in the low
        // frequency end, so this minimizes the number of iterations needed
        alpha_sp_old = alpha_sp;
        if (i > 0)
        {
          nu_lower = nu_threshold + i*deltanu;
          /// Spontaneous recombination and bf-cooling coefficient don't depend on the cutted radiation field
          gsl_integration_qag(&F_alpha_sp, nu_lower, nu_max_phixs, 0, intaccuracy, 1000, 6, wsp, &alpha_sp, &error);
          //alpha_sp *= FOURPI * sf;
          //if (zrand > alpha_sp/get_spontrecombcoeff(element,ion-1,lower,pkt_ptr->where)) break;
        }
        //printout("[debug] macroatom: zrand %g, step %d, alpha_sp %g, total_alpha_sp %g, alpha_sp/total_alpha_sp %g, nu_lower %g\n",zrand,i,alpha_sp,total_alpha_sp,alpha_sp/total_alpha_sp,nu_lower);
        if (zrand >= alpha_sp/total_alpha_sp) break;
      }
      if (i == NPHIXSPOINTS)
      {
        nu_lower = nu_max_phixs;
      }
      else if (i > 0)
      {
        nuoffset = (total_alpha_sp*zrand - alpha_sp_old) / (alpha_sp-alpha_sp_old) * deltanu;
        nu_lower = nu_threshold + (i-1)*deltanu + nuoffset;
      }
      else
        nu_lower = nu_threshold;
      pkt_ptr->nu_cmf = nu_lower;
      //printout("nu_lower %g, nu_threshold %g, nu_left %g, nu_right %g\n",nu_lower,nu_threshold,nu_threshold+(i-1)*deltanu,nu_threshold+(i)*deltanu);


      /*
      ///k-pkt emission rule
      double bfcooling_coeff,total_bfcooling_coeff,bfcooling_coeff_old;
      int ii;
      double bfcooling_integrand_gsl_2(double nu, void *paras);
      gsl_function F_bfcooling;
      F_bfcooling.function = &bfcooling_integrand_gsl_2;
      zrand = gsl_rng_uniform(rng);
      zrand = 1. - zrand;  /// Make sure that 0 < zrand <= 1
      mastate[tid].element = element;
      mastate[tid].ion = ion-1;
      mastate[tid].level = lower;
      intparas.T = T_e;
      intparas.nu_edge = nu_threshold;   /// Global variable which passes the threshold to the integrator
      F_bfcooling.params = &intparas;
      deltanu = nu_threshold * NPHIXSNUINCREMENT;
      gsl_integration_qag(&F_bfcooling, nu_threshold, nu_max_phixs, 0, intaccuracy, 1000, 6, wsp, &total_bfcooling_coeff, &error);
      bfcooling_coeff = total_bfcooling_coeff;
      for (ii= 0; ii < NPHIXSPOINTS; ii++)
      {
        bfcooling_coeff_old = bfcooling_coeff;
        if (ii > 0)
        {
          nu_lower = nu_threshold + ii*deltanu;
          /// Spontaneous recombination and bf-cooling coefficient don't depend on the cutted radiation field
          gsl_integration_qag(&F_bfcooling, nu_lower, nu_max_phixs, 0, intaccuracy, 1000, 6, wsp, &bfcooling_coeff, &error);
          //bfcooling_coeff *= FOURPI * sf;
          //if (zrand > bfcooling_coeff/get_bfcooling(element,ion,level,pkt_ptr->where)) break;
        }
        //printout("zrand %g, bfcooling_coeff %g, total_bfcooling_coeff %g, nu_lower %g\n",zrand,bfcooling_coeff,total_bfcooling_coeff,nu_lower);
        if (zrand >= bfcooling_coeff/total_bfcooling_coeff) break;
      }
      if (ii==NPHIXSPOINTS)
      {
        printout("kpkt emitts bf-photon at upper limit\n");
        nu_lower = nu_threshold * (1.0 + NPHIXSNUINCREMENT * (NPHIXSPOINTS - 1)); // + ii*deltanu;
      }
      else if (ii > 0)
      {
        nuoffset = (total_bfcooling_coeff*zrand - bfcooling_coeff_old) / (bfcooling_coeff-bfcooling_coeff_old) * deltanu;
        nu_lower = nu_threshold + (ii-1)*deltanu + nuoffset;
      }
      else
        nu_lower = nu_threshold; // + ii*deltanu;
      //printout("emitt at nu %g\n",nu_lower);
      pkt_ptr->nu_cmf = nu_lower;
        */

      #ifndef FORCE_LTE
        //mabfcount[pkt_ptr->where] += pkt_ptr->e_cmf;
        //mabfcount_thermal[pkt_ptr->where] += pkt_ptr->e_cmf*(1-nu_threshold/pkt_ptr->nu_cmf);
        //matotem[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif

/*      if (element == 6)
      {
        //printout("%g, %g, %g\n",pkt_ptr->e_cmf,nu_threshold,pkt_ptr->e_cmf/nu_threshold/H);
        cell[pkt_ptr->where].radrecomb[ion-1] += pkt_ptr->e_cmf/pkt_ptr->nu_cmf/H;
      }*/

      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   pkt_ptr->nu_cmf %g\n",pkt_ptr->nu_cmf);
        if (!isfinite(pkt_ptr->nu_cmf))
        {
          printout("[fatal] rad recombination of MA: selected frequency not finite ... abort\n");
          abort();
        }
        //if (tid == 0) ma_stat_deactivation_fb += 1;
        ma_stat_deactivation_fb += 1;
        mastate[tid].lastaction = MA_ACTION_RADRECOMB;
        pkt_ptr->interactions += 1;
        pkt_ptr->last_event = 2;
        if (debuglevel == 2) printout("[debug] do_ma: calculate_kappa_rpkt_cont after MA recombination\n");
      #endif

      /// Finally emit the packet into a randomly chosen direction, update the continuum opacity and set some flags
      emitt_rpkt(pkt_ptr,t_current);
      calculate_kappa_rpkt_cont(pkt_ptr,t_current);
      pkt_ptr->next_trans = 0;       /// continuum transition, no restrictions for further line interactions
      pkt_ptr->emissiontype = get_continuumindex(element,ion-1,lower);
      pkt_ptr->em_pos[0] = pkt_ptr->pos[0];
      pkt_ptr->em_pos[1] = pkt_ptr->pos[1];
      pkt_ptr->em_pos[2] = pkt_ptr->pos[2];
      pkt_ptr->em_time = t_current;
      pkt_ptr->nscatterings = 0;
      end_packet = 1;
    }
    else if (zrand*total_transitions < rad_deexc + col_deexc + internal_down_same + rad_recomb + col_recomb)
    {
      ///collisional recombination of macro atom => convert the packet into a k-packet
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   collisonal recombination\n");
        if (debuglevel == 2) printout("[debug] do_ma: jumps = %d\n",jumps);
        //if (tid == 0) ma_stat_deactivation_collrecomb += 1;
        ma_stat_deactivation_collrecomb += 1;
        mastate[tid].lastaction = MA_ACTION_COLRECOMB;
        pkt_ptr->interactions += 1;
        pkt_ptr->last_event = 11;
      #endif
      pkt_ptr->type = TYPE_KPKT;
      end_packet = 1;
      #ifndef FORCE_LTE
        //matotem[pkt_ptr->where] += pkt_ptr->e_cmf;
        #ifdef _OPENMP
          #pragma omp atomic
        #endif
        colheatingestimator[modelgridindex] += pkt_ptr->e_cmf;
      #endif
    }
    else if (zrand*total_transitions < rad_deexc + col_deexc + internal_down_same + rad_recomb + col_recomb + internal_down_lower)
    {
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   internal downward jump to lower ionstage\n");
        mastate[tid].lastaction = MA_ACTION_INTERNALDOWNLOWER;
        pkt_ptr->interactions += 1;
        jumps += 1;
        jump = 1;
      #endif

      /// Randomly select the occuring transition
      zrand = gsl_rng_uniform(rng);
      //zrand = 1. - 1e-14;
      rate = 0.;
      //nlevels = get_nlevels(element,ion-1);
      nlevels = get_bfcontinua(element,ion-1);
      //nlevels = get_ionisinglevels(element,ion-1);
      for (lower=0; lower < nlevels; lower++)
      {
        epsilon_target = epsilon(element,ion-1,lower);
        epsilon_trans = epsilon_current - epsilon_target;
        R = rad_recombination(modelgridindex,lower,epsilon_trans);
        C = col_recombination(modelgridindex,lower,epsilon_trans);
        rate += (R + C) * epsilon_target;
        if (zrand*internal_down_lower < rate) break;
      }
      /// and set the macroatom's new state
      mastate[tid].ion = ion-1;
      mastate[tid].level = lower;

      #ifdef DEBUG_ON
        if (lower >= nlevels)
        {
          printout("internal_down_lower  %g\n",internal_down_lower);
          printout("abort at rate %g, zrand %g\n",rate,zrand);
          abort();
        }
        if (get_ionstage(element,ion-1) == 0 && lower == 0)
	  //        if (ion-1 == 0 && lower == 0)
        {
          printout("internal downward transition to ground level occured ... abort\n");
          printout("element %d, ion %d, level %d, lower %d\n",element,ion,level,lower);
          printout("Z %d, ionstage %d, energy %g\n",elements[element].anumber,get_ionstage(element,ion-1),elements[element].ions[ion-1].levels[lower].epsilon);
          printout("[debug] do_ma:   internal downward jump to lower ionstage\n");
          abort();
        }
      #endif
    }
    else if (zrand*total_transitions < rad_deexc + col_deexc + internal_down_same + rad_recomb + col_recomb + internal_down_lower + internal_up_same)
    {
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   internal upward jump within current ionstage\n");
        mastate[tid].lastaction = MA_ACTION_INTERNALUPSAME;
        pkt_ptr->interactions += 1;
        jumps += 1;
        jump = 2;
      #endif

      ///randomly select the occuring transition
      zrand = gsl_rng_uniform(rng);
      upper = -99;
      rate = 0.;
      for (i = 1; i <= nuptrans; i++)
      {
        rate += get_individ_internal_up_same(i);
        if (zrand*internal_up_same < rate)
        {
          upper = elements[element].ions[ion].levels[level].uptrans[i].targetlevel;
          break;
        }
      }
      ///and set the macroatom's new state
      mastate[tid].ion = ion;
      mastate[tid].level = upper;
    }
    #ifdef DEBUG_ON
      else if (zrand*total_transitions < total_transitions)
    #else
      else
    #endif
    {
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_ma:   internal upward jump to next ionstage\n");
        mastate[tid].lastaction = MA_ACTION_INTERNALUPHIGHER;
        pkt_ptr->interactions += 1;
        jumps += 1;
        jump = 3;
      #endif

      /// Randomly select the occuring transition
      zrand = gsl_rng_uniform(rng);
      rate = 0.;
      for (phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
      {
        upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
        epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
        R = photoionization(modelgridindex,phixstargetindex,epsilon_trans);
        C = col_ionization(modelgridindex,phixstargetindex,epsilon_trans);
        rate += (R + C) * epsilon_current;
        if (zrand*internal_up_higher < rate) break;
      }

      /// and set the macroatom's new state
      mastate[tid].ion = ion+1;
      mastate[tid].level = upper;
    }
    #ifdef DEBUG_ON
      else
      {
        printout("[fatal] do_ma: problem with random numbers .. abort\n");
        //printout("[debug]    rad_down %g, col_down %g, internal_down %g, internal_up %g\n",rad_down,col_down,internal_down,internal_up);
        printout("[debug] do_ma: element %d, ion %d, level %d, total_transitions %g\n",element,ion,level,total_transitions);
        printout("[debug]    col_deexc %g, col_recomb %g, rad_deexc %g, rad_recomb %g\n",col_deexc,col_recomb,rad_deexc,rad_recomb);
        printout("[debug]    internal_down_same %g, internal_down_lower %g\n",internal_down_same,internal_down_lower);
        printout("[debug]    internal_up_same %g, internal_up_higher %g\n",internal_up_same,internal_up_higher);
        printout("[debug]    zrand %g\n",zrand);
        printout("[debug]    jumps %d\n",jumps);
        printout("[debug]    pkt_ptr->number %d\n",pkt_ptr->number);

        debuglevel = 777;

        if (ion > 0) ///checks only if there is a lower ion, doesn't make sure that Z(ion)=Z(ion-1)+1
        {
          printout("[debug]    check recombination\n");
          //nlevels = get_nlevels(element,ion-1);
          nlevels = get_bfcontinua(element,ion-1);
          //nlevels = get_ionisinglevels(element,ion-1);
          for (lower=0; lower < nlevels; lower++)
          {
            epsilon_target = epsilon(element,ion-1,lower);
            epsilon_trans = epsilon_current - epsilon_target;
            R = rad_recombination(modelgridindex,lower,epsilon_trans);
            C = col_recombination(modelgridindex,lower,epsilon_trans);
            printout("[debug]    recombination to ion %d, level %d, epsilon_target %g, epsilon_trans %g, R %g, C %g\n",ion-1,lower,epsilon_target,epsilon_trans,R,C);
          }
        }

        printout("[debug]    check deexcitation\n");
        printout("[debug]    ndowntrans %d %d\n",ndowntrans,elements[element].ions[ion].levels[level].downtrans[0].targetlevel);
        for (i = 1; i <= ndowntrans; i++)
        {
          lower = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
          epsilon_target = elements[element].ions[ion].levels[level].downtrans[i].epsilon;
          statweight_target = elements[element].ions[ion].levels[level].downtrans[i].stat_weight;
          lineindex = elements[element].ions[ion].levels[level].downtrans[i].lineindex;
          epsilon_trans = epsilon_current - epsilon_target;
          R = rad_deexcitation(pkt_ptr,lower,epsilon_trans,statweight_target,lineindex,t_mid);
          C = col_deexcitation(modelgridindex,lower,epsilon_trans,statweight_target,lineindex);
          printout("[debug]    deexcitation to level %d, epsilon_target %g, epsilon_trans %g, R %g, C %g\n",lower,epsilon_target,epsilon_trans,R,C);
        }

        printout("[debug]    check excitation\n");
        printout("[debug]    nuptrans %d %d\n",nuptrans,elements[element].ions[ion].levels[level].uptrans[0].targetlevel);
        for (i = 1; i <= nuptrans; i++)
        {
          upper = elements[element].ions[ion].levels[level].uptrans[i].targetlevel;
          epsilon_target = elements[element].ions[ion].levels[level].uptrans[i].epsilon;
          statweight_target = elements[element].ions[ion].levels[level].uptrans[i].stat_weight;
          lineindex = elements[element].ions[ion].levels[level].uptrans[i].lineindex;
          epsilon_trans = epsilon_target - epsilon_current;
          R = rad_excitation(pkt_ptr,upper,epsilon_trans,statweight_target,lineindex,t_mid);//,T_R,W);
          C = col_excitation(modelgridindex,upper,lineindex,epsilon_trans);
          printout("[debug]    excitation to level %d, epsilon_trans %g, R %g, C %g\n",upper,epsilon_trans,R,C);
        }

        if (ion < get_nions(element)-1)  //&& get_ionstage(element,ion) < get_element(element)+1)
        {
          printout("[debug]    check ionisation\n");
          R = 0.0;
          C = 0.0;
          for (phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
          {
            upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
            epsilon_trans = epsilon(element,ion+1,upper) - epsilon_current;
            R = photoionization(modelgridindex,phixstargetindex,epsilon_trans);
            C = col_ionization(modelgridindex,phixstargetindex,epsilon_trans);
            printout("[debug]    ionisation to ion %d, level %d, epsilon_trans %g, R %g, C %g\n",ion+1,upper,epsilon_trans,R,C);
            break;
          }
        }

        //abort();
      }
    #endif

  }///endwhile

  #ifdef DEBUG_ON
    //debuglevel = 4;
    //debuglevel = 2;
  #endif

  /// procedure ends only after a change to r or k packets has taken place and
  /// returns then the actual time, which is the same as the input t1
  /// internal transitions are carried out until a type change occurs
  gsl_integration_workspace_free(wsp);
  return t_current;
}


/// Calculation of radiative rates ///////////////////////////////////////////////////////
///****************************************************************************
double rad_deexcitation(PKT *pkt_ptr, int lower, double epsilon_trans, double statweight_target, int lineindex, double t_current)
///radiative deexcitation rate: paperII 3.5.2
/// n_1 - occupation number of ground state
{
  //double einstein_spontaneous_emission(int element, int ion, int upper, int lower);
  double einstein_spontaneous_emission(int lineindex);
  double get_levelpop(int element, int ion, int level);
  //double boltzmann(PKT *pkt_ptr, int targetlevel);
  double epsilon(int element, int ion, int level);
  double A_ul,B_ul,B_lu;
  double n_u,n_l;
  double nu_trans;
  double tau_sobolev,beta;
  double R;

  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int upper = mastate[tid].level;

  int modelgridindex = cell[pkt_ptr->where].modelgridindex;

  #ifdef DEBUG_ON
    if (upper <= lower)
    {
      printout("[fatal] rad_deexcitation: pkt %d tried to calculate upward transition ... abort\n",pkt_ptr->number);
      abort();
    }
  #endif

  nu_trans = epsilon_trans/H;
  //A_ul = einstein_spontaneous_emission(element,ion,upper,lower);
  A_ul = einstein_spontaneous_emission(lineindex);
  B_ul = CLIGHTSQUAREDOVERTWOH/pow(nu_trans,3) * A_ul;
  B_lu = mastate[tid].statweight/statweight_target * B_ul;
  //double g_ratio = mastate[tid].statweight/statweight_target;
  //B_lu = g_ratio * B_ul;

  n_u = get_levelpop(element,ion,upper);
  n_l = get_levelpop(element,ion,lower);
  //double T_R = cell[pkt_ptr->where].T_R;
  //double W = cell[pkt_ptr->where].W;
  //n_l = n_u/W / g_ratio * exp(epsilon_trans/KB/T_R);
  tau_sobolev = (B_lu*n_l - B_ul*n_u) * HCLIGHTOVERFOURPI * t_current;

  #ifdef DEBUG_ON
    if (tau_sobolev <= 0)
    {
      if (SILENT == 0) printout("[warning] rad_deexcitation: tau_sobolev %g <= 0, set beta=1\n",tau_sobolev);
      if (SILENT == 0) printout("[warning] rad_deexcitation: element %d, ion %d, upper %d, lower %d\n",element,ion,upper,lower);
      if (SILENT == 0) printout("[warning] rad_deexcitation: n_l %g, n_u %g, B_lu %g, B_ul %g\n",n_l,n_u,B_lu,B_ul);
      if (SILENT == 0) printout("[warning] rad_deexcitation: T_e %g, T_R %g, W %g in model cell %d\n",get_Te(modelgridindex),get_TR(modelgridindex),get_W(modelgridindex),modelgridindex);
      if (SILENT == 0) printout("[warning] rad_deexcitation: pkt_ptr->number %d, last_event %d\n",pkt_ptr->number,pkt_ptr->last_event);
      beta = 1.0;
      //printout("[fatal] rad_excitation: tau_sobolev <= 0 ... %g abort",tau_sobolev);
      //abort();
    }
    else
    {
      beta = 1.0/tau_sobolev * (1 - exp(-tau_sobolev));
    }
  #else
    //beta = 1.0; ///FOR DEBUGGING ONLY
    beta = 1.0/tau_sobolev * (1 - exp(-tau_sobolev));
  #endif

  R = A_ul * beta * mastate[tid].nnlevel;

  #ifdef DEBUG_ON
    if (debuglevel == 2) printout("[debug] rad_rates_down: element, ion, upper, lower %d, %d, %d, %d\n",element,ion,upper,lower);
    //printout("[debug] rad_rates_down: tau_sobolev, beta %g, %g\n",tau_sobolev,beta);
    //printout("[debug] rad_rates_down: nne %g \n",cell[pkt_ptr->where].nne);
    if (debuglevel == 2) printout("[debug] rad_deexc: A_ul %g, tau_sobolev %g, n_u %g\n",A_ul,tau_sobolev,n_u);
    if (debuglevel == 777) printout("[debug] rad_deexc: A_ul %g, tau_sobolev %g, n_u %g\n",A_ul,tau_sobolev,n_u);
    if (!isfinite(R)) {printout("fatal a1: abort\n"); abort();}
  #endif

  return R;
}



///***************************************************************************/
double rad_excitation(PKT *pkt_ptr, int upper, double epsilon_trans, double statweight_target, int lineindex, double t_current)//, double T_R, double W)
///radiative excitation rate: paperII 3.5.2
/// n_1 - occupation number of ground state
{
  //double einstein_spontaneous_emission(int element, int ion, int upper, int lower);
  double einstein_spontaneous_emission(int lineindex);
  double get_levelpop(int element, int ion, int level);
  double radfield(double nu, int modelgridindex);
  double A_ul,B_ul,B_lu;
  double n_u,n_l;//,n_u2;
  double nu_trans;
  double tau_sobolev,beta;
  double R;

  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int lower = mastate[tid].level;

  int modelgridindex = cell[pkt_ptr->where].modelgridindex;

  #ifdef DEBUG_ON
    if (upper <= lower)
    {
      printout("[fatal] rad_rates_up: tried to calculate downward transition ... abort");
      abort();
    }
  #endif

  nu_trans = epsilon_trans/H;
  //A_ul = einstein_spontaneous_emission(element,ion,upper,lower);
  A_ul = einstein_spontaneous_emission(lineindex);
  B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans,3) * A_ul;
  B_lu = statweight_target/mastate[tid].statweight * B_ul;
  //double g_ratio = statweight_target/mastate[tid].statweight;
  //B_lu = g_ratio * B_ul;

  n_u = get_levelpop(element,ion,upper);
  n_l = get_levelpop(element,ion,lower);
  //double T_R = cell[pkt_ptr->where].T_R;
  //double W = cell[pkt_ptr->where].W;
  //n_u = n_l * W * g_ratio * exp(-epsilon_trans/KB/T_R);
  tau_sobolev = (B_lu*n_l - B_ul*n_u) * HCLIGHTOVERFOURPI * t_current;

  if (tau_sobolev <= 0)
  {
    if (SILENT == 0) printout("[warning] rad_excitation: tau_sobolev %g <= 0, set beta=1\n",tau_sobolev);
    if (SILENT == 0) printout("[warning] rad_excitation: element %d, ion %d, upper %d, lower %d\n",element,ion,upper,lower);
    if (SILENT == 0) printout("[warning] rad_excitation: n_l %g, n_u %g, B_lu %g, B_ul %g\n",n_l,n_u,B_lu,B_ul);
    if (SILENT == 0) printout("[warning] rad_excitation: T_e %g, T_R %g, W %g\n",get_Te(modelgridindex),get_TR(modelgridindex),get_W(modelgridindex));
    if (SILENT == 0) printout("[warning] rad_excitation: pkt_ptr->number %d\n",pkt_ptr->number);
    beta = 1.0;
    R = 0.;

    //printout("[fatal] rad_excitation: tau_sobolev <= 0 ... %g abort",tau_sobolev);
    //abort();
  }
  else
  {
    beta = 1.0/tau_sobolev * (1. - exp(-tau_sobolev));
    //printout("[check] rad_excitation: %g, %g, %g\n",1.0/tau_sobolev,exp(-tau_sobolev),1.0/tau_sobolev * (1. - exp(-tau_sobolev)));
    //n_u2 = calculate_levelpop_fromreflevel(pkt_ptr->where,element,ion,upper,lower,mastate[tid].nnlevel);
    //R = (B_lu*mastate[tid].nnlevel - B_ul*n_u2) * beta * radfield(nu_trans,pkt_ptr->where);
    R = mastate[tid].nnlevel * (B_lu - B_ul*n_u/n_l) * beta * radfield(nu_trans,modelgridindex);
  }

  #ifdef DEBUG_ON
    //printout("tau_sobolev, beta: %g, %g\n",tau_sobolev,beta);
    if (debuglevel == 2) printout("[debug] rad_rates_up: element, ion, upper, lower, A_ul, n_u: %d, %d, %d, %d, %g, %g\n",element,ion,upper,lower,A_ul,n_l);
    if (debuglevel == 2) printout("[debug] rad_exc: A_ul %g, tau_sobolev %g, n_u %g, n_l %g, radfield %g\n",A_ul,tau_sobolev,n_u,n_l,radfield(nu_trans,modelgridindex));
    if (debuglevel == 777) printout("[debug] rad_exc: A_ul %g, tau_sobolev %g, n_u %g, n_l %g, radfield %g\n",A_ul,tau_sobolev,n_u,n_l,radfield(nu_trans,modelgridindex));
    if (!isfinite(R))
    {
      printout("[fatal] rad_excitation: abort\n");
      printout("[fatal] rad_excitation: R %g, mastate[tid].nnlevel %g, B_lu %g, B_ul %g, n_u %g, n_l %g, beta %g, radfield %g,tau_sobolev %g, t_current %g\n",R,mastate[tid].nnlevel,B_lu,B_ul,n_u,n_l,beta,radfield(nu_trans,modelgridindex),tau_sobolev,t_current);
      printout("[fatal] rad_excitation: %g, %g, %g\n",1.0/tau_sobolev,exp(-tau_sobolev),1.0/tau_sobolev * (1. - exp(-tau_sobolev)));
      abort();
    }
  #endif

  return R;
}


///****************************************************************************
double rad_recombination(int modelgridindex, int lower, double epsilon_trans)
///radiative recombination rate: paperII 3.5.2
{
  double get_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
  int get_nphixstargets(int element, int ion, int level);
  int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
  double R;

  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int upper = mastate[tid].level;
  int phixstargetindex;
  double nne = get_nne(modelgridindex);

  R = 0.0;
  for (phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion-1,lower); phixstargetindex++)
  {
    if (get_phixsupperlevel(element,ion-1,lower,phixstargetindex) == upper)
    {
      R = mastate[tid].nnlevel * nne * get_spontrecombcoeff(element,ion-1,lower,phixstargetindex,modelgridindex);// + stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+(ion-1)]);
      //printout("calculate rad_recombination: element %d, ion %d, upper %d, -> lower %d, n_u %g, nne %g, spontrecombcoeff %g\n",element,ion,upper,lower,mastate[tid].nnlevel,nne,get_spontrecombcoeff(element, ion-1, lower, pkt_ptr->where));
      break;
    }
  }

  #ifdef DEBUG_ON
    //printout("[debug]    rad_recombiantion: R %g\n",R);
    if (!isfinite(R))
    {
      printout("fatal a2: abort\n");
      abort();
    }
  #endif

  return R;
}


///***************************************************************************/
double photoionization(int modelgridindex, int phixstargetindex, double epsilon_trans)
///photoionization rate: paperII 3.5.2
/// n_1 - occupation number of ground state
{
  double get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
  int get_nphixstargets(int element, int ion, int level);
  int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
  double R;

  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int lower = mastate[tid].level;

  #ifdef DEBUG_ON
    if (phixstargetindex > get_nphixstargets(element,ion,lower))
    {
      printout("[fatal] photoionization called with phixstargetindex %g > nphixstargets %g",phixstargetindex,get_nphixstargets(element,ion,lower));
      abort();
    }
  #endif

  R = mastate[tid].nnlevel * get_corrphotoioncoeff(element,ion,lower,phixstargetindex,modelgridindex);
  //R = get_corrphotoionrate(element,ion,lower,pkt_ptr->where);

  #ifdef DEBUG_ON
    //printout("[photoionization] R %g\n",R);
    if (!isfinite(R))
    {
      printout("fatal a4: abort\n");
      abort();
    }
  #endif

  return R;
}


/// Calculation of collisional rates /////////////////////////////////////////////////////
///***************************************************************************/
double col_excitation(int modelgridindex, int upper, int lineindex, double epsilon_trans)
/// collisional excitation rate: paperII 3.5.1
{
  //double osc_strength(int element, int ion, int upper, int lower);
  double osc_strength(int lineindex);
  double coll_str(int lineindex);
  double statw_up(int lineindex);
  double statw_down(int lineindex);
  double fac1;
  double nne,T_e;
  double C;
  double g_bar,test,Gamma;

  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int lower = mastate[tid].level;
  double n_l = mastate[tid].nnlevel;

  #ifdef DEBUG_ON
    if (upper <= lower)
    {
      printout("[fatal] col_excitation: tried to calculate downward transition ... abort");
      abort();
    }
  #endif

  T_e = get_Te(modelgridindex);
  fac1 = epsilon_trans/KB/T_e;
  nne = get_nne(modelgridindex);

  if (coll_str(lineindex) > 0.0) //positive values are treated as effective collision strengths
  {
    //from Osterbrock and Ferland, p51

    C = n_l * nne * 8.629e-6 * pow(T_e,-0.5) * coll_str(lineindex) * exp(-fac1) / statw_down(lineindex);
    //test test
    //C = n_l * nne * 8.629e-6 * pow(T_e,-0.5) * 0.01 * exp(-fac1) * statw_up(lineindex);
  }
  else if ((coll_str(lineindex) < 0) && (coll_str(lineindex) > -1.5)) //i.e. to catch -1
  {
    ///permitted E1 electric dipole transitions
    ///collisional excitation: formula valid only for atoms!!!!!!!!!!!
    ///Rutten script eq. 3.32. p.50
    //C = n_l * 2.16 * pow(fac1,-1.68) * pow(T_e,-1.5) * exp(-fac1) * nne * osc_strength(element,ion,upper,lower);

    ///Van-Regemorter formula, Mihalas (1978), eq.5-75, p.133
    g_bar = 0.2; ///this should be read in from transitions data: it is 0.2 for transitions nl -> n'l' and 0.7 for transitions nl -> nl'
    //test = 0.276 * exp(fac1) * gsl_sf_expint_E1(fac1);
    /// crude approximation to the already crude Van-Regemorter formula
    test = 0.276 * exp(fac1) * (-0.5772156649 - log(fac1));
    if (g_bar >= test)
      Gamma = g_bar;
    else
      Gamma = test;
    //C = n_l * C_0 * nne * pow(T_e,0.5) * 14.5*osc_strength(element,ion,upper,lower)*pow(H_ionpot/epsilon_trans,2) * fac1 * exp(-fac1) * Gamma;
    C = n_l * C_0 * nne * pow(T_e,0.5) * 14.5*osc_strength(lineindex)*pow(H_ionpot/epsilon_trans,2) * fac1 * exp(-fac1) * Gamma;
    //C = 0.0; //TESTING ONLY DELETE THIS
  }
  else if (coll_str(lineindex) > -3.5) //to catch -2 or -3
  {
    //forbidden transitions: magnetic dipole, electric quadropole...
    C = n_l * nne * 8.629e-6 * pow(T_e,-0.5) * 0.01 * exp(-fac1) * statw_up(lineindex);
  }
  else
  {
    C = 0.0;
  }

  #ifdef DEBUG_ON
    if (debuglevel == 2)
    {
      printout("[debug] col_exc: element %d, ion %d, lower %d, upper %d\n",element,ion,lower,upper);
      printout("[debug] col_exc: n_l %g, nne %g, T_e %g, f_ul %g, epsilon_trans %g, Gamma %g\n",n_l, nne,T_e,osc_strength(lineindex),epsilon_trans,Gamma);
    }
    if (debuglevel == 777)
      printout("[debug] col_exc: n_l %g, nne %g, T_e %g, f_ul %g, epsilon_trans %g, Gamma %g\n",n_l, nne,T_e,osc_strength(lineindex),epsilon_trans,Gamma);
    if (!isfinite(C))
	{
	  printout("fatal a5: abort\n");
	  printout("[debug] col_exc: element %d, ion %d, lower %d, upper %d\n",element,ion,lower,upper);
	  printout("[debug] col_exc: n_l %g, nne %g, T_e %g, f_ul %g, epsilon_trans %g, Gamma %g\n",n_l, nne,T_e,osc_strength(lineindex),epsilon_trans,Gamma);
	  printout("[debug] col_exc: g_bar %g, fac1 %g, test %g, %g, %g, %g\n",g_bar,fac1,test,0.276 * exp(fac1),-0.5772156649 - log(fac1),0.276 * exp(fac1) * (-0.5772156649 - log(fac1)));
	  printout("[debug] col_exc: coll_str(lineindex) %g statw_up(lineindex) %g mastate[tid].statweight %g\n", coll_str(lineindex),statw_up(lineindex),mastate[tid].statweight);
	  abort();
	}
  #endif
  return C;
}


///***************************************************************************/
double col_ionization(int modelgridindex, int phixstargetindex, double epsilon_trans)
/// collisional ionization rate: paperII 3.5.1
{
  int get_ionstage(int element, int ion);
  int get_nphixstargets(int element, int ion, int level);
  float get_phixsprobability(int element, int ion, int level, int phixstargetindex);
  double epsilon(int element, int ion, int level);
  double photoionization_crosssection(double nu_edge, double nu);
  double nu_lower;
  double fac1,sigma_bf;
  double nne,T_e;
  double g,C;

  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int lower = mastate[tid].level;
  double n_l = mastate[tid].nnlevel;

  int ionstage;

  if (phixstargetindex > get_nphixstargets(element,ion,lower))
  {
    printout("[fatal] col_ionization called with phixstargetindex %g > nphixstargets %g",phixstargetindex,get_nphixstargets(element,ion,lower));
    abort();
  }

  T_e = get_Te(modelgridindex);
  nne = get_nne(modelgridindex);

  nu_lower = epsilon_trans/H;
  fac1 = epsilon_trans/KB/T_e;

  ///Seaton approximation: Mihalas (1978), eq.5-79, p.134
  ///select gaunt factor according to ionic charge
  ionstage = get_ionstage(element,ion);
  if (ionstage == 1)
    g = 0.1;
  else if (ionstage == 2)
    g = 0.2;
  else
    g = 0.3;

  sigma_bf = photoionization_crosssection(nu_lower, nu_lower) * get_phixsprobability(element,ion,lower,phixstargetindex);
  C = n_l*nne * 1.55e13 * pow(T_e,-0.5) * g * sigma_bf * exp(-fac1)/fac1; ///photoionization at the edge

  #ifdef DEBUG_ON
    if (debuglevel == 777)
    printout("[debug] col_ion: n_l %g, nne %g, T_e %g, g %g, epsilon_trans %g, sigma_bf %g\n",n_l, nne,T_e,g,epsilon_trans,photoionization_crosssection(nu_lower, nu_lower));
    if (!isfinite(C))
    {
      printout("fatal a6: abort\n");
      abort();
    }
  #endif

  return C;
}


///***************************************************************************/
double col_deexcitation(int modelgridindex, int lower, double epsilon_trans, double statweight_target, int lineindex)
/// collisional deexcitation rate: paperII 3.5.1
{
  //double osc_strength(int element, int ion, int upper, int lower);
  double osc_strength(int lineindex);
  double coll_str(int lineindex);
  double epsilon(int element, int ion, int level);
  double fac1;
  double nne,T_e;
  double C;
  double g_bar,test,Gamma,g_ratio;

  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int upper = mastate[tid].level;
  double n_u = mastate[tid].nnlevel;

  #ifdef DEBUG_ON
    if (upper <= lower)
    {
      printout("[fatal] col_deexcitation: tried to calculate upward transition ... abort");
      abort();
    }
  #endif

  T_e = get_Te(modelgridindex);
  fac1 = epsilon_trans/KB/T_e;
  nne = get_nne(modelgridindex);

  if (coll_str(lineindex) > 0.0) //positive values are treated as effective collision strengths
  {
    //from Osterbrock and Ferland, p51
    //mastate[tid].statweight is UPPER LEVEL stat weight
    //statweight_target is LOWER LEVEL stat weight
    C = n_u * nne * 8.629e-6 * pow(T_e,-0.5) * coll_str(lineindex) / mastate[tid].statweight;
    // test test
    //C = n_u * nne * 8.629e-6 * pow(T_e,-0.5) * 0.01 * statweight_target;
  }
  else if ((coll_str(lineindex) < 0) && (coll_str(lineindex) > -1.5)) //i.e. to catch -1
  {
    ///permitted E1 electric dipole transitions
    ///collisional deexcitation: formula valid only for atoms!!!!!!!!!!!
    ///Rutten script eq. 3.33. p.50
    //f = osc_strength(element,ion,upper,lower);
    //C = n_u * 2.16 * pow(fac1,-1.68) * pow(T_e,-1.5) * stat_weight(element,ion,lower)/stat_weight(element,ion,upper)  * nne * f;

    ///Van-Regemorter formula, Mihalas (1978), eq.5-75, p.133
    g_bar = 0.2; ///this should be read in from transitions data: it is 0.2 for transitions nl -> n'l' and 0.7 for transitions nl -> nl'
    //test = 0.276 * exp(fac1) * gsl_sf_expint_E1(fac1);
    /// crude approximation to the already crude Van-Regemorter formula
    test = 0.276 * exp(fac1) * (-0.5772156649 - log(fac1));
    if (g_bar >= test)
      Gamma = g_bar;
    else
      Gamma = test;
    g_ratio = statweight_target/mastate[tid].statweight;
    //C = n_u * C_0 * nne * pow(T_e,0.5) * 14.5*osc_strength(element,ion,upper,lower)*pow(H_ionpot/epsilon_trans,2) * fac1 * g_ratio * Gamma;
    C = n_u * C_0 * nne * pow(T_e,0.5) * 14.5*osc_strength(lineindex)*pow(H_ionpot/epsilon_trans,2) * fac1 * g_ratio * Gamma;
    //C = 0.0; //TESTING ONLY DELETE THIS
  }
  else if (coll_str(lineindex) > -3.5) //to catch -2 or -3
  {
    //forbidden transitions: magnetic dipole, electric quadropole...
    C = n_u * nne * 8.629e-6 * pow(T_e,-0.5) * 0.01 * statweight_target;
  }
  else
  {
    C = 0.0;
  }

  #ifdef DEBUG_ON
    if (debuglevel == 2)
    {
      printout("[debug] col_deexc: element %d, ion %d, upper %d, lower %d\n",element,ion,upper,lower);
      printout("[debug] col_deexc: n_u %g, nne %g, T_e %g, f_ul %g, epsilon_trans %g, Gamma %g, g_ratio %g\n",n_u, nne,T_e,osc_strength(lineindex),epsilon_trans,Gamma,g_ratio);
    }
    if (debuglevel == 777)
    printout("[debug] col_deexc: n_u %g, nne %g, T_e %g, f_ul %g, epsilon_trans %g, Gamma %g, g_ratio %g\n",n_u, nne,T_e,osc_strength(lineindex),epsilon_trans,Gamma,g_ratio);
    //printout("col_deexc(%d,%d,%d,%d) %g\n",element,ion,upper,lower,C);
    if (!isfinite(C))
    {
      printout("fatal a7: abort\n");
      abort();
    }
  #endif

  return C;
}


///***************************************************************************/
double col_recombination(int modelgridindex, int lower, double epsilon_trans)
/// collisional recombination rate: paperII 3.5.1
{
  double get_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold);
  double photoionization_crosssection(double nu_edge, double nu);
  int get_nphixstargets(int element, int ion, int level);
  int get_phixsupperlevel(int element, int ion, int level, int phixstargetindex);
  float get_phixsprobability(int element, int ion, int level, int phixstargetindex);
  int get_ionstage(int element, int ion);
  double nu_lower;
  double fac1;
  double nne,T_e;
  double sigma_bf,g,C;
  int phixstargetindex;

  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int upper = mastate[tid].level;
  double n_u = mastate[tid].nnlevel;
  int ionstage;

  nu_lower = epsilon_trans/H;
  T_e = get_Te(modelgridindex);
  fac1 = epsilon_trans/KB/T_e;
  nne = get_nne(modelgridindex);

  C = 0.0;
  for (phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion-1,lower); phixstargetindex++)
  {
    if (get_phixsupperlevel(element,ion-1,lower,phixstargetindex) == upper)
    {
      ///Seaton approximation: Mihalas (1978), eq.5-79, p.134
      ///select gaunt factor according to ionic charge
      ionstage = get_ionstage(element,ion);
      #ifdef DEBUG_ON
        if (debuglevel == 777) printout("ionstage %d\n",ionstage);
      #endif
      if (ionstage-1 == 1)
        g = 0.1;
      else if (ionstage-1 == 2)
        g = 0.2;
      else
        g = 0.3;
      mastate[tid].ion = ion-1;      ///the same for mastate[tid].ion
      mastate[tid].level = lower;    ///set mastate[tid].level the lower level for photoionization_crosssection
      sigma_bf = photoionization_crosssection(nu_lower, nu_lower) * get_phixsprobability(element,ion-1,lower,phixstargetindex);
      #ifdef DEBUG_ON
        if (debuglevel == 777) printout("sigma_bf %g\n",sigma_bf);
      #endif
      mastate[tid].ion = ion;
      mastate[tid].level = upper;    ///restore the old values of pkt_ptr
      C = n_u*nne*nne * get_sahafact(element,ion-1,lower,phixstargetindex,T_e,epsilon_trans) * 1.55e13 * pow(T_e,-0.5) * g * sigma_bf * exp(-fac1)/fac1;

      #ifdef DEBUG_ON
        if (debuglevel == 777)
        {
          printout("get_sahafact %g, fac1 %g, C %g\n",get_sahafact(element,ion-1,lower,phixstargetindex,T_e,epsilon_trans),fac1,C);
          ///means n_u*nne * detailed_balancing of c_ikappa
          printout("[debug] col_recomb: n_u %g, nne %g, T_e %g, g %g, epsilon_trans %g, sigma_bf %g\n",n_u, nne,T_e,g,epsilon_trans,sigma_bf);
        }
        if (!isfinite(C))
        {
          printout("fatal a8: abort\n");
          abort();
        }
      #endif
      break;
    }
  }
  return C;
}


///***************************************************************************/
double radfield(double nu, int modelgridindex)
/// calculates ambient radiation field, which is parameterised as a diluted black body
{
  double B;
  double T_R,W;

  T_R = get_TR(modelgridindex);
  W   = get_W(modelgridindex);

  B = W * TWOHOVERCLIGHTSQUARED * pow(nu,3) * 1.0/(exp(HOVERKB*nu/T_R) - 1);
  return B;
}

///***************************************************************************/
double radfield2(double nu, double T, double W)
/// calculates ambient radiation field, which is parameterised as a diluted black body
{
  double B;

  B = W * TWOHOVERCLIGHTSQUARED * pow(nu,3) * 1.0/(exp(HOVERKB*nu/T) - 1);
  return B;
}


///****************************************************************************
double get_individ_rad_deexc(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc[i];
}

///***************************************************************************/
double get_individ_internal_down_same(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[i];
}

///***************************************************************************/
double get_individ_internal_up_same(int i)
{
  int element = mastate[tid].element;
  int ion = mastate[tid].ion;
  int level = mastate[tid].level;

  return cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same[i];
}
