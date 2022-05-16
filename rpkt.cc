#include "sn3d.h"
#include "atomic.h"
#include "boundary.h"
#include "grey_emissivities.h"
#include "grid.h"
#include "kpkt.h"
#include "ltepop.h"
#include "macroatom.h"
#include "polarization.h"
#include "radfield.h"
#include "rpkt.h"
#include "stats.h"
#include "update_grid.h"
#include "vectors.h"
#include "vpkt.h"

// Material for handing r-packet propagation.

__host__ __device__
int closest_transition(const double nu_cmf, const int next_trans)
/// for the propagation through non empty cells
// find the next transition lineindex redder than nu_cmf
// return -1 if no transition can be reached
{
  int match;

  int left = next_trans;
  int right = globals::nlines - 1;
  int middle = 1;

  /// if nu_cmf is smaller than the lowest frequency in the linelist,
  /// no line interaction is possible: return negative value as a flag
  if (nu_cmf < globals::linelist[right].nu)
  {
    return -1;
  }
  if (left > right)
  {
    //printout("[debug] pp should have no line interaction anymore\n");
    return -1;
  }

  if (left > 0)
  {
    /// if left = pkt_ptr->next_trans > 0 we know the next line we should interact with, independent of the packets
    /// current nu_cmf which might be smaller than globals::linelist[left].nu due to propagation errors
    match = left;
  }
  else if (nu_cmf >= globals::linelist[0].nu)
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
    while (left <= right)  ///must be a "<=" to obtain proper search results!!!
                          ///access to negative array indices is prevented by the upper check
    {
      middle = left + ((right - left) / 2);

      //printout("[debug] middle %d, left %d, right %d, nlines %d\n",middle,left,right,nlines);
      //printout("[debug] globals::linelist[middle].nu %g, globals::linelist[middle-1].nu %g\n",linelist[middle].nu,linelist[middle-1].nu);
      if (nu_cmf >= globals::linelist[middle].nu && nu_cmf < globals::linelist[middle - 1].nu)
      {
        break;
      }
      else if (nu_cmf >= globals::linelist[middle].nu)
      {
        right = middle - 1;
      }
      else
      {
        left = middle + 1;
      }
    }
    match = middle;
  }

  /// return the transitions frequency
  return match;
}


__host__ __device__
static double get_event(
  const int modelgridindex,
  PKT *pkt_ptr,             // pointer to packet object
  const double dopplerfactor,
  int *rpkt_eventtype,
  const double tau_rnd,     // random optical depth until which the packet travels
  const double abort_dist   // maximal travel distance before packet leaves cell or time step ends
)
// returns edist, the distance to the next physical event (continuum or bound-bound)
// BE AWARE THAT THIS PROCEDURE SHOULD BE ONLY CALLED FOR NON EMPTY CELLS!!
{
  /// initialize loop variables
  double tau = 0.;        ///initial optical depth along path
  double dist = 0.;       ///initial position on path
  double edist = 0.;

  PKT dummypkt = *pkt_ptr;
  PKT *dummypkt_ptr = &dummypkt;
  bool endloop = false;
  calculate_kappa_rpkt_cont(pkt_ptr, &globals::kappa_rpkt_cont[tid]);
  const double kap_cont = globals::kappa_rpkt_cont[tid].total * dopplerfactor;
  while (!endloop)
  {
    /// calculate distance to next line encounter ldist
    /// first select the closest transition in frequency
    /// we need its frequency nu_trans, the element/ion and the corresponding levels
    /// create therefore new variables in packet, which contain next_lowerlevel, ...
    const int lineindex = closest_transition(dummypkt_ptr->nu_cmf, dummypkt_ptr->next_trans);  ///returns negative value if nu_cmf > nu_trans

    if (lineindex >= 0)
    {
      /// line interaction in principle possible (nu_cmf > nu_trans)
      // printout("[debug] get_event:   line interaction possible\n");

      const double nu_trans = globals::linelist[lineindex].nu;

      // helper variable to overcome numerical problems after line scattering
      // further scattering events should be located at lower frequencies to prevent
      // multiple scattering events of one pp in a single line
      dummypkt_ptr->next_trans = lineindex + 1;

      const int element = globals::linelist[lineindex].elementindex;
      const int ion = globals::linelist[lineindex].ionindex;
      const int upper = globals::linelist[lineindex].upperlevelindex;
      const int lower = globals::linelist[lineindex].lowerlevelindex;

      double ldist;  // distance from current position to the line interaction
      if (dummypkt_ptr->nu_cmf < nu_trans)
      {
        //printout("dummypkt_ptr->nu_cmf %g < nu_trans %g, next_trans %d, element %d, ion %d, lower%d, upper %d\n",dummypkt_ptr->nu_cmf,nu_trans,dummypkt_ptr->next_trans,element,ion,lower,upper);
        ldist = 0;  /// photon was propagated too far, make sure that we don't miss a line
      }
      else
      {
        ldist = CLIGHT * dummypkt_ptr->prop_time * (dummypkt_ptr->nu_cmf / nu_trans - 1);
      }
      //fprintf(ldist_file,"%25.16e %25.16e\n",dummypkt_ptr->nu_cmf,ldist);
      if (ldist < 0.) printout("[warning] get_event: ldist < 0 %g\n",ldist);

      // printout("[debug] get_event:     ldist %g\n",ldist);

      const double tau_cont = kap_cont * ldist;

      // printout("[debug] get_event:     tau_rnd %g, tau %g, tau_cont %g\n", tau_rnd, tau, tau_cont);

      if (tau_rnd - tau > tau_cont)
      {
        const double A_ul = einstein_spontaneous_emission(lineindex);
        const double B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans, 3) * A_ul;
        const double B_lu = stat_weight(element, ion, upper) / stat_weight(element, ion, lower) * B_ul;

        const double n_u = get_levelpop(modelgridindex, element, ion, upper);
        const double n_l = get_levelpop(modelgridindex, element, ion, lower);

        double tau_line = (B_lu * n_l - B_ul * n_u) * HCLIGHTOVERFOURPI * dummypkt_ptr->prop_time;

        if (tau_line < 0)
        {
          //printout("[warning] get_event: tau_line %g < 0, n_l %g, n_u %g, B_lu %g, B_ul %g, W %g, T_R %g, element %d, ion %d, upper %d, lower %d ... abort\n",tau_line, n_l,n_u,B_lu,B_ul,get_W(grid::get_cell_modelgridindex(pkt_ptr->where)),get_TR(grid::get_cell_modelgridindex(pkt_ptr->where)),element,ion,upper,lower);
          //printout("[warning] get_event: set tau_line = 0\n");
          tau_line = 0.;
          //printout("[fatal] get_event: tau_line < 0 ... abort\n");
          //abort();
        }

        if (false)
        {
          printout("[debug] get_event:     tau_line %g\n", tau_line);
          printout("[debug] get_event:       tau_rnd - tau > tau_cont\n");
        }

        if (tau_rnd - tau > tau_cont + tau_line)
        {
          // total optical depth still below tau_rnd: propagate to the line and continue

          if (false)
          {
            printout("[debug] get_event:         tau_rnd - tau > tau_cont + tau_line ... proceed this packets propagation\n");
            printout("[debug] get_event:         dist %g, abort_dist %g, dist-abort_dist %g\n", dist, abort_dist, dist-abort_dist);
          }

          dist = dist + ldist;
          if (dist > abort_dist)
          {
            dummypkt_ptr->next_trans -= 1;
            pkt_ptr->next_trans = dummypkt_ptr->next_trans;
            // printout("[debug] get_event:         leave propagation loop (dist %g > abort_dist %g) ... dummypkt_ptr->next_trans %d\n", dist, abort_dist, dummypkt_ptr->next_trans);
            return abort_dist + 1e20;
          }

          tau += tau_cont + tau_line;
          //dummypkt_ptr->next_trans += 1;
          dummypkt_ptr->prop_time += ldist / globals::CLIGHT_PROP;
          move_pkt(dummypkt_ptr, ldist, dummypkt_ptr->prop_time);
          radfield::increment_lineestimator(modelgridindex, lineindex, dummypkt_ptr->prop_time * CLIGHT * dummypkt_ptr->e_cmf / dummypkt_ptr->nu_cmf);

          if (false)
          {
            const int next_trans = dummypkt_ptr->next_trans;
            printout("[debug] get_event:         dummypkt_ptr->nu_cmf %g, nu(dummypkt_ptr->next_trans=%d) %g, nu(dummypkt_ptr->next_trans-1=%d) %g\n", dummypkt_ptr->nu_cmf, next_trans, globals::linelist[next_trans].nu, next_trans-1, globals::linelist[next_trans-1].nu);
            printout("[debug] get_event:         (dummypkt_ptr->nu_cmf - nu(dummypkt_ptr->next_trans-1))/dummypkt_ptr->nu_cmf %g\n", (dummypkt_ptr->nu_cmf- globals::linelist[next_trans-1].nu)/dummypkt_ptr->nu_cmf);

            if (dummypkt_ptr->nu_cmf >= globals::linelist[next_trans].nu && dummypkt_ptr->nu_cmf < globals::linelist[next_trans-1].nu)
              printout("[debug] get_event:           nu(next_trans-1) > nu_cmf >= nu(next_trans)\n");
            else if (dummypkt_ptr->nu_cmf < globals::linelist[next_trans].nu)
              printout("[debug] get_event:           nu_cmf < nu(next_trans)\n");
            else
              printout("[debug] get_event:           nu_cmf >= nu(next_trans-1)\n");
          }
        }
        else
        {
          /// bound-bound process occurs
          // printout("[debug] get_event:         tau_rnd - tau <= tau_cont + tau_line: bb-process occurs\n");

          pkt_ptr->mastate.element = element;
          pkt_ptr->mastate.ion = ion;
          pkt_ptr->mastate.level = upper;  ///if the MA will be activated it must be in the transitions upper level
          pkt_ptr->mastate.activatingline = lineindex;

          edist = dist + ldist;
          if (edist > abort_dist)
          {
            dummypkt_ptr->next_trans = dummypkt_ptr->next_trans - 1;
          }
          else if (DETAILED_LINE_ESTIMATORS_ON)
          {
            dummypkt_ptr->prop_time += ldist / globals::CLIGHT_PROP;
            move_pkt(dummypkt_ptr, ldist, dummypkt_ptr->prop_time);
            radfield::increment_lineestimator(modelgridindex, lineindex, dummypkt_ptr->prop_time * CLIGHT * dummypkt_ptr->e_cmf / dummypkt_ptr->nu_cmf);
          }

          *rpkt_eventtype = RPKT_EVENTTYPE_BB;
          /// the line and its parameters were already selected by closest_transition!
          endloop = true;
          // printout("[debug] get_event:         edist %g, abort_dist %g, edist-abort_dist %g, endloop   %d\n",edist,abort_dist,edist-abort_dist,endloop);
        }
      }
      else
      {
        /// continuum process occurs
        edist = dist + (tau_rnd - tau) / kap_cont;
        // assert_always((tau_rnd - tau) / kap_cont < ldist);
        dummypkt_ptr->next_trans -= 1;
        // printout("[debug] get_event:        distance to the occuring continuum event %g, abort_dist %g\n", edist, abort_dist);

        *rpkt_eventtype = RPKT_EVENTTYPE_CONT;
        endloop = true;
      }
    }
    else
    {
      dummypkt_ptr->next_trans = globals::nlines + 1;  ///helper variable to overcome numerical problems after line scattering
      /// no line interaction possible - check whether continuum process occurs in cell
      // printout("[debug] get_event:     line interaction impossible\n");
      const double tau_cont = kap_cont * (abort_dist - dist);
      //printout("nu_cmf %g, opticaldepths in ff %g, es %g\n",pkt_ptr->nu_cmf,kappa_rpkt_cont[tid].ff*(abort_dist-dist),kappa_rpkt_cont[tid].es*(abort_dist-dist));
      // printout("[debug] get_event:     tau_rnd %g, tau %g, tau_cont %g\n", tau_rnd, tau, tau_cont);

      if (tau_rnd - tau > tau_cont)
      {
        /// travel out of cell or time step
        // printout("[debug] get_event:       travel out of cell or time step\n");

        edist = abort_dist + 1e20;
        endloop = true;
      }
      else
      {
        /// continuum process occurs at edist
        edist = dist + (tau_rnd - tau) / kap_cont;
        // printout("[debug] get_event:       continuum process occurs at edist %g\n",edist);

        *rpkt_eventtype = RPKT_EVENTTYPE_CONT;
        endloop = true;
      }
    }
    //propagationcounter++;
  }

  pkt_ptr->next_trans = dummypkt_ptr->next_trans;
  assert_always(std::isfinite(edist));

  return edist;
}


// todo: place next to ratecoeff.cc get_index_from_cumulativesums()
static int upper_bound(double arr[], int arraylen, double searchval)
// return the index of arr such that arr[index] >= searchval, or arraylen-1
// if searchval > arr[arraylen-1]
{
  int low = 0;
  int high = arraylen;

  while (low < high)
  {
    // Find the middle index
    const int mid = low + (high - low) / 2;

    if (searchval >= arr[mid])
    {
      low = mid + 1;
    }
    else
    {
      high = mid;
    }
  }

  // if searchval is greater than arr[arraylen - 1]
  if (low < arraylen && arr[low] <= searchval)
  {
     low++;
  }

  return low;
}


__host__ __device__
static void rpkt_event_continuum(PKT *pkt_ptr, rpkt_cont_opacity_struct kappa_rpkt_cont_thisthread, const double dopplerfactor, int modelgridindex)
{
  const double nu = pkt_ptr->nu_cmf;

  const double kappa_cont = kappa_rpkt_cont_thisthread.total * dopplerfactor;
  const double sigma = kappa_rpkt_cont_thisthread.es * dopplerfactor;
  const double kappa_ff = kappa_rpkt_cont_thisthread.ff * dopplerfactor;
  const double kappa_bf = kappa_rpkt_cont_thisthread.bf * dopplerfactor;

  /// continuum process happens. select due to its probabilities sigma/kappa_cont, kappa_ff/kappa_cont, kappa_bf/kappa_cont
  const double zrand = gsl_rng_uniform(rng);
  // printout("[debug] rpkt_event:   r-pkt undergoes a continuum transition\n");
  // printout("[debug] rpkt_event:   zrand*kappa_cont %g, sigma %g, kappa_ff %g, kappa_bf %g\n", zrand * kappa_cont, sigma, kappa_ff, kappa_bf);

  if (zrand * kappa_cont < sigma)
  {
    /// electron scattering occurs
    /// in this case the packet stays a R_PKT of same nu_cmf than before (coherent scattering)
    /// but with different direction
    // printout("[debug] rpkt_event:   electron scattering\n");
    pkt_ptr->interactions += 1;
    pkt_ptr->nscatterings += 1;
    pkt_ptr->last_event = 12;
    stats::increment(stats::COUNTER_ESCOUNTER);

    // generate a virtual packet
    #ifdef VPKT_ON
      int realtype = 1;
      pkt_ptr->last_cross = NONE;
      vpkt_call_estimators(pkt_ptr, pkt_ptr->prop_time, realtype);
    #endif

    //pkt_ptr->nu_cmf = 3.7474058e+14;
    escat_rpkt(pkt_ptr);
    /// Electron scattering does not modify the last emission flag
    //pkt_ptr->emissiontype = get_continuumindex(element,ion-1,lower);
    /// but it updates the last emission position
    vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
    pkt_ptr->em_time = pkt_ptr->prop_time;

    /// Set some flags
    //pkt_ptr->next_trans = 0;   ///packet's comoving frame frequency is conserved during electron scattering
                                 ///don't touch the value of next_trans to save transition history
  }
  else if (zrand * kappa_cont < sigma + kappa_ff)
  {
    /// ff: transform to k-pkt
    // printout("[debug] rpkt_event:   free-free transition\n");
    stats::increment(stats::COUNTER_K_STAT_FROM_FF);
    pkt_ptr->interactions += 1;
    pkt_ptr->last_event = 5;
    pkt_ptr->type = TYPE_KPKT;
    pkt_ptr->absorptiontype = -1;
    #ifndef FORCE_LTE
      //kffabs[pkt_ptr->where] += pkt_ptr->e_cmf;
    #endif
  }
  else if (zrand * kappa_cont < sigma + kappa_ff + kappa_bf)
  {
    /// bf: transform to k-pkt or activate macroatom corresponding to probabilities
    // printout("[debug] rpkt_event:   bound-free transition\n");

    pkt_ptr->absorptiontype = -2;

    const double kappa_bf_inrest = kappa_rpkt_cont_thisthread.bf;

    /// Determine in which continuum the bf-absorption occurs
    const double zrand2 = gsl_rng_uniform(rng);
    const double kappa_bf_rand = zrand2 * kappa_bf_inrest;

    const int allcontindex = upper_bound(globals::phixslist[tid].kappa_bf_sum, globals::nbfcontinua, kappa_bf_rand);
    assert_always(globals::phixslist[tid].kappa_bf_sum[globals::nbfcontinua - 1] == kappa_bf_inrest);

    const double nu_edge = globals::allcont[allcontindex].nu_edge;
    const int element = globals::allcont[allcontindex].element;
    const int ion = globals::allcont[allcontindex].ion;
    const int level = globals::allcont[allcontindex].level;
    const int phixstargetindex = globals::allcont[allcontindex].phixstargetindex;

    // printout("[debug] rpkt_event:   bound-free: element %d, ion+1 %d, upper %d, ion %d, lower %d\n", element, ion + 1, 0, ion, level);
    // printout("[debug] rpkt_event:   bound-free: nu_edge %g, nu %g\n", nu_edge, nu);

    if (TRACK_ION_STATS)
    {
      stats::increment_ion_stats_contabsorption(pkt_ptr, modelgridindex, element, ion);
    }

    /// and decide whether we go to ionisation energy
    const double zrand3 = gsl_rng_uniform(rng);
    if (zrand3 < nu_edge / nu)
    {
      stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_BF);
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 3;

      if (TRACK_ION_STATS)
      {
        stats::increment_ion_stats(modelgridindex, element, ion + 1, stats::ION_MACROATOM_ENERGYIN_PHOTOION, pkt_ptr->e_cmf);
      }

      pkt_ptr->type = TYPE_MA;
      #ifndef FORCE_LTE
        //maabs[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif
      pkt_ptr->mastate.element = element;
      pkt_ptr->mastate.ion     = ion + 1;
      const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
      pkt_ptr->mastate.level   = upper;
      pkt_ptr->mastate.activatingline = -99;
    }
    /// or to the thermal pool
    else
    {
      /// transform to k-pkt
      // printout("[debug] rpkt_event:   bound-free: transform to k-pkt\n");
      stats::increment(stats::COUNTER_K_STAT_FROM_BF);
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 4;
      pkt_ptr->type = TYPE_KPKT;
    }
  }
  else
  {
    printout("ERROR: could not continuum process\n");
    abort();
  }
}


__host__ __device__
static void rpkt_event_boundbound(PKT *pkt_ptr, const int mgi)
{
  /// bound-bound transition occured
  /// activate macro-atom in corresponding upper-level. Actually all the information
  /// about the macro atoms state has already been set by closest_transition, so
  /// we need here just the activation!
  // printout("[debug] rpkt_event: bound-bound activation of macroatom\n");
  //if (tid == 0) ma_stat_activation_bb++;
  stats::increment(stats::COUNTER_MA_STAT_ACTIVATION_BB);
  pkt_ptr->interactions += 1;
  pkt_ptr->last_event = 1;

  pkt_ptr->absorptiontype = pkt_ptr->mastate.activatingline;
  pkt_ptr->absorptionfreq = pkt_ptr->nu_rf;//pkt_ptr->nu_cmf;
  pkt_ptr->absorptiondir[0] = pkt_ptr->dir[0];
  pkt_ptr->absorptiondir[1] = pkt_ptr->dir[1];
  pkt_ptr->absorptiondir[2] = pkt_ptr->dir[2];
  pkt_ptr->type = TYPE_MA;

  #if (TRACK_ION_STATS)
  const int element = pkt_ptr->mastate.element;
  const int ion = pkt_ptr->mastate.ion;
  stats::increment_ion_stats(mgi, element, ion, stats::ION_MACROATOM_ENERGYIN_RADEXC, pkt_ptr->e_cmf);

  const int et = pkt_ptr->emissiontype;
  if (et >= 0)
  {
    const int emissionelement = globals::linelist[et].elementindex;
    const int emissionion = globals::linelist[et].ionindex;
    stats::increment_ion_stats(mgi, emissionelement, emissionion, stats::ION_BOUNDBOUND_ABSORBED, pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf);
  }
  #endif

  #ifdef RECORD_LINESTAT
    if (tid == 0) globals::acounter[pkt_ptr->next_trans-1] += 1;  /// This way we will only record line statistics from OMP-thread 0
                                                         /// With an atomic pragma or a thread-private structure with subsequent
                                                         /// reduction this could be extended to all threads. However, I'm not
                                                         /// sure if this is worth the additional computational expenses.
  #endif
}


__host__ __device__
static void rpkt_event_thickcell(PKT *pkt_ptr)
/// Event handling for optically thick cells. Those cells are treated in a grey
/// approximation with electron scattering only.
/// The packet stays an R_PKT of same nu_cmf than before (coherent scattering)
/// but with different direction.
{
  // printout("[debug] rpkt_event_thickcell:   electron scattering\n");
  pkt_ptr->interactions += 1;
  pkt_ptr->nscatterings += 1;
  pkt_ptr->last_event = 12;
  stats::increment(stats::COUNTER_ESCOUNTER);

  emitt_rpkt(pkt_ptr);
  /// Electron scattering does not modify the last emission flag
  //pkt_ptr->emissiontype = get_continuumindex(element,ion-1,lower);
  /// but it updates the last emission position
  vec_copy(pkt_ptr->em_pos, pkt_ptr->pos);
  pkt_ptr->em_time = pkt_ptr->prop_time;
}


__host__ __device__
static double closest_transition_empty(PKT *pkt_ptr)
/// for the propagation through empty cells
/// here its possible that the packet jumps over several lines
{
  //int left = 0;
  int left = pkt_ptr->next_trans;
  //printout("[debug] closest_transition: initial left %d\n",left);
  int right = globals::nlines - 1;

  //printout("[debug] ___closest_transition___: initial left %d, right %d, nu_cmf %g\n",left,right,pkt_ptr->nu_cmf);
  //printout("[debug] ___closest_transition___: nu_left %g, nu_right%g\n",linelist[left].nu,linelist[right].nu);
  /// if nu_cmf is smaller than the lowest frequency in the linelist,
  /// no line interaction is possible: return negative value as a flag
  if (pkt_ptr->nu_cmf < globals::linelist[right].nu)
  {
    pkt_ptr->next_trans = globals::nlines + 1;  ///helper variable to overcome numerical problems after line scattering
    return -1;
  }
  if (left > right)
  {
    //printout("[debug] pp should have no line interaction anymore\n");
    pkt_ptr->next_trans = globals::nlines + 1;  ///helper variable to overcome numerical problems after line scattering
    return -1;
  }

  int match;
  /// no check for left > 0 in the empty case as it is possible that the packet is moved over
  /// several lines through the empty cell
  if (pkt_ptr->nu_cmf >= globals::linelist[left].nu)
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

    int middle = 1;
    while (left <= right)  // must be a "<=" to obtain proper search results!!!
                           // access to negative array indices is prevented by the upper check
    {
      middle = left + ((right-left) / 2);

      //printout("[debug] middle %d, left %d, right %d, nlines %d\n",middle,left,right,nlines);
      //printout("[debug] globals::linelist[middle].nu %g, globals::linelist[middle-1].nu %g\n",linelist[middle].nu,linelist[middle-1].nu);
      if (pkt_ptr->nu_cmf >= globals::linelist[middle].nu && pkt_ptr->nu_cmf < globals::linelist[middle-1].nu) break;

      if (pkt_ptr->nu_cmf >= globals::linelist[middle].nu)
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
  const double nu_trans = globals::linelist[match].nu;
  //pkt_ptr->mastate.element = globals::linelist[match].elementindex;
  //pkt_ptr->mastate.ion     = globals::linelist[match].ionindex;
  //pkt_ptr->mastate.level   = globals::linelist[match].upperlevelindex;  ///if the MA will be activated it must be in the transitions upper level
  //pkt_ptr->mastate.activatedfromlevel   = globals::linelist[match].lowerlevelindex;  ///helper variable for the transitions lower level

  /// For the empty case it's match not match+1: a line interaction is only possible in the next iteration
  /// of the propagation loop. We just have to make sure that the next "normal" line search knows about the
  /// current position of the photon in the frequency list.
  pkt_ptr->next_trans   = match;  /// helper variable to overcome numerical problems after line scattering
                                     ///further scattering events should be located at lower frequencies to prevent
                                     ///multiple scattering events of one pp in a single line

  /// return the transitions frequency
  return nu_trans;
}


__host__ __device__
static void update_estimators(PKT *pkt_ptr, const double distance)
/// Update the volume estimators J and nuJ
/// This is done in another routine than move, as we sometimes move dummy
/// packets which do not contribute to the radiation field.
{
  const int cellindex = pkt_ptr->where;
  const int modelgridindex = grid::get_cell_modelgridindex(cellindex);

  /// Update only non-empty cells
  if (modelgridindex == grid::get_npts_model())
  {
    return;
  }
  const double distance_e_cmf = distance * pkt_ptr->e_cmf;
  const double nu = pkt_ptr->nu_cmf;
  //double bf = exp(-HOVERKB*nu/globals::cell[modelgridindex].T_e);
  radfield::update_estimators(modelgridindex, distance_e_cmf, nu, pkt_ptr, pkt_ptr->prop_time);

#if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
  const int nelements = get_nelements();
  const int max_nions = get_max_nions();
#endif

  #ifndef FORCE_LTE
    ///ffheatingestimator does not depend on ion and element, so an array with gridsize is enough.
    ///quick and dirty solution: store info in element=ion=0, and leave the others untouched (i.e. zero)
    safeadd(globals::ffheatingestimator[modelgridindex], distance_e_cmf * globals::kappa_rpkt_cont[tid].ffheating);

    #if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
      #if (!NO_LUT_PHOTOION)
      const double distance_e_cmf_over_nu = distance_e_cmf / nu;
      #endif
      for (int i = 0; i < globals::nbfcontinua_ground; i++)
      {
        const double nu_edge = globals::groundcont[i].nu_edge;
        if (nu > nu_edge)
        {
          const int element = globals::groundcont[i].element;
          /// Cells with zero abundance for a specific element have zero contribution
          /// (set in calculate_kappa_rpkt_cont and therefore do not contribute to
          /// the estimators
          if (grid::get_elem_abundance(modelgridindex, element) > 0)
          {
            const int ion = globals::groundcont[i].ion;
            const int ionestimindex = modelgridindex * nelements * max_nions + element * max_nions + ion;
            #if (!NO_LUT_PHOTOION)
              safeadd(globals::gammaestimator[ionestimindex], globals::phixslist[tid].groundcont_gamma_contr[i] * distance_e_cmf_over_nu);

              if (!std::isfinite(globals::gammaestimator[ionestimindex]))
              {
                printout("[fatal] update_estimators: gamma estimator becomes non finite: mgi %d element %d ion %d gamma_contr %g, distance_e_cmf_over_nu %g\n", modelgridindex, element, ion, globals::phixslist[tid].groundcont_gamma_contr[i], distance_e_cmf_over_nu);
                abort();
              }
            #endif
            #if (!NO_LUT_BFHEATING)
              safeadd(globals::bfheatingestimator[ionestimindex], globals::phixslist[tid].groundcont_gamma_contr[i] * distance_e_cmf * (1. - nu_edge/nu));
            #endif
          }
        }
        else
          break; // because groundcont is sorted by nu_edge, nu < nu_edge for all remaining items
      }
    #endif

  #endif
}


__host__ __device__
static bool do_rpkt_step(PKT *pkt_ptr, const double t2)
// Routine for moving an r-packet. Similar to do_gamma in objective.
// return value - true if no mgi change, no pkttype change and not reached end of timestep, false otherwise
{
  const int cellindex = pkt_ptr->where;
  int mgi = grid::get_cell_modelgridindex(cellindex);
  const int oldmgi = mgi;

  if (pkt_ptr->next_trans > 0)
  {
    //if (globals::debuglevel == 2) printout("[debug] do_rpkt: init: pkt_ptr->nu_cmf %g, nu(pkt_ptr->next_trans=%d) %g, nu(pkt_ptr->next_trans-1=%d) %g, pkt_ptr->where %d\n", pkt_ptr->nu_cmf, pkt_ptr->next_trans, globals::linelist[pkt_ptr->next_trans].nu, pkt_ptr->next_trans-1, globals::linelist[pkt_ptr->next_trans-1].nu, pkt_ptr->where );
    if (globals::debuglevel == 2) printout("[debug] do_rpkt: init: (pkt_ptr->nu_cmf - nu(pkt_ptr->next_trans-1))/pkt_ptr->nu_cmf %g\n", (pkt_ptr->nu_cmf-globals::linelist[pkt_ptr->next_trans-1].nu)/pkt_ptr->nu_cmf);
  }

  // Assign optical depth to next physical event. And start counter of
  // optical depth for this path.
  const double zrand = gsl_rng_uniform_pos(rng);
  const double tau_next = -1. * log(zrand);

  // Start by finding the distance to the crossing of the grid cell
  // boundaries. sdist is the boundary distance and snext is the
  // grid cell into which we pass.
  int snext;
  double sdist = boundary_cross(pkt_ptr, pkt_ptr->prop_time, &snext);

  if (sdist == 0)
  {
    change_cell(pkt_ptr, snext, pkt_ptr->prop_time);
    const int cellindexnew = pkt_ptr->where;
    mgi = grid::get_cell_modelgridindex(cellindexnew);

    return (pkt_ptr->type == TYPE_RPKT && (mgi == grid::get_npts_model() || mgi == oldmgi));
  }
  else
  {
    const double maxsdist = (grid::grid_type == GRID_SPHERICAL1D) ? 2 * globals::rmax * (pkt_ptr->prop_time + sdist / globals::CLIGHT_PROP) / globals::tmin : globals::rmax * pkt_ptr->prop_time / globals::tmin;
    if (sdist > maxsdist)
    {
      printout("[fatal] do_rpkt: Unreasonably large sdist. Rpkt. Abort. %g %g %g\n", globals::rmax, pkt_ptr->prop_time/globals::tmin, sdist);
      abort();
    }

    if (sdist < 1)
    {
      const int cellindexnew = pkt_ptr->where;
      printout("[warning] r_pkt: Negative distance (sdist = %g). Abort.\n", sdist);
      printout("[warning] r_pkt: cell %d snext %d\n", cellindexnew, snext);
      printout("[warning] r_pkt: pos %g %g %g\n", pkt_ptr->pos[0], pkt_ptr->pos[1], pkt_ptr->pos[2]);
      printout("[warning] r_pkt: dir %g %g %g\n", pkt_ptr->dir[0], pkt_ptr->dir[1], pkt_ptr->dir[2]);
      printout("[warning] r_pkt: cell corner %g %g %g\n",
               grid::get_cellcoordmin(cellindexnew, 0) * pkt_ptr->prop_time / globals::tmin,
               grid::get_cellcoordmin(cellindexnew, 1) * pkt_ptr->prop_time / globals::tmin,
               grid::get_cellcoordmin(cellindexnew, 2) * pkt_ptr->prop_time / globals::tmin);
      printout("[warning] r_pkt: cell width %g\n",grid::wid_init(0)*pkt_ptr->prop_time/globals::tmin);
      //abort();
    }
    if (((snext != -99) && (snext < 0)) || (snext >= grid::ngrid))
    {
      printout("[fatal] r_pkt: Heading for inappropriate grid cell. Abort.\n");
      printout("[fatal] r_pkt: Current cell %d, target cell %d.\n", pkt_ptr->where, snext);
      abort();
    }

    if (sdist > globals::max_path_step)
    {
      sdist = globals::max_path_step;
      snext = pkt_ptr->where;
    }

    // At present there is no scattering/destruction process so all that needs to
    // happen is that we determine whether the packet reaches the boundary during the timestep.

    // Find how far it can travel during the time inverval.

    double tdist = (t2 - pkt_ptr->prop_time) * globals::CLIGHT_PROP;

    assert_always(tdist >= 0);

    const double dopplerfactor = doppler_packetpos(pkt_ptr);

    double edist;
    int rpkt_eventtype;
    bool find_nextline = false;
    if (mgi == grid::get_npts_model())
    {
      /// for empty cells no physical event occurs. The packets just propagate.
      edist = 1e99;
      find_nextline = true;
      // printout("[debug] do_rpkt: propagating through empty cell, set edist=1e99\n");
    }
    else if (grid::modelgrid[mgi].thick == 1)
    {
      /// In the case of optically thick cells, we treat the packets in grey approximation to speed up the calculation
      /// Get distance to the next physical event in this case only electron scattering
      //kappa = SIGMA_T*grid::get_nne(mgi);
      const double kappa = grid::get_kappagrey(mgi) * grid::get_rho(mgi) * dopplerfactor;
      const double tau_current = 0.0;
      edist = (tau_next - tau_current) / kappa;
      find_nextline = true;
      // printout("[debug] do_rpkt: propagating through grey cell, edist  %g\n",edist);
    }
    else
    {
      // get distance to the next physical event (continuum or bound-bound)
      edist = get_event(mgi, pkt_ptr, dopplerfactor, &rpkt_eventtype, tau_next, fmin(tdist, sdist)); //, kappacont_ptr, sigma_ptr, kappaff_ptr, kappabf_ptr);
      if (false)
      {
        const int next_trans = pkt_ptr->next_trans;
        printout("[debug] do_rpkt: after edist: pkt_ptr->nu_cmf %g, nu(pkt_ptr->next_trans=%d) %g\n", pkt_ptr->nu_cmf, next_trans, globals::linelist[next_trans].nu);
      }
    }
    assert_always(edist >= 0);
    //printout("[debug] do_rpkt: sdist, tdist, edist %g, %g, %g\n",sdist,tdist,edist);

    if ((sdist < tdist) && (sdist < edist))
    {
      // printout("[debug] do_rpkt: sdist < tdist && sdist < edist\n");
      // Move it into the new cell.
      pkt_ptr->prop_time += sdist / 2. / globals::CLIGHT_PROP;
      move_pkt(pkt_ptr, sdist / 2., pkt_ptr->prop_time);
      update_estimators(pkt_ptr, sdist);
      if (globals::do_rlc_est != 0 && globals::do_rlc_est != 3)
      {
        rlc_emiss_rpkt(pkt_ptr, sdist);
      }
      pkt_ptr->prop_time += sdist / 2. / globals::CLIGHT_PROP;
      move_pkt(pkt_ptr, sdist / 2., pkt_ptr->prop_time);

      if (snext != pkt_ptr->where)
      {
        change_cell(pkt_ptr, snext, pkt_ptr->prop_time);
        const int cellindexnew = pkt_ptr->where;
        mgi = grid::get_cell_modelgridindex(cellindexnew);
      }
      // New cell so reset the scat_counter
      pkt_ptr->scat_count = 0;
      pkt_ptr->last_event = pkt_ptr->last_event + 100;

      /// For empty or grey cells a photon can travel over several bb-lines. Thus we need to
      /// find the next possible line interaction.
      if (find_nextline)
      {
        /// However, this is only required if the new cell is non-empty or non-grey
        if (mgi != grid::get_npts_model() && grid::modelgrid[mgi].thick != 1)
        {
          closest_transition_empty(pkt_ptr);
        }
      }

      return (pkt_ptr->type == TYPE_RPKT && (mgi == grid::get_npts_model() || mgi == oldmgi));
    }
    else if ((tdist < sdist) && (tdist < edist))
    {
      // printout("[debug] do_rpkt: tdist < sdist && tdist < edist\n");
      // Doesn't reach boundary
      pkt_ptr->prop_time += tdist / 2. / globals::CLIGHT_PROP;
      move_pkt(pkt_ptr, tdist / 2., pkt_ptr->prop_time);
      update_estimators(pkt_ptr, tdist);
      if (globals::do_rlc_est != 0 && globals::do_rlc_est != 3)
      {
        rlc_emiss_rpkt(pkt_ptr, tdist);
      }
      pkt_ptr->prop_time = t2;
      move_pkt(pkt_ptr, tdist / 2., pkt_ptr->prop_time);
      pkt_ptr->last_event = pkt_ptr->last_event + 1000;

      /// For empty or grey cells a photon can travel over several bb-lines. Thus we need to
      /// find the next possible line interaction.
      if (find_nextline)
      {
        closest_transition_empty(pkt_ptr);
      }

      return false;
    }
    else if ((edist < sdist) && (edist < tdist))
    {
      // bound-bound or continuum event
      // printout("[debug] do_rpkt: edist < sdist && edist < tdist\n");
      pkt_ptr->prop_time += edist / 2. / globals::CLIGHT_PROP;
      move_pkt(pkt_ptr, edist / 2., pkt_ptr->prop_time);
      update_estimators(pkt_ptr, edist);
      if (globals::do_rlc_est != 0 && globals::do_rlc_est != 3)
      {
        rlc_emiss_rpkt(pkt_ptr, edist);
      }
      pkt_ptr->prop_time += edist / 2. / globals::CLIGHT_PROP;
      move_pkt(pkt_ptr, edist / 2., pkt_ptr->prop_time);

      // The previously selected and in pkt_ptr stored event occurs. Handling is done by rpkt_event
      if (grid::modelgrid[mgi].thick == 1)
      {
        rpkt_event_thickcell(pkt_ptr);
      }
      else if (rpkt_eventtype == RPKT_EVENTTYPE_BB)
      {
        rpkt_event_boundbound(pkt_ptr, mgi);
      }
      else if (rpkt_eventtype == RPKT_EVENTTYPE_CONT)
      {
        rpkt_event_continuum(pkt_ptr, globals::kappa_rpkt_cont[tid], dopplerfactor, mgi);
      }
      else
      {
        assert_always(false);
      }

      return (pkt_ptr->type == TYPE_RPKT && (mgi == grid::get_npts_model() || mgi == oldmgi));
    }
    else
    {
      printout("[fatal] do_rpkt: Failed to identify event . Rpkt. edist %g, sdist %g, tdist %g Abort.\n", edist, sdist, tdist);
      printout("[fatal] do_rpkt: Trouble was due to packet number %d.\n", pkt_ptr->number);
      abort();
      return false;
    }
  }
}


__host__ __device__
void do_rpkt(PKT *pkt_ptr, const double t2)
{
  while (do_rpkt_step(pkt_ptr, t2))
  {
    ;
  }
}


__host__ __device__
static double get_rpkt_escapeprob_fromdirection(const double startpos[3], double start_nu_cmf, int startcellindex, double tstart, double dirvec[3], enum cell_boundary last_cross, double *tot_tau_cont, double *tot_tau_lines)
{
  PKT vpkt;
  vpkt.type = TYPE_RPKT;
  vpkt.nu_cmf = start_nu_cmf;
  vpkt.where = startcellindex;
  vpkt.next_trans = 0;
  vpkt.last_cross = last_cross;

  vec_copy(vpkt.dir, dirvec);
  vec_copy(vpkt.pos, startpos);

  vpkt.prop_time = tstart;
  const double dopplerfactor = doppler_packetpos(&vpkt);
  vpkt.nu_rf = vpkt.nu_cmf / dopplerfactor;

  double t_future = tstart;

  int snext = -99;
  bool end_packet = false;
  while (end_packet == false)
  {
    const int cellindex = vpkt.where;
    const int mgi = grid::get_cell_modelgridindex(cellindex);
    if (grid::modelgrid[mgi].thick == 1)
    {
      return 0.;
    }

    // distance to the next cell
    const double sdist = boundary_cross(&vpkt, t_future, &snext);

    if (snext >= 0)
    {
      const int nextmgi = grid::get_cell_modelgridindex(snext);
      if (grid::modelgrid[nextmgi].thick == 1)
      {
        return 0.;
      }
    }

    vpkt.prop_time = t_future;
    calculate_kappa_rpkt_cont(&vpkt, &globals::kappa_rpkt_cont[tid]);

    const double kappa_cont = globals::kappa_rpkt_cont[tid].total * doppler_packetpos(&vpkt);

    *tot_tau_cont += kappa_cont * sdist;

    if ((*tot_tau_lines + *tot_tau_cont) > 10.)
    {
      // printout("reached tau limit of %g\n", (tot_tau_lines + tot_tau_cont));
      return 0.;
    }

    double ldist = 0.;
    while (ldist < sdist)
    {
      const int lineindex = closest_transition(vpkt.nu_cmf, vpkt.next_trans);

      if (lineindex >= 0)
      {
        const double nutrans = globals::linelist[lineindex].nu;

        vpkt.next_trans = lineindex + 1;

        if (vpkt.nu_cmf < nutrans)
        {
          ldist = 0;
        }
        else
        {
          ldist = CLIGHT * t_future * (vpkt.nu_cmf / nutrans - 1);
        }

        assert_always(ldist >= 0.);

        if (ldist > sdist)
        {
          // exit the while loop if you reach the boundary; go back to the previous transition to start next cell with the excluded line

          vpkt.next_trans -= 1;
          break;
        }

        const double t_line = t_future + ldist / CLIGHT;
        const double tau_line = get_tau_sobolev(mgi, lineindex, t_line);

        *tot_tau_lines += tau_line;
      }
      else
      {
        vpkt.next_trans = globals::nlines + 1;
        break;
      }
    }

    if (snext < 0 || grid::get_cell_modelgridindex(snext) == grid::get_npts_model())
    {
      break;
    }

    t_future += (sdist / globals::CLIGHT_PROP);
    vpkt.prop_time = t_future;
    move_pkt(&vpkt, sdist, t_future);

    if (snext != vpkt.where)
    {
      vpkt.prop_time = t_future;
      change_cell(&vpkt, snext, t_future);
      end_packet = (vpkt.type == TYPE_ESCAPE);
    }
  }

  const double tau_escape = *tot_tau_cont + *tot_tau_lines;
  const double escape_prob = exp(-tau_escape);
  // printout("  tot_tau_lines %g tot_tau_cont %g escape_prob %g\n",
  //          tot_tau_lines, tot_tau_cont, escape_prob);
  return escape_prob;
}


__host__ __device__
double get_rpkt_escape_prob(PKT *pkt_ptr, const double tstart)
{
  // return -1.; // disable this functionality and speed up the code

  const int startcellindex = pkt_ptr->where;
  double startpos[3];
  vec_copy(startpos, pkt_ptr->pos);
  const double start_nu_cmf = pkt_ptr->nu_cmf;
  const enum cell_boundary last_cross = pkt_ptr->last_cross;
  const int mgi = grid::get_cell_modelgridindex(startcellindex);
  if (grid::modelgrid[mgi].thick == 1)
  {
    // escape prob in thick cell is zero
    return 0.;
  }
  const time_t sys_time_start_escape_prob = time(NULL);

  const double pkt_radius = vec_len(startpos);
  const double rmaxnow = globals::rmax * tstart / globals::tmin;
  printout("get_rpkt_escape_prob pkt_radius %g rmax %g r/rmax %g tstart %g\n", pkt_radius, rmaxnow, pkt_radius / rmaxnow, tstart);
  // assert_always(pkt_radius <= rmaxnow);
  double escape_prob_sum = 0.;
  const int ndirs = 40; // number of random directions to sample
  for (int n = 0; n < ndirs; n++)
  {
    double dirvec[3];
    get_rand_isotropic_unitvec(dirvec);
    double tau_cont = 0.;
    double tau_lines = 0.;
    const double escape_prob = get_rpkt_escapeprob_fromdirection(startpos, start_nu_cmf, startcellindex, tstart, dirvec, last_cross, &tau_cont, &tau_lines);
    escape_prob_sum += escape_prob;

    printout("randomdir no. %d (dir dot pos) %g dir %g %g %g tau_lines %g tau_cont %g escape_prob %g escape_prob_avg %g\n",
             n, dot(startpos, dirvec), dirvec[0], dirvec[1], dirvec[2], tau_cont, tau_lines, escape_prob, escape_prob_sum / (n + 1));
  }
  const double escape_prob_avg = escape_prob_sum / ndirs;
  printout("from %d random directions, average escape probability is %g (took %ld s)\n", ndirs, escape_prob_avg, time(NULL) - sys_time_start_escape_prob);

  // reset the cell history and rpkt opacities back to values for the start point
  cellhistory_reset(mgi, false);

  return escape_prob_avg;
}


__host__ __device__
void emitt_rpkt(PKT *pkt_ptr)
{
  /// now make the packet a r-pkt and set further flags
  pkt_ptr->type = TYPE_RPKT;
  pkt_ptr->last_cross = NONE;  /// allow all further cell crossings

  /// Need to assign a new direction. Assume isotropic emission in the cmf

  double dir_cmf[3];
  get_rand_isotropic_unitvec(dir_cmf);

  double vel_vec[3];
  /// This direction is in the cmf - we want to convert it to the rest
  /// frame - use aberation of angles. We want to convert from cmf to
  /// rest so need -ve velocity.
  get_velocity(pkt_ptr->pos, vel_vec, -1. * pkt_ptr->prop_time);
  ///negative time since we want the backwards transformation here

  angle_ab(dir_cmf, vel_vec, pkt_ptr->dir);
  //printout("[debug] pkt_ptr->dir in RF: %g %g %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);

  /// Check unit vector.
  if (fabs(vec_len(pkt_ptr->dir) - 1) > 1.e-8)
  {
    printout("[fatal] emitt_rpkt: Not a unit vector. Abort.\n");
    abort();
  }

  /// Finally we want to put in the rest frame energy and frequency. And record
  /// that it's now a r-pkt.

  if (pkt_ptr->e_cmf >1e52)
  {
    printout("[fatal] emitt_rpkt: here %g\n", pkt_ptr->e_cmf);
    abort();
  }

  const double dopplerfactor = doppler_packetpos(pkt_ptr);
  pkt_ptr->nu_rf = pkt_ptr->nu_cmf / dopplerfactor;
  pkt_ptr->e_rf = pkt_ptr->e_cmf / dopplerfactor;

  // Reset polarization information
  pkt_ptr->stokes[0] = 1.;
  pkt_ptr->stokes[1] = 0.;
  pkt_ptr->stokes[2] = 0.;

  double dummy_dir[3];
  dummy_dir[0] = dummy_dir[1] = 0.0;
  dummy_dir[2] = 1.0;
  cross_prod(pkt_ptr->dir,dummy_dir,pkt_ptr->pol_dir);

  if ((dot(pkt_ptr->pol_dir,pkt_ptr->pol_dir)) < 1.e-8)
  {
    dummy_dir[0] = dummy_dir[2]=0.0;
    dummy_dir[1] = 1.0;
    cross_prod(pkt_ptr->dir,dummy_dir,pkt_ptr->pol_dir);
  }

  vec_norm(pkt_ptr->pol_dir, pkt_ptr->pol_dir);
  //printout("initialise pol state of packet %g, %g, %g, %g, %g\n",pkt_ptr->stokes_qu[0],pkt_ptr->stokes_qu[1],pkt_ptr->pol_dir[0],pkt_ptr->pol_dir[1],pkt_ptr->pol_dir[2]);
  //printout("pkt direction %g, %g, %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
}


__host__ __device__
static double calculate_kappa_ff(const int modelgridindex, const double nu)
/// free-free opacity
{
  assert_always(nu > 0.);
  const double g_ff = 1;

  const float nne = grid::get_nne(modelgridindex);
  const float T_e = grid::get_Te(modelgridindex);

  double kappa_ff = 0.;
  //kappa_ffheating = 0.;

  for (int element = 0; element < get_nelements(); element++)
  {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      ///calculate population of ionstage ...
      const double nnion = ionstagepop(modelgridindex, element, ion);
      //Z = get_element(element);  ///atomic number
      //if (get_ionstage(element,ion) > 1)
      /// Z is ionic charge in the following formula
      const int Z = get_ionstage(element,ion) - 1;
      if (Z > 0)
      {
        //kappa_ff += 3.69255e8 * pow(Z,2) / sqrt(T_e) * pow(nu,-3) * g_ff * nne * nnion * (1-exp(-HOVERKB*nu/T_e));
        kappa_ff += pow(Z, 2) * g_ff * nnion;
        //kappa_ffheating += pow(Z,2) * g_ff * nnion;
        /// heating with level dependence
        //kappa_ffheating += 3.69255e8 * pow(Z,2) / sqrt(T_e) * pow(nu,-3) * g_ff * nne * nnion * (1 - exp(-HOVERKB*nu/T_e));
        /// heating without level dependence
        //kappa_ffheating += 3.69255e8 * pow(Z,2) * pow(nu,-3) * g_ff * (1-exp(-HOVERKB*nu/T_e));
        if (!std::isfinite(kappa_ff))
        {
          printout("kappa_ff %g nne %g T_e %g mgi %d element %d ion %d nnion %g\n", kappa_ff, nne, T_e, modelgridindex, element, ion, nnion);
        }
        assert_always(std::isfinite(kappa_ff));
      }
    }
  }
  kappa_ff *= 3.69255e8 / sqrt(T_e) * pow(nu,-3) * nne * (1 - exp(-HOVERKB * nu / T_e));

  if (!std::isfinite(kappa_ff))
  {
    printout("ERRORL: kappa_ff is non-infinite mgi %d nne %g nu %g T_e %g\n", modelgridindex, nne, nu, T_e);
    abort();
  }
  //kappa_ffheating *= 3.69255e8 / sqrt(T_e) * pow(nu,-3) * nne * (1 - exp(-HOVERKB*nu/T_e));
  //kappa_ff *= 1e5;
  return kappa_ff;
}


__host__ __device__
double calculate_kappa_bf_gammacontr(const int modelgridindex, const double nu)
// bound-free opacity
{
  double kappa_bf_sum = 0.;
#if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
  for (int gphixsindex = 0; gphixsindex < globals::nbfcontinua_ground; gphixsindex++)
  {
    globals::phixslist[tid].groundcont_gamma_contr[gphixsindex] = 0.;
  }
#endif

#if (!SEPARATE_STIMRECOMB)
  const double T_e = grid::get_Te(modelgridindex);
  const double nne = grid::get_nne(modelgridindex);
#endif
  const double nnetot = grid::get_nnetot(modelgridindex);
  for (int i = 0; i < globals::nbfcontinua; i++)
  {
    const int element = globals::allcont[i].element;
    const int ion = globals::allcont[i].ion;
    const int level = globals::allcont[i].level;
    /// The bf process happens only if the current cell contains
    /// the involved atomic species

  	if ((DETAILED_BF_ESTIMATORS_ON && grid::get_elem_abundance(modelgridindex,element) > 0) ||
        (!DETAILED_BF_ESTIMATORS_ON && ((ionstagepop(modelgridindex, element, ion) / nnetot > 1.e-6) || (level == 0))))
    {
      const double nu_edge = globals::allcont[i].nu_edge;
      const int phixstargetindex = globals::allcont[i].phixstargetindex;
      const double nnlevel = get_levelpop(modelgridindex, element, ion, level);
      //printout("i %d, nu_edge %g\n",i,nu_edge);
      const double nu_max_phixs = nu_edge * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table

      if (nu >= nu_edge && nu <= nu_max_phixs && nnlevel > 0)
      {
        //printout("element %d, ion %d, level %d, nnlevel %g\n",element,ion,level,nnlevel);
        const double sigma_bf = photoionization_crosssection(element, ion, level, nu_edge, nu);
        const double probability = get_phixsprobability(element, ion, level, phixstargetindex);

#if (SEPARATE_STIMRECOMB)
        const double corrfactor = 1.; // no subtraction of stimulated recombination
#else
        const int upper = get_phixsupperlevel(element, ion, level, phixstargetindex);
        const double nnupperionlevel = get_levelpop(modelgridindex, element, ion + 1, upper);
        const double sf = calculate_sahafact(element, ion, level, upper, T_e, H * nu_edge);
        const double departure_ratio = nnupperionlevel / nnlevel * nne * sf; // put that to phixslist
        const double stimfactor = departure_ratio * exp(-HOVERKB * nu / T_e);
        double corrfactor = 1 - stimfactor; // photoionisation minus stimulated recombination
        if (corrfactor < 0)
        {
          corrfactor = 0.;
        }
        // const double corrfactor = 1.; // no subtraction of stimulated recombination
#endif

        const double kappa_bf_contr = nnlevel * sigma_bf * probability * corrfactor;

        #if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
        if (level == 0)
        {
          const int gphixsindex = globals::allcont[i].index_in_groundphixslist;
          globals::phixslist[tid].groundcont_gamma_contr[gphixsindex] += sigma_bf * probability * corrfactor;
        }
        #endif

        #if (DETAILED_BF_ESTIMATORS_ON)
        globals::phixslist[tid].gamma_contr[i] = sigma_bf * probability * corrfactor;
        #endif

        if (!std::isfinite(kappa_bf_contr))
        {
          printout("[fatal] calculate_kappa_rpkt_cont: non-finite contribution to kappa_bf_contr %g ... abort\n",kappa_bf_contr);
          printout("[fatal] phixslist index %d, element %d, ion %d, level %d\n",i,element,ion,level);
          printout("[fatal] Z=%d ionstage %d\n", get_element(element), get_ionstage(element, ion));
          printout("[fatal] globals::cell[%d].composition[%d].abundance = %g\n",modelgridindex,element,grid::get_elem_abundance(modelgridindex,element));
          printout("[fatal] nne %g, nnlevel %g, (or %g)\n", grid::get_nne(modelgridindex), nnlevel, get_levelpop(modelgridindex,element,ion,level));
          printout("[fatal] sigma_bf %g, T_e %g, nu %g, nu_edge %g\n",sigma_bf,grid::get_Te(modelgridindex),nu,nu_edge);
          abort();
        }

        kappa_bf_sum += kappa_bf_contr;
        globals::phixslist[tid].kappa_bf_sum[i] = kappa_bf_sum;
      }
      else if (nu < nu_edge) // nu < nu_edge
      {
        /// The phixslist is sorted by nu_edge in ascending order (longest to shortest wavelength)
        /// If nu < allcont[i].nu_edge no absorption in any of the following continua
        /// is possible, therefore leave the loop.

        // the rest of the list shouldn't be accessed
        // but set them to zero to be safe. This is a fast operation anyway.
        // if the packet's nu increases, this function must re-run to re-calculate kappa_bf_contr
        // a slight red-shifting is ignored
        for (int j = i; j < globals::nbfcontinua; j++)
        {
          globals::phixslist[tid].kappa_bf_sum[j] = kappa_bf_sum;
          #if (DETAILED_BF_ESTIMATORS_ON)
          globals::phixslist[tid].gamma_contr[j] = 0.;
          #endif
        }
        break; // all further processes in the list will have larger nu_edge, so stop here
      }
      else
      {
        // ignore this particular process but continue through the list
        globals::phixslist[tid].kappa_bf_sum[i] = kappa_bf_sum;
        #if (DETAILED_BF_ESTIMATORS_ON)
        globals::phixslist[tid].gamma_contr[i] = 0.;
        #endif
      }
    }
    else // no element present or not an important level
    {
      globals::phixslist[tid].kappa_bf_sum[i] = kappa_bf_sum;
      #if (DETAILED_BF_ESTIMATORS_ON)
      globals::phixslist[tid].gamma_contr[i] = 0.;
      #endif
      #if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
      if (globals::allcont[i].level == 0)
      {
        const int gphixsindex = globals::allcont[i].index_in_groundphixslist;
        globals::phixslist[tid].groundcont_gamma_contr[gphixsindex] = 0.;
      }
      #endif
    }
  }
  return kappa_bf_sum;
}


__host__ __device__
void calculate_kappa_rpkt_cont(const PKT *const pkt_ptr, rpkt_cont_opacity_struct *kappa_rpkt_cont_thisthread)
{
  const int cellindex = pkt_ptr->where;
  const int modelgridindex = grid::get_cell_modelgridindex(cellindex);
  assert_always(modelgridindex != grid::get_npts_model());
  assert_always(grid::modelgrid[modelgridindex].thick != 1);
  const double nu_cmf = pkt_ptr->nu_cmf;
  if ((modelgridindex == kappa_rpkt_cont_thisthread->modelgridindex) && (!kappa_rpkt_cont_thisthread->recalculate_required) && (fabs(kappa_rpkt_cont_thisthread->nu / nu_cmf - 1.0) < 1e-4))
  {
    // calculated values are a match already
    return;
  }

  const float nne = grid::get_nne(modelgridindex);

  double sigma = 0.0;
  double kappa_ff = 0.;
  double kappa_bf = 0.;
  double kappa_ffheating = 0.;

  if (globals::do_r_lc)
  {
    if (globals::opacity_case == 4)
    {
      /// First contribution: Thomson scattering on free electrons
      sigma = SIGMA_T * nne;
      //reduced e/s for debugging
      //sigma = 1e-30*sigma;
      //switched off e/s for debugging
      //sigma_cmf = 0. * nne;
      //sigma *= 0.1;

      /// Second contribution: free-free absorption
      kappa_ff = calculate_kappa_ff(modelgridindex, nu_cmf);
      kappa_ffheating = kappa_ff;

      /// Third contribution: bound-free absorption
      kappa_bf = calculate_kappa_bf_gammacontr(modelgridindex, nu_cmf);

      // const double pkt_lambda = 1e8 * CLIGHT / nu_cmf;
      // if (pkt_lambda < 4000)
      // {
      //   printout("lambda %7.1f kappa_bf %g \n", pkt_lambda, kappa_bf);
      // }
    }
    else
    {
      /// in the other cases kappa_grey is an mass absorption coefficient
      /// therefore use the mass density
      //sigma = globals::cell[pkt_ptr->where].kappa_grey * globals::cell[pkt_ptr->where].rho;
      //sigma = SIGMA_T * nne;

      sigma = 0.;
      // kappa_ff = 0.9*sigma;
      // sigma *= 0.1;
      // kappa_bf = 0.;

      // Second contribution: free-free absorption
      kappa_ff = 1e5 * calculate_kappa_ff(modelgridindex, nu_cmf);

      kappa_bf = 0.;
    }
  }

  kappa_rpkt_cont_thisthread->nu = nu_cmf;
  kappa_rpkt_cont_thisthread->modelgridindex = modelgridindex;
  kappa_rpkt_cont_thisthread->recalculate_required = false;
  kappa_rpkt_cont_thisthread->total = sigma + kappa_bf + kappa_ff;
  kappa_rpkt_cont_thisthread->es = sigma;
  kappa_rpkt_cont_thisthread->ff = kappa_ff;
  kappa_rpkt_cont_thisthread->bf = kappa_bf;
  kappa_rpkt_cont_thisthread->ffheating = kappa_ffheating;
  //kappa_rpkt_cont_thisthread.bfheating = kappa_bfheating;

  if (!std::isfinite(kappa_rpkt_cont_thisthread->total))
  {
    printout("[fatal] calculate_kappa_rpkt_cont: resulted in non-finite kappa_rpkt_cont.total ... abort\n");
    printout("[fatal] es %g, ff %g, bf %g\n",
    kappa_rpkt_cont_thisthread->es,kappa_rpkt_cont_thisthread->ff,kappa_rpkt_cont_thisthread->bf);
    printout("[fatal] nbfcontinua %d\n",globals::nbfcontinua);
    printout("[fatal] in cell %d with density %g\n",modelgridindex,grid::get_rho(modelgridindex));
    printout("[fatal] pkt_ptr->nu_cmf %g\n",pkt_ptr->nu_cmf);
    if (std::isfinite(kappa_rpkt_cont_thisthread->es))
    {
      kappa_rpkt_cont_thisthread->ff = 0.;
      kappa_rpkt_cont_thisthread->bf = 0.;
      kappa_rpkt_cont_thisthread->total = kappa_rpkt_cont_thisthread->es;
    }
    else
    {
      abort();
    }
  }

}


