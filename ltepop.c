#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "ltepop.h"
#include "macroatom.h"
#include "nltepop.h"
#include "ratecoeff.h"
#include "update_grid.h"

///***************************************************************************/
double nne_solution_f(double x, void *paras)
/// For libgsl bracketing type solver
/// provides the equation which has to be solved to obtain the electron number
/// density (passed by x)
{
  double innersum, abundance;
  int uppermost_ion;

  int n = ((nne_solution_paras *) paras)->cellnumber;
  //double T = cell[n].T_R;
  double rho = modelgrid[n].rho;

  //printout("n, x, T_R, T_e, W, rho %d, %g, %g, %g, %g, %g\n",n,x,T,cell[n].T_e,cell[n].W,cell[n].rho);
  double outersum = 0.;
  //printout("debug nelements %d =========================\n",nelements);
  for (int element = 0; element < nelements; element++)
  {
    abundance = modelgrid[n].composition[element].abundance;
    if (abundance > 0)
    {
      innersum = 0.;
      //printout("debug get_nions (element %d) %d =========================\n",element,get_nions(element));
      //uppermost_ion = elements[element].uppermost_ion;
      uppermost_ion = elements_uppermost_ion[tid][element];
      /*
      #ifdef FORCE_LTE
        uppermost_ion = get_nions(element)-1;
      #else
        //uppermost_ion = elements[element].uppermost_ion;
        uppermost_ion = elements_uppermost_ion[tid][element];
      #endif
      */
      int ion;
      for (ion = 0; ion <= uppermost_ion; ion++)
      {
        //printout("debug element %d, ion %d, ionfract(element,ion,T,x) %g\n",element,ion,ionfract(element,ion,T,x));
        innersum += (get_ionstage(element,ion)-1) * ionfract(element,ion,n,x);
        if (!isfinite(innersum)) abort();
      }
      outersum += abundance/elements[element].mass * innersum;
      if (!isfinite(outersum))
      {
        printout("nne_solution_f: element %d ion %d uppermostion %d abundance %g, mass %g\n",element,ion,elements_uppermost_ion[tid][element],abundance,elements[element].mass);
        printout("outersum %g\n",outersum);
        abort();
      }
    }
  }

  return rho * outersum - x;
}



///***************************************************************************/
double ionfract(int element, int ion, int modelgridindex, double nne)
/// Calculates ionization fraction for ion=ion of element=element at
/// temperature T and electron number density nne
/// modelgridindex needed to access precalculated partition functions
{
  //printout("debug nne %g\n",nne);
  //int nions = get_nions(element);

  //int uppermost_ion = elements[element].uppermost_ion;
  int uppermost_ion = elements_uppermost_ion[tid][element];
  /*#ifdef FORCE_LTE
    uppermost_ion = get_nions(element)-1;
  #else
    //uppermost_ion = elements[element].uppermost_ion;
    uppermost_ion = elements_uppermost_ion[tid][element];
  #endif*/
  //double phistorage[nions-1-ion];
  double phistorage[uppermost_ion-ion];
  double numerator = 1.;
  //for (i = ion; i < nions-1; i++)
  for (int i = ion; i < uppermost_ion; i++)
  {
    //printout("debug phi(element %d, ion %d)=%g\n",element,i,phi(element,i,modelgridindex));
    phistorage[i-ion] = phi(element,i,modelgridindex);
    numerator *= nne * phistorage[i-ion];
  }
  //printout("debug numerator %g\n",numerator);
  double denominator = 0.;
  //for (i = 0; i < nions; i++)
  for (int i = 0; i <= uppermost_ion; i++)
  {
    double factor = 1.;
    //for (ii = i; ii < nions-1; ii++)
    for (int ii = i; ii < uppermost_ion; ii++)
    {
      if (ii >= ion)
        factor *= nne * phistorage[ii-ion];
      else
        factor *= nne * phi(element,ii,modelgridindex);
    }
    denominator += factor;
  }
  //printout("debug denominator %g with nne %g and T_e %g\n",denominator, nne, get_Te(modelgridindex));

  if (!isfinite(numerator/denominator))
  {
    if (modelgridindex != MGRID)
      printout("[warning] ionfract set to zero for ionstage %d of Z=%d in cell %d with T_e %g, T_R %g\n",get_ionstage(element,ion),get_element(element),modelgridindex,get_Te(modelgridindex),get_TR(modelgridindex));
    //abort();
    printout("debug numerator %g\n",numerator);
    printout("debug denominator %g\n",denominator);
    return 0.;
  }
  return numerator / denominator;
}


///***************************************************************************/
double phi(int element, int ion, int modelgridindex)
/// Calculates population ratio (a saha factor) of two consecutive ionisation stages
/// in nebular approximation phi_j,k* = N_j,k*/(N_j+1,k* * nne)
{
  double phi;

  //double Y_nt, ionpot_in;
  //int element_in, ion_in, nions_in;
  //double rate_use;

  double ionpot = epsilon(element,ion+1,0) - epsilon(element,ion,0);
  //printout("ionpot for element %d, ion %d is %g\n",element,ion,ionpot/EV);
  float T_e = get_Te(modelgridindex);
  //float T_R = get_TR(modelgridindex);

  //float W = cell[cellnumber].W;

  /// Old ionisation formula
  //partfunct_ratio = cell[cellnumber].composition[element].partfunct[ion]/cell[cellnumber].composition[element].partfunct[ion+1];
  //phi = 1./W * sqrt(T_R/T_e) * partfunct_ratio * SAHACONST * pow(T_R,-1.5) * exp(ionpot/KB/T_R);

  /// New ionisation formula with zeta
  //zeta = interpolate_zeta(element,ion,T_e);
  //phi = 1./W * 1./(zeta+W*(1-zeta)) * sqrt(T_R/T_e) * partfunct_ratio * SAHACONST * pow(T_R,-1.5) * exp(ionpot/KB/T_R);

  /// Newest ionisation formula
  double partfunct_ratio = modelgrid[modelgridindex].composition[element].partfunct[ion]/modelgrid[modelgridindex].composition[element].partfunct[ion+1];
  #ifdef FORCE_LTE
    phi = partfunct_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e);
  #else
    if  (initial_iteration == 1 || modelgrid[modelgridindex].thick == 1)
    {
       phi = partfunct_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e);
    }
    else
#ifdef NLTE_POPS_ALL_IONS_SIMULTANEOUS
    {
      phi = ionstagepop(modelgridindex,element,ion+1) / ionstagepop(modelgridindex,element,ion);
    }
#else
    {
      //Gamma = photoionestimator[cellnumber*nelements*maxion+element*maxion+ion];
      double Gamma = gammaestimator[modelgridindex*nelements*maxion+element*maxion+ion]; //try setting to zero
      #ifdef NT_ON
        if (Gamma == 0. && rpkt_emiss[modelgridindex] == 0. && modelgrid[modelgridindex].f48cr == 0. && modelgrid[modelgridindex].fni == 0.)
      #else
        if (Gamma == 0.)
      #endif
        {
          printout("Fatal: Gamma = 0 for element %d, ion %d in phi ... abort\n",element,ion);
          exit(0);
        }

      //Alpha_st = stimrecombestimator[cellnumber*nelements*maxion+element*maxion+ion];
      double Alpha_st = 0.; ///approximate treatment neglects stimulated recombination
      //double Alpha_sp = 0.;
      double Alpha_sp = interpolate_ions_spontrecombcoeff(element,ion,T_e);
      double Col_rec = 0.;

      //     if (ion > 0)
      //	{
      int ionisinglevels = get_bfcontinua(element,ion);
      mastate[tid].element = element;
      mastate[tid].ion = ion+1;
      mastate[tid].nnlevel = 1.0;

      //nlevelsupperion = get_nlevels(element,ion+1);
      for (int level = 0; level < ionisinglevels; level++)
      {
        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element,ion,level); phixstargetindex++)
        {
          int upper = get_phixsupperlevel(element,ion,level,phixstargetindex);
          mastate[tid].level = upper;
          double epsilon_trans = epsilon(element,ion+1,upper) - epsilon(element,ion,level);
          Col_rec += col_recombination(modelgridindex,level,epsilon_trans);
          //Alpha_sp += rad_recombination(modelgridindex,level,epsilon_trans) / get_nne(modelgridindex);
        }
      }

	  //	}

      double recomb_total = Alpha_sp + Alpha_st + (Col_rec / get_nne(modelgridindex));

      /* NT TEST LINES */
      double Y_nt = 0.0;

      #ifdef NT_ON
        Y_nt = nt_ionization_rate(modelgridindex, element, ion);

      /*
      int nlevels = get_nlevels(element,ion);
      double run_tot=0.0;
      for (level = 0; level < nlevels; level++)
      {
        run_tot += calculate_exclevelpop(modelgridindex,element,ion,level);
      }
      Y_nt = Y_nt * run_tot/get_groundlevelpop(modelgridindex, element, ion);
      */

      /*
        if (elements[element].anumber > 20)
          {
            for (element_in = 0; element_in < nelements; element_in++)
              {
                if (elements[element_in].anumber > 20)
                  {
                    nions_in = get_nions(element_in);
                    for (ion_in = 0; ion_in < nions_in -1 ; ion_in++)
                      {
                        ionpot_in = epsilon(element_in,ion_in+1,0) - epsilon(element_in,ion_in,0);
                        {
                          // This is to compute the population-weighted ionization potential
                          K_nt += ionstagepop(modelgridindex,element_in,ion_in) * ionpot_in;
                        }
                      }
                  }
              }
            if (K_nt < 0.0)
              {
                printout("negative K_nt?\n");
                exit(0);
              }

            if (K_nt > 0.0)
              {
                rate_use = rpkt_emiss[modelgridindex] * 1.e20 * 4. * PI;
                // Above is the gamma-ray bit. Below is *supposed* to be the kinetic energy of positrons created by 56Co and 48V. These formulae should be checked, however.
                rate_use += (0.610*0.19*MEV)*(exp(-1.*time_step[nts_global].mid/TCOBALT) - exp(-1.*time_step[nts_global].mid/TNICKEL))/(TCOBALT-TNICKEL)*modelgrid[modelgridindex].fni*get_rho(modelgridindex)/MNI56;
                rate_use += (0.290*0.499*MEV)*(exp(-1.*time_step[nts_global].mid/T48V) - exp(-1.*time_step[nts_global].mid/T48CR))/(T48V-T48CR)*modelgrid[modelgridindex].f48cr*get_rho(modelgridindex)/MCR48;

                // This is the rate-coefficient calculation - K_nt should finally have dimensions like Gamma.
                K_nt = f_nt * rate_use / K_nt ;

                if (K_nt > 0.0)
                  {
		    //                    printout("[NT-print] element %d ion %d K_nt %g Gamma %g Alpha_sp %g Alpha_st %g ne %g cell %d gamma_rate %g total_rate %g \n", element, ion, K_nt, Gamma, Alpha_sp, Alpha_st, get_nne(modelgridindex), modelgridindex, rpkt_emiss[modelgridindex] * 1.e20 * 4. * PI, rate_use);
                  }
              }
          }
      */

      #endif
      /* END OF NT LINES */

      // || !isfinite(Gamma))
      //return phi_lte(element,ion,cellnumber);
      //gamma_lte = interpolate_photoioncoeff_below(element,ion,0,T_e) + interpolate_photoioncoeff_above(element,ion,0,T_e);
      //zeta = interpolate_zeta(element,ion,T_e);
      //alpha_sp = interpolate_spontrecombcoeff(element,ion,0,T_e);
      //phi = gamma_lte*(Alpha_sp+Apha_st)/(Gamma*alpha_sp) * partfunct_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e);

      //printout("testing: Y_nt/Gamma: %g (%d %d)\n", Y_nt/Gamma/stat_weight(element,ion,0)*modelgrid[modelgridindex].composition[element].partfunct[ion], element, ion);
      // OLD
      //phi = (Alpha_sp+Alpha_st)/(Gamma + Y_nt) * modelgrid[modelgridindex].composition[element].partfunct[ion]/
      //stat_weight(element,ion,0);
      //

      //changed July14 to include partition function to stat. weight ratio for upper ion
      // recombinations / ionizations
      //printout("[debug-luke] phi for ion %d Gamma-part %g, Y_nt %g\n",ion,(Gamma * stat_weight(element,ion,0) / modelgrid[modelgridindex].composition[element].partfunct[ion]),Y_nt);
      //Gamma = 0.0; //TODO: testing testing no gamma part
      phi = recomb_total *
          (stat_weight(element,ion+1,0) / modelgrid[modelgridindex].composition[element].partfunct[ion+1])
            /// ((Gamma * stat_weight(element,ion,0) / modelgrid[modelgridindex].composition[element].partfunct[ion]) + Y_nt);
            / Y_nt;
            //Y_nt should be higher than the Gamma term for nebular epoch

      //phi = (Alpha_sp+Alpha_st)/(Y_nt);

      /*if (element == 0)
      {
        printout("phi %g, Alpha_sp %g, Alpha_st %g, Y_nt %g, element %d, ion %d, Col_rec %g, mgi %d, get_nne %g, ratio_u %g ratio_l %g Gamma %g T_e %g\n",
              phi, Alpha_sp, Alpha_st, Y_nt, element, ion, Col_rec, modelgridindex, get_nne(modelgridindex), stat_weight(element,ion+1,0)/modelgrid[modelgridindex].composition[element].partfunct[ion+1], stat_weight(element,ion,0)/modelgrid[modelgridindex].composition[element].partfunct[ion],Gamma, T_e);
      }*/

      if (!isfinite(phi) || phi == 0.0)
      {
        printout("[fatal] phi: phi %g exceeds numerically possible range for element %d, ion %d, T_e %g, T_R %g ... remove higher or lower ionisation stages\n",phi,element,ion,T_e,T_R);
        printout("[fatal] phi: Alpha_sp %g, Alpha_st %g, Gamma %g, partfunct %g, stat_weight %g\n",Alpha_sp,Alpha_st,Gamma,modelgrid[modelgridindex].composition[element].partfunct[ion],stat_weight(element,ion,0));
        printout("[fatal] phi: recomb_total %g, upperionpartfunct %g, upperionstatweight %g\n",recomb_total,modelgrid[modelgridindex].composition[element].partfunct[ion+1],stat_weight(element,ion+1,0));
        printout("[fatal] phi: Y_nt %g Col_rec %g get_nne(modelgridindex) %g\n",Y_nt,Col_rec,get_nne(modelgridindex));
        //abort();
      }
    }
  #endif
#endif //ALL IONS SIMULTANEOUS

  return phi;
}


//***************************************************************************/
//double phi_lte(int element, int ion, int cellnumber)
/// Calculates population ratio (a saha factor) of two consecutive ionisation stages
/// in nebular approximation phi_j,k* = N_j,k*/N_j+1,k* * nne
/*{
  double partfunct_ratio;
  double phi;

  double ionpot = epsilon(element,ion+1,0) - epsilon(element,ion,0);
  float T_e = cell[cellnumber].T_e;

  partfunct_ratio = cell[cellnumber].composition[element].partfunct[ion]/cell[cellnumber].composition[element].partfunct[ion+1];
  phi = partfunct_ratio * SAHACONST * pow(T_e,-1.5) * exp(ionpot/KB/T_e);

  if (!isfinite(phi))
  {
    printout("phi_lte: phi %g exceeds numerically possible range for element %d, ion %d, T_e %g, ... remove higher or lower ionisation stages\n",phi,element,ion,T_e);
    abort();
  }
  return phi;
}*/



///***************************************************************************/
/*
double calculate_ltepartfunct(int element, int ion, double T)
/// Calculates the LTE partition function for ion=ion of element=element at
/// temperature T
{
  double U;
  double epsilon_groundlevel;
  double oneoverkbt;
  int level;
  int nlevels;

  epsilon_groundlevel = epsilon(element,ion,0);
  oneoverkbt = 1/KB/T;
  U = 0.;
  nlevels = get_nlevels(element,ion);
  for (level = 0; level < nlevels; level++)
  {
    U += stat_weight(element,ion,level) * exp(-(epsilon(element,ion,level)-epsilon_groundlevel)*oneoverkbt);
  }

  if (!isfinite(U)) abort();
  return U;
}
*/


///***************************************************************************
double calculate_partfunct_old(int element, int ion, int modelgridindex)
/// Calculates the partition function for ion=ion of element=element in
/// cell modelgridindex
{
  double T_exc = get_TJ(modelgridindex);
  double W = 1.;

  //  T_exc = MINTEMP;

/*  if (T_exc <= MINTEMP)
  {
    T_exc = get_TR(modelgridindex);
    W = get_W(modelgridindex);
  }*/
  //double T_R = modelgrid[modelgridindex].T_R;
  //double W = modelgrid[modelgridindex].W;
  //double T_J = pow(W,1./4.)*T_R;
  double oneoverkbtexc = 1/KB/T_exc;
  double epsilon_groundlevel = epsilon(element,ion,0);

  double U = stat_weight(element,ion,0);
  int nlevels = get_nlevels(element,ion);
/*  if (T_exc <= MINTEMP)
  {
    for (int level = 1; level < nlevels; level++)
    {
      if (elements[element].ions[ion].levels[level].metastable == 1)
      {
        T_exc = get_TJ(modelgridindex);
        W = 1.;
      }
      else
      {
        T_exc = get_TR(modelgridindex);
        W = get_W(modelgridindex);
      }
      U += W*stat_weight(element,ion,level) * exp(-(epsilon(element,ion,level)-epsilon_groundlevel)*oneoverkbtexc);
    }
  }
  else
  {*/
    for (int level = 1; level < nlevels; level++)
    {
      U += W * stat_weight(element,ion,level) * exp(-(epsilon(element,ion,level)-epsilon_groundlevel)*oneoverkbtexc);
      if (!isfinite(U))
      {
        printout("element %d ion %d\n",element,ion);
        printout("modelgridindex %d\n",modelgridindex);
        printout("level %d, nlevels %d\n",level,nlevels);
        printout("sw %g\n",stat_weight(element,ion,0));
        printout("T_exc %g \n",T_exc);
        abort();
      }
    }
//   }

  return U;
}

///***************************************************************************/

double calculate_partfunct(int element, int ion, int modelgridindex)
/// Calculates the partition function for ion=ion of element=element in
/// cell modelgridindex
{
  double pop_store;
  //double E_level, E_ground, test;

  //double T_exc = get_TJ(modelgridindex);
  //double W = 1.;

  //  T_exc = MINTEMP;

/*  if (T_exc <= MINTEMP)
  {
    T_exc = get_TR(modelgridindex);
    W = get_W(modelgridindex);
  }*/
  //double T_R = modelgrid[modelgridindex].T_R;
  //double W = modelgrid[modelgridindex].W;
  //double T_J = pow(W,1./4.)*T_R;
  //double oneoverkbtexc = 1/KB/T_exc;
  //double epsilon_groundlevel = epsilon(element,ion,0);

  int initial = 0;
  if (get_groundlevelpop(modelgridindex,element,ion) < MINPOP)
  {
    //either there reall is none of this ion or this is a first pass through
    //in either case, we won't have any real nlte_populations so the actual value of
    //of groundlevelpop for this calculation doesn't matter, so long as it's not zero!
    pop_store = get_groundlevelpop(modelgridindex,element,ion);
    initial = 1;
    modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = 1.0;
  }

  //printout("groundlevelpop %g\n", get_groundlevelpop(modelgridindex,element,ion));

  double U = 1.0;//stat_weight(element,ion,0);
  int nlevels = get_nlevels(element,ion);
/*  if (T_exc <= MINTEMP)
  {
    for (level = 1; level < nlevels; level++)
    {
      if (elements[element].ions[ion].levels[level].metastable == 1)
      {
        T_exc = get_TJ(modelgridindex);
        W = 1.;
      }
      else
      {
        T_exc = get_TR(modelgridindex);
        W = get_W(modelgridindex);
      }
      U += W*stat_weight(element,ion,level) * exp(-(epsilon(element,ion,level)-epsilon_groundlevel)*oneoverkbtexc);
    }
  }
  else
  {*/
    for (int level = 1; level < nlevels; level++)
    {
      double nn = calculate_exclevelpop(modelgridindex, element, ion, level) / get_groundlevelpop(modelgridindex,element,ion);//*stat_weight(element,ion,0);

      //#ifdef NLTE_POPS_ON
      //if ((is_nlte(element,ion,level) != 1) || (test = modelgrid[modelgridindex].nlte_pops[elements[element].ions[ion].first_nlte+level-1]) < -0.9)
      //{
	  /* Case for when no NLTE level information is available yet */
      //  E_level = epsilon(element,ion,level);
      //  E_ground = epsilon(element,ion,0);
      //  nn = W * stat_weight(element,ion,level) * exp(-(E_level-E_ground)/KB/T_exc);
      //	}
      //else
      //	{
	  //printout("Using an nlte population!\n");
      //	  nn = test * modelgrid[modelgridindex].rho / get_groundlevelpop(modelgridindex,element,ion)*stat_weight(element,ion,0);
      //	  if (!isfinite(nn))
      //	    {
      //
      //      printout("[fatal] NLTE population failure.\n");
      //	      printout("element %d ion %d level %d\n", element, ion, level);
      //	      printout("nn %g test %g rho %g\n", nn, test, modelgrid[modelgridindex].rho);
      //	      printout("ground level %g\n", get_groundlevelpop(modelgridindex,element,ion));
      //	      exit(0);
      //	    }

      //	}
      //#endif
	     U += nn;//W*stat_weight(element,ion,level) * exp(-(epsilon(element,ion,level)-epsilon_groundlevel)*oneoverkbtexc);
    }
//   }

  U *= stat_weight(element,ion,0);

  if (!isfinite(U))
  {
    printout("element %d ion %d\n",element,ion);
    printout("modelgridindex %d\n",modelgridindex);
    printout("nlevels %d\n",nlevels);
    printout("sw %g\n",stat_weight(element,ion,0));
    //printout("T_exc %g \n",T_exc);
    abort();
  }

  if (initial == 1)
  {
    //put back the zero, just in case it matters for something
    modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = pop_store;
  }

  return U;
}


///***************************************************************************/
/*
float calculate_groundlevelpop(int element, int ion, double T, int cellnumber, double nne, double nnnextion)
///calculates ground level population for ion=ion of element=element at
///temperature T and electron number density nne
///further the total population number nnnextion of the next higher ionisation stage is needed
{
  double partfunct(int element, int ion, double T);

  double deltaE = epsilon(element,ion+1,0) - epsilon(element,ion,0);
  double n0;

  //n0 = nnnextion * nne * stat_weight(element,ion,0)/partfunct(element,ion+1,T) * C * pow(T,-1.5) * exp(deltaE/KB/T);
  n0 = nnnextion * nne * stat_weight(element,ion,0)/cell[cellnumber].composition[element].partfunct[ion+1] * SAHACONST * pow(T,-1.5) * exp(deltaE/KB/T);

  return n0;
}
*/

///***************************************************************************/
double get_groundlevelpop(int modelgridindex, int element, int ion)
/// Returns the given ions groundlevel population for modelgridindex which was precalculated
/// during update_grid and stored to the grid.
{
  //double nn = modelgrid[modelgridindex].composition[element].groundlevelpop[ion];
  //if (nn < MINPOP) nn = MINPOP;
  //return nn;

  double nn = modelgrid[modelgridindex].composition[element].groundlevelpop[ion];
  if (nn < MINPOP)
  {
    if (get_abundance(modelgridindex,element) > 0)
      nn = MINPOP;
    else
      nn = 0.;
  }
  return nn;
}


///***************************************************************************/
double calculate_exclevelpop_old(int modelgridindex, int element, int ion, int level)
/// Calculates occupation number of level relative to the ions ground level population
/// using a modified version of the Boltzmann formula, which fulfills the diluted BB
/// approximation (or nebular approximation).
{
  double T_exc = get_TJ(modelgridindex);
  double W = 1.;

  //  T_exc = MINTEMP;


/*  if (T_exc <= MINTEMP)
  {
    if (elements[element].ions[ion].levels[level].metastable == 1)
    {
      T_exc = get_TJ(modelgridindex);
    }
    else
    {
      T_exc = get_TR(modelgridindex);
      W = get_W(modelgridindex);
    }
  }*/

  double nn;
  if (level == 0)
    nn = get_groundlevelpop(modelgridindex,element,ion);
  else
  {
    double E_level = epsilon(element,ion,level);
    double E_ground = epsilon(element,ion,0);
    nn = get_groundlevelpop(modelgridindex,element,ion) * W * stat_weight(element,ion,level) / stat_weight(element,ion,0) * exp(-(E_level-E_ground)/KB/T_exc);
  }

  if (nn < MINPOP)
  {
    if (get_abundance(modelgridindex,element) > 0)
      nn = MINPOP;
    else
      nn = 0.;
  }

  #ifdef DEBUG_ON
    if (!isfinite(nn))
    {
      printout("[fatal] calculate_exclevelpop: level %d of ion %d of element %d has infinite level population %g\n",level,ion,element,nn);
      printout("[fatal] calculate_exclevelpop: associated ground level has pop %g\n",get_groundlevelpop(modelgridindex,element,ion));
      printout("[fatal] calculate_exclevelpop: associated ion has pop %g\n",ionstagepop(modelgridindex,element,ion));
      printout("[fatal] calculate_exclevelpop: associated partition function %g\n",modelgrid[modelgridindex].composition[element].partfunct[ion]);
    }
  #endif
  return nn;
}


double calculate_levelpop_lte(int modelgridindex, int element, int ion, int level)
/// Calculates occupation population of a level assuming LTE excitation
{
  double nn;

  double T_exc = get_TJ(modelgridindex);
  double W = 1.;

  double E_level = epsilon(element,ion,level);
  double E_ground = epsilon(element,ion,0);
  nn = get_groundlevelpop(modelgridindex,element,ion) * W *
       stat_weight(element,ion,level) / stat_weight(element,ion,0) *
       exp(-(E_level - E_ground)/KB/T_exc);

  if (nn < MINPOP)
  {
    if (get_abundance(modelgridindex,element) > 0)
      nn = MINPOP;
    else
      nn = 0.;
  }

  return nn;
}


///***************************************************************************/
double calculate_exclevelpop(int modelgridindex, int element, int ion, int level)
/// Calculates the population of a level
/// (if in LTE mode) using a modified version of the Boltzmann formula, which fulfills the diluted BB
/// approximation (or nebular approximation).
{
  double nn;

  double T_exc = get_TJ(modelgridindex);
  double W = 1.;

#ifdef NLTE_POPS_ON
  double test;
  int nlte_levels;
#endif

  //  T_exc = MINTEMP;


/*  if (T_exc <= MINTEMP)
  {
    if (elements[element].ions[ion].levels[level].metastable == 1)
    {
      T_exc = get_TJ(modelgridindex);
    }
    else
    {
      T_exc = get_TR(modelgridindex);
      W = get_W(modelgridindex);
    }
  }*/

  if (level == 0)
  {
    nn = get_groundlevelpop(modelgridindex,element,ion);
  }
#ifdef NLTE_POPS_ON
  else if (is_nlte(element,ion,level) == 1)
  {
    //printout("Using an nlte population!\n");
    if ((test = modelgrid[modelgridindex].nlte_pops[elements[element].ions[ion].first_nlte+level-1]) < -0.9)
    {
      /* Case for when no NLTE level information is available yet */
      double E_level = epsilon(element,ion,level);
      double E_ground = epsilon(element,ion,0);
      nn = get_groundlevelpop(modelgridindex,element,ion) * W *
          stat_weight(element,ion,level)/stat_weight(element,ion,0) *
          exp(-(E_level-E_ground)/KB/T_exc);
    }
    else
    {
      //printout("Using an nlte population!\n");
      nn = test * modelgrid[modelgridindex].rho;
      if (!isfinite(nn))
      {
        printout("[fatal] NLTE population failure.\n");
        printout("element %d ion %d level %d\n", element, ion, level);
        printout("nn %g test %g rho %g\n", nn, test, modelgrid[modelgridindex].rho);
        printout("ground level %g\n", get_groundlevelpop(modelgridindex,element,ion));
        exit(0);
      }
      return nn;
    }
  }
  else if ((nlte_levels = get_nlevels_nlte(element,ion)) > 0)
  {
    //Case where this ion HAS nlte levels, but this isn't one of them. Then we want to use the super level to guesstimate it.
    if ((test = modelgrid[modelgridindex].nlte_pops[elements[element].ions[ion].first_nlte+nlte_levels]) < -0.9) //TODO: should change this to less than zero?
    {
      /* Case for when no NLTE level information is available yet */
      double E_level = epsilon(element,ion,level);
      double E_ground = epsilon(element,ion,0);
      nn = get_groundlevelpop(modelgridindex,element,ion) * W *
           stat_weight(element,ion,level)/stat_weight(element,ion,0) *
           exp(-(E_level-E_ground)/KB/T_exc);
    }
    else
    {
      //printout("Using a superlevel population!\n");
      nn = test * modelgrid[modelgridindex].rho * superlevel_boltzmann(modelgridindex,element,ion,level); //TOOD: fix
      if (!isfinite(nn))
      {
        printout("[fatal] NLTE population failure.\n");
        printout("element %d ion %d level %d\n", element, ion, level);
        printout("nn %g test %g rho %g\n", nn, test, modelgrid[modelgridindex].rho);
        printout("ground level %g\n", get_groundlevelpop(modelgridindex,element,ion));
        exit(0);
      }
      return nn;
    }
  }
#endif
  else
  {
    double E_level = epsilon(element,ion,level);
    double E_ground = epsilon(element,ion,0);
    nn = get_groundlevelpop(modelgridindex,element,ion) * W *
         stat_weight(element,ion,level)/stat_weight(element,ion,0) *
         exp(-(E_level-E_ground)/KB/T_exc);
  }

  if (nn < MINPOP)
  {
    if (get_abundance(modelgridindex,element) > 0)
      nn = MINPOP;
    else
      nn = 0.;
  }

  #ifdef DEBUG_ON
    if (!isfinite(nn))
    {
      printout("[fatal] calculate_exclevelpop: level %d of ion %d of element %d has infinite level population %g\n",level,ion,element,nn);
      printout("[fatal] calculate_exclevelpop: associated ground level has pop %g\n",get_groundlevelpop(modelgridindex,element,ion));
      printout("[fatal] calculate_exclevelpop: associated ion has pop %g\n",ionstagepop(modelgridindex,element,ion));
      printout("[fatal] calculate_exclevelpop: associated partition function %g\n",modelgrid[modelgridindex].composition[element].partfunct[ion]);
    }
  #endif
  return nn;
}

/*
#ifdef FORCE_LTE
{
  double nn;

  double T_R = get_TR(modelgridindex);
  //double W = cell[cellnumber].W;

  if (level == 0) nn = get_groundlevelpop(modelgridindex,element,ion);
  else nn = get_groundlevelpop(modelgridindex,element,ion) * stat_weight(element,ion,level)/stat_weight(element,ion,0) * exp(-(epsilon(element,ion,level)-epsilon(element,ion,0))/KB/T_R);

  return nn;
}
#else
{
  double E_level,E_ground;
  double nn;

  double T_J = get_TJ(modelgridindex);
  double T_R = get_TR(modelgridindex);
  double W = get_W(modelgridindex);
  //double T_J = pow(W,1./4.)*T_R;

  if (level == 0) nn = get_groundlevelpop(modelgridindex,element,ion);
  else
  {
    E_level = epsilon(element,ion,level);
    E_ground = epsilon(element,ion,0);
    nn = get_groundlevelpop(modelgridindex,element,ion) * stat_weight(element,ion,level)/stat_weight(element,ion,0) * exp(-(E_level-E_ground)/KB/T_J);
  }

  if (nn < MINPOP) nn = MINPOP;
  return nn;
}
#endif
*/


///***************************************************************************/
/*void calculate_levelpops(int modelgridindex)
/// Calculates the full level populations for a given grid cell
/// and stores them to the active entry of the cellhistory.
{
  for (int element = 0; element < nelements; element++)
  {
    int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++)
    {
      cellhistory[tid].chelements[element].chions[ion].chlevels[0].population = get_groundlevelpop(modelgridindex, element, ion);
      //printout("element %d, ion %d, level 0: population %g\n",element,ion,groundlevelpop(cellnumber, element, ion));
      int nlevels = get_nlevels(element,ion);
      for (int level = 1; level < nlevels; level++)
      {
        cellhistory[tid].chelements[element].chions[ion].chlevels[level].population = calculate_exclevelpop(modelgridindex,element,ion,level);
        //printout("element %d, ion %d, level %d: population %g\n",element,ion,level,exclevelpop(cellnumber,element,ion,level,T));
      }
    }
  }
}*/


///***************************************************************************/
double calculate_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold)
/// calculates saha factor in LTE: Phi_level,ion,element = nn_level,ion,element/(nne*nn_upper,ion+1,element)
{
  int upperionlevel = get_phixsupperlevel(element,ion,level,phixstargetindex);
  double sf = stat_weight(element,ion,level) / stat_weight(element,ion+1,upperionlevel) * SAHACONST * pow(T,-1.5) * exp(E_threshold/KB/T);
  //printout("element %d, ion %d, level %d, T, %g, E %g has sf %g (g_l %g g_u %g)\n", element, ion, level, T, E_threshold, sf,stat_weight(element,ion,level),stat_weight(element,ion+1,0) );
  if (sf < 0)
  {
    printout("[fatal] sahafact: negative saha factor");
    abort();
  }
  return sf;
}


///***************************************************************************/
double get_sahafact(int element, int ion, int level, int phixstargetindex, double T, double E_threshold)
/// retrieves or calculates saha factor in LTE: Phi_level,ion,element = nn_level,ion,element/(nne*nn_upper,ion+1,element)
{
  double sf;

  if (use_cellhist >= 0)
  {
    sf = cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].sahafact;
    if (sf < 0)
    {
      sf = calculate_sahafact(element,ion,level,phixstargetindex,T,E_threshold);
      cellhistory[tid].chelements[element].chions[ion].chlevels[level].chphixstargets[phixstargetindex].sahafact = sf;
    }
  }
  else
    sf = calculate_sahafact(element,ion,level,phixstargetindex,T,E_threshold);

  //printout("get_sahafact: sf= %g\n",sf);
  return sf;
}


/// Initialise estimator arrays which hold the last time steps values (used to damp out
/// fluctuations over timestep iterations if DO_TITER is defined) to -1.
void initialise_photoionestimators()
{
  //for (n = 0; n < ngrid; n++)
  for (int n = 0; n < npts_model; n++)
  {
    //double T_e = get_Te(n);
    #ifdef DO_TITER
      J_reduced_save[n] = -1.;
    #endif
    #ifndef FORCE_LTE
      #ifdef DO_TITER
        nuJ_reduced_save[n] = -1.;
        ffheatingestimator_save[n] = -1.;
        colheatingestimator_save[n] = -1.;
      #endif
      for (int element = 0; element < nelements; element++)
      {
        int nions = get_nions(element);
        for (int ion = 0; ion < nions-1; ion++)
        {
          //  double ionpot,Alpha_sp,sw_ratio,Gamma;
          //ionpot = epsilon(element,ion+1,0) - epsilon(element,ion,0);
          //Alpha_sp = interpolate_ions_spontrecombcoeff(element,ion,T_e);
          //sw_ratio = stat_weight(element,ion+1,0)/stat_weight(element,ion,0);
          //Gamma = Alpha_sp * sw_ratio / SAHACONST * pow(T_e,1.5) * exp(-ionpot/KB/T_e);
          ////gamma_lte = interpolate_photoioncoeff_below(element,ion,0,T_e) + interpolate_photoioncoeff_above(element,ion,0,T_e);
          ////zeta = interpolate_zeta(element,ion,T_e);
          //gammaestimator[n*nelements*maxion+element*maxion+ion] = Gamma; //gamma_lte/zeta;
          ////corrphotoionrenorm[n*nelements*maxion+element*maxion+ion] = 1.;
          ////photoionestimator[n*nelements*maxion+element*maxion+ion] = Gamma; //gamma_lte/zeta;

          #ifdef DO_TITER
            gammaestimator_save[n*nelements*maxion+element*maxion+ion] = -1.;
            bfheatingestimator_save[n*nelements*maxion+element*maxion+ion] = -1.;
            /*
            photoionestimator_save[n*nelements*maxion+element*maxion+ion] = -1.;
            stimrecombestimator_save[n*nelements*maxion+element*maxion+ion] = -1.;
            ionfluxestimator_save[n*nelements*maxion+element*maxion+ion] = -1.;
            */
          #endif
        }
      }
    #endif
  }
}
