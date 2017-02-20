#include "sn3d.h"
#include "vpkt.h"
#include <gsl/gsl_integration.h>

  
///****************************************************************************
void calculate_kpkt_rates(int modelgridindex)
/// Set up the global cooling list and determine the important entries
/// by applying a cut to the total cooling rate. Then sort the global
/// cooling list by the strength of the individual process contribution.
{
  double get_abundance(int modelgridindex, int element);
  double get_bfcooling(int element, int ion, int level, int modelgridindex);
  double col_excitation(int modelgridindex, int upper, int lineindex, double epsilon_trans);
  double col_ionization(int modelgridindex, int upper, double epsilon_trans);
  double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
  double epsilon(int element, int ion, int level);
  int get_ionstage(int element, int ion);
  double ionstagepop(int modelgridindex, int element, int ion);
  double get_groundlevelpop(int modelgridindex, int element, int ion);
  //double get_levelpop(int element, int ion, int level);
  int get_bfcontinua(int element, int ion);

  double C,C_ff,C_fb,C_exc,C_ion,C_col; //,total_k_destruct;
  double E_threshold;//,nu_threshold;
  //double alpha_sp,modified_alpha_sp;
  double nncurrention,nnnextionlevel;
  double epsilon_current,epsilon_upper,epsilon_trans;
  int element,ion,level,upper;
  int nlevels_currention,ionisinglevels;
  int i,ii,lineindex;
  int ioncharge;
  //double sf;
  int nions,nuptrans;
  double totalcooling;//partialcoolingsum;
  double contrib;

  //double interpolate_stimulated_bfcoolingcoeff(int element, int ion, int level, double T);
        
  double nne = get_nne(modelgridindex);
  double T_e = get_Te(modelgridindex);
  double T_R = get_TR(modelgridindex);
  double W = get_W(modelgridindex);
  //if (SILENT == 0) printout("[info] kpkt_cuts: sampling cell %d, T_e %g, T_R %g, W %g, nne %g\n",modelgridindex,T_e,T_R,W,nne);
  
/*  PKT dummypkt;
  dummypkt.where = cellnumber;
  PKT *pkt_ptr;
  pkt_ptr = &dummypkt;*/
  
  
  /// calculate rates for
  C = 0.;
  //C_ff = 0.;   /// free-free creation of rpkts
  //C_fb = 0.;   /// free-bound creation of rpkt
  //C_exc = 0.;  /// collisional excitation of macroatoms
  //C_ion = 0.;  /// collisional ionisation of macroatoms
  contrib = 0.;
  i = 0;
  for (element=0; element < nelements; element++)
  {
    //printout("[debug] do_kpkt: element %d\n",element);
    mastate[tid].element = element;
    nions = get_nions(element);
    //if (get_abundance(modelgridindex,element) > 0.)
    //{
      for (ion=0; ion < nions; ion++)
      {
        C_ion = 0.;
        //printout("[debug] do_kpkt: ion %d\n",ion);
        mastate[tid].ion = ion;
        nlevels_currention = get_nlevels(element,ion);
        //ionisinglevels = get_ionisinglevels(element,ion);
        ionisinglevels = get_bfcontinua(element,ion);
        nnnextionlevel = get_groundlevelpop(modelgridindex,element,ion+1);
        nncurrention = ionstagepop(modelgridindex,element,ion);
        
        /// ff creation of rpkt
        /// -------------------
        ioncharge = get_ionstage(element,ion) - 1;
        //printout("[debug] ioncharge %d, nncurrention %g, nne %g\n",ion,nncurrention,nne);
        if (ioncharge > 0)
        {
          C = 1.426e-27 * sqrt(T_e) * pow(ioncharge,2) * nncurrention * nne;
          //C_ff += C;
          C_ion += C;
          //cellhistory[tid].coolinglist[i].contribution = C;
          contrib += C;
/*          cellhistory[tid].coolinglist[i].contribution = contrib;
          cellhistory[tid].coolinglist[i].type = COOLINGTYPE_FF;
          cellhistory[tid].coolinglist[i].element = element;
          cellhistory[tid].coolinglist[i].ion = ion;
          cellhistory[tid].coolinglist[i].level = -99;
          cellhistory[tid].coolinglist[i].upperlevel = -99;*/
          i += 1;
        }
        
        for (level = 0; level < nlevels_currention; level++)
        {
          //printout("[debug] do_kpkt: element %d, ion %d, level %d\n",element,ion,level);
          epsilon_current = epsilon(element,ion,level);
          mastate[tid].level = level;
          ///Use the cellhistory populations here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          mastate[tid].nnlevel = calculate_exclevelpop(modelgridindex,element,ion,level);
          //mastate[tid].nnlevel = get_levelpop(element,ion,level);
          
          /// excitation to same ionization stage
          /// -----------------------------------
          nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
          for (ii = 1; ii <= nuptrans; ii++)
          {
            upper = elements[element].ions[ion].levels[level].uptrans[ii].targetlevel;
            lineindex = elements[element].ions[ion].levels[level].uptrans[ii].lineindex;
            //printout("    excitation to level %d possible\n",upper);
            //epsilon_trans = epsilon(element,ion,upper) - epsilon_current;
            epsilon_trans = elements[element].ions[ion].levels[level].uptrans[ii].epsilon - epsilon_current;
            C = col_excitation(modelgridindex,upper,lineindex,epsilon_trans) * epsilon_trans;
            //C = 0.;
            //C_exc += C;
            C_ion += C;
            //cellhistory[tid].coolinglist[i].contribution = C;
            contrib += C;
/*            cellhistory[tid].coolinglist[i].contribution = contrib;
            cellhistory[tid].coolinglist[i].type = COOLINGTYPE_COLLEXC;
            cellhistory[tid].coolinglist[i].element = element;
            cellhistory[tid].coolinglist[i].ion = ion;
            cellhistory[tid].coolinglist[i].level = level;
            cellhistory[tid].coolinglist[i].upperlevel = upper;
            cellhistory[tid].coolinglist[i].lineindex = lineindex;*/
            i += 1;
            //linecounter += 1;
          }
          
          if (ion < nions-1 && level < ionisinglevels) ///check whether further ionisation stage available
          //if (ion < get_nions(element)-1) ///check whether further ionisation stage available
          {
            //printout("    ionisation possible\n");
            /// ionization to higher ionization stage
            /// -------------------------------------
            /// for the moment we deal only with ionisations to the next ions groundlevel
            /// to allow an easy generalisation this was only made sure by col_ionization
            /// for speed issues (reduced number of calls to epsilon) it is now done also 
            /// here explicitly
            //nlevels_upperion = get_nlevels(element,ion+1);
            //for (upper = 0; upper < nlevels_upperion; upper++)
            //{
            upper = 0;
            epsilon_upper = epsilon(element,ion+1,upper);
            epsilon_trans = epsilon_upper - epsilon_current;
            //printout("cooling list: col_ionization\n");
            C = col_ionization(modelgridindex,upper,epsilon_trans) * epsilon_trans;
            C_ion += C;
            //C = 0.;
            //cellhistory[tid].coolinglist[i].contribution = C;
            contrib += C;
/*            cellhistory[tid].coolinglist[i].contribution = contrib;
            cellhistory[tid].coolinglist[i].type = COOLINGTYPE_COLLION;
            cellhistory[tid].coolinglist[i].element = element;
            cellhistory[tid].coolinglist[i].ion = ion;
            cellhistory[tid].coolinglist[i].level = level;
            cellhistory[tid].coolinglist[i].upperlevel = upper;*/
            i += 1;
            //}
            
            
            /// fb creation of r-pkt
            /// free bound rates are calculated from the lower ion, but associated to the higher ion
            /// --------------------
            //upper = 0;
            //epsilon_upper = epsilon(element,ion+1,0);
            //E_threshold = epsilon_upper - epsilon_current; 
            //E_threshold = epsilon_trans;
  //printout("get_bfcooling(%d,%d,%d,%d) for histindex %d\n",element,ion,level,modelgridindex,histindex);
            //epsilon_upper = epsilon(element,ion+1,0);
            //E_threshold = epsilon_upper - epsilon_current;
            C = get_bfcooling(element,ion,level,modelgridindex);
            //printout("element %d, ion %d, level %d, T_e %g, alpha_E - alpha %g, bfcoolingcoeff %g\n",element,ion,level,T_e,C/E_threshold,C);
            //double interpolate_stimulated_recomb(int element, int ion, int level, double T);
            //C = get_bfcooling(element,ion,level,modelgridindex) + (nnnextionlevel*nne * W * interpolate_stimulated_bfcoolingcoeff(element,ion,level,T_R));
            //printout("nonE %g, E %g \n",interpolate_bfcoolingcoeff(element,ion,level,T_e), interpolate_stimulated_bfcoolingcoeff(element,ion,level,T_R));
            //C = get_bfcooling(element,ion,level,modelgridindex) + (nnnextionlevel*nne * (stimrecombestimator_E_save[pkt_ptr->where*nelements*maxion+element*maxion+ion]-stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+ion]) * E_threshold);
            //printout("element %d, ion %d, modified %g, usual %g, diff %g, nonE %g\n",element,ion,stimrecombestimator_E_save[pkt_ptr->where*nelements*maxion+element*maxion+ion],stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+ion],(stimrecombestimator_E_save[pkt_ptr->where*nelements*maxion+element*maxion+ion]-stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+ion])*E_threshold,interpolate_bfcoolingcoeff(element,ion,level,T_e));
            //double interpolate_spontrecombcoeff(int element, int ion, int level, double T);
            //printout("alpha_sp %g , alpha_st %g\n",interpolate_spontrecombcoeff(element,ion,level,T_e),interpolate_stimulated_recomb(element,ion,level,T_R));
            //epsilon_upper = epsilon(element,ion+1,0);
            //E_threshold = epsilon_upper - epsilon_current;
            //C = interpolate_spontrecombcoeff(element,ion,level,T_e) * E_threshold * nnnextionlevel*nne;
            C_ion += C;
            //C_fb += C;
            //cellhistory[tid].coolinglist[i].contribution = C;
            contrib += C;
/*            cellhistory[tid].coolinglist[i].contribution = contrib;
            cellhistory[tid].coolinglist[i].type = COOLINGTYPE_FB;
            cellhistory[tid].coolinglist[i].element = element;
            cellhistory[tid].coolinglist[i].ion = ion;
            cellhistory[tid].coolinglist[i].level = level; 
            cellhistory[tid].coolinglist[i].upperlevel = 0;*/
            i += 1;
          }
        }
      
        modelgrid[modelgridindex].cooling[element].contrib[ion] = C_ion;
      }
    //}
    //else
    //{
    //  for (ion=0; ion < nions; ion++)
    //  {
    //    C_ion = 0.;
    //    modelgrid[modelgridindex].cooling[element].contrib[ion] = C_ion;
    //  }
    //}
  }
      
  //C_col = C_exc + C_ion;
  //totalcooling = C_col + C_ff + C_fb;
  //importantcoolingterms = ncoolingterms;
  //cellhistory[tid].totalcooling = contrib;
  modelgrid[modelgridindex].totalcooling = contrib;
  //printout("[debug] calculate_kpkt_rates: total cooling rate %g\n",contrib);
  
  //printout("[info] kpkt_cuts: C_ff %g, C_fb %g, C_col %g, C_exc %g, C_ion %g\n",C_ff,C_fb,C_col,C_exc,C_ion);
  /*
  if (SILENT == 0) printout("[info] kpkt_cuts: term counter %d, ioncounter %d, linecounter %d, levelcounter %d\n",iii,ioncounter,linecounter,levelcounter);
  if (SILENT == 0) printout("[info] kpkt_cuts: C_ff %g, C_fb %g, C_col %g, C_exc %g, C_ion %g\n",C_ff,C_fb,C_col,C_exc,C_ion);
  if (SILENT == 0) printout("[info] kpkt_cuts: totalcooling %g\n",totalcooling);
  
  
  /// Sort this the coolinglist
  qsort(cellhistory[tid].coolinglist,ncoolingterms,sizeof(coolinglist_entry),compare_coolinglistentry);
  //sort_coolinglist(globalcoolinglist,totalcooling);
  
  /// Now determine the number of important coolingterms cuts by applying the COOLINGCUT
  i = 0;
  partialcoolingsum = 0;
  while ((COOLINGCUT*totalcooling - partialcoolingsum) > 0 && i < ncoolingterms)
  {
    partialcoolingsum += globalcoolinglist[i].contribution;
    i++;
  }
  importantcoolingterms = i; 
  if (importantcoolingterms > ncoolingterms)
  {
    printout("[fatal] kpkt_cuts: cooling list in history not big enough to hold all important cooling contributions ... abort\n");
    printout("[fatal] kpkt_cuts: important coolingterms %d, ncoolingterms %d\n",importantcoolingterms_global,ncoolingterms);
    abort();
  }
  //printout("[info] kpkt_cuts: __important coolingterms__ %d\n",importantcoolingterms_global);
  */
}



///****************************************************************************
void calculate_kpkt_rates_ion(int modelgridindex, int element, int ion, int low, double oldcoolingsum, int high)
/// Set up the global cooling list and determine the important entries
/// by applying a cut to the total cooling rate. Then sort the global
/// cooling list by the strength of the individual process contribution.
{
  double get_bfcooling(int element, int ion, int level, int modelgridindex);
  double col_excitation(int modelgridindex, int upper, int lineindex, double epsilon_trans);
  double col_ionization(int modelgridindex, int upper, double epsilon_trans);
  //double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
  double epsilon(int element, int ion, int level);
  int get_ionstage(int element, int ion);
  double ionstagepop(int modelgridindex, int element, int ion);
  double get_groundlevelpop(int modelgridindex, int element, int ion);
  double get_levelpop(int element, int ion, int level);
  int get_bfcontinua(int element, int ion);

  double C,C_ff,C_fb,C_exc,C_ion,C_col; //,total_k_destruct;
  double E_threshold;//,nu_threshold;
  //double alpha_sp,modified_alpha_sp;
  double nncurrention,nnnextionlevel;
  double epsilon_current,epsilon_upper,epsilon_trans;
  int level,upper;
  int nlevels_currention,ionisinglevels;
  int i,ii,lineindex;
  int ioncharge;
  //double sf;
  int nions,nuptrans;
  double totalcooling;//partialcoolingsum;
  double contrib;

  //double interpolate_stimulated_bfcoolingcoeff(int element, int ion, int level, double T);
        
  double nne = get_nne(modelgridindex);
  double T_e = get_Te(modelgridindex);
  double T_R = get_TR(modelgridindex);
  double W = get_W(modelgridindex);
  
  
  /// calculate rates for
  C = 0.;
  //C_ff = 0.;   /// free-free creation of rpkts
  //C_fb = 0.;   /// free-bound creation of rpkt
  //C_exc = 0.;  /// collisional excitation of macroatoms
  //C_ion = 0.;  /// collisional ionisation of macroatoms
  contrib = oldcoolingsum;
  i = low;
  
  //printout("[debug] do_kpkt: element %d\n",element);
  nions = get_nions(element);
  mastate[tid].element = element;
  mastate[tid].ion = ion;
  nlevels_currention = get_nlevels(element,ion);
  //ionisinglevels = get_ionisinglevels(element,ion);
  ionisinglevels = get_bfcontinua(element,ion);
  nnnextionlevel = get_groundlevelpop(modelgridindex,element,ion+1);
  nncurrention = ionstagepop(modelgridindex,element,ion);
  
  /// ff creation of rpkt
  /// -------------------
  ioncharge = get_ionstage(element,ion) - 1;
  //printout("[debug] ioncharge %d, nncurrention %g, nne %g\n",ion,nncurrention,nne);
  if (ioncharge > 0)
  {
    C = 1.426e-27 * sqrt(T_e) * pow(ioncharge,2) * nncurrention * nne;
    //C_ff += C;
    //C_ion += C;
    //cellhistory[tid].coolinglist[i].contribution = C;
    contrib += C;
    cellhistory[tid].coolinglist[i].contribution = contrib;
    cellhistory[tid].coolinglist[i].type = COOLINGTYPE_FF;
    cellhistory[tid].coolinglist[i].element = element;
    cellhistory[tid].coolinglist[i].ion = ion;
    cellhistory[tid].coolinglist[i].level = -99;
    cellhistory[tid].coolinglist[i].upperlevel = -99;
    //if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d, coolingtype %d, i %d, low %d\n",contrib,oldcoolingsum,C,element,ion,-99,cellhistory[tid].coolinglist[i].type,i,low);
    i += 1;
  }
  
  for (level = 0; level < nlevels_currention; level++)
  {
    //printout("[debug] do_kpkt: element %d, ion %d, level %d\n",element,ion,level);
    epsilon_current = epsilon(element,ion,level);
    mastate[tid].level = level;
    ///Use the cellhistory populations here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //mastate[tid].nnlevel = calculate_exclevelpop(modelgridindex,element,ion,level);
    mastate[tid].nnlevel = get_levelpop(element,ion,level);
    
    /// excitation to same ionization stage
    /// -----------------------------------
    nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
    for (ii = 1; ii <= nuptrans; ii++)
    {
      upper = elements[element].ions[ion].levels[level].uptrans[ii].targetlevel;
      lineindex = elements[element].ions[ion].levels[level].uptrans[ii].lineindex;
      //printout("    excitation to level %d possible\n",upper);
      //epsilon_trans = epsilon(element,ion,upper) - epsilon_current;
      epsilon_trans = elements[element].ions[ion].levels[level].uptrans[ii].epsilon - epsilon_current;
      C = col_excitation(modelgridindex,upper,lineindex,epsilon_trans) * epsilon_trans;
      //C = 0.;
      //C_exc += C;
      //C_ion += C;
      //cellhistory[tid].coolinglist[i].contribution = C;
      contrib += C;
      cellhistory[tid].coolinglist[i].contribution = contrib;
      cellhistory[tid].coolinglist[i].type = COOLINGTYPE_COLLEXC;
      cellhistory[tid].coolinglist[i].element = element;
      cellhistory[tid].coolinglist[i].ion = ion;
      cellhistory[tid].coolinglist[i].level = level;
      cellhistory[tid].coolinglist[i].upperlevel = upper;
      cellhistory[tid].coolinglist[i].lineindex = lineindex;
      //if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d, coolingtype %d, i %d, low %d\n",contrib,oldcoolingsum,C,element,ion,level,cellhistory[tid].coolinglist[i].type,i,low);
      i += 1;
      //linecounter += 1;
    }
    
    if (ion < nions-1 && level < ionisinglevels) ///check whether further ionisation stage available
    //if (ion < get_nions(element)-1) ///check whether further ionisation stage available
    {
      //printout("    ionisation possible\n");
      /// ionization to higher ionization stage
      /// -------------------------------------
      /// for the moment we deal only with ionisations to the next ions groundlevel
      /// to allow an easy generalisation this was only made sure by col_ionization
      /// for speed issues (reduced number of calls to epsilon) it is now done also 
      /// here explicitly
      //nlevels_upperion = get_nlevels(element,ion+1);
      //for (upper = 0; upper < nlevels_upperion; upper++)
      //{
      upper = 0;
      epsilon_upper = epsilon(element,ion+1,upper);
      epsilon_trans = epsilon_upper - epsilon_current;
      //printout("cooling list: col_ionization\n");
      C = col_ionization(modelgridindex,upper,epsilon_trans) * epsilon_trans;
      //C_ion += C;
      //C = 0.;
      //cellhistory[tid].coolinglist[i].contribution = C;
      contrib += C;
      cellhistory[tid].coolinglist[i].contribution = contrib;
      cellhistory[tid].coolinglist[i].type = COOLINGTYPE_COLLION;
      cellhistory[tid].coolinglist[i].element = element;
      cellhistory[tid].coolinglist[i].ion = ion;
      cellhistory[tid].coolinglist[i].level = level;
      cellhistory[tid].coolinglist[i].upperlevel = upper;
      //if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d, coolingtype %d, i %d, low %d\n",contrib,oldcoolingsum,C,element,ion,level,cellhistory[tid].coolinglist[i].type,i,low);
      i += 1;
      //}
      
      
      /// fb creation of r-pkt
      /// free bound rates are calculated from the lower ion, but associated to the higher ion
      /// --------------------
      //upper = 0;
      //epsilon_upper = epsilon(element,ion+1,0);
      //E_threshold = epsilon_upper - epsilon_current; 
      //E_threshold = epsilon_trans;
      //printout("get_bfcooling(%d,%d,%d,%d) for histindex %d\n",element,ion,level,modelgridindex,histindex);
      //epsilon_upper = epsilon(element,ion+1,0);
      //E_threshold = epsilon_upper - epsilon_current;
      C = get_bfcooling(element,ion,level,modelgridindex);
      //printout("element %d, ion %d, level %d, T_e %g, alpha_E - alpha %g, bfcoolingcoeff %g\n",element,ion,level,T_e,C/E_threshold,C);
      //double interpolate_stimulated_recomb(int element, int ion, int level, double T);
      //C = get_bfcooling(element,ion,level,modelgridindex) + (nnnextionlevel*nne * W * interpolate_stimulated_bfcoolingcoeff(element,ion,level,T_R));
      //printout("nonE %g, E %g \n",interpolate_bfcoolingcoeff(element,ion,level,T_e), interpolate_stimulated_bfcoolingcoeff(element,ion,level,T_R));
      //C = get_bfcooling(element,ion,level,modelgridindex) + (nnnextionlevel*nne * (stimrecombestimator_E_save[pkt_ptr->where*nelements*maxion+element*maxion+ion]-stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+ion]) * E_threshold);
      //printout("element %d, ion %d, modified %g, usual %g, diff %g, nonE %g\n",element,ion,stimrecombestimator_E_save[pkt_ptr->where*nelements*maxion+element*maxion+ion],stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+ion],(stimrecombestimator_E_save[pkt_ptr->where*nelements*maxion+element*maxion+ion]-stimrecombestimator_save[pkt_ptr->where*nelements*maxion+element*maxion+ion])*E_threshold,interpolate_bfcoolingcoeff(element,ion,level,T_e));
      //double interpolate_spontrecombcoeff(int element, int ion, int level, double T);
      //printout("alpha_sp %g , alpha_st %g\n",interpolate_spontrecombcoeff(element,ion,level,T_e),interpolate_stimulated_recomb(element,ion,level,T_R));
      //epsilon_upper = epsilon(element,ion+1,0);
      //E_threshold = epsilon_upper - epsilon_current;
      //C = interpolate_spontrecombcoeff(element,ion,level,T_e) * E_threshold * nnnextionlevel*nne;
      //C_ion += C;
      //C_fb += C;
      //cellhistory[tid].coolinglist[i].contribution = C;
      contrib += C;
      cellhistory[tid].coolinglist[i].contribution = contrib;
      cellhistory[tid].coolinglist[i].type = COOLINGTYPE_FB;
      cellhistory[tid].coolinglist[i].element = element;
      cellhistory[tid].coolinglist[i].ion = ion;
      cellhistory[tid].coolinglist[i].level = level; 
      cellhistory[tid].coolinglist[i].upperlevel = 0;
      //if (contrib < oldcoolingsum) printout("contrib %g < oldcoolingsum %g: C%g, element %d, ion %d, level %d, coolingtype %d, i %d, low %d\n",contrib,oldcoolingsum,C,element,ion,level,cellhistory[tid].coolinglist[i].type,i,low);
      i += 1;
    }
  }
  
  //printout("i %d, high %d\n",i,high);

      
  //importantcoolingterms = ncoolingterms;
  //cellhistory[tid].totalcooling = contrib;
}




///****************************************************************************
double do_kpkt_bb(PKT *pkt_ptr, double t1, double t2)
/// Now routine to deal with a k-packet. Similar idea to do_gamma.
{
  double t_current;
  double zrand;
  
  double sample_planck(double T);
  void calculate_kappa_rpkt_cont(PKT *pkt_ptr, double t_current);
  double get_levelpop(int element, int ion, int level);
  void emitt_rpkt(PKT *pkt_ptr, double t_current);
  int get_element(int element);

  double coolingsum;  void update_cell(int cellnumber);

  int element,ion,upper;
  int i;
  
  

  //double nne = cell[pkt_ptr->where].nne ;
  double T_e = get_Te(cell[pkt_ptr->where].modelgridindex);
  t_current = t1;
  
  

  
  if (opacity_case == 4) 
  {
    /// Non grey opacity treatment, sample Planck in optically thick cells
    pkt_ptr->nu_cmf = sample_planck(T_e);
  }
  else
  {
    /// Grey opacity treatment, set a dummy value for nu_cmf
    pkt_ptr->nu_cmf = 1;
  }
  if (!isfinite(pkt_ptr->nu_cmf))
  {
    printout("[fatal] do_kpkt_bb: selected frequency not finite ... abort\n");
    abort();
  }
  /// and then emitt the packet randomly in the comoving frame
  emitt_rpkt(pkt_ptr,t_current);
  if (debuglevel == 2) printout("[debug] calculate_kappa_rpkt after kpkt to rpkt by ff\n");
    if (modelgrid[cell[pkt_ptr->where].modelgridindex].thick != 1){
        //printout("Abou to reset kpkt_bb\n");
        calculate_kappa_rpkt_cont(pkt_ptr,t_current);
        //printout("Abou to reset kpkt_bb\n");
    
    }
  pkt_ptr->next_trans = 0;      ///FLAG: transition history here not important, cont. process
  //if (tid == 0) k_stat_to_r_bb += 1;
  k_stat_to_r_bb += 1;
  pkt_ptr->interactions += 1;
  pkt_ptr->last_event = 6;
  pkt_ptr->emissiontype = -9999999;
  pkt_ptr->em_pos[0] = pkt_ptr->pos[0];
  pkt_ptr->em_pos[1] = pkt_ptr->pos[1];
  pkt_ptr->em_pos[2] = pkt_ptr->pos[2];
  pkt_ptr->em_time = t_current;
  pkt_ptr->nscatterings = 0;

  return(t_current);
}



///****************************************************************************
double sample_planck(double T)
/// returns a randomly chosen frequency according to the Planck 
/// distribution of temperature T
{
  double planck(double nu, double T);
  double zrand,zrand2;
  double nu,nu_peak;
  double B_peak;
  int endloop;
  int i;
  
  nu_peak = 5.879e10 * T;
  if (nu_peak > nu_max_r || nu_peak < nu_min_r) 
    printout("[warning] sample_planck: intensity peaks outside frequency range\n");
  B_peak = planck(nu_peak,T);
  endloop = 0;
  i = 0;
  while (endloop == 0)
  {
    i += 1;
    zrand = gsl_rng_uniform(rng);
    zrand2 = gsl_rng_uniform(rng);
    nu = nu_min_r + zrand*(nu_max_r - nu_min_r);
    if (zrand2*B_peak <= planck(nu,T)) endloop = 1;
    //printout("[debug] sample_planck: planck_sampling %d\n",i);
  }
  
  return nu;
}

///****************************************************************************
double planck(double nu, double T)
/// returns intensity for frequency nu and temperature T according
/// to the Planck distribution
{
  return TWOHOVERCLIGHTSQUARED*pow(nu,3) / (exp(HOVERKB*nu/T) - 1);
}






























///****************************************************************************
double do_kpkt(PKT *pkt_ptr, double t1, double t2, int nts)
/// Now routine to deal with a k-packet. Similar idea to do_gamma.
//{
//  double do_kpkt_bb(PKT *pkt_ptr, double t1, double t2);
//  return do_kpkt_bb(pkt_ptr, t1, t2);  
//}
{
  double t_current;
  double zrand;
  
  double sample_planck(double T);
  int get_continuumindex(int element, int ion, int level);
  double epsilon(int element, int ion, int level);
  //void calculate_kpkt_rates(int modelgridindex);
  void calculate_kpkt_rates_ion(int modelgridindex, int element, int ion, int low, double oldcoolingsum, int high);
  void calculate_kappa_rpkt_cont(PKT *pkt_ptr, double t_current);
  double get_levelpop(int element, int ion, int level);
  void emitt_rpkt(PKT *pkt_ptr, double t_current);
  int get_element(int element);
  int get_coolinglistoffset(int element, int ion);
  int get_ncoolingterms(int element, int ion);
  int call_estimators(PKT *pkt_ptr, double t_current, int realtype);
  void update_cell(int cellnumber);

  double coolingsum;
  int element,ion,level,upper;
  double nu_threshold;
  int i;
  int nions;
  
  
  gsl_integration_workspace *wsp;
  gslintegration_paras intparas;
  double bfcooling_integrand_gsl_2(double nu, void *paras);
  double alpha_sp_E_integrand_gsl(double nu, void *paras);
  gsl_function F_bfcooling;
  //F_bfcooling.function = &bfcooling_integrand_gsl_2;
  F_bfcooling.function = &alpha_sp_E_integrand_gsl;
  double intaccuracy = 1e-2;        /// Fractional accuracy of the integrator
  double error;
  double nu_lower,deltanu,sf,bfcooling_coeff,total_bfcooling_coeff,bfcooling_coeff_old,nuoffset;
  double calculate_sahafact(int element, int ion, int level, double T, double E_threshold);
  int ii;
  int ilow,low,high;
  double rndcool;
  double oldcoolingsum;
    
  int vflag = 0, mgi;
    
  int modelgridindex = cell[pkt_ptr->where].modelgridindex;
  
  /// Instead of doing the following, it is now made sure that kpkts in optically
  /// thick cells are treated by do_kpkt_bb in packet_prop
  /*
  /// Due to absorption of gamma-packets or decays of 52Fe and 52Mn pellets
  /// kpkts can be created in thick cells. Since in this case update_cell is 
  /// not called to calculate the level populations and reset the cell history 
  /// values it must be called here!
  void update_cell(int cellnumber);
  if (modelgrid[modelgridindex].thick == 1)
  {
    update_cell(mgi);
  }
  */
  
  /// don't calculate cooling rates after each cell crossings anylonger
  /// but only if we really get a kpkt and they hadn't been calculated already
  //if (cellhistory[tid].totalcooling == COOLING_UNDEFINED) 
  //printout("mgi %d, histindex %d\n",modelgridindex,histindex);
/*  int ondemand = 1;
  if (modelgrid[modelgridindex].totalcooling == COOLING_UNDEFINED) 
  {
    //printout("initial calculation of all cooling rates\n");
    coolingratecalccounter += 1;
    ondemand = 0;
    calculate_kpkt_rates(modelgridindex);
  }*/
  //printout("totalcooling %g\n",modelgrid[modelgridindex].totalcooling );

  //printout("[debug] do_kpkt: propagate k-pkt\n");
  //double cut = 0.99; //1.;
  

  double nne = get_nne(modelgridindex);
  double T_e = get_Te(modelgridindex);
  double deltat = 0.;
  if (nts < n_kpktdiffusion_timesteps) deltat = kpktdiffusion_timescale * time_step[nts].width;
  //double deltat = 1./(nne*1.02e-12*pow(T_e/1e4,0.843));
  //printout("kpkt diffusion time simple %g, advanced %g\n",deltat,1/(nne*1.02e-12*pow(T_e/1e4,0.843)));
  t_current = t1+deltat;

  if (t_current <= t2)
  {
    pkt_ptr->pos[0] = pkt_ptr->pos[0] * t_current / t1;
    pkt_ptr->pos[1] = pkt_ptr->pos[1] * t_current / t1;
    pkt_ptr->pos[2] = pkt_ptr->pos[2] * t_current / t1;
    wsp = gsl_integration_workspace_alloc(1000);
    
    /// Randomly select the occuring cooling process out of the important ones
    coolingsum = 0.;
    zrand = gsl_rng_uniform(rng);
    //if (debuglevel == 2) printout("do_kpkt: totalcooling %g, zrand %g, cut %g\n",cellhistory[tid].totalcooling,zrand,COOLINGCUT);
    //printout("do_kpkt: totalcooling %g, zrand %g, cut %g\n",cellhistory[tid].totalcooling,zrand,COOLINGCUT);
    /*for (i = 0; i < importantcoolingterms; i++)
    {
      coolingsum += cellhistory[tid].coolinglist[i].contribution;
      if (debuglevel == 2) printout("do_kpkt: loop i %d, coolingsum %g\n",i,coolingsum);
      //printout("do_kpkt: loop i %d, coolingsum %g, contr %g\n",i,coolingsum,cellhistory[tid].coolinglist[i].contribution);
      if (zrand*cellhistory[tid].totalcooling < coolingsum) break;
    }*/
    
    
    rndcool = zrand*modelgrid[modelgridindex].totalcooling;
    //printout("rndcool %g\n",rndcool);
    for (element = 0; element < nelements; element++)
    {
      nions = get_nions(element);
      for (ion = 0; ion < nions; ion++)
      {
        oldcoolingsum = coolingsum;
        coolingsum += modelgrid[modelgridindex].cooling[element].contrib[ion];
        //printout("element %d, ion %d, coolingsum %g\n",element,ion,coolingsum);
        if (coolingsum > rndcool) goto ENDSUM;
      }
    }
    
    ENDSUM: ;
    
  //  #ifdef DEBUG_ON
      if (element >= nelements || ion >= get_nions(element))
      {
        printout("do_kpkt: problem selecting a cooling process ... abort\n");
        printout("do_kpkt: tried to select element %d, ion %d\n",element,ion);
        printout("do_kpkt: totalcooling %g, coolingsum %g, rndcool %g\n",modelgrid[modelgridindex].totalcooling,coolingsum,rndcool);
        printout("do_kpkt: modelgridindex %d, cellno %d, nne %g\n",modelgridindex,pkt_ptr->where,get_nne(modelgridindex));
        for (element = 0; element < nelements; element++)
        {
          nions = get_nions(element);
          for (ion = 0; ion < nions; ion++)
          {
            printout("do_kpkt: element %d, ion %d, coolingcontr %g\n",element,ion,modelgrid[modelgridindex].cooling[element].contrib[ion]);
          }
        }
        abort();
      }
  //  #endif
    
    //debuglevel = 2;
    //printout("element %d, ion %d, coolingsum %g\n",element,ion,coolingsum);
    ilow = get_coolinglistoffset(element,ion);
    low = ilow;
    high = low+get_ncoolingterms(element,ion)-1;
    //printout("element %d, ion %d, low %d, high %d\n",element,ion,low,high);
    if (cellhistory[tid].coolinglist[low].contribution == COOLING_UNDEFINED)
    {
      //printout("calculate kpkt rates on demand\n");
      calculate_kpkt_rates_ion(modelgridindex,element,ion,low,oldcoolingsum,high);
    }
    //if (cellhistory[tid].coolinglist[high].contribution != coolingsum);
    //{
    //  printout("ondemand %d\n",ondemand);
    //  printout("element %d, ion %d, oldcoolingsum %g, coolingsum %g\n",element,ion,oldcoolingsum,coolingsum);
    //  printout("contrib[ilow] %g, contrib[high] %g\n",cellhistory[tid].coolinglist[ilow].contribution,cellhistory[tid].coolinglist[high].contribution);
    //}
    //low = 0;
    //high = importantcoolingterms - 1;
    //rndcool = zrand*cellhistory[tid].totalcooling;
    while (low <= high) 
    {
      i = (low + high) / 2;
      if (cellhistory[tid].coolinglist[i].contribution > rndcool)
      {
        if (i == ilow || cellhistory[tid].coolinglist[i-1].contribution < rndcool)
          break; /// found (1)
        else
          high = i - 1;
      }
      else
        low = i + 1;
      //else if (cellhistory[tid].coolinglist[i].contribution < rndcool)
      //  low = i + 1;
      //else
      //  break; /// found (2)
    }
    if (low > high) 
    {
      printout("do_kpkt: error occured while selecting a cooling channel: low %d, high %d, i %d, rndcool %g\n",low,high,i,rndcool);
      printout("element %d, ion %d, offset %d, terms %d, coolingsum %g\n",element,ion,get_coolinglistoffset(element,ion),get_ncoolingterms(element,ion),coolingsum);
      printout("oldcoolingsum %g, coolingsum %g\n",oldcoolingsum,coolingsum);
      //for (i=ilow-1; i<=ilow+get_ncoolingterms(element,ion); i++)
      //  printout("  i %d, contrib %g\n",i,cellhistory[tid].coolinglist[i].contribution);
      
      //printout("ondemand %d\n",ondemand);
        
      printout("lower %g, %g, %g\n",cellhistory[tid].coolinglist[get_coolinglistoffset(element,ion)-1].contribution,cellhistory[tid].coolinglist[get_coolinglistoffset(element,ion)].contribution,cellhistory[tid].coolinglist[get_coolinglistoffset(element,ion)+1].contribution);
      int finalpos = get_coolinglistoffset(element,ion)+get_ncoolingterms(element,ion)-1;
      printout("upper %g, %g, %g\n",cellhistory[tid].coolinglist[finalpos-1].contribution,cellhistory[tid].coolinglist[finalpos].contribution,cellhistory[tid].coolinglist[finalpos+1]);
    }
  
    
    if (debuglevel == 2) printout("do_kpkt: selected process %d, coolingsum %g, importantcoolingterms %d\n",i,coolingsum,importantcoolingterms);
    //printout("do_kpkt: selected process %d, coolingsum %g, importantcoolingterms %d, type %d\n",i,coolingsum,importantcoolingterms,cellhistory[tid].coolinglist[i].type);
    
    if (cellhistory[tid].coolinglist[i].type == COOLINGTYPE_FF)
    {
      /// The k-packet converts directly into a r-packet by free-free-emission.
      /// Need to select the r-packets frequency and a random direction in the
      /// co-moving frame.
      if (debuglevel == 2) printout("[debug] do_kpkt: k-pkt -> free-free\n");
      //kdecay.to_r += 1;
  
      /// Sample the packets comoving frame frequency according to paperII 5.4.3 eq.41
      //zrand = gsl_rng_uniform(rng);   /// delivers zrand in [0,1[
      //zrand = 1. - zrand;             /// make sure that log gets a zrand in ]0,1]
      zrand = gsl_rng_uniform_pos(rng);   /// delivers zrand in ]0,1[
      pkt_ptr->nu_cmf = -KB*T_e/H * log(zrand);
      //pkt_ptr->nu_cmf = 3.7474058e+14;
      
      if (!isfinite(pkt_ptr->nu_cmf))
      {
        printout("[fatal] ff cooling: selected frequency not finite ... abort\n");
        abort();
      }
      /// and then emitt the packet randomly in the comoving frame
      emitt_rpkt(pkt_ptr,t_current);
      if (debuglevel == 2) printout("[debug] calculate_kappa_rpkt after kpkt to rpkt by ff\n");
      
        
      pkt_ptr->next_trans = 0;      ///FLAG: transition history here not important, cont. process
      //if (tid == 0) k_stat_to_r_ff += 1;
      k_stat_to_r_ff += 1;
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 6;
      pkt_ptr->emissiontype = -9999999;
      pkt_ptr->em_pos[0] = pkt_ptr->pos[0];
      pkt_ptr->em_pos[1] = pkt_ptr->pos[1];
      pkt_ptr->em_pos[2] = pkt_ptr->pos[2];
      pkt_ptr->em_time = t_current;
      pkt_ptr->nscatterings = 0;
      #ifndef FORCE_LTE
        //kffcount[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif
       
        
        
      /* call the estimator routine - generate a virtual packet */
      #ifdef VPKT_ON
        realtype = 2 ;
       
        vflag = call_estimators(pkt_ptr, t_current, realtype);
 
        //printout("Abou to reset rpkt cell kpkt\n");
        calculate_kappa_rpkt_cont(pkt_ptr,t_current);
        //printout("Compl to reset rpkt cell kpkt\n");

    
      #else
        
        calculate_kappa_rpkt_cont(pkt_ptr,t_current);
        
      #endif
        
        
        
    }
    else if (cellhistory[tid].coolinglist[i].type == COOLINGTYPE_FB)
    {
      /// The k-packet converts directly into a r-packet by free-bound-emission.
      /// Need to select the r-packets frequency and a random direction in the
      /// co-moving frame.
      element = cellhistory[tid].coolinglist[i].element;
      ion = cellhistory[tid].coolinglist[i].ion;
      level = cellhistory[tid].coolinglist[i].level;
      upper = cellhistory[tid].coolinglist[i].upperlevel;
      nu_threshold = (epsilon(element,ion+1,upper) - epsilon(element,ion,level)) / H;
      
      #ifdef DEBUG_ON
        if (debuglevel == 2) printout("[debug] do_kpkt: k-pkt -> free-bound\n");
        if (debuglevel == 2) printout("[debug] do_kpkt: element  %d, ion %d, level %d, upper %d, nu_threshold %g\n",element,ion,level,upper,nu_threshold);
      #endif
  
      /// then randomly sample the packets frequency according to the continuums
      /// energy distribution and set some flags
      //zrand = gsl_rng_uniform(rng);   /// delivers zrand in [0,1[
      //zrand = 1. - zrand;   /// convert it to ]0,1]
      //pkt_ptr->nu_cmf = nu_threshold * (1 - KB*T_e/H/nu_threshold*log(zrand));  
      //pkt_ptr->nu_cmf = nu_threshold * (1+sqrt(1+(4*KB*T_e/H/nu_threshold)))/2 * (1 - KB*T_e/H/nu_threshold*log(zrand));  
      //pkt_ptr->nu_cmf = nu_threshold;
      
      //zrand = gsl_rng_uniform(rng);
      //if (zrand < 0.5)
      //{
      
      
        zrand = gsl_rng_uniform(rng);
        zrand = 1. - zrand;  /// Make sure that 0 < zrand <= 1
        mastate[tid].element = element;
        mastate[tid].ion = ion;
        mastate[tid].level = level;
        sf = calculate_sahafact(element,ion,level,T_e,nu_threshold*H);
        intparas.T = T_e;
        intparas.nu_edge = nu_threshold;   /// Global variable which passes the threshold to the integrator
        F_bfcooling.params = &intparas;
        deltanu = 9*nu_threshold/100;
        gsl_integration_qag(&F_bfcooling, nu_threshold, 10*nu_threshold, 0, intaccuracy, 1000, 6, wsp, &total_bfcooling_coeff, &error);
        bfcooling_coeff = total_bfcooling_coeff;
        for (ii= 0; ii < 100; ii++)
        {
          bfcooling_coeff_old = bfcooling_coeff;
          if (ii > 0)
          {
            nu_lower = nu_threshold + ii*deltanu;
            /// Spontaneous recombination and bf-cooling coefficient don't depend on the cutted radiation field
            gsl_integration_qag(&F_bfcooling, nu_lower, 10*nu_threshold, 0, intaccuracy, 1000, 6, wsp, &bfcooling_coeff, &error);
            //bfcooling_coeff *= FOURPI * sf; 
            //if (zrand > bfcooling_coeff/get_bfcooling(element,ion,level,pkt_ptr->where)) break;
          }
          //printout("zrand %g, bfcooling_coeff %g, total_bfcooling_coeff %g, nu_lower %g\n",zrand,bfcooling_coeff,total_bfcooling_coeff,nu_lower);
          if (zrand >= bfcooling_coeff/total_bfcooling_coeff) break;
        }
        if (ii==100)
        {
          printout("kpkt emitts bf-photon at upper limit\n");
          nu_lower = 10*nu_threshold; // + ii*deltanu;
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
      
      
      /*}
      else
      {
        ///Emitt like a BB
        pkt_ptr->nu_cmf = sample_planck(T_e);
      }*/
        
        //if (pkt_ptr->last_event == 4)      
        #ifndef FORCE_LTE
          //kbfcount[pkt_ptr->where] += pkt_ptr->e_cmf;
          //kbfcount_ion[pkt_ptr->where] += pkt_ptr->e_cmf*nu_threshold/pkt_ptr->nu_cmf;
        #endif
        
      //printout("nu_lower %g, nu_threshold %g, nu_left %g, nu_right %g\n",nu_lower,nu_threshold,nu_threshold+(ii-1)*deltanu,nu_threshold+(ii)*deltanu);
  
  //     if (element == 6) 
  //     {
  //       //printout("%g, %g, %g\n",pkt_ptr->e_cmf,nu_threshold,pkt_ptr->e_cmf/nu_threshold/H);
  //       cell[pkt_ptr->where].bfem[ion] += pkt_ptr->e_cmf/pkt_ptr->nu_cmf/H;
  //     }
  
      
      // /// Sample the packets comoving frame frequency according to paperII 4.2.2
      if (debuglevel == 2) printout("[debug] do_kpkt: pkt_ptr->nu_cmf %g\n",pkt_ptr->nu_cmf);
      //pkt_ptr->nu_cmf = 3.7474058e+14;
      if (!isfinite(pkt_ptr->nu_cmf))
      {
        printout("[fatal] rad deexcitation of MA: selected frequency not finite ... abort\n");
        abort();
      }
      /// and then emitt the packet randomly in the comoving frame
      emitt_rpkt(pkt_ptr,t_current);
        
      if (debuglevel == 2) printout("[debug] calculate_kappa_rpkt after kpkt to rpkt by fb\n");
        
        
      pkt_ptr->next_trans = 0;      ///FLAG: transition history here not important, cont. process
      //if (tid == 0) k_stat_to_r_fb += 1;
      k_stat_to_r_fb += 1;
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 7;
      pkt_ptr->emissiontype = get_continuumindex(element,ion,level);
      pkt_ptr->em_pos[0] = pkt_ptr->pos[0];
      pkt_ptr->em_pos[1] = pkt_ptr->pos[1];
      pkt_ptr->em_pos[2] = pkt_ptr->pos[2];
      pkt_ptr->em_time = t_current;
      pkt_ptr->nscatterings = 0;
        
        
      /* call the estimator routine - generate a virtual packet */
      #ifdef VPKT_ON
        realtype = 2 ;

        vflag = call_estimators(pkt_ptr, t_current, realtype);

        //printout("Abou to reset rpkt cell kpkt\n");
        calculate_kappa_rpkt_cont(pkt_ptr,t_current);
        //printout("Compl to reset rpkt cell kpkt\n");

      #else
        
        calculate_kappa_rpkt_cont(pkt_ptr,t_current);
      
      #endif
      
        
    }
    else if (cellhistory[tid].coolinglist[i].type == COOLINGTYPE_COLLEXC)
    {
      /// the k-packet activates a macro-atom due to collisional excitation
      if (debuglevel == 2) printout("[debug] do_kpkt: k-pkt -> collisional excitation of MA\n");
      element = cellhistory[tid].coolinglist[i].element;
      ion = cellhistory[tid].coolinglist[i].ion;
      upper = cellhistory[tid].coolinglist[i].upperlevel;
      mastate[tid].element = element;
      mastate[tid].ion = ion;
      mastate[tid].level = upper;
      mastate[tid].nnlevel = get_levelpop(element,ion,upper);
      mastate[tid].activatingline = -99;
      pkt_ptr->type = TYPE_MA;
      //if (tid == 0) ma_stat_activation_collexc += 1;
      ma_stat_activation_collexc += 1;
      //if (tid == 0) k_stat_to_ma_collexc += 1;
      k_stat_to_ma_collexc += 1;
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 8;
      #ifndef FORCE_LTE
        //maabs[pkt_ptr->where] += pkt_ptr->e_cmf;
        //kffcount[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif
    }
    else if (cellhistory[tid].coolinglist[i].type == COOLINGTYPE_COLLION)
    {
      /// the k-packet activates a macro-atom due to collisional ionisation
      if (debuglevel == 2) printout("[debug] do_kpkt: k-pkt -> collisional ionisation of MA\n");
      element = cellhistory[tid].coolinglist[i].element;
      ion = cellhistory[tid].coolinglist[i].ion+1;
      upper = cellhistory[tid].coolinglist[i].upperlevel;
      mastate[tid].element = element;
      mastate[tid].ion = ion;
      mastate[tid].level = upper;
      mastate[tid].nnlevel = get_levelpop(element,ion,upper);
      mastate[tid].activatingline = -99;
      pkt_ptr->type = TYPE_MA;
      //if (tid == 0) ma_stat_activation_collion += 1;
      ma_stat_activation_collion += 1;
      //if (tid == 0) k_stat_to_ma_collion += 1;
      k_stat_to_ma_collion += 1;
      pkt_ptr->interactions += 1;
      pkt_ptr->last_event = 9;
      #ifndef FORCE_LTE
        //maabs[pkt_ptr->where] += pkt_ptr->e_cmf;
        //kffcount[pkt_ptr->where] += pkt_ptr->e_cmf;
      #endif
    }
    else
    {
      printout("[fatal] do_kpkt: coolinglist.type mismatch\n");
      printout("[fatal] do_kpkt: zrand %g, modelgrid[modelgridindex].totalcooling %g, coolingsum %g, i %d\n",zrand,modelgrid[modelgridindex].totalcooling,  coolingsum,i);
      printout("[fatal] do_kpkt: COOLINGCUT %d, importantcoolingterms %d, tid %d\n",COOLINGCUT,importantcoolingterms,tid);
      printout("[fatal] do_kpkt: cellhistory[tid].coolinglist[i].type %d\n",cellhistory[tid].coolinglist[i].type);
      printout("[fatal] do_kpkt: pkt_ptr->where %d, mgi %d\n",pkt_ptr->where,modelgridindex);
      abort();
    }
    
    gsl_integration_workspace_free(wsp);
    return(t_current);
  }
  else
  {
    pkt_ptr->pos[0] = pkt_ptr->pos[0] * t2 / t1;
    pkt_ptr->pos[1] = pkt_ptr->pos[1] * t2 / t1;
    pkt_ptr->pos[2] = pkt_ptr->pos[2] * t2 / t1;
    return(PACKET_SAME);
  }
}





///****************************************************************************
/*int compare_coolinglistentry(const void *p1, const void *p2)
/// Helper function to sort the coolinglist by the strength of the 
/// individual cooling contributions.
{
  ionscoolinglist_t *a1, *a2;
  a1 = (ionscoolinglist_t *)(p1);
  a2 = (ionscoolinglist_t *)(p2);
  //printf("%d %d %d %d %g\n",a1->elementindex,a1->ionindex,a1->lowerlevelindex,a1->upperlevelindex,a1->nu);
  //printf("%d %d %d %d %g\n",a2->elementindex,a2->ionindex,a2->lowerlevelindex,a2->upperlevelindex,a2->nu);
  //printf("%g\n",a2->nu - a1->nu);
  if (a1->contribution - a2->contribution < 0)
    return 1;
  else if (a1->contribution - a2->contribution > 0)
    return -1;
  else
    return 0;
}*/



///***************************************************************************/
/*double get_bfcooling_direct(int element, int ion, int level, int cellnumber)
/// Returns the rate for bfheating. This can be called during packet propagation
/// or update_grid. Therefore we need to decide whether a cell history is
/// known or not.
{
  double calculate_sahafact(int element, int ion, int level, double T, double E_threshold);
  double get_groundlevelpop(int cellnumber, int element, int ion);
  double epsilon(int element, int ion, int level);
  double bfcooling;
  
  gsl_integration_workspace *wsp;
  gslintegration_paras intparas;
  double bfcooling_integrand_gsl(double nu, void *paras);
  gsl_function F_bfcooling;
  F_bfcooling.function = &bfcooling_integrand_gsl;
  double intaccuracy = 1e-2;        /// Fractional accuracy of the integrator
  double error;
  wsp = gsl_integration_workspace_alloc(1000);
  
  double T_e = cell[cellnumber].T_e;
  double nne = cell[cellnumber].nne;
  double nnionlevel = get_groundlevelpop(cellnumber,element,ion+1);
  //upper = cellhistory[tid].coolinglist[i].upperlevel;
  double nu_threshold = (epsilon(element,ion+1,0) - epsilon(element,ion,level)) / H;

  mastate[tid].element = element;
  mastate[tid].ion = ion;
  mastate[tid].level = level;
  intparas.T = T_e;
  intparas.nu_edge = nu_threshold;   /// Global variable which passes the threshold to the integrator
  F_bfcooling.params = &intparas;
  gsl_integration_qag(&F_bfcooling, nu_threshold, 10*nu_threshold, 0, intaccuracy, 1000, 6, wsp, &bfcooling, &error);
  bfcooling *= nnionlevel*nne*4*PI*calculate_sahafact(element,ion,level,T_e,nu_threshold*H);
  
  gsl_integration_workspace_free(wsp);
  return bfcooling;
}*/


int get_coolinglistoffset(int element, int ion)
{
  return elements[element].ions[ion].coolingoffset;
}

int get_ncoolingterms(int element, int ion)
{
  return elements[element].ions[ion].ncoolingterms;
}
