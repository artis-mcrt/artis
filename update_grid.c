#include "sn3d.h"
#include "vpkt.h"

/** Subroutine to update the matter quantities in the grid cells at the start
   of the new timestep. */

//int update_grid(int m, int nstart, int nblock)
int update_grid(int m, int my_rank, int nstart, int nblock, int titer)
    /// m timestep
{
  void calculate_kpkt_rates(int modelgridindex);
  int get_coolinglistoffset(int element, int ion);
  double vel;
  double get_bfheatingcoeff_ana(int element, int ion, int level, int modelgridindex);
  double hhelper;
  int jj;
  double calculate_exclevelpop(int modelgridindex, int element, int ion, int level);
  double get_corrphotoioncoeff(int element, int ion, int level, int modelgridindex);
  double get_corrphotoioncoeff_ana(int element, int ion, int level, int modelgridindex);
  double renormcoeff;
  double get_Gamma_phys(int modelgridindex, int element, int ion);
  double get_Gamma(int modelgridindex, int element, int ion);
  double bfheating_est,ffheating_est,ffcooling_ana,bfheating_ana,bfcooling_ana,photoabs_est,photoabs_twiddle_est,radrecomb_ana,thermalratio;
  int counter;
  double old,new;
  
  //void dilute_pops(int cellnumber, int nts);
  void calculate_cooling_rates(int modelgridindex);
  void calculate_oldheating_rates(int modelgridindex);
  double get_radrecomb(int element, int ion, int level, int modelgridindex);
  double get_thermalratio(int element, int ion, int level, int modelgridindex);
  double get_adcool(int modelgridindex, int nts);
  double get_ffcooling(int element, int ion, int modelgridindex);
  //double get_twiddleabs(int element, int ion, int level, int modelgridindex);
  double get_bfheating(int element, int ion, int level, int modelgridindex);
  double get_photoabs(int element, int ion, int modelgridindex);
  double get_bfcooling(int element, int ion, int level, int modelgridindex);
//  double interpolate_photoioncoeff_below(int element, int ion, int level, double T);
//  double interpolate_photoioncoeff_above(int element, int ion, int level, double T);
  //double interpolate_zeta(int element, int ion, double T);
  //double gamma_lte,zeta;
  
  double get_groundlevelpop(int modelgridindex, int element, int ion);
  void get_radfield_params(double J, double nuJ, int modelgridindex, double *T_J, double *T_R, double *W);
  double ionstagepop(int modelgridindex, int element, int ion);
  double call_T_e_finder(int modelgridindex, double t_current, int tb_info, double T_min, double T_max);
  void precalculate_partfuncts(int modelgridindex);
  double calculate_populations(int modelgridindex, int first_nonempty_cell);
  //double ionstagepop(int cellnumber, int element, int ion);
  void calculate_cooling_rates(int modelgridindex);
  void calculate_heating_rates(int modelgridindex);
  //void check_interpolation(double T_min, double T_max);
  double find_T_e(double T_e, void *paras);
  //double mintemp_solution_f(double T, void *paras);
  //double maxtemp_solution_f(double T, void *paras);
  //double calculate_ltepartfunct(int element, int ion, double T);
  //float calculate_groundlevelpop(int element, int ion, double T, int cellnumber, double nne, double nnnextion);
  void update_abundances(int modelgridindex, double t_current);
  double get_abundance(int modelgridindex, int element);

  int search_cellhistory(int cellnumber);
  void determine_kpkt_cuts(int cellnumber);

  int get_element();

  /// only needed if all level populations should be printed to the output-file
  //double pop, excitedlevelpops;
  
  int n,ncl;
  //int samplecell;
  double deltarho_old,deltaT_R_old,deltaT_e_old;
  double deltarho,deltaT_R,deltaT_e;
  int rhoindex,T_Rindex,T_eindex;
  double trat, tratmid;
  double cell_len_scale;
  double check1, check2; //MK
  int element,ion,level;
  double t_current,t_previous;
  double rho,T_R,T_e,T_J,W;//,W_D;
  double mps[MTHREADS];  /// Thread private substitution of max_path_step. Its minimum is
                         /// assigned to max_path_step after the parallel update_grid finished.
  
  double nne,nntot;
  double abundance,nnelement,nne_tot;
  double pop;
  double compton_optical_depth;
  double grey_optical_depth,grey_optical_deptha;
  
  double dt_elapsed,dt_forward;
  double Gamma,Gamma_old;
  
  int first_nonempty_cell = -1000;
  int assoc_cells;
  //int lastelement = nelements-1;
  
  //temprange_paras paras2;
  //int status;
  //int iter;
  int i;//,iter2;
  //double x_0,x_lo,x_hi;
  double T_lo,T_hi,T_min,T_max,T_e_old;
  //int kpkt_cuts_determined = 0;
  double deltat,deltaV,deltaV1;
  int nions,nlevels;
  double radial_pos;
  //char tempfilename[100],ionfractfilename[100],Gammafilename[100],Alphaspfilename[100];
  //FILE *temperature_file,*ionfract_file,*corrphotoion_file,*thermal_file,*Gamma_file,*Alpha_file,*bfcount_file;
  //FILE *gammaest_file,*gammaana_file;
  double interpolate_ions_spontrecombcoeff(int element, int ion, double T);
  
  
  /*if (&modelgrid[96].composition[0] == NULL)
  {
   printout("fatal error in ts %d abort\n",m);
   abort(); 
  }*/
  
  int tb_info = 0;
  //printout("[debug] update_grid: starting update for timestep %d...\n",m);
  check1=check2=0.0; /// MK
  trat = time_step[m].start / tmin;
  tratmid = time_step[m].mid / tmin;

  for (i=0; i < nthreads; i++)
  {
    mps[i] = 1.e35;
  }
  
  /// For debug temperature output
  /*
  sprintf(tempfilename,"d%d_thermal_%.4d.out",m,my_rank);
  //sprintf(tempfilename,"t%d_temperature_%d-%d.out",m,nstart,nstart+nblock-1);
  if ((thermal_file = fopen(tempfilename, "w")) == NULL)
  {
    printf("Cannot open %s.\n",tempfilename);
    exit(0);
  }
  setvbuf(thermal_file, NULL, _IOLBF, 1);
  */
  
  //printout("[debug] update_grid: time before initialisation of heating file %d\n",time(NULL));
  //#ifndef FORCE_LTE
  //  sprintf(filename,"h%d-%d_heating_%.4d.out",m,titer,my_rank);
  //  if ((heating_file = fopen(filename, "w")) == NULL)
  //  {
  //    printf("Cannot open %s.\n",filename);
  //    exit(0);
  //  }
  //  setvbuf(heating_file, NULL, _IOLBF, 1);
  //#endif
  //printout("[debug] update_grid: heating file initialised %d\n",time(NULL));
  
  ///Calculate the critical opacity at which opacity_case 3 switches from a
  ///regime proportional to the density to a regime independent of the density
  ///This is done by solving for tau_sobolev == 1
  ///tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/MNI56 * 3000e-8 * time_step[m].mid;
  rho_crit = ME*CLIGHT*MNI56 / (PI*QE*QE * rho_crit_para * 3000e-8 * time_step[m].mid);
  printout("update_grid: rho_crit = %g\n", rho_crit);
  //printf("time %ld\n",time(NULL));
  
  /// Needed to update abundances of radioactive isotopes.
  t_current = time_step[m].start;
  if (m == 0) 
  {
    t_previous = 0.;
    deltaV = pow(wid_init*trat,3);  /// volume of grid cell: current or previous cell size???????????????????
  }
  else 
  {
    t_previous = time_step[m-1].start;
    deltaV = pow(wid_init * time_step[m-1].mid/tmin,3);  /// volume of grid cell: current or previous cell size???????????????????
    if (m == 1)
    {
      dt_elapsed = (log(time_step[m-1].mid) - log(time_step[m-1].start));
      dt_forward = (log(time_step[m].mid) - log(time_step[m-1].mid));
    }
    else
    {
      dt_elapsed = (log(time_step[m-1].mid) - log(time_step[m-2].mid));
      dt_forward = (log(time_step[m].mid) - log(time_step[m-1].mid));
    }
  }
  
  /// and for the volume estimators
  deltat = t_current-t_previous;  /// length of previous timestep
  
  if (titer == 0)
  {
    if (m == 0)
    {
      /// Set these values, but they will not be used
      deltat = time_step[m].width; 
      deltaV = pow(wid_init*tratmid,3);
    }
    else
    {
      deltat = time_step[m-1].width; 
      deltaV = pow(wid_init * time_step[m-1].mid/tmin,3);
    }
  }
  else
  {
    deltat = time_step[m].width; 
    deltaV = pow(wid_init*tratmid,3);
  }
  
  printout("timestep %d, titer %d\n",m,titer);
  printout("deltaV %g, deltat %g\n",deltaV,deltat);
  
  
  /*
  FILE *photoion_file;
  FILE *bf_file;
  char photoion_filename[100],bf_filename[100];
  if (m != itstep)
  {
    sprintf(photoion_filename,"photoion_%.2d.out",m);
    if ((photoion_file = fopen(photoion_filename, "w")) == NULL)
    {
      printf("Cannot open %s.\n",photoion_filename);
      exit(0);
    }
    setvbuf(photoion_file, NULL, _IOLBF, 1);
    sprintf(bf_filename,"bf_%.2d.out",m);
    if ((bf_file = fopen(bf_filename, "w")) == NULL)
    {
      printf("Cannot open %s.\n",bf_filename);
      exit(0);
    }
    setvbuf(bf_file, NULL, _IOLBF, 1);

    for (ncl = 0; ncl < nblock; ncl++)
    {
      n = nonemptycells[my_rank+ncl*nprocs];
      fprintf(photoion_file,"%d ",n);
      fprintf(bf_file,"%d ",n);
      for (i=0; i < 9; i++)
      {
        fprintf(photoion_file,"%g %g ",cell[n].photoion[i],cell[n].radrecomb[i]);
        fprintf(bf_file,"%g %g ",cell[n].bfabs[i],cell[n].bfem[i]);
      }
      fprintf(photoion_file,"\n");
      fprintf(bf_file,"\n");
    }
    fclose(photoion_file);
    fclose(bf_file);
  }
  */

  
  /// Compare values of rate coefficients obtained by the tabulated values to
  /// directly integrated values.
  //check_interpolation(T_min,T_max);
  
//   #ifndef FORCE_LTE
//     //sprintf(Alphaspfilename,"bfcount%.2d-%.2d_%.4d.out",m,titer,my_rank);
//     sprintf(Alphaspfilename,"bfcount%d-%d_%.4d.out",m,titer,my_rank);
//     if ((bfcount_file = fopen(Alphaspfilename, "w")) == NULL)
//     {
//       printf("Cannot open %s.\n",Alphaspfilename);
//       exit(0);
//     }
//     setvbuf(bfcount_file, NULL, _IOLBF, 1);
//   #endif

  
  #ifdef _OPENMP
    #pragma omp parallel private(element,nions,ion,nlevels,level) 
    //copyin(nuJ,J,rhosum,T_Rsum,T_esum,Wsum,associatedcells)
    {
  #endif
      /// Do not use values which are saved in the cellhistory within update_grid
      /// and daughter routines (THREADPRIVATE VARIABLE, THEREFORE HERE!)
      use_cellhist = -1;
      
      /// All entries of the cellhistory stack must be flagged as empty at the
      /// onset of the new timestep.
      cellhistory[tid].cellnumber = -99;
      //cellhistory[tid].totalcooling = COOLING_UNDEFINED;
      for (element = 0; element < nelements; element++)
      {
        nions = get_nions(element);
        for (ion = 0; ion < nions; ion++)
        {
          cellhistory[tid].coolinglist[get_coolinglistoffset(element,ion)].contribution = COOLING_UNDEFINED;
          nlevels = get_nlevels(element,ion);
          for (level = 0; level < nlevels; level++)
          {
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].population = -99.;
            
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].sahafact = -99.;
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].spontaneousrecombrate = -99.;
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].bfcooling = -99.;
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].corrphotoioncoeff = -99.;
            
            /// This is the only flag needed for all of the following MA stuff!
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].col_deexc = -99.;
          }
        }
      }
      
      /// Reset histindex to unknown
      //printout("thread%d _ update_grid histindex %d\n",tid,histindex);
      //histindex = -99;
      //printout("thread%d _ update_grid histindex %d\n",tid,histindex);
      
      /// Updating cell information
      #ifdef _OPENMP
        #pragma omp for private(n,assoc_cells,radial_pos,tb_info,rho,element,nions,ion,nntot,nne,compton_optical_depth,grey_optical_deptha,grey_optical_depth,T_R,T_J,T_e,W,nlevels,Gamma,level,T_e_old,cell_len_scale) reduction(+:check1,check2) schedule(dynamic)
        //T_D,W_D,nne,deltarho_old,deltaT_R_old,deltaT_e_old, deltarho,deltaT_R,deltaT_e,i,rhoindex,T_Rindex,T_eindex,ncl)
      #endif
      for (n = 0; n < npts_model; n++)
      //for (n = nstart; n < nstart+nblock; n++)
      //for (ncl = 0; ncl < nblock; ncl++)
      //for (ncl = nstart; ncl < nstart+nblock; ncl++)
      //for (n = nstart; n < nstart+nblock; n++)
      {
        if (n >= nstart && n < nstart+nblock)
        {
          /// Check if this task should work on the current model grid cell.
          /// If yes, do update grid
          assoc_cells = modelgrid[n].associated_cells;
          radial_pos = modelgrid[n].initial_radial_pos*tratmid/assoc_cells;
          if (assoc_cells > 0)
          {
            //cellnumber = modelgrid[n].cellnumber;
            if (my_rank % nblock == n) tb_info = 1;
            //if (my_rank % nblock == ncl) tb_info = 1;
            else tb_info = 0;
            /// Update current mass density of cell
            //n = nonemptycells[my_rank+ncl*nprocs];
            printout("[info] update_grid: working on cell %d ...\n",n);
            //n = nonemptycells[ncl];
            //printout("[debug] update_grid: ncl %d is %d non-empty cell updating grid cell %d ... T_e %g, rho %g\n",ncl,my_rank+ncl*nprocs,n,cell[n].T_e,cell[n].rho);
            modelgrid[n].rho = modelgrid[n].rhoinit / pow(tratmid,3);
            rho = modelgrid[n].rho;
            //cell[n].rho = cell[n].rho_init / pow(tratmid,3);
            //rho = cell[n].rho;
            /// This is done outside update grid now
            //modelgrid[n].totalcooling = COOLING_UNDEFINED;
            
            
            if (opacity_case == 4)
            {
              /// Initialise shortcuts to previous values.
              /*T_R = cell[n].T_R;
              T_e = cell[n].T_e;
              T_e_old = T_e;
              W = cell[n].W;
              nne = cell[n].nne;*/
              
              /// Update abundances of radioactive isotopes
              //printout("call update abundances for timestep %d in model cell %d\n",m,n);
              update_abundances(n, time_step[m].mid);
        
              /// For timestep 0 we calculate the level populations straight forward wihout
              /// applying any temperature correction
              if (m-itstep == 0 && titer == 0)
              {
                /// Determine renormalisation factor for corrected photoionsiation cross-sections
                #ifndef FORCE_LTE
                  if (!continue_simulation)
                  {
                    for (element = 0; element < nelements; element++)
                    {
                      nions = get_nions(element);
                      for (ion = 0; ion < nions; ion++)
                      {
                        corrphotoionrenorm[n*nelements*maxion+element*maxion+ion] = 1.;
                      }
                    }
                  }
                #endif
                
		/// W == 1 indicates that this modelgrid cell was treated grey in the
		/// last timestep. Therefore it has no valid Gamma estimators and must
		/// be treated in LTE at restart. 
		if (modelgrid[n].thick != 1 && get_W(n) == 1)
		{
		  modelgrid[n].thick = 1;
		  printout("force modelgrid cell to be grey at restart\n");
		}
                printout("initial iteration %d\n",initial_iteration);
                printout("modelgrid.thick: %d\n",modelgrid[n].thick);
                precalculate_partfuncts(n);
                //printout("abundance in cell %d is %g\n",n,cell[n].composition[0].abundance);
                nntot = calculate_populations(n,0);
                nne = get_nne(n);
                compton_optical_depth = SIGMA_T*nne*wid_init*tratmid;
                        
                grey_optical_deptha = get_kappagrey(n)*get_rho(n)*wid_init*tratmid;
                printout("cell %d, compton optical depth %g, grey optical depth %g\n",n,compton_optical_depth,grey_optical_deptha);
                grey_optical_depth = get_kappagrey(n)*get_rho(n)*(rmax*tratmid-radial_pos);
                printout("pos %g, distance %g, tau_dist %g\n",radial_pos,rmax*tratmid-radial_pos,grey_optical_depth);
                modelgrid[n].grey_depth = grey_optical_depth;
                    
  //               grey_optical_depth = SIGMA_T*nne*wid_init*tratmid;
                if (grey_optical_depth > cell_is_optically_thick && m < n_grey_timesteps)
                {
                  printout("cell %d is treated in grey approximation (tau %g)\n",n,grey_optical_depth);
                  modelgrid[n].thick = 1;
                }
                else if (grey_optical_depth > cell_is_optically_thick_vpkt )
                {
                    modelgrid[n].thick = 2;
                }
                else modelgrid[n].thick = 0;
              }
              
              /// For all other timesteps temperature corrections have to be applied
              else
              {
                /// Energy deposition calculation is independent of LTE/NLTE so do it here
                energy_deposition[n] *= 1./(deltaV*deltat)/nprocs/assoc_cells;
		
                /// we have to calculate the electron density 
                /// and all the level populations
                /// Normalise estimators and make sure that they are finite.
                /// Then update T_R and W using the estimators.
                /// (This could in principle also be done for empty cells)
    
                #ifdef FORCE_LTE
                  /// LTE version of the code
                  if (!isfinite(J[n]))
                  {
                    printout("[fatal] update_grid: non finite estimator before normalisation ... abort\n");
                    abort();
                  }
                  J[n] *= ONEOVER4PI/(deltaV*deltat)/nprocs/assoc_cells;
                  
                  #ifdef DO_TITER
                    if (J_reduced_save[n] >= 0)
                    {
                      J[n] = (J[n]+J_reduced_save[n])/2.;
                    }
                    J_reduced_save[n] = J[n];
                  #endif
                  
                  T_R = pow(PI/STEBO*(J[n]),1./4.);
                  if (isfinite(T_R))
                  {
                    /// Make sure that T is in the allowed temperature range.
                    if (T_R > MAXTEMP) T_R = MAXTEMP;
                    if (T_R < MINTEMP) T_R = MINTEMP;
                  }
                  else
                  {
                    /// keep old value of T_R
                    printout("[warning] update_grid: T_R estimator infinite in cell %d, use value of last timestep\n",n);
                    T_R = modelgrid[n].TR;
                  }
                  T_J = T_R;
                  T_e = T_R;
                  W = 1.;
                  
                  modelgrid[n].TR = T_R;
                  modelgrid[n].Te = T_e;
                  modelgrid[n].TJ = T_J;
                  modelgrid[n].W = W;
                
                  /// These don't depend on T_e, therefore take them out of the T_e iteration
                  precalculate_partfuncts(n);
                  nntot = calculate_populations(n,0);
                  nne = get_nne(n);
                  compton_optical_depth = SIGMA_T*nne*wid_init*tratmid;
                  grey_optical_deptha = get_kappagrey(n)*get_rho(n)*wid_init*tratmid;
                  printout("cell %d, compton optical depth %g, grey optical depth %g\n",n,compton_optical_depth,grey_optical_deptha);
                  grey_optical_depth = get_kappagrey(n)*get_rho(n)*(rmax*tratmid-radial_pos);
                  printout("pos %g, distance %g, tau_dist %g\n",radial_pos,rmax*tratmid-radial_pos,grey_optical_depth);
                  modelgrid[n].grey_depth = grey_optical_depth;
  
  //                 grey_optical_depth = SIGMA_T*nne*wid_init*tratmid;
                  if (grey_optical_depth > cell_is_optically_thick && m < n_grey_timesteps)
                  {
                    printout("cell %d is treated in grey approximation (tau %g)\n",n,grey_optical_depth);
                    modelgrid[n].thick = 1;
                  }
                  else if (grey_optical_depth > cell_is_optically_thick_vpkt )
                  {
                      modelgrid[n].thick = 2;
                  }
                  else modelgrid[n].thick = 0;
                
                #else 
                  /// non-LTE version of the code
                  if (initial_iteration == 1 || modelgrid[n].thick == 1)
                  {
                    /// LTE version of the code
                    if (!isfinite(J[n]))
                    {
                      printout("[fatal] update_grid: non finite estimator before normalisation ... abort\n");
                      abort();
                    }
                    J[n] *= ONEOVER4PI/(deltaV*deltat)/nprocs/assoc_cells;
                    
                    #ifdef DO_TITER
                      if (J_reduced_save[n] >= 0)
                      {
                        J[n] = (J[n]+J_reduced_save[n])/2.;
                      }
                      J_reduced_save[n] = J[n];
                    #endif
                    
                    T_R = pow(PI/STEBO*(J[n]),1./4.);
                    if (isfinite(T_R))
                    {
                      /// Make sure that T is in the allowed temperature range.
                      if (T_R > MAXTEMP) T_R = MAXTEMP;
                      if (T_R < MINTEMP) T_R = MINTEMP;
                    }
                    else
                    {
                      /// keep old value of T_R
                      printout("[warning] update_grid: T_R estimator infinite in cell %d, use value of last timestep\n",n);
                      T_R = modelgrid[n].TR;
                    }
                    T_J = T_R;
                    T_e = T_R;
                    W = 1.;
                    
                    modelgrid[n].TR = T_R;
                    modelgrid[n].Te = T_e;
                    modelgrid[n].TJ = T_J;
                    modelgrid[n].W = W;
                  
                    /// These don't depend on T_e, therefore take them out of the T_e iteration
                    for (element = 0; element < nelements; element++)
                    {
                      nions = get_nions(element);
                      for (ion = 0; ion < nions-1; ion++)
                      {
                        corrphotoionrenorm[n*nelements*maxion+element*maxion+ion] = 1.;
                      }
                    }
                    precalculate_partfuncts(n);
                    nntot = calculate_populations(n,0);
                    nne = get_nne(n);
                    compton_optical_depth = SIGMA_T*nne*wid_init*tratmid;
                    grey_optical_deptha = get_kappagrey(n)*get_rho(n)*wid_init*tratmid;
                    printout("cell %d, compton optical depth %g, grey optical depth %g\n",n,compton_optical_depth,grey_optical_deptha);
                    grey_optical_depth = get_kappagrey(n)*get_rho(n)*(rmax*tratmid-radial_pos);
                    printout("pos %g, distance %g, tau_dist %g\n",radial_pos,rmax*tratmid-radial_pos,grey_optical_depth);
                    modelgrid[n].grey_depth = grey_optical_depth;
  
  //                   grey_optical_depth = SIGMA_T*nne*wid_init*tratmid;
                    if (grey_optical_depth > cell_is_optically_thick && m < n_grey_timesteps)
                    {
                      printout("cell %d is treated in grey approximation (tau %g)\n",n,grey_optical_depth);
                      modelgrid[n].thick = 1;
                    }
                    else if (grey_optical_depth > cell_is_optically_thick_vpkt )
                    {
                        modelgrid[n].thick = 2;
                    }
                    else modelgrid[n].thick = 0;
                  }
                  
                  else
                  {
                    /// Calculate estimators
                    if (!isfinite(nuJ[n]) || !isfinite(J[n]))
                    {
                      printout("[fatal] update_grid: non isfinite estimator before normalisation ... abort\n");
                      abort();
                    }
                    
                    J[n] *= ONEOVER4PI/(deltaV*deltat)/nprocs/assoc_cells;
                    nuJ[n] *= ONEOVER4PI/(deltaV*deltat)/nprocs/assoc_cells;
                    ffheatingestimator[n] *= 1/(deltaV*deltat)/nprocs/assoc_cells;
                    colheatingestimator[n] *= 1/(deltaV*deltat)/nprocs/assoc_cells;
                    
                    #ifdef DO_TITER
                      if (J_reduced_save[n] >= 0)
                      {
                        J[n] = (J[n]+J_reduced_save[n])/2.;
                      }
                      J_reduced_save[n] = J[n];
                      if (nuJ_reduced_save[n] >= 0)
                      {
                        nuJ[n] = (nuJ[n]+nuJ_reduced_save[n])/2.;
                      }
                      nuJ_reduced_save[n] = nuJ[n];
                      if (ffheatingestimator_save[n] >= 0)
                      {
                        ffheatingestimator[n] = (ffheatingestimator[n]+ffheatingestimator_save[n])/2.;
                      }
                      ffheatingestimator_save[n] = ffheatingestimator[n] ;
                      if (colheatingestimator_save[n] >= 0)
                      {
                        colheatingestimator[n] = (colheatingestimator[n]+colheatingestimator_save[n])/2.;
                      }
                      colheatingestimator_save[n] = colheatingestimator[n];
                    #endif
                    
                    for (element = 0; element < nelements; element++)
                    {
                      nions = get_nions(element);
                      for (ion = 0; ion < nions-1; ion++)
                      {
                        //printout("mgi %d, element %d, ion %d, gammaest %g\n",n,element,ion,gammaestimator[n*nelements*maxion+element*maxion+ion]);
                        gammaestimator[n*nelements*maxion+element*maxion+ion] *= 1/(deltaV*deltat)/H/nprocs/assoc_cells;
                        //printout("mgi %d, element %d, ion %d, gammaest %g\n",n,element,ion,gammaestimator[n*nelements*maxion+element*maxion+ion]);
                        #ifdef DO_TITER
                          if (gammaestimator_save[n*nelements*maxion+element*maxion+ion] >= 0)
                          {
                            gammaestimator[n*nelements*maxion+element*maxion+ion] = (gammaestimator[n*nelements*maxion+element*maxion+ion]+gammaestimator_save[n*nelements*maxion+element*maxion+ion])/2.;
                          }
                          gammaestimator_save[n*nelements*maxion+element*maxion+ion] = gammaestimator[n*nelements*maxion+element*maxion+ion];
                        #endif
                        corrphotoionrenorm[n*nelements*maxion+element*maxion+ion] = gammaestimator[n*nelements*maxion+element*maxion+ion]/get_corrphotoioncoeff_ana(element,ion,0,n);
                      
                      /// 2012-01-11. These loops should terminate here to precalculate *ALL* corrphotoionrenorm values
                      /// so that the values are known when required by the call to get_corrphotoioncoeff in the following
                      /// loops. Otherwise get_corrphotoioncoeff tries to renormalize by the closest corrphotoionrenorm 
                      /// in frequency space which can lead to zero contributions to the total photoionsation rate!
                      }
                    }
                    
                    /// Then reopen the same loops again.
                    for (element = 0; element < nelements; element++)
                    {
                      nions = get_nions(element);
                      for (ion = 0; ion < nions-1; ion++)
                      {
                        /// Reuse the gammaestimator array as temporary storage of the Gamma values during
                        /// the remaining part of the update_grid phase. Afterwards it is reset to record
                        /// the next timesteps gamma estimators.
                        //nlevels = get_nlevels(element,ion);
                        //nlevels = get_ionisinglevels(element,ion);
                        nlevels = get_bfcontinua(element,ion);
                        Gamma = 0.;
                        for (level = 0; level < nlevels; level++)
                        {
                          Gamma += calculate_exclevelpop(n,element,ion,level)*get_corrphotoioncoeff(element,ion,level,n);
                          //printout("mgi %d, element %d, ion %d, level %d, pop %g, corrphotoion %g\n",n,element,ion,level,calculate_exclevelpop(n,element,ion,level),get_corrphotoioncoeff(element,ion,level,n));
                        }
                        Gamma /= get_groundlevelpop(n, element, ion);
                        gammaestimator[n*nelements*maxion+element*maxion+ion] = Gamma;
                        //printout("mgi %d, element %d, ion %d, Gamma %g\n",n,element,ion,Gamma);
                        
                        bfheatingestimator[n*nelements*maxion+element*maxion+ion] *= 1/(deltaV*deltat)/nprocs/assoc_cells;
                        #ifdef DO_TITER
                          if (bfheatingestimator_save[n*nelements*maxion+element*maxion+ion] >= 0)
                          {
                            bfheatingestimator[n*nelements*maxion+element*maxion+ion] = (bfheatingestimator[n*nelements*maxion+element*maxion+ion]+bfheatingestimator_save[n*nelements*maxion+element*maxion+ion])/2.;
                          }
                          bfheatingestimator_save[n*nelements*maxion+element*maxion+ion] = bfheatingestimator[n*nelements*maxion+element*maxion+ion]; 
                        #endif
                        /// Now convert bfheatingestimator into the bfheating renormalisation coefficient used in get_bfheating
                        /// in the remaining part of update_grid. Later on it's reset and new contributions are added up.
                        bfheatingestimator[n*nelements*maxion+element*maxion+ion] = bfheatingestimator[n*nelements*maxion+element*maxion+ion]/get_bfheatingcoeff_ana(element,ion,0,n);
                        //printout("cell %d element %d ion %d bfheatingestimator %g\n",n,element,ion,bfheatingestimator[n*nelements*maxion+element*maxion+ion]);
                      }
                    }
                    
                    /// Get radiation field parameters out of the estimators
                    get_radfield_params(J[n],nuJ[n],n,&T_J,&T_R,&W);
                    modelgrid[n].TJ = T_J;
                    modelgrid[n].TR = T_R;
                    modelgrid[n].W = W;
                  
                    /// These don't depend on T_e, therefore take them out of the T_e iteration
                    precalculate_partfuncts(n);
                  
                    /// Find T_e as solution for thermal balance
                    T_e_old = get_Te(n);
                    if (titer == 0)
                      T_e = call_T_e_finder(n,time_step[m-1].mid,tb_info,MINTEMP,MAXTEMP);
                    else
                      T_e = call_T_e_finder(n,time_step[m].mid,tb_info,MINTEMP,MAXTEMP);
                    if (T_e > 2.*T_e_old) 
                    {
                      T_e = 2.*T_e_old;
                      printout("use T_e damping in cell %d\n",n);
                      if (T_e > MAXTEMP) T_e = MAXTEMP;
                    }
                    else if (T_e < 0.5*T_e_old) 
                    {
                      T_e = 0.5*T_e_old;
                      printout("use T_e damping in cell %d\n",n);
                      if (T_e < MINTEMP) T_e = MINTEMP;
                    }
                    //T_e = T_J;
                    set_Te(n,T_e);
      
                    /// Store population values to the grid
                    nntot = calculate_populations(n,0);
                    //calculate_cooling_rates(n);
                    //calculate_heating_rates(n);
                    nne = get_nne(n);
                    compton_optical_depth = SIGMA_T*nne*wid_init*tratmid;
                    grey_optical_deptha = get_kappagrey(n)*get_rho(n)*wid_init*tratmid;
                    printout("cell %d, compton optical depth %g, grey optical depth %g\n",n,compton_optical_depth,grey_optical_deptha);
                    grey_optical_depth = get_kappagrey(n)*get_rho(n)*(rmax*tratmid-radial_pos);
                    printout("pos %g, distance %g, tau_dist %g\n",radial_pos,rmax*tratmid-radial_pos,grey_optical_depth);
                    modelgrid[n].grey_depth = grey_optical_depth;
                    
  //                   grey_optical_depth = SIGMA_T*nne*wid_init*tratmid;
                    if (grey_optical_depth > cell_is_optically_thick && m < n_grey_timesteps)
                    {
                      printout("cell %d is treated in grey approximation (tau %g)\n",n,grey_optical_depth);
                      modelgrid[n].thick = 1;
                    }
                    else if (grey_optical_depth > cell_is_optically_thick_vpkt )
                    {
                        modelgrid[n].thick = 2;
                    }
                    else modelgrid[n].thick = 0;
                  }
                #endif
                  
                
              }
              
              
    
              
              /// Temperatures are are now written to file outside the parallel block!
              
              /// Reset the volume estimators to zero for the current timestep
              /// Moved next lines to zero_estimators routine (SS)
              /*
              nuJ[n] = 0.;
              J[n] = 0.;
              cell[n].heating_ff = 0.;
              cell[n].heating_bf = 0.;
              */
              
              
              //set_kappagrey(n,SIGMA_T);
              
              /// Cooling rates depend only on cell properties, precalculate total cooling
              /// and ion contributions inside update grid and communicate between MPI tasks
              calculate_kpkt_rates(n);
            }
        
            /// MK Begin
            else
            {
	      /// For opacity_case != 4 the opacity treatment is grey. Enforce
              /// optically thick treatment in this case (should be equivalent to grey)
	      modelgrid[n].thick = 1;

              /// Need the total number density of bound and free electrons for Compton scatterin
              update_abundances(n, time_step[m].mid);
	      nne_tot = 0.;
	      for (element = 0; element < nelements; element++)
	      {
		abundance = get_abundance(n,element);
		/// calculate number density of the current element (abundances are given by mass)
		nnelement = abundance / elements[element].mass * get_rho(n);
		nne_tot += nnelement*get_element(element); 
	      }
              set_nnetot(n,nne_tot);

	      if (opacity_case == 3)
              {
                //printout("update_grid: opacity_case 3 ... updating cell[n].kappa_grey"); //MK
                if (get_rho(n) > rho_crit)
                {
                  set_kappagrey(n, opcase3_normal * (0.9 * get_ffe(n) + 0.1) * rho_crit/get_rho(n));
                }
                else
                {
                  set_kappagrey(n, opcase3_normal * (0.9 * get_ffe(n) + 0.1));
                }
	      }
            }
            /// MK End
        
            if (do_rlc_est == 2)
            {
              if (get_nne(n) > 0)
              {
                cell_len_scale = 0.1 / get_nne(n) / SIGMA_T;
                if (cell_len_scale < mps[tid])
                {
                  mps[tid] = cell_len_scale;
                }
                cell_len_scale = 0.1 / get_rho(n) / GREY_OP;
                if (cell_len_scale <  mps[tid])
                {
                  mps[tid] = cell_len_scale;
                }
              }
            }
            
            ///Non-OpenMP output of estimator files
            #ifndef _OPENMP
              //fprintf(estimators_file,"%d %g %g %g %g %d ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),modelgrid[n].thick);
              //fprintf(estimators_file,"%d %g %g %g %g %g ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),grey_optical_depth);
              fprintf(estimators_file,"%d %g %g %g %g %g ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),modelgrid[n].grey_depth);
              //fprintf(estimators_file,"%d %g %g %g %g %g %g %g ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),grey_optical_depth,grey_optical_deptha,compton_optical_depth);
              for (element = 0; element < nelements; element++)
              {
                nions = get_nions(element);
                for (ion = 0; ion < nions; ion++)
                {
                  fprintf(estimators_file,"%g ",ionstagepop(n,element,ion));
                }
              }
              
              #ifndef FORCE_LTE
                for (element = 0; element < nelements; element++)
                {
                  nions = get_nions(element);
                  for (ion = 0; ion < nions; ion++)
                  {
                    fprintf(estimators_file,"%g ",corrphotoionrenorm[n*nelements*maxion+element*maxion+ion]);
                  }
                }
                for (element = 0; element < nelements; element++)
                {
                  nions = get_nions(element);
                  for (ion = 0; ion < nions; ion++)
                  {
                    fprintf(estimators_file,"%g ",gammaestimator[n*nelements*maxion+element*maxion+ion]);
                  }
                }
                fprintf(estimators_file,"%g %g %g %g %g %g %g %g ",heatingrates[tid].ff,heatingrates[tid].bf,heatingrates[tid].collisional, heatingrates[tid].gamma,coolingrates[tid].ff,coolingrates[tid].fb,coolingrates[tid].collisional,coolingrates[tid].adiabatic);
              #endif
	      fprintf(estimators_file,"%g\n",energy_deposition[n]);
            #endif
          }
          else 
          { 
            ///Treat modelgrid cells which are not represented in the simulation grid
            ///Set grid properties to zero
            set_TR(n,0.);
            set_TJ(n,0.);
            set_Te(n,0.);
            set_W(n,0.);
            
            #ifndef FORCE_LTE
              for (element = 0; element < nelements; element++)
              {
                nions = get_nions(element);
                for (ion = 0; ion < nions; ion++)
                {
                  corrphotoionrenorm[n*nelements*maxion+element*maxion+ion] = 0.;
                  gammaestimator[n*nelements*maxion+element*maxion+ion] = 0.;
                }
              }
            #endif
  
            ///Non-OpenMP output of estimator files
            ///Indicates that this cell is not represented in the simulation grid by printing 0
            #ifndef _OPENMP
              //fprintf(estimators_file,"%d %g %g %g %g %g ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),compton_optical_depth);
              fprintf(estimators_file,"%d %g %g %g %g %g ",n,0.,0.,0.,0.,0.);
              //fprintf(estimators_file,"%d %g %g %g %g %g %g %g ",n,0.,0.,0.,0.,0.,0.,0.);
              for (element = 0; element < nelements; element++)
              {
                nions = get_nions(element);
                for (ion = 0; ion < nions; ion++)
                {
                  fprintf(estimators_file,"%g ",0.);
                }
              }
              
              #ifndef FORCE_LTE
                for (element = 0; element < nelements; element++)
                {
                  nions = get_nions(element);
                  for (ion = 0; ion < nions; ion++)
                  {
                    fprintf(estimators_file,"%g ",0.);
                  }
                }
                for (element = 0; element < nelements; element++)
                {
                  nions = get_nions(element);
                  for (ion = 0; ion < nions; ion++)
                  {
                    fprintf(estimators_file,"%g ",0.);
                  }
                }
                fprintf(estimators_file,"%g %g %g %g %g %g %g %g ",0.,0.,0.,0.,0.,0.,0.,0.);
              #endif
	      fprintf(estimators_file,"%g\n",0.);
            #endif
          }
          
        } 
        else
        {
          /// else, only reset gammaestimator to zero. This allows us to do a global MPI
          /// communication after update_grid to synchronize gammaestimator
          /// and write a contiguous restart file with grid properties
          for (element = 0; element < nelements; element++)
          {
            nions = get_nions(element);
            for (ion = 0; ion < nions-1; ion++)
            {
              gammaestimator[n*nelements*maxion+element*maxion+ion] = 0.;
            }
          }
        }
      } ///end parallel for loop over all modelgrid cells
      /// Now atfer all the relevant taks of update_grid have been finished activate
      /// the use of the cellhistory for all OpenMP tasks, in what follows (update_packets)
      use_cellhist = 1;
    #ifdef _OPENMP
    } /// end OpenMP parallel section
  #endif
  
  
  ///Write estimators for OpenMP version
  ///compton_optical_depth and heating/cooling rates which don't live on the grid
  ///are not available so far
  #ifdef _OPENMP
    for (n = nstart; n < nstart+nblock; n++)
    {
      if (modelgrid[n].associated_cells > 0)
      {
        //fprintf(estimators_file,"%d %g %g %g %g %d ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),modelgrid[n].thick);
        fprintf(estimators_file,"%d %g %g %g %g %g ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),modelgrid[n].grey_depth);
        //fprintf(estimators_file,"%d %g %g %g %g %g ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),grey_optical_depth);
        for (element = 0; element < nelements; element++)
        {
          nions = get_nions(element);
          for (ion = 0; ion < nions; ion++)
          {
            fprintf(estimators_file,"%g ",ionstagepop(n,element,ion));
          }
        }
        
        #ifndef FORCE_LTE
          for (element = 0; element < nelements; element++)
          {
            nions = get_nions(element);
            for (ion = 0; ion < nions; ion++)
            {
              fprintf(estimators_file,"%g ",corrphotoionrenorm[n*nelements*maxion+element*maxion+ion]);
            }
          }
          for (element = 0; element < nelements; element++)
          {
            nions = get_nions(element);
            for (ion = 0; ion < nions; ion++)
            {
              fprintf(estimators_file,"%g ",gammaestimator[n*nelements*maxion+element*maxion+ion]);
            }
          }
          //fprintf(estimators_file,"%g %g %g %g %g %g %g %g ",heatingrates[tid].ff,heatingrates[tid].bf,heatingrates[tid].collisional, heatingrates[tid].gamma,coolingrates[tid].ff,coolingrates[tid].fb,coolingrates[tid].collisional,coolingrates[tid].adiabatic);
          fprintf(estimators_file,"%g %g %g %g %g %g %g %g ",0.,0.,0., 0.,0.,0.,0.,0.);
        #endif
	fprintf(estimators_file,"%g \n",energy_deposition[n]);
      }
      else
      {
        ///Write zeros for cells which are non-represented in the simulation grid
        //fprintf(estimators_file,"%d %g %g %g %g %g ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),compton_optical_depth);
        fprintf(estimators_file,"%d %g %g %g %g %g ",n,0.,0.,0.,0.,0.);
        //fprintf(estimators_file,"%d %g %g %g %g %g %g %g ",n,0.,0.,0.,0.,0.,0.,0.);
        for (element = 0; element < nelements; element++)
        {
          nions = get_nions(element);
          for (ion = 0; ion < nions; ion++)
          {
            fprintf(estimators_file,"%g ",0.);
          }
        }
        
        #ifndef FORCE_LTE
          for (element = 0; element < nelements; element++)
          {
            nions = get_nions(element);
            for (ion = 0; ion < nions; ion++)
            {
              fprintf(estimators_file,"%g ",0.);
            }
          }
          for (element = 0; element < nelements; element++)
          {
            nions = get_nions(element);
            for (ion = 0; ion < nions; ion++)
            {
              fprintf(estimators_file,"%g ",0.);
            }
          }
          fprintf(estimators_file,"%g %g %g %g %g %g %g %g ",0.,0.,0.,0.,0.,0.,0.,0.);
        #endif
	fprintf(estimators_file,"%g \n",0.);
      }
    }
  #endif
    
    
  /// Move that check to outside the update grid routine. These numbers
  /// must be the same as they were during grid_init, but on the
  /// grid as a whole. Therefore densities and kappa must be propagated.
  /*
  //fclose(thermal_file);
  check1=check2=0.;  
  int mgi;
  for (n=0; n < ngrid; n++)
  {
    mgi=cell[n].modelgridindex;
    check1 = check1 + (get_kappagrey(mgi)  * get_rho(mgi));
    check2 = check2 + get_rho(mgi);
  }
  printout("Grey normalisation check: %g\n", check1/check2);
  */

  
  /// Assign the minimum of thread private mps to the global variable max_path_step
  max_path_step = mps[0];
  for (i=1; i < nthreads; i++)
  {
    if (mps[i] < max_path_step)
    {
      max_path_step = mps[i];
    }
  }

  
  if (do_rlc_est == 2)
  {
    if (max_path_step < (wid_init * trat / 10.))
    {
      max_path_step = wid_init / 10. * trat;
    }
  }
  //printout("max_path_step %g\n", max_path_step);
  
  //#ifndef FORCE_LTE
  //  fclose(heating_file);
  //#endif
  
  //printout("[debug] update_grid: update for timestep %d finished\n",m);
  //printf("time %ld\n",time(NULL));
  printout("[debug] update_grid: process %d finished update_grid at %d\n",my_rank,time(NULL));
  return 0;
}



///****************************************************************************
/// subroutine to identify the cell index from a position and a time. */
int get_cell(x, y, z, t)
    double x, y, z, t;
{
  double trat;
  int n;
  int nx, ny, nz;

  /* Original version of this was general but very slow. Modifying to be
  faster but only work for regular grid. */

  trat = t / tmin;
  nx = (x - (cell[0].pos_init[0] * trat))/(wid_init * trat);
  ny = (y - (cell[0].pos_init[1] * trat))/(wid_init * trat);
  nz = (z - (cell[0].pos_init[2] * trat))/(wid_init * trat);

  n = nx + (nxgrid * ny) + (nxgrid * nygrid * nz);

  /* do a check */

  if (x < cell[n].pos_init[0] * trat)
  {
    printout("Problem with get_cell (1).\n");
    exit(0);
  }
  if (x > (cell[n].pos_init[0]+wid_init) * trat)
  {
    printout("Problem with get_cell (2).\n");
    exit(0);
  }
  if (y < cell[n].pos_init[1] * trat)
  {
    printout("Problem with get_cell (3).\n");
    exit(0);
  }
  if (y > (cell[n].pos_init[1]+wid_init) * trat)
  {
    printout("Problem with get_cell (4).\n");
    exit(0);
  }
  if (z < cell[n].pos_init[2] * trat)
  {
    printout("Problem with get_cell (5).\n");
    exit(0);
  }
  if (z > (cell[n].pos_init[2]+wid_init) * trat)
  {
    printout("Problem with get_cell (6).\n");
    exit(0);
  }
  return(n);



  /* OLD
  trat = t / tmin;

  for (n = 0; n < ngrid; n++)
  {
  if (
  (x > cell[n].pos_init[0] * trat) &&
  (x < (cell[n].pos_init[0] + wid_init) *trat) &&
  (y > cell[n].pos_init[1] * trat) &&
  (y < (cell[n].pos_init[1] + wid_init) *trat) &&
  (z > cell[n].pos_init[2] * trat) &&
  (z < (cell[n].pos_init[2] + wid_init) *trat))
  {
  return(n);
}
}
  END OLD */



  printout("Failed to find cell (get_cell). \n");
  printout("x %g, y %g, z %g, t %g\n", x, y, z, t);
  printout("xend %g yend %g zend %g\n", cell[ngrid-1].pos_init[0] * trat, cell[ngrid-1].pos_init[1] * trat,cell[ngrid-1].pos_init[2] * trat);
  printout("xend %g yend %g zend %g\n", cell[0].pos_init[0] * trat, cell[0].pos_init[1] * trat,cell[0].pos_init[2] * trat);
  printout("xend %g yend %g zend %g\n", (cell[0].pos_init[0]+wid_init) * trat, (cell[0].pos_init[1]+wid_init) * trat,(cell[0].pos_init[2]+wid_init) * trat);
  printout("xwid %g ywid %g zwid %g\n", (wid_init) * trat, (wid_init) * trat,(wid_init) * trat);
    
  exit(0);
}



///****************************************************************************
/*double get_abundance(int cellnumber, int element)
{
  return modelgrid[cellnumber].composition[element].abundance;
}*/


///****************************************************************************
double get_abundance(int modelgridindex, int element)
{
  return modelgrid[modelgridindex].composition[element].abundance;
}



///***************************************************************************/
void update_abundances(int modelgridindex, double t_current)
/// Updates the mass fractions of elements associated with the decay sequence
/// (56)Ni -> (56)Co -> (56)Fe at the onset of each timestep
/// Parameters: - modelgridindex: the grid cell for which to update the abundances
///             - t_current: current time (here mid of current timestep)
{
  //double get_abundance(int modelgridindex, int element);
  int get_elementindex(int anumber);
  int get_element(int element);
  double nifrac,cofrac,fefrac,mnfrac,crfrac,vfrac,tifrac;
  double ni_in,co_in,fe_in,fe52_in,cr48_in;
  double lambdani,lambdaco;
  double lambdacr,lambdav;
  double lambdafe,lambdamn;
  double deltat;
  int element;
  
  t_current -= t_model;
  lambdani = 1./TNICKEL;
  lambdaco = 1./TCOBALT;
  lambdafe = 1./T52FE;
  lambdamn = 1./T52MN;
  lambdacr = 1./T48CR;
  lambdav = 1./T48V;
  
  if (homogeneous_abundances == 1)
  {
    ni_in = elements[get_elementindex(28)].abundance;
    co_in = elements[get_elementindex(27)].abundance;
    fe_in = elements[get_elementindex(26)].abundance;
    //fe_in = cell[modelgridindex].f_fe_init;
    for (element = nelements-1; element >= 0; element--)
    {
      if (get_element(element) == 28)
      {
        nifrac = ni_in * exp(-lambdani*t_current) + modelgrid[modelgridindex].fnistable;
        modelgrid[modelgridindex].composition[element].abundance = nifrac;
      }
      if (get_element(element) == 27)
      {
        cofrac = co_in*exp(-lambdaco*t_current) + lambdani*ni_in/(lambdani-lambdaco)*(exp(-lambdaco*t_current)-exp(-lambdani*t_current)) + modelgrid[modelgridindex].fcostable;
        modelgrid[modelgridindex].composition[element].abundance = cofrac;
      }
      if (get_element(element) == 26) 
      {
        fefrac = fe_in + (co_in*lambdani - co_in*lambdaco + ni_in*lambdani - ni_in*lambdaco - co_in*lambdani*exp(-lambdaco*t_current) + co_in*lambdaco*exp(-lambdaco*t_current) - ni_in*lambdani*exp(-lambdaco*t_current) + ni_in*lambdaco*exp(-lambdani*t_current)) / (lambdani-lambdaco);
        modelgrid[modelgridindex].composition[element].abundance = fefrac;
      }
    }
  }
  else
  {
    ni_in = modelgrid[modelgridindex].fni;
    co_in = modelgrid[modelgridindex].fco;
    fe52_in = modelgrid[modelgridindex].f52fe;
    cr48_in = modelgrid[modelgridindex].f48cr;
    //printout("model cell %d, has ni_in %g, co_in %g, fe_in %g\n",modelgridindex,ni_in,co_in,fe_in);
    for (element = nelements-1; element >= 0; element--)
    {
      if (get_element(element) == 28)
      {
        nifrac = ni_in * exp(-lambdani*t_current) + modelgrid[modelgridindex].fnistable;
        modelgrid[modelgridindex].composition[element].abundance = nifrac;
        //cell[modelgridindex].composition[element].abundance = nifrac;
      }
      if (get_element(element) == 27)
      {
        cofrac = co_in*exp(-lambdaco*t_current) + lambdani*ni_in/(lambdani-lambdaco)*(exp(-lambdaco*t_current)-exp(-lambdani*t_current)) + modelgrid[modelgridindex].fcostable;
        modelgrid[modelgridindex].composition[element].abundance = cofrac;
      }
      if (get_element(element) == 26) 
      {
        fefrac = ((co_in*lambdani - co_in*lambdaco + ni_in*lambdani - ni_in*lambdaco - co_in*lambdani*exp(-lambdaco*t_current) + co_in*lambdaco*exp(-lambdaco*t_current) - ni_in*lambdani*exp(-lambdaco*t_current) + ni_in*lambdaco*exp(-lambdani*t_current)) / (lambdani-lambdaco)) + modelgrid[modelgridindex].ffestable + (fe52_in* exp(-lambdafe*t_current));
        modelgrid[modelgridindex].composition[element].abundance = fefrac;
      }
      if (get_element(element) == 25) 
      {
        mnfrac = lambdafe*fe52_in/(lambdafe-lambdamn)*(exp(-lambdamn*t_current)-exp(-lambdafe*t_current)) + modelgrid[modelgridindex].fmnstable;
        modelgrid[modelgridindex].composition[element].abundance = mnfrac;
      }
      if (get_element(element) == 24) 
      {
        crfrac = ((fe52_in*lambdafe - fe52_in*lambdamn - fe52_in*lambdafe*exp(-lambdamn*t_current) + fe52_in*lambdamn*exp(-lambdafe*t_current)) / (lambdafe-lambdamn)) + modelgrid[modelgridindex].fcrstable + (cr48_in * exp(-lambdacr*t_current));
	modelgrid[modelgridindex].composition[element].abundance = crfrac;
      }
      if (get_element(element) == 23) 
      {
        vfrac = lambdacr*cr48_in/(lambdacr-lambdav)*(exp(-lambdav*t_current)-exp(-lambdacr*t_current)) + modelgrid[modelgridindex].fvstable;
        modelgrid[modelgridindex].composition[element].abundance = vfrac;
      }
      if (get_element(element) == 22) 
      {
        tifrac = ((cr48_in*lambdacr - cr48_in*lambdav - cr48_in*lambdacr*exp(-lambdav*t_current) + cr48_in*lambdav*exp(-lambdacr*t_current)) / (lambdacr-lambdav)) + modelgrid[modelgridindex].ftistable;
        modelgrid[modelgridindex].composition[element].abundance = tifrac;
      }
 


    }
    //printout("model cell %d, has ni_in %g, co_in %g, fe_in %g, abund %g, %g, %g, stable %g,%g\n",modelgridindex,ni_in,co_in,fe_in,nifrac,cofrac,fefrac,modelgrid[modelgridindex].fnistable,modelgrid[modelgridindex].fcostable);
  }
  //printout("nifrac_old %g, cofrac_old %g, fefrac_old %g\n",nifrac_old,cofrac_old,fefrac_old);
  //printout("nifrac_new %g, cofrac_new %g, fefrac_new %g\n",nifrac_new,cofrac_new,fefrac_new);
}



///****************************************************************************
double calculate_populations(int modelgridindex, int first_nonempty_cell)
/// Determines the electron number density for a given cell using one of 
/// libgsl's root_solvers and calculates the depending level populations.
{
  double ionfract(int element, int ion, int modelgridindex, double nne);
  double phi(int element, int ion, int modelgridindex);
  double stat_weight(int element, int ion, int level);
  double nne_solution_f(double x, void *paras);
  int get_ionstage(int element, int ion);
  int get_element(int element);
  double get_abundance(int modelgridindex, int element);

  double nne,nne_lo,nne_hi,nne_check,nne_tot;
  double nnelement,nnion;
  int element,ion;
  int nions;
  int iter;
  double nntot;
  
  /// Initialise the gsl solver
  const gsl_root_fsolver_type *solvertype;
  solvertype = gsl_root_fsolver_brent;
  gsl_root_fsolver *solver;
  solver = gsl_root_fsolver_alloc(solvertype);
  double fractional_accuracy = 1e-3;
  int maxit = 100;
  int status;
  int uppermost_ion,only_neutral_ions,i,nelements_in_cell;
  double factor,abundance;
  
  /// and the solution function
  gsl_function f;
  nne_solution_paras paras;
  paras.cellnumber = modelgridindex;
  f.function = &nne_solution_f;
  f.params = &paras;
  
  neutral_flag = 0;
  
  /// Get temperatures
  double T_R = get_TR(modelgridindex);
  double T_e = get_Te(modelgridindex);
  double W = get_W(modelgridindex);
  
  //int mgi = cell[cellnumber].modelgridindex;
  
  nne_lo = 0.;  //MINPOP;
  nne_hi = get_rho(modelgridindex)/MH;
  
  /// The following selection of uppermost_ion is (so far) NOT thread safe!!!!!!!!!!!!!!!!!!!!!!!
  only_neutral_ions = 0;
  nelements_in_cell = 0;
  for (element = 0; element < nelements; element++)
  {
    nions = get_nions(element);
    //elements[element].uppermost_ion = nions-1;
    elements_uppermost_ion[tid][element] = nions-1;
    abundance = get_abundance(modelgridindex,element);
    if (abundance > 0)
    {
      #ifdef FORCE_LTE
        uppermost_ion = get_nions(element)-1;
      #else
        if (initial_iteration == 1 || modelgrid[modelgridindex].thick == 1)
        {
          uppermost_ion = get_nions(element)-1;
        }
        else
        {
          for (ion = 0; ion < nions-1; ion++)
          {
            //printout("element %d, ion %d, photoionest %g\n",element,ion,photoionestimator[modelgridindex*nelements*maxion+element*maxion+ion]);
            //if (photoionestimator[modelgridindex*nelements*maxion+element*maxion+ion] == 0) break;
            if (gammaestimator[modelgridindex*nelements*maxion+element*maxion+ion] == 0) break;
          }
          uppermost_ion = ion;
        }
      #endif
      //printout("cell %d, element %d, uppermost_ion by Gamma is %d\n",cellnumber,element,uppermost_ion);
      
      factor = 1.;
      for (i = 0; i < uppermost_ion; i++)
      {
        factor *= nne_hi * phi(element,i,modelgridindex);
        //printout("element %d, ion %d, factor %g\n",element,i,factor);
        if (!isfinite(factor)) 
        {
          printout("[info] calculate_populations: uppermost_ion limited by phi factors for element %d, ion %d in cell %d\n",element,i,modelgridindex);
          break;
        }
      }
      uppermost_ion = i;
      //printout("cell %d, element %d, final uppermost_ion is %d, factor %g\n",modelgridindex,element,uppermost_ion,factor);
      //elements[element].uppermost_ion = uppermost_ion;
      elements_uppermost_ion[tid][element] = uppermost_ion;
      if (uppermost_ion == 0) only_neutral_ions += 1;
      nelements_in_cell += 1;
    }
  }
  
  
  nne_check = 0.;
  nne_tot = 0.;   /// total number of electrons in grid cell which are possible 
                  /// targets for compton scattering of gamma rays
  nntot = 0.;
  if (only_neutral_ions == nelements_in_cell)
  {
    /// Special case of only neutral ions, set nne to some finite value that
    /// packets are not lost in kpkts
    /// Introduce a flag variable which is sent to the T_e solver so that
    /// we get this info only once when T_e is converged and not for each
    /// iteration step.
    neutral_flag = 1;
    //printout("[warning] calculate_poulations: only neutral ions in cell %d modelgridindex\n",modelgridindex);
    //exit(0);
    /// Now calculate the ground level populations in nebular approximation and store them to the grid
    
    for (element = 0; element < nelements; element++)
    {
      abundance = get_abundance(modelgridindex,element);
      nions = get_nions(element);
      /// calculate number density of the current element (abundances are given by mass)
      nnelement = abundance / elements[element].mass * get_rho(modelgridindex);
      nne_tot += nnelement*get_element(element); 
      
      /// Assign the species population to the neutral ion and set higher ions to MINPOP
      for (ion = 0; ion < nions; ion++)
      {
        if (ion == 0) nnion = nnelement;
        else 
        {
          if (abundance > 0.) nnion = MINPOP;
          else nnion = 0.;
        }
        nntot += nnion;
        nne += nnion * (get_ionstage(element,ion)-1);
        modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = nnion * stat_weight(element,ion,0) / modelgrid[modelgridindex].composition[element].partfunct[ion];
        if (!isfinite(modelgrid[modelgridindex].composition[element].groundlevelpop[ion])) 
          printout("[warning] calculate_populations: groundlevelpop infinite in connection with MINPOP\n");
      }
    }
    nntot += nne;
    if (nne < MINPOP) nne = MINPOP;
    set_nne(modelgridindex,nne);
  }
  else
  {
    /// Apply solver to get nne
    /// Search solution for nne in [nne_lo,nne_hi]
    //printout("nne@x_lo %g\n", nne_solution_f(nne_lo,f.params));
    //printout("nne@x_hi %g\n", nne_solution_f(nne_hi,f.params));
    //printout("n, x_lo, x_hi, T_R, T_e, W, rho %d, %g, %g, %g, %g, %g, %g\n",modelgridindex,x_lo,x_hi,T_R,T_e,W,cell[modelgridindex].rho);
    if (nne_solution_f(nne_lo,f.params)*nne_solution_f(nne_hi,f.params) > 0)
    {
      printout("n, nne_lo, nne_hi, T_R, T_e, W, rho %d, %g, %g, %g, %g, %g, %g\n",modelgridindex,nne_lo,nne_hi,T_R,T_e,W,get_rho(modelgridindex));
      printout("nne@x_lo %g\n", nne_solution_f(nne_lo,f.params));
      printout("nne@x_hi %g\n", nne_solution_f(nne_hi,f.params));
      #ifndef FORCE_LTE
        for (element=0; element < nelements; element++)
        {
          //printout("cell %d, element %d, uppermost_ion is %d\n",modelgridindex,element,elements[element].uppermost_ion);
          printout("cell %d, element %d, uppermost_ion is %d\n",modelgridindex,element,elements_uppermost_ion[tid][element]);
          //for (ion=0; ion <= elements[element].uppermost_ion; ion++)
          for (ion=0; ion <= elements_uppermost_ion[tid][element]; ion++)
          {
            //printout("element %d, ion %d, photoionest %g\n",element,ion,photoionestimator[modelgridindex*nelements*maxion+element*maxion+ion]);
            printout("element %d, ion %d, photoionest %g\n",element,ion,gammaestimator[modelgridindex*nelements*maxion+element*maxion+ion]);
          }
        }
      #endif
    }
    gsl_root_fsolver_set(solver, &f, nne_lo, nne_hi);
    iter = 0;
    do
    {
      iter++;
      status = gsl_root_fsolver_iterate(solver);
      nne = gsl_root_fsolver_root(solver);
      nne_lo = gsl_root_fsolver_x_lower(solver);
      nne_hi = gsl_root_fsolver_x_upper(solver);
      status = gsl_root_test_interval(nne_lo,nne_hi,0,fractional_accuracy);
      //if (first_nonempty_cell == -1000) printout("[debug] update_grid:   %d [%g, %g] %g %g\n",iter,x_lo,x_hi,x_0,x_hi-x_lo);
    }
    while (status == GSL_CONTINUE && iter < maxit);
    gsl_root_fsolver_free(solver);
    if (nne < MINPOP) nne = MINPOP;
    set_nne(modelgridindex,nne);
    //cell[modelgridindex].nne = nne;
    if (status == GSL_CONTINUE) printout("[warning] calculate_populations: nne did not converge within %d iterations\n",maxit);
    //printout("[debug] update_grid:   status = %s\n",gsl_strerror(status));
    //printout("[debug] update_grid:   converged nne %g\n",cell[modelgridindex].nne);
    
    /// Now calculate the ground level populations in nebular approximation and store them to the grid
    nne_check = 0.;
    nne_tot = 0.;   /// total number of electrons in grid cell which are possible 
                    /// targets for compton scattering of gamma rays
    
    nntot = nne;
    for (element = 0; element < nelements; element++)
    {
      abundance = get_abundance(modelgridindex,element);
      nions = get_nions(element);
      /// calculate number density of the current element (abundances are given by mass)
      nnelement = abundance / elements[element].mass * get_rho(modelgridindex);
      nne_tot += nnelement*get_element(element); 
      
      /// Use ionisationfractions to calculate the groundlevel populations
      for (ion = 0; ion < nions; ion++)
      {
        //if (ion <= elements[element].uppermost_ion) 
        if (ion <= elements_uppermost_ion[tid][element]) 
        {
          if (abundance > 0) 
          {
            nnion = nnelement * ionfract(element,ion,modelgridindex,nne);
            if (nnion < MINPOP) nnion = MINPOP;
          }
          else nnion = 0.;
        }
        else
          nnion = MINPOP;  /// uppermost_ion is only < nions-1 in cells with nonzero abundance of the given species
        nntot += nnion;
        nne_check += nnion * (get_ionstage(element,ion)-1);
        modelgrid[modelgridindex].composition[element].groundlevelpop[ion] = nnion * stat_weight(element,ion,0) / modelgrid[modelgridindex].composition[element].partfunct[ion];
        if (!isfinite(modelgrid[modelgridindex].composition[element].groundlevelpop[ion])) 
          printout("[warning] calculate_populations: groundlevelpop infinite in connection with MINPOP\n");
      }
    }
  }
  
  
  set_nnetot(modelgridindex, nne_tot);
  return nntot;
}


///****************************************************************************
void precalculate_partfuncts(int modelgridindex)
/// The partition functions depend only on T_R and W. This means they don't
/// change during any iteration on T_e. Therefore their precalculation was 
/// taken out of calculate_populations to save runtime.
{
  double calculate_partfunct(int element, int ion, int modelgridindex);
  //double calculate_ltepartfunct(int element, int ion, double T);
  int element,ion,nions;
  
  /// Precalculate partition functions for each ion in every cell
  /// this saves a factor 10 in calculation time of Saha-Boltzman populations
  for (element = 0; element < nelements; element++)
  {
    nions = get_nions(element);
    for (ion = 0; ion < nions; ion++)
    {
      //printout("precalc element %d, ion %d, mgi %d\n",element,ion,modelgridindex);
      //cell[cellnumber].composition[element].ltepartfunct[ion] = calculate_ltepartfunct(element,ion,T_R);
      modelgrid[modelgridindex].composition[element].partfunct[ion] = calculate_partfunct(element,ion,modelgridindex);
    }
  }
}







void get_radfield_params(double J, double nuJ, int modelgridindex, double *T_J, double *T_R, double *W)
{
  int i;
  double T,nubar;
  
  nubar = nuJ/J;
  if (!isfinite(nubar) || nubar == 0.)
  {
    /// Return old T_R
    printout("[warning] update_grid: T_R estimator infinite in cell %d, use value of last timestep\n",modelgridindex);
    *T_J = modelgrid[modelgridindex].TJ;
    *T_R = modelgrid[modelgridindex].TR;
    *W = modelgrid[modelgridindex].W;
  }
  else
  {
    *T_J = pow(PI/STEBO*J,1./4.);
    if (*T_J > MAXTEMP)
    {
      printout("[warning] update_grid: temperature estimator T_J=%g exceeds T_max=%g in cell %d. Set T_J = T_max!\n",*T_J,MAXTEMP,modelgridindex);
      *T_J = MAXTEMP;
    }
    if (*T_J < MINTEMP)
    {
      printout("[warning] update_grid: temperature estimator T_J=%g below T_min %g in cell %d. Set T_J = T_min!\n",*T_J,MINTEMP,modelgridindex);
      *T_J = MINTEMP;
    }
    
    *T_R = H*nubar/KB/3.832229494;
    if (*T_R > MAXTEMP)
    {
      printout("[warning] update_grid: temperature estimator T_R=%g exceeds T_max=%g in cell %d. Set T_R = T_max!\n",*T_R,MAXTEMP,modelgridindex);
      *T_R = MAXTEMP;
    }
    if (*T_R < MINTEMP)
    {
      printout("[warning] update_grid: temperature estimator T_R=%g below T_min %g in cell %d. Set T_R = T_min!\n",*T_R,MINTEMP,modelgridindex);
      *T_R = MINTEMP;
    }
    
    *W = PI*J/STEBO/pow(*T_R,4);
  }
}










///****************************************************************************
//double nuB_nu_integrand(double nu, void *paras)
//{
//  double T = ((gslintegration_paras *) paras)->T;
//  
//  return nu * TWOHOVERCLIGHTSQUARED * pow(nu,3) * 1/(exp(H*nu/KB/T)-1);
//}
//
/////****************************************************************************
//double B_nu_integrand(double nu, void *paras)
//{
//  double T = ((gslintegration_paras *) paras)->T;
//  
//  return TWOHOVERCLIGHTSQUARED * pow(nu,3) * 1/(exp(H*nu/KB/T)-1);
//}



#ifndef FORCE_LTE

/*
****************************************************************************
double get_ffcooling(int element, int ion, int cellnumber)
{
  double ionstagepop(int cellnumber, int element, int ion);
  int get_ionstage(int element, int ion);
  int ioncharge;
  double T_e,nne,nncurrention,C_ff;
  
  ion += 1;
  nncurrention = ionstagepop(cellnumber,element,ion);
  nne = cell[cellnumber].nne;
  T_e = cell[cellnumber].T_e;
      
  C_ff = 0.;
  ioncharge = get_ionstage(element,ion) - 1;
  if (ioncharge > 0)
  {
    C_ff += 1.426e-27 * sqrt(T_e) * pow(ioncharge,2) * nncurrention * nne;
  }
  
  return C_ff;
}

double get_adcool(int cellnumber, int nts)
{
  
  double t_current,nntot,T_e;
  double p,dV,V;
  
  t_current = time_step[nts].mid;
  nntot = cell[cellnumber].nn_tot;
  T_e = cell[cellnumber].T_e;
  p = nntot*KB*T_e;
  dV = 3*pow(wid_init/tmin,3)*pow(t_current,2);  
  V = pow(wid_init*t_current/tmin,3);
  
  return p*dV/V; 
}

double get_radrecomb(int element, int ion, int level, int cellnumber)
{
  double interpolate_spontrecombcoeff(int element, int ion, int level, double T);
  double get_groundlevelpop(int cellnumber, int element, int ion);
  double epsilon(int element, int ion, int level);
  double T_e,nne,nnion,E_threshold,alpha_sp;
  
  T_e = cell[cellnumber].T_e;
  nne = cell[cellnumber].nne;
  nnion = get_groundlevelpop(cellnumber,element,ion+1);
  E_threshold = epsilon(element,ion+1,0) - epsilon(element,ion,level);
  alpha_sp = interpolate_spontrecombcoeff(element,ion,level,T_e);
  ///// This is ion-wise Alpha_st is the sum of the level-wise alpha_st, so it should be added only once if mulitple levels are present
  //double Alpha_st = stimrecombestimator_save[cellnumber*nelements*maxion+element*maxion+ion];
  
  //return nnion*nne*(alpha_sp+Alpha_st)*E_threshold;
  return nnion*nne*alpha_sp*E_threshold;
}

*/

/*
double get_thermalratio(int element, int ion, int level, int cellnumber)
{
  double interpolate_spontrecombcoeff(int element, int ion, int level, double T);
  double interpolate_spontrecombcoeff_E(int element, int ion, int level, double T);
  double T_e,alpha_sp,alpha_sp_E;
  
  T_e = cell[cellnumber].T_e;
  alpha_sp = interpolate_spontrecombcoeff(element,ion,level,T_e);
  alpha_sp_E = interpolate_spontrecombcoeff_E(element,ion,level,T_e);
  
  return alpha_sp/alpha_sp_E;
}*/






/*



double get_Gamma(int cellnumber, int element, int ion)
///photoionization rate: paperII 3.5.2
/// n_1 - occupation number of ground state
{
  double get_corrphotoioncoeff(int element, int ion, int level, int cellnumber);
  double calculate_exclevelpop(int cellnumber, int element, int ion, int level);
  
  int i,level;
  int nlevels = get_nlevels(element,ion);
  double Gamma = 0.;
  
  for (i = 0; i < nlevels; i++)
  {
    Gamma += calculate_exclevelpop(cellnumber,element,ion,level)*get_corrphotoioncoeff(element,ion,level,cellnumber);
  }
  Gamma /= calculate_exclevelpop(cellnumber,element,ion,0);
  
  return Gamma;
}


double get_Gamma_phys(int cellnumber, int element, int ion)
///photoionization rate: paperII 3.5.2
/// n_1 - occupation number of ground state
{
  double get_corrphotoioncoeff_ana(int element, int ion, int level, int cellnumber);
  double get_corrphotoioncoeff(int element, int ion, int level, int cellnumber);
  double calculate_exclevelpop(int cellnumber, int element, int ion, int level);
  
  int i,level;
  int nlevels = get_nlevels(element,ion);
  double Gamma = 0.;
  
  //for (i = 0; i < nlevels; i++)
  for (i = 0; i < TAKE_N_BFCONTINUA; i++)
  {
    Gamma += get_corrphotoioncoeff_ana(element,ion,level,cellnumber);
  }
  
  return Gamma;
}

*/
#endif



void write_grid_restart_data(void)
{
  int n,element,nions,ion;
  FILE *gridsave_file;
  
  if ((gridsave_file = fopen("gridsave.dat", "w")) == NULL)
  {
    printout("Cannot open gridsave.dat.\n");
    exit(0);
  }

  for (n = 0; n < npts_model; n++)
  {
    if (modelgrid[n].associated_cells > 0)
    {
      fprintf(gridsave_file,"%d %g %g %g %g %d ",n,get_TR(n),get_Te(n),get_W(n),get_TJ(n),modelgrid[n].thick);
      
      #ifndef FORCE_LTE
        for (element = 0; element < nelements; element++)
        {
          nions = get_nions(element);
          for (ion = 0; ion < nions; ion++)
          {
            fprintf(gridsave_file,"%g ",corrphotoionrenorm[n*nelements*maxion+element*maxion+ion]);
          }
        }
        for (element = 0; element < nelements; element++)
        {
          nions = get_nions(element);
          for (ion = 0; ion < nions; ion++)
          {
            fprintf(gridsave_file,"%g ",gammaestimator[n*nelements*maxion+element*maxion+ion]);
          }
        }
      #endif
      
      fprintf(gridsave_file,"\n");
    }
    else
    {
      ///Write zeros for cells which are non-represented in the simulation grid
      fprintf(gridsave_file,"%d %g %g %g %g %d ",n,0.,0.,0.,0.,0);
      
      #ifndef FORCE_LTE
        for (element = 0; element < nelements; element++)
        {
          nions = get_nions(element);
          for (ion = 0; ion < nions; ion++)
          {
            fprintf(gridsave_file,"%g ",0.);
          }
        }
        for (element = 0; element < nelements; element++)
        {
          nions = get_nions(element);
          for (ion = 0; ion < nions; ion++)
          {
            fprintf(gridsave_file,"%g ",0.);
          }
        }
      #endif
      
      fprintf(gridsave_file,"\n");
    }
  }
  
  fclose(gridsave_file);
}
