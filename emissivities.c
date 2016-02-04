#include "sn3d.h"
#include <string.h>

int add_gam_line_emissivity(ray_ptr, nray, single_pos, single_t, lindex, dnuds)
     RAY *ray_ptr;
     int nray, lindex;
     double single_t;
     double *single_pos;
     double dnuds;
{
  
  double vel_vec[3];
  int get_velocity();
  double doppler();
  double emitt_energy;
  struct grid *grid_ptr;
  double fni(CELL *grid_ptr),f48cr(CELL *grid_ptr);
  double tfact;

  grid_ptr = &cell[ray_ptr->where];
  tfact = pow((tmin/single_t), 3);

  if (gam_line_list.type[lindex] == NI_GAM_LINE_ID)
  {
    emitt_energy = get_rhoinit(grid_ptr->modelgridindex) / MNI56 / 4. / PI
	* exp(-single_t/TNICKEL) / TNICKEL * 
	nickel_spec.probability[gam_line_list.index[lindex]]
        * nickel_spec.energy[gam_line_list.index[lindex]]
	*fni(grid_ptr)*tfact;
  }
  else if (gam_line_list.type[lindex] == CO_GAM_LINE_ID)
  {
      
    emitt_energy = get_rhoinit(grid_ptr->modelgridindex) / MNI56 / 4. / PI
      * (exp(-single_t/TNICKEL) - exp(-single_t/TCOBALT))
      / (TNICKEL - TCOBALT)
          *cobalt_spec.probability[gam_line_list.index[lindex]]
          *cobalt_spec.energy[gam_line_list.index[lindex]]
      *fni(grid_ptr)*tfact;

    if (gam_line_list.index[lindex] == 0)
    {
      emitt_energy += (compton_emiss[grid_ptr->modelgridindex][emiss_max - 1] * 1.e20 / 4. / PI);
    }
  }
  else if (gam_line_list.type[lindex] == CR48_GAM_LINE_ID)
  {
    
    emitt_energy = get_rhoinit(grid_ptr->modelgridindex) / MCR48 / 4. / PI
    * exp(-single_t/T48CR) / T48CR * 
    cr48_spec.probability[gam_line_list.index[lindex]]
    * cr48_spec.energy[gam_line_list.index[lindex]]
    *f48cr(grid_ptr)*tfact;

  }
  else if (gam_line_list.type[lindex] == V48_GAM_LINE_ID)
  {
      
      emitt_energy = get_rhoinit(grid_ptr->modelgridindex) / MCR48 / 4. / PI
	* (exp(-single_t/T48CR) - exp(-single_t/T48V))
	/ (T48CR - T48V)
        *v48_spec.probability[gam_line_list.index[lindex]]
        *v48_spec.energy[gam_line_list.index[lindex]]
	*f48cr(grid_ptr)*tfact;
  }
  else if (gam_line_list.type[lindex] == FAKE_GAM_LINE_ID)
  {
    emitt_energy = 0.0;
  }
  else
  {
    printout("unknown line??\n");
    exit(0);
  }
  

  /* I'm changing the next bit here (Jan 06) because I think what was
     here before (below) was wrong in the Doppler terms. 

     ray_ptr->e_cmf[nray]+= emitt_energy / fabs(dnuds);
     
     get_velocity(single_pos, vel_vec, single_t);
     
     ray_ptr->e_rf[nray] = ray_ptr->e_cmf[nray] / doppler(syn_dir, vel_vec);
 
  */ 
  if (emitt_energy != 0)
  {
    get_velocity(single_pos, vel_vec, single_t);
      
    ray_ptr->e_rf[nray]+= emitt_energy / fabs(dnuds)
	/doppler(syn_dir, vel_vec) / doppler(syn_dir, vel_vec) ;
  }

  
  return(0);
}

/*******************************************************/
int continuum_rt(ray_ptr, nray, ldist, single_pos, single_t, lindex)
     RAY *ray_ptr;
     int nray, lindex;
     double ldist, single_t;
     double *single_pos;
{
  /* This is called when a ray is about to be moved a distance ldist. */
  /* It should account for the changes in the ray intensity due to
     continuum processes along the path. */
  PKT dummy;
  double kap_compton, kap_photo_electric, kap_tot, tau_cont, kap_pair_prod;
  double sig_comp(), sig_photo_electric();
  double sig_pair_prod(), doppler(); 
  int get_velocity();
  double vel_vec[3], dop_fac;

  /* Make a dummy packet that carries the ray properties. */

  dummy.pos[0] = single_pos[0];
  dummy.pos[1] = single_pos[1];
  dummy.pos[2] = single_pos[2];

  dummy.dir[0] = syn_dir[0];
  dummy.dir[1] = syn_dir[1];
  dummy.dir[2] = syn_dir[2];
  
  dummy.where = ray_ptr->where;
  dummy.nu_cmf = ray_ptr->nu_cmf[nray];

  kap_compton = sig_comp(&dummy,single_t);
  kap_photo_electric = sig_photo_electric(&dummy, single_t);
  kap_pair_prod = sig_pair_prod(&dummy, single_t);
  kap_tot = kap_compton + kap_photo_electric + kap_pair_prod;
   
  /* For now no emissivity - only destruction. So very simple. */

  tau_cont = kap_tot * ldist;

  ray_ptr->e_rf[nray] = ray_ptr->e_rf[nray] * exp(-1. * tau_cont);

  /* Now adding the emissivity term. */
 
  if (lindex != RED_OF_LIST)
  {
    get_velocity(single_pos, vel_vec, single_t);
    dop_fac = 1./doppler(syn_dir, vel_vec);
      
    if (tau_cont > 1.e-6)
    {
      ray_ptr->e_rf[nray] += (dop_fac * dop_fac * compton_emiss[cell[dummy.where].modelgridindex][lindex - emiss_offset] *
				  (1. - exp(-1. * tau_cont)) / kap_tot);
    }
    else
    {
      ray_ptr->e_rf[nray] += (dop_fac * dop_fac * compton_emiss[cell[dummy.where].modelgridindex][lindex - emiss_offset] *
	  			  ldist);
    }
  }
  
  /* This MUST be followed by a call to move_one_ray() in source
     since e_cmf is NOT reset here. */

  return(0);
}

/*******************************************************/
int compton_emiss_cont(pkt_ptr, dist, t_current)
     PKT *pkt_ptr;
     double dist, t_current;
{
  /* Subroutine to add contriubtion to the MC estimator for the
compton emissivity. Called with a packet that is about to travel a 
distance dist in the lab frame. Time at start of distance is t_current.*/

  double mu_cmf;
  double vel_vec[3], cmf_dir[3], cmf_syn_dir[3];
  double dot();
  int get_velocity(), angle_ab();
  double doppler();
  double vec_len();
  double f, dsigma_domega_cmf;
  double freq_out, emiss_cont, dop_fac;
  int lindex, get_nul();
 
  /* First we need to know the scattering angle needed from the
packet's direction of motion to the desired observer. Call this angle 
mu_cmf (it's a cosine). To get it convert both the direction of
motion and the local velocity vectors to the cmf.*/

  get_velocity(pkt_ptr->pos, vel_vec, t_current);
  angle_ab(pkt_ptr->dir, vel_vec, cmf_dir);
  angle_ab(syn_dir, vel_vec, cmf_syn_dir);


  //  printout("pos %g %g %g\n", pkt_ptr->pos[0],pkt_ptr->pos[1], pkt_ptr->pos[2]);
  //  printout("dir %g %g %g\n", pkt_ptr->dir[0],pkt_ptr->dir[1], pkt_ptr->dir[2]);
  //  printout("vel %g %g %g\n", vel_vec[0], vel_vec[1], vel_vec[2]);
  //  printout("cmf_dir %g %g %g\n", cmf_dir[0], cmf_dir[1], cmf_dir[2]);
  //  printout("syn_dir %g %g %g\n", syn_dir[0], syn_dir[1], syn_dir[2]);
  //  printout("cmf_syn_dir %g %g %g\n", cmf_syn_dir[0], cmf_syn_dir[1], cmf_syn_dir[2]);
  



  mu_cmf = dot(cmf_dir,cmf_syn_dir);

  if (mu_cmf > 1 || mu_cmf < -1)
  {
    printout("problem with Compton emissivity. Abort.\n");
    exit(0);
  }

  /* Now get the factor by which the frequency will change, f, for going 
     in this direction. f = old energy / new eneregy - always should be > 1*/

  f = 1 + (H * pkt_ptr->nu_cmf / ME / CLIGHT / CLIGHT * (1. - mu_cmf));

  //  printout("compton reducion factor %g freq %g mu %g\n", f, H*pkt_ptr->nu_cmf/MEV, mu_cmf );

  /* Now work out in which frequency bin this'll happen. The scattered
     light will have frequency (nu_cmf / f) in the cmf frame. And it
     travels in direction syn_dir in the rf. */

  freq_out = pkt_ptr->nu_cmf /f; /// doppler(syn_dir, vel_vec); 
  // do we want ?/ doppler(syn_dir, vel_vec)

  lindex = get_nul(freq_out); // This is the index of the next line to 
                              // the red. The emissivity will go in this
                              // bin. However, since there's an offset
                              // in the emissivities, we shift the
                              // index by that

  /* If it's gonna be in a bin of interest carry on - otherwise leave it. */

  // printout("frequency %g\n", freq_out*H/MEV);
  // printout("lindex %d, emiss_max %d, emiss_offset %d\n", lindex, emiss_max, emiss_offset);

  if ((lindex > emiss_offset - 1) && (lindex < emiss_offset + emiss_max - 1))
  {
    
    /* Then get partial crossection dsigma_domega in cmf */
    /* Coeff is 3 / 16 / PI */
      
    dsigma_domega_cmf = 0.0596831 * SIGMA_T / f / f *
      (f + (1./f) + (mu_cmf*mu_cmf) - 1.);
    
    //speed = vec_len(vel_vec);
    //solid_angle_factor =  doppler(pkt_ptr->dir, vel_vec) * doppler(pkt_ptr->dir, vel_vec); 
    
    /* 
      pow((1 + (dot(vel_vec, syn_dir)/CLIGHT)),2)
      / (1.0 - (speed* speed / CLIGHT / CLIGHT));
    */
    
    //dsigma_domega_rf = dsigma_domega_cmf //* doppler(pkt_ptr->dir, vel_vec)
    //* solid_angle_factor;
    
    /* so now determine the contribution to the emissivity and which 
 frequency bin it should be in */
	
    dop_fac = doppler(pkt_ptr->dir, vel_vec);
	
    emiss_cont = pkt_ptr->e_rf * dsigma_domega_cmf * dist * dop_fac * dop_fac / f;
	
    /* For normalisation this needs to be
       1) divided by volume
       2) divided by frequency bin size
       3) multiplied by the cell electron number density
       4) divided by the length of the time step
       This will all be done later
    */
    
    if (lindex < emiss_offset)
    {
      printout("scarily bad error here! %d %d\n", lindex, emiss_offset);
    }
    else
    {
      #ifdef _OPENMP
        #pragma omp atomic
      #endif
      compton_emiss[cell[pkt_ptr->where].modelgridindex][lindex - emiss_offset] += emiss_cont;
    }
	
  }

  return(0);
}

/*******************************************************/
int pp_emiss_cont(pkt_ptr, dist, t_current)
     PKT *pkt_ptr;
     double dist, t_current;
{
  /* New routine for getting a pair production emissivity. Closely based on compton_emiss but simpler. The 
     emissivity itself is stored in the last row of the compton emissivity structure. Idea here is to get something 
     which, when normalised by the volume and time step, will give the energy going into the .511 MeV 
     gamma rays from pair production per unit volume per unit time in the cmf. */
  /* Called with a packet that is about to travel a 
     distance dist in the lab frame. Time at start of distance is t_current.*/


  double emiss_cont;
  double sig_pair_prod();
  
  emiss_cont = sig_pair_prod(pkt_ptr, t_current) * (2.46636e+20 / pkt_ptr->nu_cmf) * pkt_ptr->e_rf * dist;
    
  /* For normalisation this needs to be
     1) divided by volume
     2) divided by the length of the time step
     This will all be done later
  */
      
  #ifdef _OPENMP 
    #pragma omp atomic
  #endif
  compton_emiss[cell[pkt_ptr->where].modelgridindex][emiss_max - 1] += 1.e-20 * emiss_cont;
  
  //  printf("emiss_cont %g\n", emiss_cont);

  /* Note (SS May 07) - the Doppler factors are not all sorted out yet - the expression used above needs to be
     consistent with what syn_lc does. */

  return(0);
}


/***********************************************/
int zero_estimators()
{
  int i,n,m,element,ion;

  //for (n=0; n < ngrid; n++)
  for (n=0; n < npts_model; n++)
  {
    J[n] = 0.;
    #ifndef FORCE_LTE
      nuJ[n] = 0.;
      ffheatingestimator[n] = 0.;
      colheatingestimator[n] = 0.;
      /*
      mabfcount[n] = 0.;
      mabfcount_thermal[n] = 0.;
      matotem[n] = 0.;
      maabs[n] = 0.;
      kbfcount[n] = 0.;
      kbfcount_ion[n] = 0.;
      kffcount[n] = 0.;
      kffabs[n] = 0.;
      kbfabs[n] = 0.;
      kgammadep[n] = 0.;
      */
      for (element = 0; element < nelements; element++)
      {
        for (ion = 0; ion < maxion; ion++)
        {
          gammaestimator[n*nelements*maxion+element*maxion+ion] = 0.;
          bfheatingestimator[n*nelements*maxion+element*maxion+ion] = 0.;
          /*
          photoionestimator[n*nelements*maxion+element*maxion+ion] = 0.;
          stimrecombestimator[n*nelements*maxion+element*maxion+ion] = 0.;
          ionfluxestimator[n*nelements*maxion+element*maxion+ion] = 0.;
          //twiddle[n*nelements*maxion+element*maxion+ion] = 0.;
          */
        }
      }
    #endif
    // cell[n].heating_ff = 0.;
    // cell[n].heating_bf = 0.;
    for (m=0; m < emiss_max; m++)
    {
      compton_emiss[n][m] = 0.0;
    }
    
    rpkt_emiss[n] = 0.0;
  }
  
  return(0);
}


/******************************************/
int normalise_estimators(int nts)
{
  int n,m;
  int assoc_cells;
  double vol_init(CELL *pkt_ptr);
  double dfreq[EMISS_MAX], get_gam_freq();
  double time_factor;
  double volume;
  
  time_factor = 1. / pow(time_step[nts].mid / tmin, 3.0) / time_step[nts].width;

  for (m=0; m < emiss_max; m++)
  {
    dfreq[m] = get_gam_freq(&gam_line_list, m + emiss_offset+1) - 
      get_gam_freq(&gam_line_list, m + emiss_offset);
    if (dfreq[m] < 0)
    {
      printout("Problem with normalisation of estimators. Abort.\n");
      exit(0);
    }
    dfreq[m] = 1./dfreq[m]; 
  }
  
  //for (n=0; n < ngrid; n++)
  for (n=0; n < npts_model; n++)
  {
    assoc_cells = modelgrid[n].associated_cells;
    volume = 1. / vol_init(&cell[n]);  ///That's not going to work if the parameter matters!!!!!!!!!!!!!!!!!!
    for (m=0; m < emiss_max; m++)
    {
      compton_emiss[n][m] = compton_emiss[n][m] * time_factor * volume / nprocs / assoc_cells;
      
      if (m < emiss_max - 1)
      /** (emiss_max - 1) contains the pair production case so it doesn't need the nne nor the dfreq */
      {
        compton_emiss[n][m] = compton_emiss[n][m] * get_nne(n) * dfreq[m];
      }
    }
  }
  return(0);
}


/*************************************************/
int write_estimators(nts)
     int nts;
{
  int n,m;
  FILE *est_file, *dummy;
  char chch;
  char filename[100] = "est_";
  char junk[100];
  int i;
  float dum;

  if ((dummy = fopen("dummy", "w+")) == NULL)
  {
    printout("Cannot open dummy.\n");
    exit(0);
  }
  fprintf(dummy, "%d", nts);
  fclose(dummy);
  if ((dummy = fopen("dummy", "r")) == NULL)
  {
    printout("Cannot open dummy.\n");
    exit(0);
  }
  i=0;
  while ((chch=fgetc(dummy)) != EOF)
  {
    junk[i] = chch;
    i= i+1;
  }
  junk[i] = '\0';
  fclose(dummy);

  strcat(filename, junk);
  strcat(filename, ".out");
  
  if (file_set == 1)
  {
    if ((est_file = fopen(filename, "rb")) == NULL)
    {
      printout("Cannot open est_file.txt.\n");
      exit(0);
    }

      //for (n=0; n < ngrid; n++)
    for (n=0; n < npts_model; n++)
    {
      for (m=0; m < emiss_max; m++)
	    {
	      fread(&dum, sizeof(float), 1, est_file);
	      //fscanf(est_file, "%g", &dum);
	      compton_emiss[n][m] += dum;
	    }
    }
    fclose(est_file);
  }


  if ((est_file = fopen(filename, "wb+")) == NULL)
  {
    printout("Cannot open est_file.txt.\n");
    exit(0);
  }


  //for (n=0; n < ngrid; n++)
  for (n=0; n < npts_model; n++)
  {
    for (m=0; m < emiss_max; m++)
    {
      fwrite(&compton_emiss[n][m], sizeof(float), 1, est_file);
    }
  }
  fclose(est_file);
  return(0);
}

/***********************************************/
int estim_switch(nts)
     int nts; //index for time step
{
  int on_or_off;
  double tstart, tend;
  double ts_want, te_want;

  on_or_off = 1; //on

  tstart = time_step[nts].start;
  tend = time_step[nts].start + time_step[nts].width;

  ts_want = time_syn[0] * ((1. - rmax/tmin/CLIGHT_PROP));
  te_want = time_syn[nsyn_time-1] * (1. + rmax/tmin/CLIGHT_PROP);  

  if (tstart > te_want)
  {
    on_or_off = 0;
  }

  if (tend < ts_want)
  {
    on_or_off = 0;
  }

  return(on_or_off);
}
    
/*************************************************/
int emiss_load(nts)
     int nts;
{
  /* Routine to read in the stored estimators for the time step that is about to begin. */
  int n,m;
  FILE *est_file, *dummy;
  char chch;
  char filename[100] = "est_";
  char junk[100];
  int i;
  float dum;

  if ((dummy = fopen("dummy", "w+")) == NULL)
  {
    printout("Cannot open dummy.\n");
    exit(0);
  }
  fprintf(dummy, "%d", nts);
  fclose(dummy);
  if ((dummy = fopen("dummy", "r")) == NULL)
  {
    printout("Cannot open dummy.\n");
    exit(0);
  }
  i=0;
  while ((chch=fgetc(dummy)) != EOF)
  {
    junk[i] = chch;
    i= i+1;
  }
  junk[i] = '\0';
  fclose(dummy);

  strcat(filename, junk);
  strcat(filename, ".out");
  
  if ((est_file = fopen(filename, "r")) == NULL)
  {
    printout("Cannot open est_file.txt.\n");
    exit(0);
  }


  //for (n=0; n < ngrid; n++)
  for (n=0; n < npts_model; n++)
  {
    for (m=0; m < emiss_max; m++)
    {
      fscanf(est_file, "%g", &dum);
      compton_emiss[n][m] = dum;
    }
  }
  fclose(est_file);
  return(0);
}
/***********************************************/
