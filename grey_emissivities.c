#include "sn3d.h"
#include <string.h>

///In this file the call to kappa_rpkt does not fit kappa_rpkts definition any loger!!!
///This MUST BE CHANGED. But it's not only the the call of kappa_rpkt. The dummy packet
///pointer needs more information (e.g. frequency) to calculate kappa_rpkt for the non
///grey case.

/*******************************************************/
int rlc_emiss_gamma(PKT *pkt_ptr, double dist, double t_current)
{
  /* Subroutine to record the heating rate in a cell due to gamma rays.
By heating rate I mean, for now, really the rate at which the code is making
k-packets in that cell which will then convert into r-packets. This is (going 
to be) used for the new light_curve syn-style calculation. */

  /* The intention is that rpkt_emiss will contain the emissivity of r-packets 
     in the co-moving frame (which is going to be isotropic). */

  /* This is only done to order v/c for now. */

  /* Called with a packet that is about to travel a
distance dist in the lab frame. Time at start of distance is t_current.*/

  double vel_vec[3];
  double dot();
  int get_velocity();
  double heating_cont;
  double xx;
  double sig_photo_electric(), meanf_sigma(), sig_pair_prod();
  double sdist, boundary_cross();
  double kap_tot, tau_tot, calculate_kappa_rpkt_cont();
  double t_future;
  int end_packet;
  int snext;
  PKT dummy;
  int move_pkt(PKT *pkt_ptr, double distance, double time);
  int change_cell();
  double tautau, mumu;

  //printout("[debug] Execution of rlc_emiss_gamma\n");
  
  dummy.pos[0] = pkt_ptr->pos[0];
  dummy.pos[1] = pkt_ptr->pos[1];
  dummy.pos[2] = pkt_ptr->pos[2];
  dummy.dir[0] = syn_dir[0];
  dummy.dir[1] = syn_dir[1];
  dummy.dir[2] = syn_dir[2];
  dummy.where = pkt_ptr->where;
  dummy.last_cross = NONE;
  
  end_packet = 0;
  tau_tot = 0;
  sdist = 0;
  t_future = t_current;
  tautau = 1.e-8;
  
  int mgi = cell[pkt_ptr->where].modelgridindex;

  if (dist > 0)
    {
      
      /* approach for opacity weighted estimators. */
      if (do_rlc_est == 2)
	{
	  mumu = dot(pkt_ptr->dir, syn_dir);
	  tautau = dist * calculate_kappa_rpkt_cont(&dummy, t_current) * mumu;
	  while (end_packet == 0 && tau_tot < 15.)
	    {
	      sdist = boundary_cross(&dummy, t_future, &snext);
	      kap_tot = calculate_kappa_rpkt_cont(&dummy, t_future);
	      tau_tot += kap_tot * sdist * t_current * t_current * t_current
		/ (t_future * t_future * t_future);
	      t_future += (sdist / CLIGHT_PROP);
	      move_pkt(&dummy, sdist, t_future);
	      change_cell(&dummy, snext, &end_packet, t_future);
	    }
          nesc -= 1;
	}
      
      if (tau_tot < 15.)
	{
	  
	  get_velocity(pkt_ptr->pos, vel_vec, t_current);
	  xx = H * pkt_ptr->nu_cmf / ME / CLIGHT / CLIGHT;
	  
	  heating_cont = ((meanf_sigma(xx)*get_nnetot(mgi)) + sig_photo_electric(pkt_ptr,t_current) + (sig_pair_prod(pkt_ptr, t_current) * (1. - (2.46636e+20 / pkt_ptr->nu_cmf))));
	  heating_cont = heating_cont * pkt_ptr->e_rf * dist * (1. - (2.*dot(vel_vec, pkt_ptr->dir)/CLIGHT));
	  
	  /* The terms in the above are for Compton, photoelectric and pair production. The pair production  one assumes that a fraction (1. - (1.022 MeV / nu)) of the gamma's energy is thermalised. The remaining 1.022 MeV is made into gamma rays */

	  /* For normalisation this needs to be
	     1) divided by volume
	     2) divided by the length of the time step
	     3) divided by 4 pi sr
	     This will all be done later
	  */
	  
	  
	  //  printout("%g %g (1)\n",heating_cont,tautau);
	  if (tautau > 0.001)
	    {
              #ifdef _OPENMP 
                #pragma omp atomic
              #endif
	      rpkt_emiss[mgi] += 1.e-20 * heating_cont * exp( -1. * (tau_tot + (0.5*tautau))) * (exp(tautau) - 1.0) / tautau;
	    }
	  else
	    {
              #ifdef _OPENMP 
                #pragma omp atomic
              #endif
	      rpkt_emiss[mgi] += 1.e-20 * heating_cont * exp( -1. * tau_tot);
	    }
      
	}

    }
  return(0);
}

/***********************************************/
/*******************************************************/
int
rlc_emiss_rpkt(pkt_ptr, dist, t_current)
     PKT *pkt_ptr;
     double dist, t_current;
{
  /* Subroutine to record the rate of destruction (and re-creation) of
     r-packets by the grey opacity. */

  /* This is only done to order v/c for now. */

  /* Called with a packet that is about to travel a
     distance dist in the lab frame. Time at start of distance is t_current.*/

  double vel_vec[3];
  double dot();
  int get_velocity();
  double cont;
  double sdist, boundary_cross();
  double kap_tot, calculate_kappa_rpkt_cont();
  int snext;
  PKT dummy;
  int end_packet;
  double tau_tot, t_future;
  int move_pkt(PKT *pkt_ptr, double distance, double time);
  int change_cell();
  double tautau, mumu;


  dummy.pos[0] = pkt_ptr->pos[0];
  dummy.pos[1] = pkt_ptr->pos[1];
  dummy.pos[2] = pkt_ptr->pos[2];
  dummy.dir[0] = syn_dir[0];
  dummy.dir[1] = syn_dir[1];
  dummy.dir[2] = syn_dir[2];
  dummy.where = pkt_ptr->where;
  dummy.last_cross = NONE;


  end_packet = 0;
  tau_tot = 0;
  sdist = 0;
  t_future = t_current;
  tautau = 1.e-8;

  int mgi = cell[pkt_ptr->where].modelgridindex;

  if (dist > 0.0)
    {
      /* for the weighted estimators version */
      if (do_rlc_est == 2)
	{	  
	  mumu = dot(pkt_ptr->dir, syn_dir);
	  tautau = dist * calculate_kappa_rpkt_cont(&dummy, t_current) * mumu;
	  while (end_packet == 0 && tau_tot < 15.)
	    {
	      sdist = boundary_cross(&dummy, t_future, &snext);
	      kap_tot = calculate_kappa_rpkt_cont(&dummy, t_future);
	      tau_tot += kap_tot * sdist * t_current * t_current * t_current
		/ (t_future * t_future * t_future);
	      t_future += (sdist / CLIGHT_PROP);
	      move_pkt(&dummy, sdist, t_future);
	      change_cell(&dummy, snext, &end_packet, t_future);
	    }
          nesc -= 1;
	}
      
      if (tau_tot < 15.)
	{
	  get_velocity(pkt_ptr->pos, vel_vec, t_current);
	  
	  cont = (get_kappagrey(mgi) * get_rho(mgi));
	  cont = cont * pkt_ptr->e_rf * dist * (1. - (2.*dot(vel_vec, pkt_ptr->dir)/CLIGHT));
	  
	  
	  /* For normalisation this needs to be
	     1) divided by volume
	     2) divided by the length of the time step
	     3) divided by 4 pi sr
	     This will all be done later
	  */
	  
	  //printout("%g %g(2)\n",tautau,cont);
	  if (tautau > 0.001)
	    {
              #ifdef _OPENMP 
                #pragma omp atomic
              #endif
	      rpkt_emiss[mgi] += 1.e-20 * cont * exp( -1. * (tau_tot + (0.5*tautau))) * (exp(tautau) - 1.0) / tautau;
	    }
	  else
	    {
              #ifdef _OPENMP 
                #pragma omp atomic
              #endif
	      rpkt_emiss[mgi] += 1.e-20 * cont * exp( -1. * tau_tot);
	    }
	  
	  
	}
      
    }
  return(0);
}


/***********************************************/
int normalise_grey(int nts)
{
  int n;
  int assoc_cells;
  double vol_init(CELL *pkt_ptr);
  double dV,dt,helper;

  //for (n=0; n < ngrid; n++)
  dt = time_step[nts].width;
  helper = pow(time_step[nts].mid / tmin, 3.0);
  for (n=0; n < npts_model; n++)
  {
    assoc_cells = modelgrid[n].associated_cells;
    dV = vol_init(&cell[n]) * helper;
    
    rpkt_emiss[n] = rpkt_emiss[n] * ONEOVER4PI / dV / dt / nprocs / assoc_cells;
  }
  return(0);
}


/*************************************************/

int write_grey(int nts)
{
  int n;
  FILE *est_file, *dummy;
  char chch;
  char filename[100] = "grey_est_";
  char junk[100];
  int i;
  float dum;

  if ((dummy = fopen("dummy", "w+")) == NULL){
    printout("Cannot open dummy.\n");
    exit(0);
  }
  fprintf(dummy, "%d", nts);
  fclose(dummy);
  if ((dummy = fopen("dummy", "r")) == NULL){
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
      if ((est_file = fopen(filename, "r")) == NULL){
	printout("Cannot open grey_est_file.txt.\n");
	exit(0);
      }


      //for (n=0; n < ngrid; n++)
      for (n=0; n < npts_model; n++)
	{
	  fscanf(est_file, "%g", &dum);
	  rpkt_emiss[n] += dum;
	}
      fclose(est_file);
    }



  if ((est_file = fopen(filename, "w+")) == NULL){
    printout("Cannot open grey_est_file.txt.\n");
    exit(0);
  }


  //for (n=0; n < ngrid; n++)
  for (n=0; n < npts_model; n++)
    {
      fprintf(est_file, " %g\n ", rpkt_emiss[n]);
    }
  fclose(est_file);
  return(0);
}
/***********************************************/

/* Routine to compute the mean energy converted to non-thermal electrons times
the Klein-Nishina cross section. */

double
meanf_sigma(x)
     double x;
{
  double term0, term1, term2, term3, term4;
  double tot;
  double f;
  
  f=1+(2*x);

  term0 = 2/x;
  term1 = ( 1 - (2/x) - (3/(x*x)) ) * log(f);
  term2 = ( (4 / x) + (3/(x*x)) - 1) * 2 * x / f;
  term3 = ( 1 - (2/x) - (1/(x*x))) * 2 * x *(1 + x) / f / f;
  term4 =  -2. * x * ((4*x*x) + (6*x) + 3) / 3 / f / f / f;

  tot = 3 * SIGMA_T * (term0 + term1 + term2 + term3 + term4) / (8 * x);

  return(tot);

}

/**************************************************************/
/*************************************************/
int
emiss_rlc_load(nts)
     int nts;
{
  /* Routine to read in the stored estimators for the time step that is about to begin. */
  int n;
  FILE *est_file, *dummy;
  char chch;
  char filename[100] = "grey_est_";
  char junk[100];
  int i;
  float dum;

  if ((dummy = fopen("dummy", "w+")) == NULL){
    printout("Cannot open dummy.\n");
    exit(0);
  }
  fprintf(dummy, "%d", nts);
  fclose(dummy);
  if ((dummy = fopen("dummy", "r")) == NULL){
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
  
  if ((est_file = fopen(filename, "r")) == NULL){
    printout("Cannot open est_file.txt.\n");
    exit(0);
  }


  //for (n=0; n < ngrid; n++)
  for (n=0; n < npts_model; n++)
    {
      fscanf(est_file, "%g", &dum);
      rpkt_emiss[n] = dum;
    }
  fclose(est_file);
  return(0);
}
/***********************************************/
int 
grey_rt(ray_ptr, nray, ldist, single_pos, single_t, lindex)
     RAY *ray_ptr;
     int nray, lindex;
     double ldist, single_t;
     double *single_pos;
{
  /* This is called when a grey ray is about to be moved a distance ldist. */
  /* It should account for the changes in the ray intensity due to
     grey processes along the path. */
  PKT dummy;
  double kap_tot, tau_cont;
  double calculate_kappa_rpkt_cont();
  double vel_vec[3];
  int get_velocity();
  double doppler();
  double vol_init();
  double vec_len();

  /* Make a dummy packet that carries the ray properties. */

  dummy.pos[0] = single_pos[0];
  dummy.pos[1] = single_pos[1];
  dummy.pos[2] = single_pos[2];

  dummy.dir[0] = syn_dir[0];
  dummy.dir[1] = syn_dir[1];
  dummy.dir[2] = syn_dir[2];
  
  dummy.where = ray_ptr->where;
  dummy.nu_cmf = ray_ptr->nu_cmf[nray];

  get_velocity(single_pos, vel_vec, single_t);

  if (do_rlc_est == 1)
    {

      kap_tot = calculate_kappa_rpkt_cont(&dummy, single_t);
      tau_cont = kap_tot * ldist;
      ray_ptr->e_rf[nray] = ray_ptr->e_rf[nray] * exp(-1. * tau_cont);
   
    }
  else if (do_rlc_est == 2)
    {

      tau_cont = kap_tot = 0;

    }
  else
    {
      printout("Unknown rlc type. Abort.\n");
      exit(0);
    }
  

   /* Now adding the emissivity term. */
  /* I think it's a doppler^3 term because it's a doppler^2 term in Castor + there's an 
     extra doppler from the integration over frequency of the RT equation. */
  
   if (tau_cont > 1.e-6)
    {
      ray_ptr->e_rf[nray] += (rpkt_emiss[cell[dummy.where].modelgridindex] / 
			      doppler(syn_dir, vel_vec) / doppler(syn_dir, vel_vec) / doppler(syn_dir, vel_vec) 
			      *(1. - exp(-1. * tau_cont)) / kap_tot);
    }
  else
    {
      ray_ptr->e_rf[nray] += (rpkt_emiss[cell[dummy.where].modelgridindex]   /
			      doppler(syn_dir, vel_vec) / doppler(syn_dir, vel_vec) / doppler(syn_dir, vel_vec) 
			      * ldist);
    }	  


  /* This MUST be followed by a call to move_one_ray() in source
     since e_cmf is NOT reset here. */

  return(0);

}
