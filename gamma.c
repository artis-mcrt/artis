#include "sn3d.h"

/* Material for handing gamma rays - creation and propagation. */

int
pellet_decay(nts,pkt_ptr)
     int nts;
     PKT *pkt_ptr;
{

  double zrand;
  double zrand2;
  double mu, phi, sintheta;
  double vec_len();
  double dummy_dir[3], vel_vec[3];
  int get_velocity();
  int angle_ab();
  int choose_gamma_ray();
  int cross_prod();
  double dot();
  double doppler();
  void vec_norm(double x[3], double z[3]);

  /*Subroutine to convert a pellet to a gamma ray. */
  /*nts defines the time step we are in. pkt_ptr is a pointer to the packet
  that is decaying. */
  /* Record decay. */
  #ifdef _OPENMP 
    #pragma omp atomic
  #endif
  time_step[nts].pellet_decays += 1;

  /*Start by getting the position of the pellet at the point of decay. Pellet
  is moving with the matter.*/
 
  pkt_ptr->pos[0] = pkt_ptr->pos[0] * pkt_ptr->tdecay / time_step[nts].start;
  pkt_ptr->pos[1] = pkt_ptr->pos[1] * pkt_ptr->tdecay / time_step[nts].start;
  pkt_ptr->pos[2] = pkt_ptr->pos[2] * pkt_ptr->tdecay / time_step[nts].start;
  
  /*Now let's give the gamma ray a direction. */

  zrand = gsl_rng_uniform(rng);
  zrand2 = gsl_rng_uniform(rng);

  /*Assuming isotropic emission in cmf, use these two random numbers to set 
  up a cmf direction in cos(theta) and phi. */

  mu = -1 + (2.*zrand);
  phi = zrand2 * 2 * PI;
  sintheta = sqrt(1. - (mu * mu));

  pkt_ptr->dir[0] = sintheta * cos(phi);
  pkt_ptr->dir[1] = sintheta * sin(phi);
  pkt_ptr->dir[2] = mu;
  
  /* This direction is in the cmf - we want to convert it to the rest 
  frame - use aberation of angles. We want to convert from cmf to
  rest so need -ve velocity. */

  get_velocity(pkt_ptr->pos, vel_vec, (-1*(pkt_ptr->tdecay)));
  //negative time since we want the backwards transformation here
 
  angle_ab(pkt_ptr->dir, vel_vec, dummy_dir);

  pkt_ptr->dir[0] = dummy_dir[0];
  pkt_ptr->dir[1] = dummy_dir[1];
  pkt_ptr->dir[2] = dummy_dir[2];

  /*Check unit vector.*/
  if (fabs(vec_len(pkt_ptr->dir) - 1) > 1.e-8)
  {
    printout("Not a unit vector. Abort.\n");
    exit(0);
  }

  /* Now need to assign the frequency of the packet in the co-moving frame.*/

  choose_gamma_ray(pkt_ptr);


  /* Finally we want to put in the rest frame energy and frequency. And record
  that it's now a gamma ray.*/


  get_velocity(pkt_ptr->pos, vel_vec, pkt_ptr->tdecay);

  pkt_ptr->nu_rf = pkt_ptr->nu_cmf / doppler(pkt_ptr->dir, vel_vec);
  pkt_ptr->e_rf = pkt_ptr->e_cmf * pkt_ptr->nu_rf /pkt_ptr->nu_cmf; 

  pkt_ptr->type = TYPE_GAMMA;
  pkt_ptr->last_cross = NONE;

  /**initialise polarisation information */
  pkt_ptr->stokes[0]=1.0;
  pkt_ptr->stokes[1]=pkt_ptr->stokes[2]=0.0;
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
  //printout("initialise pol state of packet %g, %g, %g, %g, %g\n",pkt_ptr->stokes[1],pkt_ptr->stokes[2],pkt_ptr->pol_dir[0],pkt_ptr->pol_dir[1],pkt_ptr->pol_dir[2]);
  //printout("pkt direction %g, %g, %g\n",pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);

  return(0);
}


/******************************************************/

int
choose_gamma_ray(pkt_ptr)
    PKT *pkt_ptr;
{
  double zrand, runtot;
  int n;
  /* Routine to choose which gamma ray line it'll be. */
  
  if (pkt_ptr->type == TYPE_NICKEL_PELLET)
  {
    n=0;
    runtot=0.0;
    zrand = gsl_rng_uniform(rng);
    while (zrand > runtot)
    {
      runtot += nickel_spec.probability[n] * nickel_spec.energy[n] / ENICKEL;
      n++;
    }
    n = n -1 ;
    if (n >= nickel_spec.nlines)
    {
      printout("Failure to choose line (Ni). Abort.\n");
      abort();
    }
    else
    {
      pkt_ptr->nu_cmf = nickel_spec.energy[n] / H;
    }
  }
  else if (pkt_ptr->type == TYPE_COBALT_PELLET)
  {
    n=0;
    runtot=0.0;
    zrand = gsl_rng_uniform(rng);
    while (zrand > runtot)
    {
      runtot += cobalt_spec.probability[n] * cobalt_spec.energy[n] / ECOBALT_GAMMA;
      n++;
    }
    n = n -1 ;
    if (n >= cobalt_spec.nlines)
    {
      printout("Failure to choose line (Co). Abort.\n");
      abort();
    }
    else
    {
      pkt_ptr->nu_cmf = cobalt_spec.energy[n] / H;
    }
  }
  else if (pkt_ptr->type == TYPE_48CR_PELLET)
  {
    n=0;
    runtot=0.0;
    zrand = gsl_rng_uniform(rng);
    while (zrand > runtot)
    {
      runtot += cr48_spec.probability[n] * cr48_spec.energy[n] / E48CR;
      n++;
    }
    n = n -1 ;
    if (n >= cr48_spec.nlines)
    {
      printout("Failure to choose line (Cr48). Abort.\n");
      abort();
    }
    else
    {
      pkt_ptr->nu_cmf = cr48_spec.energy[n] / H;
    }
  }
  else if (pkt_ptr->type == TYPE_48V_PELLET)
  {
    n=0;
    runtot=0.0;
    zrand = gsl_rng_uniform(rng);
    while (zrand > runtot)
    {
      runtot += v48_spec.probability[n] * v48_spec.energy[n] / E48V;
      n++;
    }
    n = n -1 ;
    if (n >= v48_spec.nlines)
    {
      printout("Failure to choose line (V48). Abort.\n");
      abort();
    }
    else
    {
      pkt_ptr->nu_cmf = v48_spec.energy[n] / H;
    }
  }
   else 
  {
    printout("Unrecognised pellet. Abort.\n");
    exit(0);
  }

  return(0);
}
	  
      
/********************************************************/

/* Now routine for moving a gamma packet. Idea is that we have as input
a gamma packet with known properties at time t1 and we want to follow it
until time t2. */

double
      do_gamma (pkt_ptr, t1, t2)
      PKT *pkt_ptr;
     double t1, t2;
{
  double t_current;
  double sdist, tdist,edist;
  int snext;
  double zrand;
  double boundary_cross();
  double tau_next, tau_current;
  int move_pkt(PKT *pkt_ptr, double distance, double time);
  int change_cell();
  int end_packet; //tells us when to stop working on this packet
  double kap_tot, kap_compton, kap_photo_electric, kap_pair_prod;
  double sig_comp();
  double sig_photo_electric();
  double sig_pair_prod();
  int com_sca(), compton_emiss_cont(), rlc_emiss_gamma(), pair_prod();
  int pp_emiss_cont();

  end_packet = 0; //means "keep working"

  t_current = t1; //this will keep track of time in the calculation

  while (end_packet == 0)
  {
      /* Assign optical depth to next physical event. And start counter of
    optical depth for this path.*/
    zrand = gsl_rng_uniform(rng);
    tau_next = -1. * log(zrand);
    tau_current = 0.0;

      /* Start by finding the distance to the crossing of the grid cell
    boundaries. sdist is the boundary distance and snext is the
    grid cell into which we pass.*/
      
    sdist = boundary_cross(pkt_ptr, t_current, &snext);
 
    if (sdist > (rmax * t_current/tmin))
    {
      printout("Unreasonably large sdist (gamma). Abort. %g %g %g\n", rmax, t_current/tmin, sdist);
      exit(0);
    }

    if (sdist < 1)
    {
      printout("Negative distance (sdist). Abort?\n");
      sdist =0;
    }
    if (((snext != -99) && (snext < 0)) || (snext >= ngrid))
    {
      printout("Heading for inappropriate grid cell. Abort.\n");
      printout("Current cell %d, target cell %d.\n", pkt_ptr->where, snext);
      exit(0);
    }
    if (sdist > max_path_step)
    {
      sdist = max_path_step;
      snext = pkt_ptr->where;
    }
      
    /* Now consider the scattering/destruction processes. */
    /* Compton scattering - need to determine the scattering co-efficient.*/
    /* Routine returns the value in the rest frame. */

    if (gamma_grey < 0)
    {
      kap_compton = sig_comp(pkt_ptr,t_current);
    }
    else
    {
      kap_compton = 0.0;
    }      
    kap_photo_electric = sig_photo_electric(pkt_ptr, t_current);
    kap_pair_prod = sig_pair_prod(pkt_ptr, t_current);
    kap_tot = kap_compton + kap_photo_electric + kap_pair_prod;

    /* So distance before physical event is...*/

    edist = (tau_next - tau_current) / kap_tot;

    if (edist < 0)
    {
      printout("Negative distance (edist). Abort. \n");
      exit(0);
    }


    /* Find how far it can travel during the time inverval. */
      
    tdist = (t2 - t_current) * CLIGHT_PROP;
      
    if (tdist < 0)
    {
      printout("Negative distance (tdist). Abort. \n");
      exit(0);
    }

      //printout("sdist, tdist, edist %g %g %g\n",sdist, tdist, edist);
      

    if ((sdist < tdist) && (sdist < edist))
    {
      sdist = sdist / 2.;
      t_current += sdist / CLIGHT_PROP;
      move_pkt(pkt_ptr,sdist,t_current);

      /* Move it into the new cell. */
      if (kap_tot > 0)
	{
	  if (do_comp_est == 1)
	    {
	      sdist = sdist * 2.;
	      compton_emiss_cont(pkt_ptr, sdist, t_current);
	      pp_emiss_cont(pkt_ptr, sdist, t_current);
	      sdist = sdist / 2.;
	    }
	  if (do_rlc_est != 0)
	    {
	      sdist = sdist * 2.;
	      rlc_emiss_gamma(pkt_ptr, sdist, t_current);
	      sdist = sdist / 2.;
	    }
	}

      t_current += sdist / CLIGHT_PROP;
      move_pkt(pkt_ptr,sdist,t_current);
      sdist = sdist * 2.;
      if (snext != pkt_ptr->where)
      {
        change_cell(pkt_ptr, snext, &end_packet, t_current);
      }
    }
    else if ((tdist < sdist) && (tdist < edist))
    {
      /* Doesn't reach boundary. */
      tdist = tdist / 2.;
      t_current += tdist  / CLIGHT_PROP;
      move_pkt(pkt_ptr,tdist,t_current);
	  
      if (kap_tot > 0)
	{
	  if (do_comp_est == 1)
	    {
	      tdist = tdist * 2.;
	      compton_emiss_cont(pkt_ptr, tdist, t_current);
	      pp_emiss_cont(pkt_ptr, tdist, t_current);
	      tdist = tdist / 2.;
	    }
	  if (do_rlc_est != 0)
	    {
	      tdist = tdist * 2.;
	      rlc_emiss_gamma(pkt_ptr, tdist, t_current);
	      tdist = tdist / 2.;
	    }
	}
      t_current = t2;
      move_pkt(pkt_ptr,tdist,t_current);
      tdist = tdist * 2.;
      end_packet = 1;
    }
    else if ((edist < sdist) && (edist < tdist))
    {
      edist = edist / 2.;
      t_current += edist / CLIGHT_PROP;
      move_pkt(pkt_ptr,edist,t_current);
      if (kap_tot > 0)
	{
	  if (do_comp_est == 1)
	    {
	      edist = edist * 2.;
	      compton_emiss_cont(pkt_ptr, edist, t_current);  
	      pp_emiss_cont(pkt_ptr, edist, t_current);
	      edist = edist / 2.;
	    }
	  if (do_rlc_est != 0)
	    {
	      edist = edist * 2.;
	      rlc_emiss_gamma(pkt_ptr, edist, t_current);
	      edist = edist / 2.;
	    }
	}
      t_current += edist / CLIGHT_PROP;
      move_pkt(pkt_ptr,edist,t_current);
      edist = edist * 2.;


      /* event occurs. Choose which event and call the appropriate subroutine.*/
      zrand = gsl_rng_uniform(rng);
      if (kap_compton > (zrand*kap_tot))
      {
        /* Compton scattering. */
        com_sca(pkt_ptr,t_current);
        if (pkt_ptr->type != TYPE_GAMMA)
        {
          /* It's not a gamma ray any more - return.*/
          return(t_current);
        }
      }
      else if ((kap_compton + kap_photo_electric) > (zrand*kap_tot))
      {
        /* Photo electric effect - makes it a k-packet for sure.*/
        pkt_ptr->type = TYPE_KPKT;
        pkt_ptr->absorptiontype = -4;
        #ifndef FORCE_LTE
          //kgammadep[pkt_ptr->where] += pkt_ptr->e_cmf;
        #endif
        //pkt_ptr->type = TYPE_PRE_KPKT;
        //pkt_ptr->type = TYPE_GAMMA_KPKT;
        //if (tid == 0) k_stat_from_gamma += 1;
        k_stat_from_gamma += 1;
        return(t_current);
      }
      else if ((kap_compton + kap_photo_electric + kap_pair_prod) > (zrand*kap_tot))
	{
	  /*It's a pair production */
	  pair_prod(pkt_ptr,t_current);
	  if (pkt_ptr->type != TYPE_GAMMA)
	    {
	      /* It's not a gamma ray any more - return.*/
	      return(t_current);
	    }
	}
      else 
	{
        printout("Failed to identify event. Gamma (1). kap_compton %g kap_photo_electric %g kap_tot %g zrand %g Abort.\n", kap_compton, kap_photo_electric, kap_tot, zrand);
        printout(" /*cell[*/pkt_ptr->where].rho %g pkt_ptr->nu_cmf %g pkt_ptr->dir[0] %g pkt_ptr->dir[1] %g pkt_ptr->dir[2] %g pkt_ptr->pos[0] %g pkt_ptr->pos[1] %g pkt_ptr->pos[2] %g \n",get_rho(cell[pkt_ptr->where].modelgridindex), pkt_ptr->nu_cmf,pkt_ptr->dir[0],pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2],pkt_ptr->pos[1],pkt_ptr->pos[2]);


        exit(0);
      }
    }      
    else
    {
      printout("Failed to identify event. Gamma (2). edist %g, sdist %g, tdist %g Abort.\n", edist, sdist, tdist);
      exit(0);
    }
  }
  return(PACKET_SAME);
}
