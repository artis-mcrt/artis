#include "sn3d.h"

/* Master routine for moving along ray. When it called,
   it is given the time at start of inverval and at end - when it finishes,
   everything the ray passes during this time should be sorted out. */

int ray_prop(ray_ptr, t1, t2, nts)
     RAY *ray_ptr;
     double t1, t2;
     int nts; //the time step we are doing
{
  double t_current; 
  double do_gamma_ray();
  int end_packet;
  double t_change_type;

  end_packet = 0;  //means "keep working"

  t_current = t1;

  while (end_packet == 0)
    {
      /* Start by sorting out what sort of packet it is.*/
      if (ray_ptr->status == ACTIVE)
	{
	  /*It's a gamma-ray ray.*/
	  /* Call do_gamma_ray. */
	  t_change_type = do_gamma_ray( ray_ptr, t_current, t2);
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
	      printout("Error type change for a ray??\n");
	      exit(0);
	    }
	}
      else 
	{
	  printout("Trying to process an inactive ray??\n");
	  exit(0);
	}

    }

  return(0);
}

/***************************************************************/

double
do_gamma_ray(ray_ptr, t1, t2)
     RAY *ray_ptr;
     double t1, t2;
{
  double t_current;
  double sdist, tdist, ldist, trav_dist;
  int snext;
  double boundary_cross_ray();
  int end_packet; //tells us when to stop working on this ray
  RAY dum_ray;
  double stop_dist, stop_time;
  double dnuds[NSYN];
  int get_nul(), lindex;
  int copy_ray();
  int nray;
  double single_pos[3], single_t;
  double get_gam_freq();
  int move_one_ray();
  int move_ray();
  int change_cell_ray();
  int add_gam_line_emissivity();
  int continuum_rt();
  int grey_rt();

  end_packet = 0; //means "keep working"

  t_current = t1; //this will keep track of time in the calculation

  for (nray = 0; nray < NSYN; nray++)
    {
      /* need to start by identifying next spectral line in list. */
      
      lindex = get_nul(ray_ptr->nu_cmf[nray]);
      ray_ptr->lindex[nray] = lindex;
    }

  while (end_packet == 0)
    {      /* Start by finding the distance to the crossing of the grid cell
	 boundaries. sdist is the boundary distance and snext is the
	 grid cell into which we pass.*/
      
      sdist = boundary_cross_ray(ray_ptr, t_current, &snext);
 
      if (sdist > (rmax * t_current/tmin))
	{
	  printout("Unreasonably large sdist (gamma). Abort. %g %g %g\n", rmax, t_current/tmin, sdist);
	  exit(0);
	}

      if (sdist < 1)
	{
	  printout("Negative distance (sdist). Abort.\n");
	  printout("sdist %g\n", sdist);
	  exit(0);
	}
      if (((snext != -99) && (snext < 0)) || (snext >= ngrid))
	{
	  printout("Heading for inappropriate grid cell. Abort.\n");
	  printout("Current cell %d, target cell %d.\n", ray_ptr->where, snext);
	  exit(0);
	}

 
      /* Find how far it can travel during the time inverval. */
      
      tdist = (t2 - t_current) * CLIGHT_PROP;
      
      if (tdist < 0)
	{
	  printout("Negative distance (tdist). Abort. \n");
	  exit(0);
	}

      if (tdist < sdist)
	{
	  stop_dist = tdist;
	  stop_time = t2 - t_current;
	}
      else
	{
	  stop_dist = sdist;
	  stop_time = sdist / CLIGHT_PROP;
	}



      /* Now consider the gamma-ray emission processes. Need to find distance
	 to next line in the gamma-ray list. Then use emissivity of that
	 line. */

      /* At this point we need to branch to loop over all the frequencies in the ray
	 bundle. They all need to travel the same sdist and all have the same time to 
	 travel. But some will reach lines/continuum events before others. */

      /* First get an array of DNUDS for all the rays in the bundle. Copy ray and move copy to next
	 boundary. */
 
      copy_ray(ray_ptr, &dum_ray);
      move_ray(&dum_ray, stop_dist, t_current + stop_time);
      for (nray = 0; nray < NSYN; nray++)
	{
	  dnuds[nray] = (dum_ray.nu_cmf[nray] - ray_ptr->nu_cmf[nray])/stop_dist;
	  if (dnuds[nray] >= 0)
	    {
	      printout("dnuds not negative? %g %d\n", dnuds[nray], nray);
	      printout("dum_ray.nu_cmf[nray] %g, ray_ptr->nu_cmf[nray] %g, stop_dist %g\n", dum_ray.nu_cmf[nray] , ray_ptr->nu_cmf[nray], stop_dist); 
	      exit(0);
	    }

	}

      for (nray = 0; nray < NSYN; nray++)
	{
	  /* need to start by identifying next spectral line in list. */

	  lindex = ray_ptr->lindex[nray];
	  single_pos[0] = ray_ptr->pos[0];
	  single_pos[1] = ray_ptr->pos[1];
	  single_pos[2] = ray_ptr->pos[2];
	  single_t = t_current;
	  
	  trav_dist = 0.0;

	  while (trav_dist < stop_dist)
	    {
	      ldist = (get_gam_freq(&gam_line_list, lindex) - ray_ptr->nu_cmf[nray]) / dnuds[nray];
	      
	      if (ldist < 0)
		{
		  printout("negative ldist - abort\n");
		  exit(0);
		}

	      /* then find distance to line - do we reach it before stop_dist. */
	      
	      if ((ldist + trav_dist) < stop_dist)
		{
		  /*It'll reach the line */
		  single_t = single_t + (ldist / CLIGHT_PROP);
		  if (do_rlc_est != 0)
		    {
		      printout("Ermmm. What?\n");
		      exit(0);
		      grey_rt(ray_ptr, nray, ldist, single_pos, single_t, lindex);
		    }
		  else
		    {
		      continuum_rt(ray_ptr, nray, ldist, single_pos, single_t, lindex);
		    }
		  move_one_ray(ray_ptr, nray, ldist, single_pos, single_t);
		  /* now interact with the line */
		  /* for example test only */

		  add_gam_line_emissivity(ray_ptr, nray, single_pos, single_t, lindex, dnuds[nray]);

		  /* end dodgy example lines */
		  trav_dist += ldist;
		  if (lindex > 0)
		    {
		      lindex = lindex - 1;
		    }
		  else
		    {
		      lindex = RED_OF_LIST;
		    }
		  
		}
	      else
		{
		  /* Won't reach the line - just need to move it to the boundary and then do next ray
		     in bundle. */
		  single_t = single_t + ((stop_dist - trav_dist) / CLIGHT_PROP);	  
		  if (do_rlc_est != 0)
		    {
		      grey_rt(ray_ptr, nray, stop_dist - trav_dist, single_pos, single_t, lindex);
		    }
		  else
		    {
		      continuum_rt(ray_ptr, nray, stop_dist - trav_dist, single_pos, single_t, lindex);
		    }
		  move_one_ray(ray_ptr, nray, stop_dist - trav_dist, single_pos, single_t);
		  trav_dist += stop_dist; //just to force break in while loop
		  ray_ptr->lindex[nray] = lindex;  
		  if (sdist < tdist)
		    {
		      if (fabs(1. - ray_ptr->nu_cmf[nray]/dum_ray.nu_cmf[nray]) > 1.e-8)
			{
			  printout("nu_cmf doesn't match up %g %g\n", ray_ptr->nu_cmf[nray], dum_ray.nu_cmf[nray]);
			  exit(0);
			}
		    }

		}

	    }
	}

      /* At this point all rays in the bundle have reached the stop_dist. */
      /* Update the position of the ray and, if necessary, the "where" element.  */
      
      ray_ptr->pos[0] = dum_ray.pos[0];
      ray_ptr->pos[1] = dum_ray.pos[1];
      ray_ptr->pos[2] = dum_ray.pos[2];
      t_current += stop_time;

      if (tdist < sdist)
	{
	  end_packet = 1;
	}
      else
	{
	  change_cell_ray(ray_ptr, snext, &end_packet, t_current);
	}
    }

  return(PACKET_SAME);
}
 

/****************************************************/
double
boundary_cross_ray(ray_ptr, tstart, snext)
     RAY *ray_ptr;
     double tstart;
     int *snext;
{
  PKT dummy;
  double dist;
  double boundary_cross();

  /* This is just a front end to use boundary_cross with a ray. Puts relevant
     info in to a fake packet and passes. Then sets return in ray. */

  dummy.pos[0] = ray_ptr->pos[0];
  dummy.pos[1] = ray_ptr->pos[1];
  dummy.pos[2] = ray_ptr->pos[2];

  dummy.dir[0] = syn_dir[0];
  dummy.dir[1] = syn_dir[1];
  dummy.dir[2] = syn_dir[2];

  dummy.where = ray_ptr->where;
  dummy.last_cross = ray_ptr->last_cross;

  dist = boundary_cross (&dummy, tstart, snext);

  ray_ptr->last_cross = dummy.last_cross;

  return(dist);
}
  
/****************************************************/
int
change_cell_ray(ray_ptr, snext, end_packet, t_current)
     RAY *ray_ptr;
     double t_current;
     int snext;
     int *end_packet;
{
  PKT dummy;
  int change_cell();

  /* Similar to above. */

  dummy.pos[0] = ray_ptr->pos[0];
  dummy.pos[1] = ray_ptr->pos[1];
  dummy.pos[2] = ray_ptr->pos[2];

  dummy.type = TYPE_GAMMA;

  change_cell(&dummy, snext, end_packet, t_current);

  if (dummy.type == TYPE_GAMMA)
    {
      ray_ptr->where = dummy.where;
      ray_ptr->pos[0] = dummy.pos[0];
      ray_ptr->pos[1] = dummy.pos[1];
      ray_ptr->pos[2] = dummy.pos[2];
    }
  else if (dummy.type == TYPE_ESCAPE)
    {
      ray_ptr->status = FINISHED;
    }
      
  return(0);
}
  
/**************************************************/
int 
copy_ray (ray1, ray2)
     RAY *ray1, *ray2;
{
  int n;

  ray2->last_cross = ray1->last_cross;
  ray2->where = ray1->where;
  ray2->status = ray1->status;
  ray2->tstart = ray1->tstart;
  ray2->rstart[0] = ray1->rstart[0];
  ray2->rstart[1] = ray1->rstart[1];
  ray2->rstart[2] = ray1->rstart[2];
  ray2->pos[0] = ray1->pos[0];
  ray2->pos[1] = ray1->pos[1];
  ray2->pos[2] = ray1->pos[2];
  for (n = 0; n < NSYN; n++)
    { 
      ray2->nu_rf[n] = ray1->nu_rf[n];
      ray2->nu_cmf[n] = ray1->nu_cmf[n];
      ray2->e_rf[n] = ray1->e_rf[n];
      ray2->e_cmf[n] = ray1->e_cmf[n];
    }
  return(0);
}

/**********************************************/
int 
move_ray(ray_ptr, dist, time)
     RAY *ray_ptr;
     double dist;
     double time;
{
  /* just make a dummy packet and use move. */
  
  PKT dummy;
  int move_pkt(PKT *pkt_ptr, double distance, double time);
  double doppler_fac;
  int n;
  
  dummy.pos[0] = ray_ptr->pos[0];
  dummy.pos[1] = ray_ptr->pos[1];
  dummy.pos[2] = ray_ptr->pos[2];
  dummy.dir[0] = syn_dir[0];
  dummy.dir[1] = syn_dir[1];
  dummy.dir[2] = syn_dir[2];

  dummy.nu_cmf = ray_ptr->nu_cmf[0];
  dummy.nu_rf = ray_ptr->nu_rf[0];
  dummy.e_cmf = ray_ptr->e_cmf[0];
  dummy.e_rf = ray_ptr->e_rf[0];
  
  move_pkt(&dummy, dist, time);

  ray_ptr->pos[0] = dummy.pos[0];
  ray_ptr->pos[1] = dummy.pos[1];
  ray_ptr->pos[2] = dummy.pos[2];

  doppler_fac = dummy.nu_cmf / dummy.nu_rf;
  
  for (n = 0; n < NSYN; n++)
    {
      ray_ptr->nu_cmf[n] = ray_ptr->nu_rf[n] * doppler_fac;

      /* removing next line (Jan06) since e_cmf seems redundant for rays. */
      //      ray_ptr->e_cmf[n] = ray_ptr->e_rf[n] * doppler_fac;
    }
  return(0);
}
/**************************************************************/
int
move_one_ray(ray_ptr, nray, dist, single_pos, single_t)
     RAY *ray_ptr;
     int nray;
     double dist, single_t;
     double *single_pos;
{
  int get_velocity();
  double doppler();
  double vel_vec[3];
    
  if (dist < 0)
    {
      printout("Trying to move -v distance. Abort.\n");
      exit(0);
    }

  single_pos[0] += syn_dir[0] * dist;
  single_pos[1] += syn_dir[1] * dist;
  single_pos[2] += syn_dir[2] * dist;
  
  get_velocity(single_pos, vel_vec, single_t);

  ray_ptr->nu_cmf[nray] = ray_ptr->nu_rf[nray] * doppler(syn_dir, vel_vec);
  // again, rmoving next line since e_cmf seems redundant.
  //ray_ptr->e_cmf[nray] = ray_ptr->e_rf[nray] * ray_ptr->nu_cmf[nray] / ray_ptr->nu_rf[nray];

  return(0);
}



/**************************************************************/

int
get_nul(freq)
     double freq;
{
  int too_high, too_low, try;
  double freq_max, freq_min, freq_try;
  double get_gam_freq();

  freq_max = get_gam_freq(&gam_line_list, gam_line_list.total - 1);
  freq_min = get_gam_freq(&gam_line_list, 0);

  if (freq > freq_max)
    {
      return(gam_line_list.total-1);
    }
  else if (freq < freq_min)
    {
      return(RED_OF_LIST);
    }
  else
    {
      too_high = gam_line_list.total - 1;
      too_low = 0;

      while (too_high != too_low + 1)
	{
	  
	  try = (too_high + too_low)/2;
	  freq_try = get_gam_freq(&gam_line_list, try);
	  if (freq_try >= freq)
	    {
	      too_high = try;
	    }
	  else
	    {
	      too_low = try;
	    }
	}

      return(too_low);
    }
}
	
/**************************************************************/
double
get_gam_freq(line_list, n)
     LIST *line_list;
     int n;
{
  double freq;

  if (n == RED_OF_LIST)
    {
      return(0.0);
    }
  /* returns the frequency of line n */
  else if (line_list->type[n] == NI_GAM_LINE_ID)
    {
      freq = nickel_spec.energy[line_list->index[n]] / H;
    }
  else if (line_list->type[n] == CO_GAM_LINE_ID)
    {
      freq = cobalt_spec.energy[line_list->index[n]] / H;
    }
  else if (line_list->type[n] == FAKE_GAM_LINE_ID)
    {
      freq = fakeg_spec.energy[line_list->index[n]] / H;
    }
  else if (line_list->type[n] == CR48_GAM_LINE_ID)
    {
      freq = cr48_spec.energy[line_list->index[n]] / H;
    }
  else if (line_list->type[n] == V48_GAM_LINE_ID)
    {
      freq = v48_spec.energy[line_list->index[n]] / H;
    }
  else
    {
      printout("Unknown line. %d Abort.\n", n);
      printout("line_list->type[n] %d line_list->index[n] %d\n", line_list->type[n], line_list->index[n]);
      printout(" %d %d \n", gam_line_list.type[n], gam_line_list.index[n]);
      exit(0);
    }

  return(freq);

}
