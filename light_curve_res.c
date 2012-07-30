#include "sn3d.h"
#include "exspec.h"

/* Routine to make a MC light curve from the r-packets. */

int
make_light_curve_res()
{
  int gather_light_curve_res();
  int write_light_curve_res();
  
  gather_light_curve_res(0);
  write_light_curve_res();
  return(0);
}

void init_light_curve_res()
{
  int n,nn;
  
  if (ntlcbins > MTLCBINS)
  {
    printout("Too many time bins in light curve - reducing.\n");
    ntlcbins = MTLCBINS;
  }

  /* start by setting up the time bins. */ 
  /* it is all done interms of a logarithmic spacing in t - get the
     step sizes first. */
  ///Should be moved to input.c or exspec.c
  dlogtlc = (log(tmax) - log(tmin))/ntlcbins;

  for (nn = 0; nn < MALCBINS; nn++)
  {
    for (n = 0; n < ntlcbins; n++)
      {
        light_curve_res[n][nn].lower_time = exp( log(tmin) + (n * (dlogtlc)));
        light_curve_res[n][nn].delta_t = exp( log(tmin) + ((n+1) * (dlogtlc))) - light_curve[n].lower_time;
        light_curve_res[n][nn].lum = 0.0;
      }
  }
}

/************************************************************/

int
gather_light_curve_res(my_rank)
     int my_rank;
{
  void read_packets(FILE *packets_file);
  PKT *pkt_ptr;
  int add_to_lc_res();
  int i,n,p,nn;
 


  /// The grid is now set up. Now we loop over all the packets, check if they made it out or not,
  /// and if they did we add their rest frame energy to the appropriate cell.

  /// And figure out the escaping packets of those.
  for (p = 0; p < npkts; p++)
  {
    pkt_ptr = &pkt[p];
    if (pkt_ptr->type == TYPE_ESCAPE && pkt_ptr->escape_type == TYPE_RPKT)
    {
      /// It made it out.
      add_to_lc_res(pkt_ptr);
    }
  }
  return(0);
}

/*******************************************************************/
int write_light_curve_res()
{
  FILE *lc_file;
  int m,nn;
  double save[MTSTEP][MALCBINS];
  float dum1, dum2;


  
  /// Light curve is done - write it out.
  /// If needed, start by reading in existing file and storing old numbers.
  if (file_set == 1)
    {
      if ((lc_file = fopen("light_curve_res.out", "r")) == NULL){
	printout("Cannot open lc_file.txt.\n");
	exit(0);
      }
      for (nn=0; nn < MALCBINS; nn++)
	{
	  for (m=0; m < ntlcbins; m++)
	    {
	      fscanf(lc_file, "%g %g\n", &dum1, &dum2);
	      save[m][nn]=dum2;
	    }
	}

      fclose(lc_file);
    }
  else
    {
      for (m=0; m < ntlcbins; m++)
	{
	  for (nn=0; nn < MALCBINS; nn++)
	    {
	      save[m][nn]=0.0;
	    }
	}
    }



  if ((lc_file = fopen("light_curve_res.out", "w+")) == NULL){
    printout("Cannot open lc_file.txt.\n");
    exit(0);
  }

  for (nn=0; nn < MALCBINS; nn++)
    {
      for (m=0; m < ntlcbins; m++)
	{
	  fprintf(lc_file, "%g %g\n", sqrt(light_curve_res[m][nn].lower_time*(light_curve_res[m][nn].lower_time + light_curve_res[m][nn].delta_t))/DAY, ((light_curve_res[m][nn].lum/LSUN) + save[m][nn]));
	}
    }

  fclose(lc_file);

  return(0);
}

/**********************************************************************/

/*Routine to add a packet to the outcoming light-curve.*/
/*See add_to_spec.*/

int
add_to_lc_res(pkt_ptr)
     PKT *pkt_ptr;
{
  double dot(), vec_len();
  int cross_prod();
  double t_arrive;
  int nt, na;
  int thetabin, phibin;
  double vec1[3], vec2[3], xhat[3], vec3[3];
  double costheta, cosphi, testphi;

  xhat[0]=1.0;
  xhat[1]=0;
  xhat[2]=0;

  /* Formula from Leon's paper. */


  t_arrive = pkt_ptr->escape_time - (dot(pkt_ptr->pos, pkt_ptr->dir)/CLIGHT_PROP);

  /* Put this into the time grid. */
  
  if (t_arrive > tmin && t_arrive < tmax)
    {
      nt = (log(t_arrive) - log(tmin)) / dlogtlc;

      /* for angle resolved case, need to work out the correct angle bin too. */

      
      costheta = dot(pkt_ptr->dir, syn_dir);
      thetabin = ((costheta + 1.0) * sqrt(MALCBINS) / 2.0);
      cross_prod(pkt_ptr->dir, syn_dir, vec1);
      cross_prod(xhat, syn_dir, vec2);
      cosphi = dot(vec1,vec2)/vec_len(vec1)/vec_len(vec2);

      cross_prod(vec2, syn_dir, vec3);      
      testphi = dot(vec1,vec3);

      if (testphi > 0)
	{
	  phibin = (acos(cosphi) /2. / PI * sqrt(MALCBINS));
	}
      else
	{
	  phibin = ((acos(cosphi) + PI) /2. / PI * sqrt(MALCBINS));
	}

      na = (thetabin*sqrt(MALCBINS)) + phibin;

      light_curve_res[nt][na].lum += pkt_ptr->e_rf / light_curve_res[nt][na].delta_t * MALCBINS / nprocs;
    }

  return(0);

}

/**********************************************************************/

int 
cross_prod (v1, v2, v3)
     double v1[3], v2[3], v3[3];
{
  v3[0] = (v1[1]*v2[2]) - (v2[1]*v1[2]);
  v3[1] = (v1[2]*v2[0]) - (v2[2]*v1[0]);
  v3[2] = (v1[0]*v2[1]) - (v2[0]*v1[1]);

  return(0);
}
