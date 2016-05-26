#include "sn3d.h"
#include "exspec.h"
#include "gamma_light_curve.h"
#include "vectors.h"
#include "packet_init.h"

/* This is very like "light_curve" except that it now bins in angle. */

/* Routine to make a MC angle-dependent light curve for the gamma-packets. */

int make_gamma_light_curve(void)
{
  gather_gamma_light_curve(0);
  write_gamma_light_curve();

  return 0;
}


int gather_gamma_light_curve(int my_rank)
{
  PKT *pkt_ptr;
  int i,n,p,nn;
  char filename[100];               /// this must be long enough to hold "packetsxx.tmp" where xx is the number of "middle" iterations

  if (ntlcbins > MTLCBINS)
  {
    printout("Too many time bins in light curve - reducing.\n");
    ntlcbins = MTLCBINS;
  }

  /* start by setting up the time bins. */
  /* it is all done interms of a logarithmic spacing in t - get the
     step sizes first. */
  dlogtlc_angle = (log(tmax) - log(tmin))/ntlcbins;

  for (n = 0; n < ntlcbins; n++)
  {
    for (nn = 0; nn < MANGLCBINS; nn++)
    {
      light_curve_angle[n][nn].lower_time = exp( log(tmin) + (n * (dlogtlc_angle)));
      light_curve_angle[n][nn].delta_t = exp( log(tmin) + ((n+1) * (dlogtlc_angle))) - light_curve_angle[n][nn].lower_time;
      light_curve_angle[n][nn].lum = 0.0;
    }
  }

  /// The grid is now set up. Now we loop over all the packets, check if they made it out or not,
  /// and if they did we add their rest frame energy to the appropriate cell.
  int middle_iteration;
  FILE *packets_file;
  for (middle_iteration = 0; middle_iteration < n_middle_it; middle_iteration++)
  {
    for (i = 0; i < nprocs; i++)
    {
      /// Read in the next bunch of packets to work on
      //sprintf(filename,"packets%d_%d.tmp",0,i);
      sprintf(filename,"packets%.2d_%.4d.out",0,i);
      //if ((packets_file = fopen(filename, "rb")) == NULL)
      if ((packets_file = fopen(filename, "r")) == NULL)
      {
        printf("Cannot open packets file\n");
        abort();
      }
      //fread(&pkt[0], sizeof(PKT), npkts, packets_file);
      read_packets(packets_file);

      /// Close the current file.
      fclose(packets_file);

      /// And figure out the escaping packets of those.
      for (p = 0; p < npkts; p++)
      {
        pkt_ptr = &pkt[p];
        if (pkt_ptr->type == TYPE_ESCAPE && pkt_ptr->escape_type == TYPE_GAMMA)
        {
          /// It made it out.
          add_to_lc_angle(pkt_ptr);
        }
      }
    }
  }

  return 0;
}


int write_gamma_light_curve(void)
{
  double save[MTLCBINS][MANGLCBINS];

  /* Light curve is done - write it out. */

  /* If needed, start by reading in existing file and storing old numbers. */

  FILE *lc_gamma_file;

  if (file_set)
  {
    if ((lc_gamma_file = fopen("gamma_light_curve.out", "r")) == NULL)
    {
      printout("Cannot open lc_gamma_file.txt.\n");
      exit(0);
    }
    for (int m = 0; m < ntlcbins; m++)
    {
      float dum1, dum2;
      fscanf(lc_gamma_file, "%g", &dum1);
      for (int nn = 0; nn < MANGLCBINS; nn++)
      {
        fscanf(lc_gamma_file, " %g ",&dum2);
        save[m][nn]=dum2;
      }
    }
    fclose(lc_gamma_file);
  }
  else
  {
    for (int m = 0; m < ntlcbins; m++)
    {
      for (int nn = 0; nn < MANGLCBINS; nn++)
      {
        save[m][nn]=0.0;
      }
    }
  }


  if ((lc_gamma_file = fopen("gamma_light_curve.out", "w+")) == NULL)
  {
    printout("Cannot open lc_gamma_file.txt.\n");
    exit(0);
  }

  for (int m = 0; m < ntlcbins; m++)
  {
    fprintf(lc_gamma_file, "%g ", sqrt(light_curve_angle[m][0].lower_time*(light_curve_angle[m][0].lower_time + light_curve_angle[m][0].delta_t))/DAY);
    for (int nn = 0; nn < MANGLCBINS; nn++)
    {
      fprintf(lc_gamma_file, " %g ",(light_curve_angle[m][nn].lum/LSUN) + save[m][nn]);
    }
    fprintf(lc_gamma_file, "\n");
  }


  fclose(lc_gamma_file);

  return 0;
}


int add_to_lc_angle(PKT *pkt_ptr)
/*Routine to add a packet to the outcoming light-curve.*/
/*See add_to_spec.*/
{
  int nt, na;
  int thetabin, phibin;
  double vec1[3], vec2[3], xhat[3], vec3[3];
  double costheta, cosphi, testphi;

  xhat[0] = 1.0;
  xhat[1] = 0;
  xhat[2] = 0;

  /* Formula from Leon's paper. */

  double t_arrive = pkt_ptr->escape_time - (dot(pkt_ptr->pos, pkt_ptr->dir)/CLIGHT_PROP);

  /* Put this into the time grid. */

  if (t_arrive > tmin && t_arrive < tmax)
  {
    nt = (log(t_arrive) - log(tmin)) / dlogtlc_angle;
    /* Difference from light_curve is here - need to assign an angle bin. Compute the angle between
 the viewing direcation and the packet's trajectory. */

    costheta = dot(pkt_ptr->dir, syn_dir);
    thetabin = ((costheta + 1.0) * sqrt(MANGLCBINS) / 2.0);
    cross_prod(pkt_ptr->dir, syn_dir, vec1);
    cross_prod(xhat, syn_dir, vec2);
    cosphi = dot(vec1,vec2)/vec_len(vec1)/vec_len(vec2);

    cross_prod(vec2, syn_dir, vec3);
    testphi = dot(vec1,vec3);

    if (testphi > 0)
    {
      phibin = (acos(cosphi) /2. / PI * sqrt(MANGLCBINS));
    }
    else
    {
      phibin = ((acos(cosphi) + PI) /2. / PI * sqrt(MANGLCBINS));
    }

    na = (thetabin * sqrt(MANGLCBINS)) + phibin;

    light_curve_angle[nt][na].lum += pkt_ptr->e_rf / light_curve_angle[nt][na].delta_t * MANGLCBINS / nprocs;
  }

  return 0;
}
