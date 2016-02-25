#include "sn3d.h"

/* Master subroutine for computing the formal integral grey optical light curve. */

double syn_lc()
   /* Note the direction (vector pointing from origin to observer is
  externally set: syn_dir[3]. */
{
  int syn_lc_init();
  int time_init();
  int update_gamma_rays();
  int update_grid();
  FILE *syn_file;

  if ((syn_file = fopen("syn_lc.out", "w+")) == NULL)
  {
    printout("Cannot open syn.out.\n");
    exit(0);
  }

  /* main loop over the set of times for which we want spectra. */

  for (int time_loop = 0; time_loop < nsyn_time; time_loop++)
  {
    grid_init();

    syn_lc_init(time_syn[time_loop]);

    /* Now start off the calculation as in sn3d: (re)initialise the grid. */

    grid_init();

    time_init();

    for (int nts = 0; nts < ntstep; nts++)
    {

      if (do_rlc_est != 0)
      {
        double lower_test = time_syn[time_loop] / (1. + (sqrt(3.0)*vmax/CLIGHT));
        double upper_test = time_syn[time_loop] / (1. - (sqrt(3.0)*vmax/CLIGHT));
        if (time_step[nts].start > lower_test)
        {
          if ((time_step[nts].start + time_step[nts].width) < upper_test)
          {
            emiss_rlc_load(nts);
          }
        }
      }
      printout("syn time step %d\n", nts);
      update_grid(nts);
      update_gamma_rays(nts);  //because of choice of frequency, we can
                               //use the same ray tracing routines that were made for gamma rays
    }

    /* Now make the rays into a spectrum and print it out. */

    double lc_syn = 0;

    for (int m = 0; m < NRAYS_SYN; m++)
    {
      lc_syn += rays[m].e_rf[0] * DeltaA / 1.e12 / PARSEC/ PARSEC;
      if (fabs(1. - rays[m].nu_rf_init[0]/rays[m].nu_rf[0]) > 1.e-8)
      {
        printout("changes in frequencies in rest frame. Abort.\n");
        exit(0);
      }
    }
    fprintf(syn_file, "%g %g\n", time_syn[time_loop]/DAY, lc_syn*4.*PI*PARSEC*PARSEC*1.e12/LSUN);

  }

  fclose(syn_file);

  return 0;
}

/* ************************************************************** */
int syn_lc_init(time)
     double time;
{
  double max_proj;
  double r0[3];
  int get_velocity();
  double doppler();
  double vel_vec[3];
  int get_cell();
  double vec_len(), dot();
  double kick[3];

  int num_in = 0;

  /* Start by finding a bunch of rays that will reach observer with correct
     freq at correct time. */

  int n = 0;
  //for (n = 0; n < NRAYS_SYN; n++)
  while (n < NRAYS_SYN)
  {
    /* To correctly sample the rays, they should be uniformly spread over
   the observer's plane (i.e the plane perpendocular to syn_dir. */
    /* For a cuboidal grid, the maximum length from the centre in
   projection must be */

    max_proj = sqrt( (xmax*xmax) + (ymax*ymax) + (zmax*zmax)) * (time/tmin);
    max_proj = max_proj / sqrt( 1. - (vmax*vmax/CLIGHT_PROP/CLIGHT_PROP));

    /* so sample a circle (in plane of origin) with radius max_proj for
       uniform area. */

    /* The line followed by the ray is going to be
       r = r0 + lambda syn_dir
       where r, and r0 are vectors */

    double zrand = gsl_rng_uniform(rng);
    double r_proj = (sqrt( zrand)) * max_proj;
    zrand = gsl_rng_uniform(rng);
    double phi_proj = 4 * asin(1.0) * zrand;

    /* Get cartestian coordinates in a frame with zprime parallel to
       syn_dir. */

    double zprime = 0.0;
    double yprime = r_proj * cos(phi_proj);
    double xprime = r_proj * sin(phi_proj);

    /* Now rotate back to real xyz - based on syn_dir - to get r0 */

    if (syn_dir[2] > 0.9999999)
    {
      /* no need to rotate */
      r0[0] = xprime;
      r0[1] = yprime;
      r0[2] = 0;
    }
    else if (syn_dir[2] < -0.9999999)
    {
      /* no need to rotate */
      r0[0] = xprime;
      r0[1] = yprime;
      r0[2] = 0;
    }
    else
    {
      double norm1 = 1./(sqrt( (syn_dir[0]*syn_dir[0]) + (syn_dir[1]*syn_dir[1])));
      double norm2 = 1./(sqrt( (syn_dir[0]*syn_dir[0]) + (syn_dir[1]*syn_dir[1]) + (syn_dir[2]*syn_dir[2])));

      double r11 = syn_dir[1] * norm1;
      double r12 = -1 * syn_dir[0] * norm1;
      double r13 = 0.0;
      double r21 = syn_dir[0] * syn_dir[2] * norm1 * norm2;
      double r22 = syn_dir[1] * syn_dir[2] * norm1 * norm2;
      double r23 = -1 * norm2 / norm1;
      double r31 = syn_dir[0] * norm2;
      double r32 = syn_dir[1] * norm2;
      double r33 = syn_dir[2] * norm2;

      r0[0] = (r11 * xprime) + (r21 * yprime) + (r31 * zprime);
      r0[1] = (r12 * xprime) + (r22 * yprime) + (r32 * zprime);
      r0[2] = (r13 * xprime) + (r23 * yprime) + (r33 * zprime);
    }

    /***************************************************/

    //  printout("ratio %g\n",vec_len(r0)/r_proj);
    //      printout("dot %g\n", dot(r0, syn_dir)/vec_len(r0));


    /* Now, we want the observer to detect this ray at t = time so we need
       lambda = -CLIGHT * time to give us where the ray was at t = 0
       lambda = -CLIGHT * (time - tmin) to give us where it was at t= tmin */

    /* Need to answer question: does this ray ever enter the grid - if so, where and when. */

    double tin_x = ((CLIGHT_PROP * time * syn_dir[0]) - r0[0]) / ((CLIGHT_PROP * syn_dir[0]) + (xmax/tmin));
    double tou_x = ((CLIGHT_PROP * time * syn_dir[0]) - r0[0]) / ((CLIGHT_PROP * syn_dir[0]) - (xmax/tmin));

    if (syn_dir[0] < 0) //swap around if ray travelling to -ve x
    {
      double dummy = tou_x;
      tou_x = tin_x;
      tin_x = dummy;
    }
    if (tou_x < 0 && tin_x <0)
    {
      tou_x = 0.0;
      tin_x = 1.e99;
    }
    else if (tou_x < 0)
    {
      tin_x = tin_x;
      tou_x = 1.e99;
    }
    else if (tin_x < 0)
    {
      //printout("error in syn (sigh)\n");
      //printout("x in %g out %g\n",tin_x, tou_x);
      // exit(0);
      tin_x = tou_x;
      tou_x = 1.e99;
    }

    //      printout("tin_x %g tou_x %g\n", tin_x , tou_x);

    double tin_y = ((CLIGHT_PROP * time * syn_dir[1]) - r0[1]) / ((CLIGHT_PROP * syn_dir[1]) + (ymax/tmin));
    double tou_y = ((CLIGHT_PROP * time * syn_dir[1]) - r0[1]) / ((CLIGHT_PROP * syn_dir[1]) - (ymax/tmin));
    if (syn_dir[1] < 0)
    {
      double dummy = tou_y;
      tou_y = tin_y;
      tin_y = dummy;
    }
    if (tou_y < 0 && tin_y <0)
    {
      tou_y = 0.0;
      tin_y = 1.e99;
    }
    else if (tou_y < 0)
    {
      tin_y = tin_y;
      tou_y = 1.e99;
    }
    else if (tin_y < 0)
    {
      //printout("error in syn (sigh)\n");
      //printout("y in %g out %g\n",tin_y, tou_y);
      //exit(0);
      tin_y = tou_y;
      tou_y = 1.e99;
      //printout("tin_y %g tou_y %g\n", tin_y , tou_y);
    }

    double tin_z = ((CLIGHT_PROP * time * syn_dir[2]) - r0[2]) / ((CLIGHT_PROP * syn_dir[2]) + (zmax/tmin));
    double tou_z = ((CLIGHT_PROP * time * syn_dir[2]) - r0[2]) / ((CLIGHT_PROP * syn_dir[2]) - (zmax/tmin));
      //printout("tin_z %g tou_z %g r0[2] %g\n", tin_z , tou_z, r0[2]);
    if (syn_dir[2] < 0)
    {
      double dummy = tou_z;
      tou_z = tin_z;
      tin_z = dummy;
    }
    if (tou_z < 0 && tin_z <0)
    {
      tou_z = 0.0;
      tin_z = 1.e99;
    }
    else if (tou_z < 0)
    {
      tin_z = tin_z;
      tou_z = 1.e99;
    }
    else if (tin_z < 0)
    {
      //printout("error in syn (sigh)\n");
      //printout("z in %g out %g\n",tin_z, tou_z);
      //exit(0);
      tin_z = tou_z;
      tou_z = 1.e99;
    }
    //printout("tin_z %g tou_z %g\n", tin_z , tou_z);

    /* Find latest tin */
    double tin, tou;
    if (tin_x > tin_y)
    {
      if (tin_x > tin_z)
      {
        tin = tin_x;
      }
      else
      {
        tin = tin_z;
      }
    }
    else
    {
      if (tin_y > tin_z)
      {
        tin = tin_y;
      }
      else
      {
        tin = tin_z;
      }
    }

    /* Find 1st tou */
    if (tou_x < tou_y)
    {
      if (tou_x < tou_z)
      {
        tou = tou_x;
      }
      else
      {
        tou = tou_z;
      }
    }
    else
    {
      if (tou_y < tou_z)
      {
        tou = tou_y;
      }
      else
      {
        tou = tou_z;
      }
    }

    /* did it get into the grid then? */

    num_in += 1;

    if (tin < 0)
    {
      printout("negative tin (syn_lc)???\n");
      exit(0);
    }

    if (tin < tou)
    {
      /* It did enter the grid - record it. */
      // if (tin > time)
      //   {
      //    printout("hmmmmmm....\n");
      //    exit(0);
      //  }
      double tin_kick = tin * (1. + 1.e-10);
      rays[n].tstart = tin_kick;
      rays[n].rstart[0] = r0[0] - (CLIGHT_PROP*(time - tin)*syn_dir[0]);
      rays[n].rstart[1] = r0[1] - (CLIGHT_PROP*(time - tin)*syn_dir[1]);
      rays[n].rstart[2] = r0[2] - (CLIGHT_PROP*(time - tin)*syn_dir[2]);
      kick[0] = (vec_len(rays[n].rstart)*1.e-10) * syn_dir[0];
      kick[1] = (vec_len(rays[n].rstart)*1.e-10) * syn_dir[1];
      kick[2] = (vec_len(rays[n].rstart)*1.e-10) * syn_dir[2];
      rays[n].pos[0] = rays[n].rstart[0]+kick[0];
      rays[n].pos[1] = rays[n].rstart[1]+kick[1];
      rays[n].pos[2] = rays[n].rstart[2]+kick[2];
      rays[n].where = get_cell(rays[n].pos[0],rays[n].pos[1],rays[n].pos[2],tin_kick);
      /* kick is there to boof it in to the cell we want (avoid boundary). */
      rays[n].status = WAITING;

      get_velocity(rays[n].rstart, vel_vec, tin);
      if (NSYN != 1)
      {
        printout("Doing a grey calculation with coloured rays. Abort.\n");
        exit(0);
      }

      for (int i = 0; i < NSYN; i++)
      {
        rays[n].nu_rf[i] = 1.;  //grey rays - just a place holder;
        rays[n].nu_rf_init[i] = rays[n].nu_rf[i];
        rays[n].nu_cmf[i] = rays[n].nu_rf[i] * doppler(syn_dir, vel_vec);
        rays[n].e_rf[i] = 0.0;
        rays[n].e_cmf[i] = 0.0;
      }
      n = n + 1;
    }
  }

  printout("number of rays that went in %d\n", num_in);

  DeltaA = PI * max_proj * max_proj / num_in;
  printout("DeltaA %g\n", DeltaA);

  return 0;
}
