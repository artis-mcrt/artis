#include "sn3d.h"
#include "exspec.h"
#include "light_curve.h"
#include "vectors.h"


// Routine to make a MC light curve from the r-packets.

/*int make_light_curve()
{
  gather_light_curve();
  write_light_curve();

  return 0;
}
*/

struct lc
{
  float lower_time;
  float delta_t;
  double lum;
};

/// Light curve data structure
#define MTLCBINS  200
#define MALCBINS  100
#define MANGLCBINS 100


static struct lc light_curve[MTLCBINS];
static struct lc light_curve_cmf[MTLCBINS];
// static struct lc light_curve_angle[MTLCBINS][MANGLCBINS];


static void add_to_lc(const EPKT *pkt_ptr)
/**Routine to add a packet to the outcoming light-curve.*/
/**See add_to_spec.*/
{
  /// Put this into the time grid
  double t_arrive = pkt_ptr->arrive_time;
  if (t_arrive > tmin && t_arrive < tmax)
  {
    int nt = (log(t_arrive) - log(tmin)) / dlogtlc;
    light_curve[nt].lum += pkt_ptr->e_rf / light_curve[nt].delta_t / nprocs;
  }

  /// Now do the cmf light curve.
  //t_arrive = pkt_ptr->escape_time * sqrt(1. - (vmax*vmax/CLIGHTSQUARED));
  t_arrive = pkt_ptr->arrive_time_cmf;
  if (t_arrive > tmin && t_arrive < tmax)
  {
    int nt = (log(t_arrive) - log(tmin)) / dlogtlc;
    light_curve_cmf[nt].lum += pkt_ptr->e_cmf / light_curve[nt].delta_t / nprocs / sqrt(1. - (vmax*vmax/CLIGHTSQUARED));
  }
}


void init_light_curve(void)
{
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

  for (int n = 0; n < ntlcbins; n++)
  {
    light_curve[n].lower_time = exp( log(tmin) + (n * (dlogtlc)));
    light_curve[n].delta_t = exp( log(tmin) + ((n+1) * (dlogtlc))) - light_curve[n].lower_time;
    light_curve[n].lum = 0.0;
    light_curve_cmf[n].lum = 0.0;
  }
}


void write_light_curve(char lc_filename[], int current_abin)
{
  FILE *lc_file = fopen_required(lc_filename, "w");

  printout("Writing %s\n", lc_filename);

  /*
  FILE *lc_file;
  double save2[MTSTEP][2];
  double save[MTLCBINS][2];
  float dum1, dum2, dum3;

  /// Light curve is done - write it out.
  /// If needed, start by reading in existing file and storing old numbers.
  if (file_set)
  {
    if ((lc_file = fopen("light_curve.out", "r")) == NULL)
    {
      printout("Cannot open lc_file.txt.\n");
      abort();
    }
    for (int m = 0; m < ntlcbins; m++)
    {
        fscanf(lc_file, "%g %g %g\n", &dum1, &dum2, &dum3);
        save[m][0]=dum2;
        save[m][1]=dum3;
  	}

    for (int m = 0; m < ntstep; m++)
    {
      fscanf(lc_file, "%g %g %g\n", &dum1, &dum2, &dum3);
      save2[m][0]=dum2;
      save2[m][1]=dum3;
    }
    fclose(lc_file);
  }
  else
  {
    for (int m = 0; m < ntlcbins; m++)
    {
      save[m][0]=0.0;
      save[m][1]=0.0;
    }

    for (int m = 0; m < ntstep; m++)
    {
      save2[m][0]=0.0;
      save2[m][1]=0.0;
    }
  }


  if ((lc_file = fopen("light_curve.out", "w+")) == NULL){
    printout("Cannot open lc_file.txt.\n");
    abort();
  }

  for (int m = 0; m < ntlcbins; m++)
  {
    fprintf(lc_file, "%g %g %g\n", sqrt(light_curve[m].lower_time*(light_curve[m].lower_time + light_curve[m].delta_t))/DAY, ((light_curve[m].lum/LSUN) + save[m][0]), ((light_curve_cmf[m].lum/LSUN) + save[m][1]));
  }

  /// Now print out the gamma ray deposition rate in the same file.

  for (int m = 0; m < ntstep; m++)
  {
    fprintf(lc_file, "%g %g %g\n", time_step[m].mid/DAY, save2[m][0] + (time_step[m].gamma_dep/LSUN/time_step[m].width), save2[m][1] + (time_step[m].cmf_lum/time_step[m].width/LSUN));
  }

  fclose(lc_file);
  */

  /// Print out the UVOIR bolometric light curve.
  for (int m = 0; m < ntlcbins; m++)
  {
    fprintf(lc_file, "%g %g %g\n", sqrt(light_curve[m].lower_time*(light_curve[m].lower_time + light_curve[m].delta_t))/DAY, (light_curve[m].lum/LSUN), (light_curve_cmf[m].lum/LSUN));
  }

  if (current_abin == -1)
  {
    /// Now print out the gamma ray deposition rate in the same file.
    for (int m = 0; m < ntstep; m++)
    {
      fprintf(lc_file, "%g %g %g\n", time_step[m].mid/DAY, (time_step[m].gamma_dep/LSUN/time_step[m].width),  (time_step[m].cmf_lum/time_step[m].width/LSUN));
    }
  }

  fclose(lc_file);
}


void gather_light_curve(EPKT *epkts, int nepkts)
{
  //void read_packets(FILE *packets_file);
  //int i,n,p;

  /// Now add the energy of all the escaping packets to the
  /// appropriate bins.
  for (int p = 0; p < nepkts; p++)
  {
    add_to_lc(&epkts[p]);
  }
}


static void add_to_lc_res(const EPKT *pkt_ptr, int current_abin)
/**Routine to add a packet to the outcoming light-curve.*/
/**See add_to_spec.*/
{
  double vec1[3], vec2[3], xhat[3], vec3[3];

  xhat[0] = 1.0;
  xhat[1] = 0;
  xhat[2] = 0;

  /// Angle resolved case: need to work out the correct angle bin too. */
  double costheta = dot(pkt_ptr->dir, syn_dir);
  int thetabin = ((costheta + 1.0) * sqrt(MALCBINS) / 2.0);
  cross_prod(pkt_ptr->dir, syn_dir, vec1);
  cross_prod(xhat, syn_dir, vec2);
  double cosphi = dot(vec1,vec2)/vec_len(vec1)/vec_len(vec2);

  cross_prod(vec2, syn_dir, vec3);
  double testphi = dot(vec1,vec3);

  int phibin;
  if (testphi > 0)
  {
    phibin = (acos(cosphi) /2. / PI * sqrt(MALCBINS));
  }
  else
  {
    phibin = ((acos(cosphi) + PI) /2. / PI * sqrt(MALCBINS));
  }
  int na = (thetabin*sqrt(MALCBINS)) + phibin;

  /// Add only packets which escape to the current angle bin
  if (na == current_abin)
  {
    /// Put this into the time grid.
    double t_arrive = pkt_ptr->arrive_time;
    if (t_arrive > tmin && t_arrive < tmax)
    {
      int nt = (log(t_arrive) - log(tmin)) / dlogtlc;
      light_curve[nt].lum += pkt_ptr->e_rf / light_curve[nt].delta_t * MALCBINS / nprocs;
    }
  }
}


void gather_light_curve_res(EPKT *epkts, int nepkts, int current_abin)
{
  //void read_packets(FILE *packets_file);
  //int i,n,p,nn;

  /// Now add the energy of all the escaping packets to the
  /// appropriate bins.
  for (int p = 0; p < nepkts; p++)
  {
    add_to_lc_res(&epkts[p], current_abin);
  }
}
