#include "sn3d.h"
#include "exspec.h"
#include "light_curve.h"
#include "vectors.h"


// Routine to make a MC light curve from the r-packets.

/// Light curve data structure
const int MALCBINS = 100;


void write_light_curve(
  char lc_filename[], int current_abin,
  const double *light_curve_lum,
  const double *light_curve_lumcmf)
{
  FILE *lc_file = fopen_required(lc_filename, "w");

  printout("Writing %s\n", lc_filename);

  /// Print out the UVOIR bolometric light curve.
  for (int nts = 0; nts < globals::ntstep; nts++)
  {
    fprintf(lc_file, "%g %g %g\n", globals::time_step[nts].mid / DAY,
            (light_curve_lum[nts] / LSUN), (light_curve_lumcmf[nts] / LSUN));
  }

  if (current_abin == -1)
  {
    /// Now print out the gamma ray deposition rate in the same file.
    for (int m = 0; m < globals::ntstep; m++)
    {
      fprintf(lc_file, "%g %g %g\n", globals::time_step[m].mid / DAY,
              (globals::time_step[m].gamma_dep / LSUN / globals::time_step[m].width),
              (globals::time_step[m].cmf_lum / globals::time_step[m].width/LSUN));
    }
  }

  fclose(lc_file);
}


void add_to_lc_res(const PKT *pkt_ptr, int current_abin, double *light_curve_lum, double *light_curve_lumcmf)
/**Routine to add a packet to the outcoming light-curve.*/
/**See add_to_spec.*/
{
  if (current_abin == -1)
  {
    /// Put this into the time grid
    const double arrive_time = get_arrive_time(pkt_ptr);
    if (arrive_time > globals::tmin && arrive_time < globals::tmax)
    {
      const int nt = get_timestep(arrive_time);
      safeadd(light_curve_lum[nt], pkt_ptr->e_rf / globals::time_step[nt].width / globals::nprocs);
    }

    /// Now do the cmf light curve.
    //t_arrive = pkt_ptr->escape_time * sqrt(1. - (vmax*vmax/CLIGHTSQUARED));
    const double arrive_time_cmf = get_arrive_time_cmf(pkt_ptr);
    if (arrive_time_cmf > globals::tmin && arrive_time_cmf < globals::tmax)
    {
      const int nt = get_timestep(arrive_time_cmf);
      safeadd(light_curve_lumcmf[nt], pkt_ptr->e_cmf / globals::time_step[nt].width / globals::nprocs / sqrt(1. - (globals::vmax * globals::vmax/CLIGHTSQUARED)));
    }

    return;
  }

  double vec1[3], vec2[3], xhat[3], vec3[3];

  xhat[0] = 1.0;
  xhat[1] = 0;
  xhat[2] = 0;

  /// Angle resolved case: need to work out the correct angle bin too. */
  double costheta = dot(pkt_ptr->dir, globals::syn_dir);
  int thetabin = ((costheta + 1.0) * sqrt(MALCBINS) / 2.0);
  cross_prod(pkt_ptr->dir, globals::syn_dir, vec1);
  cross_prod(xhat, globals::syn_dir, vec2);
  double cosphi = dot(vec1,vec2)/vec_len(vec1)/vec_len(vec2);

  cross_prod(vec2, globals::syn_dir, vec3);
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
    double t_arrive = get_arrive_time(pkt_ptr);
    if (t_arrive > globals::tmin && t_arrive < globals::tmax)
    {
      int nt = get_timestep(t_arrive);
      safeadd(light_curve_lum[nt], pkt_ptr->e_rf / globals::time_step[nt].width * MALCBINS / globals::nprocs);
    }
  }
}


void write_partial_lightcurve(int my_rank, int nts, PKT *pkts)
{
  const time_t time_func_start = time(NULL);

  double *rpkt_light_curve_lum = (double *) calloc(globals::ntstep, sizeof(double));
  double *rpkt_light_curve_lumcmf = (double *) calloc(globals::ntstep, sizeof(double));

  for (int ii = 0; ii < globals::npkts; ii++)
  {
    // printout("packet %d escape_type %d type %d", ii, pkt[ii].escape_type, pkt[ii].type);
    if (pkts[ii].type == TYPE_ESCAPE)
    {
      if (pkts[ii].escape_type == TYPE_RPKT)
      {
        add_to_lc_res(&pkts[ii], -1, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
      }
    }
  }

  #ifdef MPI_ON
  MPI_Allreduce(MPI_IN_PLACE, rpkt_light_curve_lum, globals::ntstep, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, rpkt_light_curve_lumcmf, globals::ntstep, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  #endif

  if (my_rank == 0)
  {
    write_light_curve((char *) "light_curve.out", -1, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
  }

  free(rpkt_light_curve_lum);
  free(rpkt_light_curve_lumcmf);

  printout("Saving partial light curve took %lds\n", time(NULL) - time_func_start);
}