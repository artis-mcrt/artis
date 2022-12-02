#include "light_curve.h"

#include "exspec.h"
#include "sn3d.h"
#include "vectors.h"

// Routine to make a MC light curve from the r-packets.

void write_light_curve(const std::string &lc_filename, const int current_abin, const double *light_curve_lum,
                       const double *light_curve_lumcmf, const int numtimesteps) {
  FILE *lc_file = fopen_required(lc_filename.c_str(), "w");
  assert_always(numtimesteps <= globals::ntstep);

  printout("Writing %s\n", lc_filename.c_str());

  /// Print out the UVOIR bolometric light curve.
  for (int nts = 0; nts < numtimesteps; nts++) {
    fprintf(lc_file, "%g %g %g\n", globals::time_step[nts].mid / DAY, (light_curve_lum[nts] / LSUN),
            (light_curve_lumcmf[nts] / LSUN));
  }

  if (current_abin == -1) {
    /// Now print out the gamma ray deposition rate in the same file.
    for (int m = 0; m < numtimesteps; m++) {
      fprintf(lc_file, "%g %g %g\n", globals::time_step[m].mid / DAY,
              (globals::time_step[m].gamma_dep / LSUN / globals::time_step[m].width),
              (globals::time_step[m].cmf_lum / globals::time_step[m].width / LSUN));
    }
  }

  fclose(lc_file);
}

void add_to_lc_res(const struct packet *pkt_ptr, int current_abin, double *light_curve_lum, double *light_curve_lumcmf)
// add a packet to the outgoing light-curve.
{
  if (current_abin == -1) {
    /// Put this into the time grid
    const double arrive_time = get_arrive_time(pkt_ptr);
    if (arrive_time > globals::tmin && arrive_time < globals::tmax) {
      const int nt = get_timestep(arrive_time);
      safeadd(light_curve_lum[nt], pkt_ptr->e_rf / globals::time_step[nt].width / globals::nprocs);
    }

    /// Now do the cmf light curve.
    // t_arrive = pkt_ptr->escape_time * sqrt(1. - (vmax*vmax/CLIGHTSQUARED));
    const double arrive_time_cmf = get_arrive_time_cmf(pkt_ptr);
    if (arrive_time_cmf > globals::tmin && arrive_time_cmf < globals::tmax) {
      const int nt = get_timestep(arrive_time_cmf);
      safeadd(light_curve_lumcmf[nt], pkt_ptr->e_cmf / globals::time_step[nt].width / globals::nprocs /
                                          sqrt(1. - (globals::vmax * globals::vmax / CLIGHTSQUARED)));
    }

    return;
  } else if (get_escapedirectionbin(pkt_ptr->dir, globals::syn_dir) == current_abin) {
    // Add only packets which escape to the current angle bin
    double t_arrive = get_arrive_time(pkt_ptr);
    if (t_arrive > globals::tmin && t_arrive < globals::tmax) {
      int nt = get_timestep(t_arrive);
      safeadd(light_curve_lum[nt], pkt_ptr->e_rf / globals::time_step[nt].width * MABINS / globals::nprocs);
    }
  }
}