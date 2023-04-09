#include "light_curve.h"

#include "exspec.h"
#include "sn3d.h"
#include "vectors.h"

// Routine to make a MC light curve from the r-packets.

void write_light_curve(const std::string &lc_filename, const int current_abin, const double *light_curve_lum,
                       const double *light_curve_lumcmf, const int numtimesteps) {
  assert_always(numtimesteps <= globals::ntstep);

  std::ofstream lc_file(lc_filename);
  if (!lc_file) {
    throw std::system_error(errno, std::system_category(), "failed to open " + lc_filename);
  }

  printout("Writing %s\n", lc_filename.c_str());

  constexpr int maxlen = 1024;
  char linebuffer[maxlen];

  /// Print out the UVOIR bolometric light curve.
  for (int nts = 0; nts < numtimesteps; nts++) {
    assert_always(snprintf(linebuffer, maxlen, "%g %g %g", globals::time_step[nts].mid / DAY,
                           (light_curve_lum[nts] / LSUN), (light_curve_lumcmf[nts] / LSUN)) < maxlen);
    lc_file << linebuffer << '\n';
  }

  if (current_abin == -1) {
    /// Now print out the gamma ray deposition rate in the same file.
    for (int m = 0; m < numtimesteps; m++) {
      assert_always(snprintf(linebuffer, maxlen, "%g %g %g", globals::time_step[m].mid / DAY,
                             (globals::time_step[m].gamma_dep / LSUN / globals::time_step[m].width),
                             (globals::time_step[m].cmf_lum / globals::time_step[m].width / LSUN)) < maxlen);
      lc_file << linebuffer << '\n';
    }
  }
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
    double const t_arrive = get_arrive_time(pkt_ptr);
    if (t_arrive > globals::tmin && t_arrive < globals::tmax) {
      int const nt = get_timestep(t_arrive);
      safeadd(light_curve_lum[nt], pkt_ptr->e_rf / globals::time_step[nt].width * MABINS / globals::nprocs);
    }
  }
}