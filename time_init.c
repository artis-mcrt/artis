#include "sn3d.h"
#include "time_init.h"

/* Subroutine to define the time steps.*/

int time_init(void)
{
  /// t=tmin is the start of the calcualtion. t=tmax is the end of the calculation.
  /// ntstep is the number of time steps wanted. For now the time steps
  /// are logarithmically spaced, but the code should be written such that
  /// they don't have to be.
  if (ntstep + 1 > MTSTEP)
  {
    printout("Error: too many timesteps. (%d) Abort.\n", ntstep);
    exit(0);
  }

  /// For logarithmic steps, the logarithmic inverval will be
  double dlogt = (log(tmax) - log(tmin))/ntstep;
  //dlogt = 0.17917595;
  //dlogt = 0.03583518938456109026;
  //dlogt = (tmax - tmin)/ntstep;

  /// Now setup the individual time steps
  for (int n = 0; n < ntstep; n++)
  {
    time_step[n].start = tmin * exp(n*dlogt);
    time_step[n].width = (tmin * exp((n+1)*dlogt)) - time_step[n].start;
    time_step[n].mid = tmin * exp((n+0.5)*dlogt);
    //time_step[n].start = tmin * n*dlogt;
    //time_step[n].width = (tmin * (n+1)*dlogt) - time_step[n].start;
    //time_step[n].mid = tmin * (n+0.5)*dlogt;

    //      printout("start %g, width %g, mid %g\n",time_step[n].start,time_step[n].width,time_step[n].mid);

    time_step[n].pellet_decays = 0;
    time_step[n].gamma_dep = 0.0;
    time_step[n].cmf_lum = 0.0;
  }

  /// and add a dummy timestep which contains the endtime
  /// of the calculation
  time_step[ntstep].start = tmax;
  time_step[ntstep].mid = tmax;

  return 0;
}
