#include "sn3d.h"

/* Subroutine to move the rays for syn calculation.  */

int update_gamma_rays(nts)
     int nts;  //the time step we're doing
{
  RAY *ray_ptr;
  int ray_prop();

  /* At the start, the rays are either still at their starting points or have already been
     processed for one or more timesteps. */

  double ts = time_step[nts].start;
  double tw = time_step[nts].width;

  for (int n = 0; n < NRAYS_SYN; n++)
  {
    ray_ptr = &rays[n];
    if (ray_ptr->status == ACTIVE)
    {
      /* It should already be going - carry on! */
      ray_prop(ray_ptr, ts, ts+tw, nts);
    }
    else if (ray_ptr->tstart < ts && ray_ptr->status == WAITING)
    {
      /* It should already be going but isn't for some reason. Abort. */
      printout("Ray should be active but isn't. Abort.\n");
      exit(0);
    }
    else if (ray_ptr->tstart < (ts + tw) && ray_ptr->status == WAITING)
    {
      /* It should kick off during this step. */
      ray_ptr->status = ACTIVE;
      ray_prop(ray_ptr, ray_ptr->tstart, ts+tw, nts);
    }
    /* Note: if tstart > ts+tw it'll sit this one out - still WAITING. If it's already FINISHED
   then nothing should happen. */
  }

  return 0;
}
