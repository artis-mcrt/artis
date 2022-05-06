// #include <gsl/gsl_poly.h>
#ifndef __CUDA_ARCH__
#include <gsl/gsl_blas.h>
#endif

#include "sn3d.h"
#include "boundary.h"
#include "grid.h"
#include "rpkt.h"
#include "stats.h"
#include "update_packets.h"
#include "vectors.h"


__host__ __device__
static double get_shellcrossdist(
  const double pos[3], const double dir[3], const double shellradius, const bool isinnerboundary, const double tstart)
// find the closest forward distance to the intersection of a ray with an expanding spherical shell
// return -1 if there are no forward intersections (or if the intersection is tangential to the shell)
{
  assert_always(shellradius > 0);
  const bool debug = false;
  if (debug)
  {
    printout("get_shellcrossdist isinnerboundary %d\n", isinnerboundary);
    printout("shellradius %g tstart %g len(pos) %g\n", shellradius, tstart, vec_len(pos));
  }
  const double speed = vec_len(dir) * globals::CLIGHT_PROP;
  const double a = dot(dir, dir) - pow(shellradius / tstart / speed, 2);
  const double b = 2 * (dot(dir, pos) - pow(shellradius, 2) / tstart / speed);
  const double c = dot(pos, pos) - pow(shellradius, 2);

  const double discriminant = pow(b, 2) - 4 * a * c;

  if (discriminant < 0)
  {
    // no intersection
    assert_always(shellradius < vec_len(pos));
    if (debug)
      printout("no intersection\n");
    return -1;
  }
  else if (discriminant > 0)
  {
    // two intersections
    double d1 = (-b + sqrt(discriminant)) / 2 / a;
    double d2 = (-b - sqrt(discriminant)) / 2 / a;

    double posfinal1[3];
    double posfinal2[3];
    #ifndef __CUDA_ARCH__
    cblas_dcopy(3, pos, 1, posfinal1, 1);     // posfinal1 = pos
    cblas_daxpy(3, d1, dir, 1, posfinal1, 1); // posfinal1 += d1 * dir;

    cblas_dcopy(3, pos, 1, posfinal2, 1);
    cblas_daxpy(3, d2, dir, 1, posfinal2, 1);
    #else
    for (int d = 0; d < 3; d++)
    {
      posfinal1[d] = pos[d] + d1 * dir[d];
      posfinal2[d] = pos[d] + d2 * dir[d];
    }
    #endif

    const double shellradiusfinal1 = shellradius / tstart * (tstart + d1 / speed);
    const double shellradiusfinal2 = shellradius / tstart * (tstart + d2 / speed);
    // printout("solution1 d1 %g radiusfinal1 %g shellradiusfinal1 %g\n", d1, vec_len(posfinal1), shellradiusfinal1);
    // printout("solution2 d2 %g radiusfinal2 %g shellradiusfinal2 %g\n", d2, vec_len(posfinal2), shellradiusfinal2);
    assert_always(fabs(vec_len(posfinal1) / shellradiusfinal1 - 1.) < 1e-3);
    assert_always(fabs(vec_len(posfinal2) / shellradiusfinal2 - 1.) < 1e-3);

    // invalidate any solutions that require entering the boundary from the wrong radial direction
    if (isinnerboundary)
    {
      if (dot(posfinal1, dir) > 0.)
      {
        d1 = -1;
      }
      if (dot(posfinal2, dir) > 0.)
      {
        d2 = -1;
      }
    }
    else
    {
      if (dot(posfinal1, dir) < 0.)
      {
        d1 = -1;
      }
      if (dot(posfinal2, dir) < 0.)
      {
        d2 = -1;
      }
    }

    // negative d means in the reverse direction along the ray
    // ignore negative d values, and if two are positive then return the smaller one
    if (d1 < 0 && d2 < 0)
    {
      return -1;
    }
    else if (d2 < 0)
    {
      return d1;
    }
    else if (d1 < 0)
    {
      return d2;
    }
    else
    {
      return fmin(d1, d2);
    }
  }
  else
  {
    // exactly one intersection
    // ignore this and don't change which cell the packet is in
    assert_always(shellradius <= vec_len(pos));
    printout("single intersection\n");
    return -1.;
  }
}


__host__ __device__
double boundary_cross(PKT *const pkt_ptr, const double tstart, int *snext)
/// Basic routine to compute distance to a cell boundary.
{
  assert_always(tstart == pkt_ptr->prop_time);

  // There are six possible boundary crossings. Each of the three
  // cartesian coordinates may be taken in turn. For x, the packet
  // trajectory is
  // x = x0 + (dir.x) * c * (t - tstart)
  // the boundries follow
  // x+/- = x+/-(tmin) * (t/tmin)
  // so the crossing occurs when
  // t = (x0 - (dir.x)*c*tstart)/(x+/-(tmin)/tmin - (dir.x)c)

  // Modified so that it also returns the distance to the closest cell
  // boundary, regardless of direction.

  // d is used to loop over the coordinate indicies 0,1,2 for x,y,z

  const int cellindex = pkt_ptr->where;

  // the following four vectors are in grid coordinates, so either x,y,z or r
  const int ndim = grid::get_ngriddimensions();
  assert_testmodeonly(ndim <= 3);
  double initpos[3];       // pkt_ptr->pos converted to grid coordinates
  double cellcoordmax[3];
  double vel[3];           // pkt_ptr->dir * globals::CLIGHT_PROP converted to grid coordinates

  if (grid::grid_type == GRID_UNIFORM)
  {
    // XYZ coordinates
    for (int d = 0; d < ndim; d++)
    {
      initpos[d] = pkt_ptr->pos[d];
      cellcoordmax[d] = grid::get_cellcoordmin(cellindex, d) + grid::wid_init(0);
      vel[d] = pkt_ptr->dir[d] * globals::CLIGHT_PROP;
    }
  }
  else if (grid::grid_type == GRID_SPHERICAL1D)
  {
    // the only coordinate is radius from the origin
    initpos[0] = vec_len(pkt_ptr->pos);
    cellcoordmax[0] = grid::get_cellcoordmin(cellindex, 0) + grid::wid_init(cellindex);
    vel[0] = dot(pkt_ptr->pos, pkt_ptr->dir) / vec_len(pkt_ptr->pos) * globals::CLIGHT_PROP; // radial velocity
  }

  // for (int d = 0; d < ndim; d++)
  // {
  //   if (initpos[d] < grid::get_cellcoordmin(cellindex, d) || initpos[d] > cellcoordmax[d])
  //   {
  //     printout("WARNING: packet should have already escaped.\n");
  //     *snext = -99;
  //     return 0;
  //   }
  // }

  //printout("boundary.c: x0 %g, y0 %g, z0 %g\n", initpos[0] initpos[1] initpos[2]);

  //printout("boundary.c: vx %g, vy %g, vz %g\n",vel[0],vel[1],vel[2]);

  //printout("boundary.c: cellxmin %g, cellymin %g, cellzmin %g\n",grid::get_cellcoordmin(cellindex, 0),grid::get_cellcoordmin(cellindex, 1),grid::get_cellcoordmin(cellindex, 2));

  //printout("boundary.c: cellxmax %g, cellymax %g, cellzmax %g\n",cellcoordmax[0],cellcoordmax[1],cellcoordmax[2]);

  // enum cell_boundary not_allowed = pkt_ptr->last_cross;
  enum cell_boundary not_allowed = pkt_ptr->last_cross;
  enum cell_boundary negdirections[3] = {NEG_X, NEG_Y, NEG_Z}; // 'X' might actually be radial coordinate
  enum cell_boundary posdirections[3] = {POS_X, POS_Y, POS_Z};

  // printout("checking inside cell boundary\n");
  for (int d = 0; d < ndim; d++)
  {
    // flip is either zero or one to indicate +ve and -ve boundaries along the selected axis
    for (int flip = 0; flip < 2; flip++)
    {
      enum cell_boundary direction = flip ? posdirections[d] : negdirections[d];
      enum cell_boundary invdirection = !flip ? posdirections[d] : negdirections[d];
      const int cellindexdiff = flip ? - grid::get_coordcellindexincrement(d) : grid::get_coordcellindexincrement(d);

      const double tolerance = 1e-6 * globals::vmax * tstart;  // maximum position of out bound tolerance
      bool isoutside;
      if (flip)
      {
        isoutside = initpos[d] < (grid::get_cellcoordmin(cellindex, d) / globals::tmin * tstart - tolerance);
      }
      else
      {
        isoutside = initpos[d] > (cellcoordmax[d] / globals::tmin * tstart + tolerance);
      }

      if (isoutside && (not_allowed != direction))
      {
        for (int d2 = 0; d2 < ndim; d2++)
        {
          printout("[warning] outside coord %d %c%c boundary of cell %d. pkttype %d initpos(tmin) %g, vel %g, cellcoordmin %g, cellcoordmax %g. Abort?\n",
                   d, flip ? '-' : '+', grid::coordlabel[d], cellindex, pkt_ptr->type, initpos[d2], vel[d2], grid::get_cellcoordmin(cellindex, d2) / globals::tmin * tstart, cellcoordmax[d2] / globals::tmin * tstart);
        }
        printout("globals::tmin %g tstart %g tstart/globals::tmin %g tdecay %g\n", globals::tmin, tstart, tstart/globals::tmin, pkt_ptr->tdecay);
        // printout("[warning] pkt_ptr->number %d\n", pkt_ptr->number);
        if (flip == 0)
          printout("[warning] delta %g\n", cellcoordmax[d] - (initpos[d] * globals::tmin / tstart));
        else
          printout("[warning] delta %g\n",  (initpos[d] * globals::tmin / tstart) - grid::get_cellcoordmin(cellindex, d));

        printout("[warning] dir [%g, %g, %g]\n", pkt_ptr->dir[0], pkt_ptr->dir[1], pkt_ptr->dir[2]);
        if ((vel[d] - (initpos[d] / tstart)) > 0)
        {
          if ((grid::get_cellcoordpointnum(cellindex, d) == (grid::ncoordgrid[d] - 1) && cellindexdiff > 0) ||
              (grid::get_cellcoordpointnum(cellindex, d) == 0 && cellindexdiff < 0))
          {
            *snext = -99;
            return 0;
          }
          else
          {
            *snext = pkt_ptr->where + cellindexdiff;
            pkt_ptr->last_cross = invdirection;
            return 0;
          }
        }
        else
        {
          not_allowed = direction;
        }
      }
    }
  }


  if (false)
  {
    printout("pkt_ptr->number %d\n", pkt_ptr->number);
    printout("delta1x %g delta2x %g\n",  (initpos[0] * globals::tmin/tstart)-grid::get_cellcoordmin(cellindex, 0), cellcoordmax[0] - (initpos[0] * globals::tmin/tstart));
    printout("delta1y %g delta2y %g\n",  (initpos[1] * globals::tmin/tstart)-grid::get_cellcoordmin(cellindex, 1), cellcoordmax[1] - (initpos[1] * globals::tmin/tstart));
    printout("delta1z %g delta2z %g\n",  (initpos[2] * globals::tmin/tstart)-grid::get_cellcoordmin(cellindex, 2), cellcoordmax[2] - (initpos[2] * globals::tmin/tstart));
    printout("dir [%g, %g, %g]\n", pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
  }

  double t_plus_coordboundary[3];  // time to reach the cell's upper boundary on each coordinate
  double t_minus_coordboundary[3];  // likewise, the lower boundaries (smallest x,y,z or radius value in the cell)
  if (grid::grid_type == GRID_SPHERICAL1D)
  {
    not_allowed = NONE; // we will handle this separately by setting d_minus and d_plus negative for invalid directions
    const double r_inner = grid::get_cellcoordmin(cellindex, 0) * tstart / globals::tmin;

    const double d_minus = (r_inner > 0.) ? get_shellcrossdist(pkt_ptr->pos, pkt_ptr->dir, r_inner, true, tstart) : -1.;
    t_minus_coordboundary[0] = d_minus / globals::CLIGHT_PROP;

    const double r_outer = cellcoordmax[0] * tstart / globals::tmin;
    const double d_plus = get_shellcrossdist(pkt_ptr->pos, pkt_ptr->dir, r_outer, false, tstart);
    t_plus_coordboundary[0] = d_plus / globals::CLIGHT_PROP;

    // printout("cell %d\n", pkt_ptr->where);
    // printout("initradius %g: velrad %g\n", initpos[0], vel[0]);
    // printout("d_plus %g d_minus %g \n", d_plus, d_minus);
    // printout("t_plus %g t_minus %g \n", t_plus_coordboundary[0], t_minus_coordboundary[0]);
    // printout("cellrmin %g cellrmax %g\n",
    //          grid::get_cellcoordmin(cellindex, 0) * tstart / globals::tmin, cellcoordmax[0] * tstart / globals::tmin);
    // printout("tstart %g\n", tstart);
  }
  else
  {
    for (int d = 0; d < 3; d++)
    {
      t_plus_coordboundary[d] = ((initpos[d] - (vel[d] * tstart)) / (cellcoordmax[d] - (vel[d] * globals::tmin)) * globals::tmin) - tstart;
      t_minus_coordboundary[d] = ((initpos[d] - (vel[d] * tstart)) / (grid::get_cellcoordmin(cellindex, d) - (vel[d] * globals::tmin)) * globals::tmin) - tstart;
    }
  }

  // printout("comparing distances. not_allowed = %d\n", not_allowed);
  //We now need to identify the shortest +ve time - that's the one we want.
  int choice = 0;         ///just a control variable to
  double time = 1.e99;
  //close = 1.e99;
  //printout("bondary.c check value of not_allowed = %d\n",not_allowed);
  for (int d = 0; d < ndim; d++)
  {
    if ((t_plus_coordboundary[d] > 0) && (t_plus_coordboundary[d] < time) && (not_allowed != negdirections[d]))
    {
      choice = posdirections[d];
      time = t_plus_coordboundary[d];
      // equivalently if (nxyz[d] == (grid::ncoordgrid[d] - 1))
      // if (grid::get_cellcoordmin(cellindex, d) + 1.5 * grid::wid_init > coordmax[d])
      if (grid::get_cellcoordpointnum(cellindex, d) == (grid::ncoordgrid[d] - 1))
      {
        *snext = -99;
      }
      else
      {
        *snext = pkt_ptr->where + grid::get_coordcellindexincrement(d);
        pkt_ptr->last_cross = posdirections[d];
      }
    }

    if ((t_minus_coordboundary[d] > 0) && (t_minus_coordboundary[d] < time) && (not_allowed != posdirections[d]))
    {
      choice = negdirections[d];
      time = t_minus_coordboundary[d];
      // equivalently if (nxyz[d] == 0)
      // if (grid::get_cellcoordmin(cellindex, d) < - coordmax[d] + 0.5 * grid::wid_init)
      if (grid::get_cellcoordpointnum(cellindex, d) == 0)
      {
        *snext = -99;
      }
      else
      {
        *snext = pkt_ptr->where - grid::get_coordcellindexincrement(d);
        pkt_ptr->last_cross = negdirections[d];
      }
    }
  }

  if (choice == 0)
  // if (pkt_ptr->number == 10043)
  {
    printout("Something wrong in boundary crossing - didn't find anything.\n");
    printout("packet %d cell %d or %d\n", pkt_ptr->number, cellindex, pkt_ptr->where);
    printout("choice %d\n", choice);
    printout("globals::tmin %g tstart %g\n", globals::tmin, tstart);
    printout("not_allowed %d, type %d\n", not_allowed, pkt_ptr->type);
    for (int d2 = 0; d2 < 3; d2++)
    {
      printout("coord %d: initpos %g dir %g\n",
               d2, pkt_ptr->pos[d2], pkt_ptr->dir[d2]);
    }
    printout("|initpos| %g |dir| %g |pos.dir| %g\n",
             vec_len(pkt_ptr->pos), vec_len(pkt_ptr->dir), dot(pkt_ptr->pos, pkt_ptr->dir));
    for (int d2 = 0; d2 < ndim; d2++)
    {
      printout("coord %d: txyz_plus %g txyz_minus %g \n", d2, t_plus_coordboundary[d2], t_minus_coordboundary[d2]);
      printout("coord %d: cellcoordmin %g cellcoordmax %g\n",
               d2, grid::get_cellcoordmin(cellindex, d2) * tstart / globals::tmin, cellcoordmax[d2] * tstart / globals::tmin);
    }
    printout("tstart %g\n", tstart);

    // abort();
  }


  // Now we know what happens. The distance to crossing is....
  double distance = globals::CLIGHT_PROP * time;
  // printout("boundary_cross: time %g distance %g\n", time, distance);
  // closest = close;

  return distance;
}


__host__ __device__
void change_cell(PKT *pkt_ptr, int snext, double t_current)
/// Routine to take a packet across a boundary.
{
  assert_always(pkt_ptr->prop_time == t_current);
  if (false)
  {
    const int cellindex = pkt_ptr->where;
    printout("[debug] cellnumber %d nne %g\n",cellindex,grid::get_nne(grid::get_cell_modelgridindex(cellindex)));
    printout("[debug] snext %d\n",snext);
  }

  if (snext == -99)
  {
    // Then the packet is exiting the grid. We need to record
    // where and at what time it leaves the grid.
    pkt_ptr->escape_type = pkt_ptr->type;
    pkt_ptr->escape_time = pkt_ptr->prop_time;
    pkt_ptr->type = TYPE_ESCAPE;
    safeincrement(globals::nesc);
  }
  else
  {
    // Just need to update "where".
    // const int cellnum = pkt_ptr->where;
    // const int old_mgi = grid::get_cell_modelgridindex(cellnum);
    pkt_ptr->where = snext;
    // const int mgi = grid::get_cell_modelgridindex(snext);

    stats::increment(stats::COUNTER_CELLCROSSINGS);
  }
}


// static int locate(const PKT *pkt_ptr, double t_current)
// /// Routine to return which grid cell the packet is in.
// {
//   // Cheap and nasty version for now - assume a uniform grid.
//   int xx = (pkt_ptr->pos[0] - (globals::cell[0].pos_init[0]*t_current/globals::tmin)) / (grid::wid_init*t_current/globals::tmin);
//   int yy = (pkt_ptr->pos[1] - (globals::cell[0].pos_init[1]*t_current/globals::tmin)) / (grid::wid_init*t_current/globals::tmin);
//   int zz = (pkt_ptr->pos[2] - (globals::cell[0].pos_init[2]*t_current/globals::tmin)) / (grid::wid_init*t_current/globals::tmin);
//
//   return xx + (grid::ncoordgrid[0] * yy) + (grid::ncoordgrid[0] * grid::ncoordgrid[1] * zz);
// }
