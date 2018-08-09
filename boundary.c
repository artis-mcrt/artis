#include "sn3d.h"
#include "assert.h"
#include "boundary.h"
#include "grid_init.h"
#include "rpkt.h"
#include "update_packets.h"
#include "vectors.h"
// #include <gsl/gsl_poly.h>
#include <gsl/gsl_blas.h>



static double get_shellcrossdist(
  const double pos[3], const double dir[3], const double shellradius, const bool isinnerboundary, const double tstart, bool debug)
// find the closest forward distance to the intersection of a ray with an expanding spherical shell
// return -1 if there are no forward intersections (or if the intersection is tangential to the shell)
{
  assert(shellradius > 0);
  if (debug)
  {
    printout("get_shellcrossdist isinnerboundary %d\n", isinnerboundary);
    printout("shellradius %g tstart %g len(pos) %g\n", shellradius, tstart, vec_len(pos));
  }
  const double speed = vec_len(dir) * CLIGHT_PROP;
  const double a = dot(dir, dir) - pow(shellradius / tstart / speed, 2);
  const double b = 2 * (dot(dir, pos) - pow(shellradius, 2) / tstart / speed);
  const double c = dot(pos, pos) - pow(shellradius, 2);

  const double discriminant = pow(b, 2) - 4 * a * c;

  if (discriminant < 0)
  {
    // no intersection
    assert(shellradius < vec_len(pos));
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
    cblas_dcopy(3, pos, 1, posfinal1, 1);     // posfinal1 = pos
    cblas_daxpy(3, d1, dir, 1, posfinal1, 1); // posfinal1 += d1 * dir;

    double posfinal2[3];
    cblas_dcopy(3, pos, 1, posfinal2, 1);
    cblas_daxpy(3, d2, dir, 1, posfinal2, 1);

    const double shellradiusfinal1 = shellradius / tstart * (tstart + d1 / speed);
    const double shellradiusfinal2 = shellradius / tstart * (tstart + d2 / speed);
    // printout("solution1 d1 %g radiusfinal1 %g shellradiusfinal1 %g\n", d1, vec_len(posfinal1), shellradiusfinal1);
    // printout("solution2 d2 %g radiusfinal2 %g shellradiusfinal2 %g\n", d2, vec_len(posfinal2), shellradiusfinal2);
    assert(fabs(vec_len(posfinal1) / shellradiusfinal1 - 1.) < 1e-3);
    assert(fabs(vec_len(posfinal2) / shellradiusfinal2 - 1.) < 1e-3);

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
    assert(shellradius <= vec_len(pos));
    printout("single intersection\n");
    return -1.;
  }

}


double boundary_cross(PKT *restrict const pkt_ptr, const double tstart, int *snext)
/// Basic routine to compute distance to a cell boundary.
{
  //double close, close_try;

  /** There are six possible boundary crossings. Each of the three
  cartesian coordinates may be taken in turn. For x, the packet
  trajectory is
  x = x0 + (dir.x) * c * (t - tstart)
  the boundries follow
  x+/- = x+/-(tmin) * (t/tmin)
  so the crossing occurs when
  t = (x0 - (dir.x)*c*tstart)/(x+/-(tmin)/tmin - (dir.x)c)
   */

  /* Modified so that it also returns the distance to the closest cell
  boundary, regardless of direction. */

  // d is used to loop over the coordinate indicies 0,1,2 for x,y,z

  const int ndim = (grid_type == GRID_SPHERICAL1D) ? 1 : 3;

  const int cellindex = pkt_ptr->where;

  double initpos[ndim];  // [x0, y0, z0]
  double cellcoordmin[ndim];
  double cellcoordmax[ndim];
  double vel[ndim];

  if (grid_type == GRID_UNIFORM)
  {
    // XYZ coordinates
    for (int d = 0; d < ndim; d++)
    {
      initpos[d] = pkt_ptr->pos[d];
      cellcoordmin[d] = get_cellcoordmin(cellindex, d);
      cellcoordmax[d] = cellcoordmin[d] + wid_init(0);
      vel[d] = pkt_ptr->dir[d] * CLIGHT_PROP;
    }
  }
  else if (grid_type == GRID_SPHERICAL1D)
  {
    // the only coordinate is radius from the origin
    initpos[0] = vec_len(pkt_ptr->pos);
    cellcoordmin[0] = get_cellcoordmin(cellindex, 0);
    cellcoordmax[0] = cellcoordmin[0] + wid_init(cellindex);
    vel[0] = dot(pkt_ptr->pos, pkt_ptr->dir) / vec_len(pkt_ptr->pos) * CLIGHT_PROP; // radial velocity
  }

  // how much do we change the cellindex to move in the the x, y, z directions
  const int cellindexincrement[3] = {ncoordgrid[2] * ncoordgrid[1], ncoordgrid[2], 1};

  //printout("boundary.c: x0 %g, y0 %g, z0 %g\n", initpos[0] initpos[1] initpos[2]);

  //printout("boundary.c: vx %g, vy %g, vz %g\n",vel[0],vel[1],vel[2]);

  //printout("boundary.c: cellxmin %g, cellymin %g, cellzmin %g\n",cellcoordmin[0],cellcoordmin[1],cellcoordmin[2]);

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
      const int cellindexdiff = flip ? - cellindexincrement[d] : cellindexincrement[d];

      bool isoutside;
      if (flip)
        isoutside = initpos[d] - (cellcoordmin[d] / tmin * tstart) < -10.; // 10 cm accuracy tolerance
      else
        isoutside = initpos[d] - (cellcoordmax[d] / tmin * tstart) > -10.;

      if (isoutside && (not_allowed != direction))
      {
        for (int d2 = 0; d2 < ndim; d2++)
        {
          printout("[warning] outside coord %d boundary of cell %d. type %d initpos %g, vel %g, cellcoordmin %g, cellcoordmax %g. Abort?\n",
                   d, cellindex, pkt_ptr->type, initpos[d2] * tmin/tstart, vel[d2], cellcoordmin[d2], cellcoordmax[d2]);
        }
        printout("tmin %g tstart %g tstart/tmin %g tdecay %g\n", tmin, tstart, tstart/tmin, pkt_ptr->tdecay);
        // printout("[warning] pkt_ptr->number %d\n", pkt_ptr->number);
        if (flip == 0)
          printout("[warning] delta %g\n", cellcoordmax[d] - (initpos[d] * tmin / tstart));
        else
          printout("[warning] delta %g\n",  (initpos[d] * tmin / tstart) - cellcoordmin[d]);

        printout("[warning] dir [%g, %g, %g]\n", pkt_ptr->dir[0], pkt_ptr->dir[1], pkt_ptr->dir[2]);
        if ((vel[d] - (initpos[d] / tstart)) > 0)
        {
          if ((get_cellcoordpointnum(cellindex, d) == (ncoordgrid[d] - 1) && cellindexdiff > 0) ||
              (get_cellcoordpointnum(cellindex, d) == 0 && cellindexdiff < 0))
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


  #ifdef DEBUG_ON
    if (debuglevel == 2)
    {
      printout("pkt_ptr->number %d\n", pkt_ptr->number);
      printout("delta1x %g delta2x %g\n",  (initpos[0] * tmin/tstart)-cellcoordmin[0], cellcoordmax[0] - (initpos[0] * tmin/tstart));
      printout("delta1y %g delta2y %g\n",  (initpos[1] * tmin/tstart)-cellcoordmin[1], cellcoordmax[1] - (initpos[1] * tmin/tstart));
      printout("delta1z %g delta2z %g\n",  (initpos[2] * tmin/tstart)-cellcoordmin[2], cellcoordmax[2] - (initpos[2] * tmin/tstart));
      printout("dir [%g, %g, %g]\n", pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
    }
  #endif

  double t_plus_coordboundary[ndim];  // time to reach the cell's upper boundary on each coordinate
  double t_minus_coordboundary[ndim];  // likewise, the lower boundaries (smallest x,y,z or radius value in the cell)
  if (grid_type == GRID_SPHERICAL1D)
  {
    not_allowed = NONE; // we will handle this separately by setting d_minus and d_plus negative for invalid directions
    const double r_inner = cellcoordmin[0] * tstart / tmin;

    const double d_minus = (r_inner > 0.) ? get_shellcrossdist(pkt_ptr->pos, pkt_ptr->dir, r_inner, true, tstart, false) : -1.;
    t_minus_coordboundary[0] = d_minus / CLIGHT_PROP;

    const double r_outer = cellcoordmax[0] * tstart / tmin;
    const double d_plus = get_shellcrossdist(pkt_ptr->pos, pkt_ptr->dir, r_outer, false, tstart, false);
    t_plus_coordboundary[0] = d_plus / CLIGHT_PROP;

    // printout("cell %d\n", pkt_ptr->where);
    // printout("initradius %g: velrad %g\n", initpos[0], vel[0]);
    // printout("d_plus %g d_minus %g \n", d_plus, d_minus);
    // printout("t_plus %g t_minus %g \n", t_plus_coordboundary[0], t_minus_coordboundary[0]);
    // printout("cellrmin %g cellrmax %g\n",
    //          cellcoordmin[0] * tstart / tmin, cellcoordmax[0] * tstart / tmin);
    // printout("tstart %g\n", tstart);
  }
  else
  {
    for (int d = 0; d < 3; d++)
    {
      t_plus_coordboundary[d] = ((initpos[d] - (vel[d] * tstart)) / (cellcoordmax[d] - (vel[d] * tmin)) * tmin) - tstart;
      t_minus_coordboundary[d] = ((initpos[d] - (vel[d] * tstart)) / (cellcoordmin[d] - (vel[d] * tmin)) * tmin) - tstart;
    }
  }

  // printout("comparing distances. not_allowed = %d\n", not_allowed);
  /** We now need to identify the shortest +ve time - that's the one we want. */
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
      // equivalently if (nxyz[d] == (ncoordgrid[d] - 1))
      // if (cellcoordmin[d] + 1.5 * wid_init > coordmax[d])
      if (get_cellcoordpointnum(cellindex, d) == (ncoordgrid[d] - 1))
      {
        *snext = -99;
      }
      else
      {
        *snext = pkt_ptr->where + cellindexincrement[d];
        pkt_ptr->last_cross = posdirections[d];
      }
    }

    if ((t_minus_coordboundary[d] > 0) && (t_minus_coordboundary[d] < time) && (not_allowed != posdirections[d]))
    {
      choice = negdirections[d];
      time = t_minus_coordboundary[d];
      // equivalently if (nxyz[d] == 0)
      // if (cellcoordmin[d] < - coordmax[d] + 0.5 * wid_init)
      if (get_cellcoordpointnum(cellindex, d) == 0)
      {
        *snext = -99;
      }
      else
      {
        *snext = pkt_ptr->where - cellindexincrement[d];
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
    printout("tmin %g tstart %g\n", tmin, tstart);
    printout("not_allowed %d or %d, type %d\n", not_allowed, pkt_ptr->type);
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
               d2, cellcoordmin[d2] * tstart / tmin, cellcoordmax[d2] * tstart / tmin);
    }
    printout("tstart %g\n", tstart);

    // abort();
  }


  /** Now we know what happens. The distance to crossing is....*/
  double distance = CLIGHT_PROP * time;
  // printout("boundary_cross: time %g distance %g\n", time, distance);
  //*closest = close;

  // LJS: the cells are now ordered in z then y then x to match 3D model input
  // I think previous versions of artis swapped the z and x axes of the input models
  /** Also we've identified the cell we'll go into. The cells are ordered in x then y then z.
  So if we go up in x we go to cell + 1, up in y we go to cell + ncoordgrid[0], up in z
  we go to cell + (ncoordgrid[0] * ncoordgrid[1]).
  If's going to escape the grid this is flagged with snext = -99. */

  return distance;
}


void change_cell(PKT *restrict pkt_ptr, int snext, bool *end_packet, double t_current)
/// Routine to take a packet across a boundary.
{
  #ifdef DEBUG_ON
    if (debuglevel == 2)
    {
      const int cellindex = pkt_ptr->where;
      printout("[debug] cellnumber %d nne %g\n",cellindex,get_nne(cell[cellindex].modelgridindex));
      printout("[debug] snext %d\n",snext);
    }
  #endif

  if (snext == -99)
  {
      /* Then the packet is exiting the grid. We need to record
    where and at what time it leaves the grid. */
    pkt_ptr->escape_type = pkt_ptr->type;
    pkt_ptr->escape_time = t_current;
    pkt_ptr->type = TYPE_ESCAPE;
    nesc++;
    *end_packet = true;
  }
  else
  {
    /** Just need to update "where".*/
    const int cellnum = pkt_ptr->where;
    const int old_mgi = cell[cellnum].modelgridindex;
    pkt_ptr->where = snext;
    const int mgi = cell[snext].modelgridindex;

    cellcrossings++;

    /// and calculate the continuums opacity in the new cell if the packet is a rpkt
    /// for isothermal homogeneous grids this could be omitted if we neglect the time dependency
    /// as we do it on the rpkts way through a cell
    //if (debuglevel == 2) printout("[debug] calculate_kappa_rpkt after cell crossing\n");
    //if (pkt_ptr->type == TYPE_RPKT) calculate_kappa_rpkt_cont(pkt_ptr,t_current);

    /// check for empty cells
    if (mgi != MMODELGRID)
    {
      if (mgi != old_mgi)
      {
        /// Update the level populations and reset the precalculated rate coefficients
        //printout("change cell: cellnumber %d\n",pkt_ptr->where);
        /// This only needs to be done for non-grey cells
        if (modelgrid[mgi].thick != 1)
        {
          update_cell(mgi);
        }
        /// Using a cellhistory stack, we update the cell data only if it was not
        /// calculated already
        /*
        histindex = search_cellhistory(mgi);
        //printout("change cell: histindex found %d\n",histindex);
        if (histindex < 0)
        {
          //histindex = find_farthestcell(pkt_ptr->where);
          histindex = CELLHISTORYSIZE-1;
          //printout("change cell: histindex farthest %d\n",histindex);
          update_cell(mgi);
        }
        */
      }

      //copy_populations_to_phixslist();
      /// the rpkt's continuum opacity must be updated in any case as it depends on nu
      /// and nu changed after propagation
      if (pkt_ptr->type == TYPE_RPKT)
      {
        #ifdef DEBUG_ON
          if (debuglevel == 2) printout("[debug] calculate_kappa_rpkt after cell crossing\n");
        #endif
        /// This only needs to be done for non-grey cells
        if (modelgrid[mgi].thick != 1)
        {
          calculate_kappa_rpkt_cont(pkt_ptr, t_current);
        }
      }
    }
  }
}



void change_cell_vpkt(PKT *pkt_ptr, int snext, bool *end_packet, double t_current)
/// Routine to take a VIRTUAL packet across a boundary.
{
  //int element, ion, level;

  #ifdef DEBUG_ON
    if (debuglevel == 2)
    {
      const int cellindex = pkt_ptr->where;
      printout("[debug] cellnumber %d nne %g\n",cellindex,get_nne(cell[cellindex].modelgridindex));
      printout("[debug] snext %d\n",snext);
    }
  #endif

  if (snext == -99)
  {
    /* Then the packet is exiting the grid. We need to record
     where and at what time it leaves the grid. */
    pkt_ptr->escape_type = pkt_ptr->type;
    pkt_ptr->escape_time = t_current;
    pkt_ptr->type = TYPE_ESCAPE;
    *end_packet = true;
  }
  else
  {
    /** Just need to update "where".*/
    //int oldpos = pkt_ptr->where;
    //int old_mgi = cell[pkt_ptr->where].modelgridindex;
    pkt_ptr->where = snext;
    // mgi = cell[pkt_ptr->where].modelgridindex;

    cellcrossings++;
  }
}


// static int locate(const PKT *restrict pkt_ptr, double t_current)
// /// Routine to return which grid cell the packet is in.
// {
//   // Cheap and nasty version for now - assume a uniform grid.
//   int xx = (pkt_ptr->pos[0] - (cell[0].pos_init[0]*t_current/tmin)) / (wid_init*t_current/tmin);
//   int yy = (pkt_ptr->pos[1] - (cell[0].pos_init[1]*t_current/tmin)) / (wid_init*t_current/tmin);
//   int zz = (pkt_ptr->pos[2] - (cell[0].pos_init[2]*t_current/tmin)) / (wid_init*t_current/tmin);
//
//   return xx + (ncoordgrid[0] * yy) + (ncoordgrid[0] * ncoordgrid[1] * zz);
// }
