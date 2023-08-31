#include "boundary.h"

#include <gsl/gsl_blas.h>

#include <limits>

#include "grid.h"
#include "rpkt.h"
#include "sn3d.h"
#include "stats.h"
#include "update_packets.h"
#include "vectors.h"

static auto expanding_shell_intersection(const double pos[3], const double dir[3], const double shellradiuststart,
                                         const bool isinnerboundary, const double tstart) -> double
// find the closest forward distance to the intersection of a ray with an expanding spherical shell
// return -1 if there are no forward intersections (or if the intersection is tangential to the shell)
{
  assert_always(shellradiuststart > 0);
  constexpr bool debug = false;
  if constexpr (debug) {
    printout("expanding_shell_intersection isinnerboundary %d\n", isinnerboundary);
    printout("shellradiuststart %g tstart %g len(pos) %g\n", shellradiuststart, tstart, vec_len(pos));
  }
  const double speed = vec_len(dir) * CLIGHT_PROP;  // hopefully this is the same as CLIGHT_PROP
  const double a = dot(dir, dir) - pow(shellradiuststart / tstart / speed, 2);
  const double b = 2 * (dot(dir, pos) - pow(shellradiuststart, 2) / tstart / speed);
  const double c = dot(pos, pos) - pow(shellradiuststart, 2);

  const double discriminant = pow(b, 2) - 4 * a * c;

  if (discriminant < 0) {
    // no intersection
    assert_always(isinnerboundary);
    assert_always(shellradiuststart < vec_len(pos));
    if constexpr (debug) {
      printout("no intersection\n");
    }
    return -1;
  }
  if (discriminant > 0) {
    // two intersections
    double dist1 = (-b + sqrt(discriminant)) / 2 / a;
    double dist2 = (-b - sqrt(discriminant)) / 2 / a;

    double posfinal1[3];
    double posfinal2[3];

    for (int d = 0; d < 3; d++) {
      posfinal1[d] = pos[d] + dist1 * dir[d];
      posfinal2[d] = pos[d] + dist2 * dir[d];
    }

    const double shellradiusfinal1 = shellradiuststart / tstart * (tstart + dist1 / speed);
    const double shellradiusfinal2 = shellradiuststart / tstart * (tstart + dist2 / speed);
    // printout("solution1 dist1 %g radiusfinal1 %g shellradiusfinal1 %g\n", dist1, vec_len(posfinal1),
    // shellradiusfinal1); printout("solution2 dist2 %g radiusfinal2 %g shellradiusfinal2 %g\n", dist2,
    // vec_len(posfinal2), shellradiusfinal2);
    assert_always(fabs(vec_len(posfinal1) / shellradiusfinal1 - 1.) < 1e-3);
    assert_always(fabs(vec_len(posfinal2) / shellradiusfinal2 - 1.) < 1e-3);

    // invalidate any solutions that require entering the boundary from the wrong radial direction
    if (isinnerboundary) {
      if (dot(posfinal1, dir) > 0.) {
        dist1 = -1;
      }
      if (dot(posfinal2, dir) > 0.) {
        dist2 = -1;
      }
    } else {
      if (dot(posfinal1, dir) < 0.) {
        dist1 = -1;
      }
      if (dot(posfinal2, dir) < 0.) {
        dist2 = -1;
      }
    }

    // negative d means in the reverse direction along the ray
    // ignore negative d values, and if two are positive then return the smaller one
    if (dist1 < 0 && dist2 < 0) {
      return -1;
    }
    if (dist2 < 0) {
      return dist1;
    }
    if (dist1 < 0) {
      return dist2;
    }
    return fmin(dist1, dist2);

  }  // exactly one intersection
  // ignore this and don't change which cell the packet is in
  assert_always(shellradiuststart <= vec_len(pos));
  printout("single intersection\n");
  return -1.;
}

auto boundary_cross(struct packet *const pkt_ptr, int *snext) -> double
/// Basic routine to compute distance to a cell boundary.
{
  const double tstart = pkt_ptr->prop_time;

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
  double pktposgridcoord[3] = {0};  // pkt_ptr->pos converted to propagation grid coordinates
  double cellcoordmax[3] = {0};
  double pktvelgridcoord[3] = {0};  // pkt_ptr->dir * CLIGHT_PROP converted to grid coordinates

  if constexpr (GRID_TYPE == GRID_CARTESIAN3D) {
    // keep xyz Cartesian coordinates.
    for (int d = 0; d < ndim; d++) {
      pktposgridcoord[d] = pkt_ptr->pos[d];
      cellcoordmax[d] = grid::get_cellcoordmax(cellindex, d);
      pktvelgridcoord[d] = pkt_ptr->dir[d] * CLIGHT_PROP;
    }
  } else if constexpr (GRID_TYPE == GRID_CYLINDRICAL2D) {
    // xy plane radius and z position
    pktposgridcoord[0] = std::sqrt(std::pow(pkt_ptr->pos[0], 2) + std::pow(pkt_ptr->pos[1], 2));
    cellcoordmax[0] = grid::get_cellcoordmax(cellindex, 0);
    pktvelgridcoord[0] =
        (pkt_ptr->pos[0] * pkt_ptr->dir[0] + pkt_ptr->pos[1] * pkt_ptr->dir[1]) / pktposgridcoord[0] * CLIGHT_PROP;

    // second cylindrical coordinate is z
    pktposgridcoord[1] = pkt_ptr->pos[2];
    cellcoordmax[1] = grid::get_cellcoordmax(cellindex, 2);
    pktvelgridcoord[1] = pkt_ptr->dir[2] * CLIGHT_PROP;

  } else if constexpr (GRID_TYPE == GRID_SPHERICAL1D) {
    // the only coordinate is radius from the origin
    pktposgridcoord[0] = vec_len(pkt_ptr->pos);
    cellcoordmax[0] = grid::get_cellcoordmax(cellindex, 0);
    pktvelgridcoord[0] = dot(pkt_ptr->pos, pkt_ptr->dir) / pktposgridcoord[0] * CLIGHT_PROP;  // radial velocity
  } else {
    assert_always(false);
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

  // printout("boundary.c: x0 %g, y0 %g, z0 %g\n", initpos[0] initpos[1] initpos[2]);
  // printout("boundary.c: vx %g, vy %g, vz %g\n",vel[0],vel[1],vel[2]);
  // printout("boundary.c: cellxmin %g, cellymin %g, cellzmin %g\n",grid::get_cellcoordmin(cellindex,
  // 0),grid::get_cellcoordmin(cellindex, 1),grid::get_cellcoordmin(cellindex, 2)); printout("boundary.c: cellxmax %g,
  // cellymax %g, cellzmax %g\n",cellcoordmax[0],cellcoordmax[1],cellcoordmax[2]);

  enum cell_boundary last_cross = pkt_ptr->last_cross;
  enum cell_boundary const negdirections[3] = {COORD0_MIN, COORD1_MIN,
                                               COORD2_MIN};  // 'X' might actually be radial coordinate
  enum cell_boundary const posdirections[3] = {COORD0_MAX, COORD1_MAX, COORD2_MAX};

  // printout("checking inside cell boundary\n");
  for (int d = 0; d < ndim; d++) {
    // flip is either zero or one to indicate +ve and -ve boundaries along the selected axis
    for (int flip = 0; flip < 2; flip++) {
      enum cell_boundary const direction = flip != 0 ? posdirections[d] : negdirections[d];
      enum cell_boundary const invdirection = flip == 0 ? posdirections[d] : negdirections[d];
      const int cellindexstride =
          flip != 0 ? -grid::get_coordcellindexincrement(d) : grid::get_coordcellindexincrement(d);

      bool isoutside_thisside = false;
      if (flip != 0) {
        isoutside_thisside = pktposgridcoord[d] < (grid::get_cellcoordmin(cellindex, d) / globals::tmin * tstart -
                                                   10.);  // 10 cm accuracy tolerance
      } else {
        isoutside_thisside = pktposgridcoord[d] > (cellcoordmax[d] / globals::tmin * tstart + 10.);
      }

      if (isoutside_thisside && (last_cross != direction)) {
        // for (int d2 = 0; d2 < ndim; d2++)
        const int d2 = d;
        {
          printout(
              "[warning] packet %d outside coord %d %c%c boundary of cell %d. pkttype %d initpos(tmin) %g, vel %g, "
              "cellcoordmin %g, cellcoordmax %g\n",
              pkt_ptr->number, d, flip != 0 ? '-' : '+', grid::coordlabel[d], cellindex, pkt_ptr->type,
              pktposgridcoord[d2], pktvelgridcoord[d2], grid::get_cellcoordmin(cellindex, d2) / globals::tmin * tstart,
              cellcoordmax[d2] / globals::tmin * tstart);
        }
        printout("globals::tmin %g tstart %g tstart/globals::tmin %g tdecay %g\n", globals::tmin, tstart,
                 tstart / globals::tmin, pkt_ptr->tdecay);
        // printout("[warning] pkt_ptr->number %d\n", pkt_ptr->number);
        if (flip != 0) {
          printout("[warning] delta %g\n",
                   (pktposgridcoord[d] * globals::tmin / tstart) - grid::get_cellcoordmin(cellindex, d));
        } else {
          printout("[warning] delta %g\n", cellcoordmax[d] - (pktposgridcoord[d] * globals::tmin / tstart));
        }

        printout("[warning] dir [%g, %g, %g]\n", pkt_ptr->dir[0], pkt_ptr->dir[1], pkt_ptr->dir[2]);
        if ((pktvelgridcoord[d] - (pktposgridcoord[d] / tstart)) > 0) {
          if ((grid::get_cellcoordpointnum(cellindex, d) == (grid::ncoordgrid[d] - 1) && cellindexstride > 0) ||
              (grid::get_cellcoordpointnum(cellindex, d) == 0 && cellindexstride < 0)) {
            printout("escaping packet\n");
            *snext = -99;
            return 0;
          }
          *snext = pkt_ptr->where + cellindexstride;
          pkt_ptr->last_cross = invdirection;
          printout("[warning] swapping packet cellindex from %d to %d and setting last_cross to %d\n", pkt_ptr->where,
                   *snext, pkt_ptr->last_cross);
          return 0;
        }
        printout("pretending last_cross is %d\n", direction);
        last_cross = direction;
      }
    }
  }

  // printout("pkt_ptr->number %d\n", pkt_ptr->number);
  // printout("delta1x %g delta2x %g\n",  (initpos[0] * globals::tmin/tstart)-grid::get_cellcoordmin(cellindex, 0),
  // cellcoordmax[0] - (initpos[0] * globals::tmin/tstart)); printout("delta1y %g delta2y %g\n",  (initpos[1] *
  // globals::tmin/tstart)-grid::get_cellcoordmin(cellindex, 1), cellcoordmax[1] - (initpos[1] * globals::tmin/tstart));
  // printout("delta1z %g delta2z %g\n",  (initpos[2] * globals::tmin/tstart)-grid::get_cellcoordmin(cellindex, 2),
  // cellcoordmax[2] - (initpos[2] * globals::tmin/tstart)); printout("dir [%g, %g, %g]\n",
  // pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);

  double t_coordmaxboundary[3];  // time to reach the cell's upper boundary on each coordinate
  double t_coordminboundary[3];  // likewise, the lower boundaries (smallest x,y,z or radius value in the cell)
  if constexpr (GRID_TYPE == GRID_SPHERICAL1D) {
    last_cross = NONE;  // we will handle this separately by setting d_inner and d_outer negative for invalid directions

    const double r_outer = cellcoordmax[0] * tstart / globals::tmin;
    const double d_outer = expanding_shell_intersection(pkt_ptr->pos, pkt_ptr->dir, r_outer, false, tstart);
    t_coordmaxboundary[0] = d_outer / CLIGHT_PROP;

    const double r_inner = grid::get_cellcoordmin(cellindex, 0) * tstart / globals::tmin;
    const double d_inner =
        (r_inner > 0.) ? expanding_shell_intersection(pkt_ptr->pos, pkt_ptr->dir, r_inner, true, tstart) : -1.;
    t_coordminboundary[0] = d_inner / CLIGHT_PROP;

  } else if constexpr (GRID_TYPE == GRID_CYLINDRICAL2D) {
    // coordinate 0 is radius in x-y plane, coord 1 is z
    last_cross = NONE;  // we will handle this separately by setting d_inner and d_outer negative for invalid directions

    // to get the cylindrical intersection, get the spherical intersection when Z components are zero
    const double pktposnoz[3] = {pkt_ptr->pos[0], pkt_ptr->pos[1], 0.};
    const double pktdirnoz[3] = {pkt_ptr->dir[0], pkt_ptr->dir[1], 0.};

    const double r_inner = grid::get_cellcoordmin(cellindex, 0) * tstart / globals::tmin;
    const double d_inner =
        (r_inner > 0.) ? expanding_shell_intersection(pktposnoz, pktdirnoz, r_inner, true, tstart) : -1.;
    t_coordminboundary[0] = d_inner / CLIGHT_PROP;

    const double r_outer = cellcoordmax[0] * tstart / globals::tmin;
    const double d_outer = expanding_shell_intersection(pktposnoz, pktdirnoz, r_outer, false, tstart);
    t_coordmaxboundary[0] = d_outer / CLIGHT_PROP;

    // z boundaries same as Cartesian
    t_coordminboundary[1] =
        ((pktposgridcoord[2] - (pktvelgridcoord[2] * tstart)) /
         ((grid::get_cellcoordmin(cellindex, 1)) - (pktvelgridcoord[2] * globals::tmin)) * globals::tmin) -
        tstart;

    t_coordmaxboundary[1] = ((pktposgridcoord[2] - (pktvelgridcoord[2] * tstart)) /
                             ((cellcoordmax[1]) - (pktvelgridcoord[2] * globals::tmin)) * globals::tmin) -
                            tstart;
  } else if constexpr (GRID_TYPE == GRID_CARTESIAN3D) {
    // const double overshoot = grid::wid_init(cellindex) * 2e-7;
    constexpr double overshoot = 0.;
    for (int d = 0; d < 3; d++) {
      t_coordmaxboundary[d] = ((pktposgridcoord[d] - (pktvelgridcoord[d] * tstart)) /
                               ((cellcoordmax[d] + overshoot) - (pktvelgridcoord[d] * globals::tmin)) * globals::tmin) -
                              tstart;

      t_coordminboundary[d] =
          ((pktposgridcoord[d] - (pktvelgridcoord[d] * tstart)) /
           ((grid::get_cellcoordmin(cellindex, d) - overshoot) - (pktvelgridcoord[d] * globals::tmin)) *
           globals::tmin) -
          tstart;
    }
  } else {
    assert_always(false);
  }

  // printout("comparing distances. last_cross = %d\n", last_cross);
  // We now need to identify the shortest +ve time - that's the one we want.
  int choice = 0;  /// just a control variable to
  double time = std::numeric_limits<double>::max();
  // close = 1.e99;
  // printout("bondary.c check value of last_cross = %d\n",last_cross);
  for (int d = 0; d < ndim; d++) {
    if ((t_coordmaxboundary[d] > 0) && (t_coordmaxboundary[d] < time) && (last_cross != negdirections[d])) {
      choice = posdirections[d];
      time = t_coordmaxboundary[d];
      // equivalently if (nxyz[d] == (grid::ncoordgrid[d] - 1))
      // if (grid::get_cellcoordmin(cellindex, d) + 1.5 * grid::wid_init > coordmax[d])
      if (grid::get_cellcoordpointnum(cellindex, d) == (grid::ncoordgrid[d] - 1)) {
        *snext = -99;
      } else {
        *snext = pkt_ptr->where + grid::get_coordcellindexincrement(d);
        pkt_ptr->last_cross = posdirections[d];
      }
    }

    if ((t_coordminboundary[d] > 0) && (t_coordminboundary[d] < time) && (last_cross != posdirections[d])) {
      choice = negdirections[d];
      time = t_coordminboundary[d];
      // equivalently if (nxyz[d] == 0)
      // if (grid::get_cellcoordmin(cellindex, d) < - coordmax[d] + 0.5 * grid::wid_init)
      if (grid::get_cellcoordpointnum(cellindex, d) == 0) {
        *snext = -99;
      } else {
        *snext = pkt_ptr->where - grid::get_coordcellindexincrement(d);
        pkt_ptr->last_cross = negdirections[d];
      }
    }
  }

  if (choice == 0) {
    printout("Something wrong in boundary crossing - didn't find anything.\n");
    printout("packet %d cell %d or %d\n", pkt_ptr->number, cellindex, pkt_ptr->where);
    printout("choice %d\n", choice);
    printout("globals::tmin %g tstart %g\n", globals::tmin, tstart);
    printout("last_cross %d, type %d\n", last_cross, pkt_ptr->type);
    for (int d2 = 0; d2 < 3; d2++) {
      printout("coord %d: initpos %g dir %g\n", d2, pkt_ptr->pos[d2], pkt_ptr->dir[d2]);
    }
    printout("|initpos| %g |dir| %g |pos.dir| %g\n", vec_len(pkt_ptr->pos), vec_len(pkt_ptr->dir),
             dot(pkt_ptr->pos, pkt_ptr->dir));
    for (int d2 = 0; d2 < ndim; d2++) {
      printout("coord %d: txyz_plus %g txyz_minus %g \n", d2, t_coordmaxboundary[d2], t_coordminboundary[d2]);
      printout("coord %d: cellcoordmin %g cellcoordmax %g\n", d2,
               grid::get_cellcoordmin(cellindex, d2) * tstart / globals::tmin,
               cellcoordmax[d2] * tstart / globals::tmin);
    }
    printout("tstart %g\n", tstart);

    // abort();
  }

  // Now we know what happens. The distance to crossing is....
  double const distance = CLIGHT_PROP * time;
  // printout("boundary_cross: time %g distance %g\n", time, distance);
  // closest = close;

  return distance;
}

void change_cell(struct packet *pkt_ptr, int snext)
/// Routine to take a packet across a boundary.
{
  // const int cellindex = pkt_ptr->where;
  // printout("[debug] cellnumber %d nne %g\n", cellindex, grid::get_nne(grid::get_cell_modelgridindex(cellindex)));
  // printout("[debug] snext %d\n", snext);

  if (snext == -99) {
    // Then the packet is exiting the grid. We need to record
    // where and at what time it leaves the grid.
    pkt_ptr->escape_type = pkt_ptr->type;
    pkt_ptr->escape_time = pkt_ptr->prop_time;
    pkt_ptr->type = TYPE_ESCAPE;
    safeincrement(globals::nesc);
  } else {
    // Just need to update "where".
    pkt_ptr->where = snext;

    stats::increment(stats::COUNTER_CELLCROSSINGS);
  }
}

static auto get_cell(const double pos[3], double t) -> int
/// identify the cell index from a position and a time.
{
  assert_always(GRID_TYPE == GRID_CARTESIAN3D);  // other grid types not implemented yet

  const double trat = t / globals::tmin;
  const int nx = static_cast<int>((pos[0] - (grid::get_cellcoordmin(0, 0) * trat)) / (grid::wid_init(0, 0) * trat));
  const int ny = static_cast<int>((pos[1] - (grid::get_cellcoordmin(0, 1) * trat)) / (grid::wid_init(0, 1) * trat));
  const int nz = static_cast<int>((pos[2] - (grid::get_cellcoordmin(0, 2) * trat)) / (grid::wid_init(0, 1) * trat));

  const int cellindex = nx + (grid::ncoordgrid[0] * ny) + (grid::ncoordgrid[0] * grid::ncoordgrid[1] * nz);

  // do a check
  for (int n = 0; n < grid::get_ngriddimensions(); n++) {
    assert_always(pos[n] >= grid::get_cellcoordmin(cellindex, n));
    assert_always(pos[n] <= grid::get_cellcoordmax(cellindex, n));
  }
  return cellindex;
}