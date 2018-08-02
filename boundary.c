#include "sn3d.h"
#include "assert.h"
#include "boundary.h"
#include "grid_init.h"
#include "rpkt.h"
#include "update_packets.h"


double boundary_cross(PKT *restrict const pkt_ptr, const double tstart, int *snext)
/// Basic routine to compute distance to a call boundary.
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

  const int cellindex = pkt_ptr->where;

  double initpos[3];  // [x0, y0, z0]
  double cellxyzmin[3];
  double cellxyzmax[3];
  double vel[3];
  for (int d = 0; d < 3; d++)
  {
    initpos[d] = pkt_ptr->pos[d];
    cellxyzmin[d] = cell[cellindex].pos_init[d];
    cellxyzmax[d] = cellxyzmin[d] + wid_init;
    vel[d] = pkt_ptr->dir[d] * CLIGHT_PROP;
  }

  // how much do we change the cellindex to move in the the x, y, z directions
  const int cellindexdiffxyz[3] = {1, nxgrid, nxgrid * nygrid};
  const double xyzmax[3] = {xmax, ymax, zmax};

  //printout("boundary.c: x0 %g, y0 %g, z0 %g\n", initpos[0] initpos[1] initpos[2]);

  //printout("boundary.c: vx %g, vy %g, vz %g\n",vel[0],vel[1],vel[2]);

  //printout("boundary.c: cellxmin %g, cellymin %g, cellzmin %g\n",cellxyzmin[0],cellxyzmin[1],cellxyzmin[2]);

  //printout("boundary.c: cellxmax %g, cellymax %g, cellzmax %g\n",cellxyzmax[0],cellxyzmax[1],cellxyzmax[2]);

  // enum cell_boundary not_allowed = pkt_ptr->last_cross;
  enum cell_boundary not_allowed = pkt_ptr->last_cross;
  enum cell_boundary negdirections[3] = {NEG_X, NEG_Y, NEG_Z};
  enum cell_boundary posdirections[3] = {POS_X, POS_Y, POS_Z};

  for (int d = 0; d < 3; d++)
  {
    // flip is either zero or one to indicate one of the two boundaries along the d axis
    for (int flip = 0; flip < 2; flip++)
    {
      enum cell_boundary direction = flip ? posdirections[d] : negdirections[d];
      enum cell_boundary invdirection = !flip ? posdirections[d] : negdirections[d];
      const int cellindexdiff = flip ? - cellindexdiffxyz[d] : cellindexdiffxyz[d];

      bool isoutside;
      if (flip)
        isoutside = (initpos[d] * tmin) - (tstart * cellxyzmin[d]) < (-10. * tmin);
      else
        isoutside = (initpos[d] * tmin) - (tstart * cellxyzmax[d]) > (-10. * tmin);

      if (isoutside && (not_allowed != direction))
      {
        char dimstr[3] = {'X', 'Y', 'Z'};
        for (int d2 = 0; d2 < 3; d2++)
          printout("[warning] outside %s boundary of cell %d. x0 %g, cellxmin %g, cellxmax %g. Abort?\n", dimstr[d], pkt_ptr->where, initpos[d2] * tmin/tstart, cellxyzmin[d2], cellxyzmax[d2]);
        // printout("[warning] pkt_ptr->number %d\n", pkt_ptr->number);
        if (flip == 0)
          printout("[warning] delta %g\n", cellxyzmax[d] - (initpos[d] * tmin / tstart));
        else
          printout("[warning] delta %g\n",  (initpos[d] * tmin / tstart) - cellxyzmin[d]);

        printout("[warning] dir [%g, %g, %g]\n", pkt_ptr->dir[0], pkt_ptr->dir[1], pkt_ptr->dir[2]);
        if (((pkt_ptr->dir[d] * CLIGHT) - (initpos[d] / tstart)) > 0) //TODO: should this be CLIGHT OR CLIGHT_PROP?
        {
          *snext = pkt_ptr->where + cellindexdiff;
          pkt_ptr->last_cross = invdirection;
          return 0;
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
      printout("delta1x %g delta2x %g\n",  (initpos[0] * tmin/tstart)-cellxyzmin[0], cellxyzmax[0] - (initpos[0] * tmin/tstart));
      printout("delta1y %g delta2y %g\n",  (initpos[1] * tmin/tstart)-cellxyzmin[1], cellxyzmax[1] - (initpos[1] * tmin/tstart));
      printout("delta1z %g delta2z %g\n",  (initpos[2] * tmin/tstart)-cellxyzmin[2], cellxyzmax[2] - (initpos[2] * tmin/tstart));
      printout("dir [%g, %g, %g]\n", pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
    }
  #endif

  double t_plusxyz[3];
  double t_minusxyz[3];
  for (int d = 0; d < 3; d++)
  {
    t_plusxyz[d] = ((initpos[d] - (vel[d] * tstart)) / (cellxyzmax[d] - (vel[d] * tmin)) * tmin) - tstart;
    t_minusxyz[d] = ((initpos[d] - (vel[d] * tstart)) / (cellxyzmin[d] - (vel[d] * tmin)) * tmin) - tstart;
  }

  /** We now need to identify the shortest +ve time - that's the one we want. */
  int choice = 0;         ///just a control variable to
  double time = 1.e99;
  //close = 1.e99;
  //printout("bondary.c check value of not_allowed = %d\n",not_allowed);
  for (int d = 0; d < 3; d++)
  {
    if ((t_plusxyz[d] > 0) && (t_plusxyz[d] < time) && (not_allowed != negdirections[d]))
    {
      choice = d * 2 + 1;
      time = t_plusxyz[d];
      // equivalently if (nxyz == (nxgrid - 1))
      if (cellxyzmin[d] + 1.5 * wid_init > xyzmax[d])
      {
        *snext = -99;
      }
      else
      {
        *snext = pkt_ptr->where + cellindexdiffxyz[d];
        pkt_ptr->last_cross = posdirections[d];
      }
    }

    if ((t_minusxyz[d] > 0) && (t_minusxyz[d] < time) && (not_allowed != posdirections[d]))
    {
      choice = d * 2 + 2;
      time = t_minusxyz[d];
      // equivalently if (nxyz == 0)
      if (cellxyzmin[d] < - xyzmax[d] + 0.5 * wid_init)
      {
        *snext = -99;
      }
      else
      {
        *snext = pkt_ptr->where - cellindexdiffxyz[d];
        pkt_ptr->last_cross = negdirections[d];
      }
    }
  }

  if (choice == 0)
  {
    printout("Something wrong in boundary crossing - didn't find anything.\n");
    printout("tx_plus %g tx_minus %g \n", t_plusxyz[0], t_minusxyz[0]);
    printout("ty_plus %g ty_minus %g \n", t_plusxyz[1], t_minusxyz[1]);
    printout("tz_plus %g tz_minus %g \n", t_plusxyz[2], t_minusxyz[2]);
    printout("x0 %g y0 %g z0 %g \n", initpos[0], initpos[1], initpos[2]);
    printout("vx %g vy %g vz %g \n", vel[0], vel[1], vel[2]);
    printout("tstart %g\n", tstart);

    abort();
  }


  /** Now we know what happens. The distance to crossing is....*/
  double distance = CLIGHT_PROP * time;
  //*closest = close;

  /** Also we've identified the cell we'll go into. The cells are ordered in x then y then z.
  So if we go up in x we go to cell + 1, up in y we go to cell + nxgrid, up in z
  we go to cell + (nxgrid * nygrid).
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
    //int oldpos = pkt_ptr->where;
    const int cellnum = pkt_ptr->where;
    int old_mgi = cell[cellnum].modelgridindex;
    pkt_ptr->where = snext;
    int mgi = cell[snext].modelgridindex;
    /*
    double cellwidth = t_current/tmin * wid_init;
    if (pkt_ptr->last_cross == POS_X)
    {
      pkt_ptr->pos[0] -= cellwidth;
      pkt_ptr->last_cross = NEG_X;
    }
    else if (pkt_ptr->last_cross == NEG_X)
    {
      pkt_ptr->pos[0] += cellwidth;
      pkt_ptr->last_cross = POS_X;
    }
    else if (pkt_ptr->last_cross == POS_Y)
    {
      pkt_ptr->pos[1] -= cellwidth;
      pkt_ptr->last_cross = NEG_Y;
    }
    else if (pkt_ptr->last_cross == NEG_Y)
    {
      pkt_ptr->pos[1] += cellwidth;
      pkt_ptr->last_cross = POS_Y;
    }
    else if (pkt_ptr->last_cross == POS_Z)
    {
      pkt_ptr->pos[2] -= cellwidth;
      pkt_ptr->last_cross = NEG_Z;
    }
    else if (pkt_ptr->last_cross == NEG_Z)
    {
      pkt_ptr->pos[2] += cellwidth;
      pkt_ptr->last_cross = POS_Z;
    }
    */
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
          calculate_kappa_rpkt_cont(pkt_ptr,t_current);
        }
      }
    }
    /*
    if (pkt_ptr->last_cross == POS_X)
    {
    pkt_ptr->pos[0] = (cell[snext].pos_init[0])*t_current/tmin;
  }
    else if (pkt_ptr->last_cross == NEG_X)
    {
    pkt_ptr->pos[0] = (cell[snext].pos_init[0]+(wid_init))*t_current/tmin;
  }
    else if (pkt_ptr->last_cross == POS_Y)
    {
    pkt_ptr->pos[1] = (cell[snext].pos_init[1])*t_current/tmin;
  }
    else if (pkt_ptr->last_cross == NEG_Y)
    {
    pkt_ptr->pos[1] = (cell[snext].pos_init[1]+wid_init)*t_current/tmin;
  }
    else if (pkt_ptr->last_cross == POS_Z)
    {
    pkt_ptr->pos[2] = (cell[snext].pos_init[2])*t_current/tmin;
  }
    else if (pkt_ptr->last_cross == NEG_Z)
    {
    pkt_ptr->pos[2] = (cell[snext].pos_init[2]+wid_init)*t_current/tmin;
  }
    else
    {
    printout("Updating cell when it didn't cross a boundary. Abort?\n");
  }
    */

    /*
    if (pkt_ptr->pos[0] < (cell[snext].pos_init[0]*t_current/tmin))
    {
    pkt_ptr->pos[0] = (cell[snext].pos_init[0]+(1.e-6 *wid_init)*t_current/tmin;
  }
    if (pkt_ptr->pos[0] > ((cell[snext].pos_init[0]+wid_init)*t_current/tmin))
    {
    pkt_ptr->pos[0] = (cell[snext].pos_init[0]+(0.999999*wid_init))*t_current/tmin;
  }
    if (pkt_ptr->pos[1] < (cell[snext].pos_init[1]*t_current/tmin))
    {
    pkt_ptr->pos[1] = (cell[snext].pos_init[1]+(1.e-6 *wid_init))*t_current/tmin;
  }
    if (pkt_ptr->pos[1] > ((cell[snext].pos_init[1]+wid_init)*t_current/tmin))
    {
    pkt_ptr->pos[1] = (cell[snext].pos_init[1]+(0.999999*wid_init))*t_current/tmin;
  }
    if (pkt_ptr->pos[2] < (cell[snext].pos_init[2]*t_current/tmin))
    {
    pkt_ptr->pos[2] = (cell[snext].pos_init[2]+(1.e-6 *wid_init))*t_current/tmin;
  }
    if (pkt_ptr->pos[2] > ((cell[snext].pos_init[2]+wid_init2])*t_current/tmin))
    {
    pkt_ptr->pos[2] = (cell[snext].pos_init[2]+(0.999999*wid_init))*t_current/tmin;
  }
    */

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
//   return xx + (nxgrid * yy) + (nxgrid * nygrid * zz);
// }
