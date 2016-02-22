#include "sn3d.h"


///****************************************************************************
/// Basic routine to compute distance to a call boundary.
double boundary_cross (pkt_ptr, tstart, snext)
     PKT *pkt_ptr;
     double tstart;
     int *snext;
     //   double *closest;
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

  double x0 = pkt_ptr->pos[0];
  double y0 = pkt_ptr->pos[1];
  double z0 = pkt_ptr->pos[2];
  //printout("boundary.c: x0 %g, y0 %g, z0 %g\n",x0,y0,z0);

  double vx = pkt_ptr->dir[0] * CLIGHT_PROP;
  double vy = pkt_ptr->dir[1] * CLIGHT_PROP;
  double vz = pkt_ptr->dir[2] * CLIGHT_PROP;
  //printout("boundary.c: vx %g, vy %g, vz %g\n",vx,vy,vz);

  double cellxmin = cell[pkt_ptr->where].pos_init[0];
  double cellymin = cell[pkt_ptr->where].pos_init[1];
  double cellzmin = cell[pkt_ptr->where].pos_init[2];
  //printout("boundary.c: cellxmin %g, cellymin %g, cellzmin %g\n",cellxmin,cellymin,cellzmin);

  double cellxmax = cell[pkt_ptr->where].pos_init[0] + wid_init;
  double cellymax = cell[pkt_ptr->where].pos_init[1] + wid_init;
  double cellzmax = cell[pkt_ptr->where].pos_init[2] + wid_init;
  //printout("boundary.c: cellxmax %g, cellymax %g, cellzmax %g\n",cellxmax,cellymax,cellzmax);

  int not_allowed = pkt_ptr->last_cross;
  if ((((x0 * tmin) - (tstart * cellxmax)) > (-10.*tmin)) && (not_allowed != NEG_X))
  {
    printout("[warning] outside x boundary of cell %d. x0 %g, cellxmin %g, cellxmax %g. Abort?\n", pkt_ptr->where, x0 * tmin/tstart, cellxmin, cellxmax);
    printout("[warning] outside x boundary of cell %d. y0 %g, cellymin %g, cellymax %g. Abort?\n", pkt_ptr->where, y0 * tmin/tstart, cellymin, cellymax);
    printout("[warning] outside x boundary of cell %d. z0 %g, cellzmin %g, cellzmax %g. Abort?\n", pkt_ptr->where, z0 * tmin/tstart, cellzmin, cellzmax);
    printout("[warning] pkt_ptr->number %d\n", pkt_ptr->number);
    printout("[warning] delta %g\n", cellxmax - (x0 * tmin/tstart));
    printout("[warning] dir [%g, %g, %g]\n", pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
    if (((pkt_ptr->dir[0] * CLIGHT) - (x0/tstart)) > 0)
    {
      *snext = pkt_ptr->where + 1;
      pkt_ptr->last_cross = POS_X;
      return(0);
    }
    else
    {
      not_allowed = NEG_X;
    }
      //debuglevel = 2;
  }
  if ((((x0 * tmin) - (tstart * cellxmin)) < (10.*tmin)) && (not_allowed != POS_X))
  {
    printout("[warning] outside x boundary of cell %d. x0 %g, cellxmin %g, cellxmax %g. Abort?\n", pkt_ptr->where, x0 * tmin/tstart, cellxmin, cellxmax);
    printout("[warning] outside x boundary of cell %d. y0 %g, cellymin %g, cellymax %g. Abort?\n", pkt_ptr->where, y0 * tmin/tstart, cellymin, cellymax);
    printout("[warning] outside x boundary of cell %d. z0 %g, cellzmin %g, cellzmax %g. Abort?\n", pkt_ptr->where, z0 * tmin/tstart, cellzmin, cellzmax);
    printout("[warning] pkt_ptr->number %d\n", pkt_ptr->number);
    printout("[warning] delta %g\n",  (x0 * tmin/tstart)-cellxmin);
    printout("[warning] dir [%g, %g, %g]\n", pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
    if (((pkt_ptr->dir[0] * CLIGHT) - (x0/tstart)) < 0)
    {
      *snext = pkt_ptr->where - 1;
      pkt_ptr->last_cross = NEG_X;
      return(0);
    }
    else
    {
      not_allowed = POS_X;
    }
      //debuglevel = 2;
  }
  if ((((y0 * tmin) - (tstart * cellymax)) > (-10.*tmin)) && (not_allowed != NEG_Y))
  {
    printout("[warning] outside y boundary of cell %d. x0 %g, cellxmin %g, cellxmax %g. Abort?\n", pkt_ptr->where, x0 * tmin/tstart, cellxmin, cellxmax);
    printout("[warning] outside y boundary of cell %d. y0 %g, cellymin %g, cellymax %g. Abort?\n", pkt_ptr->where, y0 * tmin/tstart, cellymin, cellymax);
    printout("[warning] outside y boundary of cell %d. z0 %g, cellzmin %g, cellzmax %g. Abort?\n", pkt_ptr->where, z0 * tmin/tstart, cellzmin, cellzmax);
    printout("[warning] pkt_ptr->number %d\n", pkt_ptr->number);
    printout("[warning] delta %g\n",cellymax - (y0 * tmin/tstart));
    printout("[warning] dir [%g, %g, %g]\n", pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
    if (((pkt_ptr->dir[1] * CLIGHT)  - (y0/tstart)) > 0)
    {
      *snext = pkt_ptr->where + nxgrid;
      pkt_ptr->last_cross = POS_Y;
      return(0);
    }
    else
    {
      not_allowed = NEG_Y;
    }
      //debuglevel = 2;
  }
  if ((((y0 * tmin) - (tstart * cellymin)) < (10.*tmin)) && (not_allowed != POS_Y))
  {
    printout("[warning] outside y boundary of cell %d. x0 %g, cellxmin %g, cellxmax %g. Abort?\n", pkt_ptr->where, x0 * tmin/tstart, cellxmin, cellxmax);
    printout("[warning] outside y boundary of cell %d. y0 %g, cellymin %g, cellymax %g. Abort?\n", pkt_ptr->where, y0 * tmin/tstart, cellymin, cellymax);
    printout("[warning] outside y boundary of cell %d. z0 %g, cellzmin %g, cellzmax %g. Abort?\n", pkt_ptr->where, z0 * tmin/tstart, cellzmin, cellzmax);
    printout("[warning] pkt_ptr->number %d\n", pkt_ptr->number);
    printout("[warning] delta %g\n",  (y0 * tmin/tstart)-cellymin);
    printout("[warning] dir [%g, %g, %g]\n", pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
    if (((pkt_ptr->dir[1] * CLIGHT)  - (y0/tstart)) < 0)
    {
      *snext = pkt_ptr->where - nxgrid;
      pkt_ptr->last_cross = NEG_Y;
      return(0);
    }
    else
    {
      not_allowed = POS_Y;
    }
      //debuglevel = 2;
  }
  if ((((z0 * tmin) - (tstart * cellzmax)) > (-10.*tmin)) && (not_allowed != NEG_Z))
  {
    printout("[warning] outside z boundary of cell %d. x0 %g, cellxmin %g, cellxmax %g. Abort?\n", pkt_ptr->where, x0 * tmin/tstart, cellxmin, cellxmax);
    printout("[warning] outside z boundary of cell %d. y0 %g, cellymin %g, cellymax %g. Abort?\n", pkt_ptr->where, y0 * tmin/tstart, cellymin, cellymax);
    printout("[warning] outside z boundary of cell %d. z0 %g, cellzmin %g, cellzmax %g. Abort?\n", pkt_ptr->where, z0 * tmin/tstart, cellzmin, cellzmax);
    printout("[warning] pkt_ptr->number %d\n", pkt_ptr->number);
    printout("[warning] delta %g\n",  cellzmax - (z0 * tmin/tstart));
    printout("[warning] dir [%g, %g, %g]\n", pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
    if (((pkt_ptr->dir[2] * CLIGHT)  - (z0/tstart)) > 0)
    {
      *snext = pkt_ptr->where + (nxgrid*nygrid);
      pkt_ptr->last_cross = POS_Z;
      return(0);
    }
    else
    {
      not_allowed = NEG_Z;
    }
      //debuglevel = 2;
  }
  if ((((z0 * tmin) - (tstart * cellzmin)) < (10.*tmin)) && (not_allowed != POS_Z))
  {
    printout("[warning] outside z boundary of cell %d. x0 %g, cellxmin %g, cellxmax %g. Abort?\n", pkt_ptr->where, x0 * tmin/tstart, cellxmin, cellxmax);
    printout("[warning] outside z boundary of cell %d. y0 %g, cellymin %g, cellymax %g. Abort?\n", pkt_ptr->where, y0 * tmin/tstart, cellymin, cellymax);
    printout("[warning] outside z boundary of cell %d. z0 %g, cellzmin %g, cellzmax %g. Abort?\n", pkt_ptr->where, z0 * tmin/tstart, cellzmin, cellzmax);
    printout("[warning] pkt_ptr->number %d\n", pkt_ptr->number);
    printout("[warning] delta %g\n",  (z0 * tmin/tstart)-cellzmin);
    printout("[warning] dir [%g, %g, %g]\n", pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
    if (((pkt_ptr->dir[2] * CLIGHT)  - (z0/tstart)) < 0)
    {
      *snext = pkt_ptr->where - (nxgrid*nygrid);
      pkt_ptr->last_cross = NEG_Z;
      return(0);
    }
    else
    {
      not_allowed = POS_Z;
    }
      //debuglevel = 2;
  }

  #ifdef DEBUG_ON
    if (debuglevel ==2 )
    {
      printout("pkt_ptr->number %d\n", pkt_ptr->number);
      printout("delta1x %g delta2x %g\n",  (x0 * tmin/tstart)-cellxmin, cellxmax - (x0 * tmin/tstart));
      printout("delta1y %g delta2y %g\n",  (y0 * tmin/tstart)-cellymin, cellymax - (y0 * tmin/tstart));
      printout("delta1z %g delta2z %g\n",  (z0 * tmin/tstart)-cellzmin, cellzmax - (z0 * tmin/tstart));
      printout("dir [%g, %g, %g]\n", pkt_ptr->dir[0],pkt_ptr->dir[1],pkt_ptr->dir[2]);
    }
  #endif

  double tx_plus = ((x0 - (vx * tstart))/(cellxmax - (vx * tmin)) * tmin) - tstart;
  double tx_minus = ((x0 - (vx * tstart))/(cellxmin - (vx * tmin)) * tmin) - tstart;

  double ty_plus = ((y0 - (vy * tstart))/(cellymax - (vy * tmin)) * tmin) - tstart;
  double ty_minus = ((y0 - (vy * tstart))/(cellymin - (vy * tmin)) * tmin) - tstart;

  double tz_plus = ((z0 - (vz * tstart))/(cellzmax - (vz * tmin)) * tmin) - tstart;
  double tz_minus = ((z0 - (vz * tstart))/(cellzmin - (vz * tmin)) * tmin) - tstart;

  /** We now need to identify the shortest +ve time - that's the one we want. */
  int choice = 0;         ///just a control variable to
  double time = 1.e99;
  //close = 1.e99;
  //printout("bondary.c check value of not_allowed = %d\n",not_allowed);
  if ((tx_plus > 0) && (tx_plus < time) && (not_allowed != NEG_X))
  {
    choice = 1;
    time = tx_plus;
    if (cell[pkt_ptr->where].xyz[0] == (nxgrid - 1))
    {
      *snext = -99;
    }
    else
    {
      *snext = pkt_ptr->where + 1;
      pkt_ptr->last_cross = POS_X;
    }
  }

  //close_try = fabs(CLIGHT_PROP * tx_plus * pkt_ptr->dir[0]);
  //if (close_try < close)
  // {
  //    close = close_try;
  //  }

  if ((tx_minus > 0) && (tx_minus < time) && (not_allowed != POS_X))
  {
    choice = 2;
    time = tx_minus;
    if (cell[pkt_ptr->where].xyz[0] == 0)
    {
      *snext = -99;
    }
    else
    {
      *snext = pkt_ptr->where - 1;
      pkt_ptr->last_cross = NEG_X;
    }
  }

  //close_try = fabs(CLIGHT_PROP * tx_minus * pkt_ptr->dir[0]);
  //if (close_try < close)
  //  {
  //    close = close_try;
  //  }

  if ((ty_plus > 0) && (ty_plus < time) && (not_allowed != NEG_Y))
  {
    choice = 3;
    time = ty_plus;
    if (cell[pkt_ptr->where].xyz[1] == (nygrid - 1))
    {
      *snext = -99;
    }
    else
    {
      *snext = pkt_ptr->where + nxgrid;
      pkt_ptr->last_cross = POS_Y;
    }
  }

  //close_try = fabs(CLIGHT_PROP * ty_plus * pkt_ptr->dir[1]);
  //if (close_try < close)
  //  {
  //    close = close_try;
  //  }

  if ((ty_minus > 0) && (ty_minus < time) && (not_allowed != POS_Y))
  {
    choice = 4;
    time = ty_minus;
    if (cell[pkt_ptr->where].xyz[1] == 0)
    {
      *snext = -99;
    }
    else
    {
      *snext = pkt_ptr->where - nxgrid;
      pkt_ptr->last_cross = NEG_Y;
    }
  }

  //close_try = fabs(CLIGHT_PROP * ty_minus * pkt_ptr->dir[1]);
  //if (close_try < close)
   // {
  //    close = close_try;
  // }

  if ((tz_plus > 0) && (tz_plus < time) && (not_allowed != NEG_Z))
  {
    choice = 5;
    time = tz_plus;
    if (cell[pkt_ptr->where].xyz[2] == (nzgrid - 1))
    {
      *snext = -99;
    }
    else
    {
      *snext = pkt_ptr->where + (nxgrid*nygrid);
      pkt_ptr->last_cross = POS_Z;
    }
  }

  //close_try = fabs(CLIGHT_PROP * tz_plus * pkt_ptr->dir[2]);
  //if (close_try < close)
  //  {
  //    close = close_try;
  //  }

  if ((tz_minus > 0) && (tz_minus < time) && (not_allowed != POS_Z))
  {
    choice = 6;
    time = tz_minus;
    if (cell[pkt_ptr->where].xyz[2] == 0)
    {
      *snext = -99;
    }
    else
    {
      *snext = pkt_ptr->where - (nxgrid*nygrid);
      pkt_ptr->last_cross = NEG_Z;
    }
  }

  //close_try = fabs(CLIGHT_PROP * tz_minus * pkt_ptr->dir[2]);
  //if (close_try < close)
  //  {
  //    close = close_try;
  //  }


  if (choice == 0)
  {
    printout("Something wrong in boundary crossing - didn't find anything.\n");
    printout("tx_plus %g tx_minus %g \n", tx_plus, tx_minus);
    printout("ty_plus %g ty_minus %g \n", ty_plus, ty_minus);
    printout("tz_plus %g tz_minus %g \n", tz_plus, tz_minus);
    printout("x0 %g y0 %g z0 %g \n", x0, y0, z0);
    printout("vx %g vy %g vy %g \n", vx, vy, vz);
    printout("tstart %g\n", tstart);

    exit(0);
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


///****************************************************************************
/// Routine to take a packet across a boundary.
int change_cell(pkt_ptr, snext, end_packet, t_current)
    PKT *pkt_ptr;
     double t_current;
     int snext;
     int *end_packet;
{
  //void copy_populations_to_phixslist();
  int search_cellhistory(int cellnumber);
  int find_farthestcell(int cellnumber);
  void update_cell(int cellnumber);
  void calculate_kappa_rpkt_cont(PKT *pkt_ptr, double t_current);
  void determine_kpkt_cuts(int cellnumber);
  //void calculate_kpkt_rates(int cellnumber);
  void calculate_levelpops(int cellnumber);
  //int element, ion, level;
  int oldpos,mgi,old_mgi;

  #ifdef DEBUG_ON
    if (debuglevel == 2)
    {
      printout("[debug] cellnumber %d nne %g\n",pkt_ptr->where,get_nne(pkt_ptr->where));
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
    nesc += 1;
    *end_packet = 1;
    return(0);
  }
  else
  {
    /** Just need to update "where".*/
    oldpos = pkt_ptr->where;
    old_mgi = cell[pkt_ptr->where].modelgridindex;
    pkt_ptr->where = snext;
    mgi = cell[pkt_ptr->where].modelgridindex;
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
    cellcrossings += 1;

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

    return 0;
  }
}


///****************************************************************************
/// Routine to return which grid cell the packet is in.
int locate(pkt_ptr, t_current)
    PKT *pkt_ptr;
     double t_current;
{
  /* Cheap and nasty version for now - assume a uniform grid. */
  int xx = (pkt_ptr->pos[0] - (cell[0].pos_init[0]*t_current/tmin)) / (wid_init*t_current/tmin);
  int yy = (pkt_ptr->pos[1] - (cell[0].pos_init[1]*t_current/tmin)) / (wid_init*t_current/tmin);
  int zz = (pkt_ptr->pos[2] - (cell[0].pos_init[2]*t_current/tmin)) / (wid_init*t_current/tmin);

  return xx + (nxgrid * yy) + (nxgrid * nygrid * zz);
}
