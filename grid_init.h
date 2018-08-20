#ifndef GRIDINIT_H
#define GRIDINIT_H

#include "sn3d.h"

void grid_init(int my_rank);
double get_cellradialpos(int cellindex);
float get_modelinitradioabund(int modelgridindex, enum radionuclides nuclide_type);
void set_modelinitradioabund(int modelgridindex, enum radionuclides nuclide_type, float abund);
float get_stable_abund(int mgi, int anumber);
int get_numassociatedcells(int modelgridindex);


inline double wid_init(const int cellindex)
// for a uniform grid this is the extent along the x,y,z coordinate (x_2 - x_1, etc.)
// for spherical grid this is the radial extent (r_outer - r_inner)
// these values are for time tmin
{
  switch (grid_type)
  {
    case GRID_SPHERICAL1D:
    {
      const int modelgridindex = cell[cellindex].modelgridindex;
      const double v_inner = modelgridindex > 0 ? vout_model[modelgridindex - 1] : 0.;
      return (vout_model[modelgridindex] - v_inner) * tmin;
    }

    default:
      return 2 * coordmax[0] / ncoordgrid[0];
  }
}


inline double vol_init_modelcell(const int modelgridindex)
// return the model cell volume at tmin
// for a uniform cubic grid this is constant
{
  switch (grid_type)
  {
    case GRID_SPHERICAL1D:
      return 4./3. * PI * (pow(tmin * vout_model[modelgridindex], 3) - pow(tmin * (modelgridindex > 0 ? vout_model[modelgridindex - 1] : 0.), 3));

    default:
    {
      const int assoc_cells = get_numassociatedcells(modelgridindex);
      return (wid_init(0) * wid_init(0) * wid_init(0)) * assoc_cells;
    }
  }
}


inline double vol_init_gridcell(const int cellindex)
// return the propagation cell volume at tmin
// for a spherical grid, the cell index is required (and should be equivalent to a modelgridindex)
{
  switch (grid_type)
  {
    case GRID_SPHERICAL1D:
    {
      const int mgi = cell[cellindex].modelgridindex;
      return vol_init_modelcell(mgi);
    }

    default:
      return (wid_init(0) * wid_init(0) * wid_init(0));
  }
}


inline double get_cellcoordmin(const int cellindex, const int axis)
// get the minimum value of a coordinate at tmin (xyz or radial coords) of a propagation cell
// e.g., the minimum x position in xyz coords, or the minimum radius
{
  return cell[cellindex].pos_init[axis];
  // return - coordmax[axis] + (2 * get_cellcoordpointnum(cellindex, axis) * coordmax[axis] / ncoordgrid[axis]);
}


inline int get_coordcellindexincrement(const int axis)
// how much do we change the cellindex to move along a coordinately axis (e.g., the x, y, z directions, or r direction)
{
  switch (grid_type)
  {
    case GRID_SPHERICAL1D:
      return 1;

    default:
      switch (axis)
      {
        case 0:
          return ncoordgrid[2] * ncoordgrid[1];

        case 1:
          return ncoordgrid[2];

        case 2:
          return 1;

        default:
          printout("invalid coordinate index %d", axis);
          abort();
      }
  }
}

inline int get_cellcoordpointnum(const int cellindex, const int axis)
// convert a cell index number into an integer (x,y,z or r) coordinate index from 0 to ncoordgrid[axis]
{
  // return cell[cellindex].nxyz[axis];

  switch (grid_type)
  {
    case GRID_SPHERICAL1D:
      return cellindex;

    default:
      switch (axis)
      {
        // increment x first, then y, then z
        case 0:
          return cellindex % ncoordgrid[0];

        case 1:
          return (cellindex / ncoordgrid[0]) % ncoordgrid[1];

        case 2:
          return (cellindex / (ncoordgrid[0] * ncoordgrid[1])) % ncoordgrid[2];

        default:
          printout("invalid coordinate index %d", axis);
          abort();
      }
  }
}


inline int get_ngriddimensions(void)
{
  return (grid_type == GRID_SPHERICAL1D) ? 1 : 3;
}


inline float get_rhoinit(int modelgridindex)
{
  return modelgrid[modelgridindex].rhoinit;
}

inline float get_rho(int modelgridindex)
{
  return modelgrid[modelgridindex].rho;
}

inline float get_nne(int modelgridindex)
{
  return modelgrid[modelgridindex].nne;
}

inline float get_nnetot(int modelgridindex)
{
  return modelgrid[modelgridindex].nnetot;
}

// the abundances referred to below are initial abundances
inline float get_ffegrp(int modelgridindex)
{
  return modelgrid[modelgridindex].ffegrp;
}

inline float get_kappagrey(int modelgridindex)
{
  return modelgrid[modelgridindex].kappagrey;
}

inline float get_Te(int modelgridindex)
{
  return modelgrid[modelgridindex].Te;
}

inline float get_TR(int modelgridindex)
{
  return modelgrid[modelgridindex].TR;
}

inline float get_TJ(int modelgridindex)
{
  return modelgrid[modelgridindex].TJ;
}

inline float get_W(int modelgridindex)
{
  return modelgrid[modelgridindex].W;
}

inline void set_rhoinit(int modelgridindex, float x)
{
  modelgrid[modelgridindex].rhoinit = x;
}

inline void set_rho(int modelgridindex, float x)
{
  modelgrid[modelgridindex].rho = x;
}

inline void set_nne(int modelgridindex, float x)
{
  modelgrid[modelgridindex].nne = x;
}

inline void set_nnetot(int modelgridindex, float x)
{
  modelgrid[modelgridindex].nnetot = x;
}

inline void set_ffegrp(int modelgridindex, float x)
{
  modelgrid[modelgridindex].ffegrp = x;
}

inline void set_kappagrey(int modelgridindex, float kappagrey)
{
  modelgrid[modelgridindex].kappagrey = kappagrey;
}

inline void set_Te(int modelgridindex, float Te)
{
  modelgrid[modelgridindex].Te = Te;
}

inline void set_TR(int modelgridindex, float TR)
{
  modelgrid[modelgridindex].TR = TR;
}

inline void set_TJ(int modelgridindex, float TJ)
{
  modelgrid[modelgridindex].TJ = TJ;
}

inline void set_W(int modelgridindex, float W)
{
  modelgrid[modelgridindex].W = W;
}


#endif //GRIDINIT_H
