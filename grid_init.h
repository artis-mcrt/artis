#ifndef GRIDINIT_H
#define GRIDINIT_H

#include <cassert>
#include "sn3d.h"

void grid_init(int my_rank);
double get_cellradialpos(int cellindex);
float get_modelinitradioabund(int modelgridindex, enum radionuclides nuclide_type);
void set_modelinitradioabund(int modelgridindex, enum radionuclides nuclide_type, float abund);
float get_stable_abund(int mgi, int anumber);
int get_numassociatedcells(int modelgridindex);
enum model_types get_model_type(void);
void set_model_type(enum model_types model_type_value);


inline double wid_init(const int cellindex)
// for a uniform grid this is the extent along the x,y,z coordinate (x_2 - x_1, etc.)
// for spherical grid this is the radial extent (r_outer - r_inner)
// these values are for time globals::tmin
{
  switch (globals::grid_type)
  {
    case GRID_SPHERICAL1D:
    {
      const int modelgridindex = globals::cell[cellindex].modelgridindex;
      const double v_inner = modelgridindex > 0 ? globals::vout_model[modelgridindex - 1] : 0.;
      return (globals::vout_model[modelgridindex] - v_inner) * globals::tmin;
    }

    default:
      return 2 * globals::coordmax[0] / globals::ncoordgrid[0];
  }
}


inline double vol_init_modelcell(const int modelgridindex)
// return the model cell volume at globals::tmin
// for a uniform cubic grid this is constant
{
  switch (globals::grid_type)
  {
    case GRID_SPHERICAL1D:
      return 4./3. * PI * (pow(globals::tmin * globals::vout_model[modelgridindex], 3) - pow(globals::tmin * (modelgridindex > 0 ? globals::vout_model[modelgridindex - 1] : 0.), 3));

    default:
    {
      const int assoc_cells = get_numassociatedcells(modelgridindex);
      return (wid_init(0) * wid_init(0) * wid_init(0)) * assoc_cells;
    }
  }
}


inline double vol_init_gridcell(const int cellindex)
// return the propagation cell volume at globals::tmin
// for a spherical grid, the cell index is required (and should be equivalent to a modelgridindex)
{
  switch (globals::grid_type)
  {
    case GRID_SPHERICAL1D:
    {
      const int mgi = globals::cell[cellindex].modelgridindex;
      return vol_init_modelcell(mgi);
    }

    default:
      return (wid_init(0) * wid_init(0) * wid_init(0));
  }
}


inline double get_cellcoordmin(const int cellindex, const int axis)
// get the minimum value of a coordinate at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the minimum x position in xyz coords, or the minimum radius
{
  return globals::cell[cellindex].pos_init[axis];
  // return - coordmax[axis] + (2 * get_cellcoordpointnum(cellindex, axis) * coordmax[axis] / globals::ncoordgrid[axis]);
}


inline int get_coordcellindexincrement(const int axis)
// how much do we change the cellindex to move along a coordinately axis (e.g., the x, y, z directions, or r direction)
{
  switch (globals::grid_type)
  {
    case GRID_SPHERICAL1D:
      return 1;

    default:
      switch (axis)
      {
        case 0:
          return 1;

        case 1:
          return globals::ncoordgrid[0];

        case 2:
          return globals::ncoordgrid[0] * globals::ncoordgrid[1];

        default:
          printout("invalid coordinate index %d", axis);
          abort();
      }
  }
}


inline int get_cellcoordpointnum(const int cellindex, const int axis)
// convert a cell index number into an integer (x,y,z or r) coordinate index from 0 to globals::ncoordgrid[axis]
{
  // return globals::cell[cellindex].nxyz[axis];

  switch (globals::grid_type)
  {
    case GRID_SPHERICAL1D:
      return cellindex;

    default:
      switch (axis)
      {
        // increment x first, then y, then z
        case 0:
          return cellindex % globals::ncoordgrid[0];

        case 1:
          return (cellindex / globals::ncoordgrid[0]) % globals::ncoordgrid[1];

        case 2:
          return (cellindex / (globals::ncoordgrid[0] * globals::ncoordgrid[1])) % globals::ncoordgrid[2];

        default:
          printout("invalid coordinate index %d", axis);
          abort();
      }
  }
}



inline int get_ngriddimensions(void)
{
  return (globals::grid_type == GRID_SPHERICAL1D) ? 1 : 3;
}


inline float get_rhoinit(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].rhoinit;
}

inline float get_rho(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].rho;
}

inline float get_nne(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].nne;
}

inline float get_nnetot(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].nnetot;
}

// the abundances referred to below are initial abundances
inline float get_ffegrp(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].ffegrp;
}

inline float get_kappagrey(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].kappagrey;
}

inline float get_Te(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].Te;
}

inline float get_TR(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].TR;
}

inline float get_TJ(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].TJ;
}

inline float get_W(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].W;
}

inline void set_rhoinit(int modelgridindex, float x)
{
  globals::modelgrid[modelgridindex].rhoinit = x;
}

inline void set_rho(int modelgridindex, float x)
{
  globals::modelgrid[modelgridindex].rho = x;
}

inline void set_nne(int modelgridindex, float x)
{
  globals::modelgrid[modelgridindex].nne = x;
}

inline void set_nnetot(int modelgridindex, float x)
{
  globals::modelgrid[modelgridindex].nnetot = x;
}

inline void set_ffegrp(int modelgridindex, float x)
{
  assert(x >= 0);
  assert(x <= 1.);
  globals::modelgrid[modelgridindex].ffegrp = x;
}

inline void set_kappagrey(int modelgridindex, float kappagrey)
{
  globals::modelgrid[modelgridindex].kappagrey = kappagrey;
}

inline void set_Te(int modelgridindex, float Te)
{
  globals::modelgrid[modelgridindex].Te = Te;
}

inline void set_TR(int modelgridindex, float TR)
{
  globals::modelgrid[modelgridindex].TR = TR;
}

inline void set_TJ(int modelgridindex, float TJ)
{
  globals::modelgrid[modelgridindex].TJ = TJ;
}

inline void set_W(int modelgridindex, float W)
{
  globals::modelgrid[modelgridindex].W = W;
}


#endif //GRIDINIT_H
