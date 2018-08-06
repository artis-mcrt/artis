#ifndef GRIDINIT_H
#define GRIDINIT_H

#include "sn3d.h"

void grid_init(int my_rank);
int get_cellcoordpointnum(int cellindex, int axis);
double get_cellradialpos(int cellindex);
float get_modelradioabund(int modelgridindex, enum radionuclides nuclide_type);
void set_modelradioabund(int modelgridindex, enum radionuclides nuclide_type, float abund);


inline double wid_init(const int cellindex)
{
  if (grid_type == GRID_SPHERICAL1D)
  {
    const double v_inner = cellindex > 0 ? vout_model[cellindex - 1] : 0.;
    return (vout_model[cellindex] - v_inner) * tmin;
  }
  else
  {
    return 2 * coordmax[0] / ncoordgrid[0];
  }
}

/// Routine for getting the initial cell volume.
inline double vol_init(const int cellindex)//(const CELL *restrict const grid_ptr)
{
  switch (grid_type)
  {
    case GRID_UNIFORM:
      return (wid_init(0) * wid_init(0) * wid_init(0));

    case GRID_SPHERICAL1D:
      return 4./3. * PI * (pow(tmin * vout_model[cellindex], 3) - pow(tmin * (cellindex > 0 ? vout_model[cellindex - 1] : 0.), 3));
  }
  abort();
}

inline double get_cellcoordmin(const int cellindex, const int axis)
// get the minimum position along each axis at tmin (xyz or radial coords)
{
  return cell[cellindex].pos_init[axis];
  // return - coordmax[axis] + (2 * get_cellcoordpointnum(cellindex, axis) * coordmax[axis] / ncoordgrid[axis]);
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

inline float get_fnistable(int modelgridindex)
{
  return modelgrid[modelgridindex].fnistable;
}

inline float get_fcostable(int modelgridindex)
{
  return modelgrid[modelgridindex].fcostable;
}

inline float get_ffestable(int modelgridindex)
{
  return modelgrid[modelgridindex].ffestable;
}

inline float get_fmnstable(int modelgridindex)
{
  return modelgrid[modelgridindex].fmnstable;
}

inline float get_fcrstable(int modelgridindex)
{
  return modelgrid[modelgridindex].fcrstable;
}

inline float get_fvstable(int modelgridindex)
{
  return modelgrid[modelgridindex].fvstable;
}

inline float get_ftistable(int modelgridindex)
{
  return modelgrid[modelgridindex].ftistable;
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
