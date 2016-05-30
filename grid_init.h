#ifndef GRIDINIT_H
#define GRIDINIT_H

#include "sn3d.h"

int grid_init(void);
void allocate_compositiondata(int cellnumber);
void allocate_cooling(int modelgridindex);

/// Routine for getting the initial cell volume.
static inline
double vol_init() //CELL *restrict grid_ptr
{
  return (wid_init * wid_init * wid_init);
}

static inline
float get_rhoinit(int modelgridindex)
{
  return modelgrid[modelgridindex].rhoinit;
}

static inline
float get_rho(int modelgridindex)
{
  return modelgrid[modelgridindex].rho;
}

static inline
float get_nne(int modelgridindex)
{
  return modelgrid[modelgridindex].nne;
}

static inline
float get_nnetot(int modelgridindex)
{
  return modelgrid[modelgridindex].nnetot;
}

static inline
float get_fni(int modelgridindex)
{
  return modelgrid[modelgridindex].fni;
}

static inline
float get_fco(int modelgridindex)
{
  return modelgrid[modelgridindex].fco;
}

static inline
float get_f52fe(int modelgridindex)
{
  return modelgrid[modelgridindex].f52fe;
}

static inline
float get_f48cr(int modelgridindex)
{
  return modelgrid[modelgridindex].f48cr;
}

static inline
float get_ffe(int modelgridindex)
{
  return modelgrid[modelgridindex].ffe;
}

static inline
float get_fnistable(int modelgridindex)
{
  return modelgrid[modelgridindex].fnistable;
}

static inline
float get_fcostable(int modelgridindex)
{
  return modelgrid[modelgridindex].fcostable;
}

static inline
float get_ffestable(int modelgridindex)
{
  return modelgrid[modelgridindex].ffestable;
}

static inline
float get_fmnstable(int modelgridindex)
{
  return modelgrid[modelgridindex].fmnstable;
}

static inline
float get_fcrstable(int modelgridindex)
{
  return modelgrid[modelgridindex].fcrstable;
}

static inline
float get_fvstable(int modelgridindex)
{
  return modelgrid[modelgridindex].fvstable;
}

static inline
float get_ftistable(int modelgridindex)
{
  return modelgrid[modelgridindex].ftistable;
}

static inline
float get_kappagrey(int modelgridindex)
{
  return modelgrid[modelgridindex].kappagrey;
}

static inline
float get_Te(int modelgridindex)
{
  return modelgrid[modelgridindex].Te;
}

static inline
float get_TR(int modelgridindex)
{
  return modelgrid[modelgridindex].TR;
}

static inline
float get_TJ(int modelgridindex)
{
  return modelgrid[modelgridindex].TJ;
}

static inline
float get_W(int modelgridindex)
{
  return modelgrid[modelgridindex].W;
}

static inline
void set_rhoinit(int modelgridindex, float x)
{
  modelgrid[modelgridindex].rhoinit = x;
}

static inline
void set_rho(int modelgridindex, float x)
{
  modelgrid[modelgridindex].rho = x;
}

static inline
void set_nne(int modelgridindex, float x)
{
  modelgrid[modelgridindex].nne = x;
}

static inline
void set_nnetot(int modelgridindex, float x)
{
  modelgrid[modelgridindex].nnetot = x;
}

static inline
void set_fni(int modelgridindex, float x)
{
  modelgrid[modelgridindex].fni = x;
}

static inline
void set_fco(int modelgridindex, float x)
{
  modelgrid[modelgridindex].fco = x;
}

static inline
void set_f48cr(int modelgridindex, float x)
{
  modelgrid[modelgridindex].f48cr = x;
}

static inline
void set_f52fe(int modelgridindex, float x)
{
  modelgrid[modelgridindex].f52fe = x;
}

static inline
void set_ffe(int modelgridindex, float x)
{
  modelgrid[modelgridindex].ffe = x;
}

static inline
void set_fnistable(int modelgridindex, float x)
{
  if (x >= 0)
  {
    modelgrid[modelgridindex].fnistable = x;
  }
  else
  {
    //printout("Setting fnistable to 0.0 to avoid negative.\n");
    modelgrid[modelgridindex].fnistable = 0.0;
  }
}

static inline
void set_fcostable(int modelgridindex, float x)
{
  if (x >= 0)
  {
    modelgrid[modelgridindex].fcostable = x;
  }
  else
  {
    //printout("Setting fcostable to 0.0 to avoid negative.\n");
    modelgrid[modelgridindex].fcostable = 0.0;
  }
}

static inline
void set_ffestable(int modelgridindex, float x)
{
  if (x >= 0)
  {
    modelgrid[modelgridindex].ffestable = x;
  }
  else
  {
    //printout("Setting ffestable to 0.0 to avoid negative.\n");
    modelgrid[modelgridindex].ffestable = 0.0;
  }
}

static inline
void set_fmnstable(int modelgridindex, float x)
{
  if (x >= 0)
  {
    modelgrid[modelgridindex].fmnstable = x;
  }
  else
  {
    //printout("Setting ffestable to 0.0 to avoid negative.\n");
    modelgrid[modelgridindex].fmnstable = 0.0;
  }
}

static inline
void set_fcrstable(int modelgridindex, float x)
{
  if (x >= 0)
  {
    modelgrid[modelgridindex].fcrstable = x;
  }
  else
  {
    //printout("Setting fcrstable to 0.0 to avoid negative.\n");
    modelgrid[modelgridindex].fcrstable = 0.0;
  }
}

static inline
void set_fvstable(int modelgridindex, float x)
{
  if (x >= 0)
  {
    modelgrid[modelgridindex].fvstable = x;
  }
  else
  {
    //printout("Setting fcrstable to 0.0 to avoid negative.\n");
    modelgrid[modelgridindex].fvstable = 0.0;
  }
}

static inline
void set_ftistable(int modelgridindex, float x)
{
  if (x >= 0)
  {
    modelgrid[modelgridindex].ftistable = x;
  }
  else
  {
    //printout("Setting fcrstable to 0.0 to avoid negative.\n");
    modelgrid[modelgridindex].ftistable = 0.0;
  }
}

static inline
void set_kappagrey(int modelgridindex, float x)
{
  modelgrid[modelgridindex].kappagrey = x;
}

static inline
void set_Te(int modelgridindex, float x)
{
  modelgrid[modelgridindex].Te = x;
}

static inline
void set_TR(int modelgridindex, float x)
{
  modelgrid[modelgridindex].TR = x;
}

static inline
void set_TJ(int modelgridindex, float x)
{
  modelgrid[modelgridindex].TJ = x;
}

static inline
void set_W(int modelgridindex, float x)
{
  modelgrid[modelgridindex].W = x;
}


#endif //GRIDINIT_H
