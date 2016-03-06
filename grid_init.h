#ifndef GRIDINIT_H
#define GRIDINIT_H

int grid_init(void);
int uniform_grid_setup(void);
int density_1d_read(void);
int density_2d_read(void);
int density_3d_read(void);
void allocate_compositiondata(int cellnumber);
void allocate_cooling(int modelgridindex);
void abundances_setup(void);
void abundances_3d_read(void);
void abundances_1d_read(void);
void assign_temperature(void);

/// Routine for getting the initial cell volume.
inline
double vol_init(CELL *grid_ptr)
{
  return (wid_init * wid_init * wid_init);
}

inline
float get_rhoinit(int modelgridindex)
{
  return modelgrid[modelgridindex].rhoinit;
}

inline
float get_rho(int modelgridindex)
{
  return modelgrid[modelgridindex].rho;
}

inline
float get_nne(int modelgridindex)
{
  return modelgrid[modelgridindex].nne;
}

inline
float get_nnetot(int modelgridindex)
{
  return modelgrid[modelgridindex].nnetot;
}

inline
float get_fni(int modelgridindex)
{
  return modelgrid[modelgridindex].fni;
}

inline
float get_fco(int modelgridindex)
{
  return modelgrid[modelgridindex].fco;
}

inline
float get_f52fe(int modelgridindex)
{
  return modelgrid[modelgridindex].f52fe;
}

inline
float get_f48cr(int modelgridindex)
{
  return modelgrid[modelgridindex].f48cr;
}

inline
float get_ffe(int modelgridindex)
{
  return modelgrid[modelgridindex].ffe;
}

inline
float get_fnistable(int modelgridindex)
{
  return modelgrid[modelgridindex].fnistable;
}

inline
float get_fcostable(int modelgridindex)
{
  return modelgrid[modelgridindex].fcostable;
}

inline
float get_ffestable(int modelgridindex)
{
  return modelgrid[modelgridindex].ffestable;
}

inline
float get_fmnstable(int modelgridindex)
{
  return modelgrid[modelgridindex].fmnstable;
}

inline
float get_fcrstable(int modelgridindex)
{
  return modelgrid[modelgridindex].fcrstable;
}

inline
float get_fvstable(int modelgridindex)
{
  return modelgrid[modelgridindex].fvstable;
}

inline
float get_ftistable(int modelgridindex)
{
  return modelgrid[modelgridindex].ftistable;
}

inline
float get_kappagrey(int modelgridindex)
{
  return modelgrid[modelgridindex].kappagrey;
}

inline
float get_Te(int modelgridindex)
{
  return modelgrid[modelgridindex].Te;
}

inline
float get_TR(int modelgridindex)
{
  return modelgrid[modelgridindex].TR;
}

inline
float get_TJ(int modelgridindex)
{
  return modelgrid[modelgridindex].TJ;
}

inline
float get_W(int modelgridindex)
{
  return modelgrid[modelgridindex].W;
}

inline
void set_rhoinit(int modelgridindex, float x)
{
  modelgrid[modelgridindex].rhoinit = x;
}

inline
void set_rho(int modelgridindex, float x)
{
  modelgrid[modelgridindex].rho = x;
}

inline
void set_nne(int modelgridindex, float x)
{
  modelgrid[modelgridindex].nne = x;
}

inline
void set_nnetot(int modelgridindex, float x)
{
  modelgrid[modelgridindex].nnetot = x;
}

inline
void set_fni(int modelgridindex, float x)
{
  modelgrid[modelgridindex].fni = x;
}

inline
void set_fco(int modelgridindex, float x)
{
  modelgrid[modelgridindex].fco = x;
}

inline
void set_f48cr(int modelgridindex, float x)
{
  modelgrid[modelgridindex].f48cr = x;
}

inline
void set_f52fe(int modelgridindex, float x)
{
  modelgrid[modelgridindex].f52fe = x;
}

inline
void set_ffe(int modelgridindex, float x)
{
  modelgrid[modelgridindex].ffe = x;
}

inline
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

inline
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

inline
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

inline
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

inline
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

inline
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

inline
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

inline
void set_kappagrey(int modelgridindex, float x)
{
  modelgrid[modelgridindex].kappagrey = x;
}

inline
void set_Te(int modelgridindex, float x)
{
  modelgrid[modelgridindex].Te = x;
}

inline
void set_TR(int modelgridindex, float x)
{
  modelgrid[modelgridindex].TR = x;
}

inline
void set_TJ(int modelgridindex, float x)
{
  modelgrid[modelgridindex].TJ = x;
}

inline
void set_W(int modelgridindex, float x)
{
  modelgrid[modelgridindex].W = x;
}


#endif //GRIDINIT_H
