#include "assert.h"
#include "sn3d.h"
#include "atomic.h"
#include "grid_init.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "radfield.h"
#include "rpkt.h"
#include "vectors.h"


extern inline double wid_init(int cellindex);
extern inline double vol_init_model(int modelgridindex);
extern inline double vol_init_grid(int cellindex);
extern inline float get_rhoinit(int modelgridindex);
extern inline float get_rho(int modelgridindex);
extern inline float get_nne(int modelgridindex);
extern inline float get_nnetot(int modelgridindex);
extern inline float get_ffegrp(int modelgridindex);
extern inline float get_fnistable(int modelgridindex);
extern inline float get_fcostable(int modelgridindex);
extern inline float get_ffestable(int modelgridindex);
extern inline float get_fmnstable(int modelgridindex);
extern inline float get_fcrstable(int modelgridindex);
extern inline float get_fvstable(int modelgridindex);
extern inline float get_ftistable(int modelgridindex);
extern inline float get_kappagrey(int modelgridindex);
extern inline float get_Te(int modelgridindex);
extern inline float get_TR(int modelgridindex);
extern inline float get_TJ(int modelgridindex);
extern inline float get_W(int modelgridindex);
extern inline void set_rhoinit(int modelgridindex, float x);
extern inline void set_rho(int modelgridindex, float x);
extern inline void set_nne(int modelgridindex, float x);
extern inline void set_nnetot(int modelgridindex, float x);
extern inline void set_ffegrp(int modelgridindex, float x);
extern inline void set_kappagrey(int modelgridindex, float x);
extern inline void set_Te(int modelgridindex, float x);
extern inline void set_TR(int modelgridindex, float x);
extern inline void set_TJ(int modelgridindex, float x);
extern inline void set_W(int modelgridindex, float x);


static long mem_usage_nltepops = 0;

static int mg_associated_cells[MMODELGRID + 1];


int get_numassociatedcells(const int modelgridindex)
// number of propagation cells associated with each modelgrid cell
{
  return mg_associated_cells[modelgridindex];
}


float get_modelradioabund(const int modelgridindex, const enum radionuclides nuclide_type)
{
  // this function replaces get_f56ni(mgi), get_fco56(mgi), etc.
  if (model_type == RHO_UNIFORM)
  {
    // Ni 56 must be handled separate based on a grid pointer instead of a modelgridindex
    assert(nuclide_type != NUCLIDE_NI56);
    return 0.;
  }

  if (nuclide_type >= RADIONUCLIDE_COUNT)
  {
    printout("get_initfracradio called with invalid nuclide_type type %d\n", nuclide_type);
    abort();
  }

  return modelgrid[modelgridindex].fradionuclides[nuclide_type];
}


void set_modelradioabund(const int modelgridindex, const enum radionuclides nuclide_type, const float abund)
{
  // this function replaces set_f56ni(mgi), set_fco56(mgi), etc.

  if ((nuclide_type == FAKE_GAM_LINE_ID && abund > 0) || nuclide_type >= RADIONUCLIDE_COUNT)
  {
    printout("set_initfracradio called with invalid nuclide_type type %d\n", nuclide_type);
    abort();
  }

  modelgrid[modelgridindex].fradionuclides[nuclide_type] = abund;
}

/// Routine for doing a uniform density grid.
/*int uniform_density_setup ()
{
  int n;
  double vec_len();
  double dcen[3];
  double f56ni(CELL *grid_ptr);
  void allocate_compositiondata(int cellnumber);
  double opcase2_normal, opcase3_sum, check1, check2;
  int empty_cells;

  //Calculate the critical opacity at which opacity_case 3 switches from a
  //regime proportional to the density to a regime independent of the density
  //This is done by solving for tau_sobolev == 1
  //tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/MNI56 * 3000e-8 * time_step[m].mid;
  rho_crit = ME*CLIGHT*MNI56 / (PI*QE*QE * rho_crit_para * 3000e-8 * tmin);
  printout("grid_init: rho_crit = %g\n", rho_crit);

  rho_sum = 0.0;
  fe_sum = 0.0;
  empty_cells = 0;
  check1 = check2 = 0.0;
  opcase3_sum = 0;

  //MK: first loop to initialize basic variables
  for (n=0;n < ngrid; n++)
  {
    dcen[0] = cell[n].pos_init[0] + (0.5*wid_init);
    dcen[1] = cell[n].pos_init[1] + (0.5*wid_init);
    dcen[2] = cell[n].pos_init[2] + (0.5*wid_init);


    if (vec_len(dcen) < rmax)
    {
      cell[n].rho_init = cell[n].rho = 3 * mtot / 4 / PI / rmax /rmax /rmax;
      cell[n].f_ni = cell[n].f_fe = f56ni(&cell[n]);
      cell[n].f_co = 0.;
      allocate_compositiondata(n);
    //printout("[debug] uniform_density_setup: cell[%d].rho_init: %g\n",n,cell[n].rho_init);
    }
     //MK: added missing else branch for correct initialization
    else
    {
      cell[n].rho_init = cell[n].rho = 0.0;
      cell[n].kappa_grey = 0.0;
      empty_cells++;
    }
    rho_sum += cell[n].rho_init;
    fe_sum += f56ni(&cell[n]);  //MK: use for this model fni as values for iron group mass fractions
  }


  if (opacity_case == 3)
  {
    for (n = 0; n < ngrid; n++)
    {
      if (cell[n].rho_init > 0)
      {
        if (cell[n].rho_init > rho_crit)
        {
          cell[n].kappa_grey = (0.9 * f56ni(&cell[n]) + 0.1) * rho_crit/cell[n].rho_init;
        }
        else
        {
          cell[n].kappa_grey = (0.9 * f56ni(&cell[n]) + 0.1);
        }
      }
      else if (cell[n].rho_init == 0)
      {
        cell[n].kappa_grey = 0;
      }
      else if (cell[n].rho_init < 0)
      {
        printout("Error: negative density. Abort.\n");
        abort();
      }
      opcase3_sum += cell[n].kappa_grey*cell[n].rho_init;
    }
  }


  //MK: second loop to set up opacities which need integrated basic quantities for normalisation
  for (n=0;n < ngrid; n++)
  {
    if (cell[n].rho_init > 0)
    {
      if (opacity_case == 0)
      {
        cell[n].kappa_grey = GREY_OP;
      }
      else if (opacity_case == 1)
      {
        cell[n].kappa_grey = ((0.9 * f56ni(&cell[n])) + 0.1) * GREY_OP / ((0.9 *  mni56 / mtot) + 0.1);
      }
      else if (opacity_case == 2)
      {
        opcase2_normal = GREY_OP*rho_sum / ((0.9 *  fe_sum) + (0.1 * (ngrid - empty_cells)));
        cell[n].kappa_grey = opcase2_normal/cell[n].rho_init * ((0.9 * f56ni(&cell[n])) + 0.1);
      }
      else if (opacity_case == 3)
      {
        opcase3_normal = GREY_OP*rho_sum / opcase3_sum;
        cell[n].kappa_grey = cell[n].kappa_grey * opcase3_normal;
      }
      else if (opacity_case == 4)
      {
        cell[n].kappa_grey = SIGMA_T;
      }
      else
      {
        printout("Unknown opacity case. Abort.\n");
        abort();
      }
    }
    else if (cell[n].rho_init == 0)
    {
      cell[n].kappa_grey = 0.0;
    }
    else if (cell[n].rho_init < 0)
    {
      printout("Error: negative density. Abort.\n");
      abort();
    }

    check1 = check1 + (cell[n].kappa_grey  * cell[n].rho_init);
    check2 = check2 + cell[n].rho_init;
  }


  printout("Initial densities assigned uniformly.\n");
  printout("Grey normalisation check: %g\n", check1/check2);

  return 0;
}*/


static void set_fvstable(int modelgridindex, float x)
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

static void set_ftistable(int modelgridindex, float x)
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

static void set_stable_abund(const int mgi, const int anumber, const float elemabundance)
{
  // store the stable mass fraction of an element given the total element mass fraction
  // by subtracting the abundances of radioactive isotopes
  switch (anumber)
  {
    case 28:
      modelgrid[mgi].fnistable = fmax(0.,
        elemabundance - get_modelradioabund(mgi, NUCLIDE_NI56) - get_modelradioabund(mgi, NUCLIDE_NI57));
      break;

    case 27:
      modelgrid[mgi].fcostable = fmax(0.,
        elemabundance - get_modelradioabund(mgi, NUCLIDE_CO56) - get_modelradioabund(mgi, NUCLIDE_CO57));
      break;

    case 26:
      modelgrid[mgi].ffestable = fmax(0., elemabundance - get_modelradioabund(mgi, NUCLIDE_FE52));
      break;

    case 25:
      modelgrid[mgi].fmnstable = fmax(0., elemabundance);
      break;

    case 24:
      modelgrid[mgi].fcrstable = fmax(0., elemabundance - get_modelradioabund(mgi, NUCLIDE_CR48));
      break;

    case 23:
      modelgrid[mgi].fvstable = fmax(0., elemabundance);
      break;

    case 22:
      modelgrid[mgi].ftistable = fmax(0., elemabundance);
      break;
  }
}


int get_cellcoordpointnum(const int cellindex, const int axis)
// convert a cell index number into an integer x, y, or z index
{
  // return cell[cellindex].nxyz[axis];
  switch (axis)
  {
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


extern inline double get_cellcoordmin(int cellindex, int axis);


double get_cellradialpos(const int cellindex)
{
  if (grid_type == GRID_SPHERICAL1D)
  {
    return get_cellcoordmin(cellindex, 0) + (0.5 * wid_init(cellindex));
  }

  double dcen[3];
  for (int axis = 0; axis < 3; axis++)
  {
    dcen[axis] = get_cellcoordmin(cellindex, axis) + (0.5 * wid_init(0));
  }

  return vec_len(dcen);
}


static void calculate_kappagrey(void)
{
  double rho_sum = 0.0;
  double fe_sum = 0.0;
  double opcase3_sum = 0.0;
  int empty_cells = 0;

  for (int n = 0; n < ngrid; n++)
  {
    const int mgi = cell[n].modelgridindex;
    rho_sum += get_rhoinit(mgi);
    fe_sum += get_ffegrp(mgi);

    if (opacity_case == 3)
    {
      if (get_rhoinit(mgi) > 0.)
      {
        double kappagrey = (0.9 * get_ffegrp(mgi) + 0.1);

        if (get_rhoinit(mgi) > rho_crit)
          kappagrey *= rho_crit / get_rhoinit(mgi);

        set_kappagrey(mgi, kappagrey);
      }
      else if (get_rhoinit(mgi) == 0.)
      {
        set_kappagrey(mgi,0.);
      }
      else if (get_rhoinit(mgi) < 0.)
      {
        printout("Error: negative density. Abort.\n");
        abort();
      }
      opcase3_sum += get_kappagrey(mgi) * get_rhoinit(mgi);
    }
  }

  FILE *grid_file;
  if (rank_global == 0)
  {
    grid_file = fopen_required("grid.out", "w");
  }

  /// Second pass through allows calculation of normalized kappa_grey
  double check1 = 0.0;
  double check2 = 0.0;
  for (int n = 0; n < ngrid; n++)
  {
    const int mgi = cell[n].modelgridindex;
    if (rank_global == 0 && mgi != MMODELGRID)
      fprintf(grid_file,"%d %d\n", n, mgi); ///write only non-empty cells to grid file

    if (get_rhoinit(mgi) > 0)
    {
      if (opacity_case == 0)
      {
        set_kappagrey(mgi, GREY_OP);
      }
      else if (opacity_case == 1)
      {
        set_kappagrey(mgi, ((0.9 * get_ffegrp(mgi)) + 0.1) * GREY_OP / ((0.9 * mfeg / mtot) + 0.1));
      }
      else if (opacity_case == 2)
      {
        const double opcase2_normal = GREY_OP * rho_sum / ((0.9 *  fe_sum) + (0.1 * (ngrid - empty_cells)));
        set_kappagrey(mgi, opcase2_normal/get_rhoinit(mgi) * ((0.9 * get_ffegrp(mgi)) + 0.1));
      }
      else if (opacity_case == 3)
      {
        opcase3_normal = GREY_OP * rho_sum / opcase3_sum;
        set_kappagrey(mgi, get_kappagrey(mgi) * opcase3_normal);
      }
      else if (opacity_case == 4)
      {
        ///kappagrey used for initial grey approximation in this case
        set_kappagrey(mgi, ((0.9 * get_ffegrp(mgi)) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1));
        //set_kappagrey(mgi, SIGMA_T);
      }
      else
      {
        printout("Unknown opacity case. Abort.\n");
        abort();
      }
    }
    else if (get_rhoinit(mgi) == 0.)
    {
      set_kappagrey(mgi, 0.);
    }
    else if (get_rhoinit(mgi) < 0.)
    {
      printout("Error: negative density. Abort.\n");
      abort();
    }

    check1 = check1 + (get_kappagrey(mgi) * get_rhoinit(mgi));
    check2 = check2 + get_rhoinit(mgi);
  }
  if (rank_global == 0)
    fclose(grid_file);

  printout("Initial densities taken from readin.\n");
  printout("Grey normalisation check: %g\n", check1 / check2);
}


static void allocate_compositiondata(const int modelgridindex)
/// Initialise composition dependent cell data for the given cell
{
  if ((modelgrid[modelgridindex].composition = (compositionlist_entry *) malloc(nelements * sizeof(compositionlist_entry))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize compositionlist for cell %d... abort\n",modelgridindex);
    abort();
  }

  mem_usage_nltepops += total_nlte_levels * sizeof(double);

  if ((modelgrid[modelgridindex].nlte_pops = (double *) malloc(total_nlte_levels * sizeof(double))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize nlte memory for cell %d... abort\n",modelgridindex);
    abort();
  }

  for (int nlteindex = 0; nlteindex < total_nlte_levels; nlteindex++)
  {
    modelgrid[modelgridindex].nlte_pops[nlteindex] = -1.0; ///flag to indicate that there is
                                                           /// currently no information on the nlte populations
  }

  //printout("Managed to allocate memory for %d nlte levels\n", total_nlte_levels);

  for (int element = 0; element < nelements; element++)
  {
    /// Set initial abundances to zero
    modelgrid[modelgridindex].composition[element].abundance = 0.;

    /// and allocate memory to store the ground level populations for each ionisation stage
    if ((modelgrid[modelgridindex].composition[element].groundlevelpop = (float *) calloc(get_nions(element), sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize groundlevelpoplist for element %d in cell %d... abort\n",element,modelgridindex);
      abort();
    }

    if ((modelgrid[modelgridindex].composition[element].partfunct = (float *) malloc(get_nions(element) * sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize partfunctlist for element %d in cell %d... abort\n",element,modelgridindex);
      abort();
    }
    /*
    if ((modelgrid[n].composition[element].ltepartfunct = malloc(get_nions(element)*sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize lte partfunctlist for element %d in cell %d... abort\n",element,n);
      abort();
    }
    */
  }
}


static void allocate_cooling(const int modelgridindex)
/// Initialise composition dependent cell data for the given cell
{
  if ((modelgrid[modelgridindex].cooling = (mgicooling_t *) malloc(nelements * sizeof(mgicooling_t))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize coolinglist for cell %d... abort\n",modelgridindex);
    abort();
  }

  for (int element = 0; element < nelements; element++)
  {
    /// and allocate memory to store the ground level populations for each ionisation stage
    if ((modelgrid[modelgridindex].cooling[element].contrib = (double *) malloc(get_nions(element) * sizeof(double))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize coolinglist for element %d in cell %d... abort\n",element,modelgridindex);
      abort();
    }
  }
}


static void allocate_nonemptycells(void)
{
  printout("mem_usage: the modelgrid array occupies %.1f MB\n", sizeof(modelgrid) / 1024. / 1024.);
  mem_usage_nltepops = 0;
  /// This is the placeholder for empty cells. Temperatures must be positive
  /// as long as ff opacities are calculated.
  set_rhoinit(MMODELGRID, 0.);
  set_rho(MMODELGRID, 0.);
  set_nne(MMODELGRID, 0.);
  set_ffegrp(MMODELGRID, 0.);
  for (enum radionuclides iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
  {
    set_modelradioabund(MMODELGRID, iso, 0.);
  }
  set_Te(MMODELGRID, MINTEMP);
  set_TJ(MMODELGRID, MINTEMP);
  set_TR(MMODELGRID, MINTEMP);
  allocate_compositiondata(MMODELGRID);
  allocate_cooling(MMODELGRID);

  // Determine the number of simulation cells associated with the model cells
  for (int mgi = 0; mgi < npts_model; mgi++)
    mg_associated_cells[mgi] = 0;

  for (int cellindex = 0; cellindex < ngrid; cellindex++)
  {
    const int mgi = cell[cellindex].modelgridindex;
    assert(!(model_type == RHO_3D_READ) || (get_rhoinit(mgi) > 0) || (mgi == MMODELGRID));
    mg_associated_cells[mgi] += 1;
    assert(!(model_type == RHO_3D_READ) || (mg_associated_cells[mgi] == 1) || (mgi == MMODELGRID));
  }

  int numnonemptycells = 0;
  for (int mgi = 0; mgi < npts_model; mgi++)
  {
    if (get_numassociatedcells(mgi) > 0)
    {
      allocate_compositiondata(mgi);
      allocate_cooling(mgi);

      if (get_rhoinit(mgi) < 0)
      {
        printout("Error: negative density. Abort.\n");
        abort();
      }
      numnonemptycells++;
    }
    else
    {
      set_rhoinit(mgi, 0.);
      set_rho(mgi, 0.);
      for (enum radionuclides iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
      {
        set_modelradioabund(mgi, iso, 0.);
      }
    }
  }

  printout("There are %d modelgrid cells with associated propagation cells\n", numnonemptycells);

  printout("mem_usage: NLTE populations for all allocated cells occupy a total of %.1f MB\n", mem_usage_nltepops / 1024. / 1024.);
}


static void density_1d_read(void)
/// Routine for doing a density grid read from a 1-D model.
{
  for (int cellindex = 0; cellindex < ngrid; cellindex++)
  {
    const double radial_pos = get_cellradialpos(cellindex);
    const double vcell = radial_pos / tmin;
    const double vmin = 0.;
    if (radial_pos < rmax)
    {
      if (grid_type == GRID_SPHERICAL1D)
      {
        cell[cellindex].modelgridindex = cellindex;
      }
      else
      {
        cell[cellindex].modelgridindex = 0;

        for (int mgi = 0; mgi < (npts_model - 1); mgi++)
        {
          if (vout_model[mgi] < vcell)
          {
            cell[cellindex].modelgridindex = mgi + 1;
          }
        }
      }

      if (vout_model[cell[cellindex].modelgridindex] >= vmin)
      {
        modelgrid[cell[cellindex].modelgridindex].initial_radial_pos += radial_pos;
      }
      else
      {
        cell[cellindex].modelgridindex = MMODELGRID;
      }
    }
    else
    {
      cell[cellindex].modelgridindex = MMODELGRID;
    }
  }
}


/*backup copy of old version
  for (n=0;n < ngrid; n++)
  {
  dcen[0] = cell[n].pos_init[0] + (0.5*wid_init);
  dcen[1] = cell[n].pos_init[1] + (0.5*wid_init);
  dcen[2] = cell[n].pos_init[2] + (0.5*wid_init);

  radial_pos = vec_len(dcen);
  if (radial_pos < rmax)
  {
  mkeep =0;
  cell[n].rho_init = cell[n].rho = rho_model[0] * pow( (t_model/tmin), 3.);
  if (opacity_case == 0)
  {
  cell[n].kappa_grey = GREY_OP;
}
  else if (opacity_case == 1)
  {
  cell[n].kappa_grey = ((0.9 * ffegrp_model[0]) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1);
}
  else if (opacity_case == 2)
  {
  Q = GREY_OP*rho_sum / ((0.9 *  fe_sum) + (0.1 * (npts_model)));
  cell[n].kappa_grey = Q/rho_model[0] * ((0.9 * ffegrp_model[0]) + 0.1);
}
  else
  {
  printout("Unknown opacity case. Abort.\n");
  abort();
}
  for (m = 0; m < (npts_model-1); m++)
  {
  if (radial_pos > (vout_model[m] * tmin))
  {
  mkeep = m+1;
  cell[n].rho_init = cell[n].rho = rho_model[m+1] * pow( (t_model/tmin), 3.);
  if (opacity_case == 0)
  {
  cell[n].kappa_grey = GREY_OP;
}
  else if (opacity_case == 1)
  {
  cell[n].kappa_grey = ((0.9 * ffegrp_model[m+1]) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1);
}
  else if (opacity_case == 2)
  {
  Q = GREY_OP*rho_sum / ((0.9 *  fe_sum) + (0.1 * (npts_model)));
  cell[n].kappa_grey = Q/rho_model[m+1] * ((0.9 * ffegrp_model[m+1]) + 0.1);
}
  else
  {
  printout("Unknown opacity case. Abort.\n");
  abort();
}
}
}
  check1 = check1 + (cell[n].kappa_grey  * cell[n].rho_init);
  check2 = check2 + cell[n].rho_init;
  renorm[mkeep]++;

}
  else
  {
  cell[n].rho_init = cell[n].rho = 0.0;
  cell[n].kappa_grey = 0.0;
}

  if (cell[n].rho_init < 0)
  {
  printout("Error: negative density. Abort.\n");
  abort();
}

}
*/

  /* The following block of code was for a check in the 1d models - currently commented out */

  /* *******************************************************************************************

  check1 = 0.0;
  check2 = 0.0;

  if (renorm[0] > 0 )
  {
  den_norm[0] = rho_model[0] * (pow(vout_model[0],3.)) * 4. * PI * pow(tmin,3.) / 3. / (renorm[0] * wid_init * wid_init * wid_init);
}
  else
  {
  printout("No cells assigned in model shell 0\n");
}

  for (n=1;n<MMODELGRID;n++)
  {
  if (renorm[n] > 0 )
  {
  den_norm[n] = rho_model[n] * (pow(vout_model[n],3.) - pow(vout_model[n-1],3.)) * 4. * PI * pow(tmin,3.) / 3. / (renorm[n] * wid_init * wid_init * wid_init);
}
  else if (n < npts_model)
  {
  printout("No cells assigned in model shell %d\n",n);
}
}

  for (n=0;n < ngrid; n++)
  {
  dcen[0] = cell[n].pos_init[0] + (0.5*wid_init);
  dcen[1] = cell[n].pos_init[1] + (0.5*wid_init);
  dcen[2] = cell[n].pos_init[2] + (0.5*wid_init);

  radial_pos = vec_len(dcen);
  if (radial_pos < rmax)
  {
  mkeep =0;
  cell[n].rho_init = cell[n].rho = den_norm[0] * pow( (t_model/tmin), 3.);
  if (opacity_case == 1)
  {
  cell[n].kappa_grey = ((0.9 *  ffegrp_model[0]) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1);
}
  else
  {
  cell[n].kappa_grey = GREY_OP;
}
  for (m = 0; m < (npts_model-1); m++)
  {
  if (radial_pos > (vout_model[m] * tmin))
  {
  mkeep = m+1;
  cell[n].rho_init = cell[n].rho = den_norm[m+1] * pow( (t_model/tmin), 3.);
  if (opacity_case == 1)
  {
  cell[n].kappa_grey = ((0.9 *  ffegrp_model[m+1]) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1);
}
  else
  {
  cell[n].kappa_grey = GREY_OP;
}
}
}
  check1 = check1 + (cell[n].kappa_grey  * cell[n].rho_init);
  check2 = check2 + cell[n].rho_init;
  renorm[mkeep]++;

}
  else
  {
  cell[n].rho_init = cell[n].rho = 0.0;
  cell[n].kappa_grey = 0.0;
}

  if (cell[n].rho_init < 0)
  {
  printout("Error: negative density. Abort.\n");
  abort();
}

}
*/


static void density_2d_read(void)
// Routine for doing a density grid read from a 2-D model.
{
  for (int n = 0; n < ngrid; n++)
  {
    const double radial_pos = get_cellradialpos(n);

    if (radial_pos < rmax)
    {
      double dcen[3];
      for (int d = 0; d < 3; d++)
      {
        const double cellcoordmin = - coordmax[d] + (2 * get_cellcoordpointnum(n, d) * coordmax[d] / ncoordgrid[0]);
        dcen[d] = cellcoordmin + (0.5 * wid_init(0));
      }

      cell[n].modelgridindex = 0;
      const double zcylindrical = dcen[2];
      dcen[2] = 0.0;
      const double rcylindrical = vec_len(dcen);

      // Grid is uniform so only need to search in 1d to get r and z positions

      int mkeep1 = 0;
      for (int m = 0; m < ncoord1_model; m++)
      {
        if (rcylindrical > (m * dcoord1 * tmin/t_model))
        {
          mkeep1 = m;
          //cell[n].modelgridindex = m+1;
        }
      }

      int mkeep2 = 0;
      for (int m = 0; m < ncoord2_model; m++)
      {
        if (zcylindrical > (((m * dcoord2) * tmin/t_model) - rmax))
        {
          mkeep2 = m;
          //cell[n].modelgridindex = m+1;
        }
      }
      cell[n].modelgridindex = (mkeep2 * ncoord1_model) + mkeep1;
      modelgrid[cell[n].modelgridindex].initial_radial_pos += radial_pos;

      //renorm[mkeep]++;
    }
    else
    {
      cell[n].modelgridindex = MMODELGRID;
    }
  }
}


/*void abundances_setup()
{
  double vec_len(double x[3]);
  int get_elementindex(int anumber);
  double dcen[3],m_r;
  int element,n;

  if (homogeneous_abundances)
  {
    for (n = 0; n < ngrid; n++)
    {
      if (cell[n].rho >= MINDENSITY)
      {
        for (element = 0; element < nelements; element++)
        {
          /// Now set the abundances (by mass) of included elements, i.e.
          /// read out the abundances specified in the atomic data file
          cell[n].composition[element].abundance = elements[element].abundance;

          /// No stable Ni and Co isotopes in that case. Their abundances must fit
          /// the given total Ni and Co abundances and we have no separate readin
          /// for this case!
          //if (get_element(element) == 28)
          //  cell[n].f_ni =  elements[element].abundance - cell[n].f_ni;
          //if (get_element(element) == 27)
          //  cell[n].f_co =  elements[element].abundance - cell[n].f_co;
          if (get_element(element) == 26)
            cell[n].f_fe_init =  elements[element].abundance;
        }
      }
    }
  }
  else
  {
    for (n = 0; n < ngrid; n++)
    {
      if (cell[n].rho >= MINDENSITY)
      {
        dcen[0]=cell[n].pos_init[0] + (0.5*wid_init);
        dcen[1]=cell[n].pos_init[1] + (0.5*wid_init);
        dcen[2]=cell[n].pos_init[2] + (0.5*wid_init);

        m_r = vec_len(dcen) / rmax;
        m_r = pow(m_r,3) * mtot / MSUN;

        if (m_r < 0.5)
        {
          /// Inner part consists of pure Nickel
          cell[n].composition[get_elementindex(28)].abundance = 1.0;
          //cell[n].composition[get_elementindex(27)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(26)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(14)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(8)].abundance = 1e-10;
        }
        else if (m_r < 0.75)
        {
          /// In the middle we assume Nickel & intermediate mass elements
          /// which are represented by Si
          cell[n].composition[get_elementindex(28)].abundance = 1.0 - ((m_r - 0.5)*4);
          cell[n].composition[get_elementindex(14)].abundance = (m_r - 0.5)*4;
          //cell[n].composition[get_elementindex(27)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(26)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(8)].abundance = 1e-10;
        }
        else
        {
          /// For the outer shell we assume a mixture of intermediate mass
          /// elements (represented by Si) and unburned material (represented
          /// by O)
          cell[n].composition[get_elementindex(14)].abundance = 0.5;
          cell[n].composition[get_elementindex(8)].abundance = 0.5;
          //cell[n].composition[get_elementindex(28)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(27)].abundance = 1e-10;
          //cell[n].composition[get_elementindex(26)].abundance = 1e-10;
        }
      }
    }
  }
}*/


static void abundances_read(void)
{
  const bool threedimensional = (model_type == RHO_3D_READ);

  /// Open the abundances file
  FILE *abundance_file = fopen_required("abundances.txt", "r");

  /// and process through the grid to read in the abundances per cell
  /// The abundance file should only contain information for non-empty
  /// cells. Its format must be cellnumber (integer), abundance for
  /// element Z=1 (float) up to abundance for element Z=30 (float)
  /// i.e. in total one integer and 30 floats.

  // loop over propagation cells for 3D models, or modelgrid cells
  const int ncount = threedimensional ? ngrid : npts_model;
  for (int n = 0; n < ncount; n++)
  {
    const int mgi = threedimensional ? cell[n].modelgridindex : n;

    int cellnumber;
    fscanf(abundance_file, "%d", &cellnumber);

    if (cellnumber != n + 1)
    {
      printout("[fatal] %s: grid cell mismatch ... abort\n", __func__);
      printout("[fatal] n %d, cellnumber %d\n", n, cellnumber);
      abort();
    }

    double normfactor = 0.;
    float abundances_in[30];
    for (int anumber = 1; anumber <= 30; anumber++)
    {
      fscanf(abundance_file, "%g", &abundances_in[anumber - 1]);

      assert(abundances_in[anumber - 1] >= 0);
      normfactor += abundances_in[anumber - 1];
    }

    if (get_numassociatedcells(mgi) > 0)
    {
      if (threedimensional || normfactor <= 0.)
        normfactor = 1.;

      for (int element = 0; element < nelements; element++)
      {
        ///now set the abundances (by mass) of included elements, i.e.
        ///read out the abundances specified in the atomic data file
        const int anumber = get_element(element);
        const float elemabundance = abundances_in[anumber - 1] / normfactor;
        assert(elemabundance >= 0.);
        modelgrid[mgi].composition[element].abundance = elemabundance;

        // radioactive nuclide abundances should have already been set by read_model
        set_stable_abund(mgi, anumber, elemabundance);
      }
    }
  }

  fclose(abundance_file);
}


static void read_grid_restart_data(void)
{
  printout("READIN GRID SNAPSHOT\n");
  FILE *restrict gridsave_file = fopen_required("gridsave.dat", "r");

  for (int mgi = 0; mgi < npts_model; mgi++)
  {
    int mgi_in;
    float T_R;
    float T_e;
    float W;
    float T_J;
    fscanf(gridsave_file, "%d %g %g %g %g %hd %lg",
           &mgi_in, &T_R, &T_e, &W, &T_J,
           &modelgrid[mgi].thick, &rpkt_emiss[mgi]);

    if (mgi_in != mgi)
    {
      printout("[fatal] read_grid_restart_data: cell mismatch in reading input gridsave.dat ... abort\n");
      printout("[fatal] read_grid_restart_data: read cellnumber %d, expected cellnumber %d\n",mgi_in,mgi);
      abort();
    }

    assert(T_R >= 0);
    assert(T_e >= 0);
    assert(W >= 0);
    assert(T_J >= 0);
    assert(rpkt_emiss[mgi] >= 0);

    set_TR(mgi, T_R);
    set_Te(mgi, T_e);
    set_W(mgi, W);
    set_TJ(mgi, T_J);

    #ifndef FORCE_LTE
      #if (!NO_LUT_PHOTOION)
        for (int element = 0; element < nelements; element++)
        {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions; ion++)
          {
            const int estimindex = mgi * nelements * maxion + element * maxion + ion;
            fscanf(gridsave_file, " %lg %lg", &corrphotoionrenorm[estimindex], &gammaestimator[estimindex]);
          }
        }
      #endif
    #endif
  }

  // the order of these calls is very important!
  radfield_read_restart_data(gridsave_file);
  nt_read_restart_data(gridsave_file);
  nltepop_read_restart_data(gridsave_file);
  fclose(gridsave_file);
}


static void assign_temperature(void)
/// Routine for assigning temperatures to the grid cells at the start of the simulation.
{
  if (simulation_continued_from_saved)
  {
    /// For continuation of an existing simulation we read the temperatures
    /// at the end of the simulation and write them to the grid.
    read_grid_restart_data();
  }
  else
  {
    /// For a simulation started from scratch we estimate the initial temperatures

    /// We assume that for early times the material is so optically thick, that
    /// all the radiation is trapped in the cell it originates from. This
    /// means furthermore LTE, so that both temperatures can be evaluated
    /// according to the local energy density resulting from the 56Ni decay.
    /// The dilution factor is W=1 in LTE.

    const double tstart = time_step[0].mid;

    const double factor56ni = 1. / 56 / MH * (-1. / (tstart * (- T56CO + T56NI)))
      * (- E56NI * exp(- tstart / T56NI) * tstart * T56CO - E56NI * exp(- tstart / T56NI) * T56NI * T56CO
         + E56NI * exp(- tstart / T56NI) * tstart * T56NI + pow(T56NI, 2) * E56NI * exp(- tstart / T56NI)
         - T56CO * tstart * E56CO * exp(- tstart / T56CO) - pow(T56CO, 2) * E56CO * exp(- tstart / T56CO)
         + E56CO * tstart * T56NI * exp(- tstart / T56NI) + pow(T56NI, 2) * E56CO * exp(- tstart / T56NI)
         + E56NI * T56CO * T56NI - E56NI * pow(T56NI, 2) - pow(T56NI, 2) * E56CO + E56CO * pow(T56CO, 2));

    const double factor56co = 1. / 56 / MH * (1. / (tstart * T56CO))
      * (T56CO * tstart * E56CO * exp(- tstart / T56CO) + pow(T56CO, 2) * E56CO * exp(- tstart / T56CO));

    const double factor57ni = 1. / 57 / MH * (-1. / (tstart * (- T57CO + T57NI)))
      * (- E57NI * exp(- tstart / T57NI) * tstart * T57CO - E57NI * exp(- tstart / T57NI) * T57NI * T57CO
         + E57NI * exp(- tstart / T57NI) * tstart * T57NI + pow(T57NI, 2) * E57NI * exp(- tstart / T57NI)
         - T57CO * tstart * E57CO * exp(- tstart / T57CO) - pow(T57CO, 2) * E57CO * exp(- tstart / T57CO)
         + E57CO * tstart * T57NI * exp(- tstart / T57NI) + pow(T57NI, 2) * E57CO * exp(- tstart / T57NI)
         + E57NI * T57CO * T57NI - E57NI * pow(T57NI, 2) - pow(T57NI, 2) * E57CO + E57CO * pow(T57CO, 2));

    const double factor52fe = 1. / 52 / MH * (-1. / (tstart * (- T52MN + T52FE)))
      * (- E52FE * exp(- tstart / T52FE) * tstart * T52MN - E52FE * exp(- tstart / T52FE) * T52FE * T52MN
         + E52FE * exp(- tstart / T52FE) * tstart * T52FE + pow(T52FE, 2) * E52FE * exp(- tstart / T52FE)
         - T52MN * tstart * E52MN * exp(- tstart / T52MN) - pow(T52MN, 2) * E52MN * exp(- tstart / T52MN)
         + E52MN * tstart * T52FE * exp(- tstart / T52FE) + pow(T52FE, 2) * E52MN * exp(- tstart / T52FE)
         + E52FE * T52MN * T52FE - E52FE * pow(T52FE, 2) - pow(T52FE, 2) * E52MN + E52MN * pow(T52MN, 2));

    const double factor48cr = 1. / 48 / MH * (-1. / (tstart * (- T48V + T48CR)))
      * (- E48CR * exp(- tstart / T48CR) * tstart * T48V - E48CR * exp(- tstart / T48CR) * T48CR * T48V
         + E48CR * exp(- tstart / T48CR) * tstart * T48CR + pow(T48CR, 2) * E48CR * exp(- tstart / T48CR)
         - T48V * tstart * E48V * exp(- tstart / T48V) - pow(T48V, 2) * E48V * exp(- tstart / T48V)
         + E48V * tstart * T48CR * exp(- tstart / T48CR) + pow(T48CR, 2) * E48V * exp(- tstart / T48CR)
         + E48CR * T48V * T48CR - E48CR * pow(T48CR, 2) - pow(T48CR, 2) * E48V + E48V * pow(T48V, 2));

    // printout("factor56ni %g\n", factor56ni);
    // printout("factor56co %g\n", factor56co);
    // printout("factor57ni %g\n", factor57ni);
    //factor56ni = CLIGHT/4/STEBO * E56NI/56/MH;
    /// This works only for the inbuilt Lucy model
    //factor56ni = CLIGHT/4/STEBO * 3*mtot/4/PI * E56NI/56/MH  / pow(vmax,3);
    //for (n = 0; n < ngrid; n++)
    for (int n = 0; n < npts_model; n++)
    {
      //mgi = cell[n].modelgridindex;
      double T_initial = pow(CLIGHT / 4 / STEBO  * pow(tmin / tstart, 3) * get_rhoinit(n) * (
           (factor56ni * get_modelradioabund(n, NUCLIDE_NI56)) +
           (factor56co * get_modelradioabund(n, NUCLIDE_CO56)) +
           (factor57ni * get_modelradioabund(n, NUCLIDE_NI57)) +
           // (factor57co * get_modelradioabund(n, NUCLIDE_CO57)) +
           (factor52fe * get_modelradioabund(n, NUCLIDE_FE52)) +
           (factor48cr * get_modelradioabund(n, NUCLIDE_CR48))), 1. / 4.);

      //T_initial = pow(factor56ni * cell[n].f_ni * cell[n].rho_init * (1.-exp(-tmin/T56NI)), 1./4.);
      //T_initial = pow(factor56ni * cell[n].f_ni * (1.-exp(-tmin/T56NI))/pow(tmin,3), 1./4.);
      //T_initial = 30615.5;
      if (T_initial < MINTEMP)
      {
        T_initial = MINTEMP;
      }
      else if (T_initial > MAXTEMP)
      {
        T_initial = MAXTEMP;
      }

      set_Te(n, T_initial);
      set_TJ(n, T_initial);
      set_TR(n, T_initial);

      set_W(n, 1.);
    }
  }
}


static void uniform_grid_setup(void)
/// Routine for doing a uniform cuboidal grid.
{
  coordlabel[0] = 'X';
  coordlabel[1] = 'Y';
  coordlabel[2] = 'Z';
  int nxyz[3] = {0, 0, 0};
  for (int n = 0; n < ngrid; n++)
  {
    for (int axis = 0; axis < 3; axis++)
    {
      assert(nxyz[axis] == get_cellcoordpointnum(n, axis));
      cell[n].pos_init[axis] = - coordmax[axis] + (2 * nxyz[axis] * coordmax[axis] / ncoordgrid[axis]);
    }

    //cell[n].cen_init[0] = cell[n].pos_init[0] + (0.5 * wid_init);
    //cell[n].cen_init[1] = cell[n].pos_init[1] + (0.5 * wid_init);
    //cell[n].cen_init[2] = cell[n].pos_init[2] + (0.5 * wid_init);

    // cell[n].xyz[0] = nx;
    // cell[n].xyz[1] = ny;
    // cell[n].xyz[2] = nz;

    assert(n == nxyz[2] * ncoordgrid[1] * ncoordgrid[0] + nxyz[1] * ncoordgrid[0] + nxyz[0]);

    nxyz[0]++;
    if (nxyz[0] == ncoordgrid[0])
    {
      nxyz[0] = 0;
      nxyz[1]++;
    }
    if (nxyz[1] == ncoordgrid[1])
    {
      nxyz[1] = 0;
      nxyz[2]++;
    }

    ///Do we need this initialisation anywhere else (after modelgridindex was initialised) ???????????????????????????????????
    //cell[n].f_ni_stable = 0.0;
    //cell[n].f_co_stable = 0.0;
    //cell[n].f_fe_init = 0.0;
  }

  /*
  /// Finally we must also create the composition dependent data structure for
  /// the samplingcell which is located at MGRID (i.e. the MGRID+1th cell)
  n = MGRID;
  if ((cell[n].composition = malloc(nelements*sizeof(compositionlist_entry))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize compositionlist for cell %d... abort\n",n);
    abort();
  }
  for (element = 0; element < nelements; element++)
  {
    ///now set the abundances (by mass) of included elements, i.e.
    ///read out the abundances specified in the atomic data file
    ///and allocate memory to store the ground level populations for each ionisation stage
    if ((cell[n].composition[element].groundlevelpop = malloc(get_nions(element)*sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize groundlevelpoplist for element %d in cell %d... abort\n",element,n);
      abort();
    }
    if ((cell[n].composition[element].partfunct = malloc(get_nions(element)*sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize partfunctlist for element %d in cell %d... abort\n",element,n);
      abort();
    }
//     if ((cell[n].composition[element].ltepartfunct = malloc(get_nions(element)*sizeof(float))) == NULL)
//     {
//       printout("[fatal] input: not enough memory to initialize lte partfunctlist for element %d in cell %d... abort\n",element,n);
//       abort();
//     }
  }
  */
}

static void spherical1d_grid_setup(void)
{
  assert(model_type == RHO_1D_READ);
  coordlabel[0] = 'r';
  coordlabel[1] = '?';
  coordlabel[2] = '?';

  ncoordgrid[0] = npts_model;
  ncoordgrid[1] = 1;
  ncoordgrid[2] = 1;
  ngrid = ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2];
  coordmax[0] = rmax;
  coordmax[1] = 0.;
  coordmax[2] = 0.;

  // in this mode, cellindex and modelgridindex are the same thing
  for (int cellindex = 0; cellindex < npts_model; cellindex++)
  {
    const double v_inner = cellindex > 0 ? vout_model[cellindex - 1] : 0.;
    cell[cellindex].modelgridindex = cellindex;
    cell[cellindex].pos_init[0] = v_inner * tmin;
    cell[cellindex].pos_init[1] = 0.;
    cell[cellindex].pos_init[2] = 0.;
  }
}


void grid_init(int my_rank)
/// Initialises the propagation grid cells and associates them with modelgrid cells
{
  /// Start by checking that the number of grid cells is okay */
  //ngrid = ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2]; ///Moved to input.c
  //if (ngrid > MGRID)
  //{
  //  printout("[fatal] grid_init: Error: too many grid cells. Abort.");
  //  abort();
  //}
  for (int n = 0; n <= MMODELGRID; n++)
  {
    modelgrid[n].initial_radial_pos = 0;
  }

  /// The cells will be ordered by x then y, then z. Call a routine that
  /// sets up the initial positions and widths of the cells.
  if (grid_type == GRID_UNIFORM)
  {
    uniform_grid_setup();
  }
  else if (grid_type == GRID_SPHERICAL1D)
  {
    spherical1d_grid_setup();
  }
  else
  {
    printout("[fatal] grid_init: Error: Unknown grid type. Abort.");
    abort();
  }

  /// Now set up the density in each cell.
  if (model_type == RHO_UNIFORM)
  {
    //uniform_density_setup ();
    //abundances_setup();
  }
  else
  {
    // Calculate the critical opacity at which opacity_case 3 switches from a
    // regime proportional to the density to a regime independent of the density
    // This is done by solving for tau_sobolev == 1
    // tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/MNI56 * 3000e-8 * time_step[m].mid;
    rho_crit = ME * CLIGHT * MNI56 / (PI * QE * QE * rho_crit_para * 3000e-8 * tmin);
    printout("grid_init: rho_crit = %g\n", rho_crit);

    if (model_type == RHO_1D_READ)
    {
      density_1d_read();
    }
    else if (model_type == RHO_2D_READ)
    {
      density_2d_read();
    }
    else if (model_type == RHO_3D_READ)
    {
      for (int n = 0; n < ngrid; n++)
      {
        const double radial_pos = get_cellradialpos(n);
        modelgrid[cell[n].modelgridindex].initial_radial_pos = radial_pos;
      }
      // cells with rho > 0 are allocated by the above function
    }
    else
    {
      printout("[fatal] grid_init: Error: Unknown density type. Abort.");
      abort();
    }
    allocate_nonemptycells();
    calculate_kappagrey();
    abundances_read();
  }

  radfield_init(my_rank);
  if (NT_ON)
    nt_init(my_rank);

  /// and assign a temperature to the cells
  assign_temperature();

  /// Finally determine which cells are non-empty...
  /*if ((nonemptycells = malloc(ngrid*sizeof(int))) == NULL)
  {
    printout("[fatal] grid_init: not enough memory to initialize the list of non-empty cells\n");
    abort();
  }
  i = 0;
  for (n = 0; n < ngrid; n++)
  {
    if (cell[n].rho_init > MINDENSITY)
    {
      nonemptycells[i] = n;
      i++;
    }
  }
  nnonemptycells = i;*/
  /*  if ((nonemptycells = realloc(nonemptycells, nnonemptycells*sizeof(int))) == NULL)
  {
    printout("[fatal] grid_init: problem during reallocation of list of non-empty cells\n");
    abort();
  }
  */
}
