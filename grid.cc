#include "sn3d.h"
#include "atomic.h"
#include "grid.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "decay.h"
#include "radfield.h"
#include "rpkt.h"
#include "vectors.h"
#include <cstring>


__managed__ enum model_types model_type = RHO_1D_READ;
__managed__ int npts_model = 0; // number of points in 1-D input model

__managed__ double vout_model[MMODELGRID];
__managed__ int ncoord_model[3]; // the model.txt input grid dimensions
__managed__ double dcoord1;
__managed__ double dcoord2; // spacings of a 2D model grid - must be uniform grid

__managed__ double min_den; // minimum model density

__managed__ double mtot;
__managed__ double mfeg;              /// Total mass of Fe group elements in ejecta

__managed__ CELL cell[MGRID + 1];

static long mem_usage_nltepops = 0;

static __managed__ int mg_associated_cells[MMODELGRID + 1];

__managed__ double *totmassradionuclide = NULL; /// total mass of each radionuclide in the ejecta

__host__ __device__
double wid_init(const int cellindex)
// for a uniform grid this is the extent along the x,y,z coordinate (x_2 - x_1, etc.)
// for spherical grid this is the radial extent (r_outer - r_inner)
// these values are for time globals::tmin
{
  switch (globals::grid_type)
  {
    case GRID_SPHERICAL1D:
    {
      const int modelgridindex = get_cell_modelgridindex(cellindex);
      const double v_inner = modelgridindex > 0 ? vout_model[modelgridindex - 1] : 0.;
      return (vout_model[modelgridindex] - v_inner) * globals::tmin;
    }

    default:
      return 2 * globals::coordmax[0] / globals::ncoordgrid[0];
  }
}


__host__ __device__
double vol_init_modelcell(const int modelgridindex)
// return the model cell volume at globals::tmin
// for a uniform cubic grid this is constant
{
  switch (globals::grid_type)
  {
    case GRID_SPHERICAL1D:
      return 4./3. * PI * (pow(globals::tmin * vout_model[modelgridindex], 3) - pow(globals::tmin * (modelgridindex > 0 ? vout_model[modelgridindex - 1] : 0.), 3));

    default:
    {
      const int assoc_cells = get_numassociatedcells(modelgridindex);
      return (wid_init(0) * wid_init(0) * wid_init(0)) * assoc_cells;
    }
  }
}


__host__ __device__
double vol_init_gridcell(const int cellindex)
// return the propagation cell volume at globals::tmin
// for a spherical grid, the cell index is required (and should be equivalent to a modelgridindex)
{
  switch (globals::grid_type)
  {
    case GRID_SPHERICAL1D:
    {
      const int mgi = get_cell_modelgridindex(cellindex);
      return vol_init_modelcell(mgi);
    }

    default:
      return (wid_init(0) * wid_init(0) * wid_init(0));
  }
}


__host__ __device__
double get_cellcoordmin(const int cellindex, const int axis)
// get the minimum value of a coordinate at globals::tmin (xyz or radial coords) of a propagation cell
// e.g., the minimum x position in xyz coords, or the minimum radius
{
  return cell[cellindex].pos_init[axis];
  // return - coordmax[axis] + (2 * get_cellcoordpointnum(cellindex, axis) * coordmax[axis] / globals::ncoordgrid[axis]);
}


__host__ __device__
int get_coordcellindexincrement(const int axis)
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
          return -1;
      }
  }
}


__host__ __device__
int get_cellcoordpointnum(const int cellindex, const int axis)
// convert a cell index number into an integer (x,y,z or r) coordinate index from 0 to globals::ncoordgrid[axis]
{
  // return cell[cellindex].nxyz[axis];

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
          return -1;
      }
  }
}



__host__ __device__
int get_ngriddimensions(void)
{
  return (globals::grid_type == GRID_SPHERICAL1D) ? 1 : 3;
}


__host__ __device__
float get_rhoinit(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].rhoinit;
}


__host__ __device__
float get_rho(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].rho;
}


__host__ __device__
float get_nne(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].nne;
}


__host__ __device__
float get_nnetot(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].nnetot;
}


__host__ __device__
// the abundances referred to below are initial abundances
float get_ffegrp(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].ffegrp;
}


__host__ __device__
float get_elem_abundance(int modelgridindex, int element)
// mass fraction of an element (all isotopes combined)
{
  return globals::modelgrid[modelgridindex].composition[element].abundance;
}


__host__ __device__
void set_elem_abundance(int modelgridindex, int element, float newabundance)
// mass fraction of an element (all isotopes combined)
{
  globals::modelgrid[modelgridindex].composition[element].abundance = newabundance;
}


__host__ __device__
float get_kappagrey(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].kappagrey;
}


__host__ __device__
float get_Te(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].Te;
}


__host__ __device__
float get_TR(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].TR;
}


__host__ __device__
float get_TJ(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].TJ;
}


__host__ __device__
float get_W(int modelgridindex)
{
  return globals::modelgrid[modelgridindex].W;
}


__host__ __device__
void set_rhoinit(int modelgridindex, float x)
{
  globals::modelgrid[modelgridindex].rhoinit = x;
}


__host__ __device__
void set_rho(int modelgridindex, float x)
{
  globals::modelgrid[modelgridindex].rho = x;
}


__host__ __device__
void set_nne(int modelgridindex, float x)
{
  globals::modelgrid[modelgridindex].nne = x;
}


__host__ __device__
void set_nnetot(int modelgridindex, float x)
{
  globals::modelgrid[modelgridindex].nnetot = x;
}


__host__ __device__
void set_ffegrp(int modelgridindex, float x)
{
  assert_always(x >= 0);
  assert_always(x <= 1.001);
  globals::modelgrid[modelgridindex].ffegrp = x;
}


__host__ __device__
void set_kappagrey(int modelgridindex, float kappagrey)
{
  globals::modelgrid[modelgridindex].kappagrey = kappagrey;
}


__host__ __device__
void set_Te(int modelgridindex, float Te)
{
  globals::modelgrid[modelgridindex].Te = Te;
}


__host__ __device__
void set_TR(int modelgridindex, float TR)
{
  globals::modelgrid[modelgridindex].TR = TR;
}


__host__ __device__
void set_TJ(int modelgridindex, float TJ)
{
  globals::modelgrid[modelgridindex].TJ = TJ;
}


__host__ __device__
void set_W(int modelgridindex, float W)
{
  globals::modelgrid[modelgridindex].W = W;
}


__host__ __device__
enum model_types get_model_type(void)
{
  return model_type;
}


__host__ __device__
void set_model_type(enum model_types model_type_value)
{
  model_type = model_type_value;
}


__host__ __device__
int get_npts_model(void)
{
  return npts_model;
}


__host__ __device__
void set_npts_model(int new_npts_model)
{
  if (new_npts_model >= MMODELGRID)
  {
    printout("ERROR: %d = npts_model >= MMODELGRID = %d.\n", new_npts_model, MMODELGRID);
    abort();
  }
  npts_model = new_npts_model;
}


__host__ __device__
int get_cell_modelgridindex(int cellindex)
{
  assert_testmodeonly(cellindex < globals::ngrid);
  return cell[cellindex].modelgridindex;
}


__host__ __device__
static void set_cell_modelgridindex(int cellindex, int new_modelgridindex)
{
  assert_testmodeonly(cellindex < globals::ngrid);
  assert_testmodeonly(new_modelgridindex < npts_model || new_modelgridindex == MMODELGRID);
  cell[cellindex].modelgridindex = new_modelgridindex;
}


__host__ __device__
int get_numassociatedcells(const int modelgridindex)
// number of propagation cells associated with each modelgrid cell
{
  return mg_associated_cells[modelgridindex];
}


__host__ __device__
float get_modelinitradioabund(const int modelgridindex, const int z, const int a)
{
  // this function replaces get_f56ni(mgi), get_fco56(mgi), etc.

  const int nucindex = decay::get_nuc_index(z, a);
  assert_always(decay::get_nuc_z(nucindex) >= 0); // check not FAKE_GAM_LINE_ID nuclide

  assert_always(globals::modelgrid[modelgridindex].initradioabund != NULL);
  return globals::modelgrid[modelgridindex].initradioabund[nucindex];
}


__host__ __device__
void set_modelinitradioabund(const int modelgridindex, const int z, const int a, const float abund)
{
  assert_always(abund >= 0.);
  assert_always(abund <= 1.);

  const int nucindex = decay::get_nuc_index(z, a);
  assert_always(decay::get_nuc_z(nucindex) >= 0 || abund == 0.); // check not FAKE_GAM_LINE_ID nuclide

  assert_always(globals::modelgrid[modelgridindex].initradioabund != NULL);
  globals::modelgrid[modelgridindex].initradioabund[nucindex] = abund;
}


__host__ __device__
float get_stable_abund(const int mgi, const int anumber)
{
  switch (anumber)
  {
    case 28:
      return globals::modelgrid[mgi].fnistable;

    case 27:
      return globals::modelgrid[mgi].fcostable;

    case 26:
      return globals::modelgrid[mgi].ffestable;

    case 25:
      return globals::modelgrid[mgi].fmnstable;

    case 24:
      return globals::modelgrid[mgi].fcrstable;

    case 23:
      return globals::modelgrid[mgi].fvstable;

    case 22:
      return globals::modelgrid[mgi].ftistable;

    default:
      printout("ERROR: no stable abundance variable for element Z=%d\n", anumber);
      abort();
      return -1;
  }
}


__host__ __device__
static void set_elem_stable_abund_from_total(const int mgi, const int anumber, const float elemabundance)
{
  // set the stable mass fraction of an element from the total element mass fraction
  // by subtracting the abundances of radioactive isotopes.
  // if the element Z=anumber has no specific stable abundance variable then the function does nothing
  switch (anumber)
  {
    case 28:
    {
      globals::modelgrid[mgi].fnistable = elemabundance - get_modelinitradioabund(mgi, 28, 56) - get_modelinitradioabund(mgi, 28, 57);
      if (globals::modelgrid[mgi].fnistable < 0.)
      {
          printout("WARNING: cell %d Ni element abundance is less than the sum of isotopic abundances \n", mgi);
          printout("  X_Ni %g X_Ni56 %g X_Ni57 %g\n", elemabundance,
                   get_modelinitradioabund(mgi, 28, 56), get_modelinitradioabund(mgi, 28, 57));
          assert_always(globals::modelgrid[mgi].fnistable >= -1e-3);  // result is allowed to be slightly negative due to roundoff error
          globals::modelgrid[mgi].fnistable = fmax(0., globals::modelgrid[mgi].fnistable); // bring up to zero if negative
      }
      break;
    }

    case 27:
    {
      globals::modelgrid[mgi].fcostable = elemabundance - get_modelinitradioabund(mgi, 27, 56) - get_modelinitradioabund(mgi, 27, 57);
      if (globals::modelgrid[mgi].fcostable < 0)  // result can be slightly negative due to roundoff error
      {
        printout("WARNING: cell %d Co element abundance is less than the sum of isotopic abundances\n", mgi);
        printout("  X_Co %g X_Co56 %g X_Co57 %g\n", elemabundance,
                 get_modelinitradioabund(mgi, 27, 56), get_modelinitradioabund(mgi, 27, 57));
        assert_always(globals::modelgrid[mgi].fcostable >= -1e-3);  // result is allowed to be slightly negative due to roundoff error
        globals::modelgrid[mgi].fcostable = fmax(0., globals::modelgrid[mgi].fcostable); // bring up to zero if negative
      }
      break;
    }

    case 26:
    {
      globals::modelgrid[mgi].ffestable = elemabundance - get_modelinitradioabund(mgi, 26, 52);
      assert_always(globals::modelgrid[mgi].ffestable >= -2e-5);
      globals::modelgrid[mgi].ffestable = fmax(0., globals::modelgrid[mgi].ffestable);
      break;
    }

    case 25:
    {
      globals::modelgrid[mgi].fmnstable = elemabundance;
      assert_always(globals::modelgrid[mgi].fmnstable >= 0.);
      break;
    }

    case 24:
    {
      globals::modelgrid[mgi].fcrstable = elemabundance - get_modelinitradioabund(mgi, 24, 48);
      assert_always(globals::modelgrid[mgi].fcrstable >= -2e-5);
      globals::modelgrid[mgi].fcrstable = fmax(0., globals::modelgrid[mgi].fcrstable);
      break;
    }

    case 23:
    {
      globals::modelgrid[mgi].fvstable = elemabundance;
      assert_always(globals::modelgrid[mgi].fvstable >= 0.);
      break;
    }

    case 22:
    {
      globals::modelgrid[mgi].ftistable = elemabundance;
      assert_always(globals::modelgrid[mgi].ftistable >= 0.);
      break;
    }
  }
}


__host__ __device__
double get_cellradialpos(const int cellindex)
// get the radial distance from the origin to the centre of the cell
{
  // spherical coordinate case is trivial
  if (globals::grid_type == GRID_SPHERICAL1D)
  {
    return get_cellcoordmin(cellindex, 0) + (0.5 * wid_init(cellindex));
  }

  // cubic grid requires taking the length of the 3D position vector
  double dcen[3];
  for (int axis = 0; axis < 3; axis++)
  {
    dcen[axis] = get_cellcoordmin(cellindex, axis) + (0.5 * wid_init(0));
  }

  return vec_len(dcen);
}


__host__ __device__
int get_elements_uppermost_ion(const int modelgridindex, const int element)
{
  return globals::modelgrid[modelgridindex].elements_uppermost_ion[element];
}


__host__ __device__
void set_elements_uppermost_ion(const int modelgridindex, const int element, const int newvalue)
{
  globals::modelgrid[modelgridindex].elements_uppermost_ion[element] = newvalue;
}


static void calculate_kappagrey(void)
{
  double rho_sum = 0.0;
  double fe_sum = 0.0;
  double opcase3_sum = 0.0;
  int empty_cells = 0;

  for (int n = 0; n < globals::ngrid; n++)
  {
    const int mgi = get_cell_modelgridindex(n);
    rho_sum += get_rhoinit(mgi);
    fe_sum += get_ffegrp(mgi);

    if (globals::opacity_case == 3)
    {
      if (get_rhoinit(mgi) > 0.)
      {
        double kappagrey = (0.9 * get_ffegrp(mgi) + 0.1);

        if (get_rhoinit(mgi) > globals::rho_crit)
          kappagrey *= globals::rho_crit / get_rhoinit(mgi);

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
  if (globals::rank_global == 0)
  {
    grid_file = fopen_required("grid.out", "w");
  }

  /// Second pass through allows calculation of normalized kappa_grey
  double check1 = 0.0;
  double check2 = 0.0;
  for (int n = 0; n < globals::ngrid; n++)
  {
    const int mgi = get_cell_modelgridindex(n);
    if (globals::rank_global == 0 && mgi != MMODELGRID)
      fprintf(grid_file,"%d %d\n", n, mgi); ///write only non-empty cells to grid file

    if (get_rhoinit(mgi) > 0)
    {
      if (globals::opacity_case == 0)
      {
        set_kappagrey(mgi, GREY_OP);
      }
      else if (globals::opacity_case == 1)
      {
        set_kappagrey(mgi, ((0.9 * get_ffegrp(mgi)) + 0.1) * GREY_OP / ((0.9 * mfeg / mtot) + 0.1));
      }
      else if (globals::opacity_case == 2)
      {
        const double opcase2_normal = GREY_OP * rho_sum / ((0.9 *  fe_sum) + (0.1 * (globals::ngrid - empty_cells)));
        set_kappagrey(mgi, opcase2_normal/get_rhoinit(mgi) * ((0.9 * get_ffegrp(mgi)) + 0.1));
      }
      else if (globals::opacity_case == 3)
      {
        globals::opcase3_normal = GREY_OP * rho_sum / opcase3_sum;
        set_kappagrey(mgi, get_kappagrey(mgi) * globals::opcase3_normal);
      }
      else if (globals::opacity_case == 4)
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
  if (globals::rank_global == 0)
    fclose(grid_file);

  printout("Initial densities taken from readin.\n");
  printout("Grey normalisation check: %g\n", check1 / check2);
}


static void allocate_compositiondata(const int modelgridindex)
/// Initialise composition dependent cell data for the given cell
{
  globals::modelgrid[modelgridindex].elements_uppermost_ion = (int *) malloc(get_nelements() * sizeof(int));
  assert_always(globals::modelgrid[modelgridindex].elements_uppermost_ion != NULL);

  if ((globals::modelgrid[modelgridindex].composition = (compositionlist_entry *) malloc(get_nelements() * sizeof(compositionlist_entry))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize compositionlist for cell %d... abort\n",modelgridindex);
    abort();
  }

  mem_usage_nltepops += globals::total_nlte_levels * sizeof(double);

  if ((globals::modelgrid[modelgridindex].nlte_pops = (double *) malloc(globals::total_nlte_levels * sizeof(double))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize nlte memory for cell %d... abort\n",modelgridindex);
    abort();
  }

  for (int nlteindex = 0; nlteindex < globals::total_nlte_levels; nlteindex++)
  {
    globals::modelgrid[modelgridindex].nlte_pops[nlteindex] = -1.0; ///flag to indicate that there is
                                                           /// currently no information on the nlte populations
  }

  //printout("Managed to allocate memory for %d nlte levels\n", total_nlte_levels);

  for (int element = 0; element < get_nelements(); element++)
  {
    /// Set initial abundances to zero
    globals::modelgrid[modelgridindex].composition[element].abundance = 0.;

    /// and allocate memory to store the ground level populations for each ionisation stage
    if ((globals::modelgrid[modelgridindex].composition[element].groundlevelpop = (float *) calloc(get_nions(element), sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize groundlevelpoplist for element %d in cell %d... abort\n",element,modelgridindex);
      abort();
    }

    if ((globals::modelgrid[modelgridindex].composition[element].partfunct = (float *) malloc(get_nions(element) * sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize partfunctlist for element %d in cell %d... abort\n",element,modelgridindex);
      abort();
    }
    /*
    if ((globals::modelgrid[n].composition[element].ltepartfunct = malloc(get_nions(element)*sizeof(float))) == NULL)
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
  if ((globals::modelgrid[modelgridindex].cooling = (mgicooling_t *) malloc(get_nelements() * sizeof(mgicooling_t))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize coolinglist for cell %d... abort\n",modelgridindex);
    abort();
  }

  for (int element = 0; element < get_nelements(); element++)
  {
    /// and allocate memory to store the ground level populations for each ionisation stage
    if ((globals::modelgrid[modelgridindex].cooling[element].contrib = (double *) malloc(get_nions(element) * sizeof(double))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize coolinglist for element %d in cell %d... abort\n",element,modelgridindex);
      abort();
    }
  }
}


static void allocate_nonemptymodelcells(void)
{
  printout("[info] mem_usage: the modelgrid array occupies %.1f MB\n", sizeof(globals::modelgrid) / 1024. / 1024.);
  mem_usage_nltepops = 0;
  /// This is the placeholder for empty cells. Temperatures must be positive
  /// as long as ff opacities are calculated.
  set_rhoinit(MMODELGRID, 0.);
  set_rho(MMODELGRID, 0.);
  set_nne(MMODELGRID, 0.);
  set_ffegrp(MMODELGRID, 0.);
  globals::modelgrid[MMODELGRID].initradioabund = (float *) calloc(decay::get_num_nuclides(), sizeof(float));
  assert_always(globals::modelgrid[MMODELGRID].initradioabund != NULL);
  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++)
  {
    const int z = decay::get_nuc_z(nucindex);
    const int a = decay::get_nuc_a(nucindex);
    set_modelinitradioabund(MMODELGRID, z, a, 0.);
  }
  set_Te(MMODELGRID, MINTEMP);
  set_TJ(MMODELGRID, MINTEMP);
  set_TR(MMODELGRID, MINTEMP);
  allocate_compositiondata(MMODELGRID);
  allocate_cooling(MMODELGRID);

  // Determine the number of simulation cells associated with the model cells
  for (int mgi = 0; mgi < (MMODELGRID + 1); mgi++)
    mg_associated_cells[mgi] = 0;

  for (int cellindex = 0; cellindex < globals::ngrid; cellindex++)
  {
    const int mgi = get_cell_modelgridindex(cellindex);
    assert_always(!(get_model_type() == RHO_3D_READ) || (get_rhoinit(mgi) > 0) || (mgi == MMODELGRID));
    mg_associated_cells[mgi] += 1;
    assert_always(!(get_model_type() == RHO_3D_READ) || (mg_associated_cells[mgi] == 1) || (mgi == MMODELGRID));
  }

  int numnonemptycells = 0;
  for (int mgi = 0; mgi < get_npts_model(); mgi++)
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
      for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++)
      {
        const int z = decay::get_nuc_z(nucindex);
        const int a = decay::get_nuc_a(nucindex);
        set_modelinitradioabund(mgi, z, a, 0.);
      }
    }
  }

  printout("There are %d modelgrid cells with associated propagation cells\n", numnonemptycells);

  printout("[info] mem_usage: NLTE populations for all allocated cells occupy a total of %.1f MB\n", mem_usage_nltepops / 1024. / 1024.);
}


static void density_1d_read(void)
// Map 1D spherical model grid onto 3D propagation grid
{
  for (int cellindex = 0; cellindex < globals::ngrid; cellindex++)
  {
    const double radial_pos = get_cellradialpos(cellindex);
    const double vcell = radial_pos / globals::tmin;
    const double vmin = 0.;
    if (radial_pos < globals::rmax)
    {
      if (globals::grid_type == GRID_SPHERICAL1D)
      {
        set_cell_modelgridindex(cellindex, cellindex);
      }
      else
      {
        set_cell_modelgridindex(cellindex, 0);

        for (int mgi = 0; mgi < (get_npts_model() - 1); mgi++)
        {
          if (vout_model[mgi] < vcell)
          {
            set_cell_modelgridindex(cellindex, mgi + 1);
          }
        }
      }

      if (vout_model[get_cell_modelgridindex(cellindex)] >= vmin)
      {
        globals::modelgrid[get_cell_modelgridindex(cellindex)].initial_radial_pos += radial_pos;
      }
      else
      {
        set_cell_modelgridindex(cellindex, MMODELGRID);
      }
    }
    else
    {
      set_cell_modelgridindex(cellindex, MMODELGRID);
    }
  }
}


static void density_2d_read(void)
// Map 2D model grid onto 3D propagation grid
{
  for (int n = 0; n < globals::ngrid; n++)
  {
    const double radial_pos = get_cellradialpos(n);

    if (radial_pos < globals::rmax)
    {
      double dcen[3];
      for (int d = 0; d < 3; d++)
      {
        const double cellcoordmin = - globals::coordmax[d] + (2 * get_cellcoordpointnum(n, d) * globals::coordmax[d] / globals::ncoordgrid[0]);
        dcen[d] = cellcoordmin + (0.5 * wid_init(0));
      }

      set_cell_modelgridindex(n, 0);
      const double zcylindrical = dcen[2];
      dcen[2] = 0.0;
      const double rcylindrical = vec_len(dcen);

      // Grid is uniform so only need to search in 1d to get r and z positions

      int mkeep1 = 0;
      for (int m = 0; m < ncoord_model[0]; m++)
      {
        if (rcylindrical > (m * dcoord1 * globals::tmin/globals::t_model))
        {
          mkeep1 = m;
          // set_cell_modelgridindex(n, m + 1);
        }
      }

      int mkeep2 = 0;
      for (int m = 0; m < ncoord_model[1]; m++)
      {
        if (zcylindrical > (((m * dcoord2) * globals::tmin/globals::t_model) - globals::rmax))
        {
          mkeep2 = m;
          // set_cell_modelgridindex(n, m + 1);
        }
      }
      set_cell_modelgridindex(n, (mkeep2 * ncoord_model[0]) + mkeep1);
      globals::modelgrid[get_cell_modelgridindex(n)].initial_radial_pos += radial_pos;

      //renorm[mkeep]++;
    }
    else
    {
      set_cell_modelgridindex(n, MMODELGRID);
    }
  }
}


static void abundances_read(void)
{
  const bool threedimensional = (get_model_type() == RHO_3D_READ);

  /// Open the abundances file
  FILE *abundance_file = fopen_required("abundances.txt", "r");

  /// and process through the grid to read in the abundances per cell
  /// The abundance file should only contain information for non-empty
  /// cells. Its format must be cellnumber (integer), abundance for
  /// element Z=1 (float) up to abundance for element Z=30 (float)
  /// i.e. in total one integer and 30 floats.

  // loop over propagation cells for 3D models, or modelgrid cells
  const int ncount = threedimensional ? globals::ngrid : get_npts_model();
  for (int n = 0; n < ncount; n++)
  {
    const int mgi = threedimensional ? get_cell_modelgridindex(n) : n;

    assert_always(!feof(abundance_file));
    char line[2048] = "";
    assert_always(line == fgets(line, 2048, abundance_file));
    char *linepos = line;
    int offset = 0;

    int cellnumber = -1;
    assert_always(sscanf(linepos, "%d%n", &cellnumber, &offset) == 1);
    linepos += offset;

    if (cellnumber != n + 1)
    {
      printout("[fatal] %s: grid cell mismatch ... abort\n", __func__);
      printout("[fatal] n %d, cellnumber %d\n", n, cellnumber);
      abort();
    }

    // the abundances.txt file specifies the elemental mass fractions for each model cell
    // (or proportial to mass frac, e.g. element densities because they will be normalised anyway)
    // The abundances begin with hydrogen, helium, etc, going as far up the atomic numbers as required
    double normfactor = 0.;
    float abundances_in[150];
    for (int anumber = 1; anumber <= 150; anumber++)
    {
      abundances_in[anumber - 1] = 0.;
      const int itemsread = sscanf(linepos, "%g%n", &abundances_in[anumber - 1], &offset);
      linepos += offset;
      // printout("%d %d %d %g\n", cellnumber, anumber, itemsread, abundances_in[anumber - 1]);
      if (itemsread != 1)
      {
        assert_always(anumber > 1); // at least one element (hydrogen) should have been specified
        break;
      }

      assert_always(abundances_in[anumber - 1] >= 0.);
      normfactor += abundances_in[anumber - 1];
    }

    if (get_numassociatedcells(mgi) > 0)
    {
      if (threedimensional || normfactor <= 0.)
        normfactor = 1.;

      for (int element = 0; element < get_nelements(); element++)
      {
        ///now set the abundances (by mass) of included elements, i.e.
        ///read out the abundances specified in the atomic data file
        const int anumber = get_element(element);
        const float elemabundance = abundances_in[anumber - 1] / normfactor;
        assert_always(elemabundance >= 0.);
        globals::modelgrid[mgi].composition[element].abundance = elemabundance;

        // radioactive nuclide abundances should have already been set by read_??_model
        set_elem_stable_abund_from_total(mgi, anumber, elemabundance);
      }
    }
  }

  fclose(abundance_file);
}


static void read_2d3d_modelradioabundanceline(FILE * model_input, const int mgi, const bool keep)
{
  char line[2048] = "";
  if (line != fgets(line, 2048, model_input))
  {
    printout("Read failed on second line for cell %d\n", mgi);
    abort();
  }

  double f56ni_model = 0.;
  double f56co_model = 0.;
  double ffegrp_model = 0.;
  double f48cr_model = 0.;
  double f52fe_model = 0.;
  double f57ni_model = 0.;
  double f57co_model = 0.;
  const int items_read = sscanf(line, "%lg %lg %lg %lg %lg %lg %lg",
    &ffegrp_model, &f56ni_model, &f56co_model, &f52fe_model, &f48cr_model, &f57ni_model, &f57co_model);

  if (items_read == 5 || items_read == 7)
  {
    if (items_read == 10 && mgi == 0)
    {
      printout("Found Ni57 and Co57 abundance columns in model.txt\n");
    }

    // printout("mgi %d ni56 %g co56 %g fe52 %g cr48 %g ni57 %g co57 %g\n",
    //          mgi, f56ni_model, f56co_model, f52fe_model, f48cr_model, f57ni_model, f57co_model);

    if (keep)
    {
      globals::modelgrid[mgi].initradioabund = (float *) calloc(decay::get_num_nuclides(), sizeof(float));
      assert_always(globals::modelgrid[mgi].initradioabund != NULL);

      set_modelinitradioabund(mgi, 28, 56, f56ni_model);
      set_modelinitradioabund(mgi, 27, 56, f56co_model);
      set_modelinitradioabund(mgi, 28, 57, f57ni_model);
      set_modelinitradioabund(mgi, 27, 57, f57co_model);
      set_modelinitradioabund(mgi, 26, 52, f52fe_model);
      set_modelinitradioabund(mgi, 24, 48, f48cr_model);
      set_modelinitradioabund(mgi, 23, 48, 0.);

      set_ffegrp(mgi, ffegrp_model);
      //printout("mgi %d, control rho_init %g\n",mgi,get_rhoinit(mgi));
    }
  }
  else
  {
    printout("Unexpected number of values in model.txt. items_read = %d\n", items_read);
    printout("line: %s\n", line);
    abort();
  }
}


static void read_1d_model(void)
// Read in a 1D spherical model
{
  FILE *model_input = fopen_required("model.txt", "r");

  // 1st read the number of data points in the table of input model.
  int npts_model_in = 0;
  fscanf(model_input, "%d", &npts_model_in);
  set_npts_model(npts_model_in);
  ncoord_model[0] = npts_model_in;

  // Now read the time (in days) at which the model is specified.
  double t_model_days;
  fscanf(model_input, "%lg\n", &t_model_days);
  globals::t_model = t_model_days * DAY;

  // Now read in the lines of the model. Each line has 5 entries: the
  // cell number (integer) the velocity at outer boundary of cell (float),
  // the mass density in the cell (float), the abundance of Ni56 by mass
  // in the cell (float) and the total abundance of all Fe-grp elements
  // in the cell (float). For now, the last number is recorded but never
  // used.

  decay::init_nuclides();

  int mgi = 0;
  while (!feof(model_input))
  {
    char line[2048] = "";
    if (line != fgets(line, 2048, model_input))
    {
      // no more lines to read in
      break;
    }

    int cellnumberin;
    double vout_kmps;
    double log_rho;
    double f56ni_model = 0.;
    double f56co_model = 0.;
    double ffegrp_model = 0.;
    double f48cr_model = 0.;
    double f52fe_model = 0.;
    double f57ni_model = 0.;
    double f57co_model = 0.;
    const int items_read = sscanf(line, "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg",
                                   &cellnumberin, &vout_kmps, &log_rho, &ffegrp_model, &f56ni_model,
                                   &f56co_model, &f52fe_model, &f48cr_model, &f57ni_model, &f57co_model);

    if (items_read == 8 || items_read == 10)
    {
      assert_always(cellnumberin == mgi + 1);

      vout_model[mgi] = vout_kmps * 1.e5;

      const double rho_tmin = pow(10., log_rho) * pow(globals::t_model / globals::tmin, 3);
      set_rhoinit(mgi, rho_tmin);
      set_rho(mgi, rho_tmin);

      if (items_read == 10 && mgi == 0)
      {
        printout("Found Ni57 and Co57 abundance columns in model.txt\n");
      }
    }
    else
    {
      printout("Unexpected number of values in model.txt. items_read = %d\n", items_read);
      printout("line: %s\n", line);
      abort();
    }

    // printout("%d %g %g %g %g %g %g %g\n",
    //          cellnumin, vout_kmps, log_rho, ffegrp_model[n], f56ni_model[n],
    //          f56co_model[n], f52fe_model[n], f48cr_model[n]);
    // printout("   %lg %lg\n", f57ni_model[n], f57co_model[n]);
    globals::modelgrid[mgi].initradioabund = (float *) calloc(decay::get_num_nuclides(), sizeof(float));
    assert_always(globals::modelgrid[mgi].initradioabund != NULL);

    set_modelinitradioabund(mgi, 28, 56, f56ni_model);
    set_modelinitradioabund(mgi, 27, 56, f56co_model);
    set_modelinitradioabund(mgi, 28, 57, f57ni_model);
    set_modelinitradioabund(mgi, 27, 57, f57co_model);
    set_modelinitradioabund(mgi, 26, 52, f52fe_model);
    set_modelinitradioabund(mgi, 24, 48, f48cr_model);
    set_modelinitradioabund(mgi, 23, 48, 0.);
    set_ffegrp(mgi, ffegrp_model);

    mgi += 1;
    if (mgi == get_npts_model())
    {
      break;
    }
  }

  if (mgi != get_npts_model())
  {
    printout("ERROR in model.txt. Found %d only cells instead of %d expected.\n", mgi - 1, get_npts_model());
    abort();
  }

  fclose(model_input);

  globals::vmax = vout_model[get_npts_model() - 1];
}


static void read_2d_model(void)
// Read in a 2D axisymmetric spherical coordinate model
{
  FILE *model_input = fopen_required("model.txt", "r");

  // 1st read the number of data points in the table of input model.
  fscanf(model_input, "%d %d", &ncoord_model[0], &ncoord_model[1]);  // r and z (cylindrical polar)

  set_npts_model(ncoord_model[0] * ncoord_model[1]);

  // Now read the time (in days) at which the model is specified.
  double t_model_days;
  fscanf(model_input, "%lg", &t_model_days);
  globals::t_model = t_model_days * DAY;

  // Now read in globals::vmax (in cm/s)
  fscanf(model_input, "%lg\n", &globals::vmax);
  dcoord1 = globals::vmax * globals::t_model / ncoord_model[0]; //dr for input model
  dcoord2 = 2. * globals::vmax * globals::t_model / ncoord_model[1]; //dz for input model

  // Now read in the model. Each point in the model has two lines of input.
  // First is an index for the cell then its r-mid point then its z-mid point
  // then its total mass density.
  // Second is the total FeG mass, initial 56Ni mass, initial 56Co mass

  decay::init_nuclides();

  int mgi = 0;
  while (!feof(model_input))
  {
    char line[1024] = "";
    if (line != fgets(line, 1024, model_input))
    {
      // no more lines to read in
      break;
    }

    int cellnumin;
    float cell_r_in;
    float cell_z_in;
    double rho_tmodel;

    int items_read = sscanf(line, "%d %g %g %lg", &cellnumin, &cell_r_in, &cell_z_in, &rho_tmodel);
    assert_always(items_read == 4);

    const int ncoord1 = ((cellnumin - 1) % ncoord_model[0]);
    const double r_cylindrical = (ncoord1 + 0.5) * dcoord1;
    assert_always(fabs(cell_r_in / r_cylindrical - 1) < 1e-3);
    const int ncoord2 = ((cellnumin - 1) / ncoord_model[0]);
    const double z = -globals::vmax * globals::t_model + ((ncoord2 + 0.5) * dcoord2);
    assert_always(fabs(cell_z_in / z - 1) < 1e-3);

    assert_always(cellnumin == mgi + 1);

    const double rho_tmin = rho_tmodel * pow(globals::t_model / globals::tmin, 3);
    set_rhoinit(mgi, rho_tmin);
    set_rho(mgi, rho_tmin);

    read_2d3d_modelradioabundanceline(model_input, mgi, true);

    mgi++;
  }

  if (mgi != get_npts_model())
  {
    printout("ERROR in model.txt. Found %d only cells instead of %d expected.\n", mgi - 1, get_npts_model());
    abort();
  }

  fclose(model_input);
}


static void read_3d_model(void)
/// Subroutine to read in a 3-D model.
{
  FILE *model_input = fopen_required("model.txt", "r");

  /// 1st read the number of data points in the table of input model.
  /// This MUST be the same number as the maximum number of points used in the grid - if not, abort.
  int npts_model_in = 0;
  fscanf(model_input, "%d", &npts_model_in);
  if (npts_model_in > MMODELGRID)
  {
    printout("Too many points in input model. Abort. (%d > %d)\n", npts_model_in, MMODELGRID);
    abort();
  }

  ncoord_model[0] = ncoord_model[1] = ncoord_model[2] = round(pow(npts_model_in, 1/3.));
  assert_always(ncoord_model[0] * ncoord_model[1] * ncoord_model[2] == npts_model_in);

  // for a 3D input model, the progation cells will match the input cells exactly
  globals::ncoordgrid[0] = ncoord_model[0];
  globals::ncoordgrid[1] = ncoord_model[1];
  globals::ncoordgrid[2] = ncoord_model[2];
  globals::ngrid = npts_model_in;
  globals::grid_type = GRID_UNIFORM;

  /// Now read the time (in days) at which the model is specified.
  float t_model_days;
  fscanf(model_input, "%g", &t_model_days);
  globals::t_model = t_model_days * DAY;

  /// Now read in globals::vmax for the model (in cm s^-1).
  fscanf(model_input, "%lg\n", &globals::vmax);

  double xmax_tmodel = globals::vmax * globals::t_model;

  /// Now read in the lines of the model.
  min_den = 1.e99;

  // check if expected positions match in either xyz or zyx column order
  // set false if a problem is detected
  bool posmatch_xyz = true;
  bool posmatch_zyx = true;

  decay::init_nuclides();

  // mgi is the index to the model grid - empty cells are sent to MMODELGRID,
  // otherwise each input cell is one modelgrid cell
  int mgi = 0;
  int n = 0;
  while (!feof(model_input))
  {
    char line[1024] = "";
    if (line != fgets(line, 1024, model_input))
    {
      // no more lines to read in
      break;
    }

    int mgi_in;
    float cellpos_in[3];
    float rho_model;
    int items_read = sscanf(line, "%d %g %g %g %g", &mgi_in, &cellpos_in[0], &cellpos_in[1], &cellpos_in[2], &rho_model);
    assert_always(items_read == 5);
    //printout("cell %d, posz %g, posy %g, posx %g, rho %g, rho_init %g\n",dum1,dum3,dum4,dum5,rho_model,rho_model* pow( (t_model/globals::tmin), 3.));

    assert_always(mgi_in == n + 1);

    // cell coordinates in the 3D model.txt file are sometimes reordered by the scaling script
    // however, the cellindex always should increment X first, then Y, then Z

    for (int axis = 0; axis < 3; axis++)
    {
      const double cellwidth = 2 * xmax_tmodel / globals::ncoordgrid[axis];
      const double cellpos_expected = - xmax_tmodel + cellwidth * get_cellcoordpointnum(n, axis);
      // printout("n %d coord %d expected %g found %g rmax %g get_cellcoordpointnum(n, axis) %d globals::ncoordgrid %d\n",
      // n, axis, cellpos_expected, cellpos_in[axis], xmax_tmodel, get_cellcoordpointnum(n, axis), globals::ncoordgrid[axis]);
      if (fabs(cellpos_expected - cellpos_in[axis]) > 0.5 * cellwidth)
      {
        posmatch_xyz = false;
      }
      if (fabs(cellpos_expected - cellpos_in[2 - axis]) > 0.5 * cellwidth)
      {
        posmatch_zyx = false;
      }
    }

    if (rho_model < 0)
    {
      printout("negative input density %g %d\n", rho_model, n);
      abort();
    }

    if (mgi > MMODELGRID - 1)
    {
      printout("3D model wants more modelgrid cells than MMODELGRID. Abort.\n");
      abort();
    }

    const bool keepcell = (rho_model > 0);
    if (keepcell)
    {
      set_cell_modelgridindex(n, mgi);
      const double rho_tmin = rho_model * pow((globals::t_model / globals::tmin), 3.);
      //printout("mgi %d, helper %g\n",mgi,helper);
      set_rhoinit(mgi, rho_tmin);
      //printout("mgi %d, rho_init %g\n",mgi,get_rhoinit(mgi));
      set_rho(mgi, rho_tmin);

      if (min_den > rho_model)
      {
        min_den = rho_model;
      }
    }
    else
    {
      set_cell_modelgridindex(n, MMODELGRID);
    }

    read_2d3d_modelradioabundanceline(model_input, mgi, keepcell);
    if (keepcell)
    {
      mgi++;
    }

    n++;
  }
  if (n != npts_model_in)
  {
    printout("ERROR in model.txt. Found %d cells instead of %d expected.\n", n, npts_model_in);
    abort();
  }

  assert_always(posmatch_zyx ^ posmatch_xyz);  // xor because if both match then probably an infinity occurred
  if (posmatch_xyz)
  {
    printout("Cell positions in model.txt are consistent with calculated values when x-y-z column order is used.\n");
  }
  if (posmatch_zyx)
  {
    printout("Cell positions in model.txt are consistent with calculated values when z-y-x column order is used.\n");
  }

  printout("min_den %g [g/cm3]\n", min_den);
  printout("Effectively used model grid cells %d\n", mgi);

  /// Now, set actual size of the modelgrid to the number of non-empty cells.

  set_npts_model(mgi);

  fclose(model_input);
}


static void calc_totmassradionuclides(void)
{
  mtot = 0.;
  mfeg = 0.;

  assert_always(totmassradionuclide == NULL);
  totmassradionuclide = (double *) calloc(decay::get_num_nuclides(), sizeof(double));
  assert_always(totmassradionuclide != NULL);

  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++)
    totmassradionuclide[nucindex] = 0.;

  int n1 = 0;
  for (int mgi = 0; mgi < get_npts_model(); mgi++)
  {
    double cellvolume = 0.;
    if (get_model_type() == RHO_1D_READ)
    {
      const double v_inner = (mgi == 0) ? 0. : vout_model[mgi - 1];
      // mass_in_shell = rho_model[mgi] * (pow(vout_model[mgi], 3) - pow(v_inner, 3)) * 4 * PI * pow(t_model, 3) / 3.;
      cellvolume = (pow(vout_model[mgi], 3) - pow(v_inner, 3)) * 4 * PI * pow(globals::tmin, 3) / 3.;
    }
    else if (get_model_type() == RHO_2D_READ)
    {
      cellvolume = pow(globals::tmin / globals::t_model, 3) * ((2 * n1) + 1) * PI * dcoord2 * pow(dcoord1, 2.);
      n1++;
      if (n1 == ncoord_model[0])
      {
        n1 = 0;
      }
    }
    else if (get_model_type() == RHO_3D_READ)
    {
      /// Assumes cells are cubes here - all same volume.
      cellvolume = pow((2 * globals::vmax * globals::tmin), 3.) / (globals::ncoordgrid[0] * globals::ncoordgrid[1] * globals::ncoordgrid[2]);
    }
    else
    {
      printout("Unknown model type %d in function %s\n", get_model_type(), __func__);
      abort();
    }

    const double mass_in_shell = get_rhoinit(mgi) * cellvolume;

    mtot += mass_in_shell;

    for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++)
    {
      const int z = decay::get_nuc_z(nucindex);
      const int a = decay::get_nuc_a(nucindex);
      if (z > 0)  // skips FAKE_GAM_LINE_ID
      {
        totmassradionuclide[nucindex] += mass_in_shell * get_modelinitradioabund(mgi, z, a);
      }
    }

    mfeg += mass_in_shell * get_ffegrp(mgi);
  }


  printout("Masses / Msun:    Total: %9.3e  56Ni: %9.3e  56Co: %9.3e  52Fe: %9.3e  48Cr: %9.3e\n",
           mtot / MSUN, get_totmassradionuclide(28, 56) / MSUN,
           get_totmassradionuclide(27, 56) / MSUN, get_totmassradionuclide(26, 52) / MSUN,
           get_totmassradionuclide(24, 48) / MSUN);
  printout("Masses / Msun: Fe-group: %9.3e  57Ni: %9.3e  57Co: %9.3e\n",
           mfeg / MSUN, get_totmassradionuclide(28, 57) / MSUN, get_totmassradionuclide(27, 57) / MSUN);
}


void read_ejecta_model(enum model_types model_type)
{
  switch (model_type)
  {
    case RHO_UNIFORM:
    {
      assert_always(false); // needs to be reimplemented using spherical coordinate mode
      break;
    }

    case RHO_1D_READ:
    {
      printout("Read 1D model\n");
      read_1d_model();
      break;
    }

    case RHO_2D_READ:
    {
      printout("Read 2D model\n");

      read_2d_model();
      break;
    }

    case RHO_3D_READ:
    {
      printout("Read 3D model\n");

      read_3d_model();
      break;
    }

    default:
    {
      printout("Unknown model. Abort.\n");
      abort();
    }
  }

  printout("npts_model: %d\n", get_npts_model());
  globals::rmax = globals::vmax * globals::tmin;
  printout("vmax %g\n", globals::vmax);
  printout("tmin %g\n", globals::tmin);
  printout("rmax %g\n", globals::rmax);

  globals::coordmax[0] = globals::coordmax[1] = globals::coordmax[2] = globals::rmax;

  calc_totmassradionuclides();
}


static void read_grid_restart_data(const int timestep)
{
  char filename[100];
  sprintf(filename, "gridsave_ts%d.tmp", timestep);

  printout("READIN GRID SNAPSHOT from %s\n", filename);
  FILE *gridsave_file = fopen_required(filename, "r");

  int ntstep_in = -1;
  fscanf(gridsave_file, "%d ", &ntstep_in);
  assert_always(ntstep_in == globals::ntstep);

  int nprocs_in = -1;
  fscanf(gridsave_file, "%d ", &nprocs_in);
  assert_always(nprocs_in == globals::nprocs);

  int nthreads_in = -1;
  fscanf(gridsave_file, "%d ", &nthreads_in);
  assert_always(nthreads_in == get_num_threads());

  for (int nts = 0; nts < globals::ntstep; nts++)
  {
    fscanf(gridsave_file, "%lg %lg ", &globals::time_step[nts].gamma_dep, &globals::time_step[nts].positron_dep);
  }

  int timestep_in;
  fscanf(gridsave_file, "%d ", &timestep_in);
  assert_always(timestep_in = timestep);

  for (int mgi = 0; mgi < get_npts_model(); mgi++)
  {
    int mgi_in;
    float T_R;
    float T_e;
    float W;
    float T_J;
    fscanf(gridsave_file, "%d %g %g %g %g %hd %lg",
           &mgi_in, &T_R, &T_e, &W, &T_J,
           &globals::modelgrid[mgi].thick, &globals::rpkt_emiss[mgi]);

    if (mgi_in != mgi)
    {
      printout("[fatal] read_grid_restart_data: cell mismatch in reading input gridsave.dat ... abort\n");
      printout("[fatal] read_grid_restart_data: read cellnumber %d, expected cellnumber %d\n",mgi_in,mgi);
      abort();
    }

    assert_always(T_R >= 0);
    assert_always(T_e >= 0);
    assert_always(W >= 0);
    assert_always(T_J >= 0);
    assert_always(globals::rpkt_emiss[mgi] >= 0);

    set_TR(mgi, T_R);
    set_Te(mgi, T_e);
    set_W(mgi, W);
    set_TJ(mgi, T_J);

    #ifndef FORCE_LTE
      #if (!NO_LUT_PHOTOION)
        for (int element = 0; element < get_nelements(); element++)
        {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions; ion++)
          {
            const int estimindex = mgi * get_nelements() * get_max_nions() + element * get_max_nions() + ion;
            fscanf(gridsave_file, " %lg %lg", &globals::corrphotoionrenorm[estimindex], &globals::gammaestimator[estimindex]);
          }
        }
      #endif
    #endif
  }

  // the order of these calls is very important!
  radfield::read_restart_data(gridsave_file);
  nonthermal::read_restart_data(gridsave_file);
  nltepop_read_restart_data(gridsave_file);
  fclose(gridsave_file);
}


void write_grid_restart_data(const int timestep)
{
  char filename[100];
  sprintf(filename, "gridsave_ts%d.tmp", timestep);

  const time_t sys_time_start_write_restart = time(NULL);
  printout("Write grid restart data to %s...", filename);

  FILE *gridsave_file = fopen_required(filename, "w");

  fprintf(gridsave_file, "%d ", globals::ntstep);
  fprintf(gridsave_file, "%d ", globals::nprocs);
  fprintf(gridsave_file, "%d ", get_num_threads());

  for (int i = 0; i < globals::ntstep; i++)
  {
    fprintf(gridsave_file, "%lg %lg ", globals::time_step[i].gamma_dep, globals::time_step[i].positron_dep);
  }

  fprintf(gridsave_file, "%d ", timestep);

  for (int mgi = 0; mgi < get_npts_model(); mgi++)
  {
    const bool nonemptycell = (get_numassociatedcells(mgi) > 0);

    if (nonemptycell)
    {
      fprintf(gridsave_file, "%d %g %g %g %g %hd %lg",
              mgi, get_TR(mgi), get_Te(mgi), get_W(mgi), get_TJ(mgi),
              globals::modelgrid[mgi].thick, globals::rpkt_emiss[mgi]);
    }
    else
    {
      fprintf(gridsave_file, "%d %g %g %g %g %d %lg", mgi, 0., 0., 0., 0., 0, 0.);
    }

    #ifndef FORCE_LTE
      #if (!NO_LUT_PHOTOION)
        for (int element = 0; element < get_nelements(); element++)
        {
          const int nions = get_nions(element);
          for (int ion = 0; ion < nions; ion++)
          {
            const int estimindex = mgi * get_nelements() * get_max_nions() + element * get_max_nions() + ion;
            fprintf(gridsave_file, " %lg %lg",
                    (nonemptycell ? globals::corrphotoionrenorm[estimindex] : 0.),
                    (nonemptycell ? globals::gammaestimator[estimindex] : 0.));
          }
        }
      #endif
    #endif
    fprintf(gridsave_file,"\n");
  }

  // the order of these calls is very important!
  radfield::write_restart_data(gridsave_file);
  nonthermal::write_restart_data(gridsave_file);
  nltepop_write_restart_data(gridsave_file);
  fclose(gridsave_file);
  printout("done in %ld seconds.\n", time(NULL) - sys_time_start_write_restart);
}


static void assign_temperature(void)
/// Routine for assigning temperatures to the grid cells at the start of the simulation.
{
  /// For a simulation started from scratch we estimate the initial temperatures

  /// We assume that for early times the material is so optically thick, that
  /// all the radiation is trapped in the cell it originates from. This
  /// means furthermore LTE, so that both temperatures can be evaluated
  /// according to the local energy density resulting from the 56Ni decay.
  /// The dilution factor is W=1 in LTE.

  const double tstart = globals::time_step[0].mid;

  const double e_ni56 = decay::nucdecayenergy(28, 56);
  const double t_ni56 = decay::meanlife(28, 56);
  const double e_co56 = decay::nucdecayenergy(27, 56);
  const double t_co56 = decay::meanlife(27, 56);

  const double factor56ni = 1. / 56 / MH * (-1. / (tstart * (- t_co56 + t_ni56)))
    * (- e_ni56 * exp(- tstart / t_ni56) * tstart * t_co56 - e_ni56 * exp(- tstart / t_ni56) * t_ni56 * t_co56
       + e_ni56 * exp(- tstart / t_ni56) * tstart * t_ni56 + pow(t_ni56, 2) * e_ni56 * exp(- tstart / t_ni56)
       - t_co56 * tstart * e_co56 * exp(- tstart / t_co56) - pow(t_co56, 2) * e_co56 * exp(- tstart / t_co56)
       + e_co56 * tstart * t_ni56 * exp(- tstart / t_ni56) + pow(t_ni56, 2) * e_co56 * exp(- tstart / t_ni56)
       + e_ni56 * t_co56 * t_ni56 - e_ni56 * pow(t_ni56, 2) - pow(t_ni56, 2) * e_co56 + e_co56 * pow(t_co56, 2));

  const double factor56co = 1. / 56 / MH * (1. / (tstart * t_co56))
    * (t_co56 * tstart * e_co56 * exp(- tstart / t_co56) + pow(t_co56, 2) * e_co56 * exp(- tstart / t_co56));

  const double e_ni57 = decay::nucdecayenergy(28, 57);
  const double t_ni57 = decay::meanlife(28, 57);
  const double e_co57 = decay::nucdecayenergy(27, 57);
  const double t_co57 = decay::meanlife(27, 57);

  const double factor57ni = 1. / 57 / MH * (-1. / (tstart * (- t_co57 + t_ni57)))
    * (- e_ni57 * exp(- tstart / t_ni57) * tstart * t_co57 - e_ni57 * exp(- tstart / t_ni57) * t_ni57 * t_co57
       + e_ni57 * exp(- tstart / t_ni57) * tstart * t_ni57 + pow(t_ni57, 2) * e_ni57 * exp(- tstart / t_ni57)
       - t_co57 * tstart * e_co57 * exp(- tstart / t_co57) - pow(t_co57, 2) * e_co57 * exp(- tstart / t_co57)
       + e_co57 * tstart * t_ni57 * exp(- tstart / t_ni57) + pow(t_ni57, 2) * e_co57 * exp(- tstart / t_ni57)
       + e_ni57 * t_co57 * t_ni57 - e_ni57 * pow(t_ni57, 2) - pow(t_ni57, 2) * e_co57 + e_co57 * pow(t_co57, 2));

 const double e_fe52 = decay::nucdecayenergy(26, 52);
 const double t_fe52 = decay::meanlife(26, 52);
 const double e_mn52 = decay::nucdecayenergy(25, 52);
 const double t_mn52 = decay::meanlife(25, 52);

  const double factor52fe = 1. / 52 / MH * (-1. / (tstart * (- t_mn52 + t_fe52)))
    * (- e_fe52 * exp(- tstart / t_fe52) * tstart * t_mn52 - e_fe52 * exp(- tstart / t_fe52) * t_fe52 * t_mn52
       + e_fe52 * exp(- tstart / t_fe52) * tstart * t_fe52 + pow(t_fe52, 2) * e_fe52 * exp(- tstart / t_fe52)
       - t_mn52 * tstart * e_mn52 * exp(- tstart / t_mn52) - pow(t_mn52, 2) * e_mn52 * exp(- tstart / t_mn52)
       + e_mn52 * tstart * t_fe52 * exp(- tstart / t_fe52) + pow(t_fe52, 2) * e_mn52 * exp(- tstart / t_fe52)
       + e_fe52 * t_mn52 * t_fe52 - e_fe52 * pow(t_fe52, 2) - pow(t_fe52, 2) * e_mn52 + e_mn52 * pow(t_mn52, 2));

  const double e_cr48 = decay::nucdecayenergy(24, 48);
  const double t_cr48 = decay::meanlife(24, 48);
  const double e_v48 = decay::nucdecayenergy(23, 48);
  const double t_v48 = decay::meanlife(23, 48);

  const double factor48cr = 1. / 48 / MH * (-1. / (tstart * (- t_v48 + t_cr48)))
    * (- e_cr48 * exp(- tstart / t_cr48) * tstart * t_v48 - e_cr48 * exp(- tstart / t_cr48) * t_cr48 * t_v48
       + e_cr48 * exp(- tstart / t_cr48) * tstart * t_cr48 + pow(t_cr48, 2) * e_cr48 * exp(- tstart / t_cr48)
       - t_v48 * tstart * e_v48 * exp(- tstart / t_v48) - pow(t_v48, 2) * e_v48 * exp(- tstart / t_v48)
       + e_v48 * tstart * t_cr48 * exp(- tstart / t_cr48) + pow(t_cr48, 2) * e_v48 * exp(- tstart / t_cr48)
       + e_cr48 * t_v48 * t_cr48 - e_cr48 * pow(t_cr48, 2) - pow(t_cr48, 2) * e_v48 + e_v48 * pow(t_v48, 2));

  //factor56ni = CLIGHT/4/STEBO * nucdecayenergy(28, 56)/56/MH;
  /// This works only for the inbuilt Lucy model
  //factor56ni = CLIGHT/4/STEBO * 3*mtot/4/PI * nucdecayenergy(28, 56)/56/MH  / pow(vmax,3);
  for (int mgi = 0; mgi < get_npts_model(); mgi++)
  {
    // alternative, not correct because it neglects expansion factor
    // const double decayedenergy_per_mass = get_decayedenergy_per_ejectamass(mgi, tstart);
    //
    // double T_initial = pow(CLIGHT / 4 / STEBO  * pow(globals::tmin / tstart, 3) * get_rhoinit(mgi) * decayedenergy_per_mass, 1. / 4.);

    double T_initial = pow(CLIGHT / 4 / STEBO  * pow(globals::tmin / tstart, 3) * get_rhoinit(mgi) * (
         (factor56ni * get_modelinitradioabund(mgi, 28, 56)) +
         (factor56co * get_modelinitradioabund(mgi, 27, 56)) +
         (factor57ni * get_modelinitradioabund(mgi, 28, 57)) +
         // (factor57co * get_modelinitradioabund(mgi, 27, 57)) +
         (factor52fe * get_modelinitradioabund(mgi, 26, 52)) +
         (factor48cr * get_modelinitradioabund(mgi, 23, 48))), 1. / 4.);

    // printout("mgi %d: T_initial %g K tmin %g tstart %g rhoinit %g X_56Ni %g X_52Fe %g X_48cr %g\n",
    //          mgi, T_initial, globals::tmin, tstart, get_rhoinit(mgi), get_modelinitradioabund(mgi, 28, 56), get_modelinitradioabund(mgi, 26, 52), get_modelinitradioabund(mgi, 24, 48));
    if (T_initial < MINTEMP)
    {
      printout("mgi %d: T_initial of %g is below MINTEMP %g K, setting to MINTEMP.\n", mgi, T_initial, MINTEMP);
      T_initial = MINTEMP;
    }
    else if (T_initial > MAXTEMP)
    {
      printout("mgi %d: T_initial of %g is above MAXTEMP %g K, setting to MAXTEMP.\n", mgi, T_initial, MAXTEMP);
      T_initial = MAXTEMP;
    }

    set_Te(mgi, T_initial);
    set_TJ(mgi, T_initial);
    set_TR(mgi, T_initial);

    set_W(mgi, 1.);
  }
}


static void uniform_grid_setup(void)
/// Routine for doing a uniform cuboidal grid.
{

  /// Set grid size for uniform xyz grid
  if (get_model_type() == RHO_3D_READ)
  {
    // if we used in a 3D ejecta model, the propagation grid must match the input grid exactly
    globals::ncoordgrid[0] = ncoord_model[0];
    globals::ncoordgrid[1] = ncoord_model[1];
    globals::ncoordgrid[2] = ncoord_model[2];

    // in case the user specified a grid size, we should ensure that it matches
    #ifdef CUBOID_NCOORDGRID_X
    assert(globals::ncoordgrid[0] == CUBOID_NCOORDGRID_X);
    #endif
    #ifdef CUBOID_NCOORDGRID_Y
    assert(globals::ncoordgrid[1] == CUBOID_NCOORDGRID_Y);
    #endif
    #ifdef CUBOID_NCOORDGRID_Z
    assert(globals::ncoordgrid[2] == CUBOID_NCOORDGRID_Z);
    #endif
  }
  else
  {
    #ifdef CUBOID_NCOORDGRID_X
    globals::ncoordgrid[0] = CUBOID_NCOORDGRID_X;
    #else
    globals::ncoordgrid[0] = 50;
    #endif
    #ifdef CUBOID_NCOORDGRID_Y
    globals::ncoordgrid[1] = CUBOID_NCOORDGRID_Y;
    #else
    globals::ncoordgrid[1] = 50;
    #endif
    #ifdef CUBOID_NCOORDGRID_Z
    globals::ncoordgrid[2] = CUBOID_NCOORDGRID_Z;
    #else
    globals::ncoordgrid[2] = 50;
    #endif
  }

  // artis assumes in some places that the cells are cubes, not cubioids
  assert_always(globals::ncoordgrid[0] == globals::ncoordgrid[1]);
  assert_always(globals::ncoordgrid[0] == globals::ncoordgrid[2]);

  globals::ngrid = globals::ncoordgrid[0] * globals::ncoordgrid[1] * globals::ncoordgrid[2];
  assert_always(globals::ngrid <= MGRID);

  globals::coordlabel[0] = 'X';
  globals::coordlabel[1] = 'Y';
  globals::coordlabel[2] = 'Z';
  int nxyz[3] = {0, 0, 0};
  for (int n = 0; n < globals::ngrid; n++)
  {
    for (int axis = 0; axis < 3; axis++)
    {
      assert_always(nxyz[axis] == get_cellcoordpointnum(n, axis));
      cell[n].pos_init[axis] = - globals::coordmax[axis] + (2 * nxyz[axis] * globals::coordmax[axis] / globals::ncoordgrid[axis]);
      // cell[n].xyz[axis] = nxyz[axis];
    }

    assert_always(n == nxyz[2] * globals::ncoordgrid[1] * globals::ncoordgrid[2] + nxyz[1] * globals::ncoordgrid[0] + nxyz[0]);

    nxyz[0]++;  // increment x coordinate
    if (nxyz[0] == globals::ncoordgrid[0])
    {
      nxyz[0] = 0;
      nxyz[1]++;  // increment y coordinate
    }
    if (nxyz[1] == globals::ncoordgrid[1])
    {
      nxyz[1] = 0;
      nxyz[2]++;  // increment z coordinate
    }
  }
}


static void spherical1d_grid_setup(void)
{
  assert_always(get_model_type() == RHO_1D_READ);
  globals::coordlabel[0] = 'r';
  globals::coordlabel[1] = '_';
  globals::coordlabel[2] = '_';

  globals::ncoordgrid[0] = get_npts_model();
  globals::ncoordgrid[1] = 1;
  globals::ncoordgrid[2] = 1;

  globals::ngrid = globals::ncoordgrid[0] * globals::ncoordgrid[1] * globals::ncoordgrid[2];
  assert_always(globals::ngrid <= MGRID);

  globals::coordmax[0] = globals::rmax;
  globals::coordmax[1] = 0.;
  globals::coordmax[2] = 0.;

  assert_always(globals::ngrid <= MGRID);

  // in this mode, cellindex and modelgridindex are the same thing
  for (int cellindex = 0; cellindex < get_npts_model(); cellindex++)
  {
    const double v_inner = cellindex > 0 ? vout_model[cellindex - 1] : 0.;
    set_cell_modelgridindex(cellindex, cellindex);
    cell[cellindex].pos_init[0] = v_inner * globals::tmin;
    cell[cellindex].pos_init[1] = 0.;
    cell[cellindex].pos_init[2] = 0.;
  }
}


void grid_init(int my_rank)
/// Initialises the propagation grid cells and associates them with modelgrid cells
{
  /// Start by checking that the number of grid cells is okay */
  //ngrid = globals::ncoordgrid[0] * globals::ncoordgrid[1] * globals::ncoordgrid[2]; ///Moved to input.c
  //if (ngrid > MGRID)
  //{
  //  printout("[fatal] grid_init: Error: too many grid cells. Abort.");
  //  abort();
  //}
  for (int n = 0; n <= MMODELGRID; n++)
  {
    globals::modelgrid[n].initial_radial_pos = 0;
  }

  /// Select grid type
  #ifdef GRID_TYPE
  globals::grid_type = GRID_TYPE;
  #else
  globals::grid_type = GRID_UNIFORM;
  #endif

  /// The cells will be ordered by x then y, then z. Call a routine that
  /// sets up the initial positions and widths of the cells.
  char grid_type_name[256] = "";
  if (globals::grid_type == GRID_UNIFORM)
  {
    uniform_grid_setup();
    strcpy(grid_type_name, "uniform cuboidal");
  }
  else if (globals::grid_type == GRID_SPHERICAL1D)
  {
    spherical1d_grid_setup();
    strcpy(grid_type_name, "spherical");
  }
  else
  {
    printout("[fatal] grid_init: Error: Unknown grid type. Abort.");
    abort();
  }

  printout("propagation grid: %d-dimensional %s\n", get_ngriddimensions(), grid_type_name);

  for (int d = 0; d < get_ngriddimensions(); d++)
  {
    printout("    coordinate %d '%c': cells have %d position values\n", d, globals::coordlabel[d], globals::ncoordgrid[d]);
  }
  printout("    total propagration cells: %d\n", globals::ngrid);

  /// Now set up the density in each cell.

  // Calculate the critical opacity at which opacity_case 3 switches from a
  // regime proportional to the density to a regime independent of the density
  // This is done by solving for tau_sobolev == 1
  // tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/nucmass(28, 56) * 3000e-8 * globals::time_step[m].mid;
  globals::rho_crit = ME * CLIGHT * decay::nucmass(28, 56) / (PI * QE * QE * globals::rho_crit_para * 3000e-8 * globals::tmin);
  printout("grid_init: rho_crit = %g [g/cm3]\n", globals::rho_crit);

  if (get_model_type() == RHO_1D_READ)
  {
    density_1d_read();
  }
  else if (get_model_type() == RHO_2D_READ)
  {
    density_2d_read();
  }
  else if (get_model_type() == RHO_3D_READ)
  {
    assert_always(globals::grid_type == GRID_UNIFORM);
    // propagation grid must match the input model grid exactly for 3D models
    assert_always(ncoord_model[0] == globals::ncoordgrid[0]);
    assert_always(ncoord_model[1] == globals::ncoordgrid[1]);
    assert_always(ncoord_model[2] == globals::ncoordgrid[2]);

    for (int n = 0; n < globals::ngrid; n++)
    {
      globals::modelgrid[get_cell_modelgridindex(n)].initial_radial_pos = get_cellradialpos(n);
    }
    // cells with rho > 0 are allocated by the above function
  }
  else
  {
    printout("[fatal] grid_init: Error: Unknown density type. Abort.");
    abort();
  }

  printout("allocate_nonemptymodelcells\n");
  allocate_nonemptymodelcells();
  printout("calculate_kappagrey\n");
  calculate_kappagrey();
  printout("abundances_read\n");
  abundances_read();

  radfield::init(my_rank);
  nonthermal::init(my_rank);

  /// and assign a temperature to the cells
  if (globals::simulation_continued_from_saved)
  {
    /// For continuation of an existing simulation we read the temperatures
    /// at the end of the simulation and write them to the grid.
    read_grid_restart_data(globals::itstep);
  }
  else
  {
    assign_temperature();
  }

  // scale up the radioactive abundances to account for the missing masses in
  // the model cells that are not associated with any propagation cells
  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++)
  {
    const int z = decay::get_nuc_z(nucindex);
    const int a = decay::get_nuc_a(nucindex);
    if (totmassradionuclide[nucindex] <= 0)
      continue;
    double totmassradionuclide_actual = 0.;
    for (int mgi = 0; mgi < get_npts_model(); mgi++)
    {
      totmassradionuclide_actual += get_modelinitradioabund(mgi, z, a) * get_rhoinit(mgi) * vol_init_modelcell(mgi);
    }
    if (totmassradionuclide_actual >= 0.)
    {
      const double ratio = totmassradionuclide[nucindex] / totmassradionuclide_actual;
      // printout("nuclide %d ratio %g\n", nucindex, ratio);
      for (int mgi = 0; mgi < get_npts_model(); mgi++)
      {
        if (get_numassociatedcells(mgi) > 0)
        {
          const double prev_abund = get_modelinitradioabund(mgi, z, a);
          const double new_abund = prev_abund * ratio;
          set_modelinitradioabund(mgi, z, a, new_abund);
        }
      }
    }
  }
}


double get_totmassradionuclide(const int z, const int a)
{
  return totmassradionuclide[decay::get_nuc_index(z, a)];
}