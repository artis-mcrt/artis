#include "sn3d.h"
#include "atomic.h"
#include "grid.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "decay.h"
#include "radfield.h"
#include "rpkt.h"
#include "vectors.h"


__managed__ enum model_types model_type = RHO_1D_READ;
__managed__ int npts_model = 0; // number of points in 1-D input model

__managed__ CELL cell[MGRID + 1];

static long mem_usage_nltepops = 0;

static __managed__ int mg_associated_cells[MMODELGRID + 1];

__managed__ double totmassradionuclide[RADIONUCLIDE_COUNT]; /// total mass of each radionuclide in the ejecta

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
      const double v_inner = modelgridindex > 0 ? globals::vout_model[modelgridindex - 1] : 0.;
      return (globals::vout_model[modelgridindex] - v_inner) * globals::tmin;
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
      return 4./3. * PI * (pow(globals::tmin * globals::vout_model[modelgridindex], 3) - pow(globals::tmin * (modelgridindex > 0 ? globals::vout_model[modelgridindex - 1] : 0.), 3));

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
  assert(x >= 0);
  assert(x <= 1.);
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
  assert_testmodeonly(new_modelgridindex < npts_model);
  cell[cellindex].modelgridindex = new_modelgridindex;
}


__host__ __device__
int get_numassociatedcells(const int modelgridindex)
// number of propagation cells associated with each modelgrid cell
{
  return mg_associated_cells[modelgridindex];
}


__host__ __device__
float get_modelinitradioabund(const int modelgridindex, const enum radionuclides nuclide_type)
{
  // this function replaces get_f56ni(mgi), get_fco56(mgi), etc.

  assert(nuclide_type < RADIONUCLIDE_COUNT);

  return globals::modelgrid[modelgridindex].initradioabund[nuclide_type];
}


__host__ __device__
void set_modelinitradioabund(const int modelgridindex, const enum radionuclides nuclide_type, const float abund)
{
  assert(abund >= 0.);
  assert(abund <= 1.);
  // this function replaces set_f56ni(mgi), set_fco56(mgi), etc.

  if ((nuclide_type == FAKE_GAM_LINE_ID && abund > 0) || nuclide_type >= RADIONUCLIDE_COUNT)
  {
    printout("set_initfracradio called with invalid nuclide_type type %d\n", nuclide_type);
    abort();
  }

  globals::modelgrid[modelgridindex].initradioabund[nuclide_type] = abund;
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
  // store the stable mass fraction of an element given the total element mass fraction
  // by subtracting the abundances of radioactive isotopes
  // if the element Z=anumber has no specific stable abundance variable then the function does nothing
  switch (anumber)
  {
    case 28:
      globals::modelgrid[mgi].fnistable = fmax(0.,
        elemabundance - get_modelinitradioabund(mgi, NUCLIDE_NI56) - get_modelinitradioabund(mgi, NUCLIDE_NI57));
      break;

    case 27:
      globals::modelgrid[mgi].fcostable = fmax(0.,
        elemabundance - get_modelinitradioabund(mgi, NUCLIDE_CO56) - get_modelinitradioabund(mgi, NUCLIDE_CO57));
      break;

    case 26:
      globals::modelgrid[mgi].ffestable = fmax(0., elemabundance - get_modelinitradioabund(mgi, NUCLIDE_FE52));
      break;

    case 25:
      globals::modelgrid[mgi].fmnstable = fmax(0., elemabundance);
      break;

    case 24:
      globals::modelgrid[mgi].fcrstable = fmax(0., elemabundance - get_modelinitradioabund(mgi, NUCLIDE_CR48));
      break;

    case 23:
      globals::modelgrid[mgi].fvstable = fmax(0., elemabundance);
      break;

    case 22:
      globals::modelgrid[mgi].ftistable = fmax(0., elemabundance);
      break;
  }
}


__host__ __device__
double get_cellradialpos(const int cellindex)
{
  if (globals::grid_type == GRID_SPHERICAL1D)
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
        set_kappagrey(mgi, ((0.9 * get_ffegrp(mgi)) + 0.1) * GREY_OP / ((0.9 * globals::mfeg / globals::mtot) + 0.1));
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
        set_kappagrey(mgi, ((0.9 * get_ffegrp(mgi)) + 0.1) * GREY_OP / ((0.9 *  globals::mfeg / globals::mtot) + 0.1));
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
  assert(globals::modelgrid[modelgridindex].elements_uppermost_ion != NULL);

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


static void allocate_nonemptycells(void)
{
  printout("mem_usage: the modelgrid array occupies %.1f MB\n", sizeof(globals::modelgrid) / 1024. / 1024.);
  mem_usage_nltepops = 0;
  /// This is the placeholder for empty cells. Temperatures must be positive
  /// as long as ff opacities are calculated.
  set_rhoinit(MMODELGRID, 0.);
  set_rho(MMODELGRID, 0.);
  set_nne(MMODELGRID, 0.);
  set_ffegrp(MMODELGRID, 0.);
  for (int isoint = 0; isoint < RADIONUCLIDE_COUNT; isoint++)
  {
    const enum radionuclides iso = (enum radionuclides) isoint;
    set_modelinitradioabund(MMODELGRID, iso, 0.);
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
    assert(!(get_model_type() == RHO_3D_READ) || (get_rhoinit(mgi) > 0) || (mgi == MMODELGRID));
    mg_associated_cells[mgi] += 1;
    assert(!(get_model_type() == RHO_3D_READ) || (mg_associated_cells[mgi] == 1) || (mgi == MMODELGRID));
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
      for (int isoint = 0; isoint < RADIONUCLIDE_COUNT; isoint++)
      {
        const enum radionuclides iso = (enum radionuclides) isoint;
        set_modelinitradioabund(mgi, iso, 0.);
      }
    }
  }

  printout("There are %d modelgrid cells with associated propagation cells\n", numnonemptycells);

  printout("mem_usage: NLTE populations for all allocated cells occupy a total of %.1f MB\n", mem_usage_nltepops / 1024. / 1024.);
}


static void density_1d_read(void)
/// Routine for doing a density grid read from a 1-D model.
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
          if (globals::vout_model[mgi] < vcell)
          {
            set_cell_modelgridindex(cellindex, mgi + 1);
          }
        }
      }

      if (globals::vout_model[get_cell_modelgridindex(cellindex)] >= vmin)
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
// Routine for doing a density grid read from a 2-D model.
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

      cell[n].modelgridindex = 0;
      const double zcylindrical = dcen[2];
      dcen[2] = 0.0;
      const double rcylindrical = vec_len(dcen);

      // Grid is uniform so only need to search in 1d to get r and z positions

      int mkeep1 = 0;
      for (int m = 0; m < globals::ncoord1_model; m++)
      {
        if (rcylindrical > (m * globals::dcoord1 * globals::tmin/globals::t_model))
        {
          mkeep1 = m;
          //cell[n].modelgridindex = m+1;
        }
      }

      int mkeep2 = 0;
      for (int m = 0; m < globals::ncoord2_model; m++)
      {
        if (zcylindrical > (((m * globals::dcoord2) * globals::tmin/globals::t_model) - globals::rmax))
        {
          mkeep2 = m;
          //cell[n].modelgridindex = m+1;
        }
      }
      cell[n].modelgridindex = (mkeep2 * globals::ncoord1_model) + mkeep1;
      globals::modelgrid[cell[n].modelgridindex].initial_radial_pos += radial_pos;

      //renorm[mkeep]++;
    }
    else
    {
      cell[n].modelgridindex = MMODELGRID;
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

      for (int element = 0; element < get_nelements(); element++)
      {
        ///now set the abundances (by mass) of included elements, i.e.
        ///read out the abundances specified in the atomic data file
        const int anumber = get_element(element);
        const float elemabundance = abundances_in[anumber - 1] / normfactor;
        assert(elemabundance >= 0.);
        globals::modelgrid[mgi].composition[element].abundance = elemabundance;

        // radioactive nuclide abundances should have already been set by read_??_model
        set_elem_stable_abund_from_total(mgi, anumber, elemabundance);
      }
    }
  }

  fclose(abundance_file);
}


static void read_2d3d_modelabundanceline(FILE * model_input, const int mgi, const bool keep)
{
  char line[1024] = "";
  if (line != fgets(line, 1024, model_input))
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
      set_modelinitradioabund(mgi, NUCLIDE_NI56, f56ni_model);
      set_modelinitradioabund(mgi, NUCLIDE_CO56, f56co_model);
      set_modelinitradioabund(mgi, NUCLIDE_NI57, f57ni_model);
      set_modelinitradioabund(mgi, NUCLIDE_CO57, f57co_model);
      set_modelinitradioabund(mgi, NUCLIDE_FE52, f52fe_model);
      set_modelinitradioabund(mgi, NUCLIDE_CR48, f48cr_model);
      set_modelinitradioabund(mgi, NUCLIDE_V48, 0.);

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
/// Subroutine to read in a 1-D model.
{
  FILE *model_input = fopen_required("model.txt", "r");

  // 1st read the number of data points in the table of input model.
  int npts_model_in = 0;
  fscanf(model_input, "%d", &npts_model_in);
  set_npts_model(npts_model_in);

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

  int mgi = 0;
  while (!feof(model_input))
  {
    char line[1024] = "";
    if (line != fgets(line, 1024, model_input))
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
      assert(cellnumberin == mgi + 1);

      globals::vout_model[mgi] = vout_kmps * 1.e5;

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
    set_modelinitradioabund(mgi, NUCLIDE_NI56, f56ni_model);
    set_modelinitradioabund(mgi, NUCLIDE_CO56, f56co_model);
    set_modelinitradioabund(mgi, NUCLIDE_NI57, f57ni_model);
    set_modelinitradioabund(mgi, NUCLIDE_CO57, f57co_model);
    set_modelinitradioabund(mgi, NUCLIDE_FE52, f52fe_model);
    set_modelinitradioabund(mgi, NUCLIDE_CR48, f48cr_model);
    set_modelinitradioabund(mgi, NUCLIDE_V48, 0.);
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

  globals::vmax = globals::vout_model[get_npts_model() - 1];
}


static void read_2d_model(void)
/// Subroutine to read in a 2-D model.
{
  FILE *model_input = fopen_required("model.txt", "r");

  // 1st read the number of data points in the table of input model.
  fscanf(model_input, "%d %d", &globals::ncoord1_model, &globals::ncoord2_model);  // r and z (cylindrical polar)

  set_npts_model(globals::ncoord1_model * globals::ncoord2_model);

  // Now read the time (in days) at which the model is specified.
  double t_model_days;
  fscanf(model_input, "%lg", &t_model_days);
  globals::t_model = t_model_days * DAY;

  // Now read in globals::vmax (in cm/s)
  fscanf(model_input, "%lg\n", &globals::vmax);
  globals::dcoord1 = globals::vmax * globals::t_model / globals::ncoord1_model; //dr for input model
  globals::dcoord2 = 2. * globals::vmax * globals::t_model / globals::ncoord2_model; //dz for input model

  // Now read in the model. Each point in the model has two lines of input.
  // First is an index for the cell then its r-mid point then its z-mid point
  // then its total mass density.
  // Second is the total FeG mass, initial 56Ni mass, initial 56Co mass

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
    assert(items_read == 4);

    const int ncoord1 = ((cellnumin - 1) % globals::ncoord1_model);
    const double r_cylindrical = (ncoord1 + 0.5) * globals::dcoord1;
    assert(fabs(cell_r_in / r_cylindrical - 1) < 1e-3);
    const int ncoord2 = ((cellnumin - 1) / globals::ncoord1_model);
    const double z = -globals::vmax * globals::t_model + ((ncoord2 + 0.5) * globals::dcoord2);
    assert(fabs(cell_z_in / z - 1) < 1e-3);

    assert(cellnumin == mgi + 1);

    const double rho_tmin = rho_tmodel * pow(globals::t_model / globals::tmin, 3);
    set_rhoinit(mgi, rho_tmin);
    set_rho(mgi, rho_tmin);

    read_2d3d_modelabundanceline(model_input, mgi, true);

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
  if (npts_model_in != globals::ngrid)
  {
    printout("3D model/grid mismatch. Abort. %d != %d\n", npts_model_in, globals::ngrid);
    abort();
  }

  /// Now read the time (in days) at which the model is specified.
  float t_model_days;
  fscanf(model_input, "%g", &t_model_days);
  globals::t_model = t_model_days * DAY;

  /// Now read in globals::vmax for the model (in cm s^-1).
  fscanf(model_input, "%lg\n", &globals::vmax);

  // double rmax_tmodel = globals::vmax * t_model;

  /// Now read in the lines of the model.
  globals::min_den = 1.e99;

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
    int items_read = sscanf(line, "%d %g %g %g %g", &mgi_in, &cellpos_in[2], &cellpos_in[1], &cellpos_in[0], &rho_model);
    assert(items_read == 5);
    //printout("cell %d, posz %g, posy %g, posx %g, rho %g, rho_init %g\n",dum1,dum3,dum4,dum5,rho_model,rho_model* pow( (t_model/globals::tmin), 3.));

    assert(mgi_in == n + 1);

    // cell coordinates in the 3D model.txt file are sometimes reordered by the scaling script
    // however, the cellindex always should increment X first, then Y, then Z

    // for (int axis = 0; axis < 3; axis++)
    // {
    //   const double cellpos_expected = - rmax_tmodel + (2 * get_cellcoordpointnum(n, axis) * rmax_tmodel / globals::ncoordgrid[axis]);
    //   // printout("n %d axis %d expected %g found %g rmax %g get_cellcoordpointnum(n, axis) %d globals::ncoordgrid %d\n",
    //   // n, axis, cellpos_expected, cellpos_in[axis], rmax_model, get_cellcoordpointnum(n, axis), globals::ncoordgrid);
    //   assert((fabs(cellpos_in[axis] / cellpos_expected - 1) < 1e-3) || ((cellpos_in[axis] == 0) && (cellpos_expected == 0)));
    // }

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

      if (globals::min_den > rho_model)
      {
        globals::min_den = rho_model;
      }
    }
    else
    {
      set_cell_modelgridindex(n, MMODELGRID);
    }

    read_2d3d_modelabundanceline(model_input, mgi, keepcell);
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

  printout("min_den %g\n", globals::min_den);
  printout("Effectively used model grid cells %d\n", mgi);

  /// Now, set actual size of the modelgrid to the number of non-empty cells.

  set_npts_model(mgi);

  fclose(model_input);
}


void read_ejecta_model(enum model_types model_type)
{
  switch (model_type)
  {
    case RHO_UNIFORM:
    {
      assert(false); // needs to be reimplemented using spherical coordinate mode
      break;
    }

    case RHO_1D_READ:
    {
      printout("Read 1D model!\n");
      read_1d_model();
      break;
    }

    case RHO_2D_READ:
    {
      printout("Read 2D model!\n");

      read_2d_model();
      break;
    }

    case RHO_3D_READ:
    {
      printout("Read 3D model!\n");
      read_3d_model();
      break;
    }

    default:
    {
      printout("Unknown model. Abort.\n");
      abort();
    }
  }
}


static void read_grid_restart_data(const int timestep)
{
  char filename[100];
  sprintf(filename, "gridsave_ts%d.tmp", timestep);

  printout("READIN GRID SNAPSHOT from %s\n", filename);
  FILE *gridsave_file = fopen_required(filename, "r");

  int ntstep_in;
  fscanf(gridsave_file, "%d ", &ntstep_in);
  assert(ntstep_in = globals::ntstep);

  for (int nts = 0; nts < globals::ntstep; nts++)
  {
    fscanf(gridsave_file, "%lg %lg ", &globals::time_step[nts].gamma_dep, &globals::time_step[nts].positron_dep);
  }

  int timestep_in;
  fscanf(gridsave_file, "%d ", &timestep_in);
  assert(timestep_in = timestep);

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

    assert(T_R >= 0);
    assert(T_e >= 0);
    assert(W >= 0);
    assert(T_J >= 0);
    assert(globals::rpkt_emiss[mgi] >= 0);

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

  const double e_ni56 = decay::nucdecayenergy(NUCLIDE_NI56);
  const double t_ni56 = decay::meanlife(NUCLIDE_NI56);
  const double e_co56 = decay::nucdecayenergy(NUCLIDE_CO56);
  const double t_co56 = decay::meanlife(NUCLIDE_CO56);

  const double factor56ni = 1. / 56 / MH * (-1. / (tstart * (- t_co56 + t_ni56)))
    * (- e_ni56 * exp(- tstart / t_ni56) * tstart * t_co56 - e_ni56 * exp(- tstart / t_ni56) * t_ni56 * t_co56
       + e_ni56 * exp(- tstart / t_ni56) * tstart * t_ni56 + pow(t_ni56, 2) * e_ni56 * exp(- tstart / t_ni56)
       - t_co56 * tstart * e_co56 * exp(- tstart / t_co56) - pow(t_co56, 2) * e_co56 * exp(- tstart / t_co56)
       + e_co56 * tstart * t_ni56 * exp(- tstart / t_ni56) + pow(t_ni56, 2) * e_co56 * exp(- tstart / t_ni56)
       + e_ni56 * t_co56 * t_ni56 - e_ni56 * pow(t_ni56, 2) - pow(t_ni56, 2) * e_co56 + e_co56 * pow(t_co56, 2));

  const double factor56co = 1. / 56 / MH * (1. / (tstart * t_co56))
    * (t_co56 * tstart * e_co56 * exp(- tstart / t_co56) + pow(t_co56, 2) * e_co56 * exp(- tstart / t_co56));

  const double e_ni57 = decay::nucdecayenergy(NUCLIDE_NI57);
  const double t_ni57 = decay::meanlife(NUCLIDE_NI57);
  const double e_co57 = decay::nucdecayenergy(NUCLIDE_CO57);
  const double t_co57 = decay::meanlife(NUCLIDE_CO57);

  const double factor57ni = 1. / 57 / MH * (-1. / (tstart * (- t_co57 + t_ni57)))
    * (- e_ni57 * exp(- tstart / t_ni57) * tstart * t_co57 - e_ni57 * exp(- tstart / t_ni57) * t_ni57 * t_co57
       + e_ni57 * exp(- tstart / t_ni57) * tstart * t_ni57 + pow(t_ni57, 2) * e_ni57 * exp(- tstart / t_ni57)
       - t_co57 * tstart * e_co57 * exp(- tstart / t_co57) - pow(t_co57, 2) * e_co57 * exp(- tstart / t_co57)
       + e_co57 * tstart * t_ni57 * exp(- tstart / t_ni57) + pow(t_ni57, 2) * e_co57 * exp(- tstart / t_ni57)
       + e_ni57 * t_co57 * t_ni57 - e_ni57 * pow(t_ni57, 2) - pow(t_ni57, 2) * e_co57 + e_co57 * pow(t_co57, 2));

 const double e_fe52 = decay::nucdecayenergy(NUCLIDE_FE52);
 const double t_fe52 = decay::meanlife(NUCLIDE_FE52);
 const double e_mn52 = decay::nucdecayenergy(NUCLIDE_MN52);
 const double t_mn52 = decay::meanlife(NUCLIDE_MN52);

  const double factor52fe = 1. / 52 / MH * (-1. / (tstart * (- t_mn52 + t_fe52)))
    * (- e_fe52 * exp(- tstart / t_fe52) * tstart * t_mn52 - e_fe52 * exp(- tstart / t_fe52) * t_fe52 * t_mn52
       + e_fe52 * exp(- tstart / t_fe52) * tstart * t_fe52 + pow(t_fe52, 2) * e_fe52 * exp(- tstart / t_fe52)
       - t_mn52 * tstart * e_mn52 * exp(- tstart / t_mn52) - pow(t_mn52, 2) * e_mn52 * exp(- tstart / t_mn52)
       + e_mn52 * tstart * t_fe52 * exp(- tstart / t_fe52) + pow(t_fe52, 2) * e_mn52 * exp(- tstart / t_fe52)
       + e_fe52 * t_mn52 * t_fe52 - e_fe52 * pow(t_fe52, 2) - pow(t_fe52, 2) * e_mn52 + e_mn52 * pow(t_mn52, 2));

  const double e_cr48 = decay::nucdecayenergy(NUCLIDE_CR48);
  const double t_cr48 = decay::meanlife(NUCLIDE_CR48);
  const double e_v48 = decay::nucdecayenergy(NUCLIDE_V48);
  const double t_v48 = decay::meanlife(NUCLIDE_V48);

  const double factor48cr = 1. / 48 / MH * (-1. / (tstart * (- t_v48 + t_cr48)))
    * (- e_cr48 * exp(- tstart / t_cr48) * tstart * t_v48 - e_cr48 * exp(- tstart / t_cr48) * t_cr48 * t_v48
       + e_cr48 * exp(- tstart / t_cr48) * tstart * t_cr48 + pow(t_cr48, 2) * e_cr48 * exp(- tstart / t_cr48)
       - t_v48 * tstart * e_v48 * exp(- tstart / t_v48) - pow(t_v48, 2) * e_v48 * exp(- tstart / t_v48)
       + e_v48 * tstart * t_cr48 * exp(- tstart / t_cr48) + pow(t_cr48, 2) * e_v48 * exp(- tstart / t_cr48)
       + e_cr48 * t_v48 * t_cr48 - e_cr48 * pow(t_cr48, 2) - pow(t_cr48, 2) * e_v48 + e_v48 * pow(t_v48, 2));

  //factor56ni = CLIGHT/4/STEBO * nucdecayenergy(NUCLIDE_NI56)/56/MH;
  /// This works only for the inbuilt Lucy model
  //factor56ni = CLIGHT/4/STEBO * 3*mtot/4/PI * nucdecayenergy(NUCLIDE_NI56)/56/MH  / pow(vmax,3);
  for (int mgi = 0; mgi < get_npts_model(); mgi++)
  {
    // alternative, not correct because it neglects expansion factor
    // const double decayedenergy_per_mass = get_decayedenergy_per_ejectamass(mgi, tstart);
    //
    // double T_initial = pow(CLIGHT / 4 / STEBO  * pow(globals::tmin / tstart, 3) * get_rhoinit(mgi) * decayedenergy_per_mass, 1. / 4.);

    double T_initial = pow(CLIGHT / 4 / STEBO  * pow(globals::tmin / tstart, 3) * get_rhoinit(mgi) * (
         (factor56ni * get_modelinitradioabund(mgi, NUCLIDE_NI56)) +
         (factor56co * get_modelinitradioabund(mgi, NUCLIDE_CO56)) +
         (factor57ni * get_modelinitradioabund(mgi, NUCLIDE_NI57)) +
         // (factor57co * get_modelinitradioabund(mgi, NUCLIDE_CO57)) +
         (factor52fe * get_modelinitradioabund(mgi, NUCLIDE_FE52)) +
         (factor48cr * get_modelinitradioabund(mgi, NUCLIDE_CR48))), 1. / 4.);

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
  globals::coordlabel[0] = 'X';
  globals::coordlabel[1] = 'Y';
  globals::coordlabel[2] = 'Z';
  int nxyz[3] = {0, 0, 0};
  assert(globals::ngrid == globals::ncoordgrid[0] * globals::ncoordgrid[1] * globals::ncoordgrid[2]);
  for (int n = 0; n < globals::ngrid; n++)
  {
    for (int axis = 0; axis < 3; axis++)
    {
      assert(nxyz[axis] == get_cellcoordpointnum(n, axis));
      cell[n].pos_init[axis] = - globals::coordmax[axis] + (2 * nxyz[axis] * globals::coordmax[axis] / globals::ncoordgrid[axis]);
      // cell[n].xyz[axis] = nxyz[axis];
    }

    assert(n == nxyz[2] * globals::ncoordgrid[1] * globals::ncoordgrid[2] + nxyz[1] * globals::ncoordgrid[0] + nxyz[0]);

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
  assert(get_model_type() == RHO_1D_READ);
  globals::coordlabel[0] = 'r';
  globals::coordlabel[1] = '?';
  globals::coordlabel[2] = '?';

  globals::ncoordgrid[0] = get_npts_model();
  globals::ncoordgrid[1] = 1;
  globals::ncoordgrid[2] = 1;
  globals::ngrid = globals::ncoordgrid[0] * globals::ncoordgrid[1] * globals::ncoordgrid[2];
  globals::coordmax[0] = globals::rmax;
  globals::coordmax[1] = 0.;
  globals::coordmax[2] = 0.;

  assert(globals::ngrid <= MGRID);

  // in this mode, cellindex and modelgridindex are the same thing
  for (int cellindex = 0; cellindex < get_npts_model(); cellindex++)
  {
    const double v_inner = cellindex > 0 ? globals::vout_model[cellindex - 1] : 0.;
    cell[cellindex].modelgridindex = cellindex;
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

  /// The cells will be ordered by x then y, then z. Call a routine that
  /// sets up the initial positions and widths of the cells.
  if (globals::grid_type == GRID_UNIFORM)
  {
    uniform_grid_setup();
  }
  else if (globals::grid_type == GRID_SPHERICAL1D)
  {
    spherical1d_grid_setup();
  }
  else
  {
    printout("[fatal] grid_init: Error: Unknown grid type. Abort.");
    abort();
  }

  /// Now set up the density in each cell.

  // Calculate the critical opacity at which opacity_case 3 switches from a
  // regime proportional to the density to a regime independent of the density
  // This is done by solving for tau_sobolev == 1
  // tau_sobolev = PI*QE*QE/(ME*C) * rho_crit_para * rho/nucmass(NUCLIDE_NI56) * 3000e-8 * globals::time_step[m].mid;
  globals::rho_crit = ME * CLIGHT * decay::nucmass(NUCLIDE_NI56) / (PI * QE * QE * globals::rho_crit_para * 3000e-8 * globals::tmin);
  printout("grid_init: rho_crit = %g\n", globals::rho_crit);

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
    for (int n = 0; n < globals::ngrid; n++)
    {
      const double radial_pos = get_cellradialpos(n);
      globals::modelgrid[cell[n].modelgridindex].initial_radial_pos = radial_pos;
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
  for (int iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
  {
    if (totmassradionuclide[iso] <= 0)
      continue;
    double totmassradionuclide_actual = 0.;
    for (int mgi = 0; mgi < get_npts_model(); mgi++)
    {
      totmassradionuclide_actual += get_modelinitradioabund(mgi, (enum radionuclides) iso) * get_rhoinit(mgi) * vol_init_modelcell(mgi);
    }
    if (totmassradionuclide_actual >= 0.)
    {
      const double ratio = totmassradionuclide[iso] / totmassradionuclide_actual;
      // printout("nuclide %d ratio %g\n", iso, ratio);
      for (int mgi = 0; mgi < get_npts_model(); mgi++)
      {
        if (get_numassociatedcells(mgi) > 0)
        {
          const double prev_abund = get_modelinitradioabund(mgi, (enum radionuclides)(iso));
          const double new_abund = prev_abund * ratio;
          set_modelinitradioabund(mgi, (enum radionuclides)(iso), new_abund);
        }
      }
    }
  }
}


void show_totmassradionuclides(void)
{
  globals::mtot = 0.;
  globals::mfeg = 0.;

  for (int iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
    totmassradionuclide[iso] = 0.;

  int n1 = 0;
  for (int mgi = 0; mgi < get_npts_model(); mgi++)
  {
    double cellvolume = 0.;
    if (get_model_type() == RHO_1D_READ)
    {
      const double v_inner = (mgi == 0) ? 0. : globals::vout_model[mgi - 1];
      // mass_in_shell = rho_model[mgi] * (pow(globals::vout_model[mgi], 3) - pow(v_inner, 3)) * 4 * PI * pow(t_model, 3) / 3.;
      cellvolume = (pow(globals::vout_model[mgi], 3) - pow(v_inner, 3)) * 4 * PI * pow(globals::tmin, 3) / 3.;
    }
    else if (get_model_type() == RHO_2D_READ)
    {
      cellvolume = pow(globals::tmin / globals::t_model, 3) * ((2 * n1) + 1) * PI * globals::dcoord2 * pow(globals::dcoord1, 2.);
      n1++;
      if (n1 == globals::ncoord1_model)
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

    globals::mtot += mass_in_shell;

    for (int isoint = 0; isoint < RADIONUCLIDE_COUNT; isoint++)
    {
      const enum radionuclides iso = (enum radionuclides) isoint;
      totmassradionuclide[iso] += mass_in_shell * get_modelinitradioabund(mgi, iso);
    }

    globals::mfeg += mass_in_shell * get_ffegrp(mgi);
  }


  printout("Masses / Msun:    Total: %9.3e  56Ni: %9.3e  56Co: %9.3e  52Fe: %9.3e  48Cr: %9.3e\n",
           globals::mtot / MSUN, totmassradionuclide[NUCLIDE_NI56] / MSUN,
           totmassradionuclide[NUCLIDE_CO56] / MSUN, totmassradionuclide[NUCLIDE_FE52] / MSUN,
           totmassradionuclide[NUCLIDE_CR48] / MSUN);
  printout("Masses / Msun: Fe-group: %9.3e  57Ni: %9.3e  57Co: %9.3e\n",
           globals::mfeg / MSUN, totmassradionuclide[NUCLIDE_NI57] / MSUN, totmassradionuclide[NUCLIDE_CO57] / MSUN);
}


double get_totmassradionuclide(enum radionuclides nuc)
{
  return totmassradionuclide[nuc];
}