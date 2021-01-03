#include "sn3d.h"
#include "atomic.h"
#include "grid.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "decay.h"
#include "radfield.h"
#include "rpkt.h"
#include "vectors.h"


enum __managed__ model_types model_type = RHO_1D_READ;

static long mem_usage_nltepops = 0;

static __managed__ int mg_associated_cells[MMODELGRID + 1];


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
      const int modelgridindex = globals::cell[cellindex].modelgridindex;
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
      const int mgi = globals::cell[cellindex].modelgridindex;
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
  return globals::cell[cellindex].pos_init[axis];
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
      }
  }
}


__host__ __device__
int get_cellcoordpointnum(const int cellindex, const int axis)
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


__host__ __device__
static void calculate_kappagrey(void)
{
  double rho_sum = 0.0;
  double fe_sum = 0.0;
  double opcase3_sum = 0.0;
  int empty_cells = 0;

  for (int n = 0; n < globals::ngrid; n++)
  {
    const int mgi = globals::cell[n].modelgridindex;
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
    const int mgi = globals::cell[n].modelgridindex;
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


__host__ __device__
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


__host__ __device__
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
    const int mgi = globals::cell[cellindex].modelgridindex;
    assert(!(get_model_type() == RHO_3D_READ) || (get_rhoinit(mgi) > 0) || (mgi == MMODELGRID));
    mg_associated_cells[mgi] += 1;
    assert(!(get_model_type() == RHO_3D_READ) || (mg_associated_cells[mgi] == 1) || (mgi == MMODELGRID));
  }

  int numnonemptycells = 0;
  for (int mgi = 0; mgi < globals::npts_model; mgi++)
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


__host__ __device__
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
        globals::cell[cellindex].modelgridindex = cellindex;
      }
      else
      {
        globals::cell[cellindex].modelgridindex = 0;

        for (int mgi = 0; mgi < (globals::npts_model - 1); mgi++)
        {
          if (globals::vout_model[mgi] < vcell)
          {
            globals::cell[cellindex].modelgridindex = mgi + 1;
          }
        }
      }

      if (globals::vout_model[globals::cell[cellindex].modelgridindex] >= vmin)
      {
        globals::modelgrid[globals::cell[cellindex].modelgridindex].initial_radial_pos += radial_pos;
      }
      else
      {
        globals::cell[cellindex].modelgridindex = MMODELGRID;
      }
    }
    else
    {
      globals::cell[cellindex].modelgridindex = MMODELGRID;
    }
  }
}


__host__ __device__
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

      globals::cell[n].modelgridindex = 0;
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
          //globals::cell[n].modelgridindex = m+1;
        }
      }

      int mkeep2 = 0;
      for (int m = 0; m < globals::ncoord2_model; m++)
      {
        if (zcylindrical > (((m * globals::dcoord2) * globals::tmin/globals::t_model) - globals::rmax))
        {
          mkeep2 = m;
          //globals::cell[n].modelgridindex = m+1;
        }
      }
      globals::cell[n].modelgridindex = (mkeep2 * globals::ncoord1_model) + mkeep1;
      globals::modelgrid[globals::cell[n].modelgridindex].initial_radial_pos += radial_pos;

      //renorm[mkeep]++;
    }
    else
    {
      globals::cell[n].modelgridindex = MMODELGRID;
    }
  }
}


__host__ __device__
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
  const int ncount = threedimensional ? globals::ngrid : globals::npts_model;
  for (int n = 0; n < ncount; n++)
  {
    const int mgi = threedimensional ? globals::cell[n].modelgridindex : n;

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


__host__ __device__
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

  for (int mgi = 0; mgi < globals::npts_model; mgi++)
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
            const int estimindex = mgi * get_nelements() * maxion + element * maxion + ion;
            fscanf(gridsave_file, " %lg %lg", &corrphotoionrenorm[estimindex], &gammaestimator[estimindex]);
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


__host__ __device__
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
  for (int mgi = 0; mgi < globals::npts_model; mgi++)
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


__host__ __device__
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
      globals::cell[n].pos_init[axis] = - globals::coordmax[axis] + (2 * nxyz[axis] * globals::coordmax[axis] / globals::ncoordgrid[axis]);
      // globals::cell[n].xyz[axis] = nxyz[axis];
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

  /*
  /// Finally we must also create the composition dependent data structure for
  /// the samplingcell which is located at MGRID (i.e. the MGRID+1th cell)
  n = MGRID;
  if ((globals::cell[n].composition = malloc(get_nelements()*sizeof(compositionlist_entry))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize compositionlist for cell %d... abort\n",n);
    abort();
  }
  for (element = 0; element < get_nelements(); element++)
  {
    ///now set the abundances (by mass) of included elements, i.e.
    ///read out the abundances specified in the atomic data file
    ///and allocate memory to store the ground level populations for each ionisation stage
    if ((globals::cell[n].composition[element].groundlevelpop = malloc(get_nions(element)*sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize groundlevelpoplist for element %d in cell %d... abort\n",element,n);
      abort();
    }
    if ((globals::cell[n].composition[element].partfunct = malloc(get_nions(element)*sizeof(float))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize partfunctlist for element %d in cell %d... abort\n",element,n);
      abort();
    }
//     if ((globals::cell[n].composition[element].ltepartfunct = malloc(get_nions(element)*sizeof(float))) == NULL)
//     {
//       printout("[fatal] input: not enough memory to initialize lte partfunctlist for element %d in cell %d... abort\n",element,n);
//       abort();
//     }
  }
  */
}

__host__ __device__
static void spherical1d_grid_setup(void)
{
  assert(get_model_type() == RHO_1D_READ);
  globals::coordlabel[0] = 'r';
  globals::coordlabel[1] = '?';
  globals::coordlabel[2] = '?';

  globals::ncoordgrid[0] = globals::npts_model;
  globals::ncoordgrid[1] = 1;
  globals::ncoordgrid[2] = 1;
  globals::ngrid = globals::ncoordgrid[0] * globals::ncoordgrid[1] * globals::ncoordgrid[2];
  globals::coordmax[0] = globals::rmax;
  globals::coordmax[1] = 0.;
  globals::coordmax[2] = 0.;

  assert(globals::ngrid <= MGRID);

  // in this mode, cellindex and modelgridindex are the same thing
  for (int cellindex = 0; cellindex < globals::npts_model; cellindex++)
  {
    const double v_inner = cellindex > 0 ? globals::vout_model[cellindex - 1] : 0.;
    globals::cell[cellindex].modelgridindex = cellindex;
    globals::cell[cellindex].pos_init[0] = v_inner * globals::tmin;
    globals::cell[cellindex].pos_init[1] = 0.;
    globals::cell[cellindex].pos_init[2] = 0.;
  }
}


__host__ __device__
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
  if (get_model_type() == RHO_UNIFORM)
  {
    //uniform_density_setup ();
    //abundances_setup();
  }
  else
  {
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
        globals::modelgrid[globals::cell[n].modelgridindex].initial_radial_pos = radial_pos;
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
    if (globals::totmassradionuclide[iso] <= 0)
      continue;
    double totmassradionuclide_actual = 0.;
    for (int mgi = 0; mgi < globals::npts_model; mgi++)
    {
      totmassradionuclide_actual += get_modelinitradioabund(mgi, (enum radionuclides) iso) * get_rhoinit(mgi) * vol_init_modelcell(mgi);
    }
    if (totmassradionuclide_actual >= 0.)
    {
      const double ratio = globals::totmassradionuclide[iso] / totmassradionuclide_actual;
      // printout("nuclide %d ratio %g\n", iso, ratio);
      for (int mgi = 0; mgi < globals::npts_model; mgi++)
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
