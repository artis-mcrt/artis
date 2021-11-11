#include "sn3d.h"
#include "atomic.h"
#include "grid.h"
#include "input.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "decay.h"
#include "radfield.h"
#include "rpkt.h"
#include "vectors.h"
#include <cstring>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace grid
{

__managed__ modelgrid_t *modelgrid = NULL;

__managed__ int ncoordgrid[3]; /// propagration grid dimensions
__managed__ int ngrid;
__managed__ int grid_type;
__managed__ char coordlabel[3];

__managed__ enum model_types model_type = RHO_1D_READ;
__managed__ int npts_model = 0; // number of model grid cells
__managed__ int nonempty_npts_model = 0; // number of allocated non-empty model grid cells

__managed__ double t_model = -1.; // time at which densities in input model are correct.
__managed__ double *vout_model = NULL;
__managed__ int ncoord_model[3]; // the model.txt input grid dimensions
__managed__ double dcoord1;
__managed__ double dcoord2; // spacings of a 2D model grid - must be uniform grid

__managed__ double min_den; // minimum model density

__managed__ double mtot;
__managed__ double mfeg;              /// Total mass of Fe group elements in ejecta

__managed__ CELL *cell = NULL;

static long mem_usage_nltepops = 0;

static __managed__ int *mg_associated_cells = NULL;
static __managed__ int *nonemptymgi_of_mgi = NULL;

__managed__ double *totmassradionuclide = NULL; /// total mass of each radionuclide in the ejecta


__host__ __device__
double wid_init(const int cellindex)
// for a uniform grid this is the extent along the x,y,z coordinate (x_2 - x_1, etc.)
// for spherical grid this is the radial extent (r_outer - r_inner)
// these values are for time globals::tmin
{
  switch (grid_type)
  {
    case GRID_SPHERICAL1D:
    {
      const int modelgridindex = get_cell_modelgridindex(cellindex);
      const double v_inner = modelgridindex > 0 ? vout_model[modelgridindex - 1] : 0.;
      return (vout_model[modelgridindex] - v_inner) * globals::tmin;
    }

    default:
      return 2 * globals::coordmax[0] / ncoordgrid[0];
  }
}


__host__ __device__
double vol_init_modelcell(const int modelgridindex)
// return the model cell volume at globals::tmin
// for a uniform cubic grid this is constant
{
  switch (grid_type)
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
  switch (grid_type)
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
  // return - coordmax[axis] + (2 * get_cellcoordpointnum(cellindex, axis) * coordmax[axis] / ncoordgrid[axis]);
}


__host__ __device__
int get_coordcellindexincrement(const int axis)
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
          return 1;

        case 1:
          return ncoordgrid[0];

        case 2:
          return ncoordgrid[0] * ncoordgrid[1];

        default:
          printout("invalid coordinate index %d", axis);
          abort();
          return -1;
      }
  }
}


__host__ __device__
int get_cellcoordpointnum(const int cellindex, const int axis)
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
          return -1;
      }
  }
}



__host__ __device__
int get_ngriddimensions(void)
{
  return (grid_type == GRID_SPHERICAL1D) ? 1 : 3;
}


__host__ __device__
float get_rhoinit(int modelgridindex)
{
  return modelgrid[modelgridindex].rhoinit;
}


__host__ __device__
float get_rho(int modelgridindex)
{
  return modelgrid[modelgridindex].rho;
}


__host__ __device__
float get_nne(int modelgridindex)
{
  assert_testmodeonly(modelgridindex >= 0);
  assert_testmodeonly(modelgridindex < (get_npts_model() + 1));

  const double nne = modelgrid[modelgridindex].nne;
  assert_testmodeonly(std::isfinite(nne));
  return nne;
}


__host__ __device__
float get_nnetot(int modelgridindex)
{
  assert_testmodeonly(modelgridindex >= 0);
  assert_testmodeonly(modelgridindex < (get_npts_model() + 1));

  const double nnetot = modelgrid[modelgridindex].nnetot;
  assert_always(std::isfinite(nnetot));
  return nnetot;
}


__host__ __device__
// the abundances referred to below are initial abundances
float get_ffegrp(int modelgridindex)
{
  return modelgrid[modelgridindex].ffegrp;
}


__host__ __device__
float get_elem_abundance(int modelgridindex, int element)
// mass fraction of an element (all isotopes combined)
{
  return modelgrid[modelgridindex].composition[element].abundance;
}


__host__ __device__
void set_elem_abundance(int modelgridindex, int element, float newabundance)
// mass fraction of an element (all isotopes combined)
{
  modelgrid[modelgridindex].composition[element].abundance = newabundance;
}


__host__ __device__
float get_kappagrey(int modelgridindex)
{
  return modelgrid[modelgridindex].kappagrey;
}


__host__ __device__
float get_Te(int modelgridindex)
{
  return modelgrid[modelgridindex].Te;
}


__host__ __device__
float get_TR(int modelgridindex)
{
  return modelgrid[modelgridindex].TR;
}


__host__ __device__
float get_TJ(int modelgridindex)
{
  return modelgrid[modelgridindex].TJ;
}


__host__ __device__
float get_W(int modelgridindex)
{
  return modelgrid[modelgridindex].W;
}


__host__ __device__
static void set_rhoinit(int modelgridindex, float x)
{
  modelgrid[modelgridindex].rhoinit = x;
}


__host__ __device__
static void set_rho(int modelgridindex, float x)
{
  modelgrid[modelgridindex].rho = x;
}


__host__ __device__
void set_nne(int modelgridindex, float x)
{
  modelgrid[modelgridindex].nne = x;
}


__host__ __device__
void set_nnetot(int modelgridindex, float x)
{
  assert_always(x >= 0.);
  assert_always(std::isfinite(x));
  modelgrid[modelgridindex].nnetot = x;
}


__host__ __device__
static void set_ffegrp(int modelgridindex, float x)
{
  assert_always(x >= 0);
  assert_always(x <= 1.001);
  modelgrid[modelgridindex].ffegrp = x;
}


__host__ __device__
void set_kappagrey(int modelgridindex, float kappagrey)
{
  modelgrid[modelgridindex].kappagrey = kappagrey;
}


__host__ __device__
void set_Te(int modelgridindex, float Te)
{
  modelgrid[modelgridindex].Te = Te;
}


__host__ __device__
void set_TR(int modelgridindex, float TR)
{
  modelgrid[modelgridindex].TR = TR;
}


__host__ __device__
void set_TJ(int modelgridindex, float TJ)
{
  modelgrid[modelgridindex].TJ = TJ;
}


__host__ __device__
void set_W(int modelgridindex, float W)
{
  modelgrid[modelgridindex].W = W;
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
// number of model grid cells
{
  assert_always(npts_model > 0);
  return npts_model;
}


__host__ __device__
int get_nonempty_npts_model(void)
// number of model grid cells
{
  assert_always(nonempty_npts_model > 0);
  return nonempty_npts_model;
}


__host__ __device__
static void set_npts_model(int new_npts_model)
{
  npts_model = new_npts_model;

  assert_always(modelgrid == NULL);
  modelgrid = (modelgrid_t *) calloc(npts_model + 1, sizeof(modelgrid_t));
  assert_always(mg_associated_cells == NULL);
  mg_associated_cells = (int *) malloc((get_npts_model() + 1) * sizeof(int));
  assert_always(nonemptymgi_of_mgi == NULL);
  nonemptymgi_of_mgi = (int *) malloc((get_npts_model() + 1) * sizeof(int));
  assert_always(modelgrid != NULL);
}


__host__ __device__
int get_t_model(void)
// get time at which model input densities are defined
{
  assert_testmodeonly(t_model > 0.);
  return t_model;
}


__host__ __device__
int get_cell_modelgridindex(int cellindex)
{
  assert_testmodeonly(cellindex >= 0);
  assert_testmodeonly(cellindex < ngrid);
  const int mgi = cell[cellindex].modelgridindex;
  assert_testmodeonly(mgi >= 0);
  assert_testmodeonly(mgi < (get_npts_model() + 1));
  return mgi;
}


__host__ __device__
static void set_cell_modelgridindex(int cellindex, int new_modelgridindex)
{
  assert_testmodeonly(cellindex < ngrid);
  assert_testmodeonly(new_modelgridindex <= get_npts_model());
  cell[cellindex].modelgridindex = new_modelgridindex;
}


__host__ __device__
int get_numassociatedcells(const int modelgridindex)
// number of propagation cells associated with each modelgrid cell
{
  assert_testmodeonly(mg_associated_cells != NULL);
  return mg_associated_cells[modelgridindex];
}


__host__ __device__
int get_modelcell_nonemptymgi(int mgi)
{
  // get the index in the list of non-empty cells for a given model grid cell
  assert_testmodeonly(nonempty_npts_model > 0);
  assert_testmodeonly(mgi < get_npts_model());

  const int nonemptymgi = nonemptymgi_of_mgi[mgi];
  assert_testmodeonly(nonemptymgi >= 0 || get_numassociatedcells(mgi) == 0);

  return nonemptymgi;
}


__host__ __device__
float get_modelinitradioabund(const int modelgridindex, const int z, const int a)
{
  // this function replaces get_f56ni(mgi), get_fco56(mgi), etc.

  const int nucindex = decay::get_nuc_index(z, a);
  assert_always(decay::get_nuc_z(nucindex) >= 0); // check not FAKE_GAM_LINE_ID nuclide

  assert_always(modelgrid[modelgridindex].initradioabund != NULL);
  return modelgrid[modelgridindex].initradioabund[nucindex];
}


__host__ __device__
static void set_modelinitradioabund(const int modelgridindex, const int z, const int a, const float abund)
{
  assert_always(abund >= 0.);
  assert_always(abund <= 1.);

  const int nucindex = decay::get_nuc_index(z, a);
  assert_always(decay::get_nuc_z(nucindex) >= 0 || abund == 0.); // check not FAKE_GAM_LINE_ID nuclide

  assert_always(modelgrid[modelgridindex].initradioabund != NULL);
  modelgrid[modelgridindex].initradioabund[nucindex] = abund;
}


__host__ __device__
float get_stable_initabund(const int mgi, const int element)
{
  return modelgrid[mgi].initmassfracstable[element];
}


__host__ __device__
static void set_elem_stable_abund_from_total(const int mgi, const int element, const float elemabundance)
{
  // set the stable mass fraction of an element from the total element mass fraction
  // by subtracting the abundances of radioactive isotopes.
  // if the element Z=anumber has no specific stable abundance variable then the function does nothing

  const int atomic_number = get_element(element);

  double isofracsum = 0.; // mass fraction sum of radioactive isotopes
  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++)
  {
    const int a = decay::get_nuc_a(nucindex);
    if (decay::get_nuc_z(nucindex) == atomic_number)
    {
      // radioactive isotope of this element
      isofracsum += get_modelinitradioabund(mgi, atomic_number, a);
    }
  }

  double massfracstable = elemabundance - isofracsum;

  if (massfracstable < 0.)
  {
      printout("WARNING: cell %d Z=%d element abundance is less than the sum of its radioisotope abundances \n",
               mgi, atomic_number);
      printout("  massfrac(Z) %g massfrac_radioisotopes(Z) %g\n", elemabundance, isofracsum);
      assert_always(massfracstable >= -1e-3);  // result is allowed to be slightly negative due to roundoff error
      printout("  increasing elemental abundance to %g and setting stable isotopic abundance to zero\n", isofracsum);
      massfracstable = 0.; // bring up to zero if negative
  }

  modelgrid[mgi].initmassfracstable[element] = massfracstable;

  // (isofracsum + massfracstable) might not exactly match elemabundance if we had to boost it to reach isofracsum
  modelgrid[mgi].composition[element].abundance = isofracsum + massfracstable;
}


__host__ __device__
double get_cellradialpos(const int cellindex)
// get the radial distance from the origin to the centre of the cell
{
  // spherical coordinate case is trivial
  if (grid_type == GRID_SPHERICAL1D)
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
  return modelgrid[modelgridindex].elements_uppermost_ion[element];
}


__host__ __device__
void set_elements_uppermost_ion(const int modelgridindex, const int element, const int newvalue)
{
  modelgrid[modelgridindex].elements_uppermost_ion[element] = newvalue;
}


static void calculate_kappagrey(void)
{
  double rho_sum = 0.0;
  double fe_sum = 0.0;
  double opcase3_sum = 0.0;
  int empty_cells = 0;

  for (int n = 0; n < ngrid; n++)
  {
    const int mgi = get_cell_modelgridindex(n);
    rho_sum +=get_rhoinit(mgi);
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
      opcase3_sum += get_kappagrey(mgi) *get_rhoinit(mgi);
    }
  }

  FILE *grid_file = NULL;
  if (globals::rank_global == 0)
  {
    grid_file = fopen_required("grid.out", "w");
  }

  /// Second pass through allows calculation of normalized kappa_grey
  double check1 = 0.0;
  double check2 = 0.0;
  for (int n = 0; n < ngrid; n++)
  {
    const int mgi = get_cell_modelgridindex(n);
    if (globals::rank_global == 0 && mgi != get_npts_model())
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
        const double opcase2_normal = GREY_OP * rho_sum / ((0.9 * fe_sum) + (0.1 * (ngrid - empty_cells)));
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
        set_kappagrey(mgi, ((0.9 * get_ffegrp(mgi)) + 0.1) * GREY_OP / ((0.9 * mfeg / mtot) + 0.1));
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
  modelgrid[modelgridindex].elements_uppermost_ion = (int *) malloc(get_nelements() * sizeof(int));

  assert_always(modelgrid[modelgridindex].elements_uppermost_ion != NULL);

  modelgrid[modelgridindex].composition = (compositionlist_entry *) malloc(get_nelements() * sizeof(compositionlist_entry));

  if (modelgrid[modelgridindex].composition == NULL)
  {
    printout("[fatal] input: not enough memory to initialize compositionlist for cell %d... abort\n",modelgridindex);
    abort();
  }

  modelgrid[modelgridindex].initmassfracstable = (float *) malloc(get_nelements() * sizeof(float));
  assert_always(modelgrid[modelgridindex].initmassfracstable != NULL);

  mem_usage_nltepops += globals::total_nlte_levels * sizeof(double);

#ifdef MPI_ON
  MPI_Win win;
  MPI_Aint size = (globals::rank_in_node == 0) ? globals::total_nlte_levels * sizeof(double) : 0;
  MPI_Win_allocate_shared(size, sizeof(double), MPI_INFO_NULL, globals::mpi_comm_node,
                          &modelgrid[modelgridindex].nlte_pops, &win);
  if (globals::rank_in_node != 0)
  {
    int disp_unit;
    MPI_Win_shared_query(win, MPI_PROC_NULL, &size, &disp_unit, &modelgrid[modelgridindex].nlte_pops);
  }
#else
  modelgrid[modelgridindex].nlte_pops = (double *) malloc(globals::total_nlte_levels * sizeof(double));
#endif

  if (modelgrid[modelgridindex].nlte_pops == NULL)
  {
    printout("[fatal] input: not enough memory to initialize nlte memory for cell %d... abort\n",modelgridindex);
    abort();
  }

  for (int nlteindex = 0; nlteindex < globals::total_nlte_levels; nlteindex++)
  {
    modelgrid[modelgridindex].nlte_pops[nlteindex] = -1.0; ///flag to indicate that there is
                                                           /// currently no information on the nlte populations
  }

  //printout("Managed to allocate memory for %d nlte levels\n", total_nlte_levels);

  for (int element = 0; element < get_nelements(); element++)
  {
    /// Set initial abundances to zero
    modelgrid[modelgridindex].composition[element].abundance = 0.;

    /// and allocate memory to store the ground level populations for each ionisation stage
    modelgrid[modelgridindex].composition[element].groundlevelpop = (float *) calloc(get_nions(element), sizeof(float));
    if (modelgrid[modelgridindex].composition[element].groundlevelpop == NULL)
    {
      printout("[fatal] input: not enough memory to initialize groundlevelpoplist for element %d in cell %d... abort\n",element,modelgridindex);
      abort();
    }

    modelgrid[modelgridindex].composition[element].partfunct = (float *) calloc(get_nions(element), sizeof(float));

    if (modelgrid[modelgridindex].composition[element].partfunct == NULL)
    {
      printout("[fatal] input: not enough memory to initialize partfunctlist for element %d in cell %d... abort\n",element,modelgridindex);
      abort();
    }
  }
}


static void allocate_cooling(const int modelgridindex)
/// Initialise composition dependent cell data for the given cell
{
  modelgrid[modelgridindex].cooling = (mgicooling_t *) malloc(get_nelements() * sizeof(mgicooling_t));

  if (modelgrid[modelgridindex].cooling == NULL)
  {
    printout("[fatal] input: not enough memory to initialize coolinglist for cell %d... abort\n",modelgridindex);
    abort();
  }

  for (int element = 0; element < get_nelements(); element++)
  {
    /// and allocate memory to store the ground level populations for each ionisation stage

    modelgrid[modelgridindex].cooling[element].contrib = (double *) malloc(get_nions(element) * sizeof(double));
    if (modelgrid[modelgridindex].cooling[element].contrib == NULL)
    {
      printout("[fatal] input: not enough memory to initialize coolinglist for element %d in cell %d... abort\n",element,modelgridindex);
      abort();
    }
  }
}


static void allocate_nonemptymodelcells(void)
{
  mem_usage_nltepops = 0;
  /// This is the placeholder for empty cells. Temperatures must be positive
  /// as long as ff opacities are calculated.
  set_rhoinit(get_npts_model(), 0.);
  set_rho(get_npts_model(), 0.);
  set_nne(get_npts_model(), 0.);
  set_nnetot(get_npts_model(), 0.);
  set_ffegrp(get_npts_model(), 0.);
  modelgrid[get_npts_model()].initradioabund = (float *) malloc(decay::get_num_nuclides() * sizeof(float));
  assert_always(modelgrid[get_npts_model()].initradioabund != NULL);

  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++)
  {
    const int z = decay::get_nuc_z(nucindex);
    const int a = decay::get_nuc_a(nucindex);
    set_modelinitradioabund(get_npts_model(), z, a, 0.);
  }
  set_Te(get_npts_model(), MINTEMP);
  set_TJ(get_npts_model(), MINTEMP);
  set_TR(get_npts_model(), MINTEMP);
  allocate_compositiondata(get_npts_model());
  allocate_cooling(get_npts_model());

  // Determine the number of simulation cells associated with the model cells
  for (int mgi = 0; mgi < (get_npts_model() + 1); mgi++)
  {
    mg_associated_cells[mgi] = 0;
  }

  for (int cellindex = 0; cellindex < ngrid; cellindex++)
  {
    const int mgi = get_cell_modelgridindex(cellindex);
    assert_always(!(get_model_type() == RHO_3D_READ) || (get_rhoinit(mgi) > 0) || (mgi == get_npts_model()));
    mg_associated_cells[mgi] += 1;
    assert_always(!(get_model_type() == RHO_3D_READ) || (mg_associated_cells[mgi] == 1) || (mgi == get_npts_model()));
  }

  int nonemptymgi = 0;  // index within list of non-empty modelgrid cells
  for (int mgi = 0; mgi < get_npts_model(); mgi++)
  {
    if (get_numassociatedcells(mgi) > 0)
    {
      allocate_compositiondata(mgi);
      allocate_cooling(mgi);

      if (get_rhoinit(mgi) <= 0)
      {
        printout("Error: negative or zero density. Abort.\n");
        abort();
      }
      nonemptymgi_of_mgi[mgi] = nonemptymgi;
      nonemptymgi++;
    }
    else
    {
      nonemptymgi_of_mgi[mgi] = -1;
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

  nonempty_npts_model = nonemptymgi;

  printout("[info] mem_usage: the modelgrid array occupies %.3f MB\n", (get_npts_model() + 1) * sizeof(modelgrid) / 1024. / 1024.);

  printout("There are %d modelgrid cells with associated propagation cells\n", nonempty_npts_model);

  printout("[info] mem_usage: NLTE populations for all allocated cells occupy a total of %.3f MB (shared node memory)\n", mem_usage_nltepops / 1024. / 1024.);
}


static void density_1d_read(void)
// Map 1D spherical model grid onto 3D propagation grid
{
  for (int cellindex = 0; cellindex < ngrid; cellindex++)
  {
    const double radial_pos = get_cellradialpos(cellindex);
    const double vcell = radial_pos / globals::tmin;
    const double vmin = 0.;
    if (radial_pos < globals::rmax)
    {
      if (grid_type == GRID_SPHERICAL1D)
      {
        set_cell_modelgridindex(cellindex, cellindex);
      }
      else
      {
        int mgi = 0;

        for (int i = 0; i < (get_npts_model() - 1); i++)
        {
          if (vout_model[mgi] < vcell)
          {
            mgi = i + 1;
          }
        }
        set_cell_modelgridindex(cellindex, mgi);
      }

      if (vout_model[get_cell_modelgridindex(cellindex)] >= vmin)
      {
        modelgrid[get_cell_modelgridindex(cellindex)].initial_radial_pos += radial_pos;
      }
      else
      {
        set_cell_modelgridindex(cellindex, get_npts_model());
      }
    }
    else
    {
      set_cell_modelgridindex(cellindex, get_npts_model());
    }
  }
}


static void density_2d_read(void)
// Map 2D model grid onto 3D propagation grid
{
  for (int n = 0; n < ngrid; n++)
  {
    const double radial_pos = get_cellradialpos(n);

    if (radial_pos < globals::rmax)
    {
      double dcen[3];
      for (int d = 0; d < 3; d++)
      {
        const double cellcoordmin = - globals::coordmax[d] + (2 * get_cellcoordpointnum(n, d) * globals::coordmax[d] / ncoordgrid[0]);
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
        if (rcylindrical > (m * dcoord1 * globals::tmin/t_model))
        {
          mkeep1 = m;
          // set_cell_modelgridindex(n, m + 1);
        }
      }

      int mkeep2 = 0;
      for (int m = 0; m < ncoord_model[1]; m++)
      {
        if (zcylindrical > (((m * dcoord2) * globals::tmin/t_model) - globals::rmax))
        {
          mkeep2 = m;
          // set_cell_modelgridindex(n, m + 1);
        }
      }
      set_cell_modelgridindex(n, (mkeep2 * ncoord_model[0]) + mkeep1);
      modelgrid[get_cell_modelgridindex(n)].initial_radial_pos += radial_pos;

      //renorm[mkeep]++;
    }
    else
    {
      set_cell_modelgridindex(n, get_npts_model());
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
  const int ncount = threedimensional ? ngrid : get_npts_model();
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
    float abundances_in[150] = { 0. };
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

        // radioactive nuclide abundances should have already been set by read_??_model
        set_elem_stable_abund_from_total(mgi, element, elemabundance);
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
      modelgrid[mgi].initradioabund = (float *) calloc(decay::get_num_nuclides(), sizeof(float));
      assert_always(modelgrid[mgi].initradioabund != NULL);

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
  std::ifstream fmodel("model.txt");
  assert_always(fmodel.is_open());

  std::string line;

  // 1st read the number of data points in the table of input model.
  int npts_model_in = 0;
  assert_always(get_noncommentline(fmodel, line));
  std::stringstream(line) >> npts_model_in;

  set_npts_model(npts_model_in);
  ncoord_model[0] = npts_model_in;

  vout_model = (double *) malloc((get_npts_model() + 1) * sizeof(double));

  // Now read the time (in days) at which the model is specified.
  double t_model_days;
  assert_always(get_noncommentline(fmodel, line));
  std::stringstream(line) >> t_model_days;
  t_model = t_model_days * DAY;

  // Now read in the lines of the model. Each line has 5 entries: the
  // cell number (integer) the velocity at outer boundary of cell (float),
  // the mass density in the cell (float), the abundance of Ni56 by mass
  // in the cell (float) and the total abundance of all Fe-grp elements
  // in the cell (float). For now, the last number is recorded but never
  // used.

  std::vector<int> zlist;
  std::vector<int> alist;
  std::vector<std::string> customnuclidecolumns; // vector of strings containing column names
  std::streampos oldpos = fmodel.tellg();  // get position in case we need to undo getline
  std::getline(fmodel, line);
  if (lineiscommentonly(line))
  {
    // custom header line
    std::istringstream iss(line);
    std::string token;
    while (std::getline(iss, token, ' '))
    {
      if (std::all_of(token.begin(), token.end(), isspace)) // skip whitespace tokens
        continue;
      if (token.rfind("X_", 0) != 0) // skip if doesn't start with 'X_'
         continue;
      if (token == "X_Fegroup")
        continue;
      if (token == "X_Ni56")
        continue;
      if (token == "X_Co56")
        continue;
      if (token == "X_Fe52")
        continue;
      if (token == "X_Cr48")
        continue;
      if (token == "X_Ni57")
        continue;
      if (token == "X_Co57")
        continue;

      const int z = decay::get_nucstring_z(token.c_str() + 2); // + 2 skips the 'X_'
      const int a = decay::get_nucstring_a(token.c_str() + 2);
      printout("Custom nuclide column: '%s' Z %d A %d\n", token.c_str(), z, a);
      zlist.push_back(z);
      alist.push_back(a);
    }

    // alternative:
    // while(iss >> token)
    // {
    //   printout("Custom header column: %s\n", token.c_str());
    //   columns.push_back(token);
    // }

    // for(std::vector<std::string>::iterator it = columns.begin(); it != columns.end(); ++it)
    // {
    //   printout("Repeat of Custom header column: %s\n", it->c_str());
    // }
  }
  else
  {
    fmodel.seekg(oldpos); // undo getline
  }

  decay::init_nuclides(zlist, alist);

  int mgi = 0;
  while (std::getline(fmodel, line))
  {
    std::istringstream ssline(line);
    int cellnumberin;
    double vout_kmps;
    double log_rho;
    double ffegrp_model = 0.;

    double f56ni_model = 0.;
    double f56co_model = 0.;
    double f48cr_model = 0.;
    double f52fe_model = 0.;
    double f57ni_model = 0.;
    double f57co_model = 0.;

    const int items_read = sscanf(line.c_str(), "%d %lg %lg %lg %lg %lg %lg %lg %lg %lg",
                                  &cellnumberin, &vout_kmps, &log_rho, &ffegrp_model, &f56ni_model,
                                  &f56co_model, &f52fe_model, &f48cr_model, &f57ni_model, &f57co_model);

    if (items_read == 8 || items_read == 10)
    {
      assert_always(cellnumberin == mgi + 1);

      vout_model[mgi] = vout_kmps * 1.e5;

      const double rho_tmin = pow(10., log_rho) * pow(t_model / globals::tmin, 3);
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
      printout("line: %s\n", line.c_str());
      abort();
    }

    // printout("%d %g %g %g %g %g %g %g\n",
    //          cellnumin, vout_kmps, log_rho, ffegrp_model[n], f56ni_model[n],
    //          f56co_model[n], f52fe_model[n], f48cr_model[n]);
    // printout("   %lg %lg\n", f57ni_model[n], f57co_model[n]);
    modelgrid[mgi].initradioabund = (float *) calloc(decay::get_num_nuclides(), sizeof(float));
    assert_always(modelgrid[mgi].initradioabund != NULL);

    set_modelinitradioabund(mgi, 28, 56, f56ni_model);
    set_modelinitradioabund(mgi, 27, 56, f56co_model);
    set_modelinitradioabund(mgi, 28, 57, f57ni_model);
    set_modelinitradioabund(mgi, 27, 57, f57co_model);
    set_modelinitradioabund(mgi, 26, 52, f52fe_model);
    set_modelinitradioabund(mgi, 24, 48, f48cr_model);
    set_modelinitradioabund(mgi, 23, 48, 0.);
    set_ffegrp(mgi, ffegrp_model);
    for (int i = 0; i < 10; i++)
    {
      double abundin = 0.;
      ssline >> abundin; // ignore
    }
    for (int i = 0; i < (int) zlist.size(); i++)
    {
      double abundin = 0.;
      ssline >> abundin;
      set_modelinitradioabund(mgi, zlist[i], alist[i], abundin);
    }

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

  fmodel.close();

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
  t_model = t_model_days * DAY;

  // Now read in globals::vmax (in cm/s)
  fscanf(model_input, "%lg\n", &globals::vmax);
  dcoord1 = globals::vmax * t_model / ncoord_model[0]; //dr for input model
  dcoord2 = 2. * globals::vmax * t_model / ncoord_model[1]; //dz for input model

  // Now read in the model. Each point in the model has two lines of input.
  // First is an index for the cell then its r-mid point then its z-mid point
  // then its total mass density.
  // Second is the total FeG mass, initial 56Ni mass, initial 56Co mass

  decay::init_nuclides(std::vector<int>(), std::vector<int>());

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
    const double z = -globals::vmax * t_model + ((ncoord2 + 0.5) * dcoord2);
    assert_always(fabs(cell_z_in / z - 1) < 1e-3);

    assert_always(cellnumin == mgi + 1);

    const double rho_tmin = rho_tmodel * pow(t_model / globals::tmin, 3);
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

  set_npts_model(npts_model_in);

  ncoord_model[0] = ncoord_model[1] = ncoord_model[2] = round(pow(npts_model_in, 1/3.));
  assert_always(ncoord_model[0] * ncoord_model[1] * ncoord_model[2] == npts_model_in);

  // for a 3D input model, the progation cells will match the input cells exactly
  ncoordgrid[0] = ncoord_model[0];
  ncoordgrid[1] = ncoord_model[1];
  ncoordgrid[2] = ncoord_model[2];
  ngrid = npts_model_in;
  grid_type = GRID_UNIFORM;

  /// Now read the time (in days) at which the model is specified.
  float t_model_days;
  fscanf(model_input, "%g", &t_model_days);
  t_model = t_model_days * DAY;

  /// Now read in globals::vmax for the model (in cm s^-1).
  fscanf(model_input, "%lg\n", &globals::vmax);

  double xmax_tmodel = globals::vmax * t_model;

  /// Now read in the lines of the model.
  min_den = 1.e99;

  // check if expected positions match in either xyz or zyx column order
  // set false if a problem is detected
  bool posmatch_xyz = true;
  bool posmatch_zyx = true;

  decay::init_nuclides(std::vector<int>(), std::vector<int>());

  // mgi is the index to the model grid - empty cells are sent to special value get_npts_model(),
  // otherwise each input cell is one modelgrid cell
  int mgi = 0; // corresponds to model.txt index column, but zero indexed! (model.txt might be 1-indexed)
  int nonemptymgi = 0;
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

    assert_always(mgi_in == mgi + 1);

    // cell coordinates in the 3D model.txt file are sometimes reordered by the scaling script
    // however, the cellindex always should increment X first, then Y, then Z

    for (int axis = 0; axis < 3; axis++)
    {
      const double cellwidth = 2 * xmax_tmodel / ncoordgrid[axis];
      const double cellpos_expected = - xmax_tmodel + cellwidth * get_cellcoordpointnum(mgi, axis);
      // printout("n %d coord %d expected %g found %g rmax %g get_cellcoordpointnum(n, axis) %d ncoordgrid %d\n",
      // n, axis, cellpos_expected, cellpos_in[axis], xmax_tmodel, get_cellcoordpointnum(n, axis), ncoordgrid[axis]);
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
      printout("negative input density %g %d\n", rho_model, mgi);
      abort();
    }

    // in 3D cartesian, cellindex and modelgridindex are interchangeable
    const bool keepcell = (rho_model > 0);
    if (keepcell)
    {
      set_cell_modelgridindex(mgi, mgi);
      const double rho_tmin = rho_model * pow((t_model / globals::tmin), 3.);
      set_rhoinit(mgi, rho_tmin);
      set_rho(mgi, rho_tmin);

      if (min_den > rho_model)
      {
        min_den = rho_model;
      }
    }
    else
    {
      set_cell_modelgridindex(mgi, get_npts_model());
    }

    read_2d3d_modelradioabundanceline(model_input, mgi, keepcell);
    if (keepcell)
    {
      nonemptymgi++;
    }

    mgi++;
  }
  if (mgi != npts_model_in)
  {
    printout("ERROR in model.txt. Found %d cells instead of %d expected.\n", mgi, npts_model_in);
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
  printout("Effectively used model grid cells %d\n", nonemptymgi);

  /// Now, set actual size of the modelgrid to the number of non-empty cells.

  fclose(model_input);
}


static void calc_totmassradionuclides(void)
{
  mtot = 0.;
  mfeg = 0.;

  assert_always(totmassradionuclide == NULL);
  totmassradionuclide = (double *) malloc(decay::get_num_nuclides() * sizeof(double));
  assert_always(totmassradionuclide != NULL);

  for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++)
  {
    totmassradionuclide[nucindex] = 0.;
  }

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
      cellvolume = pow(globals::tmin / t_model, 3) * ((2 * n1) + 1) * PI * dcoord2 * pow(dcoord1, 2.);
      n1++;
      if (n1 == ncoord_model[0])
      {
        n1 = 0;
      }
    }
    else if (get_model_type() == RHO_3D_READ)
    {
      /// Assumes cells are cubes here - all same volume.
      cellvolume = pow((2 * globals::vmax * globals::tmin), 3.) / (ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2]);
    }
    else
    {
      printout("Unknown model type %d in function %s\n", get_model_type(), __func__);
      abort();
    }

    const double mass_in_shell =get_rhoinit(mgi) * cellvolume;

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


void read_ejecta_model(void)
{
  switch (get_model_type())
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

  globals::compton_emiss = (float *) malloc((get_npts_model() + 1) * EMISS_MAX * sizeof(float));
  globals::rpkt_emiss = (double *) malloc((get_npts_model() + 1) * sizeof(double));

  #if (!NO_LUT_PHOTOION)
  globals::corrphotoionrenorm = (double *) malloc((get_npts_model() + 1) * get_nelements() * get_max_nions() * sizeof(double));
  globals::gammaestimator = (double *) malloc((get_npts_model() + 1) * get_nelements() * get_max_nions() * sizeof(double));

  #ifdef DO_TITER
  globals::gammaestimator_save = (double *) malloc((get_npts_model() + 1) * get_nelements() * get_max_nions() * sizeof(double));
  #endif
  #endif

  #if (!NO_LUT_BFHEATING)
  globals::bfheatingestimator = (double *) malloc((get_npts_model() + 1) * get_nelements() * get_max_nions() * sizeof(double));
  #ifdef DO_TITER
  globals::bfheatingestimator_save = (double *) malloc((get_npts_model() + 1) * get_nelements() * get_max_nions() * sizeof(double));
  #endif
  #endif

  #ifndef FORCE_LTE
  globals::ffheatingestimator = (double *) malloc((get_npts_model() + 1) * sizeof(double));
  globals::colheatingestimator = (double *) malloc((get_npts_model() + 1) * sizeof(double));
  #ifdef DO_TITER
  globals::ffheatingestimator_save = (double *) malloc((get_npts_model() + 1) * sizeof(double));
  globals::colheatingestimator_save = (double *) malloc((get_npts_model() + 1) * sizeof(double));
  #endif
  #endif

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
  assert_always(timestep_in == timestep);

  for (int mgi = 0; mgi < get_npts_model(); mgi++)
  {
    int mgi_in;
    float T_R;
    float T_e;
    float W;
    float T_J;
    fscanf(gridsave_file, "%d %g %g %g %g %hd %lg",
           &mgi_in, &T_R, &T_e, &W, &T_J,
           &modelgrid[mgi].thick, &globals::rpkt_emiss[mgi]);

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
              modelgrid[mgi].thick, globals::rpkt_emiss[mgi]);
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


static void assign_initial_temperatures(void)
/// Routine for assigning temperatures to the grid cells at the start of the simulation.
{
  /// For a simulation started from scratch we estimate the initial temperatures

  /// We assume that for early times the material is so optically thick, that
  /// all the radiation is trapped in the cell it originates from. This
  /// means furthermore LTE, so that both temperatures can be evaluated
  /// according to the local energy density resulting from the 56Ni decay.
  /// The dilution factor is W=1 in LTE.

  const double tstart = globals::time_step[0].mid;

  for (int mgi = 0; mgi < get_npts_model(); mgi++)
  {
    const double decayedenergy_per_mass = decay::get_endecay_per_ejectamass_t0_to_time_withexpansion(mgi, tstart);

    double T_initial = pow(CLIGHT / 4 / STEBO * pow(globals::tmin / tstart, 3) *get_rhoinit(mgi) * decayedenergy_per_mass, 1. / 4.);

    // printout("T_initial %g T_initial_alt %g\n", T_initial, T_initial_alt);

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
    modelgrid[mgi].thick = 0;
  }
}


void get_nstart_ndo(int my_rank, int nprocesses, int *nstart, int *ndo, int *maxndo)
{
  #ifndef MPI_ON
    // no MPI, single process updates all cells
    *nstart = 0;
    *ndo = get_npts_model();
    return;
  #endif

  int n_leftover = 0;

  int nblock = get_npts_model() / nprocesses; // integer division, minimum cells for any process
  const int numtot = nblock * nprocesses; // cells counted if all processes do the minimum number of cells
  if (numtot > get_npts_model()) // LJS: should never be the case?
  {
    nblock = nblock - 1;
    *maxndo = nblock + 1;
    n_leftover = get_npts_model() - (nblock * nprocesses);
  }
  else if (numtot < get_npts_model())
  {
    *maxndo = nblock + 1;
    n_leftover = get_npts_model() - (nblock * nprocesses);
  }
  else
  {
    *maxndo = nblock;
    n_leftover = 0;
  }

  if (my_rank < n_leftover)
  {
    *ndo = nblock + 1;
    *nstart = my_rank * (nblock + 1);
  }
  else
  {
    *ndo = nblock;
    *nstart = n_leftover * (nblock + 1) + (my_rank - n_leftover) * (nblock);
  }
}


static void uniform_grid_setup(void)
/// Routine for doing a uniform cuboidal grid.
{

  /// Set grid size for uniform xyz grid
  if (get_model_type() == RHO_3D_READ)
  {
    // if we used in a 3D ejecta model, the propagation grid must match the input grid exactly
    ncoordgrid[0] = ncoord_model[0];
    ncoordgrid[1] = ncoord_model[1];
    ncoordgrid[2] = ncoord_model[2];

    // in case the user specified a grid size, we should ensure that it matches
    #ifdef CUBOID_NCOORDGRID_X
    assert(ncoordgrid[0] == CUBOID_NCOORDGRID_X);
    #endif
    #ifdef CUBOID_NCOORDGRID_Y
    assert(ncoordgrid[1] == CUBOID_NCOORDGRID_Y);
    #endif
    #ifdef CUBOID_NCOORDGRID_Z
    assert(ncoordgrid[2] == CUBOID_NCOORDGRID_Z);
    #endif
  }
  else
  {
    #ifdef CUBOID_NCOORDGRID_X
    ncoordgrid[0] = CUBOID_NCOORDGRID_X;
    #else
    ncoordgrid[0] = 50;
    #endif
    #ifdef CUBOID_NCOORDGRID_Y
    ncoordgrid[1] = CUBOID_NCOORDGRID_Y;
    #else
    ncoordgrid[1] = 50;
    #endif
    #ifdef CUBOID_NCOORDGRID_Z
    ncoordgrid[2] = CUBOID_NCOORDGRID_Z;
    #else
    ncoordgrid[2] = 50;
    #endif
  }

  // artis assumes in some places that the cells are cubes, not cubioids
  assert_always(ncoordgrid[0] == ncoordgrid[1]);
  assert_always(ncoordgrid[0] == ncoordgrid[2]);

  ngrid = ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2];
  cell = (CELL *) malloc(ngrid * sizeof(CELL));

  coordlabel[0] = 'X';
  coordlabel[1] = 'Y';
  coordlabel[2] = 'Z';
  int nxyz[3] = {0, 0, 0};
  for (int n = 0; n < ngrid; n++)
  {
    for (int axis = 0; axis < 3; axis++)
    {
      assert_always(nxyz[axis] == get_cellcoordpointnum(n, axis));
      cell[n].pos_init[axis] = - globals::coordmax[axis] + (2 * nxyz[axis] * globals::coordmax[axis] / ncoordgrid[axis]);
      // cell[n].xyz[axis] = nxyz[axis];
    }

    assert_always(n == nxyz[2] * ncoordgrid[1] * ncoordgrid[2] + nxyz[1] * ncoordgrid[0] + nxyz[0]);

    nxyz[0]++;  // increment x coordinate
    if (nxyz[0] == ncoordgrid[0])
    {
      nxyz[0] = 0;
      nxyz[1]++;  // increment y coordinate
    }
    if (nxyz[1] == ncoordgrid[1])
    {
      nxyz[1] = 0;
      nxyz[2]++;  // increment z coordinate
    }
  }
}


static void spherical1d_grid_setup(void)
{
  assert_always(get_model_type() == RHO_1D_READ);
  coordlabel[0] = 'r';
  coordlabel[1] = '_';
  coordlabel[2] = '_';

  ncoordgrid[0] = get_npts_model();
  ncoordgrid[1] = 1;
  ncoordgrid[2] = 1;

  ngrid = ncoordgrid[0] * ncoordgrid[1] * ncoordgrid[2];
  cell = (CELL *) malloc(ngrid * sizeof(CELL));

  globals::coordmax[0] = globals::rmax;
  globals::coordmax[1] = 0.;
  globals::coordmax[2] = 0.;

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
  for (int n = 0; n <= get_npts_model(); n++)
  {
    modelgrid[n].initial_radial_pos = 0;
  }

  /// Select grid type
  #ifdef GRID_TYPE
  grid_type = GRID_TYPE;
  #else
  grid_type = GRID_UNIFORM;
  #endif

  /// The cells will be ordered by x then y, then z. Call a routine that
  /// sets up the initial positions and widths of the cells.
  char grid_type_name[256] = "";
  if (grid_type == GRID_UNIFORM)
  {
    uniform_grid_setup();
    strcpy(grid_type_name, "uniform cuboidal");
  }
  else if (grid_type == GRID_SPHERICAL1D)
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
    printout("    coordinate %d '%c': cells have %d position values\n", d, coordlabel[d], ncoordgrid[d]);
  }
  printout("    total propagration cells: %d\n", ngrid);

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
    assert_always(grid_type == GRID_UNIFORM);
    // propagation grid must match the input model grid exactly for 3D models
    assert_always(ncoord_model[0] == ncoordgrid[0]);
    assert_always(ncoord_model[1] == ncoordgrid[1]);
    assert_always(ncoord_model[2] == ncoordgrid[2]);

    for (int n = 0; n < ngrid; n++)
    {
      modelgrid[get_cell_modelgridindex(n)].initial_radial_pos = get_cellradialpos(n);
    }
    // cells with rho > 0 are allocated by the above function
  }
  else
  {
    printout("[fatal] grid_init: Error: Unknown density type. Abort.");
    abort();
  }

  allocate_nonemptymodelcells();
  calculate_kappagrey();
  abundances_read();

  int nstart = 0;
  int ndo = 0;
  int maxndo = 0;
  get_nstart_ndo(my_rank, globals::nprocs, &nstart, &ndo, &maxndo);

  radfield::init(my_rank, ndo);
  nonthermal::init(my_rank, ndo);

  /// and assign a temperature to the cells
  if (globals::simulation_continued_from_saved)
  {
    /// For continuation of an existing simulation we read the temperatures
    /// at the end of the simulation and write them to the grid.
    read_grid_restart_data(globals::itstep);
  }
  else
  {
    assign_initial_temperatures();
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
      totmassradionuclide_actual += get_modelinitradioabund(mgi, z, a) *get_rhoinit(mgi) * vol_init_modelcell(mgi);
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

}