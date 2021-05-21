#include "sn3d.h"
#include "grid.h"
#include "globals.h"


__managed__ double *modelcell_energy = NULL; // energy in model grid cell read from energydistribution.txt
__managed__ double etot_fromenergyfile = 0.; // total model energy -- used to initialise pellets. Read from energydistribution.txt

__managed__ int ntimes_energydep = 0; // number of times included in energyrate.txt
__managed__ float *time_energydep = NULL; // times in seconds from energyrate.txt
__managed__ float *energy_fraction_deposited = NULL; // fraction of energy deposited by time from energyrate.txt
__managed__ float *modelcell_energydensity_init = NULL;



double get_etot_fromenergyfile(void)
{
  return etot_fromenergyfile;
}


float get_modelcell_energydensity_init(int modelgridindex)
{
  printout("mgi %d", modelgridindex);
  assert_always(modelgridindex <= MMODELGRID);
  return modelcell_energydensity_init[modelgridindex];
}


void set_modelcell_energydensity_init(int modelgridindex, float x)
{
  assert_always(modelcell_energydensity_init != NULL);
  assert_always(modelgridindex < MMODELGRID);
  if (x == 0)
  {
    assert_always(modelgridindex == MMODELGRID)
  }
  modelcell_energydensity_init[modelgridindex] = x;
}


static void read_energy_in_cells_3d(void)
{
  /// energy released in simulation during start time and end time
//  float start_time; //todo: change to model read in time
//  float end_time; //todo: do I need this?

  int number_of_cells; // number of model grid cells

  FILE *cell_energies_file = fopen_required("energydistribution.txt", "r");
  printout("read energydistribution.txt \n");

  fscanf(cell_energies_file, "%d", &number_of_cells);
//  if (number_of_cells != get_npts_model())  //Might not work for 3D - in classic it throws away empty cells for npts_model
//  {
//    printout("number of cells in energy file (%d) "
//             "does not match number of model grid cells (%d) - abort\n",
//             number_of_cells, get_npts_model());
//    abort();
//  }

  int cellnumber; // dummy value - this isn't saved
  double cell_energy;
  double energy_counter = 0.;
  int n, mgi;
  for (n = 0; n < globals::ngrid; n++)
  {
    mgi = get_cell_modelgridindex(n);
    fscanf(cell_energies_file, "%d %lf", &cellnumber, &cell_energy);
    printout("cellnumber %d cell energy %g mgi %d MMODELGRID %d npointsmodel %d\n", cellnumber, cell_energy, mgi, MMODELGRID, get_npts_model());
    if(cellnumber-1 != n)
    {
      printout("cell number in file does not match: cellnumber %d n %d \n", cellnumber, n);
      exit(0);
    }
    if (mgi != MMODELGRID)
    {
      energy_counter += cell_energy;
      modelcell_energy[mgi] = cell_energy;
    }
    else
    {
//      empty cells must have 0 energy
      modelcell_energy[mgi] = 0.;
    }
  }
  modelcell_energy[get_npts_model()] = 0.; // special empty model cell

  etot_fromenergyfile = energy_counter;
  printout("etot %g \n", etot_fromenergyfile);
}


static void read_energy_file(void)
{
  FILE *energyrate_file = fopen_required("energyrate.txt", "r");

//  int ntimes_energydep = 0;
  // number of times included in file
  fscanf(energyrate_file, "%d", &ntimes_energydep);
//  assert_always(ntimes_energydep <= globals::ntstep);  // todo: put this back when length array decided

  time_energydep = (float *) calloc(ntimes_energydep, sizeof(float));

  /// read times and fraction of energy deposited
  float time_energydep_days; //file times in days - convert to seconds
  for (int t_iter = 0; t_iter < ntimes_energydep; t_iter++)
  {
    fscanf(energyrate_file, "%g %g",
           &time_energydep_days, &energy_fraction_deposited[t_iter]);
    time_energydep[t_iter] = time_energydep_days * DAY;
  }

//  printout("%d \n", ntimes_energydep);
//  for (int t_iter = 0; t_iter < ntimes_energydep; t_iter++)
//  {
//    printout("time %g fraction deposited %g \n",
//             time_energydep[t_iter], energy_fraction_deposited[t_iter]);
//  }
}


void energy_input_init(void)
{
  printout("reading energy files npts model %d \n", get_npts_model());

//  modelcell_energy = (double *) calloc((get_npts_model() + 1), sizeof(double));
//  energy_fraction_deposited = (float *) calloc(globals::ntstep, sizeof(float));
//  modelcell_energydensity_init = (float *) calloc((get_npts_model() + 1), sizeof(float));

  modelcell_energy = (double *) calloc((MMODELGRID), sizeof(double));
  energy_fraction_deposited = (float *) calloc(300, sizeof(float));  //haven't decided how long this needs to be yet
  modelcell_energydensity_init = (float *) calloc((MMODELGRID), sizeof(float));

  read_energy_in_cells_3d();
  read_energy_file();
}


void setup_generic_pellet(const double e0, const int mgi, PKT *pkt_ptr)
{
  // I set TYPE_GENERIC_ENERGY_PELLET = 200, but let me know how you want to do this
  pkt_ptr->type = TYPE_GENERIC_ENERGY_PELLET;

  /// Choose decay time
  double zrand = gsl_rng_uniform(rng);
  // randomly sample what fraction of energy has been deposited
  // then find time by which that fraction has been deposited
  int ii = 0;
  while (energy_fraction_deposited[ii] < zrand)
  {
    ii++;
  }
  assert_always(ii < ntimes_energydep);
  pkt_ptr->tdecay = time_energydep[ii];
  //    printout("tdecay pellet %g pkt number %d\n", pkt[n].tdecay/DAY, n);
  //    printout("ii %d zrand %g energy_fraction_deposited %g time_energydep %g \n",
  //             ii, zrand, energy_fraction_deposited[ii], time_energydep[ii]/DAY);
}
