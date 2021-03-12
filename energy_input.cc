#include "sn3d.h"
#include "grid.h"
#include "globals.h"


__managed__ double *modelcell_energy; // energy in model grid cell read from energydistribution.txt
__managed__ double etot_fromenergyfile; // total model energy -- used to initialise pellets. Read from energydistribution.txt

__managed__ int ntimes_energydep; // number of times included in energyrate.txt
__managed__ float *time_energydep; // times in seconds from energyrate.txt
__managed__ float *energy_fraction_deposited; // fraction of energy deposited by time from energyrate.txt


static void read_energy_in_cells_1d(void)
{
  /// energy released in simulation during start time and end time
//  float start_time; //todo: change to model read in time
//  float end_time; //todo: do I need this?

  int number_of_cells; // number of model grid cells
  int cellnumber; // dummy value - this isn't saved

  FILE *cell_energies_file = fopen_required("energydistribution.txt", "r");

  fscanf(cell_energies_file, "%d", &number_of_cells);
  if (number_of_cells != get_npts_model())
  {
    printout("number of cells in energy file (%d) "
             "does not match number of model grid cells (%d) - abort",
             number_of_cells, get_npts_model());
    abort();
  }

  double cell_energies[number_of_cells];
  double energy_counter = 0;
  for (int mgi = 0; mgi < number_of_cells; mgi++)
  {
    fscanf(cell_energies_file, "%d %lf",
           &cellnumber, &cell_energies[mgi]);
    energy_counter += cell_energies[mgi];
    modelcell_energy[mgi] = cell_energies[mgi];
//    printout("modelcell_energy %g get_volinit_modelcell %g \n",
//             modelcell_energy[mgi], get_volinit_modelcell(mgi));
  }
  etot_fromenergyfile = energy_counter;
  //todo: broken here -- modelcell_energy and etot_fromenergyfile not defined properly.
  // trying to set them in globals.h

}


static void read_energy_file(void)
{
  FILE *energyrate_file = fopen_required("energyrate.txt", "r");

  int ntimes_energydep = 0;
  // number of times included in file
  fscanf(energyrate_file, "%d", &ntimes_energydep);

  if (ntimes_energydep > globals::ntstep)
  {
    printout("number of times in file > ntstep - abort'n");
    // Arrays time_energydep[] and energy_fraction_deposited[]
    // defined using MTSTEP - if needs to be longer? redefine
    // with new variable
    abort();
  }

  /// read times and fraction of energy deposited
  double time_energydep[ntimes_energydep];
  double energy_fraction_deposited[ntimes_energydep];
  float time_energydep_days; //file times in days - convert to seconds
  for (int t_iter = 0; t_iter < ntimes_energydep; t_iter++)
  {
    fscanf(energyrate_file, "%g %lg",
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
  printout("reading energy files \n");

  modelcell_energy = (double *) calloc(get_npts_model(), sizeof(double));
  time_energydep = (float *) calloc(globals::ntstep, sizeof(float));
  energy_fraction_deposited = (float *) calloc(globals::ntstep, sizeof(float));

  read_energy_in_cells_1d();
  read_energy_file();
}
