#include "sn3d.h"
#include "grid.h"
#include "globals.h"


static void read_energy_in_cells_1d(void)
{
  /// energy released in simulation during start time and end time
//  float start_time; //todo: change to model read in time
//  float end_time; //todo: do I need this?

  int number_of_cells; // number of model grid cells
  int cellnumber; // dummy value - this isn't saved

  FILE *cell_energies_file = fopen_required("energydistribution.txt", "r");

  fscanf(cell_energies_file, "%d", &number_of_cells);
  if (number_of_cells > MMODELGRID)
  {
    printout("too many cells in energy file (%d) > MMODELGRID %d",
             number_of_cells, MMODELGRID); //todo: is there a better check?
    exit(0);
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

//  for (int mgi = 0; mgi < get_npts_model(); mgi++)
//  {
//    printout("mgi %d cells %d \n", mgi, get_numassociatedcells(mgi));
//  }
//  exit(0);
  read_energy_in_cells_1d();
  read_energy_file();
}
