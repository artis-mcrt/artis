#include "sn3d.h"

int energy_init()
{
  void energy_in_cells_1d_read();
  void read_energy_file();

//  printout("ngrid %d \n", ngrid);
//  for (int mgi = 0; mgi < npts_model; mgi++)
//  {
//    printout("volume %g \n", get_volinit_modelcell(mgi));
//    printout("mgi %d cells %d \n", mgi, modelgrid[mgi].associated_cells);
//  }
//  exit(0);

  energy_in_cells_1d_read();
  read_energy_file();
  return(0);
}

void energy_in_cells_1d_read()
{
#ifdef USE_ENERGYINPUTFILE
  FILE *cell_energies_file;
  /// energy released in simulation during start time and end time
//  float start_time; //todo: change to model read in time
//  float end_time; //todo: do I need this?

  /// Open the file
  if ((cell_energies_file = fopen("energydistribution.txt",
                                  "r")) == NULL)
  {
    printout("Cannot open energydistribution.txt.\n");
    exit(0);
  }

  int number_of_cells; // number of model grid cells
  fscanf(cell_energies_file, "%d", &number_of_cells);
  if (number_of_cells > MMODELGRID)
  {
    printout("too many cells in energy file (%d) > MMODELGRID %d",
             number_of_cells, MMODELGRID); //todo: is there a better check?
//             "does not match number of model grid cells (%d) - abort",

    exit(0);
  }  //todo: why does this not match

  int cellnumber; // dummy value - this isn't saved
  double cell_energy;
  double energy_counter = 0.;
  //  for (int mgi = 0; mgi < number_of_cells; mgi++)
  int n, mgi;
  for (n = 0; n < ngrid; n++)
  {
//    printout("cell n %d mgi %d\n", n, cell[n].modelgridindex);
    mgi = cell[n].modelgridindex;
    fscanf(cell_energies_file, "%d %lf", &cellnumber, &cell_energy);
//    printout("cellnumber %d cell energy %g\n", cellnumber, cell_energy);
    if(cellnumber-1 != n)
    {
      printout("cell number in file does not match \n");
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
  etot_fromenergyfile = energy_counter;

//  printout("start %d end %d ncells %d \n",
//           start_time, end_time, number_of_cells);
//  for (int mgi = 0; mgi < number_of_cells; mgi++)
//  {
//    printout("%d %g \n", mgi, modelcell_energy[mgi]);
//  }
//  exit(0);
#endif
}


void read_energy_file()
{
#ifdef USE_ENERGYINPUTFILE
  FILE *energyrate_file;

  /// Open the file
  if ((energyrate_file = fopen("energyrate.txt",
                                  "r")) == NULL)
  {
    printout("Cannot open energyrate.txt.\n");
    exit(0);
  }

  // number of times included in file
  fscanf(energyrate_file, "%d", &ntimes_energydep);

  if (ntimes_energydep > MTSTEP)
  {
    printout("number of times in file > MTSTEP - abort");
    // Arrays time_energydep[] and energy_fraction_deposited[]
    // defined using MTSTEP - if needs to be longer? redefine
    // with new variable
    exit(0);
  }

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
#endif
}