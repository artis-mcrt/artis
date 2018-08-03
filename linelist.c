#include "sn3d.h"
#include "linelist.h"



// To construct an energy ordered gamma ray line list.
void get_gam_ll(void)
{
  /* Start by setting up the grid of fake lines and their energies. */
  gamma_spectra[FAKE_GAM_LINE_ID].nlines = nfake_gam;
  gamma_spectra[FAKE_GAM_LINE_ID].energy = (double *) malloc(nfake_gam * sizeof(double));
  gamma_spectra[FAKE_GAM_LINE_ID].probability = (double *) malloc(nfake_gam * sizeof(double));

  const double deltanu = (nusyn_max - nusyn_min) / (gamma_spectra[FAKE_GAM_LINE_ID].nlines - 3);
  for (int i = 0; i < gamma_spectra[FAKE_GAM_LINE_ID].nlines; i++)
  {
    gamma_spectra[FAKE_GAM_LINE_ID].energy[i] = (nusyn_min + deltanu * (i - 1)) * H;
    gamma_spectra[FAKE_GAM_LINE_ID].probability[i] = 0.0;
  }

  /* Now do the sorting. */

  int total_lines = 0;
  for (enum radionuclides iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
  {
    total_lines += gamma_spectra[iso].nlines;
  }
  printout("total gamma-ray lines %d\n", total_lines);

  gam_line_list.total = total_lines;
  gam_line_list.nuclidetype = (enum radionuclides *) malloc(total_lines * sizeof(enum radionuclides));
  gam_line_list.index = (int *) malloc(total_lines * sizeof(int));

  double energy_last = 0.0;
  int next = -99;
  enum radionuclides next_type = -99;

  for (int i = 0; i < total_lines; i++)
  {
    double energy_try = 1.e50;

    for (enum radionuclides iso = 0; iso < RADIONUCLIDE_COUNT; iso++)
    {
      // printout("iso %d nlines %d\n", iso, gamma_spectra[iso].nlines);
      for (int j = 0; j < gamma_spectra[iso].nlines; j++)
      {
        if (gamma_spectra[iso].energy[j] > energy_last && gamma_spectra[iso].energy[j] < energy_try)
        {
          // next_type = spec_type[iso];
          next_type = iso;
          next = j;
          energy_try = gamma_spectra[iso].energy[j];
        }
      }
    }

    gam_line_list.nuclidetype[i] = next_type;
    gam_line_list.index[i] = next;
    energy_last = energy_try;
  }

  FILE *const line_list = fopen_required("gammalinelist.out", "w+");

  for (int i = 0; i < total_lines; i++)
  {
    const enum radionuclides iso = gam_line_list.nuclidetype[i];
    const int index = gam_line_list.index[i];
    fprintf(line_list, "%d %d %d %g %g \n",
            i, gam_line_list.nuclidetype[i], gam_line_list.index[i],
            gamma_spectra[iso].energy[index] / MEV, gamma_spectra[iso].probability[index]);
  }
  fclose(line_list);
}
