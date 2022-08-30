#include "grey_emissivities.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <cstring>

#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "packet.h"
#include "photo_electric.h"
#include "sn3d.h"
#include "vectors.h"

static double meanf_sigma(const double x)
// Routine to compute the mean energy converted to non-thermal electrons times
// the Klein-Nishina cross section.
{
  double f = 1 + (2 * x);

  double term0 = 2 / x;
  double term1 = (1 - (2 / x) - (3 / (x * x))) * log(f);
  double term2 = ((4 / x) + (3 / (x * x)) - 1) * 2 * x / f;
  double term3 = (1 - (2 / x) - (1 / (x * x))) * 2 * x * (1 + x) / f / f;
  double term4 = -2. * x * ((4 * x * x) + (6 * x) + 3) / 3 / f / f / f;

  double tot = 3 * SIGMA_T * (term0 + term1 + term2 + term3 + term4) / (8 * x);

  return tot;
}

void rlc_emiss_gamma(const struct packet *pkt_ptr, const double dist) {
  // Subroutine to record the heating rate in a cell due to gamma rays.
  // By heating rate I mean, for now, really the rate at which the code is making
  // k-packets in that cell which will then convert into r-packets. This is (going
  // to be) used for the new light_curve syn-style calculation.

  // The intention is that rpkt_emiss will contain the emissivity of r-packets
  // in the co-moving frame (which is going to be isotropic).

  // This is only done to order v/c for now.

  // Called with a packet that is about to travel a
  // distance dist in the lab frame.

  // struct packet dummy;
  // dummy.pos[0] = pkt_ptr->pos[0];
  // dummy.pos[1] = pkt_ptr->pos[1];
  // dummy.pos[2] = pkt_ptr->pos[2];
  // dummy.dir[0] = syn_dir[0];
  // dummy.dir[1] = syn_dir[1];
  // dummy.dir[2] = syn_dir[2];
  // dummy.where = pkt_ptr->where;
  // dummy.last_cross = NONE;

  const int cellindex = pkt_ptr->where;
  const int mgi = grid::get_cell_modelgridindex(cellindex);

  if (dist > 0) {
    double vel_vec[3];
    get_velocity(pkt_ptr->pos, vel_vec, pkt_ptr->prop_time);
    const double xx = H * pkt_ptr->nu_cmf / ME / CLIGHT / CLIGHT;

    double heating_cont = ((meanf_sigma(xx) * grid::get_nnetot(mgi)) + sig_photo_electric(pkt_ptr) +
                           (sig_pair_prod(pkt_ptr) * (1. - (2.46636e+20 / pkt_ptr->nu_cmf))));
    heating_cont = heating_cont * pkt_ptr->e_rf * dist * (1. - (2. * dot(vel_vec, pkt_ptr->dir) / CLIGHT));

    // The terms in the above are for Compton, photoelectric and pair production. The pair production one
    // assumes that a fraction (1. - (1.022 MeV / nu)) of the gamma's energy is thermalised.
    // The remaining 1.022 MeV is made into gamma rays

    // For normalisation this needs to be
    //  1) divided by volume
    //  2) divided by the length of the time step
    //  3) divided by 4 pi sr
    //  This will all be done later
    assert_testmodeonly(heating_cont >= 0.);
    assert_testmodeonly(isfinite(heating_cont));
    safeadd(globals::rpkt_emiss[mgi], 1.e-20 * heating_cont);
  }
}

void rlc_emiss_rpkt(const struct packet *pkt_ptr, double dist) {
  // Subroutine to record the rate of destruction (and re-creation) of
  // r-packets by the grey opacity.

  // This is only done to order v/c for now.

  // Called with a packet that is about to travel a
  //    distance dist in the lab frame.

  /*struct packet dummy;

  dummy.pos[0] = pkt_ptr->pos[0];
  dummy.pos[1] = pkt_ptr->pos[1];
  dummy.pos[2] = pkt_ptr->pos[2];
  dummy.dir[0] = syn_dir[0];
  dummy.dir[1] = syn_dir[1];
  dummy.dir[2] = syn_dir[2];
  dummy.where = pkt_ptr->where;
  dummy.last_cross = NONE;*/

  const int cellindex = pkt_ptr->where;
  const int mgi = grid::get_cell_modelgridindex(cellindex);

  if (dist > 0.0) {
    // for the weighted estimators version

    double vel_vec[3];
    get_velocity(pkt_ptr->pos, vel_vec, pkt_ptr->prop_time);

    double cont = (grid::get_kappagrey(mgi) * grid::get_rho(mgi));
    cont = cont * pkt_ptr->e_rf * dist * (1. - (2. * dot(vel_vec, pkt_ptr->dir) / CLIGHT));

    /* For normalisation this needs to be
       1) divided by volume
       2) divided by the length of the time step
       3) divided by 4 pi sr
       This will all be done later
    */

    assert_testmodeonly(cont >= 0.);
    assert_testmodeonly(isfinite(cont));
    safeadd(globals::rpkt_emiss[mgi], 1.e-20 * cont);
  }
}

void normalise_grey(int nts) {
  const double dt = globals::time_step[nts].width;
  globals::time_step[nts].gamma_dep_pathint = 0.;
  for (int mgi = 0; mgi < grid::get_npts_model(); mgi++) {
    if (grid::get_numassociatedcells(mgi) > 0) {
      const double dV = grid::vol_init_modelcell(mgi) * pow(globals::time_step[nts].mid / globals::tmin, 3);

      globals::time_step[nts].gamma_dep_pathint += globals::rpkt_emiss[mgi] * 1.e20 / globals::nprocs;

      globals::rpkt_emiss[mgi] = globals::rpkt_emiss[mgi] * ONEOVER4PI / dV / dt / globals::nprocs;

      // assert_testmodeonly(globals::rpkt_emiss[mgi] >= 0.);
      assert_testmodeonly(isfinite(globals::rpkt_emiss[mgi]));
    }
  }
}

void write_grey(int nts) {
  FILE *est_file = NULL;
  FILE *dummy = NULL;
  char chch;
  char filename[128] = "grey_est_";
  char junk[128];

  if ((dummy = fopen("dummy", "w+")) == NULL) {
    printout("Cannot open dummy.\n");
    abort();
  }
  fprintf(dummy, "%d", nts);
  fclose(dummy);
  if ((dummy = fopen("dummy", "r")) == NULL) {
    printout("Cannot open dummy.\n");
    abort();
  }
  int i = 0;
  while ((chch = fgetc(dummy)) != EOF) {
    junk[i] = chch;
    i = i + 1;
  }
  junk[i] = '\0';
  fclose(dummy);

  strcat(filename, junk);
  strcat(filename, ".out");

  if (globals::file_set) {
    if ((est_file = fopen(filename, "r")) == NULL) {
      printout("Cannot open grey_est_file.txt.\n");
      abort();
    }

    // for (n=0; n < ngrid; n++)
    for (int n = 0; n < grid::get_npts_model(); n++) {
      float dum;
      fscanf(est_file, "%g", &dum);
      globals::rpkt_emiss[n] += dum;
    }
    fclose(est_file);
  }

  if ((est_file = fopen(filename, "w+")) == NULL) {
    printout("Cannot open grey_est_file.txt.\n");
    abort();
  }

  // for (n=0; n < ngrid; n++)
  for (int n = 0; n < grid::get_npts_model(); n++) {
    fprintf(est_file, " %g\n ", globals::rpkt_emiss[n]);
  }
  fclose(est_file);
}
