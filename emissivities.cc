#include "sn3d.h"
#include <cstring>
#include "atomic.h"
#include "grid.h"
#include "emissivities.h"
#include "packet_init.h"
#include "photo_electric.h"
#include "radfield.h"
#include "stats.h"
#include "gamma.h"
#include "vectors.h"


void compton_emiss_cont(const PKT *pkt_ptr, double dist)
{
  // Subroutine to add contribution to the MC estimator for the
  // compton emissivity. Called with a packet that is about to travel a
  // distance dist in the lab frame.

  double vel_vec[3];
  double cmf_dir[3];
  double cmf_syn_dir[3];

  // First we need to know the scattering angle needed from the
  // packet's direction of motion to the desired observer. Call this angle
  // mu_cmf (it's a cosine). To get it convert both the direction of
  // motion and the local velocity vectors to the cmf.

  get_velocity(pkt_ptr->pos, vel_vec, pkt_ptr->prop_time);
  angle_ab(pkt_ptr->dir, vel_vec, cmf_dir);
  angle_ab(globals::syn_dir, vel_vec, cmf_syn_dir);

  //  printout("pos %g %g %g\n", pkt_ptr->pos[0],pkt_ptr->pos[1], pkt_ptr->pos[2]);
  //  printout("dir %g %g %g\n", pkt_ptr->dir[0],pkt_ptr->dir[1], pkt_ptr->dir[2]);
  //  printout("vel %g %g %g\n", vel_vec[0], vel_vec[1], vel_vec[2]);
  //  printout("cmf_dir %g %g %g\n", cmf_dir[0], cmf_dir[1], cmf_dir[2]);
  //  printout("syn_dir %g %g %g\n", syn_dir[0], syn_dir[1], syn_dir[2]);
  //  printout("cmf_syn_dir %g %g %g\n", cmf_syn_dir[0], cmf_syn_dir[1], cmf_syn_dir[2]);

  const double mu_cmf = dot(cmf_dir, cmf_syn_dir);

  if (mu_cmf > 1 || mu_cmf < -1)
  {
    printout("problem with Compton emissivity. Abort.\n");
    abort();
  }

  // Now get the factor by which the frequency will change, f, for going
  // in this direction. f = old energy / new energy - always should be > 1

  const double f = 1 + (H * pkt_ptr->nu_cmf / ME / CLIGHT / CLIGHT * (1. - mu_cmf));

  //  printout("compton reducion factor %g freq %g mu %g\n", f, H*pkt_ptr->nu_cmf/MEV, mu_cmf );

  // Now work out in which frequency bin this'll happen. The scattered
  // light will have frequency (nu_cmf / f) in the cmf frame. And it
  // travels in direction syn_dir in the rf.

  const double freq_out = pkt_ptr->nu_cmf / f; /// doppler(syn_dir, vel_vec);
  // do we want ?/ doppler(syn_dir, vel_vec)

  const int lindex = get_nul(freq_out); // This is the index of the next line to
                                        // the red. The emissivity will go in this
                                        // bin. However, since there's an offset
                                        // in the emissivities, we shift the
                                        // index by that

  // If it's gonna be in a bin of interest, carry on - otherwise leave it.

  // printout("frequency %g\n", freq_out*H/MEV);
  // printout("lindex %d, emiss_max %d, emiss_offset %d\n", lindex, emiss_max, emiss_offset);

  if ((lindex > globals::emiss_offset - 1) && (lindex < globals::emiss_offset + globals::emiss_max - 1))
  {

    // Then get partial crossection dsigma_domega in cmf
    // Coeff is 3 / 16 / PI

    const double dsigma_domega_cmf = 0.0596831 * SIGMA_T / f / f * (f + (1./f) + (mu_cmf * mu_cmf) - 1.);

    //speed = vec_len(vel_vec);
    //solid_angle_factor =  doppler(pkt_ptr->dir, vel_vec) * doppler(pkt_ptr->dir, vel_vec);

    //   pow((1 + (dot(vel_vec, syn_dir)/CLIGHT)),2)
    //   / (1.0 - (speed* speed / CLIGHT / CLIGHT));

    //dsigma_domega_rf = dsigma_domega_cmf //* doppler(pkt_ptr->dir, vel_vec)
    //* solid_angle_factor;

    // so now determine the contribution to the emissivity and which
    // frequency bin it should be in

    const double dop_fac = doppler(pkt_ptr->dir, vel_vec);

    const double emiss_cont = pkt_ptr->e_rf * dsigma_domega_cmf * dist * dop_fac * dop_fac / f;

    // For normalisation this needs to be
    //    1) divided by volume
    //    2) divided by frequency bin size
    //    3) multiplied by the cell electron number density
    //    4) divided by the length of the time step
    //    This will all be done later

    if (lindex < globals::emiss_offset)
    {
      printout("scarily bad error here! %d %d\n", lindex, globals::emiss_offset);
    }
    else
    {
      const int cellindex = pkt_ptr->where;
      safeadd(globals::compton_emiss[get_cell_modelgridindex(cellindex)][lindex - globals::emiss_offset], emiss_cont);
    }

  }
}


void pp_emiss_cont(const PKT *pkt_ptr, double dist)
{
  // New routine for getting a pair production emissivity. Closely based on compton_emiss but simpler. The
  // emissivity itself is stored in the last row of the compton emissivity structure. Idea here is to get something
  // which, when normalised by the volume and time step, will give the energy going into the .511 MeV
  // gamma rays from pair production per unit volume per unit time in the cmf.

  // Called with a packet that is about to travel a
  // distance dist in the lab frame.

  const double emiss_cont = sig_pair_prod(pkt_ptr) * (2.46636e+20 / pkt_ptr->nu_cmf) * pkt_ptr->e_rf * dist;

  // For normalisation this needs to be
  //  1) divided by volume
  //  2) divided by the length of the time step
  //  This will all be done later

  const int cellindex = pkt_ptr->where;
  safeadd(globals::compton_emiss[get_cell_modelgridindex(cellindex)][globals::emiss_max - 1], 1.e-20 * emiss_cont);

  //  printf("emiss_cont %g\n", emiss_cont);

  // Note (SS May 07) - the Doppler factors are not all sorted out yet - the expression used above needs to be
  // consistent with what syn_lc does.
}


void zero_estimators(void)
{
  // for (n=0; n < ngrid; n++)
  for (int n = 0; n < get_npts_model(); n++)
  {
    radfield::zero_estimators(n);

    #ifndef FORCE_LTE
      globals::ffheatingestimator[n] = 0.;
      globals::colheatingestimator[n] = 0.;

      // mabfcount[n] = 0.;
      // mabfcount_thermal[n] = 0.;
      // matotem[n] = 0.;
      // maabs[n] = 0.;
      // kbfcount[n] = 0.;
      // kbfcount_ion[n] = 0.;
      // kffcount[n] = 0.;
      // kffabs[n] = 0.;
      // kbfabs[n] = 0.;
      // kgammadep[n] = 0.;

      if (TRACK_ION_STATS)
      {
        stats::reset_ion_stats(n);
      }

      for (int element = 0; element < get_nelements(); element++)
      {
        for (int ion = 0; ion < globals::maxion; ion++)
        {
          #if (!NO_LUT_PHOTOION)
            globals::gammaestimator[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = 0.;
          #endif
          #if (!NO_LUT_BFHEATING)
            globals::bfheatingestimator[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = 0.;
          #endif

          // photoionestimator[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = 0.;
          // stimrecombestimator[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = 0.;
          // ionfluxestimator[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = 0.;
          // twiddle[n*get_nelements()*get_max_nions()+element*get_max_nions()+ion] = 0.;
        }
      }
    #endif
    for (int m = 0; m < globals::emiss_max; m++)
    {
      globals::compton_emiss[n][m] = 0.0;
    }

    globals::rpkt_emiss[n] = 0.0;
  }
}


void normalise_compton_estimators(const int nts, struct time *time_step)
{
  double dfreq[EMISS_MAX];

  const double time_factor = 1. / pow(globals::time_step[nts].mid / globals::tmin, 3.0) / globals::time_step[nts].width;

  for (int m = 0; m < globals::emiss_max; m++)
  {
    dfreq[m] = get_gam_freq(m + globals::emiss_offset + 1) - get_gam_freq(m + globals::emiss_offset);
    if (dfreq[m] < 0)
    {
      printout("Problem with normalisation of estimators. Abort.\n");
      abort();
    }
    dfreq[m] = 1. / dfreq[m];
  }

  // for (n=0; n < ngrid; n++)
  for (int n = 0; n < get_npts_model(); n++)
  {
    const double volume = vol_init_modelcell(n);
    for (int m = 0; m < globals::emiss_max; m++)
    {
      globals::compton_emiss[n][m] = globals::compton_emiss[n][m] * time_factor / volume / globals::nprocs;

      if (m < globals::emiss_max - 1)
      // (emiss_max - 1) contains the pair production case so it doesn't need the nne nor the dfreq
      {
        globals::compton_emiss[n][m] = globals::compton_emiss[n][m] * get_nne(n) * dfreq[m];
      }
    }
  }
}


void write_compton_estimators(int nts)
{
  FILE *est_file, *dummy;
  char chch;
  char filename[100] = "est_";
  char junk[100];

  dummy = fopen_required("dummy", "w+");
  fprintf(dummy, "%d", nts);
  fclose(dummy);
  dummy = fopen_required("dummy", "r");
  int i = 0;
  while ((chch = fgetc(dummy)) != EOF)
  {
    junk[i] = chch;
    i = i+1;
  }
  junk[i] = '\0';
  fclose(dummy);

  strcat(filename, junk);
  strcat(filename, ".out");

  if (globals::file_set)
  {
    est_file = fopen_required(filename, "rb");

      //for (n=0; n < ngrid; n++)
    for (int n = 0; n < get_npts_model(); n++)
    {
      for (int m = 0; m < globals::emiss_max; m++)
      {
        float dum;
        fread(&dum, sizeof(float), 1, est_file);
        //fscanf(est_file, "%g", &dum);
        globals::compton_emiss[n][m] += dum;
      }
    }
    fclose(est_file);
  }

  est_file = fopen_required(filename, "wb+");

  for (int n = 0; n < get_npts_model(); n++)
  {
    for (int m = 0; m < globals::emiss_max; m++)
    {
      fwrite(&globals::compton_emiss[n][m], sizeof(float), 1, est_file);
    }
  }
  fclose(est_file);
}


bool estim_switch(int nts, struct time *time_step)
{
  const double tstart = globals::time_step[nts].start;
  const double tend = globals::time_step[nts].start + globals::time_step[nts].width;

  const double ts_want = globals::time_syn[0] * ((1. - globals::rmax / globals::tmin / globals::CLIGHT_PROP));
  const double te_want = globals::time_syn[globals::nsyn_time - 1] * (1. + globals::rmax / globals::tmin / globals::CLIGHT_PROP);

  return ((tstart > te_want) || (tend < ts_want));
}


// void emiss_load(int nts)
// // Routine to read in the stored estimators for the time step that is about to begin.
// {
//   FILE *est_file, *dummy;
//   char chch;
//   char filename[100] = "est_";
//   char junk[100];
//
//   dummy = fopen_required("dummy", "w+");
//   fprintf(dummy, "%d", nts);
//   fclose(dummy);
//   dummy = fopen_required("dummy", "r");
//   int i = 0;
//   while ((chch=fgetc(dummy)) != EOF)
//   {
//     junk[i] = chch;
//     i = i+1;
//   }
//   junk[i] = '\0';
//   fclose(dummy);
//
//   strcat(filename, junk);
//   strcat(filename, ".out");
//
//   est_file = fopen_required(filename, "r");
//
//
//   //for (n=0; n < ngrid; n++)
//   for (int n = 0; n < npts_model; n++)
//   {
//     for (int m = 0; m < emiss_max; m++)
//     {
//       float dum;
//       fscanf(est_file, "%g", &dum);
//       compton_emiss[n][m] = dum;
//     }
//   }
//   fclose(est_file);
// }

// static void continuum_rt(RAY *ray_ptr, int nray, double ldist, double *single_pos, double single_t, int lindex)
// // This is called when a ray is about to be moved a distance ldist.
// // It should account for the changes in the ray intensity due to
// //   continuum processes along the path.
// {
//   PKT dummy;
//
//   // Make a dummy packet that carries the ray properties.
//
//   dummy.pos[0] = single_pos[0];
//   dummy.pos[1] = single_pos[1];
//   dummy.pos[2] = single_pos[2];
//
//   dummy.dir[0] = syn_dir[0];
//   dummy.dir[1] = syn_dir[1];
//   dummy.dir[2] = syn_dir[2];
//
//   dummy.where = ray_ptr->where;
//   dummy.nu_cmf = ray_ptr->nu_cmf[nray];
//
//   double kap_compton = sig_comp(&dummy,single_t);
//   double kap_photo_electric = sig_photo_electric(&dummy, single_t);
//   double kap_pair_prod = sig_pair_prod(&dummy, single_t);
//   double kap_tot = kap_compton + kap_photo_electric + kap_pair_prod;
//
//   /* For now no emissivity - only destruction. So very simple. */
//
//   double tau_cont = kap_tot * ldist;
//
//   ray_ptr->e_rf[nray] = ray_ptr->e_rf[nray] * exp(-1. * tau_cont);
//
//   /* Now adding the emissivity term. */
//
//   if (lindex != RED_OF_LIST)
//   {
//     double vel_vec[3];
//     get_velocity(single_pos, vel_vec, single_t);
//     double dop_fac = 1./doppler(syn_dir, vel_vec);
//
//     if (tau_cont > 1.e-6)
//     {
//       ray_ptr->e_rf[nray] += (dop_fac * dop_fac * compton_emiss[get_cell_modelgridindex(dummy.where)][lindex - emiss_offset] *
//           (1. - exp(-1. * tau_cont)) / kap_tot);
//     }
//     else
//     {
//       ray_ptr->e_rf[nray] += (dop_fac * dop_fac * compton_emiss[get_cell_modelgridindex(dummy.where)][lindex - emiss_offset] *
//             ldist);
//     }
//   }
//
//   //    This MUST be followed by a call to move_one_ray() in source
//   //    since e_cmf is NOT reset here.
// }


// void add_gam_line_emissivity(RAY *ray_ptr, int nray, double *single_pos, double single_t, int lindex, double dnuds)
// {
//   double emitt_energy;
//   struct grid *grid_ptr;
//
//   grid_ptr = &globals::cell[ray_ptr->where];
//   double tfact = pow((globals::tmin/single_t), 3);
//
//   if (gam_line_list.type[lindex] == 28, 56)
//   {
//     emitt_energy = get_rhoinit(grid_ptr->modelgridindex) / nucmass(28, 56) / 4. / PI
//         * exp(-single_t/get_meanlife(28, 56)) / get_meanlife(28, 56) *
//         ni56_spec.probability[gam_line_list.index[lindex]]
//         * ni56_spec.energy[gam_line_list.index[lindex]]
//          * f56ni(grid_ptr)*tfact;
//   }
//   else if (gam_line_list.type[lindex] == NUCLIDE_CO56)
//   {
//     emitt_energy = get_rhoinit(grid_ptr->modelgridindex) / nucmass(28, 56) / 4. / PI
//       * (exp(-single_t/get_meanlife(28, 56)) - exp(-single_t/get_meanlife(NUCLIDE_CO56)))
//       / (get_meanlife(28, 56) - get_meanlife(NUCLIDE_CO56))
//       * co56_spec.probability[gam_line_list.index[lindex]]
//       * co56_spec.energy[gam_line_list.index[lindex]]
//       * f56ni(grid_ptr) * tfact;
//
//     if (gam_line_list.index[lindex] == 0)
//     {
//       emitt_energy += (compton_emiss[grid_ptr->modelgridindex][emiss_max - 1] * 1.e20 / 4. / PI);
//     }
//   }
//   else if (gam_line_list.type[lindex] == NUCLIDE_CR48)
//   {
//     emitt_energy = get_rhoinit(grid_ptr->modelgridindex) / nucmass(NUCLIDE_CR48) / 4. / PI
//         * exp(-single_t/get_meanlife(NUCLIDE_CR48)) / get_meanlife(NUCLIDE_CR48)
//         * cr48_spec.probability[gam_line_list.index[lindex]]
//         * cr48_spec.energy[gam_line_list.index[lindex]]
//         * f48cr(grid_ptr)*tfact;
//   }
//   else if (gam_line_list.type[lindex] == NUCLIDE_V48)
//   {
//
//       emitt_energy = get_rhoinit(grid_ptr->modelgridindex) / nucmass(NUCLIDE_CR48) / 4. / PI
//           * (exp(-single_t/get_meanlife(NUCLIDE_CR48)) - exp(-single_t/get_meanlife(NUCLIDE_V48)))
//           / (get_meanlife(NUCLIDE_CR48) - get_meanlife(NUCLIDE_V48))
//           * v48_spec.probability[gam_line_list.index[lindex]]
//           * v48_spec.energy[gam_line_list.index[lindex]]
//           * f48cr(grid_ptr)*tfact;
//   }
//   else if (gam_line_list.type[lindex] == FAKE_GAM_LINE_ID)
//   {
//     emitt_energy = 0.0;
//   }
//   else
//   {
//     printout("unknown line??\n");
//     abort();
//   }
//
//   /* I'm changing the next bit here (Jan 06) because I think what was
//      here before (below) was wrong in the Doppler terms.
//
//      ray_ptr->e_cmf[nray]+= emitt_energy / fabs(dnuds);
//
//      get_velocity(single_pos, vel_vec, single_t);
//
//      ray_ptr->e_rf[nray] = ray_ptr->e_cmf[nray] / doppler(syn_dir, vel_vec);
//
//   */
//   if (emitt_energy != 0)
//   {
//     double vel_vec[3];
//     get_velocity(single_pos, vel_vec, single_t);
//
//     ray_ptr->e_rf[nray]+= emitt_energy / fabs(dnuds)
//         / doppler(syn_dir, vel_vec) / doppler(syn_dir, vel_vec);
//   }
// }
