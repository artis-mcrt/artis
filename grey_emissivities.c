#include "sn3d.h"
#include "grey_emissivities.h"
#include "grid_init.h"
#include "photo_electric.h"
#include "vectors.h"
#include <string.h>

// private functions
double meanf_sigma(double x);
//int emiss_rlc_load(int nts);


///In this file the call to kappa_rpkt does not fit kappa_rpkts definition any loger!!!
///This MUST BE CHANGED. But it's not only the the call of kappa_rpkt. The dummy packet
///pointer needs more information (e.g. frequency) to calculate kappa_rpkt for the non
///grey case.

/*******************************************************/
int rlc_emiss_gamma(const PKT *pkt_ptr, double dist, double t_current)
{
  /* Subroutine to record the heating rate in a cell due to gamma rays.
By heating rate I mean, for now, really the rate at which the code is making
k-packets in that cell which will then convert into r-packets. This is (going
to be) used for the new light_curve syn-style calculation. */

  /* The intention is that rpkt_emiss will contain the emissivity of r-packets
     in the co-moving frame (which is going to be isotropic). */

  /* This is only done to order v/c for now. */

  /* Called with a packet that is about to travel a
distance dist in the lab frame. Time at start of distance is t_current.*/

  double vel_vec[3];
  PKT dummy;

  //printout("[debug] Execution of rlc_emiss_gamma\n");

  dummy.pos[0] = pkt_ptr->pos[0];
  dummy.pos[1] = pkt_ptr->pos[1];
  dummy.pos[2] = pkt_ptr->pos[2];
  dummy.dir[0] = syn_dir[0];
  dummy.dir[1] = syn_dir[1];
  dummy.dir[2] = syn_dir[2];
  dummy.where = pkt_ptr->where;
  dummy.last_cross = NONE;

  int mgi = cell[pkt_ptr->where].modelgridindex;

  if (dist > 0)
  {
    get_velocity(pkt_ptr->pos, vel_vec, t_current);
    double xx = H * pkt_ptr->nu_cmf / ME / CLIGHT / CLIGHT;

    double heating_cont = ((meanf_sigma(xx)*get_nnetot(mgi)) + sig_photo_electric(pkt_ptr,t_current) + (sig_pair_prod(pkt_ptr, t_current) * (1. - (2.46636e+20 / pkt_ptr->nu_cmf))));
    heating_cont = heating_cont * pkt_ptr->e_rf * dist * (1. - (2.*dot(vel_vec, pkt_ptr->dir)/CLIGHT));

    /* The terms in the above are for Compton, photoelectric and pair production. The pair production  one assumes that a fraction (1. - (1.022 MeV / nu)) of the gamma's energy is thermalised. The remaining 1.022 MeV is made into gamma rays */

    /* For normalisation this needs to be
       1) divided by volume
       2) divided by the length of the time step
       3) divided by 4 pi sr
       This will all be done later
    */

    #ifdef _OPENMP
      #pragma omp atomic
    #endif
    rpkt_emiss[mgi] += 1.e-20 * heating_cont;
  }
  return 0;
}

/*******************************************************/
int rlc_emiss_rpkt(const PKT *pkt_ptr, double dist, double t_current)
{
  /* Subroutine to record the rate of destruction (and re-creation) of
     r-packets by the grey opacity. */

  /* This is only done to order v/c for now. */

  /* Called with a packet that is about to travel a
     distance dist in the lab frame. Time at start of distance is t_current.*/

  double vel_vec[3];
  PKT dummy;

  dummy.pos[0] = pkt_ptr->pos[0];
  dummy.pos[1] = pkt_ptr->pos[1];
  dummy.pos[2] = pkt_ptr->pos[2];
  dummy.dir[0] = syn_dir[0];
  dummy.dir[1] = syn_dir[1];
  dummy.dir[2] = syn_dir[2];
  dummy.where = pkt_ptr->where;
  dummy.last_cross = NONE;

  int mgi = cell[pkt_ptr->where].modelgridindex;

  if (dist > 0.0)
  {
    /* for the weighted estimators version */

    get_velocity(pkt_ptr->pos, vel_vec, t_current);

    double cont = (get_kappagrey(mgi) * get_rho(mgi));
    cont = cont * pkt_ptr->e_rf * dist * (1. - (2.*dot(vel_vec, pkt_ptr->dir)/CLIGHT));

    /* For normalisation this needs to be
       1) divided by volume
       2) divided by the length of the time step
       3) divided by 4 pi sr
       This will all be done later
    */

    //printout("%g %g(2)\n",tautau,cont);
    #ifdef _OPENMP
      #pragma omp atomic
    #endif
    rpkt_emiss[mgi] += 1.e-20 * cont;
  }
  return 0;
}


/***********************************************/
int normalise_grey(int nts)
{
  //for (n=0; n < ngrid; n++)
  double dt = time_step[nts].width;
  double helper = pow(time_step[nts].mid / tmin, 3.0);
  for (int n = 0; n < npts_model; n++)
  {
    int assoc_cells = modelgrid[n].associated_cells;
    double dV = vol_init(&cell[n]) * helper;

    rpkt_emiss[n] = rpkt_emiss[n] * ONEOVER4PI / dV / dt / nprocs / assoc_cells;
  }
  return 0;
}


/*************************************************/
int write_grey(int nts)
{
  FILE *est_file, *dummy;
  char chch;
  char filename[100] = "grey_est_";
  char junk[100];

  if ((dummy = fopen("dummy", "w+")) == NULL)
  {
    printout("Cannot open dummy.\n");
    exit(0);
  }
  fprintf(dummy, "%d", nts);
  fclose(dummy);
  if ((dummy = fopen("dummy", "r")) == NULL)
  {
    printout("Cannot open dummy.\n");
    exit(0);
  }
  int i = 0;
  while ((chch=fgetc(dummy)) != EOF)
  {
    junk[i] = chch;
    i = i+1;
  }
  junk[i] = '\0';
  fclose(dummy);

  strcat(filename, junk);
  strcat(filename, ".out");

  if (file_set)
  {
    if ((est_file = fopen(filename, "r")) == NULL)
    {
      printout("Cannot open grey_est_file.txt.\n");
      exit(0);
    }

    //for (n=0; n < ngrid; n++)
    for (int n = 0; n < npts_model; n++)
    {
      float dum;
      fscanf(est_file, "%g", &dum);
      rpkt_emiss[n] += dum;
    }
    fclose(est_file);
  }

  if ((est_file = fopen(filename, "w+")) == NULL)
  {
    printout("Cannot open grey_est_file.txt.\n");
    exit(0);
  }

  //for (n=0; n < ngrid; n++)
  for (int n = 0; n < npts_model; n++)
  {
    fprintf(est_file, " %g\n ", rpkt_emiss[n]);
  }
  fclose(est_file);
  return 0;
}

/***********************************************/
// Routine to compute the mean energy converted to non-thermal electrons times
// the Klein-Nishina cross section.
double meanf_sigma(double x)
{
  double f = 1+(2*x);

  double term0 = 2/x;
  double term1 = ( 1 - (2/x) - (3/(x*x)) ) * log(f);
  double term2 = ( (4 / x) + (3/(x*x)) - 1) * 2 * x / f;
  double term3 = ( 1 - (2/x) - (1/(x*x))) * 2 * x *(1 + x) / f / f;
  double term4 =  -2. * x * ((4*x*x) + (6*x) + 3) / 3 / f / f / f;

  double tot = 3 * SIGMA_T * (term0 + term1 + term2 + term3 + term4) / (8 * x);

  return tot;
}

/**************************************************************/
/*int emiss_rlc_load(int nts)
{
  // Routine to read in the stored estimators for the time step that is about to begin.
  FILE *est_file, *dummy;
  char chch;
  char filename[100] = "grey_est_";
  char junk[100];
  float dum;

  if ((dummy = fopen("dummy", "w+")) == NULL)
  {
    printout("Cannot open dummy.\n");
    exit(0);
  }
  fprintf(dummy, "%d", nts);
  fclose(dummy);
  if ((dummy = fopen("dummy", "r")) == NULL)
  {
    printout("Cannot open dummy.\n");
    exit(0);
  }
  int i = 0;
  while ((chch=fgetc(dummy)) != EOF)
  {
    junk[i] = chch;
    i= i+1;
  }
  junk[i] = '\0';
  fclose(dummy);

  strcat(filename, junk);
  strcat(filename, ".out");

  if ((est_file = fopen(filename, "r")) == NULL)
  {
    printout("Cannot open est_file.txt.\n");
    exit(0);
  }

  //for (n=0; n < ngrid; n++)
  for (int n = 0; n < npts_model; n++)
  {
    fscanf(est_file, "%g", &dum);
    rpkt_emiss[n] = dum;
  }
  fclose(est_file);
  return 0;
}*/

/***********************************************/
/*int grey_rt(RAY *ray_ptr, int nray, double ldist, double *single_pos, double single_t, int lindex)
{
  // This is called when a grey ray is about to be moved a distance ldist.
  // It should account for the changes in the ray intensity due to
  //   grey processes along the path.
  PKT dummy;
  double kap_tot, tau_cont;
  double vel_vec[3];

  /// Make a dummy packet that carries the ray properties.

  dummy.pos[0] = single_pos[0];
  dummy.pos[1] = single_pos[1];
  dummy.pos[2] = single_pos[2];

  dummy.dir[0] = syn_dir[0];
  dummy.dir[1] = syn_dir[1];
  dummy.dir[2] = syn_dir[2];

  dummy.where = ray_ptr->where;
  dummy.nu_cmf = ray_ptr->nu_cmf[nray];

  get_velocity(single_pos, vel_vec, single_t);

  if (do_rlc_est == 1)
  {

    kap_tot = calculate_kappa_rpkt_cont(&dummy, single_t);
    tau_cont = kap_tot * ldist;
    ray_ptr->e_rf[nray] = ray_ptr->e_rf[nray] * exp(-1. * tau_cont);

  }
  else if (do_rlc_est == 2)
  {
    tau_cont = kap_tot = 0;
  }
  else
  {
    printout("Unknown rlc type. Abort.\n");
    exit(0);
  }

   // Now adding the emissivity term.
  // I think it's a doppler^3 term because it's a doppler^2 term in Castor + there's an
  //   extra doppler from the integration over frequency of the RT equation.

  if (tau_cont > 1.e-6)
  {
    ray_ptr->e_rf[nray] += (rpkt_emiss[cell[dummy.where].modelgridindex] /
      doppler(syn_dir, vel_vec) / doppler(syn_dir, vel_vec) / doppler(syn_dir, vel_vec)
      *(1. - exp(-1. * tau_cont)) / kap_tot);
  }
  else
  {
    ray_ptr->e_rf[nray] += (rpkt_emiss[cell[dummy.where].modelgridindex]   /
        doppler(syn_dir, vel_vec) / doppler(syn_dir, vel_vec) / doppler(syn_dir, vel_vec)
          * ldist);
  }

  //This MUST be followed by a call to move_one_ray() in source
  // since e_cmf is NOT reset here.

  return 0;
}*/
