/* 2007-10-30 -- MK
   Non-grey treatment of UVOIR opacity as opacity_case 4 added.
   Still not fully commented.
   Comments are marked by ///  Deactivated code by // */
/* 2007-01-17 -- MK
   Several minor modifications (some marked in the code with //MK), these include
     - global printout() routine (located in sn3d.c)
     - opacity_cases 2 and 3 added (changes in grid_init.c and update_grid.c,
       original opacity stuff was moved there from input.c) */
/* This is a code copied from Lucy 2004 paper on t-dependent supernova
   explosions. */

#include <stdio.h>
#include <stdbool.h>
#include "exspec.h"
#include "sn3d.h"
#include "input.h"
#include "light_curve.h"
#include "packet_init.h"
#include "spectrum.h"
#include "vectors.h"

int nprocs_exspec;
bool do_emission_res;

// threadprivate variables
FILE *output_file;
int tid;
bool use_cellhist;
bool neutral_flag;
gsl_rng *rng;
gsl_integration_workspace *gslworkspace;


static void get_final_packets(int rank, int nprocs, PKT pkt[])
{
  char filename[100];
  /// Read in the next bunch of packets to work on
  //sprintf(filename,"packets%d_%d.tmp",0,i);

  sprintf(filename, "packets%.2d_%.4d.out", 0, rank);
  printout("reading %s (file %d of %d)\n", filename, rank + 1, nprocs);
  if (access(filename, F_OK) != -1)
  {
    read_packets(filename, pkt);
  }
  else
  {
    printout("   WARNING %s does not exist - trying temp packets file at beginning of timestep %d...\n   ", filename, itstep);
    read_temp_packetsfile(itstep, rank, pkt);
  }
}


int main(int argc, char** argv)
{
  const int my_rank = 0;
  char filename[100];

  sprintf(filename, "exspec.txt");
  output_file = fopen_required(filename, "w");
  setvbuf(output_file, NULL, _IOLBF, 1);

  const time_t sys_time_start = time(NULL);

  printout("Begining do_exspec.\n");

  /// Get input stuff
  printout("time before input %ld\n", time(NULL));
  input(my_rank);
  printout("time after input %ld\n", time(NULL));
  nprocs = nprocs_exspec;

  PKT *pkts = (PKT *) malloc(npkts * sizeof(PKT));

  nnubins = MNUBINS; //1000;  /// frequency bins for spectrum

  init_spectrum_trace(); // needed for TRACE_EMISSION_ABSORPTION_REGION_ON

  struct spec *rpkt_spectra = alloc_spectra(do_emission_res);

  struct spec *stokes_i = NULL;
  struct spec *stokes_q = NULL;
  struct spec *stokes_u = NULL;

  #ifdef POL_ON
  stokes_i = alloc_spectra(do_emission_res);
  stokes_q = alloc_spectra(do_emission_res);
  stokes_u = alloc_spectra(do_emission_res);
  #endif

  struct spec *gamma_spectra = alloc_spectra(false);

  for (int outer_iteration = 0; outer_iteration < n_out_it; outer_iteration++)
  {
    /// Initialise the grid. Call routine that sets up the initial positions
    /// and sizes of the grid cells.
    //grid_init();
    time_init();

    const int amax = ((model_type == RHO_1D_READ)) ? 0 : MABINS;
    // a is the escape direction angle bin
    for (int a = -1; a < amax; a++)
    {
      /// Set up the light curve grid and initialise the bins to zero.
      double *rpkt_light_curve_lum = (double *) calloc(ntstep, sizeof(double));
      double *rpkt_light_curve_lumcmf = (double *) calloc(ntstep, sizeof(double));
      double *gamma_light_curve_lum = (double *) calloc(ntstep, sizeof(double));
      double *gamma_light_curve_lumcmf = (double *) calloc(ntstep, sizeof(double));
      /// Set up the spectrum grid and initialise the bins to zero.

      init_spectra(rpkt_spectra, nu_min_r, nu_max_r, do_emission_res);

      #ifdef POL_ON
      init_spectra(stokes_i, nu_min_r, nu_max_r, do_emission_res);
      init_spectra(stokes_q, nu_min_r, nu_max_r, do_emission_res);
      init_spectra(stokes_u, nu_min_r, nu_max_r, do_emission_res);
      #endif

      const double nu_min_gamma = 0.05 * MEV / H;
      const double nu_max_gamma = 4. * MEV / H;
      init_spectra(gamma_spectra, nu_min_gamma, nu_max_gamma, false);

      for (int p = 0; p < nprocs; p++)
      {
        get_final_packets(p, nprocs, pkts);
        int nesc_tot = 0;
        int nesc_gamma = 0;
        int nesc_rpkt = 0;
        for (int ii = 0; ii < npkts; ii++)
        {
          // printout("packet %d escape_type %d type %d", ii, pkts[ii].escape_type, pkts[ii].type);
          if (pkts[ii].type == TYPE_ESCAPE)
          {
            nesc_tot++;
            if (pkts[ii].escape_type == TYPE_RPKT)
            {
              nesc_rpkt++;
              add_to_lc_res(&pkts[ii], a, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
              add_to_spec_res(&pkts[ii], a, rpkt_spectra, stokes_i, stokes_q, stokes_u);
            }
            else if (pkts[ii].escape_type == TYPE_GAMMA && a == -1)
            {
              nesc_gamma++;
              add_to_lc_res(&pkts[ii], a, gamma_light_curve_lum, gamma_light_curve_lumcmf);
              add_to_spec_res(&pkts[ii], a, gamma_spectra, NULL, NULL, NULL);
            }
          }
        }
        printout("  %d of %d packets escaped (%d gamma-pkts and %d r-pkts)\n", nesc_tot, npkts, nesc_gamma, nesc_rpkt);
      }

      if (a == -1)
      {
        /// Extract angle-averaged spectra and light curves
        write_light_curve((char *) "light_curve.out", -1, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
        write_light_curve((char *) "gamma_light_curve.out", -1, gamma_light_curve_lum, gamma_light_curve_lumcmf);

        write_spectrum((char *) "spec.out", (char *) "emission.out",
                       (char *) "emissiontrue.out", (char *) "absorption.out", rpkt_spectra);
        #ifdef POL_ON
        write_specpol((char *) "specpol.out", (char *) "emissionpol.out", (char *) "absorptionpol.out",
                      stokes_i, stokes_q, stokes_u);
        #endif
        write_spectrum((char *) "gamma_spec.out", NULL, NULL, NULL, gamma_spectra);
      }
      else
      {
        /// Extract LOS dependent spectra and light curves
        char lc_filename[100] = "";
        char spec_filename[100] = "";
        char emission_filename[100] = "";
        char trueemission_filename[100] = "";
        char absorption_filename[100] = "";

        #ifdef POL_ON
        char specpol_filename[100] = "";
        sprintf(specpol_filename, "specpol_res_%.2d.out", a);
        char emissionpol_filename[100] = "";
        char absorptionpol_filename[100] = "";
        #endif

        sprintf(lc_filename, "light_curve_res_%.2d.out", a);
        sprintf(spec_filename, "spec_res_%.2d.out", a);

        if (do_emission_res)
        {
          sprintf(emission_filename, "emission_res_%.2d.out", a);
          sprintf(trueemission_filename, "emissiontrue_res_%.2d.out", a);
          sprintf(absorption_filename, "absorption_res_%.2d.out", a);
          #ifdef POL_ON
          sprintf(emissionpol_filename, "emissionpol_res_%.2d.out", a);
          sprintf(absorptionpol_filename, "absorptionpol_res_%.2d.out", a);
          #endif
        }

        write_light_curve(lc_filename, a, rpkt_light_curve_lum, rpkt_light_curve_lumcmf);
        write_spectrum(spec_filename, emission_filename, trueemission_filename, absorption_filename, rpkt_spectra);

        #ifdef POL_ON
        write_specpol(
          specpol_filename, emissionpol_filename, absorptionpol_filename,
          stokes_i, stokes_q, stokes_u);
        #endif
      }

      if (a == -1)
      {
        printout("finished angle-averaged stuff\n");
      }
      else
      {
        printout("Did %d of %d angle bins.\n", a + 1, MABINS);
      }

      free(rpkt_light_curve_lum);
      free(rpkt_light_curve_lumcmf);

      free(gamma_light_curve_lum);
      free(gamma_light_curve_lumcmf);
    }
  }

  free_spectra(rpkt_spectra);

  if (stokes_i != NULL)
    free_spectra(stokes_i);
  if (stokes_q != NULL)
    free_spectra(stokes_q);
  if (stokes_u != NULL)
    free_spectra(stokes_u);

  free_spectra(gamma_spectra);

  //fclose(ldist_file);
  //fclose(output_file);

  /* Spec syn. */
  //grid_init();
  //syn_gamma();
  free(pkts);

  printout("exspec finished at %ld (tstart + %ld seconds)\n", time(NULL), time(NULL) - sys_time_start);
  fclose(output_file);

  return 0;
}


extern inline void gsl_error_handler_printout(const char *reason, const char *file, int line, int gsl_errno);
extern inline FILE *fopen_required(const char *filename, const char *mode);


// define these functions since sn3d.c is not included when compiling exspec
void increment_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstatscounters ion_counter_type, const double increment)
{
}

double get_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstatscounters ion_counter_type)
{
  return 0.;
}

void set_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstatscounters ion_counter_type, const double newvalue)
{
}

/*void *my_malloc(size_t size)
{
  char *adr;
  #ifdef POWER6
    adr = &heap[heapcounter];
    heapcounter += size;
    if (heapcounter >= heapsize) adr = NULL;
  #else
//    adr = &heap[heapcounter];
//    heapcounter += size;
//    if (heapcounter >= heapsize) adr = NULL;
    adr = malloc(size);
  #endif
  return adr;
}*/
