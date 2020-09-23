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

#include <stdbool.h>
#include "exspec.h"
#include "sn3d.h"
#include "input.h"
#include "light_curve.h"
#include "packet_init.h"
#include "spectrum.h"
#include "vectors.h"

#ifdef DO_EXGAMMA
const bool MODE_GAMMA = true;
#else
const bool MODE_GAMMA = false;
#endif

int nprocs_exspec;
bool do_emission_res;
int nepkts;

double dlognu;


// threadprivate variables
FILE *output_file;
int tid;
bool use_cellhist;
bool neutral_flag;
gsl_rng *rng;
gsl_integration_workspace *gslworkspace;


static int get_escaped_packets(int i, int nprocs, PKT pkt[], EPKT *epkts, int npkts, enum packet_type escape_type)
{
  char filename[100];
  /// Read in the next bunch of packets to work on
  //sprintf(filename,"packets%d_%d.tmp",0,i);

  sprintf(filename, "packets%.2d_%.4d.out", 0, i);
  printout("reading %s (file %d of %d)\n", filename, i + 1, nprocs);
  if (access(filename, F_OK) != -1)
  {
    read_packets(filename, pkt);
  }
  else
  {
    printout("   WARNING %s does not exist - trying temp packets file at beginning of timestep %d...\n   ", filename, i);
    read_temp_packetsfile(itstep, i, pkt);
  }

  int nepkts = 0;
  for (int ii = 0; ii < npkts; ii++)
  {
    // printout("packet %d escape_type %d type %d", ii, pkt[ii].escape_type, pkt[ii].type);
    if (pkt[ii].type == TYPE_ESCAPE && pkt[ii].escape_type == escape_type)
    {
      /// We know that a packet escaped at "escape_time". However, we have
      /// to allow for travel time. Use the formula in Leon's paper. The extra
      /// distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
      const double arrive_time = pkt[ii].escape_time - (dot(pkt[ii].pos, pkt[ii].dir) / CLIGHT_PROP);
      epkts[nepkts].arrive_time = arrive_time;

      /// Now do the cmf time.
      const double arrive_time_cmf = pkt[ii].escape_time * sqrt(1. - (vmax * vmax / CLIGHTSQUARED));
      epkts[nepkts].arrive_time_cmf = arrive_time_cmf;
      // printout("add packet %d. arrive_time %g arrive_time_cmf %g nu_rf %g e_rf %g e_cmf %g \n",
      //          j, arrive_time, arrive_time_cmf, pkt[ii].nu_rf, pkt[ii].e_rf, pkt[ii].nu_cmf);

      epkts[nepkts].dir[0] = pkt[ii].dir[0];
      epkts[nepkts].dir[1] = pkt[ii].dir[1];
      epkts[nepkts].dir[2] = pkt[ii].dir[2];
      epkts[nepkts].nu_rf = pkt[ii].nu_rf;
      epkts[nepkts].e_rf = pkt[ii].e_rf;
      epkts[nepkts].e_cmf = pkt[ii].e_cmf;
      epkts[nepkts].em_pos[0] = pkt[ii].em_pos[0];
      epkts[nepkts].em_pos[1] = pkt[ii].em_pos[1];
      epkts[nepkts].em_pos[2] = pkt[ii].em_pos[2];
      epkts[nepkts].trueem_time = pkt[ii].trueem_time;
      epkts[nepkts].em_time = pkt[ii].em_time;
      epkts[nepkts].emissiontype = pkt[ii].emissiontype;
      epkts[nepkts].trueemissiontype = pkt[ii].trueemissiontype;
      epkts[nepkts].absorptionfreq = pkt[ii].absorptionfreq;
      epkts[nepkts].absorptiontype = pkt[ii].absorptiontype;
      epkts[nepkts].trueemissionvelocity = pkt[ii].trueemissionvelocity;

      nepkts++;
    }
  }
  return nepkts;
}


int main(int argc, char** argv)
{
  const int my_rank = 0;
  char filename[100];

  sprintf(filename, MODE_GAMMA ? "exgamma.txt" : "exspec.txt");
  output_file = fopen_required(filename, "w");
  setvbuf(output_file, NULL, _IOLBF, 1);

  const time_t sys_time_start = time(NULL);

  #ifdef EXGAMMA
  printout("Begining do_exspec for gamma rays.\n");
  #else
  printout("Begining do_exspec.\n");
  #endif

  /// Get input stuff
  printout("time before input %ld\n", time(NULL));
  input(my_rank);
  printout("time after input %ld\n", time(NULL));
  nprocs = nprocs_exspec;

  PKT *pkts = (PKT *) malloc(npkts * sizeof(PKT));
  EPKT *epkts = (EPKT *) malloc(npkts * sizeof(EPKT));

  if (epkts == NULL)
  {
    printout("[fatal] input: not enough memory to initalise escaping packets data structure ... abort\n");
    abort();
  }

  nnubins = MNUBINS; //1000;  /// frequency bins for spectrum
  if (MODE_GAMMA)
  {
    /// Spectra settings
    do_emission_res = 0;         /// We don't record information on gamma packet last interactions, thus create no emission/absorption files.
    nu_min_r = 0.05 * MEV / H;   /// Lower frequency boundary for gamma spectra (badly named variable)
    nu_max_r = 4 * MEV / H;      /// Upper frequency boundary for gamma spectra
  }

  for (int outer_iteration = 0; outer_iteration < n_out_it; outer_iteration++)
  {
    /// Initialise the grid. Call routine that sets up the initial positions
    /// and sizes of the grid cells.
    //grid_init();
    time_init();

    /// Loop over all packets in all the packets files of the simulation and check if
    /// a packet made it out as a rpkt or not. Escaping r-packets are stored in the
    /// epkts array, which is then used for the binning.
    const int amax = ((model_type == RHO_1D_READ) || MODE_GAMMA) ? 0 : MABINS;
    for (int a = -1; a < amax; a++)
    {
      /// Set up the light curve grid and initialise the bins to zero.
      init_light_curve();
      /// Set up the spectrum grid and initialise the bins to zero.
      init_spectrum();


      for (int p = 0; p < nprocs; p++)
      {
        nepkts = get_escaped_packets(p, nprocs, pkts, epkts, npkts, MODE_GAMMA ? TYPE_GAMMA : TYPE_RPKT);
        printout("  %d of %d packets escaped with type %s\n", nepkts, npkts, MODE_GAMMA ? "TYPE_GAMMA" : "TYPE_RPKT");

        if (a == -1)
        {
          printout("gather_light_curve\n");
          gather_light_curve(epkts, nepkts);
          printout("gather_spectrum\n");
          gather_spectrum(epkts, nepkts, -1, do_emission_res);
        }
        else
        {
          gather_light_curve_res(epkts, nepkts, a);
          gather_spectrum_res(epkts, nepkts, a);
        }
      }

      if (a == -1)
      {
        /// Extract angle-averaged spectra and light curves
        if (!MODE_GAMMA)
        {
          write_light_curve((char *) "light_curve.out", -1);
          write_spectrum((char *) "spec.out", do_emission_res, (char *) "emission.out",
                         (char *) "emissiontrue.out", (char *) "absorption.out");
        }
        else
        {
          write_light_curve((char *) "gamma_light_curve.out", -1);
          write_spectrum((char *) "gamma_spec.out", do_emission_res, NULL, NULL, NULL);
        }
      }
      else
      {
        /// Extract LOS dependent spectra and light curves
        char lc_filename[100] = "";
        char spec_filename[100] = "";
        char emission_filename[100] = "";
        char trueemission_filename[100] = "";
        char absorption_filename[100] = "";

        sprintf(lc_filename, "light_curve_res_%.2d.out", a);
        printout("%s \n", lc_filename);

        sprintf(spec_filename, "spec_res_%.2d.out", a);
        printout("%s \n", spec_filename);

        if (do_emission_res)
        {
          sprintf(emission_filename, "emission_res_%.2d.out", a);
          printout("%s \n", emission_filename);

          sprintf(trueemission_filename, "emissiontrue_res_%.2d.out", a);
          printout("%s \n", trueemission_filename);

          sprintf(absorption_filename, "absorption_res_%.2d.out", a);
          printout("%s \n", absorption_filename);
        }

        write_light_curve(lc_filename, a);
        write_spectrum(spec_filename, do_emission_res, emission_filename, trueemission_filename, absorption_filename);
      }
      if (a == -1)
      {
        printout("finished angle-averaged stuff\n");
      }
      else
      {
        printout("Did %d of %d angle bins.\n", a + 1, MABINS);
      }
    }
  }

  //fclose(ldist_file);
  //fclose(output_file);

  /* Spec syn. */
  //grid_init();
  //syn_gamma();
  free(epkts);
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
