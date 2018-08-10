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

#include "exspec.h"
#include "sn3d.h"
#include "threadprivate.h"
#include "input.h"
#include "time_init.h"
#include "light_curve.h"
#include "packet_init.h"
#include "spectrum.h"
#include "vectors.h"

static PKT pkt[MPKTS];


int main(int argc, char** argv)
{
  const int my_rank = 0;
  char filename[100];
  sprintf(filename,"exspec.txt");
  output_file = fopen_required(filename, "w");
  setvbuf(output_file, NULL, _IOLBF, 1);

  const time_t sys_time_start = time(NULL);
  printout("Begining do_exspec.\n");

  /// Get input stuff
  printout("time before input %d\n",time(NULL));
  input(my_rank);
  printout("time after input %d\n",time(NULL));
  nprocs = nprocs_exspec;

  /// Read binary packet files and create ASCII packets files out of them
  /*
  npkts=MPKTS;
  for (i = 0; i < nprocs; i++)
  {
    /// Read in the next bunch of packets to work on
    sprintf(filename,"packets%d_%d.tmp",0,i);
    printout("%s\n",filename);
    packets_file = fopen_required(filename, "rb")) == NULL);
    //sprintf(filename,"packets%.2d_%.4d.out",0,i);
    //packets_file = fopen_required(filename, "r");
    fread(&pkt[0], sizeof(PKT), npkts, packets_file);
    //read_packets(packets_file);
    /// Close the current file.
    fclose(packets_file);


    /// Read in the next bunch of packets to work on
    sprintf(filename,"packets%.2d_%.4d.out",0,i);
    printout("%s\n",filename);
    packets_file = fopen_required(filename, "w");
    write_packets(packets_file);
    //read_packets(packets_file);
    /// Close the current file.
    fclose(packets_file);
  }
  abort();
  */

  for (int outer_iteration = 0; outer_iteration < n_out_it; outer_iteration++)
  {
    /// Initialise the grid. Call routine that sets up the initial positions
    /// and sizes of the grid cells.
    //grid_init();
    time_init();

    if ((epkts = malloc((nprocs * npkts) * sizeof(EPKT))) == NULL)
    {
      printout("[fatal] input: not enough memory to initalise escaping packets data structure ... abort\n");
      abort();
    }

    /// Loop over all packets in all the packets files of the simulation and check if
    /// a packet made it out as a rpkt or not. Escaping r-packets are stored in the
    /// epkts array, which is then used for the binning.
    int j = 0;
    for (int i = 0; i < nprocs; i++)
    {
      /// Read in the next bunch of packets to work on
      //sprintf(filename,"packets%d_%d.tmp",0,i);
      sprintf(filename,"packets%.2d_%.4d.out", 0, i);
      printout("reading %s (file %d of %d)\n", filename, i + 1, nprocs);
      //packets_file = fopen_required(filename, "rb");
      FILE *packets_file = fopen_required(filename, "r");
      //fread(&pkt[0], sizeof(PKT), npkts, packets_file);
      read_packets(packets_file, pkt);
      fclose(packets_file);

      for (int ii = 0; ii < npkts; ii++)
      {
        PKT *pkt_ptr = &pkt[ii];
        // printout("packet %d escape_type %d type %d", ii, pkt_ptr->escape_type, pkt_ptr->type);
        if (pkt_ptr->escape_type == TYPE_RPKT && pkt_ptr->type == TYPE_ESCAPE)
        {
          /// We know that a packet escaped at "escape_time". However, we have
          /// to allow for travel time. Use the formula in Leon's paper. The extra
          /// distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
          const double arrive_time = pkt_ptr->escape_time - (dot(pkt_ptr->pos, pkt_ptr->dir) / CLIGHT_PROP);
          epkts[j].arrive_time = arrive_time;

          /// Now do the cmf time.
          const double arrive_time_cmf = pkt_ptr->escape_time * sqrt(1. - (vmax * vmax / CLIGHTSQUARED));
          epkts[j].arrive_time_cmf = arrive_time_cmf;
          // printout("add packet %d. arrive_time %g arrive_time_cmf %g nu_rf %g e_rf %g e_cmf %g \n",
          //          j, arrive_time, arrive_time_cmf, pkt_ptr->nu_rf, pkt_ptr->e_rf, pkt_ptr->nu_cmf);

          epkts[j].dir[0] = pkt_ptr->dir[0];
          epkts[j].dir[1] = pkt_ptr->dir[1];
          epkts[j].dir[2] = pkt_ptr->dir[2];
          epkts[j].nu_rf = pkt_ptr->nu_rf;
          epkts[j].e_rf = pkt_ptr->e_rf;
          epkts[j].e_cmf = pkt_ptr->e_cmf;
          epkts[j].em_pos[0] = pkt_ptr->em_pos[0];
          epkts[j].em_pos[1] = pkt_ptr->em_pos[1];
          epkts[j].em_pos[2] = pkt_ptr->em_pos[2];
          epkts[j].em_time = pkt_ptr->em_time;
          epkts[j].emissiontype = pkt_ptr->emissiontype;
          epkts[j].trueemissiontype = pkt_ptr->trueemissiontype;
          epkts[j].absorptionfreq = pkt_ptr->absorptionfreq;
          epkts[j].absorptiontype = pkt_ptr->absorptiontype;
          epkts[j].trueemissionvelocity = pkt_ptr->trueemissionvelocity;
          j++;
        }
      }
    }
    nepkts = j;


    /// Extract angle-averaged spectra and light curves
    FILE *lc_file = fopen_required("light_curve.out", "w");
    FILE *spec_file = fopen_required("spec.out", "w");
    FILE *emission_file = fopen_required("emission.out", "w");
    FILE *trueemission_file = fopen_required("emissiontrue.out", "w");
    FILE *absorption_file = fopen_required("absorption.out", "w");

    gather_spectrum(-1);
    gather_light_curve();
    write_spectrum(spec_file, emission_file, trueemission_file, absorption_file);
    write_light_curve(lc_file, -1);
    //make_gamma_light_curve();

    fclose(lc_file);
    fclose(spec_file);
    fclose(emission_file);
    fclose(trueemission_file);
    fclose(absorption_file);

    printout("finished angle-averaged stuff\n");

    /// Extract LOS dependent spectra and light curves
    if (model_type != RHO_1D_READ)
    {
      lc_file = fopen_required("light_curve_res.out", "w");
      spec_file = fopen_required("spec_res.out", "w");
      for (int i = 0; i < MABINS; i++)
      {
        if (do_emission_res == 1)
        {
          sprintf(filename, "emission_res_%.2d.out", i);
          printout("%s \n", filename);
          emission_file = fopen_required(filename, "w");

          sprintf(filename, "emissiontrue_res_%.2d.out", i);
          printout("%s \n", filename);
          trueemission_file = fopen_required(filename, "w");

          sprintf(filename, "absorption_res_%.2d.out", i);
          printout("%s \n", filename);
          absorption_file = fopen_required(filename, "w");
        }
        gather_spectrum_res(i);
        gather_light_curve_res(i);

        write_spectrum(spec_file, emission_file, trueemission_file, absorption_file);
        write_light_curve(lc_file, i);

        if (do_emission_res == 1)
        {
          fclose(emission_file);
          fclose(trueemission_file);
          fclose(absorption_file);
        }
        printout("Did %d of %d angle bins.\n",i+1,MABINS);
      }
      fclose(lc_file);
      fclose(spec_file);
    }
  }

  //fclose(ldist_file);
  //fclose(output_file);

  /* Spec syn. */
  //grid_init();
  //syn_gamma();

  printout("exspec finished at %d (tstart + %d seconds)\n", time(NULL), time(NULL) - sys_time_start);
  fclose(output_file);

  return 0;
}


extern inline int printout(const char *restrict format, ...);
extern inline void gsl_error_handler_printout(const char *reason, const char *file, int line, int gsl_errno);
extern inline FILE *fopen_required(const char *filename, const char *mode);


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
