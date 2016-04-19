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


#include "threadprivate.h"
#include "exspec.h"
#include "sn3d.h"
#include <stdarg.h>  /// MK: needed for printout()

/* Main - top level routine. */
int main(int argc, char** argv)
{
  int nts;

  FILE *emission_file,*lc_file,*spec_file,*absorption_file;
  int j,t_arrive;
  PKT *pkt_ptr;

  double syn_gamma();
  int outer_iteration;
  //int gather_spectrum(), write_spectrum(), gather_light_curve(), write_light_curve();
  //int gather_spectrum_res(), write_spectrum_res(), gather_light_curve_res(), write_light_curve_res();
  //int gather_gamma_light_curve(), write_gamma_light_curve();
  int i,ii,iii,interactions;
  double meaninteractions;
  FILE *syn_file;
  FILE *linestat_file;
  FILE *packets_file;
  FILE *temperature_file;
  int middle_iteration;
  int my_rank;
  int p;
  int nnn, nn, n;
  int element;
  #ifdef MPI_ON
    double a, b;
    int aa, bb;
    int nblock, numtot, n_leftover;
  #endif
  int nstart, nknown, ndo;
  double rho_tot, te_tot, tr_tot, w_tot, n_count;
  int position, nlp, ncl, nncl, mgi;
  double T_R_max,T_R_min,T_R_step;
  double T_e_max,T_e_min,T_e_step;
  double rho_max,rho_min,rho_step;
  char filename[100];
  int HUGEE2;
  char *buffer2;
  double nntot;
  int titer;

  double deltaV,deltat;
  int assoc_cells;

//  int HUGEE;

  #ifdef MPI_ON
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
  #else
    my_rank = 0;
    p=1;
  #endif

  nprocs = p;              /// Global variable which holds the number of MPI processes
  rank_global = my_rank;   /// Global variable which holds the rank of the active MPI process
  if (my_rank == 0)
  {


    tid = 0;
    nthreads = 1;
    sprintf(filename,"exspec_dd_%d-%d.txt",my_rank,tid);
    if ((output_file = fopen(filename, "w")) == NULL)
    {
      printf("Cannot open %s.\n",filename);
      abort();
    }
    setvbuf(output_file, NULL, _IOLBF, 1);

    printout("Begining do_exspec_dd.\n");

    /// Get input stuff
    printout("time before input %d\n",time(NULL));
    input (0);
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
      if ((packets_file = fopen(filename, "rb")) == NULL)
      //sprintf(filename,"packets%.2d_%.4d.out",0,i);
      //if ((packets_file = fopen(filename, "r")) == NULL)
      {
        printf("Cannot open packets file %s\n",filename);
        abort();
      }
      fread(&pkt[0], sizeof(PKT), npkts, packets_file);
      //read_packets(packets_file);
      /// Close the current file.
      fclose(packets_file);


      /// Read in the next bunch of packets to work on
      sprintf(filename,"packets%.2d_%.4d.out",0,i);
      printout("%s\n",filename);
      if ((packets_file = fopen(filename, "w")) == NULL)
      {
        printf("Cannot open packets file %s\n",filename);
        abort();
      }
      write_packets(packets_file);
      //read_packets(packets_file);
      /// Close the current file.
      fclose(packets_file);
    }
    exit(0);
    */

    for (outer_iteration = 0; outer_iteration < n_out_it; outer_iteration++)
    {
      /// Initialise the grid. Call routine that sets up the initial positions
      /// and sizes of the grid cells.
      //grid_init();
      time_init();

      if ((epkts = (EPKT*) malloc((nprocs*npkts)*sizeof(EPKT))) == NULL)
      {
        printout("[fatal] input: not enough memory to initalise escaping packets data structure ... abort\n");
        exit(0);
      }

      /// Loop over all packets in all the packets files of the simulation and check if
      /// a packet made it out as a rpkt or not. Escaping r-packets are stored in the
      /// epkts array, which is then used for the binning.
      j=0;
      for (i = 0; i < nprocs; i++)
      {
        /// Read in the next bunch of packets to work on
        //sprintf(filename,"packets%d_%d.tmp",0,i);
        sprintf(filename,"packets%.2d_%.4d.out",0,i);
        printout("%s, %d %d\n",filename,i,nprocs);
        //if ((packets_file = fopen(filename, "rb")) == NULL)
        if ((packets_file = fopen(filename, "r")) == NULL)
        {
          printf("Cannot open packets file %s\n",filename);
          abort();
        }
        //fread(&pkt[0], sizeof(PKT), npkts, packets_file);
        read_packets(packets_file);
        fclose(packets_file);

        for (ii = 0; ii < npkts; ii++)
        {
          pkt_ptr = &pkt[ii];
          if (pkt_ptr->type == TYPE_ESCAPE && pkt_ptr->escape_type == TYPE_RPKT)
          {
            //printout("add packet %d\n",j);
            /// We know that a packet escaped at "escape_time". However, we have
            /// to allow for travel time. Use the formula in Leon's paper. The extra
            /// distance to be travelled beyond the reference surface is ds = r_ref (1 - mu).
            t_arrive = pkt_ptr->escape_time - (dot(pkt_ptr->pos, pkt_ptr->dir)/CLIGHT_PROP);
            epkts[j].arrive_time = t_arrive;

            /// Now do the cmf time.
            t_arrive = pkt_ptr->escape_time * sqrt(1. - (vmax*vmax/CLIGHT2));
            epkts[j].arrive_time_cmf = t_arrive;

            epkts[j].dir[0] = pkt_ptr->dir[0];
            epkts[j].dir[1] = pkt_ptr->dir[1];
            epkts[j].dir[2] = pkt_ptr->dir[2];
            epkts[j].nu_rf = pkt_ptr->nu_rf;
            epkts[j].e_rf = pkt_ptr->e_rf;
            epkts[j].e_cmf = pkt_ptr->e_cmf;
            epkts[j].emissiontype = pkt_ptr->emissiontype;
            epkts[j].em_pos[0] = pkt_ptr->em_pos[0];
            epkts[j].em_pos[1] = pkt_ptr->em_pos[1];
            epkts[j].em_pos[2] = pkt_ptr->em_pos[2];
            epkts[j].em_time = pkt_ptr->em_time;
            epkts[j].absorptionfreq = pkt_ptr->absorptionfreq;
            epkts[j].absorptiontype = pkt_ptr->absorptiontype;
            j += 1;
          }
        }
      }
      nepkts = j;


      /// Extract angle-averaged spectra and light curves
      /// Start with light curves
      if ((lc_file = fopen("light_curve.out", "w")) == NULL)
      {
        printout("Cannot open light_curve.out\n");
        exit(0);
      }
      gather_light_curve();
      write_light_curve(lc_file,-1);
      //make_gamma_light_curve();
      fclose(lc_file);

      /// Now do depth-dependent spectra
      for (i = 0; i < 10; i++)
      {
        sprintf(filename,"spec%d.out",i);
        printout("%s \n",filename);

        if ((spec_file = fopen(filename, "w")) == NULL)
        {
          printout("Cannot open spec.out\n");
          exit(0);
        }

        sprintf(filename,"emission%d.out",i);
        printout("%s \n",filename);
        if ((emission_file = fopen(filename, "w")) == NULL)
        {
          printf("Cannot open emission_res.out\n");
          exit(0);
        }

        sprintf(filename,"absorption%d.out",i);
        printout("%s \n",filename);
        if ((absorption_file = fopen(filename, "w")) == NULL)
        {
          printf("Cannot open absorption_res.out\n");
          exit(0);
        }

        gather_spectrum(i);
        write_spectrum(spec_file,emission_file,absorption_file);

        fclose(spec_file);
        fclose(emission_file);
        fclose(absorption_file);
      }


      printout("finished angle-averaged stuff\n");

      /// Extract LOS dependent spectra and light curves
      if (model_type != RHO_1D_READ)
      {
        if ((lc_file = fopen("light_curve_res.out", "w")) == NULL)
        {
          printout("Cannot open light_curve_res.out\n");
          exit(0);
        }
        if ((spec_file = fopen("spec_res.out", "w")) == NULL)
        {
          printout("Cannot open spec_res.out\n");
          exit(0);
        }
        for (i = 0; i < MABINS; i++)
        {
          if (do_emission_res == 1)
          {
            sprintf(filename,"emission_res_%.2d.out",i);
            printout("%s \n",filename);
            if ((emission_file = fopen(filename, "w")) == NULL)
            {
              printf("Cannot open emission_res.out\n");
              exit(0);
            }

            sprintf(filename,"absorption_res_%.2d.out",i);
            printout("%s \n",filename);
            if ((absorption_file = fopen(filename, "w")) == NULL)
            {
              printf("Cannot open absorption_res.out\n");
              exit(0);
            }
          }
          gather_spectrum_res(i);
          gather_light_curve_res(i);

          write_spectrum(spec_file,emission_file,absorption_file);
          write_light_curve(lc_file,i);

          if (do_emission_res == 1)
          {
            fclose(emission_file);
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

    printout("simulation finished at %d\n",time(NULL));
    fclose(output_file);
  }

  #ifdef MPI_ON
    MPI_Finalize();
  #endif

  return 0;
}


// printout should be used instead of printf throughout the whole code for output messages
int printout(const char *restrict format, ...)
{
   int ret_status = 0;

   va_list args;
   va_start(args,format);
   ret_status = vfprintf(output_file,format,args);
   va_end(args);
   return ret_status;
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
