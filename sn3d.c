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
#include "sn3d.h"
#include "vpkt.h"
#include "version.h"
#include <stdarg.h>  /// MK: needed for printout()

/* Main - top level routine. */
int main(int argc, char** argv)
{
  void initialise_photoionestimators();
  void init_vspecpol();
  void init_vpkt_grid();
  int write_vspecpol(FILE *specpol_file);
  int write_vpkt_grid(FILE *vpkt_grid_file);
  int read_vspecpol(FILE *specpol_file);
  int read_vpkt_grid(FILE *vpkt_grid_file);
  FILE *vspecpol_file;
  FILE *vpkt_grid_file;
//  void precalculate_partfuncts(int cellnumber);
  void tabulate_bb_rad();
  void tabulate_ratecoefficients_gsl();
  void read_packets(FILE *packets_file);
  void write_packets(FILE *packets_file);
  void print_opticaldepth(int cellnumber, int timestep, int samplecell, int element);
  void write_grid_restart_data(void);
  void update_parameterfile(int nts);

  int get_element(int element);
  int get_ionstage(int element, int ion);
  int time_init();
  int grid_init();
  int packet_init(int middle_iteration, int my_rank);
  int input();
  int nts;
  int update_grid(int m, int my_rank, int nstart, int nblock, int titer);
  int update_packets();
  int make_spectrum_res();
  int make_spectrum();
  int make_light_curve();
  int make_light_curve_res();
  int make_gamma_light_curve();
  double syn_gamma();
  int estim_switch();
  int outer_iteration;
  int zero_estimators(), normalise_estimators(), write_estimators();
  int normalise_grey(), write_grey();
  int gather_spectrum(), write_spectrum(), gather_light_curve(), write_light_curve();
  int gather_spectrum_res(), write_spectrum_res(), gather_light_curve_res(), write_light_curve_res();
  int gather_gamma_light_curve(), write_gamma_light_curve();
  double calculate_populations(int cellnumber, int first_nonempty_cell);
  //void determine_important_bfcontinua(int cellnumber, int timestep, int samplecell);
  int i,ii,iii,interactions;
  double meaninteractions;
  FILE *syn_file;
  FILE *linestat_file;
  FILE *packets_file;
  FILE *temperature_file;
  FILE *dep_file;
  double depvalue;
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
  int last_loop;

  double deltaV,deltat;
  int assoc_cells;
  int real_time_start, do_this_full_loop;


  nvpkt=0;
  nvpkt_esc1=0;
  nvpkt_esc2=0;
  nvpkt_esc3=0;



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




  #ifdef _OPENMP
    /// Explicitly turn off dynamic threads becaucse we use the threadprivate directive!!!
    omp_set_dynamic(0);
    #pragma omp parallel private(filename)
    {
      /// Get the current threads ID, copy it to a threadprivate variable
      tid = omp_get_thread_num();
      /// and initialise the threads outputfile
      sprintf(filename,"output_%d-%d.txt",my_rank,tid);
      if ((output_file = fopen(filename, "w")) == NULL)
      {
        printf("Cannot open %s.\n",filename);
        abort();
      }
      /// Makes sure that the output_file is written line-by-line
      setvbuf(output_file, NULL, _IOLBF, 1);

      /// Get the total number of active threads
      nthreads = omp_get_num_threads();
      if (nthreads > MTHREADS)
      {
        printout("[Fatal] too many threads. Set MTHREADS (%d) > nthreads (%d). Abort.\n",MTHREADS,nthreads);
        exit(0);
      }
      printout("OpenMP parallelisation active with %d threads\n",nthreads);
    }
  #else
    tid = 0;
    nthreads = 1;
    sprintf(filename,"output_%d-%d.txt",my_rank,tid);
    if ((output_file = fopen(filename, "w")) == NULL)
    {
      printf("Cannot open %s.\n",filename);
      abort();
    }
    setvbuf(output_file, NULL, _IOLBF, 1);
  #endif


  if ((mastate = (mastate_t *) malloc(nthreads*sizeof(mastate_t))) == NULL)
  {
    printout("[fatal] input: error initializing macro atom state variables ... abort\n");
    exit(0);
  }
  if ((kappa_rpkt_cont = (rpkt_cont_opacity_struct *) malloc(nthreads*sizeof(rpkt_cont_opacity_struct))) == NULL)
  {
    printout("[fatal] input: error initializing continuum opacity communication variables ... abort\n");
    exit(0);
  }
  if ((coolingrates = (coolingrates_t *) malloc(nthreads*sizeof(coolingrates_t))) == NULL)
  {
    printout("[fatal] input: error initializing coolingrates communication variables ... abort\n");
    exit(0);
  }
  if ((heatingrates = (heatingrates_t *) malloc(nthreads*sizeof(heatingrates_t))) == NULL)
  {
    printout("[fatal] input: error initializing heatingrates communication variables ... abort\n");
    exit(0);
  }




  /*
  float *test,*test1;
  printout("Hello!\n");
  printout("sizeof float %d\n",sizeof(float));
  printout("sizeof double %d\n",sizeof(double));
  //printout("allocate 100 floats\n");
  //if ((test = calloc(100, sizeof(float))) == NULL)
  //{
  //  printout("[fatal] input: not enough memory to initialise 100 floats ... abort\n");
  //  exit(0);
  //}

  //printout("allocate 2500MB of floats\n");
  //if ((test1 = malloc(655360000*sizeof(float))) == NULL)
  //{
  //  printout("[fatal] input: not enough memory to initialise 2500MB of floats ... abort\n");
  //}

  test1 = NULL;
  printout("allocate 400MB of floats step by step\n");
  for (i=0; i < 4*26214400; i++)
  {
    if ((test = my_malloc(sizeof(float))) == NULL)
    {
      printout("[fatal] stopped initialisation at %d\n",i);
    }
    if (i % 1000 == 0)
    {
      fprintf(output_file,"i %d, diff %d, newptr %p, oldptr %p\n",i,test-test1,test,test1);
    }
    test1 = test;
  }

  printout("test finished\n");
  MPI_Finalize();
  exit(0);
  */


  /// Using this and the global variable output_file opens and closes the output_file
  /// only once, which speeds up the simulation with a lots of output switched on (debugging).
  /// The downside is that while the simulation runs, its output is only readable on that
  /// machine where the simulation is running.
  /// NB: printout also needs some changes to get this working!
  /*
  if ((output_file = fopen("output.txt", "w")) == NULL)
  {
  printf("Cannot open output.txt.\n");
  exit(0);
}
  /// Makes sure that the output_file is written line-by-line
  setvbuf(output_file, NULL, _IOLBF, 1);
  */

  /*
  if ((ldist_file = fopen("ldist.out", "w")) == NULL){
  printout("Cannot open ldist.out.\n");
  exit(0);
}
  */


  //sprintf(filename,"tb%.4d.txt",my_rank);
  //if ((tb_file = fopen(filename, "w")) == NULL)
  //{
  //  printf("Cannot open %s.\n",filename);
  //  abort();
  //}
  //setvbuf(tb_file, NULL, _IOLBF, 1);

  sprintf(filename,"estimators_%.4d.out",my_rank);
  if ((estimators_file = fopen(filename, "w")) == NULL)
  {
    printout("Cannot open %s.\n",filename);
    exit(0);
  }
  setvbuf(estimators_file, NULL, _IOLBF, 1);

  printout("Start ARTIS revision %s\n",GIT_HASH);
  printout("This binary was compiled on %s\n",COMPILETIME);
  //printout("CELLHISTORYSIZE %d\n",CELLHISTORYSIZE);


  /// Get input stuff
  real_time_start = time(NULL);
  printout("time before input %d\n",time(NULL));
  input (my_rank);


  /// Initialise linestat file
  if (my_rank == 0)
  {
    if ((linestat_file = fopen("linestat.out", "w")) == NULL){
      printout("Cannot open line_stat.out.\n");
      exit(0);
    }
    setvbuf(linestat_file, NULL, _IOLBF, 1);
    for (i = 0; i < nlines; i++) fprintf(linestat_file,"%g ", CLIGHT/linelist[i].nu);
    fprintf(linestat_file,"\n");
    for (i = 0; i < nlines; i++) fprintf(linestat_file,"%d ", get_element(linelist[i].elementindex));
    fprintf(linestat_file,"\n");
    for (i = 0; i < nlines; i++) fprintf(linestat_file,"%d ", get_ionstage(linelist[i].elementindex,linelist[i].ionindex));
    fprintf(linestat_file,"\n");
    for (i = 0; i < nlines; i++) fprintf(linestat_file,"%d ", linelist[i].upperlevelindex+1);
    fprintf(linestat_file,"\n");
    for (i = 0; i < nlines; i++) fprintf(linestat_file,"%d ", linelist[i].lowerlevelindex+1);
    fprintf(linestat_file,"\n");
  }

  printout("time after input %d\n",time(NULL));
  printout("simulation propagates %d packets through a %d x %d x %d grid\n",npkts,nxgrid,nygrid,nzgrid);
  printout("timesteps %d\n",ntstep);



  /// Precalculate the rate coefficients for spontaneous and stimulated recombination
  /// and for photoionisation. With the nebular approximation they only depend on T_e
  /// T_R and W. W is easily factored out. For stimulated recombination we must assume
  /// T_e = T_R for this precalculation.
  /// Make this parallel ?
  printout("time before tabulation of rate coefficients %d\n",time(NULL));
  tabulate_ratecoefficients_gsl();
  printout("time after tabulation of rate coefficients %d\n",time(NULL));
  //tabulate_bb_rad();
  printout("time after tabulation of radfield helpers %d\n",time(NULL));
  //abort();

  /// As a precaution, explicitly zero all the estimators here
  zero_estimators();
  printout("time after zero estimators %d\n",time(NULL));


  /// Record the chosen syn_dir
  if ((syn_file = fopen("syn_dir.txt", "w")) == NULL){
    printout("Cannot open syn_dir.txt.\n");
    exit(0);
  }
  fprintf(syn_file, "%g %g %g", syn_dir[0], syn_dir[1], syn_dir[2]);
  fclose(syn_file);
  printout("time read syn file %d\n",time(NULL));

  file_set=0;
  debuglevel = 4;  /// Selects detail level of debug output, needs still some work.
  for (outer_iteration = 0; outer_iteration < n_out_it; outer_iteration++)
  {
    /// Initialise the grid. Call routine that sets up the initial positions
    /// and sizes of the grid cells.
    time_init();
    printout("time time init %d\n",time(NULL));
    grid_init();

    /// Next we want to initialise the packets.
    /// To overcome memory limitations for large numbers of packets, which need to be
    /// propagated on the _same_ grid, this middle_iteration loop was introduced.
    for(middle_iteration = 0; middle_iteration < n_middle_it; middle_iteration++)
    {
      /// Create a bunch of npkts packets
      /// and write them to a binary file for later readin.
      packet_init(middle_iteration,my_rank);
    }

    /// For the parallelisation of update_grid, the process needs to be told which cells belong to it.
    /// The next loop is over all grid cells. For parallelisation, we want to split this loop between
    /// processes. This is done by assigning each MPI process nblock processes. The residual n_leftover
    /// cells are sent to processes 0 ... process n_leftover -1.
    #ifdef MPI_ON
      //nblock = ngrid / p;
      //nblock = nnonemptycells / p;
      nblock = npts_model / p;
      numtot = nblock * p;
      //if (numtot > ngrid)
      //if (numtot > nnonemptycells)
      if (numtot > npts_model)
      {
        nblock = nblock - 1;
        //n_leftover = ngrid - (nblock * p);
        //n_leftover = nnonemptycells - (nblock * p);
        n_leftover = npts_model - (nblock * p);
      }
      //else if (numtot < ngrid)
      //else if (numtot < nnonemptycells)
      else if (numtot < npts_model)
      {
        //n_leftover = ngrid - (nblock * p);
        //n_leftover = nnonemptycells - (nblock * p);
        n_leftover = npts_model - (nblock * p);
      }
      else
      {
        n_leftover = 0;
      }

      if (my_rank < n_leftover)
      {
        ndo = nblock + 1;
        nstart = my_rank * (nblock + 1);
      }
      else
      {
        ndo = nblock;
        nstart = n_leftover * (nblock + 1) + (my_rank - n_leftover)*(nblock);
      }

      printout("process %d doing %d cells from %d to %d\n",my_rank,ndo,nstart,nstart+ndo-1);

      /// Initialise the exchange buffer
      /// The factor 4 comes from the fact that our buffer should contain elements of 4 byte
      /// instead of 1 byte chars. But the MPI routines don't care about the buffers datatype
      //int HUGEE = 4 * ((9+2*includedions)*(nblock+1) + 1);
      //int HUGEE = 4 * ((8+2*includedions)*(nblock+1) + 1);
      int HUGEE = 4 * ((10+2*includedions)*(nblock+1) + 1) + 8 * ((1+includedions)*(nblock+1));
      printout("reserve HUGEE %d space for MPI communication buffer\n",HUGEE);
      //char buffer[HUGEE];
      char *buffer;
      if ((buffer = (char *) malloc(HUGEE*sizeof(char))) == NULL)
      {
        printout("[fatal] input: not enough memory to initialize MPI exchange buffer ... abort.\n");
        exit(0);
      }
    #else
      nstart = 0;
      //ndo = ngrid;
      //ndo = nnonemptycells;
      ndo = npts_model;
    #endif


    /** That's the end of the initialisation. */
    /** *******************************************/

    /// Now comes the loop over timesteps. Before doing this we need to
    /// define the time steps.
    //time_init();

    /// Standard approach: for loop over given number of time steps
    //for (nts = itstep; nts < ftstep; nts++)
    //{

    /// Now use while loop to allow for timed restarts
    last_loop = ftstep;
    do_this_full_loop = 1;
    nts = itstep;

    // Initialise virtual packets file and vspecpol
    #ifdef VPKT_ON
      sprintf(filename,"vspecpol_%d-%d.out",my_rank,tid);
      if ((vspecpol_file = fopen(filename, "w")) == NULL)
      {
          printout("Cannot open %s.\n",filename);
          exit(0);
      }

      if (vgrid_flag==1)
      {
          sprintf(filename,"vpkt_grid_%d-%d.out",my_rank,tid);
          if ((vpkt_grid_file = fopen(filename, "w")) == NULL)
          {
              printout("Cannot open %s.\n",filename);
              exit(0);
          }

      }

      // New simulation
      if (continue_simulation==0)
      {
          init_vspecpol();

          if (vgrid_flag == 1)
            init_vpkt_grid();
      }
      // Continue simulation: read into temporary files
      else
      {
          if (nts % 2 == 0) sprintf(filename,"vspecpol_%d_%d_odd.tmp",my_rank,tid);
          else sprintf(filename,"vspecpol_%d_%d_even.tmp",my_rank,tid);
          if ((packets_file = fopen(filename, "rb")) == NULL)
          {
              printout("Cannot read temporary packets file %s\n",filename);
              exit(0);
          }

          read_vspecpol(packets_file);


          if (vgrid_flag==1)
          {

              if (nts % 2 == 0) sprintf(filename,"vpkt_grid_%d_%d_odd.tmp",my_rank,tid);
              else sprintf(filename,"vpkt_grid_%d_%d_even.tmp",my_rank,tid);
              if ((packets_file = fopen(filename, "rb")) == NULL)
              {
                  printout("Cannot read temporary vpkt_grid file %s\n",filename);
                  exit(0);
              }

              read_vpkt_grid(packets_file);
          }
      }
    #endif

    while (nts < last_loop)
    {
      #ifdef MPI_ON
        MPI_Barrier(MPI_COMM_WORLD);
      #endif

      #ifdef TIMED_RESTARTS
        if ((time(NULL) - real_time_start) > 9000)
        {
          do_this_full_loop = 0; //This flag will make it do a write out then quit, hopefully
          printout("Going to terminate since remaining time is too short. %d\n", time(NULL) - real_time_start);
        }
        else
        {
          printout("Going to continue. Total time spent so far: %d.\n", time(NULL) - real_time_start);
        }
        #ifdef MPI_ON
          MPI_Bcast(&do_this_full_loop, 1, MPI_INT, 0, MPI_COMM_WORLD);
        #endif
      #endif

      /// The first time step must solve the ionisation balance in LTE
      if (nts == 0) initial_iteration = 1;
      else initial_iteration = 0;

      #ifndef DO_TITER
        /// Do 3 iterations on timestep 0-9
        /*if (nts == 0)
        {
          n_titer = 3;
          initial_iteration = 1;
        }
        else if (nts < 6)
        {
          n_titer = 3;
          initial_iteration = 0;
        }
        else
        {
          n_titer = 1;
          initial_iteration = 0;
        }*/
        if (nts < n_lte_timesteps)
        {
          n_titer = 1;
          initial_iteration = 1;
        }
        else
        {
          n_titer = 1;
          initial_iteration = 0;
        }
      #endif

      for (titer = 0; titer < n_titer; titer++)
      {


        /// Read the packets file for each iteration on the timestep
        if (nts % 2 == 0) sprintf(filename,"packets%d_%d_odd.tmp",tid,my_rank);
        else sprintf(filename,"packets%d_%d_even.tmp",tid,my_rank);
        //sprintf(filename,"packets%d_%d.tmp",tid,my_rank);
        if ((packets_file = fopen(filename, "rb")) == NULL)
        {
          printout("Cannot read temporary packets file %s\n",filename);
          exit(0);
        }
        fread(&pkt[0], sizeof(PKT), npkts, packets_file);
        //read_packets(packets_file);
        fclose(packets_file);

        /// Some counters on pkt-actions need to be reset to do statistics
        ma_stat_activation_collexc = 0;
        ma_stat_activation_collion = 0;
        ma_stat_activation_bb = 0;
        ma_stat_activation_bf = 0;
        ma_stat_deactivation_colldeexc = 0;
        ma_stat_deactivation_collrecomb = 0;
        ma_stat_deactivation_bb = 0;
        ma_stat_deactivation_fb = 0;
        k_stat_to_ma_collexc = 0;
        k_stat_to_ma_collion = 0;
        k_stat_to_r_ff = 0;
        k_stat_to_r_fb = 0;
        k_stat_to_r_bb = 0;
        k_stat_from_ff = 0;
        k_stat_from_bf = 0;
        k_stat_from_gamma = 0;
        k_stat_from_eminus = 0;
        k_stat_from_earlierdecay = 0;
        escounter = 0;
        cellcrossings = 0;
        updatecellcounter = 0;
        coolingratecalccounter = 0;
        resonancescatterings = 0;
        upscatter = 0;
        downscatter = 0;
        if (nts == 0) initialise_photoionestimators();
        //if (nts > 0) debuglevel = 2000;

        #ifdef RECORD_LINESTAT
          /// The same for absorption/emission of r-pkts in lines
          for (i = 0; i < nlines; i++)
          {
            acounter[i] = 0;
            ecounter[i] = 0;
          }
        #endif

        if (do_r_lc == 0)
        {
          do_comp_est = estim_switch(nts);
        }
        else
        {
          do_comp_est = 0;
        }

        nesc = 0;

        /// Update the matter quantities in the grid for the new timestep. */
        printout("time before update grid %d\n",time(NULL));

        #ifndef FORCE_LTE
          /// Initialise corrphotoionrenorm[i] to zero before update_grid is called
          /// This allows reduction after update_grid has finished
          if (continue_simulation && nts-itstep == 0 && titer == 0)
          {
            /// In this case they have been read from file and must neither be touched
            /// nor broadcasted after update_grid
          }
          else
          {
            printout("nts %d, titer %d: reset corr photoionrenorm\n",nts,titer);
            for (i = 0; i < MMODELGRID*nelements*maxion; i++)
            {
              corrphotoionrenorm[i] = 0.;
            }
            printout("after nts %d, titer %d: reset corr photoionrenorm\n",nts,titer);
          }
        #endif

        update_grid(nts,my_rank,nstart,ndo,titer);
        #ifdef DO_TITER
          /// No iterations over the zeroth timestep, set titer > n_titer
          if (nts==0) titer = n_titer+1;
        #endif
        #ifdef MPI_ON
          MPI_Barrier(MPI_COMM_WORLD);
        #endif
        printout("time after update grid %d\n",time(NULL));
        //printout("histindex %d\n",histindex);


        /// Each process has now updated its own set of cells. The results now need to be communicated between processes.
        #ifdef MPI_ON
          for (n = 0; n < p; n++)
          {
            if (my_rank == n)
            {
              position = 0;
              MPI_Pack(&ndo, 1, MPI_INT, buffer, HUGEE, &position, MPI_COMM_WORLD);
              for (mgi = nstart; mgi < (nstart+ndo); mgi++)
              //for (nncl = 0; nncl < ndo; nncl++)
              {
                //nn = nonemptycells[my_rank+nncl*nprocs];
                MPI_Pack(&mgi, 1, MPI_INT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                //if (cell[nn].rho > MINDENSITY)
                if (modelgrid[mgi].associated_cells > 0)
                {
                  MPI_Pack(&modelgrid[mgi].thick, 1, MPI_SHORT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                  MPI_Pack(&modelgrid[mgi].rho, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                  MPI_Pack(&modelgrid[mgi].nne, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                  MPI_Pack(&modelgrid[mgi].nnetot, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                  MPI_Pack(&modelgrid[mgi].kappagrey, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                  MPI_Pack(&modelgrid[mgi].Te, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                  MPI_Pack(&modelgrid[mgi].TJ, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                  MPI_Pack(&modelgrid[mgi].TR, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
    //                MPI_Pack(&cell[nn].T_D, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                  MPI_Pack(&modelgrid[mgi].W, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
    //                MPI_Pack(&cell[nn].W_D, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                    //MPI_Pack(&cell[nn].samplecell, 1, MPI_INT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                  MPI_Pack(&modelgrid[mgi].totalcooling, 1, MPI_DOUBLE, buffer, HUGEE, &position, MPI_COMM_WORLD);

                  for (element = 0; element < nelements; element++)
                  {
                    MPI_Pack(modelgrid[mgi].composition[element].groundlevelpop, get_nions(element), MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                    MPI_Pack(modelgrid[mgi].composition[element].partfunct, get_nions(element), MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
                    MPI_Pack(modelgrid[mgi].cooling[element].contrib, get_nions(element), MPI_DOUBLE, buffer, HUGEE, &position, MPI_COMM_WORLD);
                  }
                }
              }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(buffer, HUGEE, MPI_PACKED, n, MPI_COMM_WORLD);

            position = 0;
            MPI_Unpack(buffer, HUGEE, &position, &nlp, 1, MPI_INT, MPI_COMM_WORLD);
            for (nn = 0; nn < nlp; nn++)
            {
              MPI_Unpack(buffer, HUGEE, &position, &mgi, 1, MPI_INT, MPI_COMM_WORLD);
              //if (cell[ncl].rho > MINDENSITY)
              if (modelgrid[mgi].associated_cells > 0)
              {
                MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].thick, 1, MPI_SHORT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].rho, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].nne, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].nnetot, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].kappagrey, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].Te, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].TJ, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].TR, 1, MPI_FLOAT, MPI_COMM_WORLD);
                //MPI_Unpack(buffer, HUGEE, &position, &cell[ncl].T_D, 1, MPI_FLOAT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].W, 1, MPI_FLOAT, MPI_COMM_WORLD);
                //MPI_Unpack(buffer, HUGEE, &position, &cell[ncl].W_D, 1, MPI_FLOAT, MPI_COMM_WORLD);
                //MPI_Unpack(buffer, HUGEE, &position, &cell[ncl].samplecell, 1, MPI_INT, MPI_COMM_WORLD);
                MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].totalcooling, 1, MPI_DOUBLE, MPI_COMM_WORLD);

                for (element = 0; element < nelements; element++)
                {
                  MPI_Unpack(buffer, HUGEE, &position, modelgrid[mgi].composition[element].groundlevelpop, get_nions(element), MPI_FLOAT, MPI_COMM_WORLD);
                  MPI_Unpack(buffer, HUGEE, &position, modelgrid[mgi].composition[element].partfunct, get_nions(element), MPI_FLOAT, MPI_COMM_WORLD);
                  MPI_Unpack(buffer, HUGEE, &position, modelgrid[mgi].cooling[element].contrib, get_nions(element), MPI_DOUBLE, MPI_COMM_WORLD);
                }
              }
            }
          }

          #ifndef FORCE_LTE
            if (continue_simulation && nts-itstep == 0 && titer == 0)
            {
            }
            else
            {
              /// Reduce the corrphotoionrenorm array.
              printout("nts %d, titer %d: bcast corr photoionrenorm\n",nts,titer);
              MPI_Reduce(&corrphotoionrenorm, &redhelper, MMODELGRID*nelements*maxion, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID*nelements*maxion; i++)
                {
                  corrphotoionrenorm[i] = redhelper[i];
                }
              }
              MPI_Bcast(&corrphotoionrenorm, MMODELGRID*nelements*maxion, MPI_DOUBLE, 0, MPI_COMM_WORLD);

              /// Reduce the gammaestimator array. Only needed to write restart data.
              printout("nts %d, titer %d: bcast gammaestimator\n",nts,titer);
              MPI_Reduce(&gammaestimator, &redhelper, MMODELGRID*nelements*maxion, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID*nelements*maxion; i++)
                {
                  gammaestimator[i] = redhelper[i];
                }
              }
              MPI_Bcast(&gammaestimator, MMODELGRID*nelements*maxion, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
          #endif
        #endif

        /// If this is not the 0th time step of the current job step,
        /// write out a snapshot of the grid properties for further restarts
        /// and update input.txt accordingly
        if (nts-itstep != 0)
        {
          if (my_rank == 0)
          {
            printout("Write grid restart data\n");
            write_grid_restart_data();
            printout("Update input.txt for restart at time step %d\n",nts);
            update_parameterfile(nts);
            printout("input.txt successfully updated\n");
          }
        }


        /// Do this - which is only an initialisation and no need of calculation - outside update_grid to avoid communication
/*        for (n = 0; n < npts_model; n++)
        {
          modelgrid[n].totalcooling = COOLING_UNDEFINED;
        }*/
        printout("time after grid properties have been communicated %d\n",time(NULL));


        /** set all the estimators to zero before moving packets. This is now done
        after update_grid so that, if requires, the gamma-ray heating estimator is known there
        and also the photoion and stimrecomb estimators */
        zero_estimators();

        if ((nts < ftstep) && (do_this_full_loop == 1))
        {
          /// Now process the packets.
          printout("time before update packets %d\n",time(NULL));
          update_packets(nts);


          /*
          for (middle_iteration = 0; middle_iteration < n_middle_it; middle_iteration++)
          {
            /// Read in the next bunch of packets to work on.
            if (n_middle_it > 1)
            {
              sprintf(filename,"packets%d_%d.tmp",middle_iteration,my_rank);
              if ((packets_file = fopen(filename, "rb")) == NULL)
              {
                printf("Cannot open packets file\n");
                exit(0);
              }
              fread(&pkt[0], sizeof(PKT), npkts, packets_file);
              //read_packets(packets_file);
              fclose(packets_file);
            }

            /// Update those packets ...
            //sprintf(filename,"tau%d.out",nts);
            //if ((tau_file = fopen(filename, "w")) == NULL)
            //{
            //  printf("Cannot open %s.\n",filename);
            //  abort();
            //}
            update_packets(nts);
            //fclose(tau_file);

            /// And save their new state back to disc before proceeding with the next bunch of packets.
            if (n_middle_it > 1)
            {
              sprintf(filename,"packets%d_%d.tmp",middle_iteration,my_rank);
              if ((packets_file = fopen(filename, "wb")) == NULL)
              {
                printf("Cannot open packets file\n");
                exit(0);
              }
              fwrite(&pkt[0], sizeof(PKT), npkts, packets_file);
              fclose(packets_file);
            }
            //if (nts % 6 == 0 || nts == 49)
            if (nts == 49)
            {
              sprintf(filename,"packets%.2d_%.4d.out",nts,my_rank);
              //sprintf(filename,"packets%.2d_%.4d.out",middle_iteration,my_rank);
              if ((packets_file = fopen(filename, "w")) == NULL)
              {
                printf("Cannot open packets file\n");
                exit(0);
              }
              write_packets(packets_file);
              fclose(packets_file);
            }
          }
          */

          /// Calculate mean interaction per packet
          interactions = 0;
          for (i = 0; i < npkts; i++)
          {
            interactions += pkt[i].interactions;
          }
          meaninteractions = interactions / npkts;
          printout("mean number of interactions per packet = %g\n",meaninteractions);

          /// Printout packet statistics
          printout("ma_stat_activation_collexc = %d\n",ma_stat_activation_collexc);
          printout("ma_stat_activation_collion = %d\n",ma_stat_activation_collion);
          printout("ma_stat_activation_bb = %d\n",ma_stat_activation_bb);
          printout("ma_stat_activation_bf = %d\n",ma_stat_activation_bf);
          printout("ma_stat_deactivation_colldeexc = %d\n",ma_stat_deactivation_colldeexc);
          printout("ma_stat_deactivation_collrecomb = %d\n",ma_stat_deactivation_collrecomb);
          printout("ma_stat_deactivation_bb = %d\n",ma_stat_deactivation_bb);
          printout("ma_stat_deactivation_fb = %d\n",ma_stat_deactivation_fb);

          printout("k_stat_to_ma_collexc = %d\n",k_stat_to_ma_collexc);
          printout("k_stat_to_ma_collion = %d\n",k_stat_to_ma_collion);
          printout("k_stat_to_r_ff = %d\n",k_stat_to_r_ff);
          printout("k_stat_to_r_fb = %d\n",k_stat_to_r_fb);
          printout("k_stat_to_r_bb = %d\n",k_stat_to_r_bb);
          printout("k_stat_from_ff = %d\n",k_stat_from_ff);
          printout("k_stat_from_bf = %d\n",k_stat_from_bf);
          printout("k_stat_from_gamma = %d\n",k_stat_from_gamma);
          printout("k_stat_from_eminus = %d\n",k_stat_from_eminus);
          printout("k_stat_from_earlierdecay = %d\n",k_stat_from_earlierdecay);

          printout("escounter = %d\n",escounter);
          printout("cellcrossing  = %d\n",cellcrossings);
          printout("updatecellcounter  = %d\n",updatecellcounter);
          printout("coolingratecalccounter = %d\n",coolingratecalccounter);
          printout("resonancescatterings  = %d\n",resonancescatterings);

          printout("upscatterings  = %d\n",upscatter);
          printout("downscatterings  = %d\n",downscatter);


          #ifdef MPI_ON
            MPI_Barrier(MPI_COMM_WORLD); ///hold all processes once the packets are updated
          #endif
          printout("time after update packets %d\n",time(NULL));
          //exit(0);

          #ifdef MPI_ON
            /** All the processes have their own versions of the estimators for this time step now.
            Since these are going to be needed in the next time step, we will gather all the
            estimators together now, sum them, normalise on the Master thread and then pass back to the
            others*/

            /// the following blocks gather all the estimators to the zeroth (Master) thread
            MPI_Reduce(&energy_deposition, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (my_rank == 0)
            {
              for (i = 0; i < MMODELGRID; i++)
              {
                energy_deposition[i] = redhelper[i];
              }
            }
            MPI_Reduce(&J, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if (my_rank == 0)
            {
              for (i = 0; i < MMODELGRID; i++)
              {
                J[i] = redhelper[i];
              }
            }
	    MPI_Barrier(MPI_COMM_WORLD);
            #ifndef FORCE_LTE
              MPI_Reduce(&nuJ, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  nuJ[i] = redhelper[i];
                }
              }
	      MPI_Barrier(MPI_COMM_WORLD);
              MPI_Reduce(&ffheatingestimator, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  ffheatingestimator[i] = redhelper[i];
                }
              }
	      MPI_Barrier(MPI_COMM_WORLD);
              MPI_Reduce(&colheatingestimator, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  colheatingestimator[i] = redhelper[i];
                }
              }
	      MPI_Barrier(MPI_COMM_WORLD);
              MPI_Reduce(&gammaestimator, &redhelper, MMODELGRID*nelements*maxion, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID*nelements*maxion; i++)
                {
                  gammaestimator[i] = redhelper[i];
                }
              }
	      MPI_Barrier(MPI_COMM_WORLD);
              MPI_Reduce(&bfheatingestimator, &redhelper, MMODELGRID*nelements*maxion, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID*nelements*maxion; i++)
                {
                  bfheatingestimator[i] = redhelper[i];
                }
              }
    /*          MPI_Reduce(&ionfluxestimator, &redhelper, MMODELGRID*nelements*maxion, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID*nelements*maxion; i++)
                {
                  ionfluxestimator[i] = redhelper[i];
                }
              }*/
      /*        MPI_Reduce(&twiddle, &redhelper, MMODELGRID*nelements*maxion, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID*nelements*maxion; i++)
                {
                  twiddle[i] = redhelper[i];
                }
              }*/
    /*          MPI_Reduce(&stimrecombestimator, &redhelper, MMODELGRID*nelements*maxion, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID*nelements*maxion; i++)
                {
                  stimrecombestimator[i] = redhelper[i];
                }
              }*/

    /*          MPI_Reduce(&mabfcount, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  mabfcount[i] = redhelper[i]/p;
                }
              }
              MPI_Reduce(&mabfcount_thermal, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  mabfcount_thermal[i] = redhelper[i]/p;
                }
              }
              MPI_Reduce(&kbfcount, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  kbfcount[i] = redhelper[i]/p;
                }
              }
              MPI_Reduce(&kbfcount_ion, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  kbfcount_ion[i] = redhelper[i]/p;
                }
              }
              MPI_Reduce(&kffcount, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  kffcount[i] = redhelper[i]/p;
                }
              }
              MPI_Reduce(&kffabs, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  kffabs[i] = redhelper[i]/p;
                }
              }
              MPI_Reduce(&kbfabs, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  kbfabs[i] = redhelper[i]/p;
                }
              }
              MPI_Reduce(&kgammadep, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  kgammadep[i] = redhelper[i]/p;
                }
              }*/
            #endif

            #ifdef RECORD_LINESTAT
              MPI_Reduce(ecounter, linestat_reduced, nlines, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < nlines; i++)
                {
                  ecounter[i] = linestat_reduced[i];
                }
              }
              MPI_Reduce(acounter, linestat_reduced, nlines, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < nlines; i++)
                {
                  acounter[i] = linestat_reduced[i];
                }
              }
            #endif

            deltaV =  pow(wid_init * time_step[nts].mid/tmin, 3.0);
            deltat = time_step[nts].width;
            if (do_rlc_est != 0)
            {
              MPI_Reduce(&rpkt_emiss, &redhelper, MMODELGRID, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                for (i = 0; i < MMODELGRID; i++)
                {
                  rpkt_emiss[i] = redhelper[i];
                }
              }
            }
            if (do_comp_est == 1)
            {
              MPI_Reduce(&compton_emiss, &redhelper, MMODELGRID*EMISS_MAX, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
              if (my_rank == 0)
              {
                i = 0;
                for (n = 0; n < MMODELGRID; n++)
                {
                  for (nn = 0; nn < EMISS_MAX; nn++)
                  {
                    compton_emiss[n][nn] = redhelper[i];
                    i++;
                  }
                }
              }
            }
            MPI_Barrier(MPI_COMM_WORLD);
          #endif

          /** The master thread now knows the estimators (avertaged over the processors). It will now normalise them.
          Then the new values can be sent out to all threads again */
          if (my_rank == 0)
          {
            if (do_comp_est == 1)
            {
              normalise_estimators(nts);
              write_estimators(nts);
            }

            if (do_rlc_est != 0)
            {
              normalise_grey(nts);
              if (do_rlc_est != 3)
              {
                write_grey(nts);
              }
            }
          }

          #ifdef MPI_ON
            /** The master thread has normalised the rpkt and compton estimators and printed out a bunch of stuff. Now redistribute the estimators ready for the next run. */

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&J, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            #ifndef FORCE_LTE
              MPI_Bcast(&nuJ, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&ffheatingestimator, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&colheatingestimator, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&gammaestimator, MMODELGRID*nelements*maxion, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&bfheatingestimator, MMODELGRID*nelements*maxion, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              //MPI_Bcast(&photoionestimator, MMODELGRID*nelements*maxion, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              //MPI_Bcast(&stimrecombestimator, MMODELGRID*nelements*maxion, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              /*
              MPI_Bcast(&ionfluxestimator, MMODELGRID*nelements*maxion, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              //MPI_Bcast(&twiddle, MMODELGRID*nelements*maxion, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&mabfcount, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&mabfcount_thermal, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&kbfcount, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&kbfcount_ion, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&kffcount, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&kffabs, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&kbfabs, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              MPI_Bcast(&kgammadep, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              */
            #endif
            if (do_rlc_est != 0)
            {
              MPI_Bcast(&rpkt_emiss, MMODELGRID, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            }
            if (do_comp_est == 1)
            {
              MPI_Bcast(&compton_emiss, MMODELGRID*EMISS_MAX, MPI_FLOAT, 0, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);
          #endif

          /// Now printout some statistics on the current timestep
          printout("time after estimators have been communicated %d\n",time(NULL));
          printout("%d: During timestep %d on MPI process %d, %d pellets decayed and %d packets escaped. (time %g)\n",outer_iteration,nts,my_rank,time_step[nts].pellet_decays,nesc,time_step[nts].mid/DAY);

          #ifdef VPKT_ON

            printout("%d: During timestep %d on MPI process %d, %d virtual packets were generated and %d escaped. \n",outer_iteration,nts,my_rank,nvpkt,nvpkt_esc1+nvpkt_esc2+nvpkt_esc3);
            printout("%d virtual packets came from an electron scattering event, %d from a kpkt deactivation and %d from a macroatom deactivation. \n",nvpkt_esc1,nvpkt_esc2,nvpkt_esc3);

            nvpkt = 0 ;
            nvpkt_esc1 = 0 ;
            nvpkt_esc2 = 0 ;
            nvpkt_esc3 = 0 ;

          #endif


          #ifdef RECORD_LINESTAT
            if (my_rank == 0)
            {
              /// Print net absorption/emission in lines to the linestat_file
              /// Currently linestat information is only properly implemented for MPI only runs
              /// For hybrid runs only data from thread 0 is recorded
              for (i = 0; i < nlines; i++) fprintf(linestat_file,"%d ", ecounter[i]);
              fprintf(linestat_file,"\n");
              for (i = 0; i < nlines; i++) fprintf(linestat_file,"%d ", acounter[i]);
              fprintf(linestat_file,"\n");

              ///Old style
              //for (i = 0; i < nlines; i++) fprintf(linestat_file,"%g %d %d %d %d %d %d\n", CLIGHT/linelist[i].nu, get_element(linelist[i].elementindex), get_ionstage(linelist[i].elementindex,linelist[i].ionindex), linelist[i].upperlevelindex+1, linelist[i].lowerlevelindex+1,ecounter_reduced[i],acounter_reduced[i]);
      	    }
          #endif

          printout("time before write temporary packets file %d\n",time(NULL));

          if (nts % 2 == 0) sprintf(filename,"packets%d_%d_even.tmp",tid,my_rank);
          else sprintf(filename,"packets%d_%d_odd.tmp",tid,my_rank);

          if ((packets_file = fopen(filename, "wb")) == NULL)
          {
            printout("Cannot write to temporary packets file %s\n",filename);
            exit(0);
          }

          fwrite(&pkt[0], sizeof(PKT), npkts, packets_file);
          fclose(packets_file);

          #ifdef VPKT_ON
          {
            if (nts % 2 == 0) sprintf(filename,"vspecpol_%d_%d_even.tmp",my_rank,tid);
            else sprintf(filename,"vspecpol_%d_%d_odd.tmp",my_rank,tid);

            if ((packets_file = fopen(filename, "wb")) == NULL)
            {
              printout("Cannot write to temporary packets file %s\n",filename);
              exit(0);
            }

            write_vspecpol(packets_file);
            fclose(packets_file);

            // Write temporary files for vpkt_grid

            if (vgrid_flag==1) {

                if (nts % 2 == 0) sprintf(filename,"vpkt_grid_%d_%d_even.tmp",my_rank,tid);
                else sprintf(filename,"vpkt_grid_%d_%d_odd.tmp",my_rank,tid);

                if ((packets_file = fopen(filename, "wb")) == NULL)
                {
                    printout("Cannot write to vpkt_grid file %s\n",filename);
                    exit(0);
                }

                write_vpkt_grid(packets_file);
                fclose(packets_file);

            }
          }
          #endif


          printout("time after write temporary packets file %d\n",time(NULL));

          if (nts == ftstep-1)
          {
            sprintf(filename,"packets%.2d_%.4d.out",0,my_rank);
            //sprintf(filename,"packets%.2d_%.4d.out",middle_iteration,my_rank);
            if ((packets_file = fopen(filename, "w")) == NULL)
            {
              printf("Cannot write to final packets file %s\n",filename);
              exit(0);
            }
            write_packets(packets_file);
            fclose(packets_file);

            // write specpol of the virtual packets
            #ifdef VPKT_ON
              write_vspecpol(vspecpol_file);
              fclose(vspecpol_file);

              if (vgrid_flag==1)
              {
                  write_vpkt_grid(vpkt_grid_file);
                  fclose(vpkt_grid_file);
              }
            #endif

            printout("time after write final packets file %d\n",time(NULL));
          }
          /*if (nts % 6 == 0 || nts == 49)
          {
            sprintf(filename,"packets%.2d_%.4d.out",nts,my_rank);
            //sprintf(filename,"packets%.2d_%.4d.out",middle_iteration,my_rank);
            if ((packets_file = fopen(filename, "w")) == NULL)
            {
              printf("Cannot open packets file\n");
              exit(0);
            }
            write_packets(packets_file);
            fclose(packets_file);
          }*/

      }
      }

      nts++;
      if (do_this_full_loop == 0)
      {
        nts += last_loop+1; ///this will break the loop and terminate the code
      }

    }



    /// Now write a snapshot of the current model data to file to allow
    /// further continuation of the simulation.
    /// This works only for n_out_it=1 as larger n_out_it would require
    /// updating the grid over the outer iterations which is not done!
    ///------------------------------------------------------------------------

    /// The final state of the packets has already been written to
    /// the temporary packets files after the last call of
    /// update_packets for each middle-iteration by each process.
    /// As these files are binary they are not portable between
    /// different machines. To overcome this we write them here
    /// once again to formatted files.
    /*
    for (middle_iteration = 0; middle_iteration < n_middle_it; middle_iteration++)
    {
      sprintf(filename,"packets%d_%.4d.out",middle_iteration,my_rank);
      if ((packets_file = fopen(filename, "w")) == NULL)
      {
        printf("Cannot open packets file\n");
        exit(0);
      }
      write_finalpackets(packets_file);
      fclose(packets_file);
    }
    */

    /// Furthermore we need to update the grids quantities for tmax
    /// (indicated by ftstep) and write them to _one_ file using the
    /// master process. titer is here 0
    //if (ftstep < ntstep) update_grid(ftstep,my_rank,nstart,ndo,0);

    /*
    #ifdef MPI_ON
      /// Each process has now updated its own set of cells.
      /// The results now need to be communicated between processes.
      MPI_Barrier(MPI_COMM_WORLD);
      printout("time before final grid comm %d\n",time(NULL));
      for (n = 0; n < p; n++)
      {
        if (my_rank == n)
        {
          position = 0;
          MPI_Pack(&ndo, 1, MPI_INT, buffer, HUGEE, &position, MPI_COMM_WORLD);
          for (mgi = nstart; mgi < (nstart+ndo); mgi++)
          {
            MPI_Pack(&mgi, 1, MPI_INT, buffer, HUGEE, &position, MPI_COMM_WORLD);
            if (modelgrid[mgi].associated_cells > 0)
            {
              MPI_Pack(&modelgrid[mgi].Te, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
              MPI_Pack(&modelgrid[mgi].TJ, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
              MPI_Pack(&modelgrid[mgi].TR, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
              MPI_Pack(&modelgrid[mgi].W, 1, MPI_FLOAT, buffer, HUGEE, &position, MPI_COMM_WORLD);
            }
          }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(buffer, HUGEE, MPI_PACKED, n, MPI_COMM_WORLD);

        position = 0;
        MPI_Unpack(buffer, HUGEE, &position, &nlp, 1, MPI_INT, MPI_COMM_WORLD);
        for (nn = 0; nn < nlp; nn++)
        {
          MPI_Unpack(buffer, HUGEE, &position, &mgi, 1, MPI_INT, MPI_COMM_WORLD);
          if (modelgrid[mgi].associated_cells > 0)
          {
            MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].Te, 1, MPI_FLOAT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].TJ, 1, MPI_FLOAT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].TR, 1, MPI_FLOAT, MPI_COMM_WORLD);
            MPI_Unpack(buffer, HUGEE, &position, &modelgrid[mgi].W, 1, MPI_FLOAT, MPI_COMM_WORLD);
          }
        }
      }
      printout("time after final grid comm %d\n",time(NULL));
    #endif

    if (my_rank == 0)
    {
      if ((temperature_file = fopen("final_temperatures.dat", "w")) == NULL)
      {
        printf("Cannot open final_temperatures.dat\n");
        exit(0);
      }
      for (n = 0; n < ngrid; n++)
      {
        int mgi = cell[n].modelgridindex;
        if (get_rho(mgi)  > MINDENSITY)  /// plot only nonempty cells
        {
          fprintf(temperature_file,"%d %g %g %g %g\n",n,get_TR(mgi),get_Te(mgi),get_W(mgi),get_TJ(mgi));
        }
      }
      fclose(temperature_file);
    }
    */


    /// The main calculation is now over. The packets now have all stored the time, place and direction
    /// at which they left the grid. Also their rest frame energies and frequencies.
    /// Spectra and light curves are now extracted using exspec which is another make target of this
    /// code.

    #ifdef MPI_ON
      free(buffer);
    #endif
  }

  if (my_rank == 0)
  {
    fclose(linestat_file);
  }
  //fclose(ldist_file);
  //fclose(output_file);

  /* Spec syn. */
  //grid_init();
  //syn_gamma();

  if ((ntstep != ftstep) || (do_this_full_loop == 0))
  {
    printout("RESTART_NEEDED to continue model\n");
  }
  else
  {
    printout("No need for restart\n");
  }



  #ifdef MPI_ON
    /// Communicate gamma and positron deposition and write to file
    for (i=0; i < ntstep; i++)
    {
      MPI_Reduce(&time_step[i].gamma_dep, &depvalue ,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (my_rank == 0) time_step[i].gamma_dep = depvalue/nprocs;
      MPI_Reduce(&time_step[i].positron_dep, &depvalue ,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (my_rank == 0) time_step[i].positron_dep = depvalue/nprocs;
      MPI_Reduce(&time_step[i].dep, &depvalue ,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
      if (my_rank == 0) time_step[i].dep = depvalue/nprocs;
    }
  #endif
  if (my_rank == 0)
  {
    if ((dep_file = fopen("deposition.out", "w")) == NULL)
    {
      printf("Cannot open deposition.out\n");
      abort();
    }
    for (i=0; i < ntstep; i++)
    {
      fprintf(dep_file, "%g %g %g %g\n", time_step[i].mid/DAY, time_step[i].gamma_dep/time_step[i].width/LSUN, time_step[i].positron_dep/time_step[i].width/LSUN, time_step[i].dep/time_step[i].width/LSUN);
    }
    fclose(dep_file);
  }


  printout("simulation finished at %d\n",time(NULL));
  //fclose(tb_file);
  fclose(estimators_file);

  #ifdef _OPENMP
    #pragma omp parallel
    {
      fclose(output_file);
    }
  #else
      fclose(output_file);
  #endif

  #ifdef MPI_ON
    MPI_Finalize();
  #endif

  return 0;
}


/// Generalized output routine
/// following section 7.3 of "C. Programming Language."
/// by Brian W. Kernighan and Dennis Ritchie
/// As it stands it is only capable to printout floating point variables as %g
/// specifiers which determine the number of digits don't work!
void printout(char *fmt, ...)
{
  va_list ap;
  char *p, *sval;
  int ival;
  double dval;
  char filename[100];
  //FILE *output_file;

  /// To be able to follow the messages interactively the file is continuously
  /// opened and closed. As this happens always with the "a" argument the file
  /// is so far never really initialized: a existing output.txt will be continued!
/*      sprintf(filename,"output_%d-%d.txt",rank_global,tid);
      if (output_file_open == 0)
      {
        if ((output_file = fopen(filename, "a")) == NULL)
        {
          printf("Cannot open %s.\n",filename);
          abort();
  //    exit(0);
        }
        output_file_open = 1;
      }*/

      va_start(ap, fmt);
      for (p = fmt; *p; p++)
      {
        if (*p != '%')
        {
      //putchar(*p); ///for output on the screen
          fputc(*p,output_file);
          continue;
        }
        switch (*++p)
        {
          case 'd':
            ival = va_arg(ap, int);
        //printf(output_file,"%d", ival);  ///for output on the screen
            fprintf(output_file,"%d", ival);
            break;
          case 'f':
            dval = va_arg(ap, double);
            fprintf(output_file,"%f", dval);
            break;
          case 'g':
            dval = va_arg(ap, double);
            fprintf(output_file,"%g", dval);
            break;
          case 's':
            for (sval = va_arg(ap, char *); *sval; sval++)
              fputc(*sval,output_file);
            break;
          default:
            fputc(*p,output_file);
            break;
        }
      }
      va_end(ap);


      //fclose(output_file);
}











///****************************************************************************
/*void print_opticaldepth(int cellnumber, int timestep, int samplecell, int element)
{
  double epsilon(int element, int ion, int level);
  double stat_weight(int element, int ion, int level);
  double calculate_exclevelpop(int cellnumber, int element, int ion, int level);
  double einstein_spontaneous_emission(int lineindex);

  int nions,nlevels,nuptrans;
  int ion,lower,upper,lineindex;
  int i;
  double epsilon_lower,nu_trans;
  double n_u,n_l,tau_line,A_ul,B_ul,B_lu;
  FILE *tau_file;

  double t_current = time_step[timestep].mid;
  double T_e = cell[cellnumber].T_e;
  double T_R = cell[cellnumber].T_R;
  double W = cell[cellnumber].W;



  sprintf(filename,"tau%.2d_sc%.2d.out",timestep,samplecell);
  if ((tau_file = fopen(filename, "w")) == NULL)
  {
    printf("Cannot open %s.\n",filename);
    abort();
  }
  setvbuf(tau_file, NULL, _IOLBF, 1);
  fprintf(tau_file,"%d %g %g %g\n",samplecell,T_e,T_R,W);

  nions = get_nions(element);
  for (ion = 0; ion < nions; ion++)
  {
    nlevels = get_nlevels(element,ion);
    for (lower = 0; lower < nlevels; lower++)
    {
      epsilon_lower = epsilon(element,ion,lower);
      nuptrans = elements[element].ions[ion].levels[lower].uptrans[0].targetlevel;
      for (i = 1; i <= nuptrans; i++)
      {
        upper = elements[element].ions[ion].levels[lower].uptrans[i].targetlevel;
        lineindex = elements[element].ions[ion].levels[lower].uptrans[i].lineindex;
        nu_trans = (elements[element].ions[ion].levels[lower].uptrans[i].epsilon - epsilon_lower)/H;

        A_ul = einstein_spontaneous_emission(lineindex);
        B_ul = CLIGHTSQUAREDOVERTWOH / pow(nu_trans,3) * A_ul;
        B_lu = stat_weight(element,ion,upper)/stat_weight(element,ion,lower) * B_ul;

        n_u = calculate_exclevelpop(cellnumber,element,ion,upper);
        n_l = calculate_exclevelpop(cellnumber,element,ion,lower);
        tau_line = (B_lu*n_l - B_ul*n_u) * HCLIGHTOVERFOURPI * t_current;

        fprintf(tau_file,"%g %g %d\n",nu_trans,tau_line,ion);
      }
    }
  }
  fclose(tau_file);
}*/




/*void *my_malloc(size_t size)
{
  //char *adr;
  void *adr;
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

