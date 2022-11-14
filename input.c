//#include <stdio.h>
//#include <math.h>
//#include <stdlib.h>
#include "sn3d.h"
#ifdef DO_EXSPEC
  #include "exspec.h"
#endif
#ifdef VPKT_ON
  #include "vpkt.h"
#endif
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/// To govern the input. For now hardwire everything.
int input(rank)
     int rank;
{
  void read_atomicdata();
  void read_parameterfile();
  void read_parameterfile_vpkt();
  int read_1d_model();
  int read_2d_model();
  int read_3d_model();
  int energy_init();
  int get_gam_ll();
  int get_nul();

  FILE *co_lines;
  FILE *ni_lines;
  FILE *v48_lines;
  FILE *cr48_lines;

  int dum1, n, m, nn;
  float dum2, dum3;
  int lindex_min, lindex_max;

  homogeneous_abundances = 0;
  t_model = 0.0;

  /// Constants for van-Regemorter approximation
  C_0 = 5.465e-11;
  H_ionpot = 13.5979996*EV;

  /// Select grid type
  grid_type = GRID_UNIFORM;
  model_type = RHO_UNIFORM;

  maxion = MIONS;
  /// Set grid size
  nxgrid = 100; //pow(MGRID,1./3.); //10;
  nygrid = 100; //pow(MGRID,1./3.); //10;
  nzgrid = 100; //pow(MGRID,1./3.); //10;
  printout("nxgrid %d\n",nxgrid);
  /*nxgrid = 4;
  nygrid = 4;
  nzgrid = 4;*/
  ngrid = nxgrid * nygrid * nzgrid; ///Moved to input.c
  if (ngrid > MGRID)
  {
    printout("[fatal] input: Error: too many grid cells. Abort.");
    exit(0);
  }


  /// Set number of packets, outer and middle iterations
  npkts = MPKTS;
  n_out_it = 10;
  n_middle_it = 1;
/*  #ifdef FORCE_LTE
    n_titer = 1;
  #else
    n_titer = 6;
  #endif*/
  n_titer = 1;
  initial_iteration = 0;

  printout("[info] input: do n_titer %d iterations per timestep\n",n_titer);
  if (n_titer > 1)
  {
    #ifndef DO_TITER
      printout("[fatal] input: n_titer > 1, but DO_TITER not defined ... abort\n");
      exit(0);
    #endif
  }
  else if (n_titer == 1)
  {
    #ifdef DO_TITER
      printout("[warning] input: n_titer = 1 but DO_TITER defined, remove DO_TITER to save memory\n");
    #endif
  }
  else
  {
    printout("[fatal] input: no valid value for n_titer selected\n");
    exit(0);
  }


  nu_min_r = 1e14;   /// lower frequency boundary for UVOIR spectra and BB sampling
  nu_max_r = 5e15;   /// upper frequency boundary for UVOIR spectra and BB sampling
  #ifdef DO_EXSPEC
    /// Spectra settings
    nnubins = MNUBINS; //1000;  /// frequency bins for spectrum
    /// lower and upper frequency boundaries for make_spectrum_res: gamma spectra?
    nu_min = 0.05 * MEV / H;
    nu_max = 4 * MEV / H;
  #endif

  /// Lightcurve setting
  do_r_lc = 0;    /// default to no lc = gamma-ray spectrum
  do_rlc_est = 0; /// ^^


  nfake_gam = 1; ///# of fake gamma ray lines for syn




  /// Read in parameters from input.txt
  ///======================================================
  read_parameterfile(rank);
  ntbins = ntstep;   ///time bins for spectrum equal #(timesteps)
  ntlcbins = ntstep; ///time bins for light curve #(timesteps)

  /// Read in parameters from vpkt.txt
  #ifdef VPKT_ON
    read_parameterfile_vpkt();
  #endif

  /// Read in atomic data
  ///======================================================
  read_atomicdata();

  //#ifndef DO_EXSPEC
    /// Read in input model
    ///======================================================
    if (model_type == RHO_UNIFORM)
    {
      mtot = 1.39 * MSUN;
      mni56 = 0.625 * MSUN;
      vmax = 1.e9;
      rmax = vmax * tmin;
  //  rhotot = 3 * mtot / 4 / PI / rmax /rmax /rmax; //MK
      xmax = ymax = zmax = rmax;
    }
    else if (model_type == RHO_1D_READ)
    {
      printout("Read 1D model!\n");
      read_1d_model();
    }
    else if (model_type == RHO_2D_READ)
    {
      printout("Read 2D model!\n");
      read_2d_model();
    }
    else if (model_type == RHO_3D_READ)
    {
      printout("Read 3D model!\n");
      read_3d_model();
    }
    else
    {
      printout("Unknown model. Abort.\n");
      exit(0);
    }
  //#endif

  #ifdef USE_ENERGYINPUTFILE
    /// If using energy input files get cell energies
    energy_init();
  #endif

  /// Read in data for gamma ray lines.
  ///======================================================
  dum1=0;
  if ((co_lines = fopen("co_lines.txt", "r")) == NULL)
  {
    printout("Cannot open co_lines.txt.\n");
    exit(0);
  }
  fscanf(co_lines, "%d", &dum1);
  cobalt_spec.nlines = dum1;

  ECOBALT = 0.0;
  for (n = 0; n < dum1; n++)
  {
    fscanf(co_lines, "%g %g", &dum2, &dum3);
    cobalt_spec.energy[n] = dum2 * MEV;
    cobalt_spec.probability[n] = dum3;
    ECOBALT += dum2 * MEV * dum3;
  }
  /// Average energy per gamma line of Co56 decay and positron annihilation
  ECOBALT_GAMMA = ECOBALT;
  /// For total deposited energy we need to add the kinetic energy per emitted positron
  ECOBALT += 0.63*MEV * 0.19;
  //ECOBALT = ECOBALT/5; /// DEBUGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

  fclose(co_lines);

  if ((ni_lines = fopen("ni_lines.txt", "r")) == NULL)
  {
    printout("Cannot open ni_lines.txt.\n");
    exit(0);
  }
  fscanf(ni_lines, "%d", &dum1);
  nickel_spec.nlines = dum1;

  ENICKEL = 0.0;
  for (n = 0; n < dum1; n++)
  {
    fscanf(ni_lines, "%g %g", &dum2, &dum3);
    nickel_spec.energy[n] = dum2 * MEV;
    nickel_spec.probability[n] = dum3;
    ENICKEL += dum2 * MEV * dum3;
  }
  //ENICKEL = ENICKEL/5; /// DEBUGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

  fclose(ni_lines);

  if ((v48_lines = fopen("v48_lines.txt", "r")) == NULL)
  {
    printout("Cannot open v48_lines.txt.\n");
    exit(0);
  }
  fscanf(v48_lines, "%d", &dum1);
  v48_spec.nlines = dum1;

  E48V = 0.0;
  for (n = 0; n < dum1; n++)
  {
    fscanf(v48_lines, "%g %g", &dum2, &dum3);
    v48_spec.energy[n] = dum2 * MEV;
    v48_spec.probability[n] = dum3;
    E48V += dum2 * MEV * dum3;
  }

  fclose(v48_lines);

  if ((cr48_lines = fopen("cr48_lines.txt", "r")) == NULL)
  {
    printout("Cannot open cr48_lines.txt.\n");
    exit(0);
  }
  fscanf(cr48_lines, "%d", &dum1);
  cr48_spec.nlines = dum1;

  E48CR = 0.0;
  for (n = 0; n < dum1; n++)
  {
    fscanf(cr48_lines, "%g %g", &dum2, &dum3);
    cr48_spec.energy[n] = dum2 * MEV;
    cr48_spec.probability[n] = dum3;
    E48CR += dum2 * MEV * dum3;
  }

  fclose(cr48_lines);

  /// With the lines all ready in, now make a list of them in energy order.
  get_gam_ll();

  /// Now that the list exists use it to find values for spectral synthesis
  /// stuff.
  lindex_max = get_nul(nusyn_max);
  lindex_min = get_nul(nusyn_min);
  printout("lindex_max %d, lindex_min %d\n", lindex_max, lindex_min);

  emiss_offset = lindex_min;
  emiss_max = lindex_max - lindex_min + 1;
  printout("emiss_max using %d of a possible %d\n", emiss_max, EMISS_MAX);

  if (emiss_max > EMISS_MAX)
  {
    printout("Too many points needed for emissivities. Use smaller frequency range or increase EMISS_MAX. Abort.\n");
    exit(0);
  }


  #ifdef DO_EXSPEC
    /// Check if enough  memory for spectra has been assigned
    /// and allocate memory for the emission statistics
    if (nnubins > MNUBINS)
    {
      printout("Too many frequency bins in spectrum - reducing.\n");
      nnubins = MNUBINS;
    }
    if (ntbins > MTBINS)
    {
      printout("Too many time bins in spectrum - reducing.\n");
      ntbins = MTBINS;
    }
    for (n = 0; n < ntbins; n++)
    {
      for (m=0; m < nnubins; m++)
      {
        if ((spectra[n].stat[m].absorption = (double *) malloc((nelements*maxion)*sizeof(double))) == NULL)
        {
          printout("[fatal] input: not enough memory to spectrum structure ... abort\n");
          exit(0);
        }
        if ((spectra[n].stat[m].emission = (double *) malloc((2*nelements*maxion+1)*sizeof(double))) == NULL)
        {
          printout("[fatal] input: not enough memory to spectrum structure ... abort\n");
          exit(0);
        }


        #ifdef POL_ON

        if ((stokes_i[n].stat[m].absorption = (double *) malloc((nelements*maxion)*sizeof(double))) == NULL)
        {
            printout("[fatal] input: not enough memory to spectrum structure ... abort\n");
            exit(0);
        }
        if ((stokes_i[n].stat[m].emission = (double *) malloc((2*nelements*maxion+1)*sizeof(double))) ==NULL)
        {
            printout("[fatal] input: not enough memory to spectrum structure ... abort\n");
            exit(0);
        }

        if ((stokes_q[n].stat[m].absorption = (double *) malloc((nelements*maxion)*sizeof(double))) == NULL)
        {
            printout("[fatal] input: not enough memory to spectrum structure ... abort\n");
            exit(0);
        }
        if ((stokes_q[n].stat[m].emission = (double *) malloc((2*nelements*maxion+1)*sizeof(double))) ==NULL)
        {
            printout("[fatal] input: not enough memory to spectrum structure ... abort\n");
            exit(0);
        }

        if ((stokes_u[n].stat[m].absorption = (double *) malloc((nelements*maxion)*sizeof(double))) == NULL)
        {
            printout("[fatal] input: not enough memory to spectrum structure ... abort\n");
            exit(0);
        }
        if ((stokes_u[n].stat[m].emission = (double *) malloc((2*nelements*maxion+1)*sizeof(double))) == NULL)
        {
            printout("[fatal] input: not enough memory to spectrum structure ... abort\n");
            exit(0);
        }

        #endif


        /*
        if (do_emission_res == 1)
        {
          for (nn = 0; nn < MABINS; nn++)
          {
            if ((spectra_res[n][nn].emission[m].count = (int *) malloc((2*nelements*maxion+1)*sizeof(int))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialise spectra_res structure ... abort\n");
              exit(0);
            }
          }
        }
        */
      }
    }
  #endif

  return(0);
}


///****************************************************************************
/// Subroutine to read in input parameters.
void read_atomicdata()
{
  ///new atomic data scheme by readin of adata////////////////////////////////////////////////////////////////////////
  int search_groundphixslist(double nu_edge, int *index_in_groundlevelcontestimator,int element, int ion, int level);
  //int search_groundphixslist(double nu_edge, int *index_in_groundlevelcontestimator);
  int index_in_groundlevelcontestimator;
  int compare_linelistentry();
  int get_elementindex(int Z);
  double stat_weight(int element, int ion, int level);
  double epsilon(int element, int ion, int level);
  int get_element(int element);
  int get_ionstage(int element, int ion);
  int transitioncheck(int upper, int lower);
  //double einstein_spontaneous_emission(int element, int ion, int upper, int lower);
  int get_nions(int element);
  int get_nlevels(int element, int ion);
  int get_bfcontinua(int element, int ion);

  FILE *compositiondata;
  FILE *adata;
  FILE *transitiondata;
  FILE *phixsdata;
  FILE *modelatom;
  FILE *linelist_file;
  int i, ii, iii;
  double nu_edge;
  int element,ion,level;
  int Z,nions,ionstage,nlevels,nlevelsmax,nlevelsmax_readin,ionisinglevels;//,mainqn;
  int lowermost_ionstage,uppermost_ionstage;
  double abundance,mass;
  int Zcheck,Zcheck2,ionstagecheck; //,nionscheck,checknelements;
  int lineindex, tottransitions,tottransitions_all;
  int ndowntrans,nuptrans,nbfcont;
  int skip_line;
  transitiontable_entry *transitiontable;
  int add; /// Helper variable to count coolingterms per ion

  double *nutable;
  double *phixstable;

  double levelenergy,statweight,A;//,f;
  double energyoffset,ionpot;
  int upper,lower,levelindex,transitionindex,ntransitions;
  double currentlevelenergy,g;
  double A_ul,f_ul,nu_trans;
  double nu,nu_max;

  double epsilon_upper,E_threshold;
  int compare_phixslistentry_bynuedge(const void *p1, const void *p2);
  int compare_groundphixslistentry_bynuedge(const void *p1, const void *p2);

  int upperion,upperlevel,lowerion,lowerlevel,npoints;
  double energy,phixs;
  //FILE *database_file,*interpol_file;
  //char filename1[100],filename2[100];

  int T_preset;

  int totaluptrans = 0;
  int totaldowntrans = 0;
  nbfcontinua = 0;
  includedions = 0;
  int cont_index = -1;

  int *nuparr,*ndownarr;
  int targetlevel;

  //printout("start input.c\n");
  if ((modelatom = fopen("modelatom.dat", "r")) == NULL)
  {
    /// No preprocessed model atom available ==> do that now


    ///open atomic data file
    if ((compositiondata = fopen("compositiondata.txt", "r")) == NULL)
    {
      printout("Cannot open compositiondata.txt.\n");
      exit(0);
    }
    if ((adata = fopen("adata.txt", "r")) == NULL)
    {
      printout("Cannot open adata.txt.\n");
      exit(0);
    }
    //printout("adata open\n");


    /// initialize atomic data structure to number of elements
    fscanf(compositiondata,"%d",&nelements);
    if ((elements = (elementlist_entry *) calloc(nelements, sizeof(elementlist_entry))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize elementlist ... abort\n");
      exit(0);
    }
    //printout("elements initialized\n");

    /// Initialize the linelist
    if ((linelist = calloc(MLINES, sizeof(linelist_entry))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize linelist ... abort\n");
      exit(0);
    }

    /// temperature to determine relevant ionstages
    fscanf(compositiondata,"%d",&T_preset);
    fscanf(compositiondata,"%d",&homogeneous_abundances);
    if (homogeneous_abundances == 1) printout("[info] read_atomicdata: homogeneous abundances as defined in compositiondata.txt are active\n");

    /// open transition data file
    if ((transitiondata = fopen("transitiondata.txt", "r")) == NULL)
    {
      printout("Cannot open transitiondata.txt.\n");
      exit(0);
    }
    lineindex = 0;  ///counter to determine the total number of lines, initialisation

    /// readin
    int nbfcheck = 0;
    int heatingcheck = 0;
    int coolingcheck = 0;
    for (element = 0; element < nelements; element++)
    {
      /// read information about the next element which should be stored to memory
      fscanf(compositiondata,"%d %d %d %d %d %lg %lg",&Z,&nions,&lowermost_ionstage,&uppermost_ionstage,&nlevelsmax_readin,&abundance,&mass);
      printout("readin compositiondata: next element Z %d, nions %d, lowermost %d, uppermost %d, nlevelsmax %d\n",Z,nions,lowermost_ionstage,uppermost_ionstage,nlevelsmax_readin);

      /// write this elements data to memory

      elements[element].anumber = Z;
      elements[element].nions = nions;
      elements[element].abundance = abundance;       /// abundances are expected to be given by mass
      elements[element].mass = mass*MH;
      includedions += nions;

      /// Initialize the elements ionlist
      if ((elements[element].ions = (ionlist_entry *) calloc(nions, sizeof(ionlist_entry))) == NULL)
      {
          printout("[fatal] input: not enough memory to initialize ionlist ... abort\n");
          exit(0);
      }

      /// now read in data for all ions of the current element. before doing so initialize
      /// energy scale for the current element (all level energies are stored relative to
      /// the ground level of the neutral ion)
      energyoffset = 0.;
      ionpot = 0.;
      for (ion = 0; ion < nions; ion++)
      {
        nlevelsmax = nlevelsmax_readin;
        printout("ion %d\n",ion);
        /// calculate the current levels ground level energy
        energyoffset += ionpot;

        /// read information for the elements next ionstage
        fscanf(adata,"%d %d %d %lg",&Zcheck2,&ionstage,&nlevels,&ionpot);
        printout("readin level information for adata: Z %d, ionstage %d, nlevels %d\n",Zcheck2,ionstage,nlevels);
        while (Zcheck2 != Z || ionstage != lowermost_ionstage+ion)
        {
          if (Zcheck2 == Z)
          {
            printout("increasing energyoffset by ionpot %g\n",ionpot);
            energyoffset += ionpot;
          }
          for (i = 0; i < nlevels; i++)
          {
            fscanf(adata,"%d %lg %lg %d",&levelindex,&levelenergy,&statweight,&ntransitions);
          }
          fscanf(adata,"%d %d %d %lg",&Zcheck2,&ionstage,&nlevels,&ionpot);
          printout("proceed through adata: Z %d, ionstage %d, nlevels %d\n",Zcheck2,ionstage,nlevels);
        }
        if (nlevelsmax < 0)
        {
          nlevelsmax = nlevels;
        }
        else if (nlevels >= nlevelsmax)
        {
          printout("[info] read_atomicdata: reduce number of levels from %d to %d for ion %d of element %d\n",nlevels,nlevelsmax,ion,element);
        }
        else
        {
          printout("[warning] read_atomicdata: requested nlevelsmax=%d > nlevels=%d for ion %d of element %d ... reduced nlevelsmax to nlevels\n",nlevelsmax,nlevels,ion,element);
          nlevelsmax = nlevels;
          //printout("[fatal] read_atomicdata: nlevelsmax=%d > nlevels=%d for ion %d of element %d ... reduce nlevelsmax or extend atomic data \n",nlevelsmax,nlevels,ion,element);
        }


        /// and proceed through the transitionlist till we match this ionstage (if it was not the neutral one)
        fscanf(transitiondata,"%d %d %d",&Zcheck,&ionstagecheck,&tottransitions);
        printout("readin transdata: Zcheck %d, ionstagecheck %d, tottransitions %d\n",Zcheck,ionstagecheck,tottransitions);
        while (Zcheck != Z || ionstagecheck != ionstage)
        {
          for (i = 0; i < tottransitions; i++)
            fscanf(transitiondata,"%d %d %d %lg",&transitionindex,&lower,&upper,&A);
          fscanf(transitiondata,"%d %d %d",&Zcheck,&ionstagecheck,&tottransitions);
          printout("proceed through transdata: Zcheck %d, ionstagecheck %d, tottransitions %d\n",Zcheck,ionstagecheck,tottransitions);
        }

        tottransitions_all = tottransitions;
        if (ion == nions-1)
        {
          nlevelsmax = 1;
          tottransitions = 0;
        }

        /// then read in its level and transition data
        if (Zcheck == Z && ionstagecheck == ionstage)
        {
          /// load transition table for the CURRENT ion to temporary memory
          if ((transitiontable = (transitiontable_entry *) calloc(tottransitions, sizeof(transitiontable_entry))) == NULL)
          {
            if (tottransitions > 0)
              {
                printout("[fatal] input: not enough memory to initialize transitiontable ... abort\n");
                exit(0);
              }
          }
          if (tottransitions == 0)
          {
            for (i = 0; i < tottransitions_all; i++)
            {
              fscanf(transitiondata,"%d %d %d %lg",&transitionindex,&lower,&upper,&A);
              //printout("index %d, lower %d, upper %d, A %g\n",transitionindex,lower,upper,A);
            }
          }
          else
          {
            for (i = 0; i < tottransitions; i++)
            {
              fscanf(transitiondata,"%d %d %d %lg",&transitionindex,&lower,&upper,&A);
              transitiontable[i].lower = lower-1;
              transitiontable[i].upper = upper-1;
              transitiontable[i].A = A;
              //printout("index %d, lower %d, upper %d, A %g\n",transitionindex,lower,upper,A);
            }
          }


          /// store the ions data to memory and set up the ions zeta and levellist
          elements[element].ions[ion].ionstage = ionstage;
          elements[element].ions[ion].nlevels = nlevelsmax;
          elements[element].ions[ion].ionisinglevels = 0;
          elements[element].ions[ion].ionpot = ionpot*EV;
//           if ((elements[element].ions[ion].zeta = calloc(TABLESIZE, sizeof(float))) == NULL)
//           {
//             printout("[fatal] input: not enough memory to initialize zetalist for element %d, ion %d ... abort\n",element,ion);
//             exit(0);
//           }
          if ((elements[element].ions[ion].Alpha_sp = calloc(TABLESIZE, sizeof(float))) == NULL)
          {
            printout("[fatal] input: not enough memory to initialize Alpha_sp list for element %d, ion %d ... abort\n",element,ion);
            exit(0);
          }
          if ((elements[element].ions[ion].levels = (levellist_entry *) calloc(nlevelsmax, sizeof(levellist_entry))) == NULL)
          {
            printout("[fatal] input: not enough memory to initialize levelist of element %d, ion %d ... abort\n",element,ion);
            exit(0);
          }


          /// now we need to readout the data for all those levels, write them to memory
          /// and set up the list of possible transitions for each level
          if ((nuparr = calloc(nlevelsmax, sizeof(int))) == NULL)
          {
            printout("[fatal] input: not enough memory to allocate nuparr ... abort\n");
            exit(0);
          }
          if ((ndownarr = calloc(nlevelsmax, sizeof(int))) == NULL)
          {
            printout("[fatal] input: not enough memory to allocate ndownarr... abort\n");
            exit(0);
          }
          if ((transitions = calloc(nlevelsmax, sizeof(transitions_t))) == NULL)
          {
            printout("[fatal] input: not enough memory to allocate transitions ... abort\n");
            exit(0);
          }
          for (level = 0; level < nlevels; level++)
          {
            fscanf(adata,"%d %lg %lg %d",&levelindex,&levelenergy,&statweight,&ntransitions);
            //if (element == 1 && ion == 0) printf("%d %16.10f %g %d\n",levelindex,levelenergy,statweight,ntransitions);
            if (level < nlevelsmax)
            {
              ndownarr[level] = 1;
              nuparr[level] = 1;
              //elements[element].ions[ion].levels[level].epsilon = (energyoffset + levelenergy) * EV;
              currentlevelenergy = (energyoffset + levelenergy) * EV;
              //if (element == 1 && ion == 0) printf("%d %16.10e\n",levelindex,currentlevelenergy);
              //printout("energy for level %d of ionstage %d of element %d is %g\n",level,ionstage,element,currentlevelenergy/EV);
              elements[element].ions[ion].levels[level].epsilon = currentlevelenergy;

              //if (level == 0 && ion == 0) energyoffset = levelenergy;
              elements[element].ions[ion].levels[level].stat_weight = statweight;
              ///Moved to the section with ionising levels below
              //elements[element].ions[ion].levels[level].cont_index = cont_index;
              //cont_index -= 1;
              /// Initialise the metastable flag to true. Set it to false (0) if downward transition exists.
              elements[element].ions[ion].levels[level].metastable = 1;
              //elements[element].ions[ion].levels[level].main_qn = mainqn;

              /// The level contributes to the ionisinglevels if its energy
              /// is below the ionisiation potential and the level doesn't
              /// belong to the topmost ion included.
              /// Rate coefficients are only available for ionising levels.
              if (levelenergy < ionpot && ion < nions-1) ///thats only an option for pure LTE && level < TAKE_N_BFCONTINUA)
              {
                elements[element].ions[ion].ionisinglevels += 1;

                elements[element].ions[ion].levels[level].cont_index = cont_index;
                cont_index -= 1;

                if ((elements[element].ions[ion].levels[level].spontrecombcoeff = calloc(TABLESIZE, sizeof(double))) == NULL)
                {
                  printout("[fatal] input: not enough memory to initialize spontrecombcoeff table for element %d, ion %d, level %d\n",element,ion,level);
                  exit(0);
                }
                if ((elements[element].ions[ion].levels[level].corrphotoioncoeff = calloc(TABLESIZE, sizeof(double))) == NULL)
                {
                  printout("[fatal] input: not enough memory to initialize photoioncoeff table for element %d, ion %d, level %d\n",element,ion,level);
                  exit(0);
                }
                if ((elements[element].ions[ion].levels[level].bfheating_coeff = calloc(TABLESIZE, sizeof(double))) == NULL)
                {
                  printout("[fatal] input: not enough memory to initialize modified_photoioncoeff table for element %d, ion %d, level %d\n",element,ion,level);
                  exit(0);
                }
                if ((elements[element].ions[ion].levels[level].bfcooling_coeff = calloc(TABLESIZE, sizeof(double))) == NULL)
                {
                  printout("[fatal] input: not enough memory to initialize modified_photoioncoeff table for element %d, ion %d, level %d\n",element,ion,level);
                  exit(0);
                }
              }


              /// store the possible downward transitions from the current level in following order to memory
              ///     A_level,level-1; A_level,level-2; ... A_level,1
              /// entries which are not explicitly set are zero (the zero is set/initialized by calloc!)
              if ((transitions[level].to = calloc(level, sizeof(short))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize transitionlist ... abort\n");
                exit(0);
              }
              for (i = 0; i < level; i++)
              {
                transitions[level].to[i] = -99.;
              }
              if ((elements[element].ions[ion].levels[level].downtrans = malloc(sizeof(permittedtransitionlist_entry))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize downtranslist ... abort\n");
                exit(0);
              }
              /// initialize number of downward transitions to zero
              elements[element].ions[ion].levels[level].downtrans[0].targetlevel = 0;
              if ((elements[element].ions[ion].levels[level].uptrans = malloc(sizeof(permittedtransitionlist_entry))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize uptranslist ... abort\n");
                exit(0);
              }
              /// initialize number of upward transitions to zero
              elements[element].ions[ion].levels[level].uptrans[0].targetlevel = 0;
            }
          }

          for (ii = 0; ii < tottransitions; ii++)
          {
            level = transitiontable[ii].upper;
            targetlevel = transitiontable[ii].lower;

            if (level < nlevelsmax)
            {

              //if (level == transitiontable[ii].upper && level-i-1 == transitiontable[ii].lower)
              //{
              //printout("ii %d\n",ii);
              //printout("transtable upper %d, lower %d, A %g, iii %d\n",transitiontable[ii].upper,transitiontable[ii].lower, transitiontable[ii].A,iii);
              /// Make sure that we don't allow duplicate. In that case take only the lines
              /// first occurrence
              if (transitioncheck(level,targetlevel) == -99)
              {
                transitions[level].to[level-targetlevel-1] = 1;
                A_ul = transitiontable[ii].A;
                //elements[element].ions[ion].levels[level].transitions[level-targetlevel-1].einstein_A = A_ul;

                nu_trans = (epsilon(element,ion,level) - epsilon(element,ion,targetlevel)) / H;
                g = stat_weight(element,ion,level)/stat_weight(element,ion,targetlevel);
                f_ul = g * ME*pow(CLIGHT,3)/(8*pow(QE*nu_trans*PI,2)) * A_ul;
                //f_ul = g * OSCSTRENGTHCONVERSION / pow(nu_trans,2) * A_ul;
                //elements[element].ions[ion].levels[level].transitions[level-targetlevel-1].oscillator_strength = g * ME*pow(CLIGHT,3)/(8*pow(QE*nu_trans*PI,2)) * A_ul;

                //printout("lineindex %d, element %d, ion %d, lower %d, upper %d, nu %g\n",lineindex,element,ion,level-i-1,level,nu_trans);
                linelist[lineindex].elementindex = element;
                linelist[lineindex].ionindex = ion;
                linelist[lineindex].lowerlevelindex = targetlevel;
                linelist[lineindex].upperlevelindex = level;
                linelist[lineindex].nu = nu_trans;
                linelist[lineindex].einstein_A = A_ul;
                linelist[lineindex].osc_strength = f_ul;
                lineindex += 1;
                if (lineindex % MLINES == 0)
                {
                  printout("[info] read_atomicdata: increase linelistsize from %d to %d\n",lineindex,lineindex+MLINES);
                  if ((linelist = realloc(linelist, (lineindex+MLINES)*sizeof(linelist_entry))) == NULL)
                  {
                    printout("[fatal] input: not enough memory to reallocate linelist ... abort\n");
                    exit(0);
                  }
                }

                /// This is not a metastable level.
                elements[element].ions[ion].levels[level].metastable = 0;

                elements[element].ions[ion].levels[level].downtrans[0].targetlevel = ndownarr[level];
                if ((elements[element].ions[ion].levels[level].downtrans
                    = realloc(elements[element].ions[ion].levels[level].downtrans, (ndownarr[level]+1)*sizeof(permittedtransitionlist_entry))) == NULL)
                {
                  printout("[fatal] input: not enough memory to reallocate downtranslist ... abort\n");
                  exit(0);
                }
                elements[element].ions[ion].levels[level].downtrans[ndownarr[level]].targetlevel = targetlevel;
                elements[element].ions[ion].levels[level].downtrans[ndownarr[level]].epsilon = epsilon(element,ion,targetlevel);
                elements[element].ions[ion].levels[level].downtrans[ndownarr[level]].stat_weight = stat_weight(element,ion,targetlevel);
                //elements[element].ions[ion].levels[level].downtrans[ndownarr[level]].einstein_A = A_ul;
                //elements[element].ions[ion].levels[level].downtrans[ndownarr[level]].oscillator_strength = f_ul;
                ndownarr[level] += 1;

                elements[element].ions[ion].levels[targetlevel].uptrans[0].targetlevel = nuparr[targetlevel];
                if ((elements[element].ions[ion].levels[targetlevel].uptrans
                    = realloc(elements[element].ions[ion].levels[targetlevel].uptrans, (nuparr[targetlevel]+1)*sizeof(permittedtransitionlist_entry))) == NULL)
                {
                  printout("[fatal] input: not enough memory to reallocate uptranslist ... abort\n");
                  exit(0);
                }
                elements[element].ions[ion].levels[targetlevel].uptrans[nuparr[targetlevel]].targetlevel = level;
                elements[element].ions[ion].levels[targetlevel].uptrans[nuparr[targetlevel]].epsilon = epsilon(element,ion,level);
                elements[element].ions[ion].levels[targetlevel].uptrans[nuparr[targetlevel]].stat_weight = stat_weight(element,ion,level);
                nuparr[targetlevel] += 1;
              }
            }
          }
          //printf("A %g\n",elements[element].ions[ion].levels[level].transitions[i].einstein_A );
          //printout("%d -> %d has A %g\n",level,level-i-1,elements[element].ions[ion].levels[level].transitions[i].einstein_A );


          free(nuparr);
          free(ndownarr);

          for (level = 0; level < nlevelsmax; level++)
          {
            totaldowntrans += elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
            totaluptrans += elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
            free(transitions[level].to);
          }
          free(transitiontable);
          free(transitions);


          /// Also the phixslist
          if (ion < nions-1)
          {
//            elements[element].ions[ion].nbfcontinua = elements[element].ions[ion].ionisinglevels;//nlevelsmax;
            nbfcheck += elements[element].ions[ion].ionisinglevels; //nlevelsmax;
/*            if ((elements[element].ions[ion].phixslist = calloc(nlevelsmax, sizeof(ionsphixslist_t))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize phixslist for element %d, ion %d ... abort\n",element,ion);
              exit(0);
            }*/
            nbfcontinua += get_bfcontinua(element,ion);//elements[element].ions[ion].ionisinglevels;//nlevelsmax;
          }
        }
        else
        {
          printout("corrupted atomic data\n");
          exit(0);
        }
      }
    }
    fclose(adata);
    fclose(transitiondata);
    fclose(compositiondata);
    printout("nbfcheck %d\n",nbfcheck);
    printout("heatingcheck %d\n",heatingcheck);

    /// Save the linecounters value to the global variable containing the number of lines
    nlines = lineindex;
    printout("nlines %d\n",nlines);
    if (nlines > 0)
    {
      /// and release empty memory from the linelist
      if ((linelist = realloc(linelist, nlines*sizeof(linelist_entry))) == NULL)
      {
        printout("[fatal] input: not enough memory to reallocate linelist ... abort\n");
        exit(0);
      }
    }

    if (T_preset > 0) exit(0);


    /// Set up the list of allowed upward transitions for each level
    printout("total uptrans %d\n",totaluptrans);
    printout("total downtrans %d\n",totaldowntrans);
    printout("coolingcheck %d\n",coolingcheck);



    ///debug output
    /*
    FILE *linelist_file;
    if ((linelist_file = fopen("linelist_unsorted.out", "w")) == NULL){
      printout("Cannot open linelist.out.\n");
      exit(0);
    }
    for (i = 0; i < nlines; i++)
      fprintf(linelist_file,"element %d, ion %d, ionstage %d, upperindex %d, lowerindex %d, nu %g\n",linelist[i].elementindex, linelist[i].ionindex, elements[linelist[i].elementindex].ions[linelist[i].ionindex].ionstage, linelist[i].upperlevelindex, linelist[i].lowerlevelindex, linelist[i].nu);
    fclose(linelist_file);
    //exit(0);
    */

    ///then sort the linelist by decreasing frequency
    qsort(linelist,nlines,sizeof(linelist_entry),compare_linelistentry);


    /// Save sorted linelist into a file
    if (rank_global == 0)
    {
      if ((linelist_file = fopen("linelist.dat", "w")) == NULL){
        printout("Cannot open linelist.out.\n");
        exit(0);
      }
      fprintf(linelist_file,"%d\n",nlines);
      for (i = 0; i < nlines; i++)
      {
        fprintf(linelist_file,"%d %d %d %d %d %lg %lg %lg\n",i,linelist[i].elementindex, linelist[i].ionindex, linelist[i].upperlevelindex, linelist[i].lowerlevelindex, linelist[i].nu, linelist[i].einstein_A, linelist[i].osc_strength);
      }
      fclose(linelist_file);
      //exit(0);
    }



    ///Establish connection between transitions and sorted linelist
    //printout("[debug] init line counter list\n");
    printout("establish connection between transitions and sorted linelist\n");
    for(i=0; i < nlines; i++)
    {
      element = linelist[i].elementindex;
      ion = linelist[i].ionindex;
      lowerlevel = linelist[i].lowerlevelindex;
      upperlevel = linelist[i].upperlevelindex;
      for (ii=1; ii <= elements[element].ions[ion].levels[upperlevel].downtrans[0].targetlevel; ii++)
      {
        if (elements[element].ions[ion].levels[upperlevel].downtrans[ii].targetlevel == lowerlevel)
        {
          elements[element].ions[ion].levels[upperlevel].downtrans[ii].lineindex = i;
        }
      }
      for (ii=1; ii <= elements[element].ions[ion].levels[lowerlevel].uptrans[0].targetlevel; ii++)
      {
        if (elements[element].ions[ion].levels[lowerlevel].uptrans[ii].targetlevel == upperlevel)
        {
          elements[element].ions[ion].levels[lowerlevel].uptrans[ii].lineindex = i;
        }
      }
    }



    /// Photoionisation cross-sections
    ///======================================================
    ///finally read in photoionisation cross sections and store them to the atomic data structure
    printout("readin phixs datat\n");
    if ((phixsdata = fopen("phixsdata.txt", "r")) == NULL)
    {
      printout("Cannot open phixsdata.txt.\n");
      exit(0);
    }
    while (fscanf(phixsdata,"%d %d %d %d %d %d ",&Z,&upperion,&upperlevel,&lowerion,&lowerlevel,&npoints) != EOF)
    {
      //printout("[debug] Z %d, upperion %d, upperlevel %d, lowerion %d, lowerlevel, %d npoints %d\n",Z,upperion,upperlevel,lowerion,lowerlevel,npoints);
      /// translate readin anumber to element index
      element = get_elementindex(Z);

      /// store only photoionization crosssections for elements which are part of the current model atom
      if (element >= 0)
      {
        /// translate readin ionstages to ion indices
        //printout("[debug] element %d, lowermost_ionstage %d\n",element,elements[element].ions[0].ionstage);
        lowermost_ionstage = elements[element].ions[0].ionstage;
        upperion -= lowermost_ionstage;
        upperlevel -= 1;
        lowerion -= lowermost_ionstage;
        lowerlevel -= 1;
        /// store only photoionization crosssections for ions which are part of the current model atom
        /// for limited model atoms we further have to make sure that the lowerlevel is inside the limited model atom
        if (lowerion >= 0 && lowerlevel < get_nlevels(element,lowerion) && upperion < get_nions(element))
        {
          /// so far only photoionisations to upperions groundlevel are considered
          //sprintf(filename1,"database_%d_%d_%d",element,lowerion,lowerlevel);
          //sprintf(filename2,"interpol_%d_%d_%d",element,lowerion,lowerlevel);
          //if ((database_file = fopen(filename1, "w")) == NULL)
          //{
          // printf("Cannot open %s.\n",filename);
          // exit(0);
          //}
          //if ((interpol_file = fopen(filename2, "w")) == NULL)
          //{
          // printf("Cannot open %s.\n",filename);
          // exit(0);
          //}

          if (upperlevel == 0)
          {
            nu_edge = (epsilon(element,upperion,upperlevel) - epsilon(element,lowerion,lowerlevel))/H;
            //elements[element].ions[lowerion].levels[lowerlevel].photoion_xs_nu_edge = nu_edge;

            if ((nutable = calloc(npoints,sizeof(double))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize nutable... abort\n");
              exit(0);
            }
            if ((phixstable = calloc(npoints,sizeof(double))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize phixstable... abort\n");
              exit(0);
            }
            for (i = 0; i < npoints; i++)
            {
              fscanf(phixsdata,"%lg %lg",&energy,&phixs);
              nutable[i] = nu_edge + (energy*13.6*EV)/H;
              ///the photoionisation cross-sections in the database are given in Mbarn=1e6 * 1e-28m^2
              ///to convert to cgs units multiply by 1e-18
              phixstable[i] = phixs*1e-18;
              //fprintf(database_file,"%g %g\n", nutable[i], phixstable[i]);
            }
            nu_max = nutable[npoints-1];

            /// Now interpolate these cross-sections
            if ((elements[element].ions[lowerion].levels[lowerlevel].photoion_xs = calloc(NPHIXSPOINTS,sizeof(float))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize photoionxslist... abort\n");
              exit(0);
            }
            elements[element].ions[lowerion].levels[lowerlevel].photoion_xs[0] = phixstable[0];
            //fprintf(interpol_file,"%g %g\n", nu_edge, phixstable[0]);

            gsl_interp_accel *acc = gsl_interp_accel_alloc();
            gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear,npoints);
            gsl_spline_init(spline,nutable,phixstable,npoints);
            nu = nu_edge;
            for (i = 1; i < NPHIXSPOINTS; i++)
            {
              nu += 0.1*nu_edge;
              if (nu > nu_max)
              {
                phixs = phixstable[npoints-1] * pow(nu_max/nu,3);
                elements[element].ions[lowerion].levels[lowerlevel].photoion_xs[i] = phixs;
                //fprintf(interpol_file,"%g %g\n", nu, phixs);
              }
              else
              {
                phixs = gsl_spline_eval(spline,nu,acc);
                elements[element].ions[lowerion].levels[lowerlevel].photoion_xs[i] = phixs;
                //fprintf(interpol_file,"%g %g\n", nu, phixs);
              }
            }
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            free(nutable);
            free(phixstable);

            //nbfcontinua += 1;
            //printout("[debug] element %d, ion %d, level %d: phixs exists %g\n",element,lowerion,lowerlevel,phixs*1e-18);
          }
          /// for excited levels of upperion just jump over the next npoints lines
          else
          {
            for (i = 0; i < npoints; i++)
              fscanf(phixsdata,"%lg %lg",&energy,&phixs);
          }
          //fclose(database_file);
          //fclose(interpol_file);
        }
        /// for ions which are not part of the current model atom proceed through the list
        else
        {
          for (i = 0; i < npoints; i++)
            fscanf(phixsdata,"%lg %lg",&energy,&phixs);
        }

      }
      /// for elements which are not part of the current model atom proceed through the list
      else
      {
        for (i = 0; i < npoints; i++)
          fscanf(phixsdata,"%lg %lg",&energy,&phixs);
      }
    }
    fclose(phixsdata);



    printout("cont_index %d\n",cont_index);
    /// Now write the model atom to file for reuse
    if (rank_global == 0)
    {
      if ((modelatom = fopen("modelatom.dat", "w")) == NULL)
      {
        printout("Cannot write to modelatom.dat\n");
        exit(0);
      }
      int metastable,linelistindex,targetlevel;
      float einstein_A,oscillator_strength;

      fprintf(modelatom,"%d\n",nelements);
      fprintf(modelatom,"%d\n",homogeneous_abundances);
      for (element = 0; element < nelements; element++)
      {
        nions = get_nions(element);
        fprintf(modelatom,"%d %d %g %g\n",elements[element].anumber,nions,elements[element].abundance,elements[element].mass);
        for (ion = 0; ion < nions; ion++)
        {
          nlevels = get_nlevels(element,ion);
          ionisinglevels = elements[element].ions[ion].ionisinglevels;
          ionstage = elements[element].ions[ion].ionstage;
          ionpot = elements[element].ions[ion].ionpot;
//          nbfcont = elements[element].ions[ion].nbfcontinua;
          fprintf(modelatom,"%d %d %d %lg\n",nlevels,ionisinglevels,ionstage,ionpot);
          for (level = 0; level < nlevels; level++)
          {
            levelenergy = elements[element].ions[ion].levels[level].epsilon;
            statweight = elements[element].ions[ion].levels[level].stat_weight;
            cont_index = elements[element].ions[ion].levels[level].cont_index;
            metastable = elements[element].ions[ion].levels[level].metastable;
            fprintf(modelatom,"%.16e %lg %d %d\n",levelenergy,statweight,cont_index,metastable);

/*            for (i=0; i < level; i++)
            {
              einstein_A = elements[element].ions[ion].levels[level].transitions[i].einstein_A;
              oscillator_strength = elements[element].ions[ion].levels[level].transitions[i].oscillator_strength;
              linelistindex = elements[element].ions[ion].levels[level].transitions[i].linelistindex;
              fprintf(modelatom,"%g %g %d\n",einstein_A,oscillator_strength,linelistindex);
            }*/

            ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
            fprintf(modelatom,"%d\n",ndowntrans);
            for (i=1; i <= ndowntrans; i++)
            {
              targetlevel = elements[element].ions[ion].levels[level].downtrans[i].targetlevel;
              levelenergy = elements[element].ions[ion].levels[level].downtrans[i].epsilon;
              statweight = elements[element].ions[ion].levels[level].downtrans[i].stat_weight;
              lineindex = elements[element].ions[ion].levels[level].downtrans[i].lineindex;
              fprintf(modelatom,"%d %.16e %lg %d\n",targetlevel,levelenergy,statweight,lineindex);
            }

            nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
            fprintf(modelatom,"%d\n",nuptrans);
            for (i=1; i <= nuptrans; i++)
            {
              targetlevel = elements[element].ions[ion].levels[level].uptrans[i].targetlevel;
              levelenergy = elements[element].ions[ion].levels[level].uptrans[i].epsilon;
              statweight = elements[element].ions[ion].levels[level].uptrans[i].stat_weight;
              lineindex = elements[element].ions[ion].levels[level].uptrans[i].lineindex;
              fprintf(modelatom,"%d %.16e %lg %d\n",targetlevel,levelenergy,statweight,lineindex);
            }

            if (ion < nions-1)
            {
              if (level < ionisinglevels)
              {
                for (i = 0; i < NPHIXSPOINTS; i++)
                {
                  phixs = elements[element].ions[ion].levels[level].photoion_xs[i];
                  fprintf(modelatom,"%g\n",phixs);
                }
              }
            }
          }
        }
      }
      fclose(modelatom);
    }


  }

  /// Preproccessed model atom available read that in
  else
  {
    printout("[read_atomicdata] Load preprocessed model atom and linelist ...\n");
    printout("[read_atomicdata] Be aware that no consistency checks are done!!!\n");
    int metastable,linelistindex,targetlevel,dum,anumber;
    float einstein_A,oscillator_strength,abundance,mass;

    fscanf(modelatom,"%d",&nelements);
    fscanf(modelatom,"%d",&homogeneous_abundances);

    printout("nelements %d\n",nelements);

    if ((elements = (elementlist_entry *) malloc(nelements*sizeof(elementlist_entry))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize elementlist ... abort\n");
      exit(0);
    }
    printout("memory for elementlist allocated!\n");
    for (element = 0; element < nelements; element++)
    {
      printout("read in element %d\n",element);
      fscanf(modelatom,"%d %d %g %g\n",&anumber,&nions,&abundance,&mass);
      elements[element].anumber = anumber;
      elements[element].nions = nions;
      elements[element].abundance = abundance;
      elements[element].mass = mass;
      includedions += nions;

      //printout("element %d, anumber %d, nions %d, abundance %g, mass %g\n",element,anumber,nions,abundance,mass);

      if ((elements[element].ions = (ionlist_entry *) malloc(nions*sizeof(ionlist_entry))) == NULL)
      {
          printout("[fatal] input: not enough memory to initialize ionlist ... abort\n");
          exit(0);
      }
      for (ion = 0; ion < nions; ion++)
      {
        fscanf(modelatom,"%d %d %d %lg\n",&nlevels,&ionisinglevels,&ionstage,&ionpot);
        elements[element].ions[ion].nlevels =  nlevels;
        elements[element].ions[ion].ionisinglevels =  ionisinglevels;
        elements[element].ions[ion].ionstage = ionstage;
        elements[element].ions[ion].ionpot = ionpot;
//        elements[element].ions[ion].nbfcontinua = nbfcont;

        nbfcont = get_bfcontinua(element,ion);
        if ((elements[element].ions[ion].levels = (levellist_entry *) malloc(nlevels*sizeof(levellist_entry))) == NULL)
        {
          printout("[fatal] input: not enough memory to initialize levelist of element %d, ion %d ... abort\n",element,ion);
          exit(0);
        }
        for (level = 0; level < nlevels; level++)
        {
          fscanf(modelatom,"%lg %lg %d %d\n",&levelenergy,&statweight,&cont_index,&metastable);
          elements[element].ions[ion].levels[level].epsilon = levelenergy;
          elements[element].ions[ion].levels[level].stat_weight = statweight;
          elements[element].ions[ion].levels[level].cont_index = cont_index;
          elements[element].ions[ion].levels[level].metastable = metastable;

          //printout("level %d, epsilon %g, stat_weight %g, metastable %d\n",level,levelenergy,statweight,metastable);

/*          if (level > 0)
          {
            if ((elements[element].ions[ion].levels[level].transitions = malloc(level*sizeof(transitionlist_entry))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize transitionlist ... abort\n");
              exit(0);
            }
            for (i=0; i < level; i++)
            {
              fscanf(modelatom,"%g %g %d\n",&einstein_A,&oscillator_strength,&linelistindex);
              elements[element].ions[ion].levels[level].transitions[i].einstein_A = einstein_A;
              elements[element].ions[ion].levels[level].transitions[i].oscillator_strength = oscillator_strength;
              elements[element].ions[ion].levels[level].transitions[i].linelistindex = linelistindex;
            }
          }*/

          fscanf(modelatom,"%d\n",&ndowntrans);
          //printout("ndowntrans %d\n",ndowntrans);
          totaldowntrans += ndowntrans;
          if ((elements[element].ions[ion].levels[level].downtrans = (permittedtransitionlist_entry *) malloc((ndowntrans+1)*sizeof(permittedtransitionlist_entry))) == NULL)
          {
            printout("[fatal] input: not enough memory to initialize downtrans list ... abort\n");
            exit(0);
          }
          elements[element].ions[ion].levels[level].downtrans[0].targetlevel = ndowntrans;
          for (i=1; i <= ndowntrans; i++)
          {
            fscanf(modelatom,"%d %lg %lg %d\n",&targetlevel,&levelenergy,&statweight,&lineindex);
            elements[element].ions[ion].levels[level].downtrans[i].targetlevel = targetlevel;
            elements[element].ions[ion].levels[level].downtrans[i].epsilon = levelenergy;
            elements[element].ions[ion].levels[level].downtrans[i].stat_weight = statweight;
            elements[element].ions[ion].levels[level].downtrans[i].lineindex = lineindex;
          }

          fscanf(modelatom,"%d\n",&nuptrans);
          //printout("nuptrans %d\n",nuptrans);
          totaluptrans += nuptrans;
          if ((elements[element].ions[ion].levels[level].uptrans = (permittedtransitionlist_entry *) malloc((nuptrans+1)*sizeof(permittedtransitionlist_entry))) == NULL)
          {
            printout("[fatal] input: not enough memory to initialize uptrans list ... abort\n");
            exit(0);
          }
          elements[element].ions[ion].levels[level].uptrans[0].targetlevel = nuptrans;
          for (i=1; i <= nuptrans; i++)
          {
            fscanf(modelatom,"%d %lg %lg %d\n",&targetlevel,&levelenergy,&statweight,&lineindex);
            elements[element].ions[ion].levels[level].uptrans[i].targetlevel = targetlevel;
            elements[element].ions[ion].levels[level].uptrans[i].epsilon = levelenergy;
            elements[element].ions[ion].levels[level].uptrans[i].stat_weight = statweight;
            elements[element].ions[ion].levels[level].uptrans[i].lineindex = lineindex;
          }

          /// Photoionisation cross-sections are only provided up to nions-1
          if (ion < nions-1)
          {
            /// and only available for levels below the ionisation threshold
            if (level < ionisinglevels)
            {
              /// we allocate only memory for as many levels as we want to use as bf-continua
              if (level < nbfcont)
              {
                if ((elements[element].ions[ion].levels[level].photoion_xs = (float *) malloc(NPHIXSPOINTS*sizeof(float))) == NULL)
                {
                  printout("[fatal] input: not enough memory to initialize photoionxslist... abort\n");
                  exit(0);
                }
              }
              /// read all the cross sections, but use only as many as bf-continua are used
              for (i = 0; i < NPHIXSPOINTS; i++)
              {
                fscanf(modelatom,"%lg\n",&phixs);
                if (level < nbfcont) elements[element].ions[ion].levels[level].photoion_xs[i] = phixs;
              }
            }

            if (level < nbfcont)
            {
              /// For those we just need to allocate memory
              /// and we allocate only memory for as many levels as we want to use as bf-continua
              if ((elements[element].ions[ion].levels[level].spontrecombcoeff = (double *) malloc(TABLESIZE*sizeof(double))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize spontrecombcoeff table for element %d, ion %d, level %d\n",element,ion,level);
                exit(0);
              }
              if ((elements[element].ions[ion].levels[level].corrphotoioncoeff = (double *) malloc(TABLESIZE*sizeof(double))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize photoioncoeff table for element %d, ion %d, level %d\n",element,ion,level);
                exit(0);
              }
              if ((elements[element].ions[ion].levels[level].bfheating_coeff = (double *) malloc(TABLESIZE*sizeof(double))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize modified_photoioncoeff table for element %d, ion %d, level %d\n",element,ion,level);
                exit(0);
              }
              if ((elements[element].ions[ion].levels[level].bfcooling_coeff = (double *) malloc(TABLESIZE*sizeof(double))) == NULL)
              {
                printout("[fatal] input: not enough memory to initialize modified_photoioncoeff table for element %d, ion %d, level %d\n",element,ion,level);
                exit(0);
              }
            }
          }

        }

        /// For those we just need to allocate memory
//         if ((elements[element].ions[ion].zeta = (float *) malloc(TABLESIZE*sizeof(float))) == NULL)
//         {
//           printout("[fatal] input: not enough memory to initialize zetalist for element %d, ion %d ... abort\n",element,ion);
//           exit(0);
//         }
        if ((elements[element].ions[ion].Alpha_sp = (float *) malloc(TABLESIZE*sizeof(float))) == NULL)
        {
          printout("[fatal] input: not enough memory to initialize Alpha_sp list for element %d, ion %d ... abort\n",element,ion);
          exit(0);
        }
        /// Bf continua are only available up to nions-1
        if (ion < nions-1)
        {
          /*if ((elements[element].ions[ion].phixslist = malloc(nlevels*sizeof(ionsphixslist_t))) == NULL)
          {
            printout("[fatal] input: not enough memory to initialize phixslist for element %d, ion %d ... abort\n",element,ion);
            exit(0);
          }*/
          //nbfcontinua += nlevels;
          nbfcontinua += get_bfcontinua(element,ion);//ionisinglevels;
        }
      }
    }
    fclose(modelatom);

    if ((linelist_file = fopen("linelist.dat", "r")) == NULL){
      printout("Cannot open linelist.out.\n");
      exit(0);
    }
    fscanf(linelist_file,"%d\n",&nlines);
    if ((linelist = (linelist_entry *) malloc(nlines*sizeof(linelist_entry))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialize linelist ... abort\n");
      exit(0);
    }

    for (i = 0; i < nlines; i++)
    {

        fscanf(linelist_file,"%d %d %d %d %d %lg %lg %lg\n",&dum,&element,&ion,&upperlevel,&lowerlevel,&nu,&A_ul,&f_ul);
        linelist[i].elementindex = element;
        linelist[i].ionindex = ion;
        linelist[i].upperlevelindex = upperlevel;
        linelist[i].lowerlevelindex = lowerlevel;
        linelist[i].nu = nu;
        linelist[i].einstein_A = A_ul;
        linelist[i].osc_strength = f_ul;
    }

    fclose(linelist_file);

  }



  printout("included ions %d\n",includedions);



  /// INITIALISE THE ABSORPTION/EMISSION COUNTERS ARRAYS
  ///======================================================
  #ifdef RECORD_LINESTAT
    if ((ecounter  = malloc(nlines*sizeof(int))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialise ecounter array ... abort\n");
      exit(0);
    }
    if ((acounter  = malloc(nlines*sizeof(int))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialise ecounter array ... abort\n");
      exit(0);
    }
    if ((linestat_reduced  = malloc(nlines*sizeof(int))) == NULL)
    {
      printout("[fatal] input: not enough memory to initialise ecounter array ... abort\n");
      exit(0);
    }
  #endif

  /// SET UP THE COOLING LIST
  ///======================================================
  /// Determine number of processes which allow kpkts to convert to something else.
  /// This number is given by the collisional excitations (so far determined from the oscillator strengths
  /// by the van Regemorter formula, therefore totaluptrans), the number of free-bound emissions and collisional ionisations
  /// (as long as we only deal with ionisation to the ground level this means for both of these
  /// \sum_{elements,ions}get_nlevels(element,ion) and free-free which is \sum_{elements} get_nions(element)-1
  /*ncoolingterms = totaluptrans;
  for (element = 0; element < nelements; element++)
  {
    nions = get_nions(element);
    for (ion=0; ion < nions; ion++)
    {
      if (get_ionstage(element,ion) > 1) ncoolingterms += 1;
      if (ion < nions - 1) ncoolingterms += 2 * get_ionisinglevels(element,ion);
    }
  }
  printout("[info] read_atomicdata: number of coolingterms %d\n",ncoolingterms);*/


  ncoolingterms = 0;
  for (element = 0; element < nelements; element++)
  {
    nions = get_nions(element);
    for (ion=0; ion < nions; ion++)
    {
      add = 0;
      elements[element].ions[ion].coolingoffset = ncoolingterms;
      /// Ionised ions add one ff-cooling term
      if (get_ionstage(element,ion) > 1) add += 1;
      /// Ionisinglevels below the closure ion add to bf and col ionisation
      //if (ion < nions - 1) add += 2 * get_ionisinglevels(element,ion);
      if (ion < nions - 1) add += 2 * get_bfcontinua(element,ion);
      /// All the levels add number of col excitations
      nlevels = get_nlevels(element,ion);
      for (level = 0; level < nlevels; level++)
      {
        add += elements[element].ions[ion].levels[level].uptrans[0].targetlevel;
      }
      elements[element].ions[ion].ncoolingterms = add;
      ncoolingterms += add;
    }
  }
  printout("[info] read_atomicdata: number of coolingterms %d\n",ncoolingterms);

  /// And setup the global coolinglist. A shorter coolinglist with only the important contributors
  /// is part of the cellhistory.
  /*
  if ((globalcoolinglist = malloc(ncoolingterms*sizeof(coolinglist_entry))) == NULL)
  {
    printout("[fatal] read_atomicdata: not enough memory to initialize globalcoolinglist ... abort\n");
    exit(0);
  }
  */



  /// SET UP THE CELL HISTORY
  ///======================================================
  /// Stack which holds information about population and other cell specific data
  /// ===> move to update_packets
  if ((cellhistory = (cellhistory_struct *) malloc(nthreads*sizeof(cellhistory_struct))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize cellhistory of size %d... abort\n",nthreads);
    exit(0);
  }
  #ifdef _OPENMP
    #pragma omp parallel private(element,ion,level,nions,nlevels,ndowntrans,nuptrans)
    {
  #endif
      printout("[info] input: initializing cellhistory for thread %d ...\n",tid);

      cellhistory[tid].cellnumber = -99;

      cellhistory[tid].coolinglist = (cellhistorycoolinglist_t *) malloc(ncoolingterms*sizeof(cellhistorycoolinglist_t));
      if (cellhistory[tid].coolinglist == NULL)
      {
        printout("[fatal] input: not enough memory to initialize cellhistory's coolinglist ... abort\n");
        exit(0);
      }

      if ((cellhistory[tid].chelements = (chelements_struct *) malloc(nelements*sizeof(chelements_struct))) == NULL)
      {
        printout("[fatal] input: not enough memory to initialize cellhistory's elementlist ... abort\n");
        exit(0);
      }
      for (element = 0; element < nelements; element++)
      {
        nions = get_nions(element);
        if ((cellhistory[tid].chelements[element].chions = (chions_struct *) malloc(nions*sizeof(chions_struct))) == NULL)
        {
          printout("[fatal] input: not enough memory to initialize cellhistory's ionlist ... abort\n");
          exit(0);
        }
        for (ion = 0; ion < nions; ion++)
        {
          nlevels = get_nlevels(element,ion);
          if ((cellhistory[tid].chelements[element].chions[ion].chlevels = (chlevels_struct *) malloc(nlevels*sizeof(chlevels_struct))) == NULL)
          {
            printout("[fatal] input: not enough memory to initialize cellhistory's levellist ... abort\n");
            exit(0);
          }
          for (level = 0; level < nlevels; level++)
          {
            ndowntrans = elements[element].ions[ion].levels[level].downtrans[0].targetlevel;
            nuptrans = elements[element].ions[ion].levels[level].uptrans[0].targetlevel;

            if ((cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc = (double *) malloc((ndowntrans+1)*sizeof(double))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize cellhistory's individ_rad_deexc ... abort\n");
              exit(0);
            }
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_rad_deexc[0] = ndowntrans;

            if ((cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same = (double *) malloc((ndowntrans+1)*sizeof(double))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize cellhistory's individ_internal_down_same ... abort\n");
              exit(0);
            }
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_down_same[0] = ndowntrans;

            if ((cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same = (double *) malloc((nuptrans+1)*sizeof(double))) == NULL)
            {
              printout("[fatal] input: not enough memory to initialize cellhistory's individ_internal_up_same ... abort\n");
              exit(0);
            }
            cellhistory[tid].chelements[element].chions[ion].chlevels[level].individ_internal_up_same[0] = nuptrans;
          }
        }
      }
  #ifdef _OPENMP
    }
  #endif



  /// Printout some information about the read-in model atom
  ///======================================================
  //includedions = 0;
  int includedlevels = 0;
  int includedionisinglevels = 0;
  printout("[input.c] this simulation contains\n");
  printout("----------------------------------\n");
  for (element = 0; element < nelements; element++)
  {
    printout("[input.c]   element Z = %d\n",get_element(element));
    nions = get_nions(element);
    //includedions += nions;
    for (ion=0; ion < nions; ion++)
    {
      printout("[input.c]     ion %d with %d levels (%d ionising)\n",get_ionstage(element,ion),get_nlevels(element,ion),get_ionisinglevels(element,ion));
      includedlevels += get_nlevels(element,ion);
      includedionisinglevels += get_ionisinglevels(element,ion);
    }
  }
  printout("[input.c]   in total %d ions, %d levels (%d ionising) and %d lines\n",includedions,includedlevels,includedionisinglevels,nlines);


  FILE *bflist_file;
  if ((bflist = (bflist_t *) malloc(includedionisinglevels*sizeof(bflist_t))) == NULL)
  {
    printout("[fatal] input: not enough memory to initialize bflist ... abort\n");
    exit(0);
  }
  if (rank_global == 0)
  {
    if ((bflist_file = fopen("bflist.dat", "w")) == NULL){
      printout("Cannot open bflist.dat.\n");
      exit(0);
    }
    fprintf(bflist_file,"%d\n",includedionisinglevels);
  }
  i=0;
  for (element = 0; element < nelements; element++)
  {
    nions = get_nions(element);
    for (ion = 0; ion < nions; ion++)
    {
      nlevels = get_ionisinglevels(element,ion);
      for (level = 0; level < nlevels; level++)
      {
        bflist[i].elementindex = element;
        bflist[i].ionindex = ion;
        bflist[i].levelindex = level;
        if (rank_global == 0) fprintf(bflist_file,"%d %d %d %d\n",i,element,ion,level);
        i++;
      }
    }
  }
  if (rank_global == 0) fclose(bflist_file);



  /// SET UP THE PHIXSLIST
  ///======================================================
  printout("[info] read_atomicdata: number of bfcontinua %d\n",nbfcontinua);
  nbfcontinua_ground = includedions-nelements;
  printout("[info] read_atomicdata: number of ground-level bfcontinua %d\n",nbfcontinua_ground);

  int itid;
  phixslist = (phixslist_t *) malloc(nthreads*sizeof(phixslist_t));
  if (phixslist == NULL)
  {
    printout("[fatal] read_atomicdata: not enough memory to initialize phixslist... abort\n");
    exit(0);
  }

  /// MK: 2012-01-19
  /// To fix the OpenMP problem on BlueGene machines this parallel section was removed and replaced by
  /// a serial loop which intializes the phixslist data structure for all threads in a loop. I'm still
  /// not sure why this causes a problem at all and on BlueGene architectures in particular. However,
  /// it seems to fix the problem.
  //#ifdef _OPENMP
  //  #pragma omp parallel private(i,element,ion,level,nions,nlevels,epsilon_upper,E_threshold,nu_edge)
  //  {
  //#endif
  for (itid = 0; itid < nthreads; itid++)
    {
      /// Number of ground level bf-continua equals the total number of included ions minus the number
      /// of included elements, because the uppermost ionisation stages can't ionise.
      //printout("groundphixslist nbfcontinua_ground %d\n",nbfcontinua_ground);
      printout("initialising groundphixslist for itid %d\n",itid);
      phixslist[itid].groundcont = (groundphixslist_t *) malloc(nbfcontinua_ground*sizeof(groundphixslist_t));
      if (phixslist[itid].groundcont == NULL)
      {
        printout("[fatal] read_atomicdata: not enough memory to initialize phixslist[%d].groundcont... abort\n",itid);
        exit(0);
      }

      i = 0;
      for (element = 0; element < nelements; element++)
      {
        nions = get_nions(element);
        for (ion = 0; ion < nions-1; ion++)
        {
          epsilon_upper = epsilon(element,ion+1,0);
          level = 0;
          E_threshold = epsilon_upper - epsilon(element,ion,level);
          nu_edge = E_threshold/H;
          phixslist[itid].groundcont[i].element = element;
          phixslist[itid].groundcont[i].ion = ion;
          phixslist[itid].groundcont[i].level = level;
          phixslist[itid].groundcont[i].nu_edge = nu_edge;
          //printout("phixslist.groundcont nbfcontinua_ground %d, i %d, element %d, ion %d, level %d, nu_edge %g\n",nbfcontinua_ground,i,element,ion,level,nu_edge);
          i += 1;
        }
      }
      qsort(phixslist[itid].groundcont,nbfcontinua_ground,sizeof(groundphixslist_t),compare_groundphixslistentry_bynuedge);


      //if (TAKE_N_BFCONTINUA >= 0) phixslist = malloc(includedions*TAKE_N_BFCONTINUA*sizeof(phixslist_t));
      //else
      phixslist[itid].allcont = (fullphixslist_t *) malloc(nbfcontinua*sizeof(fullphixslist_t));
      if (phixslist[itid].allcont == NULL)
      {
        printout("[fatal] read_atomicdata: not enough memory to initialize phixslist... abort\n");
        exit(0);
      }

      i = 0;
      for (element = 0; element < nelements; element++)
      {
        nions = get_nions(element);
        for (ion = 0; ion < nions-1; ion++)
        {
          epsilon_upper = epsilon(element,ion+1,0);
          nlevels = get_bfcontinua(element,ion);
          //nlevels = get_ionisinglevels(element,ion);
          ///// The following line reduces the number of bf-continua per ion
          //if (nlevels > TAKE_N_BFCONTINUA) nlevels = TAKE_N_BFCONTINUA;
          for (level = 0; level < nlevels; level++)
          {
            E_threshold = epsilon_upper - epsilon(element,ion,level);
            nu_edge = E_threshold/H;
            phixslist[itid].allcont[i].element = element;
            phixslist[itid].allcont[i].ion = ion;
            phixslist[itid].allcont[i].level = level;
            phixslist[itid].allcont[i].nu_edge = nu_edge;
            phixslist[itid].allcont[i].index_in_groundphixslist = search_groundphixslist(nu_edge,&index_in_groundlevelcontestimator,element,ion,level);
            if (itid == 0) elements[element].ions[ion].levels[level].closestgroundlevelcont = index_in_groundlevelcontestimator;
            i += 1;
          }
        }
      }
      //nbfcontinua = i;
      //printout("number of bf-continua reduced to %d\n",nbfcontinua);
      qsort(phixslist[itid].allcont,nbfcontinua,sizeof(fullphixslist_t),compare_phixslistentry_bynuedge);
    }
  //#ifdef _OPENMP
  //  }
  //#endif
}



///****************************************************************************
/// Subroutine to read in input parameters from input.txt.
void read_parameterfile(rank)
     int rank;
{
  FILE *input_file;
  double rr, z1, z2, x;
  float dum2, dum3, dum4;
  int dum1, dum5, n, i;
  unsigned long int zseed; /* rnum generator seed */
  unsigned long int xseed,pre_zseed;

  if ((input_file = fopen("input.txt", "r")) == NULL)
  {
    printout("Cannot open input.txt.\n");
    exit(0);
  }

  fscanf(input_file, "%d", &dum1);
  if (dum1 > 0)
  {
    pre_zseed = dum1; ///random number seed
  }
  else
  {
    xseed = time(NULL);
    pre_zseed = xseed;
    printout("[debug] random number seed was %d\n",pre_zseed);
  }

  #ifdef _OPENMP
    #pragma omp parallel private(x,zseed,n)
    {
/*      tid = omp_get_thread_num();
      nthreads = omp_get_num_threads();
      if (nthreads > MTHREADS)
      {
        printout("[Fatal] too many threads. Set MTHREADS (%d) > nthreads (%d). Abort.\n",MTHREADS,nthreads);
        exit(0);
      }
      if (tid == 0) printout("OpenMP parallelisation active with %d threads\n",nthreads);
  #else
      tid = 0;
      nthreads = 1;*/
  #endif
      /// For MPI parallelisation, the random seed is changed based on the rank of the process
      /// For OpenMP parallelisation rng is a threadprivate variable and the seed changed according
      /// to the thread-ID tid.
      zseed = pre_zseed + (13 * rank) + (17*tid);

      printout("rank %d: thread %d has zseed %d\n",rank,tid,zseed);
      /// start by setting up the randon number generator
      rng = gsl_rng_alloc (gsl_rng_ran3);
      gsl_rng_set ( rng, zseed);
      /// call it a few times to get it in motion.
      for (n = 0; n< 100; n++)
      {
        x = gsl_rng_uniform(rng);
        //printout("zrand %g\n", x);
      }
  #ifdef _OPENMP
    }
  #endif


  fscanf(input_file, "%d", &dum1); ///number of time steps
  ntstep = dum1;
  fscanf(input_file, "%d %d", &dum1, &dum5); ///number of start and end time step
  itstep = dum1;
  ftstep = dum5;
  fscanf(input_file, "%g %g", &dum2, &dum3); ///start and end times
  tmin = dum2 * DAY;
  tmax = dum3 * DAY;
  //tlimit = dum4 * DAY;


  fscanf(input_file, "%g %g", &dum2, &dum3);
  nusyn_min = dum2 * MEV / H; ///lowest frequency to synthesise
  nusyn_max = dum3 * MEV / H; ///highest frequecnt to synthesise
  dlognusyn = (log(nusyn_max) - log(nusyn_min))/NSYN;

  fscanf(input_file, "%d", &dum1); ///number of times for synthesis
  nsyn_time = dum1;
  fscanf(input_file, "%g %g", &dum2, &dum3);///start and end times for synthesis
  for (i = 0; i < nsyn_time; i++)
  {
    time_syn[i] = exp(log(dum2) + (dum3*i)) * DAY;
  }


  fscanf(input_file, "%d", &dum1); ///model type
  if (dum1 == 1)
  {
    model_type = RHO_1D_READ;
  }
  if (dum1 == 2)
  {
    model_type = RHO_2D_READ;
  }
  if (dum1 == 3)
  {
    model_type = RHO_3D_READ;
  }

  fscanf(input_file, "%d", &dum1); ///compute the r-light curve?
  if (dum1 == 1) ///lc no estimators
  {
    do_r_lc = 1;
    do_rlc_est = 0;
  }
  else if (dum1 == 2)/// lc case with thin cells
  {
    do_r_lc = 1;
    do_rlc_est = 1;
  }
  else if (dum1 == 3)/// lc case with thick cells
  {
    do_r_lc = 1;
    do_rlc_est = 2;
  }
  else if (dum1 == 4) /// gamma-ray heating case
  {
    do_r_lc = 1;
    do_rlc_est = 3;
  }
  else if (dum1 != 0)
  {
    printout("Unknown rlc mode. Abort.\n");
    exit(0);
  }


  fscanf(input_file, "%d", &dum1); ///number of iterations
  n_out_it = dum1;


  fscanf(input_file, "%g", &dum2); ///change speed of light?
  CLIGHT_PROP = dum2 * CLIGHT;

  fscanf(input_file, "%g", &dum2); ///use grey opacity for gammas?
  gamma_grey = dum2;


  fscanf(input_file, "%g %g %g", &dum2, &dum3, &dum4); ///components of syn_dir

  if ((rr =(dum2*dum2) + (dum3*dum3) + (dum4*dum4)) > 1.e-6)
  {
    syn_dir[0] = dum2 / sqrt( rr );
    syn_dir[1] = dum3 / sqrt( rr );
    syn_dir[2] = dum4 / sqrt( rr );
  }
  else
  {
    z1 = 1. - (2.*gsl_rng_uniform(rng));
    z2 = gsl_rng_uniform(rng) * 2.0 * PI;
    syn_dir[2] = z1;
    syn_dir[0] = sqrt( (1. - (z1*z1))) * cos(z2);
    syn_dir[1] = sqrt( (1. - (z1*z1))) * sin(z2);
  }

  /// ensure that this vector is normalised.

  fscanf(input_file, "%d", &dum1); ///opacity choice
  opacity_case = dum1;

  ///MK: begin
  fscanf(input_file, "%g", &dum2); ///free parameter for calculation of rho_crit
  rho_crit_para = dum2;
  printout("input: rho_crit_para %g\n",rho_crit_para);
  ///the calculation of rho_crit itself depends on the time, therfore it happens in grid_init and update_grid

  fscanf(input_file, "%d", &dum1); /// activate debug output for packet
  debug_packet =  dum1;            /// select a negative value to deactivate
  ///MK: end

  /// Do we start a new simulation or, continue another one?
  continue_simulation = 0;          /// Preselection is to start a new simulation
  fscanf(input_file, "%d", &dum1);
  if (dum1 == 1)
  {
    continue_simulation = 1;        /// Continue simulation if dum1 = 1
    printout("input: continue simulation\n");
  }

  /// Wavelength (in Angstroems) at which the parameterisation of the radiation field
  /// switches from the nebular approximation to LTE.
  fscanf(input_file, "%g", &dum2); ///free parameter for calculation of rho_crit
  nu_rfcut = CLIGHT/(dum2*1e-8);
  printout("input: nu_rfcut %g\n",nu_rfcut);


  /// Sets the number of initial LTE timesteps for NLTE runs
  fscanf(input_file, "%d", &n_lte_timesteps);
  #ifdef FORCE_LTE
    printout("input: this is a pure LTE run\n");
  #else
    printout("input: this is a NLTE run\n");
    printout("input: do the first %d timesteps in LTE\n",n_lte_timesteps);
  #endif

  /// Set up initial grey approximation?
  fscanf(input_file, "%lg %d", &cell_is_optically_thick, &n_grey_timesteps);
  printout("input: cells with Thomson optical depth > %g are treated in grey approximation for the first %d timesteps\n",cell_is_optically_thick,n_grey_timesteps);


  /// Limit the number of bf-continua
  fscanf(input_file, "%d", &max_bf_continua);
  if (max_bf_continua == -1)
  {
    printout("input: use all bf-continua\n");
    max_bf_continua = 1e6;
  }
  else
  {
    printout("input: use only %d bf-continua per ion\n",max_bf_continua);
  }

  /// The following parameters affect the DO_EXSPEC mode only /////////////////
  /// Read number of MPI tasks for exspec
  fscanf(input_file, "%d", &dum1);
  #ifdef DO_EXSPEC
    nprocs_exspec = dum1;
    printout("input: DO_EXSPEC ... extract spectra for %d MPI tasks\n",nprocs_exspec);
    printout("input: DO_EXSPEC ... and %d packets per task\n",npkts);
  #endif

  /// Extract line-of-sight dependent information of last emission for spectrum_res
  ///   if 1, yes
  ///   if 0, no
  /// Only affects runs with DO_EXSPEC. But then it needs huge amounts of memory!!!
  /// Thus be aware to reduce MNUBINS for this mode of operation!
  fscanf(input_file, "%d", &dum1);
  #ifdef DO_EXSPEC
    do_emission_res = dum1;
    if (do_emission_res == 1) printout("input: DO_EXSPEC ... extract los dependent emission information\n");
  #endif

  /// To reduce the work imbalance between different MPI tasks I introduced a diffusion
  /// for kpkts, since it turned out that this work imbalance was largely dominated
  /// by continuous collisional interactions. By introducing a diffusion time for kpkts
  /// this loop is broken. The following two parameters control this approximation.
  /// Parameter one (a float) gives the relative fraction of a time step which individual
  /// kpkts live. Parameter two (an int) gives the number of time steps for which we
  /// want to use this approximation
  fscanf(input_file, "%g %d", &kpktdiffusion_timescale, &n_kpktdiffusion_timesteps);
  printout("input: kpkts diffuse %g of a time step's length for the first %d time steps\n",kpktdiffusion_timescale, n_kpktdiffusion_timesteps);

  fclose(input_file);
}


///****************************************************************************
/// Subroutine to read in a 1-D model.
int read_1d_model()
{
  FILE *model_input;
  int dum1, n;
  float dum2, dum3, dum4, dum5, dum6, dum7, dum8;
  double mass_in_shell;

  if ((model_input = fopen("model.txt", "r")) == NULL){
    printout("Cannot open model.txt.\n");
    exit(0);
  }

  /* 1st read the number of data points in the table of input model. */
  fscanf(model_input, "%d", &dum1);
  npts_model = dum1;
  if (npts_model > MMODELGRID)
  {
    printout("Too many points in input model. Abort.\n");
    exit(0);
  }
  /* Now read the time (in days) at which the model is specified. */
  fscanf(model_input, "%g", &dum2);
  t_model = dum2 * DAY;

  /* Now read in the lines of the model. Each line has 5 entries: the
     cell number (integer) the velocity at outer boundary of cell (float),
     the mass density in the cell (float), the abundance of Ni56 by mass
     in the cell (float) and the total abundance of all Fe-grp elements
     in the cell (float). For now, the last number is recorded but never
     used. */

  for (n = 0; n < npts_model; n++)
  {
    fscanf(model_input, "%d %g %g %g %g %g %g %g", &dum1, &dum2,  &dum3, &dum4, &dum5, &dum6, &dum7, &dum8);
    vout_model[n] = dum2 * 1.e5;
    rho_model[n] = pow(10., dum3);
    ffegrp_model[n] = dum4;
    fni_model[n] = dum5;
    fco_model[n] = dum6;
    f52fe_model[n] = dum7;
    f48cr_model[n] = dum8;
  }

  fclose(model_input);

  /* Now compute mtot and mni56 for the model. */

  mtot = 0.0;
  mni56 = 0.0;
  mco56 = 0.0;
  mfe52 = 0.0;
  mcr48 = 0.0;
  mfeg = 0.0;


  /* do inner most cell 1st */
  mass_in_shell = rho_model[0] * (pow(vout_model[0],3.)) * 4. * PI * pow(t_model,3.) / 3.;
  mtot += mass_in_shell;
  mni56 += mass_in_shell * fni_model[0];
  mco56 += mass_in_shell * fco_model[0];
  mfe52 += mass_in_shell * f52fe_model[0];
  mcr48 += mass_in_shell * f48cr_model[0];
  mfeg += mass_in_shell * ffegrp_model[0];



  /* Now do the rest. */
  for (n = 1; n < npts_model; n++)
  {
    mass_in_shell = rho_model[n] * (pow(vout_model[n],3.) - pow(vout_model[n-1],3.)) * 4. * PI * pow(t_model,3.) / 3.;
    mtot += mass_in_shell;
    mni56 += mass_in_shell * fni_model[n];
    mco56 += mass_in_shell * fco_model[n];
    mfe52 += mass_in_shell * f52fe_model[n];
    mcr48 += mass_in_shell * f48cr_model[n];
    mfeg += mass_in_shell * ffegrp_model[n];
  }

  printout("Total mass: %g. Ni mass: %g. 56Co mass: %g. 52Fe mass: %g.  48Cr mass: %g. Fe-grp mass: %g.\n", mtot/MSUN, mni56/MSUN, mco56/MSUN, mfe52/MSUN, mcr48/MSUN, mfeg/MSUN);

  vmax = vout_model[npts_model-1];
  rmax = vmax * tmin;
  xmax = ymax = zmax = rmax;
  printout("rmax %g\n", rmax);


  return(0);
}


///****************************************************************************
///****************************************************************************
/// Subroutine to read in a 2-D model.
int read_2d_model()
{
  FILE *model_input;
  int dum1, dum12, n, n1;
  float dum2, dum3, dum4, dum5, dum6, dum7, dum8, dum9;
  double mass_in_shell;

  if ((model_input = fopen("model.txt", "r")) == NULL){
    printout("Cannot open model.txt.\n");
    exit(0);
  }

  /* 1st read the number of data points in the table of input model. */
  fscanf(model_input, "%d %d", &dum1, &dum12);
  ncoord1_model = dum1; //r (cylindrical polar)
  ncoord2_model = dum12;//z (cylindrical polar)

  npts_model = dum1 * dum12;
  if (npts_model > MMODELGRID)
  {
    printout("Too many points in input model. Abort.\n");
    exit(0);
  }
  /* Now read the time (in days) at which the model is specified. */
  fscanf(model_input, "%g", &dum2);
  t_model = dum2 * DAY;

  /* Now read in vmax (in cm/s) */
  fscanf(model_input, "%g", &dum2);
  vmax = dum2;
  dcoord1 = dum2 * t_model / ncoord1_model; //dr for input model
  dcoord2 = 2. * dum2 * t_model / ncoord2_model;//dz for input model


  /* Now read in the model. Each point in the model has two lines of input.
     First is an index for the cell then its r-mid point then its z-mid point
     then its total mass density.
     Second is the total FeG mass, initial 56Ni mass, initial 56Co mass */

  for (n = 0; n < npts_model; n++)
  {
    fscanf(model_input, "%d %g %g %g", &dum1, &dum2,  &dum3, &dum4);
    fscanf(model_input, "%g %g %g %g %g",&dum5, &dum6, &dum7, &dum8, &dum9);
    rho_model[n] = dum4;
    ffegrp_model[n] = dum5;
    fni_model[n] = dum6;
    fco_model[n] = dum7;
    f52fe_model[n] = dum8;
    f48cr_model[n] = dum9;
  }

  fclose(model_input);

  /* Now compute mtot and mni56 for the model. */

  mtot = 0.0;
  mni56 = 0.0;
  mco56 = 0.0;
  mfe52 = 0.0;
  mcr48 = 0.0;
  mfeg = 0.0;

  n1=0;
  /* Now do the rest. */
  for (n = 0; n < npts_model; n++)
  {
    mass_in_shell = rho_model[n] * ((2*n1) + 1) * PI * dcoord2 * pow(dcoord1,2.);
    mtot += mass_in_shell;
    mni56 += mass_in_shell * fni_model[n];
    mco56 += mass_in_shell * fco_model[n];
    mfe52 += mass_in_shell * f52fe_model[n];
    mcr48 += mass_in_shell * f48cr_model[n];
    mfeg += mass_in_shell * ffegrp_model[n];
    n1++;
    if (n1 == ncoord1_model)
      {
	n1=0;
      }
  }

  printout("Total mass: %g. Ni mass: %g. 56Co mass: %g. 52Fe mass: %g.  48Cr mass: %g. Fe-grp mass: %g.\n", mtot/MSUN, mni56/MSUN, mco56/MSUN, mfe52/MSUN, mcr48/MSUN, mfeg/MSUN);

  rmax = vmax * tmin;
  xmax = ymax = zmax = rmax;
  printout("rmax %g\n", rmax);


  return(0);
}
///****************************************************************************
/// Subroutine to read in a 3-D model.
int read_3d_model()
{
  void allocate_compositiondata(int cellnumber);
  void allocate_cooling(int modelgridindex);
  FILE *model_input;
  int dum1, n;
  float dum2, dum3, dum4, dum5, dum6;
  float rho_model;
  double mass_in_shell, helper;
  int mgi;

  if ((model_input = fopen("model.txt", "r")) == NULL)
  {
    printout("Cannot open model.txt.\n");
    exit(0);
  }

  /// 1st read the number of data points in the table of input model.
  /// This MUST be the same number as the maximum number of points used in the grid - if not, abort.
  fscanf(model_input, "%d", &dum1);
  npts_model = dum1;
  if (npts_model > MMODELGRID)
  {
    printout("Too many points in input model. Abort.\n");
    exit(0);
  }
  if (npts_model != nxgrid * nygrid * nzgrid)
  {
    printout("3D model/grid mismatch. Abort.\n");
    exit(0);
  }

  /// Now read the time (in days) at which the model is specified.
  fscanf(model_input, "%g", &dum2);
  t_model = dum2 * DAY;

  /// Now read in vmax for the model (in cm s^-1).
  fscanf(model_input, "%g", &dum2);
  vmax = dum2;



  /// Now read in the lines of the model.
  min_den=1.e99;


  /*mgi is the index to the model grid - empty cells are sent to MMODELGRID, otherwise each input cell is one modelgrid cell */

  mgi=0;
  for (n = 0; n < npts_model; n++)
  {
    fscanf(model_input, "%d %g %g %g %g", &dum1, &dum3,  &dum4, &dum5, &rho_model);
    //printout("cell %d, vz %g, vy %g, vx %g, rho %g, rho_init %g\n",dum1,dum3,dum4,dum5,rho_model,rho_model* pow( (t_model/tmin), 3.));
    if (rho_model < 0)
    {
      printout("negative input density %g %d\n", rho_model, n);
      exit(0);
    }

    if (mgi > MMODELGRID-1)
    {
      printout("3D model wants more modelgrid cells than MMODELGRID. Abort.\n");
      exit(0);
    }

    if (rho_model > 0)
    {
      modelgrid[mgi].associated_cells = 1;
      cell[n].modelgridindex = mgi;
      helper=rho_model * pow( (t_model/tmin), 3.);
      //printout("mgi %d, helper %g\n",mgi,helper);
      set_rhoinit(mgi,helper);
      //printout("mgi %d, rho_init %g\n",mgi,get_rhoinit(mgi));
      set_rho(mgi,helper);


      if (min_den > rho_model)
      {
        min_den = rho_model;
      }
    }
    else
    {
      cell[n].modelgridindex = MMODELGRID;
    }

    fscanf(model_input, "%g %g %g %g %g", &dum2, &dum3,  &dum4, &dum5, &dum6);
    //printout("ffe %g, fni %g, fco %g, ffe52 %g, fcr48 %g\n",dum2,dum3,dum4,dum5,dum6);
    if (rho_model > 0)
      {
	set_ffe(mgi,dum2);
	set_fni(mgi,dum3);
	set_fco(mgi,dum4);
	set_f52fe(mgi,dum5);
	set_f48cr(mgi,dum6);
	allocate_compositiondata(mgi);
	allocate_cooling(mgi);
        //printout("mgi %d, control rho_init %g\n",mgi,get_rhoinit(mgi));
	mgi++;
      }
  }

  printout("min_den %g\n", min_den);
  printout("Effectively used model grid cells %d\n", mgi);

  /// Now, reduce the actual size of the modelgrid to the number of non-empty cells.
  /// Actually this doesn't reduce the memory consumption since all the related
  /// arrays are allocated statically at compile time with size MMODELGRID+1.
  /// However, it ensures that update_grid causes no idle tasks due to empty cells!
  npts_model = mgi;

  fclose(model_input);

  set_rhoinit(MMODELGRID,0.);
  set_rho(MMODELGRID,0.);
  set_nne(MMODELGRID,0.);
  set_ffe(MMODELGRID,0.);
  set_fco(MMODELGRID,0.);
  set_fni(MMODELGRID,0.);
  set_f52fe(MMODELGRID,0.);
  set_f48cr(MMODELGRID,0.);
  set_Te(MMODELGRID,MINTEMP);
  set_TJ(MMODELGRID,MINTEMP);
  set_TR(MMODELGRID,MINTEMP);
  allocate_compositiondata(MMODELGRID);

  /// Now compute mtot and mni56 for the model.
  /// Assumes cells are cubes here - all same volume.

  mtot = 0.0;
  mni56 = 0.0;
  mco56 = 0.0;
  mfe52 = 0.0;
  mcr48 = 0.0;
  mfeg = 0.0;

  //for (n = 0; n < npts_model; n++)
  for (mgi = 0; mgi < npts_model; mgi++)
  {
    //mgi = cell[n].modelgridindex;
    mass_in_shell = get_rhoinit(mgi);
    //printout("n %d, mgi %d, rho_init %g\n",n,mgi,mass_in_shell);
    mtot += mass_in_shell;
    mni56 += mass_in_shell * get_fni(mgi);
    mco56 += mass_in_shell * get_fco(mgi);
    mfe52 += mass_in_shell * get_f52fe(mgi);
    mcr48 += mass_in_shell * get_f48cr(mgi);
    mfeg += mass_in_shell * get_ffe(mgi);
  }

  double cellvolume = pow((2 * vmax * tmin),3.) / (nxgrid*nygrid*nzgrid);
  mtot = mtot * cellvolume;
  mni56 = mni56 * cellvolume;
  mco56 = mco56 * cellvolume;
  mfe52 = mfe52 * cellvolume;
  mcr48 = mcr48 * cellvolume;
  mfeg = mfeg  * cellvolume;

  printout("Total mass: %g. Ni mass: %g. 56Co mass: %g. 52Fe mass: %g.  48Cr mass: %g. Fe-grp mass: %g.\n", mtot/MSUN, mni56/MSUN, mco56/MSUN, mfe52/MSUN, mcr48/MSUN, mfeg/MSUN);

  rmax = vmax * tmin;
  xmax = ymax = zmax = rmax;
  printout("rmax %g\n", rmax);


  ///Opacity treatment moved to grid_init & update_grid
//    if (opacity_case == 0)
//      {
//        for (n = 0; n < npts_model; n++)
// 	 {
// 	   cell[n].kappa_grey = GREY_OP;
// 	 }
//      }
//    else if (opacity_case == 1)
//      {
//        for (n = 0; n < npts_model; n++)
// 	 {
// 	   cell[n].kappa_grey = ((0.9 * cell[n].f_fe) + 0.1) * GREY_OP / ((0.9 *  mfeg / mtot) + 0.1);
// 	 }
//      }
//    else if (opacity_case == 2)
//      {
//        for (n = 0; n < npts_model; n++)
// 	 {
// //           cell[n].kappa_grey = ((0.9 * cell[n].f_fe) + 0.1)/cell[n].rho_init * GREY_OP*rhotot / ((0.9 *  mfeg / mtot) + 0.1);
// //           cell[n].kappa_grey = ((0.9 * cell[n].f_fe) + 0.1)/cell[n].rho_init * GREY_OP*rhom / ((0.9 *  mfeg) + (0.1 * mtot));
//            cell[n].kappa_grey = ((0.9 * cell[n].f_fe) + 0.1)/cell[n].rho_init * GREY_OP*rho_sum / ((0.9 *  ni56_sum) + (0.1 * npts_model));
// 	 }
//      }
//    else
//      {
//        printout("Unknown opacity case. Abort.\n");
//        exit(0);
//      }

   return(0);
}


///****************************************************************************
/// Helper function to sort the linelist by frequency.

int compare_linelistentry(const void *p1, const void *p2)
{
  linelist_entry *a1, *a2;
  a1 = (linelist_entry *)(p1);
  a2 = (linelist_entry *)(p2);
  //printf("%d %d %d %d %g\n",a1->elementindex,a1->ionindex,a1->lowerlevelindex,a1->upperlevelindex,a1->nu);
  //printf("%d %d %d %d %g\n",a2->elementindex,a2->ionindex,a2->lowerlevelindex,a2->upperlevelindex,a2->nu);
  //printf("%g\n",a2->nu - a1->nu);
  if (fabs(a2->nu - a1->nu) < (1.e-10 * a1->nu))
  {
    a2->nu = a1->nu;
    if (a1->lowerlevelindex > a2->lowerlevelindex)
    {
      return -1;
    }
    else if (a1->lowerlevelindex < a2->lowerlevelindex)
    {
      return 1;
    }
    else if (a1->upperlevelindex > a2->upperlevelindex)
    {
      return -1;
    }
    else if (a1->upperlevelindex < a2->upperlevelindex)
    {
      return 1;
    }
    else
    {
      printout("Duplicate atomic line?\n");
      printout("%d %d %d %d %g\n",a1->elementindex,a1->ionindex,a1->lowerlevelindex,a1->upperlevelindex,a1->nu);
      printout("%d %d %d %d %g\n",a2->elementindex,a2->ionindex,a2->lowerlevelindex,a2->upperlevelindex,a2->nu);
      return 0;
      //exit(0);
    }
  }
  else
  {
    if ((a1->nu < a2->nu) || (a1->nu == a2->nu))
      return 1;
    else if (a1->nu > a2->nu)
      return -1;
    else
      return 0;
  }
}

/*
int compare_linelistentry(const void *p1, const void *p2)
{
  linelist_entry *a1, *a2;
  a1 = (linelist_entry *)(p1);
  a2 = (linelist_entry *)(p2);
  //printf("%d %d %d %d %g\n",a1->elementindex,a1->ionindex,a1->lowerlevelindex,a1->upperlevelindex,a1->nu);
  //printf("%d %d %d %d %g\n",a2->elementindex,a2->ionindex,a2->lowerlevelindex,a2->upperlevelindex,a2->nu);
  //printf("%g\n",a2->nu - a1->nu);
  if (a1->nu - a2->nu < 0)
    return 1;
  else if (a1->nu - a2->nu > 0)
    return -1;
  else
    return 0;
}
*/


///****************************************************************************
int search_groundphixslist(double nu_edge, int *index_in_groundlevelcontestimator, int el, int in, int ll)
/// Return the closest ground level continuum index to the given edge
/// frequency. If the given edge frequency is redder than the reddest
/// continuum return -1.
/// NB: groundphixslist must be in ascending order.
{
  double epsilon(int element, int ion, int level);
  int i,index,element,ion,level,looplevels;
  double left,right;

  if (nu_edge < phixslist[tid].groundcont[0].nu_edge)
  {
    index = -1;
    *index_in_groundlevelcontestimator = -1;
  }
  else
  {
    for (i = 1; i < nbfcontinua_ground; i++)
    {
      if (nu_edge < phixslist[tid].groundcont[i].nu_edge) break;
    }
/*    if (i == nbfcontinua_ground)
    {
      printout("[fatal] search_groundphixslist: i %d, nu_edge %g, phixslist[tid].groundcont[i-1].nu_egde %g ... abort\n",i,nu_edge,phixslist[tid].groundcont[i-1].nu_edge);
      printout("[fatal] search_groundphixslist: this is element %d, ion %d, level %d in groundphixslist at i-1\n",el,in,ll);
      //printout("[fatal] search_groundphixslist: this is element %d, ion %d, level %d in groundphixslist at i-1\n",phixslist[tid].groundcont[i-1].element,phixslist[tid].groundcont[i-1].ion,groundphixslist[i-1].level);
      abort();
    }*/
    if (i == nbfcontinua_ground)
    {
      element = phixslist[tid].groundcont[i-1].element;
      ion = phixslist[tid].groundcont[i-1].ion;
      level = phixslist[tid].groundcont[i-1].level;
      if (element == el && ion == in && level == ll)
      {
        index = i-1;
      }
      else
      {
	printout("[fatal] search_groundphixslist: element %d, ion %d, level %d has edge_frequency %g equal to the bluest ground-level continuum\n",el,in,ll,nu_edge);
	printout("[fatal] search_groundphixslist: bluest ground level continuum is element %d, ion %d, level %d at nu_edge %g\n",element,ion,level,phixslist[tid].groundcont[i-1].nu_edge);
	printout("[fatal] search_groundphixslist: i %d, nbfcontinua_ground %d\n",i,nbfcontinua_ground);
	printout("[fatal] This shouldn't happen, is hoewever possible if there are multiple levels in the adata file at energy=0\n");
	for (looplevels=0; looplevels < get_nlevels(el,in); looplevels++)
	{
	  printout("[fatal]   element %d, ion %d, level %d, energy %g\n",el,in,looplevels,epsilon(el,in,looplevels));
	}
	printout("[fatal] Abort omitted ... MAKE SURE ATOMIC DATA ARE CONSISTENT\n");
	index = i-1;
	//abort();
      }
    }
    else
    {
      left = nu_edge - phixslist[tid].groundcont[i-1].nu_edge;
      right = phixslist[tid].groundcont[i].nu_edge - nu_edge;
      if (left <= right)
      {
        element = phixslist[tid].groundcont[i-1].element;
        ion = phixslist[tid].groundcont[i-1].ion;
        index = i-1;
      }
      else
      {
        element = phixslist[tid].groundcont[i].element;
        ion = phixslist[tid].groundcont[i].ion;
        index = i;
      }
    }
    *index_in_groundlevelcontestimator = element*maxion+ion;
  }

  return index;
}


















///****************************************************************************
/// Subroutine to read in input parameters from input.txt.
void update_parameterfile(int nts)
{
  FILE *input_file;
  float dum2, dum3, dum4;
  int dum1, dum5, n, i;



  if ((input_file = fopen("input.txt", "r+")) == NULL)
  {
    printout("Cannot open input.txt.\n");
    exit(0);
  }
  //setvbuf(input_file, NULL, _IOLBF, 0);

  fscanf(input_file, "%d\n", &dum1); /// Random number seed
  fscanf(input_file, "%d\n", &dum1); /// Number of time steps

  /// Number of start and end time step, update start time step
  fseek(input_file,0,SEEK_CUR);
  fprintf(input_file, "%3.3d %3.3d\n", nts, ftstep);
  fseek(input_file,0,SEEK_CUR);

  fscanf(input_file, "%g %g\n", &dum2, &dum3); ///start and end times
  fscanf(input_file, "%g %g\n", &dum2, &dum3); ///lowest and highest frequency to synthesise in MeV
  fscanf(input_file, "%d\n", &dum1); ///number of times for synthesis
  fscanf(input_file, "%g %g\n", &dum2, &dum3);///start and end times for synthesis
  fscanf(input_file, "%d\n", &dum1); ///model type
  fscanf(input_file, "%d\n", &dum1); ///compute the r-light curve?
  fscanf(input_file, "%d\n", &dum1); ///number of iterations
  fscanf(input_file, "%g\n", &dum2); ///change speed of light?
  fscanf(input_file, "%g\n", &dum2); ///use grey opacity for gammas?
  fscanf(input_file, "%g %g %g\n", &dum2, &dum3, &dum4); ///components of syn_dir
  fscanf(input_file, "%d\n", &dum1); ///opacity choice
  fscanf(input_file, "%g\n", &dum2); ///free parameter for calculation of rho_crit
  fscanf(input_file, "%d\n", &dum1); /// activate debug output for packet

  /// Flag continuation parameter as active
  fseek(input_file,0,SEEK_CUR);
  fprintf(input_file, "%d\n", 1); /// Force continuation
  //  fseek(input_file,0,SEEK_CUR);
  //fscanf(input_file, "%d", &dum1);  /// Do we start a new simulation or, continue another one?


  fclose(input_file);
}







///****************************************************************************
/// Subroutine to read in input parameters from vpkt.txt.

#ifdef VPKT_ON

  void read_parameterfile_vpkt()
  {
      FILE *input_file;
      int dum1,dum8;
      float dum2,dum3,dum4,dum5,dum6,dum7,dum9,dum10,dum11;
      int i;

      if ((input_file = fopen("vpkt.txt", "r")) == NULL)
      {
          printout("Cannot open vpkt.txt.\n");
          exit(0);
      }

      // Nobs
      fscanf(input_file, "%d", &dum1);
      Nobs = dum1 ;

      if ( Nobs > MOBS ) {

          printout("Too many observers! Nobs > MOBS \n");
          exit(0);
      }


      // nz_obs_vpkt. Cos(theta) to the observer. A list in the case of many observers
      for (i=0;i<Nobs;i++) {

          fscanf(input_file, "%g", &dum2);
          nz_obs_vpkt[i] = dum2;

          if ( fabs(nz_obs_vpkt[i])>1 ) {

              printout("Wrong observer direction \n");
              exit(0);
          }

          else if ( nz_obs_vpkt[i]==1 ) nz_obs_vpkt[i]=0.9999 ;
          else if ( nz_obs_vpkt[i]==-1 ) nz_obs_vpkt[i]=-0.9999 ;

      }


      // phi to the observer (degrees). A list in the case of many observers
      for (i=0;i<Nobs;i++) {

          fscanf(input_file, "%g \n", &dum3);
          phiobs[i] = dum3 * PI / 180 ;

      }

      // Nspectra opacity choices (i.e. Nspectra spectra for each observer)
      fscanf(input_file, "%g ", &dum4);

      if (dum4!=1) {

          Nspectra = 1;
          exclude[0] = 0;
      }

      else {

        fscanf(input_file, "%d ", &dum1);
        Nspectra = dum1;

        if ( Nspectra > MSPECTRA ) {

            printout("Too many spectra! Nspectra > MSPECTRA \n");
            exit(0);
        }

        
        for (i=0;i<Nspectra;i++) {
    
            fscanf(input_file, "%g ", &dum2);
            exclude[i] = dum2;

            // The first number should be equal to zero!
            if (exclude[0]!=0) {

              printout("The first spectrum should allow for all opacities (exclude[i]=0) and is not \n");
              exit(0);             
            }
    
        }
      }

      // time window. If dum4=1 it restrict vpkt to time windown (dum5,dum6)
      fscanf(input_file, "%g %g %g \n", &dum4, &dum5, &dum6);
    
      if (dum4==1) {
    
          tmin_vspec_input = dum5 * DAY;
          tmax_vspec_input = dum6 * DAY;
      }
    
      else {
    
          tmin_vspec_input = tmin_vspec;
          tmax_vspec_input = tmax_vspec;
      }
    
      if (tmin_vspec_input < tmin_vspec ) {
    
          printout("tmin_vspec_input is smaller than tmin_vspec \n");
          exit(0); 
      }
    
      if (tmax_vspec_input > tmax_vspec ) {
    
          printout("tmax_vspec_input is larger than tmax_vspec \n");
          exit(0) ;
      }


      // frequency window. dum4 restrict vpkt to a frequency range, dum5 indicates the number of ranges,
      // followed by a list of ranges (dum6,dum7)
      fscanf(input_file, "%g ", &dum4);
    
      if (dum4==1) {
    
          fscanf(input_file, "%g ", &dum5);
    
          Nrange = dum5;
          if ( Nrange > MRANGE ) {
    
              printout("Too many ranges! Nrange > MRANGE \n");
              exit(0);
          }
    
          for (i=0;i<Nrange;i++) {
    
              fscanf(input_file, "%g %g", &dum6, &dum7);
    
              lmin_vspec_input[i] = dum6;
              lmax_vspec_input[i] = dum7;
    
              numin_vspec_input[i] = CLIGHT / (lmax_vspec_input[i]*1e-8) ;
              numax_vspec_input[i] = CLIGHT / (lmin_vspec_input[i]*1e-8) ;
    
          }
    
      }
    
      else {
    
          Nrange = 1;
    
          numin_vspec_input[0] = numin_vspec ;
          numax_vspec_input[0] = numax_vspec ;
    
      }

      // if dum7=1, vpkt are not created when cell optical depth is larger than cell_is_optically_thick_vpkt
      fscanf(input_file, "%g %lg \n", &dum7, &cell_is_optically_thick_vpkt);

      if (dum7!=1) cell_is_optically_thick_vpkt = cell_is_optically_thick ;

      // Maximum optical depth. If a vpkt reaches dum7 is thrown away
      fscanf(input_file, "%g \n", &dum7);
      tau_max_vpkt = dum7 ;


      // Produce velocity grid map if dum8=1
      fscanf(input_file, "%d \n", &dum8);
      vgrid_flag = dum8 ;

      if (dum8==1) {
          
          // Specify time range for velocity grid map
          fscanf(input_file, "%g %g \n", &dum9,&dum10);
          tmin_grid = dum9 * DAY;
          tmax_grid = dum10 * DAY;

          // Specify wavelength range: number of intervals (dum9) and limits (dum10,dum11)
          fscanf(input_file, "%g ", &dum9);
          Nrange_grid = dum9 ;

          if ( Nrange_grid > MRANGE_GRID ) {

              printout("Too many ranges! Nrange_grid > MRANGE_GRID \n");
              exit(0);
          }

          for (i=0;i<Nrange_grid;i++) {

              fscanf(input_file, "%g %g", &dum10, &dum11);

              nu_grid_max[i] = CLIGHT / (dum10*1e-8) ;
              nu_grid_min[i] = CLIGHT / (dum11*1e-8) ;

          }

      }

      fclose(input_file);

  }

#endif



