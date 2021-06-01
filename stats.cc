#include "atomic.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "nonthermal.h"
#include "stats.h"

namespace stats {

  static __managed__ double *ionstats = NULL;
  static __managed__ int *eventstats = NULL;


  void init(void)
  {
    if (TRACK_ION_STATS)
    {
      ionstats = (double *) malloc(grid::get_npts_model() * globals::includedions * ION_STAT_COUNT * sizeof(double));
    }
    eventstats = (int *) malloc(COUNTER_COUNT * sizeof(int));
  }


  void cleanup(void)
  {
    if (TRACK_ION_STATS)
    {
      free(ionstats);
    }
    free(eventstats);
  }


  __host__ __device__
  void increment_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstattypes ionstattype, const double increment)
  {
    if (!TRACK_ION_MASTATS && (ionstattype >= 18))
    {
      return;
    }

    assert_testmodeonly(ion < get_nions(element));
    assert_testmodeonly(ionstattype < ION_STAT_COUNT);

    const int uniqueionindex = get_uniqueionindex(element, ion);
    safeadd(ionstats[modelgridindex * globals::includedions * ION_STAT_COUNT + uniqueionindex * ION_STAT_COUNT + ionstattype], increment);
  }


  __host__ __device__
  void increment_ion_stats_contabsorption(const PKT *const pkt_ptr, const int modelgridindex, const int element, const int ion)
  {
    const double n_photons_absorbed = pkt_ptr->e_cmf / H / pkt_ptr->nu_cmf;

    stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION, n_photons_absorbed);

    const int et = pkt_ptr->emissiontype;
    if (et >= 0)  // r-packet is from bound-bound emission
    {
      stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUND, n_photons_absorbed);
      const int emissionelement = globals::linelist[et].elementindex;
      const int emissionion = globals::linelist[et].ionindex;

      stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_BOUNDBOUND_ABSORBED, n_photons_absorbed);

      if (emissionelement == element)
      {
        if (emissionion == ion + 1)
        {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSONE, n_photons_absorbed);
        }
        else if (emissionion == ion + 2)
        {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTWO, n_photons_absorbed);
        }
        else if (emissionion == ion + 3)
        {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTHREE, n_photons_absorbed);
        }
      }
    }
    else if (et != -9999999) // r-pkt is from bound-free emission (not free-free scattering)
    {
      stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBOUNDFREE, n_photons_absorbed);

      const int bfindex = -1*et - 1;
      assert_always(bfindex >= 0);
      assert_always(bfindex <= globals::nbfcontinua);
      const int emissionelement = globals::bflist[bfindex].elementindex;
      const int emissionlowerion = globals::bflist[bfindex].ionindex;
      const int emissionupperion = emissionlowerion + 1;
      const int emissionlowerlevel = globals::bflist[bfindex].levelindex;

      stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_RADRECOMB_ABSORBED, n_photons_absorbed);

      if (emissionelement == element)
      {
        stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBFSAMEELEMENT, n_photons_absorbed);
        if (emissionupperion == ion + 1)
        {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSONE, n_photons_absorbed);
        }
        else if (emissionupperion == ion + 2)
        {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSTWO, n_photons_absorbed);
        }
        else if (emissionupperion == ion + 3)
        {
          stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSTHREE, n_photons_absorbed);
        }
      }
      if (level_isinsuperlevel(emissionelement, emissionlowerion, emissionlowerlevel))
      {
        stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBFLOWERSUPERLEVEL, n_photons_absorbed);
      }
    }
  }


  double get_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstattypes ionstattype)
  {
    assert_always(ion < get_nions(element));
    assert_always(ionstattype < ION_STAT_COUNT);
    const int uniqueionindex = get_uniqueionindex(element, ion);
    return ionstats[modelgridindex * globals::includedions * ION_STAT_COUNT + uniqueionindex * ION_STAT_COUNT + ionstattype];
  }


  void set_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstattypes ionstattype, const double newvalue)
  {
    assert_always(ion < get_nions(element));
    assert_always(ionstattype < ION_STAT_COUNT);
    const int uniqueionindex = get_uniqueionindex(element, ion);
    ionstats[modelgridindex * globals::includedions * ION_STAT_COUNT + uniqueionindex * ION_STAT_COUNT + ionstattype] = newvalue;
  }


  void reset_ion_stats(int modelgridindex)
  {
    for (int element = 0; element < get_nelements(); element++)
    {
      for (int ion = 0; ion < get_nions(element); ion++)
      {
        for (int i = 0; i < ION_STAT_COUNT; i++)
        {
          set_ion_stats(modelgridindex, element, ion, (enum stats::ionstattypes)i, 0.);
        }
      }
    }
  }


  void normalise_ion_estimators(const int mgi, const double deltat, const double deltaV)
  {
    for (int element = 0; element < get_nelements(); element++)
    {
      for (int ion = 0; ion < get_nions(element); ion++)
      {
        for (int i = 0; i < ION_STAT_COUNT; i++)
        {
          // energy or event count per volume per second
          const double ratedensity = get_ion_stats(mgi, element, ion, (enum stats::ionstattypes)i) / deltaV / deltat / globals::nprocs;

          if (i < nstatcounters_ratecoeff)
          {
            // convert photon event counters into rate coefficients
            set_ion_stats(mgi, element, ion, (enum stats::ionstattypes)i, ratedensity / ionstagepop(mgi, element, ion));
          }
          else
          {
            set_ion_stats(mgi, element, ion, (enum stats::ionstattypes)i, ratedensity);
          }
        }
      }
    }
  }

  __host__ __device__
  void increment(enum eventcounters i)
  {
    assert_testmodeonly(i >= 0);
    assert_testmodeonly(i < COUNTER_COUNT);
    safeincrement(eventstats[i]);
  }

  void pkt_action_counters_reset(void)
  {
    for (int i = 0; i < COUNTER_COUNT; i++)
    {
      eventstats[i] = 0;
    }

    nonthermal::nt_reset_stats();
  }


  int get_counter(enum eventcounters i)
  {
    assert_always(i < COUNTER_COUNT);
    return eventstats[i];
  }


  void pkt_action_counters_printout(const PKT *const pkt, const int nts)
  {
    int allpktinteractions = 0;
    for (int i = 0; i < globals::npkts; i++)
    {
      allpktinteractions += pkt[i].interactions;
    }
    const double meaninteractions = (double) allpktinteractions / globals::npkts;
    printout("mean number of interactions per packet = %g\n", meaninteractions);

    const double deltat = globals::time_step[nts].width;
    double modelvolume = 0.;
    for (int mgi = 0; mgi < grid::get_npts_model(); mgi++)
    {
      modelvolume += grid::vol_init_modelcell(mgi) * pow(globals::time_step[nts].mid / globals::tmin, 3);
    }

    /// Printout packet statistics
    printout("ma_stat_activation_collexc = %d\n", get_counter(COUNTER_MA_STAT_ACTIVATION_COLLEXC));
    printout("ma_stat_activation_collion = %d\n", get_counter(COUNTER_MA_STAT_ACTIVATION_COLLION));
    printout("ma_stat_activation_ntcollexc = %d\n", get_counter(COUNTER_MA_STAT_ACTIVATION_NTCOLLEXC));
    printout("ma_stat_activation_ntcollion = %d\n", get_counter(COUNTER_MA_STAT_ACTIVATION_NTCOLLION));
    printout("ma_stat_activation_bb = %d\n", get_counter(COUNTER_MA_STAT_ACTIVATION_BB));
    printout("ma_stat_activation_bf = %d\n", get_counter(COUNTER_MA_STAT_ACTIVATION_BF));
    printout("ma_stat_activation_fb = %d\n", get_counter(COUNTER_MA_STAT_ACTIVATION_FB));
    printout("ma_stat_deactivation_colldeexc = %d\n", get_counter(COUNTER_MA_STAT_DEACTIVATION_COLLDEEXC));
    printout("ma_stat_deactivation_collrecomb = %d\n", get_counter(COUNTER_MA_STAT_DEACTIVATION_COLLRECOMB));
    printout("ma_stat_deactivation_bb = %d\n", get_counter(COUNTER_MA_STAT_DEACTIVATION_BB));
    printout("ma_stat_deactivation_fb = %d\n", get_counter(COUNTER_MA_STAT_DEACTIVATION_FB));
    printout("ma_stat_internaluphigher = %d\n", get_counter(COUNTER_MA_STAT_INTERNALUPHIGHER));
    printout("ma_stat_internaluphighernt = %d\n", get_counter(COUNTER_MA_STAT_INTERNALUPHIGHERNT));
    printout("ma_stat_internaldownlower = %d\n", get_counter(COUNTER_MA_STAT_INTERNALDOWNLOWER));

    printout("k_stat_to_ma_collexc = %d\n", get_counter(COUNTER_K_STAT_TO_MA_COLLEXC));
    printout("k_stat_to_ma_collion = %d\n", get_counter(COUNTER_K_STAT_TO_MA_COLLION));
    printout("k_stat_to_r_ff = %d\n", get_counter(COUNTER_K_STAT_TO_R_FF));
    printout("k_stat_to_r_fb = %d\n", get_counter(COUNTER_K_STAT_TO_R_FB));
    printout("k_stat_to_r_bb = %d\n", get_counter(COUNTER_K_STAT_TO_R_BB));
    printout("k_stat_from_ff = %d\n", get_counter(COUNTER_K_STAT_FROM_FF));
    printout("k_stat_from_bf = %d\n", get_counter(COUNTER_K_STAT_FROM_BF));
    printout("k_stat_from_earlierdecay = %d\n", get_counter(COUNTER_K_STAT_FROM_EARLIERDECAY));

    printout("nt_stat_from_gamma = %d\n", get_counter(COUNTER_NT_STAT_FROM_GAMMA));
    printout("nt_stat_to_ionization = %d\n", get_counter(COUNTER_NT_STAT_TO_IONIZATION));
    printout("nt_stat_to_excitation = %d\n", get_counter(COUNTER_NT_STAT_TO_EXCITATION));
    printout("nt_stat_to_kpkt = %d\n", get_counter(COUNTER_NT_STAT_TO_KPKT));
    nonthermal::nt_print_stats(nts, modelvolume, deltat);

    printout("escounter = %d\n", get_counter(COUNTER_ESCOUNTER));
    printout("cellcrossing  = %d\n", get_counter(COUNTER_CELLCROSSINGS));
    printout("updatecellcounter  = %d\n", get_counter(COUNTER_UPDATECELL));
    printout("coolingratecalccounter = %d\n", get_counter(COUNTER_COOLINGRATECALCCOUNTER));
    printout("resonancescatterings  = %d\n", get_counter(COUNTER_RESONANCESCATTERINGS));

    printout("upscatterings  = %d\n", get_counter(COUNTER_UPSCATTER));
    printout("downscatterings  = %d\n", get_counter(COUNTER_DOWNSCATTER));
  }


  void reduce_estimators(void)
  {
    #ifdef MPI_ON
    MPI_Allreduce(MPI_IN_PLACE, &stats::ionstats, grid::get_npts_model() * globals::includedions * stats::ION_STAT_COUNT, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    #endif
  }
}
