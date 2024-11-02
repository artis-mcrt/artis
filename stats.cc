#include "stats.h"

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "ltepop.h"
#include "nonthermal.h"
#include "packet.h"
#include "sn3d.h"

namespace stats {

namespace {

std::vector<double> ionstats;
std::array<ptrdiff_t, COUNTER_COUNT> eventstats{};

}  // anonymous namespace

void init() {
  if constexpr (TRACK_ION_STATS) {
    ionstats.resize(grid::get_npts_model() * get_includedions() * ION_STAT_COUNT, 0.);
  }
}

void increment_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstattypes ionstattype,
                         const double increment) {
  if (ionstattype >= 18) {
    return;
  }

  assert_testmodeonly(ion < get_nions(element));
  assert_testmodeonly(ionstattype < ION_STAT_COUNT);

  const int uniqueionindex = get_uniqueionindex(element, ion);
  atomicadd(ionstats[(modelgridindex * get_includedions() * ION_STAT_COUNT) + (uniqueionindex * ION_STAT_COUNT) +
                     ionstattype],
            increment);
}

void increment_ion_stats_contabsorption(const Packet &pkt, const int modelgridindex, const int element, const int ion) {
  const double n_photons_absorbed = pkt.e_cmf / H / pkt.nu_cmf;

  stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION, n_photons_absorbed);

  const int et = pkt.emissiontype;
  if (et >= 0) {
    // r-packet is from bound-bound emission
    stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUND, n_photons_absorbed);
    const int emissionelement = globals::linelist[et].elementindex;
    const int emissionion = globals::linelist[et].ionindex;

    stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_BOUNDBOUND_ABSORBED, n_photons_absorbed);

    if (emissionelement == element) {
      if (emissionion == ion + 1) {
        stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSONE,
                                   n_photons_absorbed);
      } else if (emissionion == ion + 2) {
        stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTWO,
                                   n_photons_absorbed);
      } else if (emissionion == ion + 3) {
        stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBOUNDBOUNDIONPLUSTHREE,
                                   n_photons_absorbed);
      }
    }
  } else if (et != EMTYPE_FREEFREE && et != EMTYPE_NOTSET) {
    // r-pkt is from bound-free emission (not free-free scattering)
    stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBOUNDFREE, n_photons_absorbed);

    const int bfindex = (-1 * et) - 1;
    assert_always(bfindex >= 0);
    assert_always(bfindex <= globals::nbfcontinua);
    const int emissionelement = globals::bflist[bfindex].elementindex;
    const int emissionlowerion = globals::bflist[bfindex].ionindex;
    const int emissionupperion = emissionlowerion + 1;
    const int emissionlowerlevel = globals::bflist[bfindex].levelindex;

    stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_RADRECOMB_ABSORBED, n_photons_absorbed);

    if (emissionelement == element) {
      stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBFSAMEELEMENT,
                                 n_photons_absorbed);
      if (emissionupperion == ion + 1) {
        stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSONE,
                                   n_photons_absorbed);
      } else if (emissionupperion == ion + 2) {
        stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSTWO,
                                   n_photons_absorbed);
      } else if (emissionupperion == ion + 3) {
        stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBFIONPLUSTHREE,
                                   n_photons_absorbed);
      }
    }
    if (level_isinsuperlevel(emissionelement, emissionlowerion, emissionlowerlevel)) {
      stats::increment_ion_stats(modelgridindex, element, ion, stats::ION_PHOTOION_FROMBFLOWERSUPERLEVEL,
                                 n_photons_absorbed);
    }
  }
}

auto get_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstattypes ionstattype)
    -> double {
  assert_always(ion < get_nions(element));
  assert_always(ionstattype < ION_STAT_COUNT);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  return ionstats[(modelgridindex * get_includedions() * ION_STAT_COUNT) + (uniqueionindex * ION_STAT_COUNT) +
                  ionstattype];
}

void set_ion_stats(const int modelgridindex, const int element, const int ion, enum ionstattypes ionstattype,
                   const double newvalue) {
  assert_always(ion < get_nions(element));
  assert_always(ionstattype < ION_STAT_COUNT);
  const int uniqueionindex = get_uniqueionindex(element, ion);
  ionstats[(modelgridindex * get_includedions() * ION_STAT_COUNT) + (uniqueionindex * ION_STAT_COUNT) + ionstattype] =
      newvalue;
}

void reset_ion_stats(int modelgridindex) {
  for (int element = 0; element < get_nelements(); element++) {
    for (int ion = 0; ion < get_nions(element); ion++) {
      for (int i = 0; i < ION_STAT_COUNT; i++) {
        set_ion_stats(modelgridindex, element, ion, static_cast<enum stats::ionstattypes>(i), 0.);
      }
    }
  }
}

void normalise_ion_estimators(const int mgi, const double deltat, const double deltaV) {
  const auto nonemptymgi = grid::get_nonemptymgi_of_mgi(mgi);
  for (int element = 0; element < get_nelements(); element++) {
    for (int ion = 0; ion < get_nions(element); ion++) {
      for (int i = 0; i < ION_STAT_COUNT; i++) {
        // energy or event count per volume per second
        const double ratedensity = get_ion_stats(mgi, element, ion, static_cast<enum stats::ionstattypes>(i)) / deltaV /
                                   deltat / globals::nprocs;

        if (i < nstatcounters_ratecoeff) {
          // convert photon event counters into rate coefficients
          set_ion_stats(mgi, element, ion, static_cast<enum stats::ionstattypes>(i),
                        ratedensity / get_nnion(nonemptymgi, element, ion));
        } else {
          set_ion_stats(mgi, element, ion, static_cast<enum stats::ionstattypes>(i), ratedensity);
        }
      }
    }
  }
}

__host__ __device__ void increment(enum eventcounters i) {
  assert_testmodeonly(i >= 0);
  assert_testmodeonly(i < COUNTER_COUNT);
  atomicadd(eventstats[i], static_cast<ptrdiff_t>(1));
}

void pkt_action_counters_reset() {
  for (int i = 0; i < COUNTER_COUNT; i++) {
    eventstats[i] = 0;
  }

  nonthermal::nt_reset_stats();
  globals::nesc = 0;
}

auto get_counter(enum eventcounters i) -> ptrdiff_t {
  assert_testmodeonly(i >= 0);
  assert_testmodeonly(i < COUNTER_COUNT);
  return eventstats[i];
}

void pkt_action_counters_printout(const int nts) {
  const double meaninteractions = static_cast<double>(get_counter(COUNTER_INTERACTIONS)) / globals::npkts;
  printout("timestep %d: mean number of interactions per packet = %g\n", nts, meaninteractions);

  const double deltat = globals::timesteps[nts].width;
  double modelvolume = 0.;
  for (int mgi = 0; mgi < grid::get_npts_model(); mgi++) {
    modelvolume += grid::get_modelcell_assocvolume_tmin(mgi) * pow(globals::timesteps[nts].mid / globals::tmin, 3);
  }

  // Printout packet statistics
  printout("timestep %d: ma_stat_activation_collexc = %td\n", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_COLLEXC));
  printout("timestep %d: ma_stat_activation_collion = %td\n", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_COLLION));
  printout("timestep %d: ma_stat_activation_ntcollexc = %td\n", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_NTCOLLEXC));
  printout("timestep %d: ma_stat_activation_ntcollion = %td\n", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_NTCOLLION));
  printout("timestep %d: ma_stat_activation_bb = %td\n", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_BB));
  printout("timestep %d: ma_stat_activation_bf = %td\n", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_BF));
  printout("timestep %d: ma_stat_activation_fb = %td\n", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_FB));
  printout("timestep %d: ma_stat_deactivation_colldeexc = %td\n", nts,
           get_counter(COUNTER_MA_STAT_DEACTIVATION_COLLDEEXC));
  printout("timestep %d: ma_stat_deactivation_collrecomb = %td\n", nts,
           get_counter(COUNTER_MA_STAT_DEACTIVATION_COLLRECOMB));
  printout("timestep %d: ma_stat_deactivation_bb = %td\n", nts, get_counter(COUNTER_MA_STAT_DEACTIVATION_BB));
  printout("timestep %d: ma_stat_deactivation_fb = %td\n", nts, get_counter(COUNTER_MA_STAT_DEACTIVATION_FB));
  printout("timestep %d: ma_stat_internaluphigher = %td\n", nts, get_counter(COUNTER_MA_STAT_INTERNALUPHIGHER));
  printout("timestep %d: ma_stat_internaluphighernt = %td\n", nts, get_counter(COUNTER_MA_STAT_INTERNALUPHIGHERNT));
  printout("timestep %d: ma_stat_internaldownlower = %td\n", nts, get_counter(COUNTER_MA_STAT_INTERNALDOWNLOWER));

  printout("timestep %d: k_stat_to_ma_collexc = %td\n", nts, get_counter(COUNTER_K_STAT_TO_MA_COLLEXC));
  printout("timestep %d: k_stat_to_ma_collion = %td\n", nts, get_counter(COUNTER_K_STAT_TO_MA_COLLION));
  printout("timestep %d: k_stat_to_r_ff = %td\n", nts, get_counter(COUNTER_K_STAT_TO_R_FF));
  printout("timestep %d: k_stat_to_r_fb = %td\n", nts, get_counter(COUNTER_K_STAT_TO_R_FB));
  printout("timestep %d: k_stat_to_r_bb = %td\n", nts, get_counter(COUNTER_K_STAT_TO_R_BB));
  printout("timestep %d: k_stat_from_ff = %td\n", nts, get_counter(COUNTER_K_STAT_FROM_FF));
  printout("timestep %d: k_stat_from_bf = %td\n", nts, get_counter(COUNTER_K_STAT_FROM_BF));
  printout("timestep %d: k_stat_from_earlierdecay = %td\n", nts, get_counter(COUNTER_K_STAT_FROM_EARLIERDECAY));

  printout("timestep %d: nt_stat_from_gamma = %td\n", nts, get_counter(COUNTER_NT_STAT_FROM_GAMMA));
  printout("timestep %d: nt_stat_to_ionization = %td\n", nts, get_counter(COUNTER_NT_STAT_TO_IONIZATION));
  printout("timestep %d: nt_stat_to_excitation = %td\n", nts, get_counter(COUNTER_NT_STAT_TO_EXCITATION));
  printout("timestep %d: nt_stat_to_kpkt = %td\n", nts, get_counter(COUNTER_NT_STAT_TO_KPKT));
  nonthermal::nt_print_stats(modelvolume, deltat);

  printout("timestep %d: escounter = %td\n", nts, get_counter(COUNTER_ESCOUNTER));
  printout("timestep %d: cellcrossing  = %td\n", nts, get_counter(COUNTER_CELLCROSSINGS));
  printout("timestep %d: updatecellcounter  = %td\n", nts, get_counter(COUNTER_UPDATECELL));
  printout("timestep %d: resonancescatterings  = %td\n", nts, get_counter(COUNTER_RESONANCESCATTERINGS));

  printout("timestep %d: upscatterings  = %td\n", nts, get_counter(COUNTER_UPSCATTER));
  printout("timestep %d: downscatterings  = %td\n", nts, get_counter(COUNTER_DOWNSCATTER));
}

void reduce_estimators() {
#ifdef MPI_ON
  MPI_Allreduce(MPI_IN_PLACE, stats::ionstats.data(),
                grid::get_npts_model() * get_includedions() * stats::ION_STAT_COUNT, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#endif
}
}  // namespace stats
