#include "stats.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>

#include "atomic.h"
#include "globals.h"
#include "grid.h"
#include "nonthermal.h"
#include "sn3d.h"

namespace stats {

namespace {

std::array<ptrdiff_t, COUNTER_COUNT> eventstats{};

}  // anonymous namespace

__host__ __device__ void increment(enum eventcounters i) {
  assert_testmodeonly(i >= 0);
  assert_testmodeonly(i < COUNTER_COUNT);
  atomicadd(eventstats[i], 1Z);
}

void pkt_action_counters_reset() {
  std::ranges::fill(eventstats, 0);

  nonthermal::reset_stats();
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
  printout("timestep %d: nt_stat_to_ionisation = %td\n", nts, get_counter(COUNTER_NT_STAT_TO_ionisation));
  printout("timestep %d: nt_stat_to_excitation = %td\n", nts, get_counter(COUNTER_NT_STAT_TO_EXCITATION));
  printout("timestep %d: nt_stat_to_kpkt = %td\n", nts, get_counter(COUNTER_NT_STAT_TO_KPKT));
  nonthermal::print_stats(modelvolume, deltat);

  printout("timestep %d: escounter = %td\n", nts, get_counter(COUNTER_ESCOUNTER));
  printout("timestep %d: cellcrossing  = %td\n", nts, get_counter(COUNTER_CELLCROSSINGS));
  printout("timestep %d: updatecellcounter  = %td\n", nts, get_counter(COUNTER_UPDATECELL));
  printout("timestep %d: resonancescatterings  = %td\n", nts, get_counter(COUNTER_RESONANCESCATTERINGS));

  printout("timestep %d: upscatterings  = %td\n", nts, get_counter(COUNTER_UPSCATTER));
  printout("timestep %d: downscatterings  = %td\n", nts, get_counter(COUNTER_DOWNSCATTER));
}

}  // namespace stats
