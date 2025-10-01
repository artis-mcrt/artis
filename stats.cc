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
  logprintlnfmt("timestep {}: mean number of interactions per packet = {:g}", nts, meaninteractions);

  const double deltat = globals::timesteps[nts].width;
  double modelvolume = 0.;
  for (int mgi = 0; mgi < grid::get_npts_model(); mgi++) {
    modelvolume += grid::get_modelcell_assocvolume_tmin(mgi) * pow(globals::timesteps[nts].mid / globals::tmin, 3);
  }

  logprintlnfmt("timestep {}: ma_stat_activation_collexc = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_COLLEXC));
  logprintlnfmt("timestep {}: ma_stat_activation_collion = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_COLLION));
  logprintlnfmt("timestep {}: ma_stat_activation_ntcollexc = {}", nts,
                get_counter(COUNTER_MA_STAT_ACTIVATION_NTCOLLEXC));
  logprintlnfmt("timestep {}: ma_stat_activation_ntcollion = {}", nts,
                get_counter(COUNTER_MA_STAT_ACTIVATION_NTCOLLION));
  logprintlnfmt("timestep {}: ma_stat_activation_bb = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_BB));
  logprintlnfmt("timestep {}: ma_stat_activation_bf = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_BF));
  logprintlnfmt("timestep {}: ma_stat_activation_fb = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_FB));
  logprintlnfmt("timestep {}: ma_stat_deactivation_colldeexc = {}", nts,
                get_counter(COUNTER_MA_STAT_DEACTIVATION_COLLDEEXC));
  logprintlnfmt("timestep {}: ma_stat_deactivation_collrecomb = {}", nts,
                get_counter(COUNTER_MA_STAT_DEACTIVATION_COLLRECOMB));
  logprintlnfmt("timestep {}: ma_stat_deactivation_bb = {}", nts, get_counter(COUNTER_MA_STAT_DEACTIVATION_BB));
  logprintlnfmt("timestep {}: ma_stat_deactivation_fb = {}", nts, get_counter(COUNTER_MA_STAT_DEACTIVATION_FB));
  logprintlnfmt("timestep {}: ma_stat_internaluphigher = {}", nts, get_counter(COUNTER_MA_STAT_INTERNALUPHIGHER));
  logprintlnfmt("timestep {}: ma_stat_internaluphighernt = {}", nts, get_counter(COUNTER_MA_STAT_INTERNALUPHIGHERNT));
  logprintlnfmt("timestep {}: ma_stat_internaldownlower = {}", nts, get_counter(COUNTER_MA_STAT_INTERNALDOWNLOWER));

  logprintlnfmt("timestep {}: k_stat_to_ma_collexc = {}", nts, get_counter(COUNTER_K_STAT_TO_MA_COLLEXC));
  logprintlnfmt("timestep {}: k_stat_to_ma_collion = {}", nts, get_counter(COUNTER_K_STAT_TO_MA_COLLION));
  logprintlnfmt("timestep {}: k_stat_to_r_ff = {}", nts, get_counter(COUNTER_K_STAT_TO_R_FF));
  logprintlnfmt("timestep {}: k_stat_to_r_fb = {}", nts, get_counter(COUNTER_K_STAT_TO_R_FB));
  logprintlnfmt("timestep {}: k_stat_to_r_bb = {}", nts, get_counter(COUNTER_K_STAT_TO_R_BB));
  logprintlnfmt("timestep {}: k_stat_from_ff = {}", nts, get_counter(COUNTER_K_STAT_FROM_FF));
  logprintlnfmt("timestep {}: k_stat_from_bf = {}", nts, get_counter(COUNTER_K_STAT_FROM_BF));
  logprintlnfmt("timestep {}: k_stat_from_earlierdecay = {}", nts, get_counter(COUNTER_K_STAT_FROM_EARLIERDECAY));

  logprintlnfmt("timestep {}: nt_stat_from_gamma = {}", nts, get_counter(COUNTER_NT_STAT_FROM_GAMMA));
  logprintlnfmt("timestep {}: nt_stat_to_ionisation = {}", nts, get_counter(COUNTER_NT_STAT_TO_ionisation));
  logprintlnfmt("timestep {}: nt_stat_to_excitation = {}", nts, get_counter(COUNTER_NT_STAT_TO_EXCITATION));
  logprintlnfmt("timestep {}: nt_stat_to_kpkt = {}", nts, get_counter(COUNTER_NT_STAT_TO_KPKT));
  nonthermal::print_stats(modelvolume, deltat);

  logprintlnfmt("timestep {}: escounter = {}", nts, get_counter(COUNTER_ESCOUNTER));
  logprintlnfmt("timestep {}: cellcrossing  = {}", nts, get_counter(COUNTER_CELLCROSSINGS));
  logprintlnfmt("timestep {}: updatecellcounter  = {}", nts, get_counter(COUNTER_UPDATECELL));
  logprintlnfmt("timestep {}: resonancescatterings  = {}", nts, get_counter(COUNTER_RESONANCESCATTERINGS));

  logprintlnfmt("timestep {}: upscatterings  = {}", nts, get_counter(COUNTER_UPSCATTER));
  logprintlnfmt("timestep {}: downscatterings  = {}", nts, get_counter(COUNTER_DOWNSCATTER));
}

}  // namespace stats
