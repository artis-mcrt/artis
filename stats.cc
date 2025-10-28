
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
  printlnlog("timestep {}: mean number of interactions per packet = {:g}", nts, meaninteractions);

  const double deltat = globals::timesteps[nts].width;
  double modelvolume = 0.;
  for (int mgi = 0; mgi < grid::get_npts_model(); mgi++) {
    modelvolume += grid::get_modelcell_assocvolume_tmin(mgi) * pow(globals::timesteps[nts].mid / globals::tmin, 3);
  }

  printlnlog("timestep {}: ma_stat_activation_collexc = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_COLLEXC));
  printlnlog("timestep {}: ma_stat_activation_collion = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_COLLION));
  printlnlog("timestep {}: ma_stat_activation_ntcollexc = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_NTCOLLEXC));
  printlnlog("timestep {}: ma_stat_activation_ntcollion = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_NTCOLLION));
  printlnlog("timestep {}: ma_stat_activation_bb = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_BB));
  printlnlog("timestep {}: ma_stat_activation_bf = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_BF));
  printlnlog("timestep {}: ma_stat_activation_fb = {}", nts, get_counter(COUNTER_MA_STAT_ACTIVATION_FB));
  printlnlog("timestep {}: ma_stat_deactivation_colldeexc = {}", nts,
             get_counter(COUNTER_MA_STAT_DEACTIVATION_COLLDEEXC));
  printlnlog("timestep {}: ma_stat_deactivation_collrecomb = {}", nts,
             get_counter(COUNTER_MA_STAT_DEACTIVATION_COLLRECOMB));
  printlnlog("timestep {}: ma_stat_deactivation_bb = {}", nts, get_counter(COUNTER_MA_STAT_DEACTIVATION_BB));
  printlnlog("timestep {}: ma_stat_deactivation_fb = {}", nts, get_counter(COUNTER_MA_STAT_DEACTIVATION_FB));
  printlnlog("timestep {}: ma_stat_internaluphigher = {}", nts, get_counter(COUNTER_MA_STAT_INTERNALUPHIGHER));
  printlnlog("timestep {}: ma_stat_internaluphighernt = {}", nts, get_counter(COUNTER_MA_STAT_INTERNALUPHIGHERNT));
  printlnlog("timestep {}: ma_stat_internaldownlower = {}", nts, get_counter(COUNTER_MA_STAT_INTERNALDOWNLOWER));

  printlnlog("timestep {}: k_stat_to_ma_collexc = {}", nts, get_counter(COUNTER_K_STAT_TO_MA_COLLEXC));
  printlnlog("timestep {}: k_stat_to_ma_collion = {}", nts, get_counter(COUNTER_K_STAT_TO_MA_COLLION));
  printlnlog("timestep {}: k_stat_to_r_ff = {}", nts, get_counter(COUNTER_K_STAT_TO_R_FF));
  printlnlog("timestep {}: k_stat_to_r_fb = {}", nts, get_counter(COUNTER_K_STAT_TO_R_FB));
  printlnlog("timestep {}: k_stat_to_r_bb = {}", nts, get_counter(COUNTER_K_STAT_TO_R_BB));
  printlnlog("timestep {}: k_stat_from_ff = {}", nts, get_counter(COUNTER_K_STAT_FROM_FF));
  printlnlog("timestep {}: k_stat_from_bf = {}", nts, get_counter(COUNTER_K_STAT_FROM_BF));
  printlnlog("timestep {}: k_stat_from_earlierdecay = {}", nts, get_counter(COUNTER_K_STAT_FROM_EARLIERDECAY));

  printlnlog("timestep {}: nt_stat_from_gamma = {}", nts, get_counter(COUNTER_NT_STAT_FROM_GAMMA));
  printlnlog("timestep {}: nt_stat_to_ionisation = {}", nts, get_counter(COUNTER_NT_STAT_TO_ionisation));
  printlnlog("timestep {}: nt_stat_to_excitation = {}", nts, get_counter(COUNTER_NT_STAT_TO_EXCITATION));
  printlnlog("timestep {}: nt_stat_to_kpkt = {}", nts, get_counter(COUNTER_NT_STAT_TO_KPKT));
  nonthermal::print_stats(modelvolume, deltat);

  printlnlog("timestep {}: escounter = {}", nts, get_counter(COUNTER_ESCOUNTER));
  printlnlog("timestep {}: cellcrossing  = {}", nts, get_counter(COUNTER_CELLCROSSINGS));
  printlnlog("timestep {}: updatecellcounter  = {}", nts, get_counter(COUNTER_UPDATECELL));
  printlnlog("timestep {}: resonancescatterings  = {}", nts, get_counter(COUNTER_RESONANCESCATTERINGS));

  printlnlog("timestep {}: upscatterings  = {}", nts, get_counter(COUNTER_UPSCATTER));
  printlnlog("timestep {}: downscatterings  = {}", nts, get_counter(COUNTER_DOWNSCATTER));
}

}  // namespace stats
