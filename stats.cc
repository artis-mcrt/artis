
#include "stats.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <utility>

#include "artisoptions.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "mpi_logging.h"
#include "nonthermal.h"

namespace stats {

namespace {

// counters are frequently incremented by all threads, so give each one its own cache line to prevent false sharing
struct EventCounter {
  ALIGNAS_AVOID_FALSE_SHARING ptrdiff_t count{0};
};

std::array<EventCounter, static_cast<size_t>(Counter::COUNT)> eventstats{};

}  // anonymous namespace

DEVICE_FUNC void increment(const Counter i) { atomicadd(eventstats[std::to_underlying(i)].count, 1Z); }

void pkt_action_counters_reset() {
  std::ranges::fill(eventstats, EventCounter{0});

  nonthermal::reset_stats();
}

auto get_counter(const Counter i) -> ptrdiff_t { return eventstats[std::to_underlying(i)].count; }

void pkt_action_counters_printout(const int nts) {
  const double meaninteractions = static_cast<double>(get_counter(Counter::INTERACTIONS)) / MPKTS;
  printlnlog("timestep {}: mean number of interactions per packet = {:g}", nts, meaninteractions);

  const double deltat = globals::timesteps[nts].width;
  double modelvolume = 0.;
  for (int mgi = 0; mgi < grid::get_npts_model(); mgi++) {
    modelvolume += grid::get_modelcell_assocvolume_tmin(mgi) * pow3(globals::timesteps[nts].mid / globals::tmin);
  }

  printlnlog("timestep {}: ma_stat_activation_collexc = {}", nts, get_counter(Counter::MA_STAT_ACTIVATION_COLLEXC));
  printlnlog("timestep {}: ma_stat_activation_collion = {}", nts, get_counter(Counter::MA_STAT_ACTIVATION_COLLION));
  printlnlog("timestep {}: ma_stat_activation_ntcollexc = {}", nts, get_counter(Counter::MA_STAT_ACTIVATION_NTCOLLEXC));
  printlnlog("timestep {}: ma_stat_activation_ntcollion = {}", nts, get_counter(Counter::MA_STAT_ACTIVATION_NTCOLLION));
  printlnlog("timestep {}: ma_stat_activation_bb = {}", nts, get_counter(Counter::MA_STAT_ACTIVATION_BB));
  printlnlog("timestep {}: ma_stat_activation_bf = {}", nts, get_counter(Counter::MA_STAT_ACTIVATION_BF));
  printlnlog("timestep {}: ma_stat_activation_fb = {}", nts, get_counter(Counter::MA_STAT_ACTIVATION_FB));
  printlnlog("timestep {}: ma_stat_deactivation_colldeexc = {}", nts,
             get_counter(Counter::MA_STAT_DEACTIVATION_COLLDEEXC));
  printlnlog("timestep {}: ma_stat_deactivation_collrecomb = {}", nts,
             get_counter(Counter::MA_STAT_DEACTIVATION_COLLRECOMB));
  printlnlog("timestep {}: ma_stat_deactivation_bb = {}", nts, get_counter(Counter::MA_STAT_DEACTIVATION_BB));
  printlnlog("timestep {}: ma_stat_deactivation_fb = {}", nts, get_counter(Counter::MA_STAT_DEACTIVATION_FB));
  printlnlog("timestep {}: ma_stat_internaluphigher = {}", nts, get_counter(Counter::MA_STAT_INTERNALUPHIGHER));
  printlnlog("timestep {}: ma_stat_internaluphighernt = {}", nts, get_counter(Counter::MA_STAT_INTERNALUPHIGHERNT));
  printlnlog("timestep {}: ma_stat_internaldownlower = {}", nts, get_counter(Counter::MA_STAT_INTERNALDOWNLOWER));

  printlnlog("timestep {}: k_stat_to_ma_collexc = {}", nts, get_counter(Counter::K_STAT_TO_MA_COLLEXC));
  printlnlog("timestep {}: k_stat_to_ma_collion = {}", nts, get_counter(Counter::K_STAT_TO_MA_COLLION));
  printlnlog("timestep {}: k_stat_to_r_ff = {}", nts, get_counter(Counter::K_STAT_TO_R_FF));
  printlnlog("timestep {}: k_stat_to_r_fb = {}", nts, get_counter(Counter::K_STAT_TO_R_FB));
  printlnlog("timestep {}: k_stat_to_r_bb = {}", nts, get_counter(Counter::K_STAT_TO_R_BB));
  printlnlog("timestep {}: k_stat_from_ff = {}", nts, get_counter(Counter::K_STAT_FROM_FF));
  printlnlog("timestep {}: k_stat_from_bf = {}", nts, get_counter(Counter::K_STAT_FROM_BF));
  printlnlog("timestep {}: k_stat_from_earlierdecay = {}", nts, get_counter(Counter::K_STAT_FROM_EARLIERDECAY));

  printlnlog("timestep {}: nt_stat_from_gamma = {}", nts, get_counter(Counter::NT_STAT_FROM_GAMMA));
  printlnlog("timestep {}: nt_stat_to_ionisation = {}", nts, get_counter(Counter::NT_STAT_TO_IONISATION));
  printlnlog("timestep {}: nt_stat_to_excitation = {}", nts, get_counter(Counter::NT_STAT_TO_EXCITATION));
  printlnlog("timestep {}: nt_stat_to_kpkt = {}", nts, get_counter(Counter::NT_STAT_TO_KPKT));
  nonthermal::print_stats(modelvolume, deltat);

  printlnlog("timestep {}: electron_scatterings = {}", nts, get_counter(Counter::ELECTRON_SCATTERINGS));
  printlnlog("timestep {}: cellcrossings = {}", nts, get_counter(Counter::CELLCROSSINGS));
  printlnlog("timestep {}: updatecellcounter = {}", nts, get_counter(Counter::UPDATECELL));
  printlnlog("timestep {}: resonancescatterings = {}", nts, get_counter(Counter::RESONANCESCATTERINGS));

  printlnlog("timestep {}: upscatterings = {}", nts, get_counter(Counter::UPSCATTER));
  printlnlog("timestep {}: downscatterings = {}", nts, get_counter(Counter::DOWNSCATTER));
}

}  // namespace stats
