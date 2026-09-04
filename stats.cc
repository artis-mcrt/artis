// Counters of Monte Carlo packet events (interactions, macro-atom transitions, estimator
// contributions), collected per timestep and printed to the log.

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

// Atomically add one to the given global event counter (safe to call from all threads).
DEVICE_FUNC void increment(const Counter i) { atomicadd(eventstats[std::to_underlying(i)].count, 1Z); }

// Reset all global event counters (and the non-thermal solver's stats) to zero, called at the start of each timestep.
void pkt_action_counters_reset() {
  std::ranges::fill(eventstats, EventCounter{0});

  nonthermal::reset_stats();
}

// Return the current value of a global event counter.
auto get_counter(const Counter i) -> ptrdiff_t { return eventstats[std::to_underlying(i)].count; }

// Log all packet-interaction statistics (macroatom/kpkt/non-thermal transitions, scatterings, cell crossings, ...)
// accumulated over timestep nts.
void pkt_action_counters_printout(const int nts) {
  const double meaninteractions = static_cast<double>(get_counter(Counter::INTERACTIONS)) / MPKTS;
  printlnlog("timestep {}: mean number of interactions per packet = {:g}", nts, meaninteractions);

  const double deltat = globals::timesteps[nts].width;
  double modelvolume = 0.;
  for (int mgi = 0; mgi < grid::get_npts_model(); mgi++) {
    modelvolume += grid::get_modelcell_assocvolume_tmin(mgi) * pow3(globals::timesteps[nts].mid / globals::tmin);
  }

  printlnlog("timestep {}: ma_stat_activation: collexc {} collion {} ntcollexc {} ntcollion {} bb {} bf {}", nts,
             get_counter(Counter::MA_STAT_ACTIVATION_COLLEXC), get_counter(Counter::MA_STAT_ACTIVATION_COLLION),
             get_counter(Counter::MA_STAT_ACTIVATION_NTCOLLEXC), get_counter(Counter::MA_STAT_ACTIVATION_NTCOLLION),
             get_counter(Counter::MA_STAT_ACTIVATION_BB), get_counter(Counter::MA_STAT_ACTIVATION_BF));
  printlnlog("timestep {}: ma_stat_deactivation: colldeexc {} collrecomb {} bb {} fb {}", nts,
             get_counter(Counter::MA_STAT_DEACTIVATION_COLLDEEXC),
             get_counter(Counter::MA_STAT_DEACTIVATION_COLLRECOMB), get_counter(Counter::MA_STAT_DEACTIVATION_BB),
             get_counter(Counter::MA_STAT_DEACTIVATION_FB));
  printlnlog("timestep {}: ma_stat_internal: uphigher {} uphighernt {} downlower {}", nts,
             get_counter(Counter::MA_STAT_INTERNALUPHIGHER), get_counter(Counter::MA_STAT_INTERNALUPHIGHERNT),
             get_counter(Counter::MA_STAT_INTERNALDOWNLOWER));

  printlnlog("timestep {}: k_stat_to: ma_collexc {} ma_collion {} r_ff {} r_fb {} r_bb {}", nts,
             get_counter(Counter::K_STAT_TO_MA_COLLEXC), get_counter(Counter::K_STAT_TO_MA_COLLION),
             get_counter(Counter::K_STAT_TO_R_FF), get_counter(Counter::K_STAT_TO_R_FB),
             get_counter(Counter::K_STAT_TO_R_BB));
  printlnlog("timestep {}: k_stat_from: ff {} bf {} earlierdecay {}", nts, get_counter(Counter::K_STAT_FROM_FF),
             get_counter(Counter::K_STAT_FROM_BF), get_counter(Counter::K_STAT_FROM_EARLIERDECAY));

  printlnlog("timestep {}: nt_stat: from_gamma {} to_ionisation {} to_excitation {} to_kpkt {}", nts,
             get_counter(Counter::NT_STAT_FROM_GAMMA), get_counter(Counter::NT_STAT_TO_IONISATION),
             get_counter(Counter::NT_STAT_TO_EXCITATION), get_counter(Counter::NT_STAT_TO_KPKT));
  nonthermal::print_stats(modelvolume, deltat);

  printlnlog("timestep {}: electron_scatterings {} resonancescatterings {} upscatterings {} downscatterings {}", nts,
             get_counter(Counter::ELECTRON_SCATTERINGS), get_counter(Counter::RESONANCESCATTERINGS),
             get_counter(Counter::UPSCATTER), get_counter(Counter::DOWNSCATTER));
  printlnlog("timestep {}: cellcrossings {} updatecellcounter {}", nts, get_counter(Counter::CELLCROSSINGS),
             get_counter(Counter::UPDATECELL));
}

}  // namespace stats
