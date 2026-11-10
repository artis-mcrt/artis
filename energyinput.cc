// The energy source of THERMAL_ENERGY_INPUT_FROM_FILES. grid.cc reads the energy of each model cell from
// energydistribution.txt. This file reads energyrate.txt, the cumulative fraction F(t) of the cell energies that the
// model cells release by the time t. F is linear in t between the rows of the file.

#include "energyinput.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <istream>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "inputfilestream.h"
#include "mpi_logging.h"
#include "packet.h"
#include "random.h"
#include "sn3d.h"

namespace energyinput {

static_assert(!THERMAL_ENERGY_INPUT_FROM_FILES || UNIFORM_PELLET_ENERGIES,
              "THERMAL_ENERGY_INPUT_FROM_FILES requires UNIFORM_PELLET_ENERGIES = true");

namespace {

// the rows of energyrate.txt: the times [s], and the cumulative fractions released by these times
std::vector<double> release_times;
std::vector<double> cumulative_fractions_released;

}  // anonymous namespace

void read_release_table() {
  int nrows = 0;
  if (globals::my_rank == 0) {
    auto filein = istream_required("energyrate.txt");
    if (!(filein >> nrows) || nrows < 1) {
      fatal_crash("energyrate.txt: the first line must give a positive number of rows");
    }
    release_times.resize(nrows);
    cumulative_fractions_released.resize(nrows);
    for (int row = 0; row < nrows; row++) {
      double time_days = 0.;
      double fraction = 0.;
      if (!(filein >> time_days >> fraction)) {
        fatal_crash("energyrate.txt: row {} of {} must give a time [d] and a cumulative fraction", row + 1, nrows);
      }
      release_times[row] = time_days * DAY;
      cumulative_fractions_released[row] = fraction;
      if (!std::isfinite(time_days) || time_days <= 0. || (row > 0 && release_times[row] <= release_times[row - 1])) {
        fatal_crash(
            "energyrate.txt: the time {:g} [d] of row {} must be positive and greater than the time of the previous "
            "row",
            time_days, row + 1);
      }
      if (!std::isfinite(fraction) || fraction < 0. || fraction > 1. ||
          (row > 0 && fraction < cumulative_fractions_released[row - 1])) {
        fatal_crash(
            "energyrate.txt: the cumulative fraction {:g} of row {} must be in [0, 1] and must not be less than the "
            "fraction of the previous row",
            fraction, row + 1);
      }
    }
    if (!(filein >> std::ws).eof()) {
      fatal_crash("energyrate.txt: the file has more data than the {} rows of the first line", nrows);
    }

    // a table with rounded values can end slightly below 1
    const double last_fraction = cumulative_fractions_released.back();
    if (std::abs(last_fraction - 1.) > 1e-3) {
      fatal_crash("energyrate.txt: the cumulative fraction of the last row is {:g}, but it must be 1", last_fraction);
    }
    for (auto& fraction : cumulative_fractions_released) {
      fraction /= last_fraction;
    }
    printlnlog("energyrate.txt: {} rows from {:g} to {:g} [d]. sn3d divides the fractions by the last fraction {:g}.",
               nrows, release_times.front() / DAY, release_times.back() / DAY, last_fraction);

    // The decay energy starts at the snapshot time of model.txt. A table that starts later has no release between
    // that time and its first row, and its first fraction is a release at the time of the first row.
    if (release_times.front() > grid::get_t_model() * (1. + 1e-6)) {
      printlnlog(
          "[warning] energyrate.txt: the first row at {:g} [d] is after the snapshot time {:g} [d] of model.txt. The "
          "table has no release between these times.",
          release_times.front() / DAY, grid::get_t_model() / DAY);
    }

    const double fraction_released_by_pellets = get_fraction_released_by_pellets();
    if (fraction_released_by_pellets <= 0.) {
      fatal_crash("energyrate.txt: the cells release no energy in the time range of the pellets");
    }
    printlnlog("energyrate.txt: the pellets release a fraction {:g} of the cell energies",
               fraction_released_by_pellets);
  }

  MPI_Bcast_safe(nrows, 0, MPI_COMM_WORLD);
  release_times.resize(nrows);
  cumulative_fractions_released.resize(nrows);
  MPI_Bcast_safe(release_times, 0, MPI_COMM_WORLD);
  MPI_Bcast_safe(cumulative_fractions_released, 0, MPI_COMM_WORLD);
}

auto get_cumulative_fraction_released(const double time_after_explosion) -> double {
  assert_always(!release_times.empty());
  // the release of the fraction of the first row occurs just after the time of that row, e.g. inside the first timestep
  if (time_after_explosion <= release_times.front()) {
    return 0.;
  }
  if (time_after_explosion >= release_times.back()) {
    return cumulative_fractions_released.back();
  }

  // the first row with a time after time_after_explosion, which lies between this row and the previous row
  const auto row = index_upperbound(release_times, time_after_explosion);
  const double fraction_of_row_interval =
      (time_after_explosion - release_times[row - 1]) / (release_times[row] - release_times[row - 1]);
  return std::lerp(cumulative_fractions_released[row - 1], cumulative_fractions_released[row],
                   fraction_of_row_interval);
}

auto get_released_fraction_after_expansion(const double time_after_explosion) -> double {
  assert_always(!release_times.empty());
  if (time_after_explosion <= release_times.front()) {
    return 0.;
  }

  // the integral of t' dF(t') from 0 to time_after_explosion. The fraction of the first row is a release at the time
  // of that row.
  double integral_release_time_dF = cumulative_fractions_released.front() * release_times.front();
  for (std::size_t row = 1; row < release_times.size() && release_times[row - 1] < time_after_explosion; row++) {
    const double t_low = release_times[row - 1];
    const double t_high = std::min(release_times[row], time_after_explosion);
    const double dF_dt = (cumulative_fractions_released[row] - cumulative_fractions_released[row - 1]) /
                         (release_times[row] - release_times[row - 1]);
    integral_release_time_dF += dF_dt * ((t_high * t_high) - (t_low * t_low)) / 2.;
  }
  return integral_release_time_dF / time_after_explosion;
}

auto get_release_rate_density(const int nonemptymgi, const int nts) -> double {
  const auto& timestep = globals::timesteps[nts];
  const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double fraction_released_in_timestep = get_cumulative_fraction_released(timestep.start + timestep.width) -
                                               get_cumulative_fraction_released(timestep.start);
  const double volume_tmid = grid::get_modelcell_assocvolume_tmin(mgi) * pow3(timestep.mid / globals::tmin);
  return grid::get_thermal_input_energy(mgi) * fraction_released_in_timestep / timestep.width / volume_tmid;
}

auto get_fraction_released_by_pellets() -> double {
  return get_cumulative_fraction_released(globals::tmax) -
         get_cumulative_fraction_released(INITIAL_PACKETS_ON ? 0. : globals::tmin);
}

void setup_thermal_input_pellet(const double e_cmf_per_packet, Packet& pkt) {
  // the inverse of F at a uniform sample of (F(first release time), F(tmax))
  const double first_release_time = INITIAL_PACKETS_ON ? 0. : globals::tmin;
  const double fraction_sample =
      std::lerp(get_cumulative_fraction_released(first_release_time), get_cumulative_fraction_released(globals::tmax),
                rng_uniform_pos(get_rngstate(pkt)));

  // the first row with a fraction at or above the sample. A rounded sample can be the last fraction, 1.
  const auto row = index_lowerbound(cumulative_fractions_released, fraction_sample);
  assert_always(row < std::ssize(cumulative_fractions_released));
  if (row == 0) {
    // the release of the fraction of the first row occurs just after the time of that row
    pkt.tdecay = std::nextafter(release_times.front(), globals::tmax);
  } else {
    const double fraction_of_row_interval =
        (fraction_sample - cumulative_fractions_released[row - 1]) /
        (cumulative_fractions_released[row] - cumulative_fractions_released[row - 1]);
    pkt.tdecay = std::lerp(release_times[row - 1], release_times[row], fraction_of_row_interval);
  }
  // A rounding error can give a time at or before the first release time. A release at the first release time would
  // take the path of a release before tmin, so the release time stays after it.
  pkt.tdecay = std::max(pkt.tdecay, std::nextafter(first_release_time, globals::tmax));

  pkt.type = TYPE_RADIOACTIVE_PELLET;
  pkt.e_cmf = e_cmf_per_packet;
  pkt.nu_cmf = e_cmf_per_packet / H;
  pkt.pellet_nucindex = -1;
  pkt.pellet_decaytype = -1;
  pkt.originated_from_particlenotgamma = false;
}

}  // namespace energyinput
