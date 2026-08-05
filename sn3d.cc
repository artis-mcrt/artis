// Main program of the sn3d simulation: time-dependent Monte Carlo radiative transfer for
// supernova and kilonova ejecta. Sets up the grid, packets, and MPI environment, then runs the
// timestep loop alternating update_grid() (solving temperatures, populations, and opacities
// from the previous timestep's estimators) with update_packets() (packet propagation), with
// periodic restart-file and output writing.

#include "sn3d.h"

#include <getopt.h>  // NOLINT(misc-include-cleaner): glibc declares getopt/optarg in bits/, no IWYU mapping to <getopt.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <ios>
#include <iterator>
#include <limits>
#include <numeric>
#include <print>
#include <string>
#include <utility>
#ifdef STDPAR_ON
#include <ranges>
#endif
#include <span>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "input.h"
#include "kpkt.h"
#include "macroatom.h"
#include "mpi_logging.h"
#include "nltepop.h"
#include "nonthermal.h"
#include "packet.h"
#include "radfield.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "spectrum_lightcurve.h"
#include "stats.h"
#include "update_grid.h"
#include "update_packets.h"
#include "version.h"
#include "vpkt.h"

namespace {

std::chrono::steady_clock::time_point real_time_start;
std::chrono::steady_clock::time_point
    time_timestep_start;  // this will be set after the first update of the grid and before packet prop
std::fstream estimators_file;

void setup_cellcache() {
  // When cellcache_singleslot is false, every non-empty cell gets its own persistent cache slot, shared by
  // all MPI ranks on the node. When it is true, each node rank gets a single reusable slot
  // (globals::cellcache[rank_in_node]) that is repopulated as its packets move between cells.
  const int num_cellcache_slots = cellcache_singleslot ? globals::node_nprocs : grid::get_nonempty_npts_model();
  reserve_resize(globals::cellcache, num_cellcache_slots);

  // per-cell sizes of each cache array
  const auto ncoolingterms = static_cast<size_t>(kpkt::ncoolingterms);
  const auto nincludedlevels = static_cast<size_t>(get_includedlevels());
  assert_always(globals::nbfcontinua >= 0);
  const auto nbfcontinua = static_cast<size_t>(globals::nbfcontinua);

  auto allphixstargetcount = 0ZU;
  auto chtransblocksize = 0ZU;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels(element, ion);

      for (int level = 0; level < nlevels; level++) {
        allphixstargetcount += get_nphixstargets(element, ion, level);

        const int ndowntrans = get_ndowntrans(element, ion, level);
        const int nuptrans = get_nuptrans(element, ion, level);
        chtransblocksize += ((2 * ndowntrans) + nuptrans);
      }
    }
  }
  assert_always(chtransblocksize <= std::numeric_limits<int>::max());

  const auto nslots = static_cast<size_t>(num_cellcache_slots);
  const auto maprocess_percell = nincludedlevels * MA_ACTION_COUNT;

  // one shared allocation per array, with a sub-range for each slot
  auto& backing = globals::cellcache_backing;
  backing.cooling_contrib.allocate(static_cast<ptrdiff_t>(ncoolingterms * nslots));
  backing.alllevels_pops.allocate(static_cast<ptrdiff_t>(nincludedlevels * nslots));
  backing.alllevels_maprocessrates.allocate(static_cast<ptrdiff_t>(maprocess_percell * nslots));
  backing.allmacroatomictransitions.allocate(static_cast<ptrdiff_t>(chtransblocksize * nslots));
  backing.allcont_modified_departureratios.allocate(static_cast<ptrdiff_t>(nbfcontinua * nslots));
  backing.allcont_nnlevel.allocate(static_cast<ptrdiff_t>(nbfcontinua * nslots));
  backing.allcont_keep.allocate(static_cast<ptrdiff_t>(nbfcontinua * nslots));
  backing.chi_ff_nnionpart.allocate(static_cast<ptrdiff_t>(nslots));
  backing.allphixstargets_corrphotoioncoeff.allocate(static_cast<ptrdiff_t>(allphixstargetcount * nslots));

  auto mem_usage_cellcache = 0ZU;
  for (int cellcachenum = 0; cellcachenum < num_cellcache_slots; cellcachenum++) {
    globals::CellCache& cacheslot = globals::cellcache[cellcachenum];

    // in multi-slot mode the slot index is the nonemptymgi (set here so that all node ranks agree on each
    // slot's cell even though the populate work is distributed); in single-slot mode the cell is not yet known
    cacheslot.nonemptymgi = cellcache_singleslot ? -1 : cellcachenum;

    const auto s = static_cast<size_t>(cellcachenum);
    cacheslot.cooling_contrib = backing.cooling_contrib.subspan(s * ncoolingterms, ncoolingterms);
    cacheslot.alllevels_pops = backing.alllevels_pops.subspan(s * nincludedlevels, nincludedlevels);
    cacheslot.alllevels_maprocessrates =
        backing.alllevels_maprocessrates.subspan(s * maprocess_percell, maprocess_percell);
    cacheslot.allmacroatomictransitions =
        backing.allmacroatomictransitions.subspan(s * chtransblocksize, chtransblocksize);
    cacheslot.allcont_modified_departureratios =
        backing.allcont_modified_departureratios.subspan(s * nbfcontinua, nbfcontinua);
    cacheslot.allcont_nnlevel = backing.allcont_nnlevel.subspan(s * nbfcontinua, nbfcontinua);
    cacheslot.allcont_keep = backing.allcont_keep.subspan(s * nbfcontinua, nbfcontinua);
    cacheslot.chi_ff_nnionpart = backing.chi_ff_nnionpart.subspan(s, 1);
    cacheslot.allphixstargets_corrphotoioncoeff =
        backing.allphixstargets_corrphotoioncoeff.subspan(s * allphixstargetcount, allphixstargetcount);

    if (cellcache_singleslot && cellcachenum == globals::rank_in_node) {
      // the lazy mutex-guarded rate calculation is only used in single-slot mode, and each rank only
      // ever accesses its own slot
      reserve_resize(cacheslot.cooling_contrib_locks, static_cast<size_t>(get_includedions()));
      std::ranges::fill(cacheslot.cooling_contrib_locks, 0);
      reserve_resize(cacheslot.allmacroatomictransitions_locks, nincludedlevels);
      std::ranges::fill(cacheslot.allmacroatomictransitions_locks, 0);

      for (size_t uniquelevelindex = 0; uniquelevelindex < nincludedlevels; uniquelevelindex++) {
        cacheslot.alllevels_maprocessrates[uniquelevelindex * MA_ACTION_COUNT] = -99.;
      }
    }

    if (cellcachenum == 0) {
      printlnlog("[info] mem_usage: cellcache coolinglist contribs for slot 0 occupies {:.3f} MB",
                 ncoolingterms * sizeof(double) / 1024. / 1024.);
      printlnlog("[info] mem_usage: cellcache macroatom stats for slot 0 occupies {:.3f} MB",
                 chtransblocksize * sizeof(double) / 1024. / 1024.);
      printlnlog("[info] mem_usage: cellcache for slot 0 occupies {:.3f} MB",
                 cacheslot.get_mem_usage() / 1024. / 1024.);
    }
    mem_usage_cellcache += cacheslot.get_mem_usage();
  }
  printlnlog("[info] mem_usage: cellcache for all {} slots occupies {:.3f} MB", num_cellcache_slots,
             mem_usage_cellcache / 1024. / 1024.);
}

void initialise_linestat_file() {
  if (globals::simulation_continued_from_saved) {
    // only write linestat.out on the first run
    return;
  }

  auto linestat_file = fstream_required("linestat.out", std::ios::out | std::ios::trunc);

  // with tens of millions of lines, per-value std::print calls to the stream are slow (each one re-checks whether
  // the stream is a terminal), so format into a buffer and write it out in large chunks
  std::string buffer;
  constexpr auto flushsize = 1UZ << 22U;
  buffer.reserve(flushsize + 64);
  const auto flush_if_full = [&linestat_file, &buffer]() {
    if (buffer.size() >= flushsize) {
      linestat_file.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
      buffer.clear();
    }
  };

  for (int i = 0; i < globals::nlines; i++) {
    std::format_to(std::back_inserter(buffer), "{:g} ", CLIGHT / globals::linelist.nu[i]);  // wavelength in cm
    flush_if_full();
  }
  buffer += '\n';

  for (int i = 0; i < globals::nlines; i++) {
    std::format_to(std::back_inserter(buffer), "{} ", get_atomicnumber(globals::linelist.elementindex[i]));
    flush_if_full();
  }
  buffer += '\n';

  for (int i = 0; i < globals::nlines; i++) {
    std::format_to(std::back_inserter(buffer), "{} ",
                   get_ionstage(globals::linelist.elementindex[i], globals::linelist.ionindex[i]));
    flush_if_full();
  }
  buffer += '\n';

  for (int i = 0; i < globals::nlines; i++) {
    const auto ionuniquelevelindexstart =
        get_ionuniquelevelindexstart(globals::linelist.elementindex[i], globals::linelist.ionindex[i]);
    const auto upper = globals::linelist.uniquelevelindex_upper[i] - ionuniquelevelindexstart;
    std::format_to(std::back_inserter(buffer), "{} ", (upper + 1));
    flush_if_full();
  }
  buffer += '\n';

  for (int i = 0; i < globals::nlines; i++) {
    const auto ionuniquelevelindexstart =
        get_ionuniquelevelindexstart(globals::linelist.elementindex[i], globals::linelist.ionindex[i]);
    const auto lower = globals::linelist.uniquelevelindex_lower[i] - ionuniquelevelindexstart;
    std::format_to(std::back_inserter(buffer), "{} ", (lower + 1));
    flush_if_full();
  }
  buffer += '\n';

  linestat_file.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
}

void write_deposition_file() {
  const int my_rank = globals::my_rank;
  const int nts = globals::timestep;
  printlog("Calculating and writing deposition.out...");
  const auto time_write_deposition_file_start = std::chrono::steady_clock::now();

  // Each analytic power [erg/s] is a sum over cells of cellmass * (rate coefficients dot the cell's initial
  // nuclide mass fractions), which factors into (rate coefficients) dot (total initial mass of each nuclide in
  // the gridded ejecta). The nuclide mass totals are time-independent, so compute them on the first call (each
  // rank sums its cell block, then a reduction) and reuse them for every timestep.
  static const auto [totmassnuclide_gridded, mtot] = [] {
    auto totmassnuclide = std::vector<double>(decay::get_num_nuclides(), 0.);
    double masstot = 0.;
    const int nstart_nonempty = grid::get_nstart_nonempty(globals::my_rank);
    const int ndo_nonempty = grid::get_ndo_nonempty(globals::my_rank);
    for (int nonemptymgi = nstart_nonempty; nonemptymgi < (nstart_nonempty + ndo_nonempty); nonemptymgi++) {
      const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
      const double cellmass = grid::get_rho_tmin(mgi) * grid::get_modelcell_assocvolume_tmin(mgi);
      masstot += cellmass;
      for (int nucindex = 0; nucindex < decay::get_num_nuclides(); nucindex++) {
        totmassnuclide[nucindex] += cellmass * grid::get_modelinitnucmassfrac(mgi, nucindex);
      }
    }
    MPI_Allreduce_safe(std::span{totmassnuclide}, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce_safe(masstot, MPI_SUM, MPI_COMM_WORLD);
    return std::pair{std::move(totmassnuclide), masstot};
  }();

  // The decay chain calculations depend only on the timestep midpoint, so the per-source-nuclide rate
  // coefficients are computed once and each power is a single dot product with the nuclide mass totals
  const double t_mid_nts = globals::timesteps[nts].mid;
  const auto nuc_massfrac_coeffs = decay::calc_nuc_massfrac_coeffs(t_mid_nts);
  const auto emissionratecoeffs = decay::calc_ana_emission_ratecoeffs(nuc_massfrac_coeffs);
  const auto qdot_ratecoeffs = decay::calc_qdot_ratecoeffs(nuc_massfrac_coeffs);

  const auto total_power = [&](const std::vector<double>& ratecoeffs_per_source) {
    return std::transform_reduce(ratecoeffs_per_source.cbegin(), ratecoeffs_per_source.cend(),
                                 totmassnuclide_gridded.cbegin(), 0.);
  };

  // power in [erg/s]
  globals::timesteps[nts].eps_positron_ana_power = total_power(emissionratecoeffs.positron);
  globals::timesteps[nts].eps_electron_ana_power = total_power(emissionratecoeffs.electron);
  globals::timesteps[nts].eps_alpha_ana_power = total_power(emissionratecoeffs.alpha);
  globals::timesteps[nts].eps_spfission_ana_power = total_power(emissionratecoeffs.spfission);
  globals::timesteps[nts].qdot_betaminus = total_power(qdot_ratecoeffs[decay::DECAYTYPE_BETAMINUS]);
  globals::timesteps[nts].qdot_alpha = total_power(qdot_ratecoeffs[decay::DECAYTYPE_ALPHA]);
  globals::timesteps[nts].qdot_spfission = total_power(qdot_ratecoeffs[decay::DECAYTYPE_SPONTFISSION]);
  globals::timesteps[nts].qdot_total = 0.;
  for (const auto decaytype : decay::all_decaytypes) {
    globals::timesteps[nts].qdot_total += total_power(qdot_ratecoeffs[decaytype]);
  }

  if (my_rank == 0) {
    const bool any_fission = decay::decaytype_is_used(decay::DECAYTYPE_SPONTFISSION);

    // each deposition.out column after the timestep number: the header name paired with the function computing
    // the value for one timestep row, so that the header and the printed values cannot fall out of alignment
    struct DepositionColumn {
      const char* header = "";
      std::function<double(const globals::TimeStep&)> value;
      bool enabled = true;
    };
    const auto columns = std::vector<DepositionColumn>{
        {.header = "tmid_days", .value = [](const auto& ts) { return ts.mid / DAY; }},
        {.header = "tmid_s", .value = [](const auto& ts) { return ts.mid; }},
        {.header = "total_dep_Lsun",
         .value =
             [](const auto& ts) {
               return (ts.gamma_dep + ts.positron_dep + ts.electron_dep + ts.alpha_dep + ts.spfission_dep_discrete) /
                      ts.width / LSUN;
             }},
        {.header = "gammadep_discrete_Lsun",
         .value = [](const auto& ts) { return ts.gamma_dep_discrete / ts.width / LSUN; }},
        {.header = "gammadep_Lsun", .value = [](const auto& ts) { return ts.gamma_dep / ts.width / LSUN; }},
        {.header = "positrondep_Lsun", .value = [](const auto& ts) { return ts.positron_dep / ts.width / LSUN; }},
        {.header = "eps_positron_ana_Lsun", .value = [](const auto& ts) { return ts.eps_positron_ana_power / LSUN; }},
        {.header = "elecdep_Lsun", .value = [](const auto& ts) { return ts.electron_dep / ts.width / LSUN; }},
        {.header = "eps_elec_Lsun", .value = [](const auto& ts) { return ts.electron_emission / ts.width / LSUN; }},
        {.header = "eps_elec_ana_Lsun", .value = [](const auto& ts) { return ts.eps_electron_ana_power / LSUN; }},
        {.header = "alphadep_Lsun", .value = [](const auto& ts) { return ts.alpha_dep / ts.width / LSUN; }},
        {.header = "eps_alpha_Lsun", .value = [](const auto& ts) { return ts.alpha_emission / ts.width / LSUN; }},
        {.header = "eps_alpha_ana_Lsun", .value = [](const auto& ts) { return ts.eps_alpha_ana_power / LSUN; }},
        {.header = "eps_spfission_ana_Lsun",
         .value = [](const auto& ts) { return ts.eps_spfission_ana_power / LSUN; },
         .enabled = any_fission},
        {.header = "eps_gamma_Lsun", .value = [](const auto& ts) { return ts.gamma_emission / ts.width / LSUN; }},
        {.header = "Qdot_betaminus_ana_erg/s/g", .value = [&](const auto& ts) { return ts.qdot_betaminus / mtot; }},
        {.header = "Qdotalpha_ana_erg/s/g", .value = [&](const auto& ts) { return ts.qdot_alpha / mtot; }},
        {.header = "Qdotspfission_ana_erg/s/g",
         .value = [&](const auto& ts) { return ts.qdot_spfission / mtot; },
         .enabled = any_fission},
        {.header = "eps_erg/s/g",
         .value =
             [&](const auto& ts) {
               return (ts.gamma_emission + ts.positron_emission + ts.electron_emission + ts.alpha_emission +
                       ts.spfission_dep_discrete) /
                      mtot / ts.width;
             }},
        {.header = "Qdot_ana_erg/s/g", .value = [&](const auto& ts) { return ts.qdot_total / mtot; }},
        {.header = "positrondep_discrete_Lsun",
         .value = [](const auto& ts) { return ts.positron_dep_discrete / ts.width / LSUN; }},
        {.header = "elecdep_discrete_Lsun",
         .value = [](const auto& ts) { return ts.electron_dep_discrete / ts.width / LSUN; }},
        {.header = "alphadep_discrete_Lsun",
         .value = [](const auto& ts) { return ts.alpha_dep_discrete / ts.width / LSUN; }},
        {.header = "spfission_dep_discrete_Lsun",
         .value = [](const auto& ts) { return ts.spfission_dep_discrete / ts.width / LSUN; },
         .enabled = any_fission},
    };

    auto dep_file = fstream_required("deposition.out.tmp", std::ios::out | std::ios::trunc);
    std::print(dep_file, "#ts");
    for (const auto& column : columns) {
      if (column.enabled) {
        std::print(dep_file, " {}", column.header);
      }
    }
    std::println(dep_file, "");

    for (int i = 0; i <= nts; i++) {
      const auto& ts = globals::timesteps[i];
      if (!any_fission) {
        assert_testmodeonly(ts.eps_spfission_ana_power == 0.);
        assert_testmodeonly(ts.qdot_spfission == 0.);
      }
      std::print(dep_file, "{}", i);
      for (const auto& column : columns) {
        if (column.enabled) {
          std::print(dep_file, " {:g}", column.value(ts));
        }
      }
      std::println(dep_file, "");
    }
    dep_file.close();

    std::remove("deposition.out");
    std::rename("deposition.out.tmp", "deposition.out");

    // energy-conservation consistency check (log only): the cumulative deposition should not exceed the
    // cumulative decay emission by more than the Monte Carlo noise of the trajectory estimators allows
    double cumulative_emission = 0.;
    double cumulative_dep = 0.;
    for (int i = 0; i <= nts; i++) {
      cumulative_emission += globals::timesteps[i].gamma_emission + globals::timesteps[i].positron_emission +
                             globals::timesteps[i].electron_emission + globals::timesteps[i].alpha_emission +
                             globals::timesteps[i].spfission_dep_discrete;
      cumulative_dep += globals::timesteps[i].gamma_dep + globals::timesteps[i].positron_dep +
                        globals::timesteps[i].electron_dep + globals::timesteps[i].alpha_dep +
                        globals::timesteps[i].spfission_dep_discrete;
    }
    if (cumulative_emission > 0. && cumulative_dep > (1.2 * cumulative_emission)) {
      printlnlog(
          "[warning] energy conservation check: cumulative deposition {:g} [erg] exceeds cumulative decay emission "
          "{:g} [erg] by more than 20 percent, which is more than Monte Carlo noise should produce",
          cumulative_dep, cumulative_emission);
    }
  }

  const auto deposition_write_duration =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - time_write_deposition_file_start).count();
  printlnlog("took {:.1f} seconds", deposition_write_duration);
}

void write_timestep_file() {
  auto timestepfile = fstream_required("timesteps.out", std::ofstream::out | std::ofstream::trunc);
  std::print(timestepfile, "#timestep tstart_days tmid_days twidth_days\n");
  for (int n = 0; n < globals::ntimesteps; n++) {
    std::println(timestepfile, "{} {:g} {:g} {:g}", n, globals::timesteps[n].start / DAY,
                 globals::timesteps[n].mid / DAY, globals::timesteps[n].width / DAY);
  }
}

void mpi_communicate_grid_properties() {
  const int nincludedions = get_includedions();
  const auto nelements = get_nelements();
  for (int root = 0; root < globals::nprocs; root++) {
    MPI_Barrier_allranks();

    int root_node_id = globals::node_id;
    MPI_Bcast_safe(root_node_id, root, MPI_COMM_WORLD);

    const ptrdiff_t root_ndo_nonempty = grid::get_ndo_nonempty(root);
    if (root_ndo_nonempty <= 0) {
      continue;
    }
    const ptrdiff_t root_nstart_nonempty = grid::get_nstart_nonempty(root);
    assert_always(root_nstart_nonempty >= 0);
    // the highest index accessed below is (nstart + ndo - 1), which must be a valid index
    assert_always((root_nstart_nonempty + root_ndo_nonempty - 1) < grid::get_nonempty_npts_model());

    if (USE_LUT_PHOTOION && globals::nbfcontinua_ground > 0) {
      MPI_Bcast_safe(std::span{globals::gammaestimator}.subspan(root_nstart_nonempty * globals::nbfcontinua_ground,
                                                                root_ndo_nonempty * globals::nbfcontinua_ground),
                     root, MPI_COMM_WORLD);
    }

    for (auto nonemptymgi = root_nstart_nonempty; nonemptymgi < (root_nstart_nonempty + root_ndo_nonempty);
         nonemptymgi++) {
      radfield::do_MPI_Bcast(nonemptymgi, root, root_node_id);

      nonthermal::nt_MPI_Bcast(nonemptymgi, root_node_id);

      MPI_Bcast_binned_opacities(nonemptymgi, root_node_id);
    }

    MPI_Barrier_allranks();
    if (globals::rank_in_node == 0) {
      MPI_Bcast_safe(grid::rho_allcells.subspan(root_nstart_nonempty, root_ndo_nonempty), root_node_id,
                     globals::mpi_comm_internode);
      MPI_Bcast_safe(grid::Te_allcells.subspan(root_nstart_nonempty, root_ndo_nonempty), root_node_id,
                     globals::mpi_comm_internode);
      MPI_Bcast_safe(grid::TJ_allcells.subspan(root_nstart_nonempty, root_ndo_nonempty), root_node_id,
                     globals::mpi_comm_internode);
      MPI_Bcast_safe(grid::TR_allcells.subspan(root_nstart_nonempty, root_ndo_nonempty), root_node_id,
                     globals::mpi_comm_internode);
      MPI_Bcast_safe(grid::W_allcells.subspan(root_nstart_nonempty, root_ndo_nonempty), root_node_id,
                     globals::mpi_comm_internode);
      MPI_Bcast_safe(grid::nne_allcells.subspan(root_nstart_nonempty, root_ndo_nonempty), root_node_id,
                     globals::mpi_comm_internode);
      MPI_Bcast_safe(grid::nnetot_allcells.subspan(root_nstart_nonempty, root_ndo_nonempty), root_node_id,
                     globals::mpi_comm_internode);
      MPI_Bcast_safe(grid::kappagrey_allcells.subspan(root_nstart_nonempty, root_ndo_nonempty), root_node_id,
                     globals::mpi_comm_internode);
      MPI_Bcast_safe(grid::grey_depth_allcells.subspan(root_nstart_nonempty, root_ndo_nonempty), root_node_id,
                     globals::mpi_comm_internode);
      MPI_Bcast_safe(grid::thick_allcells.subspan(root_nstart_nonempty, root_ndo_nonempty), root_node_id,
                     globals::mpi_comm_internode);
      if constexpr (USE_MICROCLUMPING) {
        MPI_Bcast_safe(grid::clumpfactor_allcells.subspan(root_nstart_nonempty, root_ndo_nonempty), root_node_id,
                       globals::mpi_comm_internode);
      }

      if (USE_LUT_PHOTOION && globals::nbfcontinua_ground > 0) {
        assert_always(globals::corrphotoionrenorm.data() != nullptr);
        MPI_Bcast_safe(globals::corrphotoionrenorm.subspan(root_nstart_nonempty * globals::nbfcontinua_ground,
                                                           root_ndo_nonempty * globals::nbfcontinua_ground),
                       root_node_id, globals::mpi_comm_internode);
      }

      if (globals::total_nlte_levels > 0) {
        MPI_Bcast_safe(nltepops_allcells.subspan(root_nstart_nonempty * globals::total_nlte_levels,
                                                 root_ndo_nonempty * globals::total_nlte_levels),
                       root_node_id, globals::mpi_comm_internode);
      }

      MPI_Bcast_safe(
          grid::elem_meanweight_allcells.subspan(root_nstart_nonempty * nelements, root_ndo_nonempty * nelements),
          root_node_id, globals::mpi_comm_internode);

      MPI_Bcast_safe(
          grid::elem_massfracs_allcells.subspan(root_nstart_nonempty * nelements, root_ndo_nonempty * nelements),
          root_node_id, globals::mpi_comm_internode);

      MPI_Bcast_safe(grid::ion_groundlevelpops_allcells.subspan(root_nstart_nonempty * nincludedions,
                                                                root_ndo_nonempty * nincludedions),
                     root_node_id, globals::mpi_comm_internode);

      MPI_Bcast_safe(grid::ion_partfuncts_allcells.subspan(root_nstart_nonempty * nincludedions,
                                                           root_ndo_nonempty * nincludedions),
                     root_node_id, globals::mpi_comm_internode);

      MPI_Bcast_safe(kpkt::ion_cooling_contribs_allcells.subspan(root_nstart_nonempty * nincludedions,
                                                                 root_ndo_nonempty * nincludedions),
                     root_node_id, globals::mpi_comm_internode);
    }

    MPI_Barrier_allranks();
  }

  MPI_Barrier_allranks();
}

void normalise_deposition_estimators(int nts) {
  const double dt = globals::timesteps[nts].width;
  const auto nprocs = globals::nprocs;

  globals::timesteps[nts].gamma_dep = 0.;
  globals::timesteps[nts].positron_dep = 0.;
  globals::timesteps[nts].electron_dep = 0.;
  globals::timesteps[nts].alpha_dep = 0.;

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);

    const double dV = grid::get_modelcell_assocvolume_tmin(mgi) * pow3(globals::timesteps[nts].mid / globals::tmin);

    // contribute the energy deposited (in erg) by each process in this cell to the timestep total
    globals::timesteps[nts].gamma_dep += globals::dep_estimator_gamma[nonemptymgi] / nprocs;
    globals::timesteps[nts].positron_dep += globals::dep_estimator_positron[nonemptymgi] / nprocs;
    globals::timesteps[nts].electron_dep += globals::dep_estimator_electron[nonemptymgi] / nprocs;
    globals::timesteps[nts].alpha_dep += globals::dep_estimator_alpha[nonemptymgi] / nprocs;

    // normalise the estimators to units of erg/s/cm^3
    const double estimator_normfactor = 1 / dV / dt / nprocs;

    globals::dep_estimator_gamma[nonemptymgi] *= estimator_normfactor;
    globals::dep_estimator_positron[nonemptymgi] *= estimator_normfactor;
    globals::dep_estimator_electron[nonemptymgi] *= estimator_normfactor;
    globals::dep_estimator_alpha[nonemptymgi] *= estimator_normfactor;

    assert_testmodeonly(globals::dep_estimator_gamma[nonemptymgi] >= 0.);
    assert_testmodeonly(std::isfinite(globals::dep_estimator_gamma[nonemptymgi]));
  }
}

void mpi_reduce_estimators(const int nts) {
  const int nonempty_npts_model = grid::get_nonempty_npts_model();
  radfield::reduce_estimators();
  MPI_Barrier_allranks();
  MPI_Allreduce_safe(globals::ffheatingestimator, MPI_SUM, MPI_COMM_WORLD);
  if constexpr (!DIRECT_COL_HEAT) {
    MPI_Allreduce_safe(globals::colheatingestimator, MPI_SUM, MPI_COMM_WORLD);
  }
  MPI_Barrier_allranks();

  if (globals::nbfcontinua_ground > 0) {
    if constexpr (USE_LUT_PHOTOION) {
      MPI_Barrier_allranks();
      MPI_Allreduce_safe(globals::gammaestimator, MPI_SUM, MPI_COMM_WORLD);
    }

    if constexpr (USE_ION_BFHEATING_ESTIMATORS) {
      MPI_Barrier_allranks();
      MPI_Allreduce_safe(globals::bfheatingestimator, MPI_SUM, MPI_COMM_WORLD);
    }
  }

  assert_always(std::ssize(globals::dep_estimator_gamma) == nonempty_npts_model);
  MPI_Allreduce_safe(globals::dep_estimator_gamma, MPI_SUM, MPI_COMM_WORLD);
  assert_always(std::ssize(globals::dep_estimator_positron) == nonempty_npts_model);
  MPI_Allreduce_safe(globals::dep_estimator_positron, MPI_SUM, MPI_COMM_WORLD);
  assert_always(std::ssize(globals::dep_estimator_electron) == nonempty_npts_model);
  MPI_Allreduce_safe(globals::dep_estimator_electron, MPI_SUM, MPI_COMM_WORLD);
  assert_always(std::ssize(globals::dep_estimator_alpha) == nonempty_npts_model);
  MPI_Allreduce_safe(globals::dep_estimator_alpha, MPI_SUM, MPI_COMM_WORLD);

  MPI_Barrier_allranks();

  MPI_Allreduce_safe(globals::timesteps[nts].gamma_dep_discrete, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].gamma_dep_discrete /= globals::nprocs;

  MPI_Allreduce_safe(globals::timesteps[nts].positron_dep_discrete, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].positron_dep_discrete /= globals::nprocs;

  MPI_Allreduce_safe(globals::timesteps[nts].positron_emission, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].positron_emission /= globals::nprocs;

  MPI_Allreduce_safe(globals::timesteps[nts].electron_dep_discrete, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].electron_dep_discrete /= globals::nprocs;

  MPI_Allreduce_safe(globals::timesteps[nts].electron_emission, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].electron_emission /= globals::nprocs;

  MPI_Allreduce_safe(globals::timesteps[nts].alpha_dep_discrete, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].alpha_dep_discrete /= globals::nprocs;

  MPI_Allreduce_safe(globals::timesteps[nts].alpha_emission, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].alpha_emission /= globals::nprocs;

  MPI_Allreduce_safe(globals::timesteps[nts].spfission_dep_discrete, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].spfission_dep_discrete /= globals::nprocs;

  MPI_Allreduce_safe(globals::timesteps[nts].gamma_emission, MPI_SUM, MPI_COMM_WORLD);
  globals::timesteps[nts].gamma_emission /= globals::nprocs;

  MPI_Barrier_allranks();

  // The estimators have been summed across all processes and distributed.
  // They will now be normalised independently on all processes.

  normalise_deposition_estimators(nts);
}

auto walltime_sufficient_to_continue(const int nts, const int nts_prev, const int walltimelimitseconds) -> bool {
  MPI_Barrier_allranks();
  // time is measured from just before packet propagation from one timestep to the next
  const auto estimated_time_per_timestep =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - time_timestep_start).count();
  printlnlog("TIME: time between timesteps is {:.1f} seconds (measured packet prop of ts {} and update grid of ts {})",
             estimated_time_per_timestep, nts_prev, nts);

  bool do_this_full_loop = true;
  if (walltimelimitseconds > 0 && nts < globals::timestep_finish) {
    const auto wallclock_used_seconds =
        std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - real_time_start).count();
    const auto wallclock_remaining_seconds = walltimelimitseconds - wallclock_used_seconds;
    printlnlog("TIMED_RESTARTS: Used {} of {} seconds of wall time.", wallclock_used_seconds, walltimelimitseconds);

    // This flag being false will make it update_grid, and then exit
    do_this_full_loop = (wallclock_remaining_seconds >= (1.5 * estimated_time_per_timestep));

    // communicate whatever decision the rank 0 process decided, just in case they differ
    MPI_Bcast_safe(do_this_full_loop, 0, MPI_COMM_WORLD);
    if (do_this_full_loop) {
      printlnlog("TIMED_RESTARTS: Going to continue since remaining time {} s >= 1.5 * time_per_timestep",
                 wallclock_remaining_seconds);
    } else {
      printlnlog("TIMED_RESTARTS: Going to terminate since remaining time {} s < 1.5 * time_per_timestep",
                 wallclock_remaining_seconds);
    }
  }
  return do_this_full_loop;
}

void save_grid_and_packets(const int nts, std::vector<Packet>& packets) {
  MPI_Barrier_allranks();
  const auto my_rank = globals::my_rank;

  const auto time_write_packets_file_start = std::chrono::steady_clock::now();

  if constexpr (!KEEP_ESCAPED_GAMMAS) {
    std::erase_if(packets, [](const Packet& pkt) { return pkt.type == TYPE_ESCAPE && pkt.escape_type == TYPE_GAMMA; });
  }

  // save packet state at start of current timestep (before propagation)
  write_temp_packetsfile(nts, globals::my_rank, packets);

  vpkt::write_timestep(nts, false);

  const auto time_write_packets_finished_thisrank = std::chrono::steady_clock::now();

  MPI_Barrier_allranks();

  const auto timenow = std::chrono::steady_clock::now();
  const auto packets_write_time =
      std::chrono::duration<double>(time_write_packets_finished_thisrank - time_write_packets_file_start).count();
  const auto packets_wait_time = std::chrono::duration<double>(timenow - time_write_packets_finished_thisrank).count();
  const auto packets_total_time = std::chrono::duration<double>(timenow - time_write_packets_file_start).count();

  printlnlog("timestep {}: finished writing temporary packets file (took {:.1f}s, waited {:.1f}s, total {:.1f}s)", nts,
             packets_write_time, packets_wait_time, packets_total_time);

  if (my_rank == 0) {
    grid::write_grid_restart_data(nts);
    update_parameterfile(nts);
  }

  if (!KEEP_ALL_RESTART_FILES) {
    // ensure new packets files have been written by all processes before we remove the old set
    MPI_Barrier_allranks();

    if (my_rank == 0) {
      const auto filename_prev_gridsave = std::format("gridsave_ts{}.tmp", nts - 1);
      if (std::filesystem::remove(filename_prev_gridsave)) {
        printlnlog("deleted {}", filename_prev_gridsave);
      }
    }

    // delete temp packets files from previous timestep now that all restart data for the new timestep is available
    const auto filename_prev_packetstmp = std::format("packets_{:04d}_ts{:d}.tmp", my_rank, nts - 1);
    if (std::filesystem::remove(filename_prev_packetstmp)) {
      printlnlog("deleted {}", filename_prev_packetstmp);
    }

    vpkt::remove_temp_vpkt_file(nts - 1, my_rank);
  }
}

void zero_estimators() {
  MPI_Barrier_allranks();
  radfield::zero_estimators();

  std::ranges::fill(globals::ffheatingestimator, 0.);
  std::ranges::fill(globals::colheatingestimator, 0.);
  std::ranges::fill(globals::dep_estimator_gamma, 0.);
  std::ranges::fill(globals::dep_estimator_positron, 0.);
  std::ranges::fill(globals::dep_estimator_electron, 0.);
  std::ranges::fill(globals::dep_estimator_alpha, 0.);

  if constexpr (USE_LUT_PHOTOION) {
    if (globals::nbfcontinua_ground > 0) {
      std::ranges::fill(globals::gammaestimator, 0.);
    }
  }

  if constexpr (USE_ION_BFHEATING_ESTIMATORS) {
    if (globals::nbfcontinua_ground > 0) {
      std::ranges::fill(globals::bfheatingestimator, 0.);
    }
  }

  MPI_Barrier_allranks();
}

auto do_timestep(const int nts, const int titer, std::vector<Packet>& packets, const int walltimelimitseconds) -> bool {
  bool do_this_full_loop = true;
  const int nts_prev = (titer != 0 || nts == 0) ? nts : nts - 1;
  if ((titer > 0) || (globals::simulation_continued_from_saved && (nts == globals::timestep_initial))) {
    // Read the packets file to reset before each additional iteration on the timestep
    read_temp_packetsfile(nts, globals::my_rank, packets);
  }

  // Some counters on pkt-actions need to be reset to do statistics
  stats::pkt_action_counters_reset();

  if (nts == 0) {
    radfield::initialise_prev_titer_photoionestimators();
  }

  // Update the matter quantities in the grid for the new timestep.

  update_grid(estimators_file, nts, nts_prev, titer, real_time_start);

  const auto sys_time_start_communicate_grid = std::chrono::steady_clock::now();

  // Each process has now updated its own set of cells. The results now need to be communicated between processes.
  mpi_communicate_grid_properties();

  const auto communicate_grid_duration =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - sys_time_start_communicate_grid).count();
  printlnlog("timestep {}: time after grid properties have been communicated (took {:.1f} seconds)", nts,
             communicate_grid_duration);

  // If this is not the 0th time step of the current job step,
  // write out a snapshot of the grid properties for further restarts and update input.txt accordingly
  if (((nts - globals::timestep_initial) != 0)) {
    save_grid_and_packets(nts, packets);
    do_this_full_loop = walltime_sufficient_to_continue(nts, nts_prev, walltimelimitseconds);
  }
  time_timestep_start = std::chrono::steady_clock::now();

  // set all the estimators to zero before moving packets. This is done after update_grid() so that the
  // gamma-ray heating estimator, and the photoionisation and stimulated recombination estimators, are still
  // available to it from the previous timestep.
  zero_estimators();

  MPI_Barrier_allranks();
  if ((nts < globals::timestep_finish) && do_this_full_loop) {
    // Now process the packets.

    update_packets(nts, packets);

    // All the processes have their own versions of the estimators for this time step now.
    // Since these are going to be needed in the next time step, we will gather all the
    // estimators together now, sum them, and distribute the results

    const auto time_communicate_estimators_start = std::chrono::steady_clock::now();
    mpi_reduce_estimators(nts);

    const auto communicate_estimators_duration =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - time_communicate_estimators_start).count();
    printlnlog("timestep {}: time after estimators have been communicated (took {:.1f} seconds)", nts,
               communicate_estimators_duration);

    write_deposition_file();

    write_partial_lightcurve_spectra(nts, packets);

    printlnlog("During timestep {} on MPI process {}, {} pellets decayed and {} packets escaped. (t={:g} [d])", nts,
               globals::my_rank, globals::timesteps[nts].pellet_decays, stats::get_counter(stats::Counter::PKTESCAPES),
               globals::timesteps[nts].mid / DAY);

    if (VPKT_ON) {
      printlnlog("During timestep {} on MPI process {}, {} virtual packets were generated and {} escaped.", nts,
                 globals::my_rank, vpkt::nvpkt_created,
                 vpkt::nvpkt_esc_from_rpkt + vpkt::nvpkt_esc_from_kpkt + vpkt::nvpkt_esc_from_macroatom);
      printlnlog(
          "{} virtual packets came from an electron scattering event, {} from a kpkt deactivation and {} from a "
          "macroatom deactivation.",
          vpkt::nvpkt_esc_from_rpkt, vpkt::nvpkt_esc_from_kpkt, vpkt::nvpkt_esc_from_macroatom);

      vpkt::nvpkt_created = 0;
      vpkt::nvpkt_esc_from_rpkt = 0;
      vpkt::nvpkt_esc_from_kpkt = 0;
      vpkt::nvpkt_esc_from_macroatom = 0;
    }

    if (nts == (globals::timestep_finish - 1)) {
      const auto filename = std::format("packets{:02d}_{:04d}.out", 0, globals::my_rank);
      write_text_packets(filename, packets);

      vpkt::write_timestep(nts, true);

      const auto after_final_packets_write =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - real_time_start).count();
      printlnlog("time after write final packets file (tstart + {:.1f} seconds)", after_final_packets_write);
    }
  }
  return !do_this_full_loop;
}

}  // anonymous namespace

auto main(int argc, char* argv[]) -> int {
  real_time_start = std::chrono::steady_clock::now();

  // if DETAILED_BF_ESTIMATORS_ON is true, USE_LUT_PHOTOION must be false
  assert_always(!DETAILED_BF_ESTIMATORS_ON || !USE_LUT_PHOTOION);

  if constexpr (VPKT_ON) {
    vpkt::nvpkt_created = 0;
    vpkt::nvpkt_esc_from_rpkt = 0;
    vpkt::nvpkt_esc_from_kpkt = 0;
    vpkt::nvpkt_esc_from_macroatom = 0;
  }

  MPI_Init(&argc, &argv);

  globals::setup_mpi_vars();

  check_already_running();

#ifdef STDPAR_ON
  for (int t = 1; t < get_max_threads(); t++) {
    std::filesystem::remove(std::format("output_{}-{}.txt", globals::my_rank, t));
  }
#endif

#if defined(_OPENMP) && !defined(GPU_ON)
  // Explicitly turn off dynamic threads. The per-thread log file handles in mpi_logging.h are threadprivate,
  // and OpenMP only guarantees that threadprivate data persists between parallel regions while the team size
  // stays fixed, which dynamic adjustment would break.
  omp_set_dynamic(0);
#pragma omp parallel
#endif
  {
    // initialise the thread and rank specific output file
    set_log_file(std::format("output_{}-{}.txt", globals::my_rank, get_thread_num()));

#ifdef _OPENMP
    printlnlog("OpenMP parallelisation is active with {} threads (max {})", omp_get_num_threads(), get_max_threads());
#else
    printlnlog("OpenMP parallelisation is not enabled in this build (this is normal)");
#endif
  }

#ifdef STDPAR_ON
  printlog("C++ standard parallelism (stdpar) is enabled");
#ifdef GPU_ON
  printlnlog(" with GPU support");
#else
  printlnlog(" with up to {} CPU threads", get_max_threads());
#endif
#endif

#ifdef GPU_ON
  printlnlog("GPU_ON is enabled");
#endif

  printlnlog("sn3d launched (this line marks the start time)");

  printlog("Integration method is: ");

#ifdef USE_SIMPSON_INTEGRATOR
  printlnlog("Simpson rule");
#else
  printlnlog("Boost Gauss-Kronrod quadrature");
#endif

#ifdef WALLTIMELIMITSECONDS
  int walltimelimitseconds = WALLTIMELIMITSECONDS;
#else
  int walltimelimitseconds = -1;
#endif

  int opt = 0;
  while ((opt = getopt(argc, argv, "w:")) != -1) {  // NOLINT(concurrency-mt-unsafe,misc-include-cleaner)
    if (opt == 'w') {
      char* parse_end = nullptr;
      const float walltimehours = strtof(optarg, &parse_end);  // NOLINT(misc-include-cleaner)
      if (parse_end == optarg || *parse_end != '\0' || !std::isfinite(walltimehours) || walltimehours <= 0.) {
        // silently accepting a bad value would disable the wall time limit instead of applying it
        std::println(stderr, "[error] invalid wall time hours '{}' given with -w option", optarg);
        std::abort();
      }
      walltimelimitseconds = static_cast<int>(walltimehours * HOUR);
      printlnlog("command line argument specifies wall time hours '{}', so setting walltimelimitseconds = {}", optarg,
                 walltimelimitseconds);
    } else {
      std::println(stderr, "Usage: {} [-w WALLTIMELIMITHOURS]",
                   argv[0]);  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic,)
      std::abort();
    }
  }

  std::vector<Packet> packets;
  reserve_resize(packets, MPKTS);

  printlnlog("git branch: {}", GIT_BRANCH);

  printlnlog("git version: {}", GIT_VERSION);

  printlnlog("git status: {}", GIT_STATUS);

  printlnlog("sn3d compiled at {} on {}", __TIME__, __DATE__);

#if defined TESTMODE && TESTMODE
  printlnlog("TESTMODE is ON");
#endif

  printlnlog("process id (pid): {}", getpid());
  printlnlog("MPI enabled:");
  printlnlog("  rank {} of [0..{}] in MPI_COMM_WORLD", globals::my_rank, globals::nprocs - 1);
  printlnlog("  rank {} of [0..{}] in node {} of [0..{}]", globals::rank_in_node, globals::node_nprocs - 1,
             globals::node_id, globals::node_count - 1);
#ifdef MAX_NODE_SIZE
  printlnlog(
      "[warning] Compiled with MAX_NODE_SIZE {}, which may mean that there are more nodes reported than physically "
      "present",
      MAX_NODE_SIZE);
#endif

  // Read in parameters from input.txt
  read_parameterfile(packets);

  // Read in parameters from vpkt.txt
  if constexpr (VPKT_ON) {
    vpkt::read_vpktparameterfile();
  }

  read_atomicdata();

  grid::read_ejecta_model();

  if (globals::simulation_continued_from_saved) {
    assert_always(globals::nprocs_exspec == globals::nprocs);
  } else {
    globals::nprocs_exspec = globals::nprocs;
  }

  if (globals::my_rank == 0) {
    initialise_linestat_file();
  }

  printlnlog("timesteps {}", globals::ntimesteps);

  // Precalculate the rate coefficients for spontaneous and stimulated recombination
  // and for photoionisation. With the nebular approximation they only depend on T_e
  // T_R and W. W is easily factored out. For stimulated recombination we must assume
  // T_e = T_R for this precalculation.

  ratecoefficients_init();

  MPI_Barrier_allranks();

  // Record the chosen syn_dir (only one rank writes it, since every rank would write the same file)
  if (globals::my_rank == 0) {
    auto syn_file = fstream_required("syn_dir.txt", std::ios::out | std::ios::trunc);
    std::print(syn_file, "{} {} {}", syn_dir[0], syn_dir[1], syn_dir[2]);
    syn_file.close();
  }

  bool terminate_early = false;

  setup_timesteps();

  if (globals::my_rank == 0) {
    write_timestep_file();
  }

  grid::init_grid();

  printlnlog("Simulation propagates {:g} packets per process (total {:g} with nprocs {})", 1. * MPKTS,
             1. * MPKTS * globals::nprocs, globals::nprocs);

  printlnlog("[info] mem_usage: packets occupy {:.3f} MB", MPKTS * sizeof(Packet) / 1024. / 1024.);

  if (!globals::simulation_continued_from_saved) {
    std::remove("deposition.out");
    packet_init(packets);
    zero_estimators();
  }

  // For the parallelisation of update_grid, the process needs to be told which cells belong to it.
  // The next loop is over all grid cells. For parallelisation, we want to split this loop between
  // processes. This is done by assigning each MPI process nblock cells. The residual n_leftover
  // cells are sent to processes 0 ... process n_leftover -1.
  const int nstart = grid::get_nstart(globals::my_rank);
  const int ndo = grid::get_ndo(globals::my_rank);
  const int ndo_nonempty = grid::get_ndo_nonempty(globals::my_rank);
  if (ndo > 0) {
    printlnlog(
        "process rank {} (of 0..{}) assigned {} modelgrid cells ({} nonempty): cells [{}..{}] (model "
        "has max mgi {})",
        globals::my_rank, globals::nprocs - 1, ndo, ndo_nonempty, nstart, nstart + ndo - 1, grid::get_npts_model() - 1);
  } else {
    printlnlog("process rank {} (of 0..{}) assigned {} modelgrid cells ({} nonempty)", globals::my_rank,
               globals::nprocs - 1, ndo, ndo_nonempty);
  }

  MPI_Barrier_allranks();
  globals::timestep = globals::timestep_initial;

  macroatom_open_file();
  if (ndo > 0) {
    assert_always(!estimators_file.is_open());
    estimators_file =
        fstream_required(std::format("estimators_{:04d}.out", globals::my_rank), std::ios::out | std::ios::trunc);

    if (globals::total_nlte_levels > 0 && ndo_nonempty > 0) {
      nltepop_open_file();
    }
  }

  // initialise or read in virtual packet spectra
  vpkt::init(globals::timestep, globals::simulation_continued_from_saved);

  setup_cellcache();

  while (globals::timestep < globals::timestep_finish && !terminate_early) {
    MPI_Barrier_allranks();

    // titer example: Do 3 iterations on timestep 0-6
    // globals::n_titer = (globals::timestep < 6) ? 3 : 1;
    globals::n_titer = 1;

#ifdef DO_TITER
    assert_always(globals::n_titer > 0);
#else
    assert_always(globals::n_titer == 1);
#endif
    if (globals::n_titer > 1) {
      printlnlog("Doing {} iterations on timestep {}", globals::n_titer, globals::timestep);
    }

    for (int titer = 0; titer < globals::n_titer; titer++) {
      terminate_early = do_timestep(globals::timestep, titer, packets, walltimelimitseconds);
#ifdef DO_TITER
      // No iterations over the zeroth timestep, set titer > n_titer
      if (globals::timestep == 0) titer = globals::n_titer + 1;
#endif
    }

    globals::timestep++;
  }

  // The main calculation is now over. The packets now have all stored the time, place and direction
  // at which they left the grid. Also their rest frame energies and frequencies.
  // Spectra and light curves are now extracted using exspec which is another make target of this code.

  MPI_Barrier_allranks();

  if ((globals::ntimesteps > globals::timestep_finish) || terminate_early) {
    printlnlog("RESTART_NEEDED to continue model");
  } else {
    printlnlog("No need for restart");
  }

  MPI_Barrier_allranks();
  const auto real_time_end = std::chrono::steady_clock::now();
  const auto wallclock_hours = std::chrono::duration<double>(real_time_end - real_time_start).count() / HOUR;
  printlnlog(
      "sn3d finished (job: pktprop ts {} to ts {} grid-preprop, {:.3f} wallclock hours * {} processes * {} "
      "threads = {:.3f} core hours)",
      globals::timestep_initial, globals::timestep - 1, wallclock_hours, globals::nprocs, get_max_threads(),
      wallclock_hours * globals::nprocs * get_max_threads());

  MPI_Finalize();

  std::filesystem::remove("artis.pid");

  return 0;
}
