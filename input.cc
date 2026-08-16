// Reading of the main input files: the atomic dataset (adata.txt, compositiondata.txt,
// transitiondata.txt, phixsdata_v2.txt) and the run parameters (input.txt), plus the setup of the simulation timesteps.

#include "input.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <ios>
#include <istream>
#include <iterator>
#include <limits>
#include <print>
#include <ranges>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>
#include <utility>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "kpkt.h"
#include "mpi_logging.h"
#include "packet.h"
#include "random.h"
#include "ratecoeff.h"
#include "sn3d.h"

namespace {

// level index of the ground state in the input files (0 or 1), autodetected from the first level in adata.txt
int groundstate_index_in = -1;

struct TempEnergyLevel {
  double epsilon{-1};  // Excitation energy of this level relative to the neutral ground level.
  int alltrans_startdown{};  // index into globals::alltrans for first down transition from this level
  int ndowntrans{0};  // Number of down transitions from this level
  int nuptrans{0};  // Number of up transitions to this level
  float stat_weight{0.};  // statistical weight of this level

  [[nodiscard]] constexpr auto alltrans_startup() const -> int {
    // index into globals::alltrans for first up transition from this level
    return alltrans_startdown + ndowntrans;
  }
};

struct IonTransitionsInput {
  int lower;
  int upper;
  float A;
  float coll_str;
  bool forbidden;
};

struct TempAllTransInput {
  int targetlevelindex;
  float einstein_A;
  float coll_str;
  float osc_strength;
  bool forbidden;
};

// only used temporarily during input before copy to node-shared globals::allphixstargets_levelindex and
// allphixstargets_probability
struct PhotoionTarget {
  double probability;  // fraction of phixs cross section leading to this final level
  int levelindex;  // index of upper ion level after photoionisation
};

// temporary for reading before sorting and copy to globals::linelist structure of arrays
struct TempLineTransitionInput {
  double nu;
  float einstein_A;
  int elementindex;
  int ionindex;
  int upperlevelindex;
  int lowerlevelindex;
};

// mutable views of the per-level photoionisation arrays while read_phixs_data() builds them, before
// they are published into the read-only members of globals::alllevels
struct PhixsLevelBuilders {
  std::span<int> phixsstart;
  std::span<int> nphixstargets;
  std::span<int> phixstargetstart;
};

constexpr auto inputlinecomments = std::array{
    " 0: pre_zseed: specific random number seed if > 0 or random if negative",
    " 1: ntimesteps: number of timesteps",
    " 2: timestep_start timestep_finish: timestep number range start (inclusive) and stop (not inclusive)",
    " 3: tmin_days tmax_days: start and end times [day]",
    " 4: UNUSED nusyn_min_mev nusyn_max_mev: lowest and highest frequency to synthesise [MeV]",
    " 5: UNUSED nsyn_time: number of times for synthesis",
    " 6: UNUSED start and end times for synthesis",
    " 7: UNUSED (now autodetected) model.txt number of dimensions (1, 2, or 3)",
    " 8: UNUSED compute r-light curve (1: no estimators, 2: thin cells, 3: thick cells, 4: gamma-ray heating)",
    " 9: UNUSED n_out_it: number of iterations",
    "10: UNUSED change speed of light by some factor. Change constants.h CLIGHT_PROP instead",
    "11: UNUSED gamma_kappagrey: if >0: use grey opacity for gammas, if <0: use detailed opacity",
    "12: UNUSED syn_dir: x, y, and z components of unit vector (now always 0,0,1)",
    "13: UNUSED opacity_case: opacity choice",
    "14: UNUSED rho_crit_para: free parameter for calculation of rho_crit",
    "15: UNUSED debug_packet: (>=0: activate debug output for packet id, <0: ignore)",
    "16: simulation_continued_from_saved: (0: start new simulation, 1: continue from gridsave and packets files)",
    "17: UNUSED rfcut: wavelength at which the radiation field switches from the nebular approximation to LTE.",
    "18: num_lte_timesteps",
    "19: optical_depth_is_thick num_grey_timesteps",
    "20: UNUSED max_bf_continua: (ignored; all bound-free continua are included)",
    "21: nprocs_exspec: extract spectra for n MPI tasks. sn3d will set this on start of new sim.",
    "22: UNUSED do_emission_res: this is always true for exspec, sometimes true during sn3d",
    "23: UNUSED kpktdiffusion_timescale n_kpktdiffusion_timesteps: now set in kpkt.cc",
};

// indices of the noncomment lines of input.txt that update_parameterfile() rewrites for a restart
// (the static_asserts tie each index to its description in inputlinecomments)
constexpr int inputline_timestep_range = 2;
constexpr int inputline_continue_from_saved = 16;
constexpr int inputline_nprocs_exspec = 21;
static_assert(std::string_view{inputlinecomments[inputline_timestep_range]}.starts_with(" 2:"));
static_assert(std::string_view{inputlinecomments[inputline_continue_from_saved]}.starts_with("16:"));
static_assert(std::string_view{inputlinecomments[inputline_nprocs_exspec]}.starts_with("21:"));

void read_phixs_data_table(std::istream& phixsfile, const int nphixspoints_inputtable, const int element,
                           const int lowerion, const int lowerlevel, const int upperion, int upperlevel_in,
                           std::vector<float>& tmpallphixs, std::vector<PhotoionTarget>& tmpallphixstargets,
                           size_t* mem_usage_phixs, const int phixs_file_version,
                           const PhixsLevelBuilders& levelbuilders) {
  std::string phixsline;
  const auto phixstargetstart = static_cast<int>(tmpallphixstargets.size());
  const auto lowerionlower_uniquelevelindex = get_uniquelevelindex(element, lowerion, lowerlevel);
  assert_always(levelbuilders.phixstargetstart[lowerionlower_uniquelevelindex] == -1 ||
                levelbuilders.phixstargetstart[lowerionlower_uniquelevelindex] == phixstargetstart);
  levelbuilders.phixstargetstart[lowerionlower_uniquelevelindex] = phixstargetstart;
  if (upperlevel_in >= 0) {  // file gives photoionisation to a single target state only
    int upperlevel = upperlevel_in - groundstate_index_in;
    assert_always(upperlevel >= 0);
    assert_always(levelbuilders.nphixstargets[lowerionlower_uniquelevelindex] == 0 ||
                  levelbuilders.nphixstargets[lowerionlower_uniquelevelindex] == 1);
    levelbuilders.nphixstargets[lowerionlower_uniquelevelindex] = 1;

    if (SINGLE_LEVEL_TOP_ION && (upperion == get_nions(element) - 1)) {
      // top ion has only one level, so send it to that level
      upperlevel = 0;
    }

    tmpallphixstargets.push_back({.probability = 1., .levelindex = upperlevel});
  } else {  // upperlevel < 0, indicating that a table of upper levels and their probabilities will follow
    int in_nphixstargets = 0;
    assert_always(get_noncommentline(phixsfile, phixsline));
    assert_always(std::stringstream(phixsline) >> in_nphixstargets);
    assert_always(in_nphixstargets > 0);
    // read in a table of target states and probabilities and store them
    if (!SINGLE_LEVEL_TOP_ION || upperion < get_nions(element) - 1)  // in case the top ion has nlevelsmax = 1
    {
      assert_always(levelbuilders.nphixstargets[lowerionlower_uniquelevelindex] == 0 ||
                    levelbuilders.nphixstargets[lowerionlower_uniquelevelindex] == in_nphixstargets);
      levelbuilders.nphixstargets[lowerionlower_uniquelevelindex] = in_nphixstargets;

      double probability_sum = 0.;
      for (int i = 0; i < in_nphixstargets; i++) {
        double phixstargetprobability{NAN};
        assert_always(get_noncommentline(phixsfile, phixsline));
        assert_always(std::stringstream(phixsline) >> upperlevel_in >> phixstargetprobability);
        const int upperlevel = upperlevel_in - groundstate_index_in;
        assert_always(upperlevel >= 0);
        assert_always(phixstargetprobability > 0);
        tmpallphixstargets.push_back({.probability = phixstargetprobability, .levelindex = upperlevel});

        probability_sum += phixstargetprobability;
      }
      if (fabs(probability_sum - 1.0) > 0.01) {
        printlnlog(
            "[error] photoionisation table for Z={} ionstage {} level {} has target probabilities that sum to "
            "{:g} (expected 1.0 +/- 0.01)",
            get_atomicnumber(element), get_ionstage(element, lowerion), lowerlevel, probability_sum);
        assert_always(false);
      }
    } else {  // file has table of target states and probabilities but our top ion is limited to one level
      levelbuilders.nphixstargets[lowerionlower_uniquelevelindex] = 1;

      for (int i = 0; i < in_nphixstargets; i++) {
        assert_always(get_noncommentline(phixsfile, phixsline));
      }

      // send it to the ground state of the top ion
      tmpallphixstargets.push_back({.probability = 1., .levelindex = 0});
    }
  }

  assert_always(tmpallphixs.size() % globals::NPHIXSPOINTS == 0);
  const auto tmpphixsstart = tmpallphixs.size();
  levelbuilders.phixsstart[lowerionlower_uniquelevelindex] = static_cast<int>(tmpphixsstart / globals::NPHIXSPOINTS);
  tmpallphixs.resize(tmpallphixs.size() + globals::NPHIXSPOINTS);

  const auto levelphixstable = std::span{tmpallphixs}.subspan(tmpphixsstart, globals::NPHIXSPOINTS);
  if (phixs_file_version == 1) {
    assert_always(levelbuilders.nphixstargets[lowerionlower_uniquelevelindex] == 1);
    assert_always(tmpallphixstargets[levelbuilders.phixstargetstart[lowerionlower_uniquelevelindex]].levelindex == 0);

    const double nu_edge = (epsilon(element, upperion, 0) - epsilon(element, lowerion, lowerlevel)) / H;

    auto nugrid_in = std::vector<double>(nphixspoints_inputtable);
    auto phixs_in = std::vector<double>(nphixspoints_inputtable);

    for (int i = 0; i < nphixspoints_inputtable; i++) {
      double energy = -1.;
      double phixs = -1.;
      assert_always(get_noncommentline(phixsfile, phixsline));
      assert_always(std::stringstream(phixsline) >> energy >> phixs);
      assert_always(energy >= 0);
      assert_always(phixs >= 0);
      nugrid_in[i] = nu_edge + ((energy * 13.6 * EV) / H);
      // the photoionisation cross-sections in the database are given in Mbarn=1e6 * 1e-28m^2
      // to convert to cgs units multiply by 1e-18
      phixs_in[i] = phixs * 1e-18;
    }
    const double nu_max = nugrid_in.back();

    // Now interpolate these cross-sections
    levelphixstable[0] = static_cast<float>(phixs_in[0]);

    for (int i = 1; i < globals::NPHIXSPOINTS; i++) {
      const double nu = nu_edge * (1. + (i * globals::NPHIXSNUINCREMENT));
      if (nu >= nu_max) {
        levelphixstable[i] = static_cast<float>(phixs_in[nphixspoints_inputtable - 1] * pow3(nu_max / nu));
      } else {
        assert_always(nu >= nugrid_in[0]);
        const auto index_above =
            static_cast<std::ptrdiff_t>(std::ranges::upper_bound(nugrid_in, nu) - nugrid_in.begin());
        assert_always(index_above < nphixspoints_inputtable);
        const auto index_below = index_above - 1;
        assert_always(index_below >= 0);
        const auto phixs_interp = static_cast<float>(
            std::lerp(phixs_in[index_below], phixs_in[index_above],
                      (nu - nugrid_in[index_below]) / (nugrid_in[index_above] - nugrid_in[index_below])));
        levelphixstable[i] = phixs_interp;
      }
    }
  } else {
    for (int i = 0; i < globals::NPHIXSPOINTS; i++) {
      float phixs{NAN};
      assert_always(phixsfile >> phixs);
      assert_always(phixs >= 0);

      // the photoionisation cross-sections in the database are given in Mbarn = 1e6 * 1e-28m^2
      // to convert to cgs units multiply by 1e-18
      levelphixstable[i] = static_cast<float>(phixs * 1e-18);
    }
  }

  // Every photoionisation target of this level is a level that the upper ion can recombine from, so raise the
  // upper ion's maxrecombininglevel to cover them. The topmost included ion is skipped because it has no ion
  // above it to recombine into.
  if (lowerion < get_nions(element) - 1) {
    for (int phixstargetindex = 0; phixstargetindex < levelbuilders.nphixstargets[lowerionlower_uniquelevelindex];
         phixstargetindex++) {
      const int upperlevel =
          tmpallphixstargets[levelbuilders.phixstargetstart[lowerionlower_uniquelevelindex] + phixstargetindex]
              .levelindex;
      if (upperlevel > get_maxrecombininglevel(element, lowerion + 1)) {
        globals::elements[element].ions[lowerion + 1].maxrecombininglevel = upperlevel;
      }
    }
  }

  *mem_usage_phixs += (levelbuilders.nphixstargets[lowerionlower_uniquelevelindex] * sizeof(double)) + sizeof(int);
  *mem_usage_phixs += globals::NPHIXSPOINTS * sizeof(float);
  globals::nbfcontinua += levelbuilders.nphixstargets[lowerionlower_uniquelevelindex];
  if (lowerlevel == 0 && levelbuilders.nphixstargets[lowerionlower_uniquelevelindex] > 0) {
    globals::nbfcontinua_ground++;
  }
}

void read_phixs_file(const int phixs_file_version, std::vector<float>& tmpallphixs,
                     std::vector<PhotoionTarget>& tmpallphixstargets, std::vector<int>& ion_minphixslevel_infile,
                     const PhixsLevelBuilders& levelbuilders) {
  printlnlog("reading phixs data from {}", phixsdata_filenames[phixs_file_version]);

  auto phixsfile = fstream_required(phixsdata_filenames[phixs_file_version], std::ios::in);
  std::string phixsline;
  std::istringstream ssline;
  auto mem_usage_phixs = 0ZU;

  if (phixs_file_version == 1 && phixs_file_version_exists[2]) {
    printlnlog(
        "using NPHIXSPOINTS = {} and NPHIXSNUINCREMENT = {:g} from phixsdata_v2.txt to interpolate phixsdata.txt data",
        globals::NPHIXSPOINTS, globals::NPHIXSNUINCREMENT);
    last_phixs_nuovernuedge = (1.0 + (globals::NPHIXSNUINCREMENT * (globals::NPHIXSPOINTS - 1)));
  } else if (phixs_file_version == 1) {
    globals::NPHIXSPOINTS = 100;
    globals::NPHIXSNUINCREMENT = .1;
    // not exactly where the last point is, but classic integrals go from nu_edge to 10*nu_edge
    last_phixs_nuovernuedge = 10;
    printlnlog("using NPHIXSPOINTS = {} and NPHIXSNUINCREMENT = {:g} set in input.cc", globals::NPHIXSPOINTS,
               globals::NPHIXSNUINCREMENT);
  } else {
    assert_always(phixsfile >> globals::NPHIXSPOINTS);
    assert_always(globals::NPHIXSPOINTS > 0);
    assert_always(phixsfile >> globals::NPHIXSNUINCREMENT);
    assert_always(globals::NPHIXSNUINCREMENT > 0.);
    last_phixs_nuovernuedge = (1.0 + (globals::NPHIXSNUINCREMENT * (globals::NPHIXSPOINTS - 1)));
  }

  int Z = -1;
  int upperionstage = -1;
  int upperlevel_in = -1;
  int lowerionstage = -1;
  int lowerlevel_in = -1;
  double phixs_threshold_ev = -1;  // currently just ignored, and epsilon is used instead
  while (true) {
    int nphixspoints_inputtable = 0;
    if (!get_noncommentline(phixsfile, phixsline)) {
      break;
    }
    ssline.clear();
    ssline.str(phixsline);
    if (phixs_file_version == 1) {
      assert_always(ssline >> Z >> upperionstage >> upperlevel_in >> lowerionstage >> lowerlevel_in >>
                    nphixspoints_inputtable);
    } else {
      assert_always(ssline >> Z >> upperionstage >> upperlevel_in >> lowerionstage >> lowerlevel_in >>
                    phixs_threshold_ev);
      nphixspoints_inputtable = globals::NPHIXSPOINTS;
    }
    assert_always(Z > 0);
    assert_always(upperionstage >= 2);
    assert_always(lowerionstage >= 1);

    const int element = get_elementindex(Z);

    // store only photoionisation crosssections for elements that are part of the current model atom
    bool skip_this_phixs_table = true;  // will be set to false for good data
    if (element >= 0 && get_nions(element) > 0) {
      // translate read-in ionstages to ion indices

      const int upperion = upperionstage - get_ionstage(element, 0);
      const int lowerion = lowerionstage - get_ionstage(element, 0);
      const int lowerlevel = lowerlevel_in - groundstate_index_in;
      assert_always(lowerlevel >= 0);

      // record the lowest level the file supplies a table for in each included ion, before any filtering on
      // retained levels, so that mismatched level numbering can be detected even if the misnumbered table would
      // be filtered out
      if (lowerion >= 0 && lowerion < get_nions(element)) {
        const auto lowerion_uniqueionindex = get_uniqueionindex(element, lowerion);
        ion_minphixslevel_infile[lowerion_uniqueionindex] =
            std::min(ion_minphixslevel_infile[lowerion_uniqueionindex], lowerlevel);
      }

      // store only photoionisation crosssections for ions that are part of the current model atom
      if (lowerion >= 0 && upperion < get_nions(element) && lowerlevel < get_nlevels_ionising(element, lowerion)) {
        read_phixs_data_table(phixsfile, nphixspoints_inputtable, element, lowerion, lowerlevel, upperion,
                              upperlevel_in, tmpallphixs, tmpallphixstargets, &mem_usage_phixs, phixs_file_version,
                              levelbuilders);

        skip_this_phixs_table = false;
      }
    }

    if (skip_this_phixs_table) {  // for ions or elements that are not part of the current model atom, proceed through
                                  // the table and throw away the data
      if (upperlevel_in < 0) {  // a table of target states and probabilities will follow, so read past those lines
        int nphixstargets = 0;
        assert_always(get_noncommentline(phixsfile, phixsline));
        ssline.clear();
        ssline.str(phixsline);
        assert_always(ssline >> nphixstargets);
        for (int i = 0; i < nphixstargets; i++) {
          assert_always(get_noncommentline(phixsfile, phixsline));
        }
      }
      // skip through cross section list
      for (int i = 0; i < nphixspoints_inputtable; i++) {
        if (phixs_file_version == 1) {
          assert_always(get_noncommentline(phixsfile, phixsline));
        } else {
          // one day we might want to put all of the cross section points onto a single line, so don't use getline here
          float phixs = 0;
          assert_always(phixsfile >> phixs);
        }
      }
    }
  }

  printlnlog("[info] mem_usage: photoionisation tables occupy {:.3f} MB", mem_usage_phixs / 1024. / 1024.);
}

void read_ion_levels(std::istream& adata, const int element, const int ion, const int nions, const int nlevels,
                     int nlevelsmax, const double energyoffset_ev, const double ionpot_ev,
                     std::vector<TempEnergyLevel>& temp_alllevels) {
  std::string line;
  static std::istringstream ssline;
  for (int level = 0; level < nlevels; level++) {
    int levelindex_in = 0;
    double levelenergy_ev{NAN};
    float statweight{NAN};
    int ntransitions = 0;
    assert_always(get_noncommentline(adata, line));
    ssline.clear();
    ssline.str(line);
    assert_always(ssline >> levelindex_in >> levelenergy_ev >> statweight >> ntransitions);
    if (groundstate_index_in < 0) {
      // this is the first level line of the whole file, so level == 0 and levelindex_in is the ground state index
      assert_always(levelindex_in == 0 || levelindex_in == 1);
      groundstate_index_in = levelindex_in;
      printlnlog("adata.txt: autodetected level numbering starting from {}", groundstate_index_in);
    }
    assert_always(levelindex_in == level + groundstate_index_in);

    if (level < nlevelsmax) {
      assert_always(statweight > 0.);
      const double levelepsilon = (energyoffset_ev + levelenergy_ev) * EV;
      assert_always(std::ssize(temp_alllevels) == get_ionuniquelevelindexstart(element, ion) + level);
      temp_alllevels.push_back({
          .epsilon = levelepsilon,
          .ndowntrans = 0,
          .nuptrans = 0,
          .stat_weight = statweight,
      });

      // The level contributes to the ionisinglevels if its energy
      // is below the ionisation potential and the level doesn't belong to the topmost ion included.
      // Rate coefficients are only available for ionising levels.
      if (levelenergy_ev < ionpot_ev && ion < nions - 1) {
        globals::elements[element].ions[ion].nlevels_ionising++;
      }
    }
  }
}

void read_ion_transitions(std::istream& ftransitiondata, const int ion_transition_count_in_file,
                          std::vector<IonTransitionsInput>& iontransitiontable, const int nlevels_requiretransitions,
                          const int nlevelskept, const int atomicnumber, const int ionstage,
                          const int nlevels_in_file) {
  iontransitiontable.clear();
  iontransitiontable.reserve(ion_transition_count_in_file);
  bool found_groundstate_lower = false;

  std::string line;

  // will be autodetected from first table row. old format had an index column and no collstr or forbidden columns
  bool oldtransitionformat = false;

  for (int i = 0; i < ion_transition_count_in_file; i++) {
    int lower_in = -1;
    int upper_in = -1;
    float A = 0;
    float coll_str = -1.;
    int intforbidden = 0;
    assert_always(getline(ftransitiondata, line));
    if (i == 0) {
      int column_count = 0;
      auto remainder = std::string_view{line};
      double dummy{};
      while (parse_next_token(remainder, dummy)) {
        column_count++;
      }
      assert_always(column_count == 4 || column_count == 5);
      oldtransitionformat = (column_count == 4);
    }
    auto remainder = std::string_view{line};
    if (!oldtransitionformat) {
      assert_always(parse_next_token(remainder, lower_in) && parse_next_token(remainder, upper_in) &&
                    parse_next_token(remainder, A) && parse_next_token(remainder, coll_str) &&
                    parse_next_token(remainder, intforbidden));
    } else {
      int transindex = 0;  // not used
      assert_always(parse_next_token(remainder, transindex) && parse_next_token(remainder, lower_in) &&
                    parse_next_token(remainder, upper_in) && parse_next_token(remainder, A));
    }
    const int lower = lower_in - groundstate_index_in;
    const int upper = upper_in - groundstate_index_in;
    assert_always(lower >= 0);
    assert_always(lower < upper);
    found_groundstate_lower = found_groundstate_lower || (lower == 0);
    if (upper >= nlevels_in_file) {
      printlnlog(
          "[error] transitiondata.txt: Z={} ionstage {}: transition {} -> {} refers to a level beyond the {} levels of "
          "this ion in adata.txt. The level numbering may not match the ground state index {} detected from adata.txt",
          atomicnumber, ionstage, lower_in, upper_in, nlevels_in_file, groundstate_index_in);
      assert_always(false);
    }
    if (lower >= nlevelskept || upper >= nlevelskept) {
      continue;
    }
    iontransitiontable.push_back(
        {.lower = lower, .upper = upper, .A = A, .coll_str = coll_str, .forbidden = (intforbidden == 1)});
  }

  // if the ion has any transitions then the ground state is normally the lower level of at least one of them. A
  // sparse table could legitimately omit the ground state (missing transitions are filled in below), so warn
  // without aborting
  if (ion_transition_count_in_file > 0 && !found_groundstate_lower) {
    printlnlog(
        "[warning] transitiondata.txt: Z={} ionstage {}: no transition has the ground state (level index {} in the "
        "input files) as its lower level. Check that the level numbering matches adata.txt",
        atomicnumber, ionstage, groundstate_index_in);
  }

  const auto proj_lowerupper = [](const IonTransitionsInput& t) { return std::tie(t.lower, t.upper); };
  std::ranges::sort(iontransitiontable, std::less<>{}, proj_lowerupper);

  assert_always(nlevels_requiretransitions <= nlevelskept);

  const auto old_transitioncount = std::ssize(iontransitiontable);
  for (int lower = 0; lower < nlevels_requiretransitions; lower++) {
    for (int upper = lower + 1; upper < nlevelskept; upper++) {
      // binary search the sorted table read from the file (excluding any pairs appended by this loop)
      const bool transition_exists =
          std::ranges::binary_search(iontransitiontable.begin(), iontransitiontable.begin() + old_transitioncount,
                                     std::tuple{lower, upper}, std::ranges::less{}, proj_lowerupper);
      if (!transition_exists) {
        iontransitiontable.push_back({.lower = lower, .upper = upper, .A = 0., .coll_str = -2., .forbidden = true});
      }
    }
  }
  const auto added_transitions = std::ssize(iontransitiontable) - old_transitioncount;
  if (added_transitions > 0) {
    printlnlog("[info] added {} missing transitions with A=0 to iontransitiontable", added_transitions);
    std::ranges::sort(iontransitiontable, std::less<>{}, proj_lowerupper);
  }
}

void add_transitions_to_unsorted_linelist(const int element, const int ion,
                                          const std::span<const IonTransitionsInput> iontransitiontable,
                                          std::vector<TempLineTransitionInput>& temp_linelist,
                                          std::vector<TempAllTransInput>& temp_alltranslist,
                                          std::span<TempEnergyLevel> temp_alllevels) {
  const auto nlevels = get_nlevels(element, ion);
  const auto ion_levels = std::span{temp_alllevels}.subspan(get_ionuniquelevelindexstart(element, ion), nlevels);

  const int nlines_initial = globals::nlines;
  ptrdiff_t ion_updowntranscount = 0;
  // pass 0 to get transition counts of each level
  // pass 1 to allocate and fill transition arrays
  for (int pass = 0; pass < 2; pass++) {
    globals::nlines = nlines_initial;
    if (pass == 1) {
      ptrdiff_t alltransindex = std::ssize(temp_alltranslist);
      temp_alltranslist.resize(std::ssize(temp_alltranslist) + ion_updowntranscount);
      assert_always(std::ssize(temp_alltranslist) >= std::ssize(temp_linelist));
      for (int level = 0; level < nlevels; level++) {
        ion_levels[level].alltrans_startdown = static_cast<int>(alltransindex);
        alltransindex += ion_levels[level].ndowntrans;
        alltransindex += ion_levels[level].nuptrans;

        ion_levels[level].ndowntrans = 0;
        ion_levels[level].nuptrans = 0;
      }
      assert_always(alltransindex < std::numeric_limits<int>::max());
    }

    // iontransitiontable is sorted by (lower, upper), so any duplicate transitions between the same
    // pair of levels are adjacent. Track the previous accepted transition to detect them.
    int prev_lower = -1;
    int prev_upper = -1;
    int prev_lineindex = -1;

    ion_updowntranscount = 0;
    for (const auto& transition : iontransitiontable) {
      const int level = transition.upper;
      const int lowerlevel = transition.lower;

      if ((lowerlevel >= nlevels) || (level >= nlevels)) {
        continue;
      }
      const double nu_trans = (ion_levels[level].epsilon - ion_levels[lowerlevel].epsilon) / H;
      if (!(nu_trans > 0)) {
        continue;
      }

      // Make sure that we don't allow duplicate. In that case take only the lines first occurrence
      const bool is_duplicate = (lowerlevel == prev_lower && level == prev_upper);

      if (!is_duplicate) {
        prev_lower = lowerlevel;
        prev_upper = level;
        prev_lineindex = globals::nlines++;

        const int nupperdowntrans = ion_levels[level].ndowntrans + 1;
        ion_levels[level].ndowntrans = nupperdowntrans;

        const int nloweruptrans = ion_levels[lowerlevel].nuptrans + 1;
        ion_levels[lowerlevel].nuptrans = nloweruptrans;

        ion_updowntranscount += 2;

        if (pass == 1) {
          // absorption oscillator strength f_lu from A_ul via f_lu = (g_u/g_l) * m_e c^3 / (8 pi^2 e^2 nu^2) * A_ul
          const auto g_ratio = static_cast<double>(ion_levels[level].stat_weight) / ion_levels[lowerlevel].stat_weight;
          const auto f_lu =
              static_cast<float>(g_ratio * ME * pow3(CLIGHT) / (8 * pow2(QE * nu_trans * PI)) * transition.A);
          assert_always(std::isfinite(f_lu));

          temp_linelist.push_back({
              .nu = nu_trans,
              .einstein_A = transition.A,
              .elementindex = element,
              .ionindex = ion,
              .upperlevelindex = level,
              .lowerlevelindex = lowerlevel,
          });

          // (each transition's index into the frequency-sorted linelist is set later by create_shared_alltranslist)

          temp_alltranslist[ion_levels[level].alltrans_startdown + nupperdowntrans - 1] = {
              .targetlevelindex = lowerlevel,
              .einstein_A = transition.A,
              .coll_str = transition.coll_str,
              .osc_strength = f_lu,
              .forbidden = transition.forbidden,
          };
          const auto lowerstartup = ion_levels[lowerlevel].alltrans_startup();
          temp_alltranslist[lowerstartup + nloweruptrans - 1] = {
              .targetlevelindex = level,
              .einstein_A = transition.A,
              .coll_str = transition.coll_str,
              .osc_strength = f_lu,
              .forbidden = transition.forbidden,
          };
        }

      } else if (pass == 1) {
        // A second transition between the same pair of levels (e.g. from a different physical process):
        // combine it into the existing line. The duplicate immediately follows the first occurrence
        // (the table is sorted), so the transition entries filled most recently are the ones to update.

        if ((temp_linelist[prev_lineindex].elementindex != element) ||
            (temp_linelist[prev_lineindex].ionindex != ion) ||
            (temp_linelist[prev_lineindex].upperlevelindex != level) ||
            (temp_linelist[prev_lineindex].lowerlevelindex != lowerlevel)) {
          printlnlog("[error] Failure to identify level pair for duplicate bb-transition ... going to abort now");
          printlnlog("   element {} ion {} targetlevel {} level {}", element, ion, lowerlevel, level);
          printlnlog("   duplicate of lineindex {}", prev_lineindex);
          printlnlog("   A_ul {:g}, coll_str {:g}", transition.A, transition.coll_str);
          printlnlog(
              "   globals::linelist[lineindex].elementindex {}, globals::linelist[lineindex].ionindex {}, "
              "globals::linelist[lineindex].upperlevelindex {}, globals::linelist[lineindex].lowerlevelindex {}",
              temp_linelist[prev_lineindex].elementindex, temp_linelist[prev_lineindex].ionindex,
              temp_linelist[prev_lineindex].upperlevelindex, temp_linelist[prev_lineindex].lowerlevelindex);
          std::abort();
        }

        const auto g_ratio = static_cast<double>(ion_levels[level].stat_weight) / ion_levels[lowerlevel].stat_weight;
        const auto f_lu =
            static_cast<float>(g_ratio * ME * pow3(CLIGHT) / (8 * pow2(QE * nu_trans * PI)) * transition.A);

        auto& downtransition =
            temp_alltranslist[ion_levels[level].alltrans_startdown + ion_levels[level].ndowntrans - 1];

        assert_always(downtransition.targetlevelindex == lowerlevel);

        downtransition.einstein_A += transition.A;
        downtransition.osc_strength += f_lu;
        downtransition.coll_str = std::max(downtransition.coll_str, transition.coll_str);

        const auto lowerstartup = ion_levels[lowerlevel].alltrans_startup();
        auto& uptransition = temp_alltranslist[lowerstartup + ion_levels[lowerlevel].nuptrans - 1];

        assert_always(uptransition.targetlevelindex == level);

        uptransition.einstein_A += transition.A;
        uptransition.osc_strength += f_lu;
        uptransition.coll_str = std::max(uptransition.coll_str, transition.coll_str);
      }
    }
  }
}

auto calculate_nlevels_groundterm(const int element, const int ion) -> int {
  const int nlevels = get_nlevels(element, ion);
  if (nlevels <= 2) {
    return nlevels;
  }

  int nlevels_groundterm = 1;
  // detect single-level ground term
  const double endiff10 = epsilon(element, ion, 1) - epsilon(element, ion, 0);
  const double endiff21 = epsilon(element, ion, 2) - epsilon(element, ion, 1);
  if (endiff10 > 2. * endiff21) {
    nlevels_groundterm = 1;
  } else {
    for (int level = 1; level < nlevels - 2; level++) {
      const double endiff1 = epsilon(element, ion, level) - epsilon(element, ion, level - 1);
      const double endiff2 = epsilon(element, ion, level + 1) - epsilon(element, ion, level);
      if (endiff2 > 2. * endiff1) {
        nlevels_groundterm = level + 1;
        break;
      }
    }
  }

  // there should be no duplicate stat weights within the ground term
  // limit the ground multiplet to nnnnlowest levels below the first duplicated stat weight
  for (int level_a = 1; level_a < nlevels_groundterm; level_a++) {
    const auto g_a = stat_weight(element, ion, level_a);

    for (int level_b = 0; level_b < level_a; level_b++) {
      const auto g_b = stat_weight(element, ion, level_b);
      if (fabs(g_a - g_b) < 0.4) {
        // level_a is outside the ground term because of duplicate stat weight
        // highest ground level index is level_a - 1, so nlevels_groundterm == level_a
        return level_a;
      }
    }
  }

  return nlevels_groundterm;
}

// Return the closest ground level continuum index to the given edge
// frequency. If the given edge frequency is redder than the reddest
// continuum return -1.
// groundcont_nu_edge is in ascending order (red to blue)
auto search_groundphixslist(const double nu_edge, const int element_in, const int ion_in, const int level_in) -> int {
  assert_always((USE_LUT_PHOTOION || USE_ION_BFHEATING_ESTIMATORS));

  // an atomic dataset where no ground level has a photoionisation table leaves this list empty,
  // which the rest of the code allows for, so do not index into it before checking
  if (globals::nbfcontinua_ground <= 0) {
    return -1;
  }

  if (nu_edge < globals::groundcont_nu_edge[0]) {
    return -1;
  }

  int i = 1;
  for (i = 1; i < globals::nbfcontinua_ground; i++) {
    if (nu_edge < globals::groundcont_nu_edge[i]) {
      break;
    }
  }

  if (i == globals::nbfcontinua_ground) {
    const int element = globals::groundcont_element[i - 1];
    const int ion = globals::groundcont_ion[i - 1];
    if (element == element_in && ion == ion_in && level_in == 0) {
      return i - 1;
    }

    printlnlog(
        "[error] search_groundphixslist: element {}, ion {}, level {} has edge_frequency {:g} equal to the bluest "
        "ground-level continuum",
        element_in, ion_in, level_in, nu_edge);
    printlnlog("  search_groundphixslist: bluest ground level continuum is element {}, ion {} at nu_edge {:g}", element,
               ion, globals::groundcont_nu_edge[i - 1]);
    printlnlog("  search_groundphixslist: i {}, nbfcontinua_ground {}", i, globals::nbfcontinua_ground);
    printlnlog("  This shouldn't happen, but is possible if there are multiple levels in the adata file at energy=0");
    const int nlevels_in = get_nlevels(element_in, ion_in);
    constexpr int maxshownlevels = 10;
    for (int looplevels = 0; looplevels < std::min(nlevels_in, maxshownlevels); looplevels++) {
      printlnlog("  element {}, ion {}, level {}, energy {:g}", element_in, ion_in, looplevels,
                 epsilon(element_in, ion_in, looplevels));
    }
    if (nlevels_in > maxshownlevels) {
      printlnlog("  (energies of the remaining {} levels omitted)", nlevels_in - maxshownlevels);
    }
    assert_always(false);
    return i - 1;
  }

  const double left_diff = nu_edge - globals::groundcont_nu_edge[i - 1];
  const double right_diff = globals::groundcont_nu_edge[i] - nu_edge;
  return (left_diff <= right_diff) ? i - 1 : i;
}

// set up the photoionisation transition lists
// and temporary gamma/kappa lists for each thread
void setup_phixs_list() {
  printlnlog("[info] setup_phixs_list: number of bfcontinua {}", globals::nbfcontinua);
  printlnlog("[info] setup_phixs_list: number of ground-level bfcontinua {}", globals::nbfcontinua_ground);

  struct TempPhotoionTransitionInput {
    double nu_edge;
    int element;
    int ion;
    int level;
    int phixstargetindex;
    int upperlevel;
    int uniquelevelindex;
    double probability;
    int groundcontestimindex;
  };

  auto groundcont_nu_edge = MPI_shared_array<double>(globals::nbfcontinua_ground);
  auto groundcont_element = MPI_shared_array<int>(globals::nbfcontinua_ground);
  auto groundcont_ion = MPI_shared_array<int>(globals::nbfcontinua_ground);

  // filled in by the node leaders below, then published as a read-only globals::alllevels member
  auto alllevels_closestgroundlevelcont = MPI_shared_array<int>(std::ssize(globals::alllevels.epsilon), -1);

  if (globals::rank_in_node == 0) {
    int nextgroundcontindex = 0;
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++) {
        const int level = 0;
        const int nphixstargets = get_nphixstargets(element, ion, level);
        if (nphixstargets == 0) {
          continue;
        }
        const double E_threshold = get_phixs_threshold(element, ion, level, 0);
        const double nu_edge = E_threshold / H;

        assert_testmodeonly(nextgroundcontindex < globals::nbfcontinua_ground);
        groundcont_nu_edge[nextgroundcontindex] = nu_edge;
        groundcont_element[nextgroundcontindex] = element;
        groundcont_ion[nextgroundcontindex] = ion;
        nextgroundcontindex++;
      }
    }
    assert_always(nextgroundcontindex == globals::nbfcontinua_ground);
    std::ranges::sort(std::views::zip(groundcont_nu_edge, groundcont_element, groundcont_ion),
                      [](const auto& lhs, const auto& rhs) { return std::get<0>(lhs) < std::get<0>(rhs); });
  }
  MPI_Barrier_node();
  globals::groundcont_nu_edge = std::move(groundcont_nu_edge);
  globals::groundcont_element = std::move(groundcont_element);
  globals::groundcont_ion = std::move(groundcont_ion);

  auto allcont = MPI_shared_array<TempPhotoionTransitionInput>(globals::nbfcontinua);
  printlnlog("[info] mem_usage: photoionisation list occupies {:.3f} MB",
             globals::nbfcontinua * (sizeof(TempPhotoionTransitionInput)) / 1024. / 1024.);
  const auto groundcontindices = std::ranges::iota_view{0, globals::nbfcontinua_ground};
  // the continuum list, ion ground continuum indices, and closest ground level continua are all in node-shared
  // memory, so only the node leaders fill them (synchronised by the barriers below)
  if (globals::rank_in_node == 0) {
    int allcontindex = 0;
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++) {
        int groundcontindex =
            static_cast<int>(std::ranges::find_if(groundcontindices,
                                                  [=](const auto& i) {
                                                    return (globals::groundcont_element[i] == element) &&
                                                           (globals::groundcont_ion[i] == ion);
                                                  }) -
                             groundcontindices.begin());
        if (groundcontindex >= globals::nbfcontinua_ground) {
          groundcontindex = -1;
        }
        globals::elements[element].ions[ion].groundcontindex = groundcontindex;
        const int nlevels = get_nlevels_ionising(element, ion);
        for (int level = 0; level < nlevels; level++) {
          const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
          const int nphixstargets = get_nphixstargets(uniquelevelindex);

          // nlevels_ionising comes from the level energies alone, so a level below the threshold can still
          // have no photoionisation table, and get_phixs_threshold() must not be called for those (their
          // phixstargetstart is still the -1 sentinel)
          if (nphixstargets == 0) {
            continue;
          }

          if constexpr (USE_LUT_PHOTOION || USE_ION_BFHEATING_ESTIMATORS) {
            // depends only on the level, so it is found once here rather than once per target below
            const double nu_edge_target0 = get_phixs_threshold(element, ion, level, 0) / H;
            alllevels_closestgroundlevelcont[uniquelevelindex] =
                search_groundphixslist(nu_edge_target0, element, ion, level);
          }

          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            assert_always(allcontindex < std::ssize(allcont));

            // See AllCont::groundcontestimindex. closestgroundlevelcont is left at its -1 initialiser when
            // neither estimator is enabled, so no further guard on those options is needed here.
            // TODO: the estimators are written at this nearest-edge slot but normalised in update_grid.cc
            // at get_groundcontindex(element, ion); the two disagree if two ions share a ground threshold.
            const int groundcontestimindex =
                (level == 0 && phixstargetindex == 0) ? alllevels_closestgroundlevelcont[uniquelevelindex] : -1;

            allcont[allcontindex] = {
                .nu_edge = get_phixs_threshold(element, ion, level, phixstargetindex) / H,
                .element = element,
                .ion = ion,
                .level = level,
                .phixstargetindex = phixstargetindex,
                .upperlevel = get_phixsupperlevel(uniquelevelindex, phixstargetindex),
                .uniquelevelindex = uniquelevelindex,
                .probability = get_phixsprobability(uniquelevelindex, phixstargetindex),
                .groundcontestimindex = groundcontestimindex,
            };

            allcontindex++;
          }
        }
      }
    }

    assert_always(allcontindex == globals::nbfcontinua);
  }
  MPI_Barrier_node();
  globals::alllevels.closestgroundlevelcont = std::move(alllevels_closestgroundlevelcont);

  assert_always(globals::nbfcontinua >= 0);  // was initialised as -1 before startup
  // just so that clang-tidy doesn't throw errors on the assumption that nbfcontinua is changing
  const auto nbfcontinua = globals::nbfcontinua;

  if (nbfcontinua > 0) {
    // indices above were temporary only. continuum index should be to the sorted list
    MPI_Barrier_node();
    if (globals::rank_in_node == 0) {
      std::ranges::SORT_OR_STABLE_SORT(allcont, std::ranges::less{}, &TempPhotoionTransitionInput::nu_edge);
    }
    MPI_Barrier_node();

    auto allcont_nu_edge = MPI_shared_array<double>(nbfcontinua);
    auto allcont_element = MPI_shared_array<int>(nbfcontinua);
    auto allcont_ion = MPI_shared_array<int>(nbfcontinua);
    auto allcont_level = MPI_shared_array<int>(nbfcontinua);
    auto allcont_phixstargetindex = MPI_shared_array<int>(nbfcontinua);
    auto allcont_upperlevel = MPI_shared_array<int>(nbfcontinua);
    auto allcont_uniquelevelindex = MPI_shared_array<int>(nbfcontinua);
    auto allcont_probability = MPI_shared_array<double>(nbfcontinua);
    auto allcont_groundcontestimindex = MPI_shared_array<int>(nbfcontinua);
    if (globals::rank_in_node == 0) {
      std::ranges::copy(std::views::transform(allcont, &TempPhotoionTransitionInput::nu_edge), allcont_nu_edge.begin());
      std::ranges::copy(std::views::transform(allcont, &TempPhotoionTransitionInput::element), allcont_element.begin());
      std::ranges::copy(std::views::transform(allcont, &TempPhotoionTransitionInput::ion), allcont_ion.begin());
      std::ranges::copy(std::views::transform(allcont, &TempPhotoionTransitionInput::level), allcont_level.begin());
      std::ranges::copy(std::views::transform(allcont, &TempPhotoionTransitionInput::phixstargetindex),
                        allcont_phixstargetindex.begin());
      std::ranges::copy(std::views::transform(allcont, &TempPhotoionTransitionInput::upperlevel),
                        allcont_upperlevel.begin());
      std::ranges::copy(std::views::transform(allcont, &TempPhotoionTransitionInput::uniquelevelindex),
                        allcont_uniquelevelindex.begin());
      std::ranges::copy(std::views::transform(allcont, &TempPhotoionTransitionInput::probability),
                        allcont_probability.begin());
      std::ranges::copy(std::views::transform(allcont, &TempPhotoionTransitionInput::groundcontestimindex),
                        allcont_groundcontestimindex.begin());
    }
    MPI_Barrier_node();
    globals::allcont.nu_edge = std::move(allcont_nu_edge);
    globals::allcont.element = std::move(allcont_element);
    globals::allcont.ion = std::move(allcont_ion);
    globals::allcont.level = std::move(allcont_level);
    globals::allcont.phixstargetindex = std::move(allcont_phixstargetindex);
    globals::allcont.upperlevel = std::move(allcont_upperlevel);
    globals::allcont.uniquelevelindex = std::move(allcont_uniquelevelindex);
    globals::allcont.probability = std::move(allcont_probability);
    globals::allcont.groundcontestimindex = std::move(allcont_groundcontestimindex);

    auto allcont_bfestimindex = MPI_shared_array<int>(nbfcontinua);
    std::vector<double> temp_bfestim_nu_edge;
    if (globals::rank_in_node == 0) {
      for (int i = 0; i < nbfcontinua; i++) {
        if (DETAILED_BF_ESTIMATORS_ON &&
            LEVEL_HAS_BFEST(get_atomicnumber(globals::allcont.element[i]),
                            get_ionstage(globals::allcont.element[i], globals::allcont.ion[i]),
                            globals::allcont.level[i])) {
          allcont_bfestimindex[i] = static_cast<int>(temp_bfestim_nu_edge.size());
          temp_bfestim_nu_edge.push_back(globals::allcont.nu_edge[i]);
        } else {
          allcont_bfestimindex[i] = -1;
        }
      }
    }
    globals::allcont.bfestimindex = std::move(allcont_bfestimindex);
    auto bfestimcount = std::ssize(temp_bfestim_nu_edge);
    MPI_Bcast_safe(bfestimcount, 0, globals::mpi_comm_node);
    auto bfestim_nu_edge = MPI_shared_array<double>(bfestimcount);
    if (globals::rank_in_node == 0) {
      std::ranges::copy(temp_bfestim_nu_edge, bfestim_nu_edge.begin());
    }
    MPI_Barrier_node();
    globals::bfestim_nu_edge = std::move(bfestim_nu_edge);

    setup_photoion_luts();
  }
  printlnlog("[info] bound-free estimators track bfestimcount {} photoionisation transitions",
             globals::bfestim_nu_edge.size());
}

void read_autoion_data() {
  // filled in below and published as read-only globals::alllevels members (with no autoion.txt, the
  // initial values are published as-is)
  const auto uniquelevelcount = std::ssize(globals::alllevels.epsilon);
  auto alllevels_allautoion_start = MPI_shared_array<int>(uniquelevelcount, -1);
  auto alllevels_nautoiondowntrans = MPI_shared_array<int>(uniquelevelcount, 0);
  auto alllevels_nautoionuptrans = MPI_shared_array<int>(uniquelevelcount, 0);

  // read in autoionisation rate data
  const bool have_autoion_file = std::filesystem::exists("autoion.txt");

  std::vector<globals::LevelAutoion> temp_allautoion;

  // only the node leaders parse the file, counting transitions into rank-local arrays and publishing them to the
  // node-shared arrays below
  std::vector<int> temp_nautoiondowntrans;
  std::vector<int> temp_nautoionuptrans;
  std::vector<int> temp_allautoion_start;
  ptrdiff_t nautoion_stored = 0;  // for the collective node-shared allocation below

  if (have_autoion_file && globals::rank_in_node == 0) {
    reserve_resize(temp_nautoiondowntrans, uniquelevelcount);
    std::ranges::fill(temp_nautoiondowntrans, 0);
    reserve_resize(temp_nautoionuptrans, uniquelevelcount);
    std::ranges::fill(temp_nautoionuptrans, 0);
    reserve_resize(temp_allautoion_start, uniquelevelcount);
    std::ranges::fill(temp_allautoion_start, -1);

    printlnlog("Reading autoion.txt for autoionisation data.");
    auto autoionfile = fstream_required("autoion.txt", std::ios::in);
    std::string autoionline;
    int Z = -1;
    int upperionstage = -1;
    int upperlevel_in = -1;
    int lowerionstage = -1;
    int lowerlevel_in = -1;
    double autoion_A = -1;
    bool allautoion_levels_are_not_nlte = true;
    bool allautoion_levels_are_nlte = true;
    int autoion_transitions_read = 0;

    // highest autoionising level index of each included ion as numbered in the file, before conversion
    auto ion_max_autoionlevel_in = std::vector<int>(get_includedions(), -1);

    while (get_noncommentline(autoionfile, autoionline)) {
      assert_always(std::istringstream(autoionline) >> Z >> upperionstage >> upperlevel_in >> lowerionstage >>
                    lowerlevel_in >> autoion_A);
      autoion_transitions_read++;

      assert_always(Z > 0);
      assert_always(upperionstage >= 2);
      assert_always(lowerionstage >= 1);

      const int element = get_elementindex(Z);

      if (element >= 0 && get_nions(element) > 0) {
        // translate read-in ionstages to ion indices

        const int upperion = upperionstage - get_ionstage(element, 0);
        const int lowerion = lowerionstage - get_ionstage(element, 0);
        // store only for ions that are part of the current model atom
        if (lowerion >= 0 && upperion < get_nions(element)) {
          const auto lower_uniqueionindex = get_uniqueionindex(element, lowerion);
          ion_max_autoionlevel_in[lower_uniqueionindex] =
              std::max(ion_max_autoionlevel_in[lower_uniqueionindex], lowerlevel_in);

          const int lowerlevel = lowerlevel_in - groundstate_index_in;
          const int upperlevel = upperlevel_in - groundstate_index_in;

          assert_always(upperion >= 0 && upperion < get_nions(element));
          assert_always(lowerion >= 0 && lowerion < get_nions(element));
          assert_always(lowerlevel >= 0 && lowerlevel < get_nlevels(element, lowerion));
          assert_always(upperlevel >= 0 && upperlevel < get_nlevels(element, upperion));
          assert_always(upperion > lowerion);
          const bool level_is_nlte = is_nlte(element, lowerion, lowerlevel);
          allautoion_levels_are_not_nlte = allautoion_levels_are_not_nlte && !level_is_nlte;
          allautoion_levels_are_nlte = allautoion_levels_are_nlte && level_is_nlte;

          const auto lower_uniquelevelindex = get_uniquelevelindex(element, lowerion, lowerlevel);
          const auto upper_uniquelevelindex = get_uniquelevelindex(element, upperion, upperlevel);

          temp_nautoiondowntrans[lower_uniquelevelindex] += 1;
          temp_nautoionuptrans[upper_uniquelevelindex] += 1;

          if (temp_allautoion_start[lower_uniquelevelindex] < 0) {
            assert_always(temp_nautoiondowntrans[lower_uniquelevelindex] == 1);
            //  this is the first autoionizing transition for this level, so set the start index
            temp_allautoion_start[lower_uniquelevelindex] = static_cast<int>(temp_allautoion.size());
          }

          temp_allautoion.push_back({
              .autoion_A = static_cast<float>(autoion_A),
              .elementindex = element,
              .lowerionindex = lowerion,
              .lowerlevelindex = lowerlevel,
              .upperionindex = upperion,
              .upperlevelindex = upperlevel,
          });
        }
      }
    }

    assert_always(allautoion_levels_are_not_nlte || allautoion_levels_are_nlte);

    // the autoionising levels of an ion must be a contiguous block ending at its highest level (asserted below), so
    // for each ion with autoionisation data the highest autoionising level in the file must be the ion's top level.
    // Checking the raw file indices detects an autoion.txt whose level numbering disagrees with adata.txt
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        const auto max_autoionlevel_in = ion_max_autoionlevel_in[get_uniqueionindex(element, ion)];
        const auto toplevel_in = get_nlevels(element, ion) - 1 + groundstate_index_in;
        if (max_autoionlevel_in >= 0 && max_autoionlevel_in != toplevel_in) {
          printlnlog(
              "[error] autoion.txt: Z={} ionstage {}: the highest autoionising level index {} is not the ion's highest "
              "level index {}. The autoion.txt level numbering may not match the ground state index {} detected from "
              "adata.txt",
              get_atomicnumber(element), get_ionstage(element, ion), max_autoionlevel_in, toplevel_in,
              groundstate_index_in);
          assert_always(false);
        }
      }
    }

    printlnlog("autoion.txt: read {} autoionisation transitions, stored {} for the included ions",
               autoion_transitions_read, temp_allautoion.size());

    nautoion_stored = std::ssize(temp_allautoion);
  }
  MPI_Bcast_safe(nautoion_stored, 0, globals::mpi_comm_node);

  globals::allautoion = MPI_shared_array<globals::LevelAutoion>(nautoion_stored);
  if (globals::rank_in_node == 0) {
    std::ranges::copy(temp_allautoion, globals::allautoion.begin());
    std::ranges::copy(temp_nautoiondowntrans, alllevels_nautoiondowntrans.begin());
    std::ranges::copy(temp_nautoionuptrans, alllevels_nautoionuptrans.begin());
    std::ranges::copy(temp_allautoion_start, alllevels_allautoion_start.begin());
  }
  MPI_Barrier_node();
  globals::alllevels.allautoion_start = std::move(alllevels_allautoion_start);
  globals::alllevels.nautoiondowntrans = std::move(alllevels_nautoiondowntrans);
  globals::alllevels.nautoionuptrans = std::move(alllevels_nautoionuptrans);

  // Plan is that autoionizing levels will be explicitly included in the NLTE population solver, but that their level
  // populations do not need to be accurately known - so if the ion has a superlevel already, then we will try to attach
  // the autoionizing level populations to that for all purposes outside the NLTE solver. For this, the ions need to
  // know how many autoionizing levels they have. So count those up now (only the node leaders, since the counts are
  // written to the node-shared ion data).

  if (have_autoion_file && globals::rank_in_node == 0) {
    int nlevels_autoion_sum = 0;
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        const int nlevels = get_nlevels(element, ion);
        int nlevels_autoion = 0;
        bool found_autoion_level = false;
        for (int level = 0; level < nlevels; level++) {
          const auto level_is_autoionizing = get_nautoiondowntrans(element, ion, level) > 0;

          if (level_is_autoionizing) {
            nlevels_autoion++;
            found_autoion_level = true;
          }

          // once we have found one autoionizing level, all higher levels should also be autoionizing
          if (found_autoion_level) {
            assert_always(level_is_autoionizing);
          }
        }
        globals::elements[element].ions[ion].nlevels_autoion = nlevels_autoion;
        nlevels_autoion_sum += nlevels_autoion;
      }
    }
    assert_always(nlevels_autoion_sum == nautoion_stored);
  }
  MPI_Barrier_node();
}

void read_phixs_data() {
  globals::nbfcontinua_ground = 0;
  globals::nbfcontinua = 0;
  std::vector<float> tmpallphixs;
  std::vector<PhotoionTarget> tmpallphixstargets;

  // filled in by the node-leader-only parsing below, then published as read-only globals::alllevels members
  const auto uniquelevelcount = std::ssize(globals::alllevels.epsilon);
  auto alllevels_phixsstart = MPI_shared_array<int>(uniquelevelcount, -1);
  auto alllevels_nphixstargets = MPI_shared_array<int>(uniquelevelcount, 0);
  auto alllevels_phixstargetstart = MPI_shared_array<int>(uniquelevelcount, -1);
  auto alllevels_bflist_start = MPI_shared_array<int>(uniquelevelcount, -1);
  const auto levelbuilders = PhixsLevelBuilders{.phixsstart = alllevels_phixsstart.span(),
                                                .nphixstargets = alllevels_nphixstargets.span(),
                                                .phixstargetstart = alllevels_phixstargetstart.span()};

  // read in photoionisation cross sections
  phixs_file_version_exists[0] = false;
  phixs_file_version_exists[1] = std::filesystem::exists(phixsdata_filenames[1]);
  phixs_file_version_exists[2] = std::filesystem::exists(phixsdata_filenames[2]);

  // just in case the file system was faulty and the ranks disagree on the existence of the files
  MPI_Allreduce_safe(phixs_file_version_exists, MPI_LOR, MPI_COMM_WORLD);
  assert_always(phixs_file_version_exists[1] || phixs_file_version_exists[2]);  // at least one must exist
  if (phixs_file_version_exists[1] && phixs_file_version_exists[2]) {
    printlnlog(
        "Reading two phixs files: Reading phixsdata_v2.txt first so we use NPHIXSPOINTS and NPHIXSNUINCREMENT from "
        "phixsdata_v2.txt to interpolate the phixsdata.txt data");
  }

  // per phixs file version: lowest level of each included ion that the file supplies a cross-section table for,
  // recorded before any filtering on retained levels (INT_MAX where the file has no tables for the ion)
  auto ion_minphixslevel_infile = std::array<std::vector<int>, 3>{};
  for (const int v : {1, 2}) {
    ion_minphixslevel_infile[v].assign(get_includedions(), std::numeric_limits<int>::max());
  }

  // only the node leaders parse the phixs files: the level lists they populate are in node-shared memory, and the
  // rank-local totals parsed from the files are broadcast within the node below
  if (globals::rank_in_node == 0) {
    if (phixs_file_version_exists[2]) {
      read_phixs_file(2, tmpallphixs, tmpallphixstargets, ion_minphixslevel_infile[2], levelbuilders);
    }

    if (phixs_file_version_exists[1]) {
      read_phixs_file(1, tmpallphixs, tmpallphixstargets, ion_minphixslevel_infile[1], levelbuilders);
    }

    // every ion with any photoionisation data must include a cross-section from the ground state, so an excited
    // level having one without the ground state means the phixs data numbers levels differently from adata.txt.
    // Checked on the raw file contents (recorded before filtering on retained levels) so that a misnumbered ground
    // state table cannot be filtered out silently. A file relying on its companion phixs file for the ground state
    // table is allowed with a warning, since the combined dataset satisfies the invariant
    for (const int v : {1, 2}) {
      const int othv = 3 - v;
      for (int element = 0; element < get_nelements(); element++) {
        const int nions = get_nions(element);
        for (int ion = 0; ion < nions; ion++) {
          const auto uniqueionindex = get_uniqueionindex(element, ion);
          const auto minphixslevel = ion_minphixslevel_infile[v][uniqueionindex];
          if (minphixslevel == std::numeric_limits<int>::max() || minphixslevel == 0) {
            continue;  // the file has no tables for this ion, or it includes the ground state table
          }
          if (ion_minphixslevel_infile[othv][uniqueionindex] == 0) {
            printlnlog(
                "[warning] {}: Z={} ionstage {}: the file includes a cross-section for an excited level of this ion "
                "but the ground state table comes from {}. Check that the two files use the same level numbering",
                phixsdata_filenames[v], get_atomicnumber(element), get_ionstage(element, ion),
                phixsdata_filenames[othv]);
          } else {
            printlnlog(
                "[error] {}: Z={} ionstage {}: the file includes a cross-section for an excited level of this ion but "
                "none for the ground state. The phixs level numbering may not match the ground state index {} "
                "detected from adata.txt",
                phixsdata_filenames[v], get_atomicnumber(element), get_ionstage(element, ion), groundstate_index_in);
            assert_always(false);
          }
        }
      }
    }
  }

  // share the values parsed from the phixs file headers and the continuum totals with the other ranks in the node
  // (the node-shared level lists are synchronised by the barrier)
  MPI_Bcast_safe(globals::NPHIXSPOINTS, 0, globals::mpi_comm_node);
  MPI_Bcast_safe(globals::NPHIXSNUINCREMENT, 0, globals::mpi_comm_node);
  MPI_Bcast_safe(last_phixs_nuovernuedge, 0, globals::mpi_comm_node);
  MPI_Bcast_safe(globals::nbfcontinua, 0, globals::mpi_comm_node);
  MPI_Bcast_safe(globals::nbfcontinua_ground, 0, globals::mpi_comm_node);
  MPI_Barrier_node();

  // the writes are complete, so publish as read-only
  globals::alllevels.phixsstart = std::move(alllevels_phixsstart);
  globals::alllevels.nphixstargets = std::move(alllevels_nphixstargets);
  globals::alllevels.phixstargetstart = std::move(alllevels_phixstargetstart);

  int cont_index = 0;
  ptrdiff_t nbftables = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        const int nphixstargets = get_nphixstargets(element, ion, level);
        const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
        if (globals::rank_in_node == 0) {
          alllevels_bflist_start[uniquelevelindex] = (nphixstargets > 0) ? cont_index : -1;
        }
        cont_index += nphixstargets;
        if (nphixstargets > 0) {
          nbftables++;
        }
      }
    }
  }
  if (globals::rank_in_node == 0) {
    assert_always(cont_index == std::ssize(tmpallphixstargets));
  }
  MPI_Barrier_node();
  globals::alllevels.bflist_start = std::move(alllevels_bflist_start);

  if (cont_index > 0) {
    // copy the photoionisation tables into one contiguous block of memory
    globals::allphixs = MPI_shared_array<float>(nbftables * globals::NPHIXSPOINTS);
    auto allphixstargets_levelindex = MPI_shared_array<int>(cont_index);
    auto allphixstargets_probability = MPI_shared_array<double>(cont_index);

    if (globals::rank_in_node == 0) {
      assert_always((nbftables * globals::NPHIXSPOINTS) == std::ssize(tmpallphixs));
      std::ranges::copy(tmpallphixs, globals::allphixs.begin());

      for (int i = 0; i < std::ssize(tmpallphixstargets); i++) {
        allphixstargets_levelindex[i] = tmpallphixstargets[i].levelindex;
        allphixstargets_probability[i] = tmpallphixstargets[i].probability;
      }
    }

    MPI_Barrier_node();
    globals::allphixstargets_levelindex = std::move(allphixstargets_levelindex);
    globals::allphixstargets_probability = std::move(allphixstargets_probability);
  }

  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      // below is just an extra warning consistency check
      const int nlevels_groundterm = globals::elements[element].ions[ion].nlevels_groundterm;

      // all levels in the ground term should be photoionisation targets from the lower ground state
      if (ion > 0 && ion < get_nions(element) - 1) {
        const int nphixstargets = get_nphixstargets(element, ion - 1, 0);
        if (nphixstargets > 0 && get_phixsupperlevel(element, ion - 1, 0, 0) == 0) {
          const int phixstargetlevels = get_phixsupperlevel(element, ion - 1, 0, nphixstargets - 1) + 1;

          if (nlevels_groundterm != phixstargetlevels) {
            printlnlog("[warning] Z={} ionstage {} nlevels_groundterm {} phixstargetlevels(ion-1) {}.",
                       get_atomicnumber(element), get_ionstage(element, ion), nlevels_groundterm, phixstargetlevels);
          }
        }
      }
    }
  }
}

auto read_compositiondata() -> std::vector<int> {
  auto compositiondata = fstream_required("compositiondata.txt", std::ios::in);

  int nelements_in = 0;
  assert_always(compositiondata >> nelements_in);
  globals::elements.resize(nelements_in);

  // Unused fields that the compositiondata.txt format still carries: ion stages come from the per-element
  // ranges below, and abundances from the model file.
  int T_preset = 0;
  assert_always(compositiondata >> T_preset);
  assert_always(T_preset == 0);
  int homogeneous_abundances = 0;
  assert_always(compositiondata >> homogeneous_abundances);
  assert_always(homogeneous_abundances == 0);

  auto nlevelsmax_readin = std::vector<int>(get_nelements());
  auto nions_readin = std::vector<int>(get_nelements());
  int uniqueionindex = 0;  // index into list of all ions of all elements
  for (int element = 0; element < get_nelements(); element++) {
    // read information about the next element which should be stored to memory
    int atomicnumber = 0;
    int lowermost_ionstage = 0;
    int uppermost_ionstage = 0;
    // unused column; abundances come from the model file (still range-checked below)
    double uniformabundance{NAN};
    double mass_amu{NAN};
    assert_always(compositiondata >> atomicnumber >> nions_readin[element] >> lowermost_ionstage >>
                  uppermost_ionstage >> nlevelsmax_readin[element] >> uniformabundance >> mass_amu);
    printlnlog("compositiondata.txt: Z {} nions {} lowermost {} uppermost {} nlevelsmax {}", atomicnumber,
               nions_readin[element], lowermost_ionstage, uppermost_ionstage, nlevelsmax_readin[element]);
    assert_always(atomicnumber > 0);
    assert_always(nions_readin[element] >= 0);
    assert_always(nions_readin[element] == 0 ||
                  (nions_readin[element] == (uppermost_ionstage - lowermost_ionstage + 1)));
    assert_always(uniformabundance >= 0);
    assert_always(mass_amu >= 0);

    globals::elements[element] = {
        .ions = {},  // this will be set later after the total number of ions is known for the block allocation
        .anumber = atomicnumber,
        .lowest_ionstage = lowermost_ionstage,
        .uniqueionindexstart = uniqueionindex,
        .initstablemeannucmass = static_cast<float>(mass_amu * MH),
    };
    uniqueionindex += nions_readin[element];
  }

  globals::allions = MPI_shared_array<Ion>(uniqueionindex);

  for (int element = 0; element < get_nelements(); element++) {
    globals::elements[element].ions =
        std::span{globals::allions}.subspan(globals::elements[element].uniqueionindexstart, nions_readin[element]);
  }

  return nlevelsmax_readin;
}

void read_levels_and_transitions(std::vector<TempEnergyLevel>& temp_alllevels,
                                 std::vector<TempLineTransitionInput>& temp_linelist,
                                 std::vector<TempAllTransInput>& temp_alltranslist,
                                 const std::vector<int>& nlevelsmax_readin) {
  std::string line;
  std::istringstream ssline;
  globals::nlines = 0;
  std::vector<IonTransitionsInput> iontransitiontable;
  auto adata = fstream_required("adata.txt", std::ios::in);
  auto ftransitiondata = fstream_required("transitiondata.txt", std::ios::in);
  int uniquelevelindex = 0;  // index into list of all levels of all ions of all elements
  for (int element = 0; element < get_nelements(); element++) {
    // now read in data for all ions of the current element. before doing so initialize
    // energy scale for the current element (all level energies are stored relative to
    // the ground level of the neutral ion)
    const auto atomicnumber = globals::elements[element].anumber;
    double energyoffset_ev = 0.;
    double ionpot_ev = 0.;
    const auto nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int ionstage = globals::elements[element].lowest_ionstage + ion;
      // calculate the current levels ground level energy
      assert_always(ionpot_ev >= 0);
      energyoffset_ev += ionpot_ev;

      // read information for the elements next ionstage
      int adata_Z_in = -1;
      int adata_ionstage_in = -1;
      int nlevels_in_file = 0;

      while (adata_Z_in != atomicnumber || adata_ionstage_in != ionstage) {
        // skip over this ion block
        if (adata_Z_in == atomicnumber) {
          printlnlog("increasing energyoffset by ionpot {:g}", ionpot_ev);
          energyoffset_ev += ionpot_ev;
        }
        for (int i = 0; i < nlevels_in_file; i++) {
          std::getline(adata, line);
        }

        assert_always(get_noncommentline(adata, line));
        ssline.clear();
        ssline.str(line);
        assert_always(ssline >> adata_Z_in >> adata_ionstage_in >> nlevels_in_file >> ionpot_ev);
      }

      printlnlog("adata.txt: Z {} ionstage {} nlevels_in_file {}", adata_Z_in, adata_ionstage_in, nlevels_in_file);

      // optionally limit the top ion to one level and no transitions
      const int nlevelslimit = (SINGLE_LEVEL_TOP_ION && ion == (nions - 1)) ? 1 : nlevelsmax_readin[element];
      const int nlevelskept = (nlevelslimit < 0) ? nlevels_in_file : std::min(nlevelslimit, nlevels_in_file);

      if (nlevels_in_file > nlevelskept) {
        printlnlog("[info] read_levels_and_transitions: reduce number of levels from {} to {} for Z {:2} ionstage {}",
                   nlevels_in_file, nlevelskept, adata_Z_in, adata_ionstage_in);
      }
      assert_always(nlevelskept > 0);

      // read the data for the levels and set up the list of possible transitions for each level
      // store the ions data to memory and set up the ions zeta and levellist
      globals::elements[element].ions[ion] = {
          .nlevels = nlevelskept,
          .allnltelevelsindexstart = -1,
          .nlevels_ionising = 0,
          .maxrecombininglevel = -1,
          .nlevels_groundterm = -1,
          .uniquelevelindexstart = uniquelevelindex,
          .groundcontindex = -1,
          .ionpot = ionpot_ev * EV,
      };

      assert_always(std::ssize(temp_alllevels) == uniquelevelindex);

      read_ion_levels(adata, element, ion, nions, nlevels_in_file, nlevelskept, energyoffset_ev, ionpot_ev,
                      temp_alllevels);
      uniquelevelindex += get_nlevels(element, ion);

      // and proceed through the transitionlist till we match this ionstage (if it was not the neutral one)
      int transdata_Z_in = -1;
      int transdata_ionstage_in = -1;
      int ion_transition_count_in_file = 0;
      while (transdata_Z_in != atomicnumber || transdata_ionstage_in != ionstage) {
        // skip over table
        for (int i = 0; i < ion_transition_count_in_file; i++) {
          assert_always(getline(ftransitiondata, line));
        }
        assert_always(get_noncommentline(ftransitiondata, line));  // get_noncommentline to skip over blank lines
        ssline.clear();
        ssline.str(line);
        assert_always(ssline >> transdata_Z_in >> transdata_ionstage_in >> ion_transition_count_in_file);
      }

      assert_always(transdata_Z_in == atomicnumber);
      assert_always(transdata_ionstage_in == ionstage);

      printlnlog("transitiondata.txt: Z {} ionstage {} tottransitions {}", transdata_Z_in, transdata_ionstage_in,
                 ion_transition_count_in_file);
      assert_always(ion_transition_count_in_file >= 0);

      // load transition table for the current ion to temporary memory
      if (nlevelskept <= 1) {
        // we will not read in any transitions, just skip past these lines in the file
        for (int i = 0; i < ion_transition_count_in_file; i++) {
          assert_always(getline(ftransitiondata, line));
        }
      } else {
        // first nlevels_requiretransitions levels will be fully connected with higher levels
        // by addition of forbidden collisional transitions for any missing transitions
        const int nlevels_requiretransitions =
            std::min(nlevelskept, NLEVELS_REQUIRETRANSITIONS(atomicnumber, adata_ionstage_in));

        read_ion_transitions(ftransitiondata, ion_transition_count_in_file, iontransitiontable,
                             nlevels_requiretransitions, nlevelskept, atomicnumber, ionstage, nlevels_in_file);
        add_transitions_to_unsorted_linelist(element, ion, iontransitiontable, temp_linelist, temp_alltranslist,
                                             temp_alllevels);
      }
    }
  }
}

// sort the temporary linelist by descending frequency and warn about any duplicate transitions
void sort_temp_linelist(std::vector<TempLineTransitionInput>& temp_linelist) {
  if (globals::rank_in_node == 0) {
    assert_always(globals::nlines == std::ssize(temp_linelist));
    temp_linelist.shrink_to_fit();

    // sort the linelist by frequency descending
    std::SORT_OR_STABLE_SORT(temp_linelist.begin(), temp_linelist.end(), [](const auto& a, const auto& b) {
      if (a.nu != b.nu) {
        return a.nu > b.nu;
      }
      return std::tie(a.elementindex, a.ionindex, a.lowerlevelindex, a.upperlevelindex, a.einstein_A) <
             std::tie(b.elementindex, b.ionindex, b.lowerlevelindex, b.upperlevelindex, b.einstein_A);
    });

    for (int i = 0; i < globals::nlines - 1; i++) {
      const double nu = temp_linelist[i].nu;
      const double nu_next = temp_linelist[i + 1].nu;
      if (fabs(nu_next - nu) < (1.e-10 * nu)) {
        const auto& a1 = temp_linelist[i];
        const auto& a2 = temp_linelist[i + 1];

        if ((a1.elementindex == a2.elementindex) && (a1.ionindex == a2.ionindex) &&
            (a1.lowerlevelindex == a2.lowerlevelindex) && (a1.upperlevelindex == a2.upperlevelindex)) {
          printlnlog("Duplicate transition line? {}", a1.nu == a2.nu ? "nu match exact" : "close to nu match");
          printlnlog("a: Z={} ionstage {} lower {} upper {} nu {:g} [Hz] lambda {:g} [Angstrom]",
                     get_atomicnumber(a1.elementindex), get_ionstage(a1.elementindex, a1.ionindex), a1.lowerlevelindex,
                     a1.upperlevelindex, a1.nu, 1e8 * CLIGHT / a1.nu);
          printlnlog("b: Z={} ionstage {} lower {} upper {} nu {:g} [Hz] lambda {:g} [Angstrom]",
                     get_atomicnumber(a2.elementindex), get_ionstage(a2.elementindex, a2.ionindex), a2.lowerlevelindex,
                     a2.upperlevelindex, a2.nu, 1e8 * CLIGHT / a2.nu);
        }
      }
    }
  }
}

// copy the level data into node-shared memory arrays (globals::alllevels), freeing the local copy
void create_shared_levellist(std::vector<TempEnergyLevel>& temp_alllevels) {
  ptrdiff_t nlevels = std::ssize(temp_alllevels);
  MPI_Bcast_safe(nlevels, 0, globals::mpi_comm_node);

  auto alllevels_alltrans_startdown = MPI_shared_array<int>(nlevels);
  auto alllevels_ndowntrans = MPI_shared_array<int>(nlevels);
  auto alllevels_nuptrans = MPI_shared_array<int>(nlevels);
  auto alllevels_epsilon = MPI_shared_array<double>(nlevels);
  auto alllevels_statweight = MPI_shared_array<float>(nlevels);
  auto alllevels_matransblock_start = MPI_shared_array<int>(nlevels);
  // the remaining alllevels members are built and published by read_autoion_data(), read_phixs_data(),
  // and setup_phixs_list()
  if (globals::rank_in_node == 0) {
    int chtransindex = 0;
    for (auto i = 0ZU; i < temp_alllevels.size(); i++) {
      alllevels_alltrans_startdown[i] = temp_alllevels[i].alltrans_startdown;
      alllevels_ndowntrans[i] = temp_alllevels[i].ndowntrans;
      alllevels_nuptrans[i] = temp_alllevels[i].nuptrans;
      alllevels_epsilon[i] = temp_alllevels[i].epsilon;
      alllevels_statweight[i] = temp_alllevels[i].stat_weight;
      alllevels_matransblock_start[i] = chtransindex;
      chtransindex += ((2 * alllevels_ndowntrans[i]) + alllevels_nuptrans[i]);
    }
  }
  MPI_Barrier_node();
  globals::alllevels.alltrans_startdown = std::move(alllevels_alltrans_startdown);
  globals::alllevels.ndowntrans = std::move(alllevels_ndowntrans);
  globals::alllevels.nuptrans = std::move(alllevels_nuptrans);
  globals::alllevels.epsilon = std::move(alllevels_epsilon);
  globals::alllevels.statweight = std::move(alllevels_statweight);
  globals::alllevels.matransblock_start = std::move(alllevels_matransblock_start);
  temp_alllevels.clear();
  temp_alllevels.shrink_to_fit();
}

// copy the transition data into node-shared memory arrays (globals::alltrans), freeing the
// local copy, and point every up and down transition at its index in the frequency-sorted
// linelist. Requires globals::alllevels and globals::linelist to have been created already.
void create_shared_alltranslist(std::vector<TempAllTransInput>& temp_alltranslist) {
  const int updowntranscount = []() -> int {
    const int downtranscount = std::ranges::fold_left(globals::alllevels.ndowntrans, 0, std::plus<>{});
    const int uptranscount = std::ranges::fold_left(globals::alllevels.nuptrans, 0, std::plus<>{});

    printlnlog("total uptrans {}", uptranscount);
    printlnlog("total downtrans {}", downtranscount);
    return downtranscount + uptranscount;
  }();
  printlnlog("[info] mem_usage: transition lists occupy {:.3f} MB (node shared memory)",
             updowntranscount * ((2 * sizeof(int)) + (3 * sizeof(float)) + sizeof(bool)) / 1024. / 1024.);

  MPI_Barrier_node();

  auto alltrans_lineindex = MPI_shared_array<int>(updowntranscount);
  auto alltrans_targetlevelindex = MPI_shared_array<int>(updowntranscount);
  auto alltrans_einstein_A = MPI_shared_array<float>(updowntranscount);
  auto alltrans_coll_str = MPI_shared_array<float>(updowntranscount);
  auto alltrans_osc_strength = MPI_shared_array<float>(updowntranscount);
  auto alltrans_forbidden = MPI_shared_array<bool>(updowntranscount);

  if (globals::rank_in_node == 0) {
    assert_always(std::ssize(temp_alltranslist) == updowntranscount);
    for (int t = 0; t < updowntranscount; t++) {
      alltrans_targetlevelindex[t] = temp_alltranslist[t].targetlevelindex;
      alltrans_einstein_A[t] = temp_alltranslist[t].einstein_A;
      alltrans_coll_str[t] = temp_alltranslist[t].coll_str;
      alltrans_osc_strength[t] = temp_alltranslist[t].osc_strength;
      alltrans_forbidden[t] = temp_alltranslist[t].forbidden;
    }
  }
  temp_alltranslist.clear();
  temp_alltranslist.shrink_to_fit();

  globals::alltrans.targetlevelindex = std::move(alltrans_targetlevelindex);
  globals::alltrans.einstein_A = std::move(alltrans_einstein_A);
  globals::alltrans.coll_str = std::move(alltrans_coll_str);
  globals::alltrans.osc_strength = std::move(alltrans_osc_strength);
  globals::alltrans.forbidden = std::move(alltrans_forbidden);

  // make sure that all ranks can see the filled transition arrays before the striped loop below
  MPI_Barrier_node();

  printlog("establishing connection between transitions and sorted linelist...");

  const auto time_start_establish_linelist_connections = std::chrono::steady_clock::now();

  // every transition entry belongs to exactly one line (one down and one up entry per line), so
  // the loop below writes every element of alltrans_lineindex
  assert_always(updowntranscount == 2 * globals::nlines);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int lineindex = 0; lineindex < globals::nlines; lineindex++) {
    if (lineindex % globals::node_nprocs != globals::rank_in_node) {
      continue;
    }

    const int element = globals::linelist.elementindex[lineindex];
    const int ion = globals::linelist.ionindex[lineindex];
    const auto ionuniquelevelindexstart = get_ionuniquelevelindexstart(element, ion);
    const int upper_uniquelevelindex = globals::linelist.uniquelevelindex_upper[lineindex];
    const auto lower_uniquelevelindex = globals::linelist.uniquelevelindex_lower[lineindex];
    const auto upperlevel = upper_uniquelevelindex - ionuniquelevelindexstart;
    const auto lowerlevel = lower_uniquelevelindex - ionuniquelevelindexstart;

    // there is exactly one transition per pair of levels, and each level's down and up transition target lists
    // are sorted in ascending order (add_transitions_to_unsorted_linelist fills them by iterating the transition
    // table in (lower, upper) order), so the matching entries can be found by binary search

    const auto alltrans_startdown = get_alltrans_startdown(upper_uniquelevelindex);
    const auto ndowntrans = get_ndowntrans(upper_uniquelevelindex);
    const auto downtargets = globals::alltrans.targetlevelindex.subspan(alltrans_startdown, ndowntrans);
    const auto downit = std::ranges::lower_bound(downtargets, lowerlevel);
    assert_always(downit != downtargets.end() && *downit == lowerlevel);
    const auto downtransid = alltrans_startdown + static_cast<int>(downit - downtargets.begin());
    alltrans_lineindex[downtransid] = lineindex;

    const auto alltrans_startup = get_alltrans_startup(lower_uniquelevelindex);
    const auto nuptrans = get_nuptrans(lower_uniquelevelindex);
    const auto uptargets = globals::alltrans.targetlevelindex.subspan(alltrans_startup, nuptrans);
    const auto upit = std::ranges::lower_bound(uptargets, upperlevel);
    assert_always(upit != uptargets.end() && *upit == upperlevel);
    const auto uptransid = alltrans_startup + static_cast<int>(upit - uptargets.begin());
    alltrans_lineindex[uptransid] = lineindex;
  }
  globals::alltrans.lineindex = std::move(alltrans_lineindex);

  const auto establish_linelist_connections_duration =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - time_start_establish_linelist_connections)
          .count();
  printlnlog(" (took {:.1f} seconds)", establish_linelist_connections_duration);
  MPI_Barrier_node();
}

// copy the line data into node-shared memory arrays (globals::linelist), calculating the
// Einstein B coefficients, and free the local copy
void create_shared_linelist(std::vector<TempLineTransitionInput>& temp_linelist) {
  auto linelist_nu = MPI_shared_array<double>(globals::nlines);
  auto linelist_elementindex = MPI_shared_array<int>(globals::nlines);
  auto linelist_ionindex = MPI_shared_array<int>(globals::nlines);
  auto linelist_uniquelevelindex_lower = MPI_shared_array<int>(globals::nlines);
  auto linelist_uniquelevelindex_upper = MPI_shared_array<int>(globals::nlines);
  auto linelist_B_ul = MPI_shared_array<float>(globals::nlines);
  auto linelist_B_lu = MPI_shared_array<float>(globals::nlines);

  if (globals::rank_in_node == 0) {
    assert_always(std::ssize(temp_linelist) == globals::nlines);
    for (int t = 0; t < globals::nlines; t++) {
      linelist_nu[t] = temp_linelist[t].nu;
      linelist_elementindex[t] = temp_linelist[t].elementindex;
      linelist_ionindex[t] = temp_linelist[t].ionindex;

      const auto ionuniquelevelindexstart =
          get_ionuniquelevelindexstart(temp_linelist[t].elementindex, temp_linelist[t].ionindex);
      const int uniquelevelindex_lower = ionuniquelevelindexstart + temp_linelist[t].lowerlevelindex;
      const int uniquelevelindex_upper = ionuniquelevelindexstart + temp_linelist[t].upperlevelindex;
      linelist_uniquelevelindex_lower[t] = uniquelevelindex_lower;
      linelist_uniquelevelindex_upper[t] = uniquelevelindex_upper;

      const double B_ul = CLIGHTSQUAREDOVERTWOH / pow3(temp_linelist[t].nu) * temp_linelist[t].einstein_A;
      const auto f_B_ul = static_cast<float>(B_ul);
      assert_always(std::isfinite(f_B_ul));
      linelist_B_ul[t] = f_B_ul;

      const double B_lu = stat_weight(uniquelevelindex_upper) / stat_weight(uniquelevelindex_lower) * B_ul;
      const auto f_B_lu = static_cast<float>(B_lu);
      assert_always(std::isfinite(f_B_lu));
      linelist_B_lu[t] = f_B_lu;
    }
  }
  temp_linelist.clear();
  temp_linelist.shrink_to_fit();
  MPI_Barrier_node();

  globals::linelist.nu = std::move(linelist_nu);
  globals::linelist.elementindex = std::move(linelist_elementindex);
  globals::linelist.ionindex = std::move(linelist_ionindex);
  globals::linelist.uniquelevelindex_lower = std::move(linelist_uniquelevelindex_lower);
  globals::linelist.uniquelevelindex_upper = std::move(linelist_uniquelevelindex_upper);
  globals::linelist.B_ul = std::move(linelist_B_ul);
  globals::linelist.B_lu = std::move(linelist_B_lu);

  const double linelist_mem_MB =
      ((globals::nlines * sizeof(double))  // nu
       + (globals::nlines * sizeof(int) * 4)  // elementindex, ionindex, uniquelevelindex_lower, uniquelevelindex_upper
       + (globals::nlines * sizeof(float) * 2)  // B_ul, B_lu
       ) /
      1024. / 1024;

  printlnlog("[info] mem_usage: linelist occupies {:.3f} MB (node shared memory)", linelist_mem_MB);
}

// read the full atomic dataset: the composition (compositiondata.txt), then the levels and
// transitions (adata.txt, transitiondata.txt) and photoionisation cross-sections
// (phixsdata_v2.txt) for every included ion, building the global element/ion/level/line lists
void read_atomicdata_files() {
  const auto nlevelsmax_readin = read_compositiondata();

  printlnlog("SINGLE_LEVEL_TOP_ION: {}", SINGLE_LEVEL_TOP_ION ? "true" : "false");

  std::vector<TempEnergyLevel> temp_alllevels;
  std::vector<TempLineTransitionInput> temp_linelist;
  std::vector<TempAllTransInput> temp_alltranslist;

  if (globals::rank_in_node == 0) {
    temp_linelist.reserve(1U << 22U);  // reserve initial space for 4 million lines to avoid too many reallocations
    temp_alltranslist.reserve(1U << 22U);
    read_levels_and_transitions(temp_alllevels, temp_linelist, temp_alltranslist, nlevelsmax_readin);
  }
  MPI_Barrier_node();

  // only the node leaders read adata.txt, but all ranks parse level indices in autoion.txt and the phixs files
  MPI_Bcast_safe(groundstate_index_in, 0, globals::mpi_comm_node);
  assert_always(groundstate_index_in == 0 || groundstate_index_in == 1);

  update_includedionslevels_maxnions();

  MPI_Bcast_safe(globals::nlines, 0, globals::mpi_comm_node);
  printlnlog("nlines {}", globals::nlines);

  sort_temp_linelist(temp_linelist);

  // create a shared level list and copy data across, freeing the local copy
  create_shared_levellist(temp_alllevels);

  // create a linelist shared on node and then copy data across, freeing the local copy
  create_shared_linelist(temp_linelist);

  // create the shared transitions list and point each transition at the sorted linelist
  create_shared_alltranslist(temp_alltranslist);

  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      globals::elements[element].ions[ion].nlevels_groundterm = calculate_nlevels_groundterm(element, ion);
    }
  }

  read_autoion_data();
  read_phixs_data();
  setup_phixs_list();
}

void write_bflist_file() {
  reserve_resize(globals::bflist, globals::nbfcontinua);

  int i = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels_ionising(element, ion);
      for (int level = 0; level < nlevels; level++) {
        const auto nphixstargets = get_nphixstargets(element, ion, level);
        for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
          globals::bflist[i].elementindex = element;
          globals::bflist[i].ionindex = ion;
          globals::bflist[i].levelindex = level;
          globals::bflist[i].phixstargetindex = phixstargetindex;

          const int et = -1 - i;

          assert_always(et == get_emtype_continuum(element, ion, level, phixstargetindex));

          // check the we don't overload the same packet emission type numbers
          // as the special values for free-free scattering and not set
          assert_always(et != EMTYPE_NOTSET);
          assert_always(et != EMTYPE_FREEFREE);
          i++;
        }
      }
    }
  }
  assert_always(i == globals::nbfcontinua);

  if (globals::my_rank == 0) {
    auto bflist_file = fstream_required("bflist.out", std::ios::out | std::ios::trunc);
    std::println(bflist_file, "{}", globals::nbfcontinua);
    for (i = 0; i < globals::nbfcontinua; i++) {
      const int element = globals::bflist[i].elementindex;
      const int ion = globals::bflist[i].ionindex;
      const int level = globals::bflist[i].levelindex;
      const int phixstargetindex = globals::bflist[i].phixstargetindex;
      const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
      std::println(bflist_file, "{} {} {} {} {}", i, element, ion, level, upperionlevel);
    }
  }
}

void setup_nlte_levels() {
  globals::total_nlte_levels = 0;
  int n_super_levels = 0;

  for (int element = 0; element < get_nelements(); element++) {
    globals::elements[element].has_nlte_levels = elem_has_nlte_levels_search(element);
  }

  for (int element = 0; element < get_nelements(); element++) {
    if (elem_has_nlte_levels(element)) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        globals::elements[element].ions[ion].allnltelevelsindexstart = globals::total_nlte_levels;
        const int nlevels = get_nlevels(element, ion);
        int nlevels_excited_nlte = 0;
        bool found_lte_only_level = false;
        for (int level = 1; level < nlevels; level++) {
          if (is_nlte(element, ion, level)) {
            nlevels_excited_nlte++;
            globals::total_nlte_levels++;
            assert_always(found_lte_only_level == false);  // NLTE levels must be consecutive
          } else {
            found_lte_only_level = true;
          }
        }
        globals::elements[element].ions[ion].nlevels_excited_nlte = nlevels_excited_nlte;

        // use the same definition as ion_has_superlevel(): autoionising levels get their own
        // NLTE-solver slots, so they must not count towards needing a superlevel
        const bool has_superlevel = ion_has_superlevel(element, ion);
        if (has_superlevel) {
          // If there are more levels that the ground state + the number of NLTE levels then we need an extra
          // slot to store data for the "superlevel", which is a representation of all the other levels that
          // are not treated in detail.
          globals::total_nlte_levels++;
          n_super_levels++;
        }

        printlnlog("[info]  element {:2} Z={:2} ionstage {:2} has {:5} NLTE excited levels{}. Starting at index {}",
                   element, get_atomicnumber(element), get_ionstage(element, ion),
                   get_nlevels_excited_nlte(element, ion), has_superlevel ? " plus a superlevel" : "",
                   get_allnltelevelsindexstart(element, ion));
      }
    }
  }

  printlnlog("[info] Total NLTE levels: {}, of which {} are superlevels", globals::total_nlte_levels, n_super_levels);
}

}  // anonymous namespace

// read input parameters from input.txt
void read_parameterfile(std::span<Packet> packets) {
  auto file = fstream_required("input.txt", std::ios::in);

  std::string line;
  assert_always(get_noncommentline(file, line));

  std::int64_t pre_zseed = -1;
  std::istringstream{line} >> pre_zseed;

  if (pre_zseed > 0) {
    printlnlog("input.txt specified random number seed is {}", pre_zseed);
  } else {
#if defined REPRODUCIBLE && REPRODUCIBLE
    printlnlog(
        "[error] REPRODUCIBLE mode requires a positive random number seed on the first non-comment line of input.txt "
        "(found {})",
        pre_zseed);
    std::abort();
#endif
    pre_zseed = get_rng_random_seed();
    // broadcast randomly-generated seed from rank 0 to all ranks
    MPI_Bcast_safe(pre_zseed, 0, MPI_COMM_WORLD);
    printlnlog("randomly-generated random number seed is {}", pre_zseed);
  }

  if (!packets.empty()) {
#ifdef GPU_ON
    // Give every packet its own independently-seeded generator. The ranks are spaced by the number
    // of packets that they own, so that their seed ranges do not overlap. get_max_threads() is one
    // for a GPU build, so spacing the ranks by 13 as the host generator below does would leave
    // neighbouring ranks sharing all but 13 of their seeds, and two packets given the same seed
    // follow identical histories because the grid state that they see is rank-invariant.
    // Xoshiro128PP is seeded from a 32 bit value, so distinct seeds only exist for as many packets
    // as fit in that space and the whole run has to stay within it.
    assert_always((static_cast<std::int64_t>(globals::nprocs) * std::ssize(packets)) <= (1LL << 32));
    const auto rank_seed_base =
        static_cast<std::uint32_t>(pre_zseed + (static_cast<std::int64_t>(globals::my_rank) * std::ssize(packets)));
    for (auto packetnumber = 0ZU; packetnumber < std::size(packets); packetnumber++) {
      get_rngstate(packets[packetnumber]).seed(rank_seed_base + static_cast<std::uint32_t>(packetnumber));
    }
    printlnlog("rank {}: packet rngseeds start at {}", globals::my_rank, rank_seed_base);
#else
    // For MPI parallelisation, the random seed is changed based on the rank of the process.
    // This runs on the main thread only; OpenMP/stdpar worker threads get their own randomly-seeded
    // thread_local generators on first use (see get_rngstate()), so multi-threaded runs are not
    // reproducible (they also accumulate to shared memory in a non-deterministic order)
    const auto rngseed = pre_zseed + static_cast<std::int64_t>(13 * globals::my_rank * get_max_threads());
    get_rngstate().seed(rngseed);
    for (int n = 0; n < 100; n++) {
      rng_uniform(get_rngstate());
    }
    printlnlog("rank {}: main thread has rngseed {}", globals::my_rank, rngseed);
#endif
  }

  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::ntimesteps;  // number of time steps
  assert_always(globals::ntimesteps > 0);

  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::timestep_initial >>
      globals::timestep_finish;  // number of start and end time step
  printlnlog("input: timestep_start {} timestep_finish {}", globals::timestep_initial, globals::timestep_finish);
  if (globals::timestep_finish < globals::ntimesteps) {
    printlnlog(
        "input: note that the simulation will stop after timestep {}, before the final timestep {}. Restarts advance "
        "the start timestep but never timestep_finish, so continuing past it will require editing input.txt",
        globals::timestep_finish - 1, globals::ntimesteps - 1);
  }
  assert_always(globals::timestep_initial < globals::ntimesteps);
  assert_always(globals::timestep_initial <= globals::timestep_finish);
  assert_always(globals::timestep_finish <= globals::ntimesteps);

  double tmin_days = 0.;
  double tmax_days = 0.;
  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> tmin_days >> tmax_days;  // start and end times
  assert_always(tmin_days > 0);
  assert_always(tmax_days > 0);
  assert_always(tmin_days < tmax_days);
  globals::tmin = tmin_days * DAY;
  globals::tmax = tmax_days * DAY;

  assert_always(get_noncommentline(file, line));  // UNUSED nusyn_min_mev nusyn_max_mev

  assert_always(get_noncommentline(file, line));  // UNUSED number of times for synthesis

  assert_always(get_noncommentline(file, line));  // UNUSED start and end times for synthesis

  assert_always(get_noncommentline(file, line));  // UNUSED model dimensions (now autodetected)

  assert_always(get_noncommentline(file, line));  // UNUSED compute the r-light curve?

  assert_always(get_noncommentline(file, line));  // UNUSED number of iterations

  assert_always(get_noncommentline(file, line));  // UNUSED change speed of light

  assert_always(get_noncommentline(file, line));  // UNUSED gamma_kappagrey

  assert_always(get_noncommentline(file, line));  // UNUSED components of syn_dir

  assert_always(get_noncommentline(file, line));  // UNUSED opacity choice

  assert_always(get_noncommentline(file, line));  // UNUSED free parameter for calculation of rho_crit

  assert_always(get_noncommentline(file, line));  // UNUSED activate debug output for packet

  // Do we start a new simulation or, continue another one?
  int continue_flag = 0;
  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> continue_flag;
  globals::simulation_continued_from_saved = (continue_flag == 1);
  if (globals::timestep_initial == 0) {
    // it's not possible to resume from a saved point if we start from timestep zero, so override the flag
    globals::simulation_continued_from_saved = false;
  }
  if (globals::simulation_continued_from_saved) {
    printlnlog("input: resuming simulation from saved point");
  } else {
    printlnlog("input: starting a new simulation");
  }

  // UNUSED rfcut parameter (kept for backward-compatible input parsing)
  assert_always(get_noncommentline(file, line));

  // Sets the number of initial LTE timesteps for NLTE runs
  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::num_lte_timesteps;
  printlnlog("input: doing the first {} timesteps in LTE", globals::num_lte_timesteps);

  if (NT_ON) {
    if (NT_SOLVE_SPENCERFANO) {
      printlnlog("NT_ON = true: Non-thermal ionisation with a Spencer-Fano solution is switched on for this run.");
    } else {
      printlnlog(
          "NT_ON = true: Non-thermal ionisation with the work function approximation is switched on for this run.");
    }
  } else {
    printlnlog("NT_ON = false: No non-thermal ionisation is used in this run.");
  }

  if (USE_LUT_PHOTOION) {
    printlnlog("Corrphotoioncoeff is calculated from LTE values and corrphotoionrenorm estimator.");
  } else {
    printlnlog(
        "Corrphotoioncoeff is calculated from the radiation field at each timestep in each modelgrid cell (no LUT).");
  }

  if (USE_ION_BFHEATING_ESTIMATORS) {
    printlnlog("bfheating coefficients are calculated from LTE values and bfheatingestimator.");
  } else {
    printlnlog("bfheating coefficients are calculated directly from the radiation field without bfheatingestimator.");
  }

  // Set up initial grey approximation?
  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::optical_depth_is_thick >> globals::num_grey_timesteps;
  printlnlog(
      "input: cells with Thomson optical depth > {:g} are treated in grey approximation for the first {} timesteps",
      globals::optical_depth_is_thick, globals::num_grey_timesteps);

  assert_always(get_noncommentline(file, line));  // UNUSED max_bf_continua (all bf continua are always included)

  // for exspec: read number of MPI tasks
  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::nprocs_exspec;

  // UNUSED: Extract line-of-sight dependent information of last emission for spectrum_res
  assert_always(get_noncommentline(file, line));

  // UNUSED: kpkt diffusion parameters: now set in kpkt.cc
  assert_always(get_noncommentline(file, line));

  file.close();

  if (globals::my_rank == 0 && !globals::simulation_continued_from_saved) {
    // back up original input file, adding comments to each line
    update_parameterfile(-1);
  }
}

// write out an updated input.txt to restart the simulation
void update_parameterfile(const int nts) {
  assert_always(globals::my_rank == 0);
  if (nts >= 0) {
    printlog("Update input.txt for restart at timestep {}...", nts);
  } else {
    printlog("Copying input.txt to input-newrun.txt...");
  }

  auto file = fstream_required("input.txt", std::ios::in);

  auto fileout = fstream_required("input.txt.tmp", std::ios::out | std::ios::trunc);

  std::string line;

  int noncomment_linenum = -1;
  while (std::getline(file, line)) {
    if (!lineiscommentonly(line)) {
      noncomment_linenum++;  // line number starting from 0, ignoring comment and blank lines (that start with '#')

      // overwrite particular lines to enable restarting from the current timestep
      if (nts >= 0) {
        if (noncomment_linenum == inputline_timestep_range) {
          // Number of start and end time step
          line = std::format("{:03d} {:03d}", nts, globals::timestep_finish);
        } else if (noncomment_linenum == inputline_continue_from_saved) {
          // resume from gridsave file
          line = "1";  // Force continuation
        }
      }

      // only rewrite this line when updating input.txt for a restart (sn3d), where nprocs is the
      // number of packet files being written. The nts == -1 backup path may be run by exspec
      // (nprocs == 1), which must not clobber the nprocs_exspec value it just read
      if (nts >= 0 && noncomment_linenum == inputline_nprocs_exspec) {
        // by default, exspec should use all available packet files
        globals::nprocs_exspec = globals::nprocs;
        line = std::format("{}", globals::nprocs_exspec);
      }

      if (noncomment_linenum < std::ssize(inputlinecomments)) {
        const int commentstart = 25;

        // truncate any existing comment on the line
        if (line.contains('#')) {
          line.resize(line.find('#'));
        }

        // strip trailing whitespace so that the comment column doesn't shift on repeated rewrites
        line.resize(line.find_last_not_of(" \t") + 1);

        // pad the data field out to the comment column, but never truncate it
        if (std::ssize(line) < commentstart) {
          line.resize(commentstart, ' ');
        } else {
          line.push_back(' ');
        }
        line.append("# ");
        line.append(inputlinecomments[noncomment_linenum]);
      }
    }

    std::println(fileout, "{}", line);
  }

  fileout.close();
  file.close();

  std::error_code rename_error;
  if (nts < 0) {
    // back up the original for starting a new simulation
    std::filesystem::rename("input.txt.tmp", "input-newrun.txt", rename_error);
  } else {
    std::filesystem::rename("input.txt.tmp", "input.txt", rename_error);
  }
  if (rename_error) {
    printlnlog("[error] failed to move input.txt.tmp to {}: {}", (nts < 0) ? "input-newrun.txt" : "input.txt",
               rename_error.message());
  }

  printlnlog("done");
}

void read_atomicdata() {
  read_atomicdata_files();

  kpkt::setup_coolinglist();

  // Printout some information about the read-in model atom

  int includedionisinglevels = 0;
  int includedboundboundtransitions = 0;
  int includedphotoiontransitions = 0;
  for (int element = 0; element < get_nelements(); element++) {
    printlnlog("[info]  element {} (Z={:2} {})", element, get_atomicnumber(element),
               decay::get_elname(get_atomicnumber(element)));
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      int ion_photoiontransitions = 0;
      int ion_bbtransitions = 0;
      for (int level = 0; level < get_nlevels(element, ion); level++) {
        ion_photoiontransitions += get_nphixstargets(element, ion, level);
        ion_bbtransitions += get_nuptrans(element, ion, level);
      }

      printlnlog(
          "[info]    ionstage {}: {:4} levels ({:4} ionising) {:7} lines {:6} bf transitions ("
          "epsilon_ground: {:7.2f} [eV])",
          get_ionstage(element, ion), get_nlevels(element, ion), get_nlevels_ionising(element, ion), ion_bbtransitions,
          ion_photoiontransitions, epsilon(element, ion, 0) / EV);

      includedionisinglevels += get_nlevels_ionising(element, ion);
      includedphotoiontransitions += ion_photoiontransitions;
      includedboundboundtransitions += ion_bbtransitions;
    }
  }
  assert_always(includedphotoiontransitions == globals::nbfcontinua);
  assert_always(globals::nlines == includedboundboundtransitions);

  printlnlog("[info]  in total {} ions, {} levels ({} ionising), {} lines, {} photoionisation transitions",
             get_includedions(), get_includedlevels(), includedionisinglevels, globals::nlines, globals::nbfcontinua);

  write_bflist_file();

  setup_nlte_levels();
}

// the pure timestep schemes ignore these two values, and the presets ship them as -1.
static_assert(TIMESTEP_SIZE_METHOD == TimeStepSizeMethod::LOGARITHMIC ||
                  TIMESTEP_SIZE_METHOD == TimeStepSizeMethod::CONSTANT ||
                  (FIXED_TIMESTEP_WIDTH > 0. && TIMESTEP_TRANSITION_TIME > 0.),
              "The hybrid TIMESTEP_SIZE_METHOD options require FIXED_TIMESTEP_WIDTH and TIMESTEP_TRANSITION_TIME "
              "(both in days) to be set to positive values");

// initialise the time steps
void setup_timesteps() {
  // t=globals::tmin is the start of the calculation. t=globals::tmax is the end of the calculation.
  // globals::ntimesteps is the number of time steps

  globals::timesteps.resize(globals::ntimesteps + 1);

  // Now set the individual time steps
  switch (TIMESTEP_SIZE_METHOD) {
    case TimeStepSizeMethod::LOGARITHMIC: {
      for (int n = 0; n < globals::ntimesteps; n++) {  // For logarithmic steps, the logarithmic interval will be
        const double dlogt = (log(globals::tmax) - log(globals::tmin)) / globals::ntimesteps;
        globals::timesteps[n].start = globals::tmin * exp(n * dlogt);
        globals::timesteps[n].mid = globals::tmin * exp((n + 0.5) * dlogt);
        globals::timesteps[n].width = (globals::tmin * exp((n + 1) * dlogt)) - globals::timesteps[n].start;
      }
      break;
    }

    case TimeStepSizeMethod::CONSTANT: {
      for (int n = 0; n < globals::ntimesteps; n++) {
        // for constant timesteps
        const double dt = (globals::tmax - globals::tmin) / globals::ntimesteps;
        globals::timesteps[n].start = globals::tmin + (n * dt);
        globals::timesteps[n].width = dt;
        globals::timesteps[n].mid = globals::timesteps[n].start + (0.5 * globals::timesteps[n].width);
      }
      break;
    }

    case TimeStepSizeMethod::LOGARITHMIC_THEN_CONSTANT: {
      // First part log, second part fixed timesteps
      const double t_transition = TIMESTEP_TRANSITION_TIME * DAY;  // transition from logarithmic to fixed timesteps
      const double maxtsdelta = FIXED_TIMESTEP_WIDTH * DAY;  // maximum timestep width in fixed part
      assert_always(t_transition > globals::tmin);
      assert_always(t_transition < globals::tmax);
      const int nts_fixed = ceil((globals::tmax - t_transition) / maxtsdelta);
      const double fixed_tsdelta = (globals::tmax - t_transition) / nts_fixed;
      assert_always(nts_fixed > 0);
      assert_always(nts_fixed < globals::ntimesteps);
      const int nts_log = globals::ntimesteps - nts_fixed;
      assert_always(nts_log > 0);
      assert_always(nts_log < globals::ntimesteps);
      assert_always((nts_log + nts_fixed) == globals::ntimesteps);
      for (int n = 0; n < globals::ntimesteps; n++) {
        if (n < nts_log) {
          // For logarithmic steps, the logarithmic interval will be
          const double dlogt = (log(t_transition) - log(globals::tmin)) / nts_log;
          globals::timesteps[n].start = globals::tmin * exp(n * dlogt);
          globals::timesteps[n].mid = globals::tmin * exp((n + 0.5) * dlogt);
          globals::timesteps[n].width = (globals::tmin * exp((n + 1) * dlogt)) - globals::timesteps[n].start;
        } else {
          // for constant timesteps
          const double prev_start =
              n > 0 ? (globals::timesteps[n - 1].start + globals::timesteps[n - 1].width) : globals::tmin;
          globals::timesteps[n].start = prev_start;
          globals::timesteps[n].width = fixed_tsdelta;
          globals::timesteps[n].mid = globals::timesteps[n].start + (0.5 * globals::timesteps[n].width);
        }
      }
      break;
    }

    case TimeStepSizeMethod::CONSTANT_THEN_LOGARITHMIC: {
      // First part fixed timesteps, second part log timesteps
      const double t_transition = TIMESTEP_TRANSITION_TIME * DAY;  // transition from fixed to logarithmic timesteps
      const double maxtsdelta = FIXED_TIMESTEP_WIDTH * DAY;  // timestep width of fixed timesteps
      assert_always(t_transition > globals::tmin);
      assert_always(t_transition < globals::tmax);
      const int nts_fixed = ceil((t_transition - globals::tmin) / maxtsdelta);
      const double fixed_tsdelta = (t_transition - globals::tmin) / nts_fixed;
      assert_always(nts_fixed > 0);
      assert_always(nts_fixed < globals::ntimesteps);
      const int nts_log = globals::ntimesteps - nts_fixed;
      assert_always(nts_log > 0);
      assert_always(nts_log < globals::ntimesteps);
      assert_always((nts_log + nts_fixed) == globals::ntimesteps);
      for (int n = 0; n < globals::ntimesteps; n++) {
        if (n < nts_fixed) {
          // for constant timesteps
          globals::timesteps[n].start = globals::tmin + (n * fixed_tsdelta);
          globals::timesteps[n].width = fixed_tsdelta;
          globals::timesteps[n].mid = globals::timesteps[n].start + (0.5 * globals::timesteps[n].width);
        } else {
          // For logarithmic time steps, the logarithmic interval will be
          const double dlogt = (log(globals::tmax) - log(t_transition)) / nts_log;
          const double prev_start =
              n > 0 ? (globals::timesteps[n - 1].start + globals::timesteps[n - 1].width) : globals::tmin;
          globals::timesteps[n].start = prev_start;
          globals::timesteps[n].width = (t_transition * exp((n - nts_fixed + 1) * dlogt)) - globals::timesteps[n].start;
          globals::timesteps[n].mid = globals::timesteps[n].start + (0.5 * globals::timesteps[n].width);
        }
      }
      break;
    }

    default:
      assert_always(false);
  }

  // and add a dummy timestep which contains the endtime
  // of the calculation
  globals::timesteps[globals::ntimesteps].start = globals::tmax;
  globals::timesteps[globals::ntimesteps].mid = globals::tmax;
  globals::timesteps[globals::ntimesteps].width = 0.;

  // check consistency of start + width = start_next
  for (int n = 1; n < globals::ntimesteps; n++) {
    const auto tsprev_end = globals::timesteps[n - 1].start + globals::timesteps[n - 1].width;
    assert_always(fabs((tsprev_end / globals::timesteps[n].start) - 1.) < 0.001);
  }
  const auto tsfinal_end =
      globals::timesteps[globals::ntimesteps - 1].start + globals::timesteps[globals::ntimesteps - 1].width;
  assert_always(fabs((tsfinal_end / globals::tmax) - 1.) < 0.001);

  const auto* const method_name = [] {
    switch (TIMESTEP_SIZE_METHOD) {
      case TimeStepSizeMethod::LOGARITHMIC:
        return "logarithmic";
      case TimeStepSizeMethod::CONSTANT:
        return "constant";
      case TimeStepSizeMethod::LOGARITHMIC_THEN_CONSTANT:
        return "logarithmic then constant";
      case TimeStepSizeMethod::CONSTANT_THEN_LOGARITHMIC:
        return "constant then logarithmic";
    }
    return "unknown";
  }();
  printlnlog(
      "timesteps: {} steps from tmin {:.4f} [d] to tmax {:.4f} [d] with {} sizing (first width {:.4f} [d], last "
      "width {:.4f} [d])",
      globals::ntimesteps, globals::tmin / DAY, globals::tmax / DAY, method_name, globals::timesteps[0].width / DAY,
      globals::timesteps[globals::ntimesteps - 1].width / DAY);
}
