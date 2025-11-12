#include "input.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <ios>
#include <istream>
#include <iterator>
#include <limits>
#include <ranges>
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "kpkt.h"
#include "packet.h"
#include "random.h"
#include "ratecoeff.h"
#include "sn3d.h"
#include "vpkt.h"

namespace {

const int groundstate_index_in = 1;  // starting level index in the input files

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
  int lineindex;
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

constexpr std::array<std::string_view, 24> inputlinecomments = {
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
    "10: UNUSED: change speed of light by some factor. Change constants.h CLIGHT_PROP instead",
    "11: gamma_kappagrey: if >0: use grey opacity for gammas, if <0: use detailed opacity",
    "12: UNUSED: syn_dir: x, y, and z components of unit vector (now always 0,0,1)",
    "13: opacity_case: opacity choice",
    "14: rho_crit_para: free parameter for calculation of rho_crit",
    "15: UNUSED debug_packet: (>=0: activate debug output for packet id, <0: ignore)",
    "16: simulation_continued_from_saved: (0: start new simulation, 1: continue from gridsave and packets files)",
    "17: UNUSED rfcut_angstroms: wavelength (in Angstroms) at which the parameterisation of the radiation field "
    "switches from the nebular approximation to LTE.",
    "18: num_lte_timesteps",
    "19: cell_is_optically_thick num_grey_timesteps",
    "20: UNUSED max_bf_continua: (>0: max bound-free continua per ion, <0 unlimited)",
    "21: nprocs_exspec: extract spectra for n MPI tasks. sn3d will set this on start of new sim.",
    "22: UNUSED do_emission_res: this is always true for exspec, sometimes true during sn3d",
    "23: kpktdiffusion_timescale n_kpktdiffusion_timesteps: kpkts diffuse x of a time step's length for the first y "
    "time steps"};

void read_phixs_data_table(std::istream& phixsfile, const int nphixspoints_inputtable, const int element,
                           const int lowerion, const int lowerlevel, const int upperion, int upperlevel_in,
                           std::vector<float>& tmpallphixs, std::vector<PhotoionTarget>& tmpallphixstargets,
                           size_t* mem_usage_phixs, const int phixs_file_version) {
  std::string phixsline;
  const auto phixstargetstart = static_cast<int>(tmpallphixstargets.size());
  const auto lowerionlower_uniquelevelindex = get_uniquelevelindex(element, lowerion, lowerlevel);
  assert_always(globals::alllevels.phixstargetstart[lowerionlower_uniquelevelindex] == -1 ||
                globals::alllevels.phixstargetstart[lowerionlower_uniquelevelindex] == phixstargetstart);
  globals::alllevels.phixstargetstart[lowerionlower_uniquelevelindex] = phixstargetstart;
  if (upperlevel_in >= 0) {  // file gives photoionisation to a single target state only
    int upperlevel = upperlevel_in - groundstate_index_in;
    assert_always(upperlevel >= 0);
    assert_always(globals::alllevels.nphixstargets[lowerionlower_uniquelevelindex] == 0 ||
                  globals::alllevels.nphixstargets[lowerionlower_uniquelevelindex] == 1);
    globals::alllevels.nphixstargets[lowerionlower_uniquelevelindex] = 1;

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
      assert_always(globals::alllevels.nphixstargets[lowerionlower_uniquelevelindex] == 0 ||
                    globals::alllevels.nphixstargets[lowerionlower_uniquelevelindex] == in_nphixstargets);
      globals::alllevels.nphixstargets[lowerionlower_uniquelevelindex] = in_nphixstargets;

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
        printlog("ERROR: photoionisation table for Z={} ionstage {} has probabilities that sum to {:g}",
                 get_atomicnumber(element), get_ionstage(element, lowerion), probability_sum);
        assert_always(false);
      }
    } else {  // file has table of target states and probabilities but our top ion is limited to one level
      globals::alllevels.nphixstargets[lowerionlower_uniquelevelindex] = 1;

      for (int i = 0; i < in_nphixstargets; i++) {
        assert_always(get_noncommentline(phixsfile, phixsline));
      }

      // send it to the ground state of the top ion
      tmpallphixstargets.push_back({.probability = 1., .levelindex = 0});
    }
  }

  assert_always(tmpallphixs.size() % globals::NPHIXSPOINTS == 0);
  const auto tmpphixsstart = tmpallphixs.size();
  globals::alllevels.phixsstart[lowerionlower_uniquelevelindex] =
      static_cast<int>(tmpphixsstart / globals::NPHIXSPOINTS);
  tmpallphixs.resize(tmpallphixs.size() + globals::NPHIXSPOINTS);

  const auto levelphixstable = std::span{tmpallphixs}.subspan(tmpphixsstart, globals::NPHIXSPOINTS);
  if (phixs_file_version == 1) {
    assert_always(get_nphixstargets(element, lowerion, lowerlevel) == 1);
    assert_always(tmpallphixstargets[get_allphixstargetindex(lowerionlower_uniquelevelindex, 0)].levelindex == 0);

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

    gsl_interp_accel* acc = gsl_interp_accel_alloc();
    gsl_spline* spline = gsl_spline_alloc(gsl_interp_linear, nphixspoints_inputtable);
    gsl_spline_init(spline, nugrid_in.data(), phixs_in.data(), nphixspoints_inputtable);
    for (int i = 1; i < globals::NPHIXSPOINTS; i++) {
      const double nu = nu_edge * (1. + i * globals::NPHIXSNUINCREMENT);
      if (nu > nu_max) {
        levelphixstable[i] = static_cast<float>(phixs_in[nphixspoints_inputtable - 1] * pow(nu_max / nu, 3));
      } else {
        levelphixstable[i] = static_cast<float>(gsl_spline_eval(spline, nu, acc));
      }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
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

  // The level contributes to the ionisinglevels if its energy
  // is below the ionisation potential and the level doesn't
  // belong to the topmost ion included.
  // Rate coefficients are only available for ionising levels.
  //  also need (levelenergy < ionpot && ...)?
  if (lowerion < get_nions(element) - 1) {
    for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, lowerion, lowerlevel);
         phixstargetindex++) {
      const int upperlevel =
          tmpallphixstargets[get_allphixstargetindex(lowerionlower_uniquelevelindex, phixstargetindex)].levelindex;
      if (upperlevel > get_maxrecombininglevel(element, lowerion + 1)) {
        globals::elements[element].ions[lowerion + 1].maxrecombininglevel = upperlevel;
      }
    }
  }

  *mem_usage_phixs += (get_nphixstargets(lowerionlower_uniquelevelindex) * sizeof(double)) + sizeof(int);
  *mem_usage_phixs += globals::NPHIXSPOINTS * sizeof(float);
  globals::nbfcontinua += get_nphixstargets(lowerionlower_uniquelevelindex);
  if (lowerlevel == 0 && get_nphixstargets(lowerionlower_uniquelevelindex) > 0) {
    globals::nbfcontinua_ground++;
  }
}

void read_phixs_file(const int phixs_file_version, std::vector<float>& tmpallphixs,
                     std::vector<PhotoionTarget>& tmpallphixstargets) {
  printlnlog("readin phixs data from {}", phixsdata_filenames[phixs_file_version]);

  auto phixsfile = fstream_required(phixsdata_filenames[phixs_file_version], std::ios::in);
  std::string phixsline;
  std::istringstream ssline;
  auto mem_usage_phixs = 0zU;

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
  double phixs_threshold_ev = -1;  // currently just ignored, and epilson is used instead
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
      // translate readin ionstages to ion indices

      const int upperion = upperionstage - get_ionstage(element, 0);
      const int lowerion = lowerionstage - get_ionstage(element, 0);
      const int lowerlevel = lowerlevel_in - groundstate_index_in;
      assert_always(lowerionstage >= 0);
      assert_always(lowerlevel >= 0);

      // store only photoionisation crosssections for ions that are part of the current model atom
      if (lowerion >= 0 && upperion < get_nions(element) && lowerlevel < get_nlevels_ionising(element, lowerion)) {
        read_phixs_data_table(phixsfile, nphixspoints_inputtable, element, lowerion, lowerlevel, upperion,
                              upperlevel_in, tmpallphixs, tmpallphixstargets, &mem_usage_phixs, phixs_file_version);

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
          // one day we might want to put all of the cross section points onto a single line,
          // so don't use getline here
          float phixs = 0;
          assert_always(phixsfile >> phixs);
        }
      }
    }
  }

  printlnlog("[info] mem_usage: photoionisation tables occupy {:.3f} MB", mem_usage_phixs / 1024. / 1024.);
}

constexpr auto downtranslevelstart(const int level) {
  // each level index is associated with a block of size levelindex spanning all possible down transitions.
  // so use the formula for the sum of 1 + 2 + 3 + 4 + ... + level
  return level * (level + 1) / 2;
}

void read_ion_levels(std::istream& adata, const int element, const int ion, const int nions, const int nlevels,
                     int nlevelsmax, const double energyoffset, const double ionpot,
                     std::vector<TempEnergyLevel>& temp_alllevels) {
  std::string line;
  static std::istringstream ssline;
  for (int level = 0; level < nlevels; level++) {
    int levelindex_in = 0;
    double levelenergy{NAN};
    float statweight{NAN};
    int ntransitions = 0;
    assert_always(get_noncommentline(adata, line));
    ssline.clear();
    ssline.str(line);
    assert_always(ssline >> levelindex_in >> levelenergy >> statweight >> ntransitions);
    assert_always(levelindex_in == level + groundstate_index_in);

    if (level < nlevelsmax) {
      assert_always(statweight > 0.);
      const double currentlevelenergy = (energyoffset + levelenergy) * EV;
      assert_always(std::ssize(temp_alllevels) == get_ionuniquelevelindexstart(element, ion) + level);
      temp_alllevels.push_back({
          .epsilon = currentlevelenergy,
          .ndowntrans = 0,
          .nuptrans = 0,
          .stat_weight = statweight,
      });

      // The level contributes to the ionisinglevels if its energy
      // is below the ionisation potential and the level doesn't
      // belong to the topmost ion included.
      // Rate coefficients are only available for ionising levels.
      if (levelenergy < ionpot && ion < nions - 1) {
        globals::elements[element].ions[ion].nlevels_ionising++;
      }
    }
  }
}

void read_ion_transitions(std::istream& ftransitiondata, const int ion_transition_count_in_file,
                          std::vector<IonTransitionsInput>& iontransitiontable, const int nlevels_requiretransitions,
                          const int nlevels_requiretransitions_upperlevels) {
  iontransitiontable.clear();
  iontransitiontable.reserve(ion_transition_count_in_file);

  std::string line;
  static std::istringstream ssline;

  // will be autodetected from first table row. old format had an index column and no collstr or forbidden columns
  bool oldtransitionformat = false;

  int prev_upper = -1;
  int prev_lower = 0;
  for (int i = 0; i < ion_transition_count_in_file; i++) {
    int lower_in = -1;
    int upper_in = -1;
    float A = 0;
    float coll_str = -1.;
    int intforbidden = 0;
    assert_always(getline(ftransitiondata, line));
    if (i == 0) {
      ssline.clear();
      ssline.str(line);
      std::string word;
      int column_count = 0;
      while (ssline >> word) {
        column_count++;
      }
      assert_always(column_count == 4 || column_count == 5);
      oldtransitionformat = (column_count == 4);
    }
    ssline.clear();
    ssline.str(line);
    if (!oldtransitionformat) {
      assert_always(ssline >> lower_in >> upper_in >> A >> coll_str >> intforbidden);
    } else {
      int transindex = 0;  // not used
      assert_always(ssline >> transindex >> lower_in >> upper_in >> A);
    }
    const int lower = lower_in - groundstate_index_in;
    const int upper = upper_in - groundstate_index_in;
    assert_always(lower >= 0);
    assert_always(upper > lower);
    assert_always(lower >= prev_lower);
    assert_always(upper >= prev_upper || lower > prev_lower);

    // this entire block can be removed if we don't want to add in extra collisonal
    // transitions between levels
    if (prev_lower < nlevels_requiretransitions) {
      assert_always(prev_lower >= 0);
      int stoplevel = 0;
      if (lower == prev_lower && upper > prev_upper + 1) {
        // same lower level, but some upper levels were skipped over
        stoplevel = upper - 1;
        if (stoplevel >= nlevels_requiretransitions_upperlevels) {
          stoplevel = nlevels_requiretransitions_upperlevels - 1;
        }
      } else if ((lower > prev_lower) && prev_upper < (nlevels_requiretransitions_upperlevels - 1)) {
        // we've moved onto another lower level, but the previous one was missing some required transitions
        stoplevel = nlevels_requiretransitions_upperlevels - 1;
      } else {
        stoplevel = -1;
      }

      for (int tmplevel = prev_upper + 1; tmplevel <= stoplevel; tmplevel++) {
        if (tmplevel == prev_lower) {
          continue;
        }
        assert_always(tmplevel >= 0);
        assert_always(tmplevel > prev_lower);
        iontransitiontable.push_back(
            {.lower = prev_lower, .upper = tmplevel, .A = 0., .coll_str = -2., .forbidden = true});
      }
    }

    iontransitiontable.push_back(
        {.lower = lower, .upper = upper, .A = A, .coll_str = coll_str, .forbidden = (intforbidden == 1)});
    prev_lower = lower;
    prev_upper = upper;
  }
}

void add_transitions_to_unsorted_linelist(const int element, const int ion,
                                          const std::vector<IonTransitionsInput>& iontransitiontable,
                                          std::vector<TempLineTransitionInput>& temp_linelist,
                                          std::vector<TempAllTransInput>& temp_alltranslist,
                                          std::vector<TempEnergyLevel>& temp_alllevels) {
  const auto nlevels = get_nlevels(element, ion);
  auto ion_levels = std::span{temp_alllevels}.subspan(get_ionuniquelevelindexstart(element, ion), nlevels);
  static std::vector<int> iondowntranstmplineindices;
  iondowntranstmplineindices.resize(downtranslevelstart(nlevels));

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

    std::ranges::fill(iondowntranstmplineindices, -99);

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

      // Make sure that we don't allow duplicate. In that case take only the lines
      // first occurrence
      int& downtranslineindex = iondowntranstmplineindices[downtranslevelstart(level) + lowerlevel];

      // negative means that the transition hasn't been seen yet
      if (downtranslineindex < 0) {
        downtranslineindex = globals::nlines++;

        const int nupperdowntrans = ion_levels[level].ndowntrans + 1;
        ion_levels[level].ndowntrans = nupperdowntrans;

        const int nloweruptrans = ion_levels[lowerlevel].nuptrans + 1;
        ion_levels[lowerlevel].nuptrans = nloweruptrans;

        ion_updowntranscount += 2;

        if (pass == 1) {
          const auto g_ratio = static_cast<double>(ion_levels[level].stat_weight) / ion_levels[lowerlevel].stat_weight;
          const auto f_ul =
              static_cast<float>(g_ratio * ME * pow(CLIGHT, 3) / (8 * pow(QE * nu_trans * PI, 2)) * transition.A);
          assert_always(std::isfinite(f_ul));

          temp_linelist.push_back({
              .nu = nu_trans,
              .einstein_A = transition.A,
              .elementindex = element,
              .ionindex = ion,
              .upperlevelindex = level,
              .lowerlevelindex = lowerlevel,
          });

          // the line list has not been sorted yet, so the store the level index for now and
          // the index into the sorted line list will be set later

          temp_alltranslist[ion_levels[level].alltrans_startdown + nupperdowntrans - 1] = {
              .lineindex = -1,
              .targetlevelindex = lowerlevel,
              .einstein_A = transition.A,
              .coll_str = transition.coll_str,
              .osc_strength = f_ul,
              .forbidden = transition.forbidden};
          const auto lowerstartup = ion_levels[lowerlevel].alltrans_startup();
          temp_alltranslist[lowerstartup + nloweruptrans - 1] = {.lineindex = -1,
                                                                 .targetlevelindex = level,
                                                                 .einstein_A = transition.A,
                                                                 .coll_str = transition.coll_str,
                                                                 .osc_strength = f_ul,
                                                                 .forbidden = transition.forbidden};
        }

      } else if (pass == 1) {
        // This is a new branch to deal with lines that have different types of transition. It should trip after a
        // transition is already known.

        if ((temp_linelist[downtranslineindex].elementindex != element) ||
            (temp_linelist[downtranslineindex].ionindex != ion) ||
            (temp_linelist[downtranslineindex].upperlevelindex != level) ||
            (temp_linelist[downtranslineindex].lowerlevelindex != lowerlevel)) {
          printlnlog("[input] Failure to identify level pair for duplicate bb-transition ... going to abort now");
          printlnlog("[input]   element {} ion {} targetlevel {} level {}", element, ion, lowerlevel, level);
          printlnlog("[input]   transitions[level].to[targetlevel]=lineindex {}", downtranslineindex);
          printlnlog("[input]   A_ul {:g}, coll_str {:g}", transition.A, transition.coll_str);
          printlnlog(
              "[input]   globals::linelist[lineindex].elementindex {}, globals::linelist[lineindex].ionindex {}, "
              "globals::linelist[lineindex].upperlevelindex {}, globals::linelist[lineindex].lowerlevelindex {}",
              temp_linelist[downtranslineindex].elementindex, temp_linelist[downtranslineindex].ionindex,
              temp_linelist[downtranslineindex].upperlevelindex, temp_linelist[downtranslineindex].lowerlevelindex);
          std::abort();
        }

        const auto g_ratio = static_cast<double>(ion_levels[level].stat_weight) / ion_levels[lowerlevel].stat_weight;
        const auto f_ul =
            static_cast<float>(g_ratio * ME * pow(CLIGHT, 3) / (8 * pow(QE * nu_trans * PI, 2)) * transition.A);

        auto& downtransition =
            temp_alltranslist[ion_levels[level].alltrans_startdown + ion_levels[level].ndowntrans - 1];

        assert_always(downtransition.targetlevelindex == lowerlevel);

        downtransition.einstein_A += transition.A;
        downtransition.osc_strength += f_ul;
        downtransition.coll_str = std::max(downtransition.coll_str, transition.coll_str);

        const auto lowerstartup = ion_levels[lowerlevel].alltrans_startup();
        auto& uptransition = temp_alltranslist[lowerstartup + ion_levels[lowerlevel].nuptrans - 1];

        // as above, the downtrans list should be searched to find the correct index instead of using the last one.
        // assert_always(uptransition.targetlevelindex == level);

        uptransition.einstein_A += transition.A;
        uptransition.osc_strength += f_ul;
        uptransition.coll_str = std::max(uptransition.coll_str, transition.coll_str);
      }
    }
  }
}

auto calculate_nlevels_groundterm(const int element, const int ion) -> int {
  const int nlevels = get_nlevels(element, ion);
  if (nlevels == 1) {
    return 1;
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

auto search_groundphixslist(const double nu_edge, const int element_in, const int ion_in, const int level_in) -> int
// Return the closest ground level continuum index to the given edge
// frequency. If the given edge frequency is redder than the reddest
// continuum return -1.
// NB: groundphixslist must be in ascending order.
{
  assert_always((USE_LUT_PHOTOION || USE_LUT_BFHEATING));

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
        "[fatal] search_groundphixslist: element {}, ion {}, level {} has edge_frequency {:g} equal to the bluest "
        "ground-level continuum",
        element_in, ion_in, level_in, nu_edge);
    printlnlog("[fatal] search_groundphixslist: bluest ground level continuum is element {}, ion {} at nu_edge {:g}",
               element, ion, globals::groundcont_nu_edge[i - 1]);
    printlnlog("[fatal] search_groundphixslist: i {}, nbfcontinua_ground {}", i, globals::nbfcontinua_ground);
    printlnlog(
        "[fatal] This shouldn't happen, is hoewever possible if there are multiple levels in the adata file at "
        "energy=0");
    for (int looplevels = 0; looplevels < get_nlevels(element_in, ion_in); looplevels++) {
      printlnlog("[fatal]   element {}, ion {}, level {}, energy {:g}", element_in, ion_in, looplevels,
                 epsilon(element_in, ion_in, looplevels));
    }
    printlnlog("[fatal] Abort omitted ... MAKE SURE ATOMIC DATA ARE CONSISTENT");
    return i - 1;
    // abort();
  }

  const double left_diff = nu_edge - globals::groundcont_nu_edge[i - 1];
  const double right_diff = globals::groundcont_nu_edge[i] - nu_edge;
  return (left_diff <= right_diff) ? i - 1 : i;
}

// set up the photoionisation transition lists
// and temporary gamma/kappa lists for each thread
void setup_phixs_list() {
  printlnlog("[info] read_atomicdata: number of bfcontinua {}", globals::nbfcontinua);
  printlnlog("[info] read_atomicdata: number of ground-level bfcontinua {}", globals::nbfcontinua_ground);

  struct TempGroundPhotoion {
    double nu_edge;
    int element;
    int ion;
  };

  struct TempPhotoionTransitionInput {
    double nu_edge;
    int element;
    int ion;
    int level;
    int phixstargetindex;
    int upperlevel;
    int uniquelevelindex;
    double probability;
    int index_in_groundphixslist;
    int bfestimindex;
  };

  auto groundcont_nu_edge = MPI_shared_malloc_span<double>(globals::nbfcontinua_ground);
  auto groundcont_element = MPI_shared_malloc_span<int>(globals::nbfcontinua_ground);
  auto groundcont_ion = MPI_shared_malloc_span<int>(globals::nbfcontinua_ground);

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
  MPI_Barrier(globals::mpi_comm_node);
  globals::groundcont_nu_edge = groundcont_nu_edge;
  globals::groundcont_element = groundcont_element;
  globals::groundcont_ion = groundcont_ion;

  auto allcont = MPI_shared_malloc_span<TempPhotoionTransitionInput>(globals::nbfcontinua);
  printlnlog("[info] mem_usage: photoionisation list occupies {:.3f} MB",
             globals::nbfcontinua * (sizeof(TempPhotoionTransitionInput)) / 1024. / 1024.);
  const auto groundcontindices = std::ranges::iota_view{0, globals::nbfcontinua_ground};
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

        for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
          assert_always(allcontindex < std::ssize(allcont));

          int index_in_groundphixslist = -1;
          if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
            const double nu_edge_target0 = get_phixs_threshold(uniquelevelindex, 0) / H;
            index_in_groundphixslist = search_groundphixslist(nu_edge_target0, element, ion, level);

            globals::alllevels.closestgroundlevelcont[uniquelevelindex] = index_in_groundphixslist;
          }

          allcont[allcontindex] = {
              .nu_edge = get_phixs_threshold(uniquelevelindex, phixstargetindex) / H,
              .element = element,
              .ion = ion,
              .level = level,
              .phixstargetindex = phixstargetindex,
              .upperlevel = get_phixsupperlevel(uniquelevelindex, phixstargetindex),
              .uniquelevelindex = uniquelevelindex,
              .probability = get_phixsprobability(uniquelevelindex, phixstargetindex),
              .index_in_groundphixslist = index_in_groundphixslist,
              .bfestimindex = -1,
          };

          allcontindex++;
        }
      }
    }
  }

  assert_always(allcontindex == globals::nbfcontinua);
  assert_always(globals::nbfcontinua >= 0);  // was initialised as -1 before startup
  // just so that clang-tidy doesn't throw errors on the assumption that nbfcontinua is changing
  const auto nbfcontinua = globals::nbfcontinua;

  std::vector<double> temp_bfestim_nu_edge;
  if (nbfcontinua > 0) {
    // indices above were temporary only. continuum index should be to the sorted list
    MPI_Barrier(globals::mpi_comm_node);
    if (globals::rank_in_node == 0) {
      std::ranges::SORT_OR_STABLE_SORT(allcont, std::ranges::less{}, &TempPhotoionTransitionInput::nu_edge);
    }
    MPI_Barrier(globals::mpi_comm_node);

    for (int i = 0; i < nbfcontinua; i++) {
      if (DETAILED_BF_ESTIMATORS_ON &&
          LEVEL_HAS_BFEST(get_atomicnumber(allcont[i].element), get_ionstage(allcont[i].element, allcont[i].ion),
                          allcont[i].level)) {
        allcont[i].bfestimindex = static_cast<int>(temp_bfestim_nu_edge.size());
        temp_bfestim_nu_edge.push_back(allcont[i].nu_edge);
      } else {
        allcont[i].bfestimindex = -1;
      }
    }

    auto bfestim_nu_edge = MPI_shared_malloc_span<double>(std::ssize(temp_bfestim_nu_edge));
    auto allcont_nu_edge = MPI_shared_malloc_span<double>(nbfcontinua);
    auto allcont_element = MPI_shared_malloc_span<int>(nbfcontinua);
    auto allcont_ion = MPI_shared_malloc_span<int>(nbfcontinua);
    auto allcont_level = MPI_shared_malloc_span<int>(nbfcontinua);
    auto allcont_phixstargetindex = MPI_shared_malloc_span<int>(nbfcontinua);
    auto allcont_upperlevel = MPI_shared_malloc_span<int>(nbfcontinua);
    auto allcont_uniquelevelindex = MPI_shared_malloc_span<int>(nbfcontinua);
    auto allcont_probability = MPI_shared_malloc_span<double>(nbfcontinua);
    auto allcont_index_in_groundphixslist = MPI_shared_malloc_span<int>(nbfcontinua);
    auto allcont_bfestimindex = MPI_shared_malloc_span<int>(nbfcontinua);
    if (globals::rank_in_node == 0) {
      for (int i = 0; i < std::ssize(temp_bfestim_nu_edge); i++) {
        bfestim_nu_edge[i] = temp_bfestim_nu_edge[i];
      }
      for (int i = 0; i < nbfcontinua; i++) {
        allcont_nu_edge[i] = allcont[i].nu_edge;
        allcont_element[i] = allcont[i].element;
        allcont_ion[i] = allcont[i].ion;
        allcont_level[i] = allcont[i].level;
        allcont_phixstargetindex[i] = allcont[i].phixstargetindex;
        allcont_upperlevel[i] = allcont[i].upperlevel;
        allcont_uniquelevelindex[i] = allcont[i].uniquelevelindex;
        allcont_probability[i] = allcont[i].probability;
        allcont_index_in_groundphixslist[i] = allcont[i].index_in_groundphixslist;
        allcont_bfestimindex[i] = allcont[i].bfestimindex;
      }
    }
    MPI_Barrier(globals::mpi_comm_node);
    globals::bfestim_nu_edge = bfestim_nu_edge;
    globals::allcont_nu_edge = allcont_nu_edge;
    globals::allcont_element = allcont_element;
    globals::allcont_ion = allcont_ion;
    globals::allcont_level = allcont_level;
    globals::allcont_phixstargetindex = allcont_phixstargetindex;
    globals::allcont_upperlevel = allcont_upperlevel;
    globals::allcont_uniquelevelindex = allcont_uniquelevelindex;
    globals::allcont_probability = allcont_probability;
    globals::allcont_index_in_groundphixslist = allcont_index_in_groundphixslist;
    globals::allcont_bfestimindex = allcont_bfestimindex;

    setup_photoion_luts();
  }
  printlnlog("[info] bound-free estimators track bfestimcount {} photoionisation transitions",
             globals::bfestim_nu_edge.size());
}

void read_autoion_data() {
  // read in autoionisation rate data
  if (!std::filesystem::exists("autoion.txt")) {
    return;
  }

  std::vector<globals::LevelAutoion> temp_allautoion;

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

  while (get_noncommentline(autoionfile, autoionline)) {
    assert_always(std::istringstream(autoionline) >> Z >> upperionstage >> upperlevel_in >> lowerionstage >>
                  lowerlevel_in >> autoion_A);

    assert_always(Z > 0);
    assert_always(upperionstage >= 2);
    assert_always(lowerionstage >= 1);

    const int element = get_elementindex(Z);

    if (element >= 0 && get_nions(element) > 0) {
      // translate readin ionstages to ion indices

      const int upperion = upperionstage - get_ionstage(element, 0);
      const int lowerion = lowerionstage - get_ionstage(element, 0);
      const int lowerlevel = lowerlevel_in - groundstate_index_in;
      const int upperlevel = upperlevel_in - groundstate_index_in;

      assert_always(upperion >= 0 && upperion < get_nions(element));
      assert_always(lowerion >= 0 && lowerion < get_nions(element));
      assert_always(lowerlevel >= 0 && lowerlevel < get_nlevels(element, lowerion));
      assert_always(upperlevel >= 0 && upperlevel < get_nlevels(element, upperion));
      assert_always(upperion > lowerion);
      const bool level_is_nlte = is_nlte(element, lowerion, lowerlevel - 1);
      allautoion_levels_are_not_nlte = allautoion_levels_are_not_nlte && !level_is_nlte;
      allautoion_levels_are_nlte = allautoion_levels_are_nlte && level_is_nlte;

      // store only for ions that are part of the current model atom
      if (lowerion >= 0 && upperion < get_nions(element)) {
        printlnlog("Got to noting data for Z {} upperion {} upperlvl {} lowerion {} lowerlvl {} with A {:g}", Z,
                   upperion, upperlevel, lowerion, lowerlevel, autoion_A);

        const int nautoiondowntrans = get_nautoiondowntrans(element, lowerion, lowerlevel) + 1;
        set_nautoiondowntrans(element, lowerion, lowerlevel, nautoiondowntrans);

        const int nautoionuptrans = get_nautoionuptrans(element, upperion, upperlevel) + 1;
        set_nautoionuptrans(element, upperion, upperlevel, nautoionuptrans);

        if (globals::alllevels.allautoion_start[get_uniquelevelindex(element, lowerion, lowerlevel)] < 0) {
          assert_always(nautoiondowntrans == 1);
          //  this is the first autoionizing transition for this level, so set the start index
          globals::alllevels.allautoion_start[get_uniquelevelindex(element, lowerion, lowerlevel)] =
              static_cast<int>(temp_allautoion.size());
        }

        temp_allautoion.push_back({.autoion_A = static_cast<float>(autoion_A),
                                   .elementindex = element,
                                   .lowerionindex = lowerion,
                                   .lowerlevelindex = lowerlevel,
                                   .upperionindex = upperion,
                                   .upperlevelindex = upperlevel});
      }
    }
  }

  assert_always(allautoion_levels_are_not_nlte || allautoion_levels_are_nlte);

  globals::allautoion = MPI_shared_malloc_span<globals::LevelAutoion>(temp_allautoion.size());
  if (globals::rank_in_node == 0) {
    std::copy_n(temp_allautoion.cbegin(), temp_allautoion.size(), globals::allautoion.data());
  }
  MPI_Barrier(globals::mpi_comm_node);

  // Plan is that autoionizing levels will be explicitly included in the NLTE population solver, but that their level
  // populations do not need to be accurately known - so if the ion has a superlevel already, then we will try to attach
  // the autoionizing level populations to that for all purposes outside the NLTE solber. For this, the ions need to
  // know how many autoionizing levels they have. So count those up now.

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
  assert_always(nlevels_autoion_sum == std::ssize(temp_allautoion));
}

void read_phixs_data() {
  globals::nbfcontinua_ground = 0;
  globals::nbfcontinua = 0;
  std::vector<float> tmpallphixs;
  std::vector<PhotoionTarget> tmpallphixstargets;

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

  if (phixs_file_version_exists[2]) {
    read_phixs_file(2, tmpallphixs, tmpallphixstargets);
    MPI_Barrier(globals::mpi_comm_node);
  }

  if (phixs_file_version_exists[1]) {
    read_phixs_file(1, tmpallphixs, tmpallphixstargets);
    MPI_Barrier(globals::mpi_comm_node);
  }

  int cont_index = 0;
  ptrdiff_t nbftables = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        const int nphixstargets = get_nphixstargets(element, ion, level);
        const auto uniquelevelindex = get_uniquelevelindex(element, ion, level);
        globals::alllevels.bflist_start[uniquelevelindex] = (nphixstargets > 0) ? cont_index : -1;
        cont_index += nphixstargets;
        if (nphixstargets > 0) {
          nbftables++;
        }
      }
    }
  }
  printlnlog("cont_index {}", cont_index);
  assert_always(cont_index == std::ssize(tmpallphixstargets));

  if (!tmpallphixs.empty()) {
    assert_always((nbftables * globals::NPHIXSPOINTS) == std::ssize(tmpallphixs));

    // copy the photoionisation tables into one contiguous block of memory
    globals::allphixs = MPI_shared_malloc_span<float>(std::ssize(tmpallphixs));
    auto allphixstargets_levelindex = MPI_shared_malloc_span<int>(std::ssize(tmpallphixstargets));
    auto allphixstargets_probability = MPI_shared_malloc_span<double>(std::ssize(tmpallphixstargets));

    if (globals::rank_in_node == 0) {
      std::copy_n(tmpallphixs.cbegin(), tmpallphixs.size(), globals::allphixs.data());

      for (int i = 0; i < std::ssize(tmpallphixstargets); i++) {
        allphixstargets_levelindex[i] = tmpallphixstargets[i].levelindex;
        allphixstargets_probability[i] = tmpallphixstargets[i].probability;
      }
    }

    MPI_Barrier(globals::mpi_comm_node);
    globals::allphixstargets_levelindex = allphixstargets_levelindex;
    globals::allphixstargets_probability = allphixstargets_probability;

    tmpallphixs.clear();
    tmpallphixs.shrink_to_fit();

    tmpallphixstargets.clear();
    tmpallphixstargets.shrink_to_fit();
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
            printlnlog("WARNING: Z={} ionstage {} nlevels_groundterm {} phixstargetlevels(ion-1) {}.",
                       get_atomicnumber(element), get_ionstage(element, ion), nlevels_groundterm, phixstargetlevels);
          }
        }
      }
    }
  }

  setup_phixs_list();
}

auto read_compositiondata() -> std::vector<int> {
  auto compositiondata = fstream_required("compositiondata.txt", std::ios::in);

  int nelements_in = 0;
  assert_always(compositiondata >> nelements_in);
  globals::elements.resize(nelements_in);

  // temperature to determine relevant ionstages
  int T_preset = 0;
  assert_always(compositiondata >> T_preset);
  assert_always(T_preset == 0);  // no longer in use
  int homogeneous_abundances = 0;
  assert_always(compositiondata >> homogeneous_abundances);
  assert_always(homogeneous_abundances == 0);  // no longer in use

  auto nlevelsmax_readin = std::vector<int>(get_nelements());
  auto nions_readin = std::vector<int>(get_nelements());
  int uniqueionindex = 0;  // index into list of all ions of all elements
  for (int element = 0; element < get_nelements(); element++) {
    // read information about the next element which should be stored to memory
    int atomicnumber = 0;
    int lowermost_ionstage = 0;
    int uppermost_ionstage = 0;
    double uniformabundance{NAN};  // no longer in use mode for setting uniform abundances
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

  globals::allions = MPI_shared_malloc_span<Ion>(uniqueionindex);

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
  int nbfcheck = 0;
  for (int element = 0; element < get_nelements(); element++) {
    // now read in data for all ions of the current element. before doing so initialize
    // energy scale for the current element (all level energies are stored relative to
    // the ground level of the neutral ion)
    const auto atomicnumber = globals::elements[element].anumber;
    double energyoffset = 0.;
    double ionpot = 0.;
    const auto nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int ionstage = globals::elements[element].lowest_ionstage + ion;
      // calculate the current levels ground level energy
      assert_always(ionpot >= 0);
      energyoffset += ionpot;

      // read information for the elements next ionstage
      int adata_Z_in = -1;
      int adata_ionstage_in = -1;
      int nlevels_in_file = 0;

      while (adata_Z_in != atomicnumber || adata_ionstage_in != ionstage) {
        // skip over this ion block
        if (adata_Z_in == atomicnumber) {
          printlnlog("increasing energyoffset by ionpot {:g}", ionpot);
          energyoffset += ionpot;
        }
        for (int i = 0; i < nlevels_in_file; i++) {
          std::getline(adata, line);
        }

        assert_always(get_noncommentline(adata, line));
        ssline.clear();
        ssline.str(line);
        assert_always(ssline >> adata_Z_in >> adata_ionstage_in >> nlevels_in_file >> ionpot);
      }

      printlnlog("adata.txt: Z {} ionstage {} nlevels_in_file {}", adata_Z_in, adata_ionstage_in, nlevels_in_file);

      // optionally limit the top ion to one level and no transitions
      const int nlevelslimit = (SINGLE_LEVEL_TOP_ION && ion == (nions - 1)) ? 1 : nlevelsmax_readin[element];
      const int nlevelskept = (nlevelslimit < 0) ? nlevels_in_file : std::min(nlevelslimit, nlevels_in_file);

      if (nlevels_in_file > nlevelskept) {
        printlnlog("[info] read_atomicdata: reduce number of levels from {} to {} for Z {:2} ionstage {}",
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
          .ionpot = ionpot * EV,
      };

      assert_always(std::ssize(temp_alllevels) == uniquelevelindex);

      read_ion_levels(adata, element, ion, nions, nlevels_in_file, nlevelskept, energyoffset, ionpot, temp_alllevels);
      uniquelevelindex += get_nlevels(element, ion);

      if (ion < nions - 1) {
        nbfcheck += get_nlevels_ionising(element, ion);
      }
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

      // first nlevels_requiretransitions levels will be collisionally
      // coupled to the first nlevels_requiretransitions_upperlevels levels (assumed forbidden)
      // use 0 to disable adding extra transitions

      const int nlevels_requiretransitions =
          std::min(nlevelskept, NLEVELS_REQUIRETRANSITIONS(atomicnumber, adata_ionstage_in));
      // next value with have no effect if nlevels_requiretransitions = 0
      const int nlevels_requiretransitions_upperlevels = nlevelskept;

      // load transition table for the current ion to temporary memory
      if (nlevelskept <= 1) {
        // we will not read in any transitions, just skip past these lines in the file
        for (int i = 0; i < ion_transition_count_in_file; i++) {
          assert_always(getline(ftransitiondata, line));
        }
      } else {
        read_ion_transitions(ftransitiondata, ion_transition_count_in_file, iontransitiontable,
                             nlevels_requiretransitions, nlevels_requiretransitions_upperlevels);
        add_transitions_to_unsorted_linelist(element, ion, iontransitiontable, temp_linelist, temp_alltranslist,
                                             temp_alllevels);
      }
    }
  }
  printlnlog("nbfcheck {}", nbfcheck);
}

void read_atomicdata_files() {
  auto nlevelsmax_readin = read_compositiondata();

  printlnlog("SINGLE_LEVEL_TOP_ION: {}", SINGLE_LEVEL_TOP_ION ? "true" : "false");

  std::vector<TempEnergyLevel> temp_alllevels;
  std::vector<TempLineTransitionInput> temp_linelist;
  std::vector<TempAllTransInput> temp_alltranslist;

  if (globals::rank_in_node == 0) {
    temp_linelist.reserve(1 << 22);  // reserve initial space for 4 million lines to avoid too many reallocations
    temp_alltranslist.reserve(1 << 22);
    read_levels_and_transitions(temp_alllevels, temp_linelist, temp_alltranslist, nlevelsmax_readin);
  }
  MPI_Barrier(globals::mpi_comm_node);

  update_includedionslevels_maxnions();

  MPI_Bcast_safe(globals::nlines, 0, globals::mpi_comm_node);
  printlnlog("nlines {}", globals::nlines);

  if (globals::rank_in_node == 0) {
    assert_always(globals::nlines == std::ssize(temp_linelist));
    temp_linelist.shrink_to_fit();

    // sort the lineline in descending frequency
    std::SORT_OR_STABLE_SORT(
        EXEC_PAR_UNSEQ temp_linelist.begin(), temp_linelist.end(), [](const auto& a, const auto& b) {
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
          printlnlog("a: Z={} ionstage {} lower {} upper {} nu {:g} lambda {:g}", get_atomicnumber(a1.elementindex),
                     get_ionstage(a1.elementindex, a1.ionindex), a1.lowerlevelindex, a1.upperlevelindex, a1.nu,
                     1e8 * CLIGHT / a1.nu);
          printlnlog("b: Z={} ionstage {} lower {} upper {} nu {:g} lambda {:g}", get_atomicnumber(a2.elementindex),
                     get_ionstage(a2.elementindex, a2.ionindex), a2.lowerlevelindex, a2.upperlevelindex, a2.nu,
                     1e8 * CLIGHT / a2.nu);
        }
      }
    }
  }

  // create a shared level list and copy data across, freeing the local copy
  ptrdiff_t nlevels = std::ssize(temp_alllevels);
  MPI_Bcast_safe(nlevels, 0, globals::mpi_comm_node);

  auto alllevels_alltrans_startdown = MPI_shared_malloc_span<int>(nlevels);
  auto alllevels_ndowntrans = MPI_shared_malloc_span<int>(nlevels);
  auto alllevels_nuptrans = MPI_shared_malloc_span<int>(nlevels);
  auto alllevels_epsilon = MPI_shared_malloc_span<double>(nlevels);
  auto alllevels_statweight = MPI_shared_malloc_span<float>(nlevels);
  auto alllevels_matransblock_start = MPI_shared_malloc_span<int>(nlevels);
  globals::alllevels.allautoion_start = MPI_shared_malloc_span<int>(nlevels, -1);
  globals::alllevels.nautoiondowntrans = MPI_shared_malloc_span<int>(nlevels, 0);
  globals::alllevels.nautoionuptrans = MPI_shared_malloc_span<int>(nlevels, 0);
  globals::alllevels.closestgroundlevelcont = MPI_shared_malloc_span<int>(nlevels, -1);
  globals::alllevels.phixsstart = MPI_shared_malloc_span<int>(nlevels, -1);
  globals::alllevels.nphixstargets = MPI_shared_malloc_span<int>(nlevels, 0);
  globals::alllevels.phixstargetstart = MPI_shared_malloc_span<int>(nlevels, -1);
  globals::alllevels.bflist_start = MPI_shared_malloc_span<int>(nlevels, -1);
  if (globals::rank_in_node == 0) {
    int chtransindex = 0;
    for (auto i = 0zU; i < temp_alllevels.size(); i++) {
      alllevels_alltrans_startdown[i] = temp_alllevels[i].alltrans_startdown;
      alllevels_ndowntrans[i] = temp_alllevels[i].ndowntrans;
      alllevels_nuptrans[i] = temp_alllevels[i].nuptrans;
      alllevels_epsilon[i] = temp_alllevels[i].epsilon;
      alllevels_statweight[i] = temp_alllevels[i].stat_weight;
      alllevels_matransblock_start[i] = chtransindex;
      chtransindex += ((2 * alllevels_ndowntrans[i]) + alllevels_nuptrans[i]);
    }
  }
  MPI_Barrier(globals::mpi_comm_node);
  globals::alllevels.alltrans_startdown = alllevels_alltrans_startdown;
  globals::alllevels.ndowntrans = alllevels_ndowntrans;
  globals::alllevels.nuptrans = alllevels_nuptrans;
  globals::alllevels.epsilon = alllevels_epsilon;
  globals::alllevels.statweight = alllevels_statweight;
  globals::alllevels.matransblock_start = alllevels_matransblock_start;
  temp_alllevels.clear();
  temp_alllevels.shrink_to_fit();

  const int updowntranscount = []() -> int {
    const int downtranscount = std::ranges::fold_left(globals::alllevels.ndowntrans, 0, std::plus<>{});
    const int uptranscount = std::ranges::fold_left(globals::alllevels.nuptrans, 0, std::plus<>{});

    printlnlog("total uptrans {}", uptranscount);
    printlnlog("total downtrans {}", downtranscount);
    return downtranscount + uptranscount;
  }();
  printlnlog("[info] mem_usage: transition lists occupy {:.3f} MB (node shared memory)",
             updowntranscount * (2 * sizeof(int) + 3 * sizeof(float) + sizeof(bool)) / 1024. / 1024.);

  // create a shared all transitions list and then copy data across, freeing the local copy
  MPI_Barrier(globals::mpi_comm_node);

  auto alltrans_lineindex = MPI_shared_malloc_span<int>(updowntranscount);
  auto alltrans_targetlevelindex = MPI_shared_malloc_span<int>(updowntranscount);
  auto alltrans_einstein_A = MPI_shared_malloc_span<float>(updowntranscount);
  auto alltrans_coll_str = MPI_shared_malloc_span<float>(updowntranscount);
  auto alltrans_osc_strength = MPI_shared_malloc_span<float>(updowntranscount);
  auto alltrans_forbidden = MPI_shared_malloc_span<bool>(updowntranscount);

  if (globals::rank_in_node == 0) {
    assert_always(std::ssize(temp_alltranslist) == updowntranscount);
    for (int t = 0; t < updowntranscount; t++) {
      alltrans_lineindex[t] = temp_alltranslist[t].lineindex;
      alltrans_targetlevelindex[t] = temp_alltranslist[t].targetlevelindex;
      alltrans_einstein_A[t] = temp_alltranslist[t].einstein_A;
      alltrans_coll_str[t] = temp_alltranslist[t].coll_str;
      alltrans_osc_strength[t] = temp_alltranslist[t].osc_strength;
      alltrans_forbidden[t] = temp_alltranslist[t].forbidden;
    }
  }
  temp_alltranslist.clear();
  temp_alltranslist.shrink_to_fit();

  globals::alltrans.targetlevelindex = alltrans_targetlevelindex;
  globals::alltrans.einstein_A = alltrans_einstein_A;
  globals::alltrans.coll_str = alltrans_coll_str;
  globals::alltrans.osc_strength = alltrans_osc_strength;
  globals::alltrans.forbidden = alltrans_forbidden;

  // create a linelist shared on node and then copy data across, freeing the local copy

  auto linelist_nu = MPI_shared_malloc_span<double>(globals::nlines);
  auto linelist_einstein_A = MPI_shared_malloc_span<float>(globals::nlines);
  auto linelist_elementindex = MPI_shared_malloc_span<int>(globals::nlines);
  auto linelist_ionindex = MPI_shared_malloc_span<int>(globals::nlines);
  auto linelist_upperlevelindex = MPI_shared_malloc_span<int>(globals::nlines);
  auto linelist_lowerlevelindex = MPI_shared_malloc_span<int>(globals::nlines);

  if (globals::rank_in_node == 0) {
    assert_always(std::ssize(temp_linelist) == globals::nlines);
    for (int t = 0; t < globals::nlines; t++) {
      linelist_nu[t] = temp_linelist[t].nu;
      linelist_einstein_A[t] = temp_linelist[t].einstein_A;
      linelist_elementindex[t] = temp_linelist[t].elementindex;
      linelist_ionindex[t] = temp_linelist[t].ionindex;
      linelist_upperlevelindex[t] = temp_linelist[t].upperlevelindex;
      linelist_lowerlevelindex[t] = temp_linelist[t].lowerlevelindex;
    }
  }
  temp_linelist.clear();
  temp_linelist.shrink_to_fit();
  MPI_Barrier(globals::mpi_comm_node);

  globals::linelist.nu = linelist_nu;
  globals::linelist.einstein_A = linelist_einstein_A;
  globals::linelist.elementindex = linelist_elementindex;
  globals::linelist.ionindex = linelist_ionindex;
  globals::linelist.upperlevelindex = linelist_upperlevelindex;
  globals::linelist.lowerlevelindex = linelist_lowerlevelindex;

  const double linelist_mem_MB =
      (globals::nlines * sizeof(double)  // nu
       + globals::nlines * sizeof(float)  // einstein_A
       + globals::nlines * sizeof(int) * 4  // elementindex, ionindex, upperlevelindex, lowerlevelindex
       ) /
      1024. / 1024;

  printlnlog("[info] mem_usage: linelist occupies {:.3f} MB (node shared memory)", linelist_mem_MB);

  printlnlog("establishing connection between transitions and sorted linelist...");

  auto const time_start_establish_linelist_connections = std::time(nullptr);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (int lineindex = 0; lineindex < globals::nlines; lineindex++) {
    if (lineindex % globals::node_nprocs != globals::rank_in_node) {
      continue;
    }

    const int element = globals::linelist.elementindex[lineindex];
    const int ion = globals::linelist.ionindex[lineindex];
    const int lowerlevel = globals::linelist.lowerlevelindex[lineindex];
    const int upperlevel = globals::linelist.upperlevelindex[lineindex];

    // there is never more than one transition per pair of levels,
    // so find the first up and the first down transition that match: element, ion, lowerlevel, upperlevel

    const auto upper_uniquelevelindex = get_uniquelevelindex(element, ion, upperlevel);
    const auto alltrans_startdown = get_alltrans_startdown(upper_uniquelevelindex);
    const auto ndowntrans = get_ndowntrans(upper_uniquelevelindex);
    int downtransid = -1;
    for (int i = alltrans_startdown; i < alltrans_startdown + ndowntrans; i++) {
      if (globals::alltrans.targetlevelindex[i] == lowerlevel) {
        downtransid = i;
        break;
      }
    }
    assert_always(downtransid != -1);
    alltrans_lineindex[downtransid] = lineindex;

    const auto lower_uniquelevelindex = get_uniquelevelindex(element, ion, lowerlevel);
    const auto alltrans_startup = get_alltrans_startup(lower_uniquelevelindex);
    const auto nuptrans = get_nuptrans(lower_uniquelevelindex);
    int uptransid = -1;
    for (int i = alltrans_startup; i < alltrans_startup + nuptrans; i++) {
      if (globals::alltrans.targetlevelindex[i] == upperlevel) {
        uptransid = i;
        break;
      }
    }
    assert_always(uptransid != -1);
    alltrans_lineindex[uptransid] = lineindex;
  }
  globals::alltrans.lineindex = alltrans_lineindex;

  printlnlog("  took {}s", std::time(nullptr) - time_start_establish_linelist_connections);
  MPI_Barrier(globals::mpi_comm_node);

  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      if (globals::elements[element].ions[ion].nlevels_groundterm <= 0) {
        if (single_ground_level) {
          globals::elements[element].ions[ion].nlevels_groundterm = 1;
        } else {
          globals::elements[element].ions[ion].nlevels_groundterm = calculate_nlevels_groundterm(element, ion);
        }
      }
    }
  }

  read_autoion_data();
  read_phixs_data();
}

void setup_cellcache() {
  globals::mutex_cellcachemacroatom.resize(get_includedlevels());

  // const int num_cellcache_slots = get_max_threads();
  const int num_cellcache_slots = 1;
  globals::cellcache.resize(num_cellcache_slots);

  for (int cellcachenum = 0; cellcachenum < num_cellcache_slots; cellcachenum++) {
    auto mem_usage_cellcache = 0zU;
    mem_usage_cellcache += sizeof(globals::CellCache);

    printlnlog("[info] input: initializing cellcache for thread {} ...", cellcachenum);

    globals::cellcache[cellcachenum].nonemptymgi = -1;

    const auto ncoolingterms = kpkt::ncoolingterms;
    mem_usage_cellcache += ncoolingterms * sizeof(double);
    resize_exactly(globals::cellcache[cellcachenum].cooling_contrib, ncoolingterms);
    std::ranges::fill(globals::cellcache[cellcachenum].cooling_contrib, 0.0);

    printlnlog("[info] mem_usage: cellcache coolinglist contribs for thread {} occupies {:.3f} MB", cellcachenum,
               ncoolingterms * sizeof(double) / 1024. / 1024.);

    auto allphixstargetcount = 0zU;
    auto chtransblocksize = 0zU;
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        const int nlevels = get_nlevels(element, ion);

        for (int level = 0; level < nlevels; level++) {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          allphixstargetcount += nphixstargets * sizeof(double);

          const int ndowntrans = get_ndowntrans(element, ion, level);
          const int nuptrans = get_nuptrans(element, ion, level);
          chtransblocksize += ((2 * ndowntrans) + nuptrans);
        }
      }
    }
    resize_exactly(globals::cellcache[cellcachenum].alllevels_pops, get_includedlevels());
    resize_exactly(globals::cellcache[cellcachenum].alllevels_maprocessrates, get_includedlevels());

    if (allphixstargetcount > 0) {
      resize_exactly(globals::cellcache[cellcachenum].allphixstargets_corrphotoioncoeff, allphixstargetcount);
      if constexpr (SEPARATE_STIMRECOMB) {
        resize_exactly(globals::cellcache[cellcachenum].allphixstargets_stimrecombcoeff, allphixstargetcount);
      }
    }
    mem_usage_cellcache +=
        (get_includedlevels() * (2 * sizeof(double) + sizeof(int))) + (allphixstargetcount * sizeof(double) * 2);

    assert_always(chtransblocksize <= std::numeric_limits<int>::max());
    mem_usage_cellcache += chtransblocksize * sizeof(double);
    if (chtransblocksize > 0) {
      resize_exactly(globals::cellcache[cellcachenum].allmacroatomictransitions, chtransblocksize);
    }

    for (int uniquelevelindex = 0; uniquelevelindex < get_includedlevels(); uniquelevelindex++) {
      std::ranges::fill(globals::cellcache[cellcachenum].alllevels_maprocessrates[uniquelevelindex], -99.);
    }

    assert_always(globals::nbfcontinua >= 0);
    resize_exactly(globals::cellcache[cellcachenum].allcont_departureratios, globals::nbfcontinua);
    resize_exactly(globals::cellcache[cellcachenum].allcont_nnlevel, globals::nbfcontinua);
    resize_exactly(globals::cellcache[cellcachenum].allcont_keep, globals::nbfcontinua);
    mem_usage_cellcache += 2 * globals::nbfcontinua * sizeof(double);

    printlnlog("[info] mem_usage: cellcache for thread {} occupies {:.3f} MB", cellcachenum,
               mem_usage_cellcache / 1024. / 1024.);
  }
}

void write_bflist_file() {
  resize_exactly(globals::bflist, globals::nbfcontinua);

  int i = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels_ionising(element, ion);
      for (int level = 0; level < nlevels; level++) {
        const auto nphixstargets = get_nphixstargets(element, ion, level);
        for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
          const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
          globals::bflist[i].elementindex = element;
          globals::bflist[i].ionindex = ion;
          globals::bflist[i].levelindex = level;
          globals::bflist[i].phixstargetindex = phixstargetindex;

          const int et = -1 - i;

          assert_always(et == get_emtype_continuum(element, ion, level, upperionlevel));

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
    bflist_file << globals::nbfcontinua << '\n';
    for (i = 0; i < globals::nbfcontinua; i++) {
      const int element = globals::bflist[i].elementindex;
      const int ion = globals::bflist[i].ionindex;
      const int level = globals::bflist[i].levelindex;
      const int phixstargetindex = globals::bflist[i].phixstargetindex;
      const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
      bflist_file << i << ' ' << element << ' ' << ion << ' ' << level << ' ' << upperionlevel << '\n';
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

        const bool has_superlevel = (nlevels > (nlevels_excited_nlte + 1));
        if (has_superlevel) {
          // If there are more levels that the ground state + the number of NLTE levels then we need an extra
          // slot to store data for the "superlevel", which is a representation of all the other levels that
          // are not treated in detail.
          globals::total_nlte_levels++;
          n_super_levels++;
        }

        assert_always(has_superlevel == ion_has_superlevel(element, ion));

        printlnlog("[input]  element {:2} Z={:2} ionstage {:2} has {:5} NLTE excited levels{}. Starting at index {}",
                   element, get_atomicnumber(element), get_ionstage(element, ion),
                   get_nlevels_excited_nlte(element, ion), has_superlevel ? " plus a superlevel" : "",
                   get_allnltelevelsindexstart(element, ion));
      }
    }
  }

  printlnlog("[input] Total NLTE levels: {}, of which {} are superlevels", globals::total_nlte_levels, n_super_levels);
}

void read_atomicdata() {
  read_atomicdata_files();

  // INITIALISE THE ABSORPTION/EMISSION COUNTERS ARRAYS
  if constexpr (RECORD_LINESTAT) {
    globals::ecounter.resize(globals::nlines);
    globals::acounter.resize(globals::nlines);
  }

  kpkt::setup_coolinglist();

  setup_cellcache();

  // Printout some information about the read-in model atom

  int includedionisinglevels = 0;
  int includedboundboundtransitions = 0;
  int includedphotoiontransitions = 0;
  printlnlog("[input] this simulation contains");
  printlnlog("----------------------------------");
  for (int element = 0; element < get_nelements(); element++) {
    printlnlog("[input]  element {} (Z={:2} {})", element, get_atomicnumber(element),
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
          "[input]    ionstage {}: {:4} levels ({} in groundterm, {:4} ionising) {:7} lines {:6} bf transitions "
          "(epsilon_ground: {:7.2f} eV)",
          get_ionstage(element, ion), get_nlevels(element, ion), get_nlevels_groundterm(element, ion),
          get_nlevels_ionising(element, ion), ion_bbtransitions, ion_photoiontransitions,
          epsilon(element, ion, 0) / EV);

      includedionisinglevels += get_nlevels_ionising(element, ion);
      includedphotoiontransitions += ion_photoiontransitions;
      includedboundboundtransitions += ion_bbtransitions;
    }
  }
  assert_always(includedphotoiontransitions == globals::nbfcontinua);
  assert_always(globals::nlines == includedboundboundtransitions);

  printlnlog("[input]  in total {} ions, {} levels ({} ionising), {} lines, {} photoionisation transitions",
             get_includedions(), get_includedlevels(), includedionisinglevels, globals::nlines, globals::nbfcontinua);

  write_bflist_file();

  setup_nlte_levels();
}

}  // anonymous namespace

// read input.txt, atomic data, and ejecta model
void input(int rank) {
  globals::n_titer = (globals::timestep < -1) ? 3 : 1;
  globals::lte_iteration = false;

  printlnlog("[info] input: do n_titer {} iterations per timestep", globals::n_titer);
  if (globals::n_titer > 1) {
#ifndef DO_TITER
    printlnlog("[fatal] input: n_titer > 1, but DO_TITER not defined ... abort");
    std::abort();
#endif
  } else if (globals::n_titer == 1) {
#ifdef DO_TITER
    printlnlog("[warning] input: n_titer = 1 but DO_TITER defined, remove DO_TITER to save memory");
#endif
  } else {
    printlnlog("[fatal] input: no valid value for n_titer selected");
    std::abort();
  }

  // Read in parameters from input.txt
  read_parameterfile(rank);

  // Read in parameters from vpkt.txt
  if (VPKT_ON) {
    vpkt::read_vpktparameterfile();
  }

  read_atomicdata();

  const auto time_before_barrier = std::time(nullptr);
  printlog("barrier after read_atomicdata(): time before barrier {}, ", static_cast<int>(time_before_barrier));
  MPI_Barrier(MPI_COMM_WORLD);
  printlnlog("time after barrier {} (waited {} seconds)", static_cast<int>(time(nullptr)),
             static_cast<int>(time(nullptr) - time_before_barrier));

  grid::read_ejecta_model();
}

// read the next line, skipping any comment lines beginning with '#'
auto get_noncommentline(std::istream& input, std::string& line) -> bool {
  while (true) {
    const bool linefound = !(!std::getline(input, line));
    if (!linefound) {
      return false;
    }
    if (!lineiscommentonly(line)) {
      return true;
    }
  }
}

// read input parameters from input.txt
void read_parameterfile(int rank) {
  auto file = fstream_required("input.txt", std::ios::in);

  std::string line;
  assert_always(get_noncommentline(file, line));

  std::int64_t pre_zseed = -1;
  std::istringstream{line} >> pre_zseed;

  if (pre_zseed > 0) {
    printlnlog("input.txt specified random number seed is {}\n", pre_zseed);
  } else {
    pre_zseed = get_rng_random_seed();
    // broadcast randomly-generated seed from rank 0 to all ranks
    MPI_Bcast_safe(pre_zseed, 0, MPI_COMM_WORLD);
    printlnlog("randomly-generated random number seed is {}", pre_zseed);
#if defined REPRODUCIBLE && REPRODUCIBLE
    printlnlog("ERROR: reproducible mode is on, so random number seed is required.");
    std::abort();
#endif
  }

  // For MPI parallelisation, the random seed is changed based on the rank of the process
  // Multi-threaded runs (OpenMP or stdpar) are not reproducible due to accumulation to shared memory, so the seed is
  // randomly generated
  const auto tid = get_thread_num();
  auto rngseed =
      (tid == 0) ? pre_zseed + static_cast<std::int64_t>(13 * (rank * get_max_threads() + tid)) : get_rng_random_seed();

  rng_seed(rngseed);
  printlnlog("rank {}: thread {} has rngseed {}", rank, tid, rngseed);

  // call it a few times
  for (int n = 0; n < 100; n++) {
    rng_uniform();
  }

  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::ntimesteps;  // number of time steps
  assert_always(globals::ntimesteps > 0);

  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::timestep_initial >>
      globals::timestep_finish;  // number of start and end time step
  printlnlog("input: timestep_start {} timestep_finish {}", globals::timestep_initial, globals::timestep_finish);
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

  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::gamma_kappagrey;  // use grey opacity for gammas?

  assert_always(get_noncommentline(file, line));  // UNUSED components of syn_dir

  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::opacity_case;  // opacity choice

  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::rho_crit_para;  // free parameter for calculation of rho_crit
  printlnlog("input: rho_crit_para {:g}", globals::rho_crit_para);
  // the calculation of rho_crit itself depends on the time, therefore it happens in grid_init and update_grid

  assert_always(get_noncommentline(file, line));
  int debug_packet = 0;
  std::istringstream{line} >> debug_packet;  // activate debug output for packet
  assert_always(debug_packet == -1);
  // select a negative value to deactivate

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

  // Wavelength (in Angstroms) at which the parameterisation of the radiation field
  // switches from the nebular approximation to LTE.
  float dum2{NAN};
  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> dum2;  // free parameter for calculation of rho_crit
  globals::nu_rfcut = CLIGHT / (dum2 * 1e-8);
  printlnlog("input: nu_rfcut {:g}", globals::nu_rfcut);

  // Sets the number of initial LTE timesteps for NLTE runs
  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::num_lte_timesteps;
  printlnlog("input: doing the first {} timesteps in LTE", globals::num_lte_timesteps);

  if (NT_ON) {
    if (NT_SOLVE_SPENCERFANO) {
      printlnlog("input: Non-thermal ionisation with a Spencer-Fano solution is switched on for this run.");
    } else {
      printlnlog("input: Non-thermal ionisation with the work function approximation is switched on for this run.");
    }
  } else {
    printlnlog("input: No non-thermal ionisation is used in this run.");
  }

  if (USE_LUT_PHOTOION) {
    printlnlog(
        "Corrphotoioncoeff is calculated from LTE lookup tables (ratecoeff.dat) and corrphotoionrenorm estimator.");
  } else {
    printlnlog(
        "Corrphotoioncoeff is calculated from the radiation field at each timestep in each modelgrid cell (no LUT).");
  }

  if (USE_LUT_BFHEATING) {
    printlnlog("bfheating coefficients are calculated from LTE lookup tables (ratecoeff.dat) and bfheatingestimator.");
  } else {
    printlnlog(
        "bfheating coefficients are calculated from the radiation field at each timestep in each modelgrid cell (no "
        "LUT).");
  }

  // Set up initial grey approximation?
  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::cell_is_optically_thick >> globals::num_grey_timesteps;
  printlnlog(
      "input: cells with Thomson optical depth > {:g} are treated in grey approximation for the first {} timesteps",
      globals::cell_is_optically_thick, globals::num_grey_timesteps);

  // Limit the number of bf-continua
  assert_always(get_noncommentline(file, line));
  int max_bf_continua = 0;
  std::istringstream{line} >> max_bf_continua;
  assert_always(max_bf_continua == -1);

  // for exspec: read number of MPI tasks
  assert_always(get_noncommentline(file, line));
  std::istringstream{line} >> globals::nprocs_exspec;

  // UNUSED: Extract line-of-sight dependent information of last emission for spectrum_res
  assert_always(get_noncommentline(file, line));

  // To reduce the work imbalance between different MPI tasks I introduced a diffusion
  // for kpkts, since it turned out that this work imbalance was largely dominated
  // by continuous collisional interactions. By introducing a diffusion time for kpkts
  // this loop is broken. The following two parameters control this approximation.
  // Parameter one (a float) gives the relative fraction of a time step which individual
  // kpkts live. Parameter two (an int) gives the number of time steps for which we
  // want to use this approximation
  assert_always(get_noncommentline(file, line));
  int n_kpktdiffusion_timesteps{0};
  float kpktdiffusion_timescale{0.};
  std::istringstream{line} >> kpktdiffusion_timescale >> n_kpktdiffusion_timesteps;
  kpkt::set_kpktdiffusion(kpktdiffusion_timescale, n_kpktdiffusion_timesteps);

  file.close();

  if (rank == 0 && !globals::simulation_continued_from_saved) {
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
        if (noncomment_linenum == 2) {
          // Number of start and end time step
          line = std::format("{:03d} {:03d}", nts, globals::timestep_finish);
        } else if (noncomment_linenum == 16) {
          // resume from gridsave file
          line = "1";  // Force continuation
        }
      }

      if (noncomment_linenum == 21) {
        // by default, exspec should use all available packet files
        globals::nprocs_exspec = globals::nprocs;
        line = std::format("{}", globals::nprocs_exspec);
      }

      if (noncomment_linenum < std::ssize(inputlinecomments)) {
        const int commentstart = 25;

        // truncate any existing comment on the line
        if (line.find('#') != std::string::npos) {
          line.resize(line.find('#'));
        }

        line.resize(commentstart, ' ');
        line.append("# ");
        line.append(inputlinecomments[noncomment_linenum]);
      }
    }

    fileout << line << '\n';
  }

  fileout.close();
  file.close();

  if (nts < 0) {
    std::rename("input.txt.tmp", "input-newrun.txt");  // back up the original for starting a new simulation
  } else {
    std::remove("input.txt");
    std::rename("input.txt.tmp", "input.txt");
  }

  printlnlog("done");
}

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

  // to limit the timestep durations
  // const double maxt = 0.5 * DAY;
  // for (int n = globals::ntimesteps - 1; n > 0; n--)
  // {
  //   if (globals::timesteps[n].width > maxt)
  //   {
  //     const double boundaryshift = globals::timesteps[n].width - maxt;
  //     globals::timesteps[n].width -= boundaryshift;
  //     globals::timesteps[n].start += boundaryshift;
  //     globals::timesteps[n - 1].width += boundaryshift;
  //   }
  //   else if (n < globals::ntimesteps - 1 && globals::timesteps[n + 1].width > maxt)
  //   {
  //     printout("TIME: Keeping logarithmic durations for timesteps <= %d\n", n);
  //   }
  // }
  // assert_always(globals::timesteps[0].width <= maxt); // no solution is possible with these constraints!

  // and add a dummy timestep which contains the endtime
  // of the calculation
  globals::timesteps[globals::ntimesteps].start = globals::tmax;
  globals::timesteps[globals::ntimesteps].mid = globals::tmax;
  globals::timesteps[globals::ntimesteps].width = 0.;

  // check consistency of start + width = start_next
  for (int n = 1; n < globals::ntimesteps; n++) {
    assert_always(
        fabs((globals::timesteps[n - 1].start + globals::timesteps[n - 1].width) / globals::timesteps[n].start) - 1 <
        0.001);
  }
  assert_always(
      fabs((globals::timesteps[globals::ntimesteps - 1].start + globals::timesteps[globals::ntimesteps - 1].width) /
           globals::tmax) -
          1 <
      0.001);
}

void write_timestep_file() {
  auto timestepfile = std::fstream("timesteps.out", std::ofstream::out | std::ofstream::trunc);
  assert_always(timestepfile.is_open());
  timestepfile << "#timestep tstart_days tmid_days twidth_days\n";
  for (int n = 0; n < globals::ntimesteps; n++) {
    timestepfile << n << ' ' << globals::timesteps[n].start / DAY << ' ' << globals::timesteps[n].mid / DAY << ' '
                 << globals::timesteps[n].width / DAY << '\n';
  }
}
