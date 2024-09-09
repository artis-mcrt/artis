#include "input.h"

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

#include <algorithm>
#include <array>
#include <cinttypes>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <ios>
#include <iterator>
#ifndef GPU_ON
#include <random>
#endif
#include <span>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "decay.h"
#include "globals.h"
#include "grid.h"
#include "kpkt.h"
#include "packet.h"
#include "ratecoeff.h"
#include "sn3d.h"
#include "vpkt.h"

namespace {

const int groundstate_index_in = 1;  // starting level index in the input files

struct Transition {
  int lower;
  int upper;
  float A;
  float coll_str;
  bool forbidden;
};  // only used temporarily during input

constexpr std::array<std::string_view, 24> inputlinecomments = {
    " 0: pre_zseed: specific random number seed if > 0 or random if negative",
    " 1: ntimesteps: number of timesteps",
    " 2: timestep_start timestep_finish: timestep number range start (inclusive) and stop (not inclusive)",
    " 3: tmin_days tmax_days: start and end times [day]",
    " 4: UNUSED nusyn_min_mev nusyn_max_mev: lowest and highest frequency to synthesise [MeV]",
    " 5: UNUSED nsyn_time: number of times for synthesis",
    " 6: UNUSED start and end times for synthesis",
    " 7: model_type: number of dimensions (1, 2, or 3)",
    " 8: UNUSED compute r-light curve (1: no estimators, 2: thin cells, 3: thick cells, 4: gamma-ray heating)",
    " 9: UNUSED n_out_it: number of iterations",
    "10: UNUSED: change speed of light by some factor. Change constants.h CLIGHT_PROP instead",
    "11: gamma_kappagrey: if >0: use grey opacity for gammas, if <0: use detailed opacity",
    "12: syn_dir: x, y, and z components of unit vector (will be normalised after input or randomised if zero length)",
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
    "22: do_emission_res: Extract line-of-sight dependent information of last emission for spectrum_res (1: yes, 2: "
    "no)",
    "23: kpktdiffusion_timescale n_kpktdiffusion_timesteps: kpkts diffuse x of a time step's length for the first y "
    "time steps"};

CellCachePhixsTargets *chphixstargetsblock{};

// use the temporary list, before it has been copied to node-shared memory
auto get_phixsupperlevel_tmp(const std::vector<PhotoionTarget> &temp_allphixstargets, const int element, const int ion,
                             const int level, const int phixstargetindex) {
  return temp_allphixstargets[globals::elements[element].ions[ion].levels[level].phixstargetstart + phixstargetindex]
      .levelindex;
}

void read_phixs_data_table(std::fstream &phixsfile, const int nphixspoints_inputtable, const int element,
                           const int lowerion, const int lowerlevel, const int upperion, int upperlevel_in,
                           std::vector<float> &tmpallphixs, auto &temp_allphixstargets, size_t *mem_usage_phixs,
                           const int phixs_file_version) {
  std::string phixsline;
  assert_always(globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargetstart == -1);
  globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargetstart =
      static_cast<int>(temp_allphixstargets.size());
  if (upperlevel_in >= 0) {  // file gives photoionisation to a single target state only
    int upperlevel = upperlevel_in - groundstate_index_in;
    assert_always(upperlevel >= 0);
    assert_always(globals::elements[element].ions[lowerion].levels[lowerlevel].nphixstargets == 0);
    globals::elements[element].ions[lowerion].levels[lowerlevel].nphixstargets = 1;
    *mem_usage_phixs += sizeof(PhotoionTarget);

    if (single_level_top_ion && (upperion == get_nions(element) - 1)) {
      // top ion has only one level, so send it to that level
      upperlevel = 0;
    }

    temp_allphixstargets.push_back({.probability = 1., .levelindex = upperlevel});
  } else {  // upperlevel < 0, indicating that a table of upper levels and their probabilities will follow
    int in_nphixstargets = 0;
    assert_always(get_noncommentline(phixsfile, phixsline));
    assert_always(std::stringstream(phixsline) >> in_nphixstargets);
    assert_always(in_nphixstargets > 0);
    // read in a table of target states and probabilities and store them
    if (!single_level_top_ion || upperion < get_nions(element) - 1)  // in case the top ion has nlevelsmax = 1
    {
      globals::elements[element].ions[lowerion].levels[lowerlevel].nphixstargets = in_nphixstargets;
      *mem_usage_phixs += in_nphixstargets * sizeof(PhotoionTarget);

      double probability_sum = 0.;
      for (int i = 0; i < in_nphixstargets; i++) {
        double phixstargetprobability{NAN};
        assert_always(get_noncommentline(phixsfile, phixsline));
        assert_always(std::stringstream(phixsline) >> upperlevel_in >> phixstargetprobability);
        const int upperlevel = upperlevel_in - groundstate_index_in;
        assert_always(upperlevel >= 0);
        assert_always(phixstargetprobability > 0);
        temp_allphixstargets.push_back({.probability = phixstargetprobability, .levelindex = upperlevel});

        probability_sum += phixstargetprobability;
      }
      if (fabs(probability_sum - 1.0) > 0.01) {
        printout("WARNING: photoionisation table for Z=%d ionstage %d has probabilities that sum to %g",
                 get_atomicnumber(element), get_ionstage(element, lowerion), probability_sum);
      }
    } else {  // file has table of target states and probabilities but our top ion is limited to one level
      globals::elements[element].ions[lowerion].levels[lowerlevel].nphixstargets = 1;
      *mem_usage_phixs += sizeof(PhotoionTarget);

      for (int i = 0; i < in_nphixstargets; i++) {
        assert_always(get_noncommentline(phixsfile, phixsline));
      }

      // send it to the ground state of the top ion
      temp_allphixstargets.push_back({.probability = 1., .levelindex = 0});
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
          temp_allphixstargets[globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargetstart +
                               phixstargetindex]
              .levelindex;
      if (upperlevel > get_maxrecombininglevel(element, lowerion + 1)) {
        globals::elements[element].ions[lowerion + 1].maxrecombininglevel = upperlevel;
      }
    }
  }

  *mem_usage_phixs += globals::NPHIXSPOINTS * sizeof(float);
  assert_always(tmpallphixs.size() % globals::NPHIXSPOINTS == 0);
  const auto tmpphixsstart = tmpallphixs.size();
  globals::elements[element].ions[lowerion].levels[lowerlevel].phixsstart = tmpphixsstart / globals::NPHIXSPOINTS;
  tmpallphixs.resize(tmpallphixs.size() + globals::NPHIXSPOINTS);

  auto *levelphixstable = &tmpallphixs[tmpphixsstart];
  if (phixs_file_version == 1) {
    assert_always(get_nphixstargets(element, lowerion, lowerlevel) == 1);
    assert_always(get_phixsupperlevel_tmp(temp_allphixstargets, element, lowerion, lowerlevel, 0) == 0);

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
      nugrid_in[i] = nu_edge + (energy * 13.6 * EV) / H;
      // the photoionisation cross-sections in the database are given in Mbarn=1e6 * 1e-28m^2
      // to convert to cgs units multiply by 1e-18
      phixs_in[i] = phixs * 1e-18;
    }
    const double nu_max = nugrid_in.back();

    // Now interpolate these cross-sections
    levelphixstable[0] = phixs_in[0];

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, nphixspoints_inputtable);
    gsl_spline_init(spline, nugrid_in.data(), phixs_in.data(), nphixspoints_inputtable);
    for (int i = 1; i < globals::NPHIXSPOINTS; i++) {
      const double nu = nu_edge * (1. + i * globals::NPHIXSNUINCREMENT);
      if (nu > nu_max) {
        levelphixstable[i] = phixs_in[nphixspoints_inputtable - 1] * pow(nu_max / nu, 3);
      } else {
        levelphixstable[i] = gsl_spline_eval(spline, nu, acc);
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
      levelphixstable[i] = phixs * 1e-18;
      // fprintf(database_file,"%g %g\n", nutable[i], phixstable[i]);
    }
  }

  // nbfcontinua++;
  // printout("[debug] element %d, ion %d, level %d: phixs exists %g\n",element,lowerion,lowerlevel,phixs*1e-18);
  globals::nbfcontinua += get_nphixstargets(element, lowerion, lowerlevel);
  if (lowerlevel == 0 && get_nphixstargets(element, lowerion, lowerlevel) > 0) {
    globals::nbfcontinua_ground++;
  }
}

void read_phixs_file(const int phixs_file_version, std::vector<float> &tmpallphixs,
                     std::vector<PhotoionTarget> &temp_allphixstargets) {
  size_t mem_usage_phixs = 0;

  printout("readin phixs data from %s\n", phixsdata_filenames[phixs_file_version]);

  auto phixsfile = fstream_required(phixsdata_filenames[phixs_file_version], std::ios::in);
  std::string phixsline;

  if (phixs_file_version == 1 && phixs_file_version_exists[2]) {
    printout(
        "using NPHIXSPOINTS = %d and NPHIXSNUINCREMENT = %lg from phixsdata_v2.txt to interpolate "
        "phixsdata.txt data\n",
        globals::NPHIXSPOINTS, globals::NPHIXSNUINCREMENT);
    last_phixs_nuovernuedge = (1.0 + globals::NPHIXSNUINCREMENT * (globals::NPHIXSPOINTS - 1));
  } else if (phixs_file_version == 1) {
    globals::NPHIXSPOINTS = 100;
    globals::NPHIXSNUINCREMENT = .1;
    // not exactly where the last point is, but classic integrals go from nu_edge to 10*nu_edge
    last_phixs_nuovernuedge = 10;
    printout("using NPHIXSPOINTS = %d and NPHIXSNUINCREMENT = %lg set in input.cc\n", globals::NPHIXSPOINTS,
             globals::NPHIXSNUINCREMENT);
  } else {
    assert_always(phixsfile >> globals::NPHIXSPOINTS);
    assert_always(globals::NPHIXSPOINTS > 0);
    assert_always(phixsfile >> globals::NPHIXSNUINCREMENT);
    assert_always(globals::NPHIXSNUINCREMENT > 0.);
    last_phixs_nuovernuedge = (1.0 + globals::NPHIXSNUINCREMENT * (globals::NPHIXSPOINTS - 1));
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
    if (phixs_file_version == 1) {
      assert_always(std::istringstream(phixsline) >> Z >> upperionstage >> upperlevel_in >> lowerionstage >>
                    lowerlevel_in >> nphixspoints_inputtable);
    } else {
      assert_always(std::istringstream(phixsline) >> Z >> upperionstage >> upperlevel_in >> lowerionstage >>
                    lowerlevel_in >> phixs_threshold_ev);
      nphixspoints_inputtable = globals::NPHIXSPOINTS;
    }
    assert_always(Z > 0);
    assert_always(upperionstage >= 2);
    assert_always(lowerionstage >= 1);

    const int element = get_elementindex(Z);

    // store only photoionization crosssections for elements that are part of the current model atom
    bool skip_this_phixs_table = true;  // will be set to false for good data
    if (element >= 0 && get_nions(element) > 0) {
      // translate readin ionstages to ion indices

      const int upperion = upperionstage - get_ionstage(element, 0);
      const int lowerion = lowerionstage - get_ionstage(element, 0);
      const int lowerlevel = lowerlevel_in - groundstate_index_in;
      assert_always(lowerionstage >= 0);
      assert_always(lowerlevel >= 0);

      // store only photoionization crosssections for ions that are part of the current model atom
      if (lowerion >= 0 && upperion < get_nions(element) && lowerlevel < get_ionisinglevels(element, lowerion)) {
        read_phixs_data_table(phixsfile, nphixspoints_inputtable, element, lowerion, lowerlevel, upperion,
                              upperlevel_in, tmpallphixs, temp_allphixstargets, &mem_usage_phixs, phixs_file_version);

        skip_this_phixs_table = false;
      }
    }

    if (skip_this_phixs_table) {  // for ions or elements that are not part of the current model atom, proceed through
                                  // the table and throw away the data
      if (upperlevel_in < 0) {    // a table of target states and probabilities will follow, so read past those lines
        int nphixstargets = 0;
        assert_always(get_noncommentline(phixsfile, phixsline));
        assert_always(std::stringstream(phixsline) >> nphixstargets);
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

  printout("[info] mem_usage: photoionisation tables occupy %.3f MB\n", mem_usage_phixs / 1024. / 1024.);
}

constexpr auto downtranslevelstart(const int level) {
  // each level index is associated with a block of size levelindex spanning all possible down transitions.
  // so use the formula for the sum of 1 + 2 + 3 + 4 + ... + level
  return level * (level + 1) / 2;
}

void read_ion_levels(std::fstream &adata, const int element, const int ion, const int nions, const int nlevels,
                     int nlevelsmax, const double energyoffset, const double ionpot) {
  for (int level = 0; level < nlevels; level++) {
    int levelindex_in = 0;
    double levelenergy{NAN};
    double statweight{NAN};
    int ntransitions = 0;
    std::string line;
    assert_always(get_noncommentline(adata, line));
    assert_always(std::istringstream(line) >> levelindex_in >> levelenergy >> statweight >> ntransitions);
    assert_always(levelindex_in == level + groundstate_index_in);

    if (level < nlevelsmax) {
      const double currentlevelenergy = (energyoffset + levelenergy) * EV;
      globals::elements[element].ions[ion].levels[level].nphixstargets = 0;
      globals::elements[element].ions[ion].levels[level].phixsstart = -1;
      globals::elements[element].ions[ion].levels[level].phixstargetstart = -1;
      globals::elements[element].ions[ion].levels[level].epsilon = currentlevelenergy;
      globals::elements[element].ions[ion].levels[level].stat_weight = statweight;
      assert_always(statweight > 0.);

      // The level contributes to the ionisinglevels if its energy
      // is below the ionization potential and the level doesn't
      // belong to the topmost ion included.
      // Rate coefficients are only available for ionising levels.
      if (levelenergy < ionpot && ion < nions - 1) {
        globals::elements[element].ions[ion].ionisinglevels++;
      }

      set_ndowntrans(element, ion, level, 0);
      set_nuptrans(element, ion, level, 0);
    } else {
      // globals::elements[element].ions[ion].levels[nlevelsmax - 1].stat_weight += statweight;
    }
  }
}

void read_ion_transitions(std::fstream &ftransitiondata, const int tottransitions_in_file, int &tottransitions,
                          std::vector<Transition> &iontransitiontable, const int nlevels_requiretransitions,
                          const int nlevels_requiretransitions_upperlevels) {
  iontransitiontable.clear();
  iontransitiontable.reserve(tottransitions);

  std::string line;

  if (tottransitions == 0) {
    // we will not read in any transitions, just skip past these lines in the file
    for (int i = 0; i < tottransitions_in_file; i++) {
      assert_always(getline(ftransitiondata, line));
    }
  } else {
    // will be autodetected from first table row. old format had an index column and no collstr or forbidden columns
    bool oldtransitionformat = false;

    int prev_upper = -1;
    int prev_lower = 0;
    for (int i = 0; i < tottransitions_in_file; i++) {
      int lower_in = -1;
      int upper_in = -1;
      float A = 0;
      float coll_str = -1.;
      int intforbidden = 0;
      assert_always(getline(ftransitiondata, line));
      if (i == 0) {
        std::istringstream ss(line);
        std::string word;
        int column_count = 0;
        while (ss >> word) {
          column_count++;
        }
        assert_always(column_count == 4 || column_count == 5);
        oldtransitionformat = (column_count == 4);
      }
      if (!oldtransitionformat) {
        assert_always(std::istringstream(line) >> lower_in >> upper_in >> A >> coll_str >> intforbidden);
      } else {
        int transindex = 0;  // not used
        assert_always(std::istringstream(line) >> transindex >> lower_in >> upper_in >> A);
      }
      const int lower = lower_in - groundstate_index_in;
      const int upper = upper_in - groundstate_index_in;
      assert_always(lower >= 0);
      assert_always(upper >= 0);

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
          // printout("+adding transition index %d Z=%02d ionstage %d lower %d upper %d\n", i, Z, ionstage, prev_lower,
          // tmplevel);
          tottransitions++;
          assert_always(tmplevel >= 0);
          iontransitiontable.push_back(
              {.lower = prev_lower, .upper = tmplevel, .A = 0., .coll_str = -2., .forbidden = true});
        }
      }

      iontransitiontable.push_back(
          {.lower = lower, .upper = upper, .A = A, .coll_str = coll_str, .forbidden = (intforbidden == 1)});
      // printout("index %d, lower %d, upper %d, A %g\n",transitionindex,lower,upper,A);
      //  printout("reading transition index %d lower %d upper %d\n", i, transitiontable[i].lower,
      //  transitiontable[i].upper);
      prev_lower = lower;
      prev_upper = upper;
    }
  }
}

void add_transitions_to_unsorted_linelist(const int element, const int ion, const int nlevelsmax,
                                          const std::vector<Transition> &transitiontable,
                                          std::vector<int> &iondowntranstmplineindicies, int &lineindex,
                                          std::vector<TransitionLine> &temp_linelist,
                                          std::vector<LevelTransition> &temp_alltranslist) {
  const int lineindex_initial = lineindex;
  ptrdiff_t totupdowntrans = 0;
  // pass 0 to get transition counts of each level
  // pass 1 to allocate and fill transition arrays
  for (int pass = 0; pass < 2; pass++) {
    lineindex = lineindex_initial;
    if (pass == 1) {
      int alltransindex = temp_alltranslist.size();
      temp_alltranslist.resize(temp_alltranslist.size() + totupdowntrans);
      for (int level = 0; level < nlevelsmax; level++) {
        globals::elements[element].ions[ion].levels[level].alltrans_startdown = alltransindex;
        alltransindex += get_ndowntrans(element, ion, level);
        alltransindex += get_nuptrans(element, ion, level);

        set_ndowntrans(element, ion, level, 0);
        set_nuptrans(element, ion, level, 0);
      }
    }

    std::ranges::fill(iondowntranstmplineindicies, -99);

    totupdowntrans = 0;
    for (const auto &transition : transitiontable) {
      const int level = transition.upper;
      const int lowerlevel = transition.lower;
      if (pass == 0) {
        assert_always(lowerlevel >= 0);
        assert_always(level > lowerlevel);
      }

      if ((lowerlevel >= nlevelsmax) || (level >= nlevelsmax)) {
        continue;
      }
      const double nu_trans = (epsilon(element, ion, level) - epsilon(element, ion, lowerlevel)) / H;
      if (!(nu_trans > 0)) {
        continue;
      }

      // Make sure that we don't allow duplicate. In that case take only the lines
      // first occurrence
      int &downtranslineindex = iondowntranstmplineindicies[downtranslevelstart(level) + lowerlevel];

      // negative means that the transition hasn't been seen yet
      if (downtranslineindex < 0) {
        downtranslineindex = lineindex++;

        const int nupperdowntrans = get_ndowntrans(element, ion, level) + 1;
        set_ndowntrans(element, ion, level, nupperdowntrans);

        const int nloweruptrans = get_nuptrans(element, ion, lowerlevel) + 1;
        set_nuptrans(element, ion, lowerlevel, nloweruptrans);

        totupdowntrans += 2;

        if (pass == 1 && globals::rank_in_node == 0) {
          const auto g_ratio = stat_weight(element, ion, level) / stat_weight(element, ion, lowerlevel);
          const float f_ul = g_ratio * ME * pow(CLIGHT, 3) / (8 * pow(QE * nu_trans * PI, 2)) * transition.A;
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

          temp_alltranslist[globals::elements[element].ions[ion].levels[level].alltrans_startdown + nupperdowntrans -
                            1] = {.lineindex = -1,
                                  .targetlevelindex = lowerlevel,
                                  .einstein_A = transition.A,
                                  .coll_str = transition.coll_str,
                                  .osc_strength = f_ul,
                                  .forbidden = transition.forbidden};
          temp_alltranslist[globals::elements[element].ions[ion].levels[lowerlevel].alltrans_startdown +
                            get_ndowntrans(element, ion, lowerlevel) + nloweruptrans - 1] = {
              .lineindex = -1,
              .targetlevelindex = level,
              .einstein_A = transition.A,
              .coll_str = transition.coll_str,
              .osc_strength = f_ul,
              .forbidden = transition.forbidden};
        }

      } else if (pass == 1 && globals::rank_in_node == 0) {
        // This is a new branch to deal with lines that have different types of transition. It should trip after a
        // transition is already known.

        if ((temp_linelist[downtranslineindex].elementindex != element) ||
            (temp_linelist[downtranslineindex].ionindex != ion) ||
            (temp_linelist[downtranslineindex].upperlevelindex != level) ||
            (temp_linelist[downtranslineindex].lowerlevelindex != lowerlevel)) {
          printout("[input] Failure to identify level pair for duplicate bb-transition ... going to abort now\n");
          printout("[input]   element %d ion %d targetlevel %d level %d\n", element, ion, lowerlevel, level);
          printout("[input]   transitions[level].to[targetlevel]=lineindex %d\n", downtranslineindex);
          printout("[input]   A_ul %g, coll_str %g\n", transition.A, transition.coll_str);
          printout(
              "[input]   globals::linelist[lineindex].elementindex %d, "
              "globals::linelist[lineindex].ionindex %d, globals::linelist[lineindex].upperlevelindex "
              "%d, globals::linelist[lineindex].lowerlevelindex %d\n",
              temp_linelist[downtranslineindex].elementindex, temp_linelist[downtranslineindex].ionindex,
              temp_linelist[downtranslineindex].upperlevelindex, temp_linelist[downtranslineindex].lowerlevelindex);
          std::abort();
        }

        const auto g_ratio = stat_weight(element, ion, level) / stat_weight(element, ion, lowerlevel);
        const float f_ul = g_ratio * ME * pow(CLIGHT, 3) / (8 * pow(QE * nu_trans * PI, 2)) * transition.A;

        const int nupperdowntrans = get_ndowntrans(element, ion, level);
        auto &downtransition = get_downtranslist(element, ion, level)[nupperdowntrans];

        // this is what classic did, but it is not quite correct. The downtrans list should be searched to find the
        // correct index, not just using the last one. It probably works for the case where the transitions are sorted,
        // but the assertion tripped on C IV in the classic dataset.
        // assert_always(downtransition.targetlevelindex == lowerlevel);

        downtransition.einstein_A += transition.A;
        downtransition.osc_strength += f_ul;
        downtransition.coll_str = std::max(downtransition.coll_str, transition.coll_str);

        const int nloweruptrans = get_nuptrans(element, ion, lowerlevel);
        auto &uptransition = get_uptranslist(element, ion, lowerlevel)[nloweruptrans];

        // as above, the downtrans list should be searched to find the correct index instead of using the last one.
        // assert_always(uptransition.targetlevelindex == level);

        uptransition.einstein_A += transition.A;
        uptransition.osc_strength += f_ul;
        uptransition.coll_str = std::max(uptransition.coll_str, transition.coll_str);
      }
    }
  }
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
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
    const float g_a = stat_weight(element, ion, level_a);

    for (int level_b = 0; level_b < level_a; level_b++) {
      const float g_b = stat_weight(element, ion, level_b);
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

  if (nu_edge < globals::groundcont[0].nu_edge) {
    return -1;
  }

  int i = 1;
  for (i = 1; i < globals::nbfcontinua_ground; i++) {
    if (nu_edge < globals::groundcont[i].nu_edge) {
      break;
    }
  }

  if (i == globals::nbfcontinua_ground) {
    const int element = globals::groundcont[i - 1].element;
    const int ion = globals::groundcont[i - 1].ion;
    if (element == element_in && ion == ion_in && level_in == 0) {
      return i - 1;
    }

    printout(
        "[fatal] search_groundphixslist: element %d, ion %d, level %d has edge_frequency %g equal to the "
        "bluest ground-level continuum\n",
        element_in, ion_in, level_in, nu_edge);
    printout(
        "[fatal] search_groundphixslist: bluest ground level continuum is element %d, ion %d at "
        "nu_edge %g\n",
        element, ion, globals::groundcont[i - 1].nu_edge);
    printout("[fatal] search_groundphixslist: i %d, nbfcontinua_ground %d\n", i, globals::nbfcontinua_ground);
    printout(
        "[fatal] This shouldn't happen, is hoewever possible if there are multiple levels in the adata file at "
        "energy=0\n");
    for (int looplevels = 0; looplevels < get_nlevels(element_in, ion_in); looplevels++) {
      printout("[fatal]   element %d, ion %d, level %d, energy %g\n", element_in, ion_in, looplevels,
               epsilon(element_in, ion_in, looplevels));
    }
    printout("[fatal] Abort omitted ... MAKE SURE ATOMIC DATA ARE CONSISTENT\n");
    return i - 1;
    // abort();
  }

  const double left_diff = nu_edge - globals::groundcont[i - 1].nu_edge;
  const double right_diff = globals::groundcont[i].nu_edge - nu_edge;
  return (left_diff <= right_diff) ? i - 1 : i;
}

// set up the photoionisation transition lists
// and temporary gamma/kappa lists for each thread
void setup_phixs_list() {
  printout("[info] read_atomicdata: number of bfcontinua %d\n", globals::nbfcontinua);
  printout("[info] read_atomicdata: number of ground-level bfcontinua %d\n", globals::nbfcontinua_ground);

  if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
    globals::groundcont.resize(globals::nbfcontinua_ground);

    int groundcontindex = 0;
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
        assert_always(groundcontindex < globals::nbfcontinua_ground);

        globals::groundcont[groundcontindex] = {.nu_edge = nu_edge, .element = element, .ion = ion};

        groundcontindex++;
      }
    }
    assert_always(groundcontindex == globals::nbfcontinua_ground);
    std::ranges::SORT_OR_STABLE_SORT(globals::groundcont, std::ranges::less{}, &GroundPhotoion::nu_edge);
  }

  auto *nonconstallcont =
      static_cast<FullPhotoionTransition *>(malloc(globals::nbfcontinua * sizeof(FullPhotoionTransition)));
  printout("[info] mem_usage: photoionisation list occupies %.3f MB\n",
           globals::nbfcontinua * (sizeof(FullPhotoionTransition)) / 1024. / 1024.);
  int allcontindex = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions - 1; ion++) {
      if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
        globals::elements[element].ions[ion].groundcontindex =
            static_cast<int>(std::ranges::find_if(globals::groundcont,
                                                  [=](const auto &groundcont) {
                                                    return (groundcont.element == element) && (groundcont.ion == ion);
                                                  }) -
                             globals::groundcont.begin());
        if (globals::elements[element].ions[ion].groundcontindex >= globals::nbfcontinua_ground) {
          globals::elements[element].ions[ion].groundcontindex = -1;
        }
      }
      const int nlevels = get_ionisinglevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        const int nphixstargets = get_nphixstargets(element, ion, level);

        for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
          const double nu_edge = get_phixs_threshold(element, ion, level, phixstargetindex) / H;

          assert_always(allcontindex < globals::nbfcontinua);
          nonconstallcont[allcontindex].nu_edge = nu_edge;
          nonconstallcont[allcontindex].element = element;
          nonconstallcont[allcontindex].ion = ion;
          nonconstallcont[allcontindex].level = level;
          nonconstallcont[allcontindex].phixstargetindex = phixstargetindex;
          nonconstallcont[allcontindex].probability = get_phixsprobability(element, ion, level, phixstargetindex);
          nonconstallcont[allcontindex].upperlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);

          if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
            const double nu_edge_target0 = get_phixs_threshold(element, ion, level, 0) / H;
            const auto groundcontindex = search_groundphixslist(nu_edge_target0, element, ion, level);
            nonconstallcont[allcontindex].index_in_groundphixslist = groundcontindex;

            globals::elements[element].ions[ion].levels[level].closestgroundlevelcont = groundcontindex;
          }
          allcontindex++;
        }
      }
    }
  }

  assert_always(allcontindex == globals::nbfcontinua);
  assert_always(globals::nbfcontinua >= 0);  // was initialised as -1 before startup

  globals::bfestimcount = 0;
  if (globals::nbfcontinua > 0) {
    // indicies above were temporary only. continuum index should be to the sorted list
    std::ranges::SORT_OR_STABLE_SORT(std::span(nonconstallcont, globals::nbfcontinua), std::ranges::less{},
                                     &FullPhotoionTransition::nu_edge);

    globals::bfestim_nu_edge.clear();
    for (int i = 0; i < globals::nbfcontinua; i++) {
      auto &cont = nonconstallcont[i];
      if (DETAILED_BF_ESTIMATORS_ON &&
          LEVEL_HAS_BFEST(get_atomicnumber(cont.element), get_ionstage(cont.element, cont.ion), cont.level)) {
        cont.bfestimindex = globals::bfestimcount;
        globals::bfestim_nu_edge.push_back(cont.nu_edge);
        globals::bfestimcount++;
      } else {
        cont.bfestimindex = -1;
      }
    }

    globals::allcont_nu_edge.resize(globals::nbfcontinua, 0.);
    globals::bfestim_nu_edge.shrink_to_fit();
    assert_always(globals::bfestimcount == std::ssize(globals::bfestim_nu_edge));
  }
  printout("[info] bound-free estimators track bfestimcount %d photoionisation transitions\n", globals::bfestimcount);

  if (globals::nbfcontinua > 0) {
    for (int i = 0; i < globals::nbfcontinua; i++) {
      globals::allcont_nu_edge[i] = nonconstallcont[i].nu_edge;
    }

    setup_photoion_luts();

    for (int i = 0; i < globals::nbfcontinua; i++) {
      const int element = nonconstallcont[i].element;
      const int ion = nonconstallcont[i].ion;
      const int level = nonconstallcont[i].level;
      nonconstallcont[i].photoion_xs = get_phixs_table(element, ion, level);
      assert_always(nonconstallcont[i].photoion_xs != nullptr);
    }
  }
  globals::allcont = nonconstallcont;
  nonconstallcont = nullptr;
}

void read_phixs_data() {
  globals::nbfcontinua_ground = 0;
  globals::nbfcontinua = 0;
  std::vector<float> tmpallphixs;
  std::vector<PhotoionTarget> temp_allphixstargets;

  // read in photoionisation cross sections
  phixs_file_version_exists[0] = false;
  phixs_file_version_exists[1] = std::filesystem::exists(phixsdata_filenames[1]);
  phixs_file_version_exists[2] = std::filesystem::exists(phixsdata_filenames[2]);

#ifdef MPI_ON
  // just in case the file system was faulty and the ranks disagree on the existence of the files
  MPI_Allreduce(MPI_IN_PLACE, phixs_file_version_exists.data(), 3, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);
#endif
  assert_always(phixs_file_version_exists[1] || phixs_file_version_exists[2]);  // at least one must exist
  if (phixs_file_version_exists[1] && phixs_file_version_exists[2]) {
    printout(
        "Reading two phixs files: Reading phixsdata_v2.txt first so we use NPHIXSPOINTS and NPHIXSNUINCREMENT "
        "from phixsdata_v2.txt to interpolate the phixsdata.txt data\n");
  }
  if (phixs_file_version_exists[2]) {
    read_phixs_file(2, tmpallphixs, temp_allphixstargets);
  }
  if (phixs_file_version_exists[1]) {
    read_phixs_file(1, tmpallphixs, temp_allphixstargets);
  }

  int cont_index = 0;
  ptrdiff_t nbftables = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        const int nphixstargets = get_nphixstargets(element, ion, level);
        globals::elements[element].ions[ion].levels[level].cont_index = (nphixstargets > 0) ? cont_index : -1;
        cont_index += nphixstargets;
        if (nphixstargets > 0) {
          nbftables++;
        }
      }

      // below is just an extra warning consistency check
      const int nlevels_groundterm = globals::elements[element].ions[ion].nlevels_groundterm;

      // all levels in the ground term should be photoionisation targets from the lower ground state
      if (ion > 0 && ion < get_nions(element) - 1) {
        const int nphixstargets = get_nphixstargets(element, ion - 1, 0);
        if (nphixstargets > 0 && get_phixsupperlevel_tmp(temp_allphixstargets, element, ion - 1, 0, 0) == 0) {
          const int phixstargetlevels =
              get_phixsupperlevel_tmp(temp_allphixstargets, element, ion - 1, 0, nphixstargets - 1) + 1;

          if (nlevels_groundterm != phixstargetlevels) {
            printout("WARNING: Z=%d ionstage %d nlevels_groundterm %d phixstargetlevels(ion-1) %d.\n",
                     get_atomicnumber(element), get_ionstage(element, ion), nlevels_groundterm, phixstargetlevels);
            // if (nlevels_groundterm < phixstargetlevels)
            // {
            //   printout("  -> setting to %d\n", phixstargetlevels);
            //   globals::elements[element].ions[ion].nlevels_groundterm = phixstargetlevels;
            // }
          }
        }
      }
    }
  }

  printout("cont_index %d\n", cont_index);

  if (!tmpallphixs.empty()) {
    assert_always((nbftables * globals::NPHIXSPOINTS) == std::ssize(tmpallphixs));

    // copy the photoionisation tables into one contiguous block of memory
#ifdef MPI_ON
    MPI_Win win_allphixsblock = MPI_WIN_NULL;

    const auto [_, noderank_points] =
        get_range_chunk(std::ssize(tmpallphixs), globals::node_nprocs, globals::rank_in_node);

    auto size = static_cast<MPI_Aint>(noderank_points * sizeof(float));
    int disp_unit = sizeof(float);
    MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node, &globals::allphixs,
                            &win_allphixsblock);
    MPI_Win_shared_query(win_allphixsblock, MPI_PROC_NULL, &size, &disp_unit, &globals::allphixs);

    MPI_Barrier(MPI_COMM_WORLD);
#else
    globals::allphixs = static_cast<float *>(malloc(tmpallphixs.size() * sizeof(float)));
#endif

    assert_always(globals::allphixs != nullptr);

    std::copy_n(tmpallphixs.cbegin(), tmpallphixs.size(), globals::allphixs);

#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    tmpallphixs.clear();
    tmpallphixs.shrink_to_fit();
  }

  if (!temp_allphixstargets.empty()) {
    assert_always(cont_index == std::ssize(temp_allphixstargets));

    // copy the photoionisation tables into one contiguous block of memory
#ifdef MPI_ON
    MPI_Win win_allphixstargetsblock = MPI_WIN_NULL;

    const auto [_, noderank_phixstargetcount] =
        get_range_chunk(std::ssize(temp_allphixstargets), globals::node_nprocs, globals::rank_in_node);

    auto size = static_cast<MPI_Aint>(noderank_phixstargetcount * sizeof(PhotoionTarget));
    int disp_unit = sizeof(float);
    MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node, &globals::allphixstargets,
                            &win_allphixstargetsblock);
    MPI_Win_shared_query(win_allphixstargetsblock, MPI_PROC_NULL, &size, &disp_unit, &globals::allphixstargets);

    MPI_Barrier(MPI_COMM_WORLD);
#else
    globals::allphixstargets =
        static_cast<PhotoionTarget *>(malloc(temp_allphixstargets.size() * sizeof(PhotoionTarget)));
#endif

    assert_always(globals::allphixstargets != nullptr);

    std::copy_n(temp_allphixstargets.cbegin(), temp_allphixstargets.size(), globals::allphixstargets);

#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    temp_allphixstargets.clear();
    temp_allphixstargets.shrink_to_fit();
  }

  setup_phixs_list();
}

void read_atomicdata_files() {
  int totaluptrans = 0;
  int totaldowntrans = 0;

  auto compositiondata = fstream_required("compositiondata.txt", std::ios::in);

  auto adata = fstream_required("adata.txt", std::ios::in);

  printout("single_level_top_ion: %s\n", single_level_top_ion ? "true" : "false");
  printout("single_ground_level: %s\n", single_ground_level ? "true" : "false");
  // initialize atomic data structure to number of elements
  int nelements_in = 0;
  assert_always(compositiondata >> nelements_in);
  globals::elements.resize(nelements_in);

  std::vector<TransitionLine> temp_linelist;
  std::vector<LevelTransition> temp_alltranslist;

  std::vector<Transition> iontransitiontable;
  std::vector<int> iondowntranstmplineindicies;

  // temperature to determine relevant ionstages
  int T_preset = 0;
  assert_always(compositiondata >> T_preset);
  assert_always(T_preset == 0);  // no longer in use
  int homogeneous_abundances = 0;
  assert_always(compositiondata >> homogeneous_abundances);
  assert_always(homogeneous_abundances == 0);  // no longer in use

  // open transition data file
  auto ftransitiondata = fstream_required("transitiondata.txt", std::ios::in);

  int lineindex = 0;         // counter to determine the total number of lines
  int uniqueionindex = 0;    // index into list of all ions of all elements
  int uniquelevelindex = 0;  // index into list of all levels of all ions of all elements
  int nbfcheck = 0;
  for (int element = 0; element < get_nelements(); element++) {
    // read information about the next element which should be stored to memory
    int Z = 0;
    int nions = 0;
    int lowermost_ionstage = 0;
    int uppermost_ionstage = 0;
    int nlevelsmax_readin = 0;
    double uniformabundance{NAN};  // no longer in use mode for setting uniform abundances
    double mass_amu{NAN};
    assert_always(compositiondata >> Z >> nions >> lowermost_ionstage >> uppermost_ionstage >> nlevelsmax_readin >>
                  uniformabundance >> mass_amu);
    printout("readin compositiondata: next element Z %d, nions %d, lowermost %d, uppermost %d, nlevelsmax %d\n", Z,
             nions, lowermost_ionstage, uppermost_ionstage, nlevelsmax_readin);
    assert_always(Z > 0);
    assert_always(nions >= 0);
    assert_always(nions == 0 || (nions == uppermost_ionstage - lowermost_ionstage + 1));
    assert_always(uniformabundance >= 0);
    assert_always(mass_amu >= 0);

    // write this element's data to memory
    globals::elements[element].anumber = Z;
    globals::elements[element].nions = nions;
    globals::elements[element].initstablemeannucmass = mass_amu * MH;
    globals::elements[element].uniqueionindexstart = uniqueionindex;

    // Initialize the elements ionlist
    globals::elements[element].ions = static_cast<Ion *>(malloc(nions * sizeof(Ion)));
    assert_always(globals::elements[element].ions != nullptr);

    // now read in data for all ions of the current element. before doing so initialize
    // energy scale for the current element (all level energies are stored relative to
    // the ground level of the neutral ion)
    double energyoffset = 0.;
    double ionpot = 0.;
    for (int ion = 0; ion < nions; ion++) {
      int nlevelsmax = nlevelsmax_readin;
      // printout("element %d ion %d\n", element, ion);
      // calculate the current levels ground level energy
      assert_always(ionpot >= 0);
      energyoffset += ionpot;

      // read information for the elements next ionstage
      int adata_Z_in = -1;
      int ionstage = -1;
      int nlevels = 0;

      while (adata_Z_in != Z || ionstage != lowermost_ionstage + ion)  // skip over this ion block
      {
        if (adata_Z_in == Z) {
          printout("increasing energyoffset by ionpot %g\n", ionpot);
          energyoffset += ionpot;
        }
        for (int i = 0; i < nlevels; i++) {
          double levelenergy{NAN};
          double statweight{NAN};
          int levelindex = 0;
          int ntransitions = 0;
          std::string line;
          std::getline(adata, line);

          assert_always(std::istringstream(line) >> levelindex >> levelenergy >> statweight >> ntransitions);
        }

        std::string line;
        assert_always(get_noncommentline(adata, line));
        assert_always(std::istringstream(line) >> adata_Z_in >> ionstage >> nlevels >> ionpot);
      }

      printout("adata header matched: Z %d, ionstage %d, nlevels %d\n", adata_Z_in, ionstage, nlevels);

      if (single_level_top_ion && ion == nions - 1)  // limit the top ion to one level and no transitions
      {
        nlevelsmax = 1;
      }

      // if (adata_Z_in == 26 && ionstage == 1)
      // {
      //   nlevelsmax = 5;
      // }

      if (nlevelsmax < 0) {
        nlevelsmax = nlevels;
      } else if (nlevels >= nlevelsmax) {
        printout("[info] read_atomicdata: reduce number of levels from %d to %d for Z %2d ionstage %d\n", nlevels,
                 nlevelsmax, adata_Z_in, ionstage);
      } else {
        printout(
            "[warning] read_atomicdata: requested nlevelsmax=%d > nlevels=%d for ion %d of element %d ... reduced "
            "nlevelsmax to nlevels\n",
            nlevelsmax, nlevels, ion, element);
        nlevelsmax = nlevels;
      }

      // and proceed through the transitionlist till we match this ionstage (if it was not the neutral one)
      int transdata_Z_in = -1;
      int transdata_ionstage_in = -1;
      int tottransitions_in_file = 0;
      std::string line;
      while (transdata_Z_in != Z || transdata_ionstage_in != ionstage) {
        // skip over table
        for (int i = 0; i < tottransitions_in_file; i++) {
          assert_always(getline(ftransitiondata, line));
        }
        assert_always(get_noncommentline(ftransitiondata, line));  // get_noncommentline to skip over blank lines
        assert_always(std::istringstream(line) >> transdata_Z_in >> transdata_ionstage_in >> tottransitions_in_file);
      }

      printout("transdata header matched: transdata_Z_in %d, transdata_ionstage_in %d, tottransitions %d\n",
               transdata_Z_in, transdata_ionstage_in, tottransitions_in_file);
      assert_always(tottransitions_in_file >= 0);

      // read the data for the levels and set up the list of possible transitions for each level
      // store the ions data to memory and set up the ions zeta and levellist
      globals::elements[element].ions[ion].ionstage = ionstage;
      globals::elements[element].ions[ion].nlevels = nlevelsmax;
      globals::elements[element].ions[ion].ionisinglevels = 0;
      globals::elements[element].ions[ion].maxrecombininglevel = -1;
      globals::elements[element].ions[ion].ionpot = ionpot * EV;
      globals::elements[element].ions[ion].nlevels_groundterm = -1;
      globals::elements[element].ions[ion].uniquelevelindexstart = uniquelevelindex;
      globals::elements[element].ions[ion].groundcontindex = -1;
      globals::elements[element].ions[ion].first_nlte = -1;

      globals::elements[element].ions[ion].levels = static_cast<EnergyLevel *>(calloc(nlevelsmax, sizeof(EnergyLevel)));
      assert_always(globals::elements[element].ions[ion].levels != nullptr);

      read_ion_levels(adata, element, ion, nions, nlevels, nlevelsmax, energyoffset, ionpot);

      int tottransitions = tottransitions_in_file;

      if (single_level_top_ion && ion == nions - 1)  // limit the top ion to one level and no transitions
      {
        tottransitions = 0;
      }

      assert_always(transdata_Z_in == Z);
      assert_always(transdata_ionstage_in == ionstage);

      // first nlevels_requiretransitions levels will be collisionally
      // coupled to the first nlevels_requiretransitions_upperlevels levels (assumed forbidden)
      // use 0 to disable adding extra transitions

      const int nlevels_requiretransitions = std::min(nlevelsmax, NLEVELS_REQUIRETRANSITIONS(Z, ionstage));
      // next value with have no effect if nlevels_requiretransitions = 0
      const int nlevels_requiretransitions_upperlevels = nlevelsmax;

      // load transition table for the current ion to temporary memory
      read_ion_transitions(ftransitiondata, tottransitions_in_file, tottransitions, iontransitiontable,
                           nlevels_requiretransitions, nlevels_requiretransitions_upperlevels);

      // last level index is (nlevelsmax - 1), so this is the correct size
      iondowntranstmplineindicies.resize(downtranslevelstart(nlevelsmax));

      add_transitions_to_unsorted_linelist(element, ion, nlevelsmax, iontransitiontable, iondowntranstmplineindicies,
                                           lineindex, temp_linelist, temp_alltranslist);

      for (int level = 0; level < nlevelsmax; level++) {
        uniquelevelindex++;
        totaldowntrans += get_ndowntrans(element, ion, level);
        totaluptrans += get_nuptrans(element, ion, level);
      }

      if (ion < nions - 1) {
        nbfcheck += globals::elements[element].ions[ion].ionisinglevels;  // nlevelsmax;
      }
      uniqueionindex++;
    }
  }
  printout("nbfcheck %d\n", nbfcheck);

  // Save the linecounters value to the global variable containing the number of lines
  globals::nlines = lineindex;
  printout("nlines %d\n", globals::nlines);
  if (globals::rank_in_node == 0) {
    assert_always(globals::nlines == static_cast<int>(temp_linelist.size()));
  }

  if (T_preset > 0) {
    std::abort();
  }

  // Set up the list of allowed upward transitions for each level
  printout("total uptrans %d\n", totaluptrans);
  printout("total downtrans %d\n", totaldowntrans);

  printout("[info] mem_usage: transition lists occupy %.3f MB (this rank) and %.3f MB (shared on node)\n",
           2 * uniquelevelindex * sizeof(LevelTransition *) / 1024. / 1024.,
           (totaluptrans + totaldowntrans) * sizeof(LevelTransition) / 1024. / 1024.);

  if (globals::rank_in_node == 0) {
    // sort the lineline in descending frequency
    std::SORT_OR_STABLE_SORT(EXEC_PAR_UNSEQ temp_linelist.begin(), temp_linelist.end(),
                             [](const auto &a, const auto &b) { return a.nu > b.nu; });

    for (int i = 0; i < globals::nlines - 1; i++) {
      const double nu = temp_linelist[i].nu;
      const double nu_next = temp_linelist[i + 1].nu;
      if (fabs(nu_next - nu) < (1.e-10 * nu)) {
        const auto &a1 = temp_linelist[i];
        const auto &a2 = temp_linelist[i + 1];

        if ((a1.elementindex == a2.elementindex) && (a1.ionindex == a2.ionindex) &&
            (a1.lowerlevelindex == a2.lowerlevelindex) && (a1.upperlevelindex == a2.upperlevelindex)) {
          printout("Duplicate transition line? %s\n", a1.nu == a2.nu ? "nu match exact" : "close to nu match");
          printout("a: Z=%d ionstage %d lower %d upper %d nu %g lambda %g\n", get_atomicnumber(a1.elementindex),
                   get_ionstage(a1.elementindex, a1.ionindex), a1.lowerlevelindex, a1.upperlevelindex, a1.nu,
                   1e8 * CLIGHT / a1.nu);
          printout("b: Z=%d ionstage %d lower %d upper %d nu %g lambda %g\n", get_atomicnumber(a2.elementindex),
                   get_ionstage(a2.elementindex, a2.ionindex), a2.lowerlevelindex, a2.upperlevelindex, a2.nu,
                   1e8 * CLIGHT / a2.nu);
        }
      }
    }
  }

  {
    // create a shared all transitions list and then copy data across, freeing the local copy
    const auto totupdowntrans = totaluptrans + totaldowntrans;
    assert_always(totupdowntrans == static_cast<int>(temp_alltranslist.size()));
#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win win_alltransblock = MPI_WIN_NULL;

    const auto [_, noderank_trans] = get_range_chunk(totupdowntrans, globals::node_nprocs, globals::rank_in_node);

    auto size = static_cast<MPI_Aint>(noderank_trans * sizeof(LevelTransition));
    int disp_unit = sizeof(LevelTransition);
    MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node, &globals::alltrans,
                            &win_alltransblock);

    MPI_Win_shared_query(win_alltransblock, 0, &size, &disp_unit, &globals::alltrans);
#else
    globals::alltrans = static_cast<LevelTransition *>(malloc(totupdowntrans * sizeof(LevelTransition)));
#endif
    if (globals::rank_in_node == 0) {
      std::copy_n(temp_alltranslist.data(), totupdowntrans, globals::alltrans);
      temp_alltranslist.clear();
    }
  }

  // create a linelist shared on node and then copy data across, freeing the local copy
  TransitionLine *nonconstlinelist{};
  {
#ifdef MPI_ON
    MPI_Win win_nonconstlinelist = MPI_WIN_NULL;

    const auto [_, noderank_lines] = get_range_chunk(globals::nlines, globals::node_nprocs, globals::rank_in_node);

    MPI_Aint size = noderank_lines * static_cast<MPI_Aint>(sizeof(TransitionLine));
    int disp_unit = sizeof(TransitionLine);
    MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node, &nonconstlinelist,
                            &win_nonconstlinelist);

    MPI_Win_shared_query(win_nonconstlinelist, 0, &size, &disp_unit, &nonconstlinelist);
#else
    nonconstlinelist = static_cast<TransitionLine *>(malloc(globals::nlines * sizeof(TransitionLine)));
#endif
  }

  if (globals::rank_in_node == 0) {
    memcpy(static_cast<void *>(nonconstlinelist), temp_linelist.data(), globals::nlines * sizeof(TransitionLine));
    temp_linelist.clear();
  }

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  globals::linelist = nonconstlinelist;
  nonconstlinelist = nullptr;
  printout("[info] mem_usage: linelist occupies %.3f MB (node shared memory)\n",
           globals::nlines * sizeof(TransitionLine) / 1024. / 1024);

  // Save sorted linelist into a file
  // if (rank_global == 0)
  // {
  //   FILE *linelist_file = fopen_required("linelist.dat", "w");
  //   fprintf(linelist_file,"%d\n",nlines);
  //   for (int i = 0; i < nlines; i++)
  //   {
  //     fprintf(linelist_file,"%d %d %d %d %d %lg %lg %lg %lg %d\n",
  //             i, globals::linelist[i].elementindex, globals::linelist[i].ionindex,
  //             globals::linelist[i].upperlevelindex, globals::linelist[i].lowerlevelindex,
  //             globals::linelist[i].nu, globals::linelist[i].einstein_A, globals::linelist[i].osc_strength,
  //             globals::linelist[i].coll_str, globals::linelist[i].forbidden);
  //   }
  //   fclose(linelist_file);
  // }

  printout("establishing connection between transitions and sorted linelist...\n");

  auto const time_start_establish_linelist_connections = std::time(nullptr);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (lineindex = 0; lineindex < globals::nlines; lineindex++) {
    if (lineindex % globals::node_nprocs != globals::rank_in_node) {
      continue;
    }

    const auto &line = globals::linelist[lineindex];
    const int element = line.elementindex;
    const int ion = line.ionindex;
    const int lowerlevel = line.lowerlevelindex;
    const int upperlevel = line.upperlevelindex;

    // there is never more than one transition per pair of levels,
    // so find the first matching the upper and lower transition

    const int nupperdowntrans = get_ndowntrans(element, ion, upperlevel);
    auto *downtranslist = get_downtranslist(element, ion, upperlevel);
    auto *downtransition = std::find_if(downtranslist, downtranslist + nupperdowntrans, [=](const auto &downtrans) {
      return downtrans.targetlevelindex == lowerlevel;
    });
    assert_always(downtransition != (downtranslist + nupperdowntrans));
    // assert_always(downtrans->targetlevelindex == lowerlevel);
    downtransition->lineindex = lineindex;

    const int nloweruptrans = get_nuptrans(element, ion, lowerlevel);
    auto *uptranslist = get_uptranslist(element, ion, lowerlevel);
    auto *uptransition = std::find_if(uptranslist, uptranslist + nloweruptrans,
                                      [=](const auto &uptr) { return uptr.targetlevelindex == upperlevel; });
    assert_always(uptransition != (uptranslist + nloweruptrans));
    // assert_always(uptrans->targetlevelindex == upperlevel);
    uptransition->lineindex = lineindex;
  }

  printout("  took %lds\n", std::time(nullptr) - time_start_establish_linelist_connections);
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif

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

  read_phixs_data();

  update_includedionslevels_maxnions();
}

void setup_cellcache() {
  globals::mutex_cellcachemacroatom.resize(get_includedlevels());

  // const int num_cellcache_slots = get_max_threads();
  const int num_cellcache_slots = 1;
  globals::cellcache.resize(num_cellcache_slots);

  for (int cellcachenum = 0; cellcachenum < num_cellcache_slots; cellcachenum++) {
    size_t mem_usage_cellcache = 0;
    mem_usage_cellcache += sizeof(CellCache);

    printout("[info] input: initializing cellcache for thread %d ...\n", cellcachenum);

    globals::cellcache[cellcachenum].cellnumber = -99;

    const auto ncoolingterms = kpkt::ncoolingterms;
    mem_usage_cellcache += ncoolingterms * sizeof(double);
    globals::cellcache[cellcachenum].cooling_contrib = static_cast<double *>(calloc(ncoolingterms, sizeof(double)));

    printout("[info] mem_usage: cellcache coolinglist contribs for thread %d occupies %.3f MB\n", cellcachenum,
             ncoolingterms * sizeof(double) / 1024. / 1024.);

    mem_usage_cellcache += get_nelements() * sizeof(CellCacheElements);
    globals::cellcache[cellcachenum].chelements =
        static_cast<CellCacheElements *>(malloc(get_nelements() * sizeof(CellCacheElements)));

    assert_always(globals::cellcache[cellcachenum].chelements != nullptr);

    ptrdiff_t chlevelcount = 0;
    size_t chphixsblocksize = 0;
    int chtransblocksize = 0;
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        const int nlevels = get_nlevels(element, ion);
        chlevelcount += nlevels;

        for (int level = 0; level < nlevels; level++) {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          chphixsblocksize += nphixstargets * sizeof(CellCachePhixsTargets);

          const int ndowntrans = get_ndowntrans(element, ion, level);
          const int nuptrans = get_nuptrans(element, ion, level);
          chtransblocksize += (2 * ndowntrans + nuptrans);
        }
      }
    }
    assert_always(chlevelcount > 0);
    globals::cellcache[cellcachenum].ch_all_levels.resize(chlevelcount);
    chphixstargetsblock =
        chphixsblocksize > 0 ? static_cast<CellCachePhixsTargets *>(malloc(chphixsblocksize)) : nullptr;
    mem_usage_cellcache += chlevelcount * sizeof(CellCacheLevels) + chphixsblocksize;

    mem_usage_cellcache += chtransblocksize * sizeof(double);
    double *const chtransblock =
        chtransblocksize > 0 ? static_cast<double *>(malloc(chtransblocksize * sizeof(double))) : nullptr;

    int alllevelindex = 0;
    int allphixstargetindex = 0;
    int chtransindex = 0;
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      mem_usage_cellcache += nions * sizeof(CellCacheIons);
      globals::cellcache[cellcachenum].chelements[element].chions =
          static_cast<CellCacheIons *>(malloc(nions * sizeof(CellCacheIons)));
      assert_always(globals::cellcache[cellcachenum].chelements[element].chions != nullptr);

      for (int ion = 0; ion < nions; ion++) {
        const int nlevels = get_nlevels(element, ion);
        auto &chion = globals::cellcache[cellcachenum].chelements[element].chions[ion];
        chion.chlevels = &globals::cellcache[cellcachenum].ch_all_levels[alllevelindex];

        assert_always(alllevelindex == get_uniquelevelindex(element, ion, 0));
        alllevelindex += nlevels;

        for (int level = 0; level < nlevels; level++) {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          chion.chlevels[level].chphixstargets =
              chphixsblocksize > 0 ? &chphixstargetsblock[allphixstargetindex] : nullptr;
          allphixstargetindex += nphixstargets;
        }

        for (int level = 0; level < nlevels; level++) {
          const int ndowntrans = get_ndowntrans(element, ion, level);

          chion.chlevels[level].sum_epstrans_rad_deexc = &chtransblock[chtransindex];
          chtransindex += ndowntrans;
        }

        for (int level = 0; level < nlevels; level++) {
          const int ndowntrans = get_ndowntrans(element, ion, level);
          chion.chlevels[level].sum_internal_down_same = &chtransblock[chtransindex];
          chtransindex += ndowntrans;
        }

        for (int level = 0; level < nlevels; level++) {
          const int nuptrans = get_nuptrans(element, ion, level);
          chion.chlevels[level].sum_internal_up_same = &chtransblock[chtransindex];
          chtransindex += nuptrans;
        }
      }
    }
    assert_always(chtransindex == chtransblocksize);

    assert_always(globals::nbfcontinua >= 0);
    globals::cellcache[cellcachenum].ch_allcont_departureratios.resize(globals::nbfcontinua);
    globals::cellcache[cellcachenum].ch_allcont_nnlevel.resize(globals::nbfcontinua);
    globals::cellcache[cellcachenum].ch_keep_this_cont.resize(globals::nbfcontinua);
    mem_usage_cellcache += 2 * globals::nbfcontinua * sizeof(double);

    printout("[info] mem_usage: cellcache for thread %d occupies %.3f MB\n", cellcachenum,
             mem_usage_cellcache / 1024. / 1024.);
  }
}

void write_bflist_file() {
  globals::bflist.resize(globals::nbfcontinua);

  FILE *bflist_file{};
  if (globals::rank_global == 0) {
    bflist_file = fopen_required("bflist.out", "w");
    fprintf(bflist_file, "%d\n", globals::nbfcontinua);
  }
  int i = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_ionisinglevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        const auto nphixstargets = get_nphixstargets(element, ion, level);
        for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
          const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
          globals::bflist[i].elementindex = element;
          globals::bflist[i].ionindex = ion;
          globals::bflist[i].levelindex = level;
          globals::bflist[i].phixstargetindex = phixstargetindex;

          if (globals::rank_global == 0) {
            fprintf(bflist_file, "%d %d %d %d %d\n", i, element, ion, level, upperionlevel);
          }

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
  if (globals::rank_global == 0) {
    fclose(bflist_file);
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
        globals::elements[element].ions[ion].first_nlte = globals::total_nlte_levels;
        const int nlevels = get_nlevels(element, ion);
        int fullnlteexcitedlevelcount = 0;
        bool found_lte_only_level = false;
        for (int level = 1; level < nlevels; level++) {
          if (is_nlte(element, ion, level)) {
            fullnlteexcitedlevelcount++;
            globals::total_nlte_levels++;
            assert_always(found_lte_only_level == false);  // NLTE levels must be consecutive
          } else {
            found_lte_only_level = true;
          }
        }
        globals::elements[element].ions[ion].nlevels_nlte = fullnlteexcitedlevelcount;

        const bool has_superlevel = (nlevels > (fullnlteexcitedlevelcount + 1));
        if (has_superlevel) {
          // If there are more levels that the ground state + the number of NLTE levels then we need an extra
          // slot to store data for the "superlevel", which is a representation of all the other levels that
          // are not treated in detail.
          globals::total_nlte_levels++;
          n_super_levels++;
        }

        assert_always(has_superlevel == ion_has_superlevel(element, ion));

        printout("[input]  element %2d Z=%2d ionstage %2d has %5d NLTE excited levels%s. Starting at %d\n", element,
                 get_atomicnumber(element), get_ionstage(element, ion), fullnlteexcitedlevelcount,
                 has_superlevel ? " plus a superlevel" : "", globals::elements[element].ions[ion].first_nlte);
      }
    }
  }

  printout("[input] Total NLTE levels: %d, of which %d are superlevels\n", globals::total_nlte_levels, n_super_levels);
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
  printout("[input] this simulation contains\n");
  printout("----------------------------------\n");
  for (int element = 0; element < get_nelements(); element++) {
    printout("[input]  element %d (Z=%2d %s)\n", element, get_atomicnumber(element),
             decay::get_elname(get_atomicnumber(element)).c_str());
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      int ion_photoiontransitions = 0;
      int ion_bbtransitions = 0;
      for (int level = 0; level < get_nlevels(element, ion); level++) {
        ion_photoiontransitions += get_nphixstargets(element, ion, level);
        ion_bbtransitions += get_nuptrans(element, ion, level);
      }

      printout(
          "[input]    ionstage %d: %4d levels (%d in groundterm, %4d ionising) %7d lines %6d bf transitions "
          "(epsilon_ground: %7.2f eV)\n",
          get_ionstage(element, ion), get_nlevels(element, ion), get_nlevels_groundterm(element, ion),
          get_ionisinglevels(element, ion), ion_bbtransitions, ion_photoiontransitions, epsilon(element, ion, 0) / EV);

      includedionisinglevels += get_ionisinglevels(element, ion);
      includedphotoiontransitions += ion_photoiontransitions;
      includedboundboundtransitions += ion_bbtransitions;
    }
  }
  assert_always(includedphotoiontransitions == globals::nbfcontinua);
  assert_always(globals::nlines == includedboundboundtransitions);

  printout("[input]  in total %d ions, %d levels (%d ionising), %d lines, %d photoionisation transitions\n",
           get_includedions(), get_includedlevels(), includedionisinglevels, globals::nlines, globals::nbfcontinua);

  write_bflist_file();

  setup_nlte_levels();
}

}  // anonymous namespace

// read input.txt, atomic data, and ejecta model
void input(int rank) {
  globals::n_titer = 1;
  globals::lte_iteration = false;

  printout("[info] input: do n_titer %d iterations per timestep\n", globals::n_titer);
  if (globals::n_titer > 1) {
#ifndef DO_TITER
    printout("[fatal] input: n_titer > 1, but DO_TITER not defined ... abort\n");
    std::abort();
#endif
  } else if (globals::n_titer == 1) {
#ifdef DO_TITER
    printout("[warning] input: n_titer = 1 but DO_TITER defined, remove DO_TITER to save memory\n");
#endif
  } else {
    printout("[fatal] input: no valid value for n_titer selected\n");
    std::abort();
  }

  // Read in parameters from input.txt
  read_parameterfile(rank);

  // Read in parameters from vpkt.txt
  if (VPKT_ON) {
    read_parameterfile_vpkt();
  }

  read_atomicdata();

#ifdef MPI_ON
  const auto time_before_barrier = std::time(nullptr);
  printout("barrier after read_atomicdata(): time before barrier %d, ", static_cast<int>(time_before_barrier));
  MPI_Barrier(MPI_COMM_WORLD);
  printout("time after barrier %d (waited %d seconds)\n", static_cast<int>(time(nullptr)),
           static_cast<int>(time(nullptr) - time_before_barrier));
#endif

  grid::read_ejecta_model();
}

// read the next line, skipping any comment lines beginning with '#'
auto get_noncommentline(std::fstream &input, std::string &line) -> bool {
  while (true) {
    const bool linefound = !(!std::getline(input, line));
    // printout("LINE: >%s<  linefound: %s commentonly: %s \n", line.c_str(), linefound ? "true" : "false",
    // lineiscommentonly(line) ? "true" : "false");
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
  std::istringstream(line) >> pre_zseed;

  if (pre_zseed > 0) {
    printout("using input.txt specified random number seed of %" PRId64 "\n", pre_zseed);
  } else {
#ifndef GPU_ON
    pre_zseed = std::random_device{}();
#endif
#ifdef MPI_ON
    // broadcast randomly-generated seed from rank 0 to all ranks
    MPI_Bcast(&pre_zseed, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
#endif
    printout("randomly-generated random number seed is %" PRId64 "\n", pre_zseed);
#if defined REPRODUCIBLE && REPRODUCIBLE
    printout("ERROR: reproducible mode is on, so random number seed is required.\n");
    std::abort();
#endif
  }

#if defined(_OPENMP) && !defined(GPU_ON)
#pragma omp parallel
#endif
  {
    // For MPI parallelisation, the random seed is changed based on the rank of the process
    // For OpenMP parallelisation rng is a threadprivate variable and the seed changed according
    // to the thread-ID tid.
    const auto tid = get_thread_num();
    auto rngseed = pre_zseed + static_cast<std::int64_t>(13 * (rank * get_max_threads() + tid));
#ifndef GPU_ON
    stdrng.seed(rngseed);
#endif
    printout("rank %d: thread %d has rngseed %" PRId64 "\n", rank, tid, rngseed);
    printout("rng is a std::mt19937 generator\n");

    // call it a few times
    for (int n = 0; n < 100; n++) {
      rng_uniform();
    }
  }

  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> globals::ntimesteps;  // number of time steps
  assert_always(globals::ntimesteps > 0);

  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> globals::timestep_initial >>
      globals::timestep_finish;  // number of start and end time step
  printout("input: timestep_start %d timestep_finish %d\n", globals::timestep_initial, globals::timestep_finish);
  assert_always(globals::timestep_initial < globals::ntimesteps);
  assert_always(globals::timestep_initial <= globals::timestep_finish);
  assert_always(globals::timestep_finish <= globals::ntimesteps);

  double tmin_days = 0.;
  double tmax_days = 0.;
  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> tmin_days >> tmax_days;  // start and end times
  assert_always(tmin_days > 0);
  assert_always(tmax_days > 0);
  assert_always(tmin_days < tmax_days);
  globals::tmin = tmin_days * DAY;
  globals::tmax = tmax_days * DAY;

  assert_always(get_noncommentline(file, line));  // UNUSED nusyn_min_mev nusyn_max_mev

  assert_always(get_noncommentline(file, line));  // UNUSED number of times for synthesis

  assert_always(get_noncommentline(file, line));  // UNUSED start and end times for synthesis

  assert_always(get_noncommentline(file, line));  // model dimensions
  int dum1 = 0;
  std::istringstream(line) >> dum1;
  if (dum1 == 1) {
    grid::set_model_type(GridType::SPHERICAL1D);
  } else if (dum1 == 2) {
    grid::set_model_type(GridType::CYLINDRICAL2D);
  } else if (dum1 == 3) {
    grid::set_model_type(GridType::CARTESIAN3D);
  }

  assert_always(get_noncommentline(file, line));  // UNUSED compute the r-light curve?

  assert_always(get_noncommentline(file, line));  // UNUSED number of iterations

  assert_always(get_noncommentline(file, line));  // UNUSED change speed of light

  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> globals::gamma_kappagrey;  // use grey opacity for gammas?

  std::array<float, 3> syn_dir_in{-1.};
  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> syn_dir_in[0] >> syn_dir_in[1] >> syn_dir_in[2];  // components of syn_dir

  const double rr = (syn_dir_in[0] * syn_dir_in[0]) + (syn_dir_in[1] * syn_dir_in[1]) + (syn_dir_in[2] * syn_dir_in[2]);
  // ensure that this vector is normalised.
  if (rr > 1.e-6) {
    globals::syn_dir = {syn_dir_in[0] / sqrt(rr), syn_dir_in[1] / sqrt(rr), syn_dir_in[2] / sqrt(rr)};
  } else {
    const double z1 = 1. - (2 * rng_uniform());
    const double z2 = rng_uniform() * 2.0 * PI;
    globals::syn_dir = {sqrt((1. - (z1 * z1))) * cos(z2), sqrt((1. - (z1 * z1))) * sin(z2), z1};
  }

  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> globals::opacity_case;  // opacity choice

  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> globals::rho_crit_para;  // free parameter for calculation of rho_crit
  printout("input: rho_crit_para %g\n", globals::rho_crit_para);
  // he calculation of rho_crit itself depends on the time, therefore it happens in grid_init and update_grid

  assert_always(get_noncommentline(file, line));
  int debug_packet = 0;
  std::istringstream(line) >> debug_packet;  // activate debug output for packet
  assert_always(debug_packet == -1);
  // select a negative value to deactivate

  // Do we start a new simulation or, continue another one?
  int continue_flag = 0;
  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> continue_flag;
  globals::simulation_continued_from_saved = (continue_flag == 1);
  if (globals::simulation_continued_from_saved) {
    printout("input: resuming simulation from saved point\n");
  } else {
    printout("input: starting a new simulation\n");
    assert_always(globals::timestep_initial == 0);
  }

  // Wavelength (in Angstroms) at which the parameterisation of the radiation field
  // switches from the nebular approximation to LTE.
  float dum2{NAN};
  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> dum2;  // free parameter for calculation of rho_crit
  globals::nu_rfcut = CLIGHT / (dum2 * 1e-8);
  printout("input: nu_rfcut %g\n", globals::nu_rfcut);

  // Sets the number of initial LTE timesteps for NLTE runs
  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> globals::num_lte_timesteps;
  printout("input: doing the first %d timesteps in LTE\n", globals::num_lte_timesteps);

  if (NT_ON) {
    if (NT_SOLVE_SPENCERFANO) {
      printout("input: Non-thermal ionisation with a Spencer-Fano solution is switched on for this run.\n");
    } else {
      printout("input: Non-thermal ionisation with the work function approximation is switched on for this run.\n");
    }
  } else {
    printout("input: No non-thermal ionisation is used in this run.\n");
  }

  if (USE_LUT_PHOTOION) {
    printout(
        "Corrphotoioncoeff is calculated from LTE lookup tables (ratecoeff.dat) and corrphotoionrenorm "
        "estimator.\n");
  } else {
    printout(
        "Corrphotoioncoeff is calculated from the radiation field at each timestep in each modelgrid cell (no "
        "LUT).\n");
  }

  if (USE_LUT_BFHEATING) {
    printout("bfheating coefficients are calculated from LTE lookup tables (ratecoeff.dat) and bfheatingestimator.\n");
  } else {
    printout(
        "bfheating coefficients are calculated from the radiation field at each timestep in each modelgrid cell (no "
        "LUT).\n");
  }

  // Set up initial grey approximation?
  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> globals::cell_is_optically_thick >> globals::num_grey_timesteps;
  printout(
      "input: cells with Thomson optical depth > %g are treated in grey approximation for the first %d timesteps\n",
      globals::cell_is_optically_thick, globals::num_grey_timesteps);

  // Limit the number of bf-continua
  assert_always(get_noncommentline(file, line));
  int max_bf_continua = 0;
  std::istringstream(line) >> max_bf_continua;
  assert_always(max_bf_continua == -1);

  // for exspec: read number of MPI tasks
  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> globals::nprocs_exspec;

  // Extract line-of-sight dependent information of last emission for spectrum_res
  assert_always(get_noncommentline(file, line));
  std::istringstream(line) >> dum1;
  globals::do_emission_res = (dum1 != 0);

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
  std::istringstream(line) >> kpktdiffusion_timescale >> n_kpktdiffusion_timesteps;
  kpkt::set_kpktdiffusion(kpktdiffusion_timescale, n_kpktdiffusion_timesteps);

  file.close();

  if (rank == 0 && !globals::simulation_continued_from_saved) {
    // back up original input file, adding comments to each line
    update_parameterfile(-1);
  }
}

void update_parameterfile(int nts)
// Subroutine to read in input parameters from input.txt.
{
  assert_always(globals::rank_global == 0);
  if (nts >= 0) {
    printout("Update input.txt for restart at timestep %d...", nts);
  } else {
    printout("Copying input.txt to input-newrun.txt...");
  }

  std::ifstream file("input.txt");
  assert_always(file.is_open());

  std::ofstream fileout("input.txt.tmp");
  assert_always(fileout.is_open());

  std::string line;

  // FILE *input_file = fopen_required("input.txt", "r+");
  // setvbuf(input_file,nullptr, _IOLBF, 0);

  char c_line[1024];
  int noncomment_linenum = -1;
  while (std::getline(file, line)) {
    if (!lineiscommentonly(line)) {
      noncomment_linenum++;  // line number starting from 0, ignoring comment and blank lines (that start with '#')

      // overwrite particular lines to enable restarting from the current timestep
      if (nts >= 0) {
        if (noncomment_linenum == 2) {
          // Number of start and end time step
          snprintf(c_line, sizeof(c_line), "%d %d", nts, globals::timestep_finish);
          // line.assign(c_line);
          line.replace(line.begin(), line.end(), c_line);
        } else if (noncomment_linenum == 16) {
          // resume from gridsave file
          snprintf(c_line, sizeof(c_line), "%d", 1);  // Force continuation
          line.assign(c_line);
        }
      }

      if (noncomment_linenum == 21) {
        // by default, exspec should use all available packet files
        globals::nprocs_exspec = globals::nprocs;
        snprintf(c_line, sizeof(c_line), "%d", globals::nprocs_exspec);
        line.assign(c_line);
      }

      if (noncomment_linenum < static_cast<int>(inputlinecomments.size())) {
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

  printout("done\n");
}

// initialize the time steps
void time_init() {
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
        globals::timesteps[n].start = globals::tmin + n * dt;
        globals::timesteps[n].width = dt;
        globals::timesteps[n].mid = globals::timesteps[n].start + 0.5 * globals::timesteps[n].width;
      }
      break;
    }

    case TimeStepSizeMethod::LOGARITHMIC_THEN_CONSTANT: {
      // First part log, second part fixed timesteps
      const double t_transition = TIMESTEP_TRANSITION_TIME * DAY;  // transition from logarithmic to fixed timesteps
      const double maxtsdelta = FIXED_TIMESTEP_WIDTH * DAY;        // maximum timestep width in fixed part
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
          globals::timesteps[n].mid = globals::timesteps[n].start + 0.5 * globals::timesteps[n].width;
        }
      }
      break;
    }

    case TimeStepSizeMethod::CONSTANT_THEN_LOGARITHMIC: {
      // // First part fixed timesteps, second part log timesteps
      const double t_transition = TIMESTEP_TRANSITION_TIME * DAY;  // transition from fixed to logarithmic timesteps
      const double maxtsdelta = FIXED_TIMESTEP_WIDTH * DAY;        // timestep width of fixed timesteps
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
          globals::timesteps[n].start = globals::tmin + n * fixed_tsdelta;
          globals::timesteps[n].width = fixed_tsdelta;
          globals::timesteps[n].mid = globals::timesteps[n].start + 0.5 * globals::timesteps[n].width;
        } else {
          // For logarithmic time steps, the logarithmic interval will be
          const double dlogt = (log(globals::tmax) - log(t_transition)) / nts_log;
          const double prev_start =
              n > 0 ? (globals::timesteps[n - 1].start + globals::timesteps[n - 1].width) : globals::tmin;
          globals::timesteps[n].start = prev_start;
          globals::timesteps[n].width = (t_transition * exp((n - nts_fixed + 1) * dlogt)) - globals::timesteps[n].start;
          globals::timesteps[n].mid = globals::timesteps[n].start + 0.5 * globals::timesteps[n].width;
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
  FILE *timestepfile = fopen_required("timesteps.out", "w");
  fprintf(timestepfile, "#timestep tstart_days tmid_days twidth_days\n");
  for (int n = 0; n < globals::ntimesteps; n++) {
    fprintf(timestepfile, "%d %lg %lg %lg\n", n, globals::timesteps[n].start / DAY, globals::timesteps[n].mid / DAY,
            globals::timesteps[n].width / DAY);
  }
  fclose(timestepfile);
}
