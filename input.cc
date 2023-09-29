#include "input.h"

#include <gsl/gsl_spline.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <ranges>
#include <sstream>
#include <vector>

#include "atomic.h"
#include "decay.h"
#include "exspec.h"
#include "gammapkt.h"
#include "grid.h"
#include "kpkt.h"
#include "nltepop.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "sn3d.h"
#include "vpkt.h"

const int groundstate_index_in = 1;  // starting level index in the input files
int phixs_file_version = -1;

struct transitions {
  int *to;
};

struct transitiontable_entry {
  int lower;
  int upper;
  float A;
  float coll_str;
  bool forbidden;
};  /// only used temporarily during input

constexpr std::array<std::string_view, 24> inputlinecomments = {
    " 0: pre_zseed: specific random number seed if > 0 or random if negative",
    " 1: globals::ntstep: number of timesteps",
    " 2: itstep ftstep: timestep number range start (inclusive) and stop (not inclusive)",
    " 3: tmin_days tmax_days: start and end times [day]",
    " 4: UNUSED nusyn_min_mev nusyn_max_mev: lowest and highest frequency to synthesise [MeV]",
    " 5: UNUSED nsyn_time: number of times for synthesis",
    " 6: UNUSED start and end times for synthesis",
    " 7: model_type: number of dimensions (1, 2, or 3)",
    " 8: UNUSED compute r-light curve (1: no estimators, 2: thin cells, 3: thick cells, 4: gamma-ray heating)",
    " 9: UNUSED n_out_it: number of iterations",
    "10: UNUSED: change speed of light by some factor. Change constants.h CLIGHT_PROP instead",
    "11: use grey opacity for gammas?",
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

static void read_phixs_data_table(std::fstream &phixsfile, const int nphixspoints_inputtable, const int element,
                                  const int lowerion, const int lowerlevel, const int upperion, int upperlevel_in,
                                  const double phixs_threshold_ev, size_t *mem_usage_phixs) {
  if (upperlevel_in >= 0)  // file gives photoionisation to a single target state only
  {
    int upperlevel = upperlevel_in - groundstate_index_in;
    assert_always(upperlevel >= 0);
    assert_always(globals::elements[element].ions[lowerion].levels[lowerlevel].nphixstargets == 0);
    globals::elements[element].ions[lowerion].levels[lowerlevel].nphixstargets = 1;
    *mem_usage_phixs += sizeof(phixstarget_entry);

    assert_always(globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets == nullptr);
    globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets =
        static_cast<phixstarget_entry *>(calloc(1, sizeof(phixstarget_entry)));
    assert_always(globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets != nullptr);

    if (single_level_top_ion &&
        (upperion == get_nions(element) - 1))  // top ion has only one level, so send it to that level
    {
      upperlevel = 0;
    }
    globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets[0].levelindex = upperlevel;
    globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets[0].probability = 1.0;
  } else  // upperlevel < 0, indicating that a table of upper levels and their probabilities will follow
  {
    int in_nphixstargets = 0;
    assert_always(phixsfile >> in_nphixstargets);
    assert_always(in_nphixstargets >= 0);
    // read in a table of target states and probabilities and store them
    if (!single_level_top_ion || upperion < get_nions(element) - 1)  // in case the top ion has nlevelsmax = 1
    {
      globals::elements[element].ions[lowerion].levels[lowerlevel].nphixstargets = in_nphixstargets;
      *mem_usage_phixs += in_nphixstargets * sizeof(phixstarget_entry);

      globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets =
          static_cast<phixstarget_entry *>(calloc(in_nphixstargets, sizeof(phixstarget_entry)));
      assert_always(globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets != nullptr);

      double probability_sum = 0.;
      for (int i = 0; i < in_nphixstargets; i++) {
        double phixstargetprobability = NAN;
        assert_always(phixsfile >> upperlevel_in >> phixstargetprobability);
        const int upperlevel = upperlevel_in - groundstate_index_in;
        assert_always(upperlevel >= 0);
        assert_always(phixstargetprobability > 0);
        globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets[i].levelindex = upperlevel;
        globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets[i].probability =
            phixstargetprobability;
        probability_sum += phixstargetprobability;
      }
      if (fabs(probability_sum - 1.0) > 0.01) {
        printout("WARNING: photoionisation table for Z=%d ionstage %d has probabilities that sum to %g",
                 get_atomicnumber(element), get_ionstage(element, lowerion), probability_sum);
      }
    } else  // file has table of target states and probabilities but our top ion is limited to one level
    {
      globals::elements[element].ions[lowerion].levels[lowerlevel].nphixstargets = 1;
      *mem_usage_phixs += sizeof(phixstarget_entry);
      globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets =
          static_cast<phixstarget_entry *>(calloc(1, sizeof(phixstarget_entry)));
      assert_always(globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets != nullptr);

      for (int i = 0; i < in_nphixstargets; i++) {
        double phixstargetprobability = NAN;
        assert_always(phixsfile >> upperlevel_in >> phixstargetprobability);
      }

      // send it to the ground state of the top ion
      globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets[0].levelindex = 0;
      globals::elements[element].ions[lowerion].levels[lowerlevel].phixstargets[0].probability = 1.0;
    }
  }

  /// The level contributes to the ionisinglevels if its energy
  /// is below the ionisation potential and the level doesn't
  /// belong to the topmost ion included.
  /// Rate coefficients are only available for ionising levels.
  //  also need (levelenergy < ionpot && ...)?
  if (lowerion < get_nions(element) - 1)  /// thats only an option for pure LTE && level < TAKE_N_BFCONTINUA)
  {
    for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, lowerion, lowerlevel);
         phixstargetindex++) {
      const int upperlevel = get_phixsupperlevel(element, lowerion, lowerlevel, phixstargetindex);
      if (upperlevel > get_maxrecombininglevel(element, lowerion + 1)) {
        globals::elements[element].ions[lowerion + 1].maxrecombininglevel = upperlevel;
      }
    }
  }

  *mem_usage_phixs += globals::NPHIXSPOINTS * sizeof(float);
  globals::elements[element].ions[lowerion].levels[lowerlevel].photoion_xs =
      static_cast<float *>(calloc(globals::NPHIXSPOINTS, sizeof(float)));
  assert_always(globals::elements[element].ions[lowerion].levels[lowerlevel].photoion_xs != nullptr);

  if (phixs_threshold_ev > 0) {
    globals::elements[element].ions[lowerion].levels[lowerlevel].phixs_threshold = phixs_threshold_ev * EV;
  } else {
    if (get_nphixstargets(element, lowerion, lowerlevel) > 0) {
      const int lowestupperlevel = get_phixsupperlevel(element, lowerion, lowerlevel, 0);
      const double calced_phixs_threshold =
          (epsilon(element, upperion, lowestupperlevel) - epsilon(element, lowerion, lowerlevel));
      globals::elements[element].ions[lowerion].levels[lowerlevel].phixs_threshold = calced_phixs_threshold;
    }
  }

  if (phixs_file_version == 1) {
    assert_always(get_nphixstargets(element, lowerion, lowerlevel) == 1);
    assert_always(get_phixsupperlevel(element, lowerion, lowerlevel, 0) == 0);

    double const nu_edge = (epsilon(element, upperion, 0) - epsilon(element, lowerion, lowerlevel)) / H;

    auto *nutable = static_cast<double *>(calloc(nphixspoints_inputtable, sizeof(double)));
    assert_always(nutable != nullptr);
    auto *phixstable = static_cast<double *>(calloc(nphixspoints_inputtable, sizeof(double)));
    assert_always(phixstable != nullptr);

    for (int i = 0; i < nphixspoints_inputtable; i++) {
      double energy = -1.;
      double phixs = -1.;
      assert_always(phixsfile >> energy >> phixs);
      nutable[i] = nu_edge + (energy * 13.6 * EV) / H;
      /// the photoionisation cross-sections in the database are given in Mbarn=1e6 * 1e-28m^2
      /// to convert to cgs units multiply by 1e-18
      phixstable[i] = phixs * 1e-18;
    }
    const double nu_max = nutable[nphixspoints_inputtable - 1];

    // Now interpolate these cross-sections
    globals::elements[element].ions[lowerion].levels[lowerlevel].photoion_xs[0] = phixstable[0];

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear, nphixspoints_inputtable);
    gsl_spline_init(spline, nutable, phixstable, nphixspoints_inputtable);
    for (int i = 1; i < globals::NPHIXSPOINTS; i++) {
      const double nu = nu_edge * (1. + i * globals::NPHIXSNUINCREMENT);
      if (nu > nu_max) {
        const double phixs = phixstable[nphixspoints_inputtable - 1] * pow(nu_max / nu, 3);
        globals::elements[element].ions[lowerion].levels[lowerlevel].photoion_xs[i] = phixs;
      } else {
        const double phixs = gsl_spline_eval(spline, nu, acc);
        globals::elements[element].ions[lowerion].levels[lowerlevel].photoion_xs[i] = phixs;
      }
    }
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    free(nutable);
    free(phixstable);
  } else {
    for (int i = 0; i < globals::NPHIXSPOINTS; i++) {
      float phixs = NAN;
      assert_always(phixsfile >> phixs);
      assert_always(phixs >= 0);

      /// the photoionisation cross-sections in the database are given in Mbarn = 1e6 * 1e-28m^2
      /// to convert to cgs units multiply by 1e-18
      globals::elements[element].ions[lowerion].levels[lowerlevel].photoion_xs[i] = phixs * 1e-18;
      // fprintf(database_file,"%g %g\n", nutable[i], phixstable[i]);
    }
  }

  // nbfcontinua++;
  // printout("[debug] element %d, ion %d, level %d: phixs exists %g\n",element,lowerion,lowerlevel,phixs*1e-18);
  globals::nbfcontinua += get_nphixstargets(element, lowerion, lowerlevel);
  if (lowerlevel < get_nlevels_groundterm(element, lowerion)) {
    globals::nbfcontinua_ground += get_nphixstargets(element, lowerion, lowerlevel);
  }
}

static void read_phixs_data(const int phixs_file_version) {
  size_t mem_usage_phixs = 0;

  printout("readin phixs data from %s\n", phixsdata_filenames[phixs_file_version]);

  auto phixsfile = fstream_required(phixsdata_filenames[phixs_file_version], std::ios::in);

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
  double phixs_threshold_ev = -1;
  while (true) {
    int nphixspoints_inputtable = 0;
    if (phixs_file_version == 1) {
      if (!(phixsfile >> Z >> upperionstage >> upperlevel_in >> lowerionstage >> lowerlevel_in >>
            nphixspoints_inputtable)) {
        break;
      }
    } else {
      if (!(phixsfile >> Z >> upperionstage >> upperlevel_in >> lowerionstage >> lowerlevel_in >> phixs_threshold_ev)) {
        break;
      };
      nphixspoints_inputtable = globals::NPHIXSPOINTS;
    }
    assert_always(Z > 0);
    assert_always(upperionstage >= 2);
    assert_always(lowerionstage >= 1);
    bool skip_this_phixs_table = false;
    // printout("[debug] Z %d, upperion %d, upperlevel %d, lowerion %d, lowerlevel,
    // %d\n",Z,upperion,upperlevel,lowerion,lowerlevel);
    /// translate readin anumber to element index
    const int element = get_elementindex(Z);

    /// store only photoionization crosssections for elements that are part of the current model atom
    skip_this_phixs_table = true;  // will be set to false for good data
    if (element >= 0 && get_nions(element) > 0) {
      /// translate readin ionstages to ion indices

      const int upperion = upperionstage - get_ionstage(element, 0);
      const int lowerion = lowerionstage - get_ionstage(element, 0);
      const int lowerlevel = lowerlevel_in - groundstate_index_in;
      assert_always(lowerionstage >= 0);
      assert_always(lowerlevel >= 0);
      /// store only photoionization crosssections for ions that are part of the current model atom
      if (lowerion >= 0 && upperion < get_nions(element) && lowerlevel < get_nlevels(element, lowerion)) {
        read_phixs_data_table(phixsfile, nphixspoints_inputtable, element, lowerion, lowerlevel, upperion,
                              upperlevel_in, phixs_threshold_ev, &mem_usage_phixs);

        skip_this_phixs_table = false;
      }
    }

    if (skip_this_phixs_table)  // for ions or elements that are not part of the current model atom, proceed through the
                                // lines and throw away the data
    {
      if (upperlevel_in < 0)  // a table of target states and probabilities will follow, so read past those lines
      {
        int nphixstargets = 0;
        assert_always(phixsfile >> nphixstargets);
        for (int i = 0; i < nphixstargets; i++) {
          double phixstargetprobability = NAN;
          assert_always(phixsfile >> upperlevel_in >> phixstargetprobability);
        }
      }
      for (int i = 0; i < nphixspoints_inputtable; i++)  // skip through cross section list
      {
        float phixs = 0;
        if (phixs_file_version == 1) {
          double energy = 0;
          assert_always(phixsfile >> energy >> phixs);
        } else {
          assert_always(phixsfile >> phixs);
        }
      }
    }
  }

  printout("[info] mem_usage: photoionisation tables occupy %.3f MB\n", mem_usage_phixs / 1024. / 1024.);
}

static void read_ion_levels(FILE *adata, const int element, const int ion, const int nions, const int nlevels,
                            int nlevelsmax, const double energyoffset, const double ionpot,
                            struct transitions *transitions) {
  const int nlevels_used = std::min(nlevels, nlevelsmax);
  // each level contains 0..level elements. seriess sum of 1 + 2 + 3 + 4 + ... + nlevels_used is used here
  const int transitblocksize = nlevels_used * (nlevels_used + 1) / 2;
  transitions[0].to = static_cast<int *>(malloc(transitblocksize * sizeof(int)));

  int transitionblockindex = 0;
  for (int level = 0; level < nlevels; level++) {
    int levelindex_in = 0;
    double levelenergy = NAN;
    double statweight = NAN;
    int ntransitions = 0;
    assert_always(fscanf(adata, "%d %lg %lg %d%*[^\n]\n", &levelindex_in, &levelenergy, &statweight, &ntransitions) ==
                  4);
    assert_always(levelindex_in == level + groundstate_index_in);

    if (level < nlevelsmax) {
      const double currentlevelenergy = (energyoffset + levelenergy) * EV;
      globals::elements[element].ions[ion].levels[level].epsilon = currentlevelenergy;

      // if (level == 0 && ion == 0) energyoffset = levelenergy;
      globals::elements[element].ions[ion].levels[level].stat_weight = statweight;
      assert_always(statweight > 0.);
      /// Moved to the section with ionising levels below
      // globals::elements[element].ions[ion].levels[level].cont_index = cont_index;
      // cont_index--;
      /// Initialise the metastable flag to true. Set it to false if a downward transition exists.
      globals::elements[element].ions[ion].levels[level].metastable = true;
      // globals::elements[element].ions[ion].levels[level].main_qn = mainqn;

      /// The level contributes to the ionisinglevels if its energy
      /// is below the ionization potential and the level doesn't
      /// belong to the topmost ion included.
      /// Rate coefficients are only available for ionising levels.
      if (levelenergy < ionpot && ion < nions - 1)  /// thats only an option for pure LTE && level < TAKE_N_BFCONTINUA)
      {
        globals::elements[element].ions[ion].ionisinglevels++;
      }

      /// store the possible downward transitions from the current level in following order to memory
      ///     A_level,level-1; A_level,level-2; ... A_level,1
      /// entries which are not explicitly set are zero (the zero is set/initialized by calloc!)
      transitions[level].to = &transitions[0].to[transitionblockindex];
      transitionblockindex += level;
      assert_always((transitionblockindex + level) < transitblocksize);

      /// initialize number of downward transitions to zero
      set_ndowntrans(element, ion, level, 0);

      /// initialize number of upward transitions to zero
      set_nuptrans(element, ion, level, 0);
    } else {
      // globals::elements[element].ions[ion].levels[nlevelsmax - 1].stat_weight += statweight;
    }
  }
}

static void read_ion_transitions(std::fstream &ftransitiondata, const int tottransitions_in_file, int *tottransitions,
                                 std::vector<struct transitiontable_entry> &transitiontable,
                                 const int nlevels_requiretransitions,
                                 const int nlevels_requiretransitions_upperlevels) {
  std::string line;

  if (*tottransitions == 0) {
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
        std::stringstream ss(line);
        std::string word;
        int word_count = 0;
        while (ss >> word) {
          word_count++;
        }
        assert_always(word_count == 4 || word_count == 5);
        oldtransitionformat = (word_count == 4);
      }
      if (!oldtransitionformat) {
        assert_always(sscanf(line.c_str(), "%d %d %g %g %d", &lower_in, &upper_in, &A, &coll_str, &intforbidden) == 5);
      } else {
        int transindex = 0;  // not used
        assert_always(sscanf(line.c_str(), "%d %d %d %g", &transindex, &lower_in, &upper_in, &A) == 4);
      }
      const int lower = lower_in - groundstate_index_in;
      const int upper = upper_in - groundstate_index_in;
      assert_always(lower >= 0);
      assert_always(upper >= 0);

      // this entire block can be removed if we don't want to add in extra collisonal
      // transitions between levels
      if (prev_lower < nlevels_requiretransitions) {
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
          (*tottransitions)++;
          assert_always(prev_lower >= 0);
          assert_always(tmplevel >= 0);
          transitiontable.push_back(
              {.lower = prev_lower, .upper = tmplevel, .A = 0., .coll_str = -2., .forbidden = true});
        }
      }

      transitiontable.push_back(
          {.lower = lower, .upper = upper, .A = A, .coll_str = coll_str, .forbidden = (intforbidden == 1)});
      // printout("index %d, lower %d, upper %d, A %g\n",transitionindex,lower,upper,A);
      //  printout("reading transition index %d lower %d upper %d\n", i, transitiontable[i].lower,
      //  transitiontable[i].upper);
      prev_lower = lower;
      prev_upper = upper;
    }
  }
}

static void add_transitions_to_unsorted_linelist(const int element, const int ion, const int nlevelsmax,
                                                 const std::vector<struct transitiontable_entry> &transitiontable,
                                                 struct transitions *transitions, int *lineindex,
                                                 std::vector<struct linelist_entry> &temp_linelist) {
  const int lineindex_initial = *lineindex;
  const int tottransitions = transitiontable.size();
  int totupdowntrans = 0;
  // pass 0 to get transition counts of each level
  // pass 1 to allocate and fill transition arrays
  for (int pass = 0; pass < 2; pass++) {
    *lineindex = lineindex_initial;
    if (pass == 1) {
      int alltransindex = 0;
      struct level_transition *alltransblock = nullptr;

#ifdef MPI_ON
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Win win = MPI_WIN_NULL;

      int my_rank_trans = totupdowntrans / globals::node_nprocs;
      // rank_in_node 0 gets any remainder
      if (globals::rank_in_node == 0) {
        my_rank_trans += totupdowntrans - (my_rank_trans * globals::node_nprocs);
      }

      MPI_Aint size = my_rank_trans * sizeof(linelist_entry);
      int disp_unit = sizeof(linelist_entry);
      MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node, &alltransblock, &win);

      MPI_Win_shared_query(win, 0, &size, &disp_unit, &alltransblock);
#else
      alltransblock = static_cast<struct level_transition *>(malloc(totupdowntrans * sizeof(struct level_transition)));
#endif

      for (int level = 0; level < nlevelsmax; level++) {
        globals::elements[element].ions[ion].levels[level].downtrans = &alltransblock[alltransindex];
        alltransindex += get_ndowntrans(element, ion, level);

        globals::elements[element].ions[ion].levels[level].uptrans = &alltransblock[alltransindex];
        alltransindex += get_nuptrans(element, ion, level);

        set_ndowntrans(element, ion, level, 0);
        set_nuptrans(element, ion, level, 0);
      }
    }

    for (int level = 0; level < nlevelsmax; level++) {
      std::fill_n(transitions[level].to, level, -99);
    }

    totupdowntrans = 0;
    for (int ii = 0; ii < tottransitions; ii++) {
      const int level = transitiontable[ii].upper;
      const int targetlevel = transitiontable[ii].lower;
      if (pass == 0) {
        assert_always(targetlevel >= 0);
        assert_always(level > targetlevel);
      }

      double nu_trans = -1.;
      if (targetlevel < nlevelsmax && level < nlevelsmax) {
        nu_trans = (epsilon(element, ion, level) - epsilon(element, ion, targetlevel)) / H;
      }
      if (!(nu_trans > 0)) {
        continue;
      }

      /// Make sure that we don't allow duplicate. In that case take only the lines
      /// first occurrence
      const int transitioncheck = transitions[level].to[(level - targetlevel) - 1];

      // -99 means that the transition hasn't been seen yet
      if (transitioncheck == -99) {
        transitions[level].to[level - targetlevel - 1] = *lineindex;

        const int nupperdowntrans = get_ndowntrans(element, ion, level) + 1;
        set_ndowntrans(element, ion, level, nupperdowntrans);

        const int nloweruptrans = get_nuptrans(element, ion, targetlevel) + 1;
        set_nuptrans(element, ion, targetlevel, nloweruptrans);

        totupdowntrans += 2;

        if (pass == 1 && globals::rank_in_node == 0) {
          const float A_ul = transitiontable[ii].A;
          const float coll_str = transitiontable[ii].coll_str;

          const auto g_ratio = stat_weight(element, ion, level) / stat_weight(element, ion, targetlevel);
          const float f_ul = g_ratio * ME * pow(CLIGHT, 3) / (8 * pow(QE * nu_trans * PI, 2)) * A_ul;
          assert_always(std::isfinite(f_ul));

          // printout("lineindex %d, element %d, ion %d, lower %d, upper %d, nu
          // %g\n",*lineindex,element,ion,level-i-1,level,nu_trans);

          temp_linelist.push_back({
              .nu = nu_trans,
              .einstein_A = A_ul,
              .elementindex = element,
              .ionindex = ion,
              .upperlevelindex = level,
              .lowerlevelindex = targetlevel,
          });

          // the line list has not been sorted yet, so the store the level index for now and
          // the index into the sorted line list will be set later

          globals::elements[element].ions[ion].levels[level].downtrans[nupperdowntrans - 1] = {
              .lineindex = -1,
              .targetlevelindex = targetlevel,
              .einstein_A = static_cast<float>(A_ul),
              .coll_str = coll_str,
              .osc_strength = f_ul,
              .forbidden = transitiontable[ii].forbidden};
          globals::elements[element].ions[ion].levels[targetlevel].uptrans[nloweruptrans - 1] = {
              .lineindex = -1,
              .targetlevelindex = level,
              .einstein_A = static_cast<float>(A_ul),
              .coll_str = coll_str,
              .osc_strength = f_ul,
              .forbidden = transitiontable[ii].forbidden};
        }

        /// This is not a metastable level.
        globals::elements[element].ions[ion].levels[level].metastable = false;

        (*lineindex)++;
      } else if (pass == 1 && globals::rank_in_node == 0) {
        // This is a new branch to deal with lines that have different types of transition. It should trip after a
        // transition is already known.
        const int linelistindex = transitions[level].to[level - targetlevel - 1];
        const float A_ul = transitiontable[ii].A;
        const float coll_str = transitiontable[ii].coll_str;

        const auto g_ratio = stat_weight(element, ion, level) / stat_weight(element, ion, targetlevel);
        const float f_ul = g_ratio * ME * pow(CLIGHT, 3) / (8 * pow(QE * nu_trans * PI, 2)) * A_ul;

        if ((temp_linelist[linelistindex].elementindex != element) || (temp_linelist[linelistindex].ionindex != ion) ||
            (temp_linelist[linelistindex].upperlevelindex != level) ||
            (temp_linelist[linelistindex].lowerlevelindex != targetlevel)) {
          printout("[input] Failure to identify level pair for duplicate bb-transition ... going to abort now\n");
          printout("[input]   element %d ion %d targetlevel %d level %d\n", element, ion, targetlevel, level);
          printout("[input]   transitions[level].to[level-targetlevel-1]=linelistindex %d\n",
                   transitions[level].to[level - targetlevel - 1]);
          printout("[input]   A_ul %g, coll_str %g\n", A_ul, coll_str);
          printout(
              "[input]   globals::linelist[linelistindex].elementindex %d, "
              "globals::linelist[linelistindex].ionindex %d, globals::linelist[linelistindex].upperlevelindex "
              "%d, globals::linelist[linelistindex].lowerlevelindex %d\n",
              temp_linelist[linelistindex].elementindex, temp_linelist[linelistindex].ionindex,
              temp_linelist[linelistindex].upperlevelindex, temp_linelist[linelistindex].lowerlevelindex);
          abort();
        }
        const int nupperdowntrans = get_ndowntrans(element, ion, level) + 1;

        const int nloweruptrans = get_nuptrans(element, ion, targetlevel) + 1;

        auto &downtransition = globals::elements[element].ions[ion].levels[level].downtrans[nupperdowntrans - 1];
        downtransition.einstein_A += A_ul;
        downtransition.osc_strength += f_ul;
        downtransition.coll_str = std::max(downtransition.coll_str, coll_str);

        auto &uptransition = globals::elements[element].ions[ion].levels[targetlevel].uptrans[nloweruptrans - 1];
        uptransition.einstein_A += A_ul;
        uptransition.osc_strength += f_ul;
        uptransition.coll_str = std::max(uptransition.coll_str, coll_str);
      }
    }
  }
#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

static auto calculate_nlevels_groundterm(int element, int ion) -> int {
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

static void read_atomicdata_files() {
  int totaluptrans = 0;
  int totaldowntrans = 0;

  auto compositiondata = fstream_required("compositiondata.txt", std::ios::in);

  FILE *adata = fopen_required("adata.txt", "r");

  printout("single_level_top_ion: %s\n", single_level_top_ion ? "true" : "false");
  printout("single_ground_level: %s\n", single_ground_level ? "true" : "false");
  /// initialize atomic data structure to number of elements
  int nelements_in = 0;
  assert_always(compositiondata >> nelements_in);
  set_nelements(nelements_in);
  globals::elements = static_cast<elementlist_entry *>(calloc(get_nelements(), sizeof(elementlist_entry)));
  assert_always(globals::elements != nullptr);

  /// Initialize the linelist
  std::vector<struct linelist_entry> temp_linelist;

  /// temperature to determine relevant ionstages
  int T_preset = 0;
  assert_always(compositiondata >> T_preset);
  assert_always(T_preset == 0);  // no longer in use
  int homogeneous_abundances = 0;
  assert_always(compositiondata >> homogeneous_abundances);
  assert_always(homogeneous_abundances == 0);  // no longer in use

  /// open transition data file
  auto ftransitiondata = fstream_required("transitiondata.txt", std::ios::in);

  int lineindex = 0;         /// counter to determine the total number of lines
  int uniqueionindex = 0;    // index into list of all ions of all elements
  int uniquelevelindex = 0;  // index into list of all levels of all ions of all elements
  int nbfcheck = 0;
  for (int element = 0; element < get_nelements(); element++) {
    /// read information about the next element which should be stored to memory
    int Z = 0;
    int nions = 0;
    int lowermost_ionstage = 0;
    int uppermost_ionstage = 0;
    int nlevelsmax_readin = 0;
    double abundance = NAN;
    double mass_amu = NAN;
    assert_always(compositiondata >> Z >> nions >> lowermost_ionstage >> uppermost_ionstage >> nlevelsmax_readin >>
                  abundance >> mass_amu);
    printout("readin compositiondata: next element Z %d, nions %d, lowermost %d, uppermost %d, nlevelsmax %d\n", Z,
             nions, lowermost_ionstage, uppermost_ionstage, nlevelsmax_readin);
    assert_always(Z > 0);
    assert_always(nions >= 0);
    assert_always(nions == 0 || (nions == uppermost_ionstage - lowermost_ionstage + 1));
    assert_always(abundance >= 0);
    assert_always(mass_amu >= 0);

    update_max_nions(nions);
    assert_always(nions <= get_max_nions());

    /// write this element's data to memory
    globals::elements[element].anumber = Z;
    globals::elements[element].nions = nions;
    globals::elements[element].abundance = abundance;  /// abundances are expected to be given by mass
    globals::elements[element].initstablemeannucmass = mass_amu * MH;
    increase_includedions(nions);

    /// Initialize the elements ionlist
    globals::elements[element].ions = static_cast<ionlist_entry *>(calloc(nions, sizeof(ionlist_entry)));
    assert_always(globals::elements[element].ions != nullptr);

    /// now read in data for all ions of the current element. before doing so initialize
    /// energy scale for the current element (all level energies are stored relative to
    /// the ground level of the neutral ion)
    double energyoffset = 0.;
    double ionpot = 0.;
    for (int ion = 0; ion < nions; ion++) {
      int nlevelsmax = nlevelsmax_readin;
      // printout("element %d ion %d\n", element, ion);
      /// calculate the current levels ground level energy
      assert_always(ionpot >= 0);
      energyoffset += ionpot;

      /// read information for the elements next ionstage
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
          double levelenergy = NAN;
          double statweight = NAN;
          int levelindex = 0;
          int ntransitions = 0;
          assert_always(
              fscanf(adata, "%d %lg %lg %d%*[^\n]\n", &levelindex, &levelenergy, &statweight, &ntransitions) == 4);
        }

        assert_always(fscanf(adata, "%d %d %d %lg\n", &adata_Z_in, &ionstage, &nlevels, &ionpot) == 4);
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

      /// and proceed through the transitionlist till we match this ionstage (if it was not the neutral one)
      int transdata_Z_in = -1;
      int transdata_ionstage_in = -1;
      int tottransitions_in_file = 0;
      std::string line;
      while (transdata_Z_in != Z || transdata_ionstage_in != ionstage) {
        // skip over table
        for (int i = 0; i < tottransitions_in_file; i++) {
          assert_always(getline(ftransitiondata, line));
        }
        assert_always(get_noncommentline(ftransitiondata, line));
        assert_always(
            sscanf(line.c_str(), "%d %d %d", &transdata_Z_in, &transdata_ionstage_in, &tottransitions_in_file) == 3);
      }

      printout("transdata header matched: transdata_Z_in %d, transdata_ionstage_in %d, tottransitions %d\n",
               transdata_Z_in, transdata_ionstage_in, tottransitions_in_file);
      assert_always(tottransitions_in_file >= 0);

      // read the data for the levels and set up the list of possible transitions for each level
      /// store the ions data to memory and set up the ions zeta and levellist
      globals::elements[element].ions[ion].ionstage = ionstage;
      globals::elements[element].ions[ion].nlevels = nlevelsmax;
      globals::elements[element].ions[ion].ionisinglevels = 0;
      globals::elements[element].ions[ion].maxrecombininglevel = 0;
      globals::elements[element].ions[ion].ionpot = ionpot * EV;
      globals::elements[element].ions[ion].nlevels_groundterm = -1;
      globals::elements[element].ions[ion].uniqueionindex = uniqueionindex;

      globals::elements[element].ions[ion].Alpha_sp = static_cast<float *>(calloc(TABLESIZE, sizeof(float)));
      assert_always(globals::elements[element].ions[ion].Alpha_sp != nullptr);

      globals::elements[element].ions[ion].levels =
          static_cast<struct levellist_entry *>(calloc(nlevelsmax, sizeof(struct levellist_entry)));
      assert_always(globals::elements[element].ions[ion].levels != nullptr);

      auto *transitions = static_cast<struct transitions *>(calloc(nlevelsmax, sizeof(struct transitions)));
      assert_always(transitions != nullptr);

      read_ion_levels(adata, element, ion, nions, nlevels, nlevelsmax, energyoffset, ionpot, transitions);

      int tottransitions = tottransitions_in_file;

      if (single_level_top_ion && ion == nions - 1)  // limit the top ion to one level and no transitions
      {
        tottransitions = 0;
      }

      assert_always(transdata_Z_in == Z);
      assert_always(transdata_ionstage_in == ionstage);

      /// read in the level and transition data for this ion
      std::vector<struct transitiontable_entry> transitiontable;
      transitiontable.reserve(tottransitions);

      /// load transition table for the CURRENT ion to temporary memory

      // first <nlevels_requiretransitions> levels will be collisionally
      // coupled to the first <nlevels_requiretransitions_upperlevels> levels (assumed forbidden)
      // use 0 to disable adding extra transitions

      int nlevels_requiretransitions = NLEVELS_REQUIRETRANSITIONS(Z, ionstage);
      int nlevels_requiretransitions_upperlevels = nlevelsmax;  // no effect if previous line is zero

      nlevels_requiretransitions = std::min(nlevelsmax, nlevels_requiretransitions);
      nlevels_requiretransitions_upperlevels = std::min(nlevelsmax, nlevels_requiretransitions_upperlevels);

      read_ion_transitions(ftransitiondata, tottransitions_in_file, &tottransitions, transitiontable,
                           nlevels_requiretransitions, nlevels_requiretransitions_upperlevels);

      add_transitions_to_unsorted_linelist(element, ion, nlevelsmax, transitiontable, transitions, &lineindex,
                                           temp_linelist);

      free(transitions[0].to);
      free(transitions);
      transitions = nullptr;
      transitiontable.clear();

      for (int level = 0; level < nlevelsmax; level++) {
        globals::elements[element].ions[ion].levels[level].uniquelevelindex = uniquelevelindex;
        globals::elements[element].ions[ion].levels[level].nphixstargets = 0;
        globals::elements[element].ions[ion].levels[level].phixstargets = nullptr;
        globals::elements[element].ions[ion].levels[level].photoion_xs = nullptr;
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
  fclose(adata);
  printout("nbfcheck %d\n", nbfcheck);

  /// Save the linecounters value to the global variable containing the number of lines
  globals::nlines = lineindex;
  printout("nlines %d\n", globals::nlines);
  if (globals::rank_in_node == 0) {
    assert_always(globals::nlines == static_cast<int>(temp_linelist.size()));
  }

  if (T_preset > 0) {
    abort();
  }

  /// Set up the list of allowed upward transitions for each level
  printout("total uptrans %d\n", totaluptrans);
  printout("total downtrans %d\n", totaldowntrans);

  printout("[info] mem_usage: transition lists occupy %.3f MB (this rank) and %.3f MB (shared on node)\n",
           2 * uniquelevelindex * sizeof(struct level_transition *) / 1024. / 1024.,
           (totaluptrans + totaldowntrans) * sizeof(struct level_transition) / 1024. / 1024.);

  if (globals::rank_in_node == 0) {
    // sort the lineline in descending frequency
    std::sort(temp_linelist.begin(), temp_linelist.end(),
              [](const auto &a, const auto &b) { return static_cast<bool>(a.nu > b.nu); });

    for (int i = 0; i < globals::nlines - 1; i++) {
      const double nu = temp_linelist[i].nu;
      const double nu_next = temp_linelist[i + 1].nu;
      if (fabs(nu_next - nu) < (1.e-10 * nu)) {
        auto *a1 = &temp_linelist[i];
        auto *a2 = &temp_linelist[i + 1];

        if ((a1->elementindex == a2->elementindex) && (a1->ionindex == a2->ionindex) &&
            (a1->lowerlevelindex == a2->lowerlevelindex) && (a1->upperlevelindex == a2->upperlevelindex)) {
          printout("Duplicate transition line? %s\n", a1->nu == a2->nu ? "nu match exact" : "close to nu match");
          printout("a: Z=%d ionstage %d lower %d upper %d nu %g lambda %g\n", get_atomicnumber(a1->elementindex),
                   get_ionstage(a1->elementindex, a1->ionindex), a1->lowerlevelindex, a1->upperlevelindex, a1->nu,
                   1e8 * CLIGHT / a1->nu);
          printout("b: Z=%d ionstage %d lower %d upper %d nu %g lambda %g\n", get_atomicnumber(a2->elementindex),
                   get_ionstage(a2->elementindex, a2->ionindex), a2->lowerlevelindex, a2->upperlevelindex, a2->nu,
                   1e8 * CLIGHT / a2->nu);
        }
      }
    }
  }

  // create a linelist shared on node and then copy data across, freeing the local copy
  struct linelist_entry *nonconstlinelist = nullptr;
#ifdef MPI_ON
  MPI_Win win = MPI_WIN_NULL;

  int my_rank_lines = globals::nlines / globals::node_nprocs;
  // rank_in_node 0 gets any remainder
  if (globals::rank_in_node == 0) {
    my_rank_lines += globals::nlines - (my_rank_lines * globals::node_nprocs);
  }

  MPI_Aint size = my_rank_lines * sizeof(linelist_entry);
  int disp_unit = sizeof(linelist_entry);
  MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node, &nonconstlinelist, &win);

  MPI_Win_shared_query(win, 0, &size, &disp_unit, &nonconstlinelist);
#else
  nonconstlinelist = static_cast<struct linelist_entry *>(malloc(globals::nlines * sizeof(linelist_entry)));
#endif

  if (globals::rank_in_node == 0) {
    memcpy(static_cast<void *>(nonconstlinelist), temp_linelist.data(), globals::nlines * sizeof(linelist_entry));
    temp_linelist.clear();
  }

#ifdef MPI_ON
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  globals::linelist = nonconstlinelist;
  nonconstlinelist = nullptr;
  printout("[info] mem_usage: linelist occupies %.3f MB (node shared memory)\n",
           globals::nlines * sizeof(struct linelist_entry) / 1024. / 1024);

  /// Save sorted linelist into a file
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

  time_t const time_start_establish_linelist_connections = time(nullptr);
  for (lineindex = 0; lineindex < globals::nlines; lineindex++) {
    const auto &line = globals::linelist[lineindex];
    const int element = line.elementindex;
    const int ion = line.ionindex;
    const int lowerlevel = line.lowerlevelindex;
    const int upperlevel = line.upperlevelindex;

    // there is never more than one transition per pair of levels,
    // so find the first matching the upper and lower transition

    const int nupperdowntrans = get_ndowntrans(element, ion, upperlevel);
    auto *downtranslist = globals::elements[element].ions[ion].levels[upperlevel].downtrans;
    auto *downtrans = std::find_if(downtranslist, downtranslist + nupperdowntrans,
                                   [=](const auto &downtrans) { return downtrans.targetlevelindex == lowerlevel; });
    assert_always(downtrans != (downtranslist + nupperdowntrans));
    // assert_always(downtrans->targetlevelindex == lowerlevel);
    downtrans->lineindex = lineindex;

    const int nloweruptrans = get_nuptrans(element, ion, lowerlevel);
    auto *uptranslist = globals::elements[element].ions[ion].levels[lowerlevel].uptrans;
    auto *uptrans = std::find_if(uptranslist, uptranslist + nloweruptrans,
                                 [=](const auto &uptrans) { return uptrans.targetlevelindex == upperlevel; });
    assert_always(uptrans != (uptranslist + nloweruptrans));
    // assert_always(uptrans->targetlevelindex == upperlevel);
    uptrans->lineindex = lineindex;
  }

  printout("  took %ds\n", time(nullptr) - time_start_establish_linelist_connections);

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

  globals::nbfcontinua_ground = 0;
  globals::nbfcontinua = 0;

  // read in photoionisation cross sections
  phixs_file_version_exists[1] = std::filesystem::exists(phixsdata_filenames[1]);
  phixs_file_version_exists[2] = std::filesystem::exists(phixsdata_filenames[2]);
  assert_always(phixs_file_version_exists[1] || phixs_file_version_exists[2]);  // at least one must exist
  if (phixs_file_version_exists[1] && phixs_file_version_exists[2]) {
    printout(
        "Reading two phixs files: Reading phixsdata_v2.txt first so we use NPHIXSPOINTS and NPHIXSNUINCREMENT "
        "from phixsdata_v2.txt to interpolate the phixsdata.txt data\n");
  }
  if (phixs_file_version_exists[2]) {
    read_phixs_data(2);
  }
  if (phixs_file_version_exists[1]) {
    read_phixs_data(1);
  }

  int cont_index = -1;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_nlevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        const int nphixstargets = get_nphixstargets(element, ion, level);
        globals::elements[element].ions[ion].levels[level].cont_index =
            (nphixstargets > 0) ? cont_index : std::numeric_limits<int>::max();
        cont_index -= nphixstargets;
      }

      // below is just an extra warning consistency check
      const int nlevels_groundterm = globals::elements[element].ions[ion].nlevels_groundterm;

      // all levels in the ground term should be photoionisation targets from the lower ground state
      if (ion > 0 && ion < get_nions(element) - 1) {
        const int nphixstargets = get_nphixstargets(element, ion - 1, 0);
        if (nphixstargets > 0 && get_phixsupperlevel(element, ion - 1, 0, 0) == 0) {
          const int phixstargetlevels = get_phixsupperlevel(element, ion - 1, 0, nphixstargets - 1) + 1;

          if (nlevels_groundterm != phixstargetlevels) {
            printout("WARNING: Z=%d ion_stage %d nlevels_groundterm %d phixstargetlevels(ion-1) %d.\n",
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
}

static auto search_groundphixslist(double nu_edge, int *index_in_groundlevelcontestimator, int el, int in, int ll)
    -> int
/// Return the closest ground level continuum index to the given edge
/// frequency. If the given edge frequency is redder than the reddest
/// continuum return -1.
/// NB: groundphixslist must be in ascending order.
{
  assert_always((USE_LUT_PHOTOION || USE_LUT_BFHEATING));
  int index = 0;

  if (nu_edge < globals::groundcont[0].nu_edge) {
    index = -1;
    *index_in_groundlevelcontestimator = -1;
  } else {
    int i;
    int element = -1;
    int ion = -1;
    for (i = 1; i < globals::nbfcontinua_ground; i++) {
      if (nu_edge < globals::groundcont[i].nu_edge) {
        break;
      }
    }
    if (i == globals::nbfcontinua_ground) {
      element = globals::groundcont[i - 1].element;
      ion = globals::groundcont[i - 1].ion;
      int const level = globals::groundcont[i - 1].level;
      if (element == el && ion == in && level == ll) {
        index = i - 1;
      } else {
        printout(
            "[fatal] search_groundphixslist: element %d, ion %d, level %d has edge_frequency %g equal to the "
            "bluest ground-level continuum\n",
            el, in, ll, nu_edge);
        printout(
            "[fatal] search_groundphixslist: bluest ground level continuum is element %d, ion %d, level %d at "
            "nu_edge %g\n",
            element, ion, level, globals::groundcont[i - 1].nu_edge);
        printout("[fatal] search_groundphixslist: i %d, nbfcontinua_ground %d\n", i, globals::nbfcontinua_ground);
        printout(
            "[fatal] This shouldn't happen, is hoewever possible if there are multiple levels in the adata file at "
            "energy=0\n");
        for (int looplevels = 0; looplevels < get_nlevels(el, in); looplevels++) {
          printout("[fatal]   element %d, ion %d, level %d, energy %g\n", el, in, looplevels,
                   epsilon(el, in, looplevels));
        }
        printout("[fatal] Abort omitted ... MAKE SURE ATOMIC DATA ARE CONSISTENT\n");
        index = i - 1;
        // abort();
      }
    } else {
      const double left_diff = nu_edge - globals::groundcont[i - 1].nu_edge;
      const double right_diff = globals::groundcont[i].nu_edge - nu_edge;
      index = (left_diff <= right_diff) ? i - 1 : i;
      element = globals::groundcont[index].element;
      ion = globals::groundcont[index].ion;
    }
    *index_in_groundlevelcontestimator = element * get_max_nions() + ion;
  }

  return index;
}

static void setup_cellhistory() {
  /// SET UP THE CELL HISTORY
  ///======================================================
  /// Stack which holds information about population and other cell specific data
  /// ===> move to update_packets
  globals::cellhistory = static_cast<struct cellhistory *>(malloc(get_max_threads() * sizeof(struct cellhistory)));
  assert_always(globals::cellhistory != nullptr);

#ifdef _OPENMP
#pragma omp parallel
  {
#endif
    size_t mem_usage_cellhistory = 0;
    mem_usage_cellhistory += sizeof(struct cellhistory);

    printout("[info] input: initializing cellhistory for thread %d ...\n", tid);

    globals::cellhistory[tid].cellnumber = -99;

    mem_usage_cellhistory += globals::ncoolingterms * sizeof(double);
    globals::cellhistory[tid].cooling_contrib = static_cast<double *>(calloc(globals::ncoolingterms, sizeof(double)));

    for (int element = 0; element < get_nelements(); element++) {
      for (int ion = 0; ion < get_nions(element); ion++) {
        globals::cellhistory[tid].cooling_contrib[kpkt::get_coolinglistoffset(element, ion)] = COOLING_UNDEFINED;
      }
    }

    printout("[info] mem_usage: cellhistory coolinglist contribs for thread %d occupies %.3f MB\n", tid,
             globals::ncoolingterms * sizeof(double) / 1024. / 1024.);

    mem_usage_cellhistory += get_nelements() * sizeof(struct chelements);
    globals::cellhistory[tid].chelements =
        static_cast<struct chelements *>(malloc(get_nelements() * sizeof(struct chelements)));

    assert_always(globals::cellhistory[tid].chelements != nullptr);

    size_t chlevelblocksize = 0;
    size_t chphixsblocksize = 0;
    int chtransblocksize = 0;
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        const int nlevels = get_nlevels(element, ion);
        chlevelblocksize += nlevels * sizeof(struct chlevels);

        for (int level = 0; level < nlevels; level++) {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          chphixsblocksize += nphixstargets * sizeof(chphixstargets_t);

          const int ndowntrans = get_ndowntrans(element, ion, level);
          const int nuptrans = get_nuptrans(element, ion, level);
          chtransblocksize += (2 * ndowntrans + nuptrans);
        }
      }
    }
    assert_always(chlevelblocksize > 0);
    globals::cellhistory[tid].ch_all_levels = static_cast<struct chlevels *>(malloc(chlevelblocksize));
    chphixstargets_t *chphixstargetsblock =
        chphixsblocksize > 0 ? static_cast<chphixstargets_t *>(malloc(chphixsblocksize)) : nullptr;
    mem_usage_cellhistory += chlevelblocksize + chphixsblocksize;

    mem_usage_cellhistory += chtransblocksize * sizeof(double);
    double *const chtransblock =
        chtransblocksize > 0 ? static_cast<double *>(malloc(chtransblocksize * sizeof(double))) : nullptr;

    int alllevelindex = 0;
    int allphixstargetindex = 0;
    int chtransindex = 0;
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      mem_usage_cellhistory += nions * sizeof(struct chions);
      globals::cellhistory[tid].chelements[element].chions =
          static_cast<struct chions *>(malloc(nions * sizeof(struct chions)));
      assert_always(globals::cellhistory[tid].chelements[element].chions != nullptr);

      for (int ion = 0; ion < nions; ion++) {
        const int nlevels = get_nlevels(element, ion);
        globals::cellhistory[tid].chelements[element].chions[ion].chlevels =
            &globals::cellhistory[tid].ch_all_levels[alllevelindex];

        assert_always(alllevelindex == get_uniquelevelindex(element, ion, 0));
        alllevelindex += nlevels;

        for (int level = 0; level < nlevels; level++) {
          struct chlevels *chlevel = &globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level];
          const int nphixstargets = get_nphixstargets(element, ion, level);
          chlevel->chphixstargets = chphixsblocksize > 0 ? &chphixstargetsblock[allphixstargetindex] : nullptr;
          allphixstargetindex += nphixstargets;
        }

        for (int level = 0; level < nlevels; level++) {
          struct chlevels *chlevel = &globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level];
          const int ndowntrans = get_ndowntrans(element, ion, level);

          chlevel->sum_epstrans_rad_deexc = &chtransblock[chtransindex];
          chtransindex += ndowntrans;
        }

        for (int level = 0; level < nlevels; level++) {
          struct chlevels *chlevel = &globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level];
          const int ndowntrans = get_ndowntrans(element, ion, level);
          chlevel->sum_internal_down_same = &chtransblock[chtransindex];
          chtransindex += ndowntrans;
        }

        for (int level = 0; level < nlevels; level++) {
          struct chlevels *chlevel = &globals::cellhistory[tid].chelements[element].chions[ion].chlevels[level];
          const int nuptrans = get_nuptrans(element, ion, level);
          chlevel->sum_internal_up_same = &chtransblock[chtransindex];
          chtransindex += nuptrans;
        }
      }
    }
    assert_always(chtransindex == chtransblocksize);

    assert_always(globals::nbfcontinua >= 0);
    globals::cellhistory[tid].ch_allcont_departureratios =
        static_cast<double *>(malloc(globals::nbfcontinua * sizeof(double)));
    mem_usage_cellhistory += globals::nbfcontinua * sizeof(double);

    printout("[info] mem_usage: cellhistory for thread %d occupies %.3f MB\n", tid,
             mem_usage_cellhistory / 1024. / 1024.);
#ifdef _OPENMP
  }
#endif
}

static void write_bflist_file(int includedphotoiontransitions) {
  if (includedphotoiontransitions > 0) {
    globals::bflist = static_cast<struct bflist_t *>(malloc(includedphotoiontransitions * sizeof(struct bflist_t)));
  } else {
    globals::bflist = nullptr;
  }

  FILE *bflist_file = nullptr;
  if (globals::rank_global == 0) {
    bflist_file = fopen_required("bflist.out", "w");
    fprintf(bflist_file, "%d\n", includedphotoiontransitions);
  }
  int i = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      const int nlevels = get_ionisinglevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        for (int phixstargetindex = 0; phixstargetindex < get_nphixstargets(element, ion, level); phixstargetindex++) {
          const int upperionlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);
          globals::bflist[i].elementindex = element;
          globals::bflist[i].ionindex = ion;
          globals::bflist[i].levelindex = level;
          globals::bflist[i].phixstargetindex = phixstargetindex;

          if (globals::rank_global == 0) {
            fprintf(bflist_file, "%d %d %d %d %d\n", i, element, ion, level, upperionlevel);
          }

          const int et = -1 - i;

          assert_always(et == get_continuumindex(element, ion, level, upperionlevel));

          // check the we don't overload the same packet emission type numbers
          // as the special values for free-free scattering and not set
          assert_always(et != EMTYPE_NOTSET);
          assert_always(et != EMTYPE_FREEFREE);
          i++;
        }
      }
    }
  }
  assert_always(i == includedphotoiontransitions);
  if (globals::rank_global == 0) {
    fclose(bflist_file);
  }
}

static void setup_phixs_list() {
  // set up the photoionisation transition lists
  // and temporary gamma/kappa lists for each thread

  printout("[info] read_atomicdata: number of bfcontinua %d\n", globals::nbfcontinua);
  printout("[info] read_atomicdata: number of ground-level bfcontinua %d\n", globals::nbfcontinua_ground);

  globals::phixslist = static_cast<struct phixslist *>(malloc(get_max_threads() * sizeof(struct phixslist)));
  assert_always(globals::phixslist != nullptr);

  /// MK: 2012-01-19
  /// To fix the OpenMP problem on BlueGene machines this parallel section was removed and replaced by
  /// a serial loop which intializes the phixslist data structure for all threads in a loop. I'm still
  /// not sure why this causes a problem at all and on BlueGene architectures in particular. However,
  /// it seems to fix the problem.
  // #ifdef _OPENMP
  //   #pragma omp parallel private(i,element,ion,level,nions,nlevels,epsilon_upper,E_threshold,nu_edge)
  //   {
  // #endif
  for (int itid = 0; itid < get_max_threads(); itid++) {
    /// Number of ground level bf-continua equals the total number of included ions minus the number
    /// of included elements, because the uppermost ionisation stages can't ionise.
    if ((USE_LUT_PHOTOION || USE_LUT_BFHEATING) && globals::nbfcontinua_ground > 0) {
      globals::phixslist[itid].groundcont_gamma_contr =
          static_cast<double *>(malloc(globals::nbfcontinua_ground * sizeof(double)));
      assert_always(globals::phixslist[itid].groundcont_gamma_contr != nullptr);

      for (int groundcontindex = 0; groundcontindex < globals::nbfcontinua_ground; groundcontindex++) {
        globals::phixslist[itid].groundcont_gamma_contr[groundcontindex] = 0.;
      }
    } else {
      globals::phixslist[itid].groundcont_gamma_contr = nullptr;
    }

    if (globals::nbfcontinua > 0) {
      globals::phixslist[itid].kappa_bf_sum = static_cast<double *>(malloc(globals::nbfcontinua * sizeof(double)));
      assert_always(globals::phixslist[itid].kappa_bf_sum != nullptr);

      if constexpr (DETAILED_BF_ESTIMATORS_ON) {
        globals::phixslist[itid].gamma_contr = static_cast<double *>(malloc(globals::nbfcontinua * sizeof(double)));
        assert_always(globals::phixslist[itid].gamma_contr != nullptr);
      }

      for (int allcontindex = 0; allcontindex < globals::nbfcontinua; allcontindex++) {
        globals::phixslist[itid].kappa_bf_sum[allcontindex] = 0.;

        if constexpr (DETAILED_BF_ESTIMATORS_ON) {
          globals::phixslist[itid].gamma_contr[allcontindex] = 0.;
        }
      }
    } else {
      globals::phixslist[itid].kappa_bf_sum = nullptr;
      globals::phixslist[itid].gamma_contr = nullptr;
    }

    printout("[info] mem_usage: phixslist[tid].kappa_bf_contr for thread %d occupies %.3f MB\n", itid,
             globals::nbfcontinua * sizeof(double) / 1024. / 1024.);
  }

  if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
    globals::groundcont =
        static_cast<struct groundphixslist *>(malloc(globals::nbfcontinua_ground * sizeof(struct groundphixslist)));
    assert_always(globals::groundcont != nullptr);
  }

  if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
    int groundcontindex = 0;
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++) {
        const int nlevels_groundterm = get_nlevels_groundterm(element, ion);
        for (int level = 0; level < nlevels_groundterm; level++) {
          const int nphixstargets = get_nphixstargets(element, ion, level);
          for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
            // const int upperlevel = get_phixsupperlevel(element, ion,level, 0);
            // const double E_threshold = epsilon(element, ion + 1, upperlevel) - epsilon(element, ion, level);
            const double E_threshold = get_phixs_threshold(element, ion, level, phixstargetindex);
            const double nu_edge = E_threshold / H;
            assert_always(groundcontindex < globals::nbfcontinua_ground);
            globals::groundcont[groundcontindex].element = element;
            globals::groundcont[groundcontindex].ion = ion;
            globals::groundcont[groundcontindex].level = level;
            globals::groundcont[groundcontindex].nu_edge = nu_edge;
            globals::groundcont[groundcontindex].phixstargetindex = phixstargetindex;
            groundcontindex++;
          }
        }
      }
    }
    assert_always(groundcontindex == globals::nbfcontinua_ground);
    std::sort(globals::groundcont, globals::groundcont + globals::nbfcontinua_ground,
              [](const auto &a, const auto &b) { return static_cast<bool>(a.nu_edge < b.nu_edge); });
  }

  auto *nonconstallcont =
      static_cast<struct fullphixslist *>(malloc(globals::nbfcontinua * sizeof(struct fullphixslist)));
  printout("[info] mem_usage: photoionisation list occupies %.3f MB\n",
           globals::nbfcontinua * (sizeof(fullphixslist)) / 1024. / 1024.);
  int nbftables = 0;
  int allcontindex = 0;
  for (int element = 0; element < get_nelements(); element++) {
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions - 1; ion++) {
      const int nlevels = get_ionisinglevels(element, ion);
      for (int level = 0; level < nlevels; level++) {
        const int nphixstargets = get_nphixstargets(element, ion, level);

        if (nphixstargets > 0) {
          nbftables++;
        }

        for (int phixstargetindex = 0; phixstargetindex < nphixstargets; phixstargetindex++) {
          // const int upperlevel = get_phixsupperlevel(element, ion,level, 0);
          // const double E_threshold = epsilon(element, ion + 1, upperlevel) - epsilon(element, ion, level);
          const double E_threshold = get_phixs_threshold(element, ion, level, phixstargetindex);
          const double nu_edge = E_threshold / H;

          assert_always(allcontindex < globals::nbfcontinua);
          nonconstallcont[allcontindex].nu_edge = nu_edge;
          nonconstallcont[allcontindex].element = element;
          nonconstallcont[allcontindex].ion = ion;
          nonconstallcont[allcontindex].level = level;
          nonconstallcont[allcontindex].phixstargetindex = phixstargetindex;
          nonconstallcont[allcontindex].probability = get_phixsprobability(element, ion, level, phixstargetindex);
          nonconstallcont[allcontindex].upperlevel = get_phixsupperlevel(element, ion, level, phixstargetindex);

          if constexpr (USE_LUT_PHOTOION || USE_LUT_BFHEATING) {
            int index_in_groundlevelcontestimator = 0;
            nonconstallcont[allcontindex].index_in_groundphixslist =
                search_groundphixslist(nu_edge, &index_in_groundlevelcontestimator, element, ion, level);

            globals::elements[element].ions[ion].levels[level].closestgroundlevelcont =
                index_in_groundlevelcontestimator;
          }
          allcontindex++;
        }
      }
    }
  }

  assert_always(allcontindex == globals::nbfcontinua);
  assert_always(globals::nbfcontinua >= 0);  // was initialised as -1 before startup

  if (globals::nbfcontinua > 0) {
    // indicies above were temporary only. continum index should be to the sorted list
    std::sort(nonconstallcont, nonconstallcont + globals::nbfcontinua,
              [](const auto &a, const auto &b) { return static_cast<bool>(a.nu_edge < b.nu_edge); });

    globals::allcont_nu_edge = static_cast<double *>(malloc(globals::nbfcontinua * sizeof(double)));

// copy the photoionisation tables into one contiguous block of memory
#ifdef MPI_ON
    float *allphixsblock = nullptr;
    MPI_Win win_allphixsblock = MPI_WIN_NULL;
    MPI_Aint size = (globals::rank_in_node == 0) ? nbftables * globals::NPHIXSPOINTS * sizeof(float) : 0;
    int disp_unit = sizeof(linelist_entry);

    MPI_Win_allocate_shared(size, disp_unit, MPI_INFO_NULL, globals::mpi_comm_node, &allphixsblock, &win_allphixsblock);
    MPI_Win_shared_query(win_allphixsblock, MPI_PROC_NULL, &size, &disp_unit, &allphixsblock);

    MPI_Barrier(MPI_COMM_WORLD);
#else
    auto *allphixsblock = static_cast<float *>(malloc(nbftables * globals::NPHIXSPOINTS * sizeof(float)));
#endif

    assert_always(allphixsblock != nullptr);
    int nbftableschanged = 0;
    for (int i = 0; i < globals::nbfcontinua; i++) {
      globals::allcont_nu_edge[i] = nonconstallcont[i].nu_edge;

      const int element = nonconstallcont[i].element;
      const int ion = nonconstallcont[i].ion;
      const int level = nonconstallcont[i].level;
      const int phixstargetindex = nonconstallcont[i].phixstargetindex;

      // different targets share the same cross section table, so don't repeat this process
      if (phixstargetindex == 0) {
        if (globals::rank_in_node == 0) {
          memcpy(allphixsblock, globals::elements[element].ions[ion].levels[level].photoion_xs,
                 globals::NPHIXSPOINTS * sizeof(float));
        }

        free(globals::elements[element].ions[ion].levels[level].photoion_xs);
        globals::elements[element].ions[ion].levels[level].photoion_xs = allphixsblock;

        allphixsblock += globals::NPHIXSPOINTS;
        nbftableschanged++;
      }
    }
    assert_always(nbftableschanged == nbftables);
#ifdef MPI_ON
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    for (int i = 0; i < globals::nbfcontinua; i++) {
      const int element = nonconstallcont[i].element;
      const int ion = nonconstallcont[i].ion;
      const int level = nonconstallcont[i].level;
      nonconstallcont[i].photoion_xs = globals::elements[element].ions[ion].levels[level].photoion_xs;
    }
  }
  globals::allcont = nonconstallcont;
  nonconstallcont = nullptr;

  setup_photoion_luts();
}

static void read_atomicdata() {
  read_atomicdata_files();

  printout("included ions %d\n", get_includedions());

  /// INITIALISE THE ABSORPTION/EMISSION COUNTERS ARRAYS
  if constexpr (RECORD_LINESTAT) {
    globals::ecounter = static_cast<int *>(malloc(globals::nlines * sizeof(int)));
    assert_always(globals::ecounter != nullptr);

    globals::acounter = static_cast<int *>(malloc(globals::nlines * sizeof(int)));
    assert_always(globals::acounter != nullptr);
  }

  kpkt::setup_coolinglist();

  setup_cellhistory();

  /// Printout some information about the read-in model atom

  int includedlevels = 0;
  int includedionisinglevels = 0;
  int includedboundboundtransitions = 0;
  int includedphotoiontransitions = 0;
  printout("[input] this simulation contains\n");
  printout("----------------------------------\n");
  for (int element = 0; element < get_nelements(); element++) {
    printout("[input]  element %d (Z=%2d %s)\n", element, get_atomicnumber(element),
             decay::get_elname(get_atomicnumber(element)));
    const int nions = get_nions(element);
    for (int ion = 0; ion < nions; ion++) {
      int ion_photoiontransitions = 0;
      int ion_bbtransitions = 0;
      for (int level = 0; level < get_nlevels(element, ion); level++) {
        ion_photoiontransitions += get_nphixstargets(element, ion, level);
        ion_bbtransitions += get_nuptrans(element, ion, level);
      }

      printout(
          "[input]    ion_stage %d: %4d levels (%d in groundterm, %4d ionising) %7d lines %6d bf transitions "
          "(epsilon_ground: %7.2f eV)\n",
          get_ionstage(element, ion), get_nlevels(element, ion), get_nlevels_groundterm(element, ion),
          get_ionisinglevels(element, ion), ion_bbtransitions, ion_photoiontransitions, epsilon(element, ion, 0) / EV);

      includedlevels += get_nlevels(element, ion);
      includedionisinglevels += get_ionisinglevels(element, ion);
      includedphotoiontransitions += ion_photoiontransitions;
      includedboundboundtransitions += ion_bbtransitions;
    }
  }
  assert_always(includedphotoiontransitions == globals::nbfcontinua);
  assert_always(globals::nlines == includedboundboundtransitions);

  printout("[input]  in total %d ions, %d levels (%d ionising), %d lines, %d photoionisation transitions\n",
           get_includedions(), includedlevels, includedionisinglevels, globals::nlines, globals::nbfcontinua);

  write_bflist_file(globals::nbfcontinua);

  setup_phixs_list();

  /// set-up/gather information for nlte stuff

  globals::total_nlte_levels = 0;
  int n_super_levels = 0;

  if (NLTE_POPS_ON) {
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions; ion++) {
        globals::elements[element].ions[ion].first_nlte = globals::total_nlte_levels;
        const int nlevels = get_nlevels(element, ion);
        int fullnlteexcitedlevelcount = 0;
        for (int level = 1; level < nlevels; level++) {
          if (is_nlte(element, ion, level)) {
            fullnlteexcitedlevelcount++;
            globals::total_nlte_levels++;
          }
        }

        const bool has_superlevel = (nlevels > (fullnlteexcitedlevelcount + 1));
        if (has_superlevel) {
          // If there are more levels that the ground state + the number of NLTE levels then we need an extra
          // slot to store data for the "superlevel", which is a representation of all the other levels that
          // are not treated in detail.
          globals::total_nlte_levels++;
          n_super_levels++;
        }

        globals::elements[element].ions[ion].nlevels_nlte = fullnlteexcitedlevelcount;

        assert_always(has_superlevel == ion_has_superlevel(element, ion));

        printout("[input]  element %2d Z=%2d ion_stage %2d has %5d NLTE excited levels%s. Starting at %d\n", element,
                 get_atomicnumber(element), get_ionstage(element, ion), fullnlteexcitedlevelcount,
                 has_superlevel ? " plus a superlevel" : "", globals::elements[element].ions[ion].first_nlte);
      }
    }
  }

  printout("[input] Total NLTE levels: %d, of which %d are superlevels\n", globals::total_nlte_levels, n_super_levels);
}

void input(int rank)
/// To govern the input. For now hardwire everything.
{
  globals::n_titer = 1;
  globals::initial_iteration = false;

  printout("[info] input: do n_titer %d iterations per timestep\n", globals::n_titer);
  if (globals::n_titer > 1) {
#ifndef DO_TITER
    printout("[fatal] input: n_titer > 1, but DO_TITER not defined ... abort\n");
    abort();
#endif
  } else if (globals::n_titer == 1) {
#ifdef DO_TITER
    printout("[warning] input: n_titer = 1 but DO_TITER defined, remove DO_TITER to save memory\n");
#endif
  } else {
    printout("[fatal] input: no valid value for n_titer selected\n");
    abort();
  }

  /// Read in parameters from input.txt
  read_parameterfile(rank);

  /// Read in parameters from vpkt.txt
  if (VPKT_ON) {
    read_parameterfile_vpkt();
  }

  read_atomicdata();

#ifdef MPI_ON
  const time_t time_before_barrier = time(nullptr);
  printout("barrier after read_atomicdata(): time before barrier %d, ", static_cast<int>(time_before_barrier));
  MPI_Barrier(MPI_COMM_WORLD);
  printout("time after barrier %d (waited %d seconds)\n", static_cast<int>(time(nullptr)),
           static_cast<int>(time(nullptr) - time_before_barrier));
#endif

  grid::read_ejecta_model();
}

static auto getline(std::istream &input, std::string &line) -> bool
// return true if line read, false if not (EOF)
{
  return !(!std::getline(input, line));
}

auto get_noncommentline(std::fstream &input, std::string &line) -> bool
// read the next line, skipping any comment lines beginning with '#'
{
  while (true) {
    bool const linefound = getline(input, line);
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

void read_parameterfile(int rank)
/// Subroutine to read in input parameters from input.txt.
{
  unsigned long pre_zseed = 0;

  auto file = fstream_required("input.txt", std::ios::in);

  std::string line;

  assert_always(get_noncommentline(file, line));

  uint_fast64_t zseed_input = 0;
  std::stringstream(line) >> zseed_input;

  if (zseed_input > 0) {
    pre_zseed = zseed_input;  // random number seed
    printout("using input.txt specified random number seed of %lu\n", pre_zseed);
  } else {
    // pre_zseed = time(nullptr);
    pre_zseed = std::random_device{}();
    printout("randomly-generated random number seed is %lu\n", pre_zseed);
  }

#ifdef _OPENMP
#pragma omp parallel
  {
#endif
    /// For MPI parallelisation, the random seed is changed based on the rank of the process
    /// For OpenMP parallelisation rng is a threadprivate variable and the seed changed according
    /// to the thread-ID tid.
    const uint_fast64_t zseed = pre_zseed + 13 * (rank * get_num_threads() + tid);
    printout("rank %d: thread %d has zseed %lu\n", rank, tid, zseed);
    rng_init(zseed);
    /// call it a few times
    for (int n = 0; n < 100; n++) {
      rng_uniform();
    }
#ifdef _OPENMP
  }
#endif

  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> globals::ntstep;  // number of time steps
  assert_always(globals::ntstep > 0);

  assert_always(get_noncommentline(file, line));
  // printout("line %s\n", line.c_str());
  std::stringstream(line) >> globals::itstep >> globals::ftstep;  // number of start and end time step
  printout("input: itstep %d ftstep %d\n", globals::itstep, globals::ftstep);
  assert_always(globals::itstep < globals::ntstep);
  assert_always(globals::itstep <= globals::ftstep);
  assert_always(globals::ftstep <= globals::ntstep);

  double tmin_days = 0.;
  double tmax_days = 0.;
  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> tmin_days >> tmax_days;  // start and end times
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
  std::stringstream(line) >> dum1;
  if (dum1 == 1) {
    grid::set_model_type(GRID_SPHERICAL1D);
  } else if (dum1 == 2) {
    grid::set_model_type(GRID_CYLINDRICAL2D);
  } else if (dum1 == 3) {
    grid::set_model_type(GRID_CARTESIAN3D);
  }

  assert_always(get_noncommentline(file, line));  // UNUSED compute the r-light curve?

  assert_always(get_noncommentline(file, line));  // UNUSED number of iterations

  assert_always(get_noncommentline(file, line));  // UNUSED change speed of light

  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> globals::gamma_grey;  // use grey opacity for gammas?

  float syn_dir_in[3];
  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> syn_dir_in[0] >> syn_dir_in[1] >> syn_dir_in[2];  // components of syn_dir

  const double rr = (syn_dir_in[0] * syn_dir_in[0]) + (syn_dir_in[1] * syn_dir_in[1]) + (syn_dir_in[2] * syn_dir_in[2]);
  // ensure that this vector is normalised.
  if (rr > 1.e-6) {
    globals::syn_dir[0] = syn_dir_in[0] / sqrt(rr);
    globals::syn_dir[1] = syn_dir_in[1] / sqrt(rr);
    globals::syn_dir[2] = syn_dir_in[2] / sqrt(rr);
  } else {
    const double z1 = 1. - (2 * rng_uniform());
    const double z2 = rng_uniform() * 2.0 * PI;
    globals::syn_dir[2] = z1;
    globals::syn_dir[0] = sqrt((1. - (z1 * z1))) * cos(z2);
    globals::syn_dir[1] = sqrt((1. - (z1 * z1))) * sin(z2);
  }

  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> globals::opacity_case;  // opacity choice

  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> globals::rho_crit_para;  // free parameter for calculation of rho_crit
  printout("input: rho_crit_para %g\n", globals::rho_crit_para);
  /// he calculation of rho_crit itself depends on the time, therfore it happens in grid_init and update_grid

  assert_always(get_noncommentline(file, line));
  int debug_packet = 0;
  std::stringstream(line) >> debug_packet;  // activate debug output for packet
  assert_always(debug_packet == -1);
  // select a negative value to deactivate

  // Do we start a new simulation or, continue another one?
  int continue_flag = 0;
  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> continue_flag;
  globals::simulation_continued_from_saved = (continue_flag == 1);
  if (globals::simulation_continued_from_saved) {
    printout("input: resuming simulation from saved point\n");
  } else {
    printout("input: starting a new simulation\n");
    assert_always(globals::itstep == 0);
  }

  /// Wavelength (in Angstroms) at which the parameterisation of the radiation field
  /// switches from the nebular approximation to LTE.
  float dum2 = NAN;
  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> dum2;  // free parameter for calculation of rho_crit
  globals::nu_rfcut = CLIGHT / (dum2 * 1e-8);
  printout("input: nu_rfcut %g\n", globals::nu_rfcut);

  /// Sets the number of initial LTE timesteps for NLTE runs
  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> globals::num_lte_timesteps;
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

  if (!USE_LUT_PHOTOION) {
    printout(
        "Corrphotoioncoeff is calculated from the radiation field at each timestep in each modelgrid cell (no "
        "LUT).\n");
  } else {
    printout(
        "Corrphotoioncoeff is calculated from LTE lookup tables (ratecoeff.dat) and corrphotoionrenorm "
        "estimator.\n");
  }

  if (!USE_LUT_BFHEATING) {
    printout(
        "bfheating coefficients are calculated from the radiation field at each timestep in each modelgrid cell (no "
        "LUT).\n");
  } else {
    printout("bfheating coefficients are calculated from LTE lookup tables (ratecoeff.dat) and bfheatingestimator.\n");
  }

  /// Set up initial grey approximation?
  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> globals::cell_is_optically_thick >> globals::num_grey_timesteps;
  printout(
      "input: cells with Thomson optical depth > %g are treated in grey approximation for the first %d timesteps\n",
      globals::cell_is_optically_thick, globals::num_grey_timesteps);

  /// Limit the number of bf-continua
  assert_always(get_noncommentline(file, line));
  int max_bf_continua = 0;
  std::stringstream(line) >> max_bf_continua;
  assert_always(max_bf_continua == -1);

  /// for exspec: read number of MPI tasks
  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> globals::nprocs_exspec;

  /// Extract line-of-sight dependent information of last emission for spectrum_res
  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> dum1;
  globals::do_emission_res = (dum1 != 0);

  /// To reduce the work imbalance between different MPI tasks I introduced a diffusion
  /// for kpkts, since it turned out that this work imbalance was largely dominated
  /// by continuous collisional interactions. By introducing a diffusion time for kpkts
  /// this loop is broken. The following two parameters control this approximation.
  /// Parameter one (a float) gives the relative fraction of a time step which individual
  /// kpkts live. Parameter two (an int) gives the number of time steps for which we
  /// want to use this approximation
  assert_always(get_noncommentline(file, line));
  std::stringstream(line) >> globals::kpktdiffusion_timescale >> globals::n_kpktdiffusion_timesteps;
  printout("input: kpkts diffuse %g of a time step's length for the first %d time steps\n",
           globals::kpktdiffusion_timescale, globals::n_kpktdiffusion_timesteps);

  file.close();

  if (rank == 0 && !globals::simulation_continued_from_saved) {
    // back up original input file, adding comments to each line
    update_parameterfile(-1);
  }
}

void update_parameterfile(int nts)
/// Subroutine to read in input parameters from input.txt.
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

      // if (!preceeding_comment && noncomment_linenum < inputlinecomments.size() - 1)
      // {
      //   fileout << '#' << inputlinecomments[noncomment_linenum] << '\n';
      // }

      // overwrite particular lines to enable restarting from the current timestep
      if (nts >= 0) {
        if (noncomment_linenum == 2) {
          /// Number of start and end time step
          snprintf(c_line, 1024, "%d %d", nts, globals::ftstep);
          // line.assign(c_line);
          line.replace(line.begin(), line.end(), c_line);
        } else if (noncomment_linenum == 16) {
          /// resume from gridsave file
          snprintf(c_line, 1024, "%d", 1);  /// Force continuation
          line.assign(c_line);
        }
      }

      if (noncomment_linenum == 21) {
        /// by default, exspec should use all available packet files
        globals::nprocs_exspec = globals::nprocs;
        snprintf(c_line, 1024, "%d", globals::nprocs_exspec);
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

void time_init()
// Subroutine to define the time steps.
{
  /// t=globals::tmin is the start of the calcualtion. t=globals::tmax is the end of the calculation.
  /// globals::ntstep is the number of time steps wanted.

  globals::time_step = static_cast<struct time *>(malloc((globals::ntstep + 1) * sizeof(struct time)));

  /// Now set the individual time steps
  switch (TIMESTEP_SIZE_METHOD) {
    case TIMESTEP_SIZES_LOGARITHMIC: {
      for (int n = 0; n < globals::ntstep; n++) {  // For logarithmic steps, the logarithmic inverval will be
        const double dlogt = (log(globals::tmax) - log(globals::tmin)) / globals::ntstep;
        globals::time_step[n].start = globals::tmin * exp(n * dlogt);
        globals::time_step[n].mid = globals::tmin * exp((n + 0.5) * dlogt);
        globals::time_step[n].width = (globals::tmin * exp((n + 1) * dlogt)) - globals::time_step[n].start;
      }
      break;
    }

    case TIMESTEP_SIZES_CONSTANT: {
      for (int n = 0; n < globals::ntstep; n++) {
        // for constant timesteps
        const double dt = (globals::tmax - globals::tmin) / globals::ntstep;
        globals::time_step[n].start = globals::tmin + n * dt;
        globals::time_step[n].width = dt;
        globals::time_step[n].mid = globals::time_step[n].start + 0.5 * globals::time_step[n].width;
      }
      break;
    }

    case TIMESTEP_SIZES_LOGARITHMIC_THEN_CONSTANT: {
      // First part log, second part fixed timesteps
      const double t_transition = TIMESTEP_TRANSITION_TIME * DAY;  // transition from logarithmic to fixed timesteps
      const double maxtsdelta = FIXED_TIMESTEP_WIDTH * DAY;        // maximum timestep width in fixed part
      assert_always(t_transition > globals::tmin);
      assert_always(t_transition < globals::tmax);
      const int nts_fixed = ceil((globals::tmax - t_transition) / maxtsdelta);
      const double fixed_tsdelta = (globals::tmax - t_transition) / nts_fixed;
      assert_always(nts_fixed > 0);
      assert_always(nts_fixed < globals::ntstep);
      const int nts_log = globals::ntstep - nts_fixed;
      assert_always(nts_log > 0);
      assert_always(nts_log < globals::ntstep);
      assert_always((nts_log + nts_fixed) == globals::ntstep);
      for (int n = 0; n < globals::ntstep; n++) {
        if (n < nts_log) {
          // For logarithmic steps, the logarithmic inverval will be
          const double dlogt = (log(t_transition) - log(globals::tmin)) / nts_log;
          globals::time_step[n].start = globals::tmin * exp(n * dlogt);
          globals::time_step[n].mid = globals::tmin * exp((n + 0.5) * dlogt);
          globals::time_step[n].width = (globals::tmin * exp((n + 1) * dlogt)) - globals::time_step[n].start;
        } else {
          // for constant timesteps
          const double prev_start =
              n > 0 ? (globals::time_step[n - 1].start + globals::time_step[n - 1].width) : globals::tmin;
          globals::time_step[n].start = prev_start;
          globals::time_step[n].width = fixed_tsdelta;
          globals::time_step[n].mid = globals::time_step[n].start + 0.5 * globals::time_step[n].width;
        }
      }
      break;
    }

    case TIMESTEP_SIZES_CONSTANT_THEN_LOGARITHMIC: {
      // /// First part fixed timesteps, second part log timesteps
      const double t_transition = TIMESTEP_TRANSITION_TIME * DAY;  // transition from fixed to logarithmic timesteps
      const double maxtsdelta = FIXED_TIMESTEP_WIDTH * DAY;        // timestep width of fixed timesteps
      assert_always(t_transition > globals::tmin);
      assert_always(t_transition < globals::tmax);
      const int nts_fixed = ceil((t_transition - globals::tmin) / maxtsdelta);
      const double fixed_tsdelta = (t_transition - globals::tmin) / nts_fixed;
      assert_always(nts_fixed > 0);
      assert_always(nts_fixed < globals::ntstep);
      const int nts_log = globals::ntstep - nts_fixed;
      assert_always(nts_log > 0);
      assert_always(nts_log < globals::ntstep);
      assert_always((nts_log + nts_fixed) == globals::ntstep);
      for (int n = 0; n < globals::ntstep; n++) {
        if (n < nts_fixed) {
          // for constant timesteps
          globals::time_step[n].start = globals::tmin + n * fixed_tsdelta;
          globals::time_step[n].width = fixed_tsdelta;
          globals::time_step[n].mid = globals::time_step[n].start + 0.5 * globals::time_step[n].width;
        } else {
          // For logarithmic time steps, the logarithmic interval will be
          const double dlogt = (log(globals::tmax) - log(t_transition)) / nts_log;
          const double prev_start =
              n > 0 ? (globals::time_step[n - 1].start + globals::time_step[n - 1].width) : globals::tmin;
          globals::time_step[n].start = prev_start;
          globals::time_step[n].width = (t_transition * exp((n - nts_fixed + 1) * dlogt)) - globals::time_step[n].start;
          globals::time_step[n].mid = globals::time_step[n].start + 0.5 * globals::time_step[n].width;
        }
      }
      break;
    }

    default:
      assert_always(false);
  }

  // to limit the timestep durations
  // const double maxt = 0.5 * DAY;
  // for (int n = globals::ntstep - 1; n > 0; n--)
  // {
  //   if (globals::time_step[n].width > maxt)
  //   {
  //     const double boundaryshift = globals::time_step[n].width - maxt;
  //     globals::time_step[n].width -= boundaryshift;
  //     globals::time_step[n].start += boundaryshift;
  //     globals::time_step[n - 1].width += boundaryshift;
  //   }
  //   else if (n < globals::ntstep - 1 && globals::time_step[n + 1].width > maxt)
  //   {
  //     printout("TIME: Keeping logarithmic durations for timesteps <= %d\n", n);
  //   }
  // }
  // assert_always(globals::time_step[0].width <= maxt); // no solution is possible with these constraints!

  /// and add a dummy timestep which contains the endtime
  /// of the calculation
  globals::time_step[globals::ntstep].start = globals::tmax;
  globals::time_step[globals::ntstep].mid = globals::tmax;
  globals::time_step[globals::ntstep].width = 0.;

  // check consistency of start + width = start_next
  for (int n = 1; n < globals::ntstep; n++) {
    assert_always(
        fabs((globals::time_step[n - 1].start + globals::time_step[n - 1].width) / globals::time_step[n].start) - 1 <
        0.001);
  }
  assert_always(fabs((globals::time_step[globals::ntstep - 1].start + globals::time_step[globals::ntstep - 1].width) /
                     globals::tmax) -
                    1 <
                0.001);

  for (int n = 0; n < globals::ntstep; n++) {
    globals::time_step[n].positron_dep = 0.;
    globals::time_step[n].eps_positron_ana_power = 0.;
    globals::time_step[n].electron_dep = 0.;
    globals::time_step[n].electron_emission = 0.;
    globals::time_step[n].eps_electron_ana_power = 0.;
    globals::time_step[n].alpha_dep = 0.;
    globals::time_step[n].alpha_emission = 0.;
    globals::time_step[n].eps_alpha_ana_power = 0.;
    globals::time_step[n].gamma_dep = 0.;
    globals::time_step[n].gamma_dep_pathint = 0.;
    globals::time_step[n].qdot_betaminus = 0.;
    globals::time_step[n].qdot_alpha = 0.;
    globals::time_step[n].qdot_total = 0.;
    globals::time_step[n].gamma_emission = 0.;
    globals::time_step[n].cmf_lum = 0.0;
    globals::time_step[n].pellet_decays = 0;
  }
}

void write_timestep_file() {
  FILE *timestepfile = fopen_required("timesteps.out", "w");
  fprintf(timestepfile, "#timestep tstart_days tmid_days twidth_days\n");
  for (int n = 0; n < globals::ntstep; n++) {
    fprintf(timestepfile, "%d %lg %lg %lg\n", n, globals::time_step[n].start / DAY, globals::time_step[n].mid / DAY,
            globals::time_step[n].width / DAY);
  }
  fclose(timestepfile);
}
