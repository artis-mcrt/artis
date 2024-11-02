#include "radfield.h"

#include <cstddef>

#ifdef MPI_ON
#include <mpi.h>
#endif

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_debye.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <vector>

#include "artisoptions.h"
#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "grid.h"
#include "ratecoeff.h"
#include "rpkt.h"
#include "sn3d.h"

namespace radfield {

namespace {

std::vector<double> J_normfactor;

struct RadFieldBinSolution {
  // these two parameters are used in the current timestep, but were calculated
  // from the values of J and nuJ in the previous timestep
  float W;    // dilution (scaling) factor
  float T_R;  // radiation temperature
};

struct RadFieldBin {
  double J_raw;  // value needs to be multiplied by J_normfactor to get the true value
  double nuJ_raw;
  int contribcount;
};

constexpr double radfieldbins_delta_nu =
    (nu_upper_last_initial - nu_lower_first_initial) / (RADFIELDBINCOUNT - 1);  // - 1 for the top super bin

RadFieldBin *radfieldbins{};
RadFieldBinSolution *radfieldbin_solutions{};

#ifdef MPI_ON
MPI_Win win_radfieldbin_solutions = MPI_WIN_NULL;
MPI_Win win_prev_bfrate_normed = MPI_WIN_NULL;
#endif

struct Jb_lu_estimator {
  double value = 0.;
  int contribcount = 0;
};

// reallocate the detailed line arrays in units of BLOCKSIZEJBLUE
constexpr int BLOCKSIZEJBLUE = 128;
int detailed_linecount = 0;

// array of indices into the linelist[] array for selected lines
int *detailed_lineindicies;

Jb_lu_estimator **prev_Jb_lu_normed{};  // value from the previous timestep
Jb_lu_estimator **Jb_lu_raw{};          // unnormalised estimator for the current timestep

float *prev_bfrate_normed{};     // values from the previous timestep
std::vector<double> bfrate_raw;  // unnormalised estimators for the current timestep

// expensive debugging mode to track the contributions to each bound-free rate estimator

std::vector<double> J;  // after normalisation: [ergs/s/sr/cm2/Hz]
#ifdef DO_TITER
std::vector<double> J_reduced_save;
#endif

// J and nuJ are accumulated and then normalised in-place
// i.e. be sure the normalisation has been applied (exactly once) before using the values here!
std::vector<double> nuJ;
#ifdef DO_TITER
std::vector<double> nuJ_reduced_save;
#endif

struct gsl_planck_integral_paras {
  double T_R;
  bool times_nu;
};

struct gsl_T_R_solver_paras {
  int modelgridindex;
  int binindex;
};

FILE *radfieldfile{};

constexpr auto get_bin_nu_upper(const int binindex) -> double {
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  if (binindex == RADFIELDBINCOUNT - 1) {
    return nu_upper_superbin;
  }
  return nu_lower_first_initial + ((binindex + 1) * radfieldbins_delta_nu);
}

constexpr auto get_bin_nu_lower(const int binindex) -> double {
  if (binindex > 0) {
    return get_bin_nu_upper(binindex - 1);
  }
  return nu_lower_first_initial;
}

// find the left-closed bin [nu_lower, nu_upper) that nu belongs to
#pragma omp declare simd
constexpr auto select_bin(const double nu) -> int {
  if (nu < nu_lower_first_initial) {
    return -2;  // out of range, nu lower than lowest bin's lower boundary
  }
  if (nu >= nu_upper_superbin) {
    // out of range, nu higher than highest bin's upper boundary
    return -1;
  }
  if (nu >= nu_upper_last_initial) {
    // in the superbin. separate case because the delta_nu is different to the other bins
    return RADFIELDBINCOUNT - 1;
  }

  const int binindex = static_cast<int>((nu - nu_lower_first_initial) / radfieldbins_delta_nu);

  if (nu == get_bin_nu_upper(binindex)) {
    // exactly on the upper boundary of the bin, so add 1 to ensure we get the left-closed bin
    return binindex + 1;
  }

  return binindex;
}

void realloc_detailed_lines(const int new_size) {
  auto *newptr = static_cast<int *>(realloc(detailed_lineindicies, new_size * sizeof(int)));
  if (newptr == nullptr) {
    printout("ERROR: Not enough memory to reallocate detailed Jblue estimator line list\n");
    std::abort();
  }
  assert_always(newptr != nullptr);
  detailed_lineindicies = newptr;

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    if (grid::get_numpropcells(modelgridindex) > 0) {
      prev_Jb_lu_normed[modelgridindex] = static_cast<Jb_lu_estimator *>(
          realloc(prev_Jb_lu_normed[modelgridindex], new_size * sizeof(Jb_lu_estimator)));

      Jb_lu_raw[modelgridindex] =
          static_cast<Jb_lu_estimator *>(realloc(Jb_lu_raw[modelgridindex], new_size * sizeof(Jb_lu_estimator)));

      if (prev_Jb_lu_normed[modelgridindex] == nullptr || Jb_lu_raw[modelgridindex] == nullptr) {
        printout("ERROR: Not enough memory to reallocate detailed Jblue estimator list for cell %d.\n", modelgridindex);
        std::abort();
      }
    }
  }
}

// associate a Jb_lu estimator with a particular lineindex to be used
// instead of the general radiation field model
void add_detailed_line(const int lineindex) {
  if (detailed_linecount % BLOCKSIZEJBLUE == 0) {
    const int new_size = detailed_linecount + BLOCKSIZEJBLUE;
    realloc_detailed_lines(new_size);
  }

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    if (grid::get_numpropcells(modelgridindex) > 0) {
      prev_Jb_lu_normed[modelgridindex][detailed_linecount].value = 0;
      prev_Jb_lu_normed[modelgridindex][detailed_linecount].contribcount = 0;

      // zero_estimators should do the next part anyway, but just to be sure:
      Jb_lu_raw[modelgridindex][detailed_linecount].value = 0;
      Jb_lu_raw[modelgridindex][detailed_linecount].contribcount = 0;
    }
  }
  detailed_lineindicies[detailed_linecount] = lineindex;
  detailed_linecount++;
}

// get the normalised J_nu
auto get_bin_J(const int modelgridindex, const int binindex) -> double {
  const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  assert_testmodeonly(J_normfactor[nonemptymgi] > 0.0);
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  return radfieldbins[(nonemptymgi * RADFIELDBINCOUNT) + binindex].J_raw * J_normfactor[nonemptymgi];
}

auto get_bin_nuJ(const int modelgridindex, const int binindex) -> double {
  const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  assert_testmodeonly(J_normfactor[nonemptymgi] > 0.0);
  assert_testmodeonly(modelgridindex < grid::get_npts_model());
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  return radfieldbins[(nonemptymgi * RADFIELDBINCOUNT) + binindex].nuJ_raw * J_normfactor[nonemptymgi];
}

// get <nuJ> / <J> for a bin
auto get_bin_nu_bar(const int modelgridindex, const int binindex) -> double {
  const double nuJ_sum = get_bin_nuJ(modelgridindex, binindex);
  const double J_sum = get_bin_J(modelgridindex, binindex);
  return nuJ_sum / J_sum;
}

auto get_bin_contribcount(const int modelgridindex, const int binindex) -> int {
  const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  return radfieldbins[(nonemptymgi * RADFIELDBINCOUNT) + binindex].contribcount;
}

auto get_bin_W(const int modelgridindex, const int binindex) -> float {
  const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  return radfieldbin_solutions[(nonemptymgi * RADFIELDBINCOUNT) + binindex].W;
}

auto get_bin_T_R(const int modelgridindex, const int binindex) -> float {
  const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  return radfieldbin_solutions[(nonemptymgi * RADFIELDBINCOUNT) + binindex].T_R;
}

constexpr auto gsl_integrand_planck(const double nu, void *voidparas) -> double {
  const auto *paras = static_cast<gsl_planck_integral_paras *>(voidparas);
  const auto T_R = paras->T_R;

  double integrand = TWOHOVERCLIGHTSQUARED * std::pow(nu, 3) / (std::expm1(HOVERKB * nu / T_R));

  if (paras->times_nu) {
    integrand *= nu;
  }

  return integrand;
}

void update_bfestimators(const int nonemptymgi, const double distance_e_cmf, const double nu_cmf,
                         const double doppler_nucmf_on_nurf, const Phixslist &phixslist) {
  assert_testmodeonly(DETAILED_BF_ESTIMATORS_ON);

  const double distance_e_cmf_over_nu =
      distance_e_cmf / nu_cmf * doppler_nucmf_on_nurf;  // TODO: Luke: why did I put a doppler factor here?

  // I think the nu_cmf slightly differs from when the phixslist was calculated
  // so the nu condition on this nu_cmf can truncate the list further compared to what was used in the calculation
  // of phixslist.gamma_contr

  const auto bfestimend = std::upper_bound(globals::bfestim_nu_edge.data(),
                                           globals::bfestim_nu_edge.data() + phixslist.bfestimend, nu_cmf) -
                          globals::bfestim_nu_edge.data();

  const auto bfestimbegin = std::lower_bound(globals::bfestim_nu_edge.data() + phixslist.bfestimbegin,
                                             globals::bfestim_nu_edge.data() + bfestimend, nu_cmf,
                                             [](const double nu_edge, const double find_nu_cmf) {
                                               return nu_edge * last_phixs_nuovernuedge < find_nu_cmf;
                                             }) -
                            globals::bfestim_nu_edge.data();

  const auto bfestimcount = globals::bfestimcount;
#pragma omp simd
  for (auto bfestimindex = bfestimbegin; bfestimindex < bfestimend; bfestimindex++) {
    atomicadd(bfrate_raw[(nonemptymgi * bfestimcount) + bfestimindex],
              phixslist.gamma_contr[bfestimindex] * distance_e_cmf_over_nu);
  }
}

auto planck_integral(const double T_R, const double nu_lower, const double nu_upper, const bool times_nu) -> double {
  double integral = 0.;

  double error = 0.;
  const double epsrel = 1e-10;
  const double epsabs = 0.;

  const gsl_planck_integral_paras intparas = {.T_R = T_R, .times_nu = times_nu};

  gsl_error_handler_t *previous_handler = gsl_set_error_handler(gsl_error_handler_printout);

  const int status = integrator<gsl_integrand_planck>(intparas, nu_lower, nu_upper, epsabs, epsrel, GSL_INTEG_GAUSS61,
                                                      &integral, &error);
  if (status != 0) {
    printout("planck_integral integrator status %d, GSL_FAILURE= %d. Integral value %g, setting to zero.\n", status,
             GSL_FAILURE, integral);
    integral = 0.;
  }
  gsl_set_error_handler(previous_handler);

  return integral;
}

auto delta_nu_bar(const double T_R, void *const paras) -> double
// difference between the average nu and the average nu of a Planck function
// at temperature T_R, in the frequency range corresponding to a bin
{
  const auto *params = static_cast<const gsl_T_R_solver_paras *>(paras);
  const int modelgridindex = params->modelgridindex;
  const int binindex = params->binindex;

  const double nu_lower = get_bin_nu_lower(binindex);
  const double nu_upper = get_bin_nu_upper(binindex);

  const double nu_bar_estimator = get_bin_nu_bar(modelgridindex, binindex);

  const double nu_times_planck_numerical = planck_integral(T_R, nu_lower, nu_upper, true);
  const double planck_integral_numerical = planck_integral(T_R, nu_lower, nu_upper, false);
  const double nu_bar_planck_T_R = nu_times_planck_numerical / planck_integral_numerical;

  // double nu_times_planck_integral = planck_integral_analytic(T_R, nu_lower, nu_upper, true);
  // double planck_integral_result = planck_integral_analytic(T_R, nu_lower, nu_upper, false);
  // double nu_bar_planck = nu_times_planck_integral / planck_integral_result;

  // // printout("nu_bar %g nu_bar_planck(T=%g) %g\n",nu_bar,T_R,nu_bar_planck);

  // if (!std::isfinite(nu_bar_planck)) {
  //   double nu_times_planck_numerical = planck_integral(T_R, nu_lower, nu_upper, true);
  //   double planck_integral_numerical = planck_integral(T_R, nu_lower, nu_upper, false);
  //   double nu_bar_planck_numerical = nu_times_planck_numerical / planck_integral_numerical;

  //   printout("planck_integral_analytic is %g. Replacing with numerical result of %g.\n", nu_bar_planck,
  //            nu_bar_planck_numerical);
  //   nu_bar_planck = nu_bar_planck_numerical;
  // }

  const double delta_nu_bar = nu_bar_planck_T_R - nu_bar_estimator;

  if (!std::isfinite(delta_nu_bar)) {
    printout(
        "delta_nu_bar is %g. nu_bar_planck_T_R %g nu_times_planck_numerical %g planck_integral_numerical %g "
        "nu_bar_estimator %g\n",
        delta_nu_bar, nu_bar_planck_T_R, nu_times_planck_numerical, planck_integral_numerical, nu_bar_estimator);
  }

  return delta_nu_bar;
}

auto find_T_R(const int modelgridindex, const int binindex) -> float {
  double T_R = 0.;

  gsl_T_R_solver_paras paras{};
  paras.modelgridindex = modelgridindex;
  paras.binindex = binindex;

  // Check whether the equation has a root in [T_min,T_max]
  double delta_nu_bar_min = delta_nu_bar(T_R_min, &paras);
  double delta_nu_bar_max = delta_nu_bar(T_R_max, &paras);

  // printout("find_T_R: bin %4d delta_nu_bar(T_R_min) %g, delta_nu_bar(T_R_max) %g\n",
  //          binindex, delta_nu_bar_min,delta_nu_bar_max);

  if (!std::isfinite(delta_nu_bar_min) || !std::isfinite(delta_nu_bar_max)) {
    delta_nu_bar_max = delta_nu_bar_min = -1;
  }

  if (delta_nu_bar_min * delta_nu_bar_max < 0) {
    // If there is a root in the interval, solve for T_R

    const double epsrel = 1e-4;
    const double epsabs = 0.;
    const int maxit = 100;

    gsl_function find_T_R_f = {.function = &delta_nu_bar, .params = &paras};

    // one dimensional gsl root solver, bracketing type
    gsl_root_fsolver *T_R_solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(T_R_solver, &find_T_R_f, T_R_min, T_R_max);
    int status = 0;
    for (int iteration_num = 0; iteration_num <= maxit; iteration_num++) {
      gsl_root_fsolver_iterate(T_R_solver);
      T_R = gsl_root_fsolver_root(T_R_solver);

      const double T_R_lower = gsl_root_fsolver_x_lower(T_R_solver);
      const double T_R_upper = gsl_root_fsolver_x_upper(T_R_solver);
      status = gsl_root_test_interval(T_R_lower, T_R_upper, epsabs, epsrel);

      // printout("find_T_R: bin %4d iter %d, T_R is between %7.1f and %7.1f, guess %7.1f, delta_nu_bar %g, status
      // %d\n",
      //          binindex,iteration_num,T_R_lower,T_R_upper,T_R,delta_nu_bar(T_R,&paras),status);
      if (status != GSL_CONTINUE) {
        break;
      }
    }

    if (status == GSL_CONTINUE) {
      printout("[warning] find_T_R: T_R did not converge within %d iterations\n", maxit);
    }

    gsl_root_fsolver_free(T_R_solver);
  } else if (delta_nu_bar_max < 0) {
    // Thermal balance equation always negative ===> T_R = T_min
    // Calculate the rates again at this T_e to print them to file
    T_R = T_R_max;
    printout("find_T_R: cell %d bin %4d no solution in interval, clamping to T_R_max=%g\n", modelgridindex, binindex,
             T_R_max);
  } else {
    T_R = T_R_min;
    printout("find_T_R: cell %d bin %4d no solution in interval, clamping to T_R_min=%g\n", modelgridindex, binindex,
             T_R_min);
  }

  return T_R;
}  // namespace radfield

void set_params_fullspec(const int modelgridindex, const int timestep) {
  const int nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  const double nubar = nuJ[nonemptymgi] / J[nonemptymgi];
  if (!std::isfinite(nubar) || nubar == 0.) {
    printout("[warning] T_R estimator infinite in cell %d, keep T_R, T_J, W of last timestep. J = %g. nuJ = %g\n",
             modelgridindex, J[nonemptymgi], nuJ[nonemptymgi]);
  } else {
    float T_J = pow(J[nonemptymgi] * PI / STEBO, 1 / 4.);
    if (T_J > MAXTEMP) {
      printout("[warning] temperature estimator T_J = %g exceeds T_max %g in cell %d. Setting T_J = T_max!\n", T_J,
               MAXTEMP, modelgridindex);
      T_J = MAXTEMP;
    } else if (T_J < MINTEMP) {
      printout("[warning] temperature estimator T_J = %g below T_min %g in cell %d. Setting T_J = T_min!\n", T_J,
               MINTEMP, modelgridindex);
      T_J = MINTEMP;
    }
    grid::set_TJ(modelgridindex, T_J);

    float T_R = H * nubar / KB / 3.832229494;
    if (T_R > MAXTEMP) {
      printout("[warning] temperature estimator T_R = %g exceeds T_max %g in cell %d. Setting T_R = T_max!\n", T_R,
               MAXTEMP, modelgridindex);
      T_R = MAXTEMP;
    } else if (T_R < MINTEMP) {
      printout("[warning] temperature estimator T_R = %g below T_min %g in cell %d. Setting T_R = T_min!\n", T_R,
               MINTEMP, modelgridindex);
      T_R = MINTEMP;
    }
    grid::set_TR(modelgridindex, T_R);

    const float W = J[nonemptymgi] * PI / STEBO / pow(T_R, 4);
    grid::set_W(modelgridindex, W);

    printout(
        "Full-spectrum fit radfield for cell %d at timestep %d: J %g, nubar %5.1f Angstrom, T_J %g, T_R %g, W %g\n",
        modelgridindex, timestep, J[nonemptymgi], 1e8 * CLIGHT / nubar, T_J, T_R, W);
  }
}

auto get_bfcontindex(const int element, const int lowerion, const int lower, const int phixstargetindex) -> int {
  // simple linear search seems to be faster than the binary search
  // possibly because lower frequency transitions near start of list are more likely to be called?
  const auto bfcontindex = static_cast<int>(std::find_if(globals::allcont, globals::allcont + globals::nbfcontinua,
                                                         [=](const auto &bf) {
                                                           return (bf.element == element) && (bf.ion == lowerion) &&
                                                                  (bf.level == lower) &&
                                                                  (bf.phixstargetindex == phixstargetindex);
                                                         }) -
                                            globals::allcont);

  if (bfcontindex < globals::nbfcontinua) {
    return bfcontindex;
  }
  // not found in the continua list
  return -1;
}

}  // anonymous namespace

void init(const int my_rank, const int ndo_nonempty) {
  // this should be called only after the atomic data is in memory

  const ptrdiff_t nonempty_npts_model = grid::get_nonempty_npts_model();

  J_normfactor.resize(nonempty_npts_model + 1);
  J.resize(nonempty_npts_model + 1);

#ifdef DO_TITER
  J_reduced_save.resize(nonempty_npts_model + 1);
#endif

  // J and nuJ are accumulated and then normalised in-place
  // i.e. be sure the normalisation has been applied (exactly once) before using the values here!
  nuJ.resize(nonempty_npts_model + 1);
#ifdef DO_TITER
  nuJ.resize(nonempty_npts_model + 1);
#endif

  prev_Jb_lu_normed = static_cast<Jb_lu_estimator **>(malloc((grid::get_npts_model() + 1) * sizeof(Jb_lu_estimator *)));
  Jb_lu_raw = static_cast<Jb_lu_estimator **>(malloc((grid::get_npts_model() + 1) * sizeof(Jb_lu_estimator *)));

  detailed_linecount = 0;

  detailed_lineindicies = nullptr;
  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    prev_Jb_lu_normed[modelgridindex] = nullptr;
    Jb_lu_raw[modelgridindex] = nullptr;
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    for (int i = 0; i < globals::nlines; i++) {
      const int element = globals::linelist[i].elementindex;
      const int Z = get_atomicnumber(element);
      if (Z == 26) {
        const int lowerlevel = globals::linelist[i].lowerlevelindex;
        // const int upperlevel = linelist[i].upperlevelindex;
        // const int ion = linelist[i].ionindex;
        // const int ionstage = get_ionstage(element, ion);
        const double A_ul = globals::linelist[i].einstein_A;

        bool addline = false;
        // if (ionstage == 1 && lowerlevel == 6 && upperlevel == 55)
        //   addline = true;
        // else if (ionstage == 1 && lowerlevel == 10 && upperlevel == 104)
        //   addline = true;
        // else if (ionstage == 1 && lowerlevel == 10 && upperlevel == 112)
        //   addline = true;
        // else if (ionstage == 2 && lowerlevel == 9 && upperlevel == 64)
        //   addline = true;

        if (lowerlevel <= 15 && A_ul > 0.) {  // ionstage <= 3 && A_ul > 1e3 &&
          addline = true;
        }

        if (addline) {
          // printout("Adding Jblue estimator for lineindex %d Z=%02d ionstage %d lower %d upper %d A_ul %g\n",
          //          i, Z, ionstage, lowerlevel, upperlevel, A_ul);
          add_detailed_line(i);
        }
      }
    }

    // shrink the detailed line list in case detailed_linecount isn't a multiple of BLOCKSIZEJBLUE
    // (important for saving memory if there are a lot of grid cells)
    realloc_detailed_lines(detailed_linecount);

    // these are probably sorted anyway because the previous loop goes in ascending
    // lineindex. But this sorting step is quick and makes sure that the
    // binary searching later will work correctly
    std::SORT_OR_STABLE_SORT(detailed_lineindicies, detailed_lineindicies + detailed_linecount);
  }

  printout("There are %d lines with detailed Jblue_lu estimators.\n", detailed_linecount);

  printout("DETAILED_BF_ESTIMATORS %s", DETAILED_BF_ESTIMATORS_ON ? "ON" : "OFF");
  if (DETAILED_BF_ESTIMATORS_ON) {
    printout(" from timestep %d\n", DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP);
  } else {
    printout("\n");
  }

  if (MULTIBIN_RADFIELD_MODEL_ON) {
    printout("The multibin radiation field is being used from timestep %d onwards.\n", FIRST_NLTE_RADFIELD_TIMESTEP);

    printout("Initialising multibin radiation field with %d bins from (%.2f eV, %6.1f A) to (%.2f eV, %6.1f A)\n",
             RADFIELDBINCOUNT, H * nu_lower_first_initial / EV, 1e8 * CLIGHT / nu_lower_first_initial,
             H * nu_upper_last_initial / EV, 1e8 * CLIGHT / nu_upper_last_initial);
    if (ndo_nonempty > 0) {
      char filename[MAXFILENAMELENGTH];
      snprintf(filename, MAXFILENAMELENGTH, "radfield_%.4d.out", my_rank);
      assert_always(radfieldfile == nullptr);
      radfieldfile = fopen_required(filename, "w");
      fprintf(radfieldfile, "timestep modelgridindex bin_num nu_lower nu_upper nuJ J J_nu_avg ncontrib T_R W\n");
      fflush(radfieldfile);
    }

    const size_t mem_usage_bins = nonempty_npts_model * RADFIELDBINCOUNT * sizeof(RadFieldBin);
    radfieldbins = static_cast<RadFieldBin *>(malloc(nonempty_npts_model * RADFIELDBINCOUNT * sizeof(RadFieldBin)));

    const size_t mem_usage_bin_solutions = nonempty_npts_model * RADFIELDBINCOUNT * sizeof(RadFieldBinSolution);

#ifdef MPI_ON
    std::tie(radfieldbin_solutions, win_radfieldbin_solutions) =
        MPI_shared_malloc_keepwin<RadFieldBinSolution>(nonempty_npts_model * RADFIELDBINCOUNT);
#else
    radfieldbin_solutions = static_cast<RadFieldBinSolution *>(
        malloc(nonempty_npts_model * RADFIELDBINCOUNT * sizeof(RadFieldBinSolution)));
#endif

    printout("[info] mem_usage: radiation field bin accumulators for non-empty cells occupy %.3f MB\n",
             mem_usage_bins / 1024. / 1024.);
    printout(
        "[info] mem_usage: radiation field bin solutions for non-empty cells occupy %.3f MB (node shared memory)\n",
        mem_usage_bin_solutions / 1024. / 1024.);
  } else {
    printout("The radiation field model is a full-spectrum fit to a single dilute blackbody TR & W.\n");
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    {
#ifdef MPI_ON
      std::tie(prev_bfrate_normed, win_prev_bfrate_normed) =
          MPI_shared_malloc_keepwin<float>(nonempty_npts_model * globals::bfestimcount);
#else
      prev_bfrate_normed = static_cast<float *>(malloc(nonempty_npts_model * globals::bfestimcount * sizeof(float)));
#endif
    }
    printout("[info] mem_usage: detailed bf estimators for non-empty cells occupy %.3f MB (node shared memory)\n",
             nonempty_npts_model * globals::bfestimcount * sizeof(float) / 1024. / 1024.);

    bfrate_raw.resize(nonempty_npts_model * globals::bfestimcount);

    printout("[info] mem_usage: detailed bf estimator acculumators for non-empty cells occupy %.3f MB\n",
             nonempty_npts_model * globals::bfestimcount * sizeof(double) / 1024. / 1024.);
  }

  zero_estimators();

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
#ifdef MPI_ON
    MPI_Barrier(globals::mpi_comm_node);
#endif
    if (globals::rank_in_node == 0) {
      for (ptrdiff_t nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
        for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
          const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
          radfieldbin_solutions[mgibinindex].W = -1.;
          radfieldbin_solutions[mgibinindex].T_R = -1.;
        }
      }
    }
#ifdef MPI_ON
    MPI_Barrier(globals::mpi_comm_node);
#endif
  }
}

// Initialise estimator arrays which hold the last time steps values (used to damp out
// fluctuations over timestep iterations if DO_TITER is defined) to -1.
void initialise_prev_titer_photoionestimators() {
#ifdef DO_TITER
  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    const int nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
    globals::ffheatingestimator_save[nonemptymgi] = -1.;
    globals::colheatingestimator_save[nonemptymgi] = -1.;
    J_reduced_save[nonemptymgi] = -1.;
    nuJ_reduced_save[nonemptymgi] = -1.;
    for (int element = 0; element < get_nelements(); element++) {
      const int nions = get_nions(element);
      for (int ion = 0; ion < nions - 1; ion++) {
        if constexpr (USE_LUT_PHOTOION) {
          globals::gammaestimator_save[get_ionestimindex_nonemptymgi(nonemptymgi, element, ion)] = -1.;
        }
        if constexpr (USE_LUT_BFHEATING) {
          globals::bfheatingestimator_save[get_ionestimindex_nonemptymgi(nonemptymgi, element, ion)] = -1.;
        }
      }
    }
  }
#endif
}

auto get_Jblueindex(const int lineindex) -> int {
  // returns -1 if the line does not have a Jblue estimator
  if constexpr (!DETAILED_LINE_ESTIMATORS_ON) {
    return -1;
  }

  // use a binary search, assuming the list is sorted

  int low = 0;
  int high = detailed_linecount;
  while (low <= high) {
    const int mid = low + ((high - low) / 2);
    if (detailed_lineindicies[mid] < lineindex) {
      low = mid + 1;
    } else if (detailed_lineindicies[mid] > lineindex) {
      high = mid - 1;
    } else {
      return mid;
    }
  }

  return -1;
}

auto get_Jb_lu(const int modelgridindex, const int jblueindex) -> double {
  assert_always(jblueindex >= 0);
  assert_always(jblueindex < detailed_linecount);
  return prev_Jb_lu_normed[modelgridindex][jblueindex].value;
}

auto get_Jb_lu_contribcount(const int modelgridindex, const int jblueindex) -> int {
  assert_always(jblueindex >= 0);
  assert_always(jblueindex < detailed_linecount);
  return prev_Jb_lu_normed[modelgridindex][jblueindex].contribcount;
}

void write_to_file(const int modelgridindex, const int timestep) {
  assert_always(MULTIBIN_RADFIELD_MODEL_ON);
  const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
#ifdef _OPENMP
#pragma omp critical(out_file)
  {
#endif

    int totalcontribs = 0;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      totalcontribs += get_bin_contribcount(modelgridindex, binindex);
    }

    for (int binindex = -1 - detailed_linecount; binindex < RADFIELDBINCOUNT; binindex++) {
      double nu_lower = 0.;
      double nu_upper = 0.;
      double nuJ_out = 0.;
      double J_out = 0.;
      float T_R = 0.;
      float W = 0.;
      double J_nu_bar = 0.;
      int contribcount = 0;

      const bool skipoutput = false;

      if (binindex >= 0) {
        nu_lower = get_bin_nu_lower(binindex);
        nu_upper = get_bin_nu_upper(binindex);
        nuJ_out = get_bin_nuJ(modelgridindex, binindex);
        J_out = get_bin_J(modelgridindex, binindex);
        T_R = get_bin_T_R(modelgridindex, binindex);
        W = get_bin_W(modelgridindex, binindex);
        J_nu_bar = J_out / (nu_upper - nu_lower);
        contribcount = get_bin_contribcount(modelgridindex, binindex);
      } else if (binindex == -1) {  // bin -1 is the full spectrum fit
        nuJ_out = nuJ[nonemptymgi];
        J_out = J[nonemptymgi];
        T_R = grid::get_TR(nonemptymgi);
        W = grid::get_W(nonemptymgi);
        contribcount = totalcontribs;
      } else  // use binindex < -1 for detailed line Jb_lu estimators
      {
        const int jblueindex = -2 - binindex;  // -2 is the first detailed line, -3 is the second, etc
        const int lineindex = detailed_lineindicies[jblueindex];
        const double nu_trans = globals::linelist[lineindex].nu;
        nu_lower = nu_trans;
        nu_upper = nu_trans;
        nuJ_out = -1.;
        J_out = -1.;
        T_R = -1.;
        W = -1.;
        J_nu_bar = prev_Jb_lu_normed[modelgridindex][jblueindex].value,
        contribcount = prev_Jb_lu_normed[modelgridindex][jblueindex].contribcount;
      }

      if (!skipoutput) {
        fprintf(radfieldfile, "%d %d %d %.5e %.5e %.3e %.3e %.3e %d %.1f %.5e\n", timestep, modelgridindex, binindex,
                nu_lower, nu_upper, nuJ_out, J_out, J_nu_bar, contribcount, T_R, W);
      }
    }
    fflush(radfieldfile);
#ifdef _OPENMP
  }
#endif
}

void close_file() {
  if (radfieldfile != nullptr) {
    fclose(radfieldfile);
    radfieldfile = nullptr;
  }

  if (MULTIBIN_RADFIELD_MODEL_ON) {
    free(radfieldbins);
#ifdef MPI_ON
    if (win_radfieldbin_solutions != MPI_WIN_NULL) {
      MPI_Win_free(&win_radfieldbin_solutions);
    }
#else
    if (radfieldbin_solutions != nullptr) {
      free(radfieldbin_solutions);
    }
#endif
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
#ifdef MPI_ON
    if (win_radfieldbin_solutions != MPI_WIN_NULL) {
      MPI_Win_free(&win_prev_bfrate_normed);
    }
#else
    if (prev_bfrate_normed != nullptr) {
      free(prev_bfrate_normed);
    }
#endif
  }
}

// set up the new bins and clear the estimators in preparation for a timestep
void zero_estimators() {
  std::ranges::fill(J_normfactor, -1.0);
  std::ranges::fill(J, 0.0);
  std::ranges::fill(nuJ, 0.0);
  std::ranges::fill(bfrate_raw, 0.0);

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    assert_always(radfieldbins != nullptr);
    for (ptrdiff_t nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
        const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
        radfieldbins[mgibinindex].J_raw = 0.;
        radfieldbins[mgibinindex].nuJ_raw = 0.;
        radfieldbins[mgibinindex].contribcount = 0;
      }
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
      const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
      assert_always(Jb_lu_raw != nullptr);
      assert_always(Jb_lu_raw[modelgridindex] != nullptr);
      for (int i = 0; i < detailed_linecount; i++) {
        Jb_lu_raw[modelgridindex][i].value = 0.;
        Jb_lu_raw[modelgridindex][i].contribcount = 0.;
      }
    }
  }
}

__host__ __device__ void update_estimators(const int nonemptymgi, const double distance_e_cmf, const double nu_cmf,
                                           const double doppler_nucmf_on_nurf, const Phixslist &phixslist,
                                           const bool thickcell) {
  if (distance_e_cmf == 0) {
    return;
  }

  atomicadd(J[nonemptymgi], distance_e_cmf);
  atomicadd(nuJ[nonemptymgi], distance_e_cmf * nu_cmf);

  if (thickcell) {
    return;
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    update_bfestimators(nonemptymgi, distance_e_cmf, nu_cmf, doppler_nucmf_on_nurf, phixslist);
  }

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    const int binindex = select_bin(nu_cmf);

    if (binindex >= 0) {
      const ptrdiff_t mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
      atomicadd(radfieldbins[mgibinindex].J_raw, distance_e_cmf);
      atomicadd(radfieldbins[mgibinindex].nuJ_raw, distance_e_cmf * nu_cmf);
      atomicadd(radfieldbins[mgibinindex].contribcount, 1);
    }
  }
}

__host__ __device__ void update_lineestimator(const int modelgridindex, const int lineindex, const double increment) {
  if constexpr (!DETAILED_LINE_ESTIMATORS_ON) {
    return;
  }

  const int jblueindex = get_Jblueindex(lineindex);
  if (jblueindex >= 0) {
    Jb_lu_raw[modelgridindex][jblueindex].value += increment;
    Jb_lu_raw[modelgridindex][jblueindex].contribcount += 1;
  }
}

// mean intensity J_nu [ergs/s/sr/cm2/Hz]
__host__ __device__ auto radfield(const double nu, const int nonemptymgi) -> double {
  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    if (globals::timestep >= FIRST_NLTE_RADFIELD_TIMESTEP) {
      const int binindex = select_bin(nu);
      if (binindex >= 0) {
        const auto &bin = radfieldbin_solutions[(static_cast<ptrdiff_t>(nonemptymgi) * RADFIELDBINCOUNT) + binindex];
        if (bin.W >= 0.) {
          const double J_nu = dbb(nu, bin.T_R, bin.W);
          return J_nu;
        }
      }
      return 0.;
    }
  }

  const float T_R_fullspec = grid::get_TR(nonemptymgi);
  const float W_fullspec = grid::get_W(nonemptymgi);
  const double J_nu_fullspec = dbb(nu, T_R_fullspec, W_fullspec);
  return J_nu_fullspec;
}

// return the integral of nu^3 / (exp(h nu / k T) - 1) from nu_lower to nu_upper
// or if times_nu is true, the integral of nu^4 / (exp(h nu / k T) - 1) from nu_lower to nu_upper
auto planck_integral_analytic(const double T_R, const double nu_lower, const double nu_upper, const bool times_nu)
    -> double {
  double integral = 0.;

  if (times_nu) {
    const double debye_upper = gsl_sf_debye_4(HOVERKB * nu_upper / T_R) * pow(nu_upper, 4);
    const double debye_lower = gsl_sf_debye_4(HOVERKB * nu_lower / T_R) * pow(nu_lower, 4);
    integral = TWOHOVERCLIGHTSQUARED * (debye_upper - debye_lower) * T_R / HOVERKB / 4.;
  } else {
    const double debye_upper = gsl_sf_debye_3(HOVERKB * nu_upper / T_R) * pow(nu_upper, 3);
    const double debye_lower = gsl_sf_debye_3(HOVERKB * nu_lower / T_R) * pow(nu_lower, 3);
    integral = TWOHOVERCLIGHTSQUARED * (debye_upper - debye_lower) * T_R / HOVERKB / 3.;

    if (integral == 0.) {
      // double upperexp = exp(HOVERKB * nu_upper / T_R);
      // double upperint = - pow(nu_upper,4) / 4
      //                   + pow(nu_upper,3) * log(1 - upperexp) / HOVERKB
      //                   + 3 * pow(nu_upper,2) * polylog(2,upperexp) / pow(HOVERKB,2)
      //                   - 6 * nu_upper * polylog(3,upperexp) / pow(HOVERKB,3)
      //                   + 6 * polylog(4,upperexp) / pow(HOVERKB,4);
      // double lowerexp = exp(HOVERKB * nu_lower / T_R);
      // double lowerint = - pow(nu_lower,4) / 4
      //                   + pow(nu_lower,3) * log(1 - lowerexp) / HOVERKB
      //                   + 3 * pow(nu_lower,2) * polylog(2,lowerexp) / pow(HOVERKB,2)
      //                   - 6 * nu_lower * polylog(3,lowerexp) / pow(HOVERKB,3)
      //                   + 6 * polylog(4,lowerexp) / pow(HOVERKB,4);
      // double integral2 = TWOHOVERCLIGHTSQUARED * (upperint - lowerint);

      // printout("planck_integral_analytic is zero. debye_upper %g debye_lower %g. Test alternative %g\n",
      //          debye_upper,debye_lower,integral2);
    }
  }

  return integral;
}

// finds the best fitting W and temperature parameters in each spectral bin using J and nuJ
void fit_parameters(const int modelgridindex, const int timestep) {
  set_params_fullspec(modelgridindex, timestep);

  const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    if (J_normfactor[nonemptymgi] <= 0) {
      printout("radfield: FATAL J_normfactor = %g in cell %d at call to fit_parameters", J_normfactor[nonemptymgi],
               modelgridindex);
      std::abort();
    }

    double J_bin_sum = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      J_bin_sum += get_bin_J(modelgridindex, binindex);
    }

    printout("radfield bins sum to J of %g (%.1f%% of total J).\n", J_bin_sum, 100. * J_bin_sum / J[nonemptymgi]);
    printout("radfield: Finding parameters for %d bins...\n", RADFIELDBINCOUNT);

    double J_bin_max = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      const double J_bin = get_bin_J(modelgridindex, binindex);
      J_bin_max = std::max(J_bin_max, J_bin);
    }

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      const double nu_lower = get_bin_nu_lower(binindex);
      const double nu_upper = get_bin_nu_upper(binindex);
      const double J_bin = get_bin_J(modelgridindex, binindex);
      float T_R_bin = -1.;
      double W_bin = -1.;
      const int contribcount = get_bin_contribcount(modelgridindex, binindex);

      if (contribcount > 0) {
        {
          T_R_bin = find_T_R(modelgridindex, binindex);

          if (binindex == RADFIELDBINCOUNT - 1) {
            const auto T_e = grid::get_Te(nonemptymgi);
            printout("    replacing bin %d T_R %7.1f with cell T_e = %7.1f\n", binindex,
                     get_bin_T_R(modelgridindex, binindex), T_e);
            T_R_bin = T_e;
          }

          double planck_integral_result = planck_integral(T_R_bin, nu_lower, nu_upper, false);
          //          printout("planck_integral(T_R=%g, nu_lower=%g, nu_upper=%g) = %g\n", T_R_bin, nu_lower,
          //          nu_upper, planck_integral_result);

          W_bin = J_bin / planck_integral_result;

          if (W_bin > 1e4) {
            //            printout("T_R_bin %g, nu_lower %g, nu_upper %g\n", T_R_bin, nu_lower, nu_upper);
            printout("W %g too high, trying setting T_R of bin %d to %g. J_bin %g planck_integral %g\n", W_bin,
                     binindex, T_R_max, J_bin, planck_integral_result);
            planck_integral_result = planck_integral(T_R_max, nu_lower, nu_upper, false);
            W_bin = J_bin / planck_integral_result;
            if (W_bin > 1e4) {
              printout("W still very high, W=%g. Zeroing bin...\n", W_bin);
              T_R_bin = -99.;
              W_bin = 0.;
            } else {
              printout("new W is %g. Continuing with this value\n", W_bin);
              T_R_bin = T_R_max;
            }
          }
        }
      } else {
        T_R_bin = 0.;
        W_bin = 0.;
      }

      const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
      radfieldbin_solutions[mgibinindex].T_R = T_R_bin;
      radfieldbin_solutions[mgibinindex].W = W_bin;
    }

    write_to_file(modelgridindex, timestep);
  }
}

void set_J_normfactor(const int nonemptymgi, const double normfactor) { J_normfactor[nonemptymgi] = normfactor; }

void normalise_J(const int modelgridindex, const double estimator_normfactor_over4pi) {
  const int nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  assert_always(std::isfinite(J[nonemptymgi]));
  J[nonemptymgi] *= estimator_normfactor_over4pi;
  for (int i = 0; i < detailed_linecount; i++) {
    prev_Jb_lu_normed[modelgridindex][i].value = Jb_lu_raw[modelgridindex][i].value * estimator_normfactor_over4pi;
    prev_Jb_lu_normed[modelgridindex][i].contribcount = Jb_lu_raw[modelgridindex][i].contribcount;
  }
}

void normalise_bf_estimators(const int nts, const int nts_prev, const int titer, const double deltat) {
  if (globals::rank_in_node != 0) {
    return;
  }
  if (globals::lte_iteration) {
    return;
  }
  if (nts == globals::timestep_initial && titer == 0) {
    return;
  }
  const auto bfestimcount = globals::bfestimcount;
  const ptrdiff_t nonempty_npts_model = grid::get_nonempty_npts_model();
  for (ptrdiff_t nonemptymgi = 0; nonemptymgi < nonempty_npts_model; nonemptymgi++) {
    const auto mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    if (grid::modelgrid[nonemptymgi].thick == 1) {
      continue;
    }
    const double deltaV =
        grid::get_modelcell_assocvolume_tmin(mgi) * pow(globals::timesteps[nts_prev].mid / globals::tmin, 3);
    const double estimator_normfactor = 1 / deltaV / deltat / globals::nprocs;
    for (int i = 0; i < bfestimcount; i++) {
      const auto mgibfindex = (nonemptymgi * bfestimcount) + i;
      prev_bfrate_normed[mgibfindex] = bfrate_raw[mgibfindex] * (estimator_normfactor / H);
    }
  }
}

auto get_bfrate_estimator(const int element, const int lowerion, const int lower, const int phixstargetindex,
                          const int modelgridindex) -> double {
  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    const int allcontindex = get_bfcontindex(element, lowerion, lower, phixstargetindex);
    if (allcontindex >= 0) {
      const auto bfestimindex = globals::allcont[allcontindex].bfestimindex;
      if (bfestimindex >= 0) {
        const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
        return prev_bfrate_normed[(nonemptymgi * globals::bfestimcount) + bfestimindex];
      }
    }
  }

  return -1.;
}

void normalise_nuJ(const int modelgridindex, const double estimator_normfactor_over4pi) {
  const int nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  assert_always(std::isfinite(nuJ[nonemptymgi]));
  nuJ[nonemptymgi] *= estimator_normfactor_over4pi;
}

auto get_T_J_from_J(const int modelgridindex) -> double {
  const int nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  const double T_J = pow(J[nonemptymgi] * PI / STEBO, 1. / 4.);
  if (!std::isfinite(T_J)) {
    // keep old value of T_J
    printout("[warning] get_T_J_from_J: T_J estimator infinite in cell %d, use value of last timestep\n",
             modelgridindex);
    return grid::get_TR(nonemptymgi);
  }
  // Make sure that T is in the allowed temperature range.
  if (T_J > MAXTEMP) {
    printout("[warning] get_T_J_from_J: T_J would be %.1f > MAXTEMP. Clamping to MAXTEMP = %.0f K\n", T_J, MAXTEMP);
    return MAXTEMP;
  }
  if (T_J < MINTEMP) {
    printout("[warning] get_T_J_from_J: T_J would be %.1f < MINTEMP. Clamping to MINTEMP = %.0f K\n", T_J, MINTEMP);
    return MINTEMP;
  }
  return T_J;
}

#ifdef DO_TITER
void titer_J(const int modelgridindex) {
  const int nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  if (J_reduced_save[nonemptymgi] >= 0) {
    J[nonemptymgi] = (J[nonemptymgi] + J_reduced_save[nonemptymgi]) / 2;
  }
  J_reduced_save[nonemptymgi] = J[nonemptymgi];
}

void titer_nuJ(const int modelgridindex) {
  const int nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
  if (nuJ_reduced_save[nonemptymgi] >= 0) {
    nuJ[nonemptymgi] = (nuJ[nonemptymgi] + nuJ_reduced_save[nonemptymgi]) / 2;
  }
  nuJ_reduced_save[nonemptymgi] = nuJ[nonemptymgi];
}
#endif

#ifdef MPI_ON
void reduce_estimators()
// reduce and broadcast (allreduce) the estimators for J and nuJ in all bins
{
  const ptrdiff_t nonempty_npts_model = grid::get_nonempty_npts_model();

  MPI_Allreduce(MPI_IN_PLACE, J.data(), nonempty_npts_model, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, nuJ.data(), nonempty_npts_model, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    for (ptrdiff_t nonemptymgi = 0; nonemptymgi < nonempty_npts_model; nonemptymgi++) {
      MPI_Reduce(globals::rank_in_node == 0 ? MPI_IN_PLACE : &bfrate_raw[nonemptymgi * globals::bfestimcount],
                 &bfrate_raw[nonemptymgi * globals::bfestimcount], globals::bfestimcount, MPI_DOUBLE, MPI_SUM, 0,
                 globals::mpi_comm_node);
    }
    if (globals::rank_in_node == 0) {
      MPI_Allreduce(MPI_IN_PLACE, bfrate_raw.data(), nonempty_npts_model * globals::bfestimcount, MPI_DOUBLE, MPI_SUM,
                    globals::mpi_comm_internode);
    }
    MPI_Bcast(bfrate_raw.data(), nonempty_npts_model * globals::bfestimcount, MPI_DOUBLE, 0, globals::mpi_comm_node);
  }

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    const auto sys_time_start_reduction = std::time(nullptr);
    printout("Reducing binned radiation field estimators");
    assert_always(radfieldbins != nullptr);

    for (ptrdiff_t nonemptymgi = 0; nonemptymgi < nonempty_npts_model; nonemptymgi++) {
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
        const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
        MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[mgibinindex].J_raw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[mgibinindex].nuJ_raw, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &radfieldbins[mgibinindex].contribcount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      }
    }
    const int duration_reduction = std::time(nullptr) - sys_time_start_reduction;
    printout(" (took %d s)\n", duration_reduction);
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    const auto sys_time_start_reduction = std::time(nullptr);
    printout("Reducing detailed line estimators");

    for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
      if (grid::get_numpropcells(modelgridindex) > 0) {
        for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
          MPI_Allreduce(MPI_IN_PLACE, &Jb_lu_raw[modelgridindex][jblueindex].value, 1, MPI_DOUBLE, MPI_SUM,
                        MPI_COMM_WORLD);
          MPI_Allreduce(MPI_IN_PLACE, &Jb_lu_raw[modelgridindex][jblueindex].contribcount, 1, MPI_INT, MPI_SUM,
                        MPI_COMM_WORLD);
        }
      }
    }
    const int duration_reduction = std::time(nullptr) - sys_time_start_reduction;
    printout(" (took %d s)\n", duration_reduction);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// broadcast computed radfield results including parameters
// from the cells belonging to root process to all processes
void do_MPI_Bcast(const ptrdiff_t nonemptymgi, const int root, const int root_node_id) {
  MPI_Bcast(&J_normfactor[nonemptymgi], 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
      if (globals::rank_in_node == 0) {
        MPI_Bcast(&radfieldbin_solutions[mgibinindex].W, 1, MPI_FLOAT, root_node_id, globals::mpi_comm_internode);
        MPI_Bcast(&radfieldbin_solutions[mgibinindex].T_R, 1, MPI_FLOAT, root_node_id, globals::mpi_comm_internode);
      }
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
      MPI_Bcast(&prev_Jb_lu_normed[modelgridindex][jblueindex].value, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&prev_Jb_lu_normed[modelgridindex][jblueindex].contribcount, 1, MPI_INT, root, MPI_COMM_WORLD);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
}
#endif

void write_restart_data(FILE *gridsave_file) {
  printout("binned radiation field and detailed lines, ");

  fprintf(gridsave_file, "%d\n", 30490824);  // special number marking the beginning of radfield data

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    fprintf(gridsave_file, "%d %la %la %la %la\n", RADFIELDBINCOUNT, nu_lower_first_initial, nu_upper_last_initial,
            T_R_min, T_R_max);

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      fprintf(gridsave_file, "%d %la\n", binindex, get_bin_nu_upper(binindex));
    }
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    const int nbfcontinua = globals::nbfcontinua;
    fprintf(gridsave_file, "%d\n", nbfcontinua);

    const int bfestimcount = globals::bfestimcount;
    fprintf(gridsave_file, "%d\n", bfestimcount);

    for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
      if (grid::get_numpropcells(modelgridindex) > 0) {
        const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
        fprintf(gridsave_file, "%d\n", modelgridindex);
        for (int i = 0; i < bfestimcount; i++) {
          fprintf(gridsave_file, "%a ", prev_bfrate_normed[(nonemptymgi * bfestimcount) + i]);
        }
      }
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    fprintf(gridsave_file, "%d\n", detailed_linecount);

    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
      fprintf(gridsave_file, "%d ", detailed_lineindicies[jblueindex]);
    }
  }

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    if (grid::get_numpropcells(modelgridindex) > 0) {
      const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
      assert_testmodeonly(nonemptymgi >= 0);
      fprintf(gridsave_file, "%d %la\n", modelgridindex, J_normfactor[nonemptymgi]);

      if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
        for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
          const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
          fprintf(gridsave_file, "%la %la %a %a %d\n", radfieldbins[mgibinindex].J_raw,
                  radfieldbins[mgibinindex].nuJ_raw, radfieldbin_solutions[mgibinindex].W,
                  radfieldbin_solutions[mgibinindex].T_R, radfieldbins[mgibinindex].contribcount);
        }
      }

      if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
        for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
          fprintf(gridsave_file, "%la %d\n", Jb_lu_raw[modelgridindex][jblueindex].value,
                  Jb_lu_raw[modelgridindex][jblueindex].contribcount);
        }
      }
    }
  }
  fprintf(gridsave_file, "%d\n", 42809403);  // special number marking the end of radfield data
}

void read_restart_data(FILE *gridsave_file) {
  printout("Reading restart data for radiation field\n");

  int code_check = 0;
  assert_always(fscanf(gridsave_file, "%d\n", &code_check) == 1);
  if (code_check != 30490824) {
    printout("ERROR: Beginning of radfield restart data not found! Found %d instead of 30490824\n", code_check);
    std::abort();
  }

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    double T_R_min_in{NAN};
    double T_R_max_in{NAN};
    double nu_lower_first_initial_in{NAN};
    double nu_upper_last_initial_in{NAN};
    int bincount_in = 0;
    assert_always(fscanf(gridsave_file, "%d %la %la %la %la\n", &bincount_in, &nu_lower_first_initial_in,
                         &nu_upper_last_initial_in, &T_R_min_in, &T_R_max_in) == 5);

    double nu_lower_first_ratio = nu_lower_first_initial_in / nu_lower_first_initial;
    if (nu_lower_first_ratio > 1.0) {
      nu_lower_first_ratio = 1 / nu_lower_first_ratio;
    }

    double nu_upper_last_ratio = nu_upper_last_initial_in / nu_upper_last_initial;
    if (nu_upper_last_ratio > 1.0) {
      nu_upper_last_ratio = 1 / nu_upper_last_ratio;
    }

    if (bincount_in != RADFIELDBINCOUNT || T_R_min_in != T_R_min || T_R_max_in != T_R_max ||
        nu_lower_first_ratio < 0.999 || nu_upper_last_ratio < 0.999) {
      printout(
          "ERROR: gridsave file specifies %d bins, nu_lower_first_initial %lg nu_upper_last_initial %lg T_R_min %lg "
          "T_R_max %lg\n",
          bincount_in, nu_lower_first_initial_in, nu_upper_last_initial_in, T_R_min_in, T_R_max_in);
      printout("require %d bins, nu_lower_first_initial %lg nu_upper_last_initial %lg T_R_min %lg T_R_max %lg\n",
               RADFIELDBINCOUNT, nu_lower_first_initial, nu_upper_last_initial, T_R_min, T_R_max);
      std::abort();
    }

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      int binindex_in = 0;
      double nu_upper_in = NAN;
      assert_always(fscanf(gridsave_file, "%d %la\n", &binindex_in, &nu_upper_in) == 2);
      assert_always(binindex_in == binindex);
      assert_always(nu_upper_in == get_bin_nu_upper(binindex));
    }
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    int gridsave_nbf_in = 0;
    assert_always(fscanf(gridsave_file, "%d\n", &gridsave_nbf_in) == 1);
    assert_always(gridsave_nbf_in == globals::nbfcontinua);

    int gridsave_nbfestim_in = 0;
    assert_always(fscanf(gridsave_file, "%d\n", &gridsave_nbfestim_in) == 1);
    assert_always(gridsave_nbfestim_in == globals::bfestimcount);

    for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
      if (grid::get_numpropcells(modelgridindex) > 0) {
        const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
        int mgi_in = 0;
        assert_always(fscanf(gridsave_file, "%d\n", &mgi_in) == 1);
        assert_always(mgi_in == modelgridindex);
        for (int i = 0; i < globals::bfestimcount; i++) {
          float bfrate_normed = 0;
          assert_always(fscanf(gridsave_file, "%a ", &bfrate_normed) == 1);

          if (globals::rank_in_node == 0) {
            prev_bfrate_normed[(nonemptymgi * globals::bfestimcount) + i] = bfrate_normed;
          }
        }
      }
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    int detailed_linecount_in = 0;
    assert_always(fscanf(gridsave_file, "%d\n", &detailed_linecount_in) == 1);

    if (detailed_linecount_in != detailed_linecount) {
      printout("ERROR: gridsave file specifies %d detailed lines but this simulation has %d.\n", detailed_linecount_in,
               detailed_linecount);
      std::abort();
    }

    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
      assert_always(fscanf(gridsave_file, "%d ", &detailed_lineindicies[jblueindex]) == 1);
    }
  }

  for (int modelgridindex = 0; modelgridindex < grid::get_npts_model(); modelgridindex++) {
    if (grid::get_numpropcells(modelgridindex) > 0) {
      const ptrdiff_t nonemptymgi = grid::get_nonemptymgi_of_mgi(modelgridindex);
      int mgi_in = 0;
      assert_always(fscanf(gridsave_file, "%d %la\n", &mgi_in, &J_normfactor[nonemptymgi]) == 2);
      if (mgi_in != modelgridindex) {
        printout("ERROR: expected data for cell %d but found cell %d\n", modelgridindex, mgi_in);
        std::abort();
      }

      if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
        for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
          const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
          float W = 0;
          float T_R = 0;
          assert_always(fscanf(gridsave_file, "%la %la %a %a %d\n", &radfieldbins[mgibinindex].J_raw,
                               &radfieldbins[mgibinindex].nuJ_raw, &W, &T_R,
                               &radfieldbins[mgibinindex].contribcount) == 5);
#ifdef MPI_ON
          if (globals::rank_in_node == 0)
#endif
          {
            radfieldbin_solutions[mgibinindex].W = W;
            radfieldbin_solutions[mgibinindex].T_R = T_R;
          }
        }
      }

      if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
        for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
          assert_always(fscanf(gridsave_file, "%la %d\n", &Jb_lu_raw[modelgridindex][jblueindex].value,
                               &Jb_lu_raw[modelgridindex][jblueindex].contribcount) == 2);
        }
      }
    }
  }
  assert_always(fscanf(gridsave_file, "%d\n", &code_check) == 1);
  if (code_check != 42809403) {
    printout("ERROR: End of radfield restart data not found! Found %d instead of 42809403\n", code_check);
    std::abort();
  }
}

}  // namespace radfield
