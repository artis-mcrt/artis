#include "radfield.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <format>
#include <fstream>
#include <ios>
#include <iterator>
#include <span>
#include <tuple>
#include <vector>

#pragma clang unsafe_buffer_usage begin
#if defined(USEBOOST) && USEBOOST
#include <boost/assert/source_location.hpp>
#include <boost/math/tools/toms748_solve.hpp>
#else
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#endif
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_debye.h>
#include <mpi.h>
#pragma clang unsafe_buffer_usage end

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
  float W;  // dilution (scaling) factor
  float T_R;  // radiation temperature
};

struct RadFieldBins {
  std::vector<double> J_raw;
  std::vector<double> nuJ_raw;
  std::vector<int> contribcount;

  void resize(const ptrdiff_t nonempty_npts_model) {
    resize_exactly(J_raw, nonempty_npts_model * RADFIELDBINCOUNT);
    resize_exactly(nuJ_raw, nonempty_npts_model * RADFIELDBINCOUNT);
    resize_exactly(contribcount, nonempty_npts_model * RADFIELDBINCOUNT);
  }
};

constexpr double radfieldbins_delta_nu =
    (nu_upper_last_initial - nu_lower_first_initial) / (RADFIELDBINCOUNT - 1);  // - 1 for the top super bin

RadFieldBins radfieldbins;

std::span<float> radfieldbin_solutions_W;
MPI_Win win_radfieldbin_solutions_W = MPI_WIN_NULL;

std::span<float> radfieldbin_solutions_T_R;
MPI_Win win_radfieldbin_solutions_T_R = MPI_WIN_NULL;

MPI_Win win_prev_bfrate_normed = MPI_WIN_NULL;

struct Jb_lu_estimator {
  double value = 0.;
  int contribcount = 0;
};

int detailed_linecount = 0;

// array of indices into the linelist[] array for selected lines
std::vector<int> detailed_lineindices;

std::vector<std::vector<Jb_lu_estimator>> prev_Jb_lu_normed{};  // value from the previous timestep
std::vector<std::vector<Jb_lu_estimator>> Jb_lu_raw{};  // unnormalised estimator for the current timestep

std::span<float> prev_bfrate_normed{};  // values from the previous timestep
std::vector<double> bfrate_raw;  // unnormalised estimators for the current timestep

// J and nuJ are accumulated and then normalised in-place
// i.e. be sure the normalisation has been applied (exactly once) before using the values here!

std::vector<double> J;  // after normalisation: [ergs/s/sr/cm2/Hz]
#ifdef DO_TITER
std::vector<double> J_reduced_save;
#endif

std::vector<double> nuJ;  // after normalisation: [ergs/s/sr/cm2]
#ifdef DO_TITER
std::vector<double> nuJ_reduced_save;
#endif

struct GSL_PlanckIntegralParas {
  double T_R;
  bool times_nu;
};

struct GSLTempSolverParams {
  int nonemptymgi;
  int binindex;
};

std::fstream radfieldfile;

constexpr auto get_bin_nu_upper(const int binindex) -> double {
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  if (binindex == RADFIELDBINCOUNT - 1) {
    return nu_upper_superbin;
  }
  return nu_lower_first_initial + ((binindex + 1) * radfieldbins_delta_nu);
}

constexpr auto get_bin_nu_lower(const int binindex) -> double {
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);

  if (binindex > 0) {
    return get_bin_nu_upper(binindex - 1);
  }
  return nu_lower_first_initial;
}

// find the left-closed bin [nu_lower, nu_upper) that nu belongs to
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
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < (RADFIELDBINCOUNT - 1));  // -1 because the superbin is a special case

  return binindex;
}

// associate a Jb_lu estimator with a particular lineindex to be used
// instead of the general radiation field model
void add_detailed_line(const int lineindex) {
  detailed_linecount++;
  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    prev_Jb_lu_normed[nonemptymgi].push_back({.value = 0, .contribcount = 0});
    assert_always(detailed_linecount == std::ssize(prev_Jb_lu_normed[nonemptymgi]));

    // zero_estimators should do the next part anyway, but just to be sure:
    Jb_lu_raw[nonemptymgi].push_back({.value = 0, .contribcount = 0});
    assert_always(detailed_linecount == std::ssize(Jb_lu_raw[nonemptymgi]));
  }
  detailed_lineindices.push_back(lineindex);
  assert_always(detailed_linecount == std::ssize(detailed_lineindices));
}

// get the normalised J value for a bin
auto get_bin_J(const std::ptrdiff_t nonemptymgi, const int binindex) -> double {
  assert_testmodeonly(J_normfactor[nonemptymgi] > 0.0);
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  return radfieldbins.J_raw[(nonemptymgi * RADFIELDBINCOUNT) + binindex] * J_normfactor[nonemptymgi];
}

// get the normalised nuJ value for a bin
auto get_bin_nuJ(const std::ptrdiff_t nonemptymgi, const int binindex) -> double {
  assert_testmodeonly(J_normfactor[nonemptymgi] > 0.0);
  assert_testmodeonly(binindex >= 0);
  assert_testmodeonly(binindex < RADFIELDBINCOUNT);
  return radfieldbins.nuJ_raw[(nonemptymgi * RADFIELDBINCOUNT) + binindex] * J_normfactor[nonemptymgi];
}

// get <nuJ> / <J> for a bin
auto get_bin_nu_bar(const std::ptrdiff_t nonemptymgi, const int binindex) -> double {
  const double nuJ_sum = get_bin_nuJ(nonemptymgi, binindex);
  const double J_sum = get_bin_J(nonemptymgi, binindex);
  return nuJ_sum / J_sum;
}

auto get_bin_contribcount(const std::ptrdiff_t nonemptymgi, const int binindex) -> int {
  return radfieldbins.contribcount[(nonemptymgi * RADFIELDBINCOUNT) + binindex];
}

auto get_bin_W(const std::ptrdiff_t nonemptymgi, const int binindex) -> float {
  return radfieldbin_solutions_W[(nonemptymgi * RADFIELDBINCOUNT) + binindex];
}

auto get_bin_T_R(const std::ptrdiff_t nonemptymgi, const int binindex) -> float {
  return radfieldbin_solutions_T_R[(nonemptymgi * RADFIELDBINCOUNT) + binindex];
}

constexpr auto gsl_integrand_planck(const double nu, void* const voidparas) -> double {
  const auto* paras = static_cast<const GSL_PlanckIntegralParas*>(voidparas);
  const auto T_R = paras->T_R;

  double integrand = TWOHOVERCLIGHTSQUARED * std::pow(nu, 3) / (std::expm1(HOVERKB * nu / T_R));

  if (paras->times_nu) {
    integrand *= nu;
  }

  return integrand;
}

void update_bfestimators(const ptrdiff_t nonemptymgi, const double distance_e_cmf, const double nu_cmf,
                         const double doppler_nucmf_on_nurf, const Phixslist& phixslist) {
  assert_testmodeonly(DETAILED_BF_ESTIMATORS_ON);

  const double distance_e_cmf_over_nu =
      distance_e_cmf / nu_cmf * doppler_nucmf_on_nurf;  // TODO: Luke: why did I put a doppler factor here?

  // I think the nu_cmf slightly differs from when the phixslist was calculated
  // so the nu condition on this nu_cmf can truncate the list further compared to what was used in the calculation
  // of phixslist.gamma_contr
  const auto bfestimcount = std::ssize(globals::bfestim_nu_edge);

  assert_testmodeonly(phixslist.bfestimend <= bfestimcount);
  const auto bfestimend =
      std::distance(globals::bfestim_nu_edge.begin(),
                    std::ranges::upper_bound(globals::bfestim_nu_edge.first(phixslist.bfestimend), nu_cmf));
  assert_testmodeonly(bfestimend <= bfestimcount);
  assert_testmodeonly(phixslist.bfestimbegin >= 0);
  const auto bfestimbegin =
      std::distance(globals::bfestim_nu_edge.begin(),
                    std::ranges::lower_bound(
                        globals::bfestim_nu_edge.subspan(phixslist.bfestimbegin, bfestimend - phixslist.bfestimbegin),
                        nu_cmf, [](const double nu_edge, const double find_nu_cmf) {
                          return nu_edge * last_phixs_nuovernuedge < find_nu_cmf;
                        }));

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

  const GSL_PlanckIntegralParas intparas = {.T_R = T_R, .times_nu = times_nu};

  const int status = integrator<gsl_integrand_planck>(intparas, nu_lower, nu_upper, epsabs, epsrel, GSL_INTEG_GAUSS61,
                                                      &integral, &error);
  if (status != 0) {
    printlnlog("planck_integral integrator status {}, GSL_FAILURE= {}. Integral value {:g}, setting to zero.", status,
               static_cast<int>(GSL_FAILURE), integral);
    integral = 0.;
  }

  return integral;
}

// difference between the average nu and the average nu of a Planck function
// at temperature T_R, in the frequency range corresponding to a bin
auto delta_nu_bar(const double T_R, const int nonemptymgi, const int binindex) -> double {
  const double nu_lower = get_bin_nu_lower(binindex);
  const double nu_upper = get_bin_nu_upper(binindex);
  const double nu_bar_estimator = get_bin_nu_bar(nonemptymgi, binindex);

  const double nu_times_planck_numerical = planck_integral(T_R, nu_lower, nu_upper, true);
  const double planck_integral_numerical = planck_integral(T_R, nu_lower, nu_upper, false);
  const double nu_bar_planck_T_R = nu_times_planck_numerical / planck_integral_numerical;

  // double nu_times_planck_integral = planck_integral_analytic(T_R, nu_lower, nu_upper, true);
  // double planck_integral_result = planck_integral_analytic(T_R, nu_lower, nu_upper, false);
  // double nu_bar_planck = nu_times_planck_integral / planck_integral_result;

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
    printlnlog(
        "delta_nu_bar is {:g}. nu_bar_planck_T_R {:g} nu_times_planck_numerical {:g} planck_integral_numerical {:g} "
        "nu_bar_estimator {:g}",
        delta_nu_bar, nu_bar_planck_T_R, nu_times_planck_numerical, planck_integral_numerical, nu_bar_estimator);
  }

  return delta_nu_bar;
}

auto delta_nu_bar(const double T_R, void* const voidparas) -> double {  // cppcheck-suppress constParameterPointer
  const auto* const params = static_cast<const GSLTempSolverParams*>(voidparas);
  return delta_nu_bar(T_R, params->nonemptymgi, params->binindex);
}

auto find_T_R(const int nonemptymgi, const int binindex) -> float {
  const auto f_deltanubar = [nonemptymgi, binindex](const double T_R) {
    return delta_nu_bar(T_R, nonemptymgi, binindex);
  };

  // Check whether the equation has a root in [T_min,T_max]
  double delta_nu_bar_min = f_deltanubar(T_R_min);
  double delta_nu_bar_max = f_deltanubar(T_R_max);

  if (!std::isfinite(delta_nu_bar_min) || !std::isfinite(delta_nu_bar_max)) {
    delta_nu_bar_max = delta_nu_bar_min = -1;
  }

  if (delta_nu_bar_min * delta_nu_bar_max < 0) {
    // If there is a root in the interval, solve for T_R

    constexpr double epsrel = 1e-4;
    const auto maxit = 100U;
#if defined(USEBOOST) && USEBOOST
    // use TOMS 748 solver from Boost
    boost::uintmax_t iteration_num = maxit;
    auto result = boost::math::tools::toms748_solve(f_deltanubar, T_R_min, T_R_max, ftol<epsrel>, iteration_num);
    const auto T_R_solution = static_cast<float>(0.5 * (result.first + result.second));
    if (iteration_num >= maxit) {
      printlnlog("[warning] find_T_R: T_R did not converge within {} iterations.", iteration_num);
    }
    return T_R_solution;
#else
    GSLTempSolverParams paras{.nonemptymgi = nonemptymgi, .binindex = binindex};
    gsl_function find_T_R_f = {.function = &delta_nu_bar, .params = &paras};

    // one dimensional gsl root solver, bracketing type
    gsl_root_fsolver* T_R_solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(T_R_solver, &find_T_R_f, T_R_min, T_R_max);
    int status = 0;
    float T_R_solution = 0.;
    for (auto iteration_num = 0U; iteration_num <= maxit; iteration_num++) {
      gsl_root_fsolver_iterate(T_R_solver);
      T_R_solution = static_cast<float>(gsl_root_fsolver_root(T_R_solver));

      const double T_R_lower = gsl_root_fsolver_x_lower(T_R_solver);
      const double T_R_upper = gsl_root_fsolver_x_upper(T_R_solver);
      status = gsl_root_test_interval(T_R_lower, T_R_upper, 0., epsrel);

      if (status != GSL_CONTINUE) {
        break;
      }
    }

    if (status == GSL_CONTINUE) {
      printlnlog("[warning] find_T_R: T_R did not converge within {} iterations", maxit);
    }

    gsl_root_fsolver_free(T_R_solver);
    return T_R_solution;
#endif
  } else if (delta_nu_bar_max < 0) {
    // Thermal balance equation always negative ===> T_R = T_min
    // Calculate the rates again at this T_e to print them to file
    printlnlog("find_T_R: cell {} bin {:4} no solution in interval, clamping to T_R_max={:g}",
               grid::get_mgi_of_nonemptymgi(nonemptymgi), binindex, T_R_max);
    return T_R_max;
  } else {
    printlnlog("find_T_R: cell {} bin {:4} no solution in interval, clamping to T_R_min={:g}",
               grid::get_mgi_of_nonemptymgi(nonemptymgi), binindex, T_R_min);
    return T_R_min;
  }

}  // namespace radfield

void set_params_fullspec(const int nonemptymgi, const int timestep) {
  const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
  const double nubar = nuJ[nonemptymgi] / J[nonemptymgi];
  if (!std::isfinite(nubar) || nubar == 0.) {
    printlnlog("[warning] T_R estimator infinite in cell {}, keep T_R, T_J, W of last timestep. J = {:g}. nuJ = {:g}",
               modelgridindex, J[nonemptymgi], nuJ[nonemptymgi]);
  } else {
    auto T_J = static_cast<float>(pow(J[nonemptymgi] * PI / STEBO, 1 / 4.));
    if (T_J > MAXTEMP) {
      printlnlog("[warning] temperature estimator T_J = {:g} exceeds T_max {:g} in cell {}. Setting T_J = T_max!", T_J,
                 MAXTEMP, modelgridindex);
      T_J = MAXTEMP;
    } else if (T_J < MINTEMP) {
      printlnlog("[warning] temperature estimator T_J = {:g} below T_min {:g} in cell {}. Setting T_J = T_min!", T_J,
                 MINTEMP, modelgridindex);
      T_J = MINTEMP;
    }
    grid::set_TJ(nonemptymgi, T_J);

    auto T_R = static_cast<float>(H * nubar / KB / 3.832229494);
    if (T_R > MAXTEMP) {
      printlnlog("[warning] temperature estimator T_R = {:g} exceeds T_max {:g} in cell {}. Setting T_R = T_max!", T_R,
                 MAXTEMP, modelgridindex);
      T_R = MAXTEMP;
    } else if (T_R < MINTEMP) {
      printlnlog("[warning] temperature estimator T_R = {:g} below T_min {:g} in cell {}. Setting T_R = T_min!", T_R,
                 MINTEMP, modelgridindex);
      T_R = MINTEMP;
    }
    grid::set_TR(nonemptymgi, T_R);

    const auto W = static_cast<float>(J[nonemptymgi] * PI / STEBO / pow(T_R, 4));
    grid::set_W(nonemptymgi, W);

    printlnlog(
        "Full-spectrum fit radfield for cell {} at timestep {}: J {:g}, nubar {:5.1f} Angstrom, T_J {:g}, T_R {:g}, W "
        "{:g}",
        modelgridindex, timestep, J[nonemptymgi], 1e8 * CLIGHT / nubar, T_J, T_R, W);
  }
}

auto get_bfcontindex(const int element, const int lowerion, const int lower, const int phixstargetindex) -> int {
  // simple linear search seems to be faster than the binary search
  // possibly because lower frequency transitions near start of list are more likely to be called?
  int bfcontindex = 0;
  for (; bfcontindex < globals::nbfcontinua; bfcontindex++) {
    if ((globals::allcont_element[bfcontindex] == element) && (globals::allcont_ion[bfcontindex] == lowerion) &&
        (globals::allcont_level[bfcontindex] == lower) &&
        (globals::allcont_phixstargetindex[bfcontindex] == phixstargetindex)) {
      break;
    }
  }
  if (bfcontindex < globals::nbfcontinua) {
    return bfcontindex;
  }
  // not found in the continua list
  return -1;
}

void write_to_file(const int nonemptymgi, const int timestep) {
  assert_always(MULTIBIN_RADFIELD_MODEL_ON);
  const int modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
#ifdef _OPENMP
#pragma omp critical(out_file)
  {
#endif

    int totalcontribs = 0;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      totalcontribs += get_bin_contribcount(nonemptymgi, binindex);
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
        nuJ_out = get_bin_nuJ(nonemptymgi, binindex);
        J_out = get_bin_J(nonemptymgi, binindex);
        T_R = get_bin_T_R(nonemptymgi, binindex);
        W = get_bin_W(nonemptymgi, binindex);
        J_nu_bar = J_out / (nu_upper - nu_lower);
        contribcount = get_bin_contribcount(nonemptymgi, binindex);
      } else if (binindex == -1) {  // bin -1 is the full spectrum fit
        nuJ_out = nuJ[nonemptymgi];
        J_out = J[nonemptymgi];
        T_R = grid::get_TR(nonemptymgi);
        W = grid::get_W(nonemptymgi);
        contribcount = totalcontribs;
      } else {  // use binindex < -1 for detailed line Jb_lu estimators
        const int jblueindex = -2 - binindex;  // -2 is the first detailed line, -3 is the second, etc
        const int lineindex = detailed_lineindices[jblueindex];
        const double nu_trans = globals::linelist.nu[lineindex];
        nu_lower = nu_trans;
        nu_upper = nu_trans;
        nuJ_out = -1.;
        J_out = -1.;
        T_R = -1.;
        W = -1.;
        J_nu_bar = prev_Jb_lu_normed[nonemptymgi][jblueindex].value,
        contribcount = prev_Jb_lu_normed[nonemptymgi][jblueindex].contribcount;
      }

      if (!skipoutput) {
        radfieldfile << std::format("{:d} {:d} {:d} {:.5e} {:.5e} {:.3e} {:.3e} {:.3e} {:d} {:.1f} {:.5e}\n", timestep,
                                    modelgridindex, binindex, nu_lower, nu_upper, nuJ_out, J_out, J_nu_bar,
                                    contribcount, T_R, W);
      }
    }
    radfieldfile.flush();
#ifdef _OPENMP
  }
#endif
}

}  // anonymous namespace

void init(const int my_rank, const int ndo_nonempty) {
  // this should be called only after the atomic data is in memory

  const ptrdiff_t nonempty_npts_model = grid::get_nonempty_npts_model();

  resize_exactly(J_normfactor, nonempty_npts_model + 1);
  resize_exactly(J, nonempty_npts_model + 1);

#ifdef DO_TITER
  resize_exactly(J_reduced_save, nonempty_npts_model + 1);
#endif

  // J and nuJ are accumulated and then normalised in-place
  // i.e. be sure the normalisation has been applied (exactly once) before using the values here!
  resize_exactly(nuJ, nonempty_npts_model + 1);
#ifdef DO_TITER
  resize_exactly(nuJ, nonempty_npts_model + 1);
#endif

  resize_exactly(prev_Jb_lu_normed, nonempty_npts_model);
  resize_exactly(Jb_lu_raw, nonempty_npts_model);

  detailed_linecount = 0;

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    for (int i = 0; i < globals::nlines; i++) {
      const int element = globals::linelist.elementindex[i];
      const int Z = get_atomicnumber(element);
      if (Z == 26) {
        const int lowerlevel = globals::linelist.lowerlevelindex[i];
        // const int upperlevel = linelist[i].upperlevelindex;
        // const int ion = linelist[i].ionindex;
        // const int ionstage = get_ionstage(element, ion);
        const double A_ul = globals::linelist.einstein_A[i];

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
          add_detailed_line(i);
        }
      }
    }
    printlnlog("There are {} lines with detailed Jblue_lu estimators.", detailed_linecount);
  }

  printlog("DETAILED_BF_ESTIMATORS {}", DETAILED_BF_ESTIMATORS_ON ? "ON" : "OFF");
  if (DETAILED_BF_ESTIMATORS_ON) {
    printlnlog(" from timestep {}", DETAILED_BF_ESTIMATORS_USEFROMTIMESTEP);
  } else {
    printlnlog("");
  }

  if (MULTIBIN_RADFIELD_MODEL_ON) {
    printlnlog("The multibin radiation field is being used from timestep {} onwards.", FIRST_NLTE_RADFIELD_TIMESTEP);

    printlnlog(
        "Initialising multibin radiation field with {} bins from ({:.2f} eV, {:6.1f} A) to ({:.2f} eV, {:6.1f} A)",
        RADFIELDBINCOUNT, H * nu_lower_first_initial / EV, 1e8 * CLIGHT / nu_lower_first_initial,
        H * nu_upper_last_initial / EV, 1e8 * CLIGHT / nu_upper_last_initial);
    if (ndo_nonempty > 0) {
      assert_always(!radfieldfile.is_open());
      radfieldfile = fstream_required(std::format("radfield_{:04d}.out", my_rank), std::ios::out | std::ios::trunc);
      radfieldfile << "timestep modelgridindex bin_num nu_lower nu_upper nuJ J J_nu_avg ncontrib T_R W\n";
      radfieldfile.flush();
    }

    const size_t mem_usage_bins = nonempty_npts_model * RADFIELDBINCOUNT * ((2 * sizeof(double)) + sizeof(int));
    radfieldbins.resize(nonempty_npts_model);

    printlnlog("[info] mem_usage: radiation field bin accumulators for non-empty cells occupy {:.3f} MB",
               mem_usage_bins / 1024. / 1024.);

    std::tie(radfieldbin_solutions_W, win_radfieldbin_solutions_W) =
        MPI_shared_malloc_span_keepwin<float>(nonempty_npts_model * RADFIELDBINCOUNT);

    std::tie(radfieldbin_solutions_T_R, win_radfieldbin_solutions_T_R) =
        MPI_shared_malloc_span_keepwin<float>(nonempty_npts_model * RADFIELDBINCOUNT);

    const size_t mem_usage_bin_solutions = nonempty_npts_model * RADFIELDBINCOUNT * sizeof(RadFieldBinSolution);
    printlnlog(
        "[info] mem_usage: radiation field bin solutions for non-empty cells occupy {:.3f} MB (node shared memory)",
        mem_usage_bin_solutions / 1024. / 1024.);
  } else {
    printlnlog("The radiation field model is a full-spectrum fit to a single dilute blackbody TR & W.");
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    const auto bfestimcount = std::ssize(globals::bfestim_nu_edge);
    std::tie(prev_bfrate_normed, win_prev_bfrate_normed) =
        MPI_shared_malloc_span_keepwin<float>(nonempty_npts_model * bfestimcount);
    if (globals::rank_in_node == 0) {
      std::ranges::fill(prev_bfrate_normed, 0.);
    }
    MPI_Barrier(globals::mpi_comm_node);
    printlnlog("[info] mem_usage: detailed bf estimators for non-empty cells occupy {:.3f} MB (node shared memory)",
               nonempty_npts_model * bfestimcount * sizeof(float) / 1024. / 1024.);

    resize_exactly(bfrate_raw, nonempty_npts_model * bfestimcount);

    printlnlog("[info] mem_usage: detailed bf estimator acculumators for non-empty cells occupy {:.3f} MB",
               nonempty_npts_model * bfestimcount * sizeof(double) / 1024. / 1024.);
  }

  zero_estimators();

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    MPI_Barrier(globals::mpi_comm_node);
    if (globals::rank_in_node == 0) {
      std::ranges::fill(radfieldbin_solutions_W, -1.);
      std::ranges::fill(radfieldbin_solutions_T_R, -1.);
    }
    MPI_Barrier(globals::mpi_comm_node);
  }
}

// Initialise estimator arrays which hold the last time steps values (used to damp out
// fluctuations over timestep iterations if DO_TITER is defined) to -1.
void initialise_prev_titer_photoionestimators() {
#ifdef DO_TITER
  std::ranges::fill(globals::ffheatingestimator_save, -1.);
  std::ranges::fill(globals::colheatingestimator_save, -1.);
  std::ranges::fill(J_reduced_save, -1.);
  std::ranges::fill(nuJ_reduced_save, -1.);
  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
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
  int high = detailed_linecount - 1;
  while (low <= high) {
    const int mid = low + ((high - low) / 2);
    if (detailed_lineindices[mid] < lineindex) {
      low = mid + 1;
    } else if (detailed_lineindices[mid] > lineindex) {
      high = mid - 1;
    } else {
      assert_always(mid < detailed_linecount);
      return mid;
    }
  }

  return -1;
}

auto get_Jb_lu(const int nonemptymgi, const int jblueindex) -> double {
  assert_always(jblueindex >= 0);
  assert_always(jblueindex < detailed_linecount);
  return prev_Jb_lu_normed[nonemptymgi][jblueindex].value;
}

auto get_Jb_lu_contribcount(const int nonemptymgi, const int jblueindex) -> int {
  assert_always(jblueindex >= 0);
  assert_always(jblueindex < detailed_linecount);
  return prev_Jb_lu_normed[nonemptymgi][jblueindex].contribcount;
}

void close_file() {
  if (radfieldfile.is_open()) {
    radfieldfile.close();
  }

  if (MULTIBIN_RADFIELD_MODEL_ON) {
    radfieldbins = {};
    if (win_radfieldbin_solutions_W != MPI_WIN_NULL) {
      MPI_Win_free(&win_radfieldbin_solutions_W);
      radfieldbin_solutions_W = {};
    }
    if (win_radfieldbin_solutions_T_R != MPI_WIN_NULL) {
      MPI_Win_free(&win_radfieldbin_solutions_T_R);
      radfieldbin_solutions_T_R = {};
    }
  }

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    if (win_prev_bfrate_normed != MPI_WIN_NULL) {
      MPI_Win_free(&win_prev_bfrate_normed);
      prev_bfrate_normed = {};
    }
  }
}

// set up the new bins and clear the estimators in preparation for a timestep
void zero_estimators() {
  std::ranges::fill(J_normfactor, -1.0);
  std::ranges::fill(J, 0.0);
  std::ranges::fill(nuJ, 0.0);
  std::ranges::fill(bfrate_raw, 0.0);

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    std::ranges::fill(radfieldbins.J_raw, 0.0);
    std::ranges::fill(radfieldbins.nuJ_raw, 0.0);
    std::ranges::fill(radfieldbins.contribcount, 0);
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
      std::fill_n(Jb_lu_raw[nonemptymgi].data(), detailed_linecount, Jb_lu_estimator{.value = 0., .contribcount = 0});
    }
  }
}

__host__ __device__ void update_estimators(const ptrdiff_t nonemptymgi, const double distance_e_cmf,
                                           const double nu_cmf, const double doppler_nucmf_on_nurf,
                                           const Phixslist& phixslist, const bool thickcell) {
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
      const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
      atomicadd(radfieldbins.J_raw[mgibinindex], distance_e_cmf);
      atomicadd(radfieldbins.nuJ_raw[mgibinindex], distance_e_cmf * nu_cmf);
      atomicadd(radfieldbins.contribcount[mgibinindex], 1);
    }
  }
}

__host__ __device__ void update_lineestimator(const int nonemptymgi, const int lineindex, const double increment) {
  if constexpr (!DETAILED_LINE_ESTIMATORS_ON) {
    return;
  }

  const int jblueindex = get_Jblueindex(lineindex);
  if (jblueindex >= 0) {
    Jb_lu_raw[nonemptymgi][jblueindex].value += increment;
    Jb_lu_raw[nonemptymgi][jblueindex].contribcount += 1;
  }
}

// mean intensity J_nu [ergs/s/sr/cm2/Hz]
__host__ __device__ auto radfield(const double nu, const int nonemptymgi) -> double {
  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    if (globals::timestep >= FIRST_NLTE_RADFIELD_TIMESTEP) {
      const int binindex = select_bin(nu);
      if (binindex >= 0) {
        const auto W = get_bin_W(nonemptymgi, binindex);
        if (W >= 0.) {
          return dbb(nu, get_bin_T_R(nonemptymgi, binindex), W);
        }
      }
      return 0.;
    }
  }
  // full spectrum fit to a single dilute blackbody
  return dbb(nu, grid::get_TR(nonemptymgi), grid::get_W(nonemptymgi));
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

      // printlnlog("planck_integral_analytic is zero. debye_upper {:g} debye_lower {:g}. Test alternative {:g}",
      //            debye_upper, debye_lower, integral2);
    }
  }

  return integral;
}

// finds the best fitting W and temperature parameters in each spectral bin using J and nuJ
void fit_parameters(const int nonemptymgi, const int timestep) {
  set_params_fullspec(nonemptymgi, timestep);
  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    assert_always(J_normfactor[nonemptymgi] >= 0.);

    double J_bin_sum = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      J_bin_sum += get_bin_J(nonemptymgi, binindex);
    }

    printlnlog("radfield bins sum to J of {:g} ({:.1f}% of total J).", J_bin_sum, 100. * J_bin_sum / J[nonemptymgi]);
    printlnlog("radfield: Finding parameters for {} bins...", RADFIELDBINCOUNT);

    double J_bin_max = 0.;
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      const double J_bin = get_bin_J(nonemptymgi, binindex);
      J_bin_max = std::max(J_bin_max, J_bin);
    }

    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
      const double nu_lower = get_bin_nu_lower(binindex);
      const double nu_upper = get_bin_nu_upper(binindex);
      const double J_bin = get_bin_J(nonemptymgi, binindex);
      float T_R_bin = -1.;
      float W_bin = -1.;
      const int contribcount = get_bin_contribcount(nonemptymgi, binindex);

      if (contribcount > 0) {
        {
          T_R_bin = find_T_R(nonemptymgi, binindex);

          if (binindex == RADFIELDBINCOUNT - 1) {
            const auto T_e = grid::get_Te(nonemptymgi);
            printlnlog("    replacing bin {} T_R {:7.1f} with cell T_e = {:7.1f}", binindex,
                       get_bin_T_R(nonemptymgi, binindex), T_e);
            T_R_bin = T_e;
          }

          double planck_integral_result = planck_integral(T_R_bin, nu_lower, nu_upper, false);
          //          printout("planck_integral(T_R=%g, nu_lower=%g, nu_upper=%g) = %g\n", T_R_bin, nu_lower,
          //          nu_upper, planck_integral_result);

          W_bin = static_cast<float>(J_bin / planck_integral_result);

          if (W_bin > 1e4) {
            //            printout("T_R_bin %g, nu_lower %g, nu_upper %g\n", T_R_bin, nu_lower, nu_upper);
            printlnlog("W {:g} too high, trying setting T_R of bin {} to {:g}. J_bin {:g} planck_integral {:g}", W_bin,
                       binindex, T_R_max, J_bin, planck_integral_result);
            planck_integral_result = planck_integral(T_R_max, nu_lower, nu_upper, false);
            W_bin = static_cast<float>(J_bin / planck_integral_result);
            if (W_bin > 1e4) {
              printlnlog("W still very high, W={:g}. Zeroing bin...", W_bin);
              T_R_bin = -99.;
              W_bin = 0.;
            } else {
              printlnlog("new W is {:g}. Continuing with this value", W_bin);
              T_R_bin = T_R_max;
            }
          }
        }
      } else {
        T_R_bin = 0.;
        W_bin = 0.;
      }

      const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
      radfieldbin_solutions_T_R[mgibinindex] = T_R_bin;
      radfieldbin_solutions_W[mgibinindex] = W_bin;
    }

    write_to_file(nonemptymgi, timestep);
  }
}

void set_J_normfactor(const int nonemptymgi, const double normfactor) { J_normfactor[nonemptymgi] = normfactor; }

void normalise_J(const int nonemptymgi, const double estimator_normfactor_over4pi) {
  assert_always(std::isfinite(J[nonemptymgi]));
  J[nonemptymgi] *= estimator_normfactor_over4pi;
  for (int i = 0; i < detailed_linecount; i++) {
    prev_Jb_lu_normed[nonemptymgi][i].value = Jb_lu_raw[nonemptymgi][i].value * estimator_normfactor_over4pi;
    prev_Jb_lu_normed[nonemptymgi][i].contribcount = Jb_lu_raw[nonemptymgi][i].contribcount;
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
  const auto bfestimcount = std::ssize(globals::bfestim_nu_edge);
  const ptrdiff_t nonempty_npts_model = grid::get_nonempty_npts_model();
  for (auto nonemptymgi = 0Z; nonemptymgi < nonempty_npts_model; nonemptymgi++) {
    if (grid::thick_allcells[nonemptymgi] == 1) {
      continue;
    }
    const auto mgi = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    const double deltaV =
        grid::get_modelcell_assocvolume_tmin(mgi) * pow(globals::timesteps[nts_prev].mid / globals::tmin, 3);
    const double estimator_normfactor = 1 / deltaV / deltat / globals::nprocs;
    for (int i = 0; i < bfestimcount; i++) {
      const auto mgibfindex = (nonemptymgi * bfestimcount) + i;
      prev_bfrate_normed[mgibfindex] = static_cast<float>(bfrate_raw[mgibfindex] * (estimator_normfactor / H));
    }
  }
}

__host__ __device__ auto get_bfrate_estimator(const int element, const int lowerion, const int lower,
                                              const int phixstargetindex, const int nonemptymgi) -> double {
  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    const int allcontindex = get_bfcontindex(element, lowerion, lower, phixstargetindex);
    if (allcontindex >= 0) {
      const auto bfestimindex = globals::allcont_bfestimindex[allcontindex];
      if (bfestimindex >= 0) {
        return prev_bfrate_normed[(nonemptymgi * std::ssize(globals::bfestim_nu_edge)) + bfestimindex];
      }
    }
  }

  return -1.;
}

void normalise_nuJ(const int nonemptymgi, const double estimator_normfactor_over4pi) {
  assert_always(std::isfinite(nuJ[nonemptymgi]));
  nuJ[nonemptymgi] *= estimator_normfactor_over4pi;
}

auto get_T_J_from_J(const int nonemptymgi) -> float {
  const auto T_J = static_cast<float>(pow(J[nonemptymgi] * PI / STEBO, 1. / 4.));
  if (!std::isfinite(T_J)) {
    // keep old value of T_J
    const auto modelgridindex = grid::get_mgi_of_nonemptymgi(nonemptymgi);
    printlnlog("[warning] get_T_J_from_J: T_J estimator infinite in cell {}, use value of last timestep",
               modelgridindex);
    return grid::get_TR(nonemptymgi);
  }
  // Make sure that T is in the allowed temperature range.
  if (T_J > MAXTEMP) {
    printlnlog("[warning] get_T_J_from_J: T_J would be {:.1f} > MAXTEMP. Clamping to MAXTEMP = {:.0f} K", T_J, MAXTEMP);
    return MAXTEMP;
  }
  if (T_J < MINTEMP) {
    printlnlog("[warning] get_T_J_from_J: T_J would be {:.1f} < MINTEMP. Clamping to MINTEMP = {:.0f} K", T_J, MINTEMP);
    return MINTEMP;
  }
  return T_J;
}

#ifdef DO_TITER
void titer_J(const int nonemptymgi) {
  if (J_reduced_save[nonemptymgi] >= 0) {
    J[nonemptymgi] = (J[nonemptymgi] + J_reduced_save[nonemptymgi]) / 2;
  }
  J_reduced_save[nonemptymgi] = J[nonemptymgi];
}

void titer_nuJ(const int nonemptymgi) {
  if (nuJ_reduced_save[nonemptymgi] >= 0) {
    nuJ[nonemptymgi] = (nuJ[nonemptymgi] + nuJ_reduced_save[nonemptymgi]) / 2;
  }
  nuJ_reduced_save[nonemptymgi] = nuJ[nonemptymgi];
}
#endif

void reduce_estimators()
// reduce and broadcast (allreduce) the estimators for J and nuJ in all bins
{
  MPI_Allreduce_safe(J, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce_safe(nuJ, MPI_SUM, MPI_COMM_WORLD);

  if constexpr (DETAILED_BF_ESTIMATORS_ON) {
    // reduce all ranks on each node first, then reduce the node leaders
    // then broadcast the final result to all ranks on each node
    // this seems necessary to avoid congestion compared to a single MPI_Allreduce
    // when using many ranks per node
    MPI_Reduce_safe(bfrate_raw, MPI_SUM, 0, globals::mpi_comm_node);
    if (globals::rank_in_node == 0) {
      MPI_Allreduce_safe(bfrate_raw, MPI_SUM, globals::mpi_comm_internode);
    }
    MPI_Bcast_safe(bfrate_raw, 0, globals::mpi_comm_node);
  }

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    const auto sys_time_start_reduction = std::time(nullptr);
    printlog("Reducing binned radiation field estimators");

    MPI_Allreduce_safe(radfieldbins.J_raw, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce_safe(radfieldbins.nuJ_raw, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce_safe(radfieldbins.contribcount, MPI_SUM, MPI_COMM_WORLD);

    const auto duration_reduction = std::time(nullptr) - sys_time_start_reduction;
    printlnlog(" (took {} s)", duration_reduction);
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    const auto sys_time_start_reduction = std::time(nullptr);
    printlog("Reducing detailed line estimators");

    for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
      for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
        MPI_Allreduce_safe(Jb_lu_raw[nonemptymgi][jblueindex].value, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce_safe(Jb_lu_raw[nonemptymgi][jblueindex].contribcount, MPI_SUM, MPI_COMM_WORLD);
      }
    }
    const auto duration_reduction = std::time(nullptr) - sys_time_start_reduction;
    printlnlog(" (took {} s)", duration_reduction);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// broadcast computed radfield results including parameters
// from the cells belonging to root process to all processes
void do_MPI_Bcast(const ptrdiff_t nonemptymgi, const int root, const int root_node_id) {
  MPI_Bcast_safe(J_normfactor[nonemptymgi], root, MPI_COMM_WORLD);

  if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
    if (globals::rank_in_node == 0) {
      MPI_Bcast_safe(radfieldbin_solutions_W.subspan(nonemptymgi * RADFIELDBINCOUNT, RADFIELDBINCOUNT), root_node_id,
                     globals::mpi_comm_internode);
      MPI_Bcast_safe(radfieldbin_solutions_T_R.subspan(nonemptymgi * RADFIELDBINCOUNT, RADFIELDBINCOUNT), root_node_id,
                     globals::mpi_comm_internode);
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
      MPI_Bcast_safe(prev_Jb_lu_normed[nonemptymgi][jblueindex].value, root, MPI_COMM_WORLD);
      MPI_Bcast_safe(prev_Jb_lu_normed[nonemptymgi][jblueindex].contribcount, root, MPI_COMM_WORLD);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

void write_restart_data(FILE* gridsave_file) {
  printlog("binned radiation field and detailed lines, ");

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

    const int bfestimcount = static_cast<int>(globals::bfestim_nu_edge.size());
    fprintf(gridsave_file, "%d\n", bfestimcount);

    for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
      fprintf(gridsave_file, "%d\n", nonemptymgi);
      for (int i = 0; i < bfestimcount; i++) {
        fprintf(gridsave_file, "%a ", prev_bfrate_normed[(nonemptymgi * bfestimcount) + i]);
      }
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    fprintf(gridsave_file, "%d\n", detailed_linecount);

    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
      fprintf(gridsave_file, "%d ", detailed_lineindices[jblueindex]);
    }
  }

  for (int nonemptymgi = 0; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    assert_testmodeonly(nonemptymgi >= 0);
    fprintf(gridsave_file, "%d %la\n", nonemptymgi, J_normfactor[nonemptymgi]);

    if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
        const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
        fprintf(gridsave_file, "%la %la %a %a %d\n", radfieldbins.J_raw[mgibinindex], radfieldbins.nuJ_raw[mgibinindex],
                radfieldbin_solutions_W[mgibinindex], radfieldbin_solutions_T_R[mgibinindex],
                radfieldbins.contribcount[mgibinindex]);
      }
    }

    if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
      for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
        fprintf(gridsave_file, "%la %d\n", Jb_lu_raw[nonemptymgi][jblueindex].value,
                Jb_lu_raw[nonemptymgi][jblueindex].contribcount);
      }
    }
  }
  fprintf(gridsave_file, "%d\n", 42809403);  // special number marking the end of radfield data
}

void read_restart_data(FILE* gridsave_file) {
  printlnlog("Reading restart data for radiation field");

  int code_check = 0;
  assert_always(fscanf(gridsave_file, "%d\n", &code_check) == 1);
  assert_always(code_check == 30490824);

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
      printlnlog(
          "ERROR: gridsave file specifies {} bins, nu_lower_first_initial {:g} nu_upper_last_initial {:g} T_R_min {:g} "
          "T_R_max {:g}",
          bincount_in, nu_lower_first_initial_in, nu_upper_last_initial_in, T_R_min_in, T_R_max_in);
      printlnlog("require {} bins, nu_lower_first_initial {:g} nu_upper_last_initial {:g} T_R_min {:g} T_R_max {:g}",
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

    const auto bfestimcount = std::ssize(globals::bfestim_nu_edge);
    int gridsave_nbfestim_in = 0;
    assert_always(fscanf(gridsave_file, "%d\n", &gridsave_nbfestim_in) == 1);
    assert_always(gridsave_nbfestim_in == bfestimcount);

    for (auto nonemptymgi = 0Z; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
      int nonemptymgi_in = 0;
      assert_always(fscanf(gridsave_file, "%d\n", &nonemptymgi_in) == 1);
      assert_always(nonemptymgi_in == nonemptymgi);
      for (int i = 0; i < bfestimcount; i++) {
        float bfrate_normed = 0;
        assert_always(fscanf(gridsave_file, "%a ", &bfrate_normed) == 1);

        if (globals::rank_in_node == 0) {
          prev_bfrate_normed[(nonemptymgi * bfestimcount) + i] = bfrate_normed;
        }
      }
    }
  }

  if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
    int detailed_linecount_in = 0;
    assert_always(fscanf(gridsave_file, "%d\n", &detailed_linecount_in) == 1);

    if (detailed_linecount_in != detailed_linecount) {
      printlnlog("ERROR: gridsave file specifies {} detailed lines but this simulation has {}.", detailed_linecount_in,
                 detailed_linecount);
      std::abort();
    }

    for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
      assert_always(fscanf(gridsave_file, "%d ", &detailed_lineindices[jblueindex]) == 1);
    }
  }

  for (auto nonemptymgi = 0Z; nonemptymgi < grid::get_nonempty_npts_model(); nonemptymgi++) {
    int nonemptymgi_in = 0;
    assert_always(fscanf(gridsave_file, "%d %la\n", &nonemptymgi_in, &J_normfactor[nonemptymgi]) == 2);
    assert_always(nonemptymgi_in == nonemptymgi);

    if constexpr (MULTIBIN_RADFIELD_MODEL_ON) {
      for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++) {
        const auto mgibinindex = (nonemptymgi * RADFIELDBINCOUNT) + binindex;
        float W = 0;
        float T_R = 0;
        assert_always(fscanf(gridsave_file, "%la %la %a %a %d\n", &radfieldbins.J_raw[mgibinindex],
                             &radfieldbins.nuJ_raw[mgibinindex], &W, &T_R,
                             &radfieldbins.contribcount[mgibinindex]) == 5);
        if (globals::rank_in_node == 0) {
          radfieldbin_solutions_W[mgibinindex] = W;
          radfieldbin_solutions_T_R[mgibinindex] = T_R;
        }
      }
    }

    if constexpr (DETAILED_LINE_ESTIMATORS_ON) {
      for (int jblueindex = 0; jblueindex < detailed_linecount; jblueindex++) {
        assert_always(fscanf(gridsave_file, "%la %d\n", &Jb_lu_raw[nonemptymgi][jblueindex].value,
                             &Jb_lu_raw[nonemptymgi][jblueindex].contribcount) == 2);
      }
    }
  }
  assert_always(fscanf(gridsave_file, "%d\n", &code_check) == 1);
  assert_always(code_check == 42809403);
}

}  // namespace radfield
