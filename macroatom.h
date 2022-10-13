#ifndef MACROATOM_H
#define MACROATOM_H

#include <cmath>

enum ma_action {
  /// Radiative deexcitation rate from this level.
  MA_ACTION_RADDEEXC = 0,
  /// Collisional deexcitation rate from this level.
  MA_ACTION_COLDEEXC = 1,
  /// Radiative recombination from this level.
  MA_ACTION_RADRECOMB = 2,
  /// Collisional recombination rate from this level.
  MA_ACTION_COLRECOMB = 3,
  /// Rate for internal downward transitions to same ionisation stage.
  MA_ACTION_INTERNALDOWNSAME = 4,
  /// Rate for internal upward transitions to same ionisation stage.
  MA_ACTION_INTERNALDOWNLOWER = 5,
  /// Rate for internal downward transitions to lower ionisation stage.
  MA_ACTION_INTERNALUPSAME = 6,
  /// Rate for internal upward transitions to higher ionisation stage.
  MA_ACTION_INTERNALUPHIGHER = 7,
  /// Rate for internal upward transitions to higher ionisation stage due to non-thermal collisions.
  MA_ACTION_INTERNALUPHIGHERNT = 8,
  MA_ACTION_COUNT = 9,
};

#include "atomic.h"
#include "constants.h"
#include "globals.h"
#include "packet.h"

void macroatom_open_file(const int my_rank);
void macroatom_close_file(void);

__host__ __device__ void do_macroatom(struct packet *pkt_ptr, int timestep);

__host__ __device__ double rad_deexcitation_ratecoeff(int modelgridindex, int element, int ion, int upper, int lower,
                                                      double epsilon_trans, int lineindex, double t_current);
__host__ __device__ double rad_excitation_ratecoeff(int modelgridindex, int element, int ion, int lower, int upper,
                                                    double epsilon_trans, int lineindex, double t_current);
__host__ __device__ double rad_recombination_ratecoeff(float T_e, float nne, int element, int ion, int upper, int lower,
                                                       int modelgridindex);
__host__ __device__ double stim_recombination_ratecoeff(float nne, int element, int upperion, int upper, int lower,
                                                        int modelgridindex);

__host__ __device__ double col_recombination_ratecoeff(int modelgridindex, int element, int upperion, int upper,
                                                       int lower, double epsilon_trans);
__host__ __device__ double col_ionization_ratecoeff(float T_e, float nne, int element, int ion, int lower,
                                                    int phixstargetindex, double epsilon_trans);

__host__ __device__ constexpr double col_deexcitation_ratecoeff(const float T_e, const float nne,
                                                                const double epsilon_trans,
                                                                struct linelist_entry *const line,
                                                                double lowerstatweight, double upperstatweight)
// multiply by upper level population to get a rate per second
{
  double C;
  const double coll_str_thisline = line->coll_str;
  if (coll_str_thisline < 0) {
    const bool forbidden = line->forbidden;
    if (!forbidden)  // alternative: (coll_strength > -1.5) i.e. to catch -1
    {
      /// permitted E1 electric dipole transitions
      /// collisional deexcitation: formula valid only for atoms!!!!!!!!!!!
      /// Rutten script eq. 3.33. p.50
      // f = osc_strength(element,ion,upper,lower);
      // C = n_u * 2.16 * pow(fac1,-1.68) * pow(T_e,-1.5) *
      // stat_weight(element,ion,lower)/stat_weight(element,ion,upper)  * nne * f;

      const double eoverkt = epsilon_trans / (KB * T_e);
      /// Van-Regemorter formula, Mihalas (1978), eq.5-75, p.133
      const double g_bar = 0.2;  /// this should be read in from transitions data: it is 0.2 for transitions nl -> n'l'
                                 /// and 0.7 for transitions nl -> nl'
      // test = 0.276 * exp(fac1) * gsl_sf_expint_E1(fac1);
      /// crude approximation to the already crude Van-Regemorter formula

      // double test = 0.276 * exp(fac1) * (-0.5772156649 - log(fac1));
      // double Gamma = (g_bar > test) ? g_bar : test;

      // optimisation
      const double gauntfac = (eoverkt > 0.33421) ? g_bar : 0.276 * exp(eoverkt) * (-0.5772156649 - log(eoverkt));

      const double g_ratio = lowerstatweight / upperstatweight;

      C = C_0 * 14.51039491 * nne * sqrt(T_e) * line->osc_strength * pow(H_ionpot / epsilon_trans, 2) * eoverkt *
          g_ratio * gauntfac;
    } else  // alterative: (coll_strength > -3.5) to catch -2 or -3
    {
      // forbidden transitions: magnetic dipole, electric quadropole...
      // could be Axelrod? or Maurer
      C = nne * 8.629e-6 * 0.01 * lowerstatweight / sqrt(T_e);
    }
  } else  // positive values are treated as effective collision strengths
  {
    // from Osterbrock and Ferland, p51
    // statweight_target is LOWER LEVEL stat weight
    C = nne * 8.629e-6 * coll_str_thisline / upperstatweight / sqrt(T_e);
    // test test
    // C = n_u * nne * 8.629e-6 * pow(T_e,-0.5) * 0.01 * statweight_target;
  }

  return C;
}

__host__ __device__ constexpr double col_excitation_ratecoeff(const float T_e, const float nne,
                                                              struct linelist_entry *const line,
                                                              const double epsilon_trans, const double lowerstatweight,
                                                              const double upperstatweight)
// multiply by lower level population to get a rate per second
{
  double C;
  const double coll_strength = line->coll_str;
  const double eoverkt = epsilon_trans / (KB * T_e);

  if (coll_strength < 0) {
    const bool forbidden = line->forbidden;
    if (!forbidden)  // alternative: (coll_strength > -1.5) i.e. to catch -1
    {
      /// permitted E1 electric dipole transitions
      /// collisional excitation: formula valid only for atoms!!!!!!!!!!!
      /// Rutten script eq. 3.32. p.50
      // C = n_l * 2.16 * pow(eoverkt,-1.68) * pow(T_e,-1.5) * exp(-eoverkt) * nne *
      // osc_strength(element,ion,upper,lower);

      // Van-Regemorter formula, Mihalas (1978), eq.5-75, p.133
      const double g_bar = 0.2;  // this should be read in from transitions data: it is 0.2 for transitions nl -> n'l'
                                 // and 0.7 for transitions nl -> nl'
      // test = 0.276 * exp(eoverkt) * gsl_sf_expint_E1(eoverkt);
      /// crude approximation to the already crude Van-Regemorter formula
      const double exp_eoverkt = exp(eoverkt);

      const double test = 0.276 * exp_eoverkt * (-0.5772156649 - std::log(eoverkt));
      const double Gamma = g_bar > test ? g_bar : test;
      C = C_0 * nne * std::sqrt(T_e) * 14.51039491 * line->osc_strength * pow(H_ionpot / epsilon_trans, 2) * eoverkt /
          exp_eoverkt * Gamma;
    } else  // alterative: (coll_strength > -3.5) to catch -2 or -3
    {
      // forbidden transitions: magnetic dipole, electric quadropole...
      // Axelrod's approximation (thesis 1980)
      C = nne * 8.629e-6 * 0.01 * std::exp(-eoverkt) * upperstatweight / std::sqrt(T_e);
    }
  } else {
    // from Osterbrock and Ferland, p51
    C = nne * 8.629e-6 * coll_strength * std::exp(-eoverkt) / lowerstatweight / std::sqrt(T_e);
  }

  return C;
}

#endif  // MACROATOM_H
