#ifndef MACROATOM_H
#define MACROATOM_H

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

typedef struct mastate_t
{
  int element;              /// macro atom of type element (this is an element index)
  int ion;                  /// in ionstage ion (this is an ion index)
  int level;                /// and level=level (this is a level index)
  int activatingline;       /// Linelistindex of the activating line for bb activated MAs, -99 else.
} mastate_t;


#include "cuda.h"

void macroatom_open_file(const int my_rank);
void macroatom_close_file(void);

__host__ __device__ void do_macroatom(struct packet *pkt_ptr, int timestep);

__host__ __device__ double rad_deexcitation_ratecoeff(int modelgridindex, int element, int ion, int upper, int lower, double epsilon_trans, int lineindex, double t_current);
__host__ __device__ double rad_excitation_ratecoeff(int modelgridindex, int element, int ion, int lower, int upper, double epsilon_trans, int lineindex, double t_current);
__host__ __device__ double rad_recombination_ratecoeff(float T_e, float nne, int element, int ion, int upper, int lower, int modelgridindex);
__host__ __device__ double stim_recombination_ratecoeff(float nne, int element, int upperion, int upper, int lower, int modelgridindex);

__host__ __device__ double col_deexcitation_ratecoeff(float T_e, float nne, double epsilon_trans, int lineindex);
__host__ __device__ double col_excitation_ratecoeff(float T_e, float nne, int lineindex, double epsilon_trans);
__host__ __device__ double col_recombination_ratecoeff(int modelgridindex, int element, int upperion, int upper, int lower, double epsilon_trans);
__host__ __device__ double col_ionization_ratecoeff(float T_e, float nne, int element, int ion, int lower, int phixstargetindex, double epsilon_trans);

#endif //MACROATOM_H
