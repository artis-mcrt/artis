#ifndef THREADPRIVATE_H
#define THREADPRIVATE_H

#include <stdbool.h>

int tid;
//int histindex;
bool use_cellhist;
bool neutral_flag;
gsl_rng *rng;

//int debuglevel;
//int propagationcounter;
//double J[MMODELGRID+1];                        /// Volume estimator for the frequency integrated mean intensity J  ///ATOMIC
//#ifndef FORCE_LTE
//  double nuJ[MMODELGRID+1];                      /// Volume estimator for the frequency integration of nu*J_nu       ///ATOMIC
//#endif
//struct_mastate mastate;
//rpkt_cont_opacity_struct kappa_rpkt_cont;
//coolingrates_t coolingrates;
//heatingrates_t heatingrates;
//short elements_uppermost_ion[MELEMENTS];  /// Highest ionisation stage which has a decent population for a particular element
//                                                 /// in a given cell. Be aware that this must not be used outside of the update_grid
//                                                 /// routine and their doughters.

//cellhistory_struct *cellhistory;
//phixslist_t *phixslist;
//groundphixslist_t *groundphixslist;



//float rhosum[MSAMPLEGRID];
//float T_Rsum[MSAMPLEGRID];
//float T_esum[MSAMPLEGRID];
//float Wsum[MSAMPLEGRID];
//float T_Dsum[MSAMPLEGRID];
//int associatedcells[MSAMPLEGRID];

FILE *restrict output_file;
//short output_file_open;

#endif