#ifndef LTEPOP_H
#define LTEPOP_H

#include <memory>
#include <vector>

#include "atomic.h"
#include "sn3d.h"

double get_groundlevelpop(int modelgridindex, int element, int ion);
double calculate_levelpop(int modelgridindex, int element, int ion, int level);
double calculate_levelpop_lte(int modelgridindex, int element, int ion, int level);
double get_levelpop(int modelgridindex, int element, int ion, int level);
double calculate_sahafact(int element, int ion, int level, int upperionlevel, double T, double E_threshold);
double get_nnion(int modelgridindex, int element, int ion);
void calculate_ion_balance_nne(int modelgridindex);
void calculate_cellpartfuncts(int modelgridindex, int element);
[[nodiscard]] auto calculate_ionfractions(const int element, const int modelgridindex, const double nne,
                                          const bool use_phi_lte) -> std::vector<double>;
void set_groundlevelpops(const int modelgridindex, const int element, const float nne, const bool force_lte);

#endif  // LTEPOP_H
