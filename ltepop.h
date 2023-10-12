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
double ionstagepop(int modelgridindex, int element, int ion);
void calculate_ion_balance_nne(int modelgridindex);
void calculate_cellpartfuncts(int modelgridindex, int element);
void set_groundlevelpops_if_needed(const int modelgridindex, const int element, const float nne, const bool force_lte);

#endif  // LTEPOP_H
