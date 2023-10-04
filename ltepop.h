#ifndef LTEPOP_H
#define LTEPOP_H

#include <memory>
#include <vector>

#include "atomic.h"
#include "sn3d.h"

double nne_solution_f(double x, void *paras);
std::vector<double> get_ionfractions(int element, int modelgridindex, double nne, int uppermost_ion);
double phi(int element, int ion, int modelgridindex);
double calculate_partfunct(int element, int ion, int modelgridindex);
double get_groundlevelpop(int modelgridindex, int element, int ion);
double calculate_levelpop_lte(int modelgridindex, int element, int ion, int level);
double get_levelpop(int modelgridindex, int element, int ion, int level);
double calculate_levelpop(int modelgridindex, int element, int ion, int level);
double calculate_sahafact(int element, int ion, int level, int upperionlevel, double T, double E_threshold);
double ionstagepop(int modelgridindex, int element, int ion);

#endif  // LTEPOP_H
