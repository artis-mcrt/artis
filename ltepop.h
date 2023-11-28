#ifndef LTEPOP_H
#define LTEPOP_H

#include <memory>
#include <vector>

#include "atomic.h"
#include "sn3d.h"

auto get_groundlevelpop(int modelgridindex, int element, int ion) -> double;
auto calculate_levelpop(int modelgridindex, int element, int ion, int level) -> double;
auto calculate_levelpop_lte(int modelgridindex, int element, int ion, int level) -> double;
auto get_levelpop(int modelgridindex, int element, int ion, int level) -> double;
[[nodiscard]] auto calculate_sahafact(int element, int ion, int level, int upperionlevel, double T, double E_threshold)
    -> double;
[[nodiscard]] auto get_nnion(int modelgridindex, int element, int ion) -> double;
void calculate_ion_balance_nne(int modelgridindex);
void calculate_cellpartfuncts(int modelgridindex, int element);
[[nodiscard]] auto calculate_ionfractions(int element, int modelgridindex, double nne, bool use_phi_lte)
    -> std::vector<double>;
void set_groundlevelpops(int modelgridindex, int element, float nne, bool force_lte);

#endif  // LTEPOP_H
