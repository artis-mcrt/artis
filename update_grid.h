// Declarations for the per-timestep grid update (update_grid.cc).

#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

#include <chrono>
#include <ostream>

void update_grid(std::ostream& estimators_file, int nts, int nts_prev, int titer,
                 std::chrono::steady_clock::time_point real_time_start);

#endif  // UPDATE_GRID_H
