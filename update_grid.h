// Declarations for the per-timestep grid update (update_grid.cc).

#ifndef UPDATE_GRID_H
#define UPDATE_GRID_H

#include <chrono>

class JobFolderOutputFile;

void update_grid(JobFolderOutputFile& estimators_file, int nts, int nts_prev,
                 std::chrono::steady_clock::time_point real_time_start);

#endif  // UPDATE_GRID_H
