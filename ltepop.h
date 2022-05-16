#ifndef LTEPOP_H
#define LTEPOP_H

#include "atomic.h"
#include "sn3d.h"

double nne_solution_f(double x, void *paras);
void get_ionfractions(int element, int modelgridindex, double nne, double ionfractions[], int uppermost_ion);
double phi(int element, int ion, int modelgridindex);
double calculate_partfunct(int element, int ion, int modelgridindex);
__host__ __device__ double get_groundlevelpop(int modelgridindex, int element, int ion);
__host__ __device__ double calculate_levelpop_lte(int modelgridindex, int element, int ion, int level);
__host__ __device__ double get_levelpop(int modelgridindex, int element, int ion, int level);
__host__ __device__ double calculate_levelpop(int modelgridindex, int element, int ion, int level);
__host__ __device__ double get_groundmultiplet_pop(int modelgridindex, int element, int ion);
__host__ __device__ double calculate_sahafact(int element, int ion, int level, int upperionlevel, double T, double E_threshold);
__host__ __device__ double ionstagepop(int modelgridindex, int element, int ion);



#endif //LTEPOP_H
