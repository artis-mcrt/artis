#ifndef ATOMIC_CUDA_CUH
#define ATOMIC_CUDA_CUH

__device__ double photoionization_crosssection_fromtable_gpu(float *photoion_xs, double nu_edge, double nu, int NPHIXSPOINTS, double NPHIXSNUINCREMENT);

#endif //ATOMIC_CUDA_CUH
