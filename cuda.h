#ifndef CUDA_H
#define CUDA_H

#ifndef CUDA_ENABLED
#define CUDA_ENABLED false
#endif

#ifdef __CUDACC__

#define USECUDA_BFHEATING false
#define USECUDA_PHOTOIONCOEFF false
#define USECUDA_NLTE_BOUNDBOUND false
#define USECUDA_NONTHERMAL_EXCITATION true
#define USECUDA_NONTHERMAL_IONIZATION true

#include <cuda_runtime.h>

#include <helper_cuda.h>

// copied from https://www.micc.unifi.it/bertini/download/gpu-programming-basics/2017/gpu_cuda_5.pdf
#if defined __CUDA_ARCH__ && __CUDA_ARCH__ < 600
__device__ double static inline atomicAdd(double* address, double val)
{
 unsigned long long int* address_as_ull =
 (unsigned long long int*) address;
 unsigned long long int old = *address_as_ull;
 unsigned long long int assumed;
 do {
 assumed = old;
 old = atomicCAS(address_as_ull, assumed,
 __double_as_longlong(val + __longlong_as_double(assumed)));
 // Note: uses integer comparison to avoid hang in case
 // of NaN (since NaN != NaN)
 } while (assumed != old);
 return __longlong_as_double(old);
}
#endif

#else

#include <assert.h>

assert(!CUDA_ENABLED);

#define __managed__
#define __device__
#define __host__

#endif

#endif