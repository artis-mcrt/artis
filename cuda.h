#ifndef CUDA_H
#define CUDA_H

#ifndef CUDA_ENABLED
#define CUDA_ENABLED false
#endif

#ifdef __CUDACC__

  #define USECUDA_BFHEATING false
  #define USECUDA_PHOTOIONCOEFF true // must on if USECUDA_UPDATEPACKETS is on beccause no device GSL integration
  #define USECUDA_STIMRECOMBCOEFF true
  #define USECUDA_NLTE_BOUNDBOUND false
  #define USECUDA_NONTHERMAL_EXCITATION true
  #define USECUDA_NONTHERMAL_IONIZATION true
  #define USECUDA_RPKT_CONTOPACITY false
  #define CUDA_VERIFY_CPUCONSISTENCY false

  #define USECUDA_UPDATEPACKETS true

  #include <cuda_runtime.h>
  // #include <curand.h>

  #ifndef __CUDA_ARCH__
    #include "helper_cuda.h"
  #endif
  #include "globals.h"

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

  // #ifdef __CUDA_ARCH__
  // __device__
  // inline static void checkCudaErrors(cudaError_t cudaStatus)
  // {
  //   assert(cudaStatus == cudaSuccess);
  // }
  // #endif

  #ifdef __CUDA_ARCH__
    __device__ inline static void devcheck(cudaError_t result, char const *const func, const char *const file,
               int const line) {
      if (result != cudaSuccess)
      {
        printf("CUDA error at %s:%d code=%d \"%s\" \n", file, line,
                static_cast<unsigned int>(result), func);
        assert(result == cudaSuccess);
      }
    }

    #undef checkCudaErrors
    #define checkCudaErrors(val) devcheck((val), #val, __FILE__, __LINE__)

    #define abort() assert(false)

    extern __managed__ curandState curandstates[MTHREADS];

    // #define gsl_rng_uniform(Y) curand_uniform_double(&curandstates[threadIdx.x + blockDim.x * blockIdx.x])
    __device__ inline static double gsl_rng_uniform(void *ignore)
    {
      const int idx = threadIdx.x + blockIdx.x * blockDim.x;
      const double zrand = curand_uniform_double(&curandstates[idx]);
      // printf("random %g\n", zrand);
      return zrand;
    }

    __device__ inline static double gsl_rng_uniform_pos(void *ignore)
    {
      const int idx = threadIdx.x + blockIdx.x * blockDim.x;
      while (true)
      {
        const double zrand = curand_uniform_double(&curandstates[idx]);
        if (zrand > 0.)
          return zrand;
      }
    }
    // typedef void *gsl_error_handler_t;

    // typedef struct gsl_function
    // {
    //   double *function;
    //   void *params;
    // }
    //
    // __device__ inline static int customintegration(
    //   gsl_function *gslfunc, double xlow, double xhigh, double epsabs, double epsrel,
    //   size_t wsize, int integtype, void *gslworkspace, double *integral, double *error)
    // {
    //   *integral = calculate_integral_gpu<(gslfunc->function), gsl_integral_paras_gammacorr>(
    //       gslfunc->params, xlow, xhigh);
    //   return 0;
    // }
    // #define customintegration(PARATYPE, GSLF, F, PARAS, A, B, EA, ER, WS, T, W, I, E) 0; *I = calculate_integral_gpu<F, PARATYPE>(&PARAS, A, B)

    // #define gsl_set_error_handler(x) void;
  #endif

#else

  #if CUDA_ENABLED
  #error CUDA is not available
  #endif

  #define __managed__
  #define __device__
  #define __host__
  // #define customintegration(PARATYPE, GSLF, F, PARAS, A, B, EA, ER, WS, T, W, I, E) gsl_integration_qag(GSLF, A, B, EA, ER, WS, T, W, I, E)
#endif

#endif