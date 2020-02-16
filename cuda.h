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

// #include <iostream>

// #define checkCudaErrors( call )\
// {\
//     cudaError_t result = call;\
//     if (cudaSuccess != result)\
//     {\
//         std::cerr << "CUDA error " << result << " in " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString( result ) << " (" << #call << ")" << std::endl; \
//     }\
// }

// #define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)
//
// template <typename T>
// void check(T result, char const *const func, const char *const file,
//            int const line) {
//   if (result) {
//     fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line,
//             static_cast<unsigned int>(result), _cudaGetErrorEnum(result), func);
//     // DEVICE_RESET
//     // Make sure we call CUDA Device Reset before exiting
//     exit(EXIT_FAILURE);
//   }
// }

#else

#include <assert.h>

#if CUDA_ENABLED
assert(false);
#endif

#define __managed__
#define __device__
#define __host__

#endif

#endif