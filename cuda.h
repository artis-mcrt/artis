#ifndef CUDA_H
#define CUDA_H

#ifndef CUDA_ENABLED
#define CUDA_ENABLED false
#endif

#if CUDA_ENABLED

#include <cuda_runtime.h>


#define checkCudaErrors( call )        \
{                                      \
cudaError_t result = call;             \
if (cudaSuccess != result)             \
    std::cerr << "CUDA error " << result << " in " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString( result ) << " (" << #call << ")" << std::endl; abort(); \
}


#else

#define __managed__
#define __device__
#define __host__

#endif

#endif