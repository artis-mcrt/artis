#ifndef RATECOEFF_H
#define RATECOEFF_H

void ratecoefficients_init(void);

__host__ __device__ double select_continuum_nu(int element, int ion, int level, int upperionlevel, float T_e);

__host__ __device__ double interpolate_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, double T);

__host__ __device__ double get_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, float T_e);
__host__ __device__ double get_stimrecombcoeff(int element, int lowerion, int level, int phixstargetindex, int modelgridindex);

__host__ __device__ double get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
__host__ __device__ double get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex);

__host__ __device__ double calculate_iongamma_per_gspop(int modelgridindex, int element, int ion);
__host__ __device__ double calculate_iongamma_per_ionpop(
  int modelgridindex, float T_e, int element, int lowerion,
  bool assume_lte, bool collisional_not_radiative, bool printdebug, bool use_bfest);

double calculate_ionrecombcoeff(
  int modelgridindex, float T_e,
  int element, int upperion,
  bool assume_lte, bool collisional_not_radiative, bool printdebug,
  bool lower_superlevel_only, bool per_groundmultipletpop, bool stimonly);


#if CUDA_ENABLED
template <double func_integrand(double, void *)>
__global__ void kernel_simpson_integral(void *intparas, double xlow, double deltax, int samplecount, double *integral)
/// Integrand to calculate the rate coefficient for photoionization
/// using gsl integrators. Corrected for stimulated recombination.
{
  extern __shared__ double threadcontrib[];

  const int i = blockIdx.x * blockDim.x + threadIdx.x;

  if (i < samplecount)
  {
    // Simpson's rule integral (will later be divided by 3)
    // n must be odd
    // integral = (xn - x0) / 3 * {f(x_0) + 4 * f(x_1) + 2 * f(x_2) + ... + 4 * f(x_1) + f(x_n-1)}
    // weights e.g., 1,4,2,4,2,4,1
    double weight;
    if (i == 0 || i == (samplecount - 1))
    {
      weight = 1.;
    }
    else if (i % 2 == 0)
    {
      weight = 2.;
    }
    else
    {
      weight = 4.;
    }

    const double x = xlow + deltax * i;

    threadcontrib[threadIdx.x] = weight * func_integrand(x, intparas) * deltax;
  }
  else
  {
    threadcontrib[threadIdx.x] = 0.;
  }

  __syncthreads();

  if (threadIdx.x == 0)
  {
    double blockcontrib = threadcontrib[0];
    for (int x = 1; x < blockDim.x; x++)
    {
      blockcontrib += threadcontrib[x];
    }
    atomicAdd(integral, blockcontrib / 3.); // divided by 3 for Simpson rule
  }
}


template <double func_integrand(double, void *), typename T>
__host__ __device__ double calculate_integral_gpu(T intparas, double xlow, double xhigh)
{
  T *dev_intparas;
  checkCudaErrors(cudaMalloc(&dev_intparas, sizeof(T)));

  #ifdef __CUDA_ARCH__
  memcpy(dev_intparas, &intparas, sizeof(T));
  // *dev_intparas = intparas;
  #else
  checkCudaErrors(cudaMemcpy(dev_intparas, &intparas, sizeof(T), cudaMemcpyHostToDevice));
  #endif

  double *dev_integral;
  cudaMalloc(&dev_integral, sizeof(double));

  #ifdef __CUDA_ARCH__
  *dev_integral = 0.;
  #else
  cudaMemset(dev_integral, 0, sizeof(double));
  #endif

  checkCudaErrors(cudaDeviceSynchronize());

  const int samplecount = globals::NPHIXSPOINTS * 16 + 1; // need an odd number for Simpson rule
  assert(samplecount % 2 == 1);

  dim3 threadsPerBlock(32, 1, 1);
  dim3 numBlocks((samplecount + threadsPerBlock.x - 1) / threadsPerBlock.x, 1, 1);
  size_t sharedsize = sizeof(double) * threadsPerBlock.x;

  const double deltax = (xhigh - xlow) / samplecount;

  kernel_simpson_integral<func_integrand><<<numBlocks, threadsPerBlock, sharedsize>>>(
    (void *) dev_intparas, xlow, deltax, samplecount, dev_integral);

  // Check for any errors launching the kernel
  checkCudaErrors(cudaGetLastError());

  // cudaDeviceSynchronize waits for the kernel to finish, and returns any errors encountered during the launch.
  checkCudaErrors(cudaDeviceSynchronize());

  #ifdef __CUDA_ARCH__
  const double result = *dev_integral;
  #else
  double result = 0.;
  checkCudaErrors(cudaMemcpy(&result, dev_integral, sizeof(double), cudaMemcpyDeviceToHost));
  #endif

  cudaFree(dev_integral);
  cudaFree(dev_intparas);

  return result;
}
#endif

extern __managed__ double T_step_log;

#endif //RATECOEFF_H