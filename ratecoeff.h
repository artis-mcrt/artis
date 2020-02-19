#ifndef RATECOEFF_H
#define RATECOEFF_H

#include <stdbool.h>

void ratecoefficients_init(void);

double select_continuum_nu(int element, int ion, int level, int upperionlevel, float T_e);

double interpolate_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, double T);

__host__ __device__ double get_spontrecombcoeff(int element, int ion, int level, int phixstargetindex, float T_e);
double get_stimrecombcoeff(int element, int lowerion, int level, int phixstargetindex, int modelgridindex);

double get_corrphotoioncoeff(int element, int ion, int level, int phixstargetindex, int modelgridindex);
double get_corrphotoioncoeff_ana(int element, int ion, int level, int phixstargetindex, int modelgridindex);

double calculate_iongamma_per_gspop(int modelgridindex, int element, int ion);
double calculate_iongamma_per_ionpop(
  int modelgridindex, float T_e, int element, int lowerion,
  bool assume_lte, bool collisional_not_radiative, bool printdebug, bool use_bfest);

double calculate_ionrecombcoeff(
  int modelgridindex, float T_e,
  int element, int upperion,
  bool assume_lte, bool collisional_not_radiative, bool printdebug,
  bool lower_superlevel_only, bool per_groundmultipletpop, bool stimonly);


#if CUDA_ENABLED
const int integralsamplesperxspoint = 4; // must be an even number for Simpsons rule to work

template <double func_integrand(double, void *)>
__global__ void kernel_integral(void *intparas, double nu_edge, double *integral)
/// Integrand to calculate the rate coefficient for photoionization
/// using gsl integrators. Corrected for stimulated recombination.
{
  extern __shared__ double part_integral[];
  // __shared__ double part_integral[integralsamplesperxspoint * 100];

  if (threadIdx.x < integralsamplesperxspoint && threadIdx.y < NPHIXSPOINTS)
  {
    // const double last_phixs_nuovernuedge = (1.0 + NPHIXSNUINCREMENT * (NPHIXSPOINTS - 1));

    const double nu = nu_edge * (1. + (NPHIXSNUINCREMENT * (threadIdx.y + (threadIdx.x / integralsamplesperxspoint))));

    const int sampleindex = threadIdx.y * integralsamplesperxspoint + threadIdx.x;
    const int lastsampleindex = (NPHIXSPOINTS - 1) * integralsamplesperxspoint + (integralsamplesperxspoint - 1);

    // Simpson's rule integral (will later be divided by 3)
    // n must be odd
    // integral = (xn - x0) / 3 * {f(x_0) + 4 * f(x_1) + 2 * f(x_2) + ... + 4 * f(x_1) + f(x_n-1)}
    // weights e.g., 1,4,2,4,2,4,1
    double weight = 0.;
    if (sampleindex == 0 || sampleindex == lastsampleindex)
    {
      weight = 1.;
    }
    else if (sampleindex % 2 == 0)
    {
      weight = 2.;
    }
    else
    {
      weight = 4.;
    }

    const double delta_nu = nu_edge * (NPHIXSNUINCREMENT / integralsamplesperxspoint);

    const double integrand = func_integrand(nu, intparas);

    part_integral[sampleindex] = weight * integrand * delta_nu;
  }

  __syncthreads();

  if (threadIdx.x == 0)
  {
    for (unsigned int x = 1; x < integralsamplesperxspoint; x++)
    {
      const int firstsampleindex = threadIdx.y * integralsamplesperxspoint; // first sample of the cross section point
      part_integral[firstsampleindex] += part_integral[firstsampleindex + x];
    }
  }

  __syncthreads();

  if (threadIdx.x == 0 && threadIdx.y == 0)
  {
    double total = 0.;
    for (unsigned int y = 0; y < NPHIXSPOINTS; y++)
    {
      total += part_integral[y * integralsamplesperxspoint];
    }
    atomicAdd(integral, total / 3.);
  }
}


template <double func_integrand(double, void *)>
double calculate_phixs_integral_gpu(void *dev_intparas, double nu_edge)
{
    double *dev_integral;
    checkCudaErrors(cudaMalloc(&dev_integral, sizeof(double)));
    cudaMemset(dev_integral, 0, sizeof(double));

    checkCudaErrors(cudaDeviceSynchronize());

    dim3 threadsPerBlock(integralsamplesperxspoint, NPHIXSPOINTS, 1);
    dim3 numBlocks(1, 1, 1);
    size_t sharedsize = sizeof(double) * threadsPerBlock.x * threadsPerBlock.y;

    kernel_integral<func_integrand><<<numBlocks, threadsPerBlock, sharedsize>>>(dev_intparas, nu_edge, dev_integral);

    // Check for any errors launching the kernel
    checkCudaErrors(cudaGetLastError());

    // cudaDeviceSynchronize waits for the kernel to finish, and returns any errors encountered during the launch.
    checkCudaErrors(cudaDeviceSynchronize());

    double result;

    checkCudaErrors(cudaMemcpy(&result, dev_integral, sizeof(double), cudaMemcpyDeviceToHost));

    cudaFree(dev_integral);

    return result;
}
#endif


extern __managed__ double T_step_log;

#endif //RATECOEFF_H