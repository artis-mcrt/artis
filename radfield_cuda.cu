#include <assert.h>
#include "artisoptions.h"
#include "atomic.h"
#include "atomic_cuda.cuh"
#include "radfield.h"

extern double *dev_radfieldbin_nu_upper;
extern struct radfieldbin **dev_radfieldbins;


__global__ void kernel_radfield(double nu, struct radfieldbin *radfieldbins_thiscell, double *radfieldbin_nu_upper, double *radfieldjnu)
{
    const int binindex = threadIdx.x + blockIdx.x * blockDim.x;
    const float bin_T_R = radfieldbins_thiscell[binindex].T_R;
    const float bin_W = radfieldbins_thiscell[binindex].W;
    const double bin_nu_lower = binindex == 0 ? nu_lower_first_initial : radfieldbin_nu_upper[binindex - 1];
    const double bin_nu_upper = radfieldbin_nu_upper[binindex];
    if (bin_nu_lower <= nu && bin_nu_upper > nu)
    {
      // printf("CUDAkernel: nu %lg binindex %d nu_lower %lg nu_upper %lg T_R %g W %g\n", nu, binindex, bin_nu_lower, bin_nu_upper, bin_T_R, bin_W);
      *radfieldjnu = bin_W * TWOHOVERCLIGHTSQUARED * pow(nu, 3) / expm1(HOVERKB * nu / bin_T_R);
      // printf("    radfieldjnu %g\n", *radfieldjnu);
    }
}


__device__ int select_bin_gpu(double nu, double *radfieldbin_nu_upper)
{
  // linear search one by one until found
  if (nu >= radfieldbin_nu_upper[RADFIELDBINCOUNT - 1])
    return -1; // out of range, nu higher than highest bin
  else if (nu < nu_lower_first_initial)
    return -2; // out of range, nu lower than lowest bin
  else
  {
    for (int binindex = 0; binindex < RADFIELDBINCOUNT; binindex++)
    {
      if (radfieldbin_nu_upper[binindex] > nu)
      {
        return binindex;
      }
    }

    return -3;
  }
}

const int integralsamplesperxspoint = 8; // must be an even number for Simpsons rule to work

__global__ void kernel_corrphotoion_integral(
  struct radfieldbin *radfieldbins_thiscell, double *radfieldbin_nu_upper, double nu_edge, float *photoion_xs,
  double departure_ratio, float T_e, double *integral, int NPHIXSPOINTS, double NPHIXSNUINCREMENT)
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

    const int binindex = select_bin_gpu(nu, radfieldbin_nu_upper);

    if (binindex < 0)
    {
      part_integral[sampleindex] = 0.;
    }
    else
    {
      const float bin_T_R = radfieldbins_thiscell[binindex].T_R;
      const float bin_W = radfieldbins_thiscell[binindex].W;
      // const double bin_nu_lower = binindex == 0 ? nu_lower_first_initial : radfieldbin_nu_upper[binindex - 1];
      // const double bin_nu_upper = radfieldbin_nu_upper[binindex];

      const double Jnu = bin_W * TWOHOVERCLIGHTSQUARED * pow(nu, 3) / expm1(HOVERKB * nu / bin_T_R);

      const double delta_nu = nu_edge * (NPHIXSNUINCREMENT / integralsamplesperxspoint);

      #if (SEPARATE_STIMRECOMB)
        const double corrfactor = 1.0;
      #else
        double corrfactor = 1. - departure_ratio * exp(-HOVERKB * nu / T_e);
        if (corrfactor < 0)
          corrfactor = 0.;
      #endif

      const float sigma_bf = photoionization_crosssection_fromtable_gpu(photoion_xs, nu_edge, nu, NPHIXSPOINTS, NPHIXSNUINCREMENT);

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

      part_integral[sampleindex] = weight * ONEOVERH * sigma_bf / nu * Jnu * corrfactor * delta_nu;
    }
  }

  __syncthreads();

  if (threadIdx.x == 0)
  {
    for (unsigned int x = 1; x < integralsamplesperxspoint; x++)
    {
      const int firstsampleindex = threadIdx.y * integralsamplesperxspoint;
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
    *integral = total / 3.;
  }

  __syncthreads();
}


double calculate_corrphotoioncoeff_integral_gpu(int modelgridindex, double nu_edge, float *photoion_xs, double departure_ratio, float T_e)
{
    cudaError_t cudaStatus;

    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed. CUDA-capable GPU installed?");
        abort();
    }

    void *dev_integral;

    cudaStatus = cudaMalloc(&dev_integral, sizeof(double));
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaDeviceSynchronize();
    assert(cudaStatus == cudaSuccess);

    dim3 threadsPerBlock(integralsamplesperxspoint, NPHIXSPOINTS, 1);
    dim3 numBlocks(1, 1, 1);
    size_t sharedsize = sizeof(double) * NPHIXSPOINTS * integralsamplesperxspoint;

    kernel_corrphotoion_integral<<<numBlocks, threadsPerBlock, sharedsize>>>(
      dev_radfieldbins[modelgridindex], dev_radfieldbin_nu_upper, nu_edge, photoion_xs, departure_ratio, T_e, (double *) dev_integral, NPHIXSPOINTS, NPHIXSNUINCREMENT);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    assert(cudaStatus == cudaSuccess);

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    assert(cudaStatus == cudaSuccess);

    double result;

    cudaStatus = cudaMemcpy(&result, dev_integral, sizeof(double), cudaMemcpyDeviceToHost);
    assert(cudaStatus == cudaSuccess);

    cudaFree(dev_integral);

    return result;
}
