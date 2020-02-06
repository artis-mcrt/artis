// #include <math.h>
// #include <stdint.h>
// #include <stdio.h>
// #include "cuda_runtime.h"
#include "radfield.h"
#include "artisoptions.h"

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


__device__ double photoionization_crosssection_fromtable(float *photoion_xs, double nu_edge, double nu, int NPHIXSPOINTS, double NPHIXSNUINCREMENT)
/// Calculates the photoionisation cross-section at frequency nu out of the atomic data.
/// Input: - edge frequency nu_edge of the desired bf-continuum
///        - nu
{
  float sigma_bf;
  const double ireal = (nu / nu_edge - 1.0) / NPHIXSNUINCREMENT;
  const int i = floor(ireal);

  if (i < 0)
  {
    sigma_bf = 0.0;
  }
  else if (i < NPHIXSPOINTS - 1)
  {
    // sigma_bf = elements[element].ions[ion].levels[level].photoion_xs[i];

    const double sigma_bf_a = photoion_xs[i];
    const double sigma_bf_b = photoion_xs[i + 1];
    const double factor_b = ireal - i;
    sigma_bf = ((1. - factor_b) * sigma_bf_a) + (factor_b * sigma_bf_b);
  }
  else
  {
    const double last_phixs_nuovernuedge = (1.0 + NPHIXSNUINCREMENT * (NPHIXSPOINTS - 1));
    const double nu_max_phixs = nu_edge * last_phixs_nuovernuedge; //nu of the uppermost point in the phixs table
    sigma_bf = photoion_xs[NPHIXSPOINTS-1] * pow(nu_max_phixs / nu, 3);
  }

  return sigma_bf;
}


__device__ int select_bin(double nu, double *radfieldbin_nu_upper)
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

const int integralsamplesperxspoint = 8;

__global__ void kernel_corrphotoion_integral2(
  struct radfieldbin *radfieldbins_thiscell, double *radfieldbin_nu_upper, double nu_edge, float *photoion_xs,
  double departure_ratio, float T_e, float *integral, int NPHIXSPOINTS, double NPHIXSNUINCREMENT)
/// Integrand to calculate the rate coefficient for photoionization
/// using gsl integrators. Corrected for stimulated recombination.
{
  extern __shared__ double part_integral[];
  // __shared__ double part_integral[integralsamplesperxspoint * 100];

  if (threadIdx.x >= integralsamplesperxspoint || threadIdx.y >= NPHIXSPOINTS)
    return;

  // const double last_phixs_nuovernuedge = (1.0 + NPHIXSNUINCREMENT * (NPHIXSPOINTS - 1));

  const double nu = nu_edge * (1. + (NPHIXSNUINCREMENT * (threadIdx.y + (threadIdx.x / integralsamplesperxspoint))));

  const int binindex = select_bin(nu, radfieldbin_nu_upper);
  if (binindex < 0)
    return;
  const float bin_T_R = radfieldbins_thiscell[binindex].T_R;
  const float bin_W = radfieldbins_thiscell[binindex].W;
  // const double bin_nu_lower = binindex == 0 ? nu_lower_first_initial : radfieldbin_nu_upper[binindex - 1];
  // const double bin_nu_upper = radfieldbin_nu_upper[binindex];

  const double delta_nu = nu_edge * (NPHIXSNUINCREMENT / integralsamplesperxspoint);

  #if (SEPARATE_STIMRECOMB)
    const double corrfactor = 1.0;
  #else
    double corrfactor = 1. - departure_ratio * exp(-HOVERKB * nu / T_e);
    if (corrfactor < 0)
      corrfactor = 0.;
  #endif

  // printf("kernel_corrphotoion_integral: nu %lg binindex %d nu_lower %lg nu_upper %lg T_R %g W %g\n", nu, binindex, bin_nu_lower, bin_nu_upper, bin_T_R, bin_W);

  const double Jnu = bin_W * TWOHOVERCLIGHTSQUARED * pow(nu, 3) / expm1(HOVERKB * nu / bin_T_R);

  const float sigma_bf = photoionization_crosssection_fromtable(photoion_xs, nu_edge, nu, NPHIXSPOINTS, NPHIXSNUINCREMENT);

  part_integral[threadIdx.y * integralsamplesperxspoint + threadIdx.x] = ONEOVERH * sigma_bf / nu * Jnu * corrfactor * delta_nu;

  // printf("kernel_corrphotoion_integral sigma_bf %g T_e %g corrfactor %lg part_integral %g\n", sigma_bf, T_e, corrfactor, part_integral);

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
    *integral = total;
  }

  __syncthreads();
}


double radfield_gpu(double nu, int modelgridindex)
{
    cudaError_t cudaStatus;

    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed. CUDA-capable GPU installed?");
        abort();
    }

    // Launch a kernel on the GPU with one thread for each element.
    dim3 threadsPerBlock(RADFIELDBINCOUNT, 10, 1);
    dim3 numBlocks(1, 1, 1);

    double *radfieldjnu;

    cudaMallocManaged(&radfieldjnu, sizeof(double));
    *radfieldjnu = 0;

    kernel_radfield<<<numBlocks, threadsPerBlock>>>(nu, radfieldbins[modelgridindex], radfieldbin_nu_upper, radfieldjnu);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        abort();
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        abort();
    }

    double result = *radfieldjnu;
    cudaFree(radfieldjnu);
    return result;
}


double calculate_corrphotoioncoeff_integral_gpu(int modelgridindex, double nu_edge, float *photoion_xs, double departure_ratio, float T_e)
{
    cudaError_t cudaStatus;

    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed. CUDA-capable GPU installed?");
        abort();
    }

    float *integral;

    cudaMallocManaged(&integral, sizeof(double));
    *integral = 0;

    cudaStatus = cudaDeviceSynchronize();

    dim3 threadsPerBlock(integralsamplesperxspoint, NPHIXSPOINTS, 1);
    dim3 numBlocks(1, 1, 1);
    size_t sharedsize = sizeof(double) * NPHIXSPOINTS * integralsamplesperxspoint;

    kernel_corrphotoion_integral<<<numBlocks, threadsPerBlock, sharedsize>>>(
      radfieldbins[modelgridindex], radfieldbin_nu_upper, nu_edge, photoion_xs, departure_ratio, T_e, integral, NPHIXSPOINTS, NPHIXSNUINCREMENT);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "kernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        abort();
    }

    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        abort();
    }

    double result = *integral;
    cudaFree(integral);
    return result;
}

