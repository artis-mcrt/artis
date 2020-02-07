#include <assert.h>
// #include "artisoptions.h"


__device__ double photoionization_crosssection_fromtable_gpu(float *photoion_xs, double nu_edge, double nu, int NPHIXSPOINTS, double NPHIXSNUINCREMENT)
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

