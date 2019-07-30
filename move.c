#include "sn3d.h"
#include "move.h"
#include "radfield.h"
#include "rpkt.h"
#include "update_grid.h"
#include "vectors.h"


void update_estimators(const PKT *pkt_ptr, const double distance, const double t_current)
/// Update the volume estimators J and nuJ
/// This is done in another routine than move, as we sometimes move dummy
/// packets which do not contribute to the radiation field.
{
  const int cellindex = pkt_ptr->where;
  const int modelgridindex = cell[cellindex].modelgridindex;

  /// Update only non-empty cells
  if (modelgridindex != MMODELGRID)
  {
    const double distance_e_cmf = distance * pkt_ptr->e_cmf;
    const double nu = pkt_ptr->nu_cmf;
    //double bf = exp(-HOVERKB*nu/cell[modelgridindex].T_e);

    radfield_update_estimators(modelgridindex, distance_e_cmf, nu, pkt_ptr, t_current);

    #ifndef FORCE_LTE
      ///ffheatingestimator does not depend on ion and element, so an array with gridsize is enough.
      ///quick and dirty solution: store info in element=ion=0, and leave the others untouched (i.e. zero)
      #ifdef _OPENMP
        #pragma omp atomic
      #endif
      ffheatingestimator[modelgridindex] += distance_e_cmf * kappa_rpkt_cont[tid].ffheating;

      #if (!NO_LUT_PHOTOION || !NO_LUT_BFHEATING)
        #if (!NO_LUT_PHOTOION)
        const double distance_e_cmf_over_nu = distance_e_cmf / nu;
        #endif
        for (int i = 0; i < nbfcontinua_ground; i++)
        {
          const double nu_edge = phixslist[tid].groundcont[i].nu_edge;
          if (nu > nu_edge)
          {
            const int element = phixslist[tid].groundcont[i].element;
            const int ion = phixslist[tid].groundcont[i].ion;
            /// Cells with zero abundance for a specific element have zero contribution
            /// (set in calculate_kappa_rpkt_cont and therefore do not contribute to
            /// the estimators
            if (get_abundance(modelgridindex, element) > 0)
            {
              const int ionestimindex = modelgridindex * nelements * maxion + element * maxion + ion;
              #if (!NO_LUT_PHOTOION)
                #ifdef _OPENMP
                  #pragma omp atomic
                #endif
                gammaestimator[ionestimindex] += phixslist[tid].groundcont[i].gamma_contr * distance_e_cmf_over_nu;

                #ifdef DEBUG_ON
                if (!isfinite(gammaestimator[ionestimindex]))
                {
                  printout("[fatal] update_estimators: gamma estimator becomes non finite: level %d, gamma_contr %g, distance_e_cmf_over_nu %g\n", i, phixslist[tid].groundcont[i].gamma_contr, distance_e_cmf_over_nu);
                  abort();
                }
                #endif
              #endif
              #if (!NO_LUT_BFHEATING)
                #ifdef _OPENMP
                  #pragma omp atomic
                #endif
                bfheatingestimator[ionestimindex] += phixslist[tid].groundcont[i].gamma_contr * distance_e_cmf * (1. - nu_edge/nu);
                //bfheatingestimator[ionestimindex] += phixslist[tid].groundcont[i].bfheating_contr * distance_e_cmf * (1/nu_edge - 1/nu);
              #endif
            }
          }
          else
            break; // because groundcont is sorted by nu_edge, nu < nu_edge for all remaining items
        }
      #endif

    #endif

    ///Heating estimators. These are only applicable for pure H. Other elements
    ///need advanced treatment in thermalbalance calculation.
    //cell[pkt_ptr->where].heating_ff += distance_e_cmf * kappa_rpkt_cont[tid].ffheating;
    //cell[pkt_ptr->where].heating_bf += distance_e_cmf * kappa_rpkt_cont[tid].bfheating;
  }

}


void move_pkt(PKT *restrict pkt_ptr, const double distance, const double time)
/// Subroutine to move a packet along a straight line (specified by currect
/// dir vector). The distance moved is in the rest frame. Time must be the
/// time at the end of distance travelled.
{
  /// First update pos.
  assert(distance >= 0);

  pkt_ptr->pos[0] += (pkt_ptr->dir[0] * distance);
  pkt_ptr->pos[1] += (pkt_ptr->dir[1] * distance);
  pkt_ptr->pos[2] += (pkt_ptr->dir[2] * distance);

  /// During motion, rest frame energy and frequency are conserved.
  /// But need to update the co-moving ones.
  const double dopplerfactor = doppler_packetpos(pkt_ptr, time);
  pkt_ptr->nu_cmf = pkt_ptr->nu_rf * dopplerfactor;
  pkt_ptr->e_cmf = pkt_ptr->e_rf * dopplerfactor;
}
