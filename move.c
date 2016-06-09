#include "sn3d.h"
#include "move.h"
#include "radfield.h"
#include "update_grid.h"
#include "vectors.h"


void update_estimators(const PKT *restrict pkt_ptr, double distance)
/// Update the volume estimators J and nuJ
/// This is done in another routine than move, as we sometimes move dummy
/// packets which do not contribute to the radiation field.
{
  const int cellindex = pkt_ptr->where;
  const int modelgridindex = cell[cellindex].modelgridindex;

  /// Update only non-empty cells
  if (modelgridindex != MMODELGRID)
  {
    const double helper = distance * pkt_ptr->e_cmf;
    #ifdef _OPENMP
      #pragma omp atomic
    #endif
    J[modelgridindex] += helper;

    #ifndef FORCE_LTE
      const double nu = pkt_ptr->nu_cmf;
      const double helper2 = helper / nu;
      //double bf = exp(-HOVERKB*nu/cell[modelgridindex].T_e);
      #ifdef _OPENMP
        #pragma omp atomic
      #endif
      nuJ[modelgridindex] += helper * nu;
      radfield_update_estimators(modelgridindex,distance,pkt_ptr->e_cmf,nu);

      ///ffheatingestimator does not depend on ion and element, so an array with gridsize is enough.
      ///quick and dirty solution: store info in element=ion=0, and leave the others untouched (i.e. zero)
      #ifdef _OPENMP
        #pragma omp atomic
      #endif
      ffheatingestimator[modelgridindex] += helper * kappa_rpkt_cont[tid].ffheating;
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
          if (get_abundance(modelgridindex,element) > 0)
          {
            #ifdef _OPENMP
              #pragma omp atomic
            #endif
            gammaestimator[modelgridindex*nelements*maxion+element*maxion+ion] += phixslist[tid].groundcont[i].gamma_contr * helper2;

            #ifdef DEBUG_ON
            if (!isfinite(gammaestimator[modelgridindex*nelements*maxion+element*maxion+ion]))
            {
              printout("[fatal] update_estimators: gamma estimator becomes non finite: level %d, gamma_contr %g, helper2 %g\n",i,phixslist[tid].groundcont[i].gamma_contr,helper2);
              abort();
            }
            #endif
            #ifdef _OPENMP
              #pragma omp atomic
            #endif
            bfheatingestimator[modelgridindex*nelements*maxion+element*maxion+ion] += phixslist[tid].groundcont[i].gamma_contr * helper * (1. - nu_edge/nu);
            //bfheatingestimator[modelgridindex*nelements*maxion+element*maxion+ion] += phixslist[tid].groundcont[i].bfheating_contr * helper * (1/nu_edge - 1/nu);
          }
        }
        else break;
      }

      #ifdef DEBUG_ON
        if (!isfinite(nuJ[modelgridindex]))
        {
          printout("[fatal] update_estimators: estimator becomes non finite: helper %g, nu_cmf %g ... abort\n",helper,pkt_ptr->nu_cmf);
          abort();
        }
      #endif
    #endif

    ///Heating estimators. These are only applicable for pure H. Other elements
    ///need advanced treatment in thermalbalance calculation.
    //cell[pkt_ptr->where].heating_ff += helper * kappa_rpkt_cont[tid].ffheating;
    //cell[pkt_ptr->where].heating_bf += helper * kappa_rpkt_cont[tid].bfheating;


    #ifdef DEBUG_ON
      if (!isfinite(J[modelgridindex]))
      {
        printout("[fatal] update_estimators: estimator becomes non finite: helper %g, nu_cmf %g ... abort\n",helper,pkt_ptr->nu_cmf);
        abort();
      }
    #endif
  }

}


void move_pkt(PKT *restrict pkt_ptr, double distance, double time)
/// Subroutine to move a packet along a straight line (specified by currect
/// dir vector). The distance moved is in the rest frame. Time must be the
/// time at the end of distance travelled.
{
  /// First update pos.
  if (distance < 0)
  {
    printout("Trying to move -v distance. Abort.\n");
    abort();
  }

  //printout("Move distance %g\n", distance);
  pkt_ptr->pos[0] += (pkt_ptr->dir[0] * distance);
  pkt_ptr->pos[1] += (pkt_ptr->dir[1] * distance);
  pkt_ptr->pos[2] += (pkt_ptr->dir[2] * distance);

  /// During motion, rest frame energy and frequency are conserved.
  /// But need to update the co-moving ones.
  double vel_vec[3];
  get_velocity(pkt_ptr->pos, vel_vec, time);
  double dopplerfac = doppler(pkt_ptr->dir, vel_vec);
  pkt_ptr->nu_cmf = pkt_ptr->nu_rf * dopplerfac;
  pkt_ptr->e_cmf = pkt_ptr->e_rf * dopplerfac;

  /*
  if (pkt_ptr->e_rf * pkt_ptr->nu_cmf /pkt_ptr->nu_rf > 1e46)
    {
      printout("here2 %g %g \n", pkt_ptr->e_rf, pkt_ptr->nu_cmf /pkt_ptr->nu_rf);
    }
  */
}
