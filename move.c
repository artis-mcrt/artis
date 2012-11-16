#ifndef MOVE_C
#define MOVE_C
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include "sn3d.h"
#endif


///****************************************************************************
void
update_estimators (PKT * pkt_ptr, double distance)
/// Update the volume estimators J and nuJ
/// This is done in another routine than move, as we sometimes move dummy
/// packets which do not contribute to the radiation field.
{
  double get_abundance (int modelgridindex, int element);
  int element, ion, level, i;
  double nu_edge;

  int modelgridindex = cell[pkt_ptr->where].modelgridindex;

  /// Update only non-empty cells
  if (modelgridindex != MMODELGRID)
    {
      double helper = distance * pkt_ptr->e_cmf;
      double nu = pkt_ptr->nu_cmf;

#ifdef _OPENMP
#pragma omp atomic
#endif
      J[modelgridindex] += helper;

#ifndef FORCE_LTE
      double helper2 = helper / nu;
      //double bf = exp(-HOVERKB*nu/cell[modelgridindex].T_e);
#ifdef _OPENMP
#pragma omp atomic
#endif
      nuJ[modelgridindex] += helper * nu;

      ///ffheatingestimator does not depend on ion and element, so an array with gridsize is enough.
      ///quick and dirty solution: store info in element=ion=0, and leave the others untouched (i.e. zero)
#ifdef _OPENMP
#pragma omp atomic
#endif
      ffheatingestimator[modelgridindex] +=
	helper * kappa_rpkt_cont[tid].ffheating;
      for (i = 0; i < nbfcontinua_ground; i++)
	{
	  nu_edge = phixslist[tid].groundcont[i].nu_edge;
	  if (nu > nu_edge)
	    {
	      element = phixslist[tid].groundcont[i].element;
	      ion = phixslist[tid].groundcont[i].ion;
	      /// Cells with zero abundance for a specific element have zero contribution 
	      /// (set in calculate_kappa_rpkt_cont and therefore do not contribute to
	      /// the estimators
	      if (get_abundance (modelgridindex, element) > 0)
		{
#ifdef _OPENMP
#pragma omp atomic
#endif
		  gammaestimator[modelgridindex * nelements * maxion +
				 element * maxion + ion] +=
		    phixslist[tid].groundcont[i].gamma_contr * helper2;
#ifdef _OPENMP
#pragma omp atomic
#endif
		  bfheatingestimator[modelgridindex * nelements * maxion +
				     element * maxion + ion] +=
		    phixslist[tid].groundcont[i].gamma_contr * helper * (1. -
									 nu_edge
									 /
									 nu);
		  //bfheatingestimator[modelgridindex*nelements*maxion+element*maxion+ion] += phixslist[tid].groundcont[i].bfheating_contr * helper * (1/nu_edge - 1/nu);
		}
	    }
	  else
	    break;
	}

#ifdef DEBUG_ON
      if (!finite (nuJ[modelgridindex]))
	{
	  printout
	    ("[fatal] update_estimators: estimator becomes non finite: helper %g, nu_cmf %g ... abort\n",
	     helper, pkt_ptr->nu_cmf);
	  abort ();
	}
#endif
#endif

      ///Heating estimators. These are only applicable for pure H. Other elements
      ///need advanced treatment in thermalbalance calculation.
      //cell[pkt_ptr->where].heating_ff += helper * kappa_rpkt_cont[tid].ffheating;
      //cell[pkt_ptr->where].heating_bf += helper * kappa_rpkt_cont[tid].bfheating;


#ifdef DEBUG_ON
      if (!finite (J[modelgridindex]))
	{
	  printout
	    ("[fatal] update_estimators: estimator becomes non finite: helper %g, nu_cmf %g ... abort\n",
	     helper, pkt_ptr->nu_cmf);
	  abort ();
	}
#endif
    }

}


///****************************************************************************
int
move_pkt (PKT * pkt_ptr, double distance, double time)
/// Subroutine to move a packet along a straight line (specified by currect
/// dir vector). The distance moved is in the rest frame. Time must be the
/// time at the end of distance travelled.
{
  double doppler ();
  int get_velocity ();

  double vel_vec[3];

  /// First update pos.
  if (distance < 0)
    {
      printout ("Trying to move -v distance. Abort.\n");
      abort ();
    }

#ifdef PHOTO
//This part checks if the packet crashes into the photosphere. If this happens the packet will be destroyed.
  double pos_vec[3];
  double rCenterSquared;

  pos_vec[0] = pkt_ptr->pos[0] + (pkt_ptr->dir[0] * distance);
  pos_vec[1] = pkt_ptr->pos[1] + (pkt_ptr->dir[1] * distance);
  pos_vec[2] = pkt_ptr->pos[2] + (pkt_ptr->dir[2] * distance);

  rCenterSquared =
    gsl_pow_2 (pos_vec[0]) + gsl_pow_2 (pos_vec[1]) + gsl_pow_2 (pos_vec[2]);

  if (rCenterSquared <= gsl_pow_2 (rPhotosphere))
    {
      pkt_ptr->type = TYPE_INACTIVE;
      printout
	("[info] packet %i is crashed into the photosphere and will be set to inactive.");
    }
  pkt_ptr->pos[0] = pos_vec[0];
  pkt_ptr->pos[1] = pos_vec[1];
  pkt_ptr->pos[2] = pos_vec[2];

#else
  //printout("Move distance %g\n", distance);
  pkt_ptr->pos[0] = pkt_ptr->pos[0] + (pkt_ptr->dir[0] * distance);
  pkt_ptr->pos[1] = pkt_ptr->pos[1] + (pkt_ptr->dir[1] * distance);
  pkt_ptr->pos[2] = pkt_ptr->pos[2] + (pkt_ptr->dir[2] * distance);
#endif



  /// During motion, rest frame energy and frequency are conserved.
  /// But need to update the co-moving ones.
  get_velocity (pkt_ptr->pos, vel_vec, time);
  pkt_ptr->nu_cmf = pkt_ptr->nu_rf * doppler (pkt_ptr->dir, vel_vec);
  pkt_ptr->e_cmf = pkt_ptr->e_rf * pkt_ptr->nu_cmf / pkt_ptr->nu_rf;

  /*
     if (pkt_ptr->e_rf * pkt_ptr->nu_cmf /pkt_ptr->nu_rf > 1e46)
     {
     printout("here2 %g %g \n", pkt_ptr->e_rf, pkt_ptr->nu_cmf /pkt_ptr->nu_rf);
     }
   */

  return 0;
}
