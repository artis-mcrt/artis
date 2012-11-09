//This file contain the packet initialisation for a photosphere run.

#ifdef PHOTO
#ifndef PACKET_INIT_PHOTOSPHERE_C
#define PACKET_INIT_PHOTOSPHERE_C

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_const_cgs.h>

#include "sn3d.h"

#endif
int packet_init_photosphere(int middle_iteration, int my_rank)
{
  int setup_packets_photosphere(int pktnumberoffset);
  void write_packets(FILE *packets_file); //We use the write_packets function from the  packet_init file.
  FILE *packets_file;
  int pktnumberoffset;
  char filename[100];               /// this must be long enough to hold "packetsxx.tmp" where xx is the number of "middle" iterations

  if (!continue_simulation)
  {
    pktnumberoffset = middle_iteration*npkts;
    setup_packets_photosphere(pktnumberoffset);
    sprintf(filename,"packets%d_%d_odd.tmp",middle_iteration,my_rank);
    if ((packets_file = fopen(filename, "wb")) == NULL)
    {
      printf("[fatal]: packet_init: Cannot open packets file\n");
      exit(0);
    }
    fwrite(&pkt[0], sizeof(PKT), npkts, packets_file);
    //write_packets(packets_file);
    fclose(packets_file);

  }

  return 0;
}


int setup_packets_photosphere (int pktnumberoffset)
/// Subroutine that initialises the packets if we start a new simulation.
{
  //Temporary variables. Only for debugging 
  double const LUMINOSITY =1000;
  double const rPhotosphere=1000;
  double const tPhotosphere=3600;
  double const t_current = tmin; // valid in steady stat.
  //int place_pellet(struct grid *grid_ptr, double e0, int m, int n, int pktnumberoffset); // I think we dont need this.

  double rand_nu(double _tphotosphere);
  double doppler_shift(const size_t n, double *directionVec, double *velocityVec);
  double get_cell_velocity_3(double *positionVecPrt, double *velocityVecPrt, double time);

  double etot;
  double e0;
  int n;
  double zrand;
  //CELL *grid_ptr;
  //double fni(CELL *grid_ptr),f48cr(CELL *grid_ptr),f52fe(CELL *grid_ptr);
  //double vol_init(CELL *grid_ptr);
  //float norm;
  //int mabove, mbelow;
  //int packet_reset;

  double thetaEmittDirectionLF;
  double phiEmittDirectionLF;
  double xnEmittDirection,ynEmittDirection,znEmittDirection;

  double thetaEmittPosition;
  double phiEmittPosition;
  double xnEmittPosition,ynEmittPosition,znEmittPosition;
  int xCellIndex,yCellIndex,zCellIndex;
  int cellId;

  double currentVelocity[3];

  /// The total number of pellets that we want to start with is just
  /// npkts. The total energy of the pellets is given by Luminosity L.

  /// So energy per pellet is
  e0 = LUMINOSITY / npkts / n_out_it / n_middle_it;
  printout("e0 %g\n", e0);

  /* Now we place the pellets in the ejecta. */

  /* Need to get a normalisation factor. */
  cont[ngrid] = 0.0; //I dont know if this is really necessary.

  if (npkts > MPKTS)
  {
    printout("Too many packets. Abort.\n");
    exit(EXIT_FAILURE);
  }


  //add here the photosphere code and

  for (n = 0; n < npkts; n++)
  {
    /// Get random number.
    zrand = gsl_rng_uniform(rng);

    //Position of the photon on the surface of the Photosphere 
    thetaEmittPosition = acos(1-2*gsl_rng_uniform(rng));
    phiEmittPosition = M_PI * 2. * gsl_rng_uniform(rng);
    
    xnEmittPosition = sin(thetaEmittPosition) * cos(phiEmittPosition);
    ynEmittPosition = sin(thetaEmittPosition) * cos(thetaEmittPosition);
    znEmittPosition = cos(thetaEmittPosition);


    //try to get the cell index
    //######WE HAVE TO DEFINE THE VARIABLE rPhotosphere and tPhotosphere. 
    //The current definition in this function is only temporary and should be moved to the input file.
    
    xCellIndex = (int) ((rPhotosphere * xnEmittPosition) / wid_init - 0.5 * nxgrid);
    yCellIndex = (int) ((rPhotosphere * ynEmittPosition) / wid_init - 0.5 * nygrid);
    zCellIndex = (int) ((rPhotosphere * znEmittPosition) / wid_init - 0.5 * nzgrid);

    //try to get the cell id

    cellId = zCellIndex * nygrid * nxgrid + yCellIndex * nxgrid + xCellIndex;


    


    //direction of the photon in the local frame(LF) i.e. the surface of the photosphere is assumed to be locally plain 
    thetaEmittDirectionLF = acos(1-2*gsl_rng_uniform(rng));
    phiEmittDirectionLF = M_PI * gsl_rng_uniform(rng);

    //computing the direction of the photon in the  global frame
    xnEmittDirection =   xnEmittPosition * cos(phiEmittDirectionLF) 
                       + znEmittPosition * cos(phiEmittDirectionLF);

    ynEmittDirection = - xnEmittPosition * cos(thetaEmittDirectionLF) * sin(phiEmittDirectionLF) 
                       + ynEmittPosition * cos(thetaEmittDirectionLF) 
                       - znEmittPosition * sin(thetaEmittPosition)    * sin(phiEmittPosition);

    znEmittDirection = - xnEmittPosition * cos(thetaEmittPosition)    * sin(phiEmittPosition) 
                       + ynEmittPosition * sin(thetaEmittDirectionLF) 
                       + znEmittPosition * cos(thetaEmittDirectionLF) * sin(phiEmittDirectionLF);

    

    //Transformation the photon direction from the local frame(LF) to the global frame(GF)


    //Store the properties of the packet in the corresponding structure 

    pkt[n].where = cellId;
    pkt[n].number = n + pktnumberoffset;  ///record the packets number for debugging


    /// Some fraction of the packets we reasigned because they were not going
    /// to activate in the time of interest so need to renormalise energies
    /// to account for this.


    //Position of the packet.
    pkt[n].pos[0] = xnEmittPosition * rPhotosphere;
    pkt[n].pos[1] = ynEmittPosition * rPhotosphere;
    pkt[n].pos[2] = znEmittPosition * rPhotosphere;

    //compute the current velocity.
    //I hope that t_current is a global variable which contains the current time.
    get_cell_velocity_3(pkt[n].pos, currentVelocity, t_current);
    

    //Direction of the packet.
    pkt[n].dir[0] = xnEmittDirection;
    pkt[n].dir[1] = ynEmittDirection;
    pkt[n].dir[2] = znEmittDirection;

    /// Now assign the energy and the temperature  to the pellet.
    pkt[n].e_cmf = e0; //the energy of the packet in the comoving frame 
    pkt[n].nu_cmf = rand_nu(tPhotosphere); //draw a random nu and assign it to the packet. 
    pkt[n].nu_rf = pkt[n].nu_cmf /doppler_shift(3,pkt[n].dir, currentVelocity); //compute the rest frame frequency via the doppler shift 
    pkt[n].e_rf = pkt[n].e_cmf * pkt[n].nu_cmf / pkt[n].nu_rf;  //computing the rest frame energy by using the ration between rf_nu and cmf_nu.

    pkt[n].interactions = 0;

    //Now we set the packet type to rpacket
    pkt[n].type = TYPE_RPKT;
  }

  return 0;
}


double
rand_nu(double _tphotosphere)
{
#define RAND_NU_PLANCK(v,T) ((c1 * gsl_pow_3((v)))/(gsl_sf_exp(c2*((v)/(T))) - 1 ))
  //This function draws a random nu corresponding to the temperature.
  const double c1 = 2 * GSL_CONST_CGS_PLANCKS_CONSTANT_H/(gsl_pow_2(GSL_CONST_CGS_SPEED_OF_LIGHT));
  const double c2 = GSL_CONST_CGS_PLANCKS_CONSTANT_H/ GSL_CONST_CGS_BOLTZMANN;
  const double nuPeak = 5.8789254e10 * _tphotosphere;// The constant here is taken from Wien's displacement law.
  const double bPeak = RAND_NU_PLANCK(nuPeak,_tphotosphere);


  double bRand, nuRand, bFromNu;
  double nuInterval = nu_max_r - nu_min_r;
  //check if nu range is valid 5.88e10 is derived from wien.
  if (((nuPeak) > nu_max_r)&&((nuPeak)< nu_min_r)) {
    printf("The Planck function peaks(%g < %g < %g) outside the configured frequency range. Abort.\n",nu_min_r,nuPeak,nu_max_r);
    exit(EXIT_FAILURE);
  }

  while (1)
  {
    nuRand = gsl_rng_uniform(rng) * nuInterval + nu_min_r;
    bRand = bPeak * gsl_rng_uniform(rng);
    bFromNu = RAND_NU_PLANCK(nuRand,_tphotosphere);
    if (bRand < bFromNu) break;
  }
  return nuRand;
#undef RAND_NU_PLANCK
}

double
doppler_shift(const size_t n, double *directionVec, double *velocityVec)
{
//This function computes the Doppler shift for a given velocity(wiki).

  double dot_product(double *a, double *b, const size_t n);
  double gammaCoef;

  gammaCoef = 1./(sqrt(1- (dot_product(velocityVec,velocityVec,n)/GSL_CONST_CGS_SPEED_OF_LIGHT)));

  return gammaCoef * (1 - (dot_product(directionVec, velocityVec,n)/gsl_pow_2(GSL_CONST_CGS_SPEED_OF_LIGHT)));
}

double
dot_product(double *a, double *b, const size_t n)
{
//This function computes the dot product of the two given vectors. a
//The function can be improved by using c99 feature instead of ANSI C.
        double sum = 0;
        size_t i;
 
        for (i = 0; i < n; i++) {
                sum += a[i] * b[i];
        }
        return sum;
}

double
get_cell_velocity_3(double *positionVecPrt, double *velocityVecPrt, double time)
{
  //This function computes the velocity of a particle/object at a given position assuming homologies expansion(roepke2005).

  velocityVecPrt[0] = positionVecPrt[0] / time;
  velocityVecPrt[1] = positionVecPrt[1] / time;
  velocityVecPrt[2] = positionVecPrt[2] / time;

  return 0;
}
#endif
