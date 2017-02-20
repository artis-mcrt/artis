// --------------------------------------------------------------------------------
// ---------------------------  VIRTUAL PACKETS -----------------------------------
// --------------------------------------------------------------------------------


/* FREQUENCY */
/* dlognu = (log(numax) - log(numin)) / VMNUBINS ~ 3.9e-4 (10'000 over 1e14-5e15 Hz) */ 
#define numin_vspec (CLIGHT / 10000 * 1e8)
#define numax_vspec (CLIGHT / 3500 * 1e8)
#define VMNUBINS 2500  

/* TIME */
/* dlogt = (log(tmin) - log(tmax)) / VMTBINS ~ 3.69e-2 (111 over 2-120 d) */  
#define tmin_vspec (10 * DAY)
#define tmax_vspec (30 * DAY)
#define VMTBINS 30

/* Number of spectra for each observer (total + elements switched off) */
#define MSPECTRA 12
/* Number of observers */
#define MOBS 5
/* Total number of spectra */
#define MTOT (MSPECTRA * MOBS)
/* Total number of frequency ranges */
#define MRANGE 1

struct vspecpol
{
  double flux[VMNUBINS];
  float lower_time;
  float delta_t;
} vstokes_i[VMTBINS][MTOT], vstokes_q[VMTBINS][MTOT], vstokes_u[VMTBINS][MTOT];

float lower_freq_vspec[VMNUBINS];
float delta_freq_vspec[VMNUBINS];

// --------- INPUT PARAMETERS -----------

int Nobs,Nspectra;
double nz_obs_vpkt[MOBS] ;
double phiobs[MOBS] ;
double tmin_vspec_input,tmax_vspec_input;
double Nrange;
double lmin_vspec_input[MRANGE], lmax_vspec_input[MRANGE];
double numin_vspec_input[MRANGE], numax_vspec_input[MRANGE];
double cell_is_optically_thick_vpkt;
double tau_max_vpkt;
double exclude[MSPECTRA],tau_vpkt[MSPECTRA];

// --------- VPKT GRID -----------

#define MRANGE_GRID 5
#define NY_VGRID 50
#define NZ_VGRID 50

struct vgrid
{
  double flux[MRANGE_GRID][MOBS];
  double yvel[MRANGE_GRID][MOBS];
  double zvel[MRANGE_GRID][MOBS];
} vgrid_i[NY_VGRID][NZ_VGRID], vgrid_q[NY_VGRID][NZ_VGRID], vgrid_u[NY_VGRID][NZ_VGRID];

double Nrange_grid, tmin_grid, tmax_grid, nu_grid_min[MRANGE_GRID], nu_grid_max[MRANGE_GRID] ;
int vgrid_flag;
double dlogt_vspec,dlognu_vspec;

// --------- FLAGS -----------

int realtype ;
/* number of virtual packets in a given timestep */
int nvpkt ;
/* number of escaped virtual packet in a given timestep (with tau < tau_max) */
int nvpkt_esc1 ; /* electron scattering event */
int nvpkt_esc2 ; /* kpkt deactivation */
int nvpkt_esc3 ; /* macroatom deactivation */

