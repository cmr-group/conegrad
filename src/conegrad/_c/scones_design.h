/* cones.h
 *
 * Header file for cones generating routines
 *
 * By: Jeff Stainsby
 * Last Rev: 27Nov2006
 */

/*#ifndef CONES_H
  #define CONES_H */

/* Following needed in EPIC compile */
/* extern float GAM; */
/* End */

/* Following needed in standalone compile
float GAM;
typedef struct {
  int xrt;
  int yrt;
  int zrt;
  float tx;
  float ty;
  float tz;
} LOG_GRAD;
#define max(a,b) ((a) < (b) ? (b) : (a))
#define min(a,b) ((a) < (b) ? (a) : (b))
#include "support_decl.h"

*/
/* End  */

/* Output file locations */

#ifdef SIM
#define DIR_OUT "."
#else
#define DIR_OUT "/usr/g/mrraw"
#endif
#define FILE_GXCONE "kacq_cone_gx"
#define FILE_GYCONE "kacq_cone_gy"
#define FILE_GZCONE "kacq_cone_gz"
#define FILE_KTRAJ "kacq_cone_ktraj"
#define FILE_DENS "kacq_cone_dens"
#define FILE_NPC "kacq_cone_npc"


/*#ifdef SIM
#define GXCONE_FILE "gxcone"
#define GYCONE_FILE "gycone"
#define GZCONE_FILE "gzcone"
#else
#define GXCONE_FILE "/usr/g/psddata/gxcone"
#define GYCONE_FILE "/usr/g/psddata/gycone"
#define GZCONE_FILE "/usr/g/psddata/gzcone"
#endif*/

#define MAX_THETAS 1000
#define MAX_THETAS_REDUCED 500
#define MAX_CONES 100

#ifndef GAMMA_PROTON
#define GAMMA_PROTON 4257.59
#endif
#ifndef GAMMA_FLUORINE
#define GAMMA_FLUORINE 4006.04
#endif
#ifndef GAMMA_HELIUM
#define GAMMA_HELIUM 3243.38      /* negative */
#endif
#ifndef GAMMA_PHOSPHORUS
#define GAMMA_PHOSPHORUS 1723.51
#endif
#ifndef GAMMA_LITHIUM
#define GAMMA_LITHIUM 1654.76
#endif
#ifndef GAMMA_BORON
#define GAMMA_BORON 1365.98
#endif
#ifndef GAMMA_XENON
#define GAMMA_XENON 1177.60       /* negative */
#endif
#ifndef GAMMA_SODIUM
#define GAMMA_SODIUM 1126.09
#endif
#ifndef GAMMA_CARBON
#define GAMMA_CARBON 1070.63
#endif
#ifndef GAMMA_SILICON
#define GAMMA_SILICON 845.908     /* negative */
#endif
#ifndef GAMMA_DEUTERIUM
#define GAMMA_DEUTERIUM 653.616
#endif
#ifndef GAMMA_OXYGEN
#define GAMMA_OXYGEN 577.357      /* negative */
#endif
#ifndef GAMMA_NITROGEN
#define GAMMA_NITROGEN 431.576    /* negative */
#endif


#define KTRAJ_IGNORE -99.0 /* k-space trajectory value indicating the point should be ignored */
#define DENS_IGNORE 0.0 /* density value for points to ignore */
#define MAX_LOOPS 10000 /* Max number of iterations in binary search before bail with current params */

/* Prototypes for float comparisons (DV27) */
double wcfabs(double a);
double wcfmin(double a, double b);
double wcfmax(double a, double b);
int wcfEqual(const float x, const float y);

/* 3D Cones function templates based on code from Paul Gurney's wc.c
 and additional functions based on matlab functions */
int rootof(float a, float b, float c, float *s1, float *s2);
int circcirc(float d, float e, float k, float s, float type, float *x, float *y);
int intlinecirc(float a, float b, float d, float e, float r, float *s1, float *s2);
float anglebetween(float ar, float ai, float br, float bi);
float anglebetween_new(float ar, float ai, float br, float bi);
int rto(float Gc, float kdes, float Smax, float Gmax, float gam, float GSP, float *g, float *k);
int rtz(float gcx, float gcy, float gcz, float kcx, float kcy, float kcz, float theta, int maxai, float Gmax_xy, float Gmax_z, float Smax_xy, float Smax_z, float GSP, float * pgx, float *pgy, float *pgz, float *pkx, float *pky, float *pkz, float FOV1, float FOV2, float NINT, float kmax, float dcf);
float calcn(float k, float sintheta, float FOV1, float FOV2, float NINT, float kmax, float dcf);
float calcdn(float k, float sintheta, float FOV1, float FOV2, float NINT, float kmax, float dcf);
int whirlcone(  float theta, float FOV1, float FOV2,	float FOV3, float FOV4, float dcf, float kmax, float NINT, int numpt, float GSP, float SMAX, float GMAX, float Smax_xy, float Smax_z, float Gmax_xy, float Gmax_z, float dens_min, float ***pg, float ***pk, float *pai, float *pn, int acq_mode); /*MCARL: added acq_mode, SMAX, and GMAX*/
int whirlcone_vd(  float theta, float *FOV1, float *FOV2,	float *FOV3, float *FOV4, float dcf, float kmax, float NINT, int numpt, float GSP, float SMAX, float GMAX, float Smax_xy, float Smax_z, float Gmax_xy, float Gmax_z, float dens_min, float ***pg, float ***pk, float *pai, float *pn, int acq_mode); /*MCARL: added acq_mode, SMAX, and GMAX*/
float *getthetas(float *RES, float *FOV, int *ntheta, float *dtheta_ic);
void getthetas2(float *theta, float RES[2], float FOV[2], int *ntheta, float *dtheta_ic);
int gencone(float ***g, int *len, int *totlen, float *RES, float *FOV, float NINT, float THETA[2], int MAXLEN, float GSP, float SMAX, float maxSRewind, float GMAX, float sysGMAX, int OVERSAMPLE, float DCF, float MINDENS,int nramp, int rewind_flag, int acq_mode); /*MCARL: acq_mode*/
int gencone_vd(float ***g, int *len, int *totlen, float *RES, float *FOV, float NINT, float THETA[2], int MAXLEN, float GSP, float SMAX, float maxSRewind, float GMAX, float sysGMAX, int OVERSAMPLE, float DCF, float MINDENS,int nramp, int rewind_flag, int acq_mode); /*MCARL: acq_mode*/
int findcone(float ***g, float *nint, int *leng, int *totleng, float *RES, float *FOV,int GRAD_POINTS, float THETA[2],float PRECISION,float GSP,int OVERSAMPLE,float SMAX,float maxSRewind,float GMAX,float sysGMAX,float DCF,float MINDENS,int nramp, int rewind_flag, int acq_mode); /*MCARL: acq_mode*/
int findcone_vd(float ***g, float *nint, int *leng, int *totleng, float *RES, float *FOV,int GRAD_POINTS, float THETA[2],float PRECISION,float GSP,int OVERSAMPLE,float SMAX,float maxSRewind,float GMAX,float sysGMAX,float DCF,float MINDENS,int nramp, int rewind_flag, int acq_mode); /*MCARL: acq_mode*/
float intlineellipse(float a, float b, float phi);
int mrewind(float *gx, float *gy, float *k0, float *g0, float gmax, float smax, float T);
int genktraj(float **ktmp, float *dens, float *nint, float **grad, int grad_points, int GRAD_POINTS, float xyscale, float zscale, float FOVrad, float FOVcirc, float truetheta, float dtheta_ic, float GSP,float *maxgr);
/* Python wrapper extension: added trailing `int *traj_length` (per-cone
   trajectory length in gradient samples, before rewinder) so callers can
   slice the returned waveforms at the true readout end. */
int conegrad(int *nintpc, float *rspthetas, float *RES, float *FOV, int *gxcones_host, int *gycones_host, int *gzcones_host, int NUMCONES, int GRAD_POINTS, float GSP, int OPXRES, float TSP, float PRECISION, float DCF, int OVERSAMPLE, float SMAX, float maxSRewind, float GMAX, float sysGMAX, float SysMaxRewindG, int output_grad, int ktraj_flag, int ktraj_out_flag, int endian_flag, int rewind_flag, float MINDENS, float *maxgr, float *snr_eff, int *ngradact, float t_fracx, float t_fracy, float t_fracz, int acq_mode, float rot_flag, int SymFlag, float SlabKz, int rhkacq_uid, int Cones_Plot_Flag, int *traj_length); /*MCARL: Added t_frac, acq_mode, rotflag, rhkacq_uid, Cones_Plot_Flag*/
int conegrad_vd(int *nintpc, float *rspthetas, float *RES, float *FOV, int *gxcones_host, int *gycones_host, int *gzcones_host, int NUMCONES, int GRAD_POINTS, float GSP, int OPXRES, float TSP, float PRECISION, float DCF, int OVERSAMPLE, float SMAX, float maxSRewind, float GMAX, float sysGMAX, float SysMaxRewindG, int output_grad, int ktraj_flag, int ktraj_out_flag, int endian_flag, int rewind_flag, float MINDENS, float *maxgr, float *snr_eff, int *ngradact, float t_fracx, float t_fracy, float t_fracz, int acq_mode, float rot_flag, int SymFlag, float SlabKz, int rhkacq_uid, int Cones_Plot_Flag, int *traj_length);
int amppwrewind(float targarea, float targamp, float st_amp, float en_amp, int Trise, int mindur, float *Amp, int *Attack, int *Pw, int *Decay);
double ceil2(double x);

/* JYCHENG --> */
int permuteCircGoldenCones(int[], int, int, unsigned short);
/* <-- JYCHENG */

/* MOMALAVE --> */
int permutePhyllo(int[], int[], int, int*, int[], int[], int);
/* <-- MOMALAVE */

/* #endif */
