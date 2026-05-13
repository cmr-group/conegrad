/*
 *   cones_design.c
 *
 *    By: Jeff Stainsby
 *    Last Rev: 23Nov2006
 */

/*  Functions to implement gradient and k-space trajectory design for
    3D Cones trajectory. Based on matlab code from Paul Gurney at
    Stanford. For more details see: Gurney et al., MRM, 55, 2006
*/

/* Revision history:

   epoch - 23Nov2006
      Initial implementation, converting Paul's matlab functions into
      C-code. Initial working version and successful acquisition of
      phantom and in vivo sodium images.

*/

/* Global include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*#include "cnv_endian_API.h"*/
#include "scones_design.h"

/* #include <dmalloc.h>   */

#ifndef STANDALONE
#include "pulsegen.h" /* For defn of WEOS_BIT */
#endif

/* Think this is vxWorks vs linux simulator */
#if defined(HW_IO) || (defined(IPG) && defined(PSD_HW))
#ifndef IS26X
#include <stdioLib.h> /* hardware */
#endif
#else              /* !(HW_IO || (IPG && PSD_HW)) */
#include <stdio.h> /* simulation */
#endif             /* HW_IO || (IPG && PSD_HW) */

/* Macros for min/max. Note these are previously defined on the
   target side of the EPIC compile but shouldn't really do any harm. */
#define max(a, b) ((a) < (b) ? (b) : (a))
#define min(a, b) ((a) < (b) ? (a) : (b))

#define PI 3.14159265358979
#define MAX_PG_WAMP 32766 /* This is defined in epic.h ... */
#define EQ_TOL 1.0e-09;

/* For amppwrewind fn */
#define MAX_ITER 200
#ifndef GRAD_UPDATE_TIME
#define GRAD_UPDATE_TIME 4
#endif

#ifndef RUP_GRD
#define RUP_GRD(A)                                                             \
  (((A) % GRAD_UPDATE_TIME)                                                    \
       ? (int)((A) + GRAD_UPDATE_TIME) & ~(GRAD_UPDATE_TIME - 1)               \
       : (A))
#endif

/* To allow these to be used on host and tgt side need to call the
   appropriate memory management functions */
#ifdef HOST_TGT
#define AllocMem malloc
#define FreeMem free
#else
#define AllocMem AllocNode
#define FreeMem FreeNode
#endif

/* Override these defines in standalone version */
#ifdef STANDALONE
#undef AllocMem
#undef FreeMem
#define AllocMem malloc
#define FreeMem free
#ifndef SUCCESS
#define SUCCESS 0
#endif
#ifndef FAILURE
#define FAILURE 1
#endif
#endif

/* Conegrad Python wrapper: route all printf through a runtime verbose flag.
   fprintf (used for plot/data file output) is left alone. */
int conegrad_verbose = 0;
#define printf(...) do { if (conegrad_verbose) (void)fprintf(stdout, __VA_ARGS__); } while (0)

/* Conegrad Python wrapper: upper bound on per-cone waveform sample count used
   to size the static scratch arrays in genktraj / genktraj_vd. The original
   code hard-coded 8192, which silently overflowed when grad_points exceeded
   that. Bumped to 65536 so the wrapper handles up to ~64k-sample readouts;
   the arrays are declared `static` so the larger size lives in BSS instead
   of inflating the per-call stack frame. */
#define CG_MAX_WAVEFORM_PTS 65536

double wcfabs(double a) 
{ return (a < 0) ? -a : a; }
double wcfmin(double a, double b) 
{ return (a < b) ? a : b; }
double wcfmax(double a, double b) 
{ return (a > b) ? a : b; }

/* From: realtimecollisiondetection.net/pubs/Tolerances/
   This was good too: floating-point-gui.de/errors/comparison,
   but it uses floating comparison in the algorithm which defeats the
   purpose. However, here are some good tests: 
   - floating-point-gui.de/errors/NearlyEqualsTest.java
 */
int wcfEqual(const float x, const float y) 
{
  /* In case some day I want to pass these and make them different */
  float abstol = EQ_TOL;
  float reltol = EQ_TOL;
  float maxXY = wcfmax(wcfabs(x), wcfabs(y));
  float theTol = wcfmax(abstol, reltol*maxXY);

  return (wcfabs(x - y) <= theTol) ? 1 : 0;
}

/* Function:  rootof
 *
 * To find the root of the quadratic function ax^2 + bx + c = 0.
 * Returns the roots in s1 and s2
 */
/* CHECKED */
int rootof(float a, float b, float c, float *s1, float *s2) {
  float descrim;
  if (wcfEqual(a, 0)) {
    if (wcfEqual(b, 0)) {
      *s1 = 0;
      *s2 = 0;
      return 0;
    } else {
      *s1 = -c / b;
      *s2 = -c / b;
      return 1;
    }
  } else {
    descrim = b * b - 4 * a * c;
    if (descrim < 0) {
      *s1 = 0;
      *s2 = 0;
      return 0;
    } else {
      *s1 = (-b + sqrt(descrim)) / 2 / a;
      *s2 = (-b - sqrt(descrim)) / 2 / a;
      return 2;
    } /* End if descrim<0 */
  }   /* End if a ==0, else */
}

/* Function:  circcirc
 *
 * ??? Find the intersection of two circles?
 */
/* CHECKED */
int circcirc(float d, float e, float k, float s, float type, float *x,
             float *y) {

  float s1;
  float s2;
  float t;
  float ade;
  float tx1, tx2, ty1, ty2;
  float a, b, c;
  int numsols;
  if (wcfEqual(e, 0))
    ade = 100;
  else
    ade = fabs((float)(d / e));

  if ((wcfEqual(d, 0)) && (wcfEqual(e, 0))) {
    *x = 0;
    *y = 0;
    return 0;
  } else if (ade < 1) {
    t = e;
    e = d;
    d = t;
    a = 4 * d * d + 4 * e * e;
    b = (-4 * e * e * e - 4 * k * k * e + 4 * e * s * s - 4 * d * d * e);
    c = d * d * d * d + k * k * k * k - 2 * e * e * s * s + 2 * d * d * e * e -
        2 * d * d * s * s + s * s * s * s - 2 * k * k * s * s -
        2 * k * k * d * d + e * e * e * e + 2 * k * k * e * e;
    numsols = rootof(a, b, c, &s1, &s2);
    if (numsols > 0) {
      tx1 = s1;
      ty1 = -0.5 * (2 * e * s1 - k * k - d * d - e * e + s * s) / d;
      tx2 = s2;
      ty2 = -0.5 * (2 * e * s2 - k * k - d * d - e * e + s * s) / d;
    } else {
      *y = 0;
      *x = 0;
      return 0;
    }
    t = e;
    e = d;
    d = t;
  } else {
    a = 4 * d * d + 4 * e * e;
    b = (-4 * e * e * e - 4 * k * k * e + 4 * e * s * s - 4 * d * d * e);
    c = d * d * d * d + k * k * k * k - 2 * e * e * s * s + 2 * d * d * e * e -
        2 * d * d * s * s + s * s * s * s - 2 * k * k * s * s -
        2 * k * k * d * d + e * e * e * e + 2 * k * k * e * e;
    numsols = rootof(a, b, c, &s1, &s2);
    if (numsols > 0) {
      ty1 = s1;
      tx1 = -0.5 * (2 * e * s1 - k * k - d * d - e * e + s * s) / d;
      ty2 = s2;
      tx2 = -0.5 * (2 * e * s2 - k * k - d * d - e * e + s * s) / d;
    } else {
      *y = 0;
      *x = 0;
      return 0;
    }
  }

  if ((type * (-tx1 * e + ty1 * d)) < 0) {
    *x = tx2;
    *y = ty2;
  } else {
    *x = tx1;
    *y = ty1;
  }
  return numsols;
}

/* Function:   intlinecirc
 *
 * Find the intersection between a line and a circle.
 */
int intlinecirc(float a, float b, float d, float e, float r, float *s1,
                float *s2) {
  return rootof((a * a + b * b), (-2 * e * b - 2 * d * a),
                (d * d + e * e - r * r), s1, s2);
}

/* Function:   anglebetween
 *
 * Calculates the angle between lines passing through the origin and with
 * real/imag components ar,ai and br,bi respectively.
 */
/* OLD */
float anglebetween(float ar, float ai, float br, float bi) {
  float rr, ri;
  rr = -ar * br - ai * bi;
  ri = ar * bi - ai * br;
  return atan2f(ri, rr);
}

/* NEW */
float anglebetween_new(float ar, float ai, float br, float bi) {
  float rr, ri;
  rr = ar * br + ai * bi;
  ri = -ar * bi + ai * br;
  return atan2f(ri, rr);
}

/* Function:  rto
 *
 * ???
 */
/* CHANGED - added *k to param list, replaced end of code with new */
int rto(float Gc, float kdes, float Smax, float Gmax, float gam, float GSP,
        float *g, float *k) {
  float ta, tb, tc, tba, karea, dkarea, ddkarea;
  int ai;
  float tb_leftover;
  float tc_leftover;

  /* ta:  ramp from Gc down to 0. */
  /* tb:  triangle above Gc */
  /* tc:  plateau at Gmax */

  ta = Gc / (Smax / GSP);
  karea = gam / 2 * (Smax / GSP) * ta *
          ta; /* Area achievable with an immediate ramp from Gc down to 0*/
  dkarea = kdes - karea;
  if (dkarea < 0) {
    return -1;
  }
  tb = (-gam * Gc +
        sqrt((gam * gam * Gc * Gc) - (gam * Smax / GSP) * (-dkarea))) /
       (gam * Smax / GSP) * 2;
  if (((ta + tb / 2) * Smax / GSP) > Gmax) {
    tb = 2 * (Gmax / (Smax / GSP) - ta);
    ddkarea = kdes - (karea + gam * Smax / GSP * tb * tb / 4.0 + gam * Gc * tb);
    tc = ddkarea / (gam * Gmax);
  } else {
    tc = 0;
  }

  tba = ta + tb / 2;

  /* NEW CODE */
  if (g != 0) {
    k[0] = kdes / gam;
    /* Start in the k-domain */
    for (ai = 1; ai <= floor(tb / GSP / 2); ai++) {
      k[ai] = k[0] - Gc * ai * GSP - (ai * GSP * ai * GSP) / 2 * (Smax / GSP);
    }
    tb_leftover = tb / GSP / 2 - floor(tb / GSP / 2);

    if ((tc / GSP + tb_leftover) < 1) {
      tc_leftover = tc / GSP + tb_leftover;
    } else {
      for (ai = (int)floor(tb / GSP / 2) + 1;
           ai <= (int)floor(tb / GSP / 2 + tc / GSP); ai++) {
        k[ai] = k[0] - Gc * tb / 2 - tb * tb / 8 * Smax / GSP -
                (ai * GSP - tb / 2) * Gmax;
      }
      tc_leftover = tb / GSP / 2 + tc / GSP - floor(tb / GSP / 2 + tc / GSP);
    }

    if ((tc_leftover + tba / GSP) < 1) {
      k[(int)ceil(tb / GSP / 2 + tc / GSP + tba / GSP)] = 0;
    } else {
      for (ai = (int)floor(tb / GSP / 2 + tc / GSP) + 1;
           ai <= (int)floor(tb / GSP / 2 + tc / GSP + tba / GSP); ai++) {
        k[ai] = k[0] - Gc * tb / 2 - tb * tb / 8 * Smax / GSP - (tc * Gmax) -
                (ai * GSP - tb / 2 - tc) * (Gc + tb / 2 * Smax / GSP) +
                (ai * GSP - tb / 2 - tc) * (ai * GSP - tb / 2 - tc) / 2 * Smax /
                    GSP;
      }
      k[(int)ceil(tb / GSP / 2 + tc / GSP + tba / GSP)] = 0;
    }

    for (ai = (int)ceil(tb / GSP / 2 + tc / GSP + tba / GSP); ai >= 1; ai--) {
      g[(int)ceil(tb / GSP / 2 + tc / GSP + tba / GSP) - ai + 1] =
          (k[ai - 1] - k[ai]) / GSP;
    }
    g[0] = 0;
  }

  return (int)ceil(tb / GSP / 2 + tc / GSP + tba / GSP) + 1;
}

/* Function:  rtz
 *
 * ???
 */
/* CHANGED */
int rtz(float gcx, float gcy, float gcz, float kcx, float kcy, float kcz,
        float theta, int maxai, float Gmax_xy, float Gmax_z, float Smax_xy,
        float Smax_z, float GSP, float *pgx, float *pgy, float *pgz, float *pkx,
        float *pky, float *pkz, float FOV1, float FOV2, float NINT, float kmax,
        float dcf) {

  float dang;
  float kxy, gxy;
  float kxymin, kxymax;
  int ai;
  int bi;
  int done;
  int fail;
  int lenai;
  int ttot;
  float s1, s2;
  float s1min, s2min;
  float s1max, s2max;
  float s1l, s2l;
  float Smax_c;
  float Gmax_c;
  float kr_next;
  float kr_cur;
  float n_des;
  float n_act = 0.0;
  int numsols, numsolsmin;
  float t1, tmin, tmax, tl;
  int numsolsline, numsolsmax;
  int numsolstmax;
  float range_min;
  float range_max;
  float kang, tt, tkx, tky, tkz;
  float s1kxy, s2kxy;
  int numsolskxy;
  float t1xy, t2xy;
  float s1tmax, s2tmax;
  float ttmax;

  pgx[0] = gcx;
  pgy[0] = gcy;
  pgz[0] = gcz;
  pkx[0] = kcx / GAMMA_PROTON / GSP;
  pky[0] = kcy / GAMMA_PROTON / GSP;
  pkz[0] = kcz / GAMMA_PROTON / GSP;
  ai = 0;
  done = 0;
  fail = 0;
  lenai = 1;
  /*MCARL: moved from below*/
  Smax_c = min(Smax_z / cos(theta), Smax_xy / sin(theta));
  Gmax_c = min(Gmax_z / cos(theta), Gmax_xy / sin(theta));
  kxy = sqrt(pkx[ai] * pkx[ai] + pky[ai] * pky[ai]);
  gxy = sqrt(pgx[ai] * pgx[ai] + pgy[ai] * pgy[ai]);
  dang = anglebetween(pgx[ai], pgy[ai], pkx[ai], pky[ai]);

  while ((done == 0) && (fail == 0) && (ai < (maxai))) {

    kxymin = (pkz[ai] + max(-Gmax_z, (pgz[ai] - Smax_z))) * tan(theta);
    kxymax = (pkz[ai] + min(Gmax_z, (pgz[ai] + Smax_z))) * tan(theta);

    /* Find intersection of SMAX with Gc */
    numsols = circcirc(pgx[ai], pgy[ai], gxy, Smax_xy, -1, &s1, &s2);
    pgx[ai + 1] = s1;
    pgy[ai + 1] = s2;
    s1 = s1 + pkx[ai];
    s2 = s2 + pky[ai];
    /* Find intersection of Gc with xy min (Z slew and amp constraints) */
    numsolsmin = circcirc(pkx[ai], pky[ai], kxymin, gxy, -1, &s1min, &s2min);
    /* NEW CODE */
    numsolsmax = circcirc(pkx[ai], pky[ai], kxymax, gxy, -1, &s1max, &s2max);

    /* Find intersection of Gc with line back to origin*/
    numsolsline = intlinecirc(pkx[ai] / kxy, pky[ai] / kxy, pkx[ai], pky[ai],
                              gxy, &s1l, &s2l);
    s1l = s2l * pkx[ai] / kxy;
    s2l = s2l * pky[ai] / kxy;

    /* Find intersection of kxymin and SMAX */
    numsolskxy = circcirc(pkx[ai] + pgx[ai], pky[ai] + pgy[ai], kxymin, Smax_xy,
                          1, &s1kxy, &s2kxy);

    /* Find intersection of kxy and Gc */
    numsolstmax = circcirc(pkx[ai], pky[ai], kxy, gxy, -1, &s1tmax, &s2tmax);
    ttmax = anglebetween_new(s1tmax - pkx[ai], s2tmax - pky[ai], -pkx[ai],
                             -pky[ai]);
    /* Convert to parametric form */
    if (numsols > 0) {
      t1 = anglebetween_new(s1 - pkx[ai], s2 - pky[ai], -pkx[ai], -pky[ai]);
    } else {
      t1 = 0;
    }

    if (numsolskxy > 0) {
      t1xy = anglebetween_new(pkx[ai], pky[ai], s1kxy, s2kxy);
    } else {
      t1xy = 0;
    }
    if (numsolsmin > 0) {
      t2xy = anglebetween_new(pkx[ai], pky[ai], s1min, s2min);
    } else {
      t2xy = 0;
    }

    if (t1xy < 0) {
      t1xy = 0;
    }

    tmin =
        anglebetween_new(s1min - pkx[ai], s2min - pky[ai], -pkx[ai], -pky[ai]);
    tmax =
        anglebetween_new(s1max - pkx[ai], s2max - pky[ai], -pkx[ai], -pky[ai]);
    tl = anglebetween_new(pgx[ai], pgy[ai], -pkx[ai], -pky[ai]);

    range_min = max(0, t1);
    range_min = (numsolsmin > 0) ? max(tmin, range_min) : range_min;
    range_max = ((numsolsmax > 0) && (tmax > 0)) ? min(tmax, tl) : tl;
    range_max = min(range_max, ttmax);

    kang = atan2f(-pky[ai], -pkx[ai]);

    /* Ensure a point in the range exists. If not, fail!*/
    if (range_max - range_min >= 0) {

      /* 1. Check at range_max (ensure that at least one point has n_act>=n_des.
       * If not, fail!) */
      tt = range_max;
      tkx = pkx[ai] + gxy * (cos(kang) * cos(tt) - sin(kang) * sin(tt));
      tky = pky[ai] + gxy * (cos(kang) * sin(tt) + sin(kang) * cos(tt));
      tkz = sqrt(tkx * tkx + tky * tky) / tan(theta);
      kr_cur = GAMMA_PROTON * GSP *
               sqrt(pkx[ai] * pkx[ai] + pky[ai] * pky[ai] + pkz[ai] * pkz[ai]);
      kr_next = GAMMA_PROTON * GSP * sqrt(tkx * tkx + tky * tky + tkz * tkz);
      n_des = calcn(kr_next, sin(theta), FOV1, FOV2, NINT, kmax, dcf);
      n_act = anglebetween_new(pkx[ai], pky[ai], tkx, tky) / (kr_cur - kr_next);
      pkx[ai + 1] = tkx;
      pky[ai + 1] = tky;
      pkz[ai + 1] = tkz;

      if (n_act >= n_des) {
        /* 2. Check at range_min (if range_min has n_act>=n_des, don't need to
         * search) */
        tt = range_min;
        tkx = pkx[ai] + gxy * (cos(kang) * cos(tt) - sin(kang) * sin(tt));
        tky = pky[ai] + gxy * (cos(kang) * sin(tt) + sin(kang) * cos(tt));
        tkz = sqrt(tkx * tkx + tky * tky) / tan(theta);
        kr_cur =
            GAMMA_PROTON * GSP *
            sqrt(pkx[ai] * pkx[ai] + pky[ai] * pky[ai] + pkz[ai] * pkz[ai]);
        kr_next = GAMMA_PROTON * GSP * sqrt(tkx * tkx + tky * tky + tkz * tkz);
        n_des = calcn(kr_next, sin(theta), FOV1, FOV2, NINT, kmax, dcf);
        n_act =
            anglebetween_new(pkx[ai], pky[ai], tkx, tky) / (kr_cur - kr_next);
        if (n_act >= n_des) {
          pkx[ai + 1] = tkx;
          pky[ai + 1] = tky;
          pkz[ai + 1] = tkz;
        }

        if (n_act < n_des) { /*if range_min has n_act<n_des, search for the
                                crossover point */
          /* 3. Do bissection 10 or so times to find the crossover point */
          for (bi = 0; bi < 10; bi++) {
            tt = range_min + (range_max - range_min) / 2.0;
            tkx = pkx[ai] + gxy * (cos(kang) * cos(tt) - sin(kang) * sin(tt));
            tky = pky[ai] + gxy * (cos(kang) * sin(tt) + sin(kang) * cos(tt));
            tkz = sqrt(tkx * tkx + tky * tky) / tan(theta);
            kr_cur =
                GAMMA_PROTON * GSP *
                sqrt(pkx[ai] * pkx[ai] + pky[ai] * pky[ai] + pkz[ai] * pkz[ai]);
            kr_next =
                GAMMA_PROTON * GSP * sqrt(tkx * tkx + tky * tky + tkz * tkz);
            n_des = calcn(kr_next, sin(theta), FOV1, FOV2, NINT, kmax, dcf);
            n_act = anglebetween_new(pkx[ai], pky[ai], tkx, tky) /
                    (kr_cur - kr_next);
            if (n_act >= n_des) {
              pkx[ai + 1] = tkx;
              pky[ai + 1] = tky;
              pkz[ai + 1] = tkz;
              range_max = tt;
            } else {
              range_min = tt;
            }
          }
        } else { /*range_min has n_act>=n_des*/
          if ((numsolsmin > 0) && (max(0, tmin) > max(0, t1))) {
            /* Check at t1xy */
            tkx = (pkx[ai] * cos(t1xy) + pky[ai] * sin(t1xy)) / kxy * kxymin;
            tky = (-pkx[ai] * sin(t1xy) + pky[ai] * cos(t1xy)) / kxy * kxymin;
            tkz = sqrt(tkx * tkx + tky * tky) / tan(theta);
            kr_cur =
                GAMMA_PROTON * GSP *
                sqrt(pkx[ai] * pkx[ai] + pky[ai] * pky[ai] + pkz[ai] * pkz[ai]);
            kr_next =
                GAMMA_PROTON * GSP * sqrt(tkx * tkx + tky * tky + tkz * tkz);
            n_des = calcn(kr_next, sin(theta), FOV1, FOV2, NINT, kmax, dcf);
            n_act = anglebetween_new(pkx[ai], pky[ai], tkx, tky) /
                    (kr_cur - kr_next);
            if (n_act >= n_des) {
              pkx[ai + 1] = tkx;
              pky[ai + 1] = tky;
              pkz[ai + 1] = tkz;
            }
            if (n_act < n_des) { /* find crossover point */
              for (bi = 0; bi < 10; bi++) {
                tt = t1xy + (t2xy - t1xy) / 2.0;
                tkx = (pkx[ai] * cos(tt) + pky[ai] * sin(tt)) / kxy * kxymin;
                tky = (-pkx[ai] * sin(tt) + pky[ai] * cos(tt)) / kxy * kxymin;
                tkz = sqrt(tkx * tkx + tky * tky) / tan(theta);
                kr_cur = GAMMA_PROTON * GSP *
                         sqrt(pkx[ai] * pkx[ai] + pky[ai] * pky[ai] +
                              pkz[ai] * pkz[ai]);
                kr_next = GAMMA_PROTON * GSP *
                          sqrt(tkx * tkx + tky * tky + tkz * tkz);
                n_des = calcn(kr_next, sin(theta), FOV1, FOV2, NINT, kmax, dcf);
                n_act = anglebetween_new(pkx[ai], pky[ai], tkx, tky) /
                        (kr_cur - kr_next);
                if (n_act >= n_des) {
                  pkx[ai + 1] = tkx;
                  pky[ai + 1] = tky;
                  pkz[ai + 1] = tkz;
                  t2xy = tt;
                } else {
                  t1xy = tt;
                }
              }
            }
          }
        }
        pgz[ai + 1] = pkz[ai + 1] - pkz[ai];
        pgy[ai + 1] = pky[ai + 1] - pky[ai];
        pgx[ai + 1] = pkx[ai + 1] - pkx[ai];

      } else {
        fail = 1; /* Due to there being no points in the range which have
                     n_act>n_des. */
      }
    } else {
      fail = 1; /* Due to no points existing */
    }

    if (n_act < 1e-4) {
      done = 1;
    }

    kxy = sqrt(pkx[ai + 1] * pkx[ai + 1] + pky[ai + 1] * pky[ai + 1]);
    gxy = sqrt(pgx[ai + 1] * pgx[ai + 1] + pgy[ai + 1] * pgy[ai + 1]);
    ttot = rto(sqrt(gxy * gxy + pgz[ai + 1] * pgz[ai + 1]),
               sqrt(kxy * kxy + pkz[ai + 1] * pkz[ai + 1]) * GAMMA_PROTON * GSP,
               Smax_c, Gmax_c, GAMMA_PROTON, GSP, 0, 0);

    if (ttot == -1) {
      fail = 1;
    }

    ai = ai + 1;

  } /* End while done==0 & fail==0 */
  if ((done == 1) && (fail == 0))
    return ai;
  else {
    return -1;
  }
}

/* Function:  calcn
 *
 * Calculates the "n" value from which the number of interleaves can be
 * determined?
 */
/* CHANGED - param changed from theta to sintheta, some internal logic */
float calcn(float k, float sintheta, float FOV1, float FOV2, float NINT,
            float kmax, float dcf) {
  float Gtwist;
  float n;
  if (dcf > 0) {
    Gtwist = pow(2 * PI * sintheta * k * FOV2 / NINT * k / kmax, 2) -
             FOV2 * FOV2 / (FOV1 * FOV1);
    Gtwist = sqrt((Gtwist > 0) ? Gtwist : 0);
    n = Gtwist / (k * sintheta);
  } else {
    Gtwist = pow(2 * PI * sintheta * k * FOV2 / NINT, 2) -
             FOV2 * FOV2 / (FOV1 * FOV1);
    Gtwist = sqrt((Gtwist > 0) ? Gtwist : 0);
    n = Gtwist / (k * sintheta);
  }
  return n;
}

/* Function:  calcdn
 *
 * Calculate incremental "n" value
 */
/* CHANGED - param changed from theta to sintheta, some internal logic */
float calcdn(float k, float sintheta, float FOV1, float FOV2, float NINT,
             float kmax, float dcf) {
  float Gtwist;
  float dn;
  if (dcf > 0) {
    Gtwist = pow(2 * PI * sintheta * k * FOV2 / NINT * k / kmax, 2) -
             pow(FOV2 / FOV1, 2);
    Gtwist = sqrt((Gtwist > 0) ? Gtwist : 0);
    dn = 8 * k * k * sintheta * PI * PI * FOV2 * FOV2 /
             (Gtwist * NINT * NINT * kmax * kmax) -
         Gtwist / (k * k * sintheta);
  } else {
    Gtwist = pow(2 * PI * sintheta * k * FOV2 / NINT, 2) - pow(FOV2 / FOV1, 2);
    Gtwist = sqrt((Gtwist > 0) ? Gtwist : 0);
    dn = 4 * sintheta * PI * PI * FOV2 * FOV2 / (Gtwist * NINT * NINT) -
         Gtwist / (k * k * sintheta);
  }
  return dn;
}

/* VD Start */

float calcc1(float sintheta, float FOVrad, float nintl) {
  return 4 * PI * PI * sintheta * sintheta * FOVrad * FOVrad / (nintl * nintl);
}

float calcc2(float FOVcir, float FOVrad) {
  return FOVrad * FOVrad / (FOVcir * FOVcir);
}

float calckstart(float theta, float FOVcir, float FOVrad, float nintl,
                 float kmax) {
  float c1 = calcc1(sin(theta), FOVrad, nintl);
  float c2 = calcc2(FOVcir, FOVrad);
  float kstart;
  kstart = sqrt(c1 * c2) / c1;
  return kstart;
}

void getFOVs(float kr, float *FOVcir, float *FOVrad, float *FOVkr, int numFOV,
             float *fc, float *fr) {
  int n = -1;
  float nextFrac; /* fraction of next value to use */

  /* Find first index where kr is larger than kr in fov array */
  while (((n + 1) < numFOV) && (kr > FOVkr[n + 1])) {
    n++;
  }

  if (n == -1) {
    /* smaller than all FOV values
    use first value */
    *fc = FOVcir[0];
    *fr = FOVrad[0];
  } else if (n == (numFOV - 1)) {
    /* larger than all FOV values
    use last value */
    *fc = FOVcir[n];
    *fr = FOVrad[n];
  } else {
    nextFrac = (kr - FOVkr[n]) / (FOVkr[n + 1] - FOVkr[n]);
    *fc = (1 - nextFrac) * FOVcir[n] + nextFrac * FOVcir[n + 1];
    *fr = (1 - nextFrac) * FOVrad[n] + nextFrac * FOVrad[n + 1];
  }

  return;
}

float calckstart_vd(float theta, float *FOVcir, float *FOVrad, float *FOVkr,
                    int numFOV, float kmax, float nintl) {
  float kstart, kstartNew;
  float c1, c2;
  float fr, fc;
  float sintheta = sin(theta);
  int iter = 0;

  kstart = 0;
  kstartNew = 0;
  do {
    kstart = kstartNew;
    getFOVs(kstart, FOVcir, FOVrad, FOVkr, numFOV, &fc, &fr);
    c1 = calcc1(sintheta, fr, nintl);
    c2 = calcc2(fc, fr);
    kstartNew = calckstart(theta, fc, fr, nintl, kmax);
    iter++;
  } while ((kstartNew > kstart) && (kstartNew < kmax) && (iter < 20));

  if (kstartNew > kstart)
    kstart = kstartNew;

  return kstart;
}

/* VD */

/* VD */

int whirlcone_vd(float theta, float *FOV1, float *FOV2, float *FOV3,
                 float *FOV4, float dcf, float kmax, float NINT, int numpt,
                 float GSP, float SMAX, float GMAX, float Smax_xy, float Smax_z,
                 float Gmax_xy, float Gmax_z, float dens_min, float ***pg,
                 float ***pk, float *pai, float *pn, int acq_mode) {

  float *pkx, *pky, *pkz, *pgx, *pgy, *pgz;
  float *pkfx, *pkfy, *pkfz, *pgfx, *pgfy, *pgfz;
  float *prkx, *prky, *prkz, *prgx, *prgy, *prgz;
  float *pkr;
  float *pcsndk;
  float kstart;
  float n;
  float dn;
  float dk;
  float d2k;
  float cmp;
  float ggmax;
  float a, b, c, d;
  float descrim;
  float SMS;
  float gxyend;
  float kxyend;
  float gxyzend;
  float kxyzend;
  float mg;

  float k;
  float maxs = 0.0;
  float maxg = 0.0;
  float Smax_c;
  float Gmax_c;
  float gxy_in;
  float gxy_old;
  float sintheta;
  float costheta;
  float n_kmax, kmaxg, Gc, n_cur, phi_high, phi_low, dphi, gzt;
  float Gc_top, Gc_bot;
  float gxt, gyt, kxt, kyt, kzt;
  float rgxt = 0.0, rgyt = 0.0, rgzt = 0.0;
  int rtzlen;
  int rtolen;
  int dlen;
  int ai, bi, ci, atmp;
  int rangelow, rangehigh;
  int maxrangehigh;
  float dphimax;

  /* VD */
  float FOVcir;
  float FOVrad;
  int numfkr = 100;

  /* VD */

  int WarnMeIfIEverGoHere = 0; /*MCARL*/
  int WarnMeIfIGoHere = 0;     /*MCARL*/
  int WarnMeIfIGoHereXXX = 0;  /*MCARL*/

  /* Allocate memory for desgin arrays */
  pkr = (float *)AllocMem(numpt * sizeof(float));
  pcsndk = (float *)AllocMem(numpt * sizeof(float));
  prgx = (float *)AllocMem(numpt * sizeof(float));
  prgy = (float *)AllocMem(numpt * sizeof(float));
  prgz = (float *)AllocMem(numpt * sizeof(float));
  prkx = (float *)AllocMem(numpt * sizeof(float));
  prky = (float *)AllocMem(numpt * sizeof(float));
  prkz = (float *)AllocMem(numpt * sizeof(float));
  pgx = (float *)AllocMem(numpt * sizeof(float));
  pgy = (float *)AllocMem(numpt * sizeof(float));
  pgz = (float *)AllocMem(numpt * sizeof(float));
  pkx = (float *)AllocMem(numpt * sizeof(float));
  pky = (float *)AllocMem(numpt * sizeof(float));
  pkz = (float *)AllocMem(numpt * sizeof(float));

  /* Initialize rtolen to some very large number so we default to waveform
     length of npts when rtzlen<0 */
  rtolen = 10000000;
  rtzlen = 0; /*MCARL*/

  /* Currently just going to create a down ramp to bring gradients back to
     zero at end of waveform. To preserve gradient waveform size for
     implemenation
     on scanner, calculate max. ramp points needed to ramp down from max grad */
  /*  nramp = (int)ceil(Gmax_xy/Smax_xy + 1); */
  /* Pass this in now to make sure it isn't messed up by xy/z scaling factors */

  pgfx = (*pg)[0];
  pgfy = (*pg)[1];
  pgfz = (*pg)[2];
  pkfx = (*pk)[0];
  pkfy = (*pk)[1];
  pkfz = (*pk)[2];

  sintheta = sin((double)theta);
  costheta = cos((double)theta);

  /* Find the point that the twist begins */

  /* VD */
  kstart = calckstart_vd(theta, FOV1, FOV2, FOV3, numfkr, kmax, NINT);
  /* kstart = NINT/(2*PI*sintheta*FOV1); */
  /* VD */
  SMS = Smax_xy * Smax_xy * GAMMA_PROTON * GAMMA_PROTON /
        (sintheta * sintheta * GSP * GSP);

  if (kstart > (kmax - 0.0001) || acq_mode == 2) {
    k = 0;
    mg = 0;
    ai = 0;
    pgfx[ai] = 0;
    pgfy[ai] = 0;
    pgfz[ai] = 0;
    pkfx[ai] = 0;
    pkfy[ai] = 0;
    pkfz[ai] = 0;
    if (acq_mode == 1) {
      if (theta < PI / 4.0) {
        maxs = Smax_z * sqrt(1 + tan(theta) * tan(theta));
        maxg = Gmax_z * sqrt(1 + tan(theta) * tan(theta));
      } else {
        maxs = Smax_xy * sqrt(1 + 1 / (tan(theta) * tan(theta)));
        maxg = Gmax_xy * sqrt(1 + 1 / (tan(theta) * tan(theta)));
      }
    }
    if (acq_mode == 2) {
      maxs = SMAX;
      maxg = GMAX;
    }
    while ((k <= kmax) && (ai < numpt)) {
      mg = min(maxg, mg + maxs);
      ai++;
      pgfx[ai] = mg * sintheta;
      pgfy[ai] = 0;
      pgfz[ai] = mg * costheta;
      pkfx[ai] = pkfx[ai - 1] + pgfx[ai] * GAMMA_PROTON * GSP;
      pkfy[ai] = pkfy[ai - 1] + pgfy[ai] * GAMMA_PROTON * GSP;
      pkfz[ai] = pkfz[ai - 1] + pgfz[ai] * GAMMA_PROTON * GSP;
      k = sqrt(pkfx[ai] * pkfx[ai] + pkfy[ai] * pkfy[ai] +
               pkfz[ai] * pkfz[ai]); /*MCARL: Fix FOV problem for PR*/
    }
    *pai = ai;

  } else { /* kstart>(kmax-0.0001) || acq_mode==2 */

    pkr[0] = kstart + 0.0001;
    pkx[0] = pkr[0] * sintheta;
    pky[0] = 0;
    pkz[0] = pkr[0] * costheta;
    pgx[0] = 0;
    pgy[0] = 0;
    pgz[0] = 0;
    pcsndk[0] = 0;
    dk = 0;
    ai = 0;

    while ((ai < (numpt - 1)) && (pkr[ai] < kmax)) {
      getFOVs(pkr[ai], FOV1, FOV2, FOV3, numfkr, &FOVcir, &FOVrad);
      n = calcn(pkr[ai], sintheta, FOVcir, FOVrad, NINT, kmax, dcf);
      n = (n > 0.00001) ? n : (0.00001);
      dn = calcdn(pkr[ai], sintheta, FOVcir, FOVrad, NINT, kmax, dcf);
      /* Solve differential equation for Smax in x-y */
      a = -dk * dk * n * n * pkr[ai];
      b = 1;
      c = dk * dk * (2 * n + pkr[ai] * dn);
      d = pkr[ai] * n;
      descrim = (2 * d * c * b * a - b * b * c * c + b * b * SMS -
                 d * d * a * a + d * d * SMS);
      descrim = (descrim > 0) ? descrim : 0;
      d2k = -((d * c + b * a) - sqrt(descrim)) / (b * b + d * d);
      /* Apply Smax constraint in z */
      cmp = Smax_z * GAMMA_PROTON / costheta / GSP;
      d2k = (d2k < cmp) ? d2k : cmp;
      d2k = (d2k > (-cmp)) ? d2k : (-cmp);
      dk = dk + d2k * GSP;
      /* Apply Gmax constraint in x-y */
      cmp = Gmax_xy * GAMMA_PROTON /
            (sqrt(1 + pkr[ai] * pkr[ai] * n * n) * sintheta);
      dk = (dk < cmp) ? dk : cmp;
      /* Apply Gmax constraint in z */
      cmp = Gmax_z * GAMMA_PROTON / costheta;
      dk = (dk < cmp) ? dk : cmp;
      /* Apply constant density constraint */
      if (dens_min > 0) {
        ggmax = NINT *
                sqrt((n * n * pkr[ai] * pkr[ai] * sintheta * sintheta) + 1) /
                (2 * PI * pkr[ai] * pkr[ai] * sintheta * dens_min);
        cmp = ggmax /
              sqrt(1 + pkr[ai] * pkr[ai] * n * n * sintheta * sintheta) *
              GAMMA_PROTON;
        dk = (dk < cmp) ? dk : cmp;
      }

      pkr[ai + 1] = pkr[ai] + dk * GSP;
      pcsndk[ai + 1] = pcsndk[ai] + dk * n;
      ai++;
    }
    dlen = ai;

    for (ai = 1; ai < dlen; ai++) {
      pkx[ai] = pkr[ai] * cos(pcsndk[ai] * GSP) * sintheta;
      pky[ai] = pkr[ai] * sin(pcsndk[ai] * GSP) * sintheta;
      pkz[ai] = pkr[ai] * costheta;
      pgx[ai] = (pkx[ai] - pkx[ai - 1]) / (GSP * GAMMA_PROTON);
      pgy[ai] = (pky[ai] - pky[ai - 1]) / (GSP * GAMMA_PROTON);
      pgz[ai] = (pkz[ai] - pkz[ai - 1]) / (GSP * GAMMA_PROTON);
    }

    rangelow = 1;
    rangehigh = 1;
    rtzlen = 0;
    gxy_old = 0;

    ai = 1;
    gxy_in = sqrt(pgx[ai] * pgx[ai] + pgy[ai] * pgy[ai]);
    while ((gxy_in > (gxy_old + Smax_xy / 1000)) && (ai < dlen - 1)) {
      gxy_old = gxy_in;
      gxy_in = sqrt(pgx[ai + 1] * pgx[ai + 1] + pgy[ai + 1] * pgy[ai + 1]);
      ai++;
    }
    maxrangehigh = ai - 1;

    while ((rtzlen != -1) && (rangehigh < maxrangehigh)) {
      rangehigh = min(maxrangehigh, rangehigh * 2);
      ai = rangehigh;
      gxy_in = sqrt(pgx[ai + 1] * pgx[ai + 1] + pgy[ai + 1] * pgy[ai + 1]);
      rtzlen = rtz(-pgx[ai + 1], -pgy[ai + 1], -pgz[ai + 1], pkx[ai], pky[ai],
                   pkz[ai], theta, numpt, Gmax_xy, Gmax_z, Smax_xy, Smax_z, GSP,
                   prgx, prgy, prgz, prkx, prky, prkz, FOVcir, FOVrad, NINT,
                   kmax, dcf);
    }

    if ((rangehigh == dlen - 2) && (rtzlen > 0)) {

      WarnMeIfIEverGoHere = WarnMeIfIEverGoHere + 1;
      /*printf("MCARL: WarnMeIfIEverGoHere=%d\n",WarnMeIfIEverGoHere);
      fflush(stdout);*/

      n_kmax = calcn(kmax, sin(theta), FOVcir, FOVrad, NINT, kmax, dcf);
      kmaxg = kmax / (GAMMA_PROTON * GSP) * sin(theta);
      /* We reached the end of the pre-canned waveform, which means that we */
      /* should search for the right Gc value to start with at kmax */

      /* Search over all Gc */
      Gc_top = Gmax_xy;
      Gc_bot = 0;
      for (ci = 0; ci < 10; ci++) {
        Gc = Gc_bot + (Gc_top - Gc_bot) / 2;

        /* Start by finding an appropriate starting seed */
        circcirc(kmaxg, 0, kmaxg, Gc, 1, &kxt, &kyt);
        dphimax = atan2f(kyt, kxt - kmaxg);
        gxt = Gc * cos(dphimax);
        gyt = Gc * sin(dphimax);
        kxt = kmaxg + gxt;
        kyt = gyt;
        kzt = sqrt(kxt * kxt + kyt * kyt) / tan(theta);
        gzt = kzt - kmaxg / tan(theta);
        n_cur = atan2f(kyt, kxt) /
                (GAMMA_PROTON * GSP *
                 max(1e-8, (sqrt(kxt * kxt + kyt * kyt + kzt * kzt) -
                            kmaxg / sin(theta))));
        if (n_cur > n_kmax) {
          phi_high = dphimax;
          phi_low = 0;
          for (bi = 0; bi < 20; bi++) {
            dphi = phi_low + (phi_high - phi_low) / 2;
            gxt = Gc * cos(dphi);
            gyt = Gc * sin(dphi);
            kxt = kmaxg + gxt;
            kyt = gyt;
            kzt = sqrt(kxt * kxt + kyt * kyt) / tan(theta);
            gzt = kzt - kmaxg / tan(theta);
            n_cur = atan2f(kyt, kxt) /
                    (GAMMA_PROTON * GSP *
                     max(1e-8, (sqrt(kxt * kxt + kyt * kyt + kzt * kzt) -
                                kmaxg / sin(theta))));
            if (n_cur > n_kmax) {
              phi_high = dphi;
            } else {
              phi_low = dphi;
            }
          }
          dphi = phi_high;
          gxt = Gc * cos(dphi);
          gyt = Gc * sin(dphi);
          kxt = kmaxg + gxt;
          kyt = gyt;
          kzt = sqrt(kxt * kxt + kyt * kyt) / tan(theta);
          gzt = kzt - kmaxg / tan(theta);

          rtzlen = rtz(-gxt, -gyt, -gzt, kmaxg * GAMMA_PROTON * GSP,
                       0 * GAMMA_PROTON * GSP,
                       kmaxg / tan(theta) * GAMMA_PROTON * GSP, theta, numpt,
                       Gmax_xy, Gmax_z, Smax_xy, Smax_z, GSP, prgx, prgy, prgz,
                       prkx, prky, prkz, FOVcir, FOVrad, NINT, kmax, dcf);
        } else {
          rtzlen = -1;
        }

        if (rtzlen > 0) {
          Gc_bot = Gc;
          rgxt = -gxt;
          rgyt = -gyt;
          rgzt = -gzt;
        } else {
          Gc_top = Gc;
        }
      }
      if (Gc_bot > 0) {
        rtzlen =
            rtz(rgxt, rgyt, rgzt, kmaxg * GAMMA_PROTON * GSP,
                0 * GAMMA_PROTON * GSP, kmaxg / tan(theta) * GAMMA_PROTON * GSP,
                theta, numpt, Gmax_xy, Gmax_z, Smax_xy, Smax_z, GSP, prgx, prgy,
                prgz, prkx, prky, prkz, FOVcir, FOVrad, NINT, kmax, dcf);
      } else {
        rtzlen = -1;
      }
      ai++;

    } else { // if ((rangehigh == dlen-2) && (rtzlen>0))
      WarnMeIfIGoHere = 0;
      WarnMeIfIGoHereXXX = 0;
      while ((rangehigh - rangelow) > 0) {
        ai = (int)ceil((rangehigh - rangelow) / 2.0) + rangelow;
        rtzlen = rtz(-pgx[ai + 1], -pgy[ai + 1], -pgz[ai + 1], pkx[ai], pky[ai],
                     pkz[ai], theta, numpt, Gmax_xy, Gmax_z, Smax_xy, Smax_z,
                     GSP, prgx, prgy, prgz, prkx, prky, prkz, FOVcir, FOVrad,
                     NINT, kmax, dcf);
        if (rtzlen == -1) {
          WarnMeIfIGoHere = WarnMeIfIGoHere + 1;
          rangehigh = ai - 1;
        } else {
          WarnMeIfIGoHereXXX = WarnMeIfIGoHereXXX + 1;
          rangelow = ai;
        }
      } /* End while (rangehigh-rangelow)>0 */

      ai = rangelow;
      rtzlen = rtz(-pgx[ai + 1], -pgy[ai + 1], -pgz[ai + 1], pkx[ai], pky[ai],
                   pkz[ai], theta, numpt, Gmax_xy, Gmax_z, Smax_xy, Smax_z, GSP,
                   prgx, prgy, prgz, prkx, prky, prkz, FOVcir, FOVrad, NINT,
                   kmax, dcf);
    }

    // rtzlen = -1; //If I force this -> Fig.3.7a
    // rtzlen is the location where rto starts

    /* COPY OVER SPOKE TRAJECTORY */
    if (rtzlen > 0) {
      Smax_c = min(Smax_z / cos(theta), Smax_xy / sin(theta));
      Gmax_c = min(Gmax_z / cos(theta), Gmax_xy / sin(theta));

      gxyend = sqrt(prgx[rtzlen] * prgx[rtzlen] + prgy[rtzlen] * prgy[rtzlen]);
      kxyend = sqrt(prkx[rtzlen] * prkx[rtzlen] + prky[rtzlen] * prky[rtzlen]);
      gxyzend = sqrt(prgx[rtzlen] * prgx[rtzlen] + prgy[rtzlen] * prgy[rtzlen] +
                     prgz[rtzlen] * prgz[rtzlen]);
      kxyzend = sqrt(prkx[rtzlen] * prkx[rtzlen] + prky[rtzlen] * prky[rtzlen] +
                     prkz[rtzlen] * prkz[rtzlen]);
      rtolen = rto(gxyzend, kxyzend * GAMMA_PROTON * GSP, Smax_c, Gmax_c,
                   GAMMA_PROTON, GSP, pgfx, pkfx);
      for (bi = 0; bi < rtolen; bi++) {
        pgfz[bi] = pgfx[bi] * cos(theta);
        pgfy[bi] = pgfx[bi] * sin(theta) * prky[rtzlen] / kxyend;
        pgfx[bi] = pgfx[bi] * sin(theta) * prkx[rtzlen] / kxyend;
        pkfx[bi] =
            ((bi > 0) ? pkfx[bi - 1] : 0) + pgfx[bi] * GAMMA_PROTON * GSP;
        pkfy[bi] =
            ((bi > 0) ? pkfy[bi - 1] : 0) + pgfy[bi] * GAMMA_PROTON * GSP;
        pkfz[bi] =
            ((bi > 0) ? pkfz[bi - 1] : 0) + pgfz[bi] * GAMMA_PROTON * GSP;
      }

      /* COPY OVER CIRCULAR TRAJECTORY */
      for (bi = 0; bi < min(rtzlen, numpt - rtolen); bi++) {
        pgfx[bi + rtolen] = -prgx[rtzlen - bi];
        pgfy[bi + rtolen] = -prgy[rtzlen - bi];
        pgfz[bi + rtolen] = -prgz[rtzlen - bi];
        pkfx[bi + rtolen] = prkx[rtzlen - bi - 1] * GAMMA_PROTON * GSP;
        pkfy[bi + rtolen] = prky[rtzlen - bi - 1] * GAMMA_PROTON * GSP;
        pkfz[bi + rtolen] = prkz[rtzlen - bi - 1] * GAMMA_PROTON * GSP;
      }
      /* COPY OVER DIFFERENTIAL TRAJECTORY */
      for (bi = ai + 1; bi < min(dlen, numpt - rtolen - rtzlen + (ai + 1));
           bi++) {
        pgfx[bi - (ai + 1) + rtolen + rtzlen] = pgx[bi];
        pgfy[bi - (ai + 1) + rtolen + rtzlen] = pgy[bi];
        pgfz[bi - (ai + 1) + rtolen + rtzlen] = pgz[bi];
        pkfx[bi - (ai + 1) + rtolen + rtzlen] = pkx[bi];
        pkfy[bi - (ai + 1) + rtolen + rtzlen] = pky[bi];
        pkfz[bi - (ai + 1) + rtolen + rtzlen] = pkz[bi];
      }

      *pai = min(numpt, dlen - ai + rtolen + rtzlen) - 1;

    } else { /*if (rtzlen>0) */

      Smax_c = min(Smax_z / cos(theta), Smax_xy / sin(theta));
      Gmax_c = min(Gmax_z / cos(theta), Gmax_xy / sin(theta));
      kxyend = sqrt(pkx[0] * pkx[0] + pky[0] * pky[0]);
      kxyzend = sqrt(pkx[0] * pkx[0] + pky[0] * pky[0] + pkz[0] * pkz[0]);
      rtolen = rto(0, kxyzend, Smax_c, Gmax_c, GAMMA_PROTON, GSP, pgfx, pkfx);
      for (bi = 0; bi < min(rtolen, numpt); bi++) {
        pgfz[bi] = pgfx[bi] * cos(theta);
        pgfy[bi] = pgfx[bi] * sin(theta) * pky[0] / kxyend;
        pgfx[bi] = pgfx[bi] * sin(theta) * pkx[0] / kxyend;
        pkfx[bi] =
            ((bi > 0) ? pkfx[bi - 1] : 0) + pgfx[bi] * GAMMA_PROTON * GSP;
        pkfy[bi] =
            ((bi > 0) ? pkfy[bi - 1] : 0) + pgfy[bi] * GAMMA_PROTON * GSP;
        pkfz[bi] =
            ((bi > 0) ? pkfz[bi - 1] : 0) + pgfz[bi] * GAMMA_PROTON * GSP;
      }

      for (bi = 0; bi < min(dlen, numpt - rtolen - 1); bi++) {
        pgfx[rtolen + bi] = pgx[bi];
        pgfy[rtolen + bi] = pgy[bi];
        pgfz[rtolen + bi] = pgz[bi];
        pkfx[rtolen + bi] = pkx[bi];
        pkfy[rtolen + bi] = pky[bi];
        pkfz[rtolen + bi] = pkz[bi];
      }
      *pai = min(dlen, numpt - rtolen - 1) + min(rtolen, numpt);
    } // end if (rtzlen>0)

  } /* kstart>(kmax-0.0001) || acq_mode==2 (If Twist needed)*/

  atmp = (int)(*pai - 1);

  /* Free up allocated memory */
  FreeMem(pkr);
  FreeMem(pcsndk);
  FreeMem(prgx);
  FreeMem(prgy);
  FreeMem(prgz);
  FreeMem(prkx);
  FreeMem(prky);
  FreeMem(prkz);
  FreeMem(pgx);
  FreeMem(pgy);
  FreeMem(pgz);
  FreeMem(pkx);
  FreeMem(pky);
  FreeMem(pkz);

  return 1;
}

/* VD End */

/* Function:  whirlcone
 *
 * Main generating function to calculate a cone trajectory given the following
 * input parameters:
 *  theta = polar angle of cone
 *  FOV1 = field of view
 *  FOV2 =
 *  FOV3 =
 *  FOV4 =
 *  dcf = density compensation factor
 *  kmax = maximum kspace extent
 *  NINT = number of interleaves
 *  GTWISTADJUST =
 *  numpt = number of points to design to
 *  GSP = gradient sample interval (sec)
 *  Smax = max system slew rate (T/m/s)
 *  Gmax = max gradient amplitude (T/m)
 *  pg = pointer to gradient waveforms (x, y and z)
 *  pk = pointer to k-space trajetory
 *  pai = pointer to waveform length per cone
 *  pn = pointer to "n" value per cone (Useless Variable???)
 */
/* CHANGED - arguments changed, new sections added */
int whirlcone(float theta, float FOV1, float FOV2, float FOV3, float FOV4,
              float dcf, float kmax, float NINT, int numpt, float GSP,
              float SMAX, float GMAX, float Smax_xy, float Smax_z,
              float Gmax_xy, float Gmax_z, float dens_min, float ***pg,
              float ***pk, float *pai, float *pn, int acq_mode) {

  float *pkx, *pky, *pkz, *pgx, *pgy, *pgz;
  float *pkfx, *pkfy, *pkfz, *pgfx, *pgfy, *pgfz;
  float *prkx, *prky, *prkz, *prgx, *prgy, *prgz;
  float *pkr;
  float *pcsndk;
  float kstart;
  float n;
  float dn;
  float dk;
  float d2k;
  float cmp;
  float ggmax;
  float a, b, c, d;
  float descrim;
  float SMS;
  float gxyend;
  float kxyend;
  float gxyzend;
  float kxyzend;
  float mg;

  float k;
  float maxs = 0.0;
  float maxg = 0.0;
  float Smax_c;
  float Gmax_c;
  float gxy_in;
  float gxy_old;
  float sintheta;
  float costheta;
  float n_kmax, kmaxg, Gc, n_cur, phi_high, phi_low, dphi, gzt;
  float Gc_top, Gc_bot;
  float gxt, gyt, kxt, kyt, kzt;
  float rgxt = 0.0;
  float rgyt = 0.0;
  float rgzt = 0.0;
  int rtzlen;
  int rtolen;
  int dlen;
  int ai, bi, ci, atmp;
  int rangelow, rangehigh;
  int maxrangehigh;
  float dphimax;

  int WarnMeIfIEverGoHere = 0; /*MCARL*/
  int WarnMeIfIGoHere = 0;     /*MCARL*/
  int WarnMeIfIGoHereXXX = 0;  /*MCARL*/

  /* Allocate memory for desgin arrays */
  pkr = (float *)AllocMem(numpt * sizeof(float));
  pcsndk = (float *)AllocMem(numpt * sizeof(float));
  prgx = (float *)AllocMem(numpt * sizeof(float));
  prgy = (float *)AllocMem(numpt * sizeof(float));
  prgz = (float *)AllocMem(numpt * sizeof(float));
  prkx = (float *)AllocMem(numpt * sizeof(float));
  prky = (float *)AllocMem(numpt * sizeof(float));
  prkz = (float *)AllocMem(numpt * sizeof(float));
  pgx = (float *)AllocMem(numpt * sizeof(float));
  pgy = (float *)AllocMem(numpt * sizeof(float));
  pgz = (float *)AllocMem(numpt * sizeof(float));
  pkx = (float *)AllocMem(numpt * sizeof(float));
  pky = (float *)AllocMem(numpt * sizeof(float));
  pkz = (float *)AllocMem(numpt * sizeof(float));

  /* Initialize rtolen to some very large number so we default to waveform
     length of npts when rtzlen<0 */
  rtolen = 10000000;
  rtzlen = 0; /*MCARL*/

  /* Currently just going to create a down ramp to bring gradients back to
     zero at end of waveform. To preserve gradient waveform size for
     implemenation
     on scanner, calculate max. ramp points needed to ramp down from max grad */
  /*  nramp = (int)ceil(Gmax_xy/Smax_xy + 1); */
  /* Pass this in now to make sure it isn't messed up by xy/z scaling factors */

  pgfx = (*pg)[0];
  pgfy = (*pg)[1];
  pgfz = (*pg)[2];
  pkfx = (*pk)[0];
  pkfy = (*pk)[1];
  pkfz = (*pk)[2];

  sintheta = sin((double)theta);
  costheta = cos((double)theta);

  /* Find the point that the twist begins */
  if (dcf > 0) {
    kstart = sqrt(NINT * kmax / (2 * PI * sintheta * FOV1));
  } else {
    kstart = NINT / (2 * PI * sintheta * FOV1);
  }
  SMS = Smax_xy * Smax_xy * GAMMA_PROTON * GAMMA_PROTON /
        (sintheta * sintheta * GSP * GSP);

  if (kstart > (kmax - 0.0001) || acq_mode == 2) {
    k = 0;
    mg = 0;
    ai = 0;
    pgfx[ai] = 0;
    pgfy[ai] = 0;
    pgfz[ai] = 0;
    pkfx[ai] = 0;
    pkfy[ai] = 0;
    pkfz[ai] = 0;
    if (acq_mode == 1) {
      if (theta < PI / 4.0) {
        maxs = Smax_z * sqrt(1 + tan(theta) * tan(theta));
        maxg = Gmax_z * sqrt(1 + tan(theta) * tan(theta));
      } else {
        maxs = Smax_xy * sqrt(1 + 1 / (tan(theta) * tan(theta)));
        maxg = Gmax_xy * sqrt(1 + 1 / (tan(theta) * tan(theta)));
      }
    }
    if (acq_mode == 2) {
      maxs = SMAX;
      maxg = GMAX;
    }
    while ((k <= kmax) && (ai < numpt)) {
      mg = min(maxg, mg + maxs);
      ai++;
      pgfx[ai] = mg * sintheta;
      pgfy[ai] = 0;
      pgfz[ai] = mg * costheta;
      pkfx[ai] = pkfx[ai - 1] + pgfx[ai] * GAMMA_PROTON * GSP;
      pkfy[ai] = pkfy[ai - 1] + pgfy[ai] * GAMMA_PROTON * GSP;
      pkfz[ai] = pkfz[ai - 1] + pgfz[ai] * GAMMA_PROTON * GSP;
      k = sqrt(pkfx[ai] * pkfx[ai] + pkfy[ai] * pkfy[ai] +
               pkfz[ai] * pkfz[ai]); /*MCARL: Fix FOV problem for PR*/
    }
    *pai = ai;

  } else { /* kstart>(kmax-0.0001) || acq_mode==2 */

    pkr[0] = kstart + 0.0001;
    pkx[0] = pkr[0] * sintheta;
    pky[0] = 0;
    pkz[0] = pkr[0] * costheta;
    pgx[0] = 0;
    pgy[0] = 0;
    pgz[0] = 0;
    pcsndk[0] = 0;
    dk = 0;
    ai = 0;

    while ((ai < (numpt - 1)) && (pkr[ai] < kmax)) {
      n = calcn(pkr[ai], sintheta, FOV1, FOV2, NINT, kmax, dcf);
      n = (n > 0.00001) ? n : (0.00001);
      dn = calcdn(pkr[ai], sintheta, FOV1, FOV2, NINT, kmax, dcf);
      /* Solve differential equation for Smax in x-y */
      a = -dk * dk * n * n * pkr[ai];
      b = 1;
      c = dk * dk * (2 * n + pkr[ai] * dn);
      d = pkr[ai] * n;
      descrim = (2 * d * c * b * a - b * b * c * c + b * b * SMS -
                 d * d * a * a + d * d * SMS);
      descrim = (descrim > 0) ? descrim : 0;
      d2k = -((d * c + b * a) - sqrt(descrim)) / (b * b + d * d);
      /* Apply Smax constraint in z */
      cmp = Smax_z * GAMMA_PROTON / costheta / GSP;
      d2k = (d2k < cmp) ? d2k : cmp;
      d2k = (d2k > (-cmp)) ? d2k : (-cmp);
      dk = dk + d2k * GSP;
      /* Apply Gmax constraint in x-y */
      cmp = Gmax_xy * GAMMA_PROTON /
            (sqrt(1 + pkr[ai] * pkr[ai] * n * n) * sintheta);
      dk = (dk < cmp) ? dk : cmp;
      /* Apply Gmax constraint in z */
      cmp = Gmax_z * GAMMA_PROTON / costheta;
      dk = (dk < cmp) ? dk : cmp;
      /* Apply constant density constraint */
      if (dens_min > 0) {
        ggmax = NINT *
                sqrt((n * n * pkr[ai] * pkr[ai] * sintheta * sintheta) + 1) /
                (2 * PI * pkr[ai] * pkr[ai] * sintheta * dens_min);
        cmp = ggmax /
              sqrt(1 + pkr[ai] * pkr[ai] * n * n * sintheta * sintheta) *
              GAMMA_PROTON;
        dk = (dk < cmp) ? dk : cmp;
      }

      pkr[ai + 1] = pkr[ai] + dk * GSP;
      pcsndk[ai + 1] = pcsndk[ai] + dk * n;
      ai++;
    }
    dlen = ai;

    for (ai = 1; ai < dlen; ai++) {
      pkx[ai] = pkr[ai] * cos(pcsndk[ai] * GSP) * sintheta;
      pky[ai] = pkr[ai] * sin(pcsndk[ai] * GSP) * sintheta;
      pkz[ai] = pkr[ai] * costheta;
      pgx[ai] = (pkx[ai] - pkx[ai - 1]) / (GSP * GAMMA_PROTON);
      pgy[ai] = (pky[ai] - pky[ai - 1]) / (GSP * GAMMA_PROTON);
      pgz[ai] = (pkz[ai] - pkz[ai - 1]) / (GSP * GAMMA_PROTON);
    }

    rangelow = 1;
    rangehigh = 1;
    rtzlen = 0;
    gxy_old = 0;

    ai = 1;
    gxy_in = sqrt(pgx[ai] * pgx[ai] + pgy[ai] * pgy[ai]);
    while ((gxy_in > (gxy_old + Smax_xy / 1000)) && (ai < dlen - 1)) {
      gxy_old = gxy_in;
      gxy_in = sqrt(pgx[ai + 1] * pgx[ai + 1] + pgy[ai + 1] * pgy[ai + 1]);
      ai++;
    }
    maxrangehigh = ai - 1;

    while ((rtzlen != -1) && (rangehigh < maxrangehigh)) {
      rangehigh = min(maxrangehigh, rangehigh * 2);
      ai = rangehigh;
      gxy_in = sqrt(pgx[ai + 1] * pgx[ai + 1] + pgy[ai + 1] * pgy[ai + 1]);
      rtzlen =
          rtz(-pgx[ai + 1], -pgy[ai + 1], -pgz[ai + 1], pkx[ai], pky[ai],
              pkz[ai], theta, numpt, Gmax_xy, Gmax_z, Smax_xy, Smax_z, GSP,
              prgx, prgy, prgz, prkx, prky, prkz, FOV1, FOV2, NINT, kmax, dcf);
    }

    if ((rangehigh == dlen - 2) && (rtzlen > 0)) {

      WarnMeIfIEverGoHere = WarnMeIfIEverGoHere + 1;
      /*printf("MCARL: WarnMeIfIEverGoHere=%d\n",WarnMeIfIEverGoHere);
      fflush(stdout);*/

      n_kmax = calcn(kmax, sin(theta), FOV1, FOV2, NINT, kmax, dcf);
      kmaxg = kmax / (GAMMA_PROTON * GSP) * sin(theta);
      /* We reached the end of the pre-canned waveform, which means that we */
      /* should search for the right Gc value to start with at kmax */

      /* Search over all Gc */
      Gc_top = Gmax_xy;
      Gc_bot = 0;
      for (ci = 0; ci < 10; ci++) {
        Gc = Gc_bot + (Gc_top - Gc_bot) / 2;

        /* Start by finding an appropriate starting seed */
        circcirc(kmaxg, 0, kmaxg, Gc, 1, &kxt, &kyt);
        dphimax = atan2f(kyt, kxt - kmaxg);
        gxt = Gc * cos(dphimax);
        gyt = Gc * sin(dphimax);
        kxt = kmaxg + gxt;
        kyt = gyt;
        kzt = sqrt(kxt * kxt + kyt * kyt) / tan(theta);
        gzt = kzt - kmaxg / tan(theta);
        n_cur = atan2f(kyt, kxt) /
                (GAMMA_PROTON * GSP *
                 max(1e-8, (sqrt(kxt * kxt + kyt * kyt + kzt * kzt) -
                            kmaxg / sin(theta))));
        if (n_cur > n_kmax) {
          phi_high = dphimax;
          phi_low = 0;
          for (bi = 0; bi < 20; bi++) {
            dphi = phi_low + (phi_high - phi_low) / 2;
            gxt = Gc * cos(dphi);
            gyt = Gc * sin(dphi);
            kxt = kmaxg + gxt;
            kyt = gyt;
            kzt = sqrt(kxt * kxt + kyt * kyt) / tan(theta);
            gzt = kzt - kmaxg / tan(theta);
            n_cur = atan2f(kyt, kxt) /
                    (GAMMA_PROTON * GSP *
                     max(1e-8, (sqrt(kxt * kxt + kyt * kyt + kzt * kzt) -
                                kmaxg / sin(theta))));
            if (n_cur > n_kmax) {
              phi_high = dphi;
            } else {
              phi_low = dphi;
            }
          }
          dphi = phi_high;
          gxt = Gc * cos(dphi);
          gyt = Gc * sin(dphi);
          kxt = kmaxg + gxt;
          kyt = gyt;
          kzt = sqrt(kxt * kxt + kyt * kyt) / tan(theta);
          gzt = kzt - kmaxg / tan(theta);

          rtzlen = rtz(-gxt, -gyt, -gzt, kmaxg * GAMMA_PROTON * GSP,
                       0 * GAMMA_PROTON * GSP,
                       kmaxg / tan(theta) * GAMMA_PROTON * GSP, theta, numpt,
                       Gmax_xy, Gmax_z, Smax_xy, Smax_z, GSP, prgx, prgy, prgz,
                       prkx, prky, prkz, FOV1, FOV2, NINT, kmax, dcf);
        } else {
          rtzlen = -1;
        }

        if (rtzlen > 0) {
          Gc_bot = Gc;
          rgxt = -gxt;
          rgyt = -gyt;
          rgzt = -gzt;
        } else {
          Gc_top = Gc;
        }
      }
      if (Gc_bot > 0) {
        rtzlen =
            rtz(rgxt, rgyt, rgzt, kmaxg * GAMMA_PROTON * GSP,
                0 * GAMMA_PROTON * GSP, kmaxg / tan(theta) * GAMMA_PROTON * GSP,
                theta, numpt, Gmax_xy, Gmax_z, Smax_xy, Smax_z, GSP, prgx, prgy,
                prgz, prkx, prky, prkz, FOV1, FOV2, NINT, kmax, dcf);
      } else {
        rtzlen = -1;
      }
      ai++;

    } else { // if ((rangehigh == dlen-2) && (rtzlen>0))
      WarnMeIfIGoHere = 0;
      WarnMeIfIGoHereXXX = 0;
      while ((rangehigh - rangelow) > 0) {
        ai = (int)ceil((rangehigh - rangelow) / 2.0) + rangelow;
        rtzlen = rtz(-pgx[ai + 1], -pgy[ai + 1], -pgz[ai + 1], pkx[ai], pky[ai],
                     pkz[ai], theta, numpt, Gmax_xy, Gmax_z, Smax_xy, Smax_z,
                     GSP, prgx, prgy, prgz, prkx, prky, prkz, FOV1, FOV2, NINT,
                     kmax, dcf);
        if (rtzlen == -1) {
          WarnMeIfIGoHere = WarnMeIfIGoHere + 1;
          rangehigh = ai - 1;
        } else {
          WarnMeIfIGoHereXXX = WarnMeIfIGoHereXXX + 1;
          rangelow = ai;
        }
      } /* End while (rangehigh-rangelow)>0 */

      ai = rangelow;
      rtzlen =
          rtz(-pgx[ai + 1], -pgy[ai + 1], -pgz[ai + 1], pkx[ai], pky[ai],
              pkz[ai], theta, numpt, Gmax_xy, Gmax_z, Smax_xy, Smax_z, GSP,
              prgx, prgy, prgz, prkx, prky, prkz, FOV1, FOV2, NINT, kmax, dcf);
    }

    // rtzlen = -1; //If I force this -> Fig.3.7a
    // rtzlen is the location where rto starts

    /* COPY OVER SPOKE TRAJECTORY */
    if (rtzlen > 0) {
      Smax_c = min(Smax_z / cos(theta), Smax_xy / sin(theta));
      Gmax_c = min(Gmax_z / cos(theta), Gmax_xy / sin(theta));

      gxyend = sqrt(prgx[rtzlen] * prgx[rtzlen] + prgy[rtzlen] * prgy[rtzlen]);
      kxyend = sqrt(prkx[rtzlen] * prkx[rtzlen] + prky[rtzlen] * prky[rtzlen]);
      gxyzend = sqrt(prgx[rtzlen] * prgx[rtzlen] + prgy[rtzlen] * prgy[rtzlen] +
                     prgz[rtzlen] * prgz[rtzlen]);
      kxyzend = sqrt(prkx[rtzlen] * prkx[rtzlen] + prky[rtzlen] * prky[rtzlen] +
                     prkz[rtzlen] * prkz[rtzlen]);
      rtolen = rto(gxyzend, kxyzend * GAMMA_PROTON * GSP, Smax_c, Gmax_c,
                   GAMMA_PROTON, GSP, pgfx, pkfx);
      for (bi = 0; bi < rtolen; bi++) {
        pgfz[bi] = pgfx[bi] * cos(theta);
        pgfy[bi] = pgfx[bi] * sin(theta) * prky[rtzlen] / kxyend;
        pgfx[bi] = pgfx[bi] * sin(theta) * prkx[rtzlen] / kxyend;
        pkfx[bi] =
            ((bi > 0) ? pkfx[bi - 1] : 0) + pgfx[bi] * GAMMA_PROTON * GSP;
        pkfy[bi] =
            ((bi > 0) ? pkfy[bi - 1] : 0) + pgfy[bi] * GAMMA_PROTON * GSP;
        pkfz[bi] =
            ((bi > 0) ? pkfz[bi - 1] : 0) + pgfz[bi] * GAMMA_PROTON * GSP;
      }

      /* COPY OVER CIRCULAR TRAJECTORY */
      for (bi = 0; bi < min(rtzlen, numpt - rtolen); bi++) {
        pgfx[bi + rtolen] = -prgx[rtzlen - bi];
        pgfy[bi + rtolen] = -prgy[rtzlen - bi];
        pgfz[bi + rtolen] = -prgz[rtzlen - bi];
        pkfx[bi + rtolen] = prkx[rtzlen - bi - 1] * GAMMA_PROTON * GSP;
        pkfy[bi + rtolen] = prky[rtzlen - bi - 1] * GAMMA_PROTON * GSP;
        pkfz[bi + rtolen] = prkz[rtzlen - bi - 1] * GAMMA_PROTON * GSP;
      }
      /* COPY OVER DIFFERENTIAL TRAJECTORY */
      for (bi = ai + 1; bi < min(dlen, numpt - rtolen - rtzlen + (ai + 1));
           bi++) {
        pgfx[bi - (ai + 1) + rtolen + rtzlen] = pgx[bi];
        pgfy[bi - (ai + 1) + rtolen + rtzlen] = pgy[bi];
        pgfz[bi - (ai + 1) + rtolen + rtzlen] = pgz[bi];
        pkfx[bi - (ai + 1) + rtolen + rtzlen] = pkx[bi];
        pkfy[bi - (ai + 1) + rtolen + rtzlen] = pky[bi];
        pkfz[bi - (ai + 1) + rtolen + rtzlen] = pkz[bi];
      }

      *pai = min(numpt, dlen - ai + rtolen + rtzlen) - 1;

    } else { /*if (rtzlen>0) */

      Smax_c = min(Smax_z / cos(theta), Smax_xy / sin(theta));
      Gmax_c = min(Gmax_z / cos(theta), Gmax_xy / sin(theta));
      kxyend = sqrt(pkx[0] * pkx[0] + pky[0] * pky[0]);
      kxyzend = sqrt(pkx[0] * pkx[0] + pky[0] * pky[0] + pkz[0] * pkz[0]);
      rtolen = rto(0, kxyzend, Smax_c, Gmax_c, GAMMA_PROTON, GSP, pgfx, pkfx);
      for (bi = 0; bi < min(rtolen, numpt); bi++) {
        pgfz[bi] = pgfx[bi] * cos(theta);
        pgfy[bi] = pgfx[bi] * sin(theta) * pky[0] / kxyend;
        pgfx[bi] = pgfx[bi] * sin(theta) * pkx[0] / kxyend;
        pkfx[bi] =
            ((bi > 0) ? pkfx[bi - 1] : 0) + pgfx[bi] * GAMMA_PROTON * GSP;
        pkfy[bi] =
            ((bi > 0) ? pkfy[bi - 1] : 0) + pgfy[bi] * GAMMA_PROTON * GSP;
        pkfz[bi] =
            ((bi > 0) ? pkfz[bi - 1] : 0) + pgfz[bi] * GAMMA_PROTON * GSP;
      }

      for (bi = 0; bi < min(dlen, numpt - rtolen - 1); bi++) {
        pgfx[rtolen + bi] = pgx[bi];
        pgfy[rtolen + bi] = pgy[bi];
        pgfz[rtolen + bi] = pgz[bi];
        pkfx[rtolen + bi] = pkx[bi];
        pkfy[rtolen + bi] = pky[bi];
        pkfz[rtolen + bi] = pkz[bi];
      }
      *pai = min(dlen, numpt - rtolen - 1) + min(rtolen, numpt);
    } // end if (rtzlen>0)

  } /* kstart>(kmax-0.0001) || acq_mode==2 (If Twist needed)*/

  atmp = (int)(*pai - 1);

  /* Free up allocated memory */
  FreeMem(pkr);
  FreeMem(pcsndk);
  FreeMem(prgx);
  FreeMem(prgy);
  FreeMem(prgz);
  FreeMem(prkx);
  FreeMem(prky);
  FreeMem(prkz);
  FreeMem(pgx);
  FreeMem(pgy);
  FreeMem(pgz);
  FreeMem(pkx);
  FreeMem(pky);
  FreeMem(pkz);

  return 1;
}

/* Function:  getthetas
 *
 * Calculates the necessary polar angles for cones needed to fill a sphere of
 * k-space for a given resolution and FOV (first element of RES and FOV are for
 * x/y resolution/fov and second element is for z resolution/fov
 */
float *getthetas(float *RES, float *FOV, int *ntheta, float *dtheta_ic) {
  int ai = 0;
  float *theta;
  float thetap[20 * MAX_THETAS]; /* Force fixed sized array of max reasonable
                                    number of thetas ... */
  float kmaxes[20 * MAX_THETAS];
  float KMAX[2];
  float curKMAX, curFOV;
  float dthetap;
  int i;

  thetap[ai] = 0;
  KMAX[0] = 5 / RES[0];
  KMAX[1] = 5 / RES[1];

  while (thetap[ai] <= PI / 2) {
    curFOV = (FOV[0] * FOV[1] / sqrt((FOV[1] * FOV[1] - FOV[0] * FOV[0]) *
                                         sin(thetap[ai]) * sin(thetap[ai]) +
                                     FOV[0] * FOV[0]));
    curKMAX =
        (KMAX[1] * KMAX[0] / sqrt((KMAX[0] * KMAX[0] - KMAX[1] * KMAX[1]) *
                                      sin(thetap[ai]) * sin(thetap[ai]) +
                                  KMAX[1] * KMAX[1]));
    dthetap = 1 / curFOV / curKMAX;
    thetap[ai + 1] = thetap[ai] + dthetap;
    ai++;
  } /* End while thetap */

  theta = (float *)AllocMem((ai + 1) * sizeof(float)); /* ML - not freed */
  *ntheta = ai + 1;

  for (i = 0; i <= ai; i++) {
    thetap[i] = thetap[i] / thetap[ai] * PI / 2;
    kmaxes[i] =
        (KMAX[1] * KMAX[0] / sqrt((KMAX[0] * KMAX[0] - KMAX[1] * KMAX[1]) *
                                      sin(thetap[i]) * sin(thetap[i]) +
                                  KMAX[1] * KMAX[1]));
    theta[i] = atan2f(kmaxes[i] * sin(thetap[i]) / KMAX[1],
                      kmaxes[i] * cos(thetap[i]) / KMAX[0]);
    /* Avoid divide by zero errors */
    theta[i] = min(theta[i], PI / 2 - 1e-5);
    theta[i] = max(theta[i], 1e-5);
  }

  dtheta_ic[0] = thetap[1] - thetap[0];
  for (i = 1; i <= ai - 1; i++) {
    dtheta_ic[i] = 0.5 * (thetap[i + 1] - thetap[i - 1]);
  }
  dtheta_ic[ai] = thetap[ai] - thetap[ai - 1];

  return theta;
}

/* Function:  getthetas2
 *
 * Calculates the necessary polar angles for cones needed to fill a sphere of
 * k-space for a given resolution and FOV (first element of RES and FOV are for
 * x/y resolution/fov and second element is for z resolution/fov
 */
void getthetas2(float *theta, float RES[2], float FOV[2], int *ntheta,
                float *dtheta_ic) {
  int ai = 0;
  float thetap[20 * MAX_THETAS]; /* Force fixed sized array of max reasonable
                                    number of thetas ... */
  float kmaxes[20 * MAX_THETAS];
  float KMAX[2];
  float curKMAX, curFOV;
  float dthetap;
  int i;

  thetap[ai] = 0;
  KMAX[0] = 5 / RES[0];
  KMAX[1] = 5 / RES[1];

  /*printf("MCARL: MCARL: MCARL: FOV[0] = %f\n",FOV[0]);
  printf("MCARL: MCARL: MCARL: FOV[1] = %f\n",FOV[1]);
  fflush(stdout);*/

  while (thetap[ai] <= PI / 2) {
    curFOV = (FOV[0] * FOV[1] / sqrt((FOV[1] * FOV[1] - FOV[0] * FOV[0]) *
                                         sin(thetap[ai]) * sin(thetap[ai]) +
                                     FOV[0] * FOV[0]));
    curKMAX =
        (KMAX[1] * KMAX[0] / sqrt((KMAX[0] * KMAX[0] - KMAX[1] * KMAX[1]) *
                                      sin(thetap[ai]) * sin(thetap[ai]) +
                                  KMAX[1] * KMAX[1]));
    dthetap = 1 / curFOV / curKMAX;
    thetap[ai + 1] = thetap[ai] + dthetap;
    ai++;
  } /* End while thetap */

  /* Python wrapper: removed the in-place realloc that lived here. It tried to
     shrink the caller's `theta` buffer to (ai+1) floats, but realloc may move
     the allocation -- and the original author only updated the *local* `theta`
     parameter, leaving the caller's pointer dangling (use-after-free). The
     caller in conegrad/conegrad_vd already allocates 20*nthetatmp floats which
     is far more than `ai+1`, so the shrink wasn't needed for correctness. */
  *ntheta = ai + 1;

  for (i = 0; i <= ai; i++) {
    thetap[i] = thetap[i] / thetap[ai] * PI / 2;
    kmaxes[i] =
        (KMAX[1] * KMAX[0] / sqrt((KMAX[0] * KMAX[0] - KMAX[1] * KMAX[1]) *
                                      sin(thetap[i]) * sin(thetap[i]) +
                                  KMAX[1] * KMAX[1]));
    theta[i] = atan2f(kmaxes[i] * sin(thetap[i]) / KMAX[1],
                      kmaxes[i] * cos(thetap[i]) / KMAX[0]);
    /* Avoid divide by zero errors */
    theta[i] = min(theta[i], PI / 2 - 1e-5);
    theta[i] = max(theta[i], 1e-5);
  }

  /*for (i=0; i<=ai; i++) {
  printf("MCARL: MCARL: MCARL: thetap[%d] = %f\n",i,180/PI*thetap[i]);
  fflush(stdout);
  }*/

  dtheta_ic[0] = thetap[1] - thetap[0];
  for (i = 1; i <= ai - 1; i++) {
    dtheta_ic[i] = 0.5 * (thetap[i + 1] - thetap[i - 1]);
  }
  dtheta_ic[ai] = thetap[ai] - thetap[ai - 1];

  return;
}

/* VD Start */

int gencone_vd(float ***g, int *len, int *totlen, float *RES, float *FOV,
               float NINT, float THETA[2], int MAXLEN, float GSP, float SMAX,
               float maxSRewind, float GMAX, float SysMaxRewindG,
               int OVERSAMPLE, float DCF, float MINDENS, int nramp,
               int rewind_flag, int acq_mode) {
  float THETArange, KMAXrange;
  float maxtheta, mintheta, truemaxtheta, truemintheta;
  int i, j, c, ntmp;
  float ftmp, theta;
  float *pn = NULL;
  float ai;
  int ret;
  float **k;
  float **gtmp;
  float kmax_xy, kmax_z;
  float tparam;
  float curFOV, curKMAX, dthetaCUR, dthetaBASE, densadjust;
  float xyscale, zscale;
  int glen;
  int nrampdown, nconst, nrampdowna;
  float Gstart;
  float am;
  float Kstart;
  float Garea;

  /* VD */

  float Nkr = 100.0;
  float FOVcirc_vd[100];
  float FOVrad_vd[100];
  float FOVkr_vd[100];
  float FOVrad_maxtheta_vd[100];
  float FOVrad_mintheta_vd[100];
  int idx = 0;

  /* VD */

  k = (float **)AllocMem(3 * sizeof(float *));
  k[0] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));
  k[1] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));
  k[2] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));

  gtmp = (float **)AllocMem(3 * sizeof(float *));
  gtmp[0] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));
  gtmp[1] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));
  gtmp[2] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));

  ftmp = min(fabs(THETA[0]), fabs(THETA[1]));
  theta = min(ftmp, PI / 2);
  maxtheta = max(fabs(THETA[0]), fabs(THETA[1]));
  mintheta = min(fabs(THETA[0]), fabs(THETA[1]));
  maxtheta = min(maxtheta, PI / 2 - 1e-5);
  maxtheta = max(maxtheta, 1e-5);
  mintheta = min(mintheta, PI / 2 - 1e-5);
  mintheta = max(mintheta, 1e-5);
  kmax_xy = 5.0 / RES[0];
  kmax_z = 5.0 / RES[1];

  /* Choose an appropriate angle in between the two thetas */
  THETArange = atan2f(kmax_z * sin(maxtheta), kmax_xy * cos(mintheta));
  /* Find the kmax at that angle */
  KMAXrange = intlineellipse(kmax_xy, kmax_z, THETArange);
  /* Determine the parametric theta for the chosen angle */
  tparam = atan2f((KMAXrange * sin(THETArange) * kmax_xy),
                  (KMAXrange * cos(THETArange) * kmax_z));

  /* Determine the amount by which the waveforms will be scaled to cover the
   * range of thetas */
  xyscale = cos(tparam) / cos(mintheta);
  zscale = sin(tparam) / sin(maxtheta);

  /* VD */
  /* The circumferential FOV is always just FOV_xy */
  /* FOVcirc = FOV[0]; */

  for (idx = 0; idx < ((int)Nkr); idx++) {
    FOVcirc_vd[idx] = FOV[idx];
  }

  /* VD */

  /* The radial FOV is the FOV in the radial direction */
  truemaxtheta = atan2f(sin(maxtheta) * kmax_z, cos(maxtheta) * kmax_xy);
  /*truemintheta = atan2f(sin(mintheta)*RES[1],cos(mintheta)*kmax_xy); */ /*MCARL:
   * This looks wrong*/
  truemintheta =
      atan2f(sin(mintheta) * kmax_z, cos(mintheta) * kmax_xy); /*MCARL: This?*/

  /* VD */

  /* FOVrad_maxtheta = 1/intlineellipse(1/FOV[0],1/FOV[1],truemaxtheta); */
  /* FOVrad_mintheta = 1/intlineellipse(1/FOV[0],1/FOV[1],truemintheta); */
  for (idx = 0; idx < ((int)Nkr); idx++) {
    FOVrad_maxtheta_vd[idx] =
        1 /
        intlineellipse(1 / FOV[idx], 1 / FOV[idx + ((int)Nkr)], truemaxtheta);
    FOVrad_mintheta_vd[idx] =
        1 /
        intlineellipse(1 / FOV[idx], 1 / FOV[idx + ((int)Nkr)], truemintheta);
    FOVrad_vd[idx] = max(FOVrad_maxtheta_vd[idx], FOVrad_mintheta_vd[idx]);
  }
  /*printf("MCARL: truemaxtheta=%f\n",truemaxtheta);
  printf("MCARL: truemintheta=%f\n",truemintheta);
  printf("MCARL: FOVrad_maxtheta=%f\n",FOVrad_maxtheta);
  printf("MCARL: FOVrad_mintheta=%f\n",FOVrad_mintheta);
  fflush(stdout);*/

  /* FOVrad = max(FOVrad_maxtheta, FOVrad_mintheta); */

  /* Figure out the dtheta relationship (intracones) to adjust the density */
  curFOV = (FOV[0] * FOV[((int)Nkr)] /
            sqrt((FOV[((int)Nkr)] * FOV[((int)Nkr)] - FOV[0] * FOV[0]) *
                     sin(tparam) * sin(tparam) +
                 FOV[0] * FOV[0]));
  curKMAX = (kmax_z * kmax_xy / sqrt((kmax_xy * kmax_xy - kmax_z * kmax_z) *
                                         sin(tparam) * sin(tparam) +
                                     kmax_z * kmax_z));
  dthetaCUR = 1 / curFOV / curKMAX;
  dthetaBASE = 1 / FOV[((int)Nkr)] / kmax_xy;
  densadjust = dthetaCUR / dthetaBASE;
  for (idx = 0; idx < ((int)Nkr); idx++) {
    FOVkr_vd[idx] = FOV[idx + ((int)Nkr) * 2] * curKMAX;
  }

  /* VD */

  GSP = GSP / OVERSAMPLE;

  /* To Do: get second FOVcirc/FOVrad configured properly for asymmetric fov */

  /* VD */

  ret = whirlcone_vd(PI / 2 - THETArange, FOVcirc_vd, FOVrad_vd, FOVkr_vd,
                     FOVrad_vd, DCF, KMAXrange, NINT, MAXLEN * OVERSAMPLE,
                     (float)GSP, SMAX * GSP, GMAX, (float)SMAX * GSP * xyscale,
                     (float)SMAX * GSP * zscale, (float)GMAX * xyscale,
                     (float)GMAX * zscale, MINDENS * densadjust, &gtmp, &k, &ai,
                     pn, acq_mode);

  /* VD */
  *len = (int)ai;
  *totlen = (int)ai;
  glen = min((*len), (MAXLEN - 5 -
                      1)); /* Last index holding designed gradient to use */
  ntmp = MAXLEN - 5 + nramp - glen;

  for (j = 0; j < 3; j++) {
    Garea = 0;
    /* Copy designed gradient to output array */
    for (i = 0; i < glen; i++) {
      (*g)[j][i] = gtmp[j][i];
      Garea = Garea + gtmp[j][i] * GSP * 1e6;
    }

    if (rewind_flag == 1) {
      /* Determine parameters of rewinder gradient lobe */
      Gstart = (*g)[j][glen - 1];
      Kstart = k[j][glen - 1] * 1e6 / GAMMA_PROTON;

      float TestAmp, Grew;
      int TRamp = ceil(SysMaxRewindG * 1e6 / maxSRewind);
      TRamp = (int)ceil(TRamp / GSP) * GSP;
      float Slew = fabs(SysMaxRewindG) / (float)TRamp;
      int T2Zero = ceil((fabs(Gstart) / Slew) / (GSP * 1e6)) * (GSP * 1e6);
      nrampdowna = ceil(T2Zero / (GSP * 1e6));

      float K2Zero;
      // K2Zero=Slew*T2Zero*T2Zero/2; //Initial calculation
      // if (Gstart<0) K2Zero=-K2Zero;

      // Create ramp from end of trajectory to zero
      Garea = 0;
      for (c = 1; c <= nrampdowna; c++) {
        (*g)[j][glen - 1 + c] =
            Gstart * ((float)(nrampdowna - c) / (float)nrampdowna);
        Garea = Garea + (*g)[j][glen - 1 + c] * GSP * 1e6;
      }
      K2Zero = Garea;

      float Krew = Kstart + K2Zero; // Total k-sp extend after G is back to
                                    // zero.

      float Kramp =
          Slew * TRamp * TRamp; // Max possible rewinder using a triangle.

      float DecideSign = 1;
      if (Krew >= 0)
        DecideSign = -1;
      if (Krew < 0)
        DecideSign = 1;

      float Kflat, Kdesign;

      int Tflat;
      if (fabs(Krew) > Kramp) {
        Tflat = ceil(((fabs(Krew) - Kramp) / SysMaxRewindG) / (GSP * 1e6)) *
                (GSP * 1e6);
        Kflat = (Tflat)*SysMaxRewindG;
        Kdesign = Kflat + Kramp;
        TestAmp = DecideSign * SysMaxRewindG;
      } else {
        TRamp = ceil(sqrt(fabs(Krew) / Slew) / (GSP * 1e6)) * (GSP * 1e6);
        Tflat = GSP;
        Kflat = (Tflat)*SysMaxRewindG;
        Kramp = Slew * TRamp * TRamp;
        Kdesign = Kflat + Kramp;
        TestAmp = DecideSign * Slew * TRamp;
      }

      nconst = ceil(Tflat / (GSP * 1e6));
      nrampdown = ceil(TRamp / (GSP * 1e6));

      // Test actual area
      Kdesign = 0;
      for (c = 1; c <= nrampdown; c++) {
        Grew = TestAmp * ((float)c / (float)nrampdown);
        Kdesign = Kdesign + Grew * GSP * 1e6;
      }
      for (c = 1; c <= nconst; c++) {
        Grew = TestAmp;
        Kdesign = Kdesign + Grew * GSP * 1e6;
      }
      for (c = 1; c <= nrampdown; c++) {
        Grew = TestAmp * ((float)(nrampdown - c) / (float)nrampdown);
        Kdesign = Kdesign + Grew * GSP * 1e6;
      }

      am = TestAmp * fabs(Krew / Kdesign);

      // Create ramp from zero to rewinder amplitude
      Garea = 0;
      for (c = 1; c <= nrampdown; c++) {
        (*g)[j][glen - 1 + nrampdowna + c] =
            (am) * ((float)c / (float)nrampdown);
        Garea = Garea + (*g)[j][glen - 1 + nrampdowna + c] * GSP * 1e6;
      };
      // Make constant portion of rewinder
      Garea = 0;
      for (c = 1; c <= nconst; c++) {
        (*g)[j][glen - 1 + nrampdowna + nrampdown + c] = am;
        Garea =
            Garea + (*g)[j][glen - 1 + nrampdowna + nrampdown + c] * GSP * 1e6;
      }
      // Make ramp back down to zero
      Garea = 0;
      for (c = 1; c <= nrampdown; c++) {
        (*g)[j][glen - 1 + nrampdowna + nrampdown + nconst + c] =
            (am) * ((float)(nrampdown - c) / (float)nrampdown);
        Garea =
            Garea +
            (*g)[j][glen - 1 + nrampdowna + nrampdown + nconst + c] * GSP * 1e6;
      }
      // Make sure test of waveform is zero
      Garea = 0;
      for (c = 1;
           c <= MAXLEN - 5 + nramp - glen - nrampdowna - nconst - 2 * nrampdown;
           c++) {
        (*g)[j][glen - 1 + nrampdowna + nconst + 2 * nrampdown + c] = 0.0;
        Garea = Garea +
                (*g)[j][glen - 1 + nrampdowna + nconst + 2 * nrampdown + c] *
                    GSP * 1e6;
      }

      if ((*totlen) < (ai + nrampdowna + nconst + 2 * nrampdown))
        (*totlen) = (int)(ai + nrampdowna + nconst + 2 * nrampdown);

    } else {
      /* No rewinder */

      /* Old ramp down to 0 */
      ntmp =
          (int)ceil(fabs((*g)[j][glen - 1]) * 1e6 / (maxSRewind * GSP * 1e6));

      for (c = 1; c <= ntmp; c++)
        (*g)[j][glen + c - 1] =
            (*g)[j][glen - 1] * ((float)(ntmp - (c - 1)) / (float)ntmp);

      /* Make sure rest of waveform is zero */
      for (c = 1; c <= MAXLEN - 5 + nramp - glen - ntmp; c++)
        (*g)[j][glen - 1 + ntmp + c] = 0.0;

      if ((*totlen) < (ai + ntmp))
        (*totlen) = (int)(ai + ntmp);

    } /* end if rewind_flag */

    /*printf("ai = %f\n",ai);
    printf("nramp = %d\n",nramp);
    printf("ntmp = %f\n",ntmp);
    fflush(stdout);*/

  } /* End for j (x/y/z index) */

  // printf("MCARL: totlen=%d\n",*totlen);
  // fflush(stdout);

  FreeMem(k[2]);
  FreeMem(k[1]);
  FreeMem(k[0]);
  FreeMem(k);

  FreeMem(gtmp[2]);
  FreeMem(gtmp[1]);
  FreeMem(gtmp[0]);
  FreeMem(gtmp);

  return 1;
}

/* VD End */

/* Function:  gencone
 *
 * Function responsible for handling the reuse of cones for a range of thetas
 * and
 * calling the low-level cones design function (whirlcone) appropriately.
 * NOTE: Basically assuming always calling with OVERSAMPLE=1 now.
 */
int gencone(float ***g, int *len, int *totlen, float *RES, float *FOV,
            float NINT, float THETA[2], int MAXLEN, float GSP, float SMAX,
            float maxSRewind, float GMAX, float SysMaxRewindG, int OVERSAMPLE,
            float DCF, float MINDENS, int nramp, int rewind_flag,
            int acq_mode) {
  float THETArange, KMAXrange;
  float FOVcirc, FOVrad, FOVrad_maxtheta, FOVrad_mintheta;
  float maxtheta, mintheta, truemaxtheta, truemintheta;
  int i, j, c, ntmp;
  float ftmp, theta;
  float *pn = NULL;
  float ai;
  int ret;
  float **k;
  float **gtmp;
  float kmax_xy, kmax_z;
  float tparam;
  float curFOV, curKMAX, dthetaCUR, dthetaBASE, densadjust;
  float xyscale, zscale;
  int glen;
  int nrampdown, nconst, nrampdowna;
  float Gstart;
  float am;
  float Kstart;
  float Garea;

  k = (float **)AllocMem(3 * sizeof(float *));
  k[0] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));
  k[1] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));
  k[2] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));

  gtmp = (float **)AllocMem(3 * sizeof(float *));
  gtmp[0] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));
  gtmp[1] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));
  gtmp[2] = (float *)AllocMem(MAXLEN * OVERSAMPLE * sizeof(float));

  ftmp = min(fabs(THETA[0]), fabs(THETA[1]));
  theta = min(ftmp, PI / 2);
  maxtheta = max(fabs(THETA[0]), fabs(THETA[1]));
  mintheta = min(fabs(THETA[0]), fabs(THETA[1]));
  maxtheta = min(maxtheta, PI / 2 - 1e-5);
  maxtheta = max(maxtheta, 1e-5);
  mintheta = min(mintheta, PI / 2 - 1e-5);
  mintheta = max(mintheta, 1e-5);
  kmax_xy = 5.0 / RES[0];
  kmax_z = 5.0 / RES[1];

  /* Choose an appropriate angle in between the two thetas */
  THETArange = atan2f(kmax_z * sin(maxtheta), kmax_xy * cos(mintheta));
  /* Find the kmax at that angle */
  KMAXrange = intlineellipse(kmax_xy, kmax_z, THETArange);
  /* Determine the parametric theta for the chosen angle */
  tparam = atan2f((KMAXrange * sin(THETArange) * kmax_xy),
                  (KMAXrange * cos(THETArange) * kmax_z));

  /* Determine the amount by which the waveforms will be scaled to cover the
   * range of thetas */
  xyscale = cos(tparam) / cos(mintheta);
  zscale = sin(tparam) / sin(maxtheta);

  /* The circumferential FOV is always just FOV_xy */
  FOVcirc = FOV[0];

  /* The radial FOV is the FOV in the radial direction */
  truemaxtheta = atan2f(sin(maxtheta) * kmax_z, cos(maxtheta) * kmax_xy);
  /*truemintheta = atan2f(sin(mintheta)*RES[1],cos(mintheta)*kmax_xy); */ /*MCARL:
   * This looks wrong */
  truemintheta =
      atan2f(sin(mintheta) * kmax_z, cos(mintheta) * kmax_xy); /*MCARL: This?*/

  FOVrad_maxtheta = 1 / intlineellipse(1 / FOV[0], 1 / FOV[1], truemaxtheta);
  FOVrad_mintheta = 1 / intlineellipse(1 / FOV[0], 1 / FOV[1], truemintheta);

  /*printf("MCARL: truemaxtheta=%f\n",truemaxtheta);
  printf("MCARL: truemintheta=%f\n",truemintheta);
  printf("MCARL: FOVrad_maxtheta=%f\n",FOVrad_maxtheta);
  printf("MCARL: FOVrad_mintheta=%f\n",FOVrad_mintheta);
  fflush(stdout);*/

  FOVrad = max(FOVrad_maxtheta, FOVrad_mintheta);

  /* Figure out the dtheta relationship (intracones) to adjust the density */
  curFOV = (FOV[0] * FOV[1] / sqrt((FOV[1] * FOV[1] - FOV[0] * FOV[0]) *
                                       sin(tparam) * sin(tparam) +
                                   FOV[0] * FOV[0]));
  curKMAX = (kmax_z * kmax_xy / sqrt((kmax_xy * kmax_xy - kmax_z * kmax_z) *
                                         sin(tparam) * sin(tparam) +
                                     kmax_z * kmax_z));
  dthetaCUR = 1 / curFOV / curKMAX;
  dthetaBASE = 1 / FOV[1] / kmax_xy;
  densadjust = dthetaCUR / dthetaBASE;

  GSP = GSP / OVERSAMPLE;

  /* To Do: get second FOVcirc/FOVrad configured properly for asymmetric fov */
  ret = whirlcone(PI / 2 - THETArange, FOVcirc, FOVrad, FOVcirc, FOVrad, DCF,
                  KMAXrange, NINT, MAXLEN * OVERSAMPLE, (float)GSP, SMAX * GSP,
                  GMAX, (float)SMAX * GSP * xyscale, (float)SMAX * GSP * zscale,
                  (float)GMAX * xyscale, (float)GMAX * zscale,
                  MINDENS * densadjust, &gtmp, &k, &ai, pn, acq_mode);

  *len = (int)ai;
  *totlen = (int)ai;
  glen = min((*len), (MAXLEN - 5 -
                      1)); /* Last index holding designed gradient to use */
  ntmp = MAXLEN - 5 + nramp - glen;

  for (j = 0; j < 3; j++) {
    Garea = 0;
    /* Copy designed gradient to output array */
    for (i = 0; i < glen; i++) {
      (*g)[j][i] = gtmp[j][i];
      Garea = Garea + gtmp[j][i] * GSP * 1e6;
    }

    if (rewind_flag == 1) {
      /* Determine parameters of rewinder gradient lobe */
      Gstart = (*g)[j][glen - 1];
      Kstart = k[j][glen - 1] * 1e6 / GAMMA_PROTON;

      float TestAmp, Grew;
      int TRamp = ceil(SysMaxRewindG * 1e6 / maxSRewind);
      TRamp = (int)ceil(TRamp / GSP) * GSP;
      float Slew = fabs(SysMaxRewindG) / (float)TRamp;
      int T2Zero = ceil((fabs(Gstart) / Slew) / (GSP * 1e6)) * (GSP * 1e6);
      nrampdowna = ceil(T2Zero / (GSP * 1e6));

      float K2Zero;
      // K2Zero=Slew*T2Zero*T2Zero/2; //Initial calculation
      // if (Gstart<0) K2Zero=-K2Zero;

      // Create ramp from end of trajectory to zero
      Garea = 0;
      for (c = 1; c <= nrampdowna; c++) {
        (*g)[j][glen - 1 + c] =
            Gstart * ((float)(nrampdowna - c) / (float)nrampdowna);
        Garea = Garea + (*g)[j][glen - 1 + c] * GSP * 1e6;
      }
      K2Zero = Garea;

      float Krew = Kstart + K2Zero; // Total k-sp extend after G is back to
                                    // zero.

      float Kramp =
          Slew * TRamp * TRamp; // Max possible rewinder using a triangle.

      float DecideSign = 1;
      if (Krew >= 0)
        DecideSign = -1;
      if (Krew < 0)
        DecideSign = 1;

      float Kflat, Kdesign;

      int Tflat;
      if (fabs(Krew) > Kramp) {
        Tflat = ceil(((fabs(Krew) - Kramp) / SysMaxRewindG) / (GSP * 1e6)) *
                (GSP * 1e6);
        Kflat = (Tflat)*SysMaxRewindG;
        Kdesign = Kflat + Kramp;
        TestAmp = DecideSign * SysMaxRewindG;
      } else {
        TRamp = ceil(sqrt(fabs(Krew) / Slew) / (GSP * 1e6)) * (GSP * 1e6);
        Tflat = GSP;
        Kflat = (Tflat)*SysMaxRewindG;
        Kramp = Slew * TRamp * TRamp;
        Kdesign = Kflat + Kramp;
        TestAmp = DecideSign * Slew * TRamp;
      }

      nconst = ceil(Tflat / (GSP * 1e6));
      nrampdown = ceil(TRamp / (GSP * 1e6));

      // Test actual area
      Kdesign = 0;
      for (c = 1; c <= nrampdown; c++) {
        Grew = TestAmp * ((float)c / (float)nrampdown);
        Kdesign = Kdesign + Grew * GSP * 1e6;
      }
      for (c = 1; c <= nconst; c++) {
        Grew = TestAmp;
        Kdesign = Kdesign + Grew * GSP * 1e6;
      }
      for (c = 1; c <= nrampdown; c++) {
        Grew = TestAmp * ((float)(nrampdown - c) / (float)nrampdown);
        Kdesign = Kdesign + Grew * GSP * 1e6;
      }

      am = TestAmp * fabs(Krew / Kdesign);

      // Create ramp from zero to rewinder amplitude
      Garea = 0;
      for (c = 1; c <= nrampdown; c++) {
        (*g)[j][glen - 1 + nrampdowna + c] =
            (am) * ((float)c / (float)nrampdown);
        Garea = Garea + (*g)[j][glen - 1 + nrampdowna + c] * GSP * 1e6;
      };
      // Make constant portion of rewinder
      Garea = 0;
      for (c = 1; c <= nconst; c++) {
        (*g)[j][glen - 1 + nrampdowna + nrampdown + c] = am;
        Garea =
            Garea + (*g)[j][glen - 1 + nrampdowna + nrampdown + c] * GSP * 1e6;
      }
      // Make ramp back down to zero
      Garea = 0;
      for (c = 1; c <= nrampdown; c++) {
        (*g)[j][glen - 1 + nrampdowna + nrampdown + nconst + c] =
            (am) * ((float)(nrampdown - c) / (float)nrampdown);
        Garea =
            Garea +
            (*g)[j][glen - 1 + nrampdowna + nrampdown + nconst + c] * GSP * 1e6;
      }
      // Make sure test of waveform is zero
      Garea = 0;
      for (c = 1;
           c <= MAXLEN - 5 + nramp - glen - nrampdowna - nconst - 2 * nrampdown;
           c++) {
        (*g)[j][glen - 1 + nrampdowna + nconst + 2 * nrampdown + c] = 0.0;
        Garea = Garea +
                (*g)[j][glen - 1 + nrampdowna + nconst + 2 * nrampdown + c] *
                    GSP * 1e6;
      }

      if ((*totlen) < (ai + nrampdowna + nconst + 2 * nrampdown))
        (*totlen) = (int)(ai + nrampdowna + nconst + 2 * nrampdown);

    } else {
      /* No rewinder */

      /* Old ramp down to 0 */
      ntmp =
          (int)ceil(fabs((*g)[j][glen - 1]) * 1e6 / (maxSRewind * GSP * 1e6));

      for (c = 1; c <= ntmp; c++)
        (*g)[j][glen + c - 1] =
            (*g)[j][glen - 1] * ((float)(ntmp - (c - 1)) / (float)ntmp);

      /* Make sure rest of waveform is zero */
      for (c = 1; c <= MAXLEN - 5 + nramp - glen - ntmp; c++)
        (*g)[j][glen - 1 + ntmp + c] = 0.0;

      if ((*totlen) < (ai + ntmp))
        (*totlen) = (int)(ai + ntmp);

    } /* end if rewind_flag */

    /*printf("ai = %f\n",ai);
    printf("nramp = %d\n",nramp);
    printf("ntmp = %f\n",ntmp);
    fflush(stdout);*/

  } /* End for j (x/y/z index) */

  // printf("MCARL: totlen=%d\n",*totlen);
  // fflush(stdout);

  FreeMem(k[2]);
  FreeMem(k[1]);
  FreeMem(k[0]);
  FreeMem(k);

  FreeMem(gtmp[2]);
  FreeMem(gtmp[1]);
  FreeMem(gtmp[0]);
  FreeMem(gtmp);

  return 1;
}

/* VD Start */

int findcone_vd(float ***g, float *nint, int *leng, int *totleng, float *RES,
                float *FOV, int LEN, float THETA[2], float PRECISION, float GSP,
                int OS, float SMAX, float maxSRewind, float GMAX,
                float SysMaxRewindG, float DCF, float MINDENS, int nramp,
                int rewind_flag, int acq_mode) {
  float NINTlo, NINThi;
  float LENlo, LENhi;
  float curNINT;
  int MAXLEN = 10000; /*MCARL: Not needed?*/
  int loopcnt;

  /* Initial conditions for number of interleaves used in search */
  NINTlo = min(PRECISION, 0.1);
  NINThi = 2; /*MCARL: was 100*/ /* Sri: Was 2 */
  LENlo = 10000000;
  loopcnt = 0;

  curNINT = NINThi;
  MAXLEN = LEN + 5; /*LEN is passed from findcone -> GRAD_POINTS*/

  /* Make sure the requested parameters can be achieved with a huge number
     of interleaves - sanity check */
  gencone_vd(g, leng, totleng, RES, FOV, 1e8, THETA, MAXLEN, GSP, SMAX,
             maxSRewind, GMAX, SysMaxRewindG, OS, DCF, MINDENS, nramp,
             rewind_flag, acq_mode);
  if (*leng > LEN) {
    printf("Sorry, but it is impossible to achieve that resolution in that "
           "length of time\n");
    return 0;
  }

  if (acq_mode == 1) /*Cones*/
  {

    /* Iterate over the cone generation in a binary search pattern (iterating
       over number of interleaves) */
    gencone_vd(g, leng, totleng, RES, FOV, curNINT, THETA, MAXLEN, GSP, SMAX,
               maxSRewind, GMAX, SysMaxRewindG, OS, DCF, MINDENS, nramp,
               rewind_flag, acq_mode);
    LENhi = (float)(*leng);

    /* Keep doubling the number of interleaves until we get a waveform which
       satisfies the requirements */
    while ((*leng > LEN)) {

      NINTlo = NINThi;
      NINThi = NINThi * 2;
      curNINT = NINThi;
      LENhi = (float)(*leng);
      gencone_vd(g, leng, totleng, RES, FOV, curNINT, THETA, MAXLEN, GSP, SMAX,
                 maxSRewind, GMAX, SysMaxRewindG, OS, DCF, MINDENS, nramp,
                 rewind_flag, acq_mode);
      loopcnt++;
      if (loopcnt > MAX_LOOPS)
        break; // Security against indef loops. Defined in scones_design.h
               // (MAX_LOOPS=10000)
    }

    loopcnt = 0;
    /* Do a binary search over the remaining range until we achieve the
       requested precision */
    while ((((NINThi - NINTlo) >= PRECISION) || ((*leng) > LEN)) &&
           (!wcfEqual(LENhi, (float)LEN))) {
      curNINT = (NINThi - NINTlo) / 2 + NINTlo;
      gencone_vd(g, leng, totleng, RES, FOV, curNINT, THETA, MAXLEN, GSP, SMAX,
                 maxSRewind, GMAX, SysMaxRewindG, OS, DCF, MINDENS, nramp,
                 rewind_flag, acq_mode);
      if ((*leng) > LEN) {
        NINTlo = curNINT;
        LENlo = (float)(*leng);
      } else {
        NINThi = curNINT;
        LENhi = (float)(*leng);
      }

      loopcnt++;
      if (loopcnt > MAX_LOOPS)
        break;
    }

    *nint = NINThi;

  } /*if acq_mode==1*/

  if (acq_mode == 2) {
    *nint = 0; /*MCARL: PR. Irrelevant, gets re-calculated in gentraj*/
  }

  gencone_vd(g, leng, totleng, RES, FOV, *nint, THETA, MAXLEN, GSP, SMAX,
             maxSRewind, GMAX, SysMaxRewindG, OS, DCF, MINDENS, nramp,
             rewind_flag, acq_mode);

  return 1;
}

/* VD End */

/* Function:  findcone
 *
 * This function iterates over the cone design function with varying numbers of
 * interleaves to find the cone waveform that most closely matches the requested
 * waveform duration. Uses a binary search.
 */
int findcone(float ***g, float *nint, int *leng, int *totleng, float *RES,
             float *FOV, int LEN, float THETA[2], float PRECISION, float GSP,
             int OS, float SMAX, float maxSRewind, float GMAX,
             float SysMaxRewindG, float DCF, float MINDENS, int nramp,
             int rewind_flag, int acq_mode) {
  float NINTlo, NINThi;
  float LENlo, LENhi;
  float curNINT;
  int MAXLEN = 10000; /*MCARL: Not needed?*/
  int loopcnt;

  /* Initial conditions for number of interleaves used in search */
  NINTlo = min(PRECISION, 0.1);
  NINThi = 2; /*MCARL: was 100*/
  LENlo = 10000000;
  loopcnt = 0;

  curNINT = NINThi;
  MAXLEN = LEN + 5; /*LEN is passed from findcone -> GRAD_POINTS*/

  /* Make sure the requested parameters can be achieved with a huge number
     of interleaves - sanity check */
  gencone(g, leng, totleng, RES, FOV, 1e8, THETA, MAXLEN, GSP, SMAX, maxSRewind,
          GMAX, SysMaxRewindG, OS, DCF, MINDENS, nramp, rewind_flag, acq_mode);
  if (*leng > LEN) {
    printf("Sorry, but it is impossible to achieve that resolution in that "
           "length of time\n");
    return 0;
  }

  if (acq_mode == 1) /*Cones*/
  {

    /* Iterate over the cone generation in a binary search pattern (iterating
       over number of interleaves) */
    gencone(g, leng, totleng, RES, FOV, curNINT, THETA, MAXLEN, GSP, SMAX,
            maxSRewind, GMAX, SysMaxRewindG, OS, DCF, MINDENS, nramp,
            rewind_flag, acq_mode);
    LENhi = (float)(*leng);

    /* Keep doubling the number of interleaves until we get a waveform which
       satisfies the requirements */
    while ((*leng > LEN)) {

      NINTlo = NINThi;
      NINThi = NINThi * 2;
      curNINT = NINThi;
      LENhi = (float)(*leng);
      gencone(g, leng, totleng, RES, FOV, curNINT, THETA, MAXLEN, GSP, SMAX,
              maxSRewind, GMAX, SysMaxRewindG, OS, DCF, MINDENS, nramp,
              rewind_flag, acq_mode);
      loopcnt++;
      if (loopcnt > MAX_LOOPS)
        break; // Security against indef loops. Defined in scones_design.h
               // (MAX_LOOPS=10000)
    }

    loopcnt = 0;
    /* Do a binary search over the remaining range until we achieve the
       requested precision */
    while ((((NINThi - NINTlo) >= PRECISION) || ((*leng) > LEN)) &&
           (!wcfEqual(LENhi, (float)LEN))) {
      curNINT = (NINThi - NINTlo) / 2 + NINTlo;
      gencone(g, leng, totleng, RES, FOV, curNINT, THETA, MAXLEN, GSP, SMAX,
              maxSRewind, GMAX, SysMaxRewindG, OS, DCF, MINDENS, nramp,
              rewind_flag, acq_mode);
      if ((*leng) > LEN) {
        NINTlo = curNINT;
        LENlo = (float)(*leng);
      } else {
        NINThi = curNINT;
        LENhi = (float)(*leng);
      }

      loopcnt++;
      if (loopcnt > MAX_LOOPS)
        break;
    }

    *nint = NINThi;

  } /*if acq_mode==1*/

  if (acq_mode == 2) {
    *nint = 0; /*MCARL: PR. Irrelevant, gets re-calculated in gentraj*/
  }

  gencone(g, leng, totleng, RES, FOV, *nint, THETA, MAXLEN, GSP, SMAX,
          maxSRewind, GMAX, SysMaxRewindG, OS, DCF, MINDENS, nramp, rewind_flag,
          acq_mode);

  return 1;
}

/* Function:  intlineellipse
 *
 * Find the intersection of a line and an ellipse.
 */
float intlineellipse(float a, float b, float phi) {
  float r;
  r = b * a / sqrt((a * a - b * b) * sin(phi) * sin(phi) + b * b);
  return (r);
}

/* Function:  mrewind
 *
 * Find a time optimized rewinder that preserves the zero first moment of the
 * waveform to keep cone waveform applicable to an SSFP implementation.
 *
 * NOTE: Have not been using this, likely needs some additional debugging
 */
int mrewind(float *gx, float *gy, float *k0, float *g0, float gmax, float smax,
            float T) {
  float gain = 20;
  float kacc = 0.0001;
  float gacc;
  float rotang;
  float k[2];
  int n, keepgoing;
  float kret[2];
  float ug[2];
  float gtarg[2];
  float dgvect[2];
  float dg[2];
  float gtmp[3001][2]; /* Temp array for rewinder */
  float krtmp[3001][2];
  float karr[3001][2];
  float **g = NULL;
  int i;

  gacc = 4 / 50000;
  rotang = atan2f(k0[1], k0[0]);
  k[0] = sqrt(k0[0] * k0[0] + k0[1] * k0[1]);
  k[1] = 0;

  gtmp[0][0] = g0[0] * cos(rotang) + g0[1] * sin(rotang);
  gtmp[0][1] = -g0[0] * sin(rotang) + g0[1] * cos(rotang);

  n = 0;
  keepgoing = 1;

  while (keepgoing == 1) {
    kret[0] = k[0] +
              GAMMA_PROTON * 0.5 * gtmp[n][0] *
                  sqrt(gtmp[n][0] * gtmp[n][0] + gtmp[n][1] * gtmp[n][1]) /
                  smax;
    kret[1] = k[1] +
              GAMMA_PROTON * 0.5 * gtmp[n][1] *
                  sqrt(gtmp[n][0] * gtmp[n][0] + gtmp[n][1] * gtmp[n][1]) /
                  smax;

    krtmp[n][0] = kret[0];
    krtmp[n][1] = kret[1];

    if (sqrt(kret[0] * kret[0] + kret[1] * kret[1]) < kacc) {
      ug[0] =
          gtmp[n][0] / sqrt(gtmp[n][0] * gtmp[n][0] + gtmp[n][1] * gtmp[n][1]);
      ug[1] =
          gtmp[n][1] / sqrt(gtmp[n][0] * gtmp[n][0] + gtmp[n][1] * gtmp[n][1]);

      if (sqrt(gtmp[n][0] * gtmp[n][0] + gtmp[n][1] * gtmp[n][1]) > smax * T) {
        gtmp[n + 1][0] = gtmp[n][0] - ug[0] * smax * T;
        gtmp[n + 1][1] = gtmp[n][1] - ug[1] * smax * T;
      } else {
        gtmp[n + 1][0] = 0;
        gtmp[n + 1][1] = 0;
      }
      gtarg[0] = 0;
      gtarg[1] = 0;
    } else {
      gtarg[0] = gain * (-kret[0]);
      gtarg[1] = gain * (-kret[1]);
      if (sqrt(gtarg[0] * gtarg[0] + gtarg[1] * gtarg[1]) > gmax) {
        gtarg[0] =
            gmax * gtarg[0] / sqrt(gtarg[0] * gtarg[0] + gtarg[1] * gtarg[1]);
        gtarg[1] =
            gmax * gtarg[1] / sqrt(gtarg[0] * gtarg[0] + gtarg[1] * gtarg[1]);
      }
      dgvect[0] = gtarg[0] - g[n][0];
      dgvect[1] = gtarg[1] - g[n][1];
      if (sqrt(dgvect[0] * dgvect[0] + dgvect[1] * dgvect[1]) > smax * T) {
        dg[0] = dgvect[0] /
                sqrt(dgvect[0] * dgvect[0] + dgvect[1] * dgvect[1]) * smax * T;
        dg[1] = dgvect[1] /
                sqrt(dgvect[0] * dgvect[0] + dgvect[1] * dgvect[1]) * smax * T;
      } else {
        dg[0] = dgvect[0];
        dg[1] = dgvect[1];
      }
      gtmp[n + 1][0] = gtmp[n][0] + dg[0];
      gtmp[n + 1][1] = gtmp[n][1] + dg[1];
    } /* End if abs(kret)<kacc */
    karr[n][0] = k[0];
    karr[n][1] = k[1];

    /* Iterate k */
    k[0] = k[0] + GAMMA_PROTON * (gtmp[n][0] + g[n + 1][0]) / 2 * T;
    k[1] = k[1] + GAMMA_PROTON * (gtmp[n][1] + g[n + 1][1]) / 2 * T;
    n = n + 1;

    /* Stop conditions */
    if (n > 3000)
      keepgoing = 0;

    if ((sqrt(k[0] * k[0] + k[1] * k[1]) < kacc) &&
        (sqrt(gtmp[n][0] * gtmp[n][0] + gtmp[n][1] * gtmp[n][1]) < gacc))
      keepgoing = 0;
  } /* End while keepgoing */
  gx = (float *)AllocMem(n * sizeof(float));
  gy = (float *)AllocMem(n * sizeof(float));
  for (i = 0; i < n; i++) {
    gx[i] = gtmp[i][0] * cos(rotang) - gtmp[i][1] * sin(rotang);
    gy[i] = gtmp[i][0] * sin(rotang) + gtmp[i][1] * cos(rotang);
  }

  return 1;
}

int genktraj_vd(float **ktmp, float *dens, float *nint, float **grad,
                int grad_points, int GRAD_POINTS, float xyscale, float zscale,
                float *FOVrad, float *FOVcirc, float *FOVkr, float truetheta,
                float dtheta_ic, float GSP, float *maxgr) {
  int i;
  /* Hesitantly have moved to static declarations of the following with
     "reasonable"
     maximum gradient lengths since memory management doesn't seem to work
     properly
     in the final PSD code with dynamic allocation so it looked like a big
     memory
     leak and resulted in a huge memory footprint for the PSD */
  /*MCARL: reduce from 32768 to 8192*/
  /* Python wrapper: `static` keeps this off the stack (~3 MB at the bumped
     size). Single-threaded use only -- the wrapper doesn't release the GIL. */
  static float g[CG_MAX_WAVEFORM_PTS][3];
  static float gr[CG_MAX_WAVEFORM_PTS];
  static float kr[CG_MAX_WAVEFORM_PTS];
  static float Gtwist[CG_MAX_WAVEFORM_PTS];
  static float gtangle[CG_MAX_WAVEFORM_PTS];
  static float gtdiffangle[CG_MAX_WAVEFORM_PTS];
  static float gtdiffkr[CG_MAX_WAVEFORM_PTS];
  static float abskxy[CG_MAX_WAVEFORM_PTS];
  static float nints[CG_MAX_WAVEFORM_PTS];
  static float denscomp_ic[CG_MAX_WAVEFORM_PTS]; /* inter-cone spacing */
  static float denscomp_it[CG_MAX_WAVEFORM_PTS]; /* inter-trajectory spacing */
  static float denscomp_is[CG_MAX_WAVEFORM_PTS]; /* inter-sample spacing */
  float maxkr;
  float wrap_offset = 0.0;
  float meandens = 0;

  *maxgr = -1e6;

  /* VD */

  int iNkr = 100;
  int idx = 0;
  int idx2 = 0;
  /* Conegrad Python wrapper: replaced VLAs with static size matching the
     surrounding 8192-element arrays for MSVC compatibility. */
  static float FOVrad_i[CG_MAX_WAVEFORM_PTS];
  static float FOVcirc_i[CG_MAX_WAVEFORM_PTS];
  float slope;

  /* VD */

  /* K-space trajectory is simple the integral of the gradient trajectories
     with proper xy and z-scaling */
  for (i = 0; i < grad_points; i++) {
    g[i][0] = xyscale * grad[i][0];
    g[i][1] = xyscale * grad[i][1];
    g[i][2] = zscale * grad[i][2];
    gr[i] = sqrt(g[i][0] * g[i][0] + g[i][1] * g[i][1] + g[i][2] * g[i][2]);
    if (gr[i] > *maxgr)
      *maxgr = gr[i];
  }

  maxkr = 0;

  for (i = 0; i < grad_points; i++) {
    if (i == 0) {
      (ktmp[i][0]) = (float)(g[i][0] * GAMMA_PROTON * GSP);
      (ktmp[i][1]) = (float)(g[i][1] * GAMMA_PROTON * GSP);
      (ktmp[i][2]) = (float)(g[i][2] * GAMMA_PROTON * GSP);
    } else {
      (ktmp[i][0]) = (ktmp[i - 1][0]) + (float)(g[i][0] * GAMMA_PROTON * GSP);
      (ktmp[i][1]) = (ktmp[i - 1][1]) + (float)(g[i][1] * GAMMA_PROTON * GSP);
      (ktmp[i][2]) = (ktmp[i - 1][2]) + (float)(g[i][2] * GAMMA_PROTON * GSP);
    }

    kr[i] = sqrt((ktmp[i][0]) * (ktmp[i][0]) + (ktmp[i][1]) * (ktmp[i][1]) +
                 (ktmp[i][2]) * (ktmp[i][2]));
    abskxy[i] = sqrt((ktmp[i][0]) * (ktmp[i][0]) + (ktmp[i][1]) * (ktmp[i][1]));
    maxkr = (kr[i] > maxkr) ? kr[i] : maxkr;
  }

  /* VD */

  for (idx = 0; idx < iNkr; idx++) {
    FOVkr[idx] = FOVkr[idx] * maxkr;
  }

  for (idx = 0; idx < grad_points; idx++) {
    for (idx2 = 0; idx2 < (iNkr - 1); idx2++) {
      if (FOVkr[idx2] <= kr[idx] && FOVkr[idx2 + 1] >= kr[idx]) {
        slope =
            (FOVrad[idx2 + 1] - FOVrad[idx2]) / (FOVkr[idx2 + 1] - FOVkr[idx2]);
        FOVrad_i[idx] = slope * (kr[idx] - FOVkr[idx2]) + FOVrad[idx2];
      }
    }
  }

  for (idx = 0; idx < grad_points; idx++) {
    for (idx2 = 0; idx2 < (iNkr - 1); idx2++) {
      if (FOVkr[idx2] <= kr[idx] && FOVkr[idx2 + 1] >= kr[idx]) {
        slope = (FOVcirc[idx2 + 1] - FOVcirc[idx2]) /
                (FOVkr[idx2 + 1] - FOVkr[idx2]);
        FOVcirc_i[idx] = slope * (kr[idx] - FOVkr[idx2]) + FOVcirc[idx2];
      }
    }
  }

  /* VD */

  /*printf("MCARL: MCARL: MCARL: MCARL: MCARL: maxkr = %f\n",maxkr);
          fflush(stdout);*/

  /* Not all our gradient waveforms will be the full length we requested.
     Some versions of recon will be confused if we allocate a k-space position
     of
     (0,0,0) for these last few points we haven't designed. Thus give them
     a fixed, extreme value so they will not contribute to the final image if
     considered by the version of recon we're doing */
  for (i = grad_points; i < GRAD_POINTS;
       i++) {                  /*from grad_points(nint)->GRAD_POINTS*/
    ktmp[i][0] = KTRAJ_IGNORE; /*KTRAJ_IGNORE=-99 in scones_design.h*/
    ktmp[i][1] = KTRAJ_IGNORE;
    ktmp[i][2] = KTRAJ_IGNORE;
  }

  /* Calculate Gtwist for this trajectory and use this to calculate the
     fractional
     number of interleaves needed of this trajectory */
  /* Calculate k-space angle and attempt to unwrap phase angle */
  for (i = 0; i < grad_points; i++) {
    gtangle[i] = atan2f((float)(ktmp[i][1]), (float)(ktmp[i][0])) + wrap_offset;
    if (i > 0) {
      if (gtangle[i] < (gtangle[i - 1] - PI)) {
        gtangle[i] = gtangle[i] + 2 * PI;
        wrap_offset = wrap_offset + 2 * PI;
      }
      if (gtangle[i] > (gtangle[i - 1] + PI)) {
        gtangle[i] = gtangle[i] - 2 * PI;
        wrap_offset = wrap_offset - 2 * PI;
      }
    }
  }

  *nint = 0.0;
  /* JAS - some issue with accuracy of last point between running integrated in
     PSD on scanner and running offline (e.g. in standalone recon - back off one
     point seems to help ... */
  int Mynint = 0;
  int MyCounter = 0;
  for (i = 0; i < grad_points; i++) {
    gtdiffkr[i] = max(kr[i + 1] - kr[i], 1e-8);
    gtdiffangle[i] = gtangle[i + 1] - gtangle[i];
    Gtwist[i] = gtdiffangle[i] / gtdiffkr[i] * abskxy[i];
  }
  Gtwist[grad_points - 1] = Gtwist[grad_points - 2];

  for (i = 0; i < grad_points; i++) {
    if (!wcfEqual(abskxy[i], 0)) {
      /*nints[i] =
       * 2*PI*abskxy[i]*FOVrad/sqrt(max(abs(Gtwist[i]*Gtwist[i]+(FOVrad/FOVcirc)*(FOVrad/FOVcirc)),1));*/
      nints[i] =
          2 * PI * abskxy[i] * FOVrad_i[i] /
          sqrt(Gtwist[i] * Gtwist[i] +
               (FOVrad_i[i] / FOVcirc_i[i]) * (FOVrad_i[i] / FOVcirc_i[i]));
      if ((*nint) < nints[i])
        *nint = nints[i];
      if (Mynint < nints[i])
        Mynint = nints[i];

    } /* End abskxy != 0 */
  }

  for (i = 0; i < grad_points; i++) {
    /* Special case of the NULL cone */
    if ((abskxy[grad_points - 1] / kr[grad_points - 1]) < 1e-8) {
      denscomp_ic[i] = dtheta_ic * dtheta_ic * kr[i] * kr[i] / 8;
      denscomp_it[i] = 1;
    } else {
      denscomp_ic[i] = dtheta_ic * kr[i];
      /*denscomp_it[i] = abskxy[i]/sqrt(Gtwist[i]*Gtwist[i]+1);*/

      /*denscomp_it[i] =
       * abskxy[i]*(FOVrad/FOVcirc)/sqrt(max(abs(Gtwist[i]*Gtwist[i]+(FOVrad/FOVcirc)*(FOVrad/FOVcirc)),1));*/
      denscomp_it[i] =
          abskxy[i] * (FOVrad_i[i] / FOVcirc_i[i]) /
          sqrt(Gtwist[i] * Gtwist[i] +
               (FOVrad_i[i] / FOVcirc_i[i]) * (FOVrad_i[i] / FOVcirc_i[i]));

      /*denscomp_it[i] = abskxy[i]/(maxkr*sqrt(Gtwist[i]*Gtwist[i]+1));*/ /*Confirmed
                                                                             wrong*/
      /*denscomp_it[i] = abskxy[i]/sqrt(Gtwist[i]*Gtwist[i]+(FOVrad/FOVcirc)*(FOVrad/FOVcirc));*/ /*Confirmed wrong*/
    }
    if (i == grad_points - 1)
      denscomp_is[i] = (gr[i] + gr[i]) / 2;
    else
      denscomp_is[i] = (gr[i] + gr[i + 1]) / 2;

    dens[i] = (float)(denscomp_ic[i] * denscomp_it[i] *
                      denscomp_is[i]); /*Multiply three orthogonal terms*/
  }

  for (i = grad_points; i < GRAD_POINTS; i++) /*fill in rest with zeros*/
    dens[i] = (float)DENS_IGNORE;

  for (i = 0; i < GRAD_POINTS; i++) { /*Double-check array*/
    if (abs(dens[i]) < 1e10)
      MyCounter++;
    else
      printf("MCARL: MCARL: MCARL: MCARL: MCARL: dens[%d] = %f\n", i, dens[i]);
    meandens = meandens + dens[i];
  }
  meandens = meandens / (float)GRAD_POINTS;

  /*printf("MCARL: MCARL: MCARL: genktraj: meandens = %f\n",meandens);
  fflush(stdout);*/

  return 1;
}

/* Function:  genktraj
 *
 * Generates the k-space trajectory information and the true number of
 * interleaves
 * needed for the given cones waveforms.
 */
int genktraj(float **ktmp, float *dens, float *nint, float **grad,
             int grad_points, int GRAD_POINTS, float xyscale, float zscale,
             float FOVrad, float FOVcirc, float truetheta, float dtheta_ic,
             float GSP, float *maxgr) {
  int i;
  /* Hesitantly have moved to static declarations of the following with
     "reasonable"
     maximum gradient lengths since memory management doesn't seem to work
     properly
     in the final PSD code with dynamic allocation so it looked like a big
     memory
     leak and resulted in a huge memory footprint for the PSD */
  /*MCARL: reduce from 32768 to 8192*/
  /* Python wrapper: `static` keeps this off the stack (~3 MB at the bumped
     size). Single-threaded use only -- the wrapper doesn't release the GIL. */
  static float g[CG_MAX_WAVEFORM_PTS][3];
  static float gr[CG_MAX_WAVEFORM_PTS];
  static float kr[CG_MAX_WAVEFORM_PTS];
  static float Gtwist[CG_MAX_WAVEFORM_PTS];
  static float gtangle[CG_MAX_WAVEFORM_PTS];
  static float gtdiffangle[CG_MAX_WAVEFORM_PTS];
  static float gtdiffkr[CG_MAX_WAVEFORM_PTS];
  static float abskxy[CG_MAX_WAVEFORM_PTS];
  static float nints[CG_MAX_WAVEFORM_PTS];
  static float denscomp_ic[CG_MAX_WAVEFORM_PTS]; /* inter-cone spacing */
  static float denscomp_it[CG_MAX_WAVEFORM_PTS]; /* inter-trajectory spacing */
  static float denscomp_is[CG_MAX_WAVEFORM_PTS]; /* inter-sample spacing */
  float maxkr;
  float wrap_offset = 0.0;
  float meandens = 0;

  *maxgr = -1e6;

  /* K-space trajectory is simple the integral of the gradient trajectories
     with proper xy and z-scaling */
  for (i = 0; i < grad_points; i++) {
    g[i][0] = xyscale * grad[i][0];
    g[i][1] = xyscale * grad[i][1];
    g[i][2] = zscale * grad[i][2];
    gr[i] = sqrt(g[i][0] * g[i][0] + g[i][1] * g[i][1] + g[i][2] * g[i][2]);
    if (gr[i] > *maxgr)
      *maxgr = gr[i];
  }

  maxkr = 0;

  for (i = 0; i < grad_points; i++) {
    if (i == 0) {
      (ktmp[i][0]) = (float)(g[i][0] * GAMMA_PROTON * GSP);
      (ktmp[i][1]) = (float)(g[i][1] * GAMMA_PROTON * GSP);
      (ktmp[i][2]) = (float)(g[i][2] * GAMMA_PROTON * GSP);
    } else {
      (ktmp[i][0]) = (ktmp[i - 1][0]) + (float)(g[i][0] * GAMMA_PROTON * GSP);
      (ktmp[i][1]) = (ktmp[i - 1][1]) + (float)(g[i][1] * GAMMA_PROTON * GSP);
      (ktmp[i][2]) = (ktmp[i - 1][2]) + (float)(g[i][2] * GAMMA_PROTON * GSP);
    }

    kr[i] = sqrt((ktmp[i][0]) * (ktmp[i][0]) + (ktmp[i][1]) * (ktmp[i][1]) +
                 (ktmp[i][2]) * (ktmp[i][2]));
    abskxy[i] = sqrt((ktmp[i][0]) * (ktmp[i][0]) + (ktmp[i][1]) * (ktmp[i][1]));
    maxkr = (kr[i] > maxkr) ? kr[i] : maxkr;
  }

  /*printf("MCARL: MCARL: MCARL: MCARL: MCARL: maxkr = %f\n",maxkr);
          fflush(stdout);*/

  /* Not all our gradient waveforms will be the full length we requested.
     Some versions of recon will be confused if we allocate a k-space position
     of
     (0,0,0) for these last few points we haven't designed. Thus give them
     a fixed, extreme value so they will not contribute to the final image if
     considered by the version of recon we're doing */
  for (i = grad_points; i < GRAD_POINTS;
       i++) {                  /*from grad_points(nint)->GRAD_POINTS*/
    ktmp[i][0] = KTRAJ_IGNORE; /*KTRAJ_IGNORE=-99 in scones_design.h*/
    ktmp[i][1] = KTRAJ_IGNORE;
    ktmp[i][2] = KTRAJ_IGNORE;
  }

  /* Calculate Gtwist for this trajectory and use this to calculate the
     fractional
     number of interleaves needed of this trajectory */
  /* Calculate k-space angle and attempt to unwrap phase angle */
  for (i = 0; i < grad_points; i++) {
    gtangle[i] = atan2f((float)(ktmp[i][1]), (float)(ktmp[i][0])) + wrap_offset;
    if (i > 0) {
      if (gtangle[i] < (gtangle[i - 1] - PI)) {
        gtangle[i] = gtangle[i] + 2 * PI;
        wrap_offset = wrap_offset + 2 * PI;
      }
      if (gtangle[i] > (gtangle[i - 1] + PI)) {
        gtangle[i] = gtangle[i] - 2 * PI;
        wrap_offset = wrap_offset - 2 * PI;
      }
    }
  }

  *nint = 0.0;
  /* JAS - some issue with accuracy of last point between running integrated in
     PSD on scanner and running offline (e.g. in standalone recon - back off one
     point seems to help ... */
  int Mynint = 0;
  int MyCounter = 0;
  for (i = 0; i < grad_points; i++) {
    gtdiffkr[i] = max(kr[i + 1] - kr[i], 1e-8);
    gtdiffangle[i] = gtangle[i + 1] - gtangle[i];
    Gtwist[i] = gtdiffangle[i] / gtdiffkr[i] * abskxy[i];
  }
  Gtwist[grad_points - 1] = Gtwist[grad_points - 2];

  for (i = 0; i < grad_points; i++) {
    if (!wcfEqual(abskxy[i], 0)) {
      /*nints[i] =
       * 2*PI*abskxy[i]*FOVrad/sqrt(max(abs(Gtwist[i]*Gtwist[i]+(FOVrad/FOVcirc)*(FOVrad/FOVcirc)),1));*/
      nints[i] =
          2 * PI * abskxy[i] * FOVrad /
          sqrt(Gtwist[i] * Gtwist[i] + (FOVrad / FOVcirc) * (FOVrad / FOVcirc));
      if ((*nint) < nints[i])
        *nint = nints[i];
      if (Mynint < nints[i])
        Mynint = nints[i];

    } /* End abskxy != 0 */
  }

  for (i = 0; i < grad_points; i++) {
    /* Special case of the NULL cone */
    if ((abskxy[grad_points - 1] / kr[grad_points - 1]) < 1e-8) {
      denscomp_ic[i] = dtheta_ic * dtheta_ic * kr[i] * kr[i] / 8;
      denscomp_it[i] = 1;
    } else {
      denscomp_ic[i] = dtheta_ic * kr[i];
      /*denscomp_it[i] = abskxy[i]/sqrt(Gtwist[i]*Gtwist[i]+1);*/

      /*denscomp_it[i] =
       * abskxy[i]*(FOVrad/FOVcirc)/sqrt(max(abs(Gtwist[i]*Gtwist[i]+(FOVrad/FOVcirc)*(FOVrad/FOVcirc)),1));*/
      denscomp_it[i] =
          abskxy[i] * (FOVrad / FOVcirc) /
          sqrt(Gtwist[i] * Gtwist[i] + (FOVrad / FOVcirc) * (FOVrad / FOVcirc));

      /*denscomp_it[i] = abskxy[i]/(maxkr*sqrt(Gtwist[i]*Gtwist[i]+1));*/ /*Confirmed
                                                                             wrong*/
      /*denscomp_it[i] = abskxy[i]/sqrt(Gtwist[i]*Gtwist[i]+(FOVrad/FOVcirc)*(FOVrad/FOVcirc));*/ /*Confirmed wrong*/
    }
    if (i == grad_points - 1)
      denscomp_is[i] = (gr[i] + gr[i]) / 2;
    else
      denscomp_is[i] = (gr[i] + gr[i + 1]) / 2;

    dens[i] = (float)(denscomp_ic[i] * denscomp_it[i] *
                      denscomp_is[i]); /*Multiply three orthogonal terms*/
  }

  for (i = grad_points; i < GRAD_POINTS; i++) /*fill in rest with zeros*/
    dens[i] = (float)DENS_IGNORE;

  for (i = 0; i < GRAD_POINTS; i++) { /*Double-check array*/
    if (abs(dens[i]) < 1e10)
      MyCounter++;
    else
      printf("MCARL: MCARL: MCARL: MCARL: MCARL: dens[%d] = %f\n", i, dens[i]);
    meandens = meandens + dens[i];
  }
  meandens = meandens / (float)GRAD_POINTS;

  /*printf("MCARL: MCARL: MCARL: genktraj: meandens = %f\n",meandens);
  fflush(stdout);*/

  return 1;
}

/* VD Start */

/*MCARL: Added Cones_Plot_Flag*/
int conegrad_vd(int *nintpc, float *rspthetas, float *RES, float *FOV,
                int *gxcones_host, int *gycones_host, int *gzcones_host, /* CMS */
                int NUMCONES, int GRAD_POINTS, float GSP, int READ_POINTS,
                float TSP, float PRECISION, float DCF, int OVERSAMPLE,
                float SMAX, float maxSRewind, float GMAX, float sysGMAX,
                float SysMaxRewindG, int output_grad, int ktraj_flag,
                int ktraj_out_flag, int endian_flag, int rewind_flag,
                float MINDENS, float *maxgr2, float *snr_eff, int *ngradact,
                float t_fracx, float t_fracy, float t_fracz, int acq_mode,
                float rot_flag, int SymFlag, float SlabKz,
                int rhkacq_uid, int Cones_Plot_Flag,
                int *traj_length) /*MCARL: Added t_frac, acq_mode, rot_flag,
                                     rhkacq_uid, Cones_Plot_Flag;
                                     Python wrapper added traj_length. */
{

  float **g;
  float ***k;   /* Make float as want to write to file as floats */
  float **dens; /* Make float as want to write to file as floats */
  float ***grad;
  float tmpgr;
  /* Switching to static allocation, allow up to 1000 thetas, and 100 cones */
  int conenum[MAX_THETAS];
  float xyscale[MAX_THETAS], zscale[MAX_THETAS];
  float thetas[2 * MAX_THETAS];
  float dthetas_ic[2 * MAX_THETAS];
  float thetaFinals[20 * MAX_THETAS];
  float dthetaFinals_ic[20 * MAX_THETAS];
  float dtheta_ic[20 * MAX_THETAS];
  int grad_points[MAX_CONES];
  float fnint[MAX_CONES];
  float fnints[2 * MAX_THETAS];
  float dnints[2 * MAX_THETAS];
  float *truetheta;
  float thetaup, thetadown;
  int length, totlength;
  float nint;
  float dnint;
  float *theta;
  float THETA[2];
  int ntheta1;
  int ntheta7;
  int nthetaFinal;
  int ntheta;
  float FOV1[2], FOV7[2], FOVFinal[2];
  float FOVFactor1 = 1;
  float FOVFactor7 = 7;
  float FOVFactorFinal;
  float MyRatioFinal = 4;
  float MyRatio1;
  float MyRatio7;
  float MyRatioCheck;
  float MySlope;
  float **kRot;

  int nthetatmp;
  int totints;
  int estints;
  int i, j, m;
  int nramp;
  FILE *gxf = NULL, *gyf = NULL, *gzf = NULL, *npcf = NULL;
  int tmpout;
  float tmpout2;
  float maxnints = 0.0;
  float sdens, sdens2, sdens3;
  float snr_eff_num, snr_eff_den;
  int snr_numpts;
  float snr_eff_den2;
  float dtmpx;
  float dtmpy;
  float dtmpz;
  float dtmp;
  /* For gradient rewinder calcs */
  float max_kpos;
  float rewinder_area1;
  int nramp_min;

  *maxgr2 = -1e6; /*MCARL: Need to initialize*/

  /* VD */

  float Nkr = 100.0;
  int iNkr = 100;
  float FOVrad_vd[100];
  int idx;
  float FOVcirc_vd[100];
  float FOVkr_vd[100];
  float FOVunderTemp[2];

  /* VD */

  /* Default SNR efficiency - can tell if not calculated later */
  *snr_eff = 0;
  estints = 0;
  totints = 0;

  /* The maximum length waveform actually created */
  (*ngradact) = 0; /* initialize to 0 */

  /* Set the number of points in the ramp of gradient back to zero to the max
     number that could be needed */
  nramp = (int)ceil(GMAX / (maxSRewind * GSP) + 5);
  // printf("MCARL: MCARL: MCARL: nramp1 = %d\n",nramp);
  // fflush(stdout);

  /* Calculate number of "extra" gradient points at the end to account for time
     to
     ramp back down to 0 then to apply a rewinder gradient to bring the
     trajectory
     back to k=0. Consider the worst case scenario where for each Gx or Gy
     gradient
     we could be at kx=maxK or ky=maxK. Consider we're at that point, then go
     further
     while we ramp back to zero then need time to play inverse gradient.
  */
  /* First find the furthest excursion in k-space we could be (based on encoded
     resolution and gradient ramping back down to zero
  */
  if (rewind_flag == 1) {
    max_kpos = (float)5 / min(RES[0], RES[1]) +
               0.5 * GAMMA_PROTON * GMAX * GMAX /
                   maxSRewind; /*MCARL: Second term is additional k-sp accrual
                                  during ramp-down*/
    /* Next calculate the area under a rewinder at full amplitude and min
     * duration */
    rewinder_area1 = GAMMA_PROTON * SysMaxRewindG * SysMaxRewindG / maxSRewind +
                     GAMMA_PROTON * SysMaxRewindG *
                         GSP; /*MCARL: Max area achivable during triangle,
                                 including one time-step point of flattop
                                 (second term)*/
    nramp_min = (int)(ceil(GMAX / (maxSRewind * GSP)) +
                      2 * ceil(SysMaxRewindG / (maxSRewind * GSP)) +
                      1); /*MCARL: time points to ramp down maxG plus to make
                             fast triangle to sysGmax*/
    if (rewinder_area1 > max_kpos) { /*Already have sufficent rewinder*/
      nramp = nramp_min;
    } else { /*Add flattop if needed*/
      nramp = nramp_min + (int)ceil((max_kpos - rewinder_area1) /
                                    (GAMMA_PROTON * SysMaxRewindG *
                                     GSP)); /*MCARL: GAMMA_PROTON was missing*/
    }
    nramp = nramp + 5; // For saftey
  }

  // printf("GRAD_POINTS = %d\n",GRAD_POINTS);
  // printf("nramp = %d\n",nramp);
  // fflush(stdout);

  /* Determine thetas needed for prescription */

  /* VD */

  nthetatmp = (int)((PI / 2) /
                    (1 / max(FOV[((int)(Nkr)) - 1], FOV[2 * ((int)(Nkr)) - 1]) /
                     (5 / min(RES[0], RES[1]))));

  theta = (float *)malloc(20 * nthetatmp * sizeof(float));

  /* VD */
  FOVunderTemp[0] = FOV[iNkr - 1];
  FOVunderTemp[1] = FOV[2 * iNkr - 1];
  /* VD */

  if (SymFlag == 1) {
    getthetas2(theta, RES, FOVunderTemp, &ntheta, &dtheta_ic[0]);
    // Extend to -pi/2:pi/2
    for (i = 0; i < ntheta - 1; i++) {
      thetas[i] = -1 * theta[ntheta - i - 1];
      dthetas_ic[i] = dtheta_ic[ntheta - i - 1];
      rspthetas[i] = -1 * theta[ntheta - i - 1];
    }
    for (i = 0; i < ntheta; i++) {
      thetas[ntheta - 1 + i] = theta[i];
      dthetas_ic[ntheta - 1 + i] = dtheta_ic[i];
      rspthetas[ntheta - 1 + i] = theta[i];
    }
  }

  if (SymFlag == 2) {
    FOV1[0] = FOVFactor1 * (FOVunderTemp[0]);
    FOV1[1] = FOVFactor1 * (FOVunderTemp[1]);
    FOV7[0] = FOVFactor7 * (FOVunderTemp[0]);
    FOV7[1] = FOVFactor7 * (FOVunderTemp[1]);
    getthetas(RES, FOV1, &ntheta1, &dtheta_ic[0]);
    getthetas(RES, FOV7, &ntheta7, &dtheta_ic[0]);
    ntheta = ntheta1;
    MyRatio1 = (2 * ntheta1) / (2 * ntheta1 - 1.0);
    MyRatio7 = (2 * ntheta7) / (2 * ntheta1 - 1.0);
    MySlope = (MyRatio7 - MyRatio1) / (FOVFactor7 - FOVFactor1);
    FOVFactorFinal = (MyRatioFinal - MyRatio1) / MySlope + FOVFactor1;
    FOVFinal[0] = FOVFactorFinal * (FOVunderTemp[0]);
    FOVFinal[1] = FOVFactorFinal * (FOVunderTemp[1]);
    getthetas2(theta, RES, FOVFinal, &nthetaFinal, &dtheta_ic[0]);
    MyRatioCheck = (2 * nthetaFinal) / (2 * ntheta1 - 1.0);
    printf("MCARL: MCARL: MCARL: MyRatioCheck = %f\n", MyRatioCheck);
    fflush(stdout);
    // Extend to -pi/2:pi/2
    for (i = 0; i < nthetaFinal - 1; i++) {
      thetaFinals[i] = -1 * theta[nthetaFinal - i - 1];
      dthetaFinals_ic[i] = dtheta_ic[nthetaFinal - i - 1];
    }
    for (i = 0; i < nthetaFinal; i++) {
      thetaFinals[nthetaFinal - 1 + i] = theta[i];
      dthetaFinals_ic[nthetaFinal - 1 + i] = dtheta_ic[i];
    }
    // Pick Correct Thetas
    for (i = 0; i < 2 * ntheta1 - 1; i++) {
      thetas[i] = thetaFinals[(int)(i * MyRatioFinal)];
      rspthetas[i] = thetaFinals[(int)(i * MyRatioFinal)];
      dthetas_ic[i] = dthetaFinals_ic[(int)(i * MyRatioFinal)];
    }
  }

  FreeMem(theta);

  /* VD */

  /* Allocate memory for gradient waveform arrays */
  grad = (float ***)AllocMem(NUMCONES * sizeof(float **));
  for (i = 0; i < NUMCONES; i++) {
    grad[i] = (float **)AllocMem((GRAD_POINTS + nramp) *
                                 sizeof(float *)); /* ML - Not freed */
    for (j = 0; j < GRAD_POINTS + nramp; j++)
      grad[i][j] =
          (float *)AllocMem(3 * sizeof(float)); /* ML - Not freed (lots) */
  }

  /* Python wrapper: oversize the g buffer to handle aggressive undersampling
     where the cones designer's binary search can produce per-cone lengths
     slightly past LEN+nramp. Cost is modest (~256 KB extra per axis for
     typical configs); benefit is no heap corruption when fov_under is small. */
  g = (float **)AllocMem(3 * sizeof(float *));
  g[0] = (float *)AllocMem(((GRAD_POINTS + nramp) * OVERSAMPLE + 65536) * sizeof(float));
  g[1] = (float *)AllocMem(((GRAD_POINTS + nramp) * OVERSAMPLE + 65536) * sizeof(float));
  g[2] = (float *)AllocMem(((GRAD_POINTS + nramp) * OVERSAMPLE + 65536) * sizeof(float));

  char GXCONE_FILE[128];
  char GYCONE_FILE[128];
  char GZCONE_FILE[128];

  sprintf(GXCONE_FILE, "%s/%s.%d", DIR_OUT, FILE_GXCONE, rhkacq_uid);
  sprintf(GYCONE_FILE, "%s/%s.%d", DIR_OUT, FILE_GYCONE, rhkacq_uid);
  sprintf(GZCONE_FILE, "%s/%s.%d", DIR_OUT, FILE_GZCONE, rhkacq_uid);
  /*
#ifdef SIM
  sprintf(GXCONE_FILE,"%s%d","kacq_gxcone_",rhkacq_uid);
  sprintf(GYCONE_FILE,"%s%d","kacq_gycone_",rhkacq_uid);
  sprintf(GZCONE_FILE,"%s%d","kacq_gzcone_",rhkacq_uid);
#else
  sprintf(GXCONE_FILE,"%s%d","/usr/g/psddata/kacq_gxcone_",rhkacq_uid);
  sprintf(GYCONE_FILE,"%s%d","/usr/g/psddata/kacq_gycone_",rhkacq_uid);
  sprintf(GZCONE_FILE,"%s%d","/usr/g/psddata/kacq_gzcone_",rhkacq_uid);
#endif
  */

  /* Open output files if needed */
  if (output_grad == 1) {
    if ((gxf = fopen(GXCONE_FILE, "w")) == NULL) {
      printf("Could not open gxfile for writing\n");
    }
    if ((gyf = fopen(GYCONE_FILE, "w")) == NULL) {
      printf("Could not open gyfile for writing\n");
    }
    if ((gzf = fopen(GZCONE_FILE, "w")) == NULL) {
      printf("Could not open gzfile for writing\n");
    }
  }

  /* Loop over requested number of cones to design and do the design */
  for (i = 0; i < NUMCONES; i++) {

    printf("MCARL: MCARL: MCARL: numcone = %d\n", i);
    fflush(stdout);

    /* Calculate range of thetas for current cone */
    THETA[0] = i * PI / 2 / NUMCONES;
    THETA[1] = (i + 1) * PI / 2 / NUMCONES;

    /*Call findcones for each numcone loop -> THETA = THETA[i]*/
    findcone_vd(&g, &nint, &length, &totlength, RES, FOV, GRAD_POINTS, THETA,
                PRECISION, GSP, OVERSAMPLE, SMAX, maxSRewind, GMAX,
                SysMaxRewindG, DCF, MINDENS, nramp, rewind_flag, acq_mode);
    fnint[i] = nint;

    /*printf("MCARL: MCARL: MCARL: (*ngradact)[%d] = %d\n",i,(*ngradact));
      fflush(stdout);*/

    if (totlength > (*ngradact))
      (*ngradact) = totlength;

    /* Estimate the number of interleaves for reporting to scan UI */
    estints =
        estints + (int)((float)ntheta / (float)NUMCONES) * (int)ceil(nint);

    /* Fill in full gradient waveform array based on current cone design
     * waveforms */
    for (j = 0; j < GRAD_POINTS + nramp; j++) {
      grad[i][j][0] = g[0][j];
      grad[i][j][1] = g[1][j];
      grad[i][j][2] = g[2][j];

    } /* End if j<GRAD_POINTSGTH+nramp */

    grad_points[i] = length;
    /* Python wrapper: also expose this per-cone trajectory length so callers
       can slice readouts at the true end of the spiral, before the rewinder. */
    traj_length[i] = length;
    /*printf("MCARL: MCARL: MCARL: grad_points[%d] = %d\n",i,grad_points[i]);
      fflush(stdout);*/

  } /* End for NUMCONES */

  printf("MCARL: Gradients generated\n");
  fflush(stdout);

  /* JYCHENG: re-order cones here!!!! */

  // MCARL: begin gnuplot (Gradients)
  if ((Cones_Plot_Flag == 1) | (Cones_Plot_Flag == 3)) {
    // write grad wavefrom file
    int indx;
    char temp_name1[128], temp_name2[128];
    FILE *asciiwavefile = NULL;
    sprintf(temp_name1, "grad_waveform.txt");
    asciiwavefile = fopen(temp_name1, "w");
    for (indx = 0; indx < GRAD_POINTS + nramp; indx++)
    /*for (indx = 0; indx < GRAD_POINTS; indx++) */ /*MCARL: Only plot G during DAQ*/
    {
      fprintf(asciiwavefile, "%d %g %g %g\n", indx, 10 * grad[0][indx][0],
              10 * grad[0][indx][1], 10 * grad[0][indx][2]);
    }
    fclose(asciiwavefile);
    fflush(asciiwavefile);
    // populate asciiscriptfile
    FILE *asciiscriptfile = NULL;
    sprintf(temp_name2, "plot_grad_waveform.sh");
    asciiscriptfile = fopen(temp_name2, "w");
    char cmd_line[128];
    fprintf(asciiscriptfile, "#!/usr/bin/gnuplot -persist\n");
    fprintf(asciiscriptfile, "set multiplot\n");
    fprintf(asciiscriptfile, "set size 1,1;\n");
    fprintf(asciiscriptfile, "set origin 0.0,0.0;\n");
    fprintf(asciiscriptfile, "set xlabel \"Points\"\n");
    fprintf(asciiscriptfile, "plot \"grad_waveform.txt\" using 1:2 title "
                             "\"Gx\" ,\"grad_waveform.txt\" using 1:3 title "
                             "\"Gy\",\"grad_waveform.txt\" using 1:4 title "
                             "\"Gz\"\n");
    fclose(asciiscriptfile);
    // execute asciiscriptfile
    sprintf(cmd_line, "chmod u+x %s", temp_name2);
    system(cmd_line);
    sprintf(cmd_line, "%s", temp_name2);
    system(cmd_line);
    // MCARL: end gnuplot
  } // end if (Cones_Plot_Flag==1)

  /* Loop over writing gradient files out outside above loops so we can keep
     track of the
     length of each trajectory and only write out as few total points as
     necessary
  */
  for (i = 0; i < NUMCONES; i++) {
    for (j = 0; j < (*ngradact); j++) {

      /* Write out gradient waveforms to file, convert endianness if needed
         (yes needed if reading on target side) */
      /*MCARL: Gradiens to be read by pulsegen()*/
      /* Conegrad Python wrapper: always fill the host arrays so callers get
         the gradient waveforms back regardless of file-output flag. File
         writing is still gated by output_grad. */
      tmpout = (int)(grad[i][j][0] / sysGMAX * MAX_PG_WAMP);
      gxcones_host[j + i*(*ngradact)] = tmpout;
      if (output_grad == 1 && gxf != NULL) {
        fwrite(&tmpout, sizeof(int), 1, gxf);
      }

      tmpout = (int)(grad[i][j][1] / sysGMAX * MAX_PG_WAMP);
      gycones_host[j + i*(*ngradact)] = tmpout;
      if (output_grad == 1 && gyf != NULL) {
        fwrite(&tmpout, sizeof(int), 1, gyf);
      }

      tmpout = (int)(grad[i][j][2] / sysGMAX * MAX_PG_WAMP);
      gzcones_host[j + i*(*ngradact)] = tmpout;
      if (output_grad == 1 && gzf != NULL) {
        fwrite(&tmpout, sizeof(int), 1, gzf);
      }
    }
  }

  printf("MCARL: After output_grad\n");
  fflush(stdout);

  FreeMem(g[0]);
  FreeMem(g[1]);
  FreeMem(g[2]);
  FreeMem(g);

  printf("MCARL: Gradients written\n");
  fflush(stdout);

  /* VD */

  for (idx = 0; idx < ((int)Nkr); idx++) {
    FOVcirc_vd[idx] = FOV[idx];
  }

  /* VD */

  /* VD */

  for (idx = 0; idx < ((int)Nkr); idx++) {
    FOVkr_vd[idx] = FOV[idx + 2 * ((int)Nkr)];
  }

  /* VD */

  /* ### Generate kspace trajectory (for recon) and density file ### */
  if (ktraj_flag == 1) {
    FILE *kf = NULL, *df = NULL;

    /* Allocate memory for k-space trajectories */
    k = (float ***)AllocMem((2 * ntheta - 1) * sizeof(float **));
    for (i = 0; i < 2 * ntheta - 1; i++) {
      k[i] = (float **)AllocMem(GRAD_POINTS * sizeof(float *));
      for (j = 0; j < GRAD_POINTS; j++)
        k[i][j] = (float *)AllocMem(3 * sizeof(float));
    }
    kRot = (float **)AllocMem(GRAD_POINTS * sizeof(float *));
    for (j = 0; j < GRAD_POINTS; j++)
      kRot[j] = (float *)AllocMem(3 * sizeof(float));

    /* Allocate memory for density information */
    dens = (float **)AllocMem((2 * ntheta - 1) * sizeof(float *));
    truetheta = (float *)AllocMem((2 * ntheta - 1) * sizeof(float));

    /* Loop over all actual cones (vs designed cones) and determine the kspace
       trajectory for each one */
    for (i = 0; i < 2 * ntheta - 1; i++) {

      dens[i] =
          (float *)AllocMem(GRAD_POINTS *
                            sizeof(float)); /*After this, dens is a 2D array
                                               dens[ntheta][GRAD_POINTS], used
                                               later dens[i][j].*/

      /* Which designed cone is appropriate for this actual polar angel */
      conenum[i] = min(
          NUMCONES, max(1, (int)ceil(fabs((float)thetas[i]) / (PI / 2) *
                                     NUMCONES))); /*Will feed
                                                     grad_points[conenum[i]-1]
                                                     into gentraj*/
      thetaup = (float)(conenum[i]) / (float)(NUMCONES)*PI / 2;
      thetadown = (float)(conenum[i] - 1) / (float)(NUMCONES)*PI / 2;
      /* How to scale the designed cone to achieve the current angle */
      xyscale[i] = cos(thetas[i]) / cos(thetadown);
      zscale[i] = sin(thetas[i]) / sin(thetaup);

      truetheta[i] = atan2f(sin(thetas[i]) / RES[1], cos(thetas[i]) / RES[0]);
      for (idx = 0; idx < iNkr; idx++) {
        FOVrad_vd[idx] =
            1 / intlineellipse(1 / FOV[idx], 1 / FOV[idx + iNkr], truetheta[i]);
      }

      /* Generate the k-space trajectory */
      genktraj_vd(k[i], dens[i], &dnint, grad[conenum[i] - 1],
                  grad_points[conenum[i] - 1], GRAD_POINTS, (float)xyscale[i],
                  (float)zscale[i], FOVrad_vd, FOVcirc_vd, FOVkr_vd,
                  (float)truetheta[i], (float)dthetas_ic[i], (float)GSP,
                  &tmpgr);
      if (tmpgr > *maxgr2)
        *maxgr2 = tmpgr;

      dnints[i] =
          dnint; /* Number of interleaves reqd to meet FOV constraints */
      fnints[i] = fnint[conenum[i] -
                        1]; /* Number of interleaves used by design algorithm */
    }                       /* End for i<2*ntheta-1 */

    printf("MCARL: ktraj generated\n");
    fflush(stdout);

    if (MINDENS > 0) {
      for (i = 0; i < 2 * ntheta - 1; i++) {
        if ((dnints[i] / fnints[i]) > maxnints)
          maxnints = dnints[i] / fnints[i];
      }
      for (i = 0; i < 2 * ntheta - 1; i++)
        dnints[i] = fnints[i] * maxnints;
    }

    for (i = 0; i < 2 * ntheta - 1; i++) {
      /* Save the number of interleaves */
      nintpc[i] = (int)ceil(dnints[i]);
      /*MCARL: Off Iso: Round up to nearest odd*/
      if (nintpc[i] > 2) {
        if (SymFlag == 1)
          nintpc[i] = (int)ceil(dnints[i] / 2) * 2;
        if (SymFlag == 2)
          nintpc[i] = (int)ceil(dnints[i] / 2) * 2 - 1;
      }
      /*printf("MCARL: nintpc[%d] = %d\n",i,nintpc[i]);
      fflush(stdout);*/
      totints = totints + nintpc[i];
    }

    /* ### Write out k-space trajectory and density information ### */
    /* Python wrapper: dataskip_float pulled out of the ktraj_out_flag block so
       it stays in scope for the (now-unconditional) SNR + write code below. */
    float dataskip_float = TSP / GSP;
    if (ktraj_out_flag) {

      char KTRAJ_FILE2[128];
      char DENS_FILE2[128];
      char NPC_FILE2[128];

      sprintf(KTRAJ_FILE2, "%s/%s.%d", DIR_OUT, FILE_KTRAJ, rhkacq_uid);
      sprintf(DENS_FILE2, "%s/%s.%d", DIR_OUT, FILE_DENS, rhkacq_uid);
      sprintf(NPC_FILE2, "%s/%s.%d", DIR_OUT, FILE_NPC, rhkacq_uid);
      /*
#ifdef SIM
sprintf(KTRAJ_FILE2,"%s%d","kacq_cone_ktraj_",rhkacq_uid);
sprintf(DENS_FILE2,"%s%d","kacq_cone_dens_",rhkacq_uid);
sprintf(NPC_FILE2,"%s%d","kacq_cone_npc_",rhkacq_uid);
#else
sprintf(KTRAJ_FILE2,"%s%d","/usr/g/psddata/kacq_cone_ktraj_",rhkacq_uid);
sprintf(DENS_FILE2,"%s%d","/usr/g/psddata/kacq_cone_dens_",rhkacq_uid);
sprintf(NPC_FILE2,"%s%d","/usr/g/psddata/kacq_cone_npc_",rhkacq_uid);
#endif
      */

      kf = fopen(KTRAJ_FILE2, "w");
      df = fopen(DENS_FILE2, "w");
      npcf = fopen(NPC_FILE2, "w");
      if (kf == NULL) {
        printf("Could not open kspace file for writing\n");
      }
      if (df == NULL) {
        printf("Could not open density file for writing\n");
      }
      if (npcf == NULL) {
        printf("Could not open nintpc file for writing\n");
      }
    } /* End ktraj_out_flag (Python wrapper: closed early — only fopen is gated;
         SNR + the file-write loops below are individually NULL-guarded). */

      /*write nintpc to file*/
      for (i = 0; i < 2 * ntheta - 1; i++) {
        if (npcf != NULL) {
          tmpout2 = (float)nintpc[i];
          fwrite(&tmpout2, sizeof(float), 1, npcf);
        }
      }

      printf("MCARL: nintpc written\n");
      fflush(stdout);

      snr_eff_num = 0.0;
      snr_eff_den = 0.0;
      snr_eff_den2 = 0.0;
      snr_numpts = 0;

      /* Calculate SNR efficiency and output density and k-space
         trajectory files */
      /* Loop over each actual cone */
      for (i = 0; i < 2 * ntheta - 1; i++) {
        sdens = 0.0;
        sdens2 = 0.0;
        sdens3 = 0.0;
        for (j = 0; j < GRAD_POINTS; j++) {
          dens[i][j] = dens[i][j] / (float)nintpc[i];
          sdens = sdens + (float)dens[i][j];
          sdens2 = sdens2 + (float)dens[i][j] * (float)dens[i][j];
        }
        snr_eff_num = snr_eff_num + sdens * (float)nintpc[i];
        snr_eff_den =
            snr_eff_den + sqrt(sdens2 * (float)GRAD_POINTS) * nintpc[i];
        snr_eff_den2 = snr_eff_den2 + sdens2 * nintpc[i];
        snr_numpts = snr_numpts + (int)(GRAD_POINTS * nintpc[i]);

        /* SNR Efficiency */
        *snr_eff =
            (float)(100.0 * snr_eff_num / sqrt(snr_eff_den2 * snr_numpts));

        /* Write out ktraj information for all interleaves at this angle */
        if (kf != NULL) {
          /* Skip by ratio of data sampling to gradient sampling */
          for (m = 0; m < READ_POINTS; m++) {
            /*dtmpx =
            (k[i][(int)(floor(dataskip_float*m))][0]+k[i][min(GRAD_POINTS-1,(int)(floor(dataskip_float*m+0.5)))][0])/2;
            dtmpy =
            (k[i][(int)(floor(dataskip_float*m))][1]+k[i][min(GRAD_POINTS-1,(int)(floor(dataskip_float*m+0.5)))][1])/2;
            dtmpz =
            (k[i][(int)(floor(dataskip_float*m))][2]+k[i][min(GRAD_POINTS-1,(int)(floor(dataskip_float*m+0.5)))][2])/2;*/
            dtmpx =
                k[i][(int)floor(dataskip_float * m)][0] +
                (dataskip_float * m - floor(dataskip_float * m)) *
                    (k[i][min(GRAD_POINTS - 1, (int)ceil(dataskip_float * m))]
                      [0] -
                     k[i][min(GRAD_POINTS - 1, (int)floor(dataskip_float * m))]
                      [0]);
            dtmpy =
                k[i][(int)floor(dataskip_float * m)][1] +
                (dataskip_float * m - floor(dataskip_float * m)) *
                    (k[i][min(GRAD_POINTS - 1, (int)ceil(dataskip_float * m))]
                      [1] -
                     k[i][min(GRAD_POINTS - 1, (int)floor(dataskip_float * m))]
                      [1]);
            dtmpz =
                k[i][(int)floor(dataskip_float * m)][2] +
                (dataskip_float * m - floor(dataskip_float * m)) *
                    (k[i][min(GRAD_POINTS - 1, (int)ceil(dataskip_float * m))]
                      [2] -
                     k[i][min(GRAD_POINTS - 1, (int)floor(dataskip_float * m))]
                      [2]);
            dtmpz = dtmpz + SlabKz;
            fwrite(&dtmpx, sizeof(float), 1, kf);
            fwrite(&dtmpy, sizeof(float), 1, kf);
            fwrite(&dtmpz, sizeof(float), 1, kf);
          } /* End for k<grad_points */
        }   /* End kf != NULL */
      }     /* End for i<2*ntheta-1 */

      printf("MCARL: ktraj written\n");
      fflush(stdout);

      /* Write out density information */
      float meandens;
      if (df != NULL) {
        for (i = 0; i < 2 * ntheta - 1; i++) {
          meandens = 0;
          for (m = 0; m < READ_POINTS; m++) {
            dtmp =
                dens[i][(int)floor(dataskip_float * m)] +
                (dataskip_float * m - floor(dataskip_float * m)) *
                    (dens[i]
                         [min(GRAD_POINTS - 1, (int)ceil(dataskip_float * m))] -
                     dens[i][min(GRAD_POINTS - 1,
                                 (int)floor(dataskip_float * m))]);
            /*dtmp =
             * (dens[i][(int)(floor(dataskip_float*m))]+dens[i][min(GRAD_POINTS-1,(int)(floor(dataskip_float*m+0.5)))])/2;*/
            fwrite(&dtmp, sizeof(float), 1, df);
            meandens = meandens + dtmp;
          }
          meandens = meandens / (float)READ_POINTS;
          /*printf("MCARL: MCARL: MCARL: meandens = %f\n",meandens);
          fflush(stdout);*/
        }
      }

      /*printf("MCARL: MCARL: MCARL: GRAD_POINTS = %d\n",GRAD_POINTS);
      printf("MCARL: MCARL: MCARL: READ_POINTS = %d\n",READ_POINTS);
      fflush(stdout);*/

      printf("MCARL: dens written\n");
      fflush(stdout);

      // MCARL: begin gnuplot (k-space)
      if ((Cones_Plot_Flag == 2) | (Cones_Plot_Flag == 3)) {
        // write k-space wavefrom file
        int indx;
        char temp_name1[128];
        FILE *asciiwavefile = NULL;
        sprintf(temp_name1, "kspace_waveform.txt");
        asciiwavefile = fopen(temp_name1, "w");
        for (indx = 0; indx < grad_points[conenum[ntheta] - 1]; indx++) {

          fprintf(asciiwavefile,
                  "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                  (cos(rot_flag * 0 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 0 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 0 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 0 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 1 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 1 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 1 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 1 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 2 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 2 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 2 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 2 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 3 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 3 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 3 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 3 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 4 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 4 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 4 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 4 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 5 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 5 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 5 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 5 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 6 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 6 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 6 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 6 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 7 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 7 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 7 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 7 * PI / 4.0) * k[ntheta][indx][1]));
        }
        fclose(asciiwavefile);
        fflush(asciiwavefile);

        // read and plot
        char temp_string2[1024];
        int k1, k2;
        sprintf(temp_string2, "plot ");
        // for (j=0; j<nintpc[ntheta]-1; j++)
        for (j = 0; j < 7; j++) {
          k1 = 2 * (j + 1) - 1;
          k2 = 2 * (j + 1);
          sprintf(temp_string2, "%s%s%d%s%d%s", temp_string2,
                  "\"kspace_waveform.txt\" using ", k1, ":", k2,
                  " title \"k-space-arm\" ,");
        }
        k1 = k1 + 2;
        k2 = k2 + 2;
        sprintf(temp_string2, "%s%s%d%s%d%s", temp_string2,
                "\"kspace_waveform.txt\" using ", k1, ":", k2,
                " title \"k-space-arm\"\n");
        // populate asciiscriptfile
        FILE *asciiscriptfile = NULL;
        char temp_name2[128];
        sprintf(temp_name2, "plot_kspace_waveform.sh");
        asciiscriptfile = fopen(temp_name2, "w");
        char cmd_line[128];
        fprintf(asciiscriptfile, "#!/usr/bin/gnuplot -persist\n");
        fprintf(asciiscriptfile, "set multiplot\n");
        fprintf(asciiscriptfile, "set size 1,1;\n");
        fprintf(asciiscriptfile, "set origin 0.0,0.0;\n");
        fprintf(asciiscriptfile, "set xlabel \"Points\"\n");
        fprintf(asciiscriptfile, temp_string2);

        fclose(asciiscriptfile);
        // execute asciiscriptfile
        sprintf(cmd_line, "chmod u+x %s", temp_name2);
        system(cmd_line);
        sprintf(cmd_line, "%s", temp_name2);
        system(cmd_line);
        // MCARL: end gnuplot
      } // end if (Cones_Plot_Flag==2)

      /* Close output files */
      if (df != NULL)
        fclose(df);
      if (kf != NULL)
        fclose(kf);
      if (npcf != NULL)
        fclose(npcf);
    /* End of formerly-ktraj_out_flag block (closed earlier above). */

    /* Close files and free memory as needed */
    if (output_grad == 1) {
      if (gxf != NULL)
        fclose(gxf);
      if (gyf != NULL)
        fclose(gyf);
      if (gzf != NULL)
        fclose(gzf);
    }

    /* Free memory */
    for (i = 0; i < 2 * ntheta - 1; i++) {
      for (j = 0; j < GRAD_POINTS; j++) {
        if (k[i][j] != NULL)
          /* k[i][j] = NULL; */
          FreeMem(k[i][j]);
      }
      if (k[i] != NULL)
        /*	k[i] = NULL; */
        FreeMem(k[i]);
    }
    FreeMem(k);

    for (i = 0; i < 2 * ntheta - 1; i++) {
      if (dens[i] != NULL)
        /*	dens[i] = NULL; */
        FreeMem(dens[i]);
    }
    FreeMem(dens);
  } /* End ktraj_flag */

  /* Free memory */
  for (i = 0; i < NUMCONES; i++) {
    for (j = 0; j < GRAD_POINTS + nramp; j++) {
      if (grad[i][j] != NULL)
        /* grad[i][j] = NULL; */
        FreeMem(grad[i][j]);
    }
    if (grad[i] != NULL)
      /*      grad[i] = NULL; */
      FreeMem(grad[i]);
  }
  FreeMem(grad);

  /*if (ktraj_flag == 1)
    return totints;
  else
    return 2*estints;*/

  return ntheta;
}

/* VD End */

/* Function:  conegrad
 *
 * Entry function from EPIC PSD code to generate cone gradients. Current
 * implementation is designed to be called from the host side to leverage the
 * faster CPU on the host side for the iterative design. As such it needs to
 * write the gradient waveforms to a file for reading in by the target side
 *
 * Outputs: number of interleaves needed per cone in first parameter
 * Inputs:
 * RES = resolution in xy and z orientations
 * FOV = field of view in xy and z directions
 * NUMCONES = number of cones to design (reusable for multiple angles )
 * GRAD_POINTS = maximum waveform length in sample points
 * GSP = gradient sample interval in sec
 * PRECISION = design precision of number of interleaves
 * DCF = density compensation factor
 * SMAX = max slew for design
 * maxSRewind = max slew for rewinder design
 * GMAX = max gradient amp for design
 * SysMaxRewindG = maximum system gradient amplitude
 * output_grad = flag to control writing out of gradient waveform files
 * ktraj_flag = flag to control calculating kspace trajecotry  or not (only need
 * to
 *   do once at the end since it is a more expensive calculation)
 * ktraj_out_flag = flag to control writing out of k-space trajectory file, very
 * time
 *   consuming so try to only do once just before download
 * endian_flag = flag to control endianness of output
 */

/*MCARL: Added Cones_Plot_Flag*/
int conegrad(int *nintpc, float *rspthetas, float *RES, float *FOV,
             int *gxcones_host, int *gycones_host, int *gzcones_host, /* CMS */
             int NUMCONES, int GRAD_POINTS, float GSP, int READ_POINTS,
             float TSP, float PRECISION, float DCF, int OVERSAMPLE, float SMAX,
             float maxSRewind, float GMAX, float sysGMAX, float SysMaxRewindG,
             int output_grad, int ktraj_flag, int ktraj_out_flag,
             int endian_flag, int rewind_flag, float MINDENS, float *maxgr2,
             float *snr_eff, int *ngradact, float t_fracx, float t_fracy,
             float t_fracz, int acq_mode, float rot_flag, int SymFlag,
             float SlabKz, int rhkacq_uid,
             int Cones_Plot_Flag,
             int *traj_length) /*MCARL: Added t_frac, acq_mode, rot_flag,
                                  rhkacq_uid, Cones_Plot_Flag;
                                  Python wrapper added traj_length. */
{

  float **g;
  float ***k;   /* Make float as want to write to file as floats */
  float **dens; /* Make float as want to write to file as floats */
  float ***grad;
  float tmpgr;
  /* Switching to static allocation, allow up to 1000 thetas, and 100 cones */
  int conenum[MAX_THETAS];
  float xyscale[MAX_THETAS], zscale[MAX_THETAS];
  float thetas[2 * MAX_THETAS];
  float dthetas_ic[2 * MAX_THETAS];
  float thetaFinals[20 * MAX_THETAS];
  float dthetaFinals_ic[20 * MAX_THETAS];
  float dtheta_ic[20 * MAX_THETAS];
  int grad_points[MAX_CONES];
  float fnint[MAX_CONES];
  float fnints[2 * MAX_THETAS];
  float dnints[2 * MAX_THETAS];
  float FOVrad;
  float *truetheta;
  float thetaup, thetadown;
  int length, totlength;
  float nint;
  float dnint;
  float *theta;
  float THETA[2];
  int ntheta1;
  int ntheta7;
  int nthetaFinal;
  int ntheta;
  float FOV1[2], FOV7[2], FOVFinal[2];
  float FOVFactor1 = 1;
  float FOVFactor7 = 7;
  float FOVFactorFinal;
  float MyRatioFinal = 4;
  float MyRatio1;
  float MyRatio7;
  float MyRatioCheck;
  float MySlope;
  float **kRot;

  int nthetatmp;
  int totints;
  int estints;
  int i, j, m;
  int nramp;
  FILE *gxf = NULL, *gyf = NULL, *gzf = NULL, *npcf = NULL;
  int tmpout;
  float tmpout2;
  float maxnints = 0.0;
  float sdens, sdens2, sdens3;
  float snr_eff_num, snr_eff_den;
  int snr_numpts;
  float snr_eff_den2;
  float dtmpx;
  float dtmpy;
  float dtmpz;
  float dtmp;
  /* For gradient rewinder calcs */
  float max_kpos;
  float rewinder_area1;
  int nramp_min;

  *maxgr2 = -1e6; /*MCARL: Need to initialize*/

  /* Default SNR efficiency - can tell if not calculated later */
  *snr_eff = 0;
  estints = 0;
  totints = 0;

  /* The maximum length waveform actually created */
  (*ngradact) = 0; /* initialize to 0 */

  /* Set the number of points in the ramp of gradient back to zero to the max
     number that could be needed */
  nramp = (int)ceil(GMAX / (maxSRewind * GSP) + 5);
  // printf("MCARL: MCARL: MCARL: nramp1 = %d\n",nramp);
  // fflush(stdout);

  /* Calculate number of "extra" gradient points at the end to account for time
     to
     ramp back down to 0 then to apply a rewinder gradient to bring the
     trajectory
     back to k=0. Consider the worst case scenario where for each Gx or Gy
     gradient
     we could be at kx=maxK or ky=maxK. Consider we're at that point, then go
     further
     while we ramp back to zero then need time to play inverse gradient.
  */
  /* First find the furthest excursion in k-space we could be (based on encoded
     resolution and gradient ramping back down to zero
  */
  if (rewind_flag == 1) {
    max_kpos = (float)5 / min(RES[0], RES[1]) +
               0.5 * GAMMA_PROTON * GMAX * GMAX /
                   maxSRewind; /*MCARL: Second term is additional k-sp accrual
                                  during ramp-down*/
    /* Next calculate the area under a rewinder at full amplitude and min
     * duration */
    rewinder_area1 = GAMMA_PROTON * SysMaxRewindG * SysMaxRewindG / maxSRewind +
                     GAMMA_PROTON * SysMaxRewindG *
                         GSP; /*MCARL: Max area achivable during triangle,
                                 including one time-step point of flattop
                                 (second term)*/
    nramp_min = (int)(ceil(GMAX / (maxSRewind * GSP)) +
                      2 * ceil(SysMaxRewindG / (maxSRewind * GSP)) +
                      1); /*MCARL: time points to ramp down maxG plus to make
                             fast triangle to sysGmax*/
    if (rewinder_area1 > max_kpos) { /*Already have sufficent rewinder*/
      nramp = nramp_min;
    } else { /*Add flattop if needed*/
      nramp = nramp_min + (int)ceil((max_kpos - rewinder_area1) /
                                    (GAMMA_PROTON * SysMaxRewindG *
                                     GSP)); /*MCARL: GAMMA_PROTON was missing*/
    }
    nramp = nramp + 5; // For saftey
  }

  // printf("GRAD_POINTS = %d\n",GRAD_POINTS);
  // printf("nramp = %d\n",nramp);
  // fflush(stdout);

  /* Determine thetas needed for prescription */

  nthetatmp =
      (int)((PI / 2) / (1 / max(FOV[0], FOV[1]) / (5 / min(RES[0], RES[1]))));
  theta = (float *)malloc(20 * nthetatmp * sizeof(float));

  if (SymFlag == 1) {
    getthetas2(theta, RES, FOV, &ntheta, &dtheta_ic[0]);
    // Extend to -pi/2:pi/2
    for (i = 0; i < ntheta - 1; i++) {
      thetas[i] = -1 * theta[ntheta - i - 1];
      dthetas_ic[i] = dtheta_ic[ntheta - i - 1];
      rspthetas[i] = -1 * theta[ntheta - i - 1];
    }
    for (i = 0; i < ntheta; i++) {
      thetas[ntheta - 1 + i] = theta[i];
      dthetas_ic[ntheta - 1 + i] = dtheta_ic[i];
      rspthetas[ntheta - 1 + i] = theta[i];
    }
  }

  if (SymFlag == 2) {
    FOV1[0] = FOVFactor1 * (FOV[0]);
    FOV1[1] = FOVFactor1 * (FOV[1]);
    FOV7[0] = FOVFactor7 * (FOV[0]);
    FOV7[1] = FOVFactor7 * (FOV[1]);
    getthetas(RES, FOV1, &ntheta1, &dtheta_ic[0]);
    getthetas(RES, FOV7, &ntheta7, &dtheta_ic[0]);
    ntheta = ntheta1;
    MyRatio1 = (2 * ntheta1) / (2 * ntheta1 - 1.0);
    MyRatio7 = (2 * ntheta7) / (2 * ntheta1 - 1.0);
    MySlope = (MyRatio7 - MyRatio1) / (FOVFactor7 - FOVFactor1);
    FOVFactorFinal = (MyRatioFinal - MyRatio1) / MySlope + FOVFactor1;
    FOVFinal[0] = FOVFactorFinal * (FOV[0]);
    FOVFinal[1] = FOVFactorFinal * (FOV[1]);
    getthetas2(theta, RES, FOVFinal, &nthetaFinal, &dtheta_ic[0]);
    MyRatioCheck = (2 * nthetaFinal) / (2 * ntheta1 - 1.0);
    printf("MCARL: MCARL: MCARL: MyRatioCheck = %f\n", MyRatioCheck);
    fflush(stdout);
    // Extend to -pi/2:pi/2
    for (i = 0; i < nthetaFinal - 1; i++) {
      thetaFinals[i] = -1 * theta[nthetaFinal - i - 1];
      dthetaFinals_ic[i] = dtheta_ic[nthetaFinal - i - 1];
    }
    for (i = 0; i < nthetaFinal; i++) {
      thetaFinals[nthetaFinal - 1 + i] = theta[i];
      dthetaFinals_ic[nthetaFinal - 1 + i] = dtheta_ic[i];
    }
    // Pick Correct Thetas
    for (i = 0; i < 2 * ntheta1 - 1; i++) {
      thetas[i] = thetaFinals[(int)(i * MyRatioFinal)];
      rspthetas[i] = thetaFinals[(int)(i * MyRatioFinal)];
      dthetas_ic[i] = dthetaFinals_ic[(int)(i * MyRatioFinal)];
    }
  }

  FreeMem(theta);

  /* Allocate memory for gradient waveform arrays */
  grad = (float ***)AllocMem(NUMCONES * sizeof(float **));
  for (i = 0; i < NUMCONES; i++) {
    grad[i] = (float **)AllocMem((GRAD_POINTS + nramp) *
                                 sizeof(float *)); /* ML - Not freed */
    for (j = 0; j < GRAD_POINTS + nramp; j++)
      grad[i][j] =
          (float *)AllocMem(3 * sizeof(float)); /* ML - Not freed (lots) */
  }

  /* Python wrapper: oversize the g buffer to handle aggressive undersampling
     where the cones designer's binary search can produce per-cone lengths
     slightly past LEN+nramp. Cost is modest (~256 KB extra per axis for
     typical configs); benefit is no heap corruption when fov_under is small. */
  g = (float **)AllocMem(3 * sizeof(float *));
  g[0] = (float *)AllocMem(((GRAD_POINTS + nramp) * OVERSAMPLE + 65536) * sizeof(float));
  g[1] = (float *)AllocMem(((GRAD_POINTS + nramp) * OVERSAMPLE + 65536) * sizeof(float));
  g[2] = (float *)AllocMem(((GRAD_POINTS + nramp) * OVERSAMPLE + 65536) * sizeof(float));

  char GXCONE_FILE[128];
  char GYCONE_FILE[128];
  char GZCONE_FILE[128];

  sprintf(GXCONE_FILE, "%s/%s.%d", DIR_OUT, FILE_GXCONE, rhkacq_uid);
  sprintf(GYCONE_FILE, "%s/%s.%d", DIR_OUT, FILE_GYCONE, rhkacq_uid);
  sprintf(GZCONE_FILE, "%s/%s.%d", DIR_OUT, FILE_GZCONE, rhkacq_uid);
  /*
#ifdef SIM
  sprintf(GXCONE_FILE,"%s%d","kacq_gxcone_",rhkacq_uid);
  sprintf(GYCONE_FILE,"%s%d","kacq_gycone_",rhkacq_uid);
  sprintf(GZCONE_FILE,"%s%d","kacq_gzcone_",rhkacq_uid);
#else
  sprintf(GXCONE_FILE,"%s%d","/usr/g/psddata/kacq_gxcone_",rhkacq_uid);
  sprintf(GYCONE_FILE,"%s%d","/usr/g/psddata/kacq_gycone_",rhkacq_uid);
  sprintf(GZCONE_FILE,"%s%d","/usr/g/psddata/kacq_gzcone_",rhkacq_uid);
#endif
  */

  /* Open output files if needed */
  if (output_grad == 1) {
    if ((gxf = fopen(GXCONE_FILE, "w")) == NULL) {
      printf("Could not open gxfile for writing\n");
    }
    if ((gyf = fopen(GYCONE_FILE, "w")) == NULL) {
      printf("Could not open gyfile for writing\n");
    }
    if ((gzf = fopen(GZCONE_FILE, "w")) == NULL) {
      printf("Could not open gzfile for writing\n");
    }
  }

  /* Loop over requested number of cones to design and do the design */
  for (i = 0; i < NUMCONES; i++) {

    printf("MCARL: MCARL: MCARL: numcone = %d\n", i);
    fflush(stdout);

    /* Calculate range of thetas for current cone */
    THETA[0] = i * PI / 2 / NUMCONES;
    THETA[1] = (i + 1) * PI / 2 / NUMCONES;

    /*Call findcones for each numcone loop -> THETA = THETA[i]*/
    findcone(&g, &nint, &length, &totlength, RES, FOV, GRAD_POINTS, THETA,
             PRECISION, GSP, OVERSAMPLE, SMAX, maxSRewind, GMAX, SysMaxRewindG,
             DCF, MINDENS, nramp, rewind_flag, acq_mode);
    fnint[i] = nint;

    /*printf("MCARL: MCARL: MCARL: (*ngradact)[%d] = %d\n",i,(*ngradact));
      fflush(stdout);*/

    if (totlength > (*ngradact))
      (*ngradact) = totlength;

    /* Estimate the number of interleaves for reporting to scan UI */
    estints =
        estints + (int)((float)ntheta / (float)NUMCONES) * (int)ceil(nint);

    /* Fill in full gradient waveform array based on current cone design
     * waveforms */
    for (j = 0; j < GRAD_POINTS + nramp; j++) {
      grad[i][j][0] = g[0][j];
      grad[i][j][1] = g[1][j];
      grad[i][j][2] = g[2][j];

    } /* End if j<GRAD_POINTSGTH+nramp */

    grad_points[i] = length;
    /* Python wrapper: also expose this per-cone trajectory length so callers
       can slice readouts at the true end of the spiral, before the rewinder. */
    traj_length[i] = length;
    /*printf("MCARL: MCARL: MCARL: grad_points[%d] = %d\n",i,grad_points[i]);
      fflush(stdout);*/

  } /* End for NUMCONES */

  printf("MCARL: Gradients generated\n");
  fflush(stdout);

  /* JYCHENG: re-order cones here!!!! */

  // MCARL: begin gnuplot (Gradients)
  if ((Cones_Plot_Flag == 1) | (Cones_Plot_Flag == 3)) {
    // write grad wavefrom file
    int indx;
    char temp_name1[128], temp_name2[128];
    FILE *asciiwavefile = NULL;
    sprintf(temp_name1, "grad_waveform.txt");
    asciiwavefile = fopen(temp_name1, "w");
    for (indx = 0; indx < GRAD_POINTS + nramp; indx++)
    /*for (indx = 0; indx < GRAD_POINTS; indx++) */ /*MCARL: Only plot G during DAQ*/
    {
      fprintf(asciiwavefile, "%d %g %g %g\n", indx, 10 * grad[0][indx][0],
              10 * grad[0][indx][1], 10 * grad[0][indx][2]);
    }
    fclose(asciiwavefile);
    fflush(asciiwavefile);
    // populate asciiscriptfile
    FILE *asciiscriptfile = NULL;
    sprintf(temp_name2, "plot_grad_waveform.sh");
    asciiscriptfile = fopen(temp_name2, "w");
    char cmd_line[128];
    fprintf(asciiscriptfile, "#!/usr/bin/gnuplot -persist\n");
    fprintf(asciiscriptfile, "set multiplot\n");
    fprintf(asciiscriptfile, "set size 1,1;\n");
    fprintf(asciiscriptfile, "set origin 0.0,0.0;\n");
    fprintf(asciiscriptfile, "set xlabel \"Points\"\n");
    fprintf(asciiscriptfile, "plot \"grad_waveform.txt\" using 1:2 title "
                             "\"Gx\" ,\"grad_waveform.txt\" using 1:3 title "
                             "\"Gy\",\"grad_waveform.txt\" using 1:4 title "
                             "\"Gz\"\n");
    fclose(asciiscriptfile);
    // execute asciiscriptfile
    sprintf(cmd_line, "chmod u+x %s", temp_name2);
    system(cmd_line);
    sprintf(cmd_line, "%s", temp_name2);
    system(cmd_line);
    // MCARL: end gnuplot
  } // end if (Cones_Plot_Flag==1)

  /* Loop over writing gradient files out outside above loops so we can keep
     track of the
     length of each trajectory and only write out as few total points as
     necessary
  */
  for (i = 0; i < NUMCONES; i++) {
    for (j = 0; j < (*ngradact); j++) {

      /* Write out gradient waveforms to file, convert endianness if needed
         (yes needed if reading on target side) */
      /*MCARL: Gradiens to be read by pulsegen()*/
      /* Conegrad Python wrapper: always fill the host arrays so callers get
         the gradient waveforms back regardless of file-output flag. File
         writing is still gated by output_grad. */
      tmpout = (int)(grad[i][j][0] / sysGMAX * MAX_PG_WAMP);
      gxcones_host[j + i*(*ngradact)] = tmpout;
      if (output_grad == 1 && gxf != NULL) {
        fwrite(&tmpout, sizeof(int), 1, gxf);
      }

      tmpout = (int)(grad[i][j][1] / sysGMAX * MAX_PG_WAMP);
      gycones_host[j + i*(*ngradact)] = tmpout;
      if (output_grad == 1 && gyf != NULL) {
        fwrite(&tmpout, sizeof(int), 1, gyf);
      }

      tmpout = (int)(grad[i][j][2] / sysGMAX * MAX_PG_WAMP);
      gzcones_host[j + i*(*ngradact)] = tmpout;
      if (output_grad == 1 && gzf != NULL) {
        fwrite(&tmpout, sizeof(int), 1, gzf);
      }
    }
  }

  printf("MCARL: After output_grad\n");
  fflush(stdout);

  FreeMem(g[0]);
  FreeMem(g[1]);
  FreeMem(g[2]);
  FreeMem(g);

  printf("MCARL: Gradients written\n");
  fflush(stdout);

  /* ### Generate kspace trajectory (for recon) and density file ### */
  if (ktraj_flag == 1) {
    FILE *kf = NULL, *df = NULL;

    /* Allocate memory for k-space trajectories */
    k = (float ***)AllocMem((2 * ntheta - 1) * sizeof(float **));
    for (i = 0; i < 2 * ntheta - 1; i++) {
      k[i] = (float **)AllocMem(GRAD_POINTS * sizeof(float *));
      for (j = 0; j < GRAD_POINTS; j++)
        k[i][j] = (float *)AllocMem(3 * sizeof(float));
    }
    kRot = (float **)AllocMem(GRAD_POINTS * sizeof(float *));
    for (j = 0; j < GRAD_POINTS; j++)
      kRot[j] = (float *)AllocMem(3 * sizeof(float));

    /* Allocate memory for density information */
    dens = (float **)AllocMem((2 * ntheta - 1) * sizeof(float *));
    truetheta = (float *)AllocMem((2 * ntheta - 1) * sizeof(float));

    /* Loop over all actual cones (vs designed cones) and determine the kspace
       trajectory for each one */
    for (i = 0; i < 2 * ntheta - 1; i++) {

      dens[i] =
          (float *)AllocMem(GRAD_POINTS *
                            sizeof(float)); /*After this, dens is a 2D array
                                               dens[ntheta][GRAD_POINTS], used
                                               later dens[i][j].*/

      /* Which designed cone is appropriate for this actual polar angel */
      conenum[i] = min(
          NUMCONES, max(1, (int)ceil(fabs((float)thetas[i]) / (PI / 2) *
                                     NUMCONES))); /*Will feed
                                                     grad_points[conenum[i]-1]
                                                     into gentraj*/
      thetaup = (float)(conenum[i]) / (float)(NUMCONES)*PI / 2;
      thetadown = (float)(conenum[i] - 1) / (float)(NUMCONES)*PI / 2;
      /* How to scale the designed cone to achieve the current angle */
      xyscale[i] = cos(thetas[i]) / cos(thetadown);
      zscale[i] = sin(thetas[i]) / sin(thetaup);

      truetheta[i] = atan2f(sin(thetas[i]) / RES[1], cos(thetas[i]) / RES[0]);
      FOVrad = 1 / intlineellipse(1 / FOV[0], 1 / FOV[1], truetheta[i]);

      /* Generate the k-space trajectory */
      genktraj(k[i], dens[i], &dnint, grad[conenum[i] - 1],
               grad_points[conenum[i] - 1], GRAD_POINTS, (float)xyscale[i],
               (float)zscale[i], (float)FOVrad, (float)FOV[0],
               (float)truetheta[i], (float)dthetas_ic[i], (float)GSP, &tmpgr);
      if (tmpgr > *maxgr2)
        *maxgr2 = tmpgr;

      dnints[i] =
          dnint; /* Number of interleaves reqd to meet FOV constraints */
      fnints[i] = fnint[conenum[i] -
                        1]; /* Number of interleaves used by design algorithm */
    }                       /* End for i<2*ntheta-1 */

    printf("MCARL: ktraj generated\n");
    fflush(stdout);

    if (MINDENS > 0) {
      for (i = 0; i < 2 * ntheta - 1; i++) {
        if ((dnints[i] / fnints[i]) > maxnints)
          maxnints = dnints[i] / fnints[i];
      }
      for (i = 0; i < 2 * ntheta - 1; i++)
        dnints[i] = fnints[i] * maxnints;
    }

    for (i = 0; i < 2 * ntheta - 1; i++) {
      /* Save the number of interleaves */
      nintpc[i] = (int)ceil(dnints[i]);
      /*MCARL: Off Iso: Round up to nearest odd*/
      if (nintpc[i] > 2) {
        if (SymFlag == 1)
          nintpc[i] = (int)ceil(dnints[i] / 2) * 2;
        if (SymFlag == 2)
          nintpc[i] = (int)ceil(dnints[i] / 2) * 2 - 1;
      }
      /*printf("MCARL: nintpc[%d] = %d\n",i,nintpc[i]);
      fflush(stdout);*/
      totints = totints + nintpc[i];
    }

    /* ### Write out k-space trajectory and density information ### */
    /* Python wrapper: dataskip_float pulled out of the ktraj_out_flag block so
       it stays in scope for the (now-unconditional) SNR + write code below. */
    float dataskip_float = TSP / GSP;
    if (ktraj_out_flag) {

      char KTRAJ_FILE2[128];
      char DENS_FILE2[128];
      char NPC_FILE2[128];

      sprintf(KTRAJ_FILE2, "%s/%s.%d", DIR_OUT, FILE_KTRAJ, rhkacq_uid);
      sprintf(DENS_FILE2, "%s/%s.%d", DIR_OUT, FILE_DENS, rhkacq_uid);
      sprintf(NPC_FILE2, "%s/%s.%d", DIR_OUT, FILE_NPC, rhkacq_uid);
      /*
#ifdef SIM
sprintf(KTRAJ_FILE2,"%s%d","kacq_cone_ktraj_",rhkacq_uid);
sprintf(DENS_FILE2,"%s%d","kacq_cone_dens_",rhkacq_uid);
sprintf(NPC_FILE2,"%s%d","kacq_cone_npc_",rhkacq_uid);
#else
sprintf(KTRAJ_FILE2,"%s%d","/usr/g/psddata/kacq_cone_ktraj_",rhkacq_uid);
sprintf(DENS_FILE2,"%s%d","/usr/g/psddata/kacq_cone_dens_",rhkacq_uid);
sprintf(NPC_FILE2,"%s%d","/usr/g/psddata/kacq_cone_npc_",rhkacq_uid);
#endif
      */

      kf = fopen(KTRAJ_FILE2, "w");
      df = fopen(DENS_FILE2, "w");
      npcf = fopen(NPC_FILE2, "w");
      if (kf == NULL) {
        printf("Could not open kspace file for writing\n");
      }
      if (df == NULL) {
        printf("Could not open density file for writing\n");
      }
      if (npcf == NULL) {
        printf("Could not open nintpc file for writing\n");
      }
    } /* End ktraj_out_flag (Python wrapper: closed early — only fopen is gated;
         SNR + the file-write loops below are individually NULL-guarded). */

      /*write nintpc to file*/
      for (i = 0; i < 2 * ntheta - 1; i++) {
        if (npcf != NULL) {
          tmpout2 = (float)nintpc[i];
          fwrite(&tmpout2, sizeof(float), 1, npcf);
        }
      }

      printf("MCARL: nintpc written\n");
      fflush(stdout);

      snr_eff_num = 0.0;
      snr_eff_den = 0.0;
      snr_eff_den2 = 0.0;
      snr_numpts = 0;

      /* Calculate SNR efficiency and output density and k-space
         trajectory files */
      /* Loop over each actual cone */
      for (i = 0; i < 2 * ntheta - 1; i++) {
        sdens = 0.0;
        sdens2 = 0.0;
        sdens3 = 0.0;
        for (j = 0; j < GRAD_POINTS; j++) {
          dens[i][j] = dens[i][j] / (float)nintpc[i];
          sdens = sdens + (float)dens[i][j];
          sdens2 = sdens2 + (float)dens[i][j] * (float)dens[i][j];
        }
        snr_eff_num = snr_eff_num + sdens * (float)nintpc[i];
        snr_eff_den =
            snr_eff_den + sqrt(sdens2 * (float)GRAD_POINTS) * nintpc[i];
        snr_eff_den2 = snr_eff_den2 + sdens2 * nintpc[i];
        snr_numpts = snr_numpts + (int)(GRAD_POINTS * nintpc[i]);

        /* SNR Efficiency */
        *snr_eff =
            (float)(100.0 * snr_eff_num / sqrt(snr_eff_den2 * snr_numpts));

        /* Write out ktraj information for all interleaves at this angle */
        if (kf != NULL) {
          /* Skip by ratio of data sampling to gradient sampling */
          for (m = 0; m < READ_POINTS; m++) {
            /*dtmpx =
            (k[i][(int)(floor(dataskip_float*m))][0]+k[i][min(GRAD_POINTS-1,(int)(floor(dataskip_float*m+0.5)))][0])/2;
            dtmpy =
            (k[i][(int)(floor(dataskip_float*m))][1]+k[i][min(GRAD_POINTS-1,(int)(floor(dataskip_float*m+0.5)))][1])/2;
            dtmpz =
            (k[i][(int)(floor(dataskip_float*m))][2]+k[i][min(GRAD_POINTS-1,(int)(floor(dataskip_float*m+0.5)))][2])/2;*/
            dtmpx =
                k[i][(int)floor(dataskip_float * m)][0] +
                (dataskip_float * m - floor(dataskip_float * m)) *
                    (k[i][min(GRAD_POINTS - 1, (int)ceil(dataskip_float * m))]
                      [0] -
                     k[i][min(GRAD_POINTS - 1, (int)floor(dataskip_float * m))]
                      [0]);
            dtmpy =
                k[i][(int)floor(dataskip_float * m)][1] +
                (dataskip_float * m - floor(dataskip_float * m)) *
                    (k[i][min(GRAD_POINTS - 1, (int)ceil(dataskip_float * m))]
                      [1] -
                     k[i][min(GRAD_POINTS - 1, (int)floor(dataskip_float * m))]
                      [1]);
            dtmpz =
                k[i][(int)floor(dataskip_float * m)][2] +
                (dataskip_float * m - floor(dataskip_float * m)) *
                    (k[i][min(GRAD_POINTS - 1, (int)ceil(dataskip_float * m))]
                      [2] -
                     k[i][min(GRAD_POINTS - 1, (int)floor(dataskip_float * m))]
                      [2]);
            dtmpz = dtmpz + SlabKz;
            fwrite(&dtmpx, sizeof(float), 1, kf);
            fwrite(&dtmpy, sizeof(float), 1, kf);
            fwrite(&dtmpz, sizeof(float), 1, kf);
          } /* End for k<grad_points */
        }   /* End kf != NULL */
      }     /* End for i<2*ntheta-1 */

      printf("MCARL: ktraj written\n");
      fflush(stdout);

      /* Write out density information */
      float meandens;
      if (df != NULL) {
        for (i = 0; i < 2 * ntheta - 1; i++) {
          meandens = 0;
          for (m = 0; m < READ_POINTS; m++) {
            dtmp =
                dens[i][(int)floor(dataskip_float * m)] +
                (dataskip_float * m - floor(dataskip_float * m)) *
                    (dens[i]
                         [min(GRAD_POINTS - 1, (int)ceil(dataskip_float * m))] -
                     dens[i][min(GRAD_POINTS - 1,
                                 (int)floor(dataskip_float * m))]);
            /*dtmp =
             * (dens[i][(int)(floor(dataskip_float*m))]+dens[i][min(GRAD_POINTS-1,(int)(floor(dataskip_float*m+0.5)))])/2;*/
            fwrite(&dtmp, sizeof(float), 1, df);
            meandens = meandens + dtmp;
          }
          meandens = meandens / (float)READ_POINTS;
          /*printf("MCARL: MCARL: MCARL: meandens = %f\n",meandens);
          fflush(stdout);*/
        }
      }

      /*printf("MCARL: MCARL: MCARL: GRAD_POINTS = %d\n",GRAD_POINTS);
      printf("MCARL: MCARL: MCARL: READ_POINTS = %d\n",READ_POINTS);
      fflush(stdout);*/

      printf("MCARL: dens written\n");
      fflush(stdout);

      // MCARL: begin gnuplot (k-space)
      if ((Cones_Plot_Flag == 2) | (Cones_Plot_Flag == 3)) {
        // write k-space wavefrom file
        int indx;
        char temp_name1[128];
        FILE *asciiwavefile = NULL;
        sprintf(temp_name1, "kspace_waveform.txt");
        asciiwavefile = fopen(temp_name1, "w");
        for (indx = 0; indx < grad_points[conenum[ntheta] - 1]; indx++) {

          fprintf(asciiwavefile,
                  "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
                  (cos(rot_flag * 0 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 0 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 0 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 0 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 1 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 1 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 1 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 1 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 2 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 2 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 2 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 2 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 3 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 3 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 3 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 3 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 4 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 4 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 4 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 4 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 5 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 5 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 5 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 5 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 6 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 6 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 6 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 6 * PI / 4.0) * k[ntheta][indx][1]),
                  (cos(rot_flag * 7 * PI / 4.0) * k[ntheta][indx][0] -
                   sin(rot_flag * 7 * PI / 4.0) * k[ntheta][indx][1]),
                  (sin(rot_flag * 7 * PI / 4.0) * k[ntheta][indx][0] +
                   cos(rot_flag * 7 * PI / 4.0) * k[ntheta][indx][1]));
        }
        fclose(asciiwavefile);
        fflush(asciiwavefile);

        // read and plot
        char temp_string2[1024];
        int k1, k2;
        sprintf(temp_string2, "plot ");
        // for (j=0; j<nintpc[ntheta]-1; j++)
        for (j = 0; j < 7; j++) {
          k1 = 2 * (j + 1) - 1;
          k2 = 2 * (j + 1);
          sprintf(temp_string2, "%s%s%d%s%d%s", temp_string2,
                  "\"kspace_waveform.txt\" using ", k1, ":", k2,
                  " title \"k-space-arm\" ,");
        }
        k1 = k1 + 2;
        k2 = k2 + 2;
        sprintf(temp_string2, "%s%s%d%s%d%s", temp_string2,
                "\"kspace_waveform.txt\" using ", k1, ":", k2,
                " title \"k-space-arm\"\n");
        // populate asciiscriptfile
        FILE *asciiscriptfile = NULL;
        char temp_name2[128];
        sprintf(temp_name2, "plot_kspace_waveform.sh");
        asciiscriptfile = fopen(temp_name2, "w");
        char cmd_line[128];
        fprintf(asciiscriptfile, "#!/usr/bin/gnuplot -persist\n");
        fprintf(asciiscriptfile, "set multiplot\n");
        fprintf(asciiscriptfile, "set size 1,1;\n");
        fprintf(asciiscriptfile, "set origin 0.0,0.0;\n");
        fprintf(asciiscriptfile, "set xlabel \"Points\"\n");
        fprintf(asciiscriptfile, temp_string2);

        fclose(asciiscriptfile);
        // execute asciiscriptfile
        sprintf(cmd_line, "chmod u+x %s", temp_name2);
        system(cmd_line);
        sprintf(cmd_line, "%s", temp_name2);
        system(cmd_line);
        // MCARL: end gnuplot
      } // end if (Cones_Plot_Flag==2)

      /* Close output files */
      if (df != NULL)
        fclose(df);
      if (kf != NULL)
        fclose(kf);
      if (npcf != NULL)
        fclose(npcf);
    /* End of formerly-ktraj_out_flag block (closed earlier above). */

    /* Close files and free memory as needed */
    if (output_grad == 1) {
      if (gxf != NULL)
        fclose(gxf);
      if (gyf != NULL)
        fclose(gyf);
      if (gzf != NULL)
        fclose(gzf);
    }

    /* Free memory */
    for (i = 0; i < 2 * ntheta - 1; i++) {
      for (j = 0; j < GRAD_POINTS; j++) {
        if (k[i][j] != NULL)
          /* k[i][j] = NULL; */
          FreeMem(k[i][j]);
      }
      if (k[i] != NULL)
        /*	k[i] = NULL; */
        FreeMem(k[i]);
    }
    FreeMem(k);

    for (i = 0; i < 2 * ntheta - 1; i++) {
      if (dens[i] != NULL)
        /*	dens[i] = NULL; */
        FreeMem(dens[i]);
    }
    FreeMem(dens);
  } /* End ktraj_flag */

  /* Free memory */
  for (i = 0; i < NUMCONES; i++) {
    for (j = 0; j < GRAD_POINTS + nramp; j++) {
      if (grad[i][j] != NULL)
        /* grad[i][j] = NULL; */
        FreeMem(grad[i][j]);
    }
    if (grad[i] != NULL)
      /*      grad[i] = NULL; */
      FreeMem(grad[i]);
  }
  FreeMem(grad);

  /*if (ktraj_flag == 1)
    return totints;
  else
    return 2*estints;*/

  return ntheta;
}

double ceil2(double x) {
  double y;

  y = (double)(int)x;

  if (wcfEqual((float)x-y, 0.0)) {
    return x;
  } else {
    return (double)(int)(x + 1.0);
  }
} /* end ceil2() */

/*
 * Circularly permute using the golden ratio
 *   indx_out -- output permutation array
 *   n        -- (n-1) is the max permutation index
 *   k        -- length of output array
 *   norep    -- if true, skip optimal index if already
 *               in the indx_out array
 */
int permuteCircGoldenCones(int indx_out[], int n, int k, unsigned short norep) {
  int i; /* for loop counter */

  float ai = 0.0;                               /* current angle */
  float gra = 360.0f * (sqrt(5) - 1.0f) / 2.0f; /* golden angle */
  float da = 360.0f / n;                        /* delta angle */

  int ip = 0; /* permutation element */
  int io = 0; /* current output index */

  unsigned short *exist; /* array to see if already acq */
  int existi = 0;        /* exist-array counter */
  if ((exist = AllocMem(n * sizeof(*exist))) == NULL) {
    printf("permuteCircGoldenCones: Can't allocate exist!\n");
    return FAILURE;
  }

  /* init exist array */
  for (i = 0; i < n; i++)
    exist[i] = 0;

  while (io < k) {
    ip = (int)(ai / da + 0.5);
    if (ip >= n)
      ip = 0; /* in case, we round up on last element */
    if (0) {
      char tmpstr[128];
      sprintf(tmpstr, "Current index: %d, current angle: %g, %g\n", ip, ai, da);
      fputs(tmpstr, stderr);
    }

    if (!exist[ip] || !norep) {
      indx_out[io] = ip;
      /* update exist array and counter */
      exist[ip] = 1;
      existi++;
      /* update output array index */
      io++;
    }
    ai += gra;
    if (ai >= 360)
      ai -= 360; /* unwrap angle */

    /* reset counter if finished */
    if (norep && existi >= n) {
      for (i = 0; i < n; i++)
        exist[i] = 0;
    }
  }

  FreeMem(exist);

  return SUCCESS;
} /* permuteCircGoldenCones */

/* permutePhyllo(phyPhi, sindx_p, totints, cpInter, indToTheta_i, 0); */
/*   MOMALAVE
 *   Phyllotaxis permutaiton (still needs unique elevation angles increasing
 *inter # for full phyllo)
 *	 phyTheta    -- output permutation array
 *   totints     -- (n-1) is the max permutation index
 *   cpInter     -- cones per interelave
 *   inToTheta_i -- Thetas (conenum??)
 *   phyFullFlag -- flag to use fib number to minimize distance between
 *interleaves (ToDo: increase number of interleaves for full phyllo)
 */
int permutePhyllo(int phyTheta_out[], int phyindToInt_out[], int totints,
                  int *cpInter, int indToTheta_i[], int indToInt_i[],
                  int phyFullFlag) {
  /* float gra = PI * (3 - sqrt(5)); */ /* golden angle */
  int il, cone, actcone, duflag;
  int ilNum = 0;
  int ilRemain = 0;
  int count = 0;
  int cpInterOut = *cpInter;

  int fibArray[2] = {0, 1};
  int temp = 0;
  if (phyFullFlag) {
    while (cpInterOut < totints / fibArray[1]) {
      temp = fibArray[0] + fibArray[1];
      fibArray[0] = fibArray[1];
      fibArray[1] = temp;
    }
    cpInterOut = (int)totints / fibArray[1];
  }

  ilNum = totints / cpInterOut;
  ilRemain = totints % cpInterOut; /* remaining interleaves */

  char tmpstr[80];
  sprintf(tmpstr,
          "total: %d, Rotations: %d, Remain: %d, cpInter: %d fibNum: %d \n",
          totints, ilNum, ilRemain, cpInterOut, ilNum);
  fputs(tmpstr, stderr);
  fflush(stderr);

  /* count = 0;
  for (j = 0; j < ilNum; j++){
          for (i = 0; i < cpInterOut; i++) {
                  phyTheta_out[count] = indToTheta_i[(ilNum*i) + j];
                  phyindToInt_out[count] = (ilNum*i) + j;
                  count = count + 1;
          }
  } */

  /* CMS: In the code above, each phyllotaxis interleaf resets back to the top
     cone after finishing. In this version, each interleaf is continuous with
     the next (i.e. down-up-down-up-...) */
  duflag = 0;
  count = 0;
  for (il = 0; il < ilNum; il++) {
    for (cone = 0; cone < cpInterOut; cone++) {
      actcone = (duflag) ? (cpInterOut - cone - 1) : cone;
      phyTheta_out[count] = indToTheta_i[(ilNum * actcone) + il];
      phyindToInt_out[count] = (ilNum * actcone) + il;
      count = count + 1;
    }
    duflag = (duflag + 1) % 2;
  }

  /* since we are not full phyllo we need to extract the remaining interleaves
   */
  for (cone = 0; cone < ilRemain; cone++) {
    actcone = (duflag) ? (ilRemain - cone - 1) : cone;
    phyTheta_out[count] = indToTheta_i[(ilNum * actcone) + ilNum];
    phyindToInt_out[count] = (ilNum * actcone) + ilNum;
    count = count + 1;
  }
  /* printf("count: %d \n",count); */

  *cpInter = cpInterOut;

  return SUCCESS;
} /* permutePhyllo */
