/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * convolv_2invG_nov.h
 *
 * Code generation for function 'convolv_2invG_nov'
 *
 */

#ifndef CONVOLV_2INVG_NOV_H
#define CONVOLV_2INVG_NOV_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "omp.h"
#include "IMT_analysis_April2017_types.h"

/* Function Declarations */
#ifdef _OLD_MATLAB_CODE
extern void b_convolv_2invG_nov(const emxArray_real_T *t, double m1, double s1,
  double m2, double s2, double h, emxArray_real_T *P, double *flag);
extern void c_convolv_2invG_nov(double m1, double s1, double m2, double s2,
  double P[2201], double *flag);
extern void convolv_2invG_nov(double m1, double s1, double m2, double s2, double
  P[221], double *flag);
extern void d_convolv_2invG_nov(const double t[22001], double m1, double s1,
  double m2, double s2, double P[22001], double *flag);
#endif

#endif

/* End of code generation (convolv_2invG_nov.h) */
