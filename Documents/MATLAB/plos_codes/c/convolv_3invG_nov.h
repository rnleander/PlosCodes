/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * convolv_3invG_nov.h
 *
 * Code generation for function 'convolv_3invG_nov'
 *
 */

#ifndef CONVOLV_3INVG_NOV_H
#define CONVOLV_3INVG_NOV_H

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
#include "gsl/gsl_multimin.h"


/* Function Declarations */
#ifdef _OLD_MATLAB_CODE
extern void b_convolv_3invG_nov(double m1, double s1, double m2, double s2,
  double m3, double s3, double P[266], double *h, double *flag, double *E);
extern void c_convolv_3invG_nov(double m1, double s1, double m2, double s2,
  double m3, double s3, double P[266], double *h, double *flag, double *E);
#endif

#ifdef _OLD_MATLAB_CODE
extern void convolv_3invG_nov(double m1, double s1, double m2, double s2, double
  m3, double s3, double P[266]);
#endif

extern double convolv_3invG_nov_loglikelihood(const gsl_vector *v, void *params);

extern void convolv3waldpdf(double m1, double s1, double m2, double s2, double m3, double s3, const double X[266], double Y[266], int size_XY, double h);


#endif

/* End of code generation (convolv_3invG_nov.h) */
