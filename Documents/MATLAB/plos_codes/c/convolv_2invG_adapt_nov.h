/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * convolv_2invG_adapt_nov.h
 *
 * Code generation for function 'convolv_2invG_adapt_nov'
 *
 */

#ifndef CONVOLV_2INVG_ADAPT_NOV_H
#define CONVOLV_2INVG_ADAPT_NOV_H

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
extern double convolv_2invG_adapt_nov_loglikelihood(const gsl_vector *v, void *params);

extern void conv2waldpdf(const double X[], double m1, double s1, double m2, double s2, double Y[], double h, int adaptiveMode, int size_XY);


#ifdef _OLD_MATLAB_CODE
extern void b_convolv_2invG_adapt_nov(double m1, double s1, double m2, double s2,
  double P[266], double *h, double *flag, double *E);
extern void c_convolv_2invG_adapt_nov(double m1, double s1, double m2, double s2,
  double P[266], double *h, double *flag, double *E);
extern void convolv_2invG_adapt_nov(double m1, double s1, double m2, double s2,
  double P[266]);
extern void d_convolv_2invG_adapt_nov(const double t[266], double m1, double s1,
  double m2, double s2, double P[266]);
extern void e_convolv_2invG_adapt_nov(const double t[266], double m1, double s1,
  double m2, double s2, double P[266]);
extern void f_convolv_2invG_adapt_nov(const double t[266], double m1, double s1,
  double m2, double s2, double P[266]);

#endif

#endif

/* End of code generation (convolv_2invG_adapt_nov.h) */
