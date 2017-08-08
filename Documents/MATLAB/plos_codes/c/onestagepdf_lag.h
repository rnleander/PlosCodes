/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * onestagepdf_lag.h
 *
 * Code generation for function 'onestagepdf_lag'
 *
 */

#ifndef ONESTAGEPDF_LAG_H
#define ONESTAGEPDF_LAG_H

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
extern void b_onestagepdf_lag(const double X[221], double m, double s, double l,
  double Y[221]);
extern void c_onestagepdf_lag(const double X[2201], double m, double s, double l,
  double Y[2201]);
extern void d_onestagepdf_lag(const double X[22001], double m, double s, double
  l, double Y[22001]);

extern double waldlag_loglikelihood(const gsl_vector *v, void *params);

extern void waldlagpdf(const double X[266], double mu, double s, double l, double Y[266]);


extern void onestagepdf_lag(const double X[266], double m, double s, double l,
  double Y[266]);

#endif

/* End of code generation (onestagepdf_lag.h) */
