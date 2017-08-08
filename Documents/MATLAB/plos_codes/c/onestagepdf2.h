/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * onestagepdf2.h
 *
 * Code generation for function 'onestagepdf2'
 *
 */

#ifndef ONESTAGEPDF2_H
#define ONESTAGEPDF2_H

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
extern void b_onestagepdf2(const double t[2201], double mu, double s, double Y
  [2201]);
extern double c_onestagepdf2(double t, double mu, double s);
extern void d_onestagepdf2(const emxArray_real_T *t, double mu, double s,
  emxArray_real_T *Y);
extern void e_onestagepdf2(const double t[22001], double mu, double s, double Y
  [22001]);
extern void f_onestagepdf2(const double t[221], double mu, double s, double Y
  [221]);
extern void g_onestagepdf2(const double t[221], double mu, double s, double Y
  [221]);
extern void h_onestagepdf2(const emxArray_real_T *t, double mu, double s,
  emxArray_real_T *Y);
extern void i_onestagepdf2(const double t[2201], double mu, double s, double Y
  [2201]);
extern void j_onestagepdf2(const double t[22001], double mu, double s, double Y
  [22001]);

extern void waldpdf(const double X[], double mu, double s, double Y[], int size_XY);


extern double wald_loglikelihood(const gsl_vector *v, void *params);


extern void onestagepdf2(const double t[266], double mu, double s, double Y[266]);

#endif

/* End of code generation (onestagepdf2.h) */
