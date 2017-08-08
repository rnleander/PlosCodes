/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * approxconvolv.h
 *
 * Code generation for function 'approxconvolv'
 *
 */

#ifndef APPROXCONVOLV_H
#define APPROXCONVOLV_H

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
extern void approxconvolv(const double z[2201], const double y[2201], const
  double t[266], const double x[2201], double P0[266], double *logP0);
extern void b_approxconvolv(const emxArray_real_T *z, const emxArray_real_T *y,
  double h, const emxArray_real_T *x, double P0[266], double *logP0);
extern void c_approxconvolv(const double z[22001], const double y[22001], const
  double t[266], const double x[22001], double P0[266], double *logP0);
extern void d_approxconvolv(const emxArray_real_T *z, const emxArray_real_T *y,
  const double t[266], const emxArray_real_T *x, double P0[266], double *logP0);
extern void e_approxconvolv(const emxArray_real_T *z, const emxArray_real_T *y,
  double h, const double t[266], const emxArray_real_T *x, double P0[266],
  double *logP0);
extern void f_approxconvolv(const emxArray_real_T *z, const emxArray_real_T *y,
  const double t[266], const emxArray_real_T *x, double P0[266], double *logP0);
extern void g_approxconvolv(const emxArray_real_T *z, const emxArray_real_T *y,
  const double t[266], const emxArray_real_T *x, double P0[266], double *logP0);


extern void approxconvolv_replacement(const double z[], const double y[], const double X
	[], const double x[], double Y[], double *logP0, int size_xyz, int size_XY, double h);

#endif

/* End of code generation (approxconvolv.h) */
