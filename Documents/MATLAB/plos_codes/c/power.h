/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * power.h
 *
 * Code generation for function 'power'
 *
 */

#ifndef POWER_H
#define POWER_H

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
extern void b_power(const double a[2], double y[2]);
extern void c_power(const double a[2], double y[2]);
extern void d_power(const double a[2201], double y[2201]);
extern void e_power(const emxArray_real_T *a, emxArray_real_T *y);
extern void f_power(const double a[22001], double y[22001]);
extern void g_power(const double a[3], double y[3]);
extern void h_power(const double a[3], double y[3]);
extern void i_power(const double a[3], double y[3]);
extern void j_power(const double a[221], double y[221]);
extern void k_power(const emxArray_real_T *a, emxArray_real_T *y);
extern void l_power(const emxArray_real_T *a, emxArray_real_T *y);
extern void m_power(const emxArray_real_T *a, emxArray_real_T *y);
extern void power(const double a[2], double y[2]);

#endif

/* End of code generation (power.h) */
