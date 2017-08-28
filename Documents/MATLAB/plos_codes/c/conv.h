/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * conv.h
 *
 * Code generation for function 'conv'
 *
 */

#ifndef CONV_H
#define CONV_H

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
extern void b_conv(const emxArray_real_T *A, const emxArray_real_T *B,
                   emxArray_real_T *C);
extern void c_conv(const double A[22001], const double B[22001], double C[44001]);
extern void conv(const double A[2201], const double B[2201], double C[4401]);
extern void d_conv(const double A[221], const double B[221], double C[441]);
extern void e_conv(const double A[221], const double B[221], double C[441]);
extern void f_conv(const emxArray_real_T *A, const emxArray_real_T *B,
                   emxArray_real_T *C);
extern void g_conv(const double A[2201], const double B[2201], double C[4401]);
extern void h_conv(const double A[22001], const double B[22001], double C[44001]);
#endif

#endif

/* End of code generation (conv.h) */
