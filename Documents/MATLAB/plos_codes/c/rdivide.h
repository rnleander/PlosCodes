/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * rdivide.h
 *
 * Code generation for function 'rdivide'
 *
 */

#ifndef RDIVIDE_H
#define RDIVIDE_H

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
extern void b_rdivide(const emxArray_real_T *y, emxArray_real_T *z);
extern void c_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
                      emxArray_real_T *z);
extern void rdivide(const creal_T x_data[], const int x_size[1], const creal_T
                    y_data[], creal_T z_data[], int z_size[1]);
#endif

#endif

/* End of code generation (rdivide.h) */
