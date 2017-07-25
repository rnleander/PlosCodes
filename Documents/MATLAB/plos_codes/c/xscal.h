/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xscal.h
 *
 * Code generation for function 'xscal'
 *
 */

#ifndef XSCAL_H
#define XSCAL_H

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
extern void b_xscal(int n, const creal_T a, creal_T x_data[], int ix0, int incx);
extern void xscal(int n, const creal_T a, creal_T x_data[], int ix0);

#endif

/* End of code generation (xscal.h) */
