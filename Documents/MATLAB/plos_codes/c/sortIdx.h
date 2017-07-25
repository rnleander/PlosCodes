/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sortIdx.h
 *
 * Code generation for function 'sortIdx'
 *
 */

#ifndef SORTIDX_H
#define SORTIDX_H

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
extern void b_sortIdx(const double x[3], int idx[3]);
extern void c_sortIdx(const double x[5], int idx[5]);
extern void d_sortIdx(const double x[7], int idx[7]);
extern void sortIdx(const double x[4], int idx[4]);

#endif

/* End of code generation (sortIdx.h) */
