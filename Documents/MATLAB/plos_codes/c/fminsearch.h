/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * fminsearch.h
 *
 * Code generation for function 'fminsearch'
 *
 */

#ifndef FMINSEARCH_H
#define FMINSEARCH_H

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
extern double b_fminsearch(double x[2]);
extern double c_fminsearch(double x[3]);
extern double d_fminsearch(double x[4]);
extern double e_fminsearch(double x[6]);
extern double fminsearch(double x[3]);
extern double fminsearch_generalized(double x[3], void(*pdf)(), const double data[266]);
#endif

#endif

/* End of code generation (fminsearch.h) */
