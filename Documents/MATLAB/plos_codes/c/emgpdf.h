/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * emgpdf.h
 *
 * Code generation for function 'emgpdf'
 *
 */

#ifndef EMGPDF_H
#define EMGPDF_H

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
#ifdef _OLD_MATLAB_CODE
extern void emgpdf_old(const double X[266], double l, double m, double s, double Y[266]);
#endif

extern double emgpdf_loglikelihood(const gsl_vector *v, void *params);

extern void emgpdf(const double X[266], double l, double m, double s, double Y[266]);


#endif

/* End of code generation (emgpdf.h) */
