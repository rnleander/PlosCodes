/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * onestagepdf_prime.h
 *
 * Code generation for function 'onestagepdf_prime'
 *
 */

#ifndef ONESTAGEPDF_PRIME_H
#define ONESTAGEPDF_PRIME_H

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
extern void onestagepdf_prime(const creal_T t_data[], const int t_size[1],
  double m, double s, creal_T Y_data[], int Y_size[1]);

extern void onestagepdf_prime_fixed(const double t[], int t_size, double m, double s, double Y[]);


#endif

/* End of code generation (onestagepdf_prime.h) */
