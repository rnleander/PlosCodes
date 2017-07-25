/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzlartg.h
 *
 * Code generation for function 'xzlartg'
 *
 */

#ifndef XZLARTG_H
#define XZLARTG_H

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
extern void b_xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn);
extern void xzlartg(const creal_T f, const creal_T g, double *cs, creal_T *sn,
                    creal_T *r);

#endif

/* End of code generation (xzlartg.h) */
