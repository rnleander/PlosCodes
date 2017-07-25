/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzhgeqz.h
 *
 * Code generation for function 'xzhgeqz'
 *
 */

#ifndef XZHGEQZ_H
#define XZHGEQZ_H

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
extern void xzhgeqz(const creal_T A_data[], const int A_size[2], int ilo, int
                    ihi, int *info, creal_T alpha1_data[], int alpha1_size[1],
                    creal_T beta1_data[], int beta1_size[1]);

#endif

/* End of code generation (xzhgeqz.h) */
