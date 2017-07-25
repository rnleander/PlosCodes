/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * IMT_analysis_April2017_initialize.c
 *
 * Code generation for function 'IMT_analysis_April2017_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "IMT_analysis_April2017_initialize.h"
#include "IMT_analysis_April2017_data.h"

/* Function Definitions */
void IMT_analysis_April2017_initialize(void)
{
  rt_InitInfAndNaN(8U);
  omp_init_nest_lock(&emlrtNestLockGlobal);
}

/* End of code generation (IMT_analysis_April2017_initialize.c) */
