/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * abs.c
 *
 * Code generation for function 'abs'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "abs.h"

/* Function Definitions */

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void b_abs(const double x[2], double y[2])
{
  int k;
  for (k = 0; k < 2; k++) {
    y[k] = fabs(x[k]);
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void c_abs(const double x[3], double y[3])
{
  int k;
  for (k = 0; k < 3; k++) {
    y[k] = fabs(x[k]);
  }
}
#endif
/* End of code generation (abs.c) */
