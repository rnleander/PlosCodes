/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * log.c
 *
 * Code generation for function 'log'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "log.h"

/* Function Definitions */

/*
 *
 */
void b_log(double x[266])
{
  int k;
  for (k = 0; k < 266; k++) {
    x[k] = log(x[k]);
  }
}

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void c_log(emxArray_real_T *x)
{
  int nx;
  int k;
  nx = x->size[0];
  for (k = 0; k + 1 <= nx; k++) {
    x->data[k] = log(x->data[k]);
  }
}
#endif
/* End of code generation (log.c) */
