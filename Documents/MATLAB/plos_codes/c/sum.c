/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sum.c
 *
 * Code generation for function 'sum'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "sum.h"

/* Function Definitions */

/*
 *
 */
double b_sum(const emxArray_real_T *x)
{
  double y;
  int k;
  if (x->size[1] == 0) {
    y = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= x->size[1]; k++) {
      y += x->data[k - 1];
    }
  }

  return y;
}

/*
 *
 */
double sum(const double x[266])
{
  double y;
  int k;
  y = x[0];
  for (k = 0; k < 265; k++) {
    y += x[k + 1];
  }

  return y;
}

/* End of code generation (sum.c) */
