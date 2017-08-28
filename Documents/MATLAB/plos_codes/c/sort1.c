/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sort1.c
 *
 * Code generation for function 'sort1'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "sort1.h"

/* Function Definitions */

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void b_sort(double x[3], int idx[3])
{
  boolean_T p;
  double tmp;
  if ((x[0] <= x[1]) || rtIsNaN(x[1])) {
    p = true;
  } else {
    p = false;
  }

  if (p) {
    if ((x[1] <= x[2]) || rtIsNaN(x[2])) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      idx[0] = 1;
      idx[1] = 2;
      idx[2] = 3;
    } else {
      if ((x[0] <= x[2]) || rtIsNaN(x[2])) {
        p = true;
      } else {
        p = false;
      }

      if (p) {
        idx[0] = 1;
        idx[1] = 3;
        idx[2] = 2;
        tmp = x[1];
        x[1] = x[2];
        x[2] = tmp;
      } else {
        idx[0] = 3;
        idx[1] = 1;
        idx[2] = 2;
        tmp = x[2];
        x[2] = x[1];
        x[1] = x[0];
        x[0] = tmp;
      }
    }
  } else {
    if ((x[0] <= x[2]) || rtIsNaN(x[2])) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      idx[0] = 2;
      idx[1] = 1;
      idx[2] = 3;
      tmp = x[0];
      x[0] = x[1];
      x[1] = tmp;
    } else {
      if ((x[1] <= x[2]) || rtIsNaN(x[2])) {
        p = true;
      } else {
        p = false;
      }

      if (p) {
        idx[0] = 2;
        idx[1] = 3;
        idx[2] = 1;
        tmp = x[0];
        x[0] = x[1];
        x[1] = x[2];
        x[2] = tmp;
      } else {
        idx[0] = 3;
        idx[1] = 2;
        idx[2] = 1;
        tmp = x[0];
        x[0] = x[2];
        x[2] = tmp;
      }
    }
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void sort(double x[2], int idx[2])
{
  boolean_T p;
  double tmp;
  if ((x[0] <= x[1]) || rtIsNaN(x[1])) {
    p = true;
  } else {
    p = false;
  }

  if (p) {
    idx[0] = 1;
    idx[1] = 2;
  } else {
    idx[0] = 2;
    idx[1] = 1;
    tmp = x[0];
    x[0] = x[1];
    x[1] = tmp;
  }
}
#endif

/* End of code generation (sort1.c) */
