/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sortLE.c
 *
 * Code generation for function 'sortLE'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "sortLE.h"

/* Function Definitions */

/*
 *
 */
boolean_T sortLE(const double v[162], const int col[2], int irow1, int irow2)
{
  boolean_T p;
  int k;
  boolean_T exitg1;
  int coloffset;
  p = true;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 2)) {
    coloffset = (col[k] - 1) * 81 - 1;
    if (!(v[coloffset + irow1] == v[coloffset + irow2])) {
      if ((v[coloffset + irow1] <= v[coloffset + irow2]) || rtIsNaN(v[coloffset
           + irow2])) {
        p = true;
      } else {
        p = false;
      }

      exitg1 = true;
    } else {
      k++;
    }
  }

  return p;
}

/* End of code generation (sortLE.c) */
