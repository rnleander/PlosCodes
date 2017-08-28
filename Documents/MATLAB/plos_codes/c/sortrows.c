/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sortrows.c
 *
 * Code generation for function 'sortrows'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "sortrows.h"
#include "sortLE.h"

/* Function Definitions */

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void sortrows(double y[162], double ndx[81])
{
  int k;
  int idx[81];
  int col[2];
  int i;
  int j;
  int i2;
  int pEnd;
  double ycol[81];
  int p;
  int q;
  int qEnd;
  int kEnd;
  int iwork[81];
  for (k = 0; k < 2; k++) {
    col[k] = k + 1;
  }

  memset(&idx[0], 0, 81U * sizeof(int));
  for (k = 0; k <= 78; k += 2) {
    if (sortLE(y, col, k + 1, k + 2)) {
      idx[k] = k + 1;
      idx[k + 1] = k + 2;
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
  }

  idx[80] = 81;
  i = 2;
  while (i < 81) {
    i2 = i << 1;
    j = 1;
    for (pEnd = 1 + i; pEnd < 82; pEnd = qEnd + i) {
      p = j;
      q = pEnd;
      qEnd = j + i2;
      if (qEnd > 82) {
        qEnd = 82;
      }

      k = 0;
      kEnd = qEnd - j;
      while (k + 1 <= kEnd) {
        if (sortLE(y, col, idx[p - 1], idx[q - 1])) {
          iwork[k] = idx[p - 1];
          p++;
          if (p == pEnd) {
            while (q < qEnd) {
              k++;
              iwork[k] = idx[q - 1];
              q++;
            }
          }
        } else {
          iwork[k] = idx[q - 1];
          q++;
          if (q == qEnd) {
            while (p < pEnd) {
              k++;
              iwork[k] = idx[p - 1];
              p++;
            }
          }
        }

        k++;
      }

      for (k = 0; k + 1 <= kEnd; k++) {
        idx[(j + k) - 1] = iwork[k];
      }

      j = qEnd;
    }

    i = i2;
  }

  for (j = 0; j < 2; j++) {
    for (i = 0; i < 81; i++) {
      ycol[i] = y[(idx[i] + 81 * j) - 1];
    }

    memcpy(&y[j * 81], &ycol[0], 81U * sizeof(double));
  }

  for (k = 0; k < 81; k++) {
    ndx[k] = idx[k];
  }
}
#endif

/* End of code generation (sortrows.c) */
