/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * sortIdx.c
 *
 * Code generation for function 'sortIdx'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "sortIdx.h"

/* Function Definitions */

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void b_sortIdx(const double x[3], int idx[3])
{
  boolean_T p;
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
      } else {
        idx[0] = 3;
        idx[1] = 1;
        idx[2] = 2;
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
      } else {
        idx[0] = 3;
        idx[1] = 2;
        idx[2] = 1;
      }
    }
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void c_sortIdx(const double x[5], int idx[5])
{
  int k;
  int i;
  boolean_T p;
  int i2;
  int j;
  int pEnd;
  int b_p;
  int q;
  int qEnd;
  int kEnd;
  int iwork[5];
  for (k = 0; k < 5; k++) {
    idx[k] = 0;
  }

  for (k = 0; k <= 2; k += 2) {
    if ((x[k] <= x[k + 1]) || rtIsNaN(x[k + 1])) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      idx[k] = k + 1;
      idx[k + 1] = k + 2;
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
  }

  idx[4] = 5;
  i = 2;
  while (i < 5) {
    i2 = i << 1;
    j = 1;
    for (pEnd = 1 + i; pEnd < 6; pEnd = qEnd + i) {
      b_p = j;
      q = pEnd - 1;
      qEnd = j + i2;
      if (qEnd > 6) {
        qEnd = 6;
      }

      k = 0;
      kEnd = qEnd - j;
      while (k + 1 <= kEnd) {
        if ((x[idx[b_p - 1] - 1] <= x[idx[q] - 1]) || rtIsNaN(x[idx[q] - 1])) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
          iwork[k] = idx[b_p - 1];
          b_p++;
          if (b_p == pEnd) {
            while (q + 1 < qEnd) {
              k++;
              iwork[k] = idx[q];
              q++;
            }
          }
        } else {
          iwork[k] = idx[q];
          q++;
          if (q + 1 == qEnd) {
            while (b_p < pEnd) {
              k++;
              iwork[k] = idx[b_p - 1];
              b_p++;
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
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void d_sortIdx(const double x[7], int idx[7])
{
  int k;
  int i;
  boolean_T p;
  int i2;
  int j;
  int pEnd;
  int b_p;
  int q;
  int qEnd;
  int kEnd;
  int iwork[7];
  for (k = 0; k < 7; k++) {
    idx[k] = 0;
  }

  for (k = 0; k <= 4; k += 2) {
    if ((x[k] <= x[k + 1]) || rtIsNaN(x[k + 1])) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      idx[k] = k + 1;
      idx[k + 1] = k + 2;
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
  }

  idx[6] = 7;
  i = 2;
  while (i < 7) {
    i2 = i << 1;
    j = 1;
    for (pEnd = 1 + i; pEnd < 8; pEnd = qEnd + i) {
      b_p = j;
      q = pEnd - 1;
      qEnd = j + i2;
      if (qEnd > 8) {
        qEnd = 8;
      }

      k = 0;
      kEnd = qEnd - j;
      while (k + 1 <= kEnd) {
        if ((x[idx[b_p - 1] - 1] <= x[idx[q] - 1]) || rtIsNaN(x[idx[q] - 1])) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
          iwork[k] = idx[b_p - 1];
          b_p++;
          if (b_p == pEnd) {
            while (q + 1 < qEnd) {
              k++;
              iwork[k] = idx[q];
              q++;
            }
          }
        } else {
          iwork[k] = idx[q];
          q++;
          if (q + 1 == qEnd) {
            while (b_p < pEnd) {
              k++;
              iwork[k] = idx[b_p - 1];
              b_p++;
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
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void sortIdx(const double x[4], int idx[4])
{
  int k;
  int i;
  boolean_T p;
  int i2;
  int j;
  int pEnd;
  int b_p;
  int q;
  int qEnd;
  int kEnd;
  int iwork[4];
  for (k = 0; k <= 3; k += 2) {
    if ((x[k] <= x[k + 1]) || rtIsNaN(x[k + 1])) {
      p = true;
    } else {
      p = false;
    }

    if (p) {
      idx[k] = k + 1;
      idx[k + 1] = k + 2;
    } else {
      idx[k] = k + 2;
      idx[k + 1] = k + 1;
    }
  }

  i = 2;
  while (i < 4) {
    i2 = i << 1;
    j = 1;
    for (pEnd = 1 + i; pEnd < 5; pEnd = qEnd + i) {
      b_p = j;
      q = pEnd - 1;
      qEnd = j + i2;
      k = 0;
      kEnd = qEnd - j;
      while (k + 1 <= kEnd) {
        if ((x[idx[b_p - 1] - 1] <= x[idx[q] - 1]) || rtIsNaN(x[idx[q] - 1])) {
          p = true;
        } else {
          p = false;
        }

        if (p) {
          iwork[k] = idx[b_p - 1];
          b_p++;
          if (b_p == pEnd) {
            while (q + 1 < qEnd) {
              k++;
              iwork[k] = idx[q];
              q++;
            }
          }
        } else {
          iwork[k] = idx[q];
          q++;
          if (q + 1 == qEnd) {
            while (b_p < pEnd) {
              k++;
              iwork[k] = idx[b_p - 1];
              b_p++;
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
}
#endif

/* End of code generation (sortIdx.c) */
