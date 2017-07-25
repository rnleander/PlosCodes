/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * fminsearch.c
 *
 * Code generation for function 'fminsearch'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "fminsearch.h"
#include "emgpdf.h"
#include "sortIdx.h"
#include "eps.h"
#include "onestagepdf2.h"
#include "onestagepdf_lag.h"
#include "convolv_2invG_adapt_nov.h"
#include "convolv_3invG_nov.h"

/* Function Declarations */
static void b_sort_idx_single_insert(const double fv[7], int idx[7]);
static void sort_idx_single_insert(const double fv[5], int idx[5]);

/* Function Definitions */

/*
 *
 */
static void b_sort_idx_single_insert(const double fv[7], int idx[7])
{
  int k;
  int cidx;
  for (k = 0; k < 6; k++) {
    if (fv[idx[6 - k] - 1] < fv[idx[5 - k] - 1]) {
      cidx = idx[6 - k];
      idx[6 - k] = idx[5 - k];
      idx[5 - k] = cidx;
    }
  }
}

/*
 *
 */
static void sort_idx_single_insert(const double fv[5], int idx[5])
{
  int k;
  int cidx;
  for (k = 0; k < 4; k++) {
    if (fv[idx[4 - k] - 1] < fv[idx[3 - k] - 1]) {
      cidx = idx[4 - k];
      idx[4 - k] = idx[3 - k];
      idx[3 - k] = cidx;
    }
  }
}

/*
 *
 */
double b_fminsearch(double x[2])
{
  int k;
  static const double data[266] = { 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
    12.4, 13.7, 12.4, 11.9, 10.3, 12.9, 14.7, 11.6, 13.4, 13.4, 11.6, 10.3, 9.3,
    13.7, 9.6, 10.1, 9.8, 10.9, 16.0, 9.3, 9.6, 10.3, 11.4, 10.6, 8.5, 10.3,
    11.1, 8.0, 10.6, 7.5, 12.9, 9.0, 8.5, 12.4, 11.6, 9.6, 9.6, 14.7, 9.8, 10.3,
    12.1, 8.8, 10.6, 12.1, 13.4, 12.4, 8.8, 13.2, 10.1, 11.6, 11.1, 15.8, 12.1,
    12.7, 12.7, 11.1, 13.2, 11.9, 12.4, 13.2, 14.0, 8.0, 8.8, 9.3, 16.5, 14.5,
    10.1, 14.2, 7.8, 13.2, 8.8, 8.8, 10.1, 11.9, 12.9, 14.5, 10.9, 10.6, 14.0,
    8.8, 8.8, 9.0, 10.9, 14.5, 9.6, 12.4, 11.9, 12.4, 11.1, 14.5, 10.3, 12.4,
    12.7, 11.9, 10.3, 13.7, 15.5, 14.5, 11.6, 10.6, 15.5, 14.7, 8.8, 11.6, 8.3,
    17.6, 12.4, 11.6, 15.0, 13.7, 12.7, 10.9, 7.2, 8.5, 8.3, 9.6, 11.4, 12.9,
    11.6, 13.4, 10.1, 11.6, 8.8, 12.4, 10.3, 16.3, 10.9, 10.1, 8.8, 9.3, 15.2,
    8.5, 11.1, 8.3, 11.4, 11.9, 9.3, 9.8, 16.3, 12.7, 9.0, 11.9, 9.3, 10.3, 13.4,
    11.4, 12.9, 12.4, 9.6, 10.3, 13.2, 10.6, 9.8, 11.9, 14.2, 13.4, 9.3, 9.6,
    12.1, 11.9, 10.1, 14.0, 12.9, 21.7, 11.6, 12.1, 10.3, 9.8, 14.2, 13.7, 7.2,
    10.9, 10.1, 9.6, 13.4, 13.2, 16.3, 11.6, 14.0, 10.9, 14.2, 12.4, 12.4, 13.4,
    17.6, 10.1, 10.9, 14.0, 12.9, 9.0, 13.4, 15.0, 16.0, 8.0, 9.8, 12.4, 8.5,
    9.6, 12.7, 12.1, 15.0, 16.0, 10.9, 14.2, 13.7, 11.9, 16.8, 11.4, 13.4, 12.4,
    22.0, 12.4, 16.8, 12.1, 10.3, 13.4, 11.6, 10.1, 14.5, 10.6, 11.9, 15.5, 9.8,
    12.4, 10.1, 8.0, 9.0, 9.3, 13.2, 11.1, 12.7, 12.1, 10.1, 13.2, 14.5, 10.1,
    12.7, 12.9, 11.9, 12.4, 11.1, 8.5, 14.5, 16.5, 12.4, 9.0, 11.1, 9.8, 11.1,
    11.1, 8.8, 13.2, 17.6, 16.8, 10.9, 12.4, 8.5, 14.7 };

  double b_x[266];
  int j;
  double v[6];
  double b_v;
  double fv[3];
  int idx[3];
  int colIdx;
  boolean_T doShrink;
  int itercount;
  int fun_evals;
  int lastCol;
  int firstCol;
  boolean_T exitg1;
  double xr[2];
  double maxfv;
  double d4;
  boolean_T p;
  int b_firstCol;
  double xbar[2];
  double d5;
  double xe[2];
  double fvt[3];
  int idxb[3];
  for (k = 0; k < 3; k++) {
    for (j = 0; j < 2; j++) {
      v[j + (k << 1)] = x[j];
    }
  }

  onestagepdf2(data, x[0], x[1], b_x);
  for (k = 0; k < 266; k++) {
    b_x[k] = -log(b_x[k]);
  }

  b_v = b_x[0];
  for (k = 0; k < 265; k++) {
    b_v += b_x[k + 1];
  }

  fv[0] = b_v;
  for (k = 0; k < 2; k++) {
    colIdx = (k + 1) << 1;
    if (x[k] != 0.0) {
      v[k + colIdx] = 1.05 * x[k];
    } else {
      v[k + colIdx] = 0.00025;
    }

    for (j = 0; j < 2; j++) {
      xr[j] = v[j + colIdx];
    }

    onestagepdf2(data, xr[0], xr[1], b_x);
    for (lastCol = 0; lastCol < 266; lastCol++) {
      b_x[lastCol] = -log(b_x[lastCol]);
    }

    b_v = b_x[0];
    for (lastCol = 0; lastCol < 265; lastCol++) {
      b_v += b_x[lastCol + 1];
    }

    fv[k + 1] = b_v;
  }

  b_sortIdx(fv, idx);
  doShrink = false;
  itercount = 1;
  fun_evals = 3;
  lastCol = (idx[2] - 1) << 1;
  firstCol = (idx[0] - 1) << 1;
  exitg1 = false;
  while ((!exitg1) && ((fun_evals < 10000) && (itercount < 10000))) {
    maxfv = 0.0;
    for (k = 0; k < 2; k++) {
      b_v = fabs(fv[idx[0] - 1] - fv[idx[k + 1] - 1]);
      if (b_v > maxfv) {
        maxfv = b_v;
      }
    }

    b_v = 10.0 * eps(10.0 * eps(fv[idx[0] - 1]));
    if ((0.0001 > b_v) || rtIsNaN(b_v)) {
      d4 = 0.0001;
    } else {
      d4 = b_v;
    }

    if (maxfv > d4) {
      p = false;
    } else {
      maxfv = 0.0;
      b_firstCol = (idx[0] - 1) << 1;
      for (j = 0; j < 2; j++) {
        colIdx = (idx[j + 1] - 1) << 1;
        for (k = 0; k < 2; k++) {
          b_v = fabs(v[k + b_firstCol] - v[k + colIdx]);
          if (b_v > maxfv) {
            maxfv = b_v;
          }
        }
      }

      b_v = v[b_firstCol];
      if (v[b_firstCol + 1] > v[b_firstCol]) {
        b_v = v[b_firstCol + 1];
      }

      b_v = 10.0 * eps(b_v);
      if ((0.0001 > b_v) || rtIsNaN(b_v)) {
        d5 = 0.0001;
      } else {
        d5 = b_v;
      }

      p = (maxfv <= d5);
    }

    if (!p) {
      colIdx = (idx[1] - 1) << 1;
      for (k = 0; k < 2; k++) {
        b_v = (v[k + firstCol] + v[k + colIdx]) / 2.0;
        xr[k] = 2.0 * b_v - v[k + lastCol];
        xbar[k] = b_v;
      }

      onestagepdf2(data, xr[0], xr[1], b_x);
      for (k = 0; k < 266; k++) {
        b_x[k] = -log(b_x[k]);
      }

      b_v = b_x[0];
      for (k = 0; k < 265; k++) {
        b_v += b_x[k + 1];
      }

      fun_evals++;
      if (b_v < fv[idx[0] - 1]) {
        for (k = 0; k < 2; k++) {
          xe[k] = 3.0 * xbar[k] - 2.0 * v[k + lastCol];
        }

        onestagepdf2(data, xe[0], xe[1], b_x);
        for (k = 0; k < 266; k++) {
          b_x[k] = -log(b_x[k]);
        }

        maxfv = b_x[0];
        for (k = 0; k < 265; k++) {
          maxfv += b_x[k + 1];
        }

        fun_evals++;
        if (maxfv < b_v) {
          for (k = 0; k < 2; k++) {
            v[k + lastCol] = xe[k];
          }

          fv[idx[2] - 1] = maxfv;
        } else {
          for (k = 0; k < 2; k++) {
            v[k + lastCol] = xr[k];
          }

          fv[idx[2] - 1] = b_v;
        }
      } else if (b_v < fv[idx[1] - 1]) {
        for (k = 0; k < 2; k++) {
          v[k + lastCol] = xr[k];
        }

        fv[idx[2] - 1] = b_v;
      } else {
        if (b_v < fv[idx[2] - 1]) {
          for (k = 0; k < 2; k++) {
            x[k] = 1.5 * xbar[k] - 0.5 * v[k + lastCol];
          }

          onestagepdf2(data, x[0], x[1], b_x);
          for (k = 0; k < 266; k++) {
            b_x[k] = -log(b_x[k]);
          }

          maxfv = b_x[0];
          for (k = 0; k < 265; k++) {
            maxfv += b_x[k + 1];
          }

          fun_evals++;
          if (maxfv <= b_v) {
            for (k = 0; k < 2; k++) {
              v[k + lastCol] = x[k];
            }

            fv[idx[2] - 1] = maxfv;
          } else {
            doShrink = true;
          }
        } else {
          for (k = 0; k < 2; k++) {
            x[k] = 0.5 * xbar[k] + 0.5 * v[k + lastCol];
          }

          onestagepdf2(data, x[0], x[1], b_x);
          for (k = 0; k < 266; k++) {
            b_x[k] = -log(b_x[k]);
          }

          b_v = b_x[0];
          for (k = 0; k < 265; k++) {
            b_v += b_x[k + 1];
          }

          fun_evals++;
          if (b_v < fv[idx[2] - 1]) {
            for (k = 0; k < 2; k++) {
              v[k + lastCol] = x[k];
            }

            fv[idx[2] - 1] = b_v;
          } else {
            doShrink = true;
          }
        }

        if (doShrink) {
          for (k = 0; k < 2; k++) {
            colIdx = (idx[k + 1] - 1) << 1;
            for (j = 0; j < 2; j++) {
              v[j + colIdx] = v[j + firstCol] + 0.5 * (v[j + colIdx] - v[j +
                firstCol]);
              x[j] = v[j + colIdx];
            }

            onestagepdf2(data, x[0], x[1], b_x);
            for (lastCol = 0; lastCol < 266; lastCol++) {
              b_x[lastCol] = -log(b_x[lastCol]);
            }

            b_v = b_x[0];
            for (lastCol = 0; lastCol < 265; lastCol++) {
              b_v += b_x[lastCol + 1];
            }

            fv[idx[k + 1] - 1] = b_v;
          }

          fun_evals += 2;
        }
      }

      if (doShrink) {
        for (k = 0; k < 3; k++) {
          fvt[k] = fv[idx[k] - 1];
          idxb[k] = idx[k];
        }

        b_sortIdx(fvt, idx);
        for (k = 0; k < 3; k++) {
          idx[k] = idxb[idx[k] - 1];
        }

        doShrink = false;
      } else {
        for (k = 0; k < 2; k++) {
          if (fv[idx[2 - k] - 1] < fv[idx[1 - k] - 1]) {
            lastCol = idx[2 - k];
            idx[2 - k] = idx[1 - k];
            idx[1 - k] = lastCol;
          }
        }
      }

      itercount++;
      lastCol = (idx[2] - 1) << 1;
      firstCol = (idx[0] - 1) << 1;
    } else {
      exitg1 = true;
    }
  }

  colIdx = (idx[0] - 1) << 1;
  for (k = 0; k < 2; k++) {
    x[k] = v[k + colIdx];
  }

  return fv[idx[0] - 1];
}

/*
 *
 */
double c_fminsearch(double x[3])
{
  int k;
  static const double data[266] = { 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
    12.4, 13.7, 12.4, 11.9, 10.3, 12.9, 14.7, 11.6, 13.4, 13.4, 11.6, 10.3, 9.3,
    13.7, 9.6, 10.1, 9.8, 10.9, 16.0, 9.3, 9.6, 10.3, 11.4, 10.6, 8.5, 10.3,
    11.1, 8.0, 10.6, 7.5, 12.9, 9.0, 8.5, 12.4, 11.6, 9.6, 9.6, 14.7, 9.8, 10.3,
    12.1, 8.8, 10.6, 12.1, 13.4, 12.4, 8.8, 13.2, 10.1, 11.6, 11.1, 15.8, 12.1,
    12.7, 12.7, 11.1, 13.2, 11.9, 12.4, 13.2, 14.0, 8.0, 8.8, 9.3, 16.5, 14.5,
    10.1, 14.2, 7.8, 13.2, 8.8, 8.8, 10.1, 11.9, 12.9, 14.5, 10.9, 10.6, 14.0,
    8.8, 8.8, 9.0, 10.9, 14.5, 9.6, 12.4, 11.9, 12.4, 11.1, 14.5, 10.3, 12.4,
    12.7, 11.9, 10.3, 13.7, 15.5, 14.5, 11.6, 10.6, 15.5, 14.7, 8.8, 11.6, 8.3,
    17.6, 12.4, 11.6, 15.0, 13.7, 12.7, 10.9, 7.2, 8.5, 8.3, 9.6, 11.4, 12.9,
    11.6, 13.4, 10.1, 11.6, 8.8, 12.4, 10.3, 16.3, 10.9, 10.1, 8.8, 9.3, 15.2,
    8.5, 11.1, 8.3, 11.4, 11.9, 9.3, 9.8, 16.3, 12.7, 9.0, 11.9, 9.3, 10.3, 13.4,
    11.4, 12.9, 12.4, 9.6, 10.3, 13.2, 10.6, 9.8, 11.9, 14.2, 13.4, 9.3, 9.6,
    12.1, 11.9, 10.1, 14.0, 12.9, 21.7, 11.6, 12.1, 10.3, 9.8, 14.2, 13.7, 7.2,
    10.9, 10.1, 9.6, 13.4, 13.2, 16.3, 11.6, 14.0, 10.9, 14.2, 12.4, 12.4, 13.4,
    17.6, 10.1, 10.9, 14.0, 12.9, 9.0, 13.4, 15.0, 16.0, 8.0, 9.8, 12.4, 8.5,
    9.6, 12.7, 12.1, 15.0, 16.0, 10.9, 14.2, 13.7, 11.9, 16.8, 11.4, 13.4, 12.4,
    22.0, 12.4, 16.8, 12.1, 10.3, 13.4, 11.6, 10.1, 14.5, 10.6, 11.9, 15.5, 9.8,
    12.4, 10.1, 8.0, 9.0, 9.3, 13.2, 11.1, 12.7, 12.1, 10.1, 13.2, 14.5, 10.1,
    12.7, 12.9, 11.9, 12.4, 11.1, 8.5, 14.5, 16.5, 12.4, 9.0, 11.1, 9.8, 11.1,
    11.1, 8.8, 13.2, 17.6, 16.8, 10.9, 12.4, 8.5, 14.7 };

  double b_x[266];
  int j;
  double v[12];
  double b_v;
  double fv[4];
  int idx[4];
  int colIdx;
  boolean_T doShrink;
  int itercount;
  int fun_evals;
  int lastCol;
  int firstCol;
  boolean_T exitg1;
  double xr[3];
  double maxfv;
  double d6;
  boolean_T p;
  int b_firstCol;
  double xbar[3];
  double d7;
  double xe[3];
  double fvt[4];
  int idxb[4];
  for (k = 0; k < 4; k++) {
    for (j = 0; j < 3; j++) {
      v[j + 3 * k] = x[j];
    }
  }

  onestagepdf_lag(data, x[0], x[1], x[2], b_x);
  for (k = 0; k < 266; k++) {
    b_x[k] = -log(b_x[k]);
  }

  b_v = b_x[0];
  for (k = 0; k < 265; k++) {
    b_v += b_x[k + 1];
  }

  fv[0] = b_v;
  for (k = 0; k < 3; k++) {
    colIdx = 3 * (k + 1);
    if (x[k] != 0.0) {
      v[k + colIdx] = 1.05 * x[k];
    } else {
      v[k + colIdx] = 0.00025;
    }

    for (j = 0; j < 3; j++) {
      xr[j] = v[j + colIdx];
    }

    onestagepdf_lag(data, xr[0], xr[1], xr[2], b_x);
    for (lastCol = 0; lastCol < 266; lastCol++) {
      b_x[lastCol] = -log(b_x[lastCol]);
    }

    b_v = b_x[0];
    for (lastCol = 0; lastCol < 265; lastCol++) {
      b_v += b_x[lastCol + 1];
    }

    fv[k + 1] = b_v;
  }

  sortIdx(fv, idx);
  doShrink = false;
  itercount = 1;
  fun_evals = 4;
  lastCol = 3 * (idx[3] - 1);
  firstCol = 3 * (idx[0] - 1);
  exitg1 = false;
  while ((!exitg1) && ((fun_evals < 10000) && (itercount < 10000))) {
    maxfv = 0.0;
    for (k = 0; k < 3; k++) {
      b_v = fabs(fv[idx[0] - 1] - fv[idx[k + 1] - 1]);
      if (b_v > maxfv) {
        maxfv = b_v;
      }
    }

    b_v = 10.0 * eps(10.0 * eps(fv[idx[0] - 1]));
    if ((0.0001 > b_v) || rtIsNaN(b_v)) {
      d6 = 0.0001;
    } else {
      d6 = b_v;
    }

    if (maxfv > d6) {
      p = false;
    } else {
      maxfv = 0.0;
      b_firstCol = (idx[0] - 1) * 3;
      for (j = 0; j < 3; j++) {
        colIdx = (idx[j + 1] - 1) * 3;
        for (k = 0; k < 3; k++) {
          b_v = fabs(v[k + b_firstCol] - v[k + colIdx]);
          if (b_v > maxfv) {
            maxfv = b_v;
          }
        }
      }

      b_v = v[b_firstCol];
      for (k = 0; k < 2; k++) {
        if (v[(k + b_firstCol) + 1] > b_v) {
          b_v = v[(k + b_firstCol) + 1];
        }
      }

      b_v = 10.0 * eps(b_v);
      if ((0.0001 > b_v) || rtIsNaN(b_v)) {
        d7 = 0.0001;
      } else {
        d7 = b_v;
      }

      p = (maxfv <= d7);
    }

    if (!p) {
      for (k = 0; k < 3; k++) {
        xbar[k] = v[k + firstCol];
      }

      for (k = 0; k < 2; k++) {
        colIdx = 3 * (idx[k + 1] - 1);
        for (j = 0; j < 3; j++) {
          xbar[j] += v[j + colIdx];
        }
      }

      for (k = 0; k < 3; k++) {
        b_v = xbar[k] / 3.0;
        xr[k] = 2.0 * b_v - v[k + lastCol];
        xbar[k] = b_v;
      }

      onestagepdf_lag(data, xr[0], xr[1], xr[2], b_x);
      for (k = 0; k < 266; k++) {
        b_x[k] = -log(b_x[k]);
      }

      b_v = b_x[0];
      for (k = 0; k < 265; k++) {
        b_v += b_x[k + 1];
      }

      fun_evals++;
      if (b_v < fv[idx[0] - 1]) {
        for (k = 0; k < 3; k++) {
          xe[k] = 3.0 * xbar[k] - 2.0 * v[k + lastCol];
        }

        onestagepdf_lag(data, xe[0], xe[1], xe[2], b_x);
        for (k = 0; k < 266; k++) {
          b_x[k] = -log(b_x[k]);
        }

        maxfv = b_x[0];
        for (k = 0; k < 265; k++) {
          maxfv += b_x[k + 1];
        }

        fun_evals++;
        if (maxfv < b_v) {
          for (k = 0; k < 3; k++) {
            v[k + lastCol] = xe[k];
          }

          fv[idx[3] - 1] = maxfv;
        } else {
          for (k = 0; k < 3; k++) {
            v[k + lastCol] = xr[k];
          }

          fv[idx[3] - 1] = b_v;
        }
      } else if (b_v < fv[idx[2] - 1]) {
        for (k = 0; k < 3; k++) {
          v[k + lastCol] = xr[k];
        }

        fv[idx[3] - 1] = b_v;
      } else {
        if (b_v < fv[idx[3] - 1]) {
          for (k = 0; k < 3; k++) {
            x[k] = 1.5 * xbar[k] - 0.5 * v[k + lastCol];
          }

          onestagepdf_lag(data, x[0], x[1], x[2], b_x);
          for (k = 0; k < 266; k++) {
            b_x[k] = -log(b_x[k]);
          }

          maxfv = b_x[0];
          for (k = 0; k < 265; k++) {
            maxfv += b_x[k + 1];
          }

          fun_evals++;
          if (maxfv <= b_v) {
            for (k = 0; k < 3; k++) {
              v[k + lastCol] = x[k];
            }

            fv[idx[3] - 1] = maxfv;
          } else {
            doShrink = true;
          }
        } else {
          for (k = 0; k < 3; k++) {
            x[k] = 0.5 * xbar[k] + 0.5 * v[k + lastCol];
          }

          onestagepdf_lag(data, x[0], x[1], x[2], b_x);
          for (k = 0; k < 266; k++) {
            b_x[k] = -log(b_x[k]);
          }

          b_v = b_x[0];
          for (k = 0; k < 265; k++) {
            b_v += b_x[k + 1];
          }

          fun_evals++;
          if (b_v < fv[idx[3] - 1]) {
            for (k = 0; k < 3; k++) {
              v[k + lastCol] = x[k];
            }

            fv[idx[3] - 1] = b_v;
          } else {
            doShrink = true;
          }
        }

        if (doShrink) {
          for (k = 0; k < 3; k++) {
            colIdx = 3 * (idx[k + 1] - 1);
            for (j = 0; j < 3; j++) {
              v[j + colIdx] = v[j + firstCol] + 0.5 * (v[j + colIdx] - v[j +
                firstCol]);
              x[j] = v[j + colIdx];
            }

            onestagepdf_lag(data, x[0], x[1], x[2], b_x);
            for (lastCol = 0; lastCol < 266; lastCol++) {
              b_x[lastCol] = -log(b_x[lastCol]);
            }

            b_v = b_x[0];
            for (lastCol = 0; lastCol < 265; lastCol++) {
              b_v += b_x[lastCol + 1];
            }

            fv[idx[k + 1] - 1] = b_v;
          }

          fun_evals += 3;
        }
      }

      if (doShrink) {
        for (k = 0; k < 4; k++) {
          fvt[k] = fv[idx[k] - 1];
          idxb[k] = idx[k];
        }

        sortIdx(fvt, idx);
        for (k = 0; k < 4; k++) {
          idx[k] = idxb[idx[k] - 1];
        }

        doShrink = false;
      } else {
        for (k = 0; k < 3; k++) {
          if (fv[idx[3 - k] - 1] < fv[idx[2 - k] - 1]) {
            lastCol = idx[3 - k];
            idx[3 - k] = idx[2 - k];
            idx[2 - k] = lastCol;
          }
        }
      }

      itercount++;
      lastCol = 3 * (idx[3] - 1);
      firstCol = 3 * (idx[0] - 1);
    } else {
      exitg1 = true;
    }
  }

  colIdx = 3 * (idx[0] - 1);
  for (k = 0; k < 3; k++) {
    x[k] = v[k + colIdx];
  }

  return fv[idx[0] - 1];
}

/*
 *
 */
double d_fminsearch(double x[4])
{
  double unusedExpr[266];
  double b_unusedExpr[266];
  double c_unusedExpr[266];
  double d_unusedExpr[266];
  int k;
  double b_x[266];
  int j;
  double v[20];
  double b_v;
  double fv[5];
  int idx[5];
  int colIdx;
  boolean_T doShrink;
  int itercount;
  int fun_evals;
  int lastCol;
  double xr[4];
  int firstCol;
  boolean_T exitg1;
  double maxfv;
  double d8;
  boolean_T p;
  int b_firstCol;
  double xbar[4];
  double d9;
  double xe[4];
  double fvt[5];
  int idxb[5];
  convolv_2invG_adapt_nov(x[0], x[1], x[2], x[3], unusedExpr);
  convolv_2invG_adapt_nov(x[0], x[1], x[2], x[3], b_unusedExpr);
  convolv_2invG_adapt_nov(x[0], x[1], x[2], x[3], c_unusedExpr);
  convolv_2invG_adapt_nov(x[0], x[1], x[2], x[3], d_unusedExpr);
  for (k = 0; k < 5; k++) {
    for (j = 0; j < 4; j++) {
      v[j + (k << 2)] = x[j];
    }
  }

  convolv_2invG_adapt_nov(x[0], x[1], x[2], x[3], b_x);
  for (k = 0; k < 266; k++) {
    b_x[k] = -log(b_x[k]);
  }

  b_v = b_x[0];
  for (k = 0; k < 265; k++) {
    b_v += b_x[k + 1];
  }

  fv[0] = b_v;
  for (k = 0; k < 4; k++) {
    colIdx = (k + 1) << 2;
    v[k + colIdx] = 1.05 * x[k];
    for (j = 0; j < 4; j++) {
      xr[j] = v[j + colIdx];
    }

    convolv_2invG_adapt_nov(xr[0], xr[1], xr[2], xr[3], b_x);
    for (lastCol = 0; lastCol < 266; lastCol++) {
      b_x[lastCol] = -log(b_x[lastCol]);
    }

    b_v = b_x[0];
    for (lastCol = 0; lastCol < 265; lastCol++) {
      b_v += b_x[lastCol + 1];
    }

    fv[k + 1] = b_v;
  }

  c_sortIdx(fv, idx);
  doShrink = false;
  itercount = 1;
  fun_evals = 5;
  lastCol = (idx[4] - 1) << 2;
  firstCol = (idx[0] - 1) << 2;
  exitg1 = false;
  while ((!exitg1) && ((fun_evals < 10000) && (itercount < 10000))) {
    maxfv = 0.0;
    for (k = 0; k < 4; k++) {
      b_v = fabs(fv[idx[0] - 1] - fv[idx[k + 1] - 1]);
      if (b_v > maxfv) {
        maxfv = b_v;
      }
    }

    b_v = 10.0 * eps(10.0 * eps(fv[idx[0] - 1]));
    if ((0.001 > b_v) || rtIsNaN(b_v)) {
      d8 = 0.001;
    } else {
      d8 = b_v;
    }

    if (maxfv > d8) {
      p = false;
    } else {
      maxfv = 0.0;
      b_firstCol = (idx[0] - 1) << 2;
      for (j = 0; j < 4; j++) {
        colIdx = (idx[j + 1] - 1) << 2;
        for (k = 0; k < 4; k++) {
          b_v = fabs(v[k + b_firstCol] - v[k + colIdx]);
          if (b_v > maxfv) {
            maxfv = b_v;
          }
        }
      }

      b_v = v[b_firstCol];
      for (k = 0; k < 3; k++) {
        if (v[(k + b_firstCol) + 1] > b_v) {
          b_v = v[(k + b_firstCol) + 1];
        }
      }

      b_v = 10.0 * eps(b_v);
      if ((0.001 > b_v) || rtIsNaN(b_v)) {
        d9 = 0.001;
      } else {
        d9 = b_v;
      }

      p = (maxfv <= d9);
    }

    if (!p) {
      for (k = 0; k < 4; k++) {
        xbar[k] = v[k + firstCol];
      }

      for (k = 0; k < 3; k++) {
        colIdx = (idx[k + 1] - 1) << 2;
        for (j = 0; j < 4; j++) {
          xbar[j] += v[j + colIdx];
        }
      }

      for (k = 0; k < 4; k++) {
        b_v = xbar[k] / 4.0;
        xr[k] = 2.0 * b_v - v[k + lastCol];
        xbar[k] = b_v;
      }

      convolv_2invG_adapt_nov(xr[0], xr[1], xr[2], xr[3], b_x);
      for (k = 0; k < 266; k++) {
        b_x[k] = -log(b_x[k]);
      }

      b_v = b_x[0];
      for (k = 0; k < 265; k++) {
        b_v += b_x[k + 1];
      }

      fun_evals++;
      if (b_v < fv[idx[0] - 1]) {
        for (k = 0; k < 4; k++) {
          xe[k] = 3.0 * xbar[k] - 2.0 * v[k + lastCol];
        }

        convolv_2invG_adapt_nov(xe[0], xe[1], xe[2], xe[3], b_x);
        for (k = 0; k < 266; k++) {
          b_x[k] = -log(b_x[k]);
        }

        maxfv = b_x[0];
        for (k = 0; k < 265; k++) {
          maxfv += b_x[k + 1];
        }

        fun_evals++;
        if (maxfv < b_v) {
          for (k = 0; k < 4; k++) {
            v[k + lastCol] = xe[k];
          }

          fv[idx[4] - 1] = maxfv;
        } else {
          for (k = 0; k < 4; k++) {
            v[k + lastCol] = xr[k];
          }

          fv[idx[4] - 1] = b_v;
        }
      } else if (b_v < fv[idx[3] - 1]) {
        for (k = 0; k < 4; k++) {
          v[k + lastCol] = xr[k];
        }

        fv[idx[4] - 1] = b_v;
      } else {
        if (b_v < fv[idx[4] - 1]) {
          for (k = 0; k < 4; k++) {
            x[k] = 1.5 * xbar[k] - 0.5 * v[k + lastCol];
          }

          convolv_2invG_adapt_nov(x[0], x[1], x[2], x[3], b_x);
          for (k = 0; k < 266; k++) {
            b_x[k] = -log(b_x[k]);
          }

          maxfv = b_x[0];
          for (k = 0; k < 265; k++) {
            maxfv += b_x[k + 1];
          }

          fun_evals++;
          if (maxfv <= b_v) {
            for (k = 0; k < 4; k++) {
              v[k + lastCol] = x[k];
            }

            fv[idx[4] - 1] = maxfv;
          } else {
            doShrink = true;
          }
        } else {
          for (k = 0; k < 4; k++) {
            x[k] = 0.5 * xbar[k] + 0.5 * v[k + lastCol];
          }

          convolv_2invG_adapt_nov(x[0], x[1], x[2], x[3], b_x);
          for (k = 0; k < 266; k++) {
            b_x[k] = -log(b_x[k]);
          }

          b_v = b_x[0];
          for (k = 0; k < 265; k++) {
            b_v += b_x[k + 1];
          }

          fun_evals++;
          if (b_v < fv[idx[4] - 1]) {
            for (k = 0; k < 4; k++) {
              v[k + lastCol] = x[k];
            }

            fv[idx[4] - 1] = b_v;
          } else {
            doShrink = true;
          }
        }

        if (doShrink) {
          for (k = 0; k < 4; k++) {
            colIdx = (idx[k + 1] - 1) << 2;
            for (j = 0; j < 4; j++) {
              v[j + colIdx] = v[j + firstCol] + 0.5 * (v[j + colIdx] - v[j +
                firstCol]);
              x[j] = v[j + colIdx];
            }

            convolv_2invG_adapt_nov(x[0], x[1], x[2], x[3], b_x);
            for (lastCol = 0; lastCol < 266; lastCol++) {
              b_x[lastCol] = -log(b_x[lastCol]);
            }

            b_v = b_x[0];
            for (lastCol = 0; lastCol < 265; lastCol++) {
              b_v += b_x[lastCol + 1];
            }

            fv[idx[k + 1] - 1] = b_v;
          }

          fun_evals += 4;
        }
      }

      if (doShrink) {
        for (k = 0; k < 5; k++) {
          fvt[k] = fv[idx[k] - 1];
          idxb[k] = idx[k];
        }

        c_sortIdx(fvt, idx);
        for (k = 0; k < 5; k++) {
          idx[k] = idxb[idx[k] - 1];
        }

        doShrink = false;
      } else {
        sort_idx_single_insert(fv, idx);
      }

      itercount++;
      lastCol = (idx[4] - 1) << 2;
      firstCol = (idx[0] - 1) << 2;
    } else {
      exitg1 = true;
    }
  }

  colIdx = (idx[0] - 1) << 2;
  for (k = 0; k < 4; k++) {
    x[k] = v[k + colIdx];
  }

  return fv[idx[0] - 1];
}

/*
 *
 */
double e_fminsearch(double x[6])
{
  double unusedExpr[266];
  double b_unusedExpr[266];
  double c_unusedExpr[266];
  double d_unusedExpr[266];
  int k;
  double b_x[266];
  int j;
  double v[42];
  double b_v;
  double fv[7];
  int idx[7];
  int colIdx;
  boolean_T doShrink;
  int itercount;
  int fun_evals;
  int lastCol;
  double xr[6];
  int firstCol;
  boolean_T exitg1;
  double maxfv;
  double d10;
  boolean_T p;
  int b_firstCol;
  double xbar[6];
  double d11;
  double xe[6];
  double fvt[7];
  int idxb[7];
  convolv_3invG_nov(x[0], x[1], x[2], x[3], x[4], x[5], unusedExpr);
  convolv_3invG_nov(x[0], x[1], x[2], x[3], x[4], x[5], b_unusedExpr);
  convolv_3invG_nov(x[0], x[1], x[2], x[3], x[4], x[5], c_unusedExpr);
  convolv_3invG_nov(x[0], x[1], x[2], x[3], x[4], x[5], d_unusedExpr);
  for (k = 0; k < 7; k++) {
    for (j = 0; j < 6; j++) {
      v[j + 6 * k] = x[j];
    }
  }

  convolv_3invG_nov(x[0], x[1], x[2], x[3], x[4], x[5], b_x);
  for (k = 0; k < 266; k++) {
    b_x[k] = -log(b_x[k]);
  }

  b_v = b_x[0];
  for (k = 0; k < 265; k++) {
    b_v += b_x[k + 1];
  }

  fv[0] = b_v;
  for (k = 0; k < 6; k++) {
    colIdx = 6 * (k + 1);
    v[k + colIdx] = 1.05 * x[k];
    for (j = 0; j < 6; j++) {
      xr[j] = v[j + colIdx];
    }

    convolv_3invG_nov(xr[0], xr[1], xr[2], xr[3], xr[4], xr[5], b_x);
    for (lastCol = 0; lastCol < 266; lastCol++) {
      b_x[lastCol] = -log(b_x[lastCol]);
    }

    b_v = b_x[0];
    for (lastCol = 0; lastCol < 265; lastCol++) {
      b_v += b_x[lastCol + 1];
    }

    fv[k + 1] = b_v;
  }

  d_sortIdx(fv, idx);
  doShrink = false;
  itercount = 1;
  fun_evals = 7;
  lastCol = 6 * (idx[6] - 1);
  firstCol = 6 * (idx[0] - 1);
  exitg1 = false;
  while ((!exitg1) && ((fun_evals < 10000) && (itercount < 10000))) {
    maxfv = 0.0;
    for (k = 0; k < 6; k++) {
      b_v = fabs(fv[idx[0] - 1] - fv[idx[k + 1] - 1]);
      if (b_v > maxfv) {
        maxfv = b_v;
      }
    }

    b_v = 10.0 * eps(10.0 * eps(fv[idx[0] - 1]));
    if ((0.001 > b_v) || rtIsNaN(b_v)) {
      d10 = 0.001;
    } else {
      d10 = b_v;
    }

    if (maxfv > d10) {
      p = false;
    } else {
      maxfv = 0.0;
      b_firstCol = (idx[0] - 1) * 6;
      for (j = 0; j < 6; j++) {
        colIdx = (idx[j + 1] - 1) * 6;
        for (k = 0; k < 6; k++) {
          b_v = fabs(v[k + b_firstCol] - v[k + colIdx]);
          if (b_v > maxfv) {
            maxfv = b_v;
          }
        }
      }

      b_v = v[b_firstCol];
      for (k = 0; k < 5; k++) {
        if (v[(k + b_firstCol) + 1] > b_v) {
          b_v = v[(k + b_firstCol) + 1];
        }
      }

      b_v = 10.0 * eps(b_v);
      if ((0.001 > b_v) || rtIsNaN(b_v)) {
        d11 = 0.001;
      } else {
        d11 = b_v;
      }

      p = (maxfv <= d11);
    }

    if (!p) {
      for (k = 0; k < 6; k++) {
        xbar[k] = v[k + firstCol];
      }

      for (k = 0; k < 5; k++) {
        colIdx = 6 * (idx[k + 1] - 1);
        for (j = 0; j < 6; j++) {
          xbar[j] += v[j + colIdx];
        }
      }

      for (k = 0; k < 6; k++) {
        b_v = xbar[k] / 6.0;
        xr[k] = 2.0 * b_v - v[k + lastCol];
        xbar[k] = b_v;
      }

      convolv_3invG_nov(xr[0], xr[1], xr[2], xr[3], xr[4], xr[5], b_x);
      for (k = 0; k < 266; k++) {
        b_x[k] = -log(b_x[k]);
      }

      b_v = b_x[0];
      for (k = 0; k < 265; k++) {
        b_v += b_x[k + 1];
      }

      fun_evals++;
      if (b_v < fv[idx[0] - 1]) {
        for (k = 0; k < 6; k++) {
          xe[k] = 3.0 * xbar[k] - 2.0 * v[k + lastCol];
        }

        convolv_3invG_nov(xe[0], xe[1], xe[2], xe[3], xe[4], xe[5], b_x);
        for (k = 0; k < 266; k++) {
          b_x[k] = -log(b_x[k]);
        }

        maxfv = b_x[0];
        for (k = 0; k < 265; k++) {
          maxfv += b_x[k + 1];
        }

        fun_evals++;
        if (maxfv < b_v) {
          for (k = 0; k < 6; k++) {
            v[k + lastCol] = xe[k];
          }

          fv[idx[6] - 1] = maxfv;
        } else {
          for (k = 0; k < 6; k++) {
            v[k + lastCol] = xr[k];
          }

          fv[idx[6] - 1] = b_v;
        }
      } else if (b_v < fv[idx[5] - 1]) {
        for (k = 0; k < 6; k++) {
          v[k + lastCol] = xr[k];
        }

        fv[idx[6] - 1] = b_v;
      } else {
        if (b_v < fv[idx[6] - 1]) {
          for (k = 0; k < 6; k++) {
            x[k] = 1.5 * xbar[k] - 0.5 * v[k + lastCol];
          }

          convolv_3invG_nov(x[0], x[1], x[2], x[3], x[4], x[5], b_x);
          for (k = 0; k < 266; k++) {
            b_x[k] = -log(b_x[k]);
          }

          maxfv = b_x[0];
          for (k = 0; k < 265; k++) {
            maxfv += b_x[k + 1];
          }

          fun_evals++;
          if (maxfv <= b_v) {
            for (k = 0; k < 6; k++) {
              v[k + lastCol] = x[k];
            }

            fv[idx[6] - 1] = maxfv;
          } else {
            doShrink = true;
          }
        } else {
          for (k = 0; k < 6; k++) {
            x[k] = 0.5 * xbar[k] + 0.5 * v[k + lastCol];
          }

          convolv_3invG_nov(x[0], x[1], x[2], x[3], x[4], x[5], b_x);
          for (k = 0; k < 266; k++) {
            b_x[k] = -log(b_x[k]);
          }

          b_v = b_x[0];
          for (k = 0; k < 265; k++) {
            b_v += b_x[k + 1];
          }

          fun_evals++;
          if (b_v < fv[idx[6] - 1]) {
            for (k = 0; k < 6; k++) {
              v[k + lastCol] = x[k];
            }

            fv[idx[6] - 1] = b_v;
          } else {
            doShrink = true;
          }
        }

        if (doShrink) {
          for (k = 0; k < 6; k++) {
            colIdx = 6 * (idx[k + 1] - 1);
            for (j = 0; j < 6; j++) {
              v[j + colIdx] = v[j + firstCol] + 0.5 * (v[j + colIdx] - v[j +
                firstCol]);
              x[j] = v[j + colIdx];
            }

            convolv_3invG_nov(x[0], x[1], x[2], x[3], x[4], x[5], b_x);
            for (lastCol = 0; lastCol < 266; lastCol++) {
              b_x[lastCol] = -log(b_x[lastCol]);
            }

            b_v = b_x[0];
            for (lastCol = 0; lastCol < 265; lastCol++) {
              b_v += b_x[lastCol + 1];
            }

            fv[idx[k + 1] - 1] = b_v;
          }

          fun_evals += 6;
        }
      }

      if (doShrink) {
        for (k = 0; k < 7; k++) {
          fvt[k] = fv[idx[k] - 1];
          idxb[k] = idx[k];
        }

        d_sortIdx(fvt, idx);
        for (k = 0; k < 7; k++) {
          idx[k] = idxb[idx[k] - 1];
        }

        doShrink = false;
      } else {
        b_sort_idx_single_insert(fv, idx);
      }

      itercount++;
      lastCol = 6 * (idx[6] - 1);
      firstCol = 6 * (idx[0] - 1);
    } else {
      exitg1 = true;
    }
  }

  colIdx = 6 * (idx[0] - 1);
  for (k = 0; k < 6; k++) {
    x[k] = v[k + colIdx];
  }

  return fv[idx[0] - 1];
}

/*
 *
 */
double fminsearch(double x[3])
{
  int k;
  static const double data[266] = { 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
    12.4, 13.7, 12.4, 11.9, 10.3, 12.9, 14.7, 11.6, 13.4, 13.4, 11.6, 10.3, 9.3,
    13.7, 9.6, 10.1, 9.8, 10.9, 16.0, 9.3, 9.6, 10.3, 11.4, 10.6, 8.5, 10.3,
    11.1, 8.0, 10.6, 7.5, 12.9, 9.0, 8.5, 12.4, 11.6, 9.6, 9.6, 14.7, 9.8, 10.3,
    12.1, 8.8, 10.6, 12.1, 13.4, 12.4, 8.8, 13.2, 10.1, 11.6, 11.1, 15.8, 12.1,
    12.7, 12.7, 11.1, 13.2, 11.9, 12.4, 13.2, 14.0, 8.0, 8.8, 9.3, 16.5, 14.5,
    10.1, 14.2, 7.8, 13.2, 8.8, 8.8, 10.1, 11.9, 12.9, 14.5, 10.9, 10.6, 14.0,
    8.8, 8.8, 9.0, 10.9, 14.5, 9.6, 12.4, 11.9, 12.4, 11.1, 14.5, 10.3, 12.4,
    12.7, 11.9, 10.3, 13.7, 15.5, 14.5, 11.6, 10.6, 15.5, 14.7, 8.8, 11.6, 8.3,
    17.6, 12.4, 11.6, 15.0, 13.7, 12.7, 10.9, 7.2, 8.5, 8.3, 9.6, 11.4, 12.9,
    11.6, 13.4, 10.1, 11.6, 8.8, 12.4, 10.3, 16.3, 10.9, 10.1, 8.8, 9.3, 15.2,
    8.5, 11.1, 8.3, 11.4, 11.9, 9.3, 9.8, 16.3, 12.7, 9.0, 11.9, 9.3, 10.3, 13.4,
    11.4, 12.9, 12.4, 9.6, 10.3, 13.2, 10.6, 9.8, 11.9, 14.2, 13.4, 9.3, 9.6,
    12.1, 11.9, 10.1, 14.0, 12.9, 21.7, 11.6, 12.1, 10.3, 9.8, 14.2, 13.7, 7.2,
    10.9, 10.1, 9.6, 13.4, 13.2, 16.3, 11.6, 14.0, 10.9, 14.2, 12.4, 12.4, 13.4,
    17.6, 10.1, 10.9, 14.0, 12.9, 9.0, 13.4, 15.0, 16.0, 8.0, 9.8, 12.4, 8.5,
    9.6, 12.7, 12.1, 15.0, 16.0, 10.9, 14.2, 13.7, 11.9, 16.8, 11.4, 13.4, 12.4,
    22.0, 12.4, 16.8, 12.1, 10.3, 13.4, 11.6, 10.1, 14.5, 10.6, 11.9, 15.5, 9.8,
    12.4, 10.1, 8.0, 9.0, 9.3, 13.2, 11.1, 12.7, 12.1, 10.1, 13.2, 14.5, 10.1,
    12.7, 12.9, 11.9, 12.4, 11.1, 8.5, 14.5, 16.5, 12.4, 9.0, 11.1, 9.8, 11.1,
    11.1, 8.8, 13.2, 17.6, 16.8, 10.9, 12.4, 8.5, 14.7 };

  double b_x[266];
  int j;
  double v[12];
  double b_v;
  double fv[4];
  int idx[4];
  int colIdx;
  boolean_T doShrink;
  int itercount;
  int fun_evals;
  int lastCol;
  int firstCol;
  boolean_T exitg1;
  double xr[3];
  double maxfv;
  double d2;
  boolean_T p;
  int b_firstCol;
  double xbar[3];
  double d3;
  double xe[3];
  double fvt[4];
  int idxb[4];
  for (k = 0; k < 4; k++) {
    for (j = 0; j < 3; j++) {
      v[j + 3 * k] = x[j];
    }
  }

  emgpdf(data, x[0], x[1], x[2], b_x);
  for (k = 0; k < 266; k++) {
    b_x[k] = -log(b_x[k]);
  }

  b_v = b_x[0];
  for (k = 0; k < 265; k++) {
    b_v += b_x[k + 1];
  }

  fv[0] = b_v;
  for (k = 0; k < 3; k++) {
    colIdx = 3 * (k + 1);
    if (x[k] != 0.0) {
      v[k + colIdx] = 1.05 * x[k];
    } else {
      v[k + colIdx] = 0.00025;
    }

    for (j = 0; j < 3; j++) {
      xr[j] = v[j + colIdx];
    }

    emgpdf(data, xr[0], xr[1], xr[2], b_x);
    for (lastCol = 0; lastCol < 266; lastCol++) {
      b_x[lastCol] = -log(b_x[lastCol]);
    }

    b_v = b_x[0];
    for (lastCol = 0; lastCol < 265; lastCol++) {
      b_v += b_x[lastCol + 1];
    }

    fv[k + 1] = b_v;
  }

  sortIdx(fv, idx);
  doShrink = false;
  itercount = 1;
  fun_evals = 4;
  lastCol = 3 * (idx[3] - 1);
  firstCol = 3 * (idx[0] - 1);
  exitg1 = false;
  while ((!exitg1) && ((fun_evals < 10000) && (itercount < 10000))) {
    maxfv = 0.0;
    for (k = 0; k < 3; k++) {
      b_v = fabs(fv[idx[0] - 1] - fv[idx[k + 1] - 1]);
      if (b_v > maxfv) {
        maxfv = b_v;
      }
    }

    b_v = 10.0 * eps(10.0 * eps(fv[idx[0] - 1]));
    if ((0.0001 > b_v) || rtIsNaN(b_v)) {
      d2 = 0.0001;
    } else {
      d2 = b_v;
    }

    if (maxfv > d2) {
      p = false;
    } else {
      maxfv = 0.0;
      b_firstCol = (idx[0] - 1) * 3;
      for (j = 0; j < 3; j++) {
        colIdx = (idx[j + 1] - 1) * 3;
        for (k = 0; k < 3; k++) {
          b_v = fabs(v[k + b_firstCol] - v[k + colIdx]);
          if (b_v > maxfv) {
            maxfv = b_v;
          }
        }
      }

      b_v = v[b_firstCol];
      for (k = 0; k < 2; k++) {
        if (v[(k + b_firstCol) + 1] > b_v) {
          b_v = v[(k + b_firstCol) + 1];
        }
      }

      b_v = 10.0 * eps(b_v);
      if ((0.0001 > b_v) || rtIsNaN(b_v)) {
        d3 = 0.0001;
      } else {
        d3 = b_v;
      }

      p = (maxfv <= d3);
    }

    if (!p) {
      for (k = 0; k < 3; k++) {
        xbar[k] = v[k + firstCol];
      }

      for (k = 0; k < 2; k++) {
        colIdx = 3 * (idx[k + 1] - 1);
        for (j = 0; j < 3; j++) {
          xbar[j] += v[j + colIdx];
        }
      }

      for (k = 0; k < 3; k++) {
        b_v = xbar[k] / 3.0;
        xr[k] = 2.0 * b_v - v[k + lastCol];
        xbar[k] = b_v;
      }

      emgpdf(data, xr[0], xr[1], xr[2], b_x);
      for (k = 0; k < 266; k++) {
        b_x[k] = -log(b_x[k]);
      }

      b_v = b_x[0];
      for (k = 0; k < 265; k++) {
        b_v += b_x[k + 1];
      }

      fun_evals++;
      if (b_v < fv[idx[0] - 1]) {
        for (k = 0; k < 3; k++) {
          xe[k] = 3.0 * xbar[k] - 2.0 * v[k + lastCol];
        }

        emgpdf(data, xe[0], xe[1], xe[2], b_x);
        for (k = 0; k < 266; k++) {
          b_x[k] = -log(b_x[k]);
        }

        maxfv = b_x[0];
        for (k = 0; k < 265; k++) {
          maxfv += b_x[k + 1];
        }

        fun_evals++;
        if (maxfv < b_v) {
          for (k = 0; k < 3; k++) {
            v[k + lastCol] = xe[k];
          }

          fv[idx[3] - 1] = maxfv;
        } else {
          for (k = 0; k < 3; k++) {
            v[k + lastCol] = xr[k];
          }

          fv[idx[3] - 1] = b_v;
        }
      } else if (b_v < fv[idx[2] - 1]) {
        for (k = 0; k < 3; k++) {
          v[k + lastCol] = xr[k];
        }

        fv[idx[3] - 1] = b_v;
      } else {
        if (b_v < fv[idx[3] - 1]) {
          for (k = 0; k < 3; k++) {
            x[k] = 1.5 * xbar[k] - 0.5 * v[k + lastCol];
          }

          emgpdf(data, x[0], x[1], x[2], b_x);
          for (k = 0; k < 266; k++) {
            b_x[k] = -log(b_x[k]);
          }

          maxfv = b_x[0];
          for (k = 0; k < 265; k++) {
            maxfv += b_x[k + 1];
          }

          fun_evals++;
          if (maxfv <= b_v) {
            for (k = 0; k < 3; k++) {
              v[k + lastCol] = x[k];
            }

            fv[idx[3] - 1] = maxfv;
          } else {
            doShrink = true;
          }
        } else {
          for (k = 0; k < 3; k++) {
            x[k] = 0.5 * xbar[k] + 0.5 * v[k + lastCol];
          }

          emgpdf(data, x[0], x[1], x[2], b_x);
          for (k = 0; k < 266; k++) {
            b_x[k] = -log(b_x[k]);
          }

          b_v = b_x[0];
          for (k = 0; k < 265; k++) {
            b_v += b_x[k + 1];
          }

          fun_evals++;
          if (b_v < fv[idx[3] - 1]) {
            for (k = 0; k < 3; k++) {
              v[k + lastCol] = x[k];
            }

            fv[idx[3] - 1] = b_v;
          } else {
            doShrink = true;
          }
        }

        if (doShrink) {
          for (k = 0; k < 3; k++) {
            colIdx = 3 * (idx[k + 1] - 1);
            for (j = 0; j < 3; j++) {
              v[j + colIdx] = v[j + firstCol] + 0.5 * (v[j + colIdx] - v[j +
                firstCol]);
              x[j] = v[j + colIdx];
            }

            emgpdf(data, x[0], x[1], x[2], b_x);
            for (lastCol = 0; lastCol < 266; lastCol++) {
              b_x[lastCol] = -log(b_x[lastCol]);
            }

            b_v = b_x[0];
            for (lastCol = 0; lastCol < 265; lastCol++) {
              b_v += b_x[lastCol + 1];
            }

            fv[idx[k + 1] - 1] = b_v;
          }

          fun_evals += 3;
        }
      }

      if (doShrink) {
        for (k = 0; k < 4; k++) {
          fvt[k] = fv[idx[k] - 1];
          idxb[k] = idx[k];
        }

        sortIdx(fvt, idx);
        for (k = 0; k < 4; k++) {
          idx[k] = idxb[idx[k] - 1];
        }

        doShrink = false;
      } else {
        for (k = 0; k < 3; k++) {
          if (fv[idx[3 - k] - 1] < fv[idx[2 - k] - 1]) {
            lastCol = idx[3 - k];
            idx[3 - k] = idx[2 - k];
            idx[2 - k] = lastCol;
          }
        }
      }

      itercount++;
      lastCol = 3 * (idx[3] - 1);
      firstCol = 3 * (idx[0] - 1);
    } else {
      exitg1 = true;
    }
  }

  colIdx = 3 * (idx[0] - 1);
  for (k = 0; k < 3; k++) {
    x[k] = v[k + colIdx];
  }

  return fv[idx[0] - 1];
}



double fminsearch_generalized(double x[3], void (*pdf)(), double data[266])
{
	int k;
	/*
	static const double data[266] = { 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
		12.4, 13.7, 12.4, 11.9, 10.3, 12.9, 14.7, 11.6, 13.4, 13.4, 11.6, 10.3, 9.3,
		13.7, 9.6, 10.1, 9.8, 10.9, 16.0, 9.3, 9.6, 10.3, 11.4, 10.6, 8.5, 10.3,
		11.1, 8.0, 10.6, 7.5, 12.9, 9.0, 8.5, 12.4, 11.6, 9.6, 9.6, 14.7, 9.8, 10.3,
		12.1, 8.8, 10.6, 12.1, 13.4, 12.4, 8.8, 13.2, 10.1, 11.6, 11.1, 15.8, 12.1,
		12.7, 12.7, 11.1, 13.2, 11.9, 12.4, 13.2, 14.0, 8.0, 8.8, 9.3, 16.5, 14.5,
		10.1, 14.2, 7.8, 13.2, 8.8, 8.8, 10.1, 11.9, 12.9, 14.5, 10.9, 10.6, 14.0,
		8.8, 8.8, 9.0, 10.9, 14.5, 9.6, 12.4, 11.9, 12.4, 11.1, 14.5, 10.3, 12.4,
		12.7, 11.9, 10.3, 13.7, 15.5, 14.5, 11.6, 10.6, 15.5, 14.7, 8.8, 11.6, 8.3,
		17.6, 12.4, 11.6, 15.0, 13.7, 12.7, 10.9, 7.2, 8.5, 8.3, 9.6, 11.4, 12.9,
		11.6, 13.4, 10.1, 11.6, 8.8, 12.4, 10.3, 16.3, 10.9, 10.1, 8.8, 9.3, 15.2,
		8.5, 11.1, 8.3, 11.4, 11.9, 9.3, 9.8, 16.3, 12.7, 9.0, 11.9, 9.3, 10.3, 13.4,
		11.4, 12.9, 12.4, 9.6, 10.3, 13.2, 10.6, 9.8, 11.9, 14.2, 13.4, 9.3, 9.6,
		12.1, 11.9, 10.1, 14.0, 12.9, 21.7, 11.6, 12.1, 10.3, 9.8, 14.2, 13.7, 7.2,
		10.9, 10.1, 9.6, 13.4, 13.2, 16.3, 11.6, 14.0, 10.9, 14.2, 12.4, 12.4, 13.4,
		17.6, 10.1, 10.9, 14.0, 12.9, 9.0, 13.4, 15.0, 16.0, 8.0, 9.8, 12.4, 8.5,
		9.6, 12.7, 12.1, 15.0, 16.0, 10.9, 14.2, 13.7, 11.9, 16.8, 11.4, 13.4, 12.4,
		22.0, 12.4, 16.8, 12.1, 10.3, 13.4, 11.6, 10.1, 14.5, 10.6, 11.9, 15.5, 9.8,
		12.4, 10.1, 8.0, 9.0, 9.3, 13.2, 11.1, 12.7, 12.1, 10.1, 13.2, 14.5, 10.1,
		12.7, 12.9, 11.9, 12.4, 11.1, 8.5, 14.5, 16.5, 12.4, 9.0, 11.1, 9.8, 11.1,
		11.1, 8.8, 13.2, 17.6, 16.8, 10.9, 12.4, 8.5, 14.7 };
    */

	double b_x[266];
	int j;
	double v[12];
	double b_v;
	double fv[4];
	int idx[4];
	int colIdx;
	boolean_T doShrink;
	int itercount;
	int fun_evals;
	int lastCol;
	int firstCol;
	boolean_T exitg1;
	double xr[3];
	double maxfv;
	double d2;
	boolean_T p;
	int b_firstCol;
	double xbar[3];
	double d3;
	double xe[3];
	double fvt[4];
	int idxb[4];
	for (k = 0; k < 4; k++) {
		for (j = 0; j < 3; j++) {
			v[j + 3 * k] = x[j];
		}
	}

	// PDF CALL HERE
	//emgpdf(data, x[0], x[1], x[2], b_x);
	pdf(data, x[0], x[1], x[2], b_x);
	for (k = 0; k < 266; k++) {
		b_x[k] = -log(b_x[k]);
	}

	b_v = b_x[0];
	for (k = 0; k < 265; k++) {
		b_v += b_x[k + 1];
	}

	fv[0] = b_v;
	for (k = 0; k < 3; k++) {
		colIdx = 3 * (k + 1);
		if (x[k] != 0.0) {
			v[k + colIdx] = 1.05 * x[k];
		}
		else {
			v[k + colIdx] = 0.00025;
		}

		for (j = 0; j < 3; j++) {
			xr[j] = v[j + colIdx];
		}

		// PDF CALL HERE
		//emgpdf(data, xr[0], xr[1], xr[2], b_x);
		pdf(data, xr[0], xr[1], xr[2], b_x);
		for (lastCol = 0; lastCol < 266; lastCol++) {
			b_x[lastCol] = -log(b_x[lastCol]);
		}

		b_v = b_x[0];
		for (lastCol = 0; lastCol < 265; lastCol++) {
			b_v += b_x[lastCol + 1];
		}

		fv[k + 1] = b_v;
	}

	sortIdx(fv, idx);
	doShrink = false;
	itercount = 1;
	fun_evals = 4;
	lastCol = 3 * (idx[3] - 1);
	firstCol = 3 * (idx[0] - 1);
	exitg1 = false;
	while ((!exitg1) && ((fun_evals < 10000) && (itercount < 10000))) {
		maxfv = 0.0;
		for (k = 0; k < 3; k++) {
			b_v = fabs(fv[idx[0] - 1] - fv[idx[k + 1] - 1]);
			if (b_v > maxfv) {
				maxfv = b_v;
			}
		}

		b_v = 10.0 * eps(10.0 * eps(fv[idx[0] - 1]));
		if ((0.0001 > b_v) || rtIsNaN(b_v)) {
			d2 = 0.0001;
		}
		else {
			d2 = b_v;
		}

		if (maxfv > d2) {
			p = false;
		}
		else {
			maxfv = 0.0;
			b_firstCol = (idx[0] - 1) * 3;
			for (j = 0; j < 3; j++) {
				colIdx = (idx[j + 1] - 1) * 3;
				for (k = 0; k < 3; k++) {
					b_v = fabs(v[k + b_firstCol] - v[k + colIdx]);
					if (b_v > maxfv) {
						maxfv = b_v;
					}
				}
			}

			b_v = v[b_firstCol];
			for (k = 0; k < 2; k++) {
				if (v[(k + b_firstCol) + 1] > b_v) {
					b_v = v[(k + b_firstCol) + 1];
				}
			}

			b_v = 10.0 * eps(b_v);
			if ((0.0001 > b_v) || rtIsNaN(b_v)) {
				d3 = 0.0001;
			}
			else {
				d3 = b_v;
			}

			p = (maxfv <= d3);
		}

		if (!p) {
			for (k = 0; k < 3; k++) {
				xbar[k] = v[k + firstCol];
			}

			for (k = 0; k < 2; k++) {
				colIdx = 3 * (idx[k + 1] - 1);
				for (j = 0; j < 3; j++) {
					xbar[j] += v[j + colIdx];
				}
			}

			for (k = 0; k < 3; k++) {
				b_v = xbar[k] / 3.0;
				xr[k] = 2.0 * b_v - v[k + lastCol];
				xbar[k] = b_v;
			}

			// PDF CALL HERE
			//emgpdf(data, xr[0], xr[1], xr[2], b_x);
			pdf(data, xr[0], xr[1], xr[2], b_x);
			for (k = 0; k < 266; k++) {
				b_x[k] = -log(b_x[k]);
			}

			b_v = b_x[0];
			for (k = 0; k < 265; k++) {
				b_v += b_x[k + 1];
			}

			fun_evals++;
			if (b_v < fv[idx[0] - 1]) {
				for (k = 0; k < 3; k++) {
					xe[k] = 3.0 * xbar[k] - 2.0 * v[k + lastCol];
				}

				// PDF CALL HERE
				//emgpdf(data, xe[0], xe[1], xe[2], b_x);
				pdf(data, xe[0], xe[1], xe[2], b_x);
				for (k = 0; k < 266; k++) {
					b_x[k] = -log(b_x[k]);
				}

				maxfv = b_x[0];
				for (k = 0; k < 265; k++) {
					maxfv += b_x[k + 1];
				}

				fun_evals++;
				if (maxfv < b_v) {
					for (k = 0; k < 3; k++) {
						v[k + lastCol] = xe[k];
					}

					fv[idx[3] - 1] = maxfv;
				}
				else {
					for (k = 0; k < 3; k++) {
						v[k + lastCol] = xr[k];
					}

					fv[idx[3] - 1] = b_v;
				}
			}
			else if (b_v < fv[idx[2] - 1]) {
				for (k = 0; k < 3; k++) {
					v[k + lastCol] = xr[k];
				}

				fv[idx[3] - 1] = b_v;
			}
			else {
				if (b_v < fv[idx[3] - 1]) {
					for (k = 0; k < 3; k++) {
						x[k] = 1.5 * xbar[k] - 0.5 * v[k + lastCol];
					}

					// PDF CALL HERE
					//emgpdf(data, x[0], x[1], x[2], b_x);
					pdf(data, x[0], x[1], x[2], b_x);
					for (k = 0; k < 266; k++) {
						b_x[k] = -log(b_x[k]);
					}

					maxfv = b_x[0];
					for (k = 0; k < 265; k++) {
						maxfv += b_x[k + 1];
					}

					fun_evals++;
					if (maxfv <= b_v) {
						for (k = 0; k < 3; k++) {
							v[k + lastCol] = x[k];
						}

						fv[idx[3] - 1] = maxfv;
					}
					else {
						doShrink = true;
					}
				}
				else {
					for (k = 0; k < 3; k++) {
						x[k] = 0.5 * xbar[k] + 0.5 * v[k + lastCol];
					}

					// PDF CALL HERE
					//emgpdf(data, x[0], x[1], x[2], b_x);
					pdf(data, x[0], x[1], x[2], b_x);
					for (k = 0; k < 266; k++) {
						b_x[k] = -log(b_x[k]);
					}

					b_v = b_x[0];
					for (k = 0; k < 265; k++) {
						b_v += b_x[k + 1];
					}

					fun_evals++;
					if (b_v < fv[idx[3] - 1]) {
						for (k = 0; k < 3; k++) {
							v[k + lastCol] = x[k];
						}

						fv[idx[3] - 1] = b_v;
					}
					else {
						doShrink = true;
					}
				}

				if (doShrink) {
					for (k = 0; k < 3; k++) {
						colIdx = 3 * (idx[k + 1] - 1);
						for (j = 0; j < 3; j++) {
							v[j + colIdx] = v[j + firstCol] + 0.5 * (v[j + colIdx] - v[j +
								firstCol]);
							x[j] = v[j + colIdx];
						}

						//emgpdf(data, x[0], x[1], x[2], b_x);
						pdf(data, x[0], x[1], x[2], b_x);
						for (lastCol = 0; lastCol < 266; lastCol++) {
							b_x[lastCol] = -log(b_x[lastCol]);
						}

						b_v = b_x[0];
						for (lastCol = 0; lastCol < 265; lastCol++) {
							b_v += b_x[lastCol + 1];
						}

						fv[idx[k + 1] - 1] = b_v;
					}

					fun_evals += 3;
				}
			}

			if (doShrink) {
				for (k = 0; k < 4; k++) {
					fvt[k] = fv[idx[k] - 1];
					idxb[k] = idx[k];
				}

				sortIdx(fvt, idx);
				for (k = 0; k < 4; k++) {
					idx[k] = idxb[idx[k] - 1];
				}

				doShrink = false;
			}
			else {
				for (k = 0; k < 3; k++) {
					if (fv[idx[3 - k] - 1] < fv[idx[2 - k] - 1]) {
						lastCol = idx[3 - k];
						idx[3 - k] = idx[2 - k];
						idx[2 - k] = lastCol;
					}
				}
			}

			itercount++;
			lastCol = 3 * (idx[3] - 1);
			firstCol = 3 * (idx[0] - 1);
		}
		else {
			exitg1 = true;
		}
	}

	colIdx = 3 * (idx[0] - 1);
	for (k = 0; k < 3; k++) {
		x[k] = v[k + colIdx];
	}

	return fv[idx[0] - 1];
}

/* End of code generation (fminsearch.c) */
