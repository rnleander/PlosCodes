/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * roots.c
 *
 * Code generation for function 'roots'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "roots.h"
#include "xzhseqr.h"
#include "xgehrd.h"
#include "anyNonFinite.h"
#include "rdivide.h"
#include "xzgeev.h"

/* Function Definitions */

/*
 *
 */
void roots(const double c[5], creal_T r_data[], int r_size[1])
{
  int k1;
  int k2;
  int companDim;
  boolean_T exitg1;
  int j;
  int a_size[2];
  boolean_T exitg2;
  int istart;
  double ctmp[5];
  creal_T a_data[16];
  boolean_T p;
  creal_T eiga_data[4];
  creal_T alpha1_data[4];
  int alpha1_size[1];
  creal_T beta1_data[4];
  int beta1_size[1];
  int exitg3;
  signed char iv0[2];
  int jend;
  memset(&r_data[0], 0, sizeof(creal_T) << 2);
  k1 = 1;
  while ((k1 <= 5) && (!(c[k1 - 1] != 0.0))) {
    k1++;
  }

  k2 = 5;
  while ((k2 >= k1) && (!(c[k2 - 1] != 0.0))) {
    k2--;
  }

  if (k1 < k2) {
    companDim = k2 - k1;
    exitg1 = false;
    while ((!exitg1) && (companDim > 0)) {
      j = 0;
      exitg2 = false;
      while ((!exitg2) && (j + 1 <= companDim)) {
        ctmp[j] = c[k1 + j] / c[k1 - 1];
        if (rtIsInf(fabs(ctmp[j]))) {
          exitg2 = true;
        } else {
          j++;
        }
      }

      if (j + 1 > companDim) {
        exitg1 = true;
      } else {
        k1++;
        companDim--;
      }
    }

    if (companDim < 1) {
      if (1 > 5 - k2) {
        r_size[0] = 0;
      } else {
        r_size[0] = 5 - k2;
      }
    } else {
      a_size[0] = companDim;
      a_size[1] = companDim;
      k1 = companDim * companDim;
      for (istart = 0; istart < k1; istart++) {
        a_data[istart].re = 0.0;
        a_data[istart].im = 0.0;
      }

      for (k1 = 0; k1 + 1 < companDim; k1++) {
        a_data[companDim * k1].re = -ctmp[k1];
        a_data[companDim * k1].im = 0.0;
        a_data[(k1 + companDim * k1) + 1].re = 1.0;
        a_data[(k1 + companDim * k1) + 1].im = 0.0;
      }

      a_data[companDim * (companDim - 1)].re = -ctmp[companDim - 1];
      a_data[companDim * (companDim - 1)].im = 0.0;
      for (k1 = 1; k1 <= 5 - k2; k1++) {
        r_data[k1 - 1].re = 0.0;
        r_data[k1 - 1].im = 0.0;
      }

      if (anyNonFinite(a_data, a_size)) {
        if (companDim == 1) {
          eiga_data[0].re = rtNaN;
          eiga_data[0].im = 0.0;
        } else {
          for (istart = 0; istart < companDim; istart++) {
            eiga_data[istart].re = rtNaN;
            eiga_data[istart].im = 0.0;
          }
        }
      } else if (companDim == 1) {
        eiga_data[0] = a_data[0];
      } else {
        p = true;
        j = 0;
        exitg1 = false;
        while ((!exitg1) && (j <= companDim - 1)) {
          k1 = 0;
          do {
            exitg3 = 0;
            if (k1 <= j) {
              if (!((a_data[k1 + companDim * j].re == a_data[j + companDim * k1]
                     .re) && (a_data[k1 + companDim * j].im == -a_data[j +
                              companDim * k1].im))) {
                p = false;
                exitg3 = 1;
              } else {
                k1++;
              }
            } else {
              j++;
              exitg3 = 2;
            }
          } while (exitg3 == 0);

          if (exitg3 == 1) {
            exitg1 = true;
          }
        }

        if (p) {
          if (anyNonFinite(a_data, a_size)) {
            for (istart = 0; istart < 2; istart++) {
              iv0[istart] = (signed char)a_size[istart];
            }

            a_size[0] = iv0[0];
            k1 = iv0[0] * iv0[1];
            for (istart = 0; istart < k1; istart++) {
              a_data[istart].re = rtNaN;
              a_data[istart].im = 0.0;
            }

            if (!(1 >= iv0[0])) {
              istart = 2;
              if (iv0[0] - 2 < iv0[1] - 1) {
                jend = iv0[0] - 1;
              } else {
                jend = iv0[1];
              }

              for (j = 1; j <= jend; j++) {
                for (k1 = istart; k1 <= a_size[0]; k1++) {
                  a_data[(k1 + a_size[0] * (j - 1)) - 1].re = 0.0;
                  a_data[(k1 + a_size[0] * (j - 1)) - 1].im = 0.0;
                }

                istart++;
              }
            }
          } else {
            xgehrd(a_data, a_size);
            eml_zlahqr(a_data, a_size);
            if (!(3 >= a_size[0])) {
              a_data[3].re = 0.0;
              a_data[3].im = 0.0;
            }
          }

          for (k1 = 0; k1 + 1 <= a_size[0]; k1++) {
            eiga_data[k1] = a_data[k1 + a_size[0] * k1];
          }
        } else {
          xzgeev(a_data, a_size, &k1, alpha1_data, alpha1_size, beta1_data,
                 beta1_size);
          rdivide(alpha1_data, alpha1_size, beta1_data, eiga_data, beta1_size);
        }
      }

      for (k1 = 1; k1 <= companDim; k1++) {
        r_data[(k1 - k2) + 4] = eiga_data[k1 - 1];
      }

      r_size[0] = (companDim - k2) + 5;
    }
  } else if (1 > 5 - k2) {
    r_size[0] = 0;
  } else {
    r_size[0] = 5 - k2;
  }
}

/* End of code generation (roots.c) */
