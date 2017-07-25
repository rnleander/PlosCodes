/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xgehrd.c
 *
 * Code generation for function 'xgehrd'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "xgehrd.h"
#include "recip.h"
#include "xdlapy3.h"
#include "xnrm2.h"

/* Function Definitions */

/*
 *
 */
void xgehrd(creal_T a_data[], int a_size[2])
{
  int n;
  int ia0;
  int i2;
  int i;
  creal_T work_data[4];
  int im1n;
  int in;
  creal_T alpha1;
  int c;
  double c_re;
  double c_im;
  double xnorm;
  creal_T tau_data[3];
  double beta1;
  int jy;
  int knt;
  double temp_im;
  int lastv;
  int lastc;
  int k;
  boolean_T exitg2;
  double a_data_im;
  int ix;
  int exitg1;
  creal_T b_alpha1;
  n = a_size[0];
  ia0 = (signed char)a_size[0];
  for (i2 = 0; i2 < ia0; i2++) {
    work_data[i2].re = 0.0;
    work_data[i2].im = 0.0;
  }

  for (i = 0; i + 1 < n; i++) {
    im1n = i * n + 2;
    in = (i + 1) * n;
    alpha1 = a_data[(i + a_size[0] * i) + 1];
    ia0 = i + 3;
    if (!(ia0 < n)) {
      ia0 = n;
    }

    ia0 += i * n;
    c = (n - i) - 3;
    c_re = 0.0;
    c_im = 0.0;
    if (!(c + 2 <= 0)) {
      xnorm = xnrm2(c + 1, a_data, ia0);
      if ((xnorm != 0.0) || (a_data[(i + a_size[0] * i) + 1].im != 0.0)) {
        beta1 = xdlapy3(a_data[(i + a_size[0] * i) + 1].re, a_data[(i + a_size[0]
          * i) + 1].im, xnorm);
        if (a_data[(i + a_size[0] * i) + 1].re >= 0.0) {
          beta1 = -beta1;
        }

        if (fabs(beta1) < 1.0020841800044864E-292) {
          knt = 0;
          i2 = ia0 + c;
          do {
            knt++;
            for (k = ia0; k <= i2; k++) {
              xnorm = a_data[k - 1].re;
              a_data_im = a_data[k - 1].im;
              a_data[k - 1].re = 9.9792015476736E+291 * xnorm - 0.0 * a_data_im;
              a_data[k - 1].im = 9.9792015476736E+291 * a_data_im + 0.0 * xnorm;
            }

            beta1 *= 9.9792015476736E+291;
            alpha1.re *= 9.9792015476736E+291;
            alpha1.im *= 9.9792015476736E+291;
          } while (!(fabs(beta1) >= 1.0020841800044864E-292));

          xnorm = xnrm2(c + 1, a_data, ia0);
          beta1 = xdlapy3(alpha1.re, alpha1.im, xnorm);
          if (alpha1.re >= 0.0) {
            beta1 = -beta1;
          }

          xnorm = beta1 - alpha1.re;
          if (0.0 - alpha1.im == 0.0) {
            c_re = xnorm / beta1;
            c_im = 0.0;
          } else if (xnorm == 0.0) {
            c_re = 0.0;
            c_im = (0.0 - alpha1.im) / beta1;
          } else {
            c_re = xnorm / beta1;
            c_im = (0.0 - alpha1.im) / beta1;
          }

          b_alpha1.re = alpha1.re - beta1;
          b_alpha1.im = alpha1.im;
          alpha1 = recip(b_alpha1);
          i2 = ia0 + c;
          while (ia0 <= i2) {
            xnorm = a_data[ia0 - 1].re;
            a_data_im = a_data[ia0 - 1].im;
            a_data[ia0 - 1].re = alpha1.re * xnorm - alpha1.im * a_data_im;
            a_data[ia0 - 1].im = alpha1.re * a_data_im + alpha1.im * xnorm;
            ia0++;
          }

          for (k = 1; k <= knt; k++) {
            beta1 *= 1.0020841800044864E-292;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0;
        } else {
          xnorm = beta1 - a_data[(i + a_size[0] * i) + 1].re;
          temp_im = 0.0 - a_data[(i + a_size[0] * i) + 1].im;
          if (temp_im == 0.0) {
            c_re = xnorm / beta1;
            c_im = 0.0;
          } else if (xnorm == 0.0) {
            c_re = 0.0;
            c_im = temp_im / beta1;
          } else {
            c_re = xnorm / beta1;
            c_im = temp_im / beta1;
          }

          alpha1.re = a_data[(i + a_size[0] * i) + 1].re - beta1;
          alpha1.im = a_data[(i + a_size[0] * i) + 1].im;
          alpha1 = recip(alpha1);
          i2 = ia0 + c;
          while (ia0 <= i2) {
            xnorm = a_data[ia0 - 1].re;
            a_data_im = a_data[ia0 - 1].im;
            a_data[ia0 - 1].re = alpha1.re * xnorm - alpha1.im * a_data_im;
            a_data[ia0 - 1].im = alpha1.re * a_data_im + alpha1.im * xnorm;
            ia0++;
          }

          alpha1.re = beta1;
          alpha1.im = 0.0;
        }
      }
    }

    tau_data[i].re = c_re;
    tau_data[i].im = c_im;
    a_data[(i + a_size[0] * i) + 1].re = 1.0;
    a_data[(i + a_size[0] * i) + 1].im = 0.0;
    c = (n - i) - 3;
    jy = (i + im1n) - 1;
    if ((tau_data[i].re != 0.0) || (tau_data[i].im != 0.0)) {
      lastv = c + 2;
      ia0 = jy + c;
      while ((lastv > 0) && ((a_data[ia0 + 1].re == 0.0) && (a_data[ia0 + 1].im ==
               0.0))) {
        lastv--;
        ia0--;
      }

      lastc = n;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        ia0 = in + lastc;
        c = ia0;
        do {
          exitg1 = 0;
          if (c <= ia0 + (lastv - 1) * n) {
            if ((a_data[c - 1].re != 0.0) || (a_data[c - 1].im != 0.0)) {
              exitg1 = 1;
            } else {
              c += n;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = 0;
    }

    if (lastv > 0) {
      if (lastc != 0) {
        for (ia0 = 1; ia0 <= lastc; ia0++) {
          work_data[ia0 - 1].re = 0.0;
          work_data[ia0 - 1].im = 0.0;
        }

        ix = jy;
        i2 = (in + n * (lastv - 1)) + 1;
        for (knt = in + 1; knt <= i2; knt += n) {
          c_re = a_data[ix].re - 0.0 * a_data[ix].im;
          c_im = a_data[ix].im + 0.0 * a_data[ix].re;
          ia0 = 0;
          k = (knt + lastc) - 1;
          for (c = knt; c <= k; c++) {
            work_data[ia0].re += a_data[c - 1].re * c_re - a_data[c - 1].im *
              c_im;
            work_data[ia0].im += a_data[c - 1].re * c_im + a_data[c - 1].im *
              c_re;
            ia0++;
          }

          ix++;
        }
      }

      c_re = -tau_data[i].re;
      c_im = -tau_data[i].im;
      if (!((-tau_data[i].re == 0.0) && (-tau_data[i].im == 0.0))) {
        ia0 = in;
        for (knt = 1; knt <= lastv; knt++) {
          if ((a_data[jy].re != 0.0) || (a_data[jy].im != 0.0)) {
            xnorm = a_data[jy].re * c_re + a_data[jy].im * c_im;
            temp_im = a_data[jy].re * c_im - a_data[jy].im * c_re;
            ix = 0;
            i2 = lastc + ia0;
            for (k = ia0; k + 1 <= i2; k++) {
              a_data[k].re += work_data[ix].re * xnorm - work_data[ix].im *
                temp_im;
              a_data[k].im += work_data[ix].re * temp_im + work_data[ix].im *
                xnorm;
              ix++;
            }
          }

          jy++;
          ia0 += n;
        }
      }
    }

    c = (n - i) - 3;
    im1n = (i + im1n) - 1;
    jy = (i + in) + 2;
    if ((tau_data[i].re != 0.0) || (-tau_data[i].im != 0.0)) {
      lastv = c + 2;
      ia0 = im1n + c;
      while ((lastv > 0) && ((a_data[ia0 + 1].re == 0.0) && (a_data[ia0 + 1].im ==
               0.0))) {
        lastv--;
        ia0--;
      }

      lastc = (n - i) - 1;
      exitg2 = false;
      while ((!exitg2) && (lastc > 0)) {
        ia0 = jy + (lastc - 1) * n;
        c = ia0;
        do {
          exitg1 = 0;
          if (c <= (ia0 + lastv) - 1) {
            if ((a_data[c - 1].re != 0.0) || (a_data[c - 1].im != 0.0)) {
              exitg1 = 1;
            } else {
              c++;
            }
          } else {
            lastc--;
            exitg1 = 2;
          }
        } while (exitg1 == 0);

        if (exitg1 == 1) {
          exitg2 = true;
        }
      }
    } else {
      lastv = 0;
      lastc = 0;
    }

    if (lastv > 0) {
      if (lastc != 0) {
        for (ia0 = 1; ia0 <= lastc; ia0++) {
          work_data[ia0 - 1].re = 0.0;
          work_data[ia0 - 1].im = 0.0;
        }

        ia0 = 0;
        i2 = jy + n * (lastc - 1);
        for (knt = jy; knt <= i2; knt += n) {
          ix = im1n;
          c_re = 0.0;
          c_im = 0.0;
          k = (knt + lastv) - 1;
          for (c = knt - 1; c + 1 <= k; c++) {
            c_re += a_data[c].re * a_data[ix].re + a_data[c].im * a_data[ix].im;
            c_im += a_data[c].re * a_data[ix].im - a_data[c].im * a_data[ix].re;
            ix++;
          }

          work_data[ia0].re += c_re - 0.0 * c_im;
          work_data[ia0].im += c_im + 0.0 * c_re;
          ia0++;
        }
      }

      c_re = -tau_data[i].re;
      c_im = -(-tau_data[i].im);
      if (!((-tau_data[i].re == 0.0) && (-(-tau_data[i].im) == 0.0))) {
        ia0 = jy - 1;
        jy = 0;
        for (knt = 1; knt <= lastc; knt++) {
          if ((work_data[jy].re != 0.0) || (work_data[jy].im != 0.0)) {
            xnorm = work_data[jy].re * c_re + work_data[jy].im * c_im;
            temp_im = work_data[jy].re * c_im - work_data[jy].im * c_re;
            ix = im1n;
            i2 = lastv + ia0;
            for (k = ia0; k + 1 <= i2; k++) {
              a_data_im = a_data[ix].re * temp_im + a_data[ix].im * xnorm;
              a_data[k].re += a_data[ix].re * xnorm - a_data[ix].im * temp_im;
              a_data[k].im += a_data_im;
              ix++;
            }
          }

          jy++;
          ia0 += n;
        }
      }
    }

    a_data[(i + a_size[0] * i) + 1] = alpha1;
  }
}

/* End of code generation (xgehrd.c) */
