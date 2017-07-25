/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xscal.c
 *
 * Code generation for function 'xscal'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "xscal.h"

/* Function Definitions */

/*
 *
 */
void b_xscal(int n, const creal_T a, creal_T x_data[], int ix0, int incx)
{
  int i5;
  int k;
  double x_data_re;
  double x_data_im;
  i5 = ix0 + incx * (n - 1);
  for (k = ix0; k <= i5; k += incx) {
    x_data_re = x_data[k - 1].re;
    x_data_im = x_data[k - 1].im;
    x_data[k - 1].re = a.re * x_data_re - a.im * x_data_im;
    x_data[k - 1].im = a.re * x_data_im + a.im * x_data_re;
  }
}

/*
 *
 */
void xscal(int n, const creal_T a, creal_T x_data[], int ix0)
{
  int i3;
  int k;
  double x_data_re;
  double x_data_im;
  i3 = (ix0 + n) - 1;
  for (k = ix0; k <= i3; k++) {
    x_data_re = x_data[k - 1].re;
    x_data_im = x_data[k - 1].im;
    x_data[k - 1].re = a.re * x_data_re - a.im * x_data_im;
    x_data[k - 1].im = a.re * x_data_im + a.im * x_data_re;
  }
}

/* End of code generation (xscal.c) */
