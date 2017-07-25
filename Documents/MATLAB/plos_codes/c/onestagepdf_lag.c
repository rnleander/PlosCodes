/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * onestagepdf_lag.c
 *
 * Code generation for function 'onestagepdf_lag'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "onestagepdf_lag.h"
#include "emgpdf.h"
#include "IMT_analysis_April2017_rtwutil.h"

/* Function Definitions */

/*
 * function Y=onestagepdf_lag(X,m,s,l)
 */
void b_onestagepdf_lag(const double X[221], double m, double s, double l, double
  Y[221])
{
  int k;
  double a;
  double y[221];
  double b_a[221];
  double b_y;

  /* this function evaluates a shifted inverse gaussian distribution at X */
  /* m=mu */
  /* s=sigma */
  /* l=lag */
  /* 'onestagepdf_lag:7' Y=(1./(s*(2*pi*(X-l).^3).^(1/2))).*exp(-(m*(X-l)-1).^2./(2*s^2*(X-l))); */
  for (k = 0; k < 221; k++) {
    y[k] = rt_powd_snf(6.2831853071795862 * rt_powd_snf(X[k] - l, 3.0), 0.5);
    b_a[k] = m * (X[k] - l) - 1.0;
  }

  a = 2.0 * (s * s);

  /* If some cpmponent of X is less than l, the pdf will return imaginary values,  */
  /* in this case we stipilate that the pdf returns realmin */
  /* 'onestagepdf_lag:10' Y=real(Y); */
  /* 'onestagepdf_lag:11' Y=max(Y, realmin); */
  for (k = 0; k < 221; k++) {
    b_y = 1.0 / (s * y[k]) * exp(-rt_powd_snf(b_a[k], 2.0) / (a * (X[k] - l)));
    if (!(b_y > 2.2250738585072014E-308)) {
      b_y = 2.2250738585072014E-308;
    }

    Y[k] = b_y;
  }
}

/*
 * function Y=onestagepdf_lag(X,m,s,l)
 */
void c_onestagepdf_lag(const double X[2201], double m, double s, double l,
  double Y[2201])
{
  int k;
  double a;
  double y[2201];
  double b_a[2201];
  double b_y;

  /* this function evaluates a shifted inverse gaussian distribution at X */
  /* m=mu */
  /* s=sigma */
  /* l=lag */
  /* 'onestagepdf_lag:7' Y=(1./(s*(2*pi*(X-l).^3).^(1/2))).*exp(-(m*(X-l)-1).^2./(2*s^2*(X-l))); */
  for (k = 0; k < 2201; k++) {
    y[k] = rt_powd_snf(6.2831853071795862 * rt_powd_snf(X[k] - l, 3.0), 0.5);
    b_a[k] = m * (X[k] - l) - 1.0;
  }

  a = 2.0 * (s * s);

  /* If some cpmponent of X is less than l, the pdf will return imaginary values,  */
  /* in this case we stipilate that the pdf returns realmin */
  /* 'onestagepdf_lag:10' Y=real(Y); */
  /* 'onestagepdf_lag:11' Y=max(Y, realmin); */
  for (k = 0; k < 2201; k++) {
    b_y = 1.0 / (s * y[k]) * exp(-rt_powd_snf(b_a[k], 2.0) / (a * (X[k] - l)));
    if (!(b_y > 2.2250738585072014E-308)) {
      b_y = 2.2250738585072014E-308;
    }

    Y[k] = b_y;
  }
}

/*
 * function Y=onestagepdf_lag(X,m,s,l)
 */
void d_onestagepdf_lag(const double X[22001], double m, double s, double l,
  double Y[22001])
{
  int k;
  double a;
  static double y[22001];
  static double b_a[22001];
  double b_y;

  /* this function evaluates a shifted inverse gaussian distribution at X */
  /* m=mu */
  /* s=sigma */
  /* l=lag */
  /* 'onestagepdf_lag:7' Y=(1./(s*(2*pi*(X-l).^3).^(1/2))).*exp(-(m*(X-l)-1).^2./(2*s^2*(X-l))); */
  for (k = 0; k < 22001; k++) {
    y[k] = rt_powd_snf(6.2831853071795862 * rt_powd_snf(X[k] - l, 3.0), 0.5);
    b_a[k] = m * (X[k] - l) - 1.0;
  }

  a = 2.0 * (s * s);

  /* If some cpmponent of X is less than l, the pdf will return imaginary values,  */
  /* in this case we stipilate that the pdf returns realmin */
  /* 'onestagepdf_lag:10' Y=real(Y); */
  /* 'onestagepdf_lag:11' Y=max(Y, realmin); */
  for (k = 0; k < 22001; k++) {
    b_y = 1.0 / (s * y[k]) * exp(-rt_powd_snf(b_a[k], 2.0) / (a * (X[k] - l)));
    if (!(b_y > 2.2250738585072014E-308)) {
      b_y = 2.2250738585072014E-308;
    }

    Y[k] = b_y;
  }
}

/*
 * function Y=onestagepdf_lag(X,m,s,l)
 */
void onestagepdf_lag(const double X[266], double m, double s, double l, double
                     Y[266])
{
  int k;
  double a;
  double y[266];
  double b_a[266];
  double b_y;

  /* this function evaluates a shifted inverse gaussian distribution at X */
  /* m=mu */
  /* s=sigma */
  /* l=lag */
  /* 'onestagepdf_lag:7' Y=(1./(s*(2*pi*(X-l).^3).^(1/2))).*exp(-(m*(X-l)-1).^2./(2*s^2*(X-l))); */
  for (k = 0; k < 266; k++) {
    y[k] = rt_powd_snf(6.2831853071795862 * rt_powd_snf(X[k] - l, 3.0), 0.5);
    b_a[k] = m * (X[k] - l) - 1.0;
  }

  a = 2.0 * (s * s);

  /* If some cpmponent of X is less than l, the pdf will return imaginary values,  */
  /* in this case we stipilate that the pdf returns realmin */
  /* 'onestagepdf_lag:10' Y=real(Y); */
  /* 'onestagepdf_lag:11' Y=max(Y, realmin); */
  for (k = 0; k < 266; k++) {
    b_y = 1.0 / (s * y[k]) * exp(-rt_powd_snf(b_a[k], 2.0) / (a * (X[k] - l)));
    if (!(b_y > 2.2250738585072014E-308)) {
      b_y = 2.2250738585072014E-308;
    }

    Y[k] = b_y;
  }
}

/* End of code generation (onestagepdf_lag.c) */