/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * gp_max.c
 *
 * Code generation for function 'gp_max'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "gp_max.h"
#include "onestagepdf_prime.h"
#include "roots.h"
#include "emgpdf.h"
#include "IMT_analysis_April2017_rtwutil.h"

/* Function Definitions */

/*
 * function M=gp_max(m,s)
 */
double gp_max(double m, double s)
{
  double M;
  double dv64[5];
  creal_T z_data[4];
  int z_size[1];
  int ixstart;
  int n;
  double x_data[4];
  double y_data[4];
  boolean_T tmp_data[4];
  int i;
  int partialTrueCount;
  creal_T M_data[4];
  int M_size[1];
  boolean_T exitg1;

  /* this function finds the maximum value of the first derivative of an */
  /* inverse Gaussian distribution with parameters m and s */
  /*  the second derivitive of inverse gaussian can be represented as a */
  /*  polynomial with the following coefficents, a through e */
  /* 'gp_max:8' a=m^4/(4*s^4); */
  /* 'gp_max:9' b=3*m^2/(2*s^2); */
  /* 'gp_max:10' c=15/4-m^2/(2*s^4); */
  /* 'gp_max:11' d=-5/(2*s^2); */
  /* 'gp_max:12' e=1/(4*s^4); */
  /*  find the roots of the second derivitive */
  /* 'gp_max:15' z=roots([a b c d e]); */
  dv64[0] = rt_powd_snf(m, 4.0) / (4.0 * rt_powd_snf(s, 4.0));
  dv64[1] = 3.0 * (m * m) / (2.0 * (s * s));
  dv64[2] = 3.75 - m * m / (2.0 * rt_powd_snf(s, 4.0));
  dv64[3] = -5.0 / (2.0 * (s * s));
  dv64[4] = 1.0 / (4.0 * rt_powd_snf(s, 4.0));
  roots(dv64, z_data, z_size);

  /*  we are only interested in the positive real roots */
  /* z=z(imag(z)==0); */
  /* 'gp_max:19' z=z(abs(imag(z))<=0.000000000000001); */
  ixstart = z_size[0];
  for (n = 0; n < ixstart; n++) {
    x_data[n] = z_data[n].im;
  }

  for (ixstart = 0; ixstart + 1 <= z_size[0]; ixstart++) {
    y_data[ixstart] = fabs(x_data[ixstart]);
  }

  ixstart = (signed char)z_size[0];
  for (n = 0; n < ixstart; n++) {
    tmp_data[n] = (y_data[n] <= 1.0E-15);
  }

  ixstart = (signed char)z_size[0] - 1;
  n = 0;
  for (i = 0; i <= ixstart; i++) {
    if (tmp_data[i]) {
      n++;
    }
  }

  partialTrueCount = 0;
  for (i = 0; i <= ixstart; i++) {
    if (tmp_data[i]) {
      z_data[partialTrueCount] = z_data[i];
      partialTrueCount++;
    }
  }

  /* 'gp_max:20' z=z(z>=0); */
  ixstart = n - 1;
  n = 0;
  for (i = 0; i <= ixstart; i++) {
    if (z_data[i].re >= 0.0) {
      n++;
    }
  }

  partialTrueCount = 0;
  for (i = 0; i <= ixstart; i++) {
    if (z_data[i].re >= 0.0) {
      z_data[partialTrueCount] = z_data[i];
      partialTrueCount++;
    }
  }

  z_size[0] = n;

  /*  find the first derivitive at each of these roots */
  /* 'gp_max:23' M=onestagepdf_prime(z,m,s); */
  onestagepdf_prime(z_data, z_size, m, s, M_data, M_size);

  /*  return the maximu value of the first derivative */
  /* 'gp_max:26' M=max(abs(M)); */
  for (ixstart = 0; ixstart + 1 <= M_size[0]; ixstart++) {
    y_data[ixstart] = rt_hypotd_snf(M_data[ixstart].re, M_data[ixstart].im);
  }

  ixstart = 1;
  n = (signed char)M_size[0];
  M = y_data[0];
  if ((signed char)M_size[0] > 1) {
    if (rtIsNaN(y_data[0])) {
      i = 2;
      exitg1 = false;
      while ((!exitg1) && (i <= n)) {
        ixstart = i;
        if (!rtIsNaN(y_data[i - 1])) {
          M = y_data[i - 1];
          exitg1 = true;
        } else {
          i++;
        }
      }
    }

    if (ixstart < (signed char)M_size[0]) {
      while (ixstart + 1 <= n) {
        if (y_data[ixstart] > M) {
          M = y_data[ixstart];
        }

        ixstart++;
      }
    }
  }

  return M;
}

/* End of code generation (gp_max.c) */
