/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * onestagepdf_prime.c
 *
 * Code generation for function 'onestagepdf_prime'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "onestagepdf_prime.h"
#include "emgpdf.h"
#include "gp_max.h"
#include "rdivide.h"
#include "IMT_analysis_April2017_rtwutil.h"
#include <math.h>

/* Function Declarations */
static double rt_atan2d_snf(double u0, double u1);

/* Function Definitions */
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  int b_u0;
  int b_u1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2(b_u0, b_u1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}


void onestagepdf_prime_fixed(const double t[], int t_size, double m, double s, double Y[])
{
	double a;
	int k;
	double y;
	double b_y;

	/*  find the value of the first derivitive of the inverse gaussian */
	/*  distribution with parameters m and s at each point in t */
	/*  and return a list of these values corresponding to each point t */
	/* 'onestagepdf_prime:6' Y=exp((-(1-m*t).^2)./(2*s^2*t)); */
	a = 2.0 * (s * s);
	for (k = 0; k < t_size; k++) {
		y = 1.0 - m * t[k];
		Y[k] = exp(-(y * y) / (a * t[k]));
	}

	/* 'onestagepdf_prime:7' Y=(-3/2*t.^(-5/2)+(t.^(-2)/(2*s^2)-m^2/(2*s^2)).*t.^(-3/2))*(1/(s*(2*pi)^.5)).*Y; */
	a = 2.0 * (s * s);
	y = m * m / (2.0 * (s * s));
	b_y = 1.0 / (s * 2.5066282746310002);
	for (k = 0; k < t_size; k++) {
		Y[k] *= (-1.5 * rt_powd_snf(t[k], -2.5) + (rt_powd_snf(t[k], -2.0) / a - y) *
			rt_powd_snf(t[k], -1.5)) * b_y;
	}
}


/*
 * function Y=onestagepdf_prime(t,m,s)
 */
void onestagepdf_prime(const creal_T t_data[], const int t_size[1], double m,
  double s, creal_T Y_data[], int Y_size[1])
{
  int loop_ub;
  int i1;
  creal_T a_data[4];
  double a;
  creal_T y_data[4];
  int y_size[1];
  creal_T b_y_data[4];
  double r;
  double y_data_im;
  double B;
  double y;
  double b_y;
  double re;

  /*  find the value of the first derivitive of the inverse gaussian */
  /*  distribution with parameters m and s at each point in t */
  /*  and return a list of these values corresponding to each point t */
  /* 'onestagepdf_prime:6' Y=exp((-(1-m*t).^2)./(2*s^2*t)); */
  loop_ub = t_size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    a_data[i1].re = 1.0 - m * t_data[i1].re;
    a_data[i1].im = 0.0 - m * t_data[i1].im;
  }

  for (loop_ub = 0; loop_ub + 1 <= t_size[0]; loop_ub++) {
    y_data[loop_ub].re = a_data[loop_ub].re * a_data[loop_ub].re -
      a_data[loop_ub].im * a_data[loop_ub].im;
    y_data[loop_ub].im = a_data[loop_ub].re * a_data[loop_ub].im +
      a_data[loop_ub].im * a_data[loop_ub].re;
  }

  a = 2.0 * (s * s);
  y_size[0] = (signed char)t_size[0];
  loop_ub = (signed char)t_size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    b_y_data[i1].re = -y_data[i1].re;
    b_y_data[i1].im = -y_data[i1].im;
  }

  loop_ub = t_size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    a_data[i1].re = a * t_data[i1].re;
    a_data[i1].im = a * t_data[i1].im;
  }

  rdivide(b_y_data, y_size, a_data, Y_data, Y_size);
  for (loop_ub = 0; loop_ub + 1 <= Y_size[0]; loop_ub++) {
    if (Y_data[loop_ub].im == 0.0) {
      a = exp(Y_data[loop_ub].re);
      y_data_im = 0.0;
    } else if (rtIsInf(Y_data[loop_ub].im) && rtIsInf(Y_data[loop_ub].re) &&
               (Y_data[loop_ub].re < 0.0)) {
      a = 0.0;
      y_data_im = 0.0;
    } else {
      r = exp(Y_data[loop_ub].re / 2.0);
      a = r * (r * cos(Y_data[loop_ub].im));
      y_data_im = r * (r * sin(Y_data[loop_ub].im));
    }

    Y_data[loop_ub].re = a;
    Y_data[loop_ub].im = y_data_im;
  }

  /* 'onestagepdf_prime:7' Y=(-3/2*t.^(-5/2)+(t.^(-2)/(2*s^2)-m^2/(2*s^2)).*t.^(-3/2))*(1/(s*(2*pi)^.5)).*Y; */
  for (loop_ub = 0; loop_ub + 1 <= t_size[0]; loop_ub++) {
    a = t_data[loop_ub].re;
    y_data_im = t_data[loop_ub].im;
    if ((t_data[loop_ub].im == 0.0) && (t_data[loop_ub].re >= 0.0)) {
      a = rt_powd_snf(t_data[loop_ub].re, -2.0);
      y_data_im = 0.0;
    } else if (t_data[loop_ub].re == 0.0) {
      a = -rt_powd_snf(t_data[loop_ub].im, -2.0);
      y_data_im = 0.0;
    } else {
      if ((t_data[loop_ub].im == 0.0) && rtIsNaN(t_data[loop_ub].re)) {
      } else if ((fabs(t_data[loop_ub].re) > 8.9884656743115785E+307) || (fabs
                  (t_data[loop_ub].im) > 8.9884656743115785E+307)) {
        a = log(rt_hypotd_snf(t_data[loop_ub].re / 2.0, t_data[loop_ub].im / 2.0))
          + 0.69314718055994529;
        y_data_im = rt_atan2d_snf(t_data[loop_ub].im, t_data[loop_ub].re);
      } else {
        a = log(rt_hypotd_snf(t_data[loop_ub].re, t_data[loop_ub].im));
        y_data_im = rt_atan2d_snf(t_data[loop_ub].im, t_data[loop_ub].re);
      }

      a *= -2.0;
      y_data_im *= -2.0;
      if (y_data_im == 0.0) {
        a = exp(a);
        y_data_im = 0.0;
      } else if (rtIsInf(y_data_im) && rtIsInf(a) && (a < 0.0)) {
        a = 0.0;
        y_data_im = 0.0;
      } else {
        r = exp(a / 2.0);
        a = r * (r * cos(y_data_im));
        y_data_im = r * (r * sin(y_data_im));
      }
    }

    y_data[loop_ub].re = a;
    y_data[loop_ub].im = y_data_im;
  }

  B = 2.0 * (s * s);
  for (loop_ub = 0; loop_ub + 1 <= t_size[0]; loop_ub++) {
    a = t_data[loop_ub].re;
    y_data_im = t_data[loop_ub].im;
    if ((t_data[loop_ub].im == 0.0) && (t_data[loop_ub].re >= 0.0)) {
      a = rt_powd_snf(t_data[loop_ub].re, -1.5);
      y_data_im = 0.0;
    } else {
      if ((t_data[loop_ub].im == 0.0) && rtIsNaN(t_data[loop_ub].re)) {
      } else if ((fabs(t_data[loop_ub].re) > 8.9884656743115785E+307) || (fabs
                  (t_data[loop_ub].im) > 8.9884656743115785E+307)) {
        a = log(rt_hypotd_snf(t_data[loop_ub].re / 2.0, t_data[loop_ub].im / 2.0))
          + 0.69314718055994529;
        y_data_im = rt_atan2d_snf(t_data[loop_ub].im, t_data[loop_ub].re);
      } else {
        a = log(rt_hypotd_snf(t_data[loop_ub].re, t_data[loop_ub].im));
        y_data_im = rt_atan2d_snf(t_data[loop_ub].im, t_data[loop_ub].re);
      }

      a *= -1.5;
      y_data_im *= -1.5;
      if (y_data_im == 0.0) {
        a = exp(a);
        y_data_im = 0.0;
      } else if (rtIsInf(y_data_im) && rtIsInf(a) && (a < 0.0)) {
        a = 0.0;
        y_data_im = 0.0;
      } else {
        r = exp(a / 2.0);
        a = r * (r * cos(y_data_im));
        y_data_im = r * (r * sin(y_data_im));
      }
    }

    a_data[loop_ub].re = a;
    a_data[loop_ub].im = y_data_im;
  }

  for (loop_ub = 0; loop_ub + 1 <= t_size[0]; loop_ub++) {
    a = t_data[loop_ub].re;
    y_data_im = t_data[loop_ub].im;
    if ((t_data[loop_ub].im == 0.0) && (t_data[loop_ub].re >= 0.0)) {
      a = rt_powd_snf(t_data[loop_ub].re, -2.5);
      y_data_im = 0.0;
    } else {
      if ((t_data[loop_ub].im == 0.0) && rtIsNaN(t_data[loop_ub].re)) {
      } else if ((fabs(t_data[loop_ub].re) > 8.9884656743115785E+307) || (fabs
                  (t_data[loop_ub].im) > 8.9884656743115785E+307)) {
        a = log(rt_hypotd_snf(t_data[loop_ub].re / 2.0, t_data[loop_ub].im / 2.0))
          + 0.69314718055994529;
        y_data_im = rt_atan2d_snf(t_data[loop_ub].im, t_data[loop_ub].re);
      } else {
        a = log(rt_hypotd_snf(t_data[loop_ub].re, t_data[loop_ub].im));
        y_data_im = rt_atan2d_snf(t_data[loop_ub].im, t_data[loop_ub].re);
      }

      a *= -2.5;
      y_data_im *= -2.5;
      if (y_data_im == 0.0) {
        a = exp(a);
        y_data_im = 0.0;
      } else if (rtIsInf(y_data_im) && rtIsInf(a) && (a < 0.0)) {
        a = 0.0;
        y_data_im = 0.0;
      } else {
        r = exp(a / 2.0);
        a = r * (r * cos(y_data_im));
        y_data_im = r * (r * sin(y_data_im));
      }
    }

    b_y_data[loop_ub].re = a;
    b_y_data[loop_ub].im = y_data_im;
  }

  y = m * m / (2.0 * (s * s));
  b_y = 1.0 / (s * 2.5066282746310002);
  Y_size[0] = (signed char)t_size[0];
  loop_ub = (signed char)t_size[0];
  for (i1 = 0; i1 < loop_ub; i1++) {
    if (y_data[i1].im == 0.0) {
      a = y_data[i1].re / B;
      y_data_im = 0.0;
    } else if (y_data[i1].re == 0.0) {
      a = 0.0;
      y_data_im = y_data[i1].im / B;
    } else {
      a = y_data[i1].re / B;
      y_data_im = y_data[i1].im / B;
    }

    a -= y;
    re = b_y * (-1.5 * b_y_data[i1].re + (a * a_data[i1].re - y_data_im *
      a_data[i1].im));
    a = b_y * (-1.5 * b_y_data[i1].im + (a * a_data[i1].im + y_data_im *
                a_data[i1].re));
    r = Y_data[i1].re;
    Y_data[i1].re = re * Y_data[i1].re - a * Y_data[i1].im;
    Y_data[i1].im = re * Y_data[i1].im + a * r;
  }
}

/* End of code generation (onestagepdf_prime.c) */
