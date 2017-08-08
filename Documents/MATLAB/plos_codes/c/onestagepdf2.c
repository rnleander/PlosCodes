/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * onestagepdf2.c
 *
 * Code generation for function 'onestagepdf2'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "onestagepdf2.h"
#include "emgpdf.h"
#include "IMT_analysis_April2017_emxutil.h"
#include "power.h"
#include "rdivide.h"
#include "exp.h"
#include "IMT_analysis_April2017_rtwutil.h"
#include "gsl/gsl_multimin.h"
#include "float.h"
#include "math.h"

#define _VERBOSE


/* Function Definitions */

/*
 * function Y=onestagepdf2(t,mu,s)
 */
void b_onestagepdf2(const double t[2201], double mu, double s, double Y[2201])
{
  int k;
  double a;
  double y[2201];
  double b_a[2201];
  double b_y;

  /*  find the value of the inverse gaussian with parameters mu and s */
  /*  at each point in t and return a list Y with the value cooresponding */
  /*  to the points in t */
  /*  reflect the objective fuction about the axis to allow unconstrained */
  /*  optimization */
  /* 'onestagepdf2:8' mu=abs(mu); */
  mu = fabs(mu);

  /* 'onestagepdf2:9' s=abs(s); */
  s = fabs(s);

  /* 'onestagepdf2:11' Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t))); */
  for (k = 0; k < 2201; k++) {
    y[k] = rt_powd_snf(6.2831853071795862 * rt_powd_snf(t[k], 3.0), 0.5);
    b_a[k] = mu * t[k] - 1.0;
  }

  a = 2.0 * (s * s);

  /* The pdf may return values that are zero to witin machine error */
  /* these values are also replaced by realmin */
  /* 'onestagepdf2:15' Y=max(Y, realmin); */
  for (k = 0; k < 2201; k++) {
    b_y = 1.0 / (s * y[k]) * exp(-rt_powd_snf(b_a[k], 2.0) / (a * t[k]));
    if (!(b_y > 2.2250738585072014E-308)) {
      b_y = 2.2250738585072014E-308;
    }

    Y[k] = b_y;
  }

#ifdef _VERBOSE
  printf("b_onestagepdf2 = \n");
  for (int i = 0; i < 2201; i++) {
	  if (i % 8 == 0)
		  printf("\n");
	  printf("%f ", Y[i]);
  }
  printf("\n\n");
#endif
}

/*
 * function Y=onestagepdf2(t,mu,s)
 */
double c_onestagepdf2(double t, double mu, double s)
{
  double Y;
  double a;

  /*  find the value of the inverse gaussian with parameters mu and s */
  /*  at each point in t and return a list Y with the value cooresponding */
  /*  to the points in t */
  /*  reflect the objective fuction about the axis to allow unconstrained */
  /*  optimization */
  /* 'onestagepdf2:8' mu=abs(mu); */
  mu = fabs(mu);

  /* 'onestagepdf2:9' s=abs(s); */
  s = fabs(s);

  /* 'onestagepdf2:11' Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t))); */
  a = mu * t - 1.0;
  Y = 1.0 / (s * sqrt(6.2831853071795862 * rt_powd_snf(t, 3.0))) * exp(-(a * a) /
    (2.0 * (s * s) * t));

  /* The pdf may return values that are zero to witin machine error */
  /* these values are also replaced by realmin */
  /* 'onestagepdf2:15' Y=max(Y, realmin); */
  if (!(Y > 2.2250738585072014E-308)) {
    Y = 2.2250738585072014E-308;
  }

  return Y;
}

/*
 * function Y=onestagepdf2(t,mu,s)
 */
void d_onestagepdf2(const emxArray_real_T *t, double mu, double s,
                    emxArray_real_T *Y)
{
  emxArray_real_T *y;
  int k;
  int nx;
  emxArray_real_T *b_mu;
  double a;
  int unnamed_idx_1;
  emxInit_real_T(&y, 2);

  /*  find the value of the inverse gaussian with parameters mu and s */
  /*  at each point in t and return a list Y with the value cooresponding */
  /*  to the points in t */
  /*  reflect the objective fuction about the axis to allow unconstrained */
  /*  optimization */
  /* 'onestagepdf2:8' mu=abs(mu); */
  mu = fabs(mu);

  /* 'onestagepdf2:9' s=abs(s); */
  s = fabs(s);

  /* 'onestagepdf2:11' Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t))); */
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = t->size[1];
  emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  for (k = 0; k + 1 <= t->size[1]; k++) {
    y->data[k] = rt_powd_snf(t->data[k], 3.0);
  }

  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  nx = y->size[0];
  k = y->size[1];
  nx *= k;
  for (k = 0; k < nx; k++) {
    y->data[k] *= 6.2831853071795862;
  }

  k = Y->size[0] * Y->size[1];
  Y->size[0] = 1;
  Y->size[1] = y->size[1];
  emxEnsureCapacity((emxArray__common *)Y, k, sizeof(double));
  for (k = 0; k + 1 <= y->size[1]; k++) {
    Y->data[k] = sqrt(y->data[k]);
  }

  emxInit_real_T(&b_mu, 2);
  a = 2.0 * (s * s);
  k = b_mu->size[0] * b_mu->size[1];
  b_mu->size[0] = 1;
  b_mu->size[1] = t->size[1];
  emxEnsureCapacity((emxArray__common *)b_mu, k, sizeof(double));
  nx = t->size[0] * t->size[1];
  for (k = 0; k < nx; k++) {
    b_mu->data[k] = mu * t->data[k] - 1.0;
  }

  e_power(b_mu, y);
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  nx = y->size[0];
  k = y->size[1];
  nx *= k;
  emxFree_real_T(&b_mu);
  for (k = 0; k < nx; k++) {
    y->data[k] = -y->data[k];
  }

  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  nx = y->size[0];
  k = y->size[1];
  nx *= k;
  for (k = 0; k < nx; k++) {
    y->data[k] /= a * t->data[k];
  }

  nx = y->size[1];
  for (k = 0; k + 1 <= nx; k++) {
    y->data[k] = exp(y->data[k]);
  }

  k = Y->size[0] * Y->size[1];
  Y->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)Y, k, sizeof(double));
  nx = Y->size[0];
  k = Y->size[1];
  nx *= k;
  for (k = 0; k < nx; k++) {
    Y->data[k] = 1.0 / (s * Y->data[k]) * y->data[k];
  }

  /* The pdf may return values that are zero to witin machine error */
  /* these values are also replaced by realmin */
  /* 'onestagepdf2:15' Y=max(Y, realmin); */
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = Y->size[1];
  emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  nx = Y->size[0] * Y->size[1];
  for (k = 0; k < nx; k++) {
    y->data[k] = Y->data[k];
  }

  unnamed_idx_1 = Y->size[1];
  nx = Y->size[1];
  k = Y->size[0] * Y->size[1];
  Y->size[0] = 1;
  Y->size[1] = nx;
  emxEnsureCapacity((emxArray__common *)Y, k, sizeof(double));
  for (k = 0; k + 1 <= unnamed_idx_1; k++) {
    a = y->data[k];
    if (!(a > 2.2250738585072014E-308)) {
      a = 2.2250738585072014E-308;
    }

    Y->data[k] = a;
  }

  emxFree_real_T(&y);




}

/*
 * function Y=onestagepdf2(t,mu,s)
 */
void e_onestagepdf2(const double t[22001], double mu, double s, double Y[22001])
{
  int k;
  double a;
  static double y[22001];
  double b_a[22001];
  double b_y;

  /*  find the value of the inverse gaussian with parameters mu and s */
  /*  at each point in t and return a list Y with the value cooresponding */
  /*  to the points in t */
  /*  reflect the objective fuction about the axis to allow unconstrained */
  /*  optimization */
  /* 'onestagepdf2:8' mu=abs(mu); */
  mu = fabs(mu);

  /* 'onestagepdf2:9' s=abs(s); */
  s = fabs(s);

  /* 'onestagepdf2:11' Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t))); */
  for (k = 0; k < 22001; k++) {
    y[k] = rt_powd_snf(6.2831853071795862 * rt_powd_snf(t[k], 3.0), 0.5);
    b_a[k] = mu * t[k] - 1.0;
  }

  a = 2.0 * (s * s);

  /* The pdf may return values that are zero to witin machine error */
  /* these values are also replaced by realmin */
  /* 'onestagepdf2:15' Y=max(Y, realmin); */
  for (k = 0; k < 22001; k++) {
    b_y = 1.0 / (s * y[k]) * exp(-rt_powd_snf(b_a[k], 2.0) / (a * t[k]));
    if (!(b_y > 2.2250738585072014E-308)) {
      b_y = 2.2250738585072014E-308;
    }

    Y[k] = b_y;
  }
}

/*
 * function Y=onestagepdf2(t,mu,s)
 */
void f_onestagepdf2(const double t[221], double mu, double s, double Y[221])
{
  int k;
  double a;
  double y[221];
  double b_a[221];
  double b_y;

  /*  find the value of the inverse gaussian with parameters mu and s */
  /*  at each point in t and return a list Y with the value cooresponding */
  /*  to the points in t */
  /*  reflect the objective fuction about the axis to allow unconstrained */
  /*  optimization */
  /* 'onestagepdf2:8' mu=abs(mu); */
  mu = fabs(mu);

  /* 'onestagepdf2:9' s=abs(s); */
  s = fabs(s);

  /* 'onestagepdf2:11' Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t))); */
  for (k = 0; k < 221; k++) {
    y[k] = rt_powd_snf(6.2831853071795862 * rt_powd_snf(t[k], 3.0), 0.5);
    b_a[k] = mu * t[k] - 1.0;
  }

  a = 2.0 * (s * s);

  /* The pdf may return values that are zero to witin machine error */
  /* these values are also replaced by realmin */
  /* 'onestagepdf2:15' Y=max(Y, realmin); */
  for (k = 0; k < 221; k++) {
    b_y = 1.0 / (s * y[k]) * exp(-rt_powd_snf(b_a[k], 2.0) / (a * t[k]));
    if (!(b_y > 2.2250738585072014E-308)) {
      b_y = 2.2250738585072014E-308;
    }

    Y[k] = b_y;
  }
}

/*
 * function Y=onestagepdf2(t,mu,s)
 */
void g_onestagepdf2(const double t[221], double mu, double s, double Y[221])
{
  int k;
  double a;
  double y[221];
  double b_a[221];
  double b_y;

  /*  find the value of the inverse gaussian with parameters mu and s */
  /*  at each point in t and return a list Y with the value cooresponding */
  /*  to the points in t */
  /*  reflect the objective fuction about the axis to allow unconstrained */
  /*  optimization */
  /* 'onestagepdf2:8' mu=abs(mu); */
  mu = fabs(mu);

  /* 'onestagepdf2:9' s=abs(s); */
  s = fabs(s);

  /* 'onestagepdf2:11' Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t))); */
  for (k = 0; k < 221; k++) {
    y[k] = rt_powd_snf(6.2831853071795862 * rt_powd_snf(t[k], 3.0), 0.5);
    b_a[k] = mu * t[k] - 1.0;
  }

  a = 2.0 * (s * s);

  /* The pdf may return values that are zero to witin machine error */
  /* these values are also replaced by realmin */
  /* 'onestagepdf2:15' Y=max(Y, realmin); */
  for (k = 0; k < 221; k++) {
    b_y = 1.0 / (s * y[k]) * exp(-rt_powd_snf(b_a[k], 2.0) / (a * t[k]));
    if (!(b_y > 2.2250738585072014E-308)) {
      b_y = 2.2250738585072014E-308;
    }

    Y[k] = b_y;
  }
}

/*
 * function Y=onestagepdf2(t,mu,s)
 */
void h_onestagepdf2(const emxArray_real_T *t, double mu, double s,
                    emxArray_real_T *Y)
{
  emxArray_real_T *b;
  emxArray_real_T *r6;
  int k;
  int loop_ub;
  emxArray_real_T *b_mu;
  double a;
  emxArray_real_T *r7;
  emxArray_real_T *r8;
  emxArray_real_T *b_a;
  emxArray_real_T *b_s;
  unsigned int Y_idx_0;
  unsigned int b_Y;
  emxInit_real_T1(&b, 1);
  emxInit_real_T1(&r6, 1);

  /*  find the value of the inverse gaussian with parameters mu and s */
  /*  at each point in t and return a list Y with the value cooresponding */
  /*  to the points in t */
  /*  reflect the objective fuction about the axis to allow unconstrained */
  /*  optimization */
  /* 'onestagepdf2:8' mu=abs(mu); */
  mu = fabs(mu);

  /* 'onestagepdf2:9' s=abs(s); */
  s = fabs(s);

  /* 'onestagepdf2:11' Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t))); */
  k_power(t, b);
  k = r6->size[0];
  r6->size[0] = b->size[0];
  emxEnsureCapacity((emxArray__common *)r6, k, sizeof(double));
  loop_ub = b->size[0];
  for (k = 0; k < loop_ub; k++) {
    r6->data[k] = 6.2831853071795862 * b->data[k];
  }

  emxInit_real_T1(&b_mu, 1);
  l_power(r6, b);
  a = 2.0 * (s * s);
  k = b_mu->size[0];
  b_mu->size[0] = t->size[0];
  emxEnsureCapacity((emxArray__common *)b_mu, k, sizeof(double));
  loop_ub = t->size[0];
  emxFree_real_T(&r6);
  for (k = 0; k < loop_ub; k++) {
    b_mu->data[k] = mu * t->data[k] - 1.0;
  }

  emxInit_real_T1(&r7, 1);
  emxInit_real_T1(&r8, 1);
  m_power(b_mu, r7);
  k = r8->size[0];
  r8->size[0] = r7->size[0];
  emxEnsureCapacity((emxArray__common *)r8, k, sizeof(double));
  loop_ub = r7->size[0];
  emxFree_real_T(&b_mu);
  for (k = 0; k < loop_ub; k++) {
    r8->data[k] = -r7->data[k];
  }

  emxInit_real_T1(&b_a, 1);
  k = b_a->size[0];
  b_a->size[0] = t->size[0];
  emxEnsureCapacity((emxArray__common *)b_a, k, sizeof(double));
  loop_ub = t->size[0];
  for (k = 0; k < loop_ub; k++) {
    b_a->data[k] = a * t->data[k];
  }

  emxInit_real_T1(&b_s, 1);
  c_rdivide(r8, b_a, r7);
  b_exp(r7);
  k = b_s->size[0];
  b_s->size[0] = b->size[0];
  emxEnsureCapacity((emxArray__common *)b_s, k, sizeof(double));
  loop_ub = b->size[0];
  emxFree_real_T(&b_a);
  emxFree_real_T(&r8);
  for (k = 0; k < loop_ub; k++) {
    b_s->data[k] = s * b->data[k];
  }

  b_rdivide(b_s, Y);
  k = Y->size[0];
  emxEnsureCapacity((emxArray__common *)Y, k, sizeof(double));
  loop_ub = Y->size[0];
  emxFree_real_T(&b_s);
  for (k = 0; k < loop_ub; k++) {
    Y->data[k] *= r7->data[k];
  }

  emxFree_real_T(&r7);

  /* The pdf may return values that are zero to witin machine error */
  /* these values are also replaced by realmin */
  /* 'onestagepdf2:15' Y=max(Y, realmin); */
  k = b->size[0];
  b->size[0] = Y->size[0];
  emxEnsureCapacity((emxArray__common *)b, k, sizeof(double));
  loop_ub = Y->size[0];
  for (k = 0; k < loop_ub; k++) {
    b->data[k] = Y->data[k];
  }

  Y_idx_0 = (unsigned int)Y->size[0];
  b_Y = (unsigned int)Y->size[0];
  k = Y->size[0];
  Y->size[0] = (int)b_Y;
  emxEnsureCapacity((emxArray__common *)Y, k, sizeof(double));
  for (k = 0; k + 1 <= (int)Y_idx_0; k++) {
    a = b->data[k];
    if (!(a > 2.2250738585072014E-308)) {
      a = 2.2250738585072014E-308;
    }

    Y->data[k] = a;
  }

  emxFree_real_T(&b);
}

/*
 * function Y=onestagepdf2(t,mu,s)
 */
void i_onestagepdf2(const double t[2201], double mu, double s, double Y[2201])
{
  int k;
  double a;
  double y[2201];
  double b_a[2201];
  double b_y;

  /*  find the value of the inverse gaussian with parameters mu and s */
  /*  at each point in t and return a list Y with the value cooresponding */
  /*  to the points in t */
  /*  reflect the objective fuction about the axis to allow unconstrained */
  /*  optimization */
  /* 'onestagepdf2:8' mu=abs(mu); */
  mu = fabs(mu);

  /* 'onestagepdf2:9' s=abs(s); */
  s = fabs(s);

  /* 'onestagepdf2:11' Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t))); */
  for (k = 0; k < 2201; k++) {
    y[k] = rt_powd_snf(6.2831853071795862 * rt_powd_snf(t[k], 3.0), 0.5);
    b_a[k] = mu * t[k] - 1.0;
  }

  a = 2.0 * (s * s);

  /* The pdf may return values that are zero to witin machine error */
  /* these values are also replaced by realmin */
  /* 'onestagepdf2:15' Y=max(Y, realmin); */
  for (k = 0; k < 2201; k++) {
    b_y = 1.0 / (s * y[k]) * exp(-rt_powd_snf(b_a[k], 2.0) / (a * t[k]));
    if (!(b_y > 2.2250738585072014E-308)) {
      b_y = 2.2250738585072014E-308;
    }

    Y[k] = b_y;
  }
}

/*
 * function Y=onestagepdf2(t,mu,s)
 */
void j_onestagepdf2(const double t[22001], double mu, double s, double Y[22001])
{
  int k;
  double a;
  static double y[22001];
  double b_a[22001];
  double b_y;

  /*  find the value of the inverse gaussian with parameters mu and s */
  /*  at each point in t and return a list Y with the value cooresponding */
  /*  to the points in t */
  /*  reflect the objective fuction about the axis to allow unconstrained */
  /*  optimization */
  /* 'onestagepdf2:8' mu=abs(mu); */
  mu = fabs(mu);

  /* 'onestagepdf2:9' s=abs(s); */
  s = fabs(s);

  /* 'onestagepdf2:11' Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t))); */
  for (k = 0; k < 22001; k++) {
    y[k] = rt_powd_snf(6.2831853071795862 * rt_powd_snf(t[k], 3.0), 0.5);
    b_a[k] = mu * t[k] - 1.0;
  }

  a = 2.0 * (s * s);

  /* The pdf may return values that are zero to witin machine error */
  /* these values are also replaced by realmin */
  /* 'onestagepdf2:15' Y=max(Y, realmin); */
  for (k = 0; k < 22001; k++) {
    b_y = 1.0 / (s * y[k]) * exp(-rt_powd_snf(b_a[k], 2.0) / (a * t[k]));
    if (!(b_y > 2.2250738585072014E-308)) {
      b_y = 2.2250738585072014E-308;
    }

    Y[k] = b_y;
  }
}












double wald_loglikelihood(const gsl_vector *v, void *params)
{
	double *data = (double *)params;

	double m = gsl_vector_get(v, 0);
	double s = gsl_vector_get(v, 1);

	double penalty = 0;
	if (m < 0 || s < 0)
		penalty = 1000;

	m = fabs(m);
	s = fabs(s);

	double Y[266];

	waldpdf(data, m, s, Y, 266);
	//onestagepdf2(data, m, s, Y);

	double loglikelihood = 0;
	for (int i = 0; i < 266; i++) {
		loglikelihood += log(Y[i]);
	}

	return penalty - loglikelihood;
}

void waldpdf(const double X[], double mu, double s, double Y[], int size_XY)
{
	// Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t)));
	// https://en.wikipedia.org/wiki/Inverse_Gaussian_distribution

	double a, b;
	for (int i = 0; i < size_XY; i++) {
		a = 1.0 / (s*pow(2 * M_PI  * pow(X[i], 3.0), 0.5));
		b = (pow(mu*X[i] - 1, 2)) / (2.0 * s * s * X[i]);
		Y[i] = a*exp(-b);
		if (Y[i] == 0)
			Y[i] = 2.2250738585072014E-308;
			//Y[i] = -DBL_MAX;

		if (isnan(Y[i]))
			Y[i] = 2.2250738585072014E-308;
	}

#ifdef _VERBOSE
	printf("waldpdf = \n");
	for (int i = 0; i < size_XY; i++) {
		if (i % 8 == 0)
			printf("\n");
		printf("%f ", Y[i]);
	}
	printf("\n\n");
#endif
	
}



/*
 * function Y=onestagepdf2(t,mu,s)
 */
void onestagepdf2(const double t[266], double mu, double s, double Y[266])
{
  int k;
  double a;
  double y[266];
  double b_a[266];
  double b_y;

  /*  find the value of the inverse gaussian with parameters mu and s */
  /*  at each point in t and return a list Y with the value cooresponding */
  /*  to the points in t */
  /*  reflect the objective fuction about the axis to allow unconstrained */
  /*  optimization */
  /* 'onestagepdf2:8' mu=abs(mu); */
  mu = fabs(mu);

  /* 'onestagepdf2:9' s=abs(s); */
  s = fabs(s);

  /* 'onestagepdf2:11' Y=(1./(s*(2*pi*t.^3).^(.5))).*exp(-((mu*t-1).^2)./(2*s^2*(t))); */
  for (k = 0; k < 266; k++) {
    y[k] = rt_powd_snf(6.2831853071795862 * rt_powd_snf(t[k], 3.0), 0.5);
    b_a[k] = mu * t[k] - 1.0;
  }

  a = 2.0 * (s * s);

  /* The pdf may return values that are zero to witin machine error */
  /* these values are also replaced by realmin */
  /* 'onestagepdf2:15' Y=max(Y, realmin); */
  for (k = 0; k < 266; k++) {
    b_y = 1.0 / (s * y[k]) * exp(-rt_powd_snf(b_a[k], 2.0) / (a * t[k]));
    if (!(b_y > 2.2250738585072014E-308)) {
      b_y = 2.2250738585072014E-308;
    }

    Y[k] = b_y;
  }
}

/* End of code generation (onestagepdf2.c) */
