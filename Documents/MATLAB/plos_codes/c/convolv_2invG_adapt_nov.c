/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * convolv_2invG_adapt_nov.c
 *
 * Code generation for function 'convolv_2invG_adapt_nov'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "convolv_2invG_adapt_nov.h"
#include "emgpdf.h"
#include "IMT_analysis_April2017_emxutil.h"
#include "onestagepdf_lag.h"
#include "approxconvolv.h"
#include "onestagepdf2.h"
#include "tailmass.h"
#include "sort1.h"
#include "power.h"
#include "abs.h"
#include "IMT_analysis_April2017_rtwutil.h"
#include "gsl/gsl_multimin.h"
#include "math.h"


//#define _CONV2INVG
#define _CONV2WALD
//#define _VERBOSE
//#define _PARALLEL_PDF


/* Function Definitions */

double convolv_2invG_adapt_nov_loglikelihood(const gsl_vector *v, void *params)
{
	double *data = (double *)params;

	double m1 = gsl_vector_get(v, 0);
	double s1 = gsl_vector_get(v, 1);
	double m2 = gsl_vector_get(v, 2);
	double s2 = gsl_vector_get(v, 3);

	double penalty = 0;
	if (m1 < 0 || s1 < 0 || m2 < 0 || s2 < 0)
		penalty = 1000;

	m1 = fabs(m1);
	s1 = fabs(s1);
	m2 = fabs(m2);
	s2 = fabs(s2);

	double Y[266];

	double Y_WALD[266];

	double Y_INVG[266];

#ifdef _CONV2INVG
#ifdef _VERBOSE
	printf("starting convolv_2invG_adapt_nov\n");
#endif
	convolv_2invG_adapt_nov(m1, s1, m2, s2, Y_INVG);
#ifdef _VERBOSE
	printf("convolv_2invG_adapt_nov = \n");
	for (int i = 0; i < 266; i++) {
		if (i % 8 == 0)
			printf("\n");
		printf("%f ", Y_INVG[i]);
	}
	printf("\n\n");
#endif
#endif

#ifdef _CONV2WALD
#ifdef _VERBOSE
	printf("starting conv2waldpdf\n");
#endif
	conv2waldpdf(data, m1, s1, m2, s2, Y_WALD, 0.01, 1, 266);
#ifdef _VERBOSE
	printf("conv2waldpdf = \n");
	for (int i = 0; i < 266; i++) {
		if (i % 8 == 0)
			printf("\n");
		printf("%f ", Y_WALD[i]);
	}
	printf("\n\n");
#endif
#endif

#ifdef _CONV2INVG
#ifdef _CONV2WALD
	printf("difference between convolv_2invG_adapt_nov and conv2wald = \n");
	for (int i = 0; i < 266; i++) {
		if (i % 8 == 0)
			printf("\n");
		printf("%.17f ", Y_INVG[i] - Y_WALD[i]);
	}
	printf("\n\n");
#endif
#endif

#ifdef _CONV2WALD
	for (int i = 0; i < 266; i++)
		Y[i] = Y_WALD[i];
#endif

#ifdef _CONV2INVG
	for (int i = 0; i < 266; i++)
		Y[i] = Y_INVG[i];
#endif


	double loglikelihood = 0;
	for (int i = 0; i < 266; i++) {
		loglikelihood += log(Y[i]);
		//printf("data[%d]=%.17f Y[%d]=%.17f log(%.17f)=%.17f\n", i, data[i], i, Y[i], Y[i], log(Y[i]));

	}
	//printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\nloglikelihood = %f\n\n\n", loglikelihood);

	//exit(1);

	return penalty - loglikelihood;
}

void conv2waldpdf(const double X[], double m1, double s1, double m2, double s2, double Y[], double h, int adaptiveMode, int size_XY)
{

	int flag = 0; // remember if we applied the approximation
	double eps = 0.01; // a constant that is used in determining if the Dirac approximation should be applied.

	double E;
	if (adaptiveMode)
		E = FP_INFINITE;
	else
		E = 0;

	//int N = 266; // number of points to be evaluated
	int N = size_XY;

	double m_a = m1;
	double m_b = m2;
	double s_a = s1;
	double s_b = s2;

	double m[2];
	double s[2];
	m[0] = m1; m[1] = m2;
	s[0] = s1; s[1] = s2;

	// find the variance for both sub-distributions
	double v_a = (s_a * s_a) / (m_a * m_a * m_a);
	double v_b = (s_b * s_b) / (m_b * m_b * m_b);

	// find the standard deviation for each sub-distribution
	double sd_a = pow(v_a, 0.5);
	double sd_b = pow(v_b, 0.5);

	/* reorder m and s so that the sub-distribution with the smallest
	variance comes first.  So the fist part might be approximated as a Dirac delta.*/
	if (sd_a > sd_b) {
		double tmp;
		tmp = m_a; m_a = m_b; m_b = tmp;
		tmp = s_a; s_a = s_b; s_b = tmp;
		tmp = v_a; v_a = v_b; v_b = tmp;
		tmp = sd_a; sd_a = sd_b; sd_b = tmp;
	}

	// find the largest point in t
	double maxX = 0;
	for (int i = 0; i < size_XY; i++) {
		if (X[i] > maxX)
			maxX = X[i];
	}

	// These represent our loglikelihoods
	double logP0;
	double logP1;


	/* if the first pdf is very concentrated, check to see if it can be
	approximated as a point-mass distribution */
	if (sd_a < 0.01 && checktailmass(m_a, s_a, m_b, s_b, eps, sd_a, sd_b)) {
#ifdef _VERBOSE
		printf("using dirac delta approximation\n");
#endif

		/* If there is not much probability in the tails of the convolution,
		we use a Dirac Delta for part 1.
		Flag is set to 1 to indicate this. */
		double lag = 1 / m_a;
		//onestagepdf_lag(X, m_b, s_b, lag, Y);
		waldlagpdf(X, m_b, s_b, lag, Y, size_XY);
		flag = 1;
	}
	else {

		/* produce a range of evenly spaced points to evaluate at between
		0 and Maxt with step size h. the even spacing is important when
		calculating the convolution later */
		int partitionLength = (int)(maxX / h);

		// This represents the partition
		//double * x = _mm_malloc(partitionLength * sizeof(double), 64);
		double * x = malloc(partitionLength * sizeof(double));

		// There are our two wald distributions
		//double * y = _mm_malloc(partitionLength * sizeof(double), 64);
		double * y = malloc(partitionLength * sizeof(double));


		//double * z = _mm_malloc(partitionLength * sizeof(double), 64);
		double * z = malloc(partitionLength * sizeof(double));


		// fill the partition
		/*
		for (int i = 0; i < partitionLength; i++)
		x[i] = i*h;
		*/
		double tally = 0;
		for (int i = 0; i < partitionLength; i++) {
			x[i] = tally;
			tally += h;
		}

		// evaluate the sub-distributions at each point in x

		/*
		waldpdf(x, m_a, s_a, y, partitionLength);
		waldpdf(x, m_b, s_b, z, partitionLength);
		*/
		double * w[2];

		w[0] = y;
		w[1] = z;

#ifdef _PARALLEL_PDF
#pragma omp parallel for
#endif
		for (int i = 0; i < 2; i++) {
			waldpdf(x, m[i], s[i], w[i], partitionLength);
		}

		//printf("\n\ncalling approxconv_rep from conv2waldpdf partitionLength=%d size_XY=%d h=%f\n", partitionLength, size_XY, h);
		//approxconvolv_replacement(z, y, X, x, Y, &logP0, partitionLength, 266, h);
		approxconvolv_replacement(z, y, X, x, Y, &logP0, partitionLength, size_XY, h);

		//_mm_free(x); _mm_free(y); _mm_free(z);
		free(x); free(y); free(z);

		while (E >= 0.001*fabs(logP0)) {
			h = h * 0.5; // Shrink the step size
#ifdef _VERBOSE
			printf("h=%f logP0=%f ", h, logP0);
#endif
			int partitionLength = (int)(maxX / h);
			printf("pL = %d E=%f \n", partitionLength, E);

			//double * x = _mm_malloc(partitionLength * sizeof(double), 64);
			//double * y = _mm_malloc(partitionLength * sizeof(double), 64);
			//double * z = _mm_malloc(partitionLength * sizeof(double), 64);
			double * x = malloc(partitionLength * sizeof(double));
			double * y = malloc(partitionLength * sizeof(double));
			double * z = malloc(partitionLength * sizeof(double));

			// fill the partition
			/* CAN I REPLACE THIS WITH AN ADDITIVE FILLER? x[i+1] = x[i]+h
			HOTSPIT ANALYSIS SUGGESTS THIS IS 18% OF RUNTIME */

			/*
			for (int i = 0; i < partitionLength; i++)
			x[i] = i*h;
			*/
			double tally = 0;
			for (int i = 0; i < partitionLength; i++) {
				x[i] = tally;
				tally += h;
			}


			// evaluate the sub-distributions at each point in x
			/*
			waldpdf(x, m_a, s_a, y, partitionLength);
			waldpdf(x, m_b, s_b, z, partitionLength);
			*/
			w[0] = y;
			w[1] = z;

#ifdef _PARALLEL_PDF
#pragma omp parallel for
#endif
			for (int i = 0; i < 2; i++) {
				waldpdf(x, m[i], s[i], w[i], partitionLength);
			}


			//printf("\n\ncalling approxconv_rep from conv2waldpdf partitionLength=%d size_XY=%d h=%f\n", partitionLength, size_XY, h);
			//approxconvolv_replacement(z, y, X, x, Y, &logP1, partitionLength, 266, h);
			approxconvolv_replacement(z, y, X, x, Y, &logP1, partitionLength, size_XY, h);

			//_mm_free(x); _mm_free(y); _mm_free(z);
			free(x); free(y); free(z);

			E = fabs(logP1 - logP0);

			logP0 = logP1;

			//printf("logP1=%f E=%f\n", logP1, E);

			//exit(1);
		}
	}

}

/*
 * function [P,h,flag,E]=convolv_2invG_adapt_nov(t,m1,s1,m2,s2,h)
 */
#ifdef _OLD_MATLAB_CODE
void b_convolv_2invG_adapt_nov(double m1, double s1, double m2, double s2,
  double P[266], double *h, double *flag, double *E)
{
  double x[2];
  int k;
  double m[2];
  double ndbl;
  double v[2];
  double s[2];
  double sd[2];
  int iidx[2];
  int nm1d2;
  double b_s[2];
  static const double dv66[2201] = { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
    0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
    0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32,
    0.33, 0.34, 0.35000000000000003, 0.36, 0.37, 0.38, 0.39, 0.4,
    0.41000000000000003, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47000000000000003, 0.48,
    0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57000000000000006, 0.58,
    0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68,
    0.69000000000000006, 0.70000000000000007, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76,
    0.77, 0.78, 0.79, 0.8, 0.81, 0.82000000000000006, 0.83000000000000007, 0.84,
    0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94000000000000006,
    0.95000000000000007, 0.96, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04,
    1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.1300000000000001,
    1.1400000000000001, 1.1500000000000001, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21,
    1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33, 1.34,
    1.35, 1.36, 1.37, 1.3800000000000001, 1.3900000000000001, 1.4000000000000001,
    1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53,
    1.54, 1.55, 1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.6300000000000001,
    1.6400000000000001, 1.6500000000000001, 1.6600000000000001, 1.67, 1.68, 1.69,
    1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82,
    1.83, 1.84, 1.85, 1.86, 1.87, 1.8800000000000001, 1.8900000000000001,
    1.9000000000000001, 1.9100000000000001, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97,
    1.98, 1.99, 2.0, 2.0100000000000002, 2.02, 2.0300000000000002, 2.04, 2.05,
    2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18,
    2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.2600000000000002, 2.27,
    2.2800000000000002, 2.29, 2.3000000000000003, 2.31, 2.32, 2.33, 2.34, 2.35,
    2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48,
    2.49, 2.5, 2.5100000000000002, 2.52, 2.5300000000000002, 2.54,
    2.5500000000000003, 2.56, 2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64,
    2.65, 2.66, 2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75,
    2.7600000000000002, 2.77, 2.7800000000000002, 2.79, 2.8000000000000003, 2.81,
    2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 2.89, 2.9, 2.91, 2.92, 2.93, 2.94,
    2.95, 2.96, 2.97, 2.98, 2.99, 3.0, 3.0100000000000002, 3.02,
    3.0300000000000002, 3.04, 3.0500000000000003, 3.06, 3.0700000000000003, 3.08,
    3.09, 3.1, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.2, 3.21,
    3.22, 3.23, 3.24, 3.25, 3.2600000000000002, 3.27, 3.2800000000000002, 3.29,
    3.3000000000000003, 3.31, 3.3200000000000003, 3.33, 3.34, 3.35, 3.36, 3.37,
    3.38, 3.39, 3.4, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.48, 3.49, 3.5,
    3.5100000000000002, 3.52, 3.5300000000000002, 3.54, 3.5500000000000003, 3.56,
    3.5700000000000003, 3.58, 3.59, 3.6, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66,
    3.67, 3.68, 3.69, 3.7, 3.71, 3.72, 3.73, 3.74, 3.75, 3.7600000000000002,
    3.77, 3.7800000000000002, 3.79, 3.8000000000000003, 3.81, 3.8200000000000003,
    3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95,
    3.96, 3.97, 3.98, 3.99, 4.0, 4.01, 4.0200000000000005, 4.03, 4.04, 4.05,
    4.0600000000000005, 4.07, 4.08, 4.09, 4.1, 4.11, 4.12, 4.13, 4.14, 4.15,
    4.16, 4.17, 4.18, 4.19, 4.2, 4.21, 4.22, 4.23, 4.24, 4.25, 4.26,
    4.2700000000000005, 4.28, 4.29, 4.3, 4.3100000000000005, 4.32, 4.33, 4.34,
    4.3500000000000005, 4.36, 4.37, 4.38, 4.39, 4.4, 4.41, 4.42, 4.43, 4.44,
    4.45, 4.46, 4.47, 4.48, 4.49, 4.5, 4.51, 4.5200000000000005, 4.53, 4.54,
    4.55, 4.5600000000000005, 4.57, 4.58, 4.59, 4.6000000000000005, 4.61, 4.62,
    4.63, 4.64, 4.65, 4.66, 4.67, 4.68, 4.69, 4.7, 4.71, 4.72, 4.73, 4.74, 4.75,
    4.76, 4.7700000000000005, 4.78, 4.79, 4.8, 4.8100000000000005, 4.82, 4.83,
    4.84, 4.8500000000000005, 4.86, 4.87, 4.88, 4.89, 4.9, 4.91, 4.92, 4.93,
    4.94, 4.95, 4.96, 4.97, 4.98, 4.99, 5.0, 5.01, 5.0200000000000005, 5.03,
    5.04, 5.05, 5.0600000000000005, 5.07, 5.08, 5.09, 5.1000000000000005, 5.11,
    5.12, 5.13, 5.14, 5.15, 5.16, 5.17, 5.18, 5.19, 5.2, 5.21, 5.22, 5.23, 5.24,
    5.25, 5.26, 5.2700000000000005, 5.28, 5.29, 5.3, 5.3100000000000005, 5.32,
    5.33, 5.34, 5.3500000000000005, 5.36, 5.37, 5.38, 5.39, 5.4, 5.41, 5.42,
    5.43, 5.44, 5.45, 5.46, 5.47, 5.48, 5.49, 5.5, 5.51, 5.5200000000000005,
    5.53, 5.54, 5.55, 5.5600000000000005, 5.57, 5.58, 5.59, 5.6000000000000005,
    5.61, 5.62, 5.63, 5.64, 5.65, 5.66, 5.67, 5.68, 5.69, 5.7, 5.71, 5.72, 5.73,
    5.74, 5.75, 5.76, 5.7700000000000005, 5.78, 5.79, 5.8, 5.8100000000000005,
    5.82, 5.83, 5.84, 5.8500000000000005, 5.86, 5.87, 5.88, 5.89, 5.9, 5.91,
    5.92, 5.93, 5.94, 5.95, 5.96, 5.97, 5.98, 5.99, 6.0, 6.01,
    6.0200000000000005, 6.03, 6.04, 6.05, 6.0600000000000005, 6.07, 6.08, 6.09,
    6.1000000000000005, 6.11, 6.12, 6.13, 6.1400000000000006, 6.15, 6.16, 6.17,
    6.18, 6.19, 6.2, 6.21, 6.22, 6.23, 6.24, 6.25, 6.26, 6.2700000000000005,
    6.28, 6.29, 6.3, 6.3100000000000005, 6.32, 6.33, 6.34, 6.3500000000000005,
    6.36, 6.37, 6.38, 6.3900000000000006, 6.4, 6.41, 6.42, 6.43, 6.44, 6.45,
    6.46, 6.47, 6.48, 6.49, 6.5, 6.51, 6.5200000000000005, 6.53, 6.54, 6.55,
    6.5600000000000005, 6.57, 6.58, 6.59, 6.6000000000000005, 6.61, 6.62, 6.63,
    6.6400000000000006, 6.65, 6.66, 6.67, 6.68, 6.69, 6.7, 6.71, 6.72, 6.73,
    6.74, 6.75, 6.76, 6.7700000000000005, 6.78, 6.79, 6.8, 6.8100000000000005,
    6.82, 6.83, 6.84, 6.8500000000000005, 6.86, 6.87, 6.88, 6.8900000000000006,
    6.9, 6.91, 6.92, 6.93, 6.94, 6.95, 6.96, 6.97, 6.98, 6.99, 7.0, 7.01,
    7.0200000000000005, 7.03, 7.04, 7.05, 7.0600000000000005, 7.07, 7.08, 7.09,
    7.1000000000000005, 7.11, 7.12, 7.13, 7.1400000000000006, 7.15, 7.16, 7.17,
    7.18, 7.19, 7.2, 7.21, 7.22, 7.23, 7.24, 7.25, 7.26, 7.2700000000000005,
    7.28, 7.29, 7.3, 7.3100000000000005, 7.32, 7.33, 7.34, 7.3500000000000005,
    7.36, 7.37, 7.38, 7.3900000000000006, 7.4, 7.41, 7.42, 7.43, 7.44, 7.45,
    7.46, 7.47, 7.48, 7.49, 7.5, 7.51, 7.5200000000000005, 7.53, 7.54, 7.55,
    7.5600000000000005, 7.57, 7.58, 7.59, 7.6000000000000005, 7.61, 7.62, 7.63,
    7.6400000000000006, 7.65, 7.66, 7.67, 7.68, 7.69, 7.7, 7.71, 7.72, 7.73,
    7.74, 7.75, 7.76, 7.7700000000000005, 7.78, 7.79, 7.8, 7.8100000000000005,
    7.82, 7.83, 7.84, 7.8500000000000005, 7.86, 7.87, 7.88, 7.8900000000000006,
    7.9, 7.91, 7.92, 7.9300000000000006, 7.94, 7.95, 7.96, 7.97, 7.98, 7.99, 8.0,
    8.01, 8.02, 8.03, 8.0400000000000009, 8.05, 8.06, 8.07, 8.08, 8.09, 8.1,
    8.11, 8.120000000000001, 8.13, 8.14, 8.15, 8.16, 8.17, 8.18, 8.19, 8.2, 8.21,
    8.22, 8.23, 8.24, 8.25, 8.26, 8.27, 8.28, 8.2900000000000009, 8.3, 8.31,
    8.32, 8.33, 8.34, 8.35, 8.36, 8.370000000000001, 8.38, 8.39, 8.4, 8.41, 8.42,
    8.43, 8.44, 8.45, 8.46, 8.47, 8.48, 8.49, 8.5, 8.51, 8.52, 8.53,
    8.5400000000000009, 8.55, 8.56, 8.57, 8.58, 8.59, 8.6, 8.61,
    8.620000000000001, 8.63, 8.64, 8.65, 8.66, 8.67, 8.68, 8.69,
    8.7000000000000011, 8.71, 8.72, 8.73, 8.74, 8.75, 8.76, 8.77, 8.78,
    8.7900000000000009, 8.8, 8.81, 8.82, 8.83, 8.84, 8.85, 8.86,
    8.870000000000001, 8.88, 8.89, 8.9, 8.91, 8.92, 8.93, 8.94,
    8.9500000000000011, 8.96, 8.97, 8.98, 8.99, 9.0, 9.01, 9.02, 9.03,
    9.0400000000000009, 9.05, 9.06, 9.07, 9.08, 9.09, 9.1, 9.11,
    9.120000000000001, 9.13, 9.14, 9.15, 9.16, 9.17, 9.18, 9.19,
    9.2000000000000011, 9.21, 9.22, 9.23, 9.24, 9.25, 9.26, 9.27, 9.28,
    9.2900000000000009, 9.3, 9.31, 9.32, 9.33, 9.34, 9.35, 9.36,
    9.370000000000001, 9.38, 9.39, 9.4, 9.41, 9.42, 9.43, 9.44,
    9.4500000000000011, 9.46, 9.47, 9.48, 9.49, 9.5, 9.51, 9.52, 9.53,
    9.5400000000000009, 9.55, 9.56, 9.57, 9.58, 9.59, 9.6, 9.61,
    9.620000000000001, 9.63, 9.64, 9.65, 9.66, 9.67, 9.68, 9.69,
    9.7000000000000011, 9.71, 9.72, 9.73, 9.74, 9.75, 9.76, 9.77, 9.78,
    9.7900000000000009, 9.8, 9.81, 9.82, 9.83, 9.84, 9.85, 9.86,
    9.870000000000001, 9.88, 9.89, 9.9, 9.91, 9.92, 9.93, 9.94,
    9.9500000000000011, 9.96, 9.97, 9.98, 9.99, 10.0, 10.01, 10.02, 10.03,
    10.040000000000001, 10.05, 10.06, 10.07, 10.08, 10.09, 10.1, 10.11,
    10.120000000000001, 10.13, 10.14, 10.15, 10.16, 10.17, 10.18, 10.19,
    10.200000000000001, 10.21, 10.22, 10.23, 10.24, 10.25, 10.26, 10.27, 10.28,
    10.290000000000001, 10.3, 10.31, 10.32, 10.33, 10.34, 10.35, 10.36,
    10.370000000000001, 10.38, 10.39, 10.4, 10.41, 10.42, 10.43, 10.44,
    10.450000000000001, 10.46, 10.47, 10.48, 10.49, 10.5, 10.51, 10.52, 10.53,
    10.540000000000001, 10.55, 10.56, 10.57, 10.58, 10.59, 10.6, 10.61,
    10.620000000000001, 10.63, 10.64, 10.65, 10.66, 10.67, 10.68, 10.69,
    10.700000000000001, 10.71, 10.72, 10.73, 10.74, 10.75, 10.76, 10.77, 10.78,
    10.790000000000001, 10.8, 10.81, 10.82, 10.83, 10.84, 10.85, 10.86,
    10.870000000000001, 10.88, 10.89, 10.9, 10.91, 10.92, 10.93, 10.94,
    10.950000000000001, 10.96, 10.97, 10.98, 10.99, 11.0, 11.01, 11.02, 11.03,
    11.04, 11.049999999999999, 11.06, 11.07, 11.08, 11.09, 11.1, 11.11, 11.12,
    11.129999999999999, 11.14, 11.15, 11.16, 11.17, 11.18, 11.19, 11.2,
    11.209999999999999, 11.22, 11.23, 11.24, 11.25, 11.26, 11.27, 11.28, 11.29,
    11.299999999999999, 11.31, 11.32, 11.33, 11.34, 11.35, 11.36, 11.37,
    11.379999999999999, 11.39, 11.4, 11.41, 11.42, 11.43, 11.44, 11.45,
    11.459999999999999, 11.47, 11.48, 11.49, 11.5, 11.51, 11.52, 11.53, 11.54,
    11.549999999999999, 11.56, 11.57, 11.58, 11.59, 11.6, 11.61, 11.62,
    11.629999999999999, 11.64, 11.65, 11.66, 11.67, 11.68, 11.69, 11.7,
    11.709999999999999, 11.72, 11.73, 11.74, 11.75, 11.76, 11.77, 11.78, 11.79,
    11.799999999999999, 11.81, 11.82, 11.83, 11.84, 11.85, 11.86, 11.87,
    11.879999999999999, 11.89, 11.9, 11.91, 11.92, 11.93, 11.94, 11.95,
    11.959999999999999, 11.97, 11.98, 11.99, 12.0, 12.01, 12.02, 12.03, 12.04,
    12.049999999999999, 12.06, 12.07, 12.08, 12.09, 12.1, 12.11, 12.12,
    12.129999999999999, 12.14, 12.15, 12.16, 12.17, 12.18, 12.19, 12.2,
    12.209999999999999, 12.22, 12.23, 12.24, 12.25, 12.26, 12.27, 12.28, 12.29,
    12.299999999999999, 12.31, 12.32, 12.33, 12.34, 12.35, 12.36, 12.37,
    12.379999999999999, 12.39, 12.4, 12.41, 12.42, 12.43, 12.44, 12.45,
    12.459999999999999, 12.47, 12.48, 12.49, 12.5, 12.51, 12.52, 12.53, 12.54,
    12.549999999999999, 12.56, 12.57, 12.58, 12.59, 12.6, 12.61, 12.62,
    12.629999999999999, 12.64, 12.65, 12.66, 12.67, 12.68, 12.69, 12.7,
    12.709999999999999, 12.72, 12.73, 12.74, 12.75, 12.76, 12.77, 12.78, 12.79,
    12.799999999999999, 12.81, 12.82, 12.83, 12.84, 12.85, 12.86, 12.87,
    12.879999999999999, 12.89, 12.9, 12.91, 12.92, 12.93, 12.94, 12.95,
    12.959999999999999, 12.97, 12.98, 12.99, 13.0, 13.01, 13.02, 13.03, 13.04,
    13.049999999999999, 13.06, 13.07, 13.08, 13.09, 13.1, 13.11, 13.12,
    13.129999999999999, 13.14, 13.15, 13.16, 13.17, 13.18, 13.19, 13.2,
    13.209999999999999, 13.22, 13.23, 13.24, 13.25, 13.26, 13.27, 13.28, 13.29,
    13.299999999999999, 13.31, 13.32, 13.33, 13.34, 13.35, 13.36, 13.37,
    13.379999999999999, 13.39, 13.4, 13.41, 13.42, 13.43, 13.44, 13.45,
    13.459999999999999, 13.47, 13.48, 13.49, 13.5, 13.51, 13.52, 13.53, 13.54,
    13.55, 13.56, 13.57, 13.58, 13.59, 13.6, 13.61, 13.62, 13.629999999999999,
    13.64, 13.65, 13.66, 13.67, 13.68, 13.69, 13.7, 13.709999999999999, 13.72,
    13.73, 13.74, 13.75, 13.76, 13.77, 13.78, 13.79, 13.8, 13.81, 13.82, 13.83,
    13.84, 13.85, 13.86, 13.87, 13.879999999999999, 13.89, 13.9, 13.91, 13.92,
    13.93, 13.94, 13.95, 13.959999999999999, 13.97, 13.98, 13.99, 14.0, 14.01,
    14.02, 14.030000000000001, 14.04, 14.05, 14.059999999999999, 14.07, 14.08,
    14.09, 14.1, 14.11, 14.120000000000001, 14.129999999999999, 14.14,
    14.149999999999999, 14.16, 14.17, 14.18, 14.19, 14.2, 14.21,
    14.219999999999999, 14.23, 14.24, 14.25, 14.26, 14.27, 14.280000000000001,
    14.29, 14.3, 14.309999999999999, 14.32, 14.33, 14.34, 14.35, 14.36,
    14.370000000000001, 14.379999999999999, 14.39, 14.399999999999999, 14.41,
    14.42, 14.43, 14.44, 14.45, 14.46, 14.469999999999999, 14.48, 14.49, 14.5,
    14.51, 14.52, 14.530000000000001, 14.54, 14.55, 14.559999999999999, 14.57,
    14.58, 14.59, 14.6, 14.61, 14.620000000000001, 14.629999999999999, 14.64,
    14.649999999999999, 14.66, 14.67, 14.68, 14.69, 14.7, 14.71,
    14.719999999999999, 14.73, 14.74, 14.75, 14.76, 14.77, 14.780000000000001,
    14.79, 14.8, 14.809999999999999, 14.82, 14.83, 14.84, 14.85, 14.86,
    14.870000000000001, 14.879999999999999, 14.89, 14.899999999999999, 14.91,
    14.92, 14.93, 14.94, 14.95, 14.96, 14.969999999999999, 14.98, 14.99, 15.0,
    15.01, 15.02, 15.030000000000001, 15.04, 15.05, 15.059999999999999, 15.07,
    15.08, 15.09, 15.1, 15.11, 15.120000000000001, 15.129999999999999, 15.14,
    15.149999999999999, 15.16, 15.17, 15.18, 15.19, 15.2, 15.21,
    15.219999999999999, 15.23, 15.24, 15.25, 15.26, 15.27, 15.280000000000001,
    15.29, 15.3, 15.309999999999999, 15.32, 15.33, 15.34, 15.35, 15.36,
    15.370000000000001, 15.379999999999999, 15.39, 15.399999999999999, 15.41,
    15.42, 15.43, 15.44, 15.45, 15.46, 15.469999999999999, 15.48, 15.49, 15.5,
    15.51, 15.52, 15.530000000000001, 15.54, 15.55, 15.559999999999999, 15.57,
    15.58, 15.59, 15.6, 15.61, 15.620000000000001, 15.629999999999999, 15.64,
    15.649999999999999, 15.66, 15.67, 15.68, 15.69, 15.7, 15.71,
    15.719999999999999, 15.73, 15.74, 15.75, 15.76, 15.77, 15.780000000000001,
    15.79, 15.8, 15.809999999999999, 15.82, 15.83, 15.84, 15.85, 15.86,
    15.870000000000001, 15.879999999999999, 15.89, 15.899999999999999, 15.91,
    15.92, 15.93, 15.94, 15.95, 15.96, 15.969999999999999, 15.98, 15.99, 16.0,
    16.009999999999998, 16.02, 16.03, 16.04, 16.05, 16.06, 16.07, 16.08, 16.09,
    16.1, 16.11, 16.12, 16.13, 16.14, 16.15, 16.16, 16.17, 16.18,
    16.189999999999998, 16.2, 16.21, 16.22, 16.23, 16.240000000000002, 16.25,
    16.259999999999998, 16.27, 16.28, 16.29, 16.3, 16.31, 16.32, 16.33, 16.34,
    16.35, 16.36, 16.37, 16.38, 16.39, 16.4, 16.41, 16.42, 16.43,
    16.439999999999998, 16.45, 16.46, 16.47, 16.48, 16.490000000000002, 16.5,
    16.509999999999998, 16.52, 16.53, 16.54, 16.55, 16.56, 16.57, 16.58, 16.59,
    16.6, 16.61, 16.62, 16.63, 16.64, 16.65, 16.66, 16.67, 16.68,
    16.689999999999998, 16.7, 16.71, 16.72, 16.73, 16.740000000000002, 16.75,
    16.759999999999998, 16.77, 16.78, 16.79, 16.8, 16.81, 16.82, 16.83, 16.84,
    16.85, 16.86, 16.87, 16.88, 16.89, 16.9, 16.91, 16.92, 16.93,
    16.939999999999998, 16.95, 16.96, 16.97, 16.98, 16.990000000000002, 17.0,
    17.009999999999998, 17.02, 17.03, 17.04, 17.05, 17.06, 17.07, 17.08, 17.09,
    17.1, 17.11, 17.12, 17.13, 17.14, 17.15, 17.16, 17.17, 17.18,
    17.189999999999998, 17.2, 17.21, 17.22, 17.23, 17.240000000000002, 17.25,
    17.259999999999998, 17.27, 17.28, 17.29, 17.3, 17.31, 17.32, 17.33, 17.34,
    17.35, 17.36, 17.37, 17.38, 17.39, 17.4, 17.41, 17.42, 17.43,
    17.439999999999998, 17.45, 17.46, 17.47, 17.48, 17.490000000000002, 17.5,
    17.509999999999998, 17.52, 17.53, 17.54, 17.55, 17.56, 17.57, 17.58, 17.59,
    17.6, 17.61, 17.62, 17.63, 17.64, 17.65, 17.66, 17.67, 17.68,
    17.689999999999998, 17.7, 17.71, 17.72, 17.73, 17.740000000000002, 17.75,
    17.759999999999998, 17.77, 17.78, 17.79, 17.8, 17.81, 17.82, 17.83, 17.84,
    17.85, 17.86, 17.87, 17.88, 17.89, 17.9, 17.91, 17.92, 17.93,
    17.939999999999998, 17.95, 17.96, 17.97, 17.98, 17.990000000000002, 18.0,
    18.009999999999998, 18.02, 18.03, 18.04, 18.05, 18.06, 18.07, 18.08, 18.09,
    18.1, 18.11, 18.12, 18.13, 18.14, 18.15, 18.16, 18.17, 18.18, 18.19, 18.2,
    18.21, 18.22, 18.23, 18.24, 18.25, 18.259999999999998, 18.27, 18.28, 18.29,
    18.3, 18.31, 18.32, 18.33, 18.34, 18.35, 18.36, 18.37, 18.38, 18.39, 18.4,
    18.41, 18.42, 18.43, 18.44, 18.45, 18.46, 18.47, 18.48, 18.49, 18.5,
    18.509999999999998, 18.52, 18.53, 18.54, 18.55, 18.56, 18.57, 18.58, 18.59,
    18.6, 18.61, 18.62, 18.63, 18.64, 18.65, 18.66, 18.67, 18.68, 18.69, 18.7,
    18.71, 18.72, 18.73, 18.74, 18.75, 18.759999999999998, 18.77, 18.78, 18.79,
    18.8, 18.81, 18.82, 18.83, 18.84, 18.85, 18.86, 18.87, 18.88, 18.89, 18.9,
    18.91, 18.92, 18.93, 18.94, 18.95, 18.96, 18.97, 18.98, 18.99, 19.0,
    19.009999999999998, 19.02, 19.03, 19.04, 19.05, 19.06, 19.07, 19.08, 19.09,
    19.1, 19.11, 19.12, 19.13, 19.14, 19.15, 19.16, 19.17, 19.18, 19.19, 19.2,
    19.21, 19.22, 19.23, 19.24, 19.25, 19.259999999999998, 19.27, 19.28, 19.29,
    19.3, 19.31, 19.32, 19.33, 19.34, 19.35, 19.36, 19.37, 19.38, 19.39, 19.4,
    19.41, 19.42, 19.43, 19.44, 19.45, 19.46, 19.47, 19.48, 19.49, 19.5,
    19.509999999999998, 19.52, 19.53, 19.54, 19.55, 19.56, 19.57, 19.58, 19.59,
    19.6, 19.61, 19.62, 19.63, 19.64, 19.65, 19.66, 19.67, 19.68, 19.69, 19.7,
    19.71, 19.72, 19.73, 19.74, 19.75, 19.759999999999998, 19.77, 19.78, 19.79,
    19.8, 19.81, 19.82, 19.83, 19.84, 19.85, 19.86, 19.87, 19.88, 19.89, 19.9,
    19.91, 19.92, 19.93, 19.94, 19.95, 19.96, 19.97, 19.98, 19.99, 20.0, 20.01,
    20.02, 20.03, 20.04, 20.05, 20.06, 20.07, 20.08, 20.09, 20.1, 20.11, 20.12,
    20.13, 20.14, 20.15, 20.16, 20.17, 20.18, 20.19, 20.2, 20.21, 20.22, 20.23,
    20.24, 20.25, 20.26, 20.27, 20.28, 20.29, 20.3, 20.31, 20.32, 20.33, 20.34,
    20.35, 20.36, 20.37, 20.38, 20.39, 20.4, 20.41, 20.42, 20.43, 20.44, 20.45,
    20.46, 20.47, 20.48, 20.49, 20.5, 20.51, 20.52, 20.53, 20.54, 20.55, 20.56,
    20.57, 20.58, 20.59, 20.6, 20.61, 20.62, 20.63, 20.64, 20.65, 20.66, 20.67,
    20.68, 20.69, 20.7, 20.71, 20.72, 20.73, 20.74, 20.75, 20.76, 20.77, 20.78,
    20.79, 20.8, 20.81, 20.82, 20.83, 20.84, 20.85, 20.86, 20.87, 20.88, 20.89,
    20.9, 20.91, 20.92, 20.93, 20.94, 20.95, 20.96, 20.97, 20.98, 20.99, 21.0,
    21.01, 21.02, 21.03, 21.04, 21.05, 21.06, 21.07, 21.08, 21.09, 21.1, 21.11,
    21.12, 21.13, 21.14, 21.15, 21.16, 21.17, 21.18, 21.19, 21.2, 21.21, 21.22,
    21.23, 21.24, 21.25, 21.26, 21.27, 21.28, 21.29, 21.3, 21.31, 21.32, 21.33,
    21.34, 21.35, 21.36, 21.37, 21.38, 21.39, 21.4, 21.41, 21.42, 21.43, 21.44,
    21.45, 21.46, 21.47, 21.48, 21.49, 21.5, 21.51, 21.52, 21.53, 21.54, 21.55,
    21.56, 21.57, 21.58, 21.59, 21.6, 21.61, 21.62, 21.63, 21.64, 21.65, 21.66,
    21.67, 21.68, 21.69, 21.7, 21.71, 21.72, 21.73, 21.74, 21.75, 21.76, 21.77,
    21.78, 21.79, 21.8, 21.81, 21.82, 21.83, 21.84, 21.85, 21.86, 21.87, 21.88,
    21.89, 21.9, 21.91, 21.92, 21.93, 21.94, 21.95, 21.96, 21.97, 21.98, 21.99,
    22.0 };

  double y[2201];
  double z[2201];
  emxArray_real_T *b_x;
  emxArray_real_T *b_y;
  emxArray_real_T *b_z;
  static const double dv67[266] = { 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
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

  double logP0;
  double apnd;
  double cdiff;
  int n;
  *h = 0.01;

  /*  This function evaluates the convolution of two inverse Gaussian */
  /*  distributions at vector t. */
  /*  If the variance in one of the distributions is very small so that the  */
  /*  distribution is close the a Dirac delta function, the convolution */
  /*  is approximated as a shifted inverse Gaussian, that is,  */
  /*  one part of the cell cycle is treated as deterministic in length.   */
  /*  In this case, the shift or lag is the average time to complete  */
  /*  the part of the cycle that has the smallest standard deviation. */
  /*  t is a vector of times to divide (or times to complete two parts of the */
  /*  cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1, */
  /*  s2=sigma2. */
  /*  Input parameters: */
  /*  t = the list of points at which to evaluate the distribution */
  /*  m1 = mu for the first distribution */
  /*  s1 = sigma for the first distribution */
  /*  m2 = mu for the second distribution */
  /*  s2 = sigma for the second distribution */
  /*  h = step size */
  /*  Outputs: */
  /*  P = probability at each point corresponding to points in t */
  /*  h = final step size */
  /*  flag = indicates if we used the Dirac delta */
  /*  E = is the relative error in the likelihood of the data due to the numerical integration */
  /*  log the parameters we were called with */
  /* fprintf(1,'convolv_2invG(m1=%.17g, s1=%.17g, m2=%.17g, s2=%.17g)\n', ...  */
  /*     m1, s1, m2, s2); */
  /* fprintf('    convolv_2invG_adapt(m1=%f s1=%f m2=%f s2=%f)\n',m1,s1,m2,s2); */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:37' flag=0; */
  *flag = 0.0;

  /*  a constant that is used in determining if the Dirac approximation should be applied. */
  /* 'convolv_2invG_adapt_nov:40' eps=.01; */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:43' E=Inf; */
  *E = rtInf;

  /*  number of points to be evaluated */
  /* 'convolv_2invG_adapt_nov:46' n=length(t); */
  /*  store m and s for each sub-distribution in a list */
  /* 'convolv_2invG_adapt_nov:49' m=[m1 m2]; */
  /* 'convolv_2invG_adapt_nov:50' s=[s1 s2]; */
  /*  make sure we dont have negative values for m or s by replacement with */
  /*  zero if we do */
  /* m=max(m,0); */
  /* s=max(s,0); */
  /*  reflect negative parameters into the positive space so we can use */
  /*  unconstrained optimization */
  /* 'convolv_2invG_adapt_nov:58' m=abs(m); */
  x[0] = m1;
  x[1] = m2;
  for (k = 0; k < 2; k++) {
    m[k] = fabs(x[k]);
  }

  /* 'convolv_2invG_adapt_nov:59' s=abs(s); */
  x[0] = s1;
  x[1] = s2;

  /*  find the variance for both sub-distributions */
  /* 'convolv_2invG_adapt_nov:62' v=(s.^2)./(m.^3); */
  for (k = 0; k < 2; k++) {
    ndbl = fabs(x[k]);
    v[k] = rt_powd_snf(ndbl, 2.0);
    s[k] = ndbl;
  }

  /*  find the standard deviation for each sub-distribution */
  /* 'convolv_2invG_adapt_nov:65' sd=v.^.5; */
  for (k = 0; k < 2; k++) {
    ndbl = v[k] / rt_powd_snf(m[k], 3.0);
    sd[k] = rt_powd_snf(ndbl, 0.5);
    v[k] = ndbl;
  }

  /*  reorder m and s so that the sub-distribution with the smallest */
  /*  variance comes first.  So the fist part might be approximated as a Dirac delta. */
  /* 'convolv_2invG_adapt_nov:69' [v,I]=sort(v); */
  sort(v, iidx);

  /* 'convolv_2invG_adapt_nov:70' m=m(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    v[nm1d2] = m[iidx[nm1d2] - 1];
    x[nm1d2] = iidx[nm1d2];
  }

  /* 'convolv_2invG_adapt_nov:71' s=s(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    m[nm1d2] = v[nm1d2];
    b_s[nm1d2] = s[(int)x[nm1d2] - 1];
  }

  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    s[nm1d2] = b_s[nm1d2];
  }

  /*  T1 appears to be unused, should probably remove */
  /*  T2 is only used inside FUNCTION CHECK_APPROXIMATABLE */
  /*  and should probably be moved inside that scope */
  /*  T2 is the mode of the second distribution, determined by a precomputed analytic expression. */
  /* 'convolv_2invG_adapt_nov:78' T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1))); */
  /* 'convolv_2invG_adapt_nov:79' T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2))); */
  /* When called from convolv_3invG all of the components of t may be negative. */
  /* In this case, the probability of each is set to realmin */
  /* 'convolv_2invG_adapt_nov:83' if max(t)<=0 */
  /* 'convolv_2invG_adapt_nov:87' else */
  /*  find the largest point in t */
  /* 'convolv_2invG_adapt_nov:89' Maxt=max(t); */
  /*  produce a range of evenly spaced points to evaluate at between */
  /*  0 and Maxt with step size h. the even spacing is important when */
  /*  calculating the convolution later */
  /* 'convolv_2invG_adapt_nov:94' x=0:h:Maxt; */
  /*  evaluate the sub-distributions at each point in x */
  /* 'convolv_2invG_adapt_nov:97' y=onestagepdf2(x,m(1),s(1)); */
  b_onestagepdf2(dv66, m[0], s[0], y);

  /* 'convolv_2invG_adapt_nov:98' z=onestagepdf2(x,m(2),s(2)); */
  b_onestagepdf2(dv66, m[1], s[1], z);

  /*  if the first pdf is very concentrated, check to see if it can be */
  /*  approximated as a point-mass distribution */
  /* 'convolv_2invG_adapt_nov:102' if sd(1)<.01 */
  emxInit_real_T(&b_x, 2);
  emxInit_real_T(&b_y, 2);
  emxInit_real_T(&b_z, 2);
  if (sd[0] < 0.01) {
    /* 'convolv_2invG_adapt_nov:104' check2 = tailmass(m, s, eps, T2, sd); */
    ndbl = tailmass(m, s, 1.0 / m[1] * (sqrt(1.0 + 2.25 * (rt_powd_snf(s[1], 4.0)
      / (m[1] * m[1]))) - 1.5 * (s[1] * s[1] / m[1])), sd);

    /* 'convolv_2invG_adapt_nov:106' if  check2<=eps/3 */
    if (ndbl <= 0.0033333333333333335) {
      /* If there is not much probability in the tails of the convolution, we use a Dirac Delta for part 1. */
      /* Flag is set to 1 to indicate this. */
      /* 'convolv_2invG_adapt_nov:109' flag=1; */
      *flag = 1.0;

      /* 'convolv_2invG_adapt_nov:110' sigma=s(2); */
      /* 'convolv_2invG_adapt_nov:111' mu=m(2); */
      /* l for lag. */
      /* 'convolv_2invG_adapt_nov:113' l=1/m(1); */
      /* 'convolv_2invG_adapt_nov:114' P=onestagepdf_lag(t,mu,sigma,l); */
      onestagepdf_lag(dv67, m[1], s[1], 1.0 / m[0], P);
    } else {
      /* 'convolv_2invG_adapt_nov:115' else */
      /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
      /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
      /*  Outputs: P */
      /* 'convolv_2invG_adapt_nov:120' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
      approxconvolv(z, y, dv67, dv66, P, &logP0);

      /* Keep reducing the step size in the numerical integration until we are happy with the error. */
      /* 'convolv_2invG_adapt_nov:124' while E>=.001*abs(logP0) */
      while (*E >= 0.001 * fabs(logP0)) {
        /* 'convolv_2invG_adapt_nov:125' h1=.5*h; */
        *h *= 0.5;

        /* 'convolv_2invG_adapt_nov:126' x=0:h1:Maxt; */
        if ((*h == 0.0) || (*h < 0.0)) {
          nm1d2 = b_x->size[0] * b_x->size[1];
          b_x->size[0] = 1;
          b_x->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
        } else if (0.0 == *h) {
          nm1d2 = b_x->size[0] * b_x->size[1];
          b_x->size[0] = 1;
          b_x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
          for (nm1d2 = 0; nm1d2 < 1; nm1d2++) {
            b_x->data[b_x->size[0] * nm1d2] = 0.0 * (double)nm1d2;
          }
        } else {
          ndbl = floor(22.0 / *h + 0.5);
          apnd = ndbl * *h;
          if (*h > 0.0) {
            cdiff = apnd - 22.0;
          } else {
            cdiff = 22.0 - apnd;
          }

          if (fabs(cdiff) < 9.7699626167013776E-15) {
            ndbl++;
            apnd = 22.0;
          } else if (cdiff > 0.0) {
            apnd = (ndbl - 1.0) * *h;
          } else {
            ndbl++;
          }

          if (ndbl >= 0.0) {
            n = (int)ndbl;
          } else {
            n = 0;
          }

          nm1d2 = b_x->size[0] * b_x->size[1];
          b_x->size[0] = 1;
          b_x->size[1] = n;
          emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
          if (n > 0) {
            b_x->data[0] = 0.0;
            if (n > 1) {
              b_x->data[n - 1] = apnd;
              nm1d2 = (n - 1) / 2;
              for (k = 1; k < nm1d2; k++) {
                ndbl = (double)k * *h;
                b_x->data[k] = ndbl;
                b_x->data[(n - k) - 1] = apnd - ndbl;
              }

              if (nm1d2 << 1 == n - 1) {
                b_x->data[nm1d2] = apnd / 2.0;
              } else {
                ndbl = (double)nm1d2 * *h;
                b_x->data[nm1d2] = ndbl;
                b_x->data[nm1d2 + 1] = apnd - ndbl;
              }
            }
          }
        }

        /* 'convolv_2invG_adapt_nov:127' y=onestagepdf2(x,m(1),s(1)); */
        d_onestagepdf2(b_x, m[0], s[0], b_y);

        /* 'convolv_2invG_adapt_nov:128' z=onestagepdf2(x,m(2),s(2)); */
        d_onestagepdf2(b_x, m[1], s[1], b_z);

        /* 'convolv_2invG_adapt_nov:130' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
        b_approxconvolv(b_z, b_y, *h, b_x, P, &ndbl);

        /* 'convolv_2invG_adapt_nov:132' E=abs(logP1-logP0); */
        *E = fabs(ndbl - logP0);

        /* 'convolv_2invG_adapt_nov:133' P0=P1; */
        /* 'convolv_2invG_adapt_nov:134' logP0=logP1; */
        logP0 = ndbl;

        /* 'convolv_2invG_adapt_nov:135' h=h1; */
      }

      /* 'convolv_2invG_adapt_nov:137' P=P0; */
      /*  END FUNCTION DOTHECONVOLUTION_OUTER */
    }

    /*  pdf is not very concentrated so compute the convolution directly */
  } else {
    /* 'convolv_2invG_adapt_nov:141' else */
    /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
    /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
    /*  Outputs: P */
    /* 'convolv_2invG_adapt_nov:146' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
    approxconvolv(z, y, dv67, dv66, P, &logP0);

    /* 'convolv_2invG_adapt_nov:148' while E>=.001*abs(logP0) */
    while (*E >= 0.001 * fabs(logP0)) {
      /* 'convolv_2invG_adapt_nov:149' h1=.5*h; */
      *h *= 0.5;

      /* 'convolv_2invG_adapt_nov:150' x=0:h1:Maxt; */
      if ((*h == 0.0) || (*h < 0.0)) {
        nm1d2 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
      } else if (0.0 == *h) {
        nm1d2 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
        for (nm1d2 = 0; nm1d2 < 1; nm1d2++) {
          b_x->data[b_x->size[0] * nm1d2] = 0.0 * (double)nm1d2;
        }
      } else {
        ndbl = floor(22.0 / *h + 0.5);
        apnd = ndbl * *h;
        if (*h > 0.0) {
          cdiff = apnd - 22.0;
        } else {
          cdiff = 22.0 - apnd;
        }

        if (fabs(cdiff) < 9.7699626167013776E-15) {
          ndbl++;
          apnd = 22.0;
        } else if (cdiff > 0.0) {
          apnd = (ndbl - 1.0) * *h;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          n = (int)ndbl;
        } else {
          n = 0;
        }

        nm1d2 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = n;
        emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
        if (n > 0) {
          b_x->data[0] = 0.0;
          if (n > 1) {
            b_x->data[n - 1] = apnd;
            nm1d2 = (n - 1) / 2;
            for (k = 1; k < nm1d2; k++) {
              ndbl = (double)k * *h;
              b_x->data[k] = ndbl;
              b_x->data[(n - k) - 1] = apnd - ndbl;
            }

            if (nm1d2 << 1 == n - 1) {
              b_x->data[nm1d2] = apnd / 2.0;
            } else {
              ndbl = (double)nm1d2 * *h;
              b_x->data[nm1d2] = ndbl;
              b_x->data[nm1d2 + 1] = apnd - ndbl;
            }
          }
        }
      }

      /* 'convolv_2invG_adapt_nov:151' y=onestagepdf2(x,m(1),s(1)); */
      d_onestagepdf2(b_x, m[0], s[0], b_y);

      /* 'convolv_2invG_adapt_nov:152' z=onestagepdf2(x,m(2),s(2)); */
      d_onestagepdf2(b_x, m[1], s[1], b_z);

      /* 'convolv_2invG_adapt_nov:154' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
      b_approxconvolv(b_z, b_y, *h, b_x, P, &ndbl);

      /* 'convolv_2invG_adapt_nov:157' E=abs(logP1-logP0); */
      *E = fabs(ndbl - logP0);

      /* 'convolv_2invG_adapt_nov:158' P0=P1; */
      /* 'convolv_2invG_adapt_nov:159' logP0=logP1; */
      logP0 = ndbl;

      /* 'convolv_2invG_adapt_nov:160' h=h1; */
    }

    /* 'convolv_2invG_adapt_nov:162' P=P0; */
    /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
  }

  emxFree_real_T(&b_z);
  emxFree_real_T(&b_y);
  emxFree_real_T(&b_x);
}
#endif

/*
 * function [P,h,flag,E]=convolv_2invG_adapt_nov(t,m1,s1,m2,s2,h)
 */
#ifdef _OLD_MATLAB_CODE
void c_convolv_2invG_adapt_nov(double m1, double s1, double m2, double s2,
  double P[266], double *h, double *flag, double *E)
{
  double x[2];
  int k;
  double m[2];
  double ndbl;
  double v[2];
  double s[2];
  double sd[2];
  int iidx[2];
  int nm1d2;
  double b_s[2];
  static double b_x[22001];
  static double y[22001];
  static double z[22001];
  emxArray_real_T *c_x;
  emxArray_real_T *b_y;
  emxArray_real_T *b_z;
  static const double dv68[266] = { 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
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

  double logP0;
  double apnd;
  double cdiff;
  int n;
  *h = 0.001;

  /*  This function evaluates the convolution of two inverse Gaussian */
  /*  distributions at vector t. */
  /*  If the variance in one of the distributions is very small so that the  */
  /*  distribution is close the a Dirac delta function, the convolution */
  /*  is approximated as a shifted inverse Gaussian, that is,  */
  /*  one part of the cell cycle is treated as deterministic in length.   */
  /*  In this case, the shift or lag is the average time to complete  */
  /*  the part of the cycle that has the smallest standard deviation. */
  /*  t is a vector of times to divide (or times to complete two parts of the */
  /*  cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1, */
  /*  s2=sigma2. */
  /*  Input parameters: */
  /*  t = the list of points at which to evaluate the distribution */
  /*  m1 = mu for the first distribution */
  /*  s1 = sigma for the first distribution */
  /*  m2 = mu for the second distribution */
  /*  s2 = sigma for the second distribution */
  /*  h = step size */
  /*  Outputs: */
  /*  P = probability at each point corresponding to points in t */
  /*  h = final step size */
  /*  flag = indicates if we used the Dirac delta */
  /*  E = is the relative error in the likelihood of the data due to the numerical integration */
  /*  log the parameters we were called with */
  /* fprintf(1,'convolv_2invG(m1=%.17g, s1=%.17g, m2=%.17g, s2=%.17g)\n', ...  */
  /*     m1, s1, m2, s2); */
  /* fprintf('    convolv_2invG_adapt(m1=%f s1=%f m2=%f s2=%f)\n',m1,s1,m2,s2); */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:37' flag=0; */
  *flag = 0.0;

  /*  a constant that is used in determining if the Dirac approximation should be applied. */
  /* 'convolv_2invG_adapt_nov:40' eps=.01; */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:43' E=Inf; */
  *E = rtInf;

  /*  number of points to be evaluated */
  /* 'convolv_2invG_adapt_nov:46' n=length(t); */
  /*  store m and s for each sub-distribution in a list */
  /* 'convolv_2invG_adapt_nov:49' m=[m1 m2]; */
  /* 'convolv_2invG_adapt_nov:50' s=[s1 s2]; */
  /*  make sure we dont have negative values for m or s by replacement with */
  /*  zero if we do */
  /* m=max(m,0); */
  /* s=max(s,0); */
  /*  reflect negative parameters into the positive space so we can use */
  /*  unconstrained optimization */
  /* 'convolv_2invG_adapt_nov:58' m=abs(m); */
  x[0] = m1;
  x[1] = m2;
  for (k = 0; k < 2; k++) {
    m[k] = fabs(x[k]);
  }

  /* 'convolv_2invG_adapt_nov:59' s=abs(s); */
  x[0] = s1;
  x[1] = s2;

  /*  find the variance for both sub-distributions */
  /* 'convolv_2invG_adapt_nov:62' v=(s.^2)./(m.^3); */
  for (k = 0; k < 2; k++) {
    ndbl = fabs(x[k]);
    v[k] = rt_powd_snf(ndbl, 2.0);
    s[k] = ndbl;
  }

  /*  find the standard deviation for each sub-distribution */
  /* 'convolv_2invG_adapt_nov:65' sd=v.^.5; */
  for (k = 0; k < 2; k++) {
    ndbl = v[k] / rt_powd_snf(m[k], 3.0);
    sd[k] = rt_powd_snf(ndbl, 0.5);
    v[k] = ndbl;
  }

  /*  reorder m and s so that the sub-distribution with the smallest */
  /*  variance comes first.  So the fist part might be approximated as a Dirac delta. */
  /* 'convolv_2invG_adapt_nov:69' [v,I]=sort(v); */
  sort(v, iidx);

  /* 'convolv_2invG_adapt_nov:70' m=m(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    v[nm1d2] = m[iidx[nm1d2] - 1];
    x[nm1d2] = iidx[nm1d2];
  }

  /* 'convolv_2invG_adapt_nov:71' s=s(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    m[nm1d2] = v[nm1d2];
    b_s[nm1d2] = s[(int)x[nm1d2] - 1];
  }

  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    s[nm1d2] = b_s[nm1d2];
  }

  /*  T1 appears to be unused, should probably remove */
  /*  T2 is only used inside FUNCTION CHECK_APPROXIMATABLE */
  /*  and should probably be moved inside that scope */
  /*  T2 is the mode of the second distribution, determined by a precomputed analytic expression. */
  /* 'convolv_2invG_adapt_nov:78' T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1))); */
  /* 'convolv_2invG_adapt_nov:79' T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2))); */
  /* When called from convolv_3invG all of the components of t may be negative. */
  /* In this case, the probability of each is set to realmin */
  /* 'convolv_2invG_adapt_nov:83' if max(t)<=0 */
  /* 'convolv_2invG_adapt_nov:87' else */
  /*  find the largest point in t */
  /* 'convolv_2invG_adapt_nov:89' Maxt=max(t); */
  /*  produce a range of evenly spaced points to evaluate at between */
  /*  0 and Maxt with step size h. the even spacing is important when */
  /*  calculating the convolution later */
  /* 'convolv_2invG_adapt_nov:94' x=0:h:Maxt; */
  b_x[0] = 0.0;
  b_x[22000] = 22.0;
  for (k = 0; k < 10999; k++) {
    ndbl = ((double)k + 1.0) * 0.001;
    b_x[k + 1] = ndbl;
    b_x[21999 - k] = 22.0 - ndbl;
  }

  b_x[11000] = 11.0;

  /*  evaluate the sub-distributions at each point in x */
  /* 'convolv_2invG_adapt_nov:97' y=onestagepdf2(x,m(1),s(1)); */
  e_onestagepdf2(b_x, m[0], s[0], y);

  /* 'convolv_2invG_adapt_nov:98' z=onestagepdf2(x,m(2),s(2)); */
  e_onestagepdf2(b_x, m[1], s[1], z);

  /*  if the first pdf is very concentrated, check to see if it can be */
  /*  approximated as a point-mass distribution */
  /* 'convolv_2invG_adapt_nov:102' if sd(1)<.01 */
  emxInit_real_T(&c_x, 2);
  emxInit_real_T(&b_y, 2);
  emxInit_real_T(&b_z, 2);
  if (sd[0] < 0.01) {
    /* 'convolv_2invG_adapt_nov:104' check2 = tailmass(m, s, eps, T2, sd); */
    ndbl = tailmass(m, s, 1.0 / m[1] * (sqrt(1.0 + 2.25 * (rt_powd_snf(s[1], 4.0)
      / (m[1] * m[1]))) - 1.5 * (s[1] * s[1] / m[1])), sd);

    /* 'convolv_2invG_adapt_nov:106' if  check2<=eps/3 */
    if (ndbl <= 0.0033333333333333335) {
      /* If there is not much probability in the tails of the convolution, we use a Dirac Delta for part 1. */
      /* Flag is set to 1 to indicate this. */
      /* 'convolv_2invG_adapt_nov:109' flag=1; */
      *flag = 1.0;

      /* 'convolv_2invG_adapt_nov:110' sigma=s(2); */
      /* 'convolv_2invG_adapt_nov:111' mu=m(2); */
      /* l for lag. */
      /* 'convolv_2invG_adapt_nov:113' l=1/m(1); */
      /* 'convolv_2invG_adapt_nov:114' P=onestagepdf_lag(t,mu,sigma,l); */
      onestagepdf_lag(dv68, m[1], s[1], 1.0 / m[0], P);
    } else {
      /* 'convolv_2invG_adapt_nov:115' else */
      /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
      /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
      /*  Outputs: P */
      /* 'convolv_2invG_adapt_nov:120' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
      c_approxconvolv(z, y, dv68, b_x, P, &logP0);

      /* Keep reducing the step size in the numerical integration until we are happy with the error. */
      /* 'convolv_2invG_adapt_nov:124' while E>=.001*abs(logP0) */
      while (*E >= 0.001 * fabs(logP0)) {
        /* 'convolv_2invG_adapt_nov:125' h1=.5*h; */
        *h *= 0.5;

        /* 'convolv_2invG_adapt_nov:126' x=0:h1:Maxt; */
        if ((*h == 0.0) || (*h < 0.0)) {
          nm1d2 = c_x->size[0] * c_x->size[1];
          c_x->size[0] = 1;
          c_x->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)c_x, nm1d2, sizeof(double));
        } else if (0.0 == *h) {
          nm1d2 = c_x->size[0] * c_x->size[1];
          c_x->size[0] = 1;
          c_x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)c_x, nm1d2, sizeof(double));
          for (nm1d2 = 0; nm1d2 < 1; nm1d2++) {
            c_x->data[c_x->size[0] * nm1d2] = 0.0 * (double)nm1d2;
          }
        } else {
          ndbl = floor(22.0 / *h + 0.5);
          apnd = ndbl * *h;
          if (*h > 0.0) {
            cdiff = apnd - 22.0;
          } else {
            cdiff = 22.0 - apnd;
          }

          if (fabs(cdiff) < 9.7699626167013776E-15) {
            ndbl++;
            apnd = 22.0;
          } else if (cdiff > 0.0) {
            apnd = (ndbl - 1.0) * *h;
          } else {
            ndbl++;
          }

          if (ndbl >= 0.0) {
            n = (int)ndbl;
          } else {
            n = 0;
          }

          nm1d2 = c_x->size[0] * c_x->size[1];
          c_x->size[0] = 1;
          c_x->size[1] = n;
          emxEnsureCapacity((emxArray__common *)c_x, nm1d2, sizeof(double));
          if (n > 0) {
            c_x->data[0] = 0.0;
            if (n > 1) {
              c_x->data[n - 1] = apnd;
              nm1d2 = (n - 1) / 2;
              for (k = 1; k < nm1d2; k++) {
                ndbl = (double)k * *h;
                c_x->data[k] = ndbl;
                c_x->data[(n - k) - 1] = apnd - ndbl;
              }

              if (nm1d2 << 1 == n - 1) {
                c_x->data[nm1d2] = apnd / 2.0;
              } else {
                ndbl = (double)nm1d2 * *h;
                c_x->data[nm1d2] = ndbl;
                c_x->data[nm1d2 + 1] = apnd - ndbl;
              }
            }
          }
        }

        /* 'convolv_2invG_adapt_nov:127' y=onestagepdf2(x,m(1),s(1)); */
        d_onestagepdf2(c_x, m[0], s[0], b_y);

        /* 'convolv_2invG_adapt_nov:128' z=onestagepdf2(x,m(2),s(2)); */
        d_onestagepdf2(c_x, m[1], s[1], b_z);

        /* 'convolv_2invG_adapt_nov:130' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
        b_approxconvolv(b_z, b_y, *h, c_x, P, &ndbl);

        /* 'convolv_2invG_adapt_nov:132' E=abs(logP1-logP0); */
        *E = fabs(ndbl - logP0);

        /* 'convolv_2invG_adapt_nov:133' P0=P1; */
        /* 'convolv_2invG_adapt_nov:134' logP0=logP1; */
        logP0 = ndbl;

        /* 'convolv_2invG_adapt_nov:135' h=h1; */
      }

      /* 'convolv_2invG_adapt_nov:137' P=P0; */
      /*  END FUNCTION DOTHECONVOLUTION_OUTER */
    }

    /*  pdf is not very concentrated so compute the convolution directly */
  } else {
    /* 'convolv_2invG_adapt_nov:141' else */
    /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
    /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
    /*  Outputs: P */
    /* 'convolv_2invG_adapt_nov:146' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
    c_approxconvolv(z, y, dv68, b_x, P, &logP0);

    /* 'convolv_2invG_adapt_nov:148' while E>=.001*abs(logP0) */
    while (*E >= 0.001 * fabs(logP0)) {
      /* 'convolv_2invG_adapt_nov:149' h1=.5*h; */
      *h *= 0.5;

      /* 'convolv_2invG_adapt_nov:150' x=0:h1:Maxt; */
      if ((*h == 0.0) || (*h < 0.0)) {
        nm1d2 = c_x->size[0] * c_x->size[1];
        c_x->size[0] = 1;
        c_x->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)c_x, nm1d2, sizeof(double));
      } else if (0.0 == *h) {
        nm1d2 = c_x->size[0] * c_x->size[1];
        c_x->size[0] = 1;
        c_x->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)c_x, nm1d2, sizeof(double));
        for (nm1d2 = 0; nm1d2 < 1; nm1d2++) {
          c_x->data[c_x->size[0] * nm1d2] = 0.0 * (double)nm1d2;
        }
      } else {
        ndbl = floor(22.0 / *h + 0.5);
        apnd = ndbl * *h;
        if (*h > 0.0) {
          cdiff = apnd - 22.0;
        } else {
          cdiff = 22.0 - apnd;
        }

        if (fabs(cdiff) < 9.7699626167013776E-15) {
          ndbl++;
          apnd = 22.0;
        } else if (cdiff > 0.0) {
          apnd = (ndbl - 1.0) * *h;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          n = (int)ndbl;
        } else {
          n = 0;
        }

        nm1d2 = c_x->size[0] * c_x->size[1];
        c_x->size[0] = 1;
        c_x->size[1] = n;
        emxEnsureCapacity((emxArray__common *)c_x, nm1d2, sizeof(double));
        if (n > 0) {
          c_x->data[0] = 0.0;
          if (n > 1) {
            c_x->data[n - 1] = apnd;
            nm1d2 = (n - 1) / 2;
            for (k = 1; k < nm1d2; k++) {
              ndbl = (double)k * *h;
              c_x->data[k] = ndbl;
              c_x->data[(n - k) - 1] = apnd - ndbl;
            }

            if (nm1d2 << 1 == n - 1) {
              c_x->data[nm1d2] = apnd / 2.0;
            } else {
              ndbl = (double)nm1d2 * *h;
              c_x->data[nm1d2] = ndbl;
              c_x->data[nm1d2 + 1] = apnd - ndbl;
            }
          }
        }
      }

      /* 'convolv_2invG_adapt_nov:151' y=onestagepdf2(x,m(1),s(1)); */
      d_onestagepdf2(c_x, m[0], s[0], b_y);

      /* 'convolv_2invG_adapt_nov:152' z=onestagepdf2(x,m(2),s(2)); */
      d_onestagepdf2(c_x, m[1], s[1], b_z);

      /* 'convolv_2invG_adapt_nov:154' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
      b_approxconvolv(b_z, b_y, *h, c_x, P, &ndbl);

      /* 'convolv_2invG_adapt_nov:157' E=abs(logP1-logP0); */
      *E = fabs(ndbl - logP0);

      /* 'convolv_2invG_adapt_nov:158' P0=P1; */
      /* 'convolv_2invG_adapt_nov:159' logP0=logP1; */
      logP0 = ndbl;

      /* 'convolv_2invG_adapt_nov:160' h=h1; */
    }

    /* 'convolv_2invG_adapt_nov:162' P=P0; */
    /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
  }

  emxFree_real_T(&b_z);
  emxFree_real_T(&b_y);
  emxFree_real_T(&c_x);
}
#endif

/*
 * function [P,h,flag,E]=convolv_2invG_adapt_nov(t,m1,s1,m2,s2,h)
 */
#ifdef _OLD_MATLAB_CODE
void convolv_2invG_adapt_nov(double m1, double s1, double m2, double s2, double
  P[266])
{
  double h1;
  double E;
  double x[2];
  int k;
  double m[2];
  double ndbl;
  double v[2];
  double s[2];
  double sd[2];
  int iidx[2];
  int nm1d2;
  double b_s[2];
  static const double dv62[2201] = { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
    0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
    0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32,
    0.33, 0.34, 0.35000000000000003, 0.36, 0.37, 0.38, 0.39, 0.4,
    0.41000000000000003, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47000000000000003, 0.48,
    0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57000000000000006, 0.58,
    0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68,
    0.69000000000000006, 0.70000000000000007, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76,
    0.77, 0.78, 0.79, 0.8, 0.81, 0.82000000000000006, 0.83000000000000007, 0.84,
    0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94000000000000006,
    0.95000000000000007, 0.96, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04,
    1.05, 1.06, 1.07, 1.08, 1.09, 1.1, 1.11, 1.12, 1.1300000000000001,
    1.1400000000000001, 1.1500000000000001, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21,
    1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33, 1.34,
    1.35, 1.36, 1.37, 1.3800000000000001, 1.3900000000000001, 1.4000000000000001,
    1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53,
    1.54, 1.55, 1.56, 1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.6300000000000001,
    1.6400000000000001, 1.6500000000000001, 1.6600000000000001, 1.67, 1.68, 1.69,
    1.7, 1.71, 1.72, 1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82,
    1.83, 1.84, 1.85, 1.86, 1.87, 1.8800000000000001, 1.8900000000000001,
    1.9000000000000001, 1.9100000000000001, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97,
    1.98, 1.99, 2.0, 2.0100000000000002, 2.02, 2.0300000000000002, 2.04, 2.05,
    2.06, 2.07, 2.08, 2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18,
    2.19, 2.2, 2.21, 2.22, 2.23, 2.24, 2.25, 2.2600000000000002, 2.27,
    2.2800000000000002, 2.29, 2.3000000000000003, 2.31, 2.32, 2.33, 2.34, 2.35,
    2.36, 2.37, 2.38, 2.39, 2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48,
    2.49, 2.5, 2.5100000000000002, 2.52, 2.5300000000000002, 2.54,
    2.5500000000000003, 2.56, 2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64,
    2.65, 2.66, 2.67, 2.68, 2.69, 2.7, 2.71, 2.72, 2.73, 2.74, 2.75,
    2.7600000000000002, 2.77, 2.7800000000000002, 2.79, 2.8000000000000003, 2.81,
    2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 2.89, 2.9, 2.91, 2.92, 2.93, 2.94,
    2.95, 2.96, 2.97, 2.98, 2.99, 3.0, 3.0100000000000002, 3.02,
    3.0300000000000002, 3.04, 3.0500000000000003, 3.06, 3.0700000000000003, 3.08,
    3.09, 3.1, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.2, 3.21,
    3.22, 3.23, 3.24, 3.25, 3.2600000000000002, 3.27, 3.2800000000000002, 3.29,
    3.3000000000000003, 3.31, 3.3200000000000003, 3.33, 3.34, 3.35, 3.36, 3.37,
    3.38, 3.39, 3.4, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.48, 3.49, 3.5,
    3.5100000000000002, 3.52, 3.5300000000000002, 3.54, 3.5500000000000003, 3.56,
    3.5700000000000003, 3.58, 3.59, 3.6, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66,
    3.67, 3.68, 3.69, 3.7, 3.71, 3.72, 3.73, 3.74, 3.75, 3.7600000000000002,
    3.77, 3.7800000000000002, 3.79, 3.8000000000000003, 3.81, 3.8200000000000003,
    3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95,
    3.96, 3.97, 3.98, 3.99, 4.0, 4.01, 4.0200000000000005, 4.03, 4.04, 4.05,
    4.0600000000000005, 4.07, 4.08, 4.09, 4.1, 4.11, 4.12, 4.13, 4.14, 4.15,
    4.16, 4.17, 4.18, 4.19, 4.2, 4.21, 4.22, 4.23, 4.24, 4.25, 4.26,
    4.2700000000000005, 4.28, 4.29, 4.3, 4.3100000000000005, 4.32, 4.33, 4.34,
    4.3500000000000005, 4.36, 4.37, 4.38, 4.39, 4.4, 4.41, 4.42, 4.43, 4.44,
    4.45, 4.46, 4.47, 4.48, 4.49, 4.5, 4.51, 4.5200000000000005, 4.53, 4.54,
    4.55, 4.5600000000000005, 4.57, 4.58, 4.59, 4.6000000000000005, 4.61, 4.62,
    4.63, 4.64, 4.65, 4.66, 4.67, 4.68, 4.69, 4.7, 4.71, 4.72, 4.73, 4.74, 4.75,
    4.76, 4.7700000000000005, 4.78, 4.79, 4.8, 4.8100000000000005, 4.82, 4.83,
    4.84, 4.8500000000000005, 4.86, 4.87, 4.88, 4.89, 4.9, 4.91, 4.92, 4.93,
    4.94, 4.95, 4.96, 4.97, 4.98, 4.99, 5.0, 5.01, 5.0200000000000005, 5.03,
    5.04, 5.05, 5.0600000000000005, 5.07, 5.08, 5.09, 5.1000000000000005, 5.11,
    5.12, 5.13, 5.14, 5.15, 5.16, 5.17, 5.18, 5.19, 5.2, 5.21, 5.22, 5.23, 5.24,
    5.25, 5.26, 5.2700000000000005, 5.28, 5.29, 5.3, 5.3100000000000005, 5.32,
    5.33, 5.34, 5.3500000000000005, 5.36, 5.37, 5.38, 5.39, 5.4, 5.41, 5.42,
    5.43, 5.44, 5.45, 5.46, 5.47, 5.48, 5.49, 5.5, 5.51, 5.5200000000000005,
    5.53, 5.54, 5.55, 5.5600000000000005, 5.57, 5.58, 5.59, 5.6000000000000005,
    5.61, 5.62, 5.63, 5.64, 5.65, 5.66, 5.67, 5.68, 5.69, 5.7, 5.71, 5.72, 5.73,
    5.74, 5.75, 5.76, 5.7700000000000005, 5.78, 5.79, 5.8, 5.8100000000000005,
    5.82, 5.83, 5.84, 5.8500000000000005, 5.86, 5.87, 5.88, 5.89, 5.9, 5.91,
    5.92, 5.93, 5.94, 5.95, 5.96, 5.97, 5.98, 5.99, 6.0, 6.01,
    6.0200000000000005, 6.03, 6.04, 6.05, 6.0600000000000005, 6.07, 6.08, 6.09,
    6.1000000000000005, 6.11, 6.12, 6.13, 6.1400000000000006, 6.15, 6.16, 6.17,
    6.18, 6.19, 6.2, 6.21, 6.22, 6.23, 6.24, 6.25, 6.26, 6.2700000000000005,
    6.28, 6.29, 6.3, 6.3100000000000005, 6.32, 6.33, 6.34, 6.3500000000000005,
    6.36, 6.37, 6.38, 6.3900000000000006, 6.4, 6.41, 6.42, 6.43, 6.44, 6.45,
    6.46, 6.47, 6.48, 6.49, 6.5, 6.51, 6.5200000000000005, 6.53, 6.54, 6.55,
    6.5600000000000005, 6.57, 6.58, 6.59, 6.6000000000000005, 6.61, 6.62, 6.63,
    6.6400000000000006, 6.65, 6.66, 6.67, 6.68, 6.69, 6.7, 6.71, 6.72, 6.73,
    6.74, 6.75, 6.76, 6.7700000000000005, 6.78, 6.79, 6.8, 6.8100000000000005,
    6.82, 6.83, 6.84, 6.8500000000000005, 6.86, 6.87, 6.88, 6.8900000000000006,
    6.9, 6.91, 6.92, 6.93, 6.94, 6.95, 6.96, 6.97, 6.98, 6.99, 7.0, 7.01,
    7.0200000000000005, 7.03, 7.04, 7.05, 7.0600000000000005, 7.07, 7.08, 7.09,
    7.1000000000000005, 7.11, 7.12, 7.13, 7.1400000000000006, 7.15, 7.16, 7.17,
    7.18, 7.19, 7.2, 7.21, 7.22, 7.23, 7.24, 7.25, 7.26, 7.2700000000000005,
    7.28, 7.29, 7.3, 7.3100000000000005, 7.32, 7.33, 7.34, 7.3500000000000005,
    7.36, 7.37, 7.38, 7.3900000000000006, 7.4, 7.41, 7.42, 7.43, 7.44, 7.45,
    7.46, 7.47, 7.48, 7.49, 7.5, 7.51, 7.5200000000000005, 7.53, 7.54, 7.55,
    7.5600000000000005, 7.57, 7.58, 7.59, 7.6000000000000005, 7.61, 7.62, 7.63,
    7.6400000000000006, 7.65, 7.66, 7.67, 7.68, 7.69, 7.7, 7.71, 7.72, 7.73,
    7.74, 7.75, 7.76, 7.7700000000000005, 7.78, 7.79, 7.8, 7.8100000000000005,
    7.82, 7.83, 7.84, 7.8500000000000005, 7.86, 7.87, 7.88, 7.8900000000000006,
    7.9, 7.91, 7.92, 7.9300000000000006, 7.94, 7.95, 7.96, 7.97, 7.98, 7.99, 8.0,
    8.01, 8.02, 8.03, 8.0400000000000009, 8.05, 8.06, 8.07, 8.08, 8.09, 8.1,
    8.11, 8.120000000000001, 8.13, 8.14, 8.15, 8.16, 8.17, 8.18, 8.19, 8.2, 8.21,
    8.22, 8.23, 8.24, 8.25, 8.26, 8.27, 8.28, 8.2900000000000009, 8.3, 8.31,
    8.32, 8.33, 8.34, 8.35, 8.36, 8.370000000000001, 8.38, 8.39, 8.4, 8.41, 8.42,
    8.43, 8.44, 8.45, 8.46, 8.47, 8.48, 8.49, 8.5, 8.51, 8.52, 8.53,
    8.5400000000000009, 8.55, 8.56, 8.57, 8.58, 8.59, 8.6, 8.61,
    8.620000000000001, 8.63, 8.64, 8.65, 8.66, 8.67, 8.68, 8.69,
    8.7000000000000011, 8.71, 8.72, 8.73, 8.74, 8.75, 8.76, 8.77, 8.78,
    8.7900000000000009, 8.8, 8.81, 8.82, 8.83, 8.84, 8.85, 8.86,
    8.870000000000001, 8.88, 8.89, 8.9, 8.91, 8.92, 8.93, 8.94,
    8.9500000000000011, 8.96, 8.97, 8.98, 8.99, 9.0, 9.01, 9.02, 9.03,
    9.0400000000000009, 9.05, 9.06, 9.07, 9.08, 9.09, 9.1, 9.11,
    9.120000000000001, 9.13, 9.14, 9.15, 9.16, 9.17, 9.18, 9.19,
    9.2000000000000011, 9.21, 9.22, 9.23, 9.24, 9.25, 9.26, 9.27, 9.28,
    9.2900000000000009, 9.3, 9.31, 9.32, 9.33, 9.34, 9.35, 9.36,
    9.370000000000001, 9.38, 9.39, 9.4, 9.41, 9.42, 9.43, 9.44,
    9.4500000000000011, 9.46, 9.47, 9.48, 9.49, 9.5, 9.51, 9.52, 9.53,
    9.5400000000000009, 9.55, 9.56, 9.57, 9.58, 9.59, 9.6, 9.61,
    9.620000000000001, 9.63, 9.64, 9.65, 9.66, 9.67, 9.68, 9.69,
    9.7000000000000011, 9.71, 9.72, 9.73, 9.74, 9.75, 9.76, 9.77, 9.78,
    9.7900000000000009, 9.8, 9.81, 9.82, 9.83, 9.84, 9.85, 9.86,
    9.870000000000001, 9.88, 9.89, 9.9, 9.91, 9.92, 9.93, 9.94,
    9.9500000000000011, 9.96, 9.97, 9.98, 9.99, 10.0, 10.01, 10.02, 10.03,
    10.040000000000001, 10.05, 10.06, 10.07, 10.08, 10.09, 10.1, 10.11,
    10.120000000000001, 10.13, 10.14, 10.15, 10.16, 10.17, 10.18, 10.19,
    10.200000000000001, 10.21, 10.22, 10.23, 10.24, 10.25, 10.26, 10.27, 10.28,
    10.290000000000001, 10.3, 10.31, 10.32, 10.33, 10.34, 10.35, 10.36,
    10.370000000000001, 10.38, 10.39, 10.4, 10.41, 10.42, 10.43, 10.44,
    10.450000000000001, 10.46, 10.47, 10.48, 10.49, 10.5, 10.51, 10.52, 10.53,
    10.540000000000001, 10.55, 10.56, 10.57, 10.58, 10.59, 10.6, 10.61,
    10.620000000000001, 10.63, 10.64, 10.65, 10.66, 10.67, 10.68, 10.69,
    10.700000000000001, 10.71, 10.72, 10.73, 10.74, 10.75, 10.76, 10.77, 10.78,
    10.790000000000001, 10.8, 10.81, 10.82, 10.83, 10.84, 10.85, 10.86,
    10.870000000000001, 10.88, 10.89, 10.9, 10.91, 10.92, 10.93, 10.94,
    10.950000000000001, 10.96, 10.97, 10.98, 10.99, 11.0, 11.01, 11.02, 11.03,
    11.04, 11.049999999999999, 11.06, 11.07, 11.08, 11.09, 11.1, 11.11, 11.12,
    11.129999999999999, 11.14, 11.15, 11.16, 11.17, 11.18, 11.19, 11.2,
    11.209999999999999, 11.22, 11.23, 11.24, 11.25, 11.26, 11.27, 11.28, 11.29,
    11.299999999999999, 11.31, 11.32, 11.33, 11.34, 11.35, 11.36, 11.37,
    11.379999999999999, 11.39, 11.4, 11.41, 11.42, 11.43, 11.44, 11.45,
    11.459999999999999, 11.47, 11.48, 11.49, 11.5, 11.51, 11.52, 11.53, 11.54,
    11.549999999999999, 11.56, 11.57, 11.58, 11.59, 11.6, 11.61, 11.62,
    11.629999999999999, 11.64, 11.65, 11.66, 11.67, 11.68, 11.69, 11.7,
    11.709999999999999, 11.72, 11.73, 11.74, 11.75, 11.76, 11.77, 11.78, 11.79,
    11.799999999999999, 11.81, 11.82, 11.83, 11.84, 11.85, 11.86, 11.87,
    11.879999999999999, 11.89, 11.9, 11.91, 11.92, 11.93, 11.94, 11.95,
    11.959999999999999, 11.97, 11.98, 11.99, 12.0, 12.01, 12.02, 12.03, 12.04,
    12.049999999999999, 12.06, 12.07, 12.08, 12.09, 12.1, 12.11, 12.12,
    12.129999999999999, 12.14, 12.15, 12.16, 12.17, 12.18, 12.19, 12.2,
    12.209999999999999, 12.22, 12.23, 12.24, 12.25, 12.26, 12.27, 12.28, 12.29,
    12.299999999999999, 12.31, 12.32, 12.33, 12.34, 12.35, 12.36, 12.37,
    12.379999999999999, 12.39, 12.4, 12.41, 12.42, 12.43, 12.44, 12.45,
    12.459999999999999, 12.47, 12.48, 12.49, 12.5, 12.51, 12.52, 12.53, 12.54,
    12.549999999999999, 12.56, 12.57, 12.58, 12.59, 12.6, 12.61, 12.62,
    12.629999999999999, 12.64, 12.65, 12.66, 12.67, 12.68, 12.69, 12.7,
    12.709999999999999, 12.72, 12.73, 12.74, 12.75, 12.76, 12.77, 12.78, 12.79,
    12.799999999999999, 12.81, 12.82, 12.83, 12.84, 12.85, 12.86, 12.87,
    12.879999999999999, 12.89, 12.9, 12.91, 12.92, 12.93, 12.94, 12.95,
    12.959999999999999, 12.97, 12.98, 12.99, 13.0, 13.01, 13.02, 13.03, 13.04,
    13.049999999999999, 13.06, 13.07, 13.08, 13.09, 13.1, 13.11, 13.12,
    13.129999999999999, 13.14, 13.15, 13.16, 13.17, 13.18, 13.19, 13.2,
    13.209999999999999, 13.22, 13.23, 13.24, 13.25, 13.26, 13.27, 13.28, 13.29,
    13.299999999999999, 13.31, 13.32, 13.33, 13.34, 13.35, 13.36, 13.37,
    13.379999999999999, 13.39, 13.4, 13.41, 13.42, 13.43, 13.44, 13.45,
    13.459999999999999, 13.47, 13.48, 13.49, 13.5, 13.51, 13.52, 13.53, 13.54,
    13.55, 13.56, 13.57, 13.58, 13.59, 13.6, 13.61, 13.62, 13.629999999999999,
    13.64, 13.65, 13.66, 13.67, 13.68, 13.69, 13.7, 13.709999999999999, 13.72,
    13.73, 13.74, 13.75, 13.76, 13.77, 13.78, 13.79, 13.8, 13.81, 13.82, 13.83,
    13.84, 13.85, 13.86, 13.87, 13.879999999999999, 13.89, 13.9, 13.91, 13.92,
    13.93, 13.94, 13.95, 13.959999999999999, 13.97, 13.98, 13.99, 14.0, 14.01,
    14.02, 14.030000000000001, 14.04, 14.05, 14.059999999999999, 14.07, 14.08,
    14.09, 14.1, 14.11, 14.120000000000001, 14.129999999999999, 14.14,
    14.149999999999999, 14.16, 14.17, 14.18, 14.19, 14.2, 14.21,
    14.219999999999999, 14.23, 14.24, 14.25, 14.26, 14.27, 14.280000000000001,
    14.29, 14.3, 14.309999999999999, 14.32, 14.33, 14.34, 14.35, 14.36,
    14.370000000000001, 14.379999999999999, 14.39, 14.399999999999999, 14.41,
    14.42, 14.43, 14.44, 14.45, 14.46, 14.469999999999999, 14.48, 14.49, 14.5,
    14.51, 14.52, 14.530000000000001, 14.54, 14.55, 14.559999999999999, 14.57,
    14.58, 14.59, 14.6, 14.61, 14.620000000000001, 14.629999999999999, 14.64,
    14.649999999999999, 14.66, 14.67, 14.68, 14.69, 14.7, 14.71,
    14.719999999999999, 14.73, 14.74, 14.75, 14.76, 14.77, 14.780000000000001,
    14.79, 14.8, 14.809999999999999, 14.82, 14.83, 14.84, 14.85, 14.86,
    14.870000000000001, 14.879999999999999, 14.89, 14.899999999999999, 14.91,
    14.92, 14.93, 14.94, 14.95, 14.96, 14.969999999999999, 14.98, 14.99, 15.0,
    15.01, 15.02, 15.030000000000001, 15.04, 15.05, 15.059999999999999, 15.07,
    15.08, 15.09, 15.1, 15.11, 15.120000000000001, 15.129999999999999, 15.14,
    15.149999999999999, 15.16, 15.17, 15.18, 15.19, 15.2, 15.21,
    15.219999999999999, 15.23, 15.24, 15.25, 15.26, 15.27, 15.280000000000001,
    15.29, 15.3, 15.309999999999999, 15.32, 15.33, 15.34, 15.35, 15.36,
    15.370000000000001, 15.379999999999999, 15.39, 15.399999999999999, 15.41,
    15.42, 15.43, 15.44, 15.45, 15.46, 15.469999999999999, 15.48, 15.49, 15.5,
    15.51, 15.52, 15.530000000000001, 15.54, 15.55, 15.559999999999999, 15.57,
    15.58, 15.59, 15.6, 15.61, 15.620000000000001, 15.629999999999999, 15.64,
    15.649999999999999, 15.66, 15.67, 15.68, 15.69, 15.7, 15.71,
    15.719999999999999, 15.73, 15.74, 15.75, 15.76, 15.77, 15.780000000000001,
    15.79, 15.8, 15.809999999999999, 15.82, 15.83, 15.84, 15.85, 15.86,
    15.870000000000001, 15.879999999999999, 15.89, 15.899999999999999, 15.91,
    15.92, 15.93, 15.94, 15.95, 15.96, 15.969999999999999, 15.98, 15.99, 16.0,
    16.009999999999998, 16.02, 16.03, 16.04, 16.05, 16.06, 16.07, 16.08, 16.09,
    16.1, 16.11, 16.12, 16.13, 16.14, 16.15, 16.16, 16.17, 16.18,
    16.189999999999998, 16.2, 16.21, 16.22, 16.23, 16.240000000000002, 16.25,
    16.259999999999998, 16.27, 16.28, 16.29, 16.3, 16.31, 16.32, 16.33, 16.34,
    16.35, 16.36, 16.37, 16.38, 16.39, 16.4, 16.41, 16.42, 16.43,
    16.439999999999998, 16.45, 16.46, 16.47, 16.48, 16.490000000000002, 16.5,
    16.509999999999998, 16.52, 16.53, 16.54, 16.55, 16.56, 16.57, 16.58, 16.59,
    16.6, 16.61, 16.62, 16.63, 16.64, 16.65, 16.66, 16.67, 16.68,
    16.689999999999998, 16.7, 16.71, 16.72, 16.73, 16.740000000000002, 16.75,
    16.759999999999998, 16.77, 16.78, 16.79, 16.8, 16.81, 16.82, 16.83, 16.84,
    16.85, 16.86, 16.87, 16.88, 16.89, 16.9, 16.91, 16.92, 16.93,
    16.939999999999998, 16.95, 16.96, 16.97, 16.98, 16.990000000000002, 17.0,
    17.009999999999998, 17.02, 17.03, 17.04, 17.05, 17.06, 17.07, 17.08, 17.09,
    17.1, 17.11, 17.12, 17.13, 17.14, 17.15, 17.16, 17.17, 17.18,
    17.189999999999998, 17.2, 17.21, 17.22, 17.23, 17.240000000000002, 17.25,
    17.259999999999998, 17.27, 17.28, 17.29, 17.3, 17.31, 17.32, 17.33, 17.34,
    17.35, 17.36, 17.37, 17.38, 17.39, 17.4, 17.41, 17.42, 17.43,
    17.439999999999998, 17.45, 17.46, 17.47, 17.48, 17.490000000000002, 17.5,
    17.509999999999998, 17.52, 17.53, 17.54, 17.55, 17.56, 17.57, 17.58, 17.59,
    17.6, 17.61, 17.62, 17.63, 17.64, 17.65, 17.66, 17.67, 17.68,
    17.689999999999998, 17.7, 17.71, 17.72, 17.73, 17.740000000000002, 17.75,
    17.759999999999998, 17.77, 17.78, 17.79, 17.8, 17.81, 17.82, 17.83, 17.84,
    17.85, 17.86, 17.87, 17.88, 17.89, 17.9, 17.91, 17.92, 17.93,
    17.939999999999998, 17.95, 17.96, 17.97, 17.98, 17.990000000000002, 18.0,
    18.009999999999998, 18.02, 18.03, 18.04, 18.05, 18.06, 18.07, 18.08, 18.09,
    18.1, 18.11, 18.12, 18.13, 18.14, 18.15, 18.16, 18.17, 18.18, 18.19, 18.2,
    18.21, 18.22, 18.23, 18.24, 18.25, 18.259999999999998, 18.27, 18.28, 18.29,
    18.3, 18.31, 18.32, 18.33, 18.34, 18.35, 18.36, 18.37, 18.38, 18.39, 18.4,
    18.41, 18.42, 18.43, 18.44, 18.45, 18.46, 18.47, 18.48, 18.49, 18.5,
    18.509999999999998, 18.52, 18.53, 18.54, 18.55, 18.56, 18.57, 18.58, 18.59,
    18.6, 18.61, 18.62, 18.63, 18.64, 18.65, 18.66, 18.67, 18.68, 18.69, 18.7,
    18.71, 18.72, 18.73, 18.74, 18.75, 18.759999999999998, 18.77, 18.78, 18.79,
    18.8, 18.81, 18.82, 18.83, 18.84, 18.85, 18.86, 18.87, 18.88, 18.89, 18.9,
    18.91, 18.92, 18.93, 18.94, 18.95, 18.96, 18.97, 18.98, 18.99, 19.0,
    19.009999999999998, 19.02, 19.03, 19.04, 19.05, 19.06, 19.07, 19.08, 19.09,
    19.1, 19.11, 19.12, 19.13, 19.14, 19.15, 19.16, 19.17, 19.18, 19.19, 19.2,
    19.21, 19.22, 19.23, 19.24, 19.25, 19.259999999999998, 19.27, 19.28, 19.29,
    19.3, 19.31, 19.32, 19.33, 19.34, 19.35, 19.36, 19.37, 19.38, 19.39, 19.4,
    19.41, 19.42, 19.43, 19.44, 19.45, 19.46, 19.47, 19.48, 19.49, 19.5,
    19.509999999999998, 19.52, 19.53, 19.54, 19.55, 19.56, 19.57, 19.58, 19.59,
    19.6, 19.61, 19.62, 19.63, 19.64, 19.65, 19.66, 19.67, 19.68, 19.69, 19.7,
    19.71, 19.72, 19.73, 19.74, 19.75, 19.759999999999998, 19.77, 19.78, 19.79,
    19.8, 19.81, 19.82, 19.83, 19.84, 19.85, 19.86, 19.87, 19.88, 19.89, 19.9,
    19.91, 19.92, 19.93, 19.94, 19.95, 19.96, 19.97, 19.98, 19.99, 20.0, 20.01,
    20.02, 20.03, 20.04, 20.05, 20.06, 20.07, 20.08, 20.09, 20.1, 20.11, 20.12,
    20.13, 20.14, 20.15, 20.16, 20.17, 20.18, 20.19, 20.2, 20.21, 20.22, 20.23,
    20.24, 20.25, 20.26, 20.27, 20.28, 20.29, 20.3, 20.31, 20.32, 20.33, 20.34,
    20.35, 20.36, 20.37, 20.38, 20.39, 20.4, 20.41, 20.42, 20.43, 20.44, 20.45,
    20.46, 20.47, 20.48, 20.49, 20.5, 20.51, 20.52, 20.53, 20.54, 20.55, 20.56,
    20.57, 20.58, 20.59, 20.6, 20.61, 20.62, 20.63, 20.64, 20.65, 20.66, 20.67,
    20.68, 20.69, 20.7, 20.71, 20.72, 20.73, 20.74, 20.75, 20.76, 20.77, 20.78,
    20.79, 20.8, 20.81, 20.82, 20.83, 20.84, 20.85, 20.86, 20.87, 20.88, 20.89,
    20.9, 20.91, 20.92, 20.93, 20.94, 20.95, 20.96, 20.97, 20.98, 20.99, 21.0,
    21.01, 21.02, 21.03, 21.04, 21.05, 21.06, 21.07, 21.08, 21.09, 21.1, 21.11,
    21.12, 21.13, 21.14, 21.15, 21.16, 21.17, 21.18, 21.19, 21.2, 21.21, 21.22,
    21.23, 21.24, 21.25, 21.26, 21.27, 21.28, 21.29, 21.3, 21.31, 21.32, 21.33,
    21.34, 21.35, 21.36, 21.37, 21.38, 21.39, 21.4, 21.41, 21.42, 21.43, 21.44,
    21.45, 21.46, 21.47, 21.48, 21.49, 21.5, 21.51, 21.52, 21.53, 21.54, 21.55,
    21.56, 21.57, 21.58, 21.59, 21.6, 21.61, 21.62, 21.63, 21.64, 21.65, 21.66,
    21.67, 21.68, 21.69, 21.7, 21.71, 21.72, 21.73, 21.74, 21.75, 21.76, 21.77,
    21.78, 21.79, 21.8, 21.81, 21.82, 21.83, 21.84, 21.85, 21.86, 21.87, 21.88,
    21.89, 21.9, 21.91, 21.92, 21.93, 21.94, 21.95, 21.96, 21.97, 21.98, 21.99,
    22.0 };

  double y[2201];
  double z[2201];
  emxArray_real_T *b_x;
  emxArray_real_T *b_y;
  emxArray_real_T *b_z;
  static const double dv63[266] = { 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
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

  double logP0;
  double apnd;
  int n;
  h1 = 0.01;

  /*  This function evaluates the convolution of two inverse Gaussian */
  /*  distributions at vector t. */
  /*  If the variance in one of the distributions is very small so that the  */
  /*  distribution is close the a Dirac delta function, the convolution */
  /*  is approximated as a shifted inverse Gaussian, that is,  */
  /*  one part of the cell cycle is treated as deterministic in length.   */
  /*  In this case, the shift or lag is the average time to complete  */
  /*  the part of the cycle that has the smallest standard deviation. */
  /*  t is a vector of times to divide (or times to complete two parts of the */
  /*  cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1, */
  /*  s2=sigma2. */
  /*  Input parameters: */
  /*  t = the list of points at which to evaluate the distribution */
  /*  m1 = mu for the first distribution */
  /*  s1 = sigma for the first distribution */
  /*  m2 = mu for the second distribution */
  /*  s2 = sigma for the second distribution */
  /*  h = step size */
  /*  Outputs: */
  /*  P = probability at each point corresponding to points in t */
  /*  h = final step size */
  /*  flag = indicates if we used the Dirac delta */
  /*  E = is the relative error in the likelihood of the data due to the numerical integration */
  /*  log the parameters we were called with */
  /* fprintf(1,'convolv_2invG(m1=%.17g, s1=%.17g, m2=%.17g, s2=%.17g)\n', ...  */
  /*     m1, s1, m2, s2); */
  /* fprintf('    convolv_2invG_adapt(m1=%f s1=%f m2=%f s2=%f)\n',m1,s1,m2,s2); */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:37' flag=0; */
  /*  a constant that is used in determining if the Dirac approximation should be applied. */
  /* 'convolv_2invG_adapt_nov:40' eps=.01; */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:43' E=Inf; */
  E = rtInf;

  /*  number of points to be evaluated */
  /* 'convolv_2invG_adapt_nov:46' n=length(t); */
  /*  store m and s for each sub-distribution in a list */
  /* 'convolv_2invG_adapt_nov:49' m=[m1 m2]; */
  /* 'convolv_2invG_adapt_nov:50' s=[s1 s2]; */
  /*  make sure we dont have negative values for m or s by replacement with */
  /*  zero if we do */
  /* m=max(m,0); */
  /* s=max(s,0); */
  /*  reflect negative parameters into the positive space so we can use */
  /*  unconstrained optimization */
  /* 'convolv_2invG_adapt_nov:58' m=abs(m); */
  x[0] = m1;
  x[1] = m2;
  for (k = 0; k < 2; k++) {
    m[k] = fabs(x[k]);
  }

  /* 'convolv_2invG_adapt_nov:59' s=abs(s); */
  x[0] = s1;
  x[1] = s2;

  /*  find the variance for both sub-distributions */
  /* 'convolv_2invG_adapt_nov:62' v=(s.^2)./(m.^3); */
  for (k = 0; k < 2; k++) {
    ndbl = fabs(x[k]);
    v[k] = rt_powd_snf(ndbl, 2.0);
    s[k] = ndbl;
  }

  /*  find the standard deviation for each sub-distribution */
  /* 'convolv_2invG_adapt_nov:65' sd=v.^.5; */
  for (k = 0; k < 2; k++) {
    ndbl = v[k] / rt_powd_snf(m[k], 3.0);
    sd[k] = rt_powd_snf(ndbl, 0.5);
    v[k] = ndbl;
  }

  /*  reorder m and s so that the sub-distribution with the smallest */
  /*  variance comes first.  So the fist part might be approximated as a Dirac delta. */
  /* 'convolv_2invG_adapt_nov:69' [v,I]=sort(v); */
  sort(v, iidx);

  /* 'convolv_2invG_adapt_nov:70' m=m(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    v[nm1d2] = m[iidx[nm1d2] - 1];
    x[nm1d2] = iidx[nm1d2];
  }

  /* 'convolv_2invG_adapt_nov:71' s=s(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    m[nm1d2] = v[nm1d2];
    b_s[nm1d2] = s[(int)x[nm1d2] - 1];
  }

  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    s[nm1d2] = b_s[nm1d2];
  }

  /*  T1 appears to be unused, should probably remove */
  /*  T2 is only used inside FUNCTION CHECK_APPROXIMATABLE */
  /*  and should probably be moved inside that scope */
  /*  T2 is the mode of the second distribution, determined by a precomputed analytic expression. */
  /* 'convolv_2invG_adapt_nov:78' T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1))); */
  /* 'convolv_2invG_adapt_nov:79' T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2))); */
  /* When called from convolv_3invG all of the components of t may be negative. */
  /* In this case, the probability of each is set to realmin */
  /* 'convolv_2invG_adapt_nov:83' if max(t)<=0 */
  /* 'convolv_2invG_adapt_nov:87' else */
  /*  find the largest point in t */
  /* 'convolv_2invG_adapt_nov:89' Maxt=max(t); */
  /*  produce a range of evenly spaced points to evaluate at between */
  /*  0 and Maxt with step size h. the even spacing is important when */
  /*  calculating the convolution later */
  /* 'convolv_2invG_adapt_nov:94' x=0:h:Maxt; */
  /*  evaluate the sub-distributions at each point in x */
  /* 'convolv_2invG_adapt_nov:97' y=onestagepdf2(x,m(1),s(1)); */
  b_onestagepdf2(dv62, m[0], s[0], y);

  /* 'convolv_2invG_adapt_nov:98' z=onestagepdf2(x,m(2),s(2)); */
  b_onestagepdf2(dv62, m[1], s[1], z);

  /*  if the first pdf is very concentrated, check to see if it can be */
  /*  approximated as a point-mass distribution */
  /* 'convolv_2invG_adapt_nov:102' if sd(1)<.01 */
  emxInit_real_T(&b_x, 2);
  emxInit_real_T(&b_y, 2);
  emxInit_real_T(&b_z, 2);
  if (sd[0] < 0.01) {
    /* 'convolv_2invG_adapt_nov:104' check2 = tailmass(m, s, eps, T2, sd); */
    ndbl = tailmass(m, s, 1.0 / m[1] * (sqrt(1.0 + 2.25 * (rt_powd_snf(s[1], 4.0)
      / (m[1] * m[1]))) - 1.5 * (s[1] * s[1] / m[1])), sd);

    /* 'convolv_2invG_adapt_nov:106' if  check2<=eps/3 */
    if (ndbl <= 0.0033333333333333335) {
      /* If there is not much probability in the tails of the convolution, we use a Dirac Delta for part 1. */
      /* Flag is set to 1 to indicate this. */
      /* 'convolv_2invG_adapt_nov:109' flag=1; */
      /* 'convolv_2invG_adapt_nov:110' sigma=s(2); */
      /* 'convolv_2invG_adapt_nov:111' mu=m(2); */
      /* l for lag. */
      /* 'convolv_2invG_adapt_nov:113' l=1/m(1); */
      /* 'convolv_2invG_adapt_nov:114' P=onestagepdf_lag(t,mu,sigma,l); */
      onestagepdf_lag(dv63, m[1], s[1], 1.0 / m[0], P);
    } else {
      /* 'convolv_2invG_adapt_nov:115' else */
      /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
      /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
      /*  Outputs: P */
      /* 'convolv_2invG_adapt_nov:120' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
      approxconvolv(z, y, dv63, dv62, P, &logP0);

      /* Keep reducing the step size in the numerical integration until we are happy with the error. */
      /* 'convolv_2invG_adapt_nov:124' while E>=.001*abs(logP0) */
      while (E >= 0.001 * fabs(logP0)) {
        /* 'convolv_2invG_adapt_nov:125' h1=.5*h; */
        h1 *= 0.5;

#ifdef _VERBOSE
		printf("h=%f logP0=%f ", h1, logP0);
#endif

        /* 'convolv_2invG_adapt_nov:126' x=0:h1:Maxt; */
        if ((h1 == 0.0) || (h1 < 0.0)) {
          nm1d2 = b_x->size[0] * b_x->size[1];
          b_x->size[0] = 1;
          b_x->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
        } else if (0.0 == h1) {
          nm1d2 = b_x->size[0] * b_x->size[1];
          b_x->size[0] = 1;
          b_x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
          for (nm1d2 = 0; nm1d2 < 1; nm1d2++) {
            b_x->data[b_x->size[0] * nm1d2] = 0.0 * (double)nm1d2;
          }
        } else {
          ndbl = floor(22.0 / h1 + 0.5);
          apnd = ndbl * h1;
          if (h1 > 0.0) {
            E = apnd - 22.0;
          } else {
            E = 22.0 - apnd;
          }

          if (fabs(E) < 9.7699626167013776E-15) {
            ndbl++;
            apnd = 22.0;
          } else if (E > 0.0) {
            apnd = (ndbl - 1.0) * h1;
          } else {
            ndbl++;
          }

          if (ndbl >= 0.0) {
            n = (int)ndbl;
          } else {
            n = 0;
          }

          nm1d2 = b_x->size[0] * b_x->size[1];
          b_x->size[0] = 1;
          b_x->size[1] = n;
          emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
          if (n > 0) {
            b_x->data[0] = 0.0;
            if (n > 1) {
              b_x->data[n - 1] = apnd;
              nm1d2 = (n - 1) / 2;
              for (k = 1; k < nm1d2; k++) {
                ndbl = (double)k * h1;
                b_x->data[k] = ndbl;
                b_x->data[(n - k) - 1] = apnd - ndbl;
              }

              if (nm1d2 << 1 == n - 1) {
                b_x->data[nm1d2] = apnd / 2.0;
              } else {
                ndbl = (double)nm1d2 * h1;
                b_x->data[nm1d2] = ndbl;
                b_x->data[nm1d2 + 1] = apnd - ndbl;
              }
            }
          }
        }

        /* 'convolv_2invG_adapt_nov:127' y=onestagepdf2(x,m(1),s(1)); */
        d_onestagepdf2(b_x, m[0], s[0], b_y);

        /* 'convolv_2invG_adapt_nov:128' z=onestagepdf2(x,m(2),s(2)); */
        d_onestagepdf2(b_x, m[1], s[1], b_z);

        /* 'convolv_2invG_adapt_nov:130' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
        b_approxconvolv(b_z, b_y, h1, b_x, P, &ndbl);

        /* 'convolv_2invG_adapt_nov:132' E=abs(logP1-logP0); */
        E = fabs(ndbl - logP0);

#ifdef _VERBOSE
		printf("logP1=%f E=%f\n", ndbl, E);
#endif

        /* 'convolv_2invG_adapt_nov:133' P0=P1; */
        /* 'convolv_2invG_adapt_nov:134' logP0=logP1; */
        logP0 = ndbl;

        /* 'convolv_2invG_adapt_nov:135' h=h1; */
      }

      /* 'convolv_2invG_adapt_nov:137' P=P0; */
      /*  END FUNCTION DOTHECONVOLUTION_OUTER */
    }

    /*  pdf is not very concentrated so compute the convolution directly */
  } else {
    /* 'convolv_2invG_adapt_nov:141' else */
    /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
    /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
    /*  Outputs: P */
    /* 'convolv_2invG_adapt_nov:146' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
    approxconvolv(z, y, dv63, dv62, P, &logP0);

    /* 'convolv_2invG_adapt_nov:148' while E>=.001*abs(logP0) */
    while (E >= 0.001 * fabs(logP0)) {
      /* 'convolv_2invG_adapt_nov:149' h1=.5*h; */
      h1 *= 0.5;

#ifdef _VERBOSE
	  printf("h=%f logP0=%f ", h1, logP0);
#endif


      /* 'convolv_2invG_adapt_nov:150' x=0:h1:Maxt; */
      if ((h1 == 0.0) || (h1 < 0.0)) {
        nm1d2 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
      } else if (0.0 == h1) {
        nm1d2 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
        for (nm1d2 = 0; nm1d2 < 1; nm1d2++) {
          b_x->data[b_x->size[0] * nm1d2] = 0.0 * (double)nm1d2;
        }
      } else {
        ndbl = floor(22.0 / h1 + 0.5);
        apnd = ndbl * h1;
        if (h1 > 0.0) {
          E = apnd - 22.0;
        } else {
          E = 22.0 - apnd;
        }

        if (fabs(E) < 9.7699626167013776E-15) {
          ndbl++;
          apnd = 22.0;
        } else if (E > 0.0) {
          apnd = (ndbl - 1.0) * h1;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          n = (int)ndbl;
        } else {
          n = 0;
        }

        nm1d2 = b_x->size[0] * b_x->size[1];
        b_x->size[0] = 1;
        b_x->size[1] = n;
        emxEnsureCapacity((emxArray__common *)b_x, nm1d2, sizeof(double));
        if (n > 0) {
          b_x->data[0] = 0.0;
          if (n > 1) {
            b_x->data[n - 1] = apnd;
            nm1d2 = (n - 1) / 2;
            for (k = 1; k < nm1d2; k++) {
              ndbl = (double)k * h1;
              b_x->data[k] = ndbl;
              b_x->data[(n - k) - 1] = apnd - ndbl;
            }

            if (nm1d2 << 1 == n - 1) {
              b_x->data[nm1d2] = apnd / 2.0;
            } else {
              ndbl = (double)nm1d2 * h1;
              b_x->data[nm1d2] = ndbl;
              b_x->data[nm1d2 + 1] = apnd - ndbl;
            }
          }
        }
      }

      /* 'convolv_2invG_adapt_nov:151' y=onestagepdf2(x,m(1),s(1)); */
      d_onestagepdf2(b_x, m[0], s[0], b_y);

      /* 'convolv_2invG_adapt_nov:152' z=onestagepdf2(x,m(2),s(2)); */
      d_onestagepdf2(b_x, m[1], s[1], b_z);

      /* 'convolv_2invG_adapt_nov:154' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
      b_approxconvolv(b_z, b_y, h1, b_x, P, &ndbl);

      /* 'convolv_2invG_adapt_nov:157' E=abs(logP1-logP0); */
      E = fabs(ndbl - logP0);

#ifdef _VERBOSE
	  printf("logP1=%f E=%f\n", ndbl, E);
#endif

      /* 'convolv_2invG_adapt_nov:158' P0=P1; */
      /* 'convolv_2invG_adapt_nov:159' logP0=logP1; */
      logP0 = ndbl;

      /* 'convolv_2invG_adapt_nov:160' h=h1; */
    }

    /* 'convolv_2invG_adapt_nov:162' P=P0; */
    /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
  }

  emxFree_real_T(&b_z);
  emxFree_real_T(&b_y);
  emxFree_real_T(&b_x);
}
#endif

/*
 * function [P,h,flag,E]=convolv_2invG_adapt_nov(t,m1,s1,m2,s2,h)
 */
#ifdef _OLD_MATLAB_CODE
void d_convolv_2invG_adapt_nov(const double t[266], double m1, double s1, double
  m2, double s2, double P[266])
{
  double h1;
  double E;
  double v[2];
  double m[2];
  double s[2];
  double I[2];
  int ix;
  double sd[2];
  int iidx[2];
  double b_s[2];
  int ixstart;
  double mtmp;
  boolean_T exitg1;
  emxArray_real_T *x;
  double ndbl;
  double apnd;
  double cdiff;
  emxArray_real_T *y;
  emxArray_real_T *z;
  int n;
  double logP0;
  h1 = 0.1;

  /*  This function evaluates the convolution of two inverse Gaussian */
  /*  distributions at vector t. */
  /*  If the variance in one of the distributions is very small so that the  */
  /*  distribution is close the a Dirac delta function, the convolution */
  /*  is approximated as a shifted inverse Gaussian, that is,  */
  /*  one part of the cell cycle is treated as deterministic in length.   */
  /*  In this case, the shift or lag is the average time to complete  */
  /*  the part of the cycle that has the smallest standard deviation. */
  /*  t is a vector of times to divide (or times to complete two parts of the */
  /*  cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1, */
  /*  s2=sigma2. */
  /*  Input parameters: */
  /*  t = the list of points at which to evaluate the distribution */
  /*  m1 = mu for the first distribution */
  /*  s1 = sigma for the first distribution */
  /*  m2 = mu for the second distribution */
  /*  s2 = sigma for the second distribution */
  /*  h = step size */
  /*  Outputs: */
  /*  P = probability at each point corresponding to points in t */
  /*  h = final step size */
  /*  flag = indicates if we used the Dirac delta */
  /*  E = is the relative error in the likelihood of the data due to the numerical integration */
  /*  log the parameters we were called with */
  /* fprintf(1,'convolv_2invG(m1=%.17g, s1=%.17g, m2=%.17g, s2=%.17g)\n', ...  */
  /*     m1, s1, m2, s2); */
  /* fprintf('    convolv_2invG_adapt(m1=%f s1=%f m2=%f s2=%f)\n',m1,s1,m2,s2); */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:37' flag=0; */
  /*  a constant that is used in determining if the Dirac approximation should be applied. */
  /* 'convolv_2invG_adapt_nov:40' eps=.01; */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:43' E=Inf; */
  E = rtInf;

  /*  number of points to be evaluated */
  /* 'convolv_2invG_adapt_nov:46' n=length(t); */
  /*  store m and s for each sub-distribution in a list */
  /* 'convolv_2invG_adapt_nov:49' m=[m1 m2]; */
  /* 'convolv_2invG_adapt_nov:50' s=[s1 s2]; */
  /*  make sure we dont have negative values for m or s by replacement with */
  /*  zero if we do */
  /* m=max(m,0); */
  /* s=max(s,0); */
  /*  reflect negative parameters into the positive space so we can use */
  /*  unconstrained optimization */
  /* 'convolv_2invG_adapt_nov:58' m=abs(m); */
  v[0] = m1;
  v[1] = m2;
  b_abs(v, m);

  /* 'convolv_2invG_adapt_nov:59' s=abs(s); */
  v[0] = s1;
  v[1] = s2;
  b_abs(v, s);

  /*  find the variance for both sub-distributions */
  /* 'convolv_2invG_adapt_nov:62' v=(s.^2)./(m.^3); */
  power(s, v);
  b_power(m, I);
  for (ix = 0; ix < 2; ix++) {
    v[ix] /= I[ix];
  }

  /*  find the standard deviation for each sub-distribution */
  /* 'convolv_2invG_adapt_nov:65' sd=v.^.5; */
  c_power(v, sd);

  /*  reorder m and s so that the sub-distribution with the smallest */
  /*  variance comes first.  So the fist part might be approximated as a Dirac delta. */
  /* 'convolv_2invG_adapt_nov:69' [v,I]=sort(v); */
  sort(v, iidx);

  /* 'convolv_2invG_adapt_nov:70' m=m(I); */
  for (ix = 0; ix < 2; ix++) {
    v[ix] = m[iidx[ix] - 1];
    I[ix] = iidx[ix];
  }

  /* 'convolv_2invG_adapt_nov:71' s=s(I); */
  for (ix = 0; ix < 2; ix++) {
    m[ix] = v[ix];
    b_s[ix] = s[(int)I[ix] - 1];
  }

  for (ix = 0; ix < 2; ix++) {
    s[ix] = b_s[ix];
  }

  /*  T1 appears to be unused, should probably remove */
  /*  T2 is only used inside FUNCTION CHECK_APPROXIMATABLE */
  /*  and should probably be moved inside that scope */
  /*  T2 is the mode of the second distribution, determined by a precomputed analytic expression. */
  /* 'convolv_2invG_adapt_nov:78' T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1))); */
  /* 'convolv_2invG_adapt_nov:79' T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2))); */
  /* When called from convolv_3invG all of the components of t may be negative. */
  /* In this case, the probability of each is set to realmin */
  /* 'convolv_2invG_adapt_nov:83' if max(t)<=0 */
  ixstart = 1;
  mtmp = t[0];
  if (rtIsNaN(t[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 267)) {
      ixstart = ix;
      if (!rtIsNaN(t[ix - 1])) {
        mtmp = t[ix - 1];
        exitg1 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 266) {
    while (ixstart + 1 < 267) {
      if (t[ixstart] > mtmp) {
        mtmp = t[ixstart];
      }

      ixstart++;
    }
  }

  if (mtmp <= 0.0) {
    /* 'convolv_2invG_adapt_nov:84' P=realmin.*ones(size(t)); */
    for (ixstart = 0; ixstart < 266; ixstart++) {
      P[ixstart] = 2.2250738585072014E-308;
    }

    /*  otherwise we need to calculate P for each point in t */
  } else {
    /* 'convolv_2invG_adapt_nov:87' else */
    /*  find the largest point in t */
    /* 'convolv_2invG_adapt_nov:89' Maxt=max(t); */
    ixstart = 1;
    mtmp = t[0];
    if (rtIsNaN(t[0])) {
      ix = 2;
      exitg1 = false;
      while ((!exitg1) && (ix < 267)) {
        ixstart = ix;
        if (!rtIsNaN(t[ix - 1])) {
          mtmp = t[ix - 1];
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < 266) {
      while (ixstart + 1 < 267) {
        if (t[ixstart] > mtmp) {
          mtmp = t[ixstart];
        }

        ixstart++;
      }
    }

    /*  produce a range of evenly spaced points to evaluate at between */
    /*  0 and Maxt with step size h. the even spacing is important when */
    /*  calculating the convolution later */
    /* 'convolv_2invG_adapt_nov:94' x=0:h:Maxt; */
    emxInit_real_T(&x, 2);
    if (rtIsNaN(mtmp)) {
      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
      x->data[0] = rtNaN;
    } else if (mtmp < 0.0) {
      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
    } else if (rtIsInf(mtmp) && (0.0 == mtmp)) {
      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
      x->data[0] = rtNaN;
    } else {
      ndbl = floor(mtmp / 0.1 + 0.5);
      apnd = ndbl * 0.1;
      cdiff = apnd - mtmp;
      if (fabs(cdiff) < 4.4408920985006262E-16 * mtmp) {
        ndbl++;
        apnd = mtmp;
      } else if (cdiff > 0.0) {
        apnd = (ndbl - 1.0) * 0.1;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int)ndbl;
      } else {
        n = 0;
      }

      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = n;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
      if (n > 0) {
        x->data[0] = 0.0;
        if (n > 1) {
          x->data[n - 1] = apnd;
          ixstart = (n - 1) / 2;
          for (ix = 1; ix < ixstart; ix++) {
            ndbl = (double)ix * 0.1;
            x->data[ix] = ndbl;
            x->data[(n - ix) - 1] = apnd - ndbl;
          }

          if (ixstart << 1 == n - 1) {
            x->data[ixstart] = apnd / 2.0;
          } else {
            ndbl = (double)ixstart * 0.1;
            x->data[ixstart] = ndbl;
            x->data[ixstart + 1] = apnd - ndbl;
          }
        }
      }
    }

    emxInit_real_T(&y, 2);
    emxInit_real_T(&z, 2);

    /*  evaluate the sub-distributions at each point in x */
    /* 'convolv_2invG_adapt_nov:97' y=onestagepdf2(x,m(1),s(1)); */
    d_onestagepdf2(x, m[0], s[0], y);

    /* 'convolv_2invG_adapt_nov:98' z=onestagepdf2(x,m(2),s(2)); */
    d_onestagepdf2(x, m[1], s[1], z);

    /*  if the first pdf is very concentrated, check to see if it can be */
    /*  approximated as a point-mass distribution */
    /* 'convolv_2invG_adapt_nov:102' if sd(1)<.01 */
    if (sd[0] < 0.01) {
      /* 'convolv_2invG_adapt_nov:104' check2 = tailmass(m, s, eps, T2, sd); */
      ndbl = tailmass(m, s, 1.0 / m[1] * (sqrt(1.0 + 2.25 * (rt_powd_snf(s[1],
        4.0) / (m[1] * m[1]))) - 1.5 * (s[1] * s[1] / m[1])), sd);

      /* 'convolv_2invG_adapt_nov:106' if  check2<=eps/3 */
      if (ndbl <= 0.0033333333333333335) {
        /* If there is not much probability in the tails of the convolution, we use a Dirac Delta for part 1. */
        /* Flag is set to 1 to indicate this. */
        /* 'convolv_2invG_adapt_nov:109' flag=1; */
        /* 'convolv_2invG_adapt_nov:110' sigma=s(2); */
        /* 'convolv_2invG_adapt_nov:111' mu=m(2); */
        /* l for lag. */
        /* 'convolv_2invG_adapt_nov:113' l=1/m(1); */
        /* 'convolv_2invG_adapt_nov:114' P=onestagepdf_lag(t,mu,sigma,l); */
        onestagepdf_lag(t, m[1], s[1], 1.0 / m[0], P);
      } else {
        /* 'convolv_2invG_adapt_nov:115' else */
        /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
        /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
        /*  Outputs: P */
        /* 'convolv_2invG_adapt_nov:120' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
        d_approxconvolv(z, y, t, x, P, &logP0);

        /* Keep reducing the step size in the numerical integration until we are happy with the error. */
        /* 'convolv_2invG_adapt_nov:124' while E>=.001*abs(logP0) */
        while (E >= 0.001 * fabs(logP0)) {
          /* 'convolv_2invG_adapt_nov:125' h1=.5*h; */
          h1 *= 0.5;

          /* 'convolv_2invG_adapt_nov:126' x=0:h1:Maxt; */
          if (rtIsNaN(mtmp)) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            x->data[0] = rtNaN;
          } else if ((h1 == 0.0) || ((0.0 < mtmp) && (h1 < 0.0)) || ((mtmp < 0.0)
                      && (h1 > 0.0))) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 0;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          } else if (rtIsInf(mtmp) && (rtIsInf(h1) || (0.0 == mtmp))) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            x->data[0] = rtNaN;
          } else if (rtIsInf(h1)) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            x->data[0] = 0.0;
          } else if (floor(h1) == h1) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = (int)floor(mtmp / h1) + 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            ixstart = (int)floor(mtmp / h1);
            for (ix = 0; ix <= ixstart; ix++) {
              x->data[x->size[0] * ix] = h1 * (double)ix;
            }
          } else {
            ndbl = floor(mtmp / h1 + 0.5);
            apnd = ndbl * h1;
            if (h1 > 0.0) {
              cdiff = apnd - mtmp;
            } else {
              cdiff = mtmp - apnd;
            }

            if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(mtmp)) {
              ndbl++;
              apnd = mtmp;
            } else if (cdiff > 0.0) {
              apnd = (ndbl - 1.0) * h1;
            } else {
              ndbl++;
            }

            if (ndbl >= 0.0) {
              n = (int)ndbl;
            } else {
              n = 0;
            }

            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = n;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            if (n > 0) {
              x->data[0] = 0.0;
              if (n > 1) {
                x->data[n - 1] = apnd;
                ixstart = (n - 1) / 2;
                for (ix = 1; ix < ixstart; ix++) {
                  ndbl = (double)ix * h1;
                  x->data[ix] = ndbl;
                  x->data[(n - ix) - 1] = apnd - ndbl;
                }

                if (ixstart << 1 == n - 1) {
                  x->data[ixstart] = apnd / 2.0;
                } else {
                  ndbl = (double)ixstart * h1;
                  x->data[ixstart] = ndbl;
                  x->data[ixstart + 1] = apnd - ndbl;
                }
              }
            }
          }

          /* 'convolv_2invG_adapt_nov:127' y=onestagepdf2(x,m(1),s(1)); */
          d_onestagepdf2(x, m[0], s[0], y);

          /* 'convolv_2invG_adapt_nov:128' z=onestagepdf2(x,m(2),s(2)); */
          d_onestagepdf2(x, m[1], s[1], z);

          /* 'convolv_2invG_adapt_nov:130' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
          e_approxconvolv(z, y, h1, t, x, P, &ndbl);

          /* 'convolv_2invG_adapt_nov:132' E=abs(logP1-logP0); */
          E = fabs(ndbl - logP0);

          /* 'convolv_2invG_adapt_nov:133' P0=P1; */
          /* 'convolv_2invG_adapt_nov:134' logP0=logP1; */
          logP0 = ndbl;

          /* 'convolv_2invG_adapt_nov:135' h=h1; */
        }

        /* 'convolv_2invG_adapt_nov:137' P=P0; */
        /*  END FUNCTION DOTHECONVOLUTION_OUTER */
      }

      /*  pdf is not very concentrated so compute the convolution directly */
    } else {
      /* 'convolv_2invG_adapt_nov:141' else */
      /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
      /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
      /*  Outputs: P */
      /* 'convolv_2invG_adapt_nov:146' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
      d_approxconvolv(z, y, t, x, P, &logP0);

      /* 'convolv_2invG_adapt_nov:148' while E>=.001*abs(logP0) */
      while (E >= 0.001 * fabs(logP0)) {
        /* 'convolv_2invG_adapt_nov:149' h1=.5*h; */
        h1 *= 0.5;

        /* 'convolv_2invG_adapt_nov:150' x=0:h1:Maxt; */
        if (rtIsNaN(mtmp)) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          x->data[0] = rtNaN;
        } else if ((h1 == 0.0) || ((0.0 < mtmp) && (h1 < 0.0)) || ((mtmp < 0.0) &&
                    (h1 > 0.0))) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
        } else if (rtIsInf(mtmp) && (rtIsInf(h1) || (0.0 == mtmp))) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          x->data[0] = rtNaN;
        } else if (rtIsInf(h1)) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          x->data[0] = 0.0;
        } else if (floor(h1) == h1) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = (int)floor(mtmp / h1) + 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          ixstart = (int)floor(mtmp / h1);
          for (ix = 0; ix <= ixstart; ix++) {
            x->data[x->size[0] * ix] = h1 * (double)ix;
          }
        } else {
          ndbl = floor(mtmp / h1 + 0.5);
          apnd = ndbl * h1;
          if (h1 > 0.0) {
            cdiff = apnd - mtmp;
          } else {
            cdiff = mtmp - apnd;
          }

          if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(mtmp)) {
            ndbl++;
            apnd = mtmp;
          } else if (cdiff > 0.0) {
            apnd = (ndbl - 1.0) * h1;
          } else {
            ndbl++;
          }

          if (ndbl >= 0.0) {
            n = (int)ndbl;
          } else {
            n = 0;
          }

          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = n;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          if (n > 0) {
            x->data[0] = 0.0;
            if (n > 1) {
              x->data[n - 1] = apnd;
              ixstart = (n - 1) / 2;
              for (ix = 1; ix < ixstart; ix++) {
                ndbl = (double)ix * h1;
                x->data[ix] = ndbl;
                x->data[(n - ix) - 1] = apnd - ndbl;
              }

              if (ixstart << 1 == n - 1) {
                x->data[ixstart] = apnd / 2.0;
              } else {
                ndbl = (double)ixstart * h1;
                x->data[ixstart] = ndbl;
                x->data[ixstart + 1] = apnd - ndbl;
              }
            }
          }
        }

        /* 'convolv_2invG_adapt_nov:151' y=onestagepdf2(x,m(1),s(1)); */
        d_onestagepdf2(x, m[0], s[0], y);

        /* 'convolv_2invG_adapt_nov:152' z=onestagepdf2(x,m(2),s(2)); */
        d_onestagepdf2(x, m[1], s[1], z);

        /* 'convolv_2invG_adapt_nov:154' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
        e_approxconvolv(z, y, h1, t, x, P, &ndbl);

        /* 'convolv_2invG_adapt_nov:157' E=abs(logP1-logP0); */
        E = fabs(ndbl - logP0);

        /* 'convolv_2invG_adapt_nov:158' P0=P1; */
        /* 'convolv_2invG_adapt_nov:159' logP0=logP1; */
        logP0 = ndbl;

        /* 'convolv_2invG_adapt_nov:160' h=h1; */
      }

      /* 'convolv_2invG_adapt_nov:162' P=P0; */
      /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
    }

    emxFree_real_T(&z);
    emxFree_real_T(&y);
    emxFree_real_T(&x);
  }
}
#endif

/*
 * function [P,h,flag,E]=convolv_2invG_adapt_nov(t,m1,s1,m2,s2,h)
 */
#ifdef _OLD_MATLAB_CODE
void e_convolv_2invG_adapt_nov(const double t[266], double m1, double s1, double
  m2, double s2, double P[266])
{
  double h1;
  double E;
  double v[2];
  double m[2];
  double s[2];
  double I[2];
  int ix;
  double sd[2];
  int iidx[2];
  double b_s[2];
  int ixstart;
  double mtmp;
  boolean_T exitg1;
  emxArray_real_T *x;
  double ndbl;
  double apnd;
  double cdiff;
  emxArray_real_T *y;
  emxArray_real_T *z;
  int n;
  double logP0;
  h1 = 0.01;

  /*  This function evaluates the convolution of two inverse Gaussian */
  /*  distributions at vector t. */
  /*  If the variance in one of the distributions is very small so that the  */
  /*  distribution is close the a Dirac delta function, the convolution */
  /*  is approximated as a shifted inverse Gaussian, that is,  */
  /*  one part of the cell cycle is treated as deterministic in length.   */
  /*  In this case, the shift or lag is the average time to complete  */
  /*  the part of the cycle that has the smallest standard deviation. */
  /*  t is a vector of times to divide (or times to complete two parts of the */
  /*  cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1, */
  /*  s2=sigma2. */
  /*  Input parameters: */
  /*  t = the list of points at which to evaluate the distribution */
  /*  m1 = mu for the first distribution */
  /*  s1 = sigma for the first distribution */
  /*  m2 = mu for the second distribution */
  /*  s2 = sigma for the second distribution */
  /*  h = step size */
  /*  Outputs: */
  /*  P = probability at each point corresponding to points in t */
  /*  h = final step size */
  /*  flag = indicates if we used the Dirac delta */
  /*  E = is the relative error in the likelihood of the data due to the numerical integration */
  /*  log the parameters we were called with */
  /* fprintf(1,'convolv_2invG(m1=%.17g, s1=%.17g, m2=%.17g, s2=%.17g)\n', ...  */
  /*     m1, s1, m2, s2); */
  /* fprintf('    convolv_2invG_adapt(m1=%f s1=%f m2=%f s2=%f)\n',m1,s1,m2,s2); */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:37' flag=0; */
  /*  a constant that is used in determining if the Dirac approximation should be applied. */
  /* 'convolv_2invG_adapt_nov:40' eps=.01; */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:43' E=Inf; */
  E = rtInf;

  /*  number of points to be evaluated */
  /* 'convolv_2invG_adapt_nov:46' n=length(t); */
  /*  store m and s for each sub-distribution in a list */
  /* 'convolv_2invG_adapt_nov:49' m=[m1 m2]; */
  /* 'convolv_2invG_adapt_nov:50' s=[s1 s2]; */
  /*  make sure we dont have negative values for m or s by replacement with */
  /*  zero if we do */
  /* m=max(m,0); */
  /* s=max(s,0); */
  /*  reflect negative parameters into the positive space so we can use */
  /*  unconstrained optimization */
  /* 'convolv_2invG_adapt_nov:58' m=abs(m); */
  v[0] = m1;
  v[1] = m2;
  b_abs(v, m);

  /* 'convolv_2invG_adapt_nov:59' s=abs(s); */
  v[0] = s1;
  v[1] = s2;
  b_abs(v, s);

  /*  find the variance for both sub-distributions */
  /* 'convolv_2invG_adapt_nov:62' v=(s.^2)./(m.^3); */
  power(s, v);
  b_power(m, I);
  for (ix = 0; ix < 2; ix++) {
    v[ix] /= I[ix];
  }

  /*  find the standard deviation for each sub-distribution */
  /* 'convolv_2invG_adapt_nov:65' sd=v.^.5; */
  c_power(v, sd);

  /*  reorder m and s so that the sub-distribution with the smallest */
  /*  variance comes first.  So the fist part might be approximated as a Dirac delta. */
  /* 'convolv_2invG_adapt_nov:69' [v,I]=sort(v); */
  sort(v, iidx);

  /* 'convolv_2invG_adapt_nov:70' m=m(I); */
  for (ix = 0; ix < 2; ix++) {
    v[ix] = m[iidx[ix] - 1];
    I[ix] = iidx[ix];
  }

  /* 'convolv_2invG_adapt_nov:71' s=s(I); */
  for (ix = 0; ix < 2; ix++) {
    m[ix] = v[ix];
    b_s[ix] = s[(int)I[ix] - 1];
  }

  for (ix = 0; ix < 2; ix++) {
    s[ix] = b_s[ix];
  }

  /*  T1 appears to be unused, should probably remove */
  /*  T2 is only used inside FUNCTION CHECK_APPROXIMATABLE */
  /*  and should probably be moved inside that scope */
  /*  T2 is the mode of the second distribution, determined by a precomputed analytic expression. */
  /* 'convolv_2invG_adapt_nov:78' T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1))); */
  /* 'convolv_2invG_adapt_nov:79' T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2))); */
  /* When called from convolv_3invG all of the components of t may be negative. */
  /* In this case, the probability of each is set to realmin */
  /* 'convolv_2invG_adapt_nov:83' if max(t)<=0 */
  ixstart = 1;
  mtmp = t[0];
  if (rtIsNaN(t[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 267)) {
      ixstart = ix;
      if (!rtIsNaN(t[ix - 1])) {
        mtmp = t[ix - 1];
        exitg1 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 266) {
    while (ixstart + 1 < 267) {
      if (t[ixstart] > mtmp) {
        mtmp = t[ixstart];
      }

      ixstart++;
    }
  }

  if (mtmp <= 0.0) {
    /* 'convolv_2invG_adapt_nov:84' P=realmin.*ones(size(t)); */
    for (ixstart = 0; ixstart < 266; ixstart++) {
      P[ixstart] = 2.2250738585072014E-308;
    }

    /*  otherwise we need to calculate P for each point in t */
  } else {
    /* 'convolv_2invG_adapt_nov:87' else */
    /*  find the largest point in t */
    /* 'convolv_2invG_adapt_nov:89' Maxt=max(t); */
    ixstart = 1;
    mtmp = t[0];
    if (rtIsNaN(t[0])) {
      ix = 2;
      exitg1 = false;
      while ((!exitg1) && (ix < 267)) {
        ixstart = ix;
        if (!rtIsNaN(t[ix - 1])) {
          mtmp = t[ix - 1];
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < 266) {
      while (ixstart + 1 < 267) {
        if (t[ixstart] > mtmp) {
          mtmp = t[ixstart];
        }

        ixstart++;
      }
    }

    /*  produce a range of evenly spaced points to evaluate at between */
    /*  0 and Maxt with step size h. the even spacing is important when */
    /*  calculating the convolution later */
    /* 'convolv_2invG_adapt_nov:94' x=0:h:Maxt; */
    emxInit_real_T(&x, 2);
    if (rtIsNaN(mtmp)) {
      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
      x->data[0] = rtNaN;
    } else if (mtmp < 0.0) {
      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
    } else if (rtIsInf(mtmp) && (0.0 == mtmp)) {
      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
      x->data[0] = rtNaN;
    } else {
      ndbl = floor(mtmp / 0.01 + 0.5);
      apnd = ndbl * 0.01;
      cdiff = apnd - mtmp;
      if (fabs(cdiff) < 4.4408920985006262E-16 * mtmp) {
        ndbl++;
        apnd = mtmp;
      } else if (cdiff > 0.0) {
        apnd = (ndbl - 1.0) * 0.01;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int)ndbl;
      } else {
        n = 0;
      }

      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = n;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
      if (n > 0) {
        x->data[0] = 0.0;
        if (n > 1) {
          x->data[n - 1] = apnd;
          ixstart = (n - 1) / 2;
          for (ix = 1; ix < ixstart; ix++) {
            ndbl = (double)ix * 0.01;
            x->data[ix] = ndbl;
            x->data[(n - ix) - 1] = apnd - ndbl;
          }

          if (ixstart << 1 == n - 1) {
            x->data[ixstart] = apnd / 2.0;
          } else {
            ndbl = (double)ixstart * 0.01;
            x->data[ixstart] = ndbl;
            x->data[ixstart + 1] = apnd - ndbl;
          }
        }
      }
    }

    emxInit_real_T(&y, 2);
    emxInit_real_T(&z, 2);

    /*  evaluate the sub-distributions at each point in x */
    /* 'convolv_2invG_adapt_nov:97' y=onestagepdf2(x,m(1),s(1)); */
    d_onestagepdf2(x, m[0], s[0], y);

    /* 'convolv_2invG_adapt_nov:98' z=onestagepdf2(x,m(2),s(2)); */
    d_onestagepdf2(x, m[1], s[1], z);

    /*  if the first pdf is very concentrated, check to see if it can be */
    /*  approximated as a point-mass distribution */
    /* 'convolv_2invG_adapt_nov:102' if sd(1)<.01 */
    if (sd[0] < 0.01) {
      /* 'convolv_2invG_adapt_nov:104' check2 = tailmass(m, s, eps, T2, sd); */
      ndbl = tailmass(m, s, 1.0 / m[1] * (sqrt(1.0 + 2.25 * (rt_powd_snf(s[1],
        4.0) / (m[1] * m[1]))) - 1.5 * (s[1] * s[1] / m[1])), sd);

      /* 'convolv_2invG_adapt_nov:106' if  check2<=eps/3 */
      if (ndbl <= 0.0033333333333333335) {
        /* If there is not much probability in the tails of the convolution, we use a Dirac Delta for part 1. */
        /* Flag is set to 1 to indicate this. */
        /* 'convolv_2invG_adapt_nov:109' flag=1; */
        /* 'convolv_2invG_adapt_nov:110' sigma=s(2); */
        /* 'convolv_2invG_adapt_nov:111' mu=m(2); */
        /* l for lag. */
        /* 'convolv_2invG_adapt_nov:113' l=1/m(1); */
        /* 'convolv_2invG_adapt_nov:114' P=onestagepdf_lag(t,mu,sigma,l); */
        onestagepdf_lag(t, m[1], s[1], 1.0 / m[0], P);
      } else {
        /* 'convolv_2invG_adapt_nov:115' else */
        /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
        /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
        /*  Outputs: P */
        /* 'convolv_2invG_adapt_nov:120' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
        f_approxconvolv(z, y, t, x, P, &logP0);

        /* Keep reducing the step size in the numerical integration until we are happy with the error. */
        /* 'convolv_2invG_adapt_nov:124' while E>=.001*abs(logP0) */
        while (E >= 0.001 * fabs(logP0)) {
          /* 'convolv_2invG_adapt_nov:125' h1=.5*h; */
          h1 *= 0.5;

          /* 'convolv_2invG_adapt_nov:126' x=0:h1:Maxt; */
          if (rtIsNaN(mtmp)) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            x->data[0] = rtNaN;
          } else if ((h1 == 0.0) || ((0.0 < mtmp) && (h1 < 0.0)) || ((mtmp < 0.0)
                      && (h1 > 0.0))) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 0;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          } else if (rtIsInf(mtmp) && (rtIsInf(h1) || (0.0 == mtmp))) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            x->data[0] = rtNaN;
          } else if (rtIsInf(h1)) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            x->data[0] = 0.0;
          } else if (floor(h1) == h1) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = (int)floor(mtmp / h1) + 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            ixstart = (int)floor(mtmp / h1);
            for (ix = 0; ix <= ixstart; ix++) {
              x->data[x->size[0] * ix] = h1 * (double)ix;
            }
          } else {
            ndbl = floor(mtmp / h1 + 0.5);
            apnd = ndbl * h1;
            if (h1 > 0.0) {
              cdiff = apnd - mtmp;
            } else {
              cdiff = mtmp - apnd;
            }

            if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(mtmp)) {
              ndbl++;
              apnd = mtmp;
            } else if (cdiff > 0.0) {
              apnd = (ndbl - 1.0) * h1;
            } else {
              ndbl++;
            }

            if (ndbl >= 0.0) {
              n = (int)ndbl;
            } else {
              n = 0;
            }

            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = n;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            if (n > 0) {
              x->data[0] = 0.0;
              if (n > 1) {
                x->data[n - 1] = apnd;
                ixstart = (n - 1) / 2;
                for (ix = 1; ix < ixstart; ix++) {
                  ndbl = (double)ix * h1;
                  x->data[ix] = ndbl;
                  x->data[(n - ix) - 1] = apnd - ndbl;
                }

                if (ixstart << 1 == n - 1) {
                  x->data[ixstart] = apnd / 2.0;
                } else {
                  ndbl = (double)ixstart * h1;
                  x->data[ixstart] = ndbl;
                  x->data[ixstart + 1] = apnd - ndbl;
                }
              }
            }
          }

          /* 'convolv_2invG_adapt_nov:127' y=onestagepdf2(x,m(1),s(1)); */
          d_onestagepdf2(x, m[0], s[0], y);

          /* 'convolv_2invG_adapt_nov:128' z=onestagepdf2(x,m(2),s(2)); */
          d_onestagepdf2(x, m[1], s[1], z);

          /* 'convolv_2invG_adapt_nov:130' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
          e_approxconvolv(z, y, h1, t, x, P, &ndbl);

          /* 'convolv_2invG_adapt_nov:132' E=abs(logP1-logP0); */
          E = fabs(ndbl - logP0);

          /* 'convolv_2invG_adapt_nov:133' P0=P1; */
          /* 'convolv_2invG_adapt_nov:134' logP0=logP1; */
          logP0 = ndbl;

          /* 'convolv_2invG_adapt_nov:135' h=h1; */
        }

        /* 'convolv_2invG_adapt_nov:137' P=P0; */
        /*  END FUNCTION DOTHECONVOLUTION_OUTER */
      }

      /*  pdf is not very concentrated so compute the convolution directly */
    } else {
      /* 'convolv_2invG_adapt_nov:141' else */
      /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
      /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
      /*  Outputs: P */
      /* 'convolv_2invG_adapt_nov:146' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
      f_approxconvolv(z, y, t, x, P, &logP0);

      /* 'convolv_2invG_adapt_nov:148' while E>=.001*abs(logP0) */
      while (E >= 0.001 * fabs(logP0)) {
        /* 'convolv_2invG_adapt_nov:149' h1=.5*h; */
        h1 *= 0.5;

        /* 'convolv_2invG_adapt_nov:150' x=0:h1:Maxt; */
        if (rtIsNaN(mtmp)) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          x->data[0] = rtNaN;
        } else if ((h1 == 0.0) || ((0.0 < mtmp) && (h1 < 0.0)) || ((mtmp < 0.0) &&
                    (h1 > 0.0))) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
        } else if (rtIsInf(mtmp) && (rtIsInf(h1) || (0.0 == mtmp))) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          x->data[0] = rtNaN;
        } else if (rtIsInf(h1)) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          x->data[0] = 0.0;
        } else if (floor(h1) == h1) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = (int)floor(mtmp / h1) + 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          ixstart = (int)floor(mtmp / h1);
          for (ix = 0; ix <= ixstart; ix++) {
            x->data[x->size[0] * ix] = h1 * (double)ix;
          }
        } else {
          ndbl = floor(mtmp / h1 + 0.5);
          apnd = ndbl * h1;
          if (h1 > 0.0) {
            cdiff = apnd - mtmp;
          } else {
            cdiff = mtmp - apnd;
          }

          if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(mtmp)) {
            ndbl++;
            apnd = mtmp;
          } else if (cdiff > 0.0) {
            apnd = (ndbl - 1.0) * h1;
          } else {
            ndbl++;
          }

          if (ndbl >= 0.0) {
            n = (int)ndbl;
          } else {
            n = 0;
          }

          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = n;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          if (n > 0) {
            x->data[0] = 0.0;
            if (n > 1) {
              x->data[n - 1] = apnd;
              ixstart = (n - 1) / 2;
              for (ix = 1; ix < ixstart; ix++) {
                ndbl = (double)ix * h1;
                x->data[ix] = ndbl;
                x->data[(n - ix) - 1] = apnd - ndbl;
              }

              if (ixstart << 1 == n - 1) {
                x->data[ixstart] = apnd / 2.0;
              } else {
                ndbl = (double)ixstart * h1;
                x->data[ixstart] = ndbl;
                x->data[ixstart + 1] = apnd - ndbl;
              }
            }
          }
        }

        /* 'convolv_2invG_adapt_nov:151' y=onestagepdf2(x,m(1),s(1)); */
        d_onestagepdf2(x, m[0], s[0], y);

        /* 'convolv_2invG_adapt_nov:152' z=onestagepdf2(x,m(2),s(2)); */
        d_onestagepdf2(x, m[1], s[1], z);

        /* 'convolv_2invG_adapt_nov:154' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
        e_approxconvolv(z, y, h1, t, x, P, &ndbl);

        /* 'convolv_2invG_adapt_nov:157' E=abs(logP1-logP0); */
        E = fabs(ndbl - logP0);

        /* 'convolv_2invG_adapt_nov:158' P0=P1; */
        /* 'convolv_2invG_adapt_nov:159' logP0=logP1; */
        logP0 = ndbl;

        /* 'convolv_2invG_adapt_nov:160' h=h1; */
      }

      /* 'convolv_2invG_adapt_nov:162' P=P0; */
      /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
    }

    emxFree_real_T(&z);
    emxFree_real_T(&y);
    emxFree_real_T(&x);
  }
}
#endif

/*
 * function [P,h,flag,E]=convolv_2invG_adapt_nov(t,m1,s1,m2,s2,h)
 */
#ifdef _OLD_MATLAB_CODE
void f_convolv_2invG_adapt_nov(const double t[266], double m1, double s1, double
  m2, double s2, double P[266])
{
  double h1;
  double E;
  double v[2];
  double m[2];
  double s[2];
  double I[2];
  int ix;
  double sd[2];
  int iidx[2];
  double b_s[2];
  int ixstart;
  double mtmp;
  boolean_T exitg1;
  emxArray_real_T *x;
  double ndbl;
  double apnd;
  double cdiff;
  emxArray_real_T *y;
  emxArray_real_T *z;
  int n;
  double logP0;
  h1 = 0.001;

  /*  This function evaluates the convolution of two inverse Gaussian */
  /*  distributions at vector t. */
  /*  If the variance in one of the distributions is very small so that the  */
  /*  distribution is close the a Dirac delta function, the convolution */
  /*  is approximated as a shifted inverse Gaussian, that is,  */
  /*  one part of the cell cycle is treated as deterministic in length.   */
  /*  In this case, the shift or lag is the average time to complete  */
  /*  the part of the cycle that has the smallest standard deviation. */
  /*  t is a vector of times to divide (or times to complete two parts of the */
  /*  cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1, */
  /*  s2=sigma2. */
  /*  Input parameters: */
  /*  t = the list of points at which to evaluate the distribution */
  /*  m1 = mu for the first distribution */
  /*  s1 = sigma for the first distribution */
  /*  m2 = mu for the second distribution */
  /*  s2 = sigma for the second distribution */
  /*  h = step size */
  /*  Outputs: */
  /*  P = probability at each point corresponding to points in t */
  /*  h = final step size */
  /*  flag = indicates if we used the Dirac delta */
  /*  E = is the relative error in the likelihood of the data due to the numerical integration */
  /*  log the parameters we were called with */
  /* fprintf(1,'convolv_2invG(m1=%.17g, s1=%.17g, m2=%.17g, s2=%.17g)\n', ...  */
  /*     m1, s1, m2, s2); */
  /* fprintf('    convolv_2invG_adapt(m1=%f s1=%f m2=%f s2=%f)\n',m1,s1,m2,s2); */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:37' flag=0; */
  /*  a constant that is used in determining if the Dirac approximation should be applied. */
  /* 'convolv_2invG_adapt_nov:40' eps=.01; */
  /*  ???? */
  /* 'convolv_2invG_adapt_nov:43' E=Inf; */
  E = rtInf;

  /*  number of points to be evaluated */
  /* 'convolv_2invG_adapt_nov:46' n=length(t); */
  /*  store m and s for each sub-distribution in a list */
  /* 'convolv_2invG_adapt_nov:49' m=[m1 m2]; */
  /* 'convolv_2invG_adapt_nov:50' s=[s1 s2]; */
  /*  make sure we dont have negative values for m or s by replacement with */
  /*  zero if we do */
  /* m=max(m,0); */
  /* s=max(s,0); */
  /*  reflect negative parameters into the positive space so we can use */
  /*  unconstrained optimization */
  /* 'convolv_2invG_adapt_nov:58' m=abs(m); */
  v[0] = m1;
  v[1] = m2;
  b_abs(v, m);

  /* 'convolv_2invG_adapt_nov:59' s=abs(s); */
  v[0] = s1;
  v[1] = s2;
  b_abs(v, s);

  /*  find the variance for both sub-distributions */
  /* 'convolv_2invG_adapt_nov:62' v=(s.^2)./(m.^3); */
  power(s, v);
  b_power(m, I);
  for (ix = 0; ix < 2; ix++) {
    v[ix] /= I[ix];
  }

  /*  find the standard deviation for each sub-distribution */
  /* 'convolv_2invG_adapt_nov:65' sd=v.^.5; */
  c_power(v, sd);

  /*  reorder m and s so that the sub-distribution with the smallest */
  /*  variance comes first.  So the fist part might be approximated as a Dirac delta. */
  /* 'convolv_2invG_adapt_nov:69' [v,I]=sort(v); */
  sort(v, iidx);

  /* 'convolv_2invG_adapt_nov:70' m=m(I); */
  for (ix = 0; ix < 2; ix++) {
    v[ix] = m[iidx[ix] - 1];
    I[ix] = iidx[ix];
  }

  /* 'convolv_2invG_adapt_nov:71' s=s(I); */
  for (ix = 0; ix < 2; ix++) {
    m[ix] = v[ix];
    b_s[ix] = s[(int)I[ix] - 1];
  }

  for (ix = 0; ix < 2; ix++) {
    s[ix] = b_s[ix];
  }

  /*  T1 appears to be unused, should probably remove */
  /*  T2 is only used inside FUNCTION CHECK_APPROXIMATABLE */
  /*  and should probably be moved inside that scope */
  /*  T2 is the mode of the second distribution, determined by a precomputed analytic expression. */
  /* 'convolv_2invG_adapt_nov:78' T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1))); */
  /* 'convolv_2invG_adapt_nov:79' T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2))); */
  /* When called from convolv_3invG all of the components of t may be negative. */
  /* In this case, the probability of each is set to realmin */
  /* 'convolv_2invG_adapt_nov:83' if max(t)<=0 */
  ixstart = 1;
  mtmp = t[0];
  if (rtIsNaN(t[0])) {
    ix = 2;
    exitg1 = false;
    while ((!exitg1) && (ix < 267)) {
      ixstart = ix;
      if (!rtIsNaN(t[ix - 1])) {
        mtmp = t[ix - 1];
        exitg1 = true;
      } else {
        ix++;
      }
    }
  }

  if (ixstart < 266) {
    while (ixstart + 1 < 267) {
      if (t[ixstart] > mtmp) {
        mtmp = t[ixstart];
      }

      ixstart++;
    }
  }

  if (mtmp <= 0.0) {
    /* 'convolv_2invG_adapt_nov:84' P=realmin.*ones(size(t)); */
    for (ixstart = 0; ixstart < 266; ixstart++) {
      P[ixstart] = 2.2250738585072014E-308;
    }

    /*  otherwise we need to calculate P for each point in t */
  } else {
    /* 'convolv_2invG_adapt_nov:87' else */
    /*  find the largest point in t */
    /* 'convolv_2invG_adapt_nov:89' Maxt=max(t); */
    ixstart = 1;
    mtmp = t[0];
    if (rtIsNaN(t[0])) {
      ix = 2;
      exitg1 = false;
      while ((!exitg1) && (ix < 267)) {
        ixstart = ix;
        if (!rtIsNaN(t[ix - 1])) {
          mtmp = t[ix - 1];
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < 266) {
      while (ixstart + 1 < 267) {
        if (t[ixstart] > mtmp) {
          mtmp = t[ixstart];
        }

        ixstart++;
      }
    }

    /*  produce a range of evenly spaced points to evaluate at between */
    /*  0 and Maxt with step size h. the even spacing is important when */
    /*  calculating the convolution later */
    /* 'convolv_2invG_adapt_nov:94' x=0:h:Maxt; */
    emxInit_real_T(&x, 2);
    if (rtIsNaN(mtmp)) {
      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
      x->data[0] = rtNaN;
    } else if (mtmp < 0.0) {
      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
    } else if (rtIsInf(mtmp) && (0.0 == mtmp)) {
      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
      x->data[0] = rtNaN;
    } else {
      ndbl = floor(mtmp / 0.001 + 0.5);
      apnd = ndbl * 0.001;
      cdiff = apnd - mtmp;
      if (fabs(cdiff) < 4.4408920985006262E-16 * mtmp) {
        ndbl++;
        apnd = mtmp;
      } else if (cdiff > 0.0) {
        apnd = (ndbl - 1.0) * 0.001;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int)ndbl;
      } else {
        n = 0;
      }

      ix = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = n;
      emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
      if (n > 0) {
        x->data[0] = 0.0;
        if (n > 1) {
          x->data[n - 1] = apnd;
          ixstart = (n - 1) / 2;
          for (ix = 1; ix < ixstart; ix++) {
            ndbl = (double)ix * 0.001;
            x->data[ix] = ndbl;
            x->data[(n - ix) - 1] = apnd - ndbl;
          }

          if (ixstart << 1 == n - 1) {
            x->data[ixstart] = apnd / 2.0;
          } else {
            ndbl = (double)ixstart * 0.001;
            x->data[ixstart] = ndbl;
            x->data[ixstart + 1] = apnd - ndbl;
          }
        }
      }
    }

    emxInit_real_T(&y, 2);
    emxInit_real_T(&z, 2);

    /*  evaluate the sub-distributions at each point in x */
    /* 'convolv_2invG_adapt_nov:97' y=onestagepdf2(x,m(1),s(1)); */
    d_onestagepdf2(x, m[0], s[0], y);

    /* 'convolv_2invG_adapt_nov:98' z=onestagepdf2(x,m(2),s(2)); */
    d_onestagepdf2(x, m[1], s[1], z);

    /*  if the first pdf is very concentrated, check to see if it can be */
    /*  approximated as a point-mass distribution */
    /* 'convolv_2invG_adapt_nov:102' if sd(1)<.01 */
    if (sd[0] < 0.01) {
      /* 'convolv_2invG_adapt_nov:104' check2 = tailmass(m, s, eps, T2, sd); */
      ndbl = tailmass(m, s, 1.0 / m[1] * (sqrt(1.0 + 2.25 * (rt_powd_snf(s[1],
        4.0) / (m[1] * m[1]))) - 1.5 * (s[1] * s[1] / m[1])), sd);

      /* 'convolv_2invG_adapt_nov:106' if  check2<=eps/3 */
      if (ndbl <= 0.0033333333333333335) {
        /* If there is not much probability in the tails of the convolution, we use a Dirac Delta for part 1. */
        /* Flag is set to 1 to indicate this. */
        /* 'convolv_2invG_adapt_nov:109' flag=1; */
        /* 'convolv_2invG_adapt_nov:110' sigma=s(2); */
        /* 'convolv_2invG_adapt_nov:111' mu=m(2); */
        /* l for lag. */
        /* 'convolv_2invG_adapt_nov:113' l=1/m(1); */
        /* 'convolv_2invG_adapt_nov:114' P=onestagepdf_lag(t,mu,sigma,l); */
        onestagepdf_lag(t, m[1], s[1], 1.0 / m[0], P);
      } else {
        /* 'convolv_2invG_adapt_nov:115' else */
        /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
        /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
        /*  Outputs: P */
        /* 'convolv_2invG_adapt_nov:120' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
        g_approxconvolv(z, y, t, x, P, &logP0);

        /* Keep reducing the step size in the numerical integration until we are happy with the error. */
        /* 'convolv_2invG_adapt_nov:124' while E>=.001*abs(logP0) */
        while (E >= 0.001 * fabs(logP0)) {
          /* 'convolv_2invG_adapt_nov:125' h1=.5*h; */
          h1 *= 0.5;

          /* 'convolv_2invG_adapt_nov:126' x=0:h1:Maxt; */
          if (rtIsNaN(mtmp)) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            x->data[0] = rtNaN;
          } else if ((h1 == 0.0) || ((0.0 < mtmp) && (h1 < 0.0)) || ((mtmp < 0.0)
                      && (h1 > 0.0))) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 0;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          } else if (rtIsInf(mtmp) && (rtIsInf(h1) || (0.0 == mtmp))) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            x->data[0] = rtNaN;
          } else if (rtIsInf(h1)) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            x->data[0] = 0.0;
          } else if (floor(h1) == h1) {
            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = (int)floor(mtmp / h1) + 1;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            ixstart = (int)floor(mtmp / h1);
            for (ix = 0; ix <= ixstart; ix++) {
              x->data[x->size[0] * ix] = h1 * (double)ix;
            }
          } else {
            ndbl = floor(mtmp / h1 + 0.5);
            apnd = ndbl * h1;
            if (h1 > 0.0) {
              cdiff = apnd - mtmp;
            } else {
              cdiff = mtmp - apnd;
            }

            if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(mtmp)) {
              ndbl++;
              apnd = mtmp;
            } else if (cdiff > 0.0) {
              apnd = (ndbl - 1.0) * h1;
            } else {
              ndbl++;
            }

            if (ndbl >= 0.0) {
              n = (int)ndbl;
            } else {
              n = 0;
            }

            ix = x->size[0] * x->size[1];
            x->size[0] = 1;
            x->size[1] = n;
            emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
            if (n > 0) {
              x->data[0] = 0.0;
              if (n > 1) {
                x->data[n - 1] = apnd;
                ixstart = (n - 1) / 2;
                for (ix = 1; ix < ixstart; ix++) {
                  ndbl = (double)ix * h1;
                  x->data[ix] = ndbl;
                  x->data[(n - ix) - 1] = apnd - ndbl;
                }

                if (ixstart << 1 == n - 1) {
                  x->data[ixstart] = apnd / 2.0;
                } else {
                  ndbl = (double)ixstart * h1;
                  x->data[ixstart] = ndbl;
                  x->data[ixstart + 1] = apnd - ndbl;
                }
              }
            }
          }

          /* 'convolv_2invG_adapt_nov:127' y=onestagepdf2(x,m(1),s(1)); */
          d_onestagepdf2(x, m[0], s[0], y);

          /* 'convolv_2invG_adapt_nov:128' z=onestagepdf2(x,m(2),s(2)); */
          d_onestagepdf2(x, m[1], s[1], z);

          /* 'convolv_2invG_adapt_nov:130' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
          e_approxconvolv(z, y, h1, t, x, P, &ndbl);

          /* 'convolv_2invG_adapt_nov:132' E=abs(logP1-logP0); */
          E = fabs(ndbl - logP0);

          /* 'convolv_2invG_adapt_nov:133' P0=P1; */
          /* 'convolv_2invG_adapt_nov:134' logP0=logP1; */
          logP0 = ndbl;

          /* 'convolv_2invG_adapt_nov:135' h=h1; */
        }

        /* 'convolv_2invG_adapt_nov:137' P=P0; */
        /*  END FUNCTION DOTHECONVOLUTION_OUTER */
      }

      /*  pdf is not very concentrated so compute the convolution directly */
    } else {
      /* 'convolv_2invG_adapt_nov:141' else */
      /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
      /*  Input parameters: E, m, s, Maxt, z, y, h, n, t, i, I, x */
      /*  Outputs: P */
      /* 'convolv_2invG_adapt_nov:146' [P0, logP0] = approxconvolv( z, y, h, n, t, I, x ); */
      g_approxconvolv(z, y, t, x, P, &logP0);

      /* 'convolv_2invG_adapt_nov:148' while E>=.001*abs(logP0) */
      while (E >= 0.001 * fabs(logP0)) {
        /* 'convolv_2invG_adapt_nov:149' h1=.5*h; */
        h1 *= 0.5;

        /* 'convolv_2invG_adapt_nov:150' x=0:h1:Maxt; */
        if (rtIsNaN(mtmp)) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          x->data[0] = rtNaN;
        } else if ((h1 == 0.0) || ((0.0 < mtmp) && (h1 < 0.0)) || ((mtmp < 0.0) &&
                    (h1 > 0.0))) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
        } else if (rtIsInf(mtmp) && (rtIsInf(h1) || (0.0 == mtmp))) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          x->data[0] = rtNaN;
        } else if (rtIsInf(h1)) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          x->data[0] = 0.0;
        } else if (floor(h1) == h1) {
          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = (int)floor(mtmp / h1) + 1;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          ixstart = (int)floor(mtmp / h1);
          for (ix = 0; ix <= ixstart; ix++) {
            x->data[x->size[0] * ix] = h1 * (double)ix;
          }
        } else {
          ndbl = floor(mtmp / h1 + 0.5);
          apnd = ndbl * h1;
          if (h1 > 0.0) {
            cdiff = apnd - mtmp;
          } else {
            cdiff = mtmp - apnd;
          }

          if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(mtmp)) {
            ndbl++;
            apnd = mtmp;
          } else if (cdiff > 0.0) {
            apnd = (ndbl - 1.0) * h1;
          } else {
            ndbl++;
          }

          if (ndbl >= 0.0) {
            n = (int)ndbl;
          } else {
            n = 0;
          }

          ix = x->size[0] * x->size[1];
          x->size[0] = 1;
          x->size[1] = n;
          emxEnsureCapacity((emxArray__common *)x, ix, sizeof(double));
          if (n > 0) {
            x->data[0] = 0.0;
            if (n > 1) {
              x->data[n - 1] = apnd;
              ixstart = (n - 1) / 2;
              for (ix = 1; ix < ixstart; ix++) {
                ndbl = (double)ix * h1;
                x->data[ix] = ndbl;
                x->data[(n - ix) - 1] = apnd - ndbl;
              }

              if (ixstart << 1 == n - 1) {
                x->data[ixstart] = apnd / 2.0;
              } else {
                ndbl = (double)ixstart * h1;
                x->data[ixstart] = ndbl;
                x->data[ixstart + 1] = apnd - ndbl;
              }
            }
          }
        }

        /* 'convolv_2invG_adapt_nov:151' y=onestagepdf2(x,m(1),s(1)); */
        d_onestagepdf2(x, m[0], s[0], y);

        /* 'convolv_2invG_adapt_nov:152' z=onestagepdf2(x,m(2),s(2)); */
        d_onestagepdf2(x, m[1], s[1], z);

        /* 'convolv_2invG_adapt_nov:154' [P1, logP1] = approxconvolv( z, y, h1, n, t, I, x ); */
        e_approxconvolv(z, y, h1, t, x, P, &ndbl);

        /* 'convolv_2invG_adapt_nov:157' E=abs(logP1-logP0); */
        E = fabs(ndbl - logP0);

        /* 'convolv_2invG_adapt_nov:158' P0=P1; */
        /* 'convolv_2invG_adapt_nov:159' logP0=logP1; */
        logP0 = ndbl;

        /* 'convolv_2invG_adapt_nov:160' h=h1; */
      }

      /* 'convolv_2invG_adapt_nov:162' P=P0; */
      /*  BEGIN FUNCTION DOTHECONVOLUTION_OUTER */
    }

    emxFree_real_T(&z);
    emxFree_real_T(&y);
    emxFree_real_T(&x);
  }
}
#endif
/* End of code generation (convolv_2invG_adapt_nov.c) */
