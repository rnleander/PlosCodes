/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * emgpdf.c
 *
 * Code generation for function 'emgpdf'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "emgpdf.h"
#include "IMT_analysis_April2017_rtwutil.h"
#include "sum.h"
#include "log.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_sf_erf.h"

/* Function Definitions */

double emgpdf_loglikelihood(const gsl_vector *v, void *params)
{
	double *data = (double *)params;

	double l = gsl_vector_get(v, 0);
	double m = gsl_vector_get(v, 1);
	double s = gsl_vector_get(v, 2);

	double penalty = 0;
	if (m < 0 || s < 0 || l < 0)
		penalty = 1000;

	l = fabs(l);
	m = fabs(m);
	s = fabs(s);

	double Y[266];

	// emgpdf_replacement(data, l, m, s, Y);
	emgpdf(data, l, m, s, Y);
	

	double loglikelihood = 0;
	for (int i = 0; i < 266; i++) {
		loglikelihood += log(Y[i]);
	}

	return penalty - loglikelihood;
}

void emgpdf(const double X[266], double l, double m, double s, double Y[266])
{
	// Y=(l/2)*erfc((-X+m+l*s^2)/(s*2^(1/2))).*exp((l/2)*(-2*X+2*m+l*s^2));
	// https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
	double a, b;
	for (int i = 0; i < 266; i++) {
		a = 0.5 * l * (2*m + l*pow(s, 2.0) - 2*X[i] );
		b = (m + l*pow(s, 2.0) - X[i]) / (s * pow(2, 0.5));
		Y[i] = 0.5 * l * exp(a) * gsl_sf_erfc(b);
	}
}

void emgpdf_old(const double X[266], double l, double m, double s, double Y[266])
{
  double y;
  double c;
  double B;
  int i;
  int k;
  double b_y[266];
  double c_y[266];
  double d_y;
  double absx;
  double b_s;
  double S;
  double R;
  int eint;
  int e;

  /* 'emgpdf:3' l=abs(l); */
  l = fabs(l);

  /* 'emgpdf:4' m=abs(m); */
  m = fabs(m);

  /* 'emgpdf:5' s=abs(s); */
  s = fabs(s);

  /* 'emgpdf:8' Y=(l/2)*erfc((-X+m+l*s^2)/(s*2^(1/2))).*exp((l/2)*(-2*X+2*m+l*s^2)); */
  y = l / 2.0;
  c = s * s;
  B = s * 1.4142135623730951;
  for (i = 0; i < 266; i++) {
    b_y[i] = ((-X[i] + m) + l * c) / B;
  }

#pragma omp parallel for \
 num_threads(omp_get_max_threads()) \
 private(absx,b_s,S,R,e) \
 firstprivate(eint)

  for (k = 1; k < 267; k++) {
    /* ========================== COPYRIGHT NOTICE ============================ */
    /*  The algorithms for calculating ERF(X) and ERFC(X) are derived           */
    /*  from FDLIBM, which has the following notice:                            */
    /*                                                                          */
    /*  Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.       */
    /*                                                                          */
    /*  Developed at SunSoft, a Sun Microsystems, Inc. business.                */
    /*  Permission to use, copy, modify, and distribute this                    */
    /*  software is freely granted, provided that this notice                   */
    /*  is preserved.                                                           */
    /* =============================    END    ================================ */
    absx = fabs(b_y[k - 1]);
    if (rtIsNaN(b_y[k - 1])) {
      b_s = b_y[k - 1];
    } else if (rtIsInf(b_y[k - 1])) {
      if (b_y[k - 1] < 0.0) {
        b_s = 2.0;
      } else {
        b_s = 0.0;
      }
    } else if (absx < 0.84375) {
      if (absx < 1.3877787807814457E-17) {
        b_s = 1.0 - b_y[k - 1];
      } else {
        b_s = b_y[k - 1] * b_y[k - 1];
        b_s = (0.12837916709551256 + b_s * (-0.3250421072470015 + b_s *
                (-0.02848174957559851 + b_s * (-0.0057702702964894416 + b_s *
                  -2.3763016656650163E-5)))) / (1.0 + b_s * (0.39791722395915535
          + b_s * (0.0650222499887673 + b_s * (0.0050813062818757656 + b_s *
          (0.00013249473800432164 + b_s * -3.9602282787753681E-6)))));
        if (b_y[k - 1] < 0.25) {
          b_s = 1.0 - (b_y[k - 1] + b_y[k - 1] * b_s);
        } else {
          b_s = 0.5 - (b_y[k - 1] * b_s + (b_y[k - 1] - 0.5));
        }
      }
    } else if (absx < 1.25) {
      S = -0.0023621185607526594 + (absx - 1.0) * (0.41485611868374833 + (absx -
        1.0) * (-0.37220787603570132 + (absx - 1.0) * (0.31834661990116175 +
        (absx - 1.0) * (-0.11089469428239668 + (absx - 1.0) *
                        (0.035478304325618236 + (absx - 1.0) *
                         -0.0021663755948687908)))));
      b_s = 1.0 + (absx - 1.0) * (0.10642088040084423 + (absx - 1.0) *
        (0.540397917702171 + (absx - 1.0) * (0.071828654414196266 + (absx - 1.0)
        * (0.12617121980876164 + (absx - 1.0) * (0.013637083912029051 + (absx -
        1.0) * 0.011984499846799107)))));
      if (b_y[k - 1] >= 0.0) {
        b_s = 0.15493708848953247 - S / b_s;
      } else {
        b_s = 1.0 + (0.84506291151046753 + S / b_s);
      }
    } else if (b_y[k - 1] < -6.0) {
      b_s = 2.0;
    } else if (b_y[k - 1] >= 28.0) {
      b_s = 0.0;
    } else {
      b_s = 1.0 / (absx * absx);
      if (absx < 2.8571414947509766) {
        R = -0.0098649440348471482 + b_s * (-0.69385857270718176 + b_s *
          (-10.558626225323291 + b_s * (-62.375332450326006 + b_s *
          (-162.39666946257347 + b_s * (-184.60509290671104 + b_s *
          (-81.2874355063066 + b_s * -9.8143293441691455))))));
        S = 1.0 + b_s * (19.651271667439257 + b_s * (137.65775414351904 + b_s *
          (434.56587747522923 + b_s * (645.38727173326788 + b_s *
          (429.00814002756783 + b_s * (108.63500554177944 + b_s *
          (6.5702497703192817 + b_s * -0.0604244152148581)))))));
      } else {
        R = -0.0098649429247001 + b_s * (-0.799283237680523 + b_s *
          (-17.757954917754752 + b_s * (-160.63638485582192 + b_s *
          (-637.56644336838963 + b_s * (-1025.0951316110772 + b_s *
          -483.5191916086514)))));
        S = 1.0 + b_s * (30.338060743482458 + b_s * (325.79251299657392 + b_s *
          (1536.729586084437 + b_s * (3199.8582195085955 + b_s *
          (2553.0504064331644 + b_s * (474.52854120695537 + b_s *
          -22.440952446585818))))));
      }

      if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
        b_s = frexp(absx, &eint);
        e = eint;
      } else {
        b_s = absx;
        e = 0;
      }

      b_s = floor(b_s * 2.097152E+6) / 2.097152E+6 * rt_powd_snf(2.0, e);
      b_s = exp(-b_s * b_s - 0.5625) * exp((b_s - absx) * (b_s + absx) + R / S) /
        absx;
      if (b_y[k - 1] < 0.0) {
        b_s = 2.0 - b_s;
      }
    }

    c_y[k - 1] = b_s;
  }

  B = l / 2.0;
  c = s * s;

  /* 'emgpdf:9' Y=max(realmin, Y); */
  for (i = 0; i < 266; i++) {
    d_y = y * c_y[i] * exp(B * ((-2.0 * X[i] + 2.0 * m) + l * c));
    if ((2.2250738585072014E-308 > d_y) || rtIsNaN(d_y)) {
      Y[i] = 2.2250738585072014E-308;
    } else {
      Y[i] = d_y;
    }
  }
}

/* End of code generation (emgpdf.c) */
