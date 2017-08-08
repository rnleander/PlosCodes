/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * tailmass.c
 *
 * Code generation for function 'tailmass'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "tailmass.h"
#include "onestagepdf2.h"
#include "sum.h"
#include "IMT_analysis_April2017_emxutil.h"
#include "gp_max.h"

#define _GSL_GP_MAX_FIXED

/* Function Definitions */

int checktailmass(const double m_a, const double s_a, const double m_b, const double s_b, double T2, const double sd_a, const double sd_b)
{
	double m[2];
	double s[2];
	double sd[2];
	m[0] = m_a;
	m[1] = m_b;
	s[0] = s_a;
	s[1] = s_b;
	sd[0] = sd_a;
	sd[1] = sd_b;

	if (tailmass(m, s, T2, sd) <= (T2 / 3.0))
		return 1;
	else
		return 0;
}

/*
 * function check2 = tailmass( m, s, eps, T2, sd )
 */
double tailmass(const double m[2], const double s[2], double T2, const double
                sd[2])
{
  double check2;
  double gp;
  double r;
  double nu;
  double Tu;
  emxArray_real_T *y;
  double checkerror;
  double hh;
  double nuplusr;
  int k;
  double ndbl;
  double apnd;
  double cdiff;
  emxArray_real_T *r0;
  double LeftTail;
  int n;
  int nm1d2;
  emxArray_real_T *b_y;

  /* TAILMASS Summary of this function goes here */
  /*    Detailed explanation goes here */
  /*  BEGIN FUNCTION TailMass */
  /*  Input parameters: m, s, eps, T2, nu, sd, */
  /*  Outputs: check2 */
  /*  to estimate the error (in calculating probabilities the of the */
  /*  data), that results from approximating the first pdf as a */
  /*  point-mass distribtion, find the maximum of the absolute value */
  /*  of the derivative of the second pdf. */
  /* 'tailmass:14' gp=gp_max(m(2),s(2)); */
#ifdef _GSL_GP_MAX_FIXED
  gp = gp_max_fixed(m[1], s[1]);
#else
  gp = gp_max(m[1], s[1]);
#endif

  /*  determine the radius, r, of a small interval over which the */
  /*  1. second pdf, g, is approximately constant, i.e. changes by less than eps/3 */
  /*  over any interval with that radius */
  /*  and 2. g(t) is small for t<r.  ????? */
  /* 'tailmass:20' r=min(eps/(3*gp),T2); */
  gp = 0.01 / (3.0 * gp);
  if ((gp < T2) || rtIsNaN(T2)) {
    r = gp;
  } else {
    r = T2;
  }

  /* 'tailmass:21' checkval=onestagepdf2(r,m(2),s(2)); */
  gp = c_onestagepdf2(r, m[1], s[1]);

  /* 'tailmass:22' while checkval>=eps/2 */
  while (gp >= 0.005) {
    /* 'tailmass:23' r=r/2; */
    r /= 2.0;

    /* 'tailmass:24' checkval=onestagepdf2(r,m(2),s(2)); */
    gp = c_onestagepdf2(r, m[1], s[1]);
  }

  /*  get the average value of the first pdf.  This is the point at */
  /*  which its mass is concentrated. */
  /* 'tailmass:29' nu=1/m(1); */
  nu = 1.0 / m[0];

  /*  get the maximum value of the second pdf,f. */
  /* 'tailmass:32' gm=onestagepdf2(T2,m(2),s(2)); */
  /*  ???? */
  /*  Tu is the upper limit of integral for approximating */
  /*  int_{r+nu}^{infty}g(s)ds. This is in latex. */
  /* 'tailmass:37' Tu=max(100,nu+r+1000*sd(1)); */
  gp = (nu + r) + 1000.0 * sd[0];
  if ((100.0 > gp) || rtIsNaN(gp)) {
    Tu = 100.0;
  } else {
    Tu = gp;
  }

  emxInit_real_T(&y, 2);

  /*  ???? */
  /* 'tailmass:40' checkerror=100; */
  checkerror = 100.0;

  /* 'tailmass:41' hh=.001; */
  hh = 0.001;

  /* 'tailmass:42' numinusr=nu-r; */
  gp = nu - r;

  /* 'tailmass:43' nuplusr=nu+r; */
  nuplusr = nu + r;

  /* 'tailmass:44' teeyou=Tu; */
  /* we are integrating the first pdf from zero to nu-r, and from nu+r to infinity to see if the tails */
  /* are small in probability. */
  /* 'tailmass:47' LeftTail=.001*sum(onestagepdf2((0:.001:numinusr),m(1),s(1))); */
  if (rtIsNaN(gp)) {
    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
    y->data[0] = rtNaN;
  } else if (gp < 0.0) {
    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  } else if (rtIsInf(gp) && (0.0 == gp)) {
    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
    y->data[0] = rtNaN;
  } else {
    ndbl = floor(gp / 0.001 + 0.5);
    apnd = ndbl * 0.001;
    cdiff = apnd - gp;
    if (fabs(cdiff) < 4.4408920985006262E-16 * gp) {
      ndbl++;
      apnd = gp;
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

    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n;
    emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
    if (n > 0) {
      y->data[0] = 0.0;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nm1d2 = (n - 1) / 2;
        for (k = 1; k < nm1d2; k++) {
          gp = (double)k * 0.001;
          y->data[k] = gp;
          y->data[(n - k) - 1] = apnd - gp;
        }

        if (nm1d2 << 1 == n - 1) {
          y->data[nm1d2] = apnd / 2.0;
        } else {
          gp = (double)nm1d2 * 0.001;
          y->data[nm1d2] = gp;
          y->data[nm1d2 + 1] = apnd - gp;
        }
      }
    }
  }

  emxInit_real_T(&r0, 2);
  d_onestagepdf2(y, m[0], s[0], r0);
  LeftTail = 0.001 * b_sum(r0);

  /* 'tailmass:48' RightTail=.001*sum(onestagepdf2((nuplusr:.001:teeyou),m(1),s(1))); */
  if (rtIsNaN(nuplusr) || rtIsNaN(Tu)) {
    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
    y->data[0] = rtNaN;
  } else if (Tu < nuplusr) {
    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  } else if ((rtIsInf(nuplusr) || rtIsInf(Tu)) && (nuplusr == Tu)) {
    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
    y->data[0] = rtNaN;
  } else {
    ndbl = floor((Tu - nuplusr) / 0.001 + 0.5);
    apnd = nuplusr + ndbl * 0.001;
    cdiff = apnd - Tu;
    gp = fabs(nuplusr);
    if (!((gp > Tu) || rtIsNaN(Tu))) {
      gp = Tu;
    }

    if (fabs(cdiff) < 4.4408920985006262E-16 * gp) {
      ndbl++;
      apnd = Tu;
    } else if (cdiff > 0.0) {
      apnd = nuplusr + (ndbl - 1.0) * 0.001;
    } else {
      ndbl++;
    }

    if (ndbl >= 0.0) {
      n = (int)ndbl;
    } else {
      n = 0;
    }

    k = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n;
    emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
    if (n > 0) {
      y->data[0] = nuplusr;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nm1d2 = (n - 1) / 2;
        for (k = 1; k < nm1d2; k++) {
          gp = (double)k * 0.001;
          y->data[k] = nuplusr + gp;
          y->data[(n - k) - 1] = apnd - gp;
        }

        if (nm1d2 << 1 == n - 1) {
          y->data[nm1d2] = (nuplusr + apnd) / 2.0;
        } else {
          gp = (double)nm1d2 * 0.001;
          y->data[nm1d2] = nuplusr + gp;
          y->data[nm1d2 + 1] = apnd - gp;
        }
      }
    }
  }

  d_onestagepdf2(y, m[0], s[0], r0);
  gp = 0.001 * b_sum(r0);

  /* 'tailmass:49' check1=LeftTail+RightTail; */
  LeftTail += gp;

  /* Reduce step size in above Riemann sum until the error is small, */
  /* meaning that the Riemann sum is converging. */
  /* 'tailmass:52' while checkerror>10^(-4) */
  emxInit_real_T(&b_y, 2);
  while (checkerror > 0.0001) {
    /* 'tailmass:53' hh=.5*hh; */
    hh *= 0.5;

    /* 'tailmass:54' ck1=hh*sum(onestagepdf2((0:hh:nu-r),m(1),s(1)))+hh*sum(onestagepdf2((nu+r:hh:Tu),m(1),s(1))); */
    gp = nu - r;
    if (rtIsNaN(gp)) {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
      y->data[0] = rtNaN;
    } else if ((hh == 0.0) || ((0.0 < gp) && (hh < 0.0)) || ((gp < 0.0) && (hh >
      0.0))) {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
    } else if (rtIsInf(gp) && (rtIsInf(hh) || (0.0 == gp))) {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
      y->data[0] = rtNaN;
    } else if (rtIsInf(hh)) {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
      y->data[0] = 0.0;
    } else if (floor(hh) == hh) {
      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = (int)floor(gp / hh) + 1;
      emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
      nm1d2 = (int)floor(gp / hh);
      for (k = 0; k <= nm1d2; k++) {
        y->data[y->size[0] * k] = hh * (double)k;
      }
    } else {
      ndbl = floor(gp / hh + 0.5);
      apnd = ndbl * hh;
      if (hh > 0.0) {
        cdiff = apnd - gp;
      } else {
        cdiff = gp - apnd;
      }

      if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(gp)) {
        ndbl++;
        apnd = gp;
      } else if (cdiff > 0.0) {
        apnd = (ndbl - 1.0) * hh;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int)ndbl;
      } else {
        n = 0;
      }

      k = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = n;
      emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
      if (n > 0) {
        y->data[0] = 0.0;
        if (n > 1) {
          y->data[n - 1] = apnd;
          nm1d2 = (n - 1) / 2;
          for (k = 1; k < nm1d2; k++) {
            gp = (double)k * hh;
            y->data[k] = gp;
            y->data[(n - k) - 1] = apnd - gp;
          }

          if (nm1d2 << 1 == n - 1) {
            y->data[nm1d2] = apnd / 2.0;
          } else {
            gp = (double)nm1d2 * hh;
            y->data[nm1d2] = gp;
            y->data[nm1d2 + 1] = apnd - gp;
          }
        }
      }
    }

    nuplusr = nu + r;
    if (rtIsNaN(nuplusr) || rtIsNaN(Tu)) {
      k = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)b_y, k, sizeof(double));
      b_y->data[0] = rtNaN;
    } else if ((hh == 0.0) || ((nuplusr < Tu) && (hh < 0.0)) || ((Tu < nuplusr) &&
                (hh > 0.0))) {
      k = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)b_y, k, sizeof(double));
    } else if ((rtIsInf(nuplusr) || rtIsInf(Tu)) && (rtIsInf(hh) || (nuplusr ==
                 Tu))) {
      k = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)b_y, k, sizeof(double));
      b_y->data[0] = rtNaN;
    } else if (rtIsInf(hh)) {
      k = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)b_y, k, sizeof(double));
      b_y->data[0] = nuplusr;
    } else if ((floor(nuplusr) == nuplusr) && (floor(hh) == hh)) {
      k = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = (int)floor((Tu - nuplusr) / hh) + 1;
      emxEnsureCapacity((emxArray__common *)b_y, k, sizeof(double));
      nm1d2 = (int)floor((Tu - nuplusr) / hh);
      for (k = 0; k <= nm1d2; k++) {
        b_y->data[b_y->size[0] * k] = nuplusr + hh * (double)k;
      }
    } else {
      ndbl = floor((Tu - nuplusr) / hh + 0.5);
      apnd = nuplusr + ndbl * hh;
      if (hh > 0.0) {
        cdiff = apnd - Tu;
      } else {
        cdiff = Tu - apnd;
      }

      gp = fabs(nuplusr);
      if (!((gp > Tu) || rtIsNaN(Tu))) {
        gp = Tu;
      }

      if (fabs(cdiff) < 4.4408920985006262E-16 * gp) {
        ndbl++;
        apnd = Tu;
      } else if (cdiff > 0.0) {
        apnd = nuplusr + (ndbl - 1.0) * hh;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int)ndbl;
      } else {
        n = 0;
      }

      k = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = n;
      emxEnsureCapacity((emxArray__common *)b_y, k, sizeof(double));
      if (n > 0) {
        b_y->data[0] = nuplusr;
        if (n > 1) {
          b_y->data[n - 1] = apnd;
          nm1d2 = (n - 1) / 2;
          for (k = 1; k < nm1d2; k++) {
            gp = (double)k * hh;
            b_y->data[k] = nuplusr + gp;
            b_y->data[(n - k) - 1] = apnd - gp;
          }

          if (nm1d2 << 1 == n - 1) {
            b_y->data[nm1d2] = (nuplusr + apnd) / 2.0;
          } else {
            gp = (double)nm1d2 * hh;
            b_y->data[nm1d2] = nuplusr + gp;
            b_y->data[nm1d2 + 1] = apnd - gp;
          }
        }
      }
    }

    d_onestagepdf2(y, m[0], s[0], r0);
    d_onestagepdf2(b_y, m[0], s[0], y);
    gp = hh * b_sum(r0) + hh * b_sum(y);

    /* 'tailmass:55' checkerror=abs(check1-ck1); */
    checkerror = fabs(LeftTail - gp);

    /* 'tailmass:56' check1=ck1; */
    LeftTail = gp;
  }

  emxFree_real_T(&r0);
  emxFree_real_T(&b_y);
  emxFree_real_T(&y);

  /* 'tailmass:58' check2=gm*check1; */
  check2 = c_onestagepdf2(T2, m[1], s[1]) * LeftTail;

  /*      end */
  /*  END FUNCTION TailMass */
  return check2;
}

/* End of code generation (tailmass.c) */
