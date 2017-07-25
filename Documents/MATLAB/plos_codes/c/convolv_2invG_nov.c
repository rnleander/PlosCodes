/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * convolv_2invG_nov.c
 *
 * Code generation for function 'convolv_2invG_nov'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "convolv_2invG_nov.h"
#include "IMT_analysis_April2017_emxutil.h"
#include "onestagepdf2.h"
#include "sum.h"
#include "rdivide.h"
#include "exp.h"
#include "power.h"
#include "log.h"
#include "conv.h"
#include "gp_max.h"
#include "emgpdf.h"
#include "sort1.h"
#include "onestagepdf_lag.h"
#include "IMT_analysis_April2017_rtwutil.h"

/* Function Definitions */

/*
 * function [P,flag]=convolv_2invG_nov(t,m1,s1,m2,s2,h)
 */
void b_convolv_2invG_nov(const emxArray_real_T *t, double m1, double s1, double
  m2, double s2, double h, emxArray_real_T *P, double *flag)
{
  int b_flag;
  int n;
  double y[2];
  int ixstart;
  double m[2];
  double s[2];
  double v[2];
  int nm1d2;
  int iidx[2];
  double sd[2];
  double b_s[2];
  double T2;
  int b_n;
  double kd;
  boolean_T exitg1;
  emxArray_real_T *x;
  double ndbl;
  double apnd;
  emxArray_real_T *b_y;
  double cdiff;
  emxArray_real_T *z;
  emxArray_real_T *C;
  emxArray_real_T *I;
  emxArray_real_T *b_P;
  emxArray_real_T *varargin_1;
  emxArray_real_T *r3;
  double r;
  double nu;
  double Tu;
  double checkerror;
  int i;
  double hh;
  emxArray_real_T *b_t;
  double a;
  int itmp;
  double check1;
  emxArray_real_T *b;
  emxArray_real_T *r4;
  emxArray_real_T *b_m;
  emxArray_real_T *r5;
  emxArray_real_T *b_a;
  emxArray_real_T *c_s;
  unsigned int P_idx_0;
  unsigned int c_P;

  /* This function evaluates the convolution of two inverse Gaussian */
  /* distributions at vector t. */
  /* If the variance in one of the distributions is very small so that the  */
  /* distribution is close the a Dirac delta function, the convolution */
  /* is approximated as a shifted inverse Gaussian, that is  */
  /* one part of the cell cycle is treated as deterministic in length.   */
  /* In this case, the shift or lag is the average time to complete  */
  /* the part of the cycle with the smallest sigma. */
  /* t is a vector of times to divide (or times to complete two parts of the */
  /* cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1, */
  /* s2=sigma2. */
  /* 'convolv_2invG_nov:16' flag=0; */
  b_flag = 0;

  /* 'convolv_2invG_nov:18' E=Inf; */
  /* 'convolv_2invG_nov:19' n=length(t); */
  n = t->size[0] - 1;

  /* 'convolv_2invG_nov:20' m=[m1 m2]; */
  /* 'convolv_2invG_nov:21' s=[s1 s2]; */
  /* 'convolv_2invG_nov:22' m=max(0,m); */
  y[0] = m1;
  y[1] = m2;
  for (ixstart = 0; ixstart < 2; ixstart++) {
    if ((0.0 > y[ixstart]) || rtIsNaN(y[ixstart])) {
      m[ixstart] = 0.0;
    } else {
      m[ixstart] = y[ixstart];
    }
  }

  /* 'convolv_2invG_nov:23' s=max(0,s); */
  y[0] = s1;
  y[1] = s2;
  for (ixstart = 0; ixstart < 2; ixstart++) {
    if ((0.0 > y[ixstart]) || rtIsNaN(y[ixstart])) {
      s[ixstart] = 0.0;
    } else {
      s[ixstart] = y[ixstart];
    }
  }

  /* 'convolv_2invG_nov:24' v=(s.^2)./(m.^3); */
  /* 'convolv_2invG_nov:25' [v,I]=sort(v); */
  power(s, v);
  b_power(m, y);
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    v[nm1d2] /= y[nm1d2];
  }

  sort(v, iidx);

  /* 'convolv_2invG_nov:26' sd=v.^.5; */
  c_power(v, sd);

  /* 'convolv_2invG_nov:27' m=m(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    v[nm1d2] = m[iidx[nm1d2] - 1];
    y[nm1d2] = iidx[nm1d2];
  }

  /* 'convolv_2invG_nov:28' s=s(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    m[nm1d2] = v[nm1d2];
    b_s[nm1d2] = s[(int)y[nm1d2] - 1];
  }

  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    s[nm1d2] = b_s[nm1d2];
  }

  /* 'convolv_2invG_nov:30' T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1))); */
  /* 'convolv_2invG_nov:31' T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2))); */
  T2 = 1.0 / m[1] * (sqrt(1.0 + 2.25 * (rt_powd_snf(s[1], 4.0) / (m[1] * m[1])))
                     - 1.5 * (s[1] * s[1] / m[1]));

  /* When called from convolv_3invG all of the components of t may be negative. */
  /* In this case, the probability of each is set to realmin */
  /* 'convolv_2invG_nov:36' if max(t)<=0 */
  ixstart = 1;
  b_n = t->size[0];
  kd = t->data[0];
  if (t->size[0] > 1) {
    if (rtIsNaN(t->data[0])) {
      nm1d2 = 2;
      exitg1 = false;
      while ((!exitg1) && (nm1d2 <= b_n)) {
        ixstart = nm1d2;
        if (!rtIsNaN(t->data[nm1d2 - 1])) {
          kd = t->data[nm1d2 - 1];
          exitg1 = true;
        } else {
          nm1d2++;
        }
      }
    }

    if (ixstart < t->size[0]) {
      while (ixstart + 1 <= b_n) {
        if (t->data[ixstart] > kd) {
          kd = t->data[ixstart];
        }

        ixstart++;
      }
    }
  }

  if (kd <= 0.0) {
    /* 'convolv_2invG_nov:37' P=realmin.*ones(size(t)); */
    y[0] = t->size[0];
    nm1d2 = P->size[0];
    P->size[0] = (int)y[0];
    emxEnsureCapacity((emxArray__common *)P, nm1d2, sizeof(double));
    ixstart = (int)y[0];
    for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
      P->data[nm1d2] = 2.2250738585072014E-308;
    }
  } else {
    /* 'convolv_2invG_nov:38' else */
    /* 'convolv_2invG_nov:41' Maxt=max(t); */
    ixstart = 1;
    b_n = t->size[0];
    kd = t->data[0];
    if (t->size[0] > 1) {
      if (rtIsNaN(t->data[0])) {
        nm1d2 = 2;
        exitg1 = false;
        while ((!exitg1) && (nm1d2 <= b_n)) {
          ixstart = nm1d2;
          if (!rtIsNaN(t->data[nm1d2 - 1])) {
            kd = t->data[nm1d2 - 1];
            exitg1 = true;
          } else {
            nm1d2++;
          }
        }
      }

      if (ixstart < t->size[0]) {
        while (ixstart + 1 <= b_n) {
          if (t->data[ixstart] > kd) {
            kd = t->data[ixstart];
          }

          ixstart++;
        }
      }
    }

    /* 'convolv_2invG_nov:43' x=0:h:Maxt; */
    emxInit_real_T(&x, 2);
    if (rtIsNaN(h) || rtIsNaN(kd)) {
      nm1d2 = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)x, nm1d2, sizeof(double));
      x->data[0] = rtNaN;
    } else if ((h == 0.0) || ((0.0 < kd) && (h < 0.0)) || ((kd < 0.0) && (h >
                 0.0))) {
      nm1d2 = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)x, nm1d2, sizeof(double));
    } else if (rtIsInf(kd) && (rtIsInf(h) || (0.0 == kd))) {
      nm1d2 = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)x, nm1d2, sizeof(double));
      x->data[0] = rtNaN;
    } else if (rtIsInf(h)) {
      nm1d2 = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)x, nm1d2, sizeof(double));
      x->data[0] = 0.0;
    } else if (floor(h) == h) {
      nm1d2 = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = (int)floor(kd / h) + 1;
      emxEnsureCapacity((emxArray__common *)x, nm1d2, sizeof(double));
      ixstart = (int)floor(kd / h);
      for (nm1d2 = 0; nm1d2 <= ixstart; nm1d2++) {
        x->data[x->size[0] * nm1d2] = h * (double)nm1d2;
      }
    } else {
      ndbl = floor(kd / h + 0.5);
      apnd = ndbl * h;
      if (h > 0.0) {
        cdiff = apnd - kd;
      } else {
        cdiff = kd - apnd;
      }

      if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(kd)) {
        ndbl++;
        apnd = kd;
      } else if (cdiff > 0.0) {
        apnd = (ndbl - 1.0) * h;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        b_n = (int)ndbl;
      } else {
        b_n = 0;
      }

      nm1d2 = x->size[0] * x->size[1];
      x->size[0] = 1;
      x->size[1] = b_n;
      emxEnsureCapacity((emxArray__common *)x, nm1d2, sizeof(double));
      if (b_n > 0) {
        x->data[0] = 0.0;
        if (b_n > 1) {
          x->data[b_n - 1] = apnd;
          nm1d2 = (b_n - 1) / 2;
          for (ixstart = 1; ixstart < nm1d2; ixstart++) {
            kd = (double)ixstart * h;
            x->data[ixstart] = kd;
            x->data[(b_n - ixstart) - 1] = apnd - kd;
          }

          if (nm1d2 << 1 == b_n - 1) {
            x->data[nm1d2] = apnd / 2.0;
          } else {
            kd = (double)nm1d2 * h;
            x->data[nm1d2] = kd;
            x->data[nm1d2 + 1] = apnd - kd;
          }
        }
      }
    }

    emxInit_real_T(&b_y, 2);
    emxInit_real_T(&z, 2);

    /* 'convolv_2invG_nov:45' y=onestagepdf2(x,m(1),s(1)); */
    d_onestagepdf2(x, m[0], s[0], b_y);

    /* 'convolv_2invG_nov:46' z=onestagepdf2(x,m(2),s(2)); */
    d_onestagepdf2(x, m[1], s[1], z);

    /* 'convolv_2invG_nov:48' if sd(1)<.01 */
    emxInit_real_T(&C, 2);
    emxInit_real_T(&I, 2);
    emxInit_real_T(&b_P, 2);
    emxInit_real_T(&varargin_1, 2);
    emxInit_real_T1(&r3, 1);
    if (sd[0] < 0.01) {
      /* to estimate the error (in calculating probabilities the of the data),  */
      /* that results from approximating the first pdf as a point-mass */
      /* distribtion, find the maximum of the absolute value of the  */
      /* derivative of the second pdf. */
      /* 'convolv_2invG_nov:55' gp=gp_max(m(2),s(2)); */
      kd = gp_max(m[1], s[1]);

      /* determine the radius, r, of a small interval over which the second pdf, g, is */
      /* approximately constant and so that g(t) is small for t<r.   */
      /* 'convolv_2invG_nov:60' r=min(eps/(3*gp),T2); */
      kd = 2.2204460492503131E-16 / (3.0 * kd);
      if ((kd < T2) || rtIsNaN(T2)) {
        r = kd;
      } else {
        r = T2;
      }

      /* 'convolv_2invG_nov:61' checkval=onestagepdf2(r,m(2),s(2)); */
      kd = c_onestagepdf2(r, m[1], s[1]);

      /* 'convolv_2invG_nov:62' while checkval>=eps/2 */
      while (kd >= 1.1102230246251565E-16) {
        /* 'convolv_2invG_nov:63' r=r/2; */
        r /= 2.0;

        /* 'convolv_2invG_nov:64' checkval=onestagepdf2(r,m(2),s(2)); */
        kd = c_onestagepdf2(r, m[1], s[1]);
      }

      /* get the average value of the first pdf.  This is the point at which */
      /* its mass is concentrated. */
      /* 'convolv_2invG_nov:69' nu=1/m(1); */
      nu = 1.0 / m[0];

      /* get the maximum value f the second pdf. */
      /* 'convolv_2invG_nov:72' gm=onestagepdf2(T2,m(2),s(2)); */
      /* get upper limit of integral for approximating */
      /* int_{r+nu}^{infty}f(s)ds. */
      /* 'convolv_2invG_nov:75' Tu=max(100,nu+r+1000*sd(1)); */
      kd = (nu + r) + 1000.0 * sd[0];
      if ((100.0 > kd) || rtIsNaN(kd)) {
        Tu = 100.0;
      } else {
        Tu = kd;
      }

      /* 'convolv_2invG_nov:78' checkerror=100; */
      checkerror = 100.0;

      /* 'convolv_2invG_nov:79' hh=.001; */
      hh = 0.001;

      /* 'convolv_2invG_nov:80' check1=.001*sum(onestagepdf2((0:.001:nu-r),m(1),s(1)))+.001*sum(onestagepdf2((nu+r:.001:Tu),m(1),s(1))); */
      kd = nu - r;
      if (rtIsNaN(kd)) {
        nm1d2 = I->size[0] * I->size[1];
        I->size[0] = 1;
        I->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
        I->data[0] = rtNaN;
      } else if (kd < 0.0) {
        nm1d2 = I->size[0] * I->size[1];
        I->size[0] = 1;
        I->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
      } else if (rtIsInf(kd) && (0.0 == kd)) {
        nm1d2 = I->size[0] * I->size[1];
        I->size[0] = 1;
        I->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
        I->data[0] = rtNaN;
      } else {
        ndbl = floor(kd / 0.001 + 0.5);
        apnd = ndbl * 0.001;
        cdiff = apnd - kd;
        if (fabs(cdiff) < 4.4408920985006262E-16 * kd) {
          ndbl++;
          apnd = kd;
        } else if (cdiff > 0.0) {
          apnd = (ndbl - 1.0) * 0.001;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          b_n = (int)ndbl;
        } else {
          b_n = 0;
        }

        nm1d2 = I->size[0] * I->size[1];
        I->size[0] = 1;
        I->size[1] = b_n;
        emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
        if (b_n > 0) {
          I->data[0] = 0.0;
          if (b_n > 1) {
            I->data[b_n - 1] = apnd;
            nm1d2 = (b_n - 1) / 2;
            for (ixstart = 1; ixstart < nm1d2; ixstart++) {
              kd = (double)ixstart * 0.001;
              I->data[ixstart] = kd;
              I->data[(b_n - ixstart) - 1] = apnd - kd;
            }

            if (nm1d2 << 1 == b_n - 1) {
              I->data[nm1d2] = apnd / 2.0;
            } else {
              kd = (double)nm1d2 * 0.001;
              I->data[nm1d2] = kd;
              I->data[nm1d2 + 1] = apnd - kd;
            }
          }
        }
      }

      a = nu + r;
      if (rtIsNaN(a) || rtIsNaN(Tu)) {
        nm1d2 = b_P->size[0] * b_P->size[1];
        b_P->size[0] = 1;
        b_P->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
        b_P->data[0] = rtNaN;
      } else if (Tu < a) {
        nm1d2 = b_P->size[0] * b_P->size[1];
        b_P->size[0] = 1;
        b_P->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
      } else if ((rtIsInf(a) || rtIsInf(Tu)) && (a == Tu)) {
        nm1d2 = b_P->size[0] * b_P->size[1];
        b_P->size[0] = 1;
        b_P->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
        b_P->data[0] = rtNaN;
      } else {
        ndbl = floor((Tu - a) / 0.001 + 0.5);
        apnd = a + ndbl * 0.001;
        cdiff = apnd - Tu;
        kd = fabs(a);
        if (!((kd > Tu) || rtIsNaN(Tu))) {
          kd = Tu;
        }

        if (fabs(cdiff) < 4.4408920985006262E-16 * kd) {
          ndbl++;
          apnd = Tu;
        } else if (cdiff > 0.0) {
          apnd = a + (ndbl - 1.0) * 0.001;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          b_n = (int)ndbl;
        } else {
          b_n = 0;
        }

        nm1d2 = b_P->size[0] * b_P->size[1];
        b_P->size[0] = 1;
        b_P->size[1] = b_n;
        emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
        if (b_n > 0) {
          b_P->data[0] = a;
          if (b_n > 1) {
            b_P->data[b_n - 1] = apnd;
            nm1d2 = (b_n - 1) / 2;
            for (ixstart = 1; ixstart < nm1d2; ixstart++) {
              kd = (double)ixstart * 0.001;
              b_P->data[ixstart] = a + kd;
              b_P->data[(b_n - ixstart) - 1] = apnd - kd;
            }

            if (nm1d2 << 1 == b_n - 1) {
              b_P->data[nm1d2] = (a + apnd) / 2.0;
            } else {
              kd = (double)nm1d2 * 0.001;
              b_P->data[nm1d2] = a + kd;
              b_P->data[nm1d2 + 1] = apnd - kd;
            }
          }
        }
      }

      d_onestagepdf2(I, m[0], s[0], varargin_1);
      d_onestagepdf2(b_P, m[0], s[0], I);
      check1 = 0.001 * b_sum(varargin_1) + 0.001 * b_sum(I);

      /* 'convolv_2invG_nov:81' while checkerror>10^(-4) */
      while (checkerror > 0.0001) {
        /* 'convolv_2invG_nov:82' hh=.5*hh; */
        hh *= 0.5;

        /* 'convolv_2invG_nov:83' ck1=hh*sum(onestagepdf2((0:hh:nu-r),m(1),s(1)))+hh*sum(onestagepdf2((nu+r:hh:Tu),m(1),s(1))); */
        kd = nu - r;
        if (rtIsNaN(kd)) {
          nm1d2 = I->size[0] * I->size[1];
          I->size[0] = 1;
          I->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
          I->data[0] = rtNaN;
        } else if ((hh == 0.0) || ((0.0 < kd) && (hh < 0.0)) || ((kd < 0.0) &&
                    (hh > 0.0))) {
          nm1d2 = I->size[0] * I->size[1];
          I->size[0] = 1;
          I->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
        } else if (rtIsInf(kd) && (rtIsInf(hh) || (0.0 == kd))) {
          nm1d2 = I->size[0] * I->size[1];
          I->size[0] = 1;
          I->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
          I->data[0] = rtNaN;
        } else if (rtIsInf(hh)) {
          nm1d2 = I->size[0] * I->size[1];
          I->size[0] = 1;
          I->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
          I->data[0] = 0.0;
        } else if (floor(hh) == hh) {
          nm1d2 = I->size[0] * I->size[1];
          I->size[0] = 1;
          I->size[1] = (int)floor(kd / hh) + 1;
          emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
          ixstart = (int)floor(kd / hh);
          for (nm1d2 = 0; nm1d2 <= ixstart; nm1d2++) {
            I->data[I->size[0] * nm1d2] = hh * (double)nm1d2;
          }
        } else {
          ndbl = floor(kd / hh + 0.5);
          apnd = ndbl * hh;
          if (hh > 0.0) {
            cdiff = apnd - kd;
          } else {
            cdiff = kd - apnd;
          }

          if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(kd)) {
            ndbl++;
            apnd = kd;
          } else if (cdiff > 0.0) {
            apnd = (ndbl - 1.0) * hh;
          } else {
            ndbl++;
          }

          if (ndbl >= 0.0) {
            b_n = (int)ndbl;
          } else {
            b_n = 0;
          }

          nm1d2 = I->size[0] * I->size[1];
          I->size[0] = 1;
          I->size[1] = b_n;
          emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
          if (b_n > 0) {
            I->data[0] = 0.0;
            if (b_n > 1) {
              I->data[b_n - 1] = apnd;
              nm1d2 = (b_n - 1) / 2;
              for (ixstart = 1; ixstart < nm1d2; ixstart++) {
                kd = (double)ixstart * hh;
                I->data[ixstart] = kd;
                I->data[(b_n - ixstart) - 1] = apnd - kd;
              }

              if (nm1d2 << 1 == b_n - 1) {
                I->data[nm1d2] = apnd / 2.0;
              } else {
                kd = (double)nm1d2 * hh;
                I->data[nm1d2] = kd;
                I->data[nm1d2 + 1] = apnd - kd;
              }
            }
          }
        }

        a = nu + r;
        if (rtIsNaN(a) || rtIsNaN(Tu)) {
          nm1d2 = b_P->size[0] * b_P->size[1];
          b_P->size[0] = 1;
          b_P->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
          b_P->data[0] = rtNaN;
        } else if ((hh == 0.0) || ((a < Tu) && (hh < 0.0)) || ((Tu < a) && (hh >
          0.0))) {
          nm1d2 = b_P->size[0] * b_P->size[1];
          b_P->size[0] = 1;
          b_P->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
        } else if ((rtIsInf(a) || rtIsInf(Tu)) && (rtIsInf(hh) || (a == Tu))) {
          nm1d2 = b_P->size[0] * b_P->size[1];
          b_P->size[0] = 1;
          b_P->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
          b_P->data[0] = rtNaN;
        } else if (rtIsInf(hh)) {
          nm1d2 = b_P->size[0] * b_P->size[1];
          b_P->size[0] = 1;
          b_P->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
          b_P->data[0] = a;
        } else if ((floor(a) == a) && (floor(hh) == hh)) {
          nm1d2 = b_P->size[0] * b_P->size[1];
          b_P->size[0] = 1;
          b_P->size[1] = (int)floor((Tu - a) / hh) + 1;
          emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
          ixstart = (int)floor((Tu - a) / hh);
          for (nm1d2 = 0; nm1d2 <= ixstart; nm1d2++) {
            b_P->data[b_P->size[0] * nm1d2] = a + hh * (double)nm1d2;
          }
        } else {
          ndbl = floor((Tu - a) / hh + 0.5);
          apnd = a + ndbl * hh;
          if (hh > 0.0) {
            cdiff = apnd - Tu;
          } else {
            cdiff = Tu - apnd;
          }

          kd = fabs(a);
          if (!((kd > Tu) || rtIsNaN(Tu))) {
            kd = Tu;
          }

          if (fabs(cdiff) < 4.4408920985006262E-16 * kd) {
            ndbl++;
            apnd = Tu;
          } else if (cdiff > 0.0) {
            apnd = a + (ndbl - 1.0) * hh;
          } else {
            ndbl++;
          }

          if (ndbl >= 0.0) {
            b_n = (int)ndbl;
          } else {
            b_n = 0;
          }

          nm1d2 = b_P->size[0] * b_P->size[1];
          b_P->size[0] = 1;
          b_P->size[1] = b_n;
          emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
          if (b_n > 0) {
            b_P->data[0] = a;
            if (b_n > 1) {
              b_P->data[b_n - 1] = apnd;
              nm1d2 = (b_n - 1) / 2;
              for (ixstart = 1; ixstart < nm1d2; ixstart++) {
                kd = (double)ixstart * hh;
                b_P->data[ixstart] = a + kd;
                b_P->data[(b_n - ixstart) - 1] = apnd - kd;
              }

              if (nm1d2 << 1 == b_n - 1) {
                b_P->data[nm1d2] = (a + apnd) / 2.0;
              } else {
                kd = (double)nm1d2 * hh;
                b_P->data[nm1d2] = a + kd;
                b_P->data[nm1d2 + 1] = apnd - kd;
              }
            }
          }
        }

        d_onestagepdf2(I, m[0], s[0], varargin_1);
        d_onestagepdf2(b_P, m[0], s[0], I);
        kd = hh * b_sum(varargin_1) + hh * b_sum(I);

        /* 'convolv_2invG_nov:84' checkerror=abs(check1-ck1); */
        checkerror = fabs(check1 - kd);

        /* 'convolv_2invG_nov:85' check1=ck1; */
        check1 = kd;
      }

      /* 'convolv_2invG_nov:87' check2=gm*check1; */
      /* 'convolv_2invG_nov:88' if  check2<=eps/3 */
      if (c_onestagepdf2(T2, m[1], s[1]) * check1 <= 7.4014868308343765E-17) {
        emxInit_real_T1(&b_t, 1);

        /* 'convolv_2invG_nov:90' flag=1; */
        b_flag = 1;

        /* 'convolv_2invG_nov:91' sigma=s(2); */
        /* 'convolv_2invG_nov:92' mu=m(2); */
        /* 'convolv_2invG_nov:93' l=1/m(1); */
        kd = 1.0 / m[0];

        /* 'convolv_2invG_nov:95' P=onestagepdf_lag(t,mu,sigma,l); */
        /* this function evaluates a shifted inverse gaussian distribution at X */
        /* m=mu */
        /* s=sigma */
        /* l=lag */
        /* 'onestagepdf_lag:7' Y=(1./(s*(2*pi*(X-l).^3).^(1/2))).*exp(-(m*(X-l)-1).^2./(2*s^2*(X-l))); */
        nm1d2 = b_t->size[0];
        b_t->size[0] = t->size[0];
        emxEnsureCapacity((emxArray__common *)b_t, nm1d2, sizeof(double));
        ixstart = t->size[0];
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          b_t->data[nm1d2] = t->data[nm1d2] - kd;
        }

        emxInit_real_T1(&b, 1);
        emxInit_real_T1(&r4, 1);
        k_power(b_t, b);
        nm1d2 = r4->size[0];
        r4->size[0] = b->size[0];
        emxEnsureCapacity((emxArray__common *)r4, nm1d2, sizeof(double));
        ixstart = b->size[0];
        emxFree_real_T(&b_t);
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          r4->data[nm1d2] = 6.2831853071795862 * b->data[nm1d2];
        }

        emxInit_real_T1(&b_m, 1);
        l_power(r4, b);
        a = 2.0 * (s[1] * s[1]);
        nm1d2 = b_m->size[0];
        b_m->size[0] = t->size[0];
        emxEnsureCapacity((emxArray__common *)b_m, nm1d2, sizeof(double));
        ixstart = t->size[0];
        emxFree_real_T(&r4);
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          b_m->data[nm1d2] = m[1] * (t->data[nm1d2] - kd) - 1.0;
        }

        emxInit_real_T1(&r5, 1);
        m_power(b_m, r3);
        nm1d2 = r5->size[0];
        r5->size[0] = r3->size[0];
        emxEnsureCapacity((emxArray__common *)r5, nm1d2, sizeof(double));
        ixstart = r3->size[0];
        emxFree_real_T(&b_m);
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          r5->data[nm1d2] = -r3->data[nm1d2];
        }

        emxInit_real_T1(&b_a, 1);
        nm1d2 = b_a->size[0];
        b_a->size[0] = t->size[0];
        emxEnsureCapacity((emxArray__common *)b_a, nm1d2, sizeof(double));
        ixstart = t->size[0];
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          b_a->data[nm1d2] = a * (t->data[nm1d2] - kd);
        }

        emxInit_real_T1(&c_s, 1);
        c_rdivide(r5, b_a, r3);
        b_exp(r3);
        nm1d2 = c_s->size[0];
        c_s->size[0] = b->size[0];
        emxEnsureCapacity((emxArray__common *)c_s, nm1d2, sizeof(double));
        ixstart = b->size[0];
        emxFree_real_T(&b_a);
        emxFree_real_T(&r5);
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          c_s->data[nm1d2] = s[1] * b->data[nm1d2];
        }

        b_rdivide(c_s, P);
        nm1d2 = P->size[0];
        emxEnsureCapacity((emxArray__common *)P, nm1d2, sizeof(double));
        ixstart = P->size[0];
        emxFree_real_T(&c_s);
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          P->data[nm1d2] *= r3->data[nm1d2];
        }

        /* If some cpmponent of X is less than l, the pdf will return imaginary values,  */
        /* in this case we stipilate that the pdf returns realmin */
        /* 'onestagepdf_lag:10' Y=real(Y); */
        /* 'onestagepdf_lag:11' Y=max(Y, realmin); */
        nm1d2 = b->size[0];
        b->size[0] = P->size[0];
        emxEnsureCapacity((emxArray__common *)b, nm1d2, sizeof(double));
        ixstart = P->size[0];
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          b->data[nm1d2] = P->data[nm1d2];
        }

        P_idx_0 = (unsigned int)P->size[0];
        c_P = (unsigned int)P->size[0];
        nm1d2 = P->size[0];
        P->size[0] = (int)c_P;
        emxEnsureCapacity((emxArray__common *)P, nm1d2, sizeof(double));
        for (ixstart = 0; ixstart + 1 <= (int)P_idx_0; ixstart++) {
          kd = b->data[ixstart];
          if (!(kd > 2.2250738585072014E-308)) {
            kd = 2.2250738585072014E-308;
          }

          P->data[ixstart] = kd;
        }

        emxFree_real_T(&b);
      } else {
        /* 'convolv_2invG_nov:97' else */
        /*  find the discrete convolution of the vectors y and z */
        /*  the (i-1)th element of v approximates the convolution of the pdfs  */
        /*  over [.001, x(i)] as a left-hand Riemann sum. */
        /* 'convolv_2invG_nov:102' C=conv(z,y)*h; */
        b_conv(z, b_y, C);
        nm1d2 = C->size[0] * C->size[1];
        C->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)C, nm1d2, sizeof(double));
        ixstart = C->size[0];
        nm1d2 = C->size[1];
        ixstart *= nm1d2;
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          C->data[nm1d2] *= h;
        }

        /* 'convolv_2invG_nov:103' N=length(y); */
        /*  only the first N elements of the convolution are valid */
        /* 'convolv_2invG_nov:105' C=C(1:N); */
        if (1 > b_y->size[1]) {
          nm1d2 = 0;
        } else {
          nm1d2 = b_y->size[1];
        }

        ixstart = C->size[0] * C->size[1];
        C->size[1] = nm1d2;
        emxEnsureCapacity((emxArray__common *)C, ixstart, sizeof(double));

        /* 'convolv_2invG_nov:106' I=zeros(1,n); */
        /* 'convolv_2invG_nov:107' P=zeros(1,n); */
        /* 'convolv_2invG_nov:108' for i=1:n */
        nm1d2 = I->size[0] * I->size[1];
        I->size[0] = 1;
        I->size[1] = t->size[0];
        emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
        nm1d2 = b_P->size[0] * b_P->size[1];
        b_P->size[0] = 1;
        b_P->size[1] = t->size[0];
        emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
        i = 0;
        emxInit_real_T(&b_t, 2);
        while (i <= n) {
          /* find element of x that is closest to t(i) */
          /* 'convolv_2invG_nov:110' [~,I(i)]=min((t(i)-x).^2); */
          nm1d2 = b_t->size[0] * b_t->size[1];
          b_t->size[0] = 1;
          b_t->size[1] = x->size[1];
          emxEnsureCapacity((emxArray__common *)b_t, nm1d2, sizeof(double));
          kd = t->data[i];
          ixstart = x->size[0] * x->size[1];
          for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
            b_t->data[nm1d2] = kd - x->data[nm1d2];
          }

          e_power(b_t, varargin_1);
          ixstart = 1;
          b_n = varargin_1->size[1];
          kd = varargin_1->data[0];
          itmp = 1;
          if (varargin_1->size[1] > 1) {
            if (rtIsNaN(varargin_1->data[0])) {
              nm1d2 = 2;
              exitg1 = false;
              while ((!exitg1) && (nm1d2 <= b_n)) {
                ixstart = nm1d2;
                if (!rtIsNaN(varargin_1->data[nm1d2 - 1])) {
                  kd = varargin_1->data[nm1d2 - 1];
                  itmp = nm1d2;
                  exitg1 = true;
                } else {
                  nm1d2++;
                }
              }
            }

            if (ixstart < varargin_1->size[1]) {
              while (ixstart + 1 <= b_n) {
                if (varargin_1->data[ixstart] < kd) {
                  kd = varargin_1->data[ixstart];
                  itmp = ixstart + 1;
                }

                ixstart++;
              }
            }
          }

          I->data[i] = itmp;

          /* 'convolv_2invG_nov:110' ~ */
          /* If t(i)<0 the probability is set to zero, otherwise the */
          /* probability is approxiated as a value from the vector x. */
          /* 'convolv_2invG_nov:113' if t(i)>0 && I(i)>1 */
          if ((t->data[i] > 0.0) && ((int)(unsigned int)I->data[i] > 1)) {
            /* 'convolv_2invG_nov:114' P(i)=C(I(i)-1); */
            b_P->data[i] = C->data[(int)(unsigned int)I->data[i] - 2];
          } else {
            /* 'convolv_2invG_nov:115' else */
            /* 'convolv_2invG_nov:116' P(i)=realmin; */
            b_P->data[i] = 2.2250738585072014E-308;
          }

          i++;
        }

        emxFree_real_T(&b_t);

        /* toc */
        /* 'convolv_2invG_nov:121' P=P'; */
        nm1d2 = P->size[0];
        P->size[0] = b_P->size[1];
        emxEnsureCapacity((emxArray__common *)P, nm1d2, sizeof(double));
        ixstart = b_P->size[1];
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          P->data[nm1d2] = b_P->data[b_P->size[0] * nm1d2];
        }

        /* 'convolv_2invG_nov:123' logP=sum(log(P)); */
        nm1d2 = r3->size[0];
        r3->size[0] = P->size[0];
        emxEnsureCapacity((emxArray__common *)r3, nm1d2, sizeof(double));
        ixstart = P->size[0];
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          r3->data[nm1d2] = P->data[nm1d2];
        }

        c_log(r3);
      }
    } else {
      /* 'convolv_2invG_nov:126' else */
      /*  find the discrete convolution of the vectors y and z */
      /*  the (i-1)th element of v approximates the convolution of the pdfs  */
      /*  over [.001, x(i)] as a left-hand Riemann sum. */
      /* 'convolv_2invG_nov:131' C=conv(z,y)*h; */
      b_conv(z, b_y, C);
      nm1d2 = C->size[0] * C->size[1];
      C->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)C, nm1d2, sizeof(double));
      ixstart = C->size[0];
      nm1d2 = C->size[1];
      ixstart *= nm1d2;
      for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
        C->data[nm1d2] *= h;
      }

      /* 'convolv_2invG_nov:132' N=length(y); */
      /*  only the first N elements of the convolution are valid */
      /* 'convolv_2invG_nov:134' C=C(1:N); */
      if (1 > b_y->size[1]) {
        nm1d2 = 0;
      } else {
        nm1d2 = b_y->size[1];
      }

      ixstart = C->size[0] * C->size[1];
      C->size[1] = nm1d2;
      emxEnsureCapacity((emxArray__common *)C, ixstart, sizeof(double));

      /* 'convolv_2invG_nov:135' I=zeros(1,n); */
      /* 'convolv_2invG_nov:136' P=zeros(1,n); */
      /* 'convolv_2invG_nov:137' for i=1:n */
      nm1d2 = I->size[0] * I->size[1];
      I->size[0] = 1;
      I->size[1] = t->size[0];
      emxEnsureCapacity((emxArray__common *)I, nm1d2, sizeof(double));
      nm1d2 = b_P->size[0] * b_P->size[1];
      b_P->size[0] = 1;
      b_P->size[1] = t->size[0];
      emxEnsureCapacity((emxArray__common *)b_P, nm1d2, sizeof(double));
      i = 0;
      emxInit_real_T(&b_t, 2);
      while (i <= n) {
        /* find element of x that is closest to t(i) */
        /* 'convolv_2invG_nov:139' [~,I(i)]=min((t(i)-x).^2); */
        nm1d2 = b_t->size[0] * b_t->size[1];
        b_t->size[0] = 1;
        b_t->size[1] = x->size[1];
        emxEnsureCapacity((emxArray__common *)b_t, nm1d2, sizeof(double));
        kd = t->data[i];
        ixstart = x->size[0] * x->size[1];
        for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
          b_t->data[nm1d2] = kd - x->data[nm1d2];
        }

        e_power(b_t, varargin_1);
        ixstart = 1;
        b_n = varargin_1->size[1];
        kd = varargin_1->data[0];
        itmp = 1;
        if (varargin_1->size[1] > 1) {
          if (rtIsNaN(varargin_1->data[0])) {
            nm1d2 = 2;
            exitg1 = false;
            while ((!exitg1) && (nm1d2 <= b_n)) {
              ixstart = nm1d2;
              if (!rtIsNaN(varargin_1->data[nm1d2 - 1])) {
                kd = varargin_1->data[nm1d2 - 1];
                itmp = nm1d2;
                exitg1 = true;
              } else {
                nm1d2++;
              }
            }
          }

          if (ixstart < varargin_1->size[1]) {
            while (ixstart + 1 <= b_n) {
              if (varargin_1->data[ixstart] < kd) {
                kd = varargin_1->data[ixstart];
                itmp = ixstart + 1;
              }

              ixstart++;
            }
          }
        }

        I->data[i] = itmp;

        /* 'convolv_2invG_nov:139' ~ */
        /* If t(i)<0 the probability is set to zero, otherwise the */
        /* probability is approxiated as a value from the vector x. */
        /* 'convolv_2invG_nov:142' if t(i)>0 && I(i)>1 */
        if ((t->data[i] > 0.0) && ((int)(unsigned int)I->data[i] > 1)) {
          /* 'convolv_2invG_nov:143' P(i)=C(I(i)-1); */
          b_P->data[i] = C->data[(int)(unsigned int)I->data[i] - 2];
        } else {
          /* 'convolv_2invG_nov:144' else */
          /* 'convolv_2invG_nov:145' P(i)=realmin; */
          b_P->data[i] = 2.2250738585072014E-308;
        }

        i++;
      }

      emxFree_real_T(&b_t);

      /* toc */
      /* 'convolv_2invG_nov:150' P=P'; */
      nm1d2 = P->size[0];
      P->size[0] = b_P->size[1];
      emxEnsureCapacity((emxArray__common *)P, nm1d2, sizeof(double));
      ixstart = b_P->size[1];
      for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
        P->data[nm1d2] = b_P->data[b_P->size[0] * nm1d2];
      }

      /* 'convolv_2invG_nov:152' logP=sum(log(P)); */
      nm1d2 = r3->size[0];
      r3->size[0] = P->size[0];
      emxEnsureCapacity((emxArray__common *)r3, nm1d2, sizeof(double));
      ixstart = P->size[0];
      for (nm1d2 = 0; nm1d2 < ixstart; nm1d2++) {
        r3->data[nm1d2] = P->data[nm1d2];
      }

      c_log(r3);
    }

    emxFree_real_T(&r3);
    emxFree_real_T(&varargin_1);
    emxFree_real_T(&b_P);
    emxFree_real_T(&I);
    emxFree_real_T(&C);
    emxFree_real_T(&z);
    emxFree_real_T(&b_y);
    emxFree_real_T(&x);
  }

  *flag = b_flag;
}

/*
 * function [P,flag]=convolv_2invG_nov(t,m1,s1,m2,s2,h)
 */
void c_convolv_2invG_nov(double m1, double s1, double m2, double s2, double P
  [2201], double *flag)
{
  double y[2];
  int ixstart;
  double m[2];
  double s[2];
  double v[2];
  int nm1d2;
  int iidx[2];
  double sd[2];
  double b_s[2];
  double T2;
  static const double x[2201] = { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07,
    0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2,
    0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33,
    0.34, 0.35000000000000003, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41000000000000003,
    0.42, 0.43, 0.44, 0.45, 0.46, 0.47000000000000003, 0.48, 0.49, 0.5, 0.51,
    0.52, 0.53, 0.54, 0.55, 0.56, 0.57000000000000006, 0.58, 0.59, 0.6, 0.61,
    0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69000000000000006,
    0.70000000000000007, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79,
    0.8, 0.81, 0.82000000000000006, 0.83000000000000007, 0.84, 0.85, 0.86, 0.87,
    0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94000000000000006, 0.95000000000000007,
    0.96, 0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08,
    1.09, 1.1, 1.11, 1.12, 1.1300000000000001, 1.1400000000000001,
    1.1500000000000001, 1.16, 1.17, 1.18, 1.19, 1.2, 1.21, 1.22, 1.23, 1.24,
    1.25, 1.26, 1.27, 1.28, 1.29, 1.3, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37,
    1.3800000000000001, 1.3900000000000001, 1.4000000000000001, 1.41, 1.42, 1.43,
    1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.5, 1.51, 1.52, 1.53, 1.54, 1.55, 1.56,
    1.57, 1.58, 1.59, 1.6, 1.61, 1.62, 1.6300000000000001, 1.6400000000000001,
    1.6500000000000001, 1.6600000000000001, 1.67, 1.68, 1.69, 1.7, 1.71, 1.72,
    1.73, 1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.8, 1.81, 1.82, 1.83, 1.84, 1.85,
    1.86, 1.87, 1.8800000000000001, 1.8900000000000001, 1.9000000000000001,
    1.9100000000000001, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 2.0,
    2.0100000000000002, 2.02, 2.0300000000000002, 2.04, 2.05, 2.06, 2.07, 2.08,
    2.09, 2.1, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.2, 2.21,
    2.22, 2.23, 2.24, 2.25, 2.2600000000000002, 2.27, 2.2800000000000002, 2.29,
    2.3000000000000003, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39,
    2.4, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, 2.5,
    2.5100000000000002, 2.52, 2.5300000000000002, 2.54, 2.5500000000000003, 2.56,
    2.57, 2.58, 2.59, 2.6, 2.61, 2.62, 2.63, 2.64, 2.65, 2.66, 2.67, 2.68, 2.69,
    2.7, 2.71, 2.72, 2.73, 2.74, 2.75, 2.7600000000000002, 2.77,
    2.7800000000000002, 2.79, 2.8000000000000003, 2.81, 2.82, 2.83, 2.84, 2.85,
    2.86, 2.87, 2.88, 2.89, 2.9, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98,
    2.99, 3.0, 3.0100000000000002, 3.02, 3.0300000000000002, 3.04,
    3.0500000000000003, 3.06, 3.0700000000000003, 3.08, 3.09, 3.1, 3.11, 3.12,
    3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.2, 3.21, 3.22, 3.23, 3.24, 3.25,
    3.2600000000000002, 3.27, 3.2800000000000002, 3.29, 3.3000000000000003, 3.31,
    3.3200000000000003, 3.33, 3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 3.4, 3.41,
    3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.48, 3.49, 3.5, 3.5100000000000002,
    3.52, 3.5300000000000002, 3.54, 3.5500000000000003, 3.56, 3.5700000000000003,
    3.58, 3.59, 3.6, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.7,
    3.71, 3.72, 3.73, 3.74, 3.75, 3.7600000000000002, 3.77, 3.7800000000000002,
    3.79, 3.8000000000000003, 3.81, 3.8200000000000003, 3.83, 3.84, 3.85, 3.86,
    3.87, 3.88, 3.89, 3.9, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99,
    4.0, 4.01, 4.0200000000000005, 4.03, 4.04, 4.05, 4.0600000000000005, 4.07,
    4.08, 4.09, 4.1, 4.11, 4.12, 4.13, 4.14, 4.15, 4.16, 4.17, 4.18, 4.19, 4.2,
    4.21, 4.22, 4.23, 4.24, 4.25, 4.26, 4.2700000000000005, 4.28, 4.29, 4.3,
    4.3100000000000005, 4.32, 4.33, 4.34, 4.3500000000000005, 4.36, 4.37, 4.38,
    4.39, 4.4, 4.41, 4.42, 4.43, 4.44, 4.45, 4.46, 4.47, 4.48, 4.49, 4.5, 4.51,
    4.5200000000000005, 4.53, 4.54, 4.55, 4.5600000000000005, 4.57, 4.58, 4.59,
    4.6000000000000005, 4.61, 4.62, 4.63, 4.64, 4.65, 4.66, 4.67, 4.68, 4.69,
    4.7, 4.71, 4.72, 4.73, 4.74, 4.75, 4.76, 4.7700000000000005, 4.78, 4.79, 4.8,
    4.8100000000000005, 4.82, 4.83, 4.84, 4.8500000000000005, 4.86, 4.87, 4.88,
    4.89, 4.9, 4.91, 4.92, 4.93, 4.94, 4.95, 4.96, 4.97, 4.98, 4.99, 5.0, 5.01,
    5.0200000000000005, 5.03, 5.04, 5.05, 5.0600000000000005, 5.07, 5.08, 5.09,
    5.1000000000000005, 5.11, 5.12, 5.13, 5.14, 5.15, 5.16, 5.17, 5.18, 5.19,
    5.2, 5.21, 5.22, 5.23, 5.24, 5.25, 5.26, 5.2700000000000005, 5.28, 5.29, 5.3,
    5.3100000000000005, 5.32, 5.33, 5.34, 5.3500000000000005, 5.36, 5.37, 5.38,
    5.39, 5.4, 5.41, 5.42, 5.43, 5.44, 5.45, 5.46, 5.47, 5.48, 5.49, 5.5, 5.51,
    5.5200000000000005, 5.53, 5.54, 5.55, 5.5600000000000005, 5.57, 5.58, 5.59,
    5.6000000000000005, 5.61, 5.62, 5.63, 5.64, 5.65, 5.66, 5.67, 5.68, 5.69,
    5.7, 5.71, 5.72, 5.73, 5.74, 5.75, 5.76, 5.7700000000000005, 5.78, 5.79, 5.8,
    5.8100000000000005, 5.82, 5.83, 5.84, 5.8500000000000005, 5.86, 5.87, 5.88,
    5.89, 5.9, 5.91, 5.92, 5.93, 5.94, 5.95, 5.96, 5.97, 5.98, 5.99, 6.0, 6.01,
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

  double b_y[2201];
  double z[2201];
  static double C[4401];
  double gp;
  double r;
  int i;
  double nu;
  double t[2201];
  double varargin_1[2201];
  static const double b_t[2201] = { 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
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

  int ix;
  boolean_T exitg1;
  double Tu;
  emxArray_real_T *c_y;
  double checkerror;
  double hh;
  double ndbl;
  double apnd;
  double cdiff;
  emxArray_real_T *d_y;
  double a;
  emxArray_real_T *r10;
  double check1;

  /* This function evaluates the convolution of two inverse Gaussian */
  /* distributions at vector t. */
  /* If the variance in one of the distributions is very small so that the  */
  /* distribution is close the a Dirac delta function, the convolution */
  /* is approximated as a shifted inverse Gaussian, that is  */
  /* one part of the cell cycle is treated as deterministic in length.   */
  /* In this case, the shift or lag is the average time to complete  */
  /* the part of the cycle with the smallest sigma. */
  /* t is a vector of times to divide (or times to complete two parts of the */
  /* cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1, */
  /* s2=sigma2. */
  /* 'convolv_2invG_nov:16' flag=0; */
  *flag = 0.0;

  /* 'convolv_2invG_nov:18' E=Inf; */
  /* 'convolv_2invG_nov:19' n=length(t); */
  /* 'convolv_2invG_nov:20' m=[m1 m2]; */
  /* 'convolv_2invG_nov:21' s=[s1 s2]; */
  /* 'convolv_2invG_nov:22' m=max(0,m); */
  y[0] = m1;
  y[1] = m2;
  for (ixstart = 0; ixstart < 2; ixstart++) {
    if ((0.0 > y[ixstart]) || rtIsNaN(y[ixstart])) {
      m[ixstart] = 0.0;
    } else {
      m[ixstart] = y[ixstart];
    }
  }

  /* 'convolv_2invG_nov:23' s=max(0,s); */
  y[0] = s1;
  y[1] = s2;
  for (ixstart = 0; ixstart < 2; ixstart++) {
    if ((0.0 > y[ixstart]) || rtIsNaN(y[ixstart])) {
      s[ixstart] = 0.0;
    } else {
      s[ixstart] = y[ixstart];
    }
  }

  /* 'convolv_2invG_nov:24' v=(s.^2)./(m.^3); */
  /* 'convolv_2invG_nov:25' [v,I]=sort(v); */
  power(s, v);
  b_power(m, y);
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    v[nm1d2] /= y[nm1d2];
  }

  sort(v, iidx);

  /* 'convolv_2invG_nov:26' sd=v.^.5; */
  c_power(v, sd);

  /* 'convolv_2invG_nov:27' m=m(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    v[nm1d2] = m[iidx[nm1d2] - 1];
    y[nm1d2] = iidx[nm1d2];
  }

  /* 'convolv_2invG_nov:28' s=s(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    m[nm1d2] = v[nm1d2];
    b_s[nm1d2] = s[(int)y[nm1d2] - 1];
  }

  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    s[nm1d2] = b_s[nm1d2];
  }

  /* 'convolv_2invG_nov:30' T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1))); */
  /* 'convolv_2invG_nov:31' T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2))); */
  T2 = 1.0 / m[1] * (sqrt(1.0 + 2.25 * (rt_powd_snf(s[1], 4.0) / (m[1] * m[1])))
                     - 1.5 * (s[1] * s[1] / m[1]));

  /* When called from convolv_3invG all of the components of t may be negative. */
  /* In this case, the probability of each is set to realmin */
  /* 'convolv_2invG_nov:36' if max(t)<=0 */
  /* 'convolv_2invG_nov:38' else */
  /* 'convolv_2invG_nov:41' Maxt=max(t); */
  /* 'convolv_2invG_nov:43' x=0:h:Maxt; */
  /* 'convolv_2invG_nov:45' y=onestagepdf2(x,m(1),s(1)); */
  b_onestagepdf2(x, m[0], s[0], b_y);

  /* 'convolv_2invG_nov:46' z=onestagepdf2(x,m(2),s(2)); */
  b_onestagepdf2(x, m[1], s[1], z);

  /* 'convolv_2invG_nov:48' if sd(1)<.01 */
  if (sd[0] < 0.01) {
    /* to estimate the error (in calculating probabilities the of the data),  */
    /* that results from approximating the first pdf as a point-mass */
    /* distribtion, find the maximum of the absolute value of the  */
    /* derivative of the second pdf. */
    /* 'convolv_2invG_nov:55' gp=gp_max(m(2),s(2)); */
    gp = gp_max(m[1], s[1]);

    /* determine the radius, r, of a small interval over which the second pdf, g, is */
    /* approximately constant and so that g(t) is small for t<r.   */
    /* 'convolv_2invG_nov:60' r=min(eps/(3*gp),T2); */
    gp = 2.2204460492503131E-16 / (3.0 * gp);
    if ((gp < T2) || rtIsNaN(T2)) {
      r = gp;
    } else {
      r = T2;
    }

    /* 'convolv_2invG_nov:61' checkval=onestagepdf2(r,m(2),s(2)); */
    gp = c_onestagepdf2(r, m[1], s[1]);

    /* 'convolv_2invG_nov:62' while checkval>=eps/2 */
    while (gp >= 1.1102230246251565E-16) {
      /* 'convolv_2invG_nov:63' r=r/2; */
      r /= 2.0;

      /* 'convolv_2invG_nov:64' checkval=onestagepdf2(r,m(2),s(2)); */
      gp = c_onestagepdf2(r, m[1], s[1]);
    }

    /* get the average value of the first pdf.  This is the point at which */
    /* its mass is concentrated. */
    /* 'convolv_2invG_nov:69' nu=1/m(1); */
    nu = 1.0 / m[0];

    /* get the maximum value f the second pdf. */
    /* 'convolv_2invG_nov:72' gm=onestagepdf2(T2,m(2),s(2)); */
    /* get upper limit of integral for approximating */
    /* int_{r+nu}^{infty}f(s)ds. */
    /* 'convolv_2invG_nov:75' Tu=max(100,nu+r+1000*sd(1)); */
    gp = (nu + r) + 1000.0 * sd[0];
    if ((100.0 > gp) || rtIsNaN(gp)) {
      Tu = 100.0;
    } else {
      Tu = gp;
    }

    emxInit_real_T(&c_y, 2);

    /* 'convolv_2invG_nov:78' checkerror=100; */
    checkerror = 100.0;

    /* 'convolv_2invG_nov:79' hh=.001; */
    hh = 0.001;

    /* 'convolv_2invG_nov:80' check1=.001*sum(onestagepdf2((0:.001:nu-r),m(1),s(1)))+.001*sum(onestagepdf2((nu+r:.001:Tu),m(1),s(1))); */
    gp = nu - r;
    if (rtIsNaN(gp)) {
      nm1d2 = c_y->size[0] * c_y->size[1];
      c_y->size[0] = 1;
      c_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
      c_y->data[0] = rtNaN;
    } else if (gp < 0.0) {
      nm1d2 = c_y->size[0] * c_y->size[1];
      c_y->size[0] = 1;
      c_y->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
    } else if (rtIsInf(gp) && (0.0 == gp)) {
      nm1d2 = c_y->size[0] * c_y->size[1];
      c_y->size[0] = 1;
      c_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
      c_y->data[0] = rtNaN;
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
        ix = (int)ndbl;
      } else {
        ix = 0;
      }

      nm1d2 = c_y->size[0] * c_y->size[1];
      c_y->size[0] = 1;
      c_y->size[1] = ix;
      emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
      if (ix > 0) {
        c_y->data[0] = 0.0;
        if (ix > 1) {
          c_y->data[ix - 1] = apnd;
          nm1d2 = (ix - 1) / 2;
          for (ixstart = 1; ixstart < nm1d2; ixstart++) {
            gp = (double)ixstart * 0.001;
            c_y->data[ixstart] = gp;
            c_y->data[(ix - ixstart) - 1] = apnd - gp;
          }

          if (nm1d2 << 1 == ix - 1) {
            c_y->data[nm1d2] = apnd / 2.0;
          } else {
            gp = (double)nm1d2 * 0.001;
            c_y->data[nm1d2] = gp;
            c_y->data[nm1d2 + 1] = apnd - gp;
          }
        }
      }
    }

    emxInit_real_T(&d_y, 2);
    a = nu + r;
    if (rtIsNaN(a) || rtIsNaN(Tu)) {
      nm1d2 = d_y->size[0] * d_y->size[1];
      d_y->size[0] = 1;
      d_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
      d_y->data[0] = rtNaN;
    } else if (Tu < a) {
      nm1d2 = d_y->size[0] * d_y->size[1];
      d_y->size[0] = 1;
      d_y->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
    } else if ((rtIsInf(a) || rtIsInf(Tu)) && (a == Tu)) {
      nm1d2 = d_y->size[0] * d_y->size[1];
      d_y->size[0] = 1;
      d_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
      d_y->data[0] = rtNaN;
    } else {
      ndbl = floor((Tu - a) / 0.001 + 0.5);
      apnd = a + ndbl * 0.001;
      cdiff = apnd - Tu;
      gp = fabs(a);
      if (!((gp > Tu) || rtIsNaN(Tu))) {
        gp = Tu;
      }

      if (fabs(cdiff) < 4.4408920985006262E-16 * gp) {
        ndbl++;
        apnd = Tu;
      } else if (cdiff > 0.0) {
        apnd = a + (ndbl - 1.0) * 0.001;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        ix = (int)ndbl;
      } else {
        ix = 0;
      }

      nm1d2 = d_y->size[0] * d_y->size[1];
      d_y->size[0] = 1;
      d_y->size[1] = ix;
      emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
      if (ix > 0) {
        d_y->data[0] = a;
        if (ix > 1) {
          d_y->data[ix - 1] = apnd;
          nm1d2 = (ix - 1) / 2;
          for (ixstart = 1; ixstart < nm1d2; ixstart++) {
            gp = (double)ixstart * 0.001;
            d_y->data[ixstart] = a + gp;
            d_y->data[(ix - ixstart) - 1] = apnd - gp;
          }

          if (nm1d2 << 1 == ix - 1) {
            d_y->data[nm1d2] = (a + apnd) / 2.0;
          } else {
            gp = (double)nm1d2 * 0.001;
            d_y->data[nm1d2] = a + gp;
            d_y->data[nm1d2 + 1] = apnd - gp;
          }
        }
      }
    }

    emxInit_real_T(&r10, 2);
    d_onestagepdf2(c_y, m[0], s[0], r10);
    d_onestagepdf2(d_y, m[0], s[0], c_y);
    check1 = 0.001 * b_sum(r10) + 0.001 * b_sum(c_y);

    /* 'convolv_2invG_nov:81' while checkerror>10^(-4) */
    while (checkerror > 0.0001) {
      /* 'convolv_2invG_nov:82' hh=.5*hh; */
      hh *= 0.5;

      /* 'convolv_2invG_nov:83' ck1=hh*sum(onestagepdf2((0:hh:nu-r),m(1),s(1)))+hh*sum(onestagepdf2((nu+r:hh:Tu),m(1),s(1))); */
      gp = nu - r;
      if (rtIsNaN(gp)) {
        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
        c_y->data[0] = rtNaN;
      } else if ((hh == 0.0) || ((0.0 < gp) && (hh < 0.0)) || ((gp < 0.0) && (hh
        > 0.0))) {
        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
      } else if (rtIsInf(gp) && (rtIsInf(hh) || (0.0 == gp))) {
        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
        c_y->data[0] = rtNaN;
      } else if (rtIsInf(hh)) {
        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
        c_y->data[0] = 0.0;
      } else if (floor(hh) == hh) {
        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = (int)floor(gp / hh) + 1;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
        ixstart = (int)floor(gp / hh);
        for (nm1d2 = 0; nm1d2 <= ixstart; nm1d2++) {
          c_y->data[c_y->size[0] * nm1d2] = hh * (double)nm1d2;
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
          ix = (int)ndbl;
        } else {
          ix = 0;
        }

        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = ix;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
        if (ix > 0) {
          c_y->data[0] = 0.0;
          if (ix > 1) {
            c_y->data[ix - 1] = apnd;
            nm1d2 = (ix - 1) / 2;
            for (ixstart = 1; ixstart < nm1d2; ixstart++) {
              gp = (double)ixstart * hh;
              c_y->data[ixstart] = gp;
              c_y->data[(ix - ixstart) - 1] = apnd - gp;
            }

            if (nm1d2 << 1 == ix - 1) {
              c_y->data[nm1d2] = apnd / 2.0;
            } else {
              gp = (double)nm1d2 * hh;
              c_y->data[nm1d2] = gp;
              c_y->data[nm1d2 + 1] = apnd - gp;
            }
          }
        }
      }

      a = nu + r;
      if (rtIsNaN(a) || rtIsNaN(Tu)) {
        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
        d_y->data[0] = rtNaN;
      } else if ((hh == 0.0) || ((a < Tu) && (hh < 0.0)) || ((Tu < a) && (hh >
                   0.0))) {
        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
      } else if ((rtIsInf(a) || rtIsInf(Tu)) && (rtIsInf(hh) || (a == Tu))) {
        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
        d_y->data[0] = rtNaN;
      } else if (rtIsInf(hh)) {
        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
        d_y->data[0] = a;
      } else if ((floor(a) == a) && (floor(hh) == hh)) {
        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = (int)floor((Tu - a) / hh) + 1;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
        ixstart = (int)floor((Tu - a) / hh);
        for (nm1d2 = 0; nm1d2 <= ixstart; nm1d2++) {
          d_y->data[d_y->size[0] * nm1d2] = a + hh * (double)nm1d2;
        }
      } else {
        ndbl = floor((Tu - a) / hh + 0.5);
        apnd = a + ndbl * hh;
        if (hh > 0.0) {
          cdiff = apnd - Tu;
        } else {
          cdiff = Tu - apnd;
        }

        gp = fabs(a);
        if (!((gp > Tu) || rtIsNaN(Tu))) {
          gp = Tu;
        }

        if (fabs(cdiff) < 4.4408920985006262E-16 * gp) {
          ndbl++;
          apnd = Tu;
        } else if (cdiff > 0.0) {
          apnd = a + (ndbl - 1.0) * hh;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          ix = (int)ndbl;
        } else {
          ix = 0;
        }

        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = ix;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
        if (ix > 0) {
          d_y->data[0] = a;
          if (ix > 1) {
            d_y->data[ix - 1] = apnd;
            nm1d2 = (ix - 1) / 2;
            for (ixstart = 1; ixstart < nm1d2; ixstart++) {
              gp = (double)ixstart * hh;
              d_y->data[ixstart] = a + gp;
              d_y->data[(ix - ixstart) - 1] = apnd - gp;
            }

            if (nm1d2 << 1 == ix - 1) {
              d_y->data[nm1d2] = (a + apnd) / 2.0;
            } else {
              gp = (double)nm1d2 * hh;
              d_y->data[nm1d2] = a + gp;
              d_y->data[nm1d2 + 1] = apnd - gp;
            }
          }
        }
      }

      d_onestagepdf2(c_y, m[0], s[0], r10);
      d_onestagepdf2(d_y, m[0], s[0], c_y);
      gp = hh * b_sum(r10) + hh * b_sum(c_y);

      /* 'convolv_2invG_nov:84' checkerror=abs(check1-ck1); */
      checkerror = fabs(check1 - gp);

      /* 'convolv_2invG_nov:85' check1=ck1; */
      check1 = gp;
    }

    emxFree_real_T(&r10);
    emxFree_real_T(&d_y);
    emxFree_real_T(&c_y);

    /* 'convolv_2invG_nov:87' check2=gm*check1; */
    /* 'convolv_2invG_nov:88' if  check2<=eps/3 */
    if (c_onestagepdf2(T2, m[1], s[1]) * check1 <= 7.4014868308343765E-17) {
      /* 'convolv_2invG_nov:90' flag=1; */
      *flag = 1.0;

      /* 'convolv_2invG_nov:91' sigma=s(2); */
      /* 'convolv_2invG_nov:92' mu=m(2); */
      /* 'convolv_2invG_nov:93' l=1/m(1); */
      /* 'convolv_2invG_nov:95' P=onestagepdf_lag(t,mu,sigma,l); */
      c_onestagepdf_lag(b_t, m[1], s[1], 1.0 / m[0], P);
    } else {
      /* 'convolv_2invG_nov:97' else */
      /*  find the discrete convolution of the vectors y and z */
      /*  the (i-1)th element of v approximates the convolution of the pdfs  */
      /*  over [.001, x(i)] as a left-hand Riemann sum. */
      /* 'convolv_2invG_nov:102' C=conv(z,y)*h; */
      conv(z, b_y, C);
      for (nm1d2 = 0; nm1d2 < 4401; nm1d2++) {
        C[nm1d2] *= 0.01;
      }

      /* 'convolv_2invG_nov:103' N=length(y); */
      /*  only the first N elements of the convolution are valid */
      /* 'convolv_2invG_nov:105' C=C(1:N); */
      /* 'convolv_2invG_nov:106' I=zeros(1,n); */
      /* 'convolv_2invG_nov:107' P=zeros(1,n); */
      /* 'convolv_2invG_nov:108' for i=1:n */
      /* toc */
      /* 'convolv_2invG_nov:121' P=P'; */
      for (i = 0; i < 2201; i++) {
        /* find element of x that is closest to t(i) */
        /* 'convolv_2invG_nov:110' [~,I(i)]=min((t(i)-x).^2); */
        for (nm1d2 = 0; nm1d2 < 2201; nm1d2++) {
          t[nm1d2] = b_t[i] - x[nm1d2];
        }

        d_power(t, varargin_1);
        ixstart = 1;
        gp = varargin_1[0];
        nm1d2 = 1;
        if (rtIsNaN(varargin_1[0])) {
          ix = 2;
          exitg1 = false;
          while ((!exitg1) && (ix < 2202)) {
            ixstart = ix;
            if (!rtIsNaN(varargin_1[ix - 1])) {
              gp = varargin_1[ix - 1];
              nm1d2 = ix;
              exitg1 = true;
            } else {
              ix++;
            }
          }
        }

        if (ixstart < 2201) {
          while (ixstart + 1 < 2202) {
            if (varargin_1[ixstart] < gp) {
              gp = varargin_1[ixstart];
              nm1d2 = ixstart + 1;
            }

            ixstart++;
          }
        }

        b_y[i] = nm1d2;

        /* 'convolv_2invG_nov:110' ~ */
        /* If t(i)<0 the probability is set to zero, otherwise the */
        /* probability is approxiated as a value from the vector x. */
        /* 'convolv_2invG_nov:113' if t(i)>0 && I(i)>1 */
        if ((b_t[i] > 0.0) && ((short)b_y[i] > 1)) {
          /* 'convolv_2invG_nov:114' P(i)=C(I(i)-1); */
          z[i] = C[(short)b_y[i] - 2];
        } else {
          /* 'convolv_2invG_nov:115' else */
          /* 'convolv_2invG_nov:116' P(i)=realmin; */
          z[i] = 2.2250738585072014E-308;
        }

        P[i] = z[i];
      }

      /* 'convolv_2invG_nov:123' logP=sum(log(P)); */
    }
  } else {
    /* 'convolv_2invG_nov:126' else */
    /*  find the discrete convolution of the vectors y and z */
    /*  the (i-1)th element of v approximates the convolution of the pdfs  */
    /*  over [.001, x(i)] as a left-hand Riemann sum. */
    /* 'convolv_2invG_nov:131' C=conv(z,y)*h; */
    conv(z, b_y, C);
    for (nm1d2 = 0; nm1d2 < 4401; nm1d2++) {
      C[nm1d2] *= 0.01;
    }

    /* 'convolv_2invG_nov:132' N=length(y); */
    /*  only the first N elements of the convolution are valid */
    /* 'convolv_2invG_nov:134' C=C(1:N); */
    /* 'convolv_2invG_nov:135' I=zeros(1,n); */
    /* 'convolv_2invG_nov:136' P=zeros(1,n); */
    /* 'convolv_2invG_nov:137' for i=1:n */
    /* toc */
    /* 'convolv_2invG_nov:150' P=P'; */
    for (i = 0; i < 2201; i++) {
      /* find element of x that is closest to t(i) */
      /* 'convolv_2invG_nov:139' [~,I(i)]=min((t(i)-x).^2); */
      for (nm1d2 = 0; nm1d2 < 2201; nm1d2++) {
        t[nm1d2] = b_t[i] - x[nm1d2];
      }

      d_power(t, varargin_1);
      ixstart = 1;
      gp = varargin_1[0];
      nm1d2 = 1;
      if (rtIsNaN(varargin_1[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix < 2202)) {
          ixstart = ix;
          if (!rtIsNaN(varargin_1[ix - 1])) {
            gp = varargin_1[ix - 1];
            nm1d2 = ix;
            exitg1 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < 2201) {
        while (ixstart + 1 < 2202) {
          if (varargin_1[ixstart] < gp) {
            gp = varargin_1[ixstart];
            nm1d2 = ixstart + 1;
          }

          ixstart++;
        }
      }

      b_y[i] = nm1d2;

      /* 'convolv_2invG_nov:139' ~ */
      /* If t(i)<0 the probability is set to zero, otherwise the */
      /* probability is approxiated as a value from the vector x. */
      /* 'convolv_2invG_nov:142' if t(i)>0 && I(i)>1 */
      if ((b_t[i] > 0.0) && ((short)b_y[i] > 1)) {
        /* 'convolv_2invG_nov:143' P(i)=C(I(i)-1); */
        z[i] = C[(short)b_y[i] - 2];
      } else {
        /* 'convolv_2invG_nov:144' else */
        /* 'convolv_2invG_nov:145' P(i)=realmin; */
        z[i] = 2.2250738585072014E-308;
      }

      P[i] = z[i];
    }

    /* 'convolv_2invG_nov:152' logP=sum(log(P)); */
  }
}

/*
 * function [P,flag]=convolv_2invG_nov(t,m1,s1,m2,s2,h)
 */
void convolv_2invG_nov(double m1, double s1, double m2, double s2, double P[221],
  double *flag)
{
  double y[2];
  int ixstart;
  double m[2];
  double s[2];
  double v[2];
  int nm1d2;
  int iidx[2];
  double sd[2];
  double b_s[2];
  double T2;
  static const double x[221] = { 0.0, 0.1, 0.2, 0.30000000000000004, 0.4, 0.5,
    0.60000000000000009, 0.70000000000000007, 0.8, 0.9, 1.0, 1.1,
    1.2000000000000002, 1.3, 1.4000000000000001, 1.5, 1.6, 1.7000000000000002,
    1.8, 1.9000000000000001, 2.0, 2.1, 2.2, 2.3000000000000003,
    2.4000000000000004, 2.5, 2.6, 2.7, 2.8000000000000003, 2.9000000000000004,
    3.0, 3.1, 3.2, 3.3000000000000003, 3.4000000000000004, 3.5, 3.6, 3.7,
    3.8000000000000003, 3.9000000000000004, 4.0, 4.1000000000000005, 4.2, 4.3,
    4.4, 4.5, 4.6000000000000005, 4.7, 4.8000000000000007, 4.9, 5.0,
    5.1000000000000005, 5.2, 5.3000000000000007, 5.4, 5.5, 5.6000000000000005,
    5.7, 5.8000000000000007, 5.9, 6.0, 6.1000000000000005, 6.2,
    6.3000000000000007, 6.4, 6.5, 6.6000000000000005, 6.7, 6.8000000000000007,
    6.9, 7.0, 7.1000000000000005, 7.2, 7.3000000000000007, 7.4, 7.5,
    7.6000000000000005, 7.7, 7.8000000000000007, 7.9, 8.0, 8.1,
    8.2000000000000011, 8.3, 8.4, 8.5, 8.6, 8.7000000000000011, 8.8, 8.9, 9.0,
    9.1, 9.2000000000000011, 9.3, 9.4, 9.5, 9.6000000000000014,
    9.7000000000000011, 9.8, 9.9, 10.0, 10.100000000000001, 10.200000000000001,
    10.3, 10.4, 10.5, 10.600000000000001, 10.700000000000001, 10.8, 10.9, 11.0,
    11.1, 11.2, 11.299999999999999, 11.399999999999999, 11.5, 11.6, 11.7,
    11.799999999999999, 11.899999999999999, 12.0, 12.1, 12.2, 12.299999999999999,
    12.399999999999999, 12.5, 12.6, 12.7, 12.799999999999999, 12.9, 13.0, 13.1,
    13.2, 13.299999999999999, 13.4, 13.5, 13.6, 13.7, 13.799999999999999, 13.9,
    14.0, 14.1, 14.2, 14.3, 14.399999999999999, 14.5, 14.6, 14.7, 14.8,
    14.899999999999999, 15.0, 15.1, 15.2, 15.3, 15.399999999999999, 15.5, 15.6,
    15.7, 15.8, 15.899999999999999, 16.0, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6,
    16.7, 16.8, 16.9, 17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9,
    18.0, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19.0, 19.1, 19.2,
    19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20.0, 20.1, 20.2, 20.3, 20.4, 20.5,
    20.6, 20.7, 20.8, 20.9, 21.0, 21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8,
    21.9, 22.0 };

  double b_y[221];
  double z[221];
  double C[441];
  double gp;
  double r;
  int i;
  double nu;
  double t[221];
  double varargin_1[221];
  static const double b_t[221] = { 0.0, 0.1, 0.2, 0.30000000000000004, 0.4, 0.5,
    0.60000000000000009, 0.70000000000000007, 0.8, 0.9, 1.0, 1.1,
    1.2000000000000002, 1.3, 1.4000000000000001, 1.5, 1.6, 1.7000000000000002,
    1.8, 1.9000000000000001, 2.0, 2.1, 2.2, 2.3000000000000003,
    2.4000000000000004, 2.5, 2.6, 2.7, 2.8000000000000003, 2.9000000000000004,
    3.0, 3.1, 3.2, 3.3000000000000003, 3.4000000000000004, 3.5, 3.6, 3.7,
    3.8000000000000003, 3.9000000000000004, 4.0, 4.1000000000000005, 4.2, 4.3,
    4.4, 4.5, 4.6000000000000005, 4.7, 4.8000000000000007, 4.9, 5.0,
    5.1000000000000005, 5.2, 5.3000000000000007, 5.4, 5.5, 5.6000000000000005,
    5.7, 5.8000000000000007, 5.9, 6.0, 6.1000000000000005, 6.2,
    6.3000000000000007, 6.4, 6.5, 6.6000000000000005, 6.7, 6.8000000000000007,
    6.9, 7.0, 7.1000000000000005, 7.2, 7.3000000000000007, 7.4, 7.5,
    7.6000000000000005, 7.7, 7.8000000000000007, 7.9, 8.0, 8.1,
    8.2000000000000011, 8.3, 8.4, 8.5, 8.6, 8.7000000000000011, 8.8, 8.9, 9.0,
    9.1, 9.2000000000000011, 9.3, 9.4, 9.5, 9.6000000000000014,
    9.7000000000000011, 9.8, 9.9, 10.0, 10.100000000000001, 10.200000000000001,
    10.3, 10.4, 10.5, 10.600000000000001, 10.700000000000001, 10.8, 10.9, 11.0,
    11.1, 11.2, 11.299999999999999, 11.399999999999999, 11.5, 11.6, 11.7,
    11.799999999999999, 11.899999999999999, 12.0, 12.1, 12.2, 12.299999999999999,
    12.399999999999999, 12.5, 12.6, 12.7, 12.799999999999999, 12.9, 13.0, 13.1,
    13.2, 13.299999999999999, 13.4, 13.5, 13.6, 13.7, 13.799999999999999, 13.9,
    14.0, 14.1, 14.2, 14.3, 14.399999999999999, 14.5, 14.6, 14.7, 14.8,
    14.899999999999999, 15.0, 15.1, 15.2, 15.3, 15.399999999999999, 15.5, 15.6,
    15.7, 15.8, 15.899999999999999, 16.0, 16.1, 16.2, 16.3, 16.4, 16.5, 16.6,
    16.7, 16.8, 16.9, 17.0, 17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9,
    18.0, 18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19.0, 19.1, 19.2,
    19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20.0, 20.1, 20.2, 20.3, 20.4, 20.5,
    20.6, 20.7, 20.8, 20.9, 21.0, 21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8,
    21.9, 22.0 };

  int ix;
  boolean_T exitg1;
  double Tu;
  emxArray_real_T *c_y;
  double checkerror;
  double hh;
  double ndbl;
  double apnd;
  double cdiff;
  emxArray_real_T *d_y;
  double a;
  emxArray_real_T *r2;
  double check1;

  /* This function evaluates the convolution of two inverse Gaussian */
  /* distributions at vector t. */
  /* If the variance in one of the distributions is very small so that the  */
  /* distribution is close the a Dirac delta function, the convolution */
  /* is approximated as a shifted inverse Gaussian, that is  */
  /* one part of the cell cycle is treated as deterministic in length.   */
  /* In this case, the shift or lag is the average time to complete  */
  /* the part of the cycle with the smallest sigma. */
  /* t is a vector of times to divide (or times to complete two parts of the */
  /* cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1, */
  /* s2=sigma2. */
  /* 'convolv_2invG_nov:16' flag=0; */
  *flag = 0.0;

  /* 'convolv_2invG_nov:18' E=Inf; */
  /* 'convolv_2invG_nov:19' n=length(t); */
  /* 'convolv_2invG_nov:20' m=[m1 m2]; */
  /* 'convolv_2invG_nov:21' s=[s1 s2]; */
  /* 'convolv_2invG_nov:22' m=max(0,m); */
  y[0] = m1;
  y[1] = m2;
  for (ixstart = 0; ixstart < 2; ixstart++) {
    if ((0.0 > y[ixstart]) || rtIsNaN(y[ixstart])) {
      m[ixstart] = 0.0;
    } else {
      m[ixstart] = y[ixstart];
    }
  }

  /* 'convolv_2invG_nov:23' s=max(0,s); */
  y[0] = s1;
  y[1] = s2;
  for (ixstart = 0; ixstart < 2; ixstart++) {
    if ((0.0 > y[ixstart]) || rtIsNaN(y[ixstart])) {
      s[ixstart] = 0.0;
    } else {
      s[ixstart] = y[ixstart];
    }
  }

  /* 'convolv_2invG_nov:24' v=(s.^2)./(m.^3); */
  /* 'convolv_2invG_nov:25' [v,I]=sort(v); */
  power(s, v);
  b_power(m, y);
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    v[nm1d2] /= y[nm1d2];
  }

  sort(v, iidx);

  /* 'convolv_2invG_nov:26' sd=v.^.5; */
  c_power(v, sd);

  /* 'convolv_2invG_nov:27' m=m(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    v[nm1d2] = m[iidx[nm1d2] - 1];
    y[nm1d2] = iidx[nm1d2];
  }

  /* 'convolv_2invG_nov:28' s=s(I); */
  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    m[nm1d2] = v[nm1d2];
    b_s[nm1d2] = s[(int)y[nm1d2] - 1];
  }

  for (nm1d2 = 0; nm1d2 < 2; nm1d2++) {
    s[nm1d2] = b_s[nm1d2];
  }

  /* 'convolv_2invG_nov:30' T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1))); */
  /* 'convolv_2invG_nov:31' T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2))); */
  T2 = 1.0 / m[1] * (sqrt(1.0 + 2.25 * (rt_powd_snf(s[1], 4.0) / (m[1] * m[1])))
                     - 1.5 * (s[1] * s[1] / m[1]));

  /* When called from convolv_3invG all of the components of t may be negative. */
  /* In this case, the probability of each is set to realmin */
  /* 'convolv_2invG_nov:36' if max(t)<=0 */
  /* 'convolv_2invG_nov:38' else */
  /* 'convolv_2invG_nov:41' Maxt=max(t); */
  /* 'convolv_2invG_nov:43' x=0:h:Maxt; */
  /* 'convolv_2invG_nov:45' y=onestagepdf2(x,m(1),s(1)); */
  f_onestagepdf2(x, m[0], s[0], b_y);

  /* 'convolv_2invG_nov:46' z=onestagepdf2(x,m(2),s(2)); */
  f_onestagepdf2(x, m[1], s[1], z);

  /* 'convolv_2invG_nov:48' if sd(1)<.01 */
  if (sd[0] < 0.01) {
    /* to estimate the error (in calculating probabilities the of the data),  */
    /* that results from approximating the first pdf as a point-mass */
    /* distribtion, find the maximum of the absolute value of the  */
    /* derivative of the second pdf. */
    /* 'convolv_2invG_nov:55' gp=gp_max(m(2),s(2)); */
    gp = gp_max(m[1], s[1]);

    /* determine the radius, r, of a small interval over which the second pdf, g, is */
    /* approximately constant and so that g(t) is small for t<r.   */
    /* 'convolv_2invG_nov:60' r=min(eps/(3*gp),T2); */
    gp = 2.2204460492503131E-16 / (3.0 * gp);
    if ((gp < T2) || rtIsNaN(T2)) {
      r = gp;
    } else {
      r = T2;
    }

    /* 'convolv_2invG_nov:61' checkval=onestagepdf2(r,m(2),s(2)); */
    gp = c_onestagepdf2(r, m[1], s[1]);

    /* 'convolv_2invG_nov:62' while checkval>=eps/2 */
    while (gp >= 1.1102230246251565E-16) {
      /* 'convolv_2invG_nov:63' r=r/2; */
      r /= 2.0;

      /* 'convolv_2invG_nov:64' checkval=onestagepdf2(r,m(2),s(2)); */
      gp = c_onestagepdf2(r, m[1], s[1]);
    }

    /* get the average value of the first pdf.  This is the point at which */
    /* its mass is concentrated. */
    /* 'convolv_2invG_nov:69' nu=1/m(1); */
    nu = 1.0 / m[0];

    /* get the maximum value f the second pdf. */
    /* 'convolv_2invG_nov:72' gm=onestagepdf2(T2,m(2),s(2)); */
    /* get upper limit of integral for approximating */
    /* int_{r+nu}^{infty}f(s)ds. */
    /* 'convolv_2invG_nov:75' Tu=max(100,nu+r+1000*sd(1)); */
    gp = (nu + r) + 1000.0 * sd[0];
    if ((100.0 > gp) || rtIsNaN(gp)) {
      Tu = 100.0;
    } else {
      Tu = gp;
    }

    emxInit_real_T(&c_y, 2);

    /* 'convolv_2invG_nov:78' checkerror=100; */
    checkerror = 100.0;

    /* 'convolv_2invG_nov:79' hh=.001; */
    hh = 0.001;

    /* 'convolv_2invG_nov:80' check1=.001*sum(onestagepdf2((0:.001:nu-r),m(1),s(1)))+.001*sum(onestagepdf2((nu+r:.001:Tu),m(1),s(1))); */
    gp = nu - r;
    if (rtIsNaN(gp)) {
      nm1d2 = c_y->size[0] * c_y->size[1];
      c_y->size[0] = 1;
      c_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
      c_y->data[0] = rtNaN;
    } else if (gp < 0.0) {
      nm1d2 = c_y->size[0] * c_y->size[1];
      c_y->size[0] = 1;
      c_y->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
    } else if (rtIsInf(gp) && (0.0 == gp)) {
      nm1d2 = c_y->size[0] * c_y->size[1];
      c_y->size[0] = 1;
      c_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
      c_y->data[0] = rtNaN;
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
        ix = (int)ndbl;
      } else {
        ix = 0;
      }

      nm1d2 = c_y->size[0] * c_y->size[1];
      c_y->size[0] = 1;
      c_y->size[1] = ix;
      emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
      if (ix > 0) {
        c_y->data[0] = 0.0;
        if (ix > 1) {
          c_y->data[ix - 1] = apnd;
          nm1d2 = (ix - 1) / 2;
          for (ixstart = 1; ixstart < nm1d2; ixstart++) {
            gp = (double)ixstart * 0.001;
            c_y->data[ixstart] = gp;
            c_y->data[(ix - ixstart) - 1] = apnd - gp;
          }

          if (nm1d2 << 1 == ix - 1) {
            c_y->data[nm1d2] = apnd / 2.0;
          } else {
            gp = (double)nm1d2 * 0.001;
            c_y->data[nm1d2] = gp;
            c_y->data[nm1d2 + 1] = apnd - gp;
          }
        }
      }
    }

    emxInit_real_T(&d_y, 2);
    a = nu + r;
    if (rtIsNaN(a) || rtIsNaN(Tu)) {
      nm1d2 = d_y->size[0] * d_y->size[1];
      d_y->size[0] = 1;
      d_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
      d_y->data[0] = rtNaN;
    } else if (Tu < a) {
      nm1d2 = d_y->size[0] * d_y->size[1];
      d_y->size[0] = 1;
      d_y->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
    } else if ((rtIsInf(a) || rtIsInf(Tu)) && (a == Tu)) {
      nm1d2 = d_y->size[0] * d_y->size[1];
      d_y->size[0] = 1;
      d_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
      d_y->data[0] = rtNaN;
    } else {
      ndbl = floor((Tu - a) / 0.001 + 0.5);
      apnd = a + ndbl * 0.001;
      cdiff = apnd - Tu;
      gp = fabs(a);
      if (!((gp > Tu) || rtIsNaN(Tu))) {
        gp = Tu;
      }

      if (fabs(cdiff) < 4.4408920985006262E-16 * gp) {
        ndbl++;
        apnd = Tu;
      } else if (cdiff > 0.0) {
        apnd = a + (ndbl - 1.0) * 0.001;
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        ix = (int)ndbl;
      } else {
        ix = 0;
      }

      nm1d2 = d_y->size[0] * d_y->size[1];
      d_y->size[0] = 1;
      d_y->size[1] = ix;
      emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
      if (ix > 0) {
        d_y->data[0] = a;
        if (ix > 1) {
          d_y->data[ix - 1] = apnd;
          nm1d2 = (ix - 1) / 2;
          for (ixstart = 1; ixstart < nm1d2; ixstart++) {
            gp = (double)ixstart * 0.001;
            d_y->data[ixstart] = a + gp;
            d_y->data[(ix - ixstart) - 1] = apnd - gp;
          }

          if (nm1d2 << 1 == ix - 1) {
            d_y->data[nm1d2] = (a + apnd) / 2.0;
          } else {
            gp = (double)nm1d2 * 0.001;
            d_y->data[nm1d2] = a + gp;
            d_y->data[nm1d2 + 1] = apnd - gp;
          }
        }
      }
    }

    emxInit_real_T(&r2, 2);
    d_onestagepdf2(c_y, m[0], s[0], r2);
    d_onestagepdf2(d_y, m[0], s[0], c_y);
    check1 = 0.001 * b_sum(r2) + 0.001 * b_sum(c_y);

    /* 'convolv_2invG_nov:81' while checkerror>10^(-4) */
    while (checkerror > 0.0001) {
      /* 'convolv_2invG_nov:82' hh=.5*hh; */
      hh *= 0.5;

      /* 'convolv_2invG_nov:83' ck1=hh*sum(onestagepdf2((0:hh:nu-r),m(1),s(1)))+hh*sum(onestagepdf2((nu+r:hh:Tu),m(1),s(1))); */
      gp = nu - r;
      if (rtIsNaN(gp)) {
        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
        c_y->data[0] = rtNaN;
      } else if ((hh == 0.0) || ((0.0 < gp) && (hh < 0.0)) || ((gp < 0.0) && (hh
        > 0.0))) {
        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
      } else if (rtIsInf(gp) && (rtIsInf(hh) || (0.0 == gp))) {
        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
        c_y->data[0] = rtNaN;
      } else if (rtIsInf(hh)) {
        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
        c_y->data[0] = 0.0;
      } else if (floor(hh) == hh) {
        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = (int)floor(gp / hh) + 1;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
        ixstart = (int)floor(gp / hh);
        for (nm1d2 = 0; nm1d2 <= ixstart; nm1d2++) {
          c_y->data[c_y->size[0] * nm1d2] = hh * (double)nm1d2;
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
          ix = (int)ndbl;
        } else {
          ix = 0;
        }

        nm1d2 = c_y->size[0] * c_y->size[1];
        c_y->size[0] = 1;
        c_y->size[1] = ix;
        emxEnsureCapacity((emxArray__common *)c_y, nm1d2, sizeof(double));
        if (ix > 0) {
          c_y->data[0] = 0.0;
          if (ix > 1) {
            c_y->data[ix - 1] = apnd;
            nm1d2 = (ix - 1) / 2;
            for (ixstart = 1; ixstart < nm1d2; ixstart++) {
              gp = (double)ixstart * hh;
              c_y->data[ixstart] = gp;
              c_y->data[(ix - ixstart) - 1] = apnd - gp;
            }

            if (nm1d2 << 1 == ix - 1) {
              c_y->data[nm1d2] = apnd / 2.0;
            } else {
              gp = (double)nm1d2 * hh;
              c_y->data[nm1d2] = gp;
              c_y->data[nm1d2 + 1] = apnd - gp;
            }
          }
        }
      }

      a = nu + r;
      if (rtIsNaN(a) || rtIsNaN(Tu)) {
        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
        d_y->data[0] = rtNaN;
      } else if ((hh == 0.0) || ((a < Tu) && (hh < 0.0)) || ((Tu < a) && (hh >
                   0.0))) {
        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
      } else if ((rtIsInf(a) || rtIsInf(Tu)) && (rtIsInf(hh) || (a == Tu))) {
        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
        d_y->data[0] = rtNaN;
      } else if (rtIsInf(hh)) {
        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
        d_y->data[0] = a;
      } else if ((floor(a) == a) && (floor(hh) == hh)) {
        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = (int)floor((Tu - a) / hh) + 1;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
        ixstart = (int)floor((Tu - a) / hh);
        for (nm1d2 = 0; nm1d2 <= ixstart; nm1d2++) {
          d_y->data[d_y->size[0] * nm1d2] = a + hh * (double)nm1d2;
        }
      } else {
        ndbl = floor((Tu - a) / hh + 0.5);
        apnd = a + ndbl * hh;
        if (hh > 0.0) {
          cdiff = apnd - Tu;
        } else {
          cdiff = Tu - apnd;
        }

        gp = fabs(a);
        if (!((gp > Tu) || rtIsNaN(Tu))) {
          gp = Tu;
        }

        if (fabs(cdiff) < 4.4408920985006262E-16 * gp) {
          ndbl++;
          apnd = Tu;
        } else if (cdiff > 0.0) {
          apnd = a + (ndbl - 1.0) * hh;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          ix = (int)ndbl;
        } else {
          ix = 0;
        }

        nm1d2 = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = ix;
        emxEnsureCapacity((emxArray__common *)d_y, nm1d2, sizeof(double));
        if (ix > 0) {
          d_y->data[0] = a;
          if (ix > 1) {
            d_y->data[ix - 1] = apnd;
            nm1d2 = (ix - 1) / 2;
            for (ixstart = 1; ixstart < nm1d2; ixstart++) {
              gp = (double)ixstart * hh;
              d_y->data[ixstart] = a + gp;
              d_y->data[(ix - ixstart) - 1] = apnd - gp;
            }

            if (nm1d2 << 1 == ix - 1) {
              d_y->data[nm1d2] = (a + apnd) / 2.0;
            } else {
              gp = (double)nm1d2 * hh;
              d_y->data[nm1d2] = a + gp;
              d_y->data[nm1d2 + 1] = apnd - gp;
            }
          }
        }
      }

      d_onestagepdf2(c_y, m[0], s[0], r2);
      d_onestagepdf2(d_y, m[0], s[0], c_y);
      gp = hh * b_sum(r2) + hh * b_sum(c_y);

      /* 'convolv_2invG_nov:84' checkerror=abs(check1-ck1); */
      checkerror = fabs(check1 - gp);

      /* 'convolv_2invG_nov:85' check1=ck1; */
      check1 = gp;
    }

    emxFree_real_T(&r2);
    emxFree_real_T(&d_y);
    emxFree_real_T(&c_y);

    /* 'convolv_2invG_nov:87' check2=gm*check1; */
    /* 'convolv_2invG_nov:88' if  check2<=eps/3 */
    if (c_onestagepdf2(T2, m[1], s[1]) * check1 <= 7.4014868308343765E-17) {
      /* 'convolv_2invG_nov:90' flag=1; */
      *flag = 1.0;

      /* 'convolv_2invG_nov:91' sigma=s(2); */
      /* 'convolv_2invG_nov:92' mu=m(2); */
      /* 'convolv_2invG_nov:93' l=1/m(1); */
      /* 'convolv_2invG_nov:95' P=onestagepdf_lag(t,mu,sigma,l); */
      b_onestagepdf_lag(b_t, m[1], s[1], 1.0 / m[0], P);
    } else {
      /* 'convolv_2invG_nov:97' else */
      /*  find the discrete convolution of the vectors y and z */
      /*  the (i-1)th element of v approximates the convolution of the pdfs  */
      /*  over [.001, x(i)] as a left-hand Riemann sum. */
      /* 'convolv_2invG_nov:102' C=conv(z,y)*h; */
      d_conv(z, b_y, C);
      for (nm1d2 = 0; nm1d2 < 441; nm1d2++) {
        C[nm1d2] *= 0.1;
      }

      /* 'convolv_2invG_nov:103' N=length(y); */
      /*  only the first N elements of the convolution are valid */
      /* 'convolv_2invG_nov:105' C=C(1:N); */
      /* 'convolv_2invG_nov:106' I=zeros(1,n); */
      /* 'convolv_2invG_nov:107' P=zeros(1,n); */
      /* 'convolv_2invG_nov:108' for i=1:n */
      /* toc */
      /* 'convolv_2invG_nov:121' P=P'; */
      for (i = 0; i < 221; i++) {
        /* find element of x that is closest to t(i) */
        /* 'convolv_2invG_nov:110' [~,I(i)]=min((t(i)-x).^2); */
        for (nm1d2 = 0; nm1d2 < 221; nm1d2++) {
          t[nm1d2] = b_t[i] - x[nm1d2];
        }

        j_power(t, varargin_1);
        ixstart = 1;
        gp = varargin_1[0];
        nm1d2 = 1;
        if (rtIsNaN(varargin_1[0])) {
          ix = 2;
          exitg1 = false;
          while ((!exitg1) && (ix < 222)) {
            ixstart = ix;
            if (!rtIsNaN(varargin_1[ix - 1])) {
              gp = varargin_1[ix - 1];
              nm1d2 = ix;
              exitg1 = true;
            } else {
              ix++;
            }
          }
        }

        if (ixstart < 221) {
          while (ixstart + 1 < 222) {
            if (varargin_1[ixstart] < gp) {
              gp = varargin_1[ixstart];
              nm1d2 = ixstart + 1;
            }

            ixstart++;
          }
        }

        b_y[i] = nm1d2;

        /* 'convolv_2invG_nov:110' ~ */
        /* If t(i)<0 the probability is set to zero, otherwise the */
        /* probability is approxiated as a value from the vector x. */
        /* 'convolv_2invG_nov:113' if t(i)>0 && I(i)>1 */
        if ((b_t[i] > 0.0) && ((unsigned char)b_y[i] > 1)) {
          /* 'convolv_2invG_nov:114' P(i)=C(I(i)-1); */
          z[i] = C[(unsigned char)b_y[i] - 2];
        } else {
          /* 'convolv_2invG_nov:115' else */
          /* 'convolv_2invG_nov:116' P(i)=realmin; */
          z[i] = 2.2250738585072014E-308;
        }

        P[i] = z[i];
      }

      /* 'convolv_2invG_nov:123' logP=sum(log(P)); */
    }
  } else {
    /* 'convolv_2invG_nov:126' else */
    /*  find the discrete convolution of the vectors y and z */
    /*  the (i-1)th element of v approximates the convolution of the pdfs  */
    /*  over [.001, x(i)] as a left-hand Riemann sum. */
    /* 'convolv_2invG_nov:131' C=conv(z,y)*h; */
    d_conv(z, b_y, C);
    for (nm1d2 = 0; nm1d2 < 441; nm1d2++) {
      C[nm1d2] *= 0.1;
    }

    /* 'convolv_2invG_nov:132' N=length(y); */
    /*  only the first N elements of the convolution are valid */
    /* 'convolv_2invG_nov:134' C=C(1:N); */
    /* 'convolv_2invG_nov:135' I=zeros(1,n); */
    /* 'convolv_2invG_nov:136' P=zeros(1,n); */
    /* 'convolv_2invG_nov:137' for i=1:n */
    /* toc */
    /* 'convolv_2invG_nov:150' P=P'; */
    for (i = 0; i < 221; i++) {
      /* find element of x that is closest to t(i) */
      /* 'convolv_2invG_nov:139' [~,I(i)]=min((t(i)-x).^2); */
      for (nm1d2 = 0; nm1d2 < 221; nm1d2++) {
        t[nm1d2] = b_t[i] - x[nm1d2];
      }

      j_power(t, varargin_1);
      ixstart = 1;
      gp = varargin_1[0];
      nm1d2 = 1;
      if (rtIsNaN(varargin_1[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix < 222)) {
          ixstart = ix;
          if (!rtIsNaN(varargin_1[ix - 1])) {
            gp = varargin_1[ix - 1];
            nm1d2 = ix;
            exitg1 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < 221) {
        while (ixstart + 1 < 222) {
          if (varargin_1[ixstart] < gp) {
            gp = varargin_1[ixstart];
            nm1d2 = ixstart + 1;
          }

          ixstart++;
        }
      }

      b_y[i] = nm1d2;

      /* 'convolv_2invG_nov:139' ~ */
      /* If t(i)<0 the probability is set to zero, otherwise the */
      /* probability is approxiated as a value from the vector x. */
      /* 'convolv_2invG_nov:142' if t(i)>0 && I(i)>1 */
      if ((b_t[i] > 0.0) && ((unsigned char)b_y[i] > 1)) {
        /* 'convolv_2invG_nov:143' P(i)=C(I(i)-1); */
        z[i] = C[(unsigned char)b_y[i] - 2];
      } else {
        /* 'convolv_2invG_nov:144' else */
        /* 'convolv_2invG_nov:145' P(i)=realmin; */
        z[i] = 2.2250738585072014E-308;
      }

      P[i] = z[i];
    }

    /* 'convolv_2invG_nov:152' logP=sum(log(P)); */
  }
}

/*
 * function [P,flag]=convolv_2invG_nov(t,m1,s1,m2,s2,h)
 */
void d_convolv_2invG_nov(const double t[22001], double m1, double s1, double m2,
  double s2, double P[22001], double *flag)
{
  double y[2];
  int ixstart;
  double m[2];
  double s[2];
  double v[2];
  int itmp;
  int iidx[2];
  double sd[2];
  double b_s[2];
  double T2;
  double kd;
  int nm1d2;
  boolean_T exitg1;
  int i;
  emxArray_real_T *b_y;
  double ndbl;
  double apnd;
  double cdiff;
  emxArray_real_T *c_y;
  emxArray_real_T *r12;
  int n;
  emxArray_real_T *z;
  emxArray_real_T *C;
  emxArray_real_T *a;
  double r;
  emxArray_real_T *b_a;
  double nu;
  double Tu;
  double checkerror;
  double hh;
  emxArray_real_T *b_t;
  emxArray_real_T *d_y;
  double c_a;
  double check1;
  unsigned int I[22001];

  /* This function evaluates the convolution of two inverse Gaussian */
  /* distributions at vector t. */
  /* If the variance in one of the distributions is very small so that the  */
  /* distribution is close the a Dirac delta function, the convolution */
  /* is approximated as a shifted inverse Gaussian, that is  */
  /* one part of the cell cycle is treated as deterministic in length.   */
  /* In this case, the shift or lag is the average time to complete  */
  /* the part of the cycle with the smallest sigma. */
  /* t is a vector of times to divide (or times to complete two parts of the */
  /* cell cycle when called by convolv_3invG), m1=mu1, m2=mu2, s1=sigma1, */
  /* s2=sigma2. */
  /* 'convolv_2invG_nov:16' flag=0; */
  *flag = 0.0;

  /* 'convolv_2invG_nov:18' E=Inf; */
  /* 'convolv_2invG_nov:19' n=length(t); */
  /* 'convolv_2invG_nov:20' m=[m1 m2]; */
  /* 'convolv_2invG_nov:21' s=[s1 s2]; */
  /* 'convolv_2invG_nov:22' m=max(0,m); */
  y[0] = m1;
  y[1] = m2;
  for (ixstart = 0; ixstart < 2; ixstart++) {
    if ((0.0 > y[ixstart]) || rtIsNaN(y[ixstart])) {
      m[ixstart] = 0.0;
    } else {
      m[ixstart] = y[ixstart];
    }
  }

  /* 'convolv_2invG_nov:23' s=max(0,s); */
  y[0] = s1;
  y[1] = s2;
  for (ixstart = 0; ixstart < 2; ixstart++) {
    if ((0.0 > y[ixstart]) || rtIsNaN(y[ixstart])) {
      s[ixstart] = 0.0;
    } else {
      s[ixstart] = y[ixstart];
    }
  }

  /* 'convolv_2invG_nov:24' v=(s.^2)./(m.^3); */
  /* 'convolv_2invG_nov:25' [v,I]=sort(v); */
  power(s, v);
  b_power(m, y);
  for (itmp = 0; itmp < 2; itmp++) {
    v[itmp] /= y[itmp];
  }

  sort(v, iidx);

  /* 'convolv_2invG_nov:26' sd=v.^.5; */
  c_power(v, sd);

  /* 'convolv_2invG_nov:27' m=m(I); */
  for (itmp = 0; itmp < 2; itmp++) {
    v[itmp] = m[iidx[itmp] - 1];
    y[itmp] = iidx[itmp];
  }

  /* 'convolv_2invG_nov:28' s=s(I); */
  for (itmp = 0; itmp < 2; itmp++) {
    m[itmp] = v[itmp];
    b_s[itmp] = s[(int)y[itmp] - 1];
  }

  for (itmp = 0; itmp < 2; itmp++) {
    s[itmp] = b_s[itmp];
  }

  /* 'convolv_2invG_nov:30' T1=(1/m(1))*((1+(9/4)*(s(1)^4/m(1)^2))^.5-(3/2)*(s(1)^2/m(1))); */
  /* 'convolv_2invG_nov:31' T2=(1/m(2))*((1+(9/4)*(s(2)^4/m(2)^2))^.5-(3/2)*(s(2)^2/m(2))); */
  T2 = 1.0 / m[1] * (sqrt(1.0 + 2.25 * (rt_powd_snf(s[1], 4.0) / (m[1] * m[1])))
                     - 1.5 * (s[1] * s[1] / m[1]));

  /* When called from convolv_3invG all of the components of t may be negative. */
  /* In this case, the probability of each is set to realmin */
  /* 'convolv_2invG_nov:36' if max(t)<=0 */
  ixstart = 1;
  kd = t[0];
  if (rtIsNaN(t[0])) {
    nm1d2 = 2;
    exitg1 = false;
    while ((!exitg1) && (nm1d2 < 22002)) {
      ixstart = nm1d2;
      if (!rtIsNaN(t[nm1d2 - 1])) {
        kd = t[nm1d2 - 1];
        exitg1 = true;
      } else {
        nm1d2++;
      }
    }
  }

  if (ixstart < 22001) {
    while (ixstart + 1 < 22002) {
      if (t[ixstart] > kd) {
        kd = t[ixstart];
      }

      ixstart++;
    }
  }

  if (kd <= 0.0) {
    /* 'convolv_2invG_nov:37' P=realmin.*ones(size(t)); */
    for (i = 0; i < 22001; i++) {
      P[i] = 2.2250738585072014E-308;
    }
  } else {
    /* 'convolv_2invG_nov:38' else */
    /* 'convolv_2invG_nov:41' Maxt=max(t); */
    ixstart = 1;
    kd = t[0];
    if (rtIsNaN(t[0])) {
      nm1d2 = 2;
      exitg1 = false;
      while ((!exitg1) && (nm1d2 < 22002)) {
        ixstart = nm1d2;
        if (!rtIsNaN(t[nm1d2 - 1])) {
          kd = t[nm1d2 - 1];
          exitg1 = true;
        } else {
          nm1d2++;
        }
      }
    }

    if (ixstart < 22001) {
      while (ixstart + 1 < 22002) {
        if (t[ixstart] > kd) {
          kd = t[ixstart];
        }

        ixstart++;
      }
    }

    /* 'convolv_2invG_nov:43' x=0:h:Maxt; */
    emxInit_real_T(&b_y, 2);
    if (rtIsNaN(kd)) {
      itmp = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)b_y, itmp, sizeof(double));
      b_y->data[0] = rtNaN;
    } else if (kd < 0.0) {
      itmp = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = 0;
      emxEnsureCapacity((emxArray__common *)b_y, itmp, sizeof(double));
    } else if (rtIsInf(kd) && (0.0 == kd)) {
      itmp = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = 1;
      emxEnsureCapacity((emxArray__common *)b_y, itmp, sizeof(double));
      b_y->data[0] = rtNaN;
    } else {
      ndbl = floor(kd / 0.001 + 0.5);
      apnd = ndbl * 0.001;
      cdiff = apnd - kd;
      if (fabs(cdiff) < 4.4408920985006262E-16 * kd) {
        ndbl++;
        apnd = kd;
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

      itmp = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = n;
      emxEnsureCapacity((emxArray__common *)b_y, itmp, sizeof(double));
      if (n > 0) {
        b_y->data[0] = 0.0;
        if (n > 1) {
          b_y->data[n - 1] = apnd;
          nm1d2 = (n - 1) / 2;
          for (ixstart = 1; ixstart < nm1d2; ixstart++) {
            kd = (double)ixstart * 0.001;
            b_y->data[ixstart] = kd;
            b_y->data[(n - ixstart) - 1] = apnd - kd;
          }

          if (nm1d2 << 1 == n - 1) {
            b_y->data[nm1d2] = apnd / 2.0;
          } else {
            kd = (double)nm1d2 * 0.001;
            b_y->data[nm1d2] = kd;
            b_y->data[nm1d2 + 1] = apnd - kd;
          }
        }
      }
    }

    emxInit_real_T(&c_y, 2);
    emxInit_real_T(&r12, 2);

    /* 'convolv_2invG_nov:45' y=onestagepdf2(x,m(1),s(1)); */
    d_onestagepdf2(b_y, m[0], s[0], r12);
    itmp = c_y->size[0] * c_y->size[1];
    c_y->size[0] = 1;
    c_y->size[1] = r12->size[1];
    emxEnsureCapacity((emxArray__common *)c_y, itmp, sizeof(double));
    nm1d2 = r12->size[0] * r12->size[1];
    for (itmp = 0; itmp < nm1d2; itmp++) {
      c_y->data[itmp] = r12->data[itmp];
    }

    emxInit_real_T(&z, 2);

    /* 'convolv_2invG_nov:46' z=onestagepdf2(x,m(2),s(2)); */
    d_onestagepdf2(b_y, m[1], s[1], r12);
    itmp = z->size[0] * z->size[1];
    z->size[0] = 1;
    z->size[1] = r12->size[1];
    emxEnsureCapacity((emxArray__common *)z, itmp, sizeof(double));
    nm1d2 = r12->size[0] * r12->size[1];
    for (itmp = 0; itmp < nm1d2; itmp++) {
      z->data[itmp] = r12->data[itmp];
    }

    /* 'convolv_2invG_nov:48' if sd(1)<.01 */
    emxInit_real_T(&C, 2);
    emxInit_real_T(&a, 2);
    if (sd[0] < 0.01) {
      /* to estimate the error (in calculating probabilities the of the data),  */
      /* that results from approximating the first pdf as a point-mass */
      /* distribtion, find the maximum of the absolute value of the  */
      /* derivative of the second pdf. */
      /* 'convolv_2invG_nov:55' gp=gp_max(m(2),s(2)); */
      kd = gp_max(m[1], s[1]);

      /* determine the radius, r, of a small interval over which the second pdf, g, is */
      /* approximately constant and so that g(t) is small for t<r.   */
      /* 'convolv_2invG_nov:60' r=min(eps/(3*gp),T2); */
      kd = 2.2204460492503131E-16 / (3.0 * kd);
      if ((kd < T2) || rtIsNaN(T2)) {
        r = kd;
      } else {
        r = T2;
      }

      /* 'convolv_2invG_nov:61' checkval=onestagepdf2(r,m(2),s(2)); */
      kd = c_onestagepdf2(r, m[1], s[1]);

      /* 'convolv_2invG_nov:62' while checkval>=eps/2 */
      while (kd >= 1.1102230246251565E-16) {
        /* 'convolv_2invG_nov:63' r=r/2; */
        r /= 2.0;

        /* 'convolv_2invG_nov:64' checkval=onestagepdf2(r,m(2),s(2)); */
        kd = c_onestagepdf2(r, m[1], s[1]);
      }

      /* get the average value of the first pdf.  This is the point at which */
      /* its mass is concentrated. */
      /* 'convolv_2invG_nov:69' nu=1/m(1); */
      nu = 1.0 / m[0];

      /* get the maximum value f the second pdf. */
      /* 'convolv_2invG_nov:72' gm=onestagepdf2(T2,m(2),s(2)); */
      /* get upper limit of integral for approximating */
      /* int_{r+nu}^{infty}f(s)ds. */
      /* 'convolv_2invG_nov:75' Tu=max(100,nu+r+1000*sd(1)); */
      kd = (nu + r) + 1000.0 * sd[0];
      if ((100.0 > kd) || rtIsNaN(kd)) {
        Tu = 100.0;
      } else {
        Tu = kd;
      }

      /* 'convolv_2invG_nov:78' checkerror=100; */
      checkerror = 100.0;

      /* 'convolv_2invG_nov:79' hh=.001; */
      hh = 0.001;

      /* 'convolv_2invG_nov:80' check1=.001*sum(onestagepdf2((0:.001:nu-r),m(1),s(1)))+.001*sum(onestagepdf2((nu+r:.001:Tu),m(1),s(1))); */
      kd = nu - r;
      if (rtIsNaN(kd)) {
        itmp = a->size[0] * a->size[1];
        a->size[0] = 1;
        a->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
        a->data[0] = rtNaN;
      } else if (kd < 0.0) {
        itmp = a->size[0] * a->size[1];
        a->size[0] = 1;
        a->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
      } else if (rtIsInf(kd) && (0.0 == kd)) {
        itmp = a->size[0] * a->size[1];
        a->size[0] = 1;
        a->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
        a->data[0] = rtNaN;
      } else {
        ndbl = floor(kd / 0.001 + 0.5);
        apnd = ndbl * 0.001;
        cdiff = apnd - kd;
        if (fabs(cdiff) < 4.4408920985006262E-16 * kd) {
          ndbl++;
          apnd = kd;
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

        itmp = a->size[0] * a->size[1];
        a->size[0] = 1;
        a->size[1] = n;
        emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
        if (n > 0) {
          a->data[0] = 0.0;
          if (n > 1) {
            a->data[n - 1] = apnd;
            nm1d2 = (n - 1) / 2;
            for (ixstart = 1; ixstart < nm1d2; ixstart++) {
              kd = (double)ixstart * 0.001;
              a->data[ixstart] = kd;
              a->data[(n - ixstart) - 1] = apnd - kd;
            }

            if (nm1d2 << 1 == n - 1) {
              a->data[nm1d2] = apnd / 2.0;
            } else {
              kd = (double)nm1d2 * 0.001;
              a->data[nm1d2] = kd;
              a->data[nm1d2 + 1] = apnd - kd;
            }
          }
        }
      }

      emxInit_real_T(&d_y, 2);
      c_a = nu + r;
      if (rtIsNaN(c_a) || rtIsNaN(Tu)) {
        itmp = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)d_y, itmp, sizeof(double));
        d_y->data[0] = rtNaN;
      } else if (Tu < c_a) {
        itmp = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = 0;
        emxEnsureCapacity((emxArray__common *)d_y, itmp, sizeof(double));
      } else if ((rtIsInf(c_a) || rtIsInf(Tu)) && (c_a == Tu)) {
        itmp = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)d_y, itmp, sizeof(double));
        d_y->data[0] = rtNaN;
      } else {
        ndbl = floor((Tu - c_a) / 0.001 + 0.5);
        apnd = c_a + ndbl * 0.001;
        cdiff = apnd - Tu;
        kd = fabs(c_a);
        if (!((kd > Tu) || rtIsNaN(Tu))) {
          kd = Tu;
        }

        if (fabs(cdiff) < 4.4408920985006262E-16 * kd) {
          ndbl++;
          apnd = Tu;
        } else if (cdiff > 0.0) {
          apnd = c_a + (ndbl - 1.0) * 0.001;
        } else {
          ndbl++;
        }

        if (ndbl >= 0.0) {
          n = (int)ndbl;
        } else {
          n = 0;
        }

        itmp = d_y->size[0] * d_y->size[1];
        d_y->size[0] = 1;
        d_y->size[1] = n;
        emxEnsureCapacity((emxArray__common *)d_y, itmp, sizeof(double));
        if (n > 0) {
          d_y->data[0] = c_a;
          if (n > 1) {
            d_y->data[n - 1] = apnd;
            nm1d2 = (n - 1) / 2;
            for (ixstart = 1; ixstart < nm1d2; ixstart++) {
              kd = (double)ixstart * 0.001;
              d_y->data[ixstart] = c_a + kd;
              d_y->data[(n - ixstart) - 1] = apnd - kd;
            }

            if (nm1d2 << 1 == n - 1) {
              d_y->data[nm1d2] = (c_a + apnd) / 2.0;
            } else {
              kd = (double)nm1d2 * 0.001;
              d_y->data[nm1d2] = c_a + kd;
              d_y->data[nm1d2 + 1] = apnd - kd;
            }
          }
        }
      }

      d_onestagepdf2(a, m[0], s[0], r12);
      d_onestagepdf2(d_y, m[0], s[0], a);
      check1 = 0.001 * b_sum(r12) + 0.001 * b_sum(a);

      /* 'convolv_2invG_nov:81' while checkerror>10^(-4) */
      while (checkerror > 0.0001) {
        /* 'convolv_2invG_nov:82' hh=.5*hh; */
        hh *= 0.5;

        /* 'convolv_2invG_nov:83' ck1=hh*sum(onestagepdf2((0:hh:nu-r),m(1),s(1)))+hh*sum(onestagepdf2((nu+r:hh:Tu),m(1),s(1))); */
        kd = nu - r;
        if (rtIsNaN(kd)) {
          itmp = a->size[0] * a->size[1];
          a->size[0] = 1;
          a->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
          a->data[0] = rtNaN;
        } else if ((hh == 0.0) || ((0.0 < kd) && (hh < 0.0)) || ((kd < 0.0) &&
                    (hh > 0.0))) {
          itmp = a->size[0] * a->size[1];
          a->size[0] = 1;
          a->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
        } else if (rtIsInf(kd) && (rtIsInf(hh) || (0.0 == kd))) {
          itmp = a->size[0] * a->size[1];
          a->size[0] = 1;
          a->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
          a->data[0] = rtNaN;
        } else if (rtIsInf(hh)) {
          itmp = a->size[0] * a->size[1];
          a->size[0] = 1;
          a->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
          a->data[0] = 0.0;
        } else if (floor(hh) == hh) {
          itmp = a->size[0] * a->size[1];
          a->size[0] = 1;
          a->size[1] = (int)floor(kd / hh) + 1;
          emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
          nm1d2 = (int)floor(kd / hh);
          for (itmp = 0; itmp <= nm1d2; itmp++) {
            a->data[a->size[0] * itmp] = hh * (double)itmp;
          }
        } else {
          ndbl = floor(kd / hh + 0.5);
          apnd = ndbl * hh;
          if (hh > 0.0) {
            cdiff = apnd - kd;
          } else {
            cdiff = kd - apnd;
          }

          if (fabs(cdiff) < 4.4408920985006262E-16 * fabs(kd)) {
            ndbl++;
            apnd = kd;
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

          itmp = a->size[0] * a->size[1];
          a->size[0] = 1;
          a->size[1] = n;
          emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
          if (n > 0) {
            a->data[0] = 0.0;
            if (n > 1) {
              a->data[n - 1] = apnd;
              nm1d2 = (n - 1) / 2;
              for (ixstart = 1; ixstart < nm1d2; ixstart++) {
                kd = (double)ixstart * hh;
                a->data[ixstart] = kd;
                a->data[(n - ixstart) - 1] = apnd - kd;
              }

              if (nm1d2 << 1 == n - 1) {
                a->data[nm1d2] = apnd / 2.0;
              } else {
                kd = (double)nm1d2 * hh;
                a->data[nm1d2] = kd;
                a->data[nm1d2 + 1] = apnd - kd;
              }
            }
          }
        }

        c_a = nu + r;
        if (rtIsNaN(c_a) || rtIsNaN(Tu)) {
          itmp = d_y->size[0] * d_y->size[1];
          d_y->size[0] = 1;
          d_y->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)d_y, itmp, sizeof(double));
          d_y->data[0] = rtNaN;
        } else if ((hh == 0.0) || ((c_a < Tu) && (hh < 0.0)) || ((Tu < c_a) &&
                    (hh > 0.0))) {
          itmp = d_y->size[0] * d_y->size[1];
          d_y->size[0] = 1;
          d_y->size[1] = 0;
          emxEnsureCapacity((emxArray__common *)d_y, itmp, sizeof(double));
        } else if ((rtIsInf(c_a) || rtIsInf(Tu)) && (rtIsInf(hh) || (c_a == Tu)))
        {
          itmp = d_y->size[0] * d_y->size[1];
          d_y->size[0] = 1;
          d_y->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)d_y, itmp, sizeof(double));
          d_y->data[0] = rtNaN;
        } else if (rtIsInf(hh)) {
          itmp = d_y->size[0] * d_y->size[1];
          d_y->size[0] = 1;
          d_y->size[1] = 1;
          emxEnsureCapacity((emxArray__common *)d_y, itmp, sizeof(double));
          d_y->data[0] = c_a;
        } else if ((floor(c_a) == c_a) && (floor(hh) == hh)) {
          itmp = d_y->size[0] * d_y->size[1];
          d_y->size[0] = 1;
          d_y->size[1] = (int)floor((Tu - c_a) / hh) + 1;
          emxEnsureCapacity((emxArray__common *)d_y, itmp, sizeof(double));
          nm1d2 = (int)floor((Tu - c_a) / hh);
          for (itmp = 0; itmp <= nm1d2; itmp++) {
            d_y->data[d_y->size[0] * itmp] = c_a + hh * (double)itmp;
          }
        } else {
          ndbl = floor((Tu - c_a) / hh + 0.5);
          apnd = c_a + ndbl * hh;
          if (hh > 0.0) {
            cdiff = apnd - Tu;
          } else {
            cdiff = Tu - apnd;
          }

          kd = fabs(c_a);
          if (!((kd > Tu) || rtIsNaN(Tu))) {
            kd = Tu;
          }

          if (fabs(cdiff) < 4.4408920985006262E-16 * kd) {
            ndbl++;
            apnd = Tu;
          } else if (cdiff > 0.0) {
            apnd = c_a + (ndbl - 1.0) * hh;
          } else {
            ndbl++;
          }

          if (ndbl >= 0.0) {
            n = (int)ndbl;
          } else {
            n = 0;
          }

          itmp = d_y->size[0] * d_y->size[1];
          d_y->size[0] = 1;
          d_y->size[1] = n;
          emxEnsureCapacity((emxArray__common *)d_y, itmp, sizeof(double));
          if (n > 0) {
            d_y->data[0] = c_a;
            if (n > 1) {
              d_y->data[n - 1] = apnd;
              nm1d2 = (n - 1) / 2;
              for (ixstart = 1; ixstart < nm1d2; ixstart++) {
                kd = (double)ixstart * hh;
                d_y->data[ixstart] = c_a + kd;
                d_y->data[(n - ixstart) - 1] = apnd - kd;
              }

              if (nm1d2 << 1 == n - 1) {
                d_y->data[nm1d2] = (c_a + apnd) / 2.0;
              } else {
                kd = (double)nm1d2 * hh;
                d_y->data[nm1d2] = c_a + kd;
                d_y->data[nm1d2 + 1] = apnd - kd;
              }
            }
          }
        }

        d_onestagepdf2(a, m[0], s[0], r12);
        d_onestagepdf2(d_y, m[0], s[0], a);
        kd = hh * b_sum(r12) + hh * b_sum(a);

        /* 'convolv_2invG_nov:84' checkerror=abs(check1-ck1); */
        checkerror = fabs(check1 - kd);

        /* 'convolv_2invG_nov:85' check1=ck1; */
        check1 = kd;
      }

      emxFree_real_T(&d_y);

      /* 'convolv_2invG_nov:87' check2=gm*check1; */
      /* 'convolv_2invG_nov:88' if  check2<=eps/3 */
      if (c_onestagepdf2(T2, m[1], s[1]) * check1 <= 7.4014868308343765E-17) {
        /* 'convolv_2invG_nov:90' flag=1; */
        *flag = 1.0;

        /* 'convolv_2invG_nov:91' sigma=s(2); */
        /* 'convolv_2invG_nov:92' mu=m(2); */
        /* 'convolv_2invG_nov:93' l=1/m(1); */
        /* 'convolv_2invG_nov:95' P=onestagepdf_lag(t,mu,sigma,l); */
        d_onestagepdf_lag(t, m[1], s[1], 1.0 / m[0], P);
      } else {
        /* 'convolv_2invG_nov:97' else */
        /*  find the discrete convolution of the vectors y and z */
        /*  the (i-1)th element of v approximates the convolution of the pdfs  */
        /*  over [.001, x(i)] as a left-hand Riemann sum. */
        /* 'convolv_2invG_nov:102' C=conv(z,y)*h; */
        b_conv(z, c_y, a);
        itmp = a->size[0] * a->size[1];
        a->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
        ixstart = a->size[0];
        nm1d2 = a->size[1];
        nm1d2 *= ixstart;
        for (itmp = 0; itmp < nm1d2; itmp++) {
          a->data[itmp] *= 0.001;
        }

        /* 'convolv_2invG_nov:103' N=length(y); */
        /*  only the first N elements of the convolution are valid */
        /* 'convolv_2invG_nov:105' C=C(1:N); */
        if (1 > c_y->size[1]) {
          nm1d2 = 0;
        } else {
          nm1d2 = c_y->size[1];
        }

        emxInit_real_T(&b_a, 2);
        itmp = b_a->size[0] * b_a->size[1];
        b_a->size[0] = 1;
        b_a->size[1] = a->size[1];
        emxEnsureCapacity((emxArray__common *)b_a, itmp, sizeof(double));
        ixstart = a->size[1];
        for (itmp = 0; itmp < ixstart; itmp++) {
          b_a->data[b_a->size[0] * itmp] = a->data[a->size[0] * itmp];
        }

        itmp = C->size[0] * C->size[1];
        C->size[0] = 1;
        C->size[1] = nm1d2;
        emxEnsureCapacity((emxArray__common *)C, itmp, sizeof(double));
        for (itmp = 0; itmp < nm1d2; itmp++) {
          C->data[C->size[0] * itmp] = b_a->data[itmp];
        }

        emxFree_real_T(&b_a);

        /* 'convolv_2invG_nov:106' I=zeros(1,n); */
        /* 'convolv_2invG_nov:107' P=zeros(1,n); */
        /* 'convolv_2invG_nov:108' for i=1:n */
        emxInit_real_T(&b_t, 2);
        for (i = 0; i < 22001; i++) {
          /* find element of x that is closest to t(i) */
          /* 'convolv_2invG_nov:110' [~,I(i)]=min((t(i)-x).^2); */
          itmp = b_t->size[0] * b_t->size[1];
          b_t->size[0] = 1;
          b_t->size[1] = b_y->size[1];
          emxEnsureCapacity((emxArray__common *)b_t, itmp, sizeof(double));
          nm1d2 = b_y->size[0] * b_y->size[1];
          for (itmp = 0; itmp < nm1d2; itmp++) {
            b_t->data[itmp] = t[i] - b_y->data[itmp];
          }

          e_power(b_t, a);
          ixstart = 1;
          n = a->size[1];
          kd = a->data[0];
          itmp = 1;
          if (a->size[1] > 1) {
            if (rtIsNaN(a->data[0])) {
              nm1d2 = 2;
              exitg1 = false;
              while ((!exitg1) && (nm1d2 <= n)) {
                ixstart = nm1d2;
                if (!rtIsNaN(a->data[nm1d2 - 1])) {
                  kd = a->data[nm1d2 - 1];
                  itmp = nm1d2;
                  exitg1 = true;
                } else {
                  nm1d2++;
                }
              }
            }

            if (ixstart < a->size[1]) {
              while (ixstart + 1 <= n) {
                if (a->data[ixstart] < kd) {
                  kd = a->data[ixstart];
                  itmp = ixstart + 1;
                }

                ixstart++;
              }
            }
          }

          I[i] = (unsigned int)itmp;

          /* 'convolv_2invG_nov:110' ~ */
          /* If t(i)<0 the probability is set to zero, otherwise the */
          /* probability is approxiated as a value from the vector x. */
          /* 'convolv_2invG_nov:113' if t(i)>0 && I(i)>1 */
          if ((t[i] > 0.0) && ((int)I[i] > 1)) {
            /* 'convolv_2invG_nov:114' P(i)=C(I(i)-1); */
            P[i] = C->data[(int)I[i] - 2];
          } else {
            /* 'convolv_2invG_nov:115' else */
            /* 'convolv_2invG_nov:116' P(i)=realmin; */
            P[i] = 2.2250738585072014E-308;
          }
        }

        emxFree_real_T(&b_t);

        /* toc */
        /* 'convolv_2invG_nov:121' P=P'; */
        /* 'convolv_2invG_nov:123' logP=sum(log(P)); */
      }
    } else {
      /* 'convolv_2invG_nov:126' else */
      /*  find the discrete convolution of the vectors y and z */
      /*  the (i-1)th element of v approximates the convolution of the pdfs  */
      /*  over [.001, x(i)] as a left-hand Riemann sum. */
      /* 'convolv_2invG_nov:131' C=conv(z,y)*h; */
      b_conv(z, c_y, a);
      itmp = a->size[0] * a->size[1];
      a->size[0] = 1;
      emxEnsureCapacity((emxArray__common *)a, itmp, sizeof(double));
      ixstart = a->size[0];
      nm1d2 = a->size[1];
      nm1d2 *= ixstart;
      for (itmp = 0; itmp < nm1d2; itmp++) {
        a->data[itmp] *= 0.001;
      }

      /* 'convolv_2invG_nov:132' N=length(y); */
      /*  only the first N elements of the convolution are valid */
      /* 'convolv_2invG_nov:134' C=C(1:N); */
      if (1 > c_y->size[1]) {
        nm1d2 = 0;
      } else {
        nm1d2 = c_y->size[1];
      }

      emxInit_real_T(&b_a, 2);
      itmp = b_a->size[0] * b_a->size[1];
      b_a->size[0] = 1;
      b_a->size[1] = a->size[1];
      emxEnsureCapacity((emxArray__common *)b_a, itmp, sizeof(double));
      ixstart = a->size[1];
      for (itmp = 0; itmp < ixstart; itmp++) {
        b_a->data[b_a->size[0] * itmp] = a->data[a->size[0] * itmp];
      }

      itmp = C->size[0] * C->size[1];
      C->size[0] = 1;
      C->size[1] = nm1d2;
      emxEnsureCapacity((emxArray__common *)C, itmp, sizeof(double));
      for (itmp = 0; itmp < nm1d2; itmp++) {
        C->data[C->size[0] * itmp] = b_a->data[itmp];
      }

      emxFree_real_T(&b_a);

      /* 'convolv_2invG_nov:135' I=zeros(1,n); */
      /* 'convolv_2invG_nov:136' P=zeros(1,n); */
      /* 'convolv_2invG_nov:137' for i=1:n */
      emxInit_real_T(&b_t, 2);
      for (i = 0; i < 22001; i++) {
        /* find element of x that is closest to t(i) */
        /* 'convolv_2invG_nov:139' [~,I(i)]=min((t(i)-x).^2); */
        itmp = b_t->size[0] * b_t->size[1];
        b_t->size[0] = 1;
        b_t->size[1] = b_y->size[1];
        emxEnsureCapacity((emxArray__common *)b_t, itmp, sizeof(double));
        nm1d2 = b_y->size[0] * b_y->size[1];
        for (itmp = 0; itmp < nm1d2; itmp++) {
          b_t->data[itmp] = t[i] - b_y->data[itmp];
        }

        e_power(b_t, a);
        ixstart = 1;
        n = a->size[1];
        kd = a->data[0];
        itmp = 1;
        if (a->size[1] > 1) {
          if (rtIsNaN(a->data[0])) {
            nm1d2 = 2;
            exitg1 = false;
            while ((!exitg1) && (nm1d2 <= n)) {
              ixstart = nm1d2;
              if (!rtIsNaN(a->data[nm1d2 - 1])) {
                kd = a->data[nm1d2 - 1];
                itmp = nm1d2;
                exitg1 = true;
              } else {
                nm1d2++;
              }
            }
          }

          if (ixstart < a->size[1]) {
            while (ixstart + 1 <= n) {
              if (a->data[ixstart] < kd) {
                kd = a->data[ixstart];
                itmp = ixstart + 1;
              }

              ixstart++;
            }
          }
        }

        I[i] = (unsigned int)itmp;

        /* 'convolv_2invG_nov:139' ~ */
        /* If t(i)<0 the probability is set to zero, otherwise the */
        /* probability is approxiated as a value from the vector x. */
        /* 'convolv_2invG_nov:142' if t(i)>0 && I(i)>1 */
        if ((t[i] > 0.0) && ((int)I[i] > 1)) {
          /* 'convolv_2invG_nov:143' P(i)=C(I(i)-1); */
          P[i] = C->data[(int)I[i] - 2];
        } else {
          /* 'convolv_2invG_nov:144' else */
          /* 'convolv_2invG_nov:145' P(i)=realmin; */
          P[i] = 2.2250738585072014E-308;
        }
      }

      emxFree_real_T(&b_t);

      /* toc */
      /* 'convolv_2invG_nov:150' P=P'; */
      /* 'convolv_2invG_nov:152' logP=sum(log(P)); */
    }

    emxFree_real_T(&r12);
    emxFree_real_T(&a);
    emxFree_real_T(&b_y);
    emxFree_real_T(&C);
    emxFree_real_T(&z);
    emxFree_real_T(&c_y);
  }
}

/* End of code generation (convolv_2invG_nov.c) */
