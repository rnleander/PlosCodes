/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * approxconvolv.c
 *
 * Code generation for function 'approxconvolv'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "approxconvolv.h"
#include "power.h"
#include "IMT_analysis_April2017_emxutil.h"
#include "conv.h"
#include "sum.h"
#include "log.h"

/* Function Definitions */

/*
 * function [P0, logP0] = approxconvolv( z, y, h, n, t, I, x )
 */
void approxconvolv(const double z[2201], const double y[2201], const double t
                   [266], const double x[2201], double P0[266], double *logP0)
{
  double C[4401];
  int ixstart;
  int i;
  double dv65[266];
  double b_t[2201];
  double varargin_1[2201];
  double mtmp;
  int itmp;
  int ix;
  boolean_T exitg1;
  short I[266];
  double P[266];

  /* APPROXCONVOLV Summary of this function goes here */
  /*    Detailed explanation goes here */
  /*  BEGIN FUNCTION ApproxConvolv */
  /*  Input parameters: z, y, h, n, t, I, x */
  /*  Outputs: logP0 */
  /*  find the discrete convolution of the vectors y and z */
  /*  the (i-1)th element of v approximates the convolution of the pdfs  */
  /*  over [0, x(i)] as a left-hand Riemann sum. */
  /* 'approxconvolv:14' C=conv(z,y)*h; */
  conv(z, y, C);
  for (ixstart = 0; ixstart < 4401; ixstart++) {
    C[ixstart] *= 0.01;
  }

  /* 'approxconvolv:15' N=length(y); */
  /*  only the first N elements of the convolution are valid */
  /* 'approxconvolv:17' C=C(1:N); */
  /* 'approxconvolv:18' I=zeros(1,n); */
  /* 'approxconvolv:19' P=zeros(1,n); */
  /* 'approxconvolv:20' for i=1:n */
  /* toc */
  /* 'approxconvolv:32' P=P'; */
  /* 'approxconvolv:33' P0=max(realmin,P); */
  /* 'approxconvolv:34' logP0=sum(log(P0)); */
  for (i = 0; i < 266; i++) {
    /* find element of x that is closest to t(i) */
    /* 'approxconvolv:22' [~,I(i)]=min((t(i)-x).^2); */
    for (ixstart = 0; ixstart < 2201; ixstart++) {
      b_t[ixstart] = t[i] - x[ixstart];
    }

    d_power(b_t, varargin_1);
    ixstart = 1;
    mtmp = varargin_1[0];
    itmp = 1;
    if (rtIsNaN(varargin_1[0])) {
      ix = 2;
      exitg1 = false;
      while ((!exitg1) && (ix < 2202)) {
        ixstart = ix;
        if (!rtIsNaN(varargin_1[ix - 1])) {
          mtmp = varargin_1[ix - 1];
          itmp = ix;
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < 2201) {
      while (ixstart + 1 < 2202) {
        if (varargin_1[ixstart] < mtmp) {
          mtmp = varargin_1[ixstart];
          itmp = ixstart + 1;
        }

        ixstart++;
      }
    }

    I[i] = (short)itmp;

    /* 'approxconvolv:22' ~ */
    /* If t(i)<0 the probability is set to zero, otherwise the */
    /* probability is approxiated as a value from the vector x. */
    /* 'approxconvolv:25' if t(i)>0 && I(i)>1 */
    if (I[i] > 1) {
      /* 'approxconvolv:26' P(i)=C(I(i)-1); */
      P[i] = C[I[i] - 2];
    } else {
      /* 'approxconvolv:27' else */
      /* 'approxconvolv:28' P(i)=realmin; */
      P[i] = 2.2250738585072014E-308;
    }

    if ((2.2250738585072014E-308 > P[i]) || rtIsNaN(P[i])) {
      P0[i] = 2.2250738585072014E-308;
    } else {
      P0[i] = P[i];
    }

    dv65[i] = P0[i];
  }

  b_log(dv65);
  *logP0 = sum(dv65);

  /*  END FUNCTION DOTHECONVOLUTION_INNER */
}

/*
 * function [P0, logP0] = approxconvolv( z, y, h, n, t, I, x )
 */
void b_approxconvolv(const emxArray_real_T *z, const emxArray_real_T *y, double
                     h, const emxArray_real_T *x, double P0[266], double *logP0)
{
  emxArray_real_T *C;
  int n;
  int ixstart;
  emxArray_real_T *varargin_1;
  emxArray_real_T *t;
  int i;
  double b_x[266];
  double P[266];
  static const double b_t[266] = { 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
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

  double extremum;
  int itmp;
  unsigned int I[266];
  int ix;
  boolean_T exitg1;
  emxInit_real_T(&C, 2);

  /* APPROXCONVOLV Summary of this function goes here */
  /*    Detailed explanation goes here */
  /*  BEGIN FUNCTION ApproxConvolv */
  /*  Input parameters: z, y, h, n, t, I, x */
  /*  Outputs: logP0 */
  /*  find the discrete convolution of the vectors y and z */
  /*  the (i-1)th element of v approximates the convolution of the pdfs  */
  /*  over [0, x(i)] as a left-hand Riemann sum. */
  /* 'approxconvolv:14' C=conv(z,y)*h; */
  b_conv(z, y, C);
  n = C->size[0] * C->size[1];
  C->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)C, n, sizeof(double));
  ixstart = C->size[0];
  n = C->size[1];
  ixstart *= n;
  for (n = 0; n < ixstart; n++) {
    C->data[n] *= h;
  }

  /* 'approxconvolv:15' N=length(y); */
  /*  only the first N elements of the convolution are valid */
  /* 'approxconvolv:17' C=C(1:N); */
  if (1 > y->size[1]) {
    n = 0;
  } else {
    n = y->size[1];
  }

  ixstart = C->size[0] * C->size[1];
  C->size[1] = n;
  emxEnsureCapacity((emxArray__common *)C, ixstart, sizeof(double));

  /* 'approxconvolv:18' I=zeros(1,n); */
  /* 'approxconvolv:19' P=zeros(1,n); */
  /* 'approxconvolv:20' for i=1:n */
  emxInit_real_T(&varargin_1, 2);
  emxInit_real_T(&t, 2);
  for (i = 0; i < 266; i++) {
    /* find element of x that is closest to t(i) */
    /* 'approxconvolv:22' [~,I(i)]=min((t(i)-x).^2); */
    n = t->size[0] * t->size[1];
    t->size[0] = 1;
    t->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)t, n, sizeof(double));
    ixstart = x->size[0] * x->size[1];
    for (n = 0; n < ixstart; n++) {
      t->data[n] = b_t[i] - x->data[n];
    }

    e_power(t, varargin_1);
    ixstart = 1;
    n = varargin_1->size[1];
    extremum = varargin_1->data[0];
    itmp = 1;
    if (varargin_1->size[1] > 1) {
      if (rtIsNaN(varargin_1->data[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix <= n)) {
          ixstart = ix;
          if (!rtIsNaN(varargin_1->data[ix - 1])) {
            extremum = varargin_1->data[ix - 1];
            itmp = ix;
            exitg1 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < varargin_1->size[1]) {
        while (ixstart + 1 <= n) {
          if (varargin_1->data[ixstart] < extremum) {
            extremum = varargin_1->data[ixstart];
            itmp = ixstart + 1;
          }

          ixstart++;
        }
      }
    }

    I[i] = (unsigned int)itmp;

    /* 'approxconvolv:22' ~ */
    /* If t(i)<0 the probability is set to zero, otherwise the */
    /* probability is approxiated as a value from the vector x. */
    /* 'approxconvolv:25' if t(i)>0 && I(i)>1 */
    if ((int)I[i] > 1) {
      /* 'approxconvolv:26' P(i)=C(I(i)-1); */
      P[i] = C->data[(int)I[i] - 2];
    } else {
      /* 'approxconvolv:27' else */
      /* 'approxconvolv:28' P(i)=realmin; */
      P[i] = 2.2250738585072014E-308;
    }
  }

  emxFree_real_T(&t);
  emxFree_real_T(&varargin_1);
  emxFree_real_T(&C);

  /* toc */
  /* 'approxconvolv:32' P=P'; */
  /* 'approxconvolv:33' P0=max(realmin,P); */
  /* 'approxconvolv:34' logP0=sum(log(P0)); */
  for (ixstart = 0; ixstart < 266; ixstart++) {
    if ((2.2250738585072014E-308 > P[ixstart]) || rtIsNaN(P[ixstart])) {
      extremum = 2.2250738585072014E-308;
    } else {
      extremum = P[ixstart];
    }

    P0[ixstart] = extremum;
    b_x[ixstart] = log(extremum);
  }

  *logP0 = b_x[0];
  for (ixstart = 0; ixstart < 265; ixstart++) {
    *logP0 += b_x[ixstart + 1];
  }

  /*  END FUNCTION DOTHECONVOLUTION_INNER */
}

/*
 * function [P0, logP0] = approxconvolv( z, y, h, n, t, I, x )
 */
void c_approxconvolv(const double z[22001], const double y[22001], const double
                     t[266], const double x[22001], double P0[266], double
                     *logP0)
{
  static double C[44001];
  int ixstart;
  int i;
  double dv69[266];
  static double b_t[22001];
  double varargin_1[22001];
  double mtmp;
  int itmp;
  int ix;
  boolean_T exitg1;
  short I[266];
  double P[266];

  /* APPROXCONVOLV Summary of this function goes here */
  /*    Detailed explanation goes here */
  /*  BEGIN FUNCTION ApproxConvolv */
  /*  Input parameters: z, y, h, n, t, I, x */
  /*  Outputs: logP0 */
  /*  find the discrete convolution of the vectors y and z */
  /*  the (i-1)th element of v approximates the convolution of the pdfs  */
  /*  over [0, x(i)] as a left-hand Riemann sum. */
  /* 'approxconvolv:14' C=conv(z,y)*h; */
  c_conv(z, y, C);
  for (ixstart = 0; ixstart < 44001; ixstart++) {
    C[ixstart] *= 0.001;
  }

  /* 'approxconvolv:15' N=length(y); */
  /*  only the first N elements of the convolution are valid */
  /* 'approxconvolv:17' C=C(1:N); */
  /* 'approxconvolv:18' I=zeros(1,n); */
  /* 'approxconvolv:19' P=zeros(1,n); */
  /* 'approxconvolv:20' for i=1:n */
  /* toc */
  /* 'approxconvolv:32' P=P'; */
  /* 'approxconvolv:33' P0=max(realmin,P); */
  /* 'approxconvolv:34' logP0=sum(log(P0)); */
  for (i = 0; i < 266; i++) {
    /* find element of x that is closest to t(i) */
    /* 'approxconvolv:22' [~,I(i)]=min((t(i)-x).^2); */
    for (ixstart = 0; ixstart < 22001; ixstart++) {
      b_t[ixstart] = t[i] - x[ixstart];
    }

    f_power(b_t, varargin_1);
    ixstart = 1;
    mtmp = varargin_1[0];
    itmp = 1;
    if (rtIsNaN(varargin_1[0])) {
      ix = 2;
      exitg1 = false;
      while ((!exitg1) && (ix < 22002)) {
        ixstart = ix;
        if (!rtIsNaN(varargin_1[ix - 1])) {
          mtmp = varargin_1[ix - 1];
          itmp = ix;
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < 22001) {
      while (ixstart + 1 < 22002) {
        if (varargin_1[ixstart] < mtmp) {
          mtmp = varargin_1[ixstart];
          itmp = ixstart + 1;
        }

        ixstart++;
      }
    }

    I[i] = (short)itmp;

    /* 'approxconvolv:22' ~ */
    /* If t(i)<0 the probability is set to zero, otherwise the */
    /* probability is approxiated as a value from the vector x. */
    /* 'approxconvolv:25' if t(i)>0 && I(i)>1 */
    if (I[i] > 1) {
      /* 'approxconvolv:26' P(i)=C(I(i)-1); */
      P[i] = C[I[i] - 2];
    } else {
      /* 'approxconvolv:27' else */
      /* 'approxconvolv:28' P(i)=realmin; */
      P[i] = 2.2250738585072014E-308;
    }

    if ((2.2250738585072014E-308 > P[i]) || rtIsNaN(P[i])) {
      P0[i] = 2.2250738585072014E-308;
    } else {
      P0[i] = P[i];
    }

    dv69[i] = P0[i];
  }

  b_log(dv69);
  *logP0 = sum(dv69);

  /*  END FUNCTION DOTHECONVOLUTION_INNER */
}

/*
 * function [P0, logP0] = approxconvolv( z, y, h, n, t, I, x )
 */
void d_approxconvolv(const emxArray_real_T *z, const emxArray_real_T *y, const
                     double t[266], const emxArray_real_T *x, double P0[266],
                     double *logP0)
{
  emxArray_real_T *C;
  int n;
  int ixstart;
  emxArray_real_T *varargin_1;
  emxArray_real_T *b_t;
  int i;
  double b_x[266];
  double P[266];
  double extremum;
  int itmp;
  unsigned int I[266];
  int ix;
  boolean_T exitg1;
  emxInit_real_T(&C, 2);

  /* APPROXCONVOLV Summary of this function goes here */
  /*    Detailed explanation goes here */
  /*  BEGIN FUNCTION ApproxConvolv */
  /*  Input parameters: z, y, h, n, t, I, x */
  /*  Outputs: logP0 */
  /*  find the discrete convolution of the vectors y and z */
  /*  the (i-1)th element of v approximates the convolution of the pdfs  */
  /*  over [0, x(i)] as a left-hand Riemann sum. */
  /* 'approxconvolv:14' C=conv(z,y)*h; */
  b_conv(z, y, C);
  n = C->size[0] * C->size[1];
  C->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)C, n, sizeof(double));
  ixstart = C->size[0];
  n = C->size[1];
  ixstart *= n;
  for (n = 0; n < ixstart; n++) {
    C->data[n] *= 0.1;
  }

  /* 'approxconvolv:15' N=length(y); */
  /*  only the first N elements of the convolution are valid */
  /* 'approxconvolv:17' C=C(1:N); */
  if (1 > y->size[1]) {
    n = 0;
  } else {
    n = y->size[1];
  }

  ixstart = C->size[0] * C->size[1];
  C->size[1] = n;
  emxEnsureCapacity((emxArray__common *)C, ixstart, sizeof(double));

  /* 'approxconvolv:18' I=zeros(1,n); */
  /* 'approxconvolv:19' P=zeros(1,n); */
  /* 'approxconvolv:20' for i=1:n */
  emxInit_real_T(&varargin_1, 2);
  emxInit_real_T(&b_t, 2);
  for (i = 0; i < 266; i++) {
    /* find element of x that is closest to t(i) */
    /* 'approxconvolv:22' [~,I(i)]=min((t(i)-x).^2); */
    n = b_t->size[0] * b_t->size[1];
    b_t->size[0] = 1;
    b_t->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)b_t, n, sizeof(double));
    ixstart = x->size[0] * x->size[1];
    for (n = 0; n < ixstart; n++) {
      b_t->data[n] = t[i] - x->data[n];
    }

    e_power(b_t, varargin_1);
    ixstart = 1;
    n = varargin_1->size[1];
    extremum = varargin_1->data[0];
    itmp = 1;
    if (varargin_1->size[1] > 1) {
      if (rtIsNaN(varargin_1->data[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix <= n)) {
          ixstart = ix;
          if (!rtIsNaN(varargin_1->data[ix - 1])) {
            extremum = varargin_1->data[ix - 1];
            itmp = ix;
            exitg1 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < varargin_1->size[1]) {
        while (ixstart + 1 <= n) {
          if (varargin_1->data[ixstart] < extremum) {
            extremum = varargin_1->data[ixstart];
            itmp = ixstart + 1;
          }

          ixstart++;
        }
      }
    }

    I[i] = (unsigned int)itmp;

    /* 'approxconvolv:22' ~ */
    /* If t(i)<0 the probability is set to zero, otherwise the */
    /* probability is approxiated as a value from the vector x. */
    /* 'approxconvolv:25' if t(i)>0 && I(i)>1 */
    if ((t[i] > 0.0) && ((int)I[i] > 1)) {
      /* 'approxconvolv:26' P(i)=C(I(i)-1); */
      P[i] = C->data[(int)I[i] - 2];
    } else {
      /* 'approxconvolv:27' else */
      /* 'approxconvolv:28' P(i)=realmin; */
      P[i] = 2.2250738585072014E-308;
    }
  }

  emxFree_real_T(&b_t);
  emxFree_real_T(&varargin_1);
  emxFree_real_T(&C);

  /* toc */
  /* 'approxconvolv:32' P=P'; */
  /* 'approxconvolv:33' P0=max(realmin,P); */
  /* 'approxconvolv:34' logP0=sum(log(P0)); */
  for (ixstart = 0; ixstart < 266; ixstart++) {
    if ((2.2250738585072014E-308 > P[ixstart]) || rtIsNaN(P[ixstart])) {
      extremum = 2.2250738585072014E-308;
    } else {
      extremum = P[ixstart];
    }

    P0[ixstart] = extremum;
    b_x[ixstart] = log(extremum);
  }

  *logP0 = b_x[0];
  for (ixstart = 0; ixstart < 265; ixstart++) {
    *logP0 += b_x[ixstart + 1];
  }

  /*  END FUNCTION DOTHECONVOLUTION_INNER */
}

/*
 * function [P0, logP0] = approxconvolv( z, y, h, n, t, I, x )
 */
void e_approxconvolv(const emxArray_real_T *z, const emxArray_real_T *y, double
                     h, const double t[266], const emxArray_real_T *x, double
                     P0[266], double *logP0)
{
  emxArray_real_T *C;
  int n;
  int ixstart;
  emxArray_real_T *varargin_1;
  emxArray_real_T *b_t;
  int i;
  double b_x[266];
  double P[266];
  double extremum;
  int itmp;
  unsigned int I[266];
  int ix;
  boolean_T exitg1;
  emxInit_real_T(&C, 2);

  /* APPROXCONVOLV Summary of this function goes here */
  /*    Detailed explanation goes here */
  /*  BEGIN FUNCTION ApproxConvolv */
  /*  Input parameters: z, y, h, n, t, I, x */
  /*  Outputs: logP0 */
  /*  find the discrete convolution of the vectors y and z */
  /*  the (i-1)th element of v approximates the convolution of the pdfs  */
  /*  over [0, x(i)] as a left-hand Riemann sum. */
  /* 'approxconvolv:14' C=conv(z,y)*h; */
  b_conv(z, y, C);
  n = C->size[0] * C->size[1];
  C->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)C, n, sizeof(double));
  ixstart = C->size[0];
  n = C->size[1];
  ixstart *= n;
  for (n = 0; n < ixstart; n++) {
    C->data[n] *= h;
  }

  /* 'approxconvolv:15' N=length(y); */
  /*  only the first N elements of the convolution are valid */
  /* 'approxconvolv:17' C=C(1:N); */
  if (1 > y->size[1]) {
    n = 0;
  } else {
    n = y->size[1];
  }

  ixstart = C->size[0] * C->size[1];
  C->size[1] = n;
  emxEnsureCapacity((emxArray__common *)C, ixstart, sizeof(double));

  /* 'approxconvolv:18' I=zeros(1,n); */
  /* 'approxconvolv:19' P=zeros(1,n); */
  /* 'approxconvolv:20' for i=1:n */
  emxInit_real_T(&varargin_1, 2);
  emxInit_real_T(&b_t, 2);
  for (i = 0; i < 266; i++) {
    /* find element of x that is closest to t(i) */
    /* 'approxconvolv:22' [~,I(i)]=min((t(i)-x).^2); */
    n = b_t->size[0] * b_t->size[1];
    b_t->size[0] = 1;
    b_t->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)b_t, n, sizeof(double));
    ixstart = x->size[0] * x->size[1];
    for (n = 0; n < ixstart; n++) {
      b_t->data[n] = t[i] - x->data[n];
    }

    e_power(b_t, varargin_1);
    ixstart = 1;
    n = varargin_1->size[1];
    extremum = varargin_1->data[0];
    itmp = 1;
    if (varargin_1->size[1] > 1) {
      if (rtIsNaN(varargin_1->data[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix <= n)) {
          ixstart = ix;
          if (!rtIsNaN(varargin_1->data[ix - 1])) {
            extremum = varargin_1->data[ix - 1];
            itmp = ix;
            exitg1 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < varargin_1->size[1]) {
        while (ixstart + 1 <= n) {
          if (varargin_1->data[ixstart] < extremum) {
            extremum = varargin_1->data[ixstart];
            itmp = ixstart + 1;
          }

          ixstart++;
        }
      }
    }

    I[i] = (unsigned int)itmp;

    /* 'approxconvolv:22' ~ */
    /* If t(i)<0 the probability is set to zero, otherwise the */
    /* probability is approxiated as a value from the vector x. */
    /* 'approxconvolv:25' if t(i)>0 && I(i)>1 */
    if ((t[i] > 0.0) && ((int)I[i] > 1)) {
      /* 'approxconvolv:26' P(i)=C(I(i)-1); */
      P[i] = C->data[(int)I[i] - 2];
    } else {
      /* 'approxconvolv:27' else */
      /* 'approxconvolv:28' P(i)=realmin; */
      P[i] = 2.2250738585072014E-308;
    }
  }

  emxFree_real_T(&b_t);
  emxFree_real_T(&varargin_1);
  emxFree_real_T(&C);

  /* toc */
  /* 'approxconvolv:32' P=P'; */
  /* 'approxconvolv:33' P0=max(realmin,P); */
  /* 'approxconvolv:34' logP0=sum(log(P0)); */
  for (ixstart = 0; ixstart < 266; ixstart++) {
    if ((2.2250738585072014E-308 > P[ixstart]) || rtIsNaN(P[ixstart])) {
      extremum = 2.2250738585072014E-308;
    } else {
      extremum = P[ixstart];
    }

    P0[ixstart] = extremum;
    b_x[ixstart] = log(extremum);
  }

  *logP0 = b_x[0];
  for (ixstart = 0; ixstart < 265; ixstart++) {
    *logP0 += b_x[ixstart + 1];
  }

  /*  END FUNCTION DOTHECONVOLUTION_INNER */
}

/*
 * function [P0, logP0] = approxconvolv( z, y, h, n, t, I, x )
 */
void f_approxconvolv(const emxArray_real_T *z, const emxArray_real_T *y, const
                     double t[266], const emxArray_real_T *x, double P0[266],
                     double *logP0)
{
  emxArray_real_T *C;
  int n;
  int ixstart;
  emxArray_real_T *varargin_1;
  emxArray_real_T *b_t;
  int i;
  double b_x[266];
  double P[266];
  double extremum;
  int itmp;
  unsigned int I[266];
  int ix;
  boolean_T exitg1;
  emxInit_real_T(&C, 2);

  /* APPROXCONVOLV Summary of this function goes here */
  /*    Detailed explanation goes here */
  /*  BEGIN FUNCTION ApproxConvolv */
  /*  Input parameters: z, y, h, n, t, I, x */
  /*  Outputs: logP0 */
  /*  find the discrete convolution of the vectors y and z */
  /*  the (i-1)th element of v approximates the convolution of the pdfs  */
  /*  over [0, x(i)] as a left-hand Riemann sum. */
  /* 'approxconvolv:14' C=conv(z,y)*h; */
  b_conv(z, y, C);
  n = C->size[0] * C->size[1];
  C->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)C, n, sizeof(double));
  ixstart = C->size[0];
  n = C->size[1];
  ixstart *= n;
  for (n = 0; n < ixstart; n++) {
    C->data[n] *= 0.01;
  }

  /* 'approxconvolv:15' N=length(y); */
  /*  only the first N elements of the convolution are valid */
  /* 'approxconvolv:17' C=C(1:N); */
  if (1 > y->size[1]) {
    n = 0;
  } else {
    n = y->size[1];
  }

  ixstart = C->size[0] * C->size[1];
  C->size[1] = n;
  emxEnsureCapacity((emxArray__common *)C, ixstart, sizeof(double));

  /* 'approxconvolv:18' I=zeros(1,n); */
  /* 'approxconvolv:19' P=zeros(1,n); */
  /* 'approxconvolv:20' for i=1:n */
  emxInit_real_T(&varargin_1, 2);
  emxInit_real_T(&b_t, 2);
  for (i = 0; i < 266; i++) {
    /* find element of x that is closest to t(i) */
    /* 'approxconvolv:22' [~,I(i)]=min((t(i)-x).^2); */
    n = b_t->size[0] * b_t->size[1];
    b_t->size[0] = 1;
    b_t->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)b_t, n, sizeof(double));
    ixstart = x->size[0] * x->size[1];
    for (n = 0; n < ixstart; n++) {
      b_t->data[n] = t[i] - x->data[n];
    }

    e_power(b_t, varargin_1);
    ixstart = 1;
    n = varargin_1->size[1];
    extremum = varargin_1->data[0];
    itmp = 1;
    if (varargin_1->size[1] > 1) {
      if (rtIsNaN(varargin_1->data[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix <= n)) {
          ixstart = ix;
          if (!rtIsNaN(varargin_1->data[ix - 1])) {
            extremum = varargin_1->data[ix - 1];
            itmp = ix;
            exitg1 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < varargin_1->size[1]) {
        while (ixstart + 1 <= n) {
          if (varargin_1->data[ixstart] < extremum) {
            extremum = varargin_1->data[ixstart];
            itmp = ixstart + 1;
          }

          ixstart++;
        }
      }
    }

    I[i] = (unsigned int)itmp;

    /* 'approxconvolv:22' ~ */
    /* If t(i)<0 the probability is set to zero, otherwise the */
    /* probability is approxiated as a value from the vector x. */
    /* 'approxconvolv:25' if t(i)>0 && I(i)>1 */
    if ((t[i] > 0.0) && ((int)I[i] > 1)) {
      /* 'approxconvolv:26' P(i)=C(I(i)-1); */
      P[i] = C->data[(int)I[i] - 2];
    } else {
      /* 'approxconvolv:27' else */
      /* 'approxconvolv:28' P(i)=realmin; */
      P[i] = 2.2250738585072014E-308;
    }
  }

  emxFree_real_T(&b_t);
  emxFree_real_T(&varargin_1);
  emxFree_real_T(&C);

  /* toc */
  /* 'approxconvolv:32' P=P'; */
  /* 'approxconvolv:33' P0=max(realmin,P); */
  /* 'approxconvolv:34' logP0=sum(log(P0)); */
  for (ixstart = 0; ixstart < 266; ixstart++) {
    if ((2.2250738585072014E-308 > P[ixstart]) || rtIsNaN(P[ixstart])) {
      extremum = 2.2250738585072014E-308;
    } else {
      extremum = P[ixstart];
    }

    P0[ixstart] = extremum;
    b_x[ixstart] = log(extremum);
  }

  *logP0 = b_x[0];
  for (ixstart = 0; ixstart < 265; ixstart++) {
    *logP0 += b_x[ixstart + 1];
  }

  /*  END FUNCTION DOTHECONVOLUTION_INNER */
}

/*
 * function [P0, logP0] = approxconvolv( z, y, h, n, t, I, x )
 */
void g_approxconvolv(const emxArray_real_T *z, const emxArray_real_T *y, const
                     double t[266], const emxArray_real_T *x, double P0[266],
                     double *logP0)
{
  emxArray_real_T *C;
  int n;
  int ixstart;
  emxArray_real_T *varargin_1;
  emxArray_real_T *b_t;
  int i;
  double b_x[266];
  double P[266];
  double extremum;
  int itmp;
  unsigned int I[266];
  int ix;
  boolean_T exitg1;
  emxInit_real_T(&C, 2);

  /* APPROXCONVOLV Summary of this function goes here */
  /*    Detailed explanation goes here */
  /*  BEGIN FUNCTION ApproxConvolv */
  /*  Input parameters: z, y, h, n, t, I, x */
  /*  Outputs: logP0 */
  /*  find the discrete convolution of the vectors y and z */
  /*  the (i-1)th element of v approximates the convolution of the pdfs  */
  /*  over [0, x(i)] as a left-hand Riemann sum. */
  /* 'approxconvolv:14' C=conv(z,y)*h; */
  b_conv(z, y, C);
  n = C->size[0] * C->size[1];
  C->size[0] = 1;
  emxEnsureCapacity((emxArray__common *)C, n, sizeof(double));
  ixstart = C->size[0];
  n = C->size[1];
  ixstart *= n;
  for (n = 0; n < ixstart; n++) {
    C->data[n] *= 0.001;
  }

  /* 'approxconvolv:15' N=length(y); */
  /*  only the first N elements of the convolution are valid */
  /* 'approxconvolv:17' C=C(1:N); */
  if (1 > y->size[1]) {
    n = 0;
  } else {
    n = y->size[1];
  }

  ixstart = C->size[0] * C->size[1];
  C->size[1] = n;
  emxEnsureCapacity((emxArray__common *)C, ixstart, sizeof(double));

  /* 'approxconvolv:18' I=zeros(1,n); */
  /* 'approxconvolv:19' P=zeros(1,n); */
  /* 'approxconvolv:20' for i=1:n */
  emxInit_real_T(&varargin_1, 2);
  emxInit_real_T(&b_t, 2);
  for (i = 0; i < 266; i++) {
    /* find element of x that is closest to t(i) */
    /* 'approxconvolv:22' [~,I(i)]=min((t(i)-x).^2); */
    n = b_t->size[0] * b_t->size[1];
    b_t->size[0] = 1;
    b_t->size[1] = x->size[1];
    emxEnsureCapacity((emxArray__common *)b_t, n, sizeof(double));
    ixstart = x->size[0] * x->size[1];
    for (n = 0; n < ixstart; n++) {
      b_t->data[n] = t[i] - x->data[n];
    }

    e_power(b_t, varargin_1);
    ixstart = 1;
    n = varargin_1->size[1];
    extremum = varargin_1->data[0];
    itmp = 1;
    if (varargin_1->size[1] > 1) {
      if (rtIsNaN(varargin_1->data[0])) {
        ix = 2;
        exitg1 = false;
        while ((!exitg1) && (ix <= n)) {
          ixstart = ix;
          if (!rtIsNaN(varargin_1->data[ix - 1])) {
            extremum = varargin_1->data[ix - 1];
            itmp = ix;
            exitg1 = true;
          } else {
            ix++;
          }
        }
      }

      if (ixstart < varargin_1->size[1]) {
        while (ixstart + 1 <= n) {
          if (varargin_1->data[ixstart] < extremum) {
            extremum = varargin_1->data[ixstart];
            itmp = ixstart + 1;
          }

          ixstart++;
        }
      }
    }

    I[i] = (unsigned int)itmp;

    /* 'approxconvolv:22' ~ */
    /* If t(i)<0 the probability is set to zero, otherwise the */
    /* probability is approxiated as a value from the vector x. */
    /* 'approxconvolv:25' if t(i)>0 && I(i)>1 */
    if ((t[i] > 0.0) && ((int)I[i] > 1)) {
      /* 'approxconvolv:26' P(i)=C(I(i)-1); */
      P[i] = C->data[(int)I[i] - 2];
    } else {
      /* 'approxconvolv:27' else */
      /* 'approxconvolv:28' P(i)=realmin; */
      P[i] = 2.2250738585072014E-308;
    }
  }

  emxFree_real_T(&b_t);
  emxFree_real_T(&varargin_1);
  emxFree_real_T(&C);

  /* toc */
  /* 'approxconvolv:32' P=P'; */
  /* 'approxconvolv:33' P0=max(realmin,P); */
  /* 'approxconvolv:34' logP0=sum(log(P0)); */
  for (ixstart = 0; ixstart < 266; ixstart++) {
    if ((2.2250738585072014E-308 > P[ixstart]) || rtIsNaN(P[ixstart])) {
      extremum = 2.2250738585072014E-308;
    } else {
      extremum = P[ixstart];
    }

    P0[ixstart] = extremum;
    b_x[ixstart] = log(extremum);
  }

  *logP0 = b_x[0];
  for (ixstart = 0; ixstart < 265; ixstart++) {
    *logP0 += b_x[ixstart + 1];
  }

  /*  END FUNCTION DOTHECONVOLUTION_INNER */
}

/* End of code generation (approxconvolv.c) */
