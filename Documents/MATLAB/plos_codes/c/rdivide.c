/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * rdivide.c
 *
 * Code generation for function 'rdivide'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "rdivide.h"
#include "IMT_analysis_April2017_emxutil.h"

/* Function Definitions */

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void b_rdivide(const emxArray_real_T *y, emxArray_real_T *z)
{
  int i2;
  int loop_ub;
  i2 = z->size[0];
  z->size[0] = y->size[0];
  emxEnsureCapacity((emxArray__common *)z, i2, sizeof(double));
  loop_ub = y->size[0];
  for (i2 = 0; i2 < loop_ub; i2++) {
    z->data[i2] = 1.0 / y->data[i2];
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void c_rdivide(const emxArray_real_T *x, const emxArray_real_T *y,
               emxArray_real_T *z)
{
  int i3;
  int loop_ub;
  i3 = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)z, i3, sizeof(double));
  loop_ub = x->size[0];
  for (i3 = 0; i3 < loop_ub; i3++) {
    z->data[i3] = x->data[i3] / y->data[i3];
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void rdivide(const creal_T x_data[], const int x_size[1], const creal_T y_data[],
             creal_T z_data[], int z_size[1])
{
  int loop_ub;
  int i0;
  double brm;
  double bim;
  double d;
  z_size[0] = x_size[0];
  loop_ub = x_size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    if (y_data[i0].im == 0.0) {
      if (x_data[i0].im == 0.0) {
        z_data[i0].re = x_data[i0].re / y_data[i0].re;
        z_data[i0].im = 0.0;
      } else if (x_data[i0].re == 0.0) {
        z_data[i0].re = 0.0;
        z_data[i0].im = x_data[i0].im / y_data[i0].re;
      } else {
        z_data[i0].re = x_data[i0].re / y_data[i0].re;
        z_data[i0].im = x_data[i0].im / y_data[i0].re;
      }
    } else if (y_data[i0].re == 0.0) {
      if (x_data[i0].re == 0.0) {
        z_data[i0].re = x_data[i0].im / y_data[i0].im;
        z_data[i0].im = 0.0;
      } else if (x_data[i0].im == 0.0) {
        z_data[i0].re = 0.0;
        z_data[i0].im = -(x_data[i0].re / y_data[i0].im);
      } else {
        z_data[i0].re = x_data[i0].im / y_data[i0].im;
        z_data[i0].im = -(x_data[i0].re / y_data[i0].im);
      }
    } else {
      brm = fabs(y_data[i0].re);
      bim = fabs(y_data[i0].im);
      if (brm > bim) {
        bim = y_data[i0].im / y_data[i0].re;
        d = y_data[i0].re + bim * y_data[i0].im;
        z_data[i0].re = (x_data[i0].re + bim * x_data[i0].im) / d;
        z_data[i0].im = (x_data[i0].im - bim * x_data[i0].re) / d;
      } else if (bim == brm) {
        if (y_data[i0].re > 0.0) {
          bim = 0.5;
        } else {
          bim = -0.5;
        }

        if (y_data[i0].im > 0.0) {
          d = 0.5;
        } else {
          d = -0.5;
        }

        z_data[i0].re = (x_data[i0].re * bim + x_data[i0].im * d) / brm;
        z_data[i0].im = (x_data[i0].im * bim - x_data[i0].re * d) / brm;
      } else {
        bim = y_data[i0].re / y_data[i0].im;
        d = y_data[i0].im + bim * y_data[i0].re;
        z_data[i0].re = (bim * x_data[i0].re + x_data[i0].im) / d;
        z_data[i0].im = (bim * x_data[i0].im - x_data[i0].re) / d;
      }
    }
  }
}
#endif

/* End of code generation (rdivide.c) */
