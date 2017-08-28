/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * power.c
 *
 * Code generation for function 'power'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "power.h"
#include "IMT_analysis_April2017_emxutil.h"
#include "emgpdf.h"
#include "IMT_analysis_April2017_rtwutil.h"

/* Function Definitions */

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void b_power(const double a[2], double y[2])
{
  int k;
  for (k = 0; k < 2; k++) {
    y[k] = rt_powd_snf(a[k], 3.0);
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void c_power(const double a[2], double y[2])
{
  int k;
  for (k = 0; k < 2; k++) {
    y[k] = rt_powd_snf(a[k], 0.5);
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void d_power(const double a[2201], double y[2201])
{
  int k;
  for (k = 0; k < 2201; k++) {
    y[k] = rt_powd_snf(a[k], 2.0);
  }
}
#endif

/*
 *
 */
void e_power(const emxArray_real_T *a, emxArray_real_T *y)
{
  int k;
  k = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = a->size[1];
  emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  for (k = 0; k + 1 <= a->size[1]; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void f_power(const double a[22001], double y[22001])
{
  int k;
  for (k = 0; k < 22001; k++) {
    y[k] = rt_powd_snf(a[k], 2.0);
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void g_power(const double a[3], double y[3])
{
  int k;
  for (k = 0; k < 3; k++) {
    y[k] = rt_powd_snf(a[k], 2.0);
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void h_power(const double a[3], double y[3])
{
  int k;
  for (k = 0; k < 3; k++) {
    y[k] = rt_powd_snf(a[k], 3.0);
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void i_power(const double a[3], double y[3])
{
  int k;
  for (k = 0; k < 3; k++) {
    y[k] = rt_powd_snf(a[k], 0.5);
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void j_power(const double a[221], double y[221])
{
  int k;
  for (k = 0; k < 221; k++) {
    y[k] = rt_powd_snf(a[k], 2.0);
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void k_power(const emxArray_real_T *a, emxArray_real_T *y)
{
  unsigned int a_idx_0;
  int k;
  a_idx_0 = (unsigned int)a->size[0];
  k = y->size[0];
  y->size[0] = (int)a_idx_0;
  emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  for (k = 0; k + 1 <= a->size[0]; k++) {
    y->data[k] = rt_powd_snf(a->data[k], 3.0);
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void l_power(const emxArray_real_T *a, emxArray_real_T *y)
{
  unsigned int a_idx_0;
  int k;
  a_idx_0 = (unsigned int)a->size[0];
  k = y->size[0];
  y->size[0] = (int)a_idx_0;
  emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  for (k = 0; k + 1 <= a->size[0]; k++) {
    y->data[k] = sqrt(a->data[k]);
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void m_power(const emxArray_real_T *a, emxArray_real_T *y)
{
  unsigned int a_idx_0;
  int k;
  a_idx_0 = (unsigned int)a->size[0];
  k = y->size[0];
  y->size[0] = (int)a_idx_0;
  emxEnsureCapacity((emxArray__common *)y, k, sizeof(double));
  for (k = 0; k + 1 <= a->size[0]; k++) {
    y->data[k] = a->data[k] * a->data[k];
  }
}
#endif

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
void power(const double a[2], double y[2])
{
  int k;
  for (k = 0; k < 2; k++) {
    y[k] = rt_powd_snf(a[k], 2.0);
  }
}
#endif

/* End of code generation (power.c) */
