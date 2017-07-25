/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * conv.c
 *
 * Code generation for function 'conv'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "conv.h"
#include "IMT_analysis_April2017_emxutil.h"

/* Function Definitions */

/*
 *
 */
void b_conv(const emxArray_real_T *A, const emxArray_real_T *B, emxArray_real_T *
            C)
{
  int nA;
  int nB;
  int nApnB;
  int k;
  nA = A->size[1];
  nB = B->size[1];
  nApnB = A->size[1] + B->size[1];
  if ((A->size[1] == 0) || (B->size[1] == 0)) {
  } else {
    nApnB--;
  }

  k = C->size[0] * C->size[1];
  C->size[0] = 1;
  C->size[1] = nApnB;
  emxEnsureCapacity((emxArray__common *)C, k, sizeof(double));
  for (k = 0; k < nApnB; k++) {
    C->data[k] = 0.0;
  }

  if ((A->size[1] > 0) && (B->size[1] > 0)) {
    if (B->size[1] > A->size[1]) {
      for (nApnB = 0; nApnB + 1 <= nA; nApnB++) {
        for (k = 0; k < nB; k++) {
          C->data[nApnB + k] += A->data[nApnB] * B->data[k];
        }
      }
    } else {
      for (nApnB = 0; nApnB + 1 <= nB; nApnB++) {
        for (k = 0; k < nA; k++) {
          C->data[nApnB + k] += B->data[nApnB] * A->data[k];
        }
      }
    }
  }
}

/*
 *
 */
void c_conv(const double A[22001], const double B[22001], double C[44001])
{
  int k;
  int b_k;
  memset(&C[0], 0, 44001U * sizeof(double));
  for (k = 0; k < 22001; k++) {
    for (b_k = 0; b_k < 22001; b_k++) {
      C[k + b_k] += B[k] * A[b_k];
    }
  }
}

/*
 *
 */
void conv(const double A[2201], const double B[2201], double C[4401])
{
  int k;
  int b_k;
  memset(&C[0], 0, 4401U * sizeof(double));
  for (k = 0; k < 2201; k++) {
    for (b_k = 0; b_k < 2201; b_k++) {
      C[k + b_k] += B[k] * A[b_k];
    }
  }
}

/*
 *
 */
void d_conv(const double A[221], const double B[221], double C[441])
{
  int k;
  int b_k;
  memset(&C[0], 0, 441U * sizeof(double));
  for (k = 0; k < 221; k++) {
    for (b_k = 0; b_k < 221; b_k++) {
      C[k + b_k] += B[k] * A[b_k];
    }
  }
}

/*
 *
 */
void e_conv(const double A[221], const double B[221], double C[441])
{
  int k;
  int b_k;
  memset(&C[0], 0, 441U * sizeof(double));
  for (k = 0; k < 221; k++) {
    for (b_k = 0; b_k < 221; b_k++) {
      C[k + b_k] += B[k] * A[b_k];
    }
  }
}

/*
 *
 */
void f_conv(const emxArray_real_T *A, const emxArray_real_T *B, emxArray_real_T *
            C)
{
  int nA;
  int nB;
  int nApnB;
  int k;
  nA = A->size[0];
  nB = B->size[0];
  nApnB = A->size[0] + B->size[0];
  if ((A->size[0] == 0) || (B->size[0] == 0)) {
  } else {
    nApnB--;
  }

  k = C->size[0];
  C->size[0] = nApnB;
  emxEnsureCapacity((emxArray__common *)C, k, sizeof(double));
  for (k = 0; k < nApnB; k++) {
    C->data[k] = 0.0;
  }

  if ((A->size[0] > 0) && (B->size[0] > 0)) {
    if (B->size[0] > A->size[0]) {
      for (nApnB = 0; nApnB + 1 <= nA; nApnB++) {
        for (k = 0; k < nB; k++) {
          C->data[nApnB + k] += A->data[nApnB] * B->data[k];
        }
      }
    } else {
      for (nApnB = 0; nApnB + 1 <= nB; nApnB++) {
        for (k = 0; k < nA; k++) {
          C->data[nApnB + k] += B->data[nApnB] * A->data[k];
        }
      }
    }
  }
}

/*
 *
 */
void g_conv(const double A[2201], const double B[2201], double C[4401])
{
  int k;
  int b_k;
  memset(&C[0], 0, 4401U * sizeof(double));
  for (k = 0; k < 2201; k++) {
    for (b_k = 0; b_k < 2201; b_k++) {
      C[k + b_k] += B[k] * A[b_k];
    }
  }
}

/*
 *
 */
void h_conv(const double A[22001], const double B[22001], double C[44001])
{
  int k;
  int b_k;
  memset(&C[0], 0, 44001U * sizeof(double));
  for (k = 0; k < 22001; k++) {
    for (b_k = 0; b_k < 22001; b_k++) {
      C[k + b_k] += B[k] * A[b_k];
    }
  }
}

/* End of code generation (conv.c) */
