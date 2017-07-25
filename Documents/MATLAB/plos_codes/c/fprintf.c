/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * fprintf.c
 *
 * Code generation for function 'fprintf'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "fprintf.h"
#include "fileManager.h"

/* Function Definitions */

/*
 *
 */
int ab_cfprintf(double varargin_1, double varargin_2, double varargin_3, double
                varargin_4, double varargin_5, double varargin_6)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[25] = { ' ', ' ', 'p', '=', '[', '%', 'f', ' ', '%',
    'f', ' ', '%', 'f', ' ', '%', 'f', ' ', '%', 'f', ' ', '%', 'f', ']', '\x0a',
    '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  ab_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3,
                        varargin_4, varargin_5, varargin_6);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int b_cfprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[40] = { 'N', 'o', ' ', 'v', 'a', 'l', 'i', 'd', ' ',
    'm', 'o', 'd', 'e', 'l', ' ', 'w', 'a', 's', ' ', 's', 'e', 'l', 'e', 'c',
    't', 'e', 'd', ',', ' ', 'q', 'u', 'i', 't', 't', 'i', 'n', 'g', '.', '\x0a',
    '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  b_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int bb_cfprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[56] = { 'R', 'e', 'c', 'a', 'l', 'c', 'u', 'l', 'a',
    't', 'i', 'n', 'g', ' ', 'c', 'a', 'n', 'i', 'd', 'a', 't', 'e', ' ', 's',
    'o', 'l', 'u', 't', 'i', 'o', 'n', 's', ' ', 'w', 'i', 't', 'h', ' ', 's',
    'm', 'a', 'l', 'l', 'e', 'r', ' ', 's', 't', 'e', 'p', 's', 'i', 'z', 'e',
    '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  bb_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int c_cfprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[14] = { 'e', 'm', 'g', 'f', 'i', 't', 'n', 'o', 'm',
    'l', 'e', '\x0a', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  c_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int cb_cfprintf(double varargin_1, double varargin_2, double varargin_3, double
                varargin_4, double varargin_5, double varargin_6)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[29] = { 'p', 'd', '_', 'm', 'a', 'x', '=', '[', '%',
    'f', ' ', '%', 'f', ' ', '%', 'f', ' ', '%', 'f', ' ', '%', 'f', ' ', '%',
    'f', ']', '\x0a', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  cb_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3,
                        varargin_4, varargin_5, varargin_6);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int cfprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[13] = { 'F', 'U', 'C', 'C', 'I', ' ', 'D', 'a', 't',
    'a', '\x0a', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int d_cfprintf(double varargin_1)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[6] = { 'i', '=', '%', 'f', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  d_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int e_cfprintf(double varargin_1, double varargin_2, double varargin_3)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[16] = { ' ', ' ', 'P', '=', '[', '%', 'f', ' ', '%',
    'f', ' ', '%', 'f', ']', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  e_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int f_cfprintf(double varargin_1, double varargin_2, double varargin_3)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[17] = { ' ', ' ', 'e', 'p', '=', '[', '%', 'f', ' ',
    '%', 'f', ' ', '%', 'f', ']', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  f_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int g_cfprintf(double varargin_1)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[9] = { ' ', ' ', 'l', '=', '%', 'f', '\x0a', '\x0a',
    '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  g_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int h_cfprintf(double varargin_1, double varargin_2)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[21] = { 'm', 'a', 'x', '_', 'l', 'e', '=', '%', 'f',
    ' ', 'i', 'n', 'd', '_', 'l', 'e', '=', '%', 'f', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  h_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int i_cfprintf(double varargin_1, double varargin_2, double varargin_3)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[20] = { 'e', 'p', '_', 'm', 'a', 'x', '=', '[', '%',
    'f', ' ', '%', 'f', ' ', '%', 'f', ']', '\x0a', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  i_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int j_cfprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[19] = { 'o', 'n', 'e', 's', 't', 'a', 'g', 'e', 'f',
    'i', 't', 'n', 'o', 'm', 'l', 'e', '\x0a', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  j_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int k_cfprintf(double varargin_1, double varargin_2)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[13] = { ' ', ' ', 'P', '=', '[', '%', 'f', ' ', '%',
    'f', ']', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  k_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int l_cfprintf(double varargin_1, double varargin_2)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[13] = { ' ', ' ', 'p', '=', '[', '%', 'f', ' ', '%',
    'f', ']', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  l_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int m_cfprintf(double varargin_1, double varargin_2)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[21] = { 'm', 'a', 'x', '_', 'l', 'd', '=', '%', 'f',
    ' ', 'r', 'o', 'w', '_', 'l', 'd', '=', '%', 'f', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  m_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int n_cfprintf(double varargin_1, double varargin_2)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[17] = { 'p', 'd', '_', 'm', 'a', 'x', '=', '[', '%',
    'f', ' ', '%', 'f', ']', '\x0a', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  n_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int o_cfprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[22] = { 'o', 'n', 'e', 's', 't', 'a', 'g', 'e', 'f',
    'i', 't', 'l', 'a', 'g', 'n', 'o', 'm', 'l', 'e', '\x0a', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  o_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int p_cfprintf(double varargin_1, double varargin_2, double varargin_3)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[16] = { ' ', ' ', 'p', '=', '[', '%', 'f', ' ', '%',
    'f', ' ', '%', 'f', ']', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  p_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int q_cfprintf(double varargin_1, double varargin_2, double varargin_3)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[20] = { 'p', 'd', '_', 'm', 'a', 'x', '=', '[', '%',
    'f', ' ', '%', 'f', ' ', '%', 'f', ']', '\x0a', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  q_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int r_cfprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[18] = { 't', 'w', 'o', 's', 't', 'a', 'g', 'e', 'f',
    'i', 't', 'n', 'o', 'm', 'l', 'e', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  r_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int s_cfprintf(double varargin_1, double varargin_2, double varargin_3, double
               varargin_4)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[19] = { ' ', ' ', 'P', '=', '[', '%', 'f', ' ', '%',
    'f', ' ', '%', 'f', ' ', '%', 'f', ']', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  s_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3,
                        varargin_4);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int t_cfprintf(double varargin_1, double varargin_2, double varargin_3, double
               varargin_4)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[19] = { ' ', ' ', 'p', '=', '[', '%', 'f', ' ', '%',
    'f', ' ', '%', 'f', ' ', '%', 'f', ']', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  t_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3,
                        varargin_4);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int u_cfprintf(double varargin_1, double varargin_2, double varargin_3, double
               varargin_4)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[28] = { ' ', ' ', 'l', '=', '%', 'f', ' ', 'h', 'p',
    '=', '%', 'f', ' ', 'f', 'l', 'a', 'g', '=', '%', 'f', ' ', 'E', '=', '%',
    'f', '\x0a', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  u_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3,
                        varargin_4);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int v_cfprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[56] = { 'r', 'e', 'c', 'a', 'l', 'c', 'u', 'l', 'a',
    't', 'i', 'n', 'g', ' ', 'c', 'a', 'n', 'i', 'd', 'a', 't', 'e', ' ', 's',
    'o', 'l', 'u', 't', 'i', 'o', 'n', 's', ' ', 'w', 'i', 't', 'h', ' ', 's',
    'm', 'a', 'l', 'l', 'e', 'r', ' ', 's', 't', 'e', 'p', 's', 'i', 'z', 'e',
    '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  v_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int w_cfprintf(double varargin_1, double varargin_2, double varargin_3, double
               varargin_4)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[23] = { 'p', 'd', '_', 'm', 'a', 'x', '=', '[', '%',
    'f', ' ', '%', 'f', ' ', '%', 'f', ' ', '%', 'f', ']', '\x0a', '\x0a',
    '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  w_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3,
                        varargin_4);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int x_cfprintf(void)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[20] = { 't', 'h', 'r', 'e', 'e', 's', 't', 'a', 'g',
    'e', 'f', 'i', 't', 'n', 'o', 'm', 'l', 'e', '\x0a', '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  x_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt);
    fflush(filestar);
  }

  return nbytesint;
}

/*
 *
 */
int y_cfprintf(double varargin_1, double varargin_2, double varargin_3, double
               varargin_4, double varargin_5, double varargin_6)
{
  int nbytesint;
  FILE * b_NULL;
  FILE * filestar;
  boolean_T autoflush;
  static const char cfmt[25] = { ' ', ' ', 'P', '=', '[', '%', 'f', ' ', '%',
    'f', ' ', '%', 'f', ' ', '%', 'f', ' ', '%', 'f', ' ', '%', 'f', ']', '\x0a',
    '\x00' };

  b_NULL = NULL;
  nbytesint = 0;
  y_fileManager(&filestar, &autoflush);
  if (!(filestar == b_NULL)) {
    nbytesint = fprintf(filestar, cfmt, varargin_1, varargin_2, varargin_3,
                        varargin_4, varargin_5, varargin_6);
    fflush(filestar);
  }

  return nbytesint;
}

/* End of code generation (fprintf.c) */
