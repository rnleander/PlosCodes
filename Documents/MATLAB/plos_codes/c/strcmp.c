/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * strcmp.c
 *
 * Code generation for function 'strcmp'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "strcmp.h"

/* Function Definitions */

/*
 *
 */
#ifdef _OLD_MATLAB_CODE
boolean_T b_strcmp(const char a[3])
{
  boolean_T b_bool;
  int kstr;
  int exitg1;
  static const char cv0[3] = { 'a', 'l', 'l' };

  b_bool = false;
  kstr = 0;
  do {
    exitg1 = 0;
    if (kstr + 1 < 4) {
      if (a[kstr] != cv0[kstr]) {
        exitg1 = 1;
      } else {
        kstr++;
      }
    } else {
      b_bool = true;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return b_bool;
}
#endif

/* End of code generation (strcmp.c) */
