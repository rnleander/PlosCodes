/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * tailmass.h
 *
 * Code generation for function 'tailmass'
 *
 */

#ifndef TAILMASS_H
#define TAILMASS_H

/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_defines.h"
#include "rt_nonfinite.h"
#include "rtwtypes.h"
#include "omp.h"
#include "IMT_analysis_April2017_types.h"

/* Function Declarations */
extern double tailmass(const double m[2], const double s[2], double T2, const double sd[2]);

extern int checktailmass(const double m_a, const double s_a, const double m_b, const double s_b, double T2, const double sd_a, const double sd_b);

#endif

/* End of code generation (tailmass.h) */
