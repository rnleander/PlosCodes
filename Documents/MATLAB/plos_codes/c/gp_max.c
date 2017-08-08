/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * gp_max.c
 *
 * Code generation for function 'gp_max'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "gp_max.h"
#include "onestagepdf_prime.h"
#include "emgpdf.h"
#include "IMT_analysis_April2017_rtwutil.h"
#include <gsl/gsl_poly.h>

/* Function Definitions */

/*
 * function M=gp_max(m,s)
 */

double gp_max_fixed(double m, double s)
{
	double coeff[5];
	coeff[4] = rt_powd_snf(m, 4.0) / (4.0 * rt_powd_snf(s, 4.0));
	coeff[3] = 3.0 * (m * m) / (2.0 * (s * s));
	coeff[2] = 3.75 - m * m / (2.0 * rt_powd_snf(s, 4.0));
	coeff[1] = -5.0 / (2.0 * (s * s));
	coeff[0] = 1.0 / (4.0 * rt_powd_snf(s, 4.0));

	double foundRoots[10];

	gsl_poly_complex_workspace * workspace = gsl_poly_complex_workspace_alloc(5);
	gsl_poly_complex_solve(coeff, 5, workspace, foundRoots);
	gsl_poly_complex_workspace_free(workspace);

	int numRealRoots = 0;
	double realRoots[10];
	double re;
	double im;
	for (int i = 0; i < 4; i++)
	{
		re = foundRoots[2 * i];
		im = foundRoots[2 * i + 1];

		if (im <= 0.000000000000001) {
			if (re >= 0) {
				realRoots[numRealRoots] = re;
				numRealRoots++;
			}
		}
	}
	double y[10];
	onestagepdf_prime_fixed(realRoots, numRealRoots, m, s, y);

	double largest = 0;
	for (int i = 0; i < numRealRoots; i++) {
		if (y[i] > largest) {
			largest = y[i];
		}
	}

	return largest;
}

/* End of code generation (gp_max.c) */
