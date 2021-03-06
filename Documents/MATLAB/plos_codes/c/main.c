/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * main.c
 *
 * Code generation for function 'main'
 *
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/
/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "main.h"
#include "IMT_analysis_April2017_terminate.h"
#include "IMT_analysis_April2017_initialize.h"


#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_poly.h>





/* Function Declarations */
static void main_IMT_analysis_April2017(const char *);

/* Function Definitions */
static void main_IMT_analysis_April2017(const char *model)
{
  /* Call the entry-point 'IMT_analysis_April2017'. */
  IMT_analysis_April2017(model);
}

int main(int argc, const char * const argv[])
{
  (void)argc;
  (void)argv;

  

	  

  /* Initialize the application.
     You do not need to do this more than one time. */
  IMT_analysis_April2017_initialize();

  /* Invoke the entry-point functions.
     You can call entry-point functions multiple times. */
  if (argc == 2) {
	  main_IMT_analysis_April2017(argv[1]);
  }
  else {
	  printf("invalid number of arguments, please specify model\n");
  }

  IMT_analysis_April2017_terminate();
  return 0;
}

/* End of code generation (main.c) */
