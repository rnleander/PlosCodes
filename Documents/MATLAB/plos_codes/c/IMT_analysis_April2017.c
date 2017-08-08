/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * IMT_analysis_April2017.c
 *
 * Code generation for function 'IMT_analysis_April2017'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "sum.h"
#include "log.h"
#include "emgpdf.h"
#include "fminsearch.h"
#include "onestagepdf2.h"
#include "onestagepdf_lag.h"
#include "convolv_2invG_adapt_nov.h"
#include "convolv_3invG_nov.h"
#include "strcmp.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_statistics_double.h"

#define _GSL_MINIMIZE
#define _VERBOSE

/* Type Definitions */
#ifndef typedef_cell_wrap_3
#define typedef_cell_wrap_3

typedef struct {
  double f1[2];
} cell_wrap_3;

#endif                                 /*typedef_cell_wrap_3*/

/* Function Definitions */

/*
 * function IMT_analysis_April2017(model)
 */
void IMT_analysis_April2017(const char *model)
{

#ifdef _GSL_MINIMIZE
	printf("Using GSL minimizer\n");
#else
	printf("Using Matlab's fminsearch\n");
#endif

  int twostagefitnomle;
  int onestagelagnomle;
  int onestagefitnomle;
  int emgfitnomle;
  int threestagefitnomle;
  int itmp;
  double P[81];



  double mtmp;

  int ix;
  boolean_T exitg1;
  double p[3];
  double ep[81];

  /*
  static const double dv27[266] = { 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
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
	*/


  double l[266];
  double b_P[32];

  /*
  static const double dv28[2] = { 0.04246759052302191, 0.11796527321068413 };
  static const double dv29[2] = { 0.08493518104604382, 0.029491318302671033 };
  static const double dv30[2] = { 0.12740277156906574, 0.029491318302671033 };
  static const double dv31[2] = { 0.12740277156906574, 0.058982636605342066 };
  static const double dv32[2] = { 0.16987036209208764, 0.029491318302671033 };
  static const double dv33[2] = { 0.16987036209208764, 0.058982636605342066 };
  */

  double ld[16];


  double b_p[2];

  double pd[32];
  static const double dv34[3] = { 0.027593515034192811, 0.043688113214883209, 0.67839610301742326 };
  static const double dv35[3] = { 0.027593515034192811, 0.043688113214883209, 1.3567922060348465 };
  static const double dv36[3] = { 0.027593515034192811, 0.043688113214883209, 2.0351883090522698 };
  static const double dv37[3] = { 0.027593515034192811, 0.087376226429766418, 0.67839610301742326 };
  static const double dv38[3] = { 0.027593515034192811, 0.087376226429766418, 1.3567922060348465 };
  static const double dv39[3] = { 0.027593515034192811, 0.087376226429766418, 2.0351883090522698 };
  static const double dv40[3] = { 0.027593515034192811, 0.17475245285953284, 0.67839610301742326 };
  static const double dv41[3] = { 0.027593515034192811, 0.17475245285953284, 1.3567922060348465 };
  static const double dv42[3] = { 0.027593515034192811, 0.17475245285953284, 2.0351883090522698 };
  static const double dv43[3] = { 0.055187030068385622, 0.043688113214883209, 0.67839610301742326 };
  static const double dv44[3] = { 0.055187030068385622, 0.043688113214883209, 1.3567922060348465 };
  static const double dv45[3] = { 0.055187030068385622, 0.043688113214883209, 2.0351883090522698 };
  static const double dv46[3] = { 0.055187030068385622, 0.087376226429766418, 0.67839610301742326 };
  static const double dv47[3] = { 0.055187030068385622, 0.087376226429766418, 1.3567922060348465 };
  static const double dv48[3] = { 0.055187030068385622, 0.087376226429766418, 2.0351883090522698 };
  static const double dv49[3] = { 0.055187030068385622, 0.17475245285953284, 0.67839610301742326 };
  static const double dv50[3] = { 0.055187030068385622, 0.17475245285953284, 1.3567922060348465 };
  static const double dv51[3] = { 0.055187030068385622, 0.17475245285953284, 2.0351883090522698 };
  static const double dv52[3] = { 0.082780545102578429, 0.043688113214883209, 0.67839610301742326 };
  static const double dv53[3] = { 0.082780545102578429, 0.043688113214883209, 1.3567922060348465 };
  static const double dv54[3] = { 0.082780545102578429, 0.043688113214883209, 2.0351883090522698 };
  static const double dv55[3] = { 0.082780545102578429, 0.087376226429766418, 0.67839610301742326 };
  static const double dv56[3] = { 0.082780545102578429, 0.087376226429766418, 1.3567922060348465 };
  static const double dv57[3] = { 0.082780545102578429, 0.087376226429766418, 2.0351883090522698 };
  static const double dv58[3] = { 0.082780545102578429, 0.17475245285953284, 0.67839610301742326 };
  static const double dv59[3] = { 0.082780545102578429, 0.17475245285953284, 1.3567922060348465 };
  static const double dv60[3] = { 0.082780545102578429, 0.17475245285953284, 2.0351883090522698 };

  cell_wrap_3 pcell[9];
  double c_P[180];
  static const signed char id[90] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
    2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 6, 6, 7,
    7, 7, 8, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 2, 3, 4, 5, 6, 7, 8, 9, 3, 4, 5, 6,
    7, 8, 9, 4, 5, 6, 7, 8, 9, 5, 6, 7, 8, 9, 6, 7, 8, 9, 7, 8, 9, 8, 9, 9 };

  double b_pd[180];
  double flag;
  double E;
  double b_flag[45];
  double c_p[4];
  cell_wrap_3 b_pcell[4];
  static const double dv61[2] = { 0.8494, 0.0843 };

  double d_P[120];
  static const signed char b_id[60] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2,
    2, 2, 3, 3, 3, 4, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 2, 2, 2, 3, 3, 4, 3, 3, 4, 4,
    1, 2, 3, 4, 2, 3, 4, 3, 4, 4, 2, 3, 4, 3, 4, 4, 3, 4, 4, 4 };

  double c_pd[120];
  double c_flag[20];
  double d_p[6];


  static const double data[266] = { 11.9, 10.6, 11.6, 9.8, 9.3, 9.0, 11.6, 11.1,
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





  /*  This file fits EMG, one-, two-, and three-stage stochasttic models to IMT data. */
  /* clear all */
  /* This section of code gets the IMT data */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* for erlotinib */
  /* load('erlot_imts_April2017.mat') */
  /* data=imt_b; */
  /* for AT1  */
  /* load('AT1_imts_April2017.mat') */
  /* data=imt_b; */
  /* for MCF */
  /* load('MCF_imts_April2017.mat') */
  /* data=imt_b; */
  /* for DMSO  */
  /* load('DMSO_imts_April2017.mat') */
  /* data=imt_b; */
  /* for CHX */
  /* load('CHX_imts_April2017.mat') */
  /* data=imt_b; */
  /* for FUCCI data */
  /* 'IMT_analysis_April2017:31' fprintf('FUCCI Data\n\n'); */
  //cfprintf();
  printf("FUCCI Data\n\n");

  /* GetProcessedDataParts */
  /* old_enviroment = load('FUCCI_April2017.mat'); */
  /* data=imt_b; */
  /* data=old_enviroment.G2Time_b; */
  /* data = csvread('fucci.csv'); */
  /* 'IMT_analysis_April2017:37' data = [11.9000;10.6000;11.6000;9.8000;9.3000;9.0000;11.6000;11.1000;12.4000;13.7000;12.4000;11.9000;10.3000;12.9000;14.7000;11.6000;13.4000;13.4000;11.6000;10.3000;9.3000;13.7000;9.6000;10.1000;9.8000;10.9000;16.0000;9.3000;9.6000;10.3000;11.4000;10.6000;8.5000;10.3000;11.1000;8.0000;10.6000;7.5000;12.9000;9.0000;8.5000;12.4000;11.6000;9.6000;9.6000;14.7000;9.8000;10.3000;12.1000;8.8000;10.6000;12.1000;13.4000;12.4000;8.8000;13.2000;10.1000;11.6000;11.1000;15.8000;12.1000;12.7000;12.7000;11.1000;13.2000;11.9000;12.4000;13.2000;14.0000;8.0000;8.8000;9.3000;16.5000;14.5000;10.1000;14.2000;7.8000;13.2000;8.8000;8.8000;10.1000;11.9000;12.9000;14.5000;10.9000;10.6000;14.0000;8.8000;8.8000;9.0000;10.9000;14.5000;9.6000;12.4000;11.9000;12.4000;11.1000;14.5000;10.3000;12.4000;12.7000;11.9000;10.3000;13.7000;15.5000;14.5000;11.6000;10.6000;15.5000;14.7000;8.8000;11.6000;8.3000;17.6000;12.4000;11.6000;15.0000;13.7000;12.7000;10.9000;7.2000;8.5000;8.3000;9.6000;11.4000;12.9000;11.6000;13.4000;10.1000;11.6000;8.8000;12.4000;10.3000;16.3000;10.9000;10.1000;8.8000;9.3000;15.2000;8.5000;11.1000;8.3000;11.4000;11.9000;9.3000;9.8000;16.3000;12.7000;9.0000;11.9000;9.3000;10.3000;13.4000;11.4000;12.9000;12.4000;9.6000;10.3000;13.2000;10.6000;9.8000;11.9000;14.2000;13.4000;9.3000;9.6000;12.1000;11.9000;10.1000;14.0000;12.9000;21.7000;11.6000;12.1000;10.3000;9.8000;14.2000;13.7000;7.2000;10.9000;10.1000;9.6000;13.4000;13.2000;16.3000;11.6000;14.0000;10.9000;14.2000;12.4000;12.4000;13.4000;17.6000;10.1000;10.9000;14.0000;12.9000;9.0000;13.4000;15.0000;16.0000;8.0000;9.8000;12.4000;8.5000;9.6000;12.7000;12.1000;15.0000;16.0000;10.9000;14.2000;13.7000;11.9000;16.8000;11.4000;13.4000;12.4000;22.0000;12.4000;16.8000;12.1000;10.3000;13.4000;11.6000;10.1000;14.5000;10.6000;11.9000;15.5000;9.8000;12.4000;10.1000;8.0000;9.0000;9.3000;13.2000;11.1000;12.7000;12.1000;10.1000;13.2000;14.5000;10.1000;12.7000;12.9000;11.9000;12.4000;11.1000;8.5000;14.5000;16.5000;12.4000;9.0000;11.1000;9.8000;11.1000;11.1000;8.8000;13.2000;17.6000;16.8000;10.9000;12.4000;8.5000;14.7000]; */
  /* data=G1Time_b; */
  /* for PC9 cells */
  /* load('PC9_April2017.mat') */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */


  /* choose model to fit */
  twostagefitnomle = 0;
  const char * twostagefitnomle_str = "twostage";
  onestagelagnomle = 0;
  const char * onestagelagnomle_str = "onestagelag";
  onestagefitnomle = 0;
  const char * onestagefitnomle_str = "onestage";
  emgfitnomle = 0;
  const char * emgfitnomle_str = "emg";
  threestagefitnomle = 0;
  const char * threestagefitnomle_str = "threestage";
  const char * all_str = "all";

  if (strcmp(model, twostagefitnomle_str) == 0) {
    twostagefitnomle = 1;
  } else if (strcmp(model, onestagelagnomle_str) == 0) {
    onestagelagnomle = 1;
  } else if (strcmp(model, onestagefitnomle_str) == 0) {
    onestagefitnomle = 1;
  } else if (strcmp(model, emgfitnomle_str) == 0) {
    emgfitnomle = 1;
  } else if (strcmp(model, threestagefitnomle_str) == 0) {
    threestagefitnomle = 1;
  } else if (strcmp(model, all_str) == 0) {
    emgfitnomle = 1;
    onestagefitnomle = 1;
    onestagelagnomle = 1;
    twostagefitnomle = 1;
    threestagefitnomle = 1;
  } else {
    /* 'IMT_analysis_April2017:85' fprintf('No valid model was selected, quitting.\n'); */
    //b_cfprintf();
	  printf("No valid model was selected, quitting.\n");
  }

  /* get sample statistics for fitting initializing the model parameters */
  /* 'IMT_analysis_April2017:97' num = length(data); */
  /* 'IMT_analysis_April2017:98' C1 = mean(data); */
  /* 'IMT_analysis_April2017:99' C2 = var(data); */
  /* 'IMT_analysis_April2017:100' C3 = sum((data-C1).^3)/(length(data)); */
  /*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */


























  /* Fit the EMG model */
  if (emgfitnomle == 1) {

	printf("emgfitnomle\n\n");

	/* prepare statistical variables */
	double mean = gsl_stats_mean(data, 1, 266);
	double variance = gsl_stats_variance(data, 1, 266);
	double vry[3] = { 0.25, 0.5, 0.75 };

	/*we vary the parameters so the the Gaussian and exponential parts of
	the cell cycle are responsible for a fraction of the total mean and
	variance in the IMT.
	*/
	double lambda[3];
	for (int i = 0; i < 3; i++) 
		lambda[i] = 1.0 / (mean*vry[i]);

	double mu[3];
	for (int i = 0; i < 3; i++) 
		mu[i] = mean * vry[i];

	double sigma[3];
	for (int i = 0; i < 3; i++) 
		sigma[i] = pow(variance*vry[i], 0.5);

	/* prepare parameter seeds */
	int seedIdx = 0;
	double paramSeeds[27][3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				paramSeeds[seedIdx][0] = lambda[i];
				paramSeeds[seedIdx][1] = sigma[j];
				paramSeeds[seedIdx][2] = mu[k];
				seedIdx++;
			}
		}
	}

	/* optimize parameters */
	double optimizedParams[27][3];
	double le[27];
    for (seedIdx = 0; seedIdx < 27; seedIdx++) {

	  printf("i = %d\n", 1 + seedIdx);
	  printf("  P=[%.17f %.17f %.17f]\n", paramSeeds[seedIdx][0], paramSeeds[seedIdx][1], paramSeeds[seedIdx][2]);

	  /* select our working seed */
	  double paramSeed[3];
      for (int i = 0; i < 3; i++) 
		  paramSeed[i] = paramSeeds[seedIdx][i];

#ifdef _GSL_MINIMIZE
	  // https://www.gnu.org/software/gsl/doc/html/multimin.html#algorithms-without-derivatives

	  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	  gsl_multimin_fminimizer *s = NULL;
	  gsl_vector *ss, *x;
	  gsl_multimin_function minex_func;

	  size_t iter = 0;
	  int status;
	  double size;

	  /* Starting point */
	  x = gsl_vector_alloc(3);
	  gsl_vector_set(x, 0, paramSeed[0]);
	  gsl_vector_set(x, 1, paramSeed[1]);
	  gsl_vector_set(x, 2, paramSeed[2]);

	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc(3);
	  gsl_vector_set_all(ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 3;
	  minex_func.f = emgpdf_loglikelihood;
	  minex_func.params = (void *)data;

	  s = gsl_multimin_fminimizer_alloc(T, 3);
	  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	  do
	  {
		  iter++;
		  status = gsl_multimin_fminimizer_iterate(s);

		  if (status)
			  break;

		  size = gsl_multimin_fminimizer_size(s);
		  status = gsl_multimin_test_size(size, 1e-4);

#ifdef _VERBOSE
		  if (status == GSL_SUCCESS)
		  {
			  printf("converged to minimum at\n");
		  }

		  printf("%5d %.17f %.17f %.17f f() = %.17f size = %.3f\n",
			  iter,
			  gsl_vector_get(s->x, 0),
			  gsl_vector_get(s->x, 1),
			  gsl_vector_get(s->x, 2),
			  s->fval, size);
#endif
	  } while (status == GSL_CONTINUE && iter < 10000);

	  optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
	  optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));
	  optimizedParams[seedIdx][2] = fabs(gsl_vector_get(s->x, 2));

	  gsl_vector_free(x);
	  gsl_vector_free(ss);
	  gsl_multimin_fminimizer_free(s);

#else
	  void(*emgpdf_ptr)(const double X[266], double l, double m, double s, double Y[266]) = &emgpdf;
	  fminsearch_generalized(paramSeed, emgpdf_ptr, data);
	  for (int i = 0; i < 3; i++) 
		  optimizedParams[seedIdx][i] = paramSeed[i];
#endif

	  printf("  ep=[%.17f %.17f %.17f]\n", optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2]);

	  /* recalculate our best fit */
      emgpdf(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2], l);

	  /* calculate the log likelihood for our best fit */
	  double loglikelihood = 0;
	  for (int i = 0; i < 266; i++) {
		  loglikelihood += log(l[i]);
	  }

	  printf("  l=%.17f\n\n", loglikelihood);
      le[seedIdx] = loglikelihood;
    }

	/* find the best optimized parameter set for all starting seeds tried */
	double maxLikelihood = 0;
	int ind_le = 0;
	for (int i = 0; i < 27; i++) {
		if (le[i] < maxLikelihood){
			maxLikelihood = le[i];
			ind_le = i;
		}
	}

	printf("max_le=%.17f ind_le=%d\n", maxLikelihood, ind_le+1);
	printf("ep_max=[%.17f %.17f %.17f]\n\n", optimizedParams[ind_le][0], optimizedParams[ind_le][1], optimizedParams[ind_le][2]);
  }

 




























  /* Fit the one-stage model */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* 'IMT_analysis_April2017:283' if onestagefitnomle == 1 */
  if (onestagefitnomle == 1) {

	  printf("onestagefitnomle\n\n");

	  /* prepare statistical variables */
	  double mean = gsl_stats_mean(data, 1, 266);
	  double variance = gsl_stats_variance(data, 1, 266);
	  double vry[4] = { 0.5, 1.0, 1.5, 2.0 };

	  double mu = 1.0 / mean;
	  double sigma = pow((variance / pow(mean, 3)), 0.5);

	  double m[4];
	  for (int i = 0; i < 4; i++) 
		  m[i] = mu*vry[i];
	  
	  double s[4];
	  for (int i = 0; i < 4; i++)
		  s[i] = sigma*vry[i];

	  /* prepare parameter seeds */
	  int seedIdx = 0;
	  double paramSeeds[16][2];
	  for (int i = 0; i < 4; i++) {
		  for (int j = 0; j < 4; j++) {
				  paramSeeds[seedIdx][0] = m[i];
				  paramSeeds[seedIdx][1] = s[j];
				  seedIdx++;
		  }
	  }

	/* optimize parameters */
	  double optimizedParams[27][3];
	  double le[27];
    for (seedIdx = 0; seedIdx < 16; seedIdx++) {
		printf("i=%d\n", 1.0 + seedIdx);
		printf("  P=[%f %f]\n", paramSeeds[seedIdx][0], paramSeeds[seedIdx][1]);

		

#ifdef _GSL_MINIMIZE
	  // https://www.gnu.org/software/gsl/doc/html/multimin.html#algorithms-without-derivatives

	  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	  gsl_multimin_fminimizer *s = NULL;
	  gsl_vector *ss, *x;
	  gsl_multimin_function minex_func;

	  size_t iter = 0;
	  int status;
	  double size;

	  /* Starting point */
	  x = gsl_vector_alloc(2);
	  gsl_vector_set(x, 0, paramSeeds[seedIdx][0]);
	  gsl_vector_set(x, 1, paramSeeds[seedIdx][1]);

	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc(2);
	  gsl_vector_set_all(ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 2;
	  minex_func.f = wald_loglikelihood;
	  minex_func.params = (void *)data;

	  s = gsl_multimin_fminimizer_alloc(T, 2);
	  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	  do
	  {
		  iter++;
		  status = gsl_multimin_fminimizer_iterate(s);

		  if (status)
			  break;

		  size = gsl_multimin_fminimizer_size(s);
		  status = gsl_multimin_test_size(size, 1e-4);

#ifdef _VERBOSE
		  if (status == GSL_SUCCESS)
		  {
			  printf("converged to minimum at\n");
		  }

		  printf("%5d %.17f %.17f f() = %.17f size = %.3f\n",
			  iter,
			  gsl_vector_get(s->x, 0),
			  gsl_vector_get(s->x, 1),
			  s->fval, size);
#endif
	  } while (status == GSL_CONTINUE && iter < 10000);

	  optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
	  optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));

	  gsl_vector_free(x);
	  gsl_vector_free(ss);
	  gsl_multimin_fminimizer_free(s);
#else
		for (int i = 0; i < 2; i++) {
			b_p[i] = paramSeeds[seedIdx][i];
		}
      b_fminsearch(b_p);
	  optimizedParams[seedIdx][0] = b_p[0];
	  optimizedParams[seedIdx][1] = b_p[1];
#endif


	  printf("  p=[%f %f]\n", optimizedParams[seedIdx][0], optimizedParams[seedIdx][1]);

      //onestagepdf2(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], l);
	  waldpdf(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], l, 266);


	  /* calculate the log likelihood for our best fit */
	  double loglikelihood = 0;
	  for (int i = 0; i < 266; i++) {
		  loglikelihood += log(l[i]);
	  }

	  printf("  l=%.17f\n\n", loglikelihood);
	  ld[seedIdx] = loglikelihood;
    }

	/* find the best optimized parameter set for all starting seeds tried */
	double maxLikelihood = 0;
	int ind_ld = 0;
	for (int i = 0; i < 27; i++) {
		if (ld[i] < maxLikelihood) {
			maxLikelihood = ld[i];
			ind_ld = i;
		}
	}

	printf("max_ld=%f row_ld=%f\n", maxLikelihood, ind_ld);
	printf("pd_max=[%f %f]\n\n", optimizedParams[ind_ld][0], optimizedParams[ind_ld][1]);
  }

  

































  /* Fit one-stage model with lag with fminsearch */
  if (onestagelagnomle == 1) {

	  printf("onestagefitlagnomle\n\n");

	  /* prepare statistical variables */
	  double mean = gsl_stats_mean(data, 1, 266);
	  double variance = gsl_stats_variance(data, 1, 266);

	  // C3 = sum((data-C1).^3)/(length(data));
	  double c3sum = 0;
	  for (int i = 0; i < 266; i++) {
		  c3sum += pow(data[i] - mean, 3.0);
	  }
	  c3sum = c3sum / 266;

	  double mu = c3sum / (3 * pow(variance, 2.0));
	  double sigma = pow((pow(c3sum, 3.0) / (27*pow(variance, 5.0))), 0.5);
	  double lag = mean - 3 * pow(variance, 2.0) / c3sum;

	  double vryv[3] = { 0.5, 1, 2 };
	  double vrym[3] = { 0.25, 0.5, 0.75 };

	  double m[3];
	  for (int i = 0; i < 3; i++)
		  m[i] = mu*vrym[i];

	  double s[3];
	  for (int i = 0; i < 3; i++)
		  s[i] = sigma*vryv[i];

	  double L[3];
	  for (int i = 0; i < 3; i++)
		  L[i] = lag*vrym[i];


	  /* prepare parameter seeds */
	  int seedIdx = 0;
	  double paramSeeds[64][3];
	  for (int i = 0; i < 3; i++) {
		  for (int j = 0; j < 3; j++) {
			  for (int k = 0; k < 3; k++) {
				  paramSeeds[seedIdx][0] = m[i];
				  paramSeeds[seedIdx][1] = s[j];
				  paramSeeds[seedIdx][2] = L[k];
				  seedIdx++;
			  }
		  }
	  }

	/* optimize parameters */
	double optimizedParams[27][3];
	double le[27];

    for (seedIdx = 0; seedIdx < 27; seedIdx++) {

		printf("i=%d\n", 1 + seedIdx);
		printf("  P=[%f %f %f]\n", paramSeeds[seedIdx][0], paramSeeds[seedIdx][1], paramSeeds[seedIdx][2]);


      /* 'IMT_analysis_April2017:434' x0=P(i,:); */
      /* [p,conf1]=mle(data,'pdf',@onestagepdf_lag,'start',x0, 'upperbound', [Inf Inf Inf],'lowerbound',[0 0 0],'options',options); */
      /* 'IMT_analysis_April2017:436' f=@(x)(sum(-log(onestagepdf_lag(data,x(1),x(2),x(3))))); */
      /* 'IMT_analysis_April2017:437' [p,pval]=fminsearch(f, x0, options); */

#ifdef _GSL_MINIMIZE
	  // https://www.gnu.org/software/gsl/doc/html/multimin.html#algorithms-without-derivatives

	  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	  gsl_multimin_fminimizer *s = NULL;
	  gsl_vector *ss, *x;
	  gsl_multimin_function minex_func;

	  size_t iter = 0;
	  int status;
	  double size;

	  /* Starting point */
	  x = gsl_vector_alloc(3);
	  gsl_vector_set(x, 0, paramSeeds[seedIdx][0]);
	  gsl_vector_set(x, 1, paramSeeds[seedIdx][1]);
	  gsl_vector_set(x, 2, paramSeeds[seedIdx][2]);

	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc(3);
	  gsl_vector_set_all(ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 3;
	  minex_func.f = waldlag_loglikelihood;
	  minex_func.params = (void *)data;

	  s = gsl_multimin_fminimizer_alloc(T, 3);
	  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	  do
	  {
		  iter++;
		  status = gsl_multimin_fminimizer_iterate(s);

		  if (status)
			  break;

		  size = gsl_multimin_fminimizer_size(s);
		  status = gsl_multimin_test_size(size, 1e-4);

#ifdef _VERBOSE
		  if (status == GSL_SUCCESS)
		  {
			  printf("converged to minimum at\n");
		  }

		  printf("%5d %.17f %.17f %.17f f() = %.17f size = %.3f\n",
			  iter,
			  gsl_vector_get(s->x, 0),
			  gsl_vector_get(s->x, 1),
			  gsl_vector_get(s->x, 2),
			  s->fval, size);
#endif
	  } while (status == GSL_CONTINUE && iter < 10000);

	  optimizedParams[seedIdx][0] = fabs(gsl_vector_get(s->x, 0));
	  optimizedParams[seedIdx][1] = fabs(gsl_vector_get(s->x, 1));
	  optimizedParams[seedIdx][2] = fabs(gsl_vector_get(s->x, 2));

	  gsl_vector_free(x);
	  gsl_vector_free(ss);
	  gsl_multimin_fminimizer_free(s);

#else
      for (int i = 0; i < 3; i++) {
			p[i] = paramSeeds[seedIdx][i];
	  }
      c_fminsearch(p);
	  for (int i = 0; i < 3; i++) {
		  optimizedParams[seedIdx][i] = p[i];
	  }

#endif

	  printf("  p=[%f %f %f]\n", optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2]);

      //onestagepdf_lag(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2], l);
	  waldlagpdf(data, optimizedParams[seedIdx][0], optimizedParams[seedIdx][1], optimizedParams[seedIdx][2], l);


	  /* calculate the log likelihood for our best fit */
	  double loglikelihood = 0;
	  for (int i = 0; i < 266; i++) {
		  loglikelihood += log(l[i]);
	  }

	  printf("  l=%.17f\n\n", loglikelihood);
	  le[seedIdx] = loglikelihood;
    }

	/* find the best optimized parameter set for all starting seeds tried */
	double maxLikelihood = 0;
	int ind_ld = 0;
	for (int i = 0; i < 27; i++) {
		if (le[i] < maxLikelihood) {
			maxLikelihood = le[i];
			ind_ld = i;
		}
	}

	printf("max_ld=%f row_ld=%d\n", maxLikelihood, ind_ld + 1);
	printf("pd_max=[%f %f %f]\n\n", optimizedParams[ind_ld][0], optimizedParams[ind_ld][1], optimizedParams[ind_ld][2]);
  }
































  
  /*  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* Fit two-stage model without mle */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* 'IMT_analysis_April2017:547' if twostagefitnomle == 1 */
  if (twostagefitnomle == 1) {
    /*  BEGIN FUNCTION FIT_TWOSTAGE */
    /* 'IMT_analysis_April2017:549' fprintf('twostagefitnomle\n'); */
    //r_cfprintf();
	  printf("twostagefitnomle\n");

	  /* prepare statistical variables */
	  double mean = gsl_stats_mean(data, 1, 266);
	  double variance = gsl_stats_variance(data, 1, 266);

	  double vry[3] = { 0.25, 0.5, 0.75 };
	  double vrys[3] = { 0.01, 1, 10 };

	  double m[3];
	  double s[3];
	  for (int i = 0; i < 3; i++) {
		  m[i] = 1 / (mean*vry[i]);
		  s[i] = pow((variance*vry[i])/(pow(mean*vry[i],3.0)), 0.5);
	  }

	  double pcomb[9][2];
	  int pcomb_idx = 0;
	  for (int i = 0; i < 3; i++) {
		  for (int j = 0; j < 3; j++) {
			  pcomb[pcomb_idx][0] = m[i];
			  pcomb[pcomb_idx][1] = s[j];
		  }
	  }

    /*  prepare parameter seeds */
    /* get all pairs of the form [m(i),s(j)] */
    /* these pairs represent all possible unique  */
    /* parameter choices for each part of the cell */
    /* cycle.   */
    /* pcomb = allcomb(m,s) */
    /* 'IMT_analysis_April2017:572' pcomb = [0.3397, 0.2359; 0.3397, 0.1180; 0.3397, 0.0786; 0.1699, 0.2359; 0.1699, 0.1180; 0.1699, 0.0786; 0.1132, 0.2359; 0.1132, 0.1180; 0.1132, 0.0786]; */
    /* place paramter pairs into a cell.  The parameters choices for each part */
    /* are now indexed */
    /* 'IMT_analysis_April2017:577' pcell = cell(length(pcomb),1); */
    /* 'IMT_analysis_April2017:578' for i = 1:length(pcomb) */
    /* 'IMT_analysis_April2017:579' pcell{i} = pcomb(i,:); */
    /* 'IMT_analysis_April2017:579' pcell{i} = pcomb(i,:); */
    /* 'IMT_analysis_April2017:579' pcell{i} = pcomb(i,:); */
    /* 'IMT_analysis_April2017:579' pcell{i} = pcomb(i,:); */
    /* 'IMT_analysis_April2017:579' pcell{i} = pcomb(i,:); */
    /* 'IMT_analysis_April2017:579' pcell{i} = pcomb(i,:); */
    /* 'IMT_analysis_April2017:579' pcell{i} = pcomb(i,:); */
    /* 'IMT_analysis_April2017:579' pcell{i} = pcomb(i,:); */
    /* 'IMT_analysis_April2017:579' pcell{i} = pcomb(i,:); */
    for (itmp = 0; itmp < 2; itmp++) {
      pcell[0].f1[itmp] = 0.3397 + -0.1038 * (double)itmp;
      pcell[1].f1[itmp] = 0.3397 + -0.2217 * (double)itmp;
      pcell[2].f1[itmp] = 0.3397 + -0.2611 * (double)itmp;
      pcell[3].f1[itmp] = 0.1699 + 0.066 * (double)itmp;
      pcell[4].f1[itmp] = 0.1699 + -0.0519 * (double)itmp;
      pcell[5].f1[itmp] = 0.1699 + -0.091299999999999992 * (double)itmp;
      pcell[6].f1[itmp] = 0.1132 + 0.1227 * (double)itmp;
      pcell[7].f1[itmp] = 0.1132 + 0.0047999999999999987 * (double)itmp;
      pcell[8].f1[itmp] = 0.1132 + -0.034599999999999992 * (double)itmp;
    }

    /* get all pairs of indices for the parameter  */
    /* choices for each part of the cycle to get all  */
    /* parameter choices for the entire cycle */
    /* id = allcomb(1:length(pcomb),1:length(pcomb)) */
    /* id = [1,1;1,2;1,3;1,4;1,5;1,6;1,7;1,8;1,9;2,1;2,2;2,3;2,4;2,5;2,6;2,7;2,8;2,9;3,1;3,2;3,3;3,4;3,5;3,6;3,7;3,8;3,9;4,1;4,2;4,3;4,4;4,5;4,6;4,7;4,8;4,9;5,1;5,2;5,3;5,4;5,5;5,6;5,7;5,8;5,9;6,1;6,2;6,3;6,4;6,5;6,6;6,7;6,8;6,9;7,1;7,2;7,3;7,4;7,5;7,6;7,7;7,8;7,9;8,1;8,2;8,3;8,4;8,5;8,6;8,7;8,8;8,9;9,1;9,2;9,3;9,4;9,5;9,6;9,7;9,8;9,9]; */
    /* sort the pairs in ascending order.   */
    /* This equates choices of the form [i,j] and [j,i]. */
    /* id = sort(id,2); */
    /* remove repeats */
    /* id = unique(id,'rows') */
    /* 'IMT_analysis_April2017:592' id = [1,1;1,2;1,3;1,4;1,5;1,6;1,7;1,8;1,9;2,2;2,3;2,4;2,5;2,6;2,7;2,8;2,9;3,3;3,4;3,5;3,6;3,7;3,8;3,9;4,4;4,5;4,6;4,7;4,8;4,9;5,5;5,6;5,7;5,8;5,9;6,6;6,7;6,8;6,9;7,7;7,8;7,9;8,8;8,9;9,9]; */
    /* create a matrix of unique parameter choices for the cell cycle */
    /* 'IMT_analysis_April2017:594' P = zeros(length(id),4); */
    /* 'IMT_analysis_April2017:595' for ii = 1:length(id) */
    for (emgfitnomle = 0; emgfitnomle < 45; emgfitnomle++) {
      /* 'IMT_analysis_April2017:596' P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)}]; */
      for (itmp = 0; itmp < 2; itmp++) {
        c_P[emgfitnomle + 45 * itmp] = pcell[id[emgfitnomle] - 1].f1[itmp];
        c_P[emgfitnomle + 45 * (itmp + 2)] = pcell[id[45 + emgfitnomle] - 1]
          .f1[itmp];
      }
    }

    /* 'IMT_analysis_April2017:599' confint=nan(2,4,length(id)); */
    /* 'IMT_analysis_April2017:600' hp=nan(1,length(id)); */
    /* 'IMT_analysis_April2017:601' E=nan(1,length(id)); */
    /* 'IMT_analysis_April2017:602' hp_true=nan(1,length(id)); */
    /* 'IMT_analysis_April2017:603' flag_true=nan(1,length(id)); */
    /* 'IMT_analysis_April2017:604' E_true=nan(1,length(id)); */
    /*  optimize parameters */
    /* 'IMT_analysis_April2017:607' pd=zeros(length(P),4); */
    /* 'IMT_analysis_April2017:608' ld = NaN*ones(length(P),1); */
    /* 'IMT_analysis_April2017:609' options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs'); */
    /* 'IMT_analysis_April2017:610' flag=zeros(length(id),1); */
    /* 'IMT_analysis_April2017:611' for i=1:length(id) */
    for (emgfitnomle = 0; emgfitnomle < 45; emgfitnomle++) {
      /* 'IMT_analysis_April2017:612' fprintf('i=%f\n',i); */
      //d_cfprintf(1.0 + (double)emgfitnomle);
		printf("i=%f\n", 1.0 + (double)emgfitnomle);

      /* 'IMT_analysis_April2017:613' fprintf('  P=[%f %f %f %f]\n',P(i,1),P(i,2),P(i,3),P(i,4)); */
      //s_cfprintf(c_P[emgfitnomle], c_P[45 + emgfitnomle], c_P[90 + emgfitnomle],
      //           c_P[135 + emgfitnomle]);
	  printf("  P=[%f %f %f %f]\n", c_P[emgfitnomle], c_P[45 + emgfitnomle], c_P[90 + emgfitnomle], c_P[135 + emgfitnomle]);

      /* 'IMT_analysis_April2017:614' x0 = P(i,:); */
      /* f=@(x,m1,s1,m2,s2)convolv_2invG_adapt_nov(x,m1,s1,m2,s2,.01); */
      /* 'IMT_analysis_April2017:618' f=@(x)(sum(-log(convolv_2invG_adapt_nov(data,x(1),x(2),x(3),x(4),.01)))); */
      /*  p found best pdf parameters */
      /*  conf1  95% confidence intervals for the parameters */
      /* [p,conf1]=mle(data,'pdf',f,'start',x0, 'upperbound', [Inf Inf Inf Inf],'lowerbound',[0 0 0 0],'options',options) */
      /* [p,pval]=fmincon(f,x0,[],[],[],[],[0 0 0 0],[.2 4 .2 4],[],options); */
      /* 'IMT_analysis_April2017:628' [p,pval]=fminsearch(f, x0, options); */
      for (itmp = 0; itmp < 4; itmp++) {
        c_p[itmp] = c_P[emgfitnomle + 45 * itmp];
      }
#ifdef _GSL_MINIMIZE
	  // https://www.gnu.org/software/gsl/doc/html/multimin.html#algorithms-without-derivatives

	  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	  gsl_multimin_fminimizer *s = NULL;
	  gsl_vector *ss, *x;
	  gsl_multimin_function minex_func;

	  size_t iter = 0;
	  int status;
	  double size;

	  /* Starting point */
	  x = gsl_vector_alloc(4);
	  gsl_vector_set(x, 0, c_p[0]);
	  gsl_vector_set(x, 1, c_p[1]);
	  gsl_vector_set(x, 2, c_p[2]);
	  gsl_vector_set(x, 3, c_p[3]);


	  // Store previous state
	  gsl_vector *prev;
	  prev = gsl_vector_alloc(4);
	  gsl_vector_set_all(prev, 0);
	  int repeated = 0;


	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc(4);
	  gsl_vector_set_all(ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 4;
	  minex_func.f = convolv_2invG_adapt_nov_loglikelihood;
	  minex_func.params = (void *)data;

	  s = gsl_multimin_fminimizer_alloc(T, 4);
	  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	  do
	  {
		  iter++;
		  status = gsl_multimin_fminimizer_iterate(s);

		  if (status)
			  break;

		  size = gsl_multimin_fminimizer_size(s);
		  status = gsl_multimin_test_size(size, 1e-4);
		  
		  /*
		  // Detect repeated states
		  gsl_vector_sub(prev, s->x);
		  if (gsl_vector_isnull(prev))
			  repeated++;
		  else
			  repeated = 0;
		  if (repeated > 10){
			  printf("ERROR - QUITTING BECAUSE EXCESSIVE REPEATS!!!\n");
			  status = GSL_SUCCESS;
		  }
		  gsl_vector_set(prev, 0, gsl_vector_get(s->x, 0));
		  gsl_vector_set(prev, 1, gsl_vector_get(s->x, 1));
		  gsl_vector_set(prev, 2, gsl_vector_get(s->x, 2));
		  gsl_vector_set(prev, 3, gsl_vector_get(s->x, 3));
		   */ 


#ifdef _VERBOSE
		  if (status == GSL_SUCCESS)
		  {
			  printf("converged to minimum at\n");
		  }

		  printf("%5d %.17f %.17f %.17f %.17f f() = %.17f size = %.3f\n",
			  iter,
			  gsl_vector_get(s->x, 0),
			  gsl_vector_get(s->x, 1),
			  gsl_vector_get(s->x, 2),
			  gsl_vector_get(s->x, 3),
			  s->fval, size);
#endif
	  } while (status == GSL_CONTINUE && iter < 10000);

	  c_p[0] = fabs(gsl_vector_get(s->x, 0));
	  c_p[1] = fabs(gsl_vector_get(s->x, 1));
	  c_p[2] = fabs(gsl_vector_get(s->x, 2));
	  c_p[3] = fabs(gsl_vector_get(s->x, 3));

	  gsl_vector_free(x);
	  gsl_vector_free(ss);
	  gsl_multimin_fminimizer_free(s);
#else
      d_fminsearch(c_p);
#endif


	  printf("  p=[%f %f %f %f]\n", c_p[0], c_p[1], c_p[2], c_p[3]);

      /* 'IMT_analysis_April2017:631' pd(i,:)=p; */
      for (itmp = 0; itmp < 4; itmp++) {
        b_pd[emgfitnomle + 45 * itmp] = c_p[itmp];
      }


      /* 'IMT_analysis_April2017:639' [l,hp(i),flag(i),E(i)]=convolv_2invG_adapt_nov(data,p(1),p(2),p(3),p(4),.01); */
      b_convolv_2invG_adapt_nov(c_p[0], c_p[1], c_p[2], c_p[3], l, &mtmp, &flag, &E);

      /* 'IMT_analysis_April2017:640' l=sum(log(l)); */
      /* 'IMT_analysis_April2017:641' ld(i)=l; */
      /* 'IMT_analysis_April2017:642' fprintf('  l=%f hp=%f flag=%f E=%f\n\n',l,hp(i),flag(i),E(i)); */
      b_log(l);
	  printf("  l=%f hp=%f flag=%f E=%f\n\n", sum(l), mtmp, flag, E);
    }

    /*  we previously optimized with a larger step size, recalculate with */
    /*  a smaller stepsize after the fact */
	printf("recalculating canidate solutions with smaller stepsize\n");

    /* 'IMT_analysis_April2017:648' ld_true=zeros(length(ld),1); */
    /* 'IMT_analysis_April2017:649' for i=1:length(ld) */
    for (emgfitnomle = 0; emgfitnomle < 45; emgfitnomle++) {
      /* 'IMT_analysis_April2017:650' [l,hp_true(i),flag_true(i),E_true(i)]=convolv_2invG_adapt_nov(data,pd(i,1),pd(i,2),pd(i,3),pd(i,4),.001); */
      c_convolv_2invG_adapt_nov(b_pd[emgfitnomle], b_pd[45 + emgfitnomle], b_pd
        [90 + emgfitnomle], b_pd[135 + emgfitnomle], l, &mtmp, &flag, &E);

      /* 'IMT_analysis_April2017:651' ld_true(i)=sum(log(l)); */
      b_log(l);
      b_flag[emgfitnomle] = sum(l);
    }

    /*  common to each fit, consider factoring out */
    /* 'IMT_analysis_April2017:655' [max_ld,row_ld]=max(ld_true); */
    emgfitnomle = 1;
    mtmp = b_flag[0];
    itmp = 0;
    if (rtIsNaN(b_flag[0])) {
      ix = 1;
      exitg1 = false;
      while ((!exitg1) && (ix + 1 < 46)) {
        emgfitnomle = ix + 1;
        if (!rtIsNaN(b_flag[ix])) {
          mtmp = b_flag[ix];
          itmp = ix;
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (emgfitnomle < 45) {
      while (emgfitnomle + 1 < 46) {
        if (b_flag[emgfitnomle] > mtmp) {
          mtmp = b_flag[emgfitnomle];
          itmp = emgfitnomle;
        }

        emgfitnomle++;
      }
    }


	printf("max_ld=%f row_ld=%f\n", mtmp, (double) itmp + 1);
	printf("pd_max=[%f %f %f %f]\n\n", b_pd[itmp], b_pd[45 + itmp], b_pd[90 + itmp], b_pd[135 + itmp]);

    /*  END FUNCTION FIT_TWOSTAGE */
  }





































 
  /* Fit three-stage model with fminsearch */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* 'IMT_analysis_April2017:819' if threestagefitnomle == 1 */
  if (threestagefitnomle == 1) {
    /*  BEGIN FUNCTION FIT_THREESTAGE */
    /* 'IMT_analysis_April2017:821' fprintf('threestagefitnomle\n'); */
    //x_cfprintf();
	  printf("threestagefitnomle\n");

    /*  prepare statistical parameters */
    /* 'IMT_analysis_April2017:823' vry = [0.1 0.7]'; */
    /* 'IMT_analysis_April2017:824' c1 = C1*vry; */
    /* 'IMT_analysis_April2017:825' c2 = C2*vry; */
    /* 'IMT_analysis_April2017:826' m = 1./c1; */
    /* 'IMT_analysis_April2017:827' s = (c2./c1.^3).^0.5; */
    /* 'IMT_analysis_April2017:828' N = length(vry); */
    /*  prepare parameter seeds */
    /* pcomb = allcomb(m,s) */
    /* 'IMT_analysis_April2017:832' pcomb = [0.8494,0.5898;0.8494,0.0843;0.1213,0.5898;0.1213,0.0843]; */
    /* 'IMT_analysis_April2017:834' pcell = cell(length(pcomb),1); */
    /* 'IMT_analysis_April2017:835' for i = 1:length(pcomb) */
    /* 'IMT_analysis_April2017:836' pcell{i} = pcomb(i,:); */
    /* 'IMT_analysis_April2017:836' pcell{i} = pcomb(i,:); */
    /* 'IMT_analysis_April2017:836' pcell{i} = pcomb(i,:); */
    /* 'IMT_analysis_April2017:836' pcell{i} = pcomb(i,:); */
    for (itmp = 0; itmp < 2; itmp++) {
      b_pcell[0].f1[itmp] = 0.8494 + -0.25960000000000005 * (double)itmp;
      b_pcell[1].f1[itmp] = dv61[itmp];
      b_pcell[2].f1[itmp] = 0.1213 + 0.46849999999999997 * (double)itmp;
      b_pcell[3].f1[itmp] = 0.1213 + -0.037000000000000005 * (double)itmp;
    }

    /* id = allcomb(1:length(pcomb),1:length(pcomb),1:length(pcomb)); */
    /* id = sort(id,2); */
    /* id = unique(id,'rows') */
    /* 'IMT_analysis_April2017:841' id = [1,1,1;1,1,2;1,1,3;1,1,4;1,2,2;1,2,3;1,2,4;1,3,3;1,3,4;1,4,4;2,2,2;2,2,3;2,2,4;2,3,3;2,3,4;2,4,4;3,3,3;3,3,4;3,4,4;4,4,4]; */
    /* 'IMT_analysis_April2017:842' P = zeros(length(id),6); */
    /* 'IMT_analysis_April2017:843' for ii = 1:length(id) */
    /* 'IMT_analysis_April2017:847' hp=nan(1,length(id)); */
    /* 'IMT_analysis_April2017:848' E=nan(1,length(id)); */
    for (emgfitnomle = 0; emgfitnomle < 20; emgfitnomle++) {
      /* 'IMT_analysis_April2017:844' P(ii,:) = [pcell{id(ii,1)},pcell{id(ii,2)},pcell{id(ii,3)}]; */
      for (itmp = 0; itmp < 2; itmp++) {
        d_P[emgfitnomle + 20 * itmp] = b_pcell[b_id[emgfitnomle] - 1].f1[itmp];
        d_P[emgfitnomle + 20 * (itmp + 2)] = b_pcell[b_id[20 + emgfitnomle] - 1]
          .f1[itmp];
        d_P[emgfitnomle + 20 * (itmp + 4)] = b_pcell[b_id[40 + emgfitnomle] - 1]
          .f1[itmp];
      }
    }

    /* 'IMT_analysis_April2017:849' hp_true=nan(1,length(id)); */
    /* 'IMT_analysis_April2017:850' flag_true=nan(1,length(id)); */
    /* 'IMT_analysis_April2017:851' E_true=nan(1,length(id)); */
    /*  optimize parameters */
    /* 'IMT_analysis_April2017:854' pd=zeros(length(P),6); */
    memset(&c_pd[0], 0, 120U * sizeof(double));

    /* 'IMT_analysis_April2017:855' ld = NaN*ones(length(P),1); */
    /* 'IMT_analysis_April2017:856' flag=zeros(length(P),1); */
    /* options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs','Display','iter'); */
    /* 'IMT_analysis_April2017:858' options = statset('MaxIter',10000, 'MaxFunEvals',10000,'TolFun',1e-3,'TolX',1e-3,'TolTypeFun','rel', 'TolTypeX', 'abs'); */
    /*  WARNING WARNING WARNING */
    /*  I've skipped over P(18),P(19) and P(20) for the moment b/c P(18) */
    /*  is pathological. this needs fixed!!!! */
    /*  WARNING WARNING WARNING */
    /* 'IMT_analysis_April2017:863' for i=1:length(P)-3 */
    for (emgfitnomle = 0; emgfitnomle < 17; emgfitnomle++) {
      /* 'IMT_analysis_April2017:864' fprintf('i=%f\n',i); */
      //d_cfprintf(1.0 + (double)emgfitnomle);
		printf("i=%f\n", 1.0 + (double)emgfitnomle);

      /* 'IMT_analysis_April2017:865' fprintf('  P=[%f %f %f %f %f %f]\n',P(i,1),P(i,2),P(i,3),P(i,4),P(i,5),P(i,6)); */
      /*y_cfprintf(d_P[emgfitnomle], d_P[20 + emgfitnomle], d_P[40 + emgfitnomle],
                 d_P[60 + emgfitnomle], d_P[80 + emgfitnomle], d_P[100 +
                 emgfitnomle]);*/
		printf("  P=[%f %f %f %f %f %f]\n", d_P[emgfitnomle], d_P[20 + emgfitnomle], d_P[40 + emgfitnomle],
			d_P[60 + emgfitnomle], d_P[80 + emgfitnomle], d_P[100 +
			emgfitnomle]);

      /* 'IMT_analysis_April2017:866' x0=P(i,:); */
      /* g=@(x,m1,s1,m2,s2,m3,s3)convolv_3invG_nov(x,m1,s1,m2,s2,m3,s3,.01); */
      /* [p,conf1]=mle(data,'pdf',g,'start',x0, 'upperbound', [Inf Inf Inf Inf Inf Inf],'lowerbound',[0 0 0 0 0 0],'options',options); */
      /* 'IMT_analysis_April2017:870' f=@(x)(sum(-log(convolv_3invG_nov(data,x(1),x(2),x(3),x(4),x(5),x(6),0.1)))); */
      /* 'IMT_analysis_April2017:871' [p,pval]=fminsearch(f, x0, options); */
      for (itmp = 0; itmp < 6; itmp++) {
        d_p[itmp] = d_P[emgfitnomle + 20 * itmp];
      }


#ifdef _GSL_MINIMIZE
	  // https://www.gnu.org/software/gsl/doc/html/multimin.html#algorithms-without-derivatives

	  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	  gsl_multimin_fminimizer *s = NULL;
	  gsl_vector *ss, *x;
	  gsl_multimin_function minex_func;

	  size_t iter = 0;
	  int status;
	  double size;

	  /* Starting point */
	  x = gsl_vector_alloc(6);
	  gsl_vector_set(x, 0, d_p[0]);
	  gsl_vector_set(x, 1, d_p[1]);
	  gsl_vector_set(x, 2, d_p[2]);
	  gsl_vector_set(x, 3, d_p[3]);
	  gsl_vector_set(x, 4, d_p[4]);
	  gsl_vector_set(x, 5, d_p[5]);

	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc(6);
	  gsl_vector_set_all(ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 6;
	  minex_func.f = convolv_3invG_nov_loglikelihood;
	  minex_func.params = (void *)data;

	  s = gsl_multimin_fminimizer_alloc(T, 6);
	  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

	  do
	  {
		  iter++;
		  status = gsl_multimin_fminimizer_iterate(s);

		  if (status)
			  break;

		  size = gsl_multimin_fminimizer_size(s);
		  status = gsl_multimin_test_size(size, 1e-4);

#ifdef _VERBOSE
		  if (status == GSL_SUCCESS)
		  {
			  printf("converged to minimum at\n");
		  }

		  printf("%5d %.17f %.17f %.17f %.17f %.17f %.17f f() = %.17f size = %.3f\n",
			  iter,
			  gsl_vector_get(s->x, 0),
			  gsl_vector_get(s->x, 1),
			  gsl_vector_get(s->x, 2),
			  gsl_vector_get(s->x, 3),
			  gsl_vector_get(s->x, 4),
			  gsl_vector_get(s->x, 5),
			  s->fval, size);
#endif
	  } while (status == GSL_CONTINUE && iter < 10000);

	  d_p[0] = fabs(gsl_vector_get(s->x, 0));
	  d_p[1] = fabs(gsl_vector_get(s->x, 1));
	  d_p[2] = fabs(gsl_vector_get(s->x, 2));
	  d_p[3] = fabs(gsl_vector_get(s->x, 3));
	  d_p[4] = fabs(gsl_vector_get(s->x, 4));
	  d_p[5] = fabs(gsl_vector_get(s->x, 5));

	  gsl_vector_free(x);
	  gsl_vector_free(ss);
	  gsl_multimin_fminimizer_free(s);
#else
      e_fminsearch(d_p);
#endif

      /* 'IMT_analysis_April2017:873' fprintf('  p=[%f %f %f %f %f %f]\n',p(1),p(2),p(3),p(4),p(5),p(6)); */
      //ab_cfprintf(d_p[0], d_p[1], d_p[2], d_p[3], d_p[4], d_p[5]);
	  printf("  p=[%f %f %f %f %f %f]\n", d_p[0], d_p[1], d_p[2], d_p[3], d_p[4], d_p[5]);


      /* 'IMT_analysis_April2017:874' pd(i,:)=p; */
      for (itmp = 0; itmp < 6; itmp++) {
        c_pd[emgfitnomle + 20 * itmp] = d_p[itmp];
      }

      /* confint(:,:,i)=conf1(:); */
      /* 'IMT_analysis_April2017:876' [l,hp(i),flag(i),E(i)]=convolv_3invG_nov(data,p(1),p(2),p(3),p(4),p(5),p(6),.01); */
      b_convolv_3invG_nov(d_p[0], d_p[1], d_p[2], d_p[3], d_p[4], d_p[5], l,
                          &mtmp, &flag, &E);

      /* 'IMT_analysis_April2017:877' l=sum(log(l)); */
      /* 'IMT_analysis_April2017:878' ld(i)=l; */
      /* 'IMT_analysis_April2017:879' fprintf('  l=%f hp=%f flag=%f E=%f\n\n',l,hp(i),flag(i),E(i)); */
      b_log(l);
      //u_cfprintf(sum(l), mtmp, flag, E);
	  printf("  l=%f hp=%f flag=%f E=%f\n\n", sum(l), mtmp, flag, E);
    }

    /*  we previously optimized with a larger step size, recalculate with */
    /*  a smaller stepsize after the fact */
    /* 'IMT_analysis_April2017:884' fprintf('Recalculating canidate solutions with smaller stepsize\n'); */
    //bb_cfprintf();
	printf("Recalculating canidate solutions with smaller stepsize\n");

    /* 'IMT_analysis_April2017:885' ld_true=zeros(length(ld),1); */
    /* 'IMT_analysis_April2017:886' for i=1:length(ld) */
    for (emgfitnomle = 0; emgfitnomle < 20; emgfitnomle++) {
      /* 'IMT_analysis_April2017:887' [l,hp_true(i),flag_true(i),E_true(i)]=convolv_3invG_nov(data,pd(i,1),pd(i,2),pd(i,3),pd(i,4),pd(i,5),pd(i,6),.001); */
      c_convolv_3invG_nov(c_pd[emgfitnomle], c_pd[20 + emgfitnomle], c_pd[40 +
                          emgfitnomle], c_pd[60 + emgfitnomle], c_pd[80 +
                          emgfitnomle], c_pd[100 + emgfitnomle], l, &mtmp, &flag,
                          &E);

      /* 'IMT_analysis_April2017:888' ld_true(i)=sum(log(l)); */
      b_log(l);
      c_flag[emgfitnomle] = sum(l);
    }

    /*  common to each fit, consider factoring out */
    /* 'IMT_analysis_April2017:892' [max_ld,row_ld]=max(ld_true); */
    emgfitnomle = 1;
    mtmp = c_flag[0];
    itmp = 0;
    if (rtIsNaN(c_flag[0])) {
      ix = 1;
      exitg1 = false;
      while ((!exitg1) && (ix + 1 < 21)) {
        emgfitnomle = ix + 1;
        if (!rtIsNaN(c_flag[ix])) {
          mtmp = c_flag[ix];
          itmp = ix;
          exitg1 = true;
        } else {
          ix++;
        }
      }
    }

    if (emgfitnomle < 20) {
      while (emgfitnomle + 1 < 21) {
        if (c_flag[emgfitnomle] > mtmp) {
          mtmp = c_flag[emgfitnomle];
          itmp = emgfitnomle;
        }

        emgfitnomle++;
      }
    }

    /* 'IMT_analysis_April2017:893' pd_max = pd(row_ld,:); */
    /* confint_max=confint(:,:,row_ld); */
    /* 'IMT_analysis_April2017:895' fprintf('max_ld=%f row_ld=%f\n',max_ld,row_ld); */
    //m_cfprintf(mtmp, itmp + 1);
	printf("max_ld=%f row_ld=%f\n", mtmp, (double) itmp + 1);

    /* 'IMT_analysis_April2017:896' fprintf('pd_max=[%f %f %f %f %f %f]\n\n',pd_max(1),pd_max(2),pd_max(3),pd_max(4),pd_max(5),pd_max(6)); */
    /*cb_cfprintf(c_pd[itmp], c_pd[20 + itmp], c_pd[40 + itmp], c_pd[60 + itmp],
                c_pd[80 + itmp], c_pd[100 + itmp]);
	printf("pd_max=[%f %f %f %f %f %f]\n\n", c_pd[itmp], c_pd[20 + itmp], c_pd[40 + itmp], c_pd[60 + itmp],
		c_pd[80 + itmp], c_pd[100 + itmp]);

    /*  END FUNCTION FIT_THREESTAGE */
  }
}

/* End of code generation (IMT_analysis_April2017.c) */
