/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzgeev.c
 *
 * Code generation for function 'xzgeev'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "xzgeev.h"
#include "xzlartg.h"
#include "xzhgeqz.h"
#include "gp_max.h"
#include "IMT_analysis_April2017_rtwutil.h"

/* Function Definitions */

/*
 *
 */
void xzgeev(const creal_T A_data[], const int A_size[2], int *info, creal_T
            alpha1_data[], int alpha1_size[1], creal_T beta1_data[], int
            beta1_size[1])
{
  int At_size[2];
  int ii;
  int nzcount;
  creal_T At_data[16];
  double anrm;
  boolean_T exitg1;
  double absxk;
  boolean_T ilascl;
  double anrmto;
  int ilo;
  double ctoc;
  int ihi;
  boolean_T notdone;
  int exitg3;
  double stemp_im;
  int i;
  double cto1;
  int j;
  double mul;
  creal_T atmp;
  int jj;
  boolean_T exitg4;
  int exitg2;
  At_size[0] = A_size[0];
  At_size[1] = A_size[1];
  ii = A_size[0] * A_size[1];
  for (nzcount = 0; nzcount < ii; nzcount++) {
    At_data[nzcount] = A_data[nzcount];
  }

  *info = 0;
  anrm = 0.0;
  nzcount = A_size[0] * A_size[1];
  ii = 0;
  exitg1 = false;
  while ((!exitg1) && (ii <= nzcount - 1)) {
    absxk = rt_hypotd_snf(At_data[ii].re, At_data[ii].im);
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg1 = true;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      ii++;
    }
  }

  if (!((!rtIsInf(anrm)) && (!rtIsNaN(anrm)))) {
    alpha1_size[0] = A_size[0];
    ii = A_size[0];
    for (nzcount = 0; nzcount < ii; nzcount++) {
      alpha1_data[nzcount].re = rtNaN;
      alpha1_data[nzcount].im = 0.0;
    }

    beta1_size[0] = A_size[0];
    ii = A_size[0];
    for (nzcount = 0; nzcount < ii; nzcount++) {
      beta1_data[nzcount].re = rtNaN;
      beta1_data[nzcount].im = 0.0;
    }
  } else {
    ilascl = false;
    anrmto = anrm;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = true;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = true;
      }
    }

    if (ilascl) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = true;
      while (notdone) {
        stemp_im = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((stemp_im > ctoc) && (ctoc != 0.0)) {
          mul = 2.0041683600089728E-292;
          absxk = stemp_im;
        } else if (cto1 > absxk) {
          mul = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          mul = ctoc / absxk;
          notdone = false;
        }

        ii = At_size[0] * At_size[1];
        for (nzcount = 0; nzcount < ii; nzcount++) {
          At_data[nzcount].re *= mul;
          At_data[nzcount].im *= mul;
        }
      }
    }

    ilo = 0;
    ihi = A_size[0];
    if (A_size[0] <= 1) {
      ihi = 1;
    } else {
      do {
        exitg3 = 0;
        i = 0;
        j = 0;
        notdone = false;
        ii = ihi;
        exitg1 = false;
        while ((!exitg1) && (ii > 0)) {
          nzcount = 0;
          i = ii;
          j = ihi;
          jj = 1;
          exitg4 = false;
          while ((!exitg4) && (jj <= ihi)) {
            if ((At_data[(ii + At_size[0] * (jj - 1)) - 1].re != 0.0) ||
                (At_data[(ii + At_size[0] * (jj - 1)) - 1].im != 0.0) || (ii ==
                 jj)) {
              if (nzcount == 0) {
                j = jj;
                nzcount = 1;
                jj++;
              } else {
                nzcount = 2;
                exitg4 = true;
              }
            } else {
              jj++;
            }
          }

          if (nzcount < 2) {
            notdone = true;
            exitg1 = true;
          } else {
            ii--;
          }
        }

        if (!notdone) {
          exitg3 = 2;
        } else {
          if (i != ihi) {
            for (ii = 0; ii + 1 <= At_size[0]; ii++) {
              atmp = At_data[(i + At_size[0] * ii) - 1];
              At_data[(i + At_size[0] * ii) - 1] = At_data[(ihi + At_size[0] *
                ii) - 1];
              At_data[(ihi + At_size[0] * ii) - 1] = atmp;
            }
          }

          if (j != ihi) {
            for (ii = 0; ii + 1 <= ihi; ii++) {
              atmp = At_data[ii + At_size[0] * (j - 1)];
              At_data[ii + At_size[0] * (j - 1)] = At_data[ii + At_size[0] *
                (ihi - 1)];
              At_data[ii + At_size[0] * (ihi - 1)] = atmp;
            }
          }

          ihi--;
          if (ihi == 1) {
            exitg3 = 1;
          }
        }
      } while (exitg3 == 0);

      if (exitg3 == 1) {
      } else {
        do {
          exitg2 = 0;
          i = 0;
          j = 0;
          notdone = false;
          jj = ilo + 1;
          exitg1 = false;
          while ((!exitg1) && (jj <= ihi)) {
            nzcount = 0;
            i = ihi;
            j = jj;
            ii = ilo + 1;
            exitg4 = false;
            while ((!exitg4) && (ii <= ihi)) {
              if ((At_data[(ii + At_size[0] * (jj - 1)) - 1].re != 0.0) ||
                  (At_data[(ii + At_size[0] * (jj - 1)) - 1].im != 0.0) || (ii ==
                   jj)) {
                if (nzcount == 0) {
                  i = ii;
                  nzcount = 1;
                  ii++;
                } else {
                  nzcount = 2;
                  exitg4 = true;
                }
              } else {
                ii++;
              }
            }

            if (nzcount < 2) {
              notdone = true;
              exitg1 = true;
            } else {
              jj++;
            }
          }

          if (!notdone) {
            exitg2 = 1;
          } else {
            if (i != ilo + 1) {
              for (ii = ilo; ii + 1 <= At_size[0]; ii++) {
                atmp = At_data[(i + At_size[0] * ii) - 1];
                At_data[(i + At_size[0] * ii) - 1] = At_data[ilo + At_size[0] *
                  ii];
                At_data[ilo + At_size[0] * ii] = atmp;
              }
            }

            if (j != ilo + 1) {
              for (ii = 0; ii + 1 <= ihi; ii++) {
                atmp = At_data[ii + At_size[0] * (j - 1)];
                At_data[ii + At_size[0] * (j - 1)] = At_data[ii + At_size[0] *
                  ilo];
                At_data[ii + At_size[0] * ilo] = atmp;
              }
            }

            ilo++;
            if (ilo + 1 == ihi) {
              exitg2 = 1;
            }
          }
        } while (exitg2 == 0);
      }
    }

    if ((!(A_size[0] <= 1)) && (!(ihi < ilo + 3))) {
      for (ii = ilo; ii + 1 < ihi - 1; ii++) {
        for (nzcount = ihi - 1; nzcount + 1 > ii + 2; nzcount--) {
          xzlartg(At_data[(nzcount + At_size[0] * ii) - 1], At_data[nzcount +
                  At_size[0] * ii], &absxk, &atmp, &At_data[(nzcount + At_size[0]
                   * ii) - 1]);
          At_data[nzcount + At_size[0] * ii].re = 0.0;
          At_data[nzcount + At_size[0] * ii].im = 0.0;
          for (j = ii + 1; j + 1 <= At_size[0]; j++) {
            ctoc = absxk * At_data[(nzcount + At_size[0] * j) - 1].re + (atmp.re
              * At_data[nzcount + At_size[0] * j].re - atmp.im * At_data[nzcount
              + At_size[0] * j].im);
            stemp_im = absxk * At_data[(nzcount + At_size[0] * j) - 1].im +
              (atmp.re * At_data[nzcount + At_size[0] * j].im + atmp.im *
               At_data[nzcount + At_size[0] * j].re);
            cto1 = At_data[(nzcount + At_size[0] * j) - 1].re;
            At_data[nzcount + At_size[0] * j].re = absxk * At_data[nzcount +
              At_size[0] * j].re - (atmp.re * At_data[(nzcount + At_size[0] * j)
              - 1].re + atmp.im * At_data[(nzcount + At_size[0] * j) - 1].im);
            At_data[nzcount + At_size[0] * j].im = absxk * At_data[nzcount +
              At_size[0] * j].im - (atmp.re * At_data[(nzcount + At_size[0] * j)
              - 1].im - atmp.im * cto1);
            At_data[(nzcount + At_size[0] * j) - 1].re = ctoc;
            At_data[(nzcount + At_size[0] * j) - 1].im = stemp_im;
          }

          atmp.re = -atmp.re;
          atmp.im = -atmp.im;
          for (i = 0; i + 1 <= ihi; i++) {
            ctoc = absxk * At_data[i + At_size[0] * nzcount].re + (atmp.re *
              At_data[i + At_size[0] * (nzcount - 1)].re - atmp.im * At_data[i +
              At_size[0] * (nzcount - 1)].im);
            stemp_im = absxk * At_data[i + At_size[0] * nzcount].im + (atmp.re *
              At_data[i + At_size[0] * (nzcount - 1)].im + atmp.im * At_data[i +
              At_size[0] * (nzcount - 1)].re);
            cto1 = At_data[i + At_size[0] * nzcount].re;
            At_data[i + At_size[0] * (nzcount - 1)].re = absxk * At_data[i +
              At_size[0] * (nzcount - 1)].re - (atmp.re * At_data[i + At_size[0]
              * nzcount].re + atmp.im * At_data[i + At_size[0] * nzcount].im);
            At_data[i + At_size[0] * (nzcount - 1)].im = absxk * At_data[i +
              At_size[0] * (nzcount - 1)].im - (atmp.re * At_data[i + At_size[0]
              * nzcount].im - atmp.im * cto1);
            At_data[i + At_size[0] * nzcount].re = ctoc;
            At_data[i + At_size[0] * nzcount].im = stemp_im;
          }
        }
      }
    }

    xzhgeqz(At_data, At_size, ilo + 1, ihi, info, alpha1_data, alpha1_size,
            beta1_data, beta1_size);
    if ((*info == 0) && ilascl) {
      notdone = true;
      while (notdone) {
        stemp_im = anrmto * 2.0041683600089728E-292;
        cto1 = anrm / 4.9896007738368E+291;
        if ((stemp_im > anrm) && (anrm != 0.0)) {
          mul = 2.0041683600089728E-292;
          anrmto = stemp_im;
        } else if (cto1 > anrmto) {
          mul = 4.9896007738368E+291;
          anrm = cto1;
        } else {
          mul = anrm / anrmto;
          notdone = false;
        }

        ii = alpha1_size[0];
        for (nzcount = 0; nzcount < ii; nzcount++) {
          alpha1_data[nzcount].re *= mul;
          alpha1_data[nzcount].im *= mul;
        }
      }
    }
  }
}

/* End of code generation (xzgeev.c) */
