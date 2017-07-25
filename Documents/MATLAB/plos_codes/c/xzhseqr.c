/*
 * Student License - for use by students to meet course requirements and
 * perform academic research at degree granting institutions only.  Not
 * for government, commercial, or other organizational use.
 *
 * xzhseqr.c
 *
 * Code generation for function 'xzhseqr'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "IMT_analysis_April2017.h"
#include "xzhseqr.h"
#include "gp_max.h"
#include "xscal.h"
#include "xzlarfg.h"
#include "sqrt.h"
#include "IMT_analysis_April2017_rtwutil.h"

/* Function Definitions */

/*
 *
 */
int eml_zlahqr(creal_T h_data[], int h_size[2])
{
  int info;
  int n;
  int ldh;
  int j;
  int i;
  double SMLNUM;
  double tst;
  boolean_T exitg1;
  double htmp1;
  double ba;
  int L;
  creal_T u2;
  boolean_T goto140;
  int its;
  boolean_T exitg2;
  int k;
  boolean_T exitg3;
  int i4;
  double t_re;
  double ab;
  creal_T y;
  boolean_T goto70;
  int m;
  double u_re;
  double u_im;
  double aa;
  double s;
  int b_k;
  double b_SMLNUM;
  creal_T v[2];
  double b_u_re;
  n = h_size[0];
  ldh = h_size[0];
  info = 0;
  if (1 != h_size[0]) {
    j = 1;
    while (j <= n - 3) {
      h_data[2].re = 0.0;
      h_data[2].im = 0.0;
      h_data[3].re = 0.0;
      h_data[3].im = 0.0;
      j = 2;
    }

    if (1 <= n - 2) {
      h_data[(n + h_size[0] * (n - 3)) - 1].re = 0.0;
      h_data[(n + h_size[0] * (n - 3)) - 1].im = 0.0;
    }

    for (i = 1; i + 1 <= n; i++) {
      if (h_data[i + h_size[0] * (i - 1)].im != 0.0) {
        tst = h_data[i + h_size[0] * (i - 1)].re;
        htmp1 = h_data[i + h_size[0] * (i - 1)].im;
        ba = fabs(h_data[i + h_size[0] * (i - 1)].re) + fabs(h_data[i + h_size[0]
          * (i - 1)].im);
        if (htmp1 == 0.0) {
          u2.re = tst / ba;
          u2.im = 0.0;
        } else if (tst == 0.0) {
          u2.re = 0.0;
          u2.im = htmp1 / ba;
        } else {
          u2.re = tst / ba;
          u2.im = htmp1 / ba;
        }

        ba = rt_hypotd_snf(u2.re, u2.im);
        if (-u2.im == 0.0) {
          u2.re /= ba;
          u2.im = 0.0;
        } else if (u2.re == 0.0) {
          u2.re = 0.0;
          u2.im = -u2.im / ba;
        } else {
          u2.re /= ba;
          u2.im = -u2.im / ba;
        }

        h_data[i + h_size[0] * (i - 1)].re = rt_hypotd_snf(h_data[i + h_size[0] *
          (i - 1)].re, h_data[i + h_size[0] * (i - 1)].im);
        h_data[i + h_size[0] * (i - 1)].im = 0.0;
        L = i + i * ldh;
        i4 = (L + ldh * ((n - i) - 1)) + 1;
        while (L + 1 <= i4) {
          tst = h_data[L].re;
          htmp1 = h_data[L].im;
          h_data[L].re = u2.re * tst - u2.im * htmp1;
          h_data[L].im = u2.re * htmp1 + u2.im * tst;
          L += ldh;
        }

        L = i * ldh;
        u2.im = -u2.im;
        j = i + 2;
        if (n < j) {
          j = n;
        }

        i4 = L + j;
        while (L + 1 <= i4) {
          tst = h_data[L].re;
          htmp1 = h_data[L].im;
          h_data[L].re = u2.re * tst - u2.im * htmp1;
          h_data[L].im = u2.re * htmp1 + u2.im * tst;
          L++;
        }
      }
    }

    SMLNUM = 2.2250738585072014E-308 * ((double)n / 2.2204460492503131E-16);
    i = n - 1;
    exitg1 = false;
    while ((!exitg1) && (i + 1 >= 1)) {
      L = -1;
      goto140 = false;
      its = 0;
      exitg2 = false;
      while ((!exitg2) && (its < 31)) {
        k = i;
        exitg3 = false;
        while ((!exitg3) && ((k + 1 > L + 2) && (!(fabs(h_data[k + h_size[0] *
                   (k - 1)].re) + fabs(h_data[k + h_size[0] * (k - 1)].im) <=
                  SMLNUM)))) {
          tst = (fabs(h_data[(k + h_size[0] * (k - 1)) - 1].re) + fabs(h_data[(k
                   + h_size[0] * (k - 1)) - 1].im)) + (fabs(h_data[k + h_size[0]
            * k].re) + fabs(h_data[k + h_size[0] * k].im));
          if (tst == 0.0) {
            if (k - 1 >= 1) {
              tst = fabs(h_data[(k + h_size[0] * (k - 2)) - 1].re);
            }

            if (k + 2 <= n) {
              tst += fabs(h_data[(k + h_size[0] * k) + 1].re);
            }
          }

          if (fabs(h_data[k + h_size[0] * (k - 1)].re) <= 2.2204460492503131E-16
              * tst) {
            htmp1 = fabs(h_data[k + h_size[0] * (k - 1)].re) + fabs(h_data[k +
              h_size[0] * (k - 1)].im);
            tst = fabs(h_data[(k + h_size[0] * k) - 1].re) + fabs(h_data[(k +
              h_size[0] * k) - 1].im);
            if (htmp1 > tst) {
              ab = htmp1;
              ba = tst;
            } else {
              ab = tst;
              ba = htmp1;
            }

            htmp1 = fabs(h_data[k + h_size[0] * k].re) + fabs(h_data[k + h_size
              [0] * k].im);
            tst = fabs(h_data[(k + h_size[0] * (k - 1)) - 1].re - h_data[k +
                       h_size[0] * k].re) + fabs(h_data[(k + h_size[0] * (k - 1))
              - 1].im - h_data[k + h_size[0] * k].im);
            if (htmp1 > tst) {
              aa = htmp1;
              htmp1 = tst;
            } else {
              aa = tst;
            }

            s = aa + ab;
            tst = 2.2204460492503131E-16 * (htmp1 * (aa / s));
            if ((SMLNUM > tst) || rtIsNaN(tst)) {
              b_SMLNUM = SMLNUM;
            } else {
              b_SMLNUM = tst;
            }

            if (ba * (ab / s) <= b_SMLNUM) {
              exitg3 = true;
            } else {
              k--;
            }
          } else {
            k--;
          }
        }

        L = k - 1;
        if (k + 1 > 1) {
          h_data[k + h_size[0] * (k - 1)].re = 0.0;
          h_data[k + h_size[0] * (k - 1)].im = 0.0;
        }

        if (k + 1 >= i + 1) {
          goto140 = true;
          exitg2 = true;
        } else {
          if (its == 10) {
            t_re = 0.75 * fabs(h_data[(k + h_size[0] * k) + 1].re) + h_data[k +
              h_size[0] * k].re;
            ba = h_data[k + h_size[0] * k].im;
          } else if (its == 20) {
            t_re = 0.75 * fabs(h_data[i + h_size[0] * (i - 1)].re) + h_data[i +
              h_size[0] * i].re;
            ba = h_data[i + h_size[0] * i].im;
          } else {
            t_re = h_data[i + h_size[0] * i].re;
            ba = h_data[i + h_size[0] * i].im;
            y = h_data[(i + h_size[0] * i) - 1];
            b_sqrt(&y);
            u2 = h_data[i + h_size[0] * (i - 1)];
            b_sqrt(&u2);
            u_re = y.re * u2.re - y.im * u2.im;
            u_im = y.re * u2.im + y.im * u2.re;
            s = fabs(u_re) + fabs(u_im);
            if (s != 0.0) {
              htmp1 = 0.5 * (h_data[(i + h_size[0] * (i - 1)) - 1].re - h_data[i
                             + h_size[0] * i].re);
              ab = 0.5 * (h_data[(i + h_size[0] * (i - 1)) - 1].im - h_data[i +
                          h_size[0] * i].im);
              aa = fabs(htmp1) + fabs(ab);
              tst = fabs(htmp1) + fabs(ab);
              if (!((s > tst) || rtIsNaN(tst))) {
                s = tst;
              }

              if (ab == 0.0) {
                t_re = htmp1 / s;
                ba = 0.0;
              } else if (htmp1 == 0.0) {
                t_re = 0.0;
                ba = ab / s;
              } else {
                t_re = htmp1 / s;
                ba = ab / s;
              }

              tst = t_re;
              t_re = t_re * t_re - ba * ba;
              ba = tst * ba + ba * tst;
              if (u_im == 0.0) {
                u2.re = u_re / s;
                u2.im = 0.0;
              } else if (u_re == 0.0) {
                u2.re = 0.0;
                u2.im = u_im / s;
              } else {
                u2.re = u_re / s;
                u2.im = u_im / s;
              }

              y.re = t_re + (u2.re * u2.re - u2.im * u2.im);
              y.im = ba + (u2.re * u2.im + u2.im * u2.re);
              b_sqrt(&y);
              y.re *= s;
              y.im *= s;
              if (aa > 0.0) {
                if (ab == 0.0) {
                  t_re = htmp1 / aa;
                  ba = 0.0;
                } else if (htmp1 == 0.0) {
                  t_re = 0.0;
                  ba = ab / aa;
                } else {
                  t_re = htmp1 / aa;
                  ba = ab / aa;
                }

                if (t_re * y.re + ba * y.im < 0.0) {
                  y.re = -y.re;
                  y.im = -y.im;
                }
              }

              ba = htmp1 + y.re;
              aa = ab + y.im;
              if (aa == 0.0) {
                if (u_im == 0.0) {
                  b_u_re = u_re / ba;
                  tst = 0.0;
                } else if (u_re == 0.0) {
                  b_u_re = 0.0;
                  tst = u_im / ba;
                } else {
                  b_u_re = u_re / ba;
                  tst = u_im / ba;
                }
              } else if (ba == 0.0) {
                if (u_re == 0.0) {
                  b_u_re = u_im / aa;
                  tst = 0.0;
                } else if (u_im == 0.0) {
                  b_u_re = 0.0;
                  tst = -(u_re / aa);
                } else {
                  b_u_re = u_im / aa;
                  tst = -(u_re / aa);
                }
              } else {
                ab = fabs(ba);
                tst = fabs(aa);
                if (ab > tst) {
                  s = aa / ba;
                  tst = ba + s * aa;
                  b_u_re = (u_re + s * u_im) / tst;
                  tst = (u_im - s * u_re) / tst;
                } else if (tst == ab) {
                  if (ba > 0.0) {
                    htmp1 = 0.5;
                  } else {
                    htmp1 = -0.5;
                  }

                  if (aa > 0.0) {
                    tst = 0.5;
                  } else {
                    tst = -0.5;
                  }

                  b_u_re = (u_re * htmp1 + u_im * tst) / ab;
                  tst = (u_im * htmp1 - u_re * tst) / ab;
                } else {
                  s = ba / aa;
                  tst = aa + s * ba;
                  b_u_re = (s * u_re + u_im) / tst;
                  tst = (s * u_im - u_re) / tst;
                }
              }

              t_re = h_data[i + h_size[0] * i].re - (u_re * b_u_re - u_im * tst);
              ba = h_data[i + h_size[0] * i].im - (u_re * tst + u_im * b_u_re);
            }
          }

          goto70 = false;
          m = i;
          exitg3 = false;
          while ((!exitg3) && (m > k + 1)) {
            u2.re = h_data[(m + h_size[0] * (m - 1)) - 1].re - t_re;
            u2.im = h_data[(m + h_size[0] * (m - 1)) - 1].im - ba;
            tst = h_data[m + h_size[0] * (m - 1)].re;
            s = (fabs(u2.re) + fabs(u2.im)) + fabs(tst);
            if (u2.im == 0.0) {
              u2.re /= s;
              u2.im = 0.0;
            } else if (u2.re == 0.0) {
              u2.re = 0.0;
              u2.im /= s;
            } else {
              u2.re /= s;
              u2.im /= s;
            }

            tst /= s;
            v[0] = u2;
            v[1].re = tst;
            v[1].im = 0.0;
            if (fabs(h_data[(m + h_size[0] * (m - 2)) - 1].re) * fabs(tst) <=
                2.2204460492503131E-16 * ((fabs(u2.re) + fabs(u2.im)) * ((fabs
                   (h_data[(m + h_size[0] * (m - 1)) - 1].re) + fabs(h_data[(m +
                     h_size[0] * (m - 1)) - 1].im)) + (fabs(h_data[m + h_size[0]
                    * m].re) + fabs(h_data[m + h_size[0] * m].im))))) {
              goto70 = true;
              exitg3 = true;
            } else {
              m--;
            }
          }

          if (!goto70) {
            u2.re = h_data[k + h_size[0] * k].re - t_re;
            u2.im = h_data[k + h_size[0] * k].im - ba;
            tst = h_data[(k + h_size[0] * k) + 1].re;
            s = (fabs(u2.re) + fabs(u2.im)) + fabs(tst);
            if (u2.im == 0.0) {
              u2.re /= s;
              u2.im = 0.0;
            } else if (u2.re == 0.0) {
              u2.re = 0.0;
              u2.im /= s;
            } else {
              u2.re /= s;
              u2.im /= s;
            }

            tst /= s;
            v[0] = u2;
            v[1].re = tst;
            v[1].im = 0.0;
          }

          for (b_k = m; b_k <= i; b_k++) {
            if (b_k > m) {
              v[0] = h_data[(b_k + h_size[0] * (b_k - 2)) - 1];
              v[1] = h_data[b_k + h_size[0] * (b_k - 2)];
            }

            u2 = xzlarfg(&v[0], &v[1]);
            if (b_k > m) {
              h_data[(b_k + h_size[0] * (b_k - 2)) - 1] = v[0];
              h_data[b_k + h_size[0] * (b_k - 2)].re = 0.0;
              h_data[b_k + h_size[0] * (b_k - 2)].im = 0.0;
            }

            tst = u2.re * v[1].re - u2.im * v[1].im;
            for (j = b_k - 1; j + 1 <= n; j++) {
              t_re = (u2.re * h_data[(b_k + h_size[0] * j) - 1].re - -u2.im *
                      h_data[(b_k + h_size[0] * j) - 1].im) + tst * h_data[b_k +
                h_size[0] * j].re;
              ba = (u2.re * h_data[(b_k + h_size[0] * j) - 1].im + -u2.im *
                    h_data[(b_k + h_size[0] * j) - 1].re) + tst * h_data[b_k +
                h_size[0] * j].im;
              h_data[(b_k + h_size[0] * j) - 1].re -= t_re;
              h_data[(b_k + h_size[0] * j) - 1].im -= ba;
              h_data[b_k + h_size[0] * j].re -= t_re * v[1].re - ba * v[1].im;
              h_data[b_k + h_size[0] * j].im -= t_re * v[1].im + ba * v[1].re;
            }

            if (b_k + 2 < i + 1) {
              i4 = b_k;
            } else {
              i4 = i - 1;
            }

            for (j = 0; j + 1 <= i4 + 2; j++) {
              t_re = (u2.re * h_data[j + h_size[0] * (b_k - 1)].re - u2.im *
                      h_data[j + h_size[0] * (b_k - 1)].im) + tst * h_data[j +
                h_size[0] * b_k].re;
              ba = (u2.re * h_data[j + h_size[0] * (b_k - 1)].im + u2.im *
                    h_data[j + h_size[0] * (b_k - 1)].re) + tst * h_data[j +
                h_size[0] * b_k].im;
              h_data[j + h_size[0] * (b_k - 1)].re -= t_re;
              h_data[j + h_size[0] * (b_k - 1)].im -= ba;
              h_data[j + h_size[0] * b_k].re -= t_re * v[1].re - ba * -v[1].im;
              h_data[j + h_size[0] * b_k].im -= t_re * -v[1].im + ba * v[1].re;
            }

            if ((b_k == m) && (m > k + 1)) {
              u2.re = 1.0 - u2.re;
              u2.im = 0.0 - u2.im;
              ba = rt_hypotd_snf(u2.re, u2.im);
              if (u2.im == 0.0) {
                u2.re /= ba;
                u2.im = 0.0;
              } else if (u2.re == 0.0) {
                u2.re = 0.0;
                u2.im /= ba;
              } else {
                u2.re /= ba;
                u2.im /= ba;
              }

              tst = h_data[m + h_size[0] * (m - 1)].re;
              htmp1 = h_data[m + h_size[0] * (m - 1)].im;
              h_data[m + h_size[0] * (m - 1)].re = tst * u2.re - htmp1 * -u2.im;
              h_data[m + h_size[0] * (m - 1)].im = tst * -u2.im + htmp1 * u2.re;
              if (m + 2 <= i + 1) {
                tst = h_data[3 + h_size[0] * m].re;
                htmp1 = h_data[3 + h_size[0] * m].im;
                h_data[3 + h_size[0] * m].re = tst * u2.re - htmp1 * u2.im;
                h_data[3 + h_size[0] * m].im = tst * u2.im + htmp1 * u2.re;
              }

              for (j = m; j <= i + 1; j++) {
                if (j != m + 1) {
                  if (n > j) {
                    b_xscal(n - j, u2, h_data, j + j * ldh, ldh);
                  }

                  y.re = u2.re;
                  y.im = -u2.im;
                  xscal(j - 1, y, h_data, 1 + (j - 1) * ldh);
                }
              }
            }
          }

          u2 = h_data[i + h_size[0] * (i - 1)];
          if (h_data[i + h_size[0] * (i - 1)].im != 0.0) {
            tst = rt_hypotd_snf(h_data[i + h_size[0] * (i - 1)].re, h_data[i +
                                h_size[0] * (i - 1)].im);
            h_data[i + h_size[0] * (i - 1)].re = tst;
            h_data[i + h_size[0] * (i - 1)].im = 0.0;
            if (u2.im == 0.0) {
              u2.re /= tst;
              u2.im = 0.0;
            } else if (u2.re == 0.0) {
              u2.re = 0.0;
              u2.im /= tst;
            } else {
              u2.re /= tst;
              u2.im /= tst;
            }

            if (n > i + 1) {
              y.re = u2.re;
              y.im = -u2.im;
              b_xscal((n - i) - 1, y, h_data, (i + (i + 1) * ldh) + 1, ldh);
            }

            xscal(i, u2, h_data, 1 + i * ldh);
          }

          its++;
        }
      }

      if (!goto140) {
        info = i + 1;
        exitg1 = true;
      } else {
        i = L;
      }
    }
  }

  return info;
}

/* End of code generation (xzhseqr.c) */
