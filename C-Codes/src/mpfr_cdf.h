#ifndef _MPFR_CDF_H_INCLUDED_
#define _MPFR_CDF_H_INCLUDED_

#include <stdio.h>

#include <gsl/gsl_machine.h>

#include <mpfr_array.h>

#include <count_temp_decl.h>


#define MPFR_CDF_MAX_ITER 50000
/**** MPF FUNCTIONS ****/
void mpf_gamma_pdf(mpf_ptr pdf, mpf_ptr x, mpf_ptr a, mpf_ptr b);
NTEMP_DECL(mpf_gamma_pdf);

void mpf_gamma_cdf_Q(mpf_t Q, mpf_t x, mpf_t a, mpf_t b);
NTEMP_DECL(mpf_gamma_cdf_Q);

/**** MPFR FUNCTIONS ****/
void mpfr_gamma_pdf(mpfr_ptr pdf, mpfr_ptr x, mpfr_ptr a, mpfr_ptr b);
NTEMP_DECL(mpfr_gamma_pdf);

void mpfr_gamma_cdf_Q(mpfr_t Q, mpfr_t x, mpfr_t a, mpfr_t b);
NTEMP_DECL(mpfr_gamma_cdf_Q);

void mpfr_pdf_gaussian(mpfr_ptr pdf, mpfr_ptr x, mpfr_ptr mean, mpfr_ptr stdev);
NTEMP_DECL(mpfr_pdf_gaussian);
void mpfr_pdf_ugaussian(mpfr_ptr pdf, mpfr_ptr x);
NTEMP_DECL(mpfr_pdf_ugaussian);
void mpfr_cdf_gaussian_P(mpfr_ptr P, mpfr_ptr x, mpfr_ptr mean, mpfr_ptr stdev);
NTEMP_DECL(mpfr_cdf_gaussian_P);
void mpfr_cdf_ugaussian_P(mpfr_ptr P, mpfr_ptr x);
NTEMP_DECL(mpfr_cdf_ugaussian_P);
void mpfr_cdf_ugaussian_Q(mpfr_ptr Q, mpfr_ptr x);
NTEMP_DECL(mpfr_cdf_ugaussian_Q);

void mpfr_gamma_cdf_Q_routine(mpfr_t Q, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp);
NTEMP_DECL(mpfr_gamma_cdf_Q_routine);

void mpfr_gamma_cdf_D(mpfr_ptr D, mpfr_ptr x, mpfr_ptr a, mpfr_t *temp);
NTEMP_DECL(mpfr_gamma_cdf_D);

void mpfr_gamma_cdf_F_CF(mpfr_t F_CF, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp);
NTEMP_DECL(mpfr_gamma_cdf_F_CF);

void mpfr_gamma_cdf_Q_CF(mpfr_t Q_CF, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp);
NTEMP_DECL(mpfr_gamma_cdf_Q_CF);

void mpfr_gamma_cdf_P_series(mpfr_t P, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp);
NTEMP_DECL(mpfr_gamma_cdf_P_series);

void mpfr_gamma_cdf_Q_series(mpfr_t Q, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp);
NTEMP_DECL(mpfr_gamma_cdf_Q_series);

void mpfr_gamma_cdf_Q_asymp_unif(mpfr_t Q, mpfr_t x, mpfr_t a, mpfr_t eps, mpfr_t *temp);
NTEMP_DECL(mpfr_gamma_cdf_Q_asymp_unif);

void mpfr_exprel_n_CF(mpfr_t CF, mpfr_t x, mpfr_t N, mpfr_t eps, mpfr_t *temp);
NTEMP_DECL(mpfr_exprel_n_CF);

void mpfr_erf_approx(mpfr_ptr efr, mpfr_ptr x);
NTEMP_DECL(mpfr_erf_approx);

#endif