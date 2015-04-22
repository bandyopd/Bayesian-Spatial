#ifndef _POLY_GAUSSIAN_H_INCLUDED_
#define _POLY_GAUSSIAN_H_INCLUDED_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <mpfr_optimization.h>
#include <poly_eval.h>
#include <mpfr_cdf.h>

#include <count_temp_decl.h>

typedef struct {
	double mean;
	double stdev;
} poly_gaussian_proposal_param_t;

void mpfr_poly_gaussian_gradient_log(mpfr_ptr grad, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr mean, mpfr_ptr stdev);
NTEMP_DECL(mpfr_poly_gaussian_gradient_log);
void mpfr_poly_gaussian_gradient_log_v(mpfr_ptr grad, mpfr_ptr x, va_list pass_args);
NTEMP_DECL(mpfr_poly_gaussian_gradient_log_v);

void mpfr_poly_gaussian_hessian_log(mpfr_ptr hess, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr mean, mpfr_ptr stdev);
NTEMP_DECL(mpfr_poly_gaussian_hessian_log);
void mpfr_poly_gaussian_hessian_log_v(mpfr_ptr grad, mpfr_ptr x, va_list pass_args);
NTEMP_DECL(mpfr_poly_gaussian_hessian_log_v);

poly_gaussian_proposal_param_t 
	poly_gaussian_proposal_param_max_using_mpfr(size_t nterms, double* coef_mult, double constr1, double constr2, double mean, double stdev);
NTEMP_DECL(poly_gaussian_proposal_param_max_using_mpfr);

void mpfr_poly_gaussian_proposal_param_max(mpfr_ptr mean_par, mpfr_ptr stdev_par, size_t nterms, mpfr_t* coef, mpfr_ptr constr1, mpfr_ptr constr2, mpfr_ptr mean, mpfr_ptr stdev);
NTEMP_DECL(mpfr_poly_gaussian_proposal_param_max);

void mpfr_poly_gaussian_pdf_pass_temp(mpfr_t dist, mpfr_t x, size_t nterms, mpfr_t *coef, mpfr_t c1, mpfr_t c2, mpfr_t mean, mpfr_t stdev, mpfr_t* temp);
NTEMP_DECL(mpfr_poly_gaussian_pdf_pass_temp);
double poly_gaussian_pdf_using_mpfr(double x, size_t nterms, double *terms, double c1, double c2, double mean, double stdev);
NTEMP_DECL(poly_gaussian_pdf_using_mpfr);
double poly_gaussian_pdf_ratio_using_mpfr(double x1, double x2, size_t nterms, double *terms, double c1, double c2, double mean, double stdev);
NTEMP_DECL(poly_gaussian_pdf_ratio_using_mpfr);
double poly_gaussian_proposal_pdf(double x, poly_gaussian_proposal_param_t par);
double poly_gaussian_proposal_samp(poly_gaussian_proposal_param_t par, gsl_rng* ran_gen);

#endif
