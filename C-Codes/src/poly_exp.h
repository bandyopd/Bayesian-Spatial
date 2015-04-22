#ifndef _POLY_EXP_H_INCLUDED_
#define _POLY_EXP_H_INCLUDED_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <mpfr_optimization.h>
#include <poly_eval.h>

#include <count_temp_decl.h>

typedef struct {
	char gamma;
	double offset;
	int sign;
	double shape;
	double scale;
	double mean;
	double stdev;
} poly_exp_proposal_param_t;

void mpfr_poly_exp_gradient_log(mpfr_ptr grad, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr beta);
NTEMP_DECL(mpfr_poly_exp_gradient_log);
void mpfr_poly_exp_gradient_log_v(mpfr_ptr grad, mpfr_ptr x, va_list args);
NTEMP_DECL(mpfr_poly_exp_gradient_log_v);

void mpfr_poly_exp_hessian_log(mpfr_ptr hess, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr beta);
NTEMP_DECL(mpfr_poly_exp_hessian_log);
void mpfr_poly_exp_hessian_log_v(mpfr_ptr grad, mpfr_ptr x, va_list args);
NTEMP_DECL(mpfr_poly_exp_hessian_log_v);

poly_exp_proposal_param_t
	poly_exp_proposal_param_max_using_mpfr(size_t nterms, double* coef, double c1, double c2, double beta);
NTEMP_DECL(poly_exp_proposal_param_max_using_mpfr);

void mpfr_poly_exp_gaussian_proposal_param_max(mpfr_ptr mean_par, mpfr_ptr stdev_par, size_t nterms, mpfr_t* coef, mpfr_ptr constr1, mpfr_ptr constr2, mpfr_ptr beta);
NTEMP_DECL(mpfr_poly_exp_gaussian_proposal_param_max);
void mpfr_poly_exp_gamma_proposal_param_max(mpfr_ptr offset_par, mpfr_ptr shape_par, mpfr_ptr scale_par, int *sign, size_t nterms, mpfr_t* coef, mpfr_ptr c1, mpfr_ptr c2, mpfr_ptr beta);
NTEMP_DECL(mpfr_poly_exp_gamma_proposal_param_max);

void mpfr_poly_exp_find_max(mpfr_ptr xmax, size_t nterms, mpfr_t *coef, mpfr_ptr c1, mpfr_ptr c2, mpfr_ptr beta);
NTEMP_DECL(mpfr_poly_exp_find_max);

void mpfr_poly_exp_pdf(mpfr_ptr dist, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr c1, mpfr_ptr c2, mpfr_ptr beta);
NTEMP_DECL(mpfr_poly_exp_pdf);
void mpfr_poly_exp_pdf_pass_temp(mpfr_ptr dist, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr c1, mpfr_ptr c2, mpfr_ptr beta, mpfr_t *temp);
NTEMP_DECL(mpfr_poly_exp_pdf_pass_temp);

double poly_exp_pdf_ratio_using_mpfr(double x1, double x2, size_t nterms, double *terms, double c1, double c2, double beta);
NTEMP_DECL(poly_exp_pdf_ratio_using_mpfr);
double poly_exp_pdf_using_mpfr(double x, size_t nterms, double *terms, double c1, double c2, double beta);
NTEMP_DECL(poly_exp_pdf_using_mpfr);
double poly_exp_proposal_pdf(double x, poly_exp_proposal_param_t par);
double poly_exp_proposal_samp(poly_exp_proposal_param_t par, gsl_rng* ran_gen);

#endif
