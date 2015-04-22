#ifndef _POLY_GAMMA_H_INCLUDED_
#define _POLY_GAMMA_H_INCLUDED_

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <mpfr_optimization.h>
#include <poly_eval.h>
#include <mpfr_cdf.h>

#include <count_temp_decl.h>

typedef struct {
	double theta;
	double shape;
	double scale;
} poly_gamma_proposal_param_t;

void mpfr_poly_gamma_gradient_log(mpfr_ptr grad, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr shape, mpfr_ptr scale);
NTEMP_DECL(mpfr_poly_gamma_gradient_log);
void mpfr_poly_gamma_gradient_log_v(mpfr_ptr grad, mpfr_ptr x, va_list args);
NTEMP_DECL(mpfr_poly_gamma_gradient_log_v);

void mpfr_poly_gamma_hessian_log(mpfr_ptr hess, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr shape, mpfr_ptr scale);
NTEMP_DECL(mpfr_poly_gamma_hessian_log);
void mpfr_poly_gamma_hessian_log_v(mpfr_ptr grad, mpfr_ptr x, va_list args);
NTEMP_DECL(mpfr_poly_gamma_hessian_log_v);

void mpfr_poly_gamma_third_deriv_log(mpfr_ptr third_deriv, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr shape, mpfr_ptr scale);
NTEMP_DECL(mpfr_poly_gamma_third_deriv_log);
void mpfr_poly_gamma_third_deriv_log_v(mpfr_ptr third_deriv, mpfr_ptr x, va_list args);
NTEMP_DECL(mpfr_poly_gamma_third_deriv_log_v);

poly_gamma_proposal_param_t
		poly_gamma_proposal_param_max_using_mpfr(size_t nterms, double *terms, double constr, double shape, double scale);
NTEMP_DECL(poly_gamma_proposal_param_max_using_mpfr);

void mpfr_poly_gamma_proposal_param_max(mpfr_ptr theta_par, mpfr_ptr shape_par, mpfr_ptr scale_par, size_t nterms, mpfr_t* coef, mpfr_ptr constr, mpfr_ptr shape, mpfr_ptr scale);
NTEMP_DECL(mpfr_poly_gamma_proposal_param_max);

double poly_gamma_pdf_ratio_using_mpfr(double x1, double x2, size_t nterms, double* terms, double constr, double shape, double scale);
NTEMP_DECL(poly_gamma_pdf_ratio_using_mpfr);
double poly_gamma_pdf_using_mpfr(double x, size_t nterms, double* terms, double constr, double shape, double scale);
NTEMP_DECL(poly_gamma_pdf_using_mpfr);

void mpfr_poly_gamma_pdf(mpfr_ptr pdf, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr constr, mpfr_ptr shape, mpfr_ptr scale);
NTEMP_DECL(mpfr_poly_gamma_pdf);

void mpfr_poly_gamma_pdf_pass_temp(mpfr_ptr pdf, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr constr, mpfr_ptr shape, mpfr_ptr scale, mpfr_t *temp);
NTEMP_DECL(mpfr_poly_gamma_pdf_pass_temp);

double poly_gamma_proposal_pdf_using_mpfr(double x, poly_gamma_proposal_param_t par);
NTEMP_DECL(poly_gamma_proposal_pdf_using_mpfr);

/******* REGULAR FUNCTIONS *******/
double poly_gamma_proposal_samp(gsl_rng * r, poly_gamma_proposal_param_t par);
double poly_gamma_proposal_pdf(double x, poly_gamma_proposal_param_t par);
double poly_gamma_pdf(double x, double constr, size_t nterms, double* terms, double shape, double scale);

#endif
