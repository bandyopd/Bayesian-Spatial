#ifndef _SAMP_ALL_INCLUDED_
#define _SAMP_ALL_INCLUDED_

#include <stdarg.h>

#include <mcmc.h>
#include <poly_gamma.h>
#include <poly_gaussian.h>
#include <poly_exp.h>

double samp_all_car(size_t iter, size_t pos, size_t nparam, double* param, gsl_rng* ran_gen, size_t nobs, size_t nfrail, 
			size_t ncov, size_t ntimes, double* time_int, size_t* frail_tbl, size_t* nevents, 
			size_t** event_obs, double* atrisk_propsum, double *cov, double* cov_lim, double r0, double c0, 
			size_t* nevcov, size_t** evcov_obs, double* atrisk_sumcov,
			size_t *nevfrail, size_t** evfrail_obs, double *atrisk_sumprop_frail, size_t *nadj, size_t **adj_ind,
			double theta_shape, double theta_scale, double var_proportion);

double samp_all_car_v(size_t iter, size_t pos, size_t nparam, double* param, gsl_rng* ran_gen, va_list args);



double baseline_fullcond_samp(size_t iter, size_t pos, double current, size_t ncov, size_t nfrail, size_t ntimes, size_t nevents, 
			size_t *evobs, double dt, double psum, double *regcoef, double *frail, size_t *frail_tbl, 
			double (*cov)[ncov][ntimes], double (*cov_lim)[2], double r0, double c0, gsl_rng *ran_gen);
		
double regfun_fullcond_samp(size_t iter, size_t cov_pos, size_t time_pos, double current, 
			size_t ncov, size_t nfrail, size_t ntimes, size_t nevcov,
			size_t *evcov_obs, double dt, double sum_cv_ar, double baseline, double *regcoef, 
			double *frail, size_t *frail_tbl, double (*cov)[ncov][ntimes], double (*covlim)[2],
			gsl_rng *ran_gen);
		
double frail_fullcond_samp(size_t iter, size_t lpos, size_t jpos, double current, size_t ncov, size_t nfrail, size_t ntimes, 
			size_t nevfrail, size_t *evfrail_obs, double dt,
			double psumfr, size_t nadj, size_t *adj_ind, 
			double baseline, double *regcoef, double *frail, size_t *frail_tbl, double (*cov)[ncov][ntimes], double (*cov_lim)[2], 
			double theta, double var_proportion, gsl_rng *ran_gen);
				
double theta_fullcond_samp(size_t iter, size_t jpos, double theta, size_t nfrail, size_t ntimes, double *frail, 
			size_t *nadj, size_t **adj_ind, double theta_shape, double theta_scale, double var_proportion, gsl_rng *ran_gen);

double baslin_samp(double current, gsl_rng * ran_gen, poly_gamma_proposal_param_t par);
double baslin_samp_v(double current, gsl_rng * ran_gen, va_list args);
double baslin_samp_args(double current, gsl_rng * ran_gen, void** args);
double baslin_samp1_v(double current, gsl_rng * ran_gen, va_list args);
double baslin_pdf_ratio(double x1, double x2, double constr, size_t nterms, double *terms, double shape, double scale);
double baslin_pdf_ratio_v(double x1, double x2, va_list args);
double baslin_pdf(double x, double constr, size_t nterms, double *terms, double shape, double scale);
double baslin_pdf_v(double x, va_list args);
double baslin_pdf1_v(double x, va_list args);
double baslin_pdf_args(double x, void** args);
double baslin_proposal_pdf(double x, double x_old, poly_gamma_proposal_param_t par);
double baslin_proposal_pdf_v(double x, double x_old, va_list args);
double baslin_proposal_pdf1_v(double x, double x_old, va_list args);
double baslin_proposal_pdf_args(double x, double x_old, void** args);

double regfun_gaussian_samp(double current, gsl_rng* ran_gen, poly_gaussian_proposal_param_t par);
double regfun_gaussian_samp_v(double current, gsl_rng* ran_gen, va_list pass_args);
double regfun_gaussian_pdf_ratio(double x1, double x2, size_t nterms, double *terms, double c1, double c2, double mean, double stdev);
double regfun_gaussian_pdf_ratio_v(double x1, double x2, va_list pass_args);
double regfun_gaussian_pdf(double x, size_t nterms, double *terms, double c1, double c2, double mean, double stdev);
double regfun_gaussian_pdf_v(double x, va_list pass_args);
double regfun_gaussian_proposal_pdf(double x, double x_old, poly_gaussian_proposal_param_t par);
double regfun_gaussian_proposal_pdf_v(double x, double x_old, va_list pass_args);

double regfun_exp_samp_direct(gsl_rng* ran_gen, double constr1, double constr2, poly_exp_proposal_param_t par);
double regfun_exp_samp(double current, gsl_rng* ran_gen, poly_exp_proposal_param_t par);
double regfun_exp_samp_v(double current, gsl_rng* ran_gen, va_list pass_args);
double regfun_exp_pdf_ratio(double x1, double x2, size_t nterms, double *terms, double c1, double c2, double beta);
double regfun_exp_pdf_ratio_v(double x1, double x2, va_list pass_args);
double regfun_exp_pdf(double x, size_t nterms, double *terms, double c1, double c2, double beta);
double regfun_exp_pdf_v(double x, va_list pass_args);
double regfun_exp_proposal_pdf(double x, double x_old, poly_exp_proposal_param_t par);
double regfun_exp_proposal_pdf_v(double x, double x_old, va_list pass_args);


#endif
