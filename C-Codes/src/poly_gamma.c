#include <poly_gamma.h>

#include <count_temp_start.h>

/**** MPFR FUNCTIONS ****/

void mpfr_poly_gamma_gradient_log(mpfr_ptr grad, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr shape, mpfr_ptr scale) {
	size_t ntemp = NTEMP(mpfr_poly_gamma_gradient_log);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr sum = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(set_ui, sum, 0);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(add, tmp, x, terms[i]);
		MPFR_FUN(ui_div, tmp, 1, tmp);
		MPFR_FUN(add, sum, sum, tmp);
	}
	MPFR_FUN(sub_ui, tmp, shape, 1);
	MPFR_FUN(div, tmp, tmp, x);
	MPFR_FUN(add, sum, sum, tmp);
	MPFR_FUN(ui_div, tmp, 1, scale);
	MPFR_FUN(sub, sum, sum, tmp);
	MPFR_FUN(set, grad, sum);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_poly_gamma_gradient_log
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gamma_gradient_log

void mpfr_poly_gamma_gradient_log_v(mpfr_ptr grad, mpfr_ptr x, va_list args) {
	size_t nterms = va_arg(args, size_t);
	mpfr_t *terms = va_arg(args, mpfr_t*);
	mpfr_ptr shape = va_arg(args, mpfr_ptr);
	mpfr_ptr scale = va_arg(args, mpfr_ptr);
	mpfr_poly_gamma_gradient_log(grad, x, nterms, terms, shape, scale);
	return;
}
#define FUN_NAME mpfr_poly_gamma_gradient_log_v
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gamma_gradient_log_v

void mpfr_poly_gamma_hessian_log(mpfr_ptr hess, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr shape, mpfr_ptr scale) {
	size_t ntemp = NTEMP(mpfr_poly_gamma_hessian_log);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr sum = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(set_ui, sum, 0);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(add, tmp, x, terms[i]);
		MPFR_FUN(pow_ui, tmp, tmp, 2);
		MPFR_FUN(ui_div, tmp, 1, tmp);
		MPFR_FUN(sub, sum, sum, tmp);
	}
	MPFR_FUN(sub_ui, tmp, shape, 1);
	MPFR_FUN(div, tmp, tmp, x);
	MPFR_FUN(div, tmp, tmp, x);
	MPFR_FUN(sub, sum, sum, tmp);
	MPFR_FUN(set, hess, sum);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_poly_gamma_hessian_log
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gamma_hessian_log

void mpfr_poly_gamma_hessian_log_v(mpfr_ptr grad, mpfr_ptr x, va_list args) {
	size_t nterms = va_arg(args, size_t);
	mpfr_t *terms = va_arg(args, mpfr_t*);
	mpfr_ptr shape = va_arg(args, mpfr_ptr);
	mpfr_ptr scale = va_arg(args, mpfr_ptr);
	mpfr_poly_gamma_hessian_log(grad, x, nterms, terms, shape, scale);
	return;
}
#define FUN_NAME mpfr_poly_gamma_hessian_log_v
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gamma_hessian_log_v

void mpfr_poly_gamma_third_deriv_log(mpfr_ptr third_deriv, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr shape, mpfr_ptr scale) {
	size_t ntemp = NTEMP(mpfr_poly_gamma_third_deriv_log);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr sum = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(set_ui, sum, 0);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(add, tmp, x, terms[i]);
		MPFR_FUN(pow_ui, tmp, tmp, 3);
		MPFR_FUN(ui_div, tmp, 2, tmp);
		MPFR_FUN(add, sum, sum, tmp);
	}
	MPFR_FUN(sub_ui, tmp, shape, 1);
	MPFR_FUN(div, tmp, tmp, x);
	MPFR_FUN(div, tmp, tmp, x);
	MPFR_FUN(div, tmp, tmp, x);
	MPFR_FUN(mul_ui, tmp, tmp, 2);
	MPFR_FUN(add, sum, sum, tmp);
	MPFR_FUN(set, third_deriv, sum);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_poly_gamma_third_deriv_log
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gamma_third_deriv_log

void mpfr_poly_gamma_third_deriv_log_v(mpfr_ptr third_deriv, mpfr_ptr x, va_list args) {
	size_t nterms = va_arg(args, size_t);
	mpfr_t *terms = va_arg(args, mpfr_t*);
	mpfr_ptr shape = va_arg(args, mpfr_ptr);
	mpfr_ptr scale = va_arg(args, mpfr_ptr);
	mpfr_poly_gamma_third_deriv_log(third_deriv, x, nterms, terms, shape, scale);
	return;
}
#define FUN_NAME mpfr_poly_gamma_third_deriv_log_v
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gamma_third_deriv_log_v

poly_gamma_proposal_param_t
		poly_gamma_proposal_param_max_using_mpfr(size_t nterms, double *terms, double constr, double shape, double scale) {
	size_t ntemp = NTEMP(poly_gamma_proposal_param_max_using_mpfr);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr constr_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr shape_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr scale_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr theta_par = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr shape_par = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr scale_par = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *terms_m = mpfr_array_alloc_init(nterms);
	MPFR_FUN(set_d, constr_m, constr);
	MPFR_FUN(set_d, shape_m, shape);
	MPFR_FUN(set_d, scale_m, scale);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(set_d, terms_m[i], terms[i]);
	}
	mpfr_poly_gamma_proposal_param_max(theta_par, shape_par, scale_par, nterms, terms_m, constr_m, shape_m, scale_m);
	poly_gamma_proposal_param_t par;
	par.theta = MPFR_FUN(get_d, theta_par);
	par.shape = MPFR_FUN(get_d, shape_par);
	par.scale = MPFR_FUN(get_d, scale_par);
	mpfr_array_clear_free(terms_m, nterms);
	mpfr_array_clear_free(temp, ntemp);
	return par;
}
#define FUN_NAME poly_gamma_proposal_param_max_using_mpfr
#include <count_temp.h>
#define FUN_NAME1 poly_gamma_proposal_param_max_using_mpfr

void mpfr_poly_gamma_proposal_param_max(mpfr_ptr theta_par, mpfr_ptr shape_par, mpfr_ptr scale_par, size_t nterms, mpfr_t* coef, mpfr_ptr constr, mpfr_ptr shape, mpfr_ptr scale) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_gamma_pdf_pass_temp)
	if (mpfr_sgn(constr)>0) {
		MPFR_FUN(set, theta_par, constr);
	} else {
		MPFR_FUN(set_ui, theta_par, 0);
	}
	MPFR_FUN(set, scale_par, scale);
	if (nterms==0) {
		MPFR_FUN(set, shape_par, shape);
		return;
	}
	size_t ntemp = NTEMP(mpfr_poly_gamma_proposal_param_max);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	//MPFR_PRINT_VECTOR("coef", coef, nterms);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr tmp1 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr tmp2 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr xmax = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr init = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr c1 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr c2 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr term_max = coef[0];
	mpfr_ptr term_min = coef[0];
	for (size_t i=1; i<nterms; i++) {
		mpfr_ptr term = coef[i];
		if (mpfr_cmp(term,term_max)>0) {
			term_max = term;
		}
		if (mpfr_cmp(term,term_min)<0) {
			term_min = term;
		}
	}
	
	/* computing upper bound */
	MPFR_FUN(add_ui, tmp1, shape, nterms-1);
	MPFR_FUN(mul, tmp1, tmp1, scale);
	MPFR_FUN(sub, tmp1, tmp1, term_min); // scale*(shape+nterms-1)-term_min
	
	MPFR_FUN(pow_ui, tmp2, tmp1, 2); // the same squared
	
	MPFR_FUN(sub_ui, tmp, shape, 1);
	MPFR_FUN(mul, tmp, tmp, scale);
	MPFR_FUN(mul, tmp, tmp, term_min);
	MPFR_FUN(mul_ui, tmp, tmp, 4); // 4*scale*(shape-1)*term_min
	
	MPFR_FUN(add, tmp, tmp, tmp2);
	
	if (mpfr_sgn(tmp)<0) { // no max at all
		MPFR_FUN(set, xmax, theta_par);
	} else {
		MPFR_FUN(sqrt, tmp, tmp);
		MPFR_FUN(add, tmp, tmp, tmp1);
		MPFR_FUN(div_ui, c2, tmp, 2); // upper bound
		if (mpfr_cmp(c2, theta_par)<=0) { // no max in the region
			MPFR_FUN(set, xmax, theta_par);
		} else {
			/* computing lower bound */
			MPFR_FUN(add_ui, tmp1, shape, nterms-1);
			MPFR_FUN(mul, tmp1, tmp1, scale);
			MPFR_FUN(sub, tmp1, tmp1, term_max); // scale*(shape+nterms-1)-term_max
			
			MPFR_FUN(pow_ui, tmp2, tmp1, 2); // the same squared
	
			MPFR_FUN(sub_ui, tmp, shape, 1);
			MPFR_FUN(mul, tmp, tmp, scale);
			MPFR_FUN(mul, tmp, tmp, term_max);
			MPFR_FUN(mul_ui, tmp, tmp, 4); // 4*scale*(shape-1)*term_max
			
			MPFR_FUN(add, tmp, tmp, tmp2);
			if (mpfr_sgn(tmp)<0) { // no lower bound
				MPFR_FUN(set, c1, theta_par);
			} else {
				MPFR_FUN(sqrt, tmp, tmp);
				MPFR_FUN(add, tmp, tmp, tmp1);
				MPFR_FUN(div_ui, c1, tmp, 2);
				if (mpfr_cmp(c1, theta_par)<0) {
					MPFR_FUN(set, c1, theta_par);
				}
			}
	
			MPFR_FUN(add, init, c1, c2);
			MPFR_FUN(div_ui, init, init, 2);
			mpfr_optimization_newton(xmax, init, c1, c2, mpfr_poly_gamma_gradient_log_v, 
					mpfr_poly_gamma_hessian_log_v, nterms, coef, shape, scale);
		}
	}
	MPFR_FUN(sub, tmp, xmax, theta_par);
	MPFR_FUN(div, tmp, tmp, scale);
	MPFR_FUN(add_ui, shape_par, tmp, 1);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_poly_gamma_proposal_param_max
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gamma_proposal_param_max

double poly_gamma_pdf_ratio_using_mpfr(double x1, double x2, size_t nterms, double* terms, double constr, double shape, double scale) {
	size_t ntemp = NTEMP(poly_gamma_pdf_ratio_using_mpfr);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr pdf1_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr pdf2_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr pdf_ratio_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr x1_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr x2_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr constr_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr shape_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr scale_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *terms_m = mpfr_array_alloc_init(nterms);
	MPFR_FUN(set_d, x1_m, x1);
	MPFR_FUN(set_d, x2_m, x2);
	MPFR_FUN(set_d, constr_m, constr);
	MPFR_FUN(set_d, shape_m, shape);
	MPFR_FUN(set_d, scale_m, scale);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(set_d, terms_m[i], terms[i]);
	}
	mpfr_poly_gamma_pdf(pdf1_m, x1_m, nterms, terms_m, constr_m, shape_m, scale_m);
	mpfr_poly_gamma_pdf(pdf2_m, x2_m, nterms, terms_m, constr_m, shape_m, scale_m);
	MPFR_FUN(div, pdf_ratio_m, pdf1_m, pdf2_m);
	double ratio = MPFR_FUN(get_d, pdf_ratio_m);
	mpfr_array_clear_free(temp, ntemp);
	mpfr_array_clear_free(terms_m, nterms);
	return ratio;
}
#define FUN_NAME poly_gamma_pdf_ratio_using_mpfr
#include <count_temp.h>
#define FUN_NAME1 poly_gamma_pdf_ratio_using_mpfr

double poly_gamma_pdf_using_mpfr(double x, size_t nterms, double* terms, double constr, double shape, double scale) {
	size_t ntemp = NTEMP(poly_gamma_pdf_using_mpfr);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr pdf_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr x_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr constr_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr shape_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr scale_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *terms_m = mpfr_array_alloc_init(nterms);
	MPFR_FUN(set_d, x_m, x);
	MPFR_FUN(set_d, constr_m, constr);
	MPFR_FUN(set_d, shape_m, shape);
	MPFR_FUN(set_d, scale_m, scale);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(set_d, terms_m[i], terms[i]);
	}
	mpfr_poly_gamma_pdf(pdf_m, x_m, nterms, terms_m, constr_m, shape_m, scale_m);
	double pdf = MPFR_FUN(get_d, pdf_m);
	mpfr_array_clear_free(temp, ntemp);
	mpfr_array_clear_free(terms_m, nterms);
	return pdf;
}
#define FUN_NAME poly_gamma_pdf_using_mpfr
#include <count_temp.h>
#define FUN_NAME1 poly_gamma_pdf_using_mpfr

void mpfr_poly_gamma_pdf(mpfr_ptr pdf, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr constr, mpfr_ptr shape, mpfr_ptr scale) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_gamma_pdf_pass_temp)
	size_t ntemp = NTEMP(mpfr_poly_gamma_pdf);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_poly_gamma_pdf_pass_temp(pdf, x, nterms, terms, constr, shape, scale, pass_temp);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_poly_gamma_pdf
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gamma_pdf

void mpfr_poly_gamma_pdf_pass_temp(mpfr_ptr pdf, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr constr, mpfr_ptr shape, mpfr_ptr scale, mpfr_t *temp) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_eval_mult_pass_temp)
	if (mpfr_sgn(x)<0 || mpfr_cmp(x,constr)<0 || mpfr_inf_p(x)) {
		MPFR_FUN(set_ui, pdf, 0);
		return;
	}
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr val = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_poly_eval_mult_pass_temp(val, x, nterms, terms, pass_temp);
	mpfr_gamma_pdf(tmp, x, shape, scale);
	MPFR_FUN(mul, val, val, tmp);
	MPFR_FUN(set, pdf, val);
	return;
}
#define FUN_NAME mpfr_poly_gamma_pdf_pass_temp
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gamma_pdf_pass_temp

double poly_gamma_proposal_pdf_using_mpfr(double x, poly_gamma_proposal_param_t par) {
	size_t ntemp = NTEMP(poly_gamma_proposal_pdf_using_mpfr);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr pdf_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr x_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr theta = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr shape = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr scale = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(set_d, x_m, x);
	MPFR_FUN(set_d, theta, par.theta);
	MPFR_FUN(set_d, shape, par.shape);
	MPFR_FUN(set_d, scale, par.scale);
	MPFR_FUN(sub, x_m, x_m, theta);
	mpfr_gamma_pdf(pdf_m, x_m, shape, scale);
	double pdf = MPFR_FUN(get_d, pdf_m);
	mpfr_array_clear_free(temp, ntemp);
	return pdf;
}
#define FUN_NAME poly_gamma_proposal_pdf_using_mpfr
#include <count_temp.h>
#define FUN_NAME1 poly_gamma_proposal_pdf_using_mpfr

/**** Regular functions ****/

double poly_gamma_proposal_samp(gsl_rng * r, poly_gamma_proposal_param_t par) {
	double samp = par.theta+gsl_ran_gamma_mt(r, par.shape, par.scale);
	return samp;
}

double poly_gamma_proposal_pdf(double x, poly_gamma_proposal_param_t par) {
	double dist = gsl_ran_gamma_pdf(x-par.theta,par.shape,par.scale);
	return dist;
}

double poly_gamma_pdf(double x, double constr, size_t nterms, double* terms, double shape, double scale) {
	return (x>constr) ? poly_eval_mult_regular(x,nterms,terms) * gsl_ran_gamma_pdf(x,shape,scale) : 0;
}
