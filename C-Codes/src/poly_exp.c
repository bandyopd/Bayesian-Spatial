#include <poly_exp.h>

#include <count_temp_start.h>

void mpfr_poly_exp_gradient_log(mpfr_ptr grad, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr beta) {
	size_t ntemp = NTEMP(mpfr_poly_exp_gradient_log);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr sum = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(neg, sum, beta);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(add, tmp, x, terms[i]);
		MPFR_FUN(ui_div, tmp, 1, tmp);
		MPFR_FUN(add, sum, sum, tmp);
	}
	MPFR_FUN(set, grad, sum);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_poly_exp_gradient_log
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_exp_gradient_log

void mpfr_poly_exp_gradient_log_v(mpfr_ptr grad, mpfr_ptr x, va_list args) {
	size_t nterms = va_arg(args, size_t);
	mpfr_t *terms = va_arg(args, mpfr_t*);
	mpfr_ptr beta = va_arg(args, mpfr_ptr);
	mpfr_poly_exp_gradient_log(grad, x, nterms, terms, beta);
	return;
}
#define FUN_NAME mpfr_poly_exp_gradient_log_v
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_exp_gradient_log_v

void mpfr_poly_exp_hessian_log(mpfr_ptr hess, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr beta) {
	size_t ntemp = NTEMP(mpfr_poly_exp_hessian_log);
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
	MPFR_FUN(set, hess, sum);
	mpfr_array_clear_free(temp, ntemp);
}
#define FUN_NAME mpfr_poly_exp_hessian_log
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_exp_hessian_log

void mpfr_poly_exp_hessian_log_v(mpfr_ptr grad, mpfr_ptr x, va_list args) {
	size_t nterms = va_arg(args, size_t);
	mpfr_t *terms = va_arg(args, mpfr_t*);
	mpfr_ptr beta = va_arg(args, mpfr_ptr);
	mpfr_poly_exp_hessian_log(grad, x, nterms, terms, beta);
	return;
}
#define FUN_NAME mpfr_poly_exp_hessian_log_v
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_exp_hessian_log_v

poly_exp_proposal_param_t
			poly_exp_proposal_param_max_using_mpfr(size_t nterms, double* coef, double c1, double c2, double beta) {
	size_t ntemp = NTEMP(poly_exp_proposal_param_max_using_mpfr);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr c1_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr c2_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr beta_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *coef_m = mpfr_array_alloc_init(nterms);
	
	MPFR_FUN(set_d, c1_m, c1);
	MPFR_FUN(set_d, c2_m, c2);
	MPFR_FUN(set_d, beta_m, beta);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(set_d, coef_m[i], coef[i]);
	}
	
	mpfr_ptr first_par = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr second_par = GET_NEXT_TEMP_PTR(temp);
	
	poly_exp_proposal_param_t par;
	if (beta==0.0) {
		mpfr_ptr mean_par = first_par;
		mpfr_ptr stdev_par = second_par;
		
		mpfr_poly_exp_gaussian_proposal_param_max(mean_par, stdev_par, nterms, coef_m, c1_m, c2_m, beta_m);
		
		par.gamma = 0;
		par.mean = MPFR_FUN(get_d, mean_par);
		par.stdev = MPFR_FUN(get_d, stdev_par);
	} else {
		mpfr_ptr shape_par = first_par;
		mpfr_ptr scale_par = second_par;
		mpfr_ptr offset_par = GET_NEXT_TEMP_PTR(temp);
		int sign_par;
		
		mpfr_poly_exp_gamma_proposal_param_max(offset_par, shape_par, scale_par, &sign_par, nterms, coef_m, c1_m, c2_m, beta_m);
		
		par.gamma = 1;
		par.offset = MPFR_FUN(get_d, offset_par);
		par.shape = MPFR_FUN(get_d, shape_par);
		par.scale = MPFR_FUN(get_d, scale_par);
		par.sign = sign_par;
	}
	mpfr_array_clear_free(coef_m, nterms);
	mpfr_array_clear_free(temp, ntemp);
	return par;
}
#define FUN_NAME poly_exp_proposal_param_max_using_mpfr
#include <count_temp.h>
#define FUN_NAME1 poly_exp_proposal_param_max_using_mpfr

void mpfr_poly_exp_find_max(mpfr_ptr xmax, size_t nterms, mpfr_t *coef, mpfr_ptr c1, mpfr_ptr c2, mpfr_ptr beta) {
	size_t ntemp = NTEMP(mpfr_poly_exp_find_max);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr init = GET_NEXT_TEMP_PTR(temp);
	if (mpfr_inf_p(c1)) {
		MPFR_FUN(set, init, c2);
		MPFR_FUN(sub_ui, init, init, 1); // can not use c2 because the function can be zero here
	} else if (mpfr_inf_p(c2)) {
		MPFR_FUN(set, init, c1);
		MPFR_FUN(add_ui, init, init, 1); // can not use c1 because the function can be zero here
	} else {
		MPFR_FUN(add, init, c1, c2);
		MPFR_FUN(div_ui, init, init, 2);
	}
	mpfr_optimization_newton(xmax, init, c1, c2, mpfr_poly_exp_gradient_log_v, mpfr_poly_exp_hessian_log_v, nterms, coef, beta);
	//mpfr_printf("init: %Rg max: %Rg\n", init, xmax);
	mpfr_array_clear_free(temp, ntemp);
}
#define FUN_NAME mpfr_poly_exp_find_max
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_exp_find_max

void mpfr_poly_exp_gaussian_proposal_param_max(mpfr_ptr mean_par, mpfr_ptr stdev_par, size_t nterms, mpfr_t* coef, mpfr_ptr c1, mpfr_ptr c2, mpfr_ptr beta) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_exp_pdf_pass_temp)
	size_t ntemp = NTEMP(mpfr_poly_exp_gaussian_proposal_param_max);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	//MPFR_PRINT_VECTOR("coef", coef, nterms);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr xmax = GET_NEXT_TEMP_PTR(temp);
	
	/* finding the maximum */
	mpfr_poly_exp_find_max(xmax, nterms, coef, c1, c2, beta);
	
	mpfr_ptr xc1 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr xc2 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr mid = GET_NEXT_TEMP_PTR(temp);
	/* set xc1 to almost c1 (99%) */
	MPFR_FUN(mul_ui, xc1, c1, 99);
	MPFR_FUN(add, xc1, xmax, xc1);
	MPFR_FUN(div_ui, xc1, xc1, 100);
	/* set xc2 to almost c2 (99%) */
	MPFR_FUN(mul_ui, xc2, c2, 99);
	MPFR_FUN(add, xc2, xmax, xc2);
	MPFR_FUN(div_ui, xc2, xc2, 100);
	/* set mid to the middle between c1 and c2 (assume that they are both finite, which is always the case when we need to use this function) */
	MPFR_FUN(add, mid, c1, c2);
	MPFR_FUN(div_ui, mid, mid, 2);
	
	mpfr_ptr fc1 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr fc2 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr fxc1 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr fxc2 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr fmax = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr fmid = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	
	mpfr_poly_exp_pdf_pass_temp(fc1, c1, nterms, coef, c1, c2, beta, pass_temp);
	mpfr_poly_exp_pdf_pass_temp(fc2, c2, nterms, coef, c1, c2, beta, pass_temp);
	mpfr_poly_exp_pdf_pass_temp(fxc1, xc1, nterms, coef, c1, c2, beta, pass_temp);
	mpfr_poly_exp_pdf_pass_temp(fxc2, xc2, nterms, coef, c1, c2, beta, pass_temp);
	mpfr_poly_exp_pdf_pass_temp(fmax, xmax, nterms, coef, c1, c2, beta, pass_temp);
	mpfr_poly_exp_pdf_pass_temp(fmid, mid, nterms, coef, c1, c2, beta, pass_temp);
	/*
	TEST_IN_FUNCTION_START
		mpfr_printf("xc1: %Rg xc2: %Rg xinit: %Rg xmax: %Rg\n", xc1, xc2, init, xmax);
		mpfr_printf("f(xc1): %Rg f(xc2): %Rg f_init: %Rg f_max: %Rg\n", fxc1, fxc2, finit, fmax);
	TEST_IN_FUNCTION_END
	*/
	/* checking for correct maximization */
	if (mpfr_cmpabs(fc1,fmax)>=0) {
		MPFR_FUN(set, fmax, fc1);
		MPFR_FUN(set, xmax, c1);
	} else if (mpfr_cmpabs(fc2,fmax)>=0) {
		MPFR_FUN(set, fmax, fc2);
		MPFR_FUN(set, xmax, c2);
	}
	
	mpfr_ptr xcmp;
	mpfr_ptr fcmp;
	if (mpfr_cmpabs(fxc1,fxc2)>=0 && mpfr_sgn(fxc1)!=0 && !mpfr_equal_p(xmax,xc1) && !mpfr_inf_p(xc1) && !mpfr_equal_p(fmax, fxc1)) {
		fcmp = fxc1;
		xcmp = xc1;
	} else if (mpfr_cmpabs(fxc2,fxc1)>=0 && mpfr_sgn(fxc2)!=0 && !mpfr_equal_p(xmax,xc2) && !mpfr_inf_p(xc2) && !mpfr_equal_p(fmax, fxc2)) {
		fcmp = fxc2;
		xcmp = xc2;
	} else if (mpfr_sgn(fxc1)!=0 && !mpfr_equal_p(xmax,xc1) && !mpfr_inf_p(xc1) && !mpfr_equal_p(fmax, fxc1)){
		fcmp = fxc1;
		xcmp = xc1;
	} else if (mpfr_sgn(fxc2)!=0 && !mpfr_equal_p(xmax,xc2) && !mpfr_inf_p(xc2) && !mpfr_equal_p(fmax, fxc2)) {
		fcmp = fxc2;
		xcmp = xc2;
	} else {
		fcmp = fmid;
		xcmp = mid;
	}
	
	MPFR_FUN(set, mean_par, xmax);
	
	MPFR_FUN(div, tmp, fmax, fcmp);
	MPFR_FUN(log, tmp, tmp);
	MPFR_FUN(mul_ui, tmp, tmp, 2);
	MPFR_FUN(sqrt, tmp, tmp);
	MPFR_FUN(sub, stdev_par, xcmp, xmax);
	MPFR_FUN(abs, stdev_par, stdev_par);
	MPFR_FUN(div, stdev_par, stdev_par, tmp);
	/*
	TEST_IN_FUNCTION_START
		mpfr_printf("mean: %Rg stdev: %Rg\n", mean_par, stdev_par);
	TEST_IN_FUNCTION_END
	*/
	mpfr_array_clear_free(temp, ntemp);
}
#define FUN_NAME mpfr_poly_exp_gaussian_proposal_param_max
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_exp_gaussian_proposal_param_max

void mpfr_poly_exp_gamma_proposal_param_max(mpfr_ptr offset_par, mpfr_ptr shape_par, mpfr_ptr scale_par, int *sign, size_t nterms, mpfr_t* coef, mpfr_ptr c1, mpfr_ptr c2, mpfr_ptr beta) {
	size_t ntemp = NTEMP(mpfr_poly_exp_gamma_proposal_param_max);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	
	mpfr_ptr xmax = GET_NEXT_TEMP_PTR(temp);
	
	if (mpfr_sgn(beta)>0) {
		*sign = 1;
	} else {
		*sign = -1;
	}
	MPFR_FUN(si_div, scale_par, *sign, beta);
	if (*sign>0) {
		MPFR_FUN(set, offset_par, c1);
	} else {
		MPFR_FUN(set, offset_par, c2);
	}
	/* finding the maximum */
	if (nterms>0) {
		mpfr_poly_exp_find_max(xmax, nterms, coef, c1, c2, beta);
	} else {
		MPFR_FUN(set, xmax, (*sign>0)? c1: c2);
	}

	MPFR_FUN(sub, shape_par, xmax, offset_par);
	MPFR_FUN(mul, shape_par, shape_par, beta);
	MPFR_FUN(add_ui, shape_par, shape_par, 1);
	//mpfr_printf("c1: %Rg c2: %Rg nterms: %zd beta: %Rg\n",c1,c2,nterms, beta);
	//mpfr_printf("sign: %d offset: %Rg shape: %Rg scale: %Rg\n",*sign, offset_par, shape_par, scale_par);
	mpfr_array_clear_free(temp, ntemp);
}
#define FUN_NAME mpfr_poly_exp_gamma_proposal_param_max
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_exp_gamma_proposal_param_max

double poly_exp_pdf_ratio_using_mpfr(double x1, double x2, size_t nterms, double *terms, double c1, double c2, double beta) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_exp_pdf_pass_temp)
	size_t ntemp = NTEMP(poly_exp_pdf_ratio_using_mpfr);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr x1_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr x2_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr c1_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr c2_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr beta_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr dist1_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr dist2_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr ratio_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_t *terms_mpfr = mpfr_array_alloc_init(nterms);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(set_d,terms_mpfr[i], terms[i]);
	}
	MPFR_FUN(set_d,x1_mpfr, x1);
	MPFR_FUN(set_d,x2_mpfr, x2);
	MPFR_FUN(set_d,c1_mpfr, c1);
	MPFR_FUN(set_d,c2_mpfr, c2);
	MPFR_FUN(set_d,beta_mpfr, beta);
	
	mpfr_poly_exp_pdf_pass_temp(dist1_mpfr, x1_mpfr, nterms, terms_mpfr, c1_mpfr, c2_mpfr, beta_mpfr, pass_temp);
	mpfr_poly_exp_pdf_pass_temp(dist2_mpfr, x2_mpfr, nterms, terms_mpfr, c1_mpfr, c2_mpfr, beta_mpfr, pass_temp);

	MPFR_FUN(div, ratio_mpfr, dist1_mpfr, dist2_mpfr);
	
	double ratio = MPFR_FUN(get_d,ratio_mpfr);

	mpfr_array_clear_free(temp, ntemp);
	mpfr_array_clear_free(terms_mpfr, nterms);
	
	return ratio;
}
#define FUN_NAME poly_exp_pdf_ratio_using_mpfr
#include <count_temp.h>
#define FUN_NAME1 poly_exp_pdf_ratio_using_mpfr

double poly_exp_pdf_using_mpfr(double x, size_t nterms, double *terms, double c1, double c2, double beta) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_exp_pdf_pass_temp)
	size_t ntemp = NTEMP(poly_exp_pdf_using_mpfr);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr x_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr c1_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr c2_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr beta_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr dist_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_t *terms_mpfr = mpfr_array_alloc_init(nterms);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(set_d,terms_mpfr[i], terms[i]);
	}
	MPFR_FUN(set_d,x_mpfr, x);
	MPFR_FUN(set_d,c1_mpfr, c1);
	MPFR_FUN(set_d,c2_mpfr, c2);
	MPFR_FUN(set_d,beta_mpfr, beta);
	
	mpfr_poly_exp_pdf_pass_temp(dist_mpfr, x_mpfr, nterms, terms_mpfr, c1_mpfr, c2_mpfr, beta_mpfr, pass_temp);
	
	double dist = MPFR_FUN(get_d,dist_mpfr);
	
	mpfr_array_clear_free(temp, ntemp);
	mpfr_array_clear_free(terms_mpfr, nterms);
	
	return dist;
}
#define FUN_NAME poly_exp_pdf_using_mpfr
#include <count_temp.h>
#define FUN_NAME1 poly_exp_pdf_using_mpfr

void mpfr_poly_exp_pdf(mpfr_ptr dist, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr c1, mpfr_ptr c2, mpfr_ptr beta) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_exp_pdf_pass_temp)
	size_t ntemp = NTEMP(mpfr_poly_exp_pdf);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_poly_exp_pdf_pass_temp(dist, x, nterms, terms, c1, c2, beta, pass_temp);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_poly_exp_pdf
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_exp_pdf

void mpfr_poly_exp_pdf_pass_temp(mpfr_ptr dist, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr c1, mpfr_ptr c2, mpfr_ptr beta, mpfr_t *temp) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_eval_mult_pass_temp)
	if (mpfr_cmp(x,c1)<0 || mpfr_cmp(x,c2)>0 || mpfr_inf_p(x)) {
		MPFR_FUN(set_ui, dist, 0);
		return;
	}
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr val = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_poly_eval_mult_pass_temp(val, x, nterms, terms, pass_temp);
	MPFR_FUN(mul, tmp, x, beta);
	MPFR_FUN(neg, tmp, tmp);
	MPFR_FUN(exp, tmp, tmp);
	MPFR_FUN(mul, dist, tmp, val);
	return;
}
#define FUN_NAME mpfr_poly_exp_pdf_pass_temp
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_exp_pdf_pass_temp

double poly_exp_proposal_pdf(double x, poly_exp_proposal_param_t par) {
	double dist;
	if (par.gamma) {
		dist = gsl_ran_gamma_pdf(par.sign*x-par.offset, par.shape, par.scale);
	} else {
		dist = gsl_ran_gaussian_pdf(x-par.mean, par.stdev);
	}
	return dist;
}

double poly_exp_proposal_samp(poly_exp_proposal_param_t par, gsl_rng* ran_gen) {
	double samp;
	if (par.gamma) {
		samp = par.offset + par.sign*gsl_ran_gamma(ran_gen, par.shape, par.scale);
	} else {
		samp = par.mean+gsl_ran_gaussian(ran_gen, par.stdev);
	}
	return samp;
}
