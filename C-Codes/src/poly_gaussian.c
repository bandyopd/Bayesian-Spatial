#include <poly_gaussian.h>

#include <count_temp_start.h>

void mpfr_poly_gaussian_gradient_log(mpfr_ptr grad, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr mean, mpfr_ptr stdev) {
	size_t ntemp = NTEMP(mpfr_poly_gaussian_gradient_log);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr sum = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(set_ui, sum, 0);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(add, tmp, x, terms[i]);
		MPFR_FUN(ui_div, tmp, 1, tmp);
		MPFR_FUN(add, sum, sum, tmp);
	}
	MPFR_FUN(sub, tmp, x, mean);
	MPFR_FUN(div, tmp, tmp, stdev);
	MPFR_FUN(div, tmp, tmp, stdev);
	MPFR_FUN(sub, sum, sum, tmp);
	MPFR_FUN(set, grad, sum);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_poly_gaussian_gradient_log
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gaussian_gradient_log

void mpfr_poly_gaussian_gradient_log_v(mpfr_ptr grad, mpfr_ptr x, va_list pass_args) {
	size_t nterms = va_arg(pass_args, size_t);
	mpfr_t *terms = va_arg(pass_args, mpfr_t*);
	mpfr_ptr mean = va_arg(pass_args, mpfr_ptr);
	mpfr_ptr stdev = va_arg(pass_args, mpfr_ptr);
	mpfr_poly_gaussian_gradient_log(grad, x, nterms, terms, mean, stdev);
	return;
}
#define FUN_NAME mpfr_poly_gaussian_gradient_log_v
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gaussian_gradient_log_v

void mpfr_poly_gaussian_hessian_log(mpfr_ptr hess, mpfr_ptr x, size_t nterms, mpfr_t *terms, mpfr_ptr mean, mpfr_ptr stdev) {
	size_t ntemp = NTEMP(mpfr_poly_gaussian_hessian_log);
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
	MPFR_FUN(pow_ui, tmp, stdev, 2);
	MPFR_FUN(ui_div, tmp, 1, tmp);
	MPFR_FUN(sub, sum, sum, tmp);
	MPFR_FUN(set, hess, sum);
	mpfr_array_clear_free(temp, ntemp);
}
#define FUN_NAME mpfr_poly_gaussian_hessian_log
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gaussian_hessian_log

void mpfr_poly_gaussian_hessian_log_v(mpfr_ptr grad, mpfr_ptr x, va_list pass_args) {
	size_t nterms = va_arg(pass_args, size_t);
	mpfr_t *terms = va_arg(pass_args, mpfr_t*);
	mpfr_ptr mean = va_arg(pass_args, mpfr_ptr);
	mpfr_ptr stdev = va_arg(pass_args, mpfr_ptr);
	mpfr_poly_gaussian_hessian_log(grad, x, nterms, terms, mean, stdev);
	return;
}
#define FUN_NAME mpfr_poly_gaussian_hessian_log_v
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gaussian_hessian_log_v

poly_gaussian_proposal_param_t
			poly_gaussian_proposal_param_max_using_mpfr(size_t nterms, double* coef_mult, double constr1, double constr2, double mean, double stdev) {
	size_t ntemp = NTEMP(poly_gaussian_proposal_param_max_using_mpfr);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr constr1_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr constr2_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr mean_mpfr = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr stdev_mpfr = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(set_d,constr1_mpfr, constr1);
	MPFR_FUN(set_d,constr2_mpfr, constr2);
	MPFR_FUN(set_d,mean_mpfr, mean);
	MPFR_FUN(set_d,stdev_mpfr, stdev);
	mpfr_t *coef_mpfr = mpfr_array_alloc_init(nterms);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(set_d,coef_mpfr[i], coef_mult[i]);
	}
	mpfr_ptr mean_par = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr stdev_par = GET_NEXT_TEMP_PTR(temp);
	mpfr_poly_gaussian_proposal_param_max(mean_par, stdev_par, nterms, coef_mpfr, constr1_mpfr, constr2_mpfr, mean_mpfr, stdev_mpfr);
	poly_gaussian_proposal_param_t par;
	par.mean = MPFR_FUN(get_d,mean_par);
	par.stdev = MPFR_FUN(get_d,stdev_par);
	mpfr_array_clear_free(temp, ntemp);
	mpfr_array_clear_free(coef_mpfr, nterms);
	return par;
}
#define FUN_NAME poly_gaussian_proposal_param_max_using_mpfr
#include <count_temp.h>
#define FUN_NAME1 poly_gaussian_proposal_param_max_using_mpfr

void mpfr_poly_gaussian_proposal_param_max(mpfr_ptr mean_par, mpfr_ptr stdev_par, size_t nterms, mpfr_t* coef, mpfr_ptr c1, mpfr_ptr c2, mpfr_ptr mean, mpfr_ptr stdev) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_gaussian_pdf_pass_temp)
	size_t ntemp = NTEMP(mpfr_poly_gaussian_proposal_param_max);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr xmax = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr init = GET_NEXT_TEMP_PTR(temp);
	if (mpfr_inf_p(c1)) {
		MPFR_FUN(sub_ui, init, c2, 1); // can not use c2 because the function can be zero here
	} else if (mpfr_inf_p(c2)) {
		MPFR_FUN(add_ui, init, c1, 1); // can not use c1 because the function can be zero here
	} else {
		MPFR_FUN(add, init, c1, c2);
		MPFR_FUN(div_ui, init, init, 2);
	}
	mpfr_optimization_newton(xmax, init, c1, c2, mpfr_poly_gaussian_gradient_log_v, mpfr_poly_gaussian_hessian_log_v, nterms, coef, mean, stdev);
	
	//mpfr_ptr xc1 = GET_NEXT_TEMP_PTR(temp);
	//mpfr_ptr xc2 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr xcmp = GET_NEXT_TEMP_PTR(temp);
	/* set xc1 to almost c1 (99%) */
	//MPFR_FUN(mul_ui, xc1, c1, 99);
	//MPFR_FUN(add, xc1, xmax, xc1);
	//MPFR_FUN(div_ui, xc1, xc1, 100);
	/* set xc2 to almost c2 (99%) */
	//MPFR_FUN(mul_ui, xc2, c2, 99);
	//MPFR_FUN(add, xc2, xmax, xc2);
	//MPFR_FUN(div_ui, xc2, xc2, 100);
	
	//mpfr_ptr fc1 = GET_NEXT_TEMP_PTR(temp);
	//mpfr_ptr fc2 = GET_NEXT_TEMP_PTR(temp);
	//mpfr_ptr fxc1 = GET_NEXT_TEMP_PTR(temp);
	//mpfr_ptr fxc2 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr fcmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr fmax = GET_NEXT_TEMP_PTR(temp);
	//mpfr_ptr finit = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	
	//mpfr_poly_gaussian_pdf_pass_temp(fc1, c1, nterms, coef, c1, c2, mean, stdev, pass_temp);
	//mpfr_poly_gaussian_pdf_pass_temp(fc2, c2, nterms, coef, c1, c2, mean, stdev, pass_temp);
	//mpfr_poly_gaussian_pdf_pass_temp(fxc1, xc1, nterms, coef, c1, c2, mean, stdev, pass_temp);
	//mpfr_poly_gaussian_pdf_pass_temp(fxc2, xc2, nterms, coef, c1, c2, mean, stdev, pass_temp);
	mpfr_poly_gaussian_pdf_pass_temp(fmax, xmax, nterms, coef, c1, c2, mean, stdev, pass_temp);
	//mpfr_poly_gaussian_pdf_pass_temp(finit, init, nterms, coef, c1, c2, mean, stdev, pass_temp);
	/*
	TEST_IN_FUNCTION_START
		mpfr_printf("xc1: %Rg xc2: %Rg xinit: %Rg xmax: %Rg\n", xc1, xc2, init, xmax);
		mpfr_printf("f(xc1): %Rg f(xc2): %Rg f_init: %Rg f_max: %Rg\n", fxc1, fxc2, finit, fmax);
	TEST_IN_FUNCTION_END
	*/
	/* checking for correct maximization */
	/*if (mpfr_cmpabs(fc1,fmax)>=0) {
		MPFR_FUN(set, fmax, fc1);
		MPFR_FUN(set, xmax, c1);
	} else if (mpfr_cmpabs(fc2,fmax)>=0) {
		MPFR_FUN(set, fmax, fc2);
		MPFR_FUN(set, xmax, c2);
	}*/
/*	
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
		fcmp = finit;
		xcmp = init;
	}
*/
	if (mpfr_cmpabs(xmax, c1)>0) {
		MPFR_FUN(abs, xcmp, xmax);
	} else {
		MPFR_FUN(abs, xcmp, c1);
	}
	if (mpfr_sgn(xcmp)==0) { // shouldn't be so but just in case
		MPFR_FUN(set_ui, xcmp, 1);
	}
	MPFR_FUN(add, xcmp, xcmp, xmax);
	mpfr_poly_gaussian_pdf_pass_temp(fcmp, xcmp, nterms, coef, c1, c2, mean, stdev, pass_temp);
	
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
	return;
}
#define FUN_NAME mpfr_poly_gaussian_proposal_param_max
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gaussian_proposal_param_max

void mpfr_poly_gaussian_pdf_pass_temp(mpfr_t dist, mpfr_t x, size_t nterms, mpfr_t *coef, mpfr_t c1, mpfr_t c2, mpfr_t mean, mpfr_t stdev, mpfr_t* temp) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_eval_mult_pass_temp)
	if (mpfr_cmp(x,c1)<0 || mpfr_cmp(x,c2)>0 || mpfr_inf_p(x)) {
		MPFR_FUN(set_ui,dist,0);
		return;
	}
	mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr tmp1 = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_poly_eval_mult_pass_temp(tmp, x, nterms, coef, pass_temp);
	mpfr_pdf_gaussian(tmp1, x, mean, stdev);
	MPFR_FUN(mul,dist,tmp,tmp1);
	return;
}
#define FUN_NAME mpfr_poly_gaussian_pdf_pass_temp
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_gaussian_pdf_pass_temp

double poly_gaussian_pdf_using_mpfr(double x, size_t nterms, double *terms, double c1, double c2, double mean, double stdev) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_gaussian_pdf_pass_temp)
	size_t ntemp = NTEMP(poly_gaussian_pdf_using_mpfr);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_t *x_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *c1_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *c2_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *mean_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *stdev_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *dist_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_t *terms_mpfr = mpfr_array_alloc_init(nterms);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(set_d,terms_mpfr[i], terms[i]);
	}
	MPFR_FUN(set_d,*x_mpfr, x);
	MPFR_FUN(set_d,*c1_mpfr, c1);
	MPFR_FUN(set_d,*c2_mpfr, c2);
	MPFR_FUN(set_d,*mean_mpfr, mean);
	MPFR_FUN(set_d,*stdev_mpfr, stdev);
	
	mpfr_poly_gaussian_pdf_pass_temp(*dist_mpfr, *x_mpfr, nterms, terms_mpfr, *c1_mpfr, *c2_mpfr, *mean_mpfr, *stdev_mpfr, pass_temp);
	
	double dist = MPFR_FUN(get_d,*dist_mpfr);
	
	mpfr_array_clear_free(temp, ntemp);
	mpfr_array_clear_free(terms_mpfr, nterms);
	
	return dist;
}
#define FUN_NAME poly_gaussian_pdf_using_mpfr
#include <count_temp.h>
#define FUN_NAME1 poly_gaussian_pdf_using_mpfr

double poly_gaussian_pdf_ratio_using_mpfr(double x1, double x2, size_t nterms, double *terms, double c1, double c2, double mean, double stdev) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_poly_gaussian_pdf_pass_temp)
	size_t ntemp = NTEMP(poly_gaussian_pdf_ratio_using_mpfr);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_t *x1_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *x2_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *c1_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *c2_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *mean_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *stdev_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *dist1_mpfr = GET_NEXT_TEMP(temp);
	mpfr_t *dist2_mpfr = GET_NEXT_TEMP(temp);
	mpfr_ptr ratio_m = GET_NEXT_TEMP_PTR(temp);
	mpfr_t *pass_temp = GET_NEXT_TEMP(temp);
	mpfr_t *terms_mpfr = mpfr_array_alloc_init(nterms);
	for (size_t i=0; i<nterms; i++) {
		MPFR_FUN(set_d,terms_mpfr[i], terms[i]);
	}
	MPFR_FUN(set_d,*x1_mpfr, x1);
	MPFR_FUN(set_d,*x2_mpfr, x2);
	MPFR_FUN(set_d,*c1_mpfr, c1);
	MPFR_FUN(set_d,*c2_mpfr, c2);
	MPFR_FUN(set_d,*mean_mpfr, mean);
	MPFR_FUN(set_d,*stdev_mpfr, stdev);
	
	mpfr_poly_gaussian_pdf_pass_temp(*dist1_mpfr, *x1_mpfr, nterms, terms_mpfr, *c1_mpfr, *c2_mpfr, *mean_mpfr, *stdev_mpfr, pass_temp);
	mpfr_poly_gaussian_pdf_pass_temp(*dist2_mpfr, *x2_mpfr, nterms, terms_mpfr, *c1_mpfr, *c2_mpfr, *mean_mpfr, *stdev_mpfr, pass_temp);

	MPFR_FUN(div, ratio_m, *dist1_mpfr, *dist2_mpfr);
	
	double ratio = MPFR_FUN(get_d,ratio_m);
	
	mpfr_array_clear_free(temp, ntemp);
	mpfr_array_clear_free(terms_mpfr, nterms);
	
	return ratio;
}
#define FUN_NAME poly_gaussian_pdf_ratio_using_mpfr
#include <count_temp.h>
#define FUN_NAME1 poly_gaussian_pdf_ratio_using_mpfr

double poly_gaussian_proposal_pdf(double x, poly_gaussian_proposal_param_t par) {
	double dist = gsl_ran_gaussian_pdf(x-par.mean, par.stdev);
	return dist;
}

double poly_gaussian_proposal_samp(poly_gaussian_proposal_param_t par, gsl_rng* ran_gen) {
	double samp = par.mean+gsl_ran_gaussian(ran_gen, par.stdev);
	return samp;
}
