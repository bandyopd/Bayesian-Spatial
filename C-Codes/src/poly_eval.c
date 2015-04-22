#include <poly_eval.h>

#include <count_temp_start.h>

double poly_eval_mult_regular(double x, size_t nterms, double* coef) {
	if (nterms==0 || coef==NULL)
		return 1;
	double px=1;
	for (size_t i=0; i<nterms; i++) {
		px *= (x+coef[i]);
	}
	return px;
}

double poly_eval_mult_using_mpf(double x, size_t nterms, double* coef) {
	if (nterms==0 || coef==NULL)
		return 1.0;
	size_t ntemp = NTEMP(poly_eval_mult_using_mpf);
	mpf_t *temp = mpf_array_alloc(ntemp);
	mpf_array_init(temp,ntemp);
	
	mpf_t *px = GET_NEXT_TEMP(temp);
	mpf_t *tmp = GET_NEXT_TEMP(temp);
	
	mpf_set_d(*px,1.0);
	for (size_t i=0; i<nterms; i++) {
		mpf_set_d(*tmp,x+coef[i]);
		mpf_mul(*px,*px,*tmp);
	}
	double pval = mpf_get_d(*px);
	mpf_array_clear(temp,ntemp);
	free(temp);
	return pval;
}
#define FUN_NAME poly_eval_mult_using_mpf
#include <count_temp.h>
#define FUN_NAME1 poly_eval_mult_using_mpf

void poly_eval_mult_mpf(mpf_t val, mpf_t x, size_t nterms, mpf_t* coef, mpf_t* temp) {
	mpf_set_d(val,1.0);
	if (nterms==0 || coef==NULL) {
		return;
	}
	for (size_t i=0; i<nterms; i++) {
		mpf_t *tmp = GET_NEXT_TEMP(temp);
		mpf_add(*tmp, x, coef[i]);
		mpf_mul(val, val, *tmp);
	}
	return;
}
#define FUN_NAME poly_eval_mult_mpf
#include <count_temp.h>
#define FUN_NAME1 poly_eval_mult_mpf

void mpfr_poly_eval_mult_pass_temp(mpfr_t val, mpfr_t x, size_t nterms, mpfr_t* coef, mpfr_t* temp) {
	MPFR_FUN(set_ui,val,1);
	if (nterms==0 || coef==NULL) {
		return;
	}
	for (size_t i=0; i<nterms; i++) {
		mpfr_t *tmp = GET_NEXT_TEMP(temp);
		MPFR_FUN(add, *tmp, x, coef[i]);
		MPFR_FUN(mul, val, val, *tmp);
	}
	return;
}
#define FUN_NAME mpfr_poly_eval_mult_pass_temp
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_eval_mult_pass_temp

double poly_eval(double x, size_t nterms, double* coef) {
	if (nterms<1 || coef==NULL)
		return NAN;
	double px=0.0;
	double x_pow=1.0;
	for (size_t i=0; i<nterms; i++) {
		px += x_pow*coef[i];
		x_pow *= x;
	}
	return px;
}

void poly_eval_mpf(mpf_t val, mpf_t x, size_t nterms, mpf_t* coef, mpf_t *temp) {
	if (nterms<1 || coef==NULL) {
		mpf_set_d(val,1.0);
		return;
	}
	mpf_t* x_pow = GET_NEXT_TEMP(temp);
	mpf_set_d(*x_pow,1.0);
	mpf_set_d(val,0.0);
	for (size_t i=0; i<nterms; i++) {
		mpf_t* tmp = GET_NEXT_TEMP(temp);
		mpf_mul(*tmp,*x_pow, coef[i]);
		mpf_add(val,val,*tmp);
		mpf_mul(*x_pow,*x_pow,x);
	}
	return;
}
#define FUN_NAME poly_eval_mpf
#include <count_temp.h>
#define FUN_NAME1 poly_eval_mpf

void mpfr_poly_eval(mpfr_ptr val, mpfr_ptr x, size_t nterms, mpfr_t* coef, mpfr_t *temp) {
	if (nterms<1 || coef==NULL) {
		MPFR_FUN(set_ui,val,1);
		return;
	}
	mpfr_ptr x_pow = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr sum = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(set_ui,x_pow,1);
	MPFR_FUN(set_ui,sum,0);
	for (size_t i=0; i<nterms; i++) {
		mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
		MPFR_FUN(mul,tmp,x_pow, coef[i]);
		MPFR_FUN(add,sum,sum,tmp);
		MPFR_FUN(mul,x_pow,x_pow,x);
	}
	MPFR_FUN(set, val, sum);
	return;
}
#define FUN_NAME mpfr_poly_eval
#include <count_temp.h>
#define FUN_NAME1 mpfr_poly_eval

#include <count_temp_end.h>