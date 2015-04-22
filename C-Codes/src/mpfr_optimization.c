#include <mpfr_optimization.h>

#include <count_temp_start.h>

void mpfr_optimization_gradient_descent_max(mpfr_ptr xmax, mpfr_ptr xinit, mpfr_ptr lbound, mpfr_ptr ubound, void (*gradient_fun)(mpfr_ptr grad, mpfr_ptr x, va_list pass_args), mpfr_ptr gamma, ...) {
	va_list pass_args;
	va_start(pass_args, gamma);
	mpfr_optimization_gradient_descent_max_v(xmax, xinit, lbound, ubound, gradient_fun, gamma, pass_args);
	va_end(pass_args);
	return;
}
#define FUN_NAME mpfr_optimization_gradient_descent_max
#include <count_temp.h>
#define FUN_NAME1 mpfr_optimization_gradient_descent_max

void mpfr_optimization_gradient_descent_max_v(mpfr_ptr xmax, mpfr_ptr xinit, mpfr_ptr lbound, mpfr_ptr ubound, void (*gradient_fun)(mpfr_ptr grad, mpfr_ptr x, va_list pass_args), mpfr_ptr gamma, va_list pass_args) {
	size_t niter = MPFR_OPTIMIZATION_MAX_ITER;
	size_t ntemp = NTEMP(mpfr_optimization_gradient_descent_max_v);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr xn = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr eps = GET_NEXT_TEMP_PTR(temp);
	mpfr_array_get_epsilon(eps);
	MPFR_FUN(set, xn, xinit);
	for (size_t i=0; i<niter; i++) {
		va_list p_args;
		va_copy(p_args, pass_args);
		mpfr_ptr grad = GET_NEXT_TEMP_PTR(temp);
		mpfr_ptr dx = GET_NEXT_TEMP_PTR(temp);
		mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
		(*gradient_fun)(grad, xn, p_args);
		va_end(p_args);
		MPFR_FUN(mul, dx, grad, gamma);
		MPFR_FUN(add, tmp, xn, dx);
		if (mpfr_cmp(tmp, lbound)<0) {			
			MPFR_FUN(sub, dx, lbound, xn);
			MPFR_FUN(set, xn, lbound);
			MPFR_FUN(mul_d, gamma, gamma, 0.9);
		} else if (mpfr_cmp(tmp, ubound)>0) {
			MPFR_FUN(sub, dx, ubound, xn);
			MPFR_FUN(set, xn, ubound);
			MPFR_FUN(mul_d, gamma, gamma, 0.9);
		} else {
			MPFR_FUN(set, xn, tmp);
		}
		if (mpfr_cmpabs(dx, eps)<0) {
			break;
		}
	}
	MPFR_FUN(set, xmax, xn);
	mpfr_array_clear_free(temp, ntemp);
	return;
}
#define FUN_NAME mpfr_optimization_gradient_descent_max_v
#include <count_temp.h>
#define FUN_NAME1 mpfr_optimization_gradient_descent_max_v

void mpfr_optimization_newton(mpfr_ptr xmax, mpfr_ptr xinit, mpfr_ptr lbound, mpfr_ptr ubound, void (*gradient_fun)(mpfr_ptr grad, mpfr_ptr x, va_list pass_args), void (*hessian_fun)(mpfr_ptr hess, mpfr_ptr x, va_list pass_args), ...) {
	va_list pass_args;
	va_start(pass_args, hessian_fun);
	mpfr_optimization_newton_v(xmax, xinit, lbound, ubound, gradient_fun, hessian_fun, pass_args);
	va_end(pass_args);
	return;
}
#define FUN_NAME mpfr_optimization_newton
#include <count_temp.h>
#define FUN_NAME1 mpfr_optimization_newton

void mpfr_optimization_newton_v(mpfr_ptr xmax, mpfr_ptr xinit, mpfr_ptr lbound, mpfr_ptr ubound, void (*gradient_fun)(mpfr_ptr grad, mpfr_ptr x, va_list pass_args), void (*hessian_fun)(mpfr_ptr hess, mpfr_ptr x, va_list pass_args), va_list pass_args) {
	size_t niter = MPFR_OPTIMIZATION_MAX_ITER;
	size_t ntemp = NTEMP(mpfr_optimization_newton_v);
	mpfr_t *temp = mpfr_array_alloc_init(ntemp);
	mpfr_ptr xn = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr eps = GET_NEXT_TEMP_PTR(temp);
	mpfr_array_get_epsilon(eps);
	MPFR_FUN(sqrt, eps, eps);
	MPFR_FUN(set, xn, xinit);
	for (size_t i=0; i<niter; i++) {
		mpfr_ptr dx = GET_NEXT_TEMP_PTR(temp);
		mpfr_ptr grad = GET_NEXT_TEMP_PTR(temp);
		mpfr_ptr hess = GET_NEXT_TEMP_PTR(temp);
		mpfr_ptr tmp = GET_NEXT_TEMP_PTR(temp);
		va_list p_args1, p_args2;
		va_copy(p_args1, pass_args);
		va_copy(p_args2, pass_args);
		(*gradient_fun)(grad, xn, p_args1);
		(*hessian_fun)(hess, xn, p_args2);
		va_end(p_args1);
		va_end(p_args2);
		MPFR_FUN(div, dx, grad, hess);
		MPFR_FUN(neg, dx, dx);
		if (mpfr_nan_p(dx)||mpfr_inf_p(dx)) {
			//mpfr_set_nan(xn);
			break;
		}
		MPFR_FUN(add, tmp, xn, dx);
		if (mpfr_cmp(tmp, lbound)<=0) {
			MPFR_FUN(sub, dx, lbound, xn);
			/* do not let xn reach the bound */
			MPFR_FUN(mul_ui, dx, dx, 99);
			MPFR_FUN(div_ui, dx, dx, 100);
		} else if (mpfr_cmp(tmp, ubound)>=0) {
			MPFR_FUN(sub, dx, ubound, xn);
			/* do not let xn reach the bound */
			MPFR_FUN(mul_ui, dx, dx, 99);
			MPFR_FUN(div_ui, dx, dx, 100);
		}
		MPFR_FUN(set, tmp, xn);
		MPFR_FUN(add, xn, xn, dx);
		MPFR_FUN(div, dx, dx, xn);
		if (mpfr_cmpabs(dx, eps)<0) {
			break;
		}
	}
	MPFR_FUN(set, xmax, xn);
	mpfr_array_clear_free(temp, ntemp);
}
#define FUN_NAME mpfr_optimization_newton_v
#include <count_temp.h>
#define FUN_NAME1 mpfr_optimization_newton_v
