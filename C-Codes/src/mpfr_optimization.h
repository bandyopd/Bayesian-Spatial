#ifndef _MPFR_OPTIMIZATION_H_INCLUDED_
#define _MPFR_OPTIMIZATION_H_INCLUDED_

#include <stdarg.h>

#include <mpfr_array.h>

#include <test_in_function.h>

#include <count_temp_decl.h>

#define MPFR_OPTIMIZATION_MAX_ITER 500

void mpfr_optimization_gradient_descent_max(mpfr_ptr xmax, mpfr_ptr xinit, mpfr_ptr lbound, mpfr_ptr ubound, void (*gradient_fun)(mpfr_ptr grad, mpfr_ptr x, va_list pass_args), mpfr_ptr gamma, ...);
NTEMP_DECL(mpfr_optimization_gradient_descent_max);

void mpfr_optimization_gradient_descent_max_v(mpfr_ptr xmax, mpfr_ptr xinit, mpfr_ptr lbound, mpfr_ptr ubound, void (*gradient_fun)(mpfr_ptr grad, mpfr_ptr x, va_list pass_args), mpfr_ptr gamma, va_list pass_args);
NTEMP_DECL(mpfr_optimization_gradient_descent_max_v);

void mpfr_optimization_newton(mpfr_ptr xmax, mpfr_ptr xinit, mpfr_ptr lbound, mpfr_ptr ubound, void (*gradient_fun)(mpfr_ptr grad, mpfr_ptr x, va_list pass_args), void (*hessian_fun)(mpfr_ptr hess, mpfr_ptr x, va_list pass_args), ...);
NTEMP_DECL(mpfr_optimization_newton);

void mpfr_optimization_newton_v(mpfr_ptr xmax, mpfr_ptr xinit, mpfr_ptr lbound, mpfr_ptr ubound, void (*gradient_fun)(mpfr_ptr grad, mpfr_ptr x, va_list pass_args), void (*hessian_fun)(mpfr_ptr hess, mpfr_ptr x, va_list pass_args), va_list pass_args);
NTEMP_DECL(mpfr_optimization_newton_v);

#endif