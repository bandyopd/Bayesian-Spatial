#ifndef _POLY_EVAL_H_INCLUDED_
#define _POLY_EVAL_H_INCLUDED_

#include <math.h>
#include <stdlib.h>

#include <mpfr_array.h>

#include <count_temp_decl.h>

double poly_eval_mult_regular(double x, size_t nterms, double* coef);
double poly_eval_mult_using_mpf(double x, size_t nterms, double* coef);
NTEMP_DECL(poly_eval_mult_using_mpf);
void poly_eval_mult_mpf(mpf_t val, mpf_t x, size_t nterms, mpf_t* coef, mpf_t* temp);
NTEMP_DECL(poly_eval_mult_mpf);
void mpfr_poly_eval_mult_pass_temp(mpfr_t val, mpfr_t x, size_t nterms, mpfr_t* coef, mpfr_t* temp);
NTEMP_DECL(mpfr_poly_eval_mult_pass_temp);
double poly_eval(double x, size_t nterms, double* coef);
void poly_eval_mpf(mpf_t val, mpf_t x, size_t nterms, mpf_t* coef, mpf_t *temp);
NTEMP_DECL(poly_eval_mpf);
void mpfr_poly_eval(mpfr_ptr val, mpfr_ptr x, size_t nterms, mpfr_t* coef, mpfr_t *temp);
NTEMP_DECL(mpfr_poly_eval);

#endif