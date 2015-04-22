#ifndef _MPF_ARRAY_H_INCLUDED_
#define _MPF_ARRAY_H_INCLUDED_

#include <stdlib.h>

#include <gmp.h>
#include <mpfr.h>

#include <count_temp_decl.h>

#ifndef DEFAULT_ARRAY_PREC
	#define DEFAULT_ARRAY_PREC 1<<7
#endif

#define MPFR_FUN(NAME,ARGS...) mpfr_##NAME(ARGS,mpfr_get_default_rounding_mode())
#define MPF_FUN(NAME,ARGS...) mpf_##NAME(ARGS)

#define MAX(a,b) ( ((a)>(b)) ? (a) : (b) )
#define MIN(a,b) ( ((a)<(b)) ? (a) : (b) )

mpf_t *mpf_array_alloc(size_t size);
void mpf_array_init(mpf_t *arr, size_t size);
void mpf_array_init_prec(mpf_t *arr, size_t size, mp_bitcnt_t prec);
mpf_t *mpf_array_alloc_init_prec(size_t size, mp_bitcnt_t prec);
mpf_t *mpf_array_alloc_init(size_t size);
void mpf_array_init_one(mpf_t x);
void mpf_array_init_one_prec(mpf_t x, mp_bitcnt_t prec);
void mpf_array_clear(mpf_t *arr, size_t size);
void mpf_array_clear_free(mpf_t *arr, size_t size);

void mpf_array_get_epsilon(mpf_t eps);
NTEMP_DECL(mpf_array_get_epsilon);

void mpf_array_get_epsilon_routine(mpf_t eps, mpf_t* temp);
NTEMP_DECL(mpf_array_get_epsilon_routine);

mpfr_t *mpfr_array_alloc(size_t size);
mpfr_t *mpfr_array_alloc_init_prec(size_t size, mpfr_prec_t prec);
mpfr_t *mpfr_array_alloc_init(size_t size);
void mpfr_array_init(mpfr_t *arr, size_t size);
void mpfr_array_init_prec(mpfr_t *arr, size_t size, mpfr_prec_t prec);
void mpfr_array_clear(mpfr_t *arr, size_t size);
void mpfr_array_clear_free(mpfr_t *arr, size_t size);

mpfr_ptr *mpfr_ptr_array_alloc(size_t size);
mpfr_ptr *mpfr_ptr_array_alloc_init_prec(size_t size, mpfr_prec_t prec);
mpfr_ptr *mpfr_ptr_array_alloc_init(size_t size);
void mpfr_ptr_array_init(mpfr_ptr *arr, size_t size);
void mpfr_ptr_array_init_prec(mpfr_ptr *arr, size_t size, mpfr_prec_t prec);
void mpfr_ptr_array_clear(mpfr_ptr *arr, size_t size);
void mpfr_ptr_array_clear_free(mpfr_ptr *arr, size_t size);

void mpfr_array_init_one(mpfr_ptr x);
void mpfr_array_init_one_prec(mpfr_ptr x, mpfr_prec_t prec);

void mpfr_array_get_epsilon(mpfr_ptr eps);
NTEMP_DECL(mpfr_array_get_epsilon);

void mpfr_array_get_epsilon_routine(mpfr_ptr eps, mpfr_t* temp);
NTEMP_DECL(mpfr_array_get_epsilon_routine);

#endif
