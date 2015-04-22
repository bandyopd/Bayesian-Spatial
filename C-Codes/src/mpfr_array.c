#include <mpfr_array.h>

#include <count_temp_start.h>

mpf_t *mpf_array_alloc(size_t size) {
	mpf_t *arr = malloc(size*sizeof(*arr));
	return arr;
}

void mpf_array_init(mpf_t *arr, size_t size) {
	mpf_array_init_prec(arr, size, DEFAULT_ARRAY_PREC);
}

mpf_t *mpf_array_alloc_init_prec(size_t size, mp_bitcnt_t prec) {
	mpf_t *arr = mpf_array_alloc(size);
	if (arr!=0) {
		mpf_array_init_prec(arr, size, prec);
	}
	return arr;
}

mpf_t *mpf_array_alloc_init(size_t size) {
	return mpf_array_alloc_init_prec(size, DEFAULT_ARRAY_PREC);
}

void mpf_array_init_prec(mpf_t *arr, size_t size, mp_bitcnt_t prec) {
	for (size_t i=0; i<size; i++) {
		mpf_init2(arr[i],prec);
	}
	return;
}

void mpf_array_init_one(mpf_t x) {
	mpf_array_init_one_prec(x, DEFAULT_ARRAY_PREC);
}

void mpf_array_init_one_prec(mpf_t x, mp_bitcnt_t prec) {
	mpf_init2(x,prec);
	return;
}

void mpf_array_clear(mpf_t *arr, size_t size) {
	for (size_t i=0; i<size; i++) {
		mpf_clear(arr[i]);
	}
	return;
}

void mpf_array_clear_free(mpf_t *arr, size_t size) {
	mpf_array_clear(arr,size);
	free(arr);
	return;
}

void mpf_array_get_epsilon(mpf_t eps) {
	#undef NPASS
	#define NPASS NTEMP(mpf_array_get_epsilon_routine)
	mp_bitcnt_t prec = mpf_get_prec(eps);
	size_t ntemp = NTEMP(mpf_array_get_epsilon);
	mpf_t *temp = mpf_array_alloc(ntemp);
	mpf_array_init_prec(temp,ntemp,prec);
	mpf_t *pass_temp = GET_NEXT_TEMP(temp);
	mpf_array_get_epsilon_routine(eps,pass_temp);
	mpf_array_clear(temp,ntemp);
	free(temp);
}
#define FUN_NAME mpf_array_get_epsilon
#include <count_temp.h>
#define FUN_NAME1 mpf_array_get_epsilon

void mpf_array_get_epsilon_routine(mpf_t eps, mpf_t* temp) {
	MPF_FUN(set_ui,eps,1);
	mpf_t* eps_p_1 = GET_NEXT_TEMP(temp);
	mpf_t *one = GET_NEXT_TEMP(temp);
	MPF_FUN(set_ui,*one,1);
	do {
		MPF_FUN(div_ui,eps,eps,2);
		MPF_FUN(add,*eps_p_1, eps, *one);
	} while (mpf_cmp(*eps_p_1,*one)!=0);
	MPF_FUN(mul_ui, eps, eps, 2);
	return;
}
#define FUN_NAME mpf_array_get_epsilon_routine
#include <count_temp.h>
#define FUN_NAME1 mpf_array_get_epsilon_routine

mpfr_t *mpfr_array_alloc(size_t size) {
	mpfr_t *arr = malloc(size*sizeof(*arr));
	return arr;
}

mpfr_t *mpfr_array_alloc_init_prec(size_t size, mpfr_prec_t prec) {
	mpfr_t *arr = mpfr_array_alloc(size);
	if (arr!=0) {
		mpfr_array_init_prec(arr, size, prec);
	}
	return arr;
}

mpfr_t *mpfr_array_alloc_init(size_t size) {
	return mpfr_array_alloc_init_prec(size,DEFAULT_ARRAY_PREC);
}

void mpfr_array_init(mpfr_t *arr, size_t size) {
	mpfr_array_init_prec(arr, size, DEFAULT_ARRAY_PREC);
}

void mpfr_array_init_prec(mpfr_t *arr, size_t size, mpfr_prec_t prec) {
	for (size_t i=0; i<size; i++) {
		mpfr_init2(arr[i],prec);
	}
	return;
}

void mpfr_array_clear(mpfr_t *arr, size_t size) {
	if (arr==0) {
		return;
	}
	for (size_t i=0; i<size; i++) {
		mpfr_clear(arr[i]);
	}
	return;
}

void mpfr_array_clear_free(mpfr_t *arr, size_t size) {
	mpfr_array_clear(arr, size);
	free(arr);
}

mpfr_ptr *mpfr_ptr_array_alloc(size_t size) {
	mpfr_ptr *arr = malloc(size*sizeof(*arr));
	return arr;
}

mpfr_ptr *mpfr_ptr_array_alloc_init_prec(size_t size, mpfr_prec_t prec) {
	mpfr_ptr *arr = mpfr_ptr_array_alloc(size);
	if (arr!=0) {
		mpfr_ptr_array_init_prec(arr, size, prec);
	}
	return arr;
}

mpfr_ptr *mpfr_ptr_array_alloc_init(size_t size) {
	return mpfr_ptr_array_alloc_init_prec(size,DEFAULT_ARRAY_PREC);
}

void mpfr_ptr_array_init(mpfr_ptr *arr, size_t size) {
	mpfr_ptr_array_init_prec(arr, size, DEFAULT_ARRAY_PREC);
}

void mpfr_ptr_array_init_prec(mpfr_ptr *arr, size_t size, mpfr_prec_t prec) {
	for (size_t i=0; i<size; i++) {
		arr[i] = malloc(sizeof(*(arr[i])));
		mpfr_array_init_one_prec(arr[i],prec);
	}
	return;
}

void mpfr_ptr_array_clear(mpfr_ptr *arr, size_t size) {
	if (arr==0) {
		return;
	}
	for (size_t i=0; i<size; i++) {
		mpfr_clear(arr[i]);
		free(arr[i]);
	}
	return;
}

void mpfr_ptr_array_clear_free(mpfr_ptr *arr, size_t size) {
	mpfr_ptr_array_clear(arr, size);
	free(arr);
}

void mpfr_array_init_one(mpfr_ptr x) {
	mpfr_array_init_one_prec(x, DEFAULT_ARRAY_PREC);
}

void mpfr_array_init_one_prec(mpfr_ptr x, mpfr_prec_t prec) {
	mpfr_init2(x,prec);
	return;
}

void mpfr_array_get_epsilon(mpfr_ptr eps) {
	#undef NPASS
	#define NPASS NTEMP(mpfr_array_get_epsilon_routine)
	mpfr_prec_t prec = mpfr_get_prec(eps);
	size_t ntemp = NTEMP(mpfr_array_get_epsilon);
	mpfr_t *temp = mpfr_array_alloc(ntemp);
	mpfr_array_init_prec(temp,ntemp,prec);
	size_t pass_temp_ind = NEXT_TEMP_INDEX;
	mpfr_array_get_epsilon_routine(eps,temp+pass_temp_ind);
	mpfr_array_clear(temp,ntemp);
	free(temp);
}
#define FUN_NAME mpfr_array_get_epsilon
#include <count_temp.h>
#define FUN_NAME1 mpfr_array_get_epsilon

void mpfr_array_get_epsilon_routine(mpfr_ptr eps, mpfr_t* temp) {
	MPFR_FUN(set_ui,eps,1);
	mpfr_ptr eps_p_1 = GET_NEXT_TEMP_PTR(temp);
	mpfr_ptr one = GET_NEXT_TEMP_PTR(temp);
	MPFR_FUN(set_ui,one,1);
	do {
		MPFR_FUN(div_ui,eps,eps,2);
		MPFR_FUN(add,eps_p_1, eps, one);
	} while (mpfr_cmp(eps_p_1,one)!=0);
	MPFR_FUN(mul_ui, eps, eps, 2);
	return;
}
#define FUN_NAME mpfr_array_get_epsilon_routine
#include <count_temp.h>
#define FUN_NAME1 mpfr_array_get_epsilon_routine


#include <count_temp_end.h>