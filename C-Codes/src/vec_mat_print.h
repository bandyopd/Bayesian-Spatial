#ifndef _VEC_MAT_PRINT_H_INCLUDED_
#define _VEC_MAT_PRINT_H_INCLUDED_

#include <stdio.h>
#include <gmp.h>

#define PRINT_VECTOR(NAME,V,SIZE) {												\
	const char* _v_print_format_= sizeof((V)[0])==sizeof(double)?"%g ":"%d ";	\
	PRINT_VECTOR_FORMAT(NAME,V,SIZE,_v_print_format_)							\
}

#define PRINT_MATRIX(NAME,M,SIZE1,SIZE2) {		 								\
	const char* _v_print_format_= sizeof((M)[0])==sizeof(double)?"%g ":"%d ";	\
	PRINT_MATRIX_FORMAT(NAME,M,SIZE1,SIZE2,_v_print_format_)					\
}

#define PRINT_VECTOR_FORMAT(NAME,V,SIZE,FORMAT) {							\
	printf(NAME ":\n");														\
	for (size_t i=0; i<SIZE; i++) {											\
		printf(FORMAT,(V)[i]);												\
	}																		\
	printf("\n");															\
}

#define PRINT_MATRIX_FORMAT(NAME,M,SIZE1,SIZE2,FORMAT) {					\
	printf(NAME ":\n");														\
	for (size_t i=0; i<SIZE1; i++) {										\
		for (size_t _v_print_j_=0;_v_print_j_<SIZE2;_v_print_j_++){ 		\
			printf(FORMAT,(M)[i*SIZE2+_v_print_j_]);						\
		}																	\
		printf("\n");														\
	}																		\
}

#define GMP_PRINT_VECTOR_MPF(NAME,V,SIZE) GMP_PRINT_VECTOR_FORMAT(NAME,V,SIZE,"%Fg ")

#define GMP_PRINT_VECTOR_FORMAT(NAME,V,SIZE,FORMAT) {						\
	printf(NAME ":\n");														\
	for (size_t i=0; i<SIZE; i++) {											\
		gmp_printf(FORMAT,(V)[i]);											\
	}																		\
	printf("\n");															\
}

#define GMP_PRINT_MATRIX_MPF(NAME,M,SIZE1,SIZE2) GMP_PRINT_MATRIX_FORMAT(NAME,M,SIZE1,SIZE2,"%Fg ")

#define GMP_PRINT_MATRIX_FORMAT(NAME,M,SIZE1,SIZE2,FORMAT) {				\
	printf(NAME ":\n");														\
	for (size_t i=0; i<SIZE1; i++) {										\
		for (size_t _v_print_j_=0;_v_print_j_<SIZE2;_v_print_j_++){ 		\
			gmp_printf(FORMAT,(M)[i*SIZE2+_v_print_j_]);					\
		}																	\
		printf("\n");														\
	}																		\
}

#define MPFR_PRINT_VECTOR_FORMAT(NAME,V,SIZE,FORMAT) {						\
	printf(NAME ":\n");														\
	for (size_t i=0; i<SIZE; i++) {											\
		mpfr_printf(FORMAT,(V)[i]);											\
	}																		\
	printf("\n");															\
}

#define MPFR_PRINT_VECTOR(NAME,V,SIZE) MPFR_PRINT_VECTOR_FORMAT(NAME,V,SIZE,"%Rg ")

#endif