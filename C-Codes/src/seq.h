#ifndef _SEQ_H_INCLUDED_
#define _SEQ_H_INCLUDED_

#include <stdlib.h>

double* seq_double(double* dest, double start, double end, size_t n);
size_t seq_double_step(double* dest, double start, double end, double step, size_t max_size);

int* seq_int(int* dest, int start, int end, size_t n);
size_t seq_int_step(int* dest, int start, int end, int step, size_t max_size);

#endif