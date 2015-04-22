#include <seq.h>

double* seq_double(double* dest, double start, double end, size_t n) {
	if (dest==NULL || n<=0)
		return NULL;
	double step = (n>1) ? (end-start)/(n-1) : 0.0;
	for (size_t i=0; i<n; i++)
		dest[i] = start+i*step;
	return dest;
}

size_t seq_double_step(double* dest, double start, double end, double step, size_t max_size) {
	if (dest==NULL || max_size==0 || (end-start)*step<0)
		return 0;
	size_t n;
	if (step==0)
		n=1;
	else
		n=(end-start)/step;
	size_t n_real = (n<max_size)? n : max_size;
	for (size_t i=0; i<n_real; i++)
		dest[i] = start+i*step;
	return n_real;
}

int* seq_int(int* dest, int start, int end, size_t n) {
	if (dest==NULL || n<=0)
		return NULL;
	int step = (n>1) ? (end-start)/(n-1) : 0;
	for (int i=0; i<n; i++)
		dest[i] = start+i*step;
	return dest;
}

size_t seq_int_step(int* dest, int start, int end, int step, size_t max_size) {
	if (dest==NULL || max_size==0 || (end-start)*step<0)
		return 0;
	size_t n;
	if (step==0)
		n=1;
	else
		n=(end-start)/step;
	size_t n_real = (n<max_size)? n : max_size;
	for (size_t i=0; i<n_real; i++)
		dest[i] = start+i*step;
	return n_real;
}