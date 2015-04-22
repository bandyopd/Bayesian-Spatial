#include <frailty_car.h>

void frailty_car_nadj(size_t *nadj, size_t nfrail, char *adj_matr) {
	char (*am)[nfrail] = (char (*)[nfrail])adj_matr;
	for (size_t i=0; i<nfrail; i++) {
		size_t cnt=0;
		for (size_t j=0; j<nfrail; j++) {
			if (j!=i && am[i][j]) {
				cnt++;
			}
		}
		nadj[i] = cnt;
	}
	return;
}

void frailty_car_adj_ind(size_t **adj_ind, size_t nfrail, size_t *nadj, char *adj_matr) {
	char (*am)[nfrail] = (char (*)[nfrail])adj_matr;
	for (size_t i=0; i<nfrail; i++) {
		size_t pos=0;
		for (size_t j=0; j<nfrail; j++) {
			if (j!=i && am[i][j]) {
				adj_ind[i][pos] = j;
				pos++;
			}
		}
	}
	return;
}

void frailty_car_adj_matr_diag_dist(char *adj_matr, size_t diag_dist, size_t nfrail) {
	char (*am)[nfrail] = (char (*)[nfrail])adj_matr;
	for (size_t i=0; i<nfrail; i++) {
		for (size_t j=0; j<nfrail; j++) {
			am[i][j] = (i-j<=diag_dist || j-i<=diag_dist); // 0 or 1, respectively
		}
	}
	return;
}