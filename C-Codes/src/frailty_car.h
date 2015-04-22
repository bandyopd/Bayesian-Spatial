#ifndef _FRAILTY_CAR_H_INCLUDED_
#define _FRAILTY_CAR_H_INCLUDED_

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

void frailty_car_adj_matr_diag_dist(char *adj_matr, size_t diag_dist, size_t nfrail);

void frailty_car_nadj(size_t *nadj, size_t nfrail, char *adj_matr);
void frailty_car_adj_ind(size_t **adj_ind, size_t nfrail, size_t *nadj, char *adj_matr);

#endif