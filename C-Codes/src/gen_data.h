#ifndef _GEN_DATA_H_INCLUDED_
#define _GEN_DATA_H_INCLUDED_

#include <stdio.h>
#include <math.h>

#include <gsl/gsl_randist.h>

double* gen_data_event_times(double* event_times, size_t nobs, size_t ntimes, double* time, double* hazard, gsl_rng* ran_gen);
void gen_data_event(char* event, double* time, size_t nobs, double* event_times, double tau, double prob_cens_tau, gsl_rng* ran_gen);
size_t* gen_data_region(size_t* reg, size_t nobs, size_t nfrail, gsl_rng* ran_gen);

void gen_data_stepwise_approx(double *val, size_t npoints, double *x, double (*fun)(double x));
void gen_data_rand_walk_cov(double *cov, size_t nobs, size_t ncov, size_t ntimes, double var, double (*cov_lim)[2], gsl_rng *ran_gen);

#endif