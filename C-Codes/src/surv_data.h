#ifndef _SURV_DATA_H_INCLUDED_
#define _SURV_DATA_H_INCLUDED_

#include <stdlib.h>
#include <string.h>

typedef struct {
	size_t nobs; 
	double* enter;
	double* time;
	char* event;
	size_t* region; 
} surv_data;

surv_data* surv_data_alloc(size_t nobs);
surv_data* surv_data_calloc(size_t nobs);
void surv_data_free(surv_data* dat);

char* surv_data_atrisk(char* atrisk_matr, size_t ntimes, double* time_vec, surv_data* dat);
char* surv_data_event(char* event_matr, size_t ntimes, double* time_vec, surv_data* dat);
double* surv_data_hazard(double* result, size_t nobs, size_t nfrail, size_t ncov, size_t ntimes, double* baseline, double* cov, double* regfun, double* frail, size_t* frail_tbl);

size_t* surv_data_nevents(size_t* nevents, size_t nobs, size_t ntimes, char* events);
size_t** surv_data_event_obs(size_t** event_obs, size_t nobs, size_t ntimes, char* events);

double* surv_data_atrisk_prop(double* atrisk_prop, size_t nobs, size_t ntimes, double* time_int, double* enter_time, double* event_time);
double* surv_data_atrisk_sumprop(double* atrisk_sumprop, size_t nobs, size_t ntimes, double* atrisk_prop);
 
size_t* surv_data_nevcov(size_t* nevcov, size_t ncov, size_t ntimes, size_t* nevents, size_t** event_obs, double* cov);
size_t** surv_data_evcov_obs(size_t** evcov_obs, size_t ncov, size_t ntimes, size_t* nevents, size_t** event_obs, double* cov);
double* surv_data_atrisk_sumcov(double* atrisk_sumcov, size_t nobs, size_t ncov, size_t ntimes, double* cov, double* atrisk_prop);

void surv_data_nevfrail(size_t *nevfrail, size_t nfrail, size_t ntimes, size_t *nevents, size_t **event_obs, size_t *frail_tbl);
void surv_data_evfrail_obs(size_t **evfrail_obs, size_t nfrail, size_t ntimes, size_t *nevents, size_t **event_obs, size_t *frail_tbl);
void surv_data_atrisk_sumprop_frail(double *atrisk_sumprop_frail, size_t nobs, size_t nfrail, size_t ntimes, double *atrisk_prop, size_t *frail_tbl);

#endif