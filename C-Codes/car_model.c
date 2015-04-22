#include<R.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <vec_mat_print.h>
#include <mpfr_array.h>

#include <mcmc.h>
#include <samp_all.h>
#include <surv_data.h>
#include <gen_data.h>
#include <frailty_car.h>
#include <seq.h>

gsl_rng *ran_gen_init(){
	static gsl_rng* static_ran_gen=NULL;
	static const gsl_rng_type * T;
	if (static_ran_gen==NULL){
		gsl_rng_env_setup();
		T = gsl_rng_default;
		static_ran_gen = gsl_rng_alloc (T);
		gsl_rng_set(static_ran_gen,time(NULL));
	}
	return static_ran_gen;
}

void compute_hazard(double *hazard, int *n_obs, int *n_cov, int *n_frail, int *n_times, double *baseline, double *cov, double *regfun, double *frail, int *frail_tbl) {
	size_t nobs = *n_obs;
	size_t ncov = *n_cov;
	size_t nfrail = *n_frail;
	size_t ntimes = *n_times;
	size_t *frtbl = malloc(nobs*sizeof(*frtbl));
	//char str[100];
	for (size_t i=0; i<nobs; i++) {
		frtbl[i] = frail_tbl[i];
	}
	surv_data_hazard(hazard, nobs, nfrail, ncov, ntimes, baseline, cov, regfun, frail, frtbl);
	free(frtbl);
	return;
}

void gen_stepwise_events(double *event_times, int *n_obs, int *n_times, double *time_int, double *hazard) {
	size_t nobs = *n_obs;
	size_t ntimes = *n_times;
	gsl_rng * ran_gen = ran_gen_init();
	gen_data_event_times(event_times, nobs, ntimes, time_int, hazard, ran_gen);
	return;
}

void gen_cens_events(double *time, int *event, int *n_obs, double *event_times, double *end_of_study, double *prob_of_cens) {
	gsl_rng * ran_gen = ran_gen_init();
	size_t nobs = *n_obs;
	double tau = *end_of_study;
	double prob = *prob_of_cens;
	char *ev = malloc(nobs*sizeof(*ev));
	gen_data_event(ev, time, nobs, event_times, tau, prob, ran_gen);
	for (size_t i=0; i<nobs; i++) {
		event[i] = ev[i];
	}
	free(ev);
	return;
}

void fit_car_model(double *sample, // result
			int *n_iter, int *n_obs, int *n_cov, int *n_frail, int *n_times, // sizes
			double *enter_times, double *event_times, int *event_ind, int *regions, // data
			double *time_int, double *covariates, double *cov_lim, int *adj_matrix, // additional parameters from data
			double *p_r0, double *p_c0, double *alpha, double *beta // parameters of priors) 
			){
	size_t niter = *n_iter;
	size_t nobs = *n_obs;
	size_t ncov = *n_cov;
	size_t nfrail = *n_frail;
	size_t ntimes = *n_times;
	double r0 = *p_r0;
	double c0 = *p_c0;
	double (*cov)[ncov][ntimes] = (double (*)[ncov][ntimes])covariates;
	double theta_shape = *alpha;
	double theta_scale = 1.0/(*beta);
	
	gsl_set_error_handler_off ();
	
	gsl_rng * ran_gen = ran_gen_init();
	
	surv_data *data = surv_data_alloc(nobs);
	for (size_t i=0; i<nobs; i++) {
		data->enter[i] = enter_times[i];
		data->time[i] = event_times[i];
		data->event[i] = event_ind[i];
		data->region[i] = regions[i];
	}
	
	char (*adj_matr)[nfrail] = malloc(nfrail*sizeof(*adj_matr));
	for (size_t l=0; l<nfrail; l++) {
		for (size_t ll=0; ll<nfrail;ll++) {
			adj_matr[l][ll] = adj_matrix[nfrail*l+ll];
		}
	}
	
	char (*event_matr)[ntimes] = (char (*)[ntimes])malloc(nobs*sizeof(*event_matr));
	char (*atrisk_matr)[ntimes] = (char (*)[ntimes])malloc(nobs*sizeof(*atrisk_matr));
	surv_data_event((char*)event_matr, ntimes, time_int, data);
	surv_data_atrisk((char*)atrisk_matr, ntimes, time_int, data);
	
	size_t *nevents = malloc(ntimes*sizeof(*nevents));
	size_t *natrisk = malloc(ntimes*sizeof(*natrisk));
	surv_data_nevents(nevents, nobs, ntimes, (char*)event_matr);
	surv_data_nevents(natrisk, nobs, ntimes, (char*)atrisk_matr);
	size_t nevents_total = 0;
	for (size_t j=0; j<ntimes; j++) {
		nevents_total += nevents[j];
	}
	double (*atrisk_prop)[ntimes] = malloc(nobs*sizeof(*atrisk_prop));
	surv_data_atrisk_prop((double*)atrisk_prop, nobs, ntimes, time_int, data->enter, data->time);
	double *atrisk_sumprop = malloc(ntimes*sizeof(*atrisk_sumprop));
	surv_data_atrisk_sumprop(atrisk_sumprop, nobs, ntimes, (double*)atrisk_prop);
	double (*atrisk_sumprop_frail)[ntimes] = malloc(nfrail*sizeof(*atrisk_sumprop_frail));
	surv_data_atrisk_sumprop_frail((double*)atrisk_sumprop_frail, nobs, nfrail, ntimes, (double*)atrisk_prop, data->region);
	
	size_t* *event_obs = malloc(ntimes*sizeof(*event_obs));
	for (size_t j=0; j<ntimes; j++) {
		event_obs[j] = malloc(nevents[j]*sizeof(event_obs[j]));
	}
	surv_data_event_obs(event_obs, nobs, ntimes, (char*)event_matr);
	
	size_t (*nevcov)[ntimes] = malloc(ncov*sizeof(*nevcov));
	surv_data_nevcov((size_t*)nevcov, ncov, ntimes, nevents, event_obs, (double*)cov);
	
	size_t* (*evcov_obs)[ntimes] = malloc(ncov*sizeof(*evcov_obs));
	for (size_t k=0; k<ncov; k++) {
		for (size_t j=0; j<ntimes; j++) {
			evcov_obs[k][j] = malloc(nevcov[k][j]*sizeof(*(evcov_obs[k][j])));
		}
	}
	surv_data_evcov_obs((size_t**)evcov_obs, ncov, ntimes, nevents, event_obs, (double*)cov);
	
	double (*atrisk_sumcov)[ntimes] = malloc(ncov*sizeof(*atrisk_sumcov));
	surv_data_atrisk_sumcov((double*)atrisk_sumcov, nobs, ncov, ntimes, (double*)cov, (double*)atrisk_prop);
	
	size_t (*nevfrail)[ntimes] = malloc(nfrail*sizeof(*nevfrail));
	surv_data_nevfrail((size_t*)nevfrail, nfrail, ntimes, nevents, event_obs, data->region);
	
	size_t* (*evfrail_obs)[ntimes] = malloc(nfrail*sizeof(*evfrail_obs));
	for (size_t l=0; l<nfrail; l++) {
		for (size_t j=0; j<ntimes; j++) {
			evfrail_obs[l][j] = malloc(nevfrail[l][j]*sizeof(*(evfrail_obs[l][j])));
		}
	}
	surv_data_evfrail_obs((size_t**)evfrail_obs, nfrail, ntimes, nevents, event_obs, data->region);
	size_t *nadj = malloc(nfrail*sizeof(*nadj));
	frailty_car_nadj(nadj, nfrail, (char*)adj_matr);
	size_t* (*adj_ind) = malloc(nfrail*sizeof(*adj_ind));
	for (size_t l=0; l<nfrail; l++) {
		adj_ind[l] = malloc(nadj[l]*sizeof(*(adj_ind[l])));
	}
	frailty_car_adj_ind(adj_ind, nfrail, nadj, (char*)adj_matr);
	
	size_t nparam = ntimes*(1+ncov+nfrail+1);
	double *init = malloc(nparam*sizeof(*init));
	for (size_t i=0; i<nparam; i++) {
		if (i<ntimes) { //baseline
			init[i] = (double)nevents[i]/natrisk[i]/(time_int[i]-(i>0?time_int[i-1]:0));
		} else if (i<ntimes*(1+ncov)) { // regcoef
			init[i] = 0;
		} else if (i<ntimes*(1+ncov+nfrail)) { // frail
			init[i] = 0;
		} else { // theta
			init[i] = theta_shape*theta_scale;
		}
	}
	
	mcmc_gibbs_samp(sample, niter, nparam, init, samp_all_car_v, ran_gen, 
			nobs, nfrail, ncov, ntimes, time_int, data->region, nevents, event_obs, atrisk_sumprop, covariates, cov_lim, r0, c0, 
			(size_t*)nevcov, (size_t**)evcov_obs, (double*)atrisk_sumcov,
			(size_t*)nevfrail, (size_t**)evfrail_obs, (double*)atrisk_sumprop_frail, nadj, adj_ind,
			theta_shape, theta_scale, 0);
	free(adj_matr);
	free(event_matr);
	free(atrisk_matr);
	free(nevents);
	free(natrisk);
	free(atrisk_prop);
	free(atrisk_sumprop);
	free(atrisk_sumcov);
	free(nevcov);
	for (size_t j=0; j<ntimes; j++) {
		free(event_obs[j]);
	}
	free(event_obs);
	for (size_t k=0; k<ncov; k++) {
		for (size_t j=0; j<ntimes; j++) {
			free(evcov_obs[k][j]);
		}
	}
	free(evcov_obs);
	free(nevfrail);
	for (size_t l=0; l<nfrail; l++) {
		for (size_t j=0; j<ntimes; j++) {
			free(evfrail_obs[l][j]);
		}
	}
	free(evfrail_obs);
	free(nadj);
	for (size_t l=0; l<nfrail; l++) {
		free(adj_ind[l]);
	}
	free(adj_ind);
	surv_data_free(data);
	free(init);
	return;
}
