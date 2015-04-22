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

gsl_rng * ran_gen;

void generate_data_car() {
	const char data_file_name[] = "out/data.dat";
	const char cov_file_name[] = "out/cov.dat";
	const char covlim_file_name[] = "out/cov_lim.dat";
	double tau = 1.0;
	size_t nobs = 10000;
	size_t nfrail = 4;
	size_t ncov = 2;
	size_t ntimes = 10;
	
	FILE *fp;
	
	double baseline_fun(double x) {
		return x;
	}
	
	double reg_fun_1(double x) {
		return -0.3*x;
	}
	double reg_fun_2(double x) {
		return -0.4*x*x;
	}
	
	double (**reg_fun)(double x) = malloc(ncov*sizeof(*reg_fun));
	reg_fun[0] = reg_fun_1;
	reg_fun[1] = reg_fun_2;
	double tmin = tau/ntimes;
	double (*time) = malloc(ntimes*sizeof(*time));
	seq_double(time,tmin,tau,ntimes);
	double (*baseline) = malloc(ntimes*sizeof(*baseline));
	gen_data_stepwise_approx(baseline, ntimes, time, baseline_fun);
	double (*regfun)[ntimes] = malloc(ncov*sizeof(*regfun));
	for (size_t k=0; k<ncov; k++) {
		gen_data_stepwise_approx((double*)regfun[k], ntimes, time, reg_fun[k]);
	}
	double (*cov_lim)[2] = malloc(ncov*sizeof(*cov_lim));
	for (size_t k=0; k<ncov; k++) {
		cov_lim[k][0] = -1;
		cov_lim[k][1] = +1;
	}
	fp = fopen(covlim_file_name, "w");
	if (!fp) {
		fprintf(stderr, "Can not create file: %s\n", covlim_file_name);
	} else {
		fprintf(fp, "%zd\n", ncov);
		for (size_t k=0; k<ncov; k++) {
			fprintf(fp, "%.5g\t%.5g\n", cov_lim[k][0], cov_lim[k][1]);
		}
		fclose(fp);
	}
	
	double (*cov)[ncov][ntimes] = malloc(nobs*sizeof(*cov));
	gen_data_rand_walk_cov((double*)cov, nobs, ncov, ntimes, tmin, cov_lim, ran_gen);
	fp = fopen(cov_file_name, "w");
	if (!fp) {
		fprintf(stderr, "Can not create file: %s\n", cov_file_name);
	} else {
		fprintf(fp, "%zd\t%zd\t%zd\n", nobs, ncov, ntimes);
		for (size_t i=0; i<nobs; i++) {
			for (size_t k=0; k<ncov; k++) {
				for (size_t j=0; j<ntimes; j++) {
					fprintf(fp, "%.5g\t", cov[i][k][j]);
				}
				fprintf(fp, "\n");
			}
		}
		fclose(fp);
	}
	
	char (*adj_matr)[nfrail] = malloc(nfrail*sizeof(*adj_matr));
	frailty_car_adj_matr_diag_dist((char*)adj_matr, 1, nfrail);
	
	free(time);
	free(baseline);
	free(regfun);
	free(cov_lim);
	free(cov);
	free(reg_fun);
	return;
}

void gen_data_test() {
	size_t nobs = 10000;
	size_t nfrail = 5;
	size_t ncov = 2;
	size_t ntimes = 10;
	
	double tau = 1;
	double prob_gtau = 0.5;
	double prob_cens_tau = 0.7;
	double FRAIL_PROP = -0.1;
	
	size_t niter = 100;
	size_t burnin = niter/4;
	
	double rate = -log(prob_gtau)/tau;
	
	printf("Observations: %zu\n",nobs);
	
	double time_int[ntimes];
	double tmin = tau/ntimes;
	seq_double(time_int,tmin,tau,ntimes);
	PRINT_VECTOR("time_int points",time_int,ntimes);
	
	double *baseline = malloc(ntimes*sizeof(*baseline));
	for (size_t j=0; j<ntimes; j++) {
		//baseline[j] = rate*(1+0.5*j/ntimes);
		double tt = (j>0)? time_int[j-1] : 0;
		double t = time_int[j];
		baseline[j] = 1 + (t+tt)/2;
	}
	
	PRINT_VECTOR("baseline",baseline,ntimes);
	
	double cov_lim[ncov][2];
	for (size_t k=0; k<ncov; k++) {
		cov_lim[k][0] = -1;
		cov_lim[k][1] = 1;
	}
	
	double regfun[ncov][ntimes];
	for (size_t k=0; k<ncov; k++) {
		for (size_t j=0; j<ntimes; j++) {
			/*
			double tt = (j>0)? time_int[j-1] : 0;
			double t = time_int[j];
			if (k==0) {
				regfun[k][j] = -0.3*(t+tt)/2;
			} else if (k==1){
				regfun[k][j] = -0.4*(t+tt)/2*(t+tt)/2;
			} else {
				regfun[k][j] = (-1.0+gsl_rng_uniform(ran_gen)*2.0)*baseline[j]*(1.0-FRAIL_PROP)/(MAX(cov_lim[k][1],fabs(cov_lim[k][0])))/(double)ncov;
			}
			*/
			regfun[k][j] = 0;
		}
	}
	PRINT_MATRIX("Regression functions", (double*)regfun, ncov, ntimes);
	
	double (*frail)[ntimes] = malloc(nfrail*sizeof(*frail));
	for (size_t j=0; j<ntimes; j++) {
		//double sum = 0;
		/*
		for (size_t l=0; l<nfrail; l++) {
			if (j>0) {
				frail[l][j] = frail[l][0];
				continue;
			}
			double res;
			do {
				if (l==0) {
					res = gsl_rng_uniform(ran_gen)*baseline[j]*FRAIL_PROP*2;
				} else {
					res = frail[0][j]*(0.5+gsl_rng_uniform(ran_gen));
				}
			} while (res<=-0.3);
			frail[l][j] = res;
			//sum += frail[l][j];
		}
		//double avg = sum / nfrail;
		for (size_t l=1; l<nfrail; l++) {
			//frail[l][j] -= avg;
			frail[l][j] -= frail[0][j];
		}
		frail[0][j] = 0;
		*/
		for (size_t l=0; l<nfrail; l++) {
			frail[l][j] = 0;
		}
	}
	PRINT_MATRIX("Frailty", (double*)frail, nfrail, ntimes);
	
	double (*cov)[ncov][ntimes] = malloc(nobs*sizeof(*cov));
	for (size_t i=0; i<nobs; i++) {
		for (size_t k=0; k<ncov; k++) {
			for (size_t j=0; j<ntimes; j++) {
				cov[i][k][j] = cov_lim[k][0]+gsl_rng_uniform_pos(ran_gen)*(cov_lim[k][1]-cov_lim[k][0]);
			}
		}
	}
	
	surv_data *data = surv_data_calloc(nobs);
	gen_data_region(data->region, nobs, nfrail, ran_gen);
	
	double (*hazard)[ntimes] = malloc(nobs*sizeof(*hazard));
	surv_data_hazard((double*)hazard, nobs, nfrail, ncov, ntimes, baseline, (double*)cov, (double*)regfun, (double*)frail, data->region);
	
	double *event_times = malloc(nobs*sizeof(*event_times));
	gen_data_event_times(event_times, nobs, ntimes,time_int, (double*)hazard, ran_gen);
	gen_data_event(data->event, data->time, nobs, event_times, tau, prob_cens_tau, ran_gen);
	
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
	printf("Total events: %zu\n",nevents_total);
	PRINT_VECTOR_FORMAT("nevents",nevents,ntimes,"%zu ");
	
	
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
	
	double r0 = (double)nevents_total/nobs/tau;
	double c0 = 0.001;
	
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
	
	size_t diag_dist = 1;
	char (*adj_matr)[nfrail] = malloc(nfrail*sizeof(*adj_matr));
	frailty_car_adj_matr_diag_dist((char*)adj_matr, diag_dist, nfrail);
	PRINT_MATRIX("adj_matr", (char*)adj_matr, nfrail, nfrail);
	size_t *nadj = malloc(nfrail*sizeof(*nadj));
	frailty_car_nadj(nadj, nfrail, (char*)adj_matr);
	size_t* (*adj_ind) = malloc(nfrail*sizeof(*adj_ind));
	for (size_t l=0; l<nfrail; l++) {
		adj_ind[l] = malloc(nadj[l]*sizeof(*(adj_ind[l])));
	}
	frailty_car_adj_ind(adj_ind, nfrail, nadj, (char*)adj_matr);
	
	double frail_stdev = (double)nevents_total/nobs/tau/3;
	double theta_mean = 1.0/frail_stdev/frail_stdev;
	//double theta_mean = 1.0e+04;
	double theta_scale = 1000; // big variance
	double theta_shape = theta_mean/theta_scale;
	double var_proportion = 0.01;
	
	size_t nparam = ntimes*(1+ncov+nfrail+1);
	double *init = malloc(nparam*sizeof(*init));
	for (size_t i=0; i<nparam; i++) {
		if (i<ntimes) { //baseline
			init[i] = (double)nevents[i]/natrisk[i]/(time_int[i]-(i>0?time_int[i-1]:0));
			//init[i] = baseline[i];
		} else if (i<ntimes*(1+ncov)) { // regcoef
			//init[i] = regfun[(i-ntimes)%ncov][(i-ntimes)/ncov];
			init[i] = 0;
		} else if (i<ntimes*(1+ncov+nfrail)) { // frail
			//init[i] = frail[(i-ntimes*(1+ncov))%nfrail][(i-ntimes*(1+ncov))/nfrail];
			init[i] = 0;
		} else { // theta
			init[i] = theta_mean;
		}
	}
	
	PRINT_VECTOR("Init",init,nparam);
	
	double (*samp)[nparam] = malloc(niter*sizeof(*samp));
	printf("Sampling...\n");
	clock_t clock_start, clock_end;
	clock_start = clock();
	mcmc_gibbs_samp((double*)samp, niter, nparam, init, samp_all_car_v, ran_gen, 
			nobs, nfrail, ncov, ntimes, time_int, data->region, nevents, event_obs, atrisk_sumprop, (double*)cov, (double*)cov_lim, r0, c0, 
			(size_t*)nevcov, (size_t**)evcov_obs, (double*)atrisk_sumcov,
			(size_t*)nevfrail, (size_t**)evfrail_obs, (double*)atrisk_sumprop_frail, nadj, adj_ind,
			theta_shape, theta_scale, var_proportion);
	clock_end = clock();
	printf("Iterations: %zu\n", niter);
	printf("Sampling time: %g\n",(double)(clock_end-clock_start)/CLOCKS_PER_SEC);
	
	double *avg = malloc(nparam*sizeof(*avg));
	for (size_t j=0; j<nparam; j++) {
		avg[j] = 0;
	}
	double cnt = 0;
	for (size_t i=burnin; i<niter; i++) {
		for (size_t j=0; j<nparam; j++) {
			avg[j] += samp[i][j];
		}
		cnt++;
	}
	for (size_t j=0; j<nparam; j++) {
		avg[j] /= cnt;
	}
	PRINT_VECTOR("Average", avg, nparam);
	
	FILE *fp;
	fp = fopen("out/data", "w");
	fprintf(fp, "Number of time points: %zu\n", ntimes);
	fprintf(fp, "Time points:\n");
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp, "%g ", time_int[j]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "Number of observations: %zu\n", nobs);
	fprintf(fp, "Total number of events: %zu\n", nevents_total);
	fprintf(fp, "Number of events at time intervals:\n");
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp, "[t%zu]\t",j);
	}
	fprintf(fp, "\n");
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp, "%zu\t", nevents[j]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "The sum of proportions of time the individuals are at risk:\n");
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp, "[t%zu]\t",j);
	}
	fprintf(fp, "\n");
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp, "%g\t", atrisk_sumprop[j]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "Frailty adjustment matrix:\n");
	fprintf(fp, "  \t");
	for (size_t l=0; l<nfrail; l++) {
		fprintf(fp, "[%zu]\t", l);
	}
	fprintf(fp, "\n");
	for (size_t l=0; l<nfrail; l++) {
		fprintf(fp, "[%zu]\t", l);
		for (size_t ll=0; ll<nfrail; ll++) {
			fprintf(fp, " %d \t", adj_matr[l][ll]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"Observations:\n");
	fprintf(fp,"[obs#]\tenter\ttime\tevent\tregion\n");
	for (size_t i=0; i<nobs; i++) {
		fprintf(fp,"[%zu]\t%.3g\t%.3g\t%d\t%zu\n", i, data->enter[i], data->time[i], data->event[i], data->region[i]);
	}
	fprintf(fp, "Covariate limits:\n");
	fprintf(fp, "[cov#]\tlower\tupper\n");
	for (size_t k=0; k<ncov; k++) {
		fprintf(fp, "[%zu]\t%g\t%g\n", k, cov_lim[k][0], cov_lim[k][1]);
	}
	fprintf(fp, "Covariates:\n");
	fprintf(fp, "[cov#]\t");
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp, "[t%zu]\t", j);
	}
	fprintf(fp, "\n");
	for (size_t i=0; i<nobs; i++) {
		fprintf(fp, "obs=%zu\n", i);
		for (size_t k=0; k<ncov; k++) {
			fprintf(fp, "[%zu]\t", k);
			for (size_t j=0; j<ntimes; j++) {
				fprintf(fp, "%g\t", cov[i][k][j]);
			}
			fprintf(fp, "\n");
		}
	}
	fclose(fp);
	
	fp = fopen("out/samp", "w");
	fprintf(fp, "Baseline:\n");
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp,"[t%zu]\t", j);
	}
	fprintf(fp,"\n");
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp, "%g\t", baseline[j]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "Regression functions:\n");
	fprintf(fp, "[cov#]\t");
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp,"[t%zu]\t", j);
	}
	fprintf(fp,"\n");
	for (size_t k=0; k<ncov; k++) {
		fprintf(fp, "[%zu]\t", k);
		for (size_t j=0; j<ntimes; j++) {
			fprintf(fp, "%g\t", regfun[k][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "Frailties:\n");
	fprintf(fp, "[fr#]\t");
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp,"[t%zu]\t", j);
	}
	fprintf(fp,"\n");
	for (size_t l=0; l<nfrail; l++) {
		fprintf(fp, "[%zu]\t", l);
		for (size_t j=0; j<ntimes; j++) {
			fprintf(fp, "%g\t", frail[l][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "Number of iterations: %zu\n", niter);
	fprintf(fp, "Burn-in: %zu\n", burnin);
	fprintf(fp, "Number of parameters: %zu\n", nparam);
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp, "bl(t%zu)\t", j);
	}
	for (size_t j=0; j<ntimes; j++) {
		for (size_t k=0; k<ncov; k++) {
			fprintf(fp, "rf[%zu](t%zu)\t", k, j);
		}
	}
	for (size_t j=0; j<ntimes; j++) {
		for (size_t l=0; l<nfrail; l++) {
			fprintf(fp, "fr[%zu](t%zu)\t", l, j);
		}
	}
	for (size_t j=0; j<ntimes; j++) {
		fprintf(fp, "theta[%zu]\t", j);
	}
	fprintf(fp, "\n");
	fprintf(fp, "Initial values:\n");
	for (size_t i=0; i<nparam; i++) {
		fprintf(fp, "%g\t", init[i]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "Real values:\n");
	for (size_t i=0; i<nparam; i++) {
		if (i<ntimes) {
			size_t jpos = i;
			fprintf(fp, "%g\t", baseline[jpos]);
		} else if (i<ntimes*(1+ncov)) {
			size_t jpos = (i-ntimes)/ncov;
			size_t kpos = (i-ntimes)%ncov;
			fprintf(fp, "%g\t", regfun[kpos][jpos]);
		} else if (i<ntimes*(1+ncov+nfrail)) {
			size_t jpos = (i-ntimes*(1+ncov))/nfrail;
			size_t lpos = (i-ntimes*(1+ncov))%nfrail;
			fprintf(fp, "%g\t", frail[lpos][jpos]);
		} else {
			//size_t jpos = i-ntimes*(1+ncov+nfrail);
			fprintf(fp, "---\t");
		}
	}
	fprintf(fp, "\n");
	fprintf(fp, "Mean:\n");
	for (size_t i=0; i<nparam; i++) {
		fprintf(fp, "%g\t", avg[i]);
	}
	fprintf(fp, "\n");
	fprintf(fp, "Bias:\n");
	for (size_t i=0; i<nparam; i++) {
		if (i<ntimes) {
			size_t jpos = i;
			fprintf(fp, "%g\t", avg[i] - baseline[jpos]);
		} else if (i<ntimes*(1+ncov)) {
			size_t jpos = (i-ntimes)/ncov;
			size_t kpos = (i-ntimes)%ncov;
			fprintf(fp, "%g\t", avg[i] - regfun[kpos][jpos]);
		} else if (i<ntimes*(1+ncov+nfrail)) {
			size_t jpos = (i-ntimes*(1+ncov))/nfrail;
			size_t lpos = (i-ntimes*(1+ncov))%nfrail;
			fprintf(fp, "%g\t", avg[i] - frail[lpos][jpos]);
		} else {
			//size_t jpos = i-ntimes*(1+ncov+nfrail);
			fprintf(fp, "---\t");
		}
	}
	fprintf(fp, "\n");
	fprintf(fp, "Sample:\n");
	fprintf(fp, "[iter#]\tValues of parameters\n");
	for (size_t i=0; i<niter; i++) {
		fprintf(fp, "[%zu]\t", i);
		for (size_t p=0; p<nparam; p++) {
			fprintf(fp, "%g\t", samp[i][p]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	
	/* FREEING THE MEMORY */
	free(baseline);
	free(cov);
	free(frail);
	free(hazard);
	free(event_times);
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
	free(adj_matr);
	free(nadj);
	for (size_t l=0; l<nfrail; l++) {
		free(adj_ind[l]);
	}
	free(adj_ind);
	surv_data_free(data);
	free(samp);
	free(init);
	free(avg);
	return;
}

void ran_gen_init(){
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	ran_gen = gsl_rng_alloc (T);
	gsl_rng_set(ran_gen,time(NULL));
	return;
}


int main() {
	ran_gen_init();
	gen_data_test();
	//generate_data_car();
	//prop_dist_norm_test();
	//cdf_gaussian_test();
	gsl_rng_free(ran_gen);
	return 0;
}
