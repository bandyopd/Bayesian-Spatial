#include <gen_data.h>

double* gen_data_event_times(double* event_times, size_t nobs, size_t ntimes, double* time, double* hazard, gsl_rng* ran_gen) {
	if( event_times==NULL || time==NULL || hazard==NULL || ran_gen==NULL)
		return NULL;
	double (*haz)[ntimes] = (double (*)[ntimes]) hazard;
	double surv;
	double surv_eval;
	double dt;
	size_t jq;
	double cum_haz;
	for (size_t i=0; i<nobs; i++) {
		surv = gsl_rng_uniform_pos(ran_gen);
		surv_eval = 1;
		jq = ntimes;
		cum_haz = 0;
		for (size_t j=0; j<ntimes; j++) {
			dt = time[j] - ((j>0) ? time[j-1] : 0);
			surv_eval *= exp(-haz[i][j]*dt);
			if (surv_eval<=surv) {
				jq = j;
				break;
			}
			cum_haz += haz[i][j]*dt;
		}
		if (surv_eval>surv) {
			event_times[i] = INFINITY;
		} else {
			event_times[i] = 1/haz[i][jq]*(-log(surv)-cum_haz)+((jq>0)?time[jq-1]:0);
		}
	}
	return event_times;
}

void gen_data_event(char* event, double* time, size_t nobs, double* event_times, double tau, double prob_cens_tau, gsl_rng* ran_gen) {
	double rate = -log(prob_cens_tau)/tau;
	double cens_time;
	for (size_t i=0; i<nobs; i++) {
		cens_time = gsl_ran_exponential(ran_gen, 1/rate);
		cens_time = (cens_time<tau) ? cens_time : tau;
		event[i] = (cens_time>event_times[i]);
		time[i] = (event[i]) ? event_times[i] : cens_time;
	}
	return;
}

size_t* gen_data_region(size_t* reg, size_t nobs, size_t nfrail, gsl_rng* ran_gen) {
	for (size_t i=0; i<nobs; i++) {
		reg[i] = gsl_rng_uniform_int(ran_gen, nfrail);
	}
	return reg;
}

void gen_data_stepwise_approx(double *val, size_t npoints, double *x, double (*fun)(double x)) {
	for (size_t i=0; i<npoints; i++) {
		double xprev = (i>0)?x[i-1]:0;
		val[i] = fun((x[i]+xprev)/2);
	}
	return;
}

void gen_data_rand_walk_cov(double *cov, size_t nobs, size_t ncov, size_t ntimes, double var, double (*cov_lim)[2], gsl_rng *ran_gen) {
	double (*cv)[ncov][ntimes] = (double (*)[ncov][ntimes])cov;
	for (size_t i=0; i<nobs; i++) {
		for (size_t k=0; k<ncov; k++) {
			double low = cov_lim[k][0];
			double up = cov_lim[k][1];
			cv[i][k][0] = low+gsl_rng_uniform(ran_gen)*(up-low);
			for (size_t j=1; j<ntimes; j++) {
				double prev = cv[i][k][j-1];
				double next;
				do {
					double step = gsl_ran_gaussian(ran_gen,sqrt(var));
					next = prev + step;
				} while (next<low || next>up);
				cv[i][k][j] = next;
			}
		}
	}
	return;
}
