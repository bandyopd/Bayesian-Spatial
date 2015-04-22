#include <surv_data.h>

surv_data* surv_data_alloc(size_t nobs) {
	surv_data* dat = (surv_data*)malloc(sizeof(surv_data));
	if (dat==NULL)
		return NULL;
	dat->nobs = nobs;
	dat->enter = (double*) malloc(nobs*sizeof(double));
	dat->time = (double*) malloc(nobs*sizeof(double));
	dat->event = (char*) malloc(nobs*sizeof(char));
	dat->region = (size_t*) malloc(nobs*sizeof(size_t));
	if (dat->enter==NULL || dat->time==NULL || dat->event==NULL || dat->region==NULL) {
		free(dat->enter);
		free(dat->time);
		free(dat->event);
		free(dat->region);
		free(dat);
		return NULL;
	}
	return dat;
}

surv_data* surv_data_calloc(size_t nobs) {
	surv_data* dat = surv_data_alloc(nobs);
	if (dat==NULL)
		return NULL;
	memset(dat->enter,0,nobs*sizeof(dat->enter[0]));
	memset(dat->time,0,nobs*sizeof(dat->time[0]));
	memset(dat->event,0,nobs*sizeof(dat->event[0]));
	memset(dat->region,0,nobs*sizeof(dat->region[0]));
	return dat;
}

void surv_data_free(surv_data* dat) {
	if (dat==NULL)
		return;
	free(dat->enter);
	dat->enter=NULL;
	free(dat->time);
	dat->time=NULL;
	free(dat->event);
	dat->event=NULL;
	free(dat->region);
	dat->region=NULL;
	free(dat);
}

char* surv_data_atrisk(char* atrisk_matr, size_t ntimes, double* time_vec, surv_data* dat) {
	if (atrisk_matr==NULL || dat==NULL || time_vec== NULL)
		return NULL;
	char (*ar)[ntimes] = (char (*)[ntimes])atrisk_matr;
	double tt;
	for (size_t i=0; i<dat->nobs; i++) {
		for (size_t j=0; j<ntimes; j++) {
			tt = (j>0) ? time_vec[j-1] : 0;
			ar[i][j] = (dat->enter[i]<time_vec[j]&&dat->time[i]>tt);
		}
	}
	return atrisk_matr;
}

char* surv_data_event(char* event_matr, size_t ntimes, double* time_vec, surv_data* dat) {
	if (event_matr==NULL || dat==NULL || time_vec==NULL)
		return NULL;
	char (*ev)[ntimes] = (char (*)[ntimes])event_matr;
	double tt;
	for (size_t i=0; i<dat->nobs; i++) {
		for (size_t j=0; j<ntimes; j++) {
			tt = (j>0) ? time_vec[j-1] : 0;
			ev[i][j] = (dat->time[i]>tt && dat->time[i]<=time_vec[j] && dat->event[i]) ? 1 : 0;;
		}
	}
	return event_matr;
}

double* surv_data_hazard(double* result, size_t nobs, size_t nfrail, size_t ncov, size_t ntimes, double* baseline, double* cov, double* regfun, double* frail, size_t* frail_tbl) {
	if (result==NULL || baseline==NULL || cov==NULL || regfun==NULL || frail==NULL || frail_tbl==NULL)
		return NULL;
	size_t N = nobs;
	size_t p = ncov;
	size_t m = ntimes;
	
	double (*hz)[ntimes] = (double (*)[ntimes])result;
	double (*cv)[ncov][ntimes] = (double (*)[ncov][ntimes])cov;
	double (*rf)[ntimes] = (double (*)[ntimes])regfun;
	double (*fr)[ntimes] = (double (*)[ntimes])frail;
	
	double cov_term;
	double haz;
	for (size_t i=0; i<N; i++) {
		for (size_t j=0; j<m; j++) {
			cov_term=0;
			for (size_t k=0; k<p; k++) {
				cov_term+=cv[i][k][j]*rf[k][j];
			}
			haz = baseline[j]+cov_term+fr[frail_tbl[i]][j];
			hz[i][j] = haz;
		}
	}
	return result;
}

size_t* surv_data_nevents(size_t* nevents, size_t nobs, size_t ntimes, char* events) {
	if (nevents==NULL || events==NULL || ntimes==0)
		return NULL;
	char (*cp)[ntimes] = (char (*)[ntimes])events;
	size_t ncp;
	for (size_t j=0; j<ntimes; j++) {
		ncp=0;
		for (size_t i=0; i<nobs; i++) {
			if (cp[i][j])
				ncp++;
		}
		nevents[j] = ncp;
	}
	return nevents;
}

size_t** surv_data_event_obs(size_t** event_obs, size_t nobs, size_t ntimes, char* events) {
	if (event_obs==NULL || events==NULL)
		return NULL;
	char (*cp)[ntimes] = (char (*)[ntimes])events;
	size_t pos;
	for (size_t j=0; j<ntimes; j++) {
		pos=0;
		for (size_t i=0; i<nobs; i++) {
			if (cp[i][j]) {
				event_obs[j][pos] = i;
				pos++;
			}
		}
	}
	return event_obs;
}

double* surv_data_atrisk_prop(double* atrisk_prop, size_t nobs, size_t ntimes, double* time_int, double* enter_time, double* event_time) {
	double (*arp)[ntimes] = (double (*)[ntimes])atrisk_prop;
	double start;
	double end;
	double t, tt;
	for (size_t i=0; i<nobs; i++) {
		for (size_t j=0; j<ntimes; j++) {
			t = time_int[j];
			tt = (j>0) ? time_int[j-1] : 0;
			start = (tt>enter_time[i]) ? tt : enter_time[i];
			end = (t<event_time[i]) ? t : event_time[i];
			arp[i][j] = (end>start) ? (end-start)/(t-tt) : 0;
		}
	}
	return atrisk_prop;
}

double* surv_data_atrisk_sumprop(double* atrisk_sumprop, size_t nobs, size_t ntimes, double* atrisk_prop) {
	double (*arp)[ntimes] = (double (*)[ntimes])atrisk_prop;
	double sum;
	for (size_t j=0; j<ntimes; j++) {
		sum = 0;
		for (size_t i=0; i<nobs; i++) {
			sum += arp[i][j];
		}
		atrisk_sumprop[j] = sum;
	}
	return atrisk_sumprop;
}

size_t* surv_data_nevcov(size_t* nevcov, size_t ncov, size_t ntimes, size_t* nevents, size_t** event_obs, double* cov) {
	double (*cv)[ncov][ntimes] = (double (*)[ncov][ntimes])cov;
	size_t (*nevcv)[ntimes] = (size_t (*)[ntimes])nevcov;
	size_t cnt;
	for (size_t k=0; k<ncov; k++) {
		for (size_t j=0; j<ntimes; j++) {
			cnt = 0;
			for (size_t i=0; i<nevents[j]; i++) {
				if (cv[event_obs[j][i]][k][j]!=0)
					cnt++;
			}
			nevcv[k][j] = cnt;
		}
	}
	return nevcov;
}

size_t** surv_data_evcov_obs(size_t** evcov_obs, size_t ncov, size_t ntimes, size_t* nevents, size_t** event_obs, double* cov) {
	double (*cv)[ncov][ntimes] = (double (*)[ncov][ntimes])cov;
	size_t* (*evcvo)[ntimes] = (size_t* (*)[ntimes])evcov_obs;
	size_t pos;
	for (size_t k=0; k<ncov; k++) {
		for (size_t j=0; j<ntimes; j++) {
			pos = 0;
			for (size_t i=0; i<nevents[j]; i++) {
				if (cv[event_obs[j][i]][k][j]!=0) {
					evcvo[k][j][pos] = event_obs[j][i];
					pos++;
				}
			}
		}
	}
	return evcov_obs;
}

double* surv_data_atrisk_sumcov(double* atrisk_sumcov, size_t nobs, size_t ncov, size_t ntimes, double* cov, double* atrisk_prop) {
	double (*cv)[ncov][ntimes] = (double (*)[ncov][ntimes])cov;
	double (*arp)[ntimes] = (double (*)[ntimes])atrisk_prop;
	double (*arsc)[ntimes] = (double (*)[ntimes])atrisk_sumcov;
	double sum;
	for (size_t k=0; k<ncov; k++) {
		for (size_t j=0; j<ntimes; j++) {
			sum = 0;
			for (size_t i=0; i<nobs; i++) {
				sum += cv[i][k][j]*arp[i][j];
			}
			arsc[k][j] = sum;
		}
	}
	return atrisk_sumcov;
}

void surv_data_nevfrail(size_t *nevfrail, size_t nfrail, size_t ntimes, size_t *nevents, size_t **event_obs, size_t *frail_tbl) {
	size_t (*nevfr)[ntimes] = (size_t (*)[ntimes])nevfrail;
	for (size_t l=0; l<nfrail; l++) {
		for (size_t j=0; j<ntimes; j++) {
			nevfr[l][j] = 0;
		}
	}
	for (size_t j=0; j<ntimes; j++) {
		for (size_t i=0; i<nevents[j]; i++) {
			size_t l=frail_tbl[event_obs[j][i]];
			nevfr[l][j]++;
		}
	}
	return;
}

void surv_data_evfrail_obs(size_t **evfrail_obs, size_t nfrail, size_t ntimes, size_t *nevents, size_t **event_obs, size_t *frail_tbl) {
	size_t* (*evfrobs)[ntimes] = (size_t* (*)[ntimes])evfrail_obs;
	size_t *frpos = malloc(nfrail*sizeof(*frpos));
	for (size_t j=0; j<ntimes; j++) {
		for (size_t l=0; l<nfrail; l++) {
			frpos[l] = 0;
		}
		for (size_t i=0; i<nevents[j]; i++) {
			size_t obs = event_obs[j][i];
			size_t l=frail_tbl[obs];
			evfrobs[l][j][frpos[l]] = obs;
			frpos[l]++;
		}
	}
	free(frpos);
	return;
}

void surv_data_atrisk_sumprop_frail(double *atrisk_sumprop_frail, size_t nobs, size_t nfrail, size_t ntimes, double *atrisk_prop, size_t *frail_tbl) {
	double (*ar_spfr)[ntimes] = (double (*)[ntimes])atrisk_sumprop_frail;
	double (*arp)[ntimes] = (double (*)[ntimes])atrisk_prop;
	for (size_t l=0; l<nfrail; l++) {
		for (size_t j=0; j<ntimes; j++) {
			ar_spfr[l][j] = 0;
		}
	}
	for (size_t i=0; i<nobs; i++) {
		for (size_t j=0; j<ntimes; j++) {
			ar_spfr[frail_tbl[i]][j] += arp[i][j];
		}
	}
	return;
}