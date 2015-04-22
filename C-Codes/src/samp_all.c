#include <samp_all.h>

double samp_all_car(size_t iter, size_t pos, size_t nparam, double* param, gsl_rng* ran_gen, size_t nobs, size_t nfrail, 
			size_t ncov, size_t ntimes, double* time_int, size_t* frail_tbl, size_t* nevents, 
			size_t** event_obs, double* atrisk_propsum, double *cov, double* cov_lim, double r0, double c0, 
			size_t* nevcov, size_t** evcov_obs, double* atrisk_sumcov,
			size_t *nevfrail, size_t** evfrail_obs, double *atrisk_sumprop_frail, size_t *nadj, size_t **adj_ind,
			double theta_shape, double theta_scale, double var_proportion) {
	double (*cv)[ncov][ntimes] = (double (*)[ncov][ntimes])cov;
	double (*cvl)[2] = (double (*)[2])cov_lim;
	size_t (*nevcv)[ntimes] = (size_t (*)[ntimes])nevcov;
	size_t* (*evcvo)[ntimes] = (size_t* (*)[ntimes])evcov_obs;
	double (*arsc)[ntimes] = (double (*)[ntimes])atrisk_sumcov;
	size_t (*nevfr)[ntimes] = (size_t (*)[ntimes])nevfrail;
	size_t*  (*evfro)[ntimes] = (size_t* (*)[ntimes])evfrail_obs;
	double (*arspfr)[ntimes] = (double (*)[ntimes])atrisk_sumprop_frail;

	if (pos<ntimes) {
		size_t* evobs = event_obs[pos];
		size_t nev = nevents[pos];
		double psum = atrisk_propsum[pos];
		double* regcoef = param+(ntimes+ncov*pos);
		double* frail = param+(ntimes+ncov*ntimes+nfrail*pos);
		double dt = time_int[pos] - ((pos>0) ? time_int[pos-1] : 0);
		
		double samp = baseline_fullcond_samp(iter, pos, param[pos], ncov, nfrail, ntimes, nev,
				evobs, dt, psum, regcoef, frail, frail_tbl, cv, cvl, r0, c0, ran_gen);
		
		//double samp = param[pos];
		return samp;
	} else if (pos<ntimes*(1+ncov)) { /* regression function */
		size_t jpos = (pos-ntimes) / ncov;
		size_t kpos = (pos-ntimes) % ncov;
		size_t nevc = nevcv[kpos][jpos];
		size_t* evcvobs = evcvo[kpos][jpos];
		double sum_cv_ar = arsc[kpos][jpos];
		double baseline = param[jpos];
		double* regcoef = param+(ntimes+ncov*jpos);
		double* frail = param+(ntimes+ncov*ntimes+nfrail*jpos);
		double dt = time_int[jpos] - ((jpos>0) ? time_int[jpos-1] : 0);
		
		double samp = regfun_fullcond_samp(iter, kpos, jpos, param[pos], ncov, nfrail, ntimes, nevc, evcvobs, dt,
				sum_cv_ar, baseline, regcoef, frail, frail_tbl, cv, cvl, ran_gen);
		
		//double samp = param[pos];
		return samp;
	} else if (pos<ntimes*(1+ncov+nfrail)) { /* frailty */
		size_t jpos = (pos-ntimes*(1+ncov)) / nfrail;
		size_t lpos = (pos-ntimes*(1+ncov)) % nfrail;
		size_t nevf = nevfr[lpos][jpos];
		size_t *evfrobs = evfro[lpos][jpos];
		double psumfr = arspfr[lpos][jpos];
		double baseline = param[jpos];
		double *regcoef = param+(ntimes+ncov*jpos);
		double *frail = param+(ntimes+ncov*ntimes+nfrail*jpos);
		double theta = param[ntimes*(1+ncov+nfrail)+jpos];
		double dt = time_int[jpos] - ((jpos>0) ? time_int[jpos-1] : 0);
		size_t num_adj = nadj[lpos];
		size_t *adj_i = adj_ind[lpos];
		
		double samp;
		
		if (lpos==0) {
			samp = 0;
		} else{
			samp = frail_fullcond_samp(iter, lpos, jpos, param[pos], ncov, nfrail, ntimes, nevf, evfrobs, dt,
				psumfr, num_adj, adj_i, baseline, regcoef, frail, frail_tbl, cv, cvl, theta, var_proportion, ran_gen);
		}
			
		//samp = param[pos];
		return samp;
	} else { /* theta */
		size_t jpos = pos-ntimes*(1+ncov+nfrail);
		double *frail = param+ntimes*(1+ncov)+nfrail*jpos;
		double samp = theta_fullcond_samp(iter, jpos, param[pos], nfrail, ntimes, frail, nadj, adj_ind, theta_shape, theta_scale, var_proportion, ran_gen);
		//double samp = param[pos];
		return samp;
	}
}

double samp_all_car_v(size_t iter, size_t pos, size_t nparam, double* param, gsl_rng* ran_gen, va_list args) {
	size_t nobs = va_arg(args,size_t);
	size_t nfrail = va_arg(args,size_t);
	size_t ncov = va_arg(args,size_t);
	size_t ntimes = va_arg(args,size_t);
	double* time_int = va_arg(args,double*);
	size_t* frail_tbl = va_arg(args,size_t*);
	size_t* nevents = va_arg(args,size_t*);
	size_t** event_obs = va_arg(args,size_t**);
	double* atrisk_propsum = va_arg(args,double*);
	double* cov = va_arg(args,double*);
	double* cov_lim = va_arg(args,double*);
	double r0 = va_arg(args,double);
	double c0 = va_arg(args,double);
	size_t* nevcov = va_arg(args,size_t*);
	size_t** evcov_obs = va_arg(args,size_t**);
	double* atrisk_sumcov = va_arg(args,double*);
	size_t *nevfrail = va_arg(args,size_t*);
	size_t **evfrail_obs = va_arg(args, size_t**);
	double *atrisk_sumprop_frail = va_arg(args, double*);
	size_t *nadj = va_arg(args,size_t*);
	size_t **adj_ind = va_arg(args,size_t**);
	double theta_shape = va_arg(args, double);
	double theta_scale = va_arg(args, double);
	double var_proportion = va_arg(args, double);
	return samp_all_car(iter,pos,nparam,param,ran_gen,nobs,nfrail,ncov, ntimes, time_int, frail_tbl, 
			nevents, event_obs, atrisk_propsum, cov, cov_lim, r0, c0, nevcov, evcov_obs, atrisk_sumcov, 
			nevfrail, evfrail_obs, atrisk_sumprop_frail, nadj, adj_ind,
			theta_shape, theta_scale, var_proportion);
}

double baseline_fullcond_samp(size_t iter, size_t pos, double current, size_t ncov, size_t nfrail, size_t ntimes, size_t nevents, 
		size_t *evobs, double dt, double psum, double *regcoef, double *frail, size_t *frail_tbl, 
		double (*cov)[ncov][ntimes], double (*cov_lim)[2], double r0, double c0, gsl_rng *ran_gen) {
	double* terms;
	if (nevents==0) {
		terms = NULL;
	} else {
		terms = malloc(nevents*sizeof(*terms));
	}
	double covsum;
	for (size_t i=0; i<nevents; i++) {
		covsum=0;
		for (size_t k=0; k<ncov; k++) {
			covsum+=cov[evobs[i]][k][pos]*regcoef[k];
		}
		terms[i] = (covsum+frail[frail_tbl[evobs[i]]]);
	}
	double covmin = 0;
	double lim1;
	double lim2;
	for (size_t k=0; k<ncov; k++) {
		lim1 = regcoef[k]*cov_lim[k][0];
		lim2 = regcoef[k]*cov_lim[k][1];
		covmin += (lim1<lim2) ? lim1 : lim2;
	}
	double frailmin = frail[0];
	for (size_t l=1; l<nfrail; l++) {
		frailmin = (frailmin<frail[l]) ? frailmin : frail[l];
	}
	double constr = -(covmin+frailmin);
	double shape = c0*r0*dt;
	double scale = 1.0/(c0+psum)/dt;
	/*
	if (iter==0) {
		printf("shape: %g\n",shape);
		printf("psum: %g\n",psum);
	}
	*/
	poly_gamma_proposal_param_t par = poly_gamma_proposal_param_max_using_mpfr(nevents, terms, constr, shape, scale);
	/*
	if (iter==0) {
		printf("theta=%g shape=%g scale=%g mean=%g\n",par.theta, par.shape, par.scale, par.shape*par.scale+par.theta);
	}
	*/
	/*
	poly_gamma_proposal_param_t par1 = poly_gamma_proposal_param_integral_using_mpf(nevents, terms, constr, shape, scale);
	if (iter==0) {
		printf("theta1=%g shape1=%g scale1=%g mean1=%g\n",par1.theta, par1.shape, par1.scale, par1.shape*par1.scale+par1.theta);
	}
	*/
	double samp;
	if (iter==0) {
		samp = baslin_samp(current, ran_gen, par);
	} else {
		samp = mcmc_metr_hast_step_ratio(current, baslin_pdf_ratio_v, baslin_proposal_pdf_v, baslin_samp_v, ran_gen, 
		par, constr, nevents, terms, shape, scale);
	}
	free(terms);
	return samp;
}

double regfun_fullcond_samp(size_t iter, size_t cov_pos, size_t time_pos, double current, 
		size_t ncov, size_t nfrail, size_t ntimes, size_t nevcov,
		size_t *evcov_obs, double dt, double sum_cv_ar, double baseline, double *regcoef, 
		double *frail, size_t *frail_tbl, double (*cov)[ncov][ntimes], double (*covlim)[2],
		gsl_rng *ran_gen) {	
	double* terms;
	if (nevcov==0) {
		terms = NULL;
	} else {
		terms = malloc(nevcov*sizeof(*terms));
	}
	double covsum;
	for (size_t i=0; i<nevcov; i++) {
		covsum=0;
		for (size_t k=0; k<ncov; k++) {
			if (k!=cov_pos)
				covsum+=cov[evcov_obs[i]][k][time_pos]*regcoef[k];
		}
		terms[i] = (baseline+covsum+frail[frail_tbl[evcov_obs[i]]])/cov[evcov_obs[i]][cov_pos][time_pos];
	}
	double covmin = 0;
	double lim1;
	double lim2;
	for (size_t k=0; k<ncov; k++) {
		if (k!=cov_pos) {
			lim1 = regcoef[k]*covlim[k][0];
			lim2 = regcoef[k]*covlim[k][1];
			covmin += (lim1<lim2) ? lim1 : lim2;
		}
	}
	double frailmin = frail[0];
	for (size_t l=1; l<nfrail; l++) {
		frailmin = (frailmin<frail[l]) ? frailmin : frail[l];
	}
	double constr = - (baseline+covmin+frailmin);
	double constr1 = (covlim[cov_pos][1]>0) ? constr/covlim[cov_pos][1] : -INFINITY;
	double constr2 = (covlim[cov_pos][0]<0) ? constr/covlim[cov_pos][0] : INFINITY;
	double beta = sum_cv_ar*dt;
	poly_exp_proposal_param_t par = poly_exp_proposal_param_max_using_mpfr(nevcov, terms, constr1, constr2, beta);
	/*
	if (iter==0) {
		printf("sum_cv_ar: %g cov_lim: (%g %g) constr: (%g %g)\n", sum_cv_ar, covlim[cov_pos][0], covlim[cov_pos][1], constr1, constr2);
		if (par.gamma) {
			printf("par.offset: %g par.shape: %g par.scale: %g par.sign %d mean: %g max: %g\n", par.offset, par.shape, par.scale, par.sign, par.shape*par.scale*par.sign+par.offset, (par.shape-1)*par.scale*par.sign+par.offset);
		}
	}
	*/
	double samp;
//	if (iter==0) {
//		samp = regfun_exp_samp_direct(ran_gen, constr1, constr2, par);
//		if (isnan(samp)){
//			samp = mcmc_metr_hast_step(current, regfun_exp_pdf_v, regfun_exp_proposal_pdf_v, regfun_exp_samp_v, ran_gen, 
//				par, nevcov, terms, constr1, constr2, beta);
//		}
//	} else {
		samp = mcmc_metr_hast_step_ratio(current, regfun_exp_pdf_ratio_v, regfun_exp_proposal_pdf_v, regfun_exp_samp_v, ran_gen, 
				par, nevcov, terms, constr1, constr2, beta);
//	}
	free(terms);
	return samp;
}

double frail_fullcond_samp(size_t iter, size_t frail_pos, size_t time_pos, double current, size_t ncov, size_t nfrail, size_t ntimes, 
				size_t nevfrail, size_t *evfrail_obs, double dt,
				double psumfr, size_t nadj, size_t *adj_ind, 
				double baseline, double *regcoef, double *frail, size_t *frail_tbl, double (*cov)[ncov][ntimes], double (*covlim)[2], 
				double theta, double var_proportion, gsl_rng *ran_gen) {
	double *terms;
	if (nevfrail==0) {
		terms = NULL;
	} else {
		terms = malloc(nevfrail*sizeof(*terms));
	}
	double covsum;
	for (size_t i=0; i<nevfrail; i++) {
		covsum=0;
		for (size_t k=0; k<ncov; k++) {
			size_t obs = evfrail_obs[i];
			covsum+=cov[obs][k][time_pos]*regcoef[k];
		}
		terms[i] = baseline+covsum;
	}
	double covmin = 0;
	double lim1;
	double lim2;
	for (size_t k=0; k<ncov; k++) {
		lim1 = regcoef[k]*covlim[k][0];
		lim2 = regcoef[k]*covlim[k][1];
		covmin += (lim1<lim2) ? lim1 : lim2;
	}
	double constr = -(baseline+covmin);
	double adj_sum=0;
	for (size_t i=0; i<nadj; i++) {
		size_t l = adj_ind[i];
		adj_sum += frail[l];
	}
	//PRINT_VECTOR("adj_ind", adj_ind, nadj);
	//printf("psumfr: %g theta: %g\n", psumfr, theta);
	//double adj_mean = adj_sum / nadj;
	//double mean = adj_mean - psumfr*dt/theta/nadj;
	double mean = (adj_sum-psumfr*dt/theta)/(nadj+var_proportion);
	double var = 1/theta/(nadj+var_proportion);
	double stdev = sqrt(var);
	/*
	if (iter==0) {
		printf("constr: %g mean: %g stdev: %g\n", constr, mean, stdev);
	}
	*/
	poly_gaussian_proposal_param_t par = poly_gaussian_proposal_param_max_using_mpfr(nevfrail, terms, constr, INFINITY, mean, stdev);
	/*
	if (iter==0) {
		printf("par.mean: %g par.stdev: %g\n", par.mean, par.stdev);
	}
	*/
	double samp = mcmc_metr_hast_step_ratio(current, regfun_gaussian_pdf_ratio_v, regfun_gaussian_proposal_pdf_v, regfun_gaussian_samp_v, ran_gen,
			par, nevfrail, terms, constr, INFINITY, mean, stdev);
	//double samp = current;
	free(terms);
	return samp;
}

double theta_fullcond_samp(size_t iter, size_t jpos, double theta, size_t nfrail, size_t ntimes, double *frail, 
			size_t *nadj, size_t **adj_ind, double theta_shape, double theta_scale, double var_proportion, gsl_rng *ran_gen) {
	double sumsq = 0;	
	for (size_t l=0; l<nfrail; l++) {
		double frail_l = frail[l];
		for (size_t li=0; li<nadj[l]; li++) {
			size_t ll = adj_ind[l][li];
			if (ll>l) {
				double frail_ll = frail[ll];
				sumsq += (frail_ll-frail_l)*(frail_ll-frail_l);
			}
		}
		sumsq += frail_l*frail_l*var_proportion;
	}
	double shape = theta_shape + (double)nfrail/2;
	double scale = 1.0/(1.0/theta_scale+sumsq/2);
	double samp = gsl_ran_gamma(ran_gen, shape, scale);
	return samp;
}

double baslin_samp(double current, gsl_rng *ran_gen, poly_gamma_proposal_param_t par) {
	double samp = poly_gamma_proposal_samp(ran_gen,par);
	return samp;
}

double baslin_samp_v(double current, gsl_rng *ran_gen, va_list args) {
	poly_gamma_proposal_param_t par = va_arg(args,poly_gamma_proposal_param_t);
	return baslin_samp(current, ran_gen, par);
}

double baslin_samp_args(double current, gsl_rng * ran_gen, void** args) {
	poly_gamma_proposal_param_t* par = (poly_gamma_proposal_param_t*)(args[0]);
	return baslin_samp(current,ran_gen,*par);
}

double baslin_pdf_ratio(double x1, double x2, double constr, size_t nterms, double *terms, double shape, double scale) {
	double pdf = poly_gamma_pdf_ratio_using_mpfr(x1, x2, nterms, terms, constr, shape, scale);
	return pdf;
}
double baslin_pdf_ratio_v(double x1, double x2, va_list args) {
	va_arg(args, poly_gamma_proposal_param_t);
	double constr = va_arg(args,double);
	size_t nterms = va_arg(args,size_t);
	double* terms = va_arg(args,double*);
	double shape = va_arg(args,double);
	double scale = va_arg(args,double);
	return baslin_pdf_ratio(x1, x2, constr, nterms, terms, shape, scale);
}

double baslin_pdf(double x, double constr, size_t nterms, double *terms, double shape, double scale) {
	double pdf = poly_gamma_pdf_using_mpfr(x, nterms, terms, constr, shape, scale);
	return pdf;
}

double baslin_pdf_v(double x, va_list args) {
	va_arg(args, poly_gamma_proposal_param_t);
	double constr = va_arg(args,double);
	size_t nterms = va_arg(args,size_t);
	double* terms = va_arg(args,double*);
	double shape = va_arg(args,double);
	double scale = va_arg(args,double);
	return baslin_pdf(x, constr, nterms, terms, shape, scale);
}

double baslin_pdf1_v(double x, va_list args) {
	double constr = va_arg(args,double);
	size_t nterms = va_arg(args,size_t);
	double* terms = va_arg(args,double*);
	double shape = va_arg(args,double);
	double scale = va_arg(args,double);
	return baslin_pdf(x, constr, nterms, terms, shape, scale);
}

double baslin_pdf_args(double x, void** args) {
	double* constr = (double*)(args[1]);
	size_t* nterms = (size_t*)(args[2]);
	double* terms = (double*)(args[3]);
	double* shape = (double*)(args[4]);
	double* scale = (double*)(args[5]);
	return baslin_pdf(x, *constr, *nterms, terms, *shape, *scale);
}

double baslin_proposal_pdf(double x, double x_old, poly_gamma_proposal_param_t par) {
	return poly_gamma_proposal_pdf_using_mpfr(x,par);
}

double baslin_proposal_pdf_v(double x, double x_old, va_list args) {
	poly_gamma_proposal_param_t par = va_arg(args, poly_gamma_proposal_param_t);
	return baslin_proposal_pdf(x, x_old, par);
}

double baslin_proposal_pdf_args(double x, double x_old, void** args) {
	poly_gamma_proposal_param_t* par = (poly_gamma_proposal_param_t*)(args[0]);
	return baslin_proposal_pdf(x,x_old,*par);
}

double baslin_proposal_pdf1_v(double x, double x_old, va_list args) {
	va_arg(args,double); // constr
	va_arg(args,size_t); // nterms
	va_arg(args,double*); // terms
	double shape = va_arg(args,double);
	double scale = va_arg(args,double);
	
	double sd = sqrt(shape*scale*scale);
	return gsl_ran_gaussian_pdf(x-x_old, sd);
}

double regfun_gaussian_samp(double current, gsl_rng* ran_gen, poly_gaussian_proposal_param_t par) {
	return poly_gaussian_proposal_samp(par, ran_gen);
}

double regfun_gaussian_samp_v(double current, gsl_rng* ran_gen, va_list pass_args) {
	poly_gaussian_proposal_param_t par = va_arg(pass_args, poly_gaussian_proposal_param_t);
	//size_t nterms = va_arg(pass_args, size_t);
	//double *terms = va_arg(pass_args, double*);
	//double c1 = va_arg(pass_args, double);
	//double c2 = va_arg(pass_args, double);
	//double mean = va_arg(pass_args, double);
	//double stdev = va_arg(pass_args, double);
	double samp = regfun_gaussian_samp(current, ran_gen, par);
	return samp;
}

double regfun_gaussian_pdf_ratio(double x1, double x2, size_t nterms, double *terms, double c1, double c2, double mean, double stdev) {
	return poly_gaussian_pdf_ratio_using_mpfr(x1, x2, nterms, terms, c1, c2, mean, stdev);	
}
double regfun_gaussian_pdf_ratio_v(double x1, double x2, va_list pass_args) {
	va_arg(pass_args, poly_gaussian_proposal_param_t); // par
	size_t nterms = va_arg(pass_args, size_t);
	double *terms = va_arg(pass_args, double*);
	double c1 = va_arg(pass_args, double);
	double c2 = va_arg(pass_args, double);
	double mean = va_arg(pass_args, double);
	double stdev = va_arg(pass_args, double);
	double dist = regfun_gaussian_pdf_ratio(x1, x2, nterms, terms, c1, c2, mean, stdev);
	return dist;
}

double regfun_gaussian_pdf(double x, size_t nterms, double *terms, double c1, double c2, double mean, double stdev) {
	return poly_gaussian_pdf_using_mpfr(x, nterms, terms, c1, c2, mean, stdev);	
}

double regfun_gaussian_pdf_v(double x, va_list pass_args) {
	va_arg(pass_args, poly_gaussian_proposal_param_t); // par
	size_t nterms = va_arg(pass_args, size_t);
	double *terms = va_arg(pass_args, double*);
	double c1 = va_arg(pass_args, double);
	double c2 = va_arg(pass_args, double);
	double mean = va_arg(pass_args, double);
	double stdev = va_arg(pass_args, double);
	double dist = regfun_gaussian_pdf(x, nterms, terms, c1, c2, mean, stdev);
	return dist;
}

double regfun_gaussian_proposal_pdf(double x, double x_old, poly_gaussian_proposal_param_t par) {
	return poly_gaussian_proposal_pdf(x, par);
}

double regfun_gaussian_proposal_pdf_v(double x, double x_old, va_list pass_args) {
	poly_gaussian_proposal_param_t par = va_arg(pass_args, poly_gaussian_proposal_param_t);
	//size_t nterms = va_arg(pass_args, size_t);
	//double *terms = va_arg(pass_args, double*);
	//double c1 = va_arg(pass_args, double);
	//double c2 = va_arg(pass_args, double);
	//double mean = va_arg(pass_args, double);
	//double stdev = va_arg(pass_args, double);
	double dist = regfun_gaussian_proposal_pdf(x, x_old, par);
	return dist;
}

double regfun_exp_samp_direct(gsl_rng* ran_gen, double constr1, double constr2, poly_exp_proposal_param_t par) {
	double samp;
	const size_t MAX_TRIALS = 100;
	size_t cnt = 0;
	do {
		if (cnt>=MAX_TRIALS) {
			samp = NAN;
			break;
		} else {
			samp = poly_exp_proposal_samp(par, ran_gen);
			cnt++;
		}
	} while (samp<constr1 || samp>constr2);
	return samp;
}

double regfun_exp_samp(double current, gsl_rng* ran_gen, poly_exp_proposal_param_t par) {
		return poly_exp_proposal_samp(par, ran_gen);
}

double regfun_exp_samp_v(double current, gsl_rng* ran_gen, va_list pass_args) {
	poly_exp_proposal_param_t par = va_arg(pass_args, poly_exp_proposal_param_t);
	//size_t nterms = va_arg(pass_args, size_t);
	//double *terms = va_arg(pass_args, double*);
	//double c1 = va_arg(pass_args, double);
	//double c2 = va_arg(pass_args, double);
	//double beta = va_arg(pass_args, double);
	double samp = regfun_exp_samp(current, ran_gen, par);
	return samp;
}

double regfun_exp_pdf_ratio(double x1, double x2, size_t nterms, double *terms, double c1, double c2, double beta) {
	return poly_exp_pdf_ratio_using_mpfr(x1, x2, nterms, terms, c1, c2, beta);	
}

double regfun_exp_pdf_ratio_v(double x1, double x2, va_list pass_args) {
	va_arg(pass_args, poly_exp_proposal_param_t); // par
	size_t nterms = va_arg(pass_args, size_t);
	double *terms = va_arg(pass_args, double*);
	double c1 = va_arg(pass_args, double);
	double c2 = va_arg(pass_args, double);
	double beta = va_arg(pass_args, double);
	double dist = regfun_exp_pdf_ratio(x1, x2, nterms, terms, c1, c2, beta);
	return dist;
}

double regfun_exp_pdf(double x, size_t nterms, double *terms, double c1, double c2, double beta) {
	return poly_exp_pdf_using_mpfr(x, nterms, terms, c1, c2, beta);	
}

double regfun_exp_pdf_v(double x, va_list pass_args) {
	va_arg(pass_args, poly_exp_proposal_param_t); // par
	size_t nterms = va_arg(pass_args, size_t);
	double *terms = va_arg(pass_args, double*);
	double c1 = va_arg(pass_args, double);
	double c2 = va_arg(pass_args, double);
	double beta = va_arg(pass_args, double);
	double dist = regfun_exp_pdf(x, nterms, terms, c1, c2, beta);
	return dist;
}

double regfun_exp_proposal_pdf(double x, double x_old, poly_exp_proposal_param_t par) {
	return poly_exp_proposal_pdf(x, par);
}

double regfun_exp_proposal_pdf_v(double x, double x_old, va_list pass_args) {
	poly_exp_proposal_param_t par = va_arg(pass_args, poly_exp_proposal_param_t);
	//size_t nterms = va_arg(pass_args, size_t);
	//double *terms = va_arg(pass_args, double*);
	//double c1 = va_arg(pass_args, double);
	//double c2 = va_arg(pass_args, double);
	//double beta = va_arg(pass_args, double);
	double dist = regfun_exp_proposal_pdf(x, x_old, par);
	return dist;
}
