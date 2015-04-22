#include<mcmc.h>

double* mcmc_gibbs_samp(double* samp, size_t niter, size_t nparam, double* init,  double (*samp_fun)(size_t iter,size_t pos,size_t nparam, double* param, gsl_rng* ran_gen, va_list pass_args), gsl_rng* ran_gen,...) {
	va_list pass_args;
	va_start(pass_args,ran_gen);
	double* result = mcmc_gibbs_samp_v(samp,niter,nparam,init,samp_fun,ran_gen, pass_args);
	va_end(pass_args);
	return result;
}

//#include<R.h>
static inline void loadBar(size_t x, size_t n, size_t r, size_t w) {
	// Only update r times.
	if (n/r != 0) {
		if ( x % (n/r) != 0 ) return;
	}
	
	// Calculuate the ratio of complete-to-incomplete.
	float ratio = x/(float)n;
	size_t c = ratio * w;
	
	// Show the percentage complete.
	fprintf(stderr,"%3zu%% [", (size_t)(ratio*100) );
	
	// Show the load bar.
	for (size_t x=0; x<c; x++) {
	    	fprintf(stderr,"=");
	    }
	    
	for (size_t x=c; x<w; x++) {
		fprintf(stderr," ");
	}
	fprintf(stderr,"]\r"); // Move to the first column
	fflush(stderr);
	//fprintf(stderr,"]\n\033[F");
}

double* mcmc_gibbs_samp_v(double* samp, size_t niter, size_t nparam, double* init,  double (*samp_fun)(size_t iter,size_t pos,size_t nparam,double* param,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen,va_list pass_args) {
	if (samp==NULL || init==NULL || samp_fun==NULL || ran_gen==NULL)
		return NULL;
	double (*smp)[nparam] = (double (*)[nparam])samp;
	double current[nparam];
	memcpy(current,init,nparam*sizeof(double));
	loadBar(0, niter, 100, 50);
	for (size_t i=0; i<niter; i++) {
		for (size_t j=0; j<nparam; j++) {
			va_list args;
			va_copy (args, pass_args);
			current[j] = samp_fun(i,j,nparam,current,ran_gen,args);
			va_end(args);
		}
		memcpy(smp[i],current,nparam*sizeof(double));
		loadBar(i+1, niter, 100, 50);
	}
	fprintf(stderr, "\n");
	return samp;
}

double* mcmc_gibbs_samp_args(double* samp, size_t niter, size_t nparam, double* init, double (*samp_fun)(size_t iter,size_t pos,size_t nparam,double* param,gsl_rng* ran_gen,void** pass_args), gsl_rng* ran_gen, void** pass_args) {
	if (samp==NULL || init==NULL || samp_fun==NULL || ran_gen==NULL)
		return NULL;
	double (*smp)[nparam] = (double (*)[nparam])samp;
	double current[nparam];
	memcpy(current,init,nparam*sizeof(double));
	for (size_t i=0; i<niter; i++) {
		for (size_t j=0; j<nparam; j++) {
			current[j] = samp_fun(i,j,nparam,(double*)current,ran_gen,pass_args);
		}
		memcpy((double*)(smp[i]),(double*)current,nparam*sizeof(double));
	}
	return samp;
}

double* mcmc_metr_hast_args(double* samp, size_t niter, size_t nparam, double* init, double (*dist_fun)(size_t nparam,double* x,void**),
		double (*propdist_fun)(size_t nparam, double* x, double* x_old, void** pass_args), void (*samp_fun)(double* samp,size_t nparam,double* current,gsl_rng* ran_gen,void** pass_args), gsl_rng* ran_gen, void** pass_args){
	if (samp==NULL || init==NULL || dist_fun==NULL || propdist_fun==NULL || samp_fun==NULL || ran_gen==NULL)
		return NULL;
	double (*smp)[nparam] = (double (*)[nparam])samp;
	double current[nparam];
	double new[nparam];
	memcpy(current,init,nparam*sizeof(double));
	for (size_t i=0; i<niter; i++) {
		
		samp_fun(new,nparam,current,ran_gen,pass_args);
		double acc_ratio = dist_fun(nparam,new,pass_args)/dist_fun(nparam,current,pass_args)*
			propdist_fun(nparam,current,new,pass_args)/propdist_fun(nparam,new,current,pass_args);
		double acc_prob = (acc_ratio<1)? acc_ratio : 1;
		double accept = gsl_rng_uniform(ran_gen);
		if (accept < acc_prob) {
			memcpy(current,new,nparam*sizeof(double));
		}
		memcpy(smp[i],current,nparam*sizeof(double));
	}
	return samp;
}

double* mcmc_metr_hast_v(double* samp, size_t niter, size_t nparam, double* init, double (*dist_fun)(size_t nparam,double*x,va_list pass_args),
		double (*propdist_fun)(size_t nparam, double* x, double* x_old, va_list pass_args), void (*samp_fun)(double* samp,size_t nparam,double*,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen, va_list pass_args){
	if (samp==NULL || init==NULL || dist_fun==NULL || propdist_fun==NULL || samp_fun==NULL || ran_gen==NULL)
		return NULL;
	double (*smp)[nparam] = (double (*)[nparam])samp;
	double *current = malloc(nparam*sizeof(*current));
	double *new = malloc(nparam*sizeof(*new));
	va_list args0,args1,args2,args3,args4;
	memcpy(current,init,nparam*sizeof(double));
	for (size_t i=0; i<niter; i++) {
		va_copy(args0,pass_args);
		va_copy(args1,pass_args);
		va_copy(args2,pass_args);
		va_copy(args3,pass_args);
		va_copy(args4,pass_args);
		samp_fun(new,nparam,current,ran_gen,args0);
		double acc_ratio = dist_fun(nparam,new,args1)/dist_fun(nparam,current,args2)*
			propdist_fun(nparam,current,new,args3)/propdist_fun(nparam,new,current,args4);
		double acc_prob = (acc_ratio<1)? acc_ratio : 1;
		double accept = gsl_rng_uniform(ran_gen);
		if (accept < acc_prob) {
			memcpy(current,new,nparam*sizeof(double));
		}
		memcpy(smp[i],current,nparam*sizeof(double));
		va_end(args0);
		va_end(args1);
		va_end(args2);
		va_end(args3);
		va_end(args4);
	}
	free(current);
	free(new);
	return samp;
}

double mcmc_metr_hast_step(double current, double (*dist_fun)(double x,va_list pass_args), double (*propdist_fun)(double x, double x_old, va_list pass_args), 
		double (*samp_fun)(double current,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen, ...) {
	va_list pass_args;
	va_start(pass_args, ran_gen);
	double samp = mcmc_metr_hast_step_v(current,dist_fun,propdist_fun,samp_fun, ran_gen, pass_args);
	va_end(pass_args);
	return samp;
}
double mcmc_metr_hast_step_v(double current, double (*dist_fun)(double x,va_list pass_args), double (*propdist_fun)(double x, double x_old, va_list pass_args), 
		double (*samp_fun)(double current,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen, va_list pass_args) {
	
	va_list args0,args1,args2,args3,args4;
	va_copy(args0,pass_args);
	va_copy(args1,pass_args);
	va_copy(args2,pass_args);
	va_copy(args3,pass_args);
	va_copy(args4,pass_args);
	double new = samp_fun(current,ran_gen,args0);
	double dist_new = dist_fun(new,args1);
	double dist_cur = dist_fun(current,args2);
	double prop_new = propdist_fun(new,current,args4);
	double prop_cur = propdist_fun(current,new,args3);
//	double acc_ratio = dist_fun(new,args1)/dist_fun(current,args2)*
//			propdist_fun(current,new,args3)/propdist_fun(new,current,args4);
	double acc_ratio = dist_new/dist_cur*prop_cur/prop_new;
	double acc_prob = (acc_ratio<1)? acc_ratio : 1;
	double accept = gsl_rng_uniform(ran_gen);
	double samp;
	if (accept<acc_prob)
		samp = new;
	else
		samp = current;
	va_end(args0);
	va_end(args1);
	va_end(args2);
	va_end(args3);
	va_end(args4);
	
	return samp;
}

double mcmc_metr_hast_step_ratio(double current, double (*dist_ratio_fun)(double new, double current,va_list pass_args), double (*propdist_fun)(double x, double x_old, va_list pass_args), 
		double (*samp_fun)(double current,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen, ...) {
	va_list pass_args;
	va_start(pass_args, ran_gen);
	double samp = mcmc_metr_hast_step_ratio_v(current,dist_ratio_fun,propdist_fun,samp_fun, ran_gen, pass_args);
	va_end(pass_args);
	return samp;
}
double mcmc_metr_hast_step_ratio_v(double current, double (*dist_ratio_fun)(double new, double current,va_list pass_args), double (*propdist_fun)(double x, double x_old, va_list pass_args), 
		double (*samp_fun)(double current,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen, va_list pass_args) {
	
	va_list args0,args1,args2,args3;
	va_copy(args0,pass_args);
	va_copy(args1,pass_args);
	va_copy(args2,pass_args);
	va_copy(args3,pass_args);
	double new = samp_fun(current,ran_gen,args0);
	double dist_ratio = dist_ratio_fun(new, current,args1);
	double prop_new = propdist_fun(new,current,args2);
	double prop_cur = propdist_fun(current,new,args3);

	double acc_ratio = dist_ratio*prop_cur/prop_new;
	double acc_prob = (acc_ratio<1)? acc_ratio : 1;
	double accept = gsl_rng_uniform(ran_gen);
	double samp;
	if (accept<acc_prob)
		samp = new;
	else
		samp = current;
	va_end(args0);
	va_end(args1);
	va_end(args2);
	va_end(args3);
	
	return samp;
}

double mcmc_metr_hast_step_args(double current, double (*dist_fun)(double x,void** pass_args), double (*propdist_fun)(double x, double x_old, void** pass_args), 
		double (*samp_fun)(double current,gsl_rng* ran_gen,void** pass_args), gsl_rng* ran_gen, void** pass_args) {
	double new = samp_fun(current,ran_gen,pass_args);
	double prob_new = dist_fun(new,pass_args);
	double prob_cur = dist_fun(current,pass_args);
	double prob_p_cur = propdist_fun(current,new,pass_args);
	double prob_p_new = propdist_fun(new,current,pass_args);
	double acc_ratio = prob_new/prob_cur*prob_p_cur/prob_p_new;
	double acc_prob = (acc_ratio<1)? acc_ratio : 1;
	double accept = gsl_rng_uniform(ran_gen);
	double samp;
	if (accept<acc_prob)
		samp = new;
	else
		samp = current;
	return samp;
}
