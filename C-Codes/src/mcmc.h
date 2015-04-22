#ifndef _MCMC_H_INCLUDED_
#define _MCMC_H_INCLUDED_

#include<stdlib.h>
#include<string.h>
#include<stdarg.h>

#include <gsl/gsl_rng.h>

#include <test_in_function.h>

double* mcmc_gibbs_samp(double* samp, size_t niter, size_t nparam, double* init,  double (*samp_fun)(size_t iter,size_t pos,size_t nparam,double* param,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen,...);
double* mcmc_gibbs_samp_v(double* samp, size_t niter, size_t nparam, double* init,  double (*samp_fun)(size_t iter,size_t pos,size_t nparam,double* param,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen, va_list pass_args);
double* mcmc_gibbs_samp_args(double* samp, size_t niter, size_t nparam, double* init, double (*samp_fun)(size_t iter,size_t pos,size_t nparam,double* param,gsl_rng* ran_gen,void** pass_args), gsl_rng* ran_gen, void** pass_args);

double* mcmc_metr_hast_args(double* samp, size_t niter, size_t nparam, double* init, double (*dist_fun)(size_t nparam,double* x,void** pass_args),
		double (*propdist_fun)(size_t nparam,double* x, double* x_old,void**), void (*samp_fun)(double* samp, size_t nparam, double* current,gsl_rng* ran_gen,void** pass_args), gsl_rng* ran_gen, void** pass_args);

double* mcmc_metr_hast_v(double* samp, size_t niter, size_t nparam, double* init, double (*dist_fun)(size_t nparam,double* x,va_list pass_args),
		double (*propdist_fun)(size_t nparam,double* x, double* x_old,va_list pass_args), void (*samp_fun)(double* samp, size_t nparam, double* current,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen, va_list pass_args);

double mcmc_metr_hast_step_ratio(double current, double (*dist_ratio_fun)(double new, double current,va_list pass_args), double (*propdist_fun)(double x, double x_old, va_list pass_args), 
		double (*samp_fun)(double current,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen, ...);
double mcmc_metr_hast_step_ratio_v(double current, double (*dist_ratio_fun)(double new, double current,va_list pass_args), double (*propdist_fun)(double x, double x_old, va_list pass_args), 
		double (*samp_fun)(double current,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen, va_list pass_args);

double mcmc_metr_hast_step(double current, double (*dist_fun)(double x, va_list pass_args), double (*propdist_fun)(double x, double x_old, va_list pass_args), 
		double (*samp_fun)(double current, gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen, ...);

double mcmc_metr_hast_step_v(double current, double (*dist_fun)(double x,va_list pass_args), double (*propdist_fun)(double x, double x_old, va_list pass_args), 
		double (*samp_fun)(double current,gsl_rng* ran_gen,va_list pass_args), gsl_rng* ran_gen, va_list pass_args);

double mcmc_metr_hast_step_args(double current, double (*dist_fun)(double x,void** pass_args), double (*propdist_fun)(double x, double x_old, void** pass_args), 
		double (*samp_fun)(double current,gsl_rng* ran_gen,void** pass_args), gsl_rng* ran_gen, void** pass_args);
		
#endif
