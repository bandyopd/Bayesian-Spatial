# Calculate stepwise approximations of covariate functions 
gen_data_stepwise_approx = function(x, fun) {
	n = length(x);
	xprev = c(0,x[-n]);
	val = fun((x+xprev)/2);
	return (val);
}

# Calculate covariance matrix of CAR Frailties with first frailty equal to zero
frailty_car_cov_matr_firstzero = function(adj_matr,var_param) {
	n = dim(adj_matr)[1];
	prec_matr = adj_matr[-1,-1]
	nadj = apply(adj_matr,1,sum)[-1]
	prec_matr = diag(nadj)-prec_matr
	cov_matr = var_param*chol2inv(chol(prec_matr))
	return(cov_matr)
}

# Calculate hazard based on input data (a wrapper for a C function)
compute_hazard = function(baseline, covariates, regfun, frail, frail_tbl) {
	nobs = dim(covariates)[1];
	ncov = dim(covariates)[2];
	nfrail = dim(frail)[1];
	ntimes = length(baseline);
	res = .C('compute_hazard',haz=double(nobs*ntimes), as.integer(nobs), as.integer(ncov),as.integer(nfrail), as.integer(ntimes), as.double(baseline), as.double(aperm(covariates)), as.double(aperm(regfun)),as.double(aperm(frail)), as.integer(frail_tbl))
	haz = matrix(res[['haz']],nobs, ntimes, byrow=TRUE)
	return(haz)
}

# Fit the model with CAR frailties (a wrapper for a C function)
fit_car_model = function(niter, dat, times, covar, covlim, adj_matr, r0=sum(dat$event)/dim(dat)[1]/max(times), c0=0.001, alpha=0.01, beta=0.001) {
	nobs = dim(covar)[1]
	ncov = dim(covar)[2]
	ntimes = dim(covar)[3]
	nfrail = dim(adj_matr)[1]
	npar = ntimes*(ncov+nfrail+2)
	res = .C('fit_car_model', samp=double(npar*niter), 
			as.integer(niter), as.integer(nobs), as.integer(ncov), as.integer(nfrail), as.integer(ntimes), 
			as.double(dat$enter), as.double(dat$time), as.integer(dat$event), as.integer(dat$region), 
			as.double(times), as.double(aperm(covar)), as.double(aperm(covlim)), as.integer(adj_matr),
			as.double(r0), as.double(c0), as.double(alpha), as.double(beta))
	samp = matrix(res[['samp']],niter,npar,byrow=TRUE)
	samp[,npar] = 1/samp[,npar]
	return(samp)
}

# Create a PDF file containing a Regression Function plot. If "real function" is available (as in case of simulations), it is plotted as well
plot_to_pdf <- function(prefix='',filename, real_fun, quant, label, ylim, ...){
	pdf(paste0('img/',prefix,filename,'_',as.integer(nobs),'-',as.integer(ntimes_samp),'-',as.integer(niter),'.pdf'))
	if (missing(ylim)){
		ylim <- range(quant,na.rm=TRUE, finite=TRUE)
		if (!missing(real_fun)) {
			ylim <- range(c(real_fun(times),ylim),na.rm=TRUE, finite=TRUE)
		}
	}
	if (!missing(real_fun)) {
		plot(real_fun,0,tau,xlab='t',ylab=label,ylim=ylim,main=label,lwd=2, ...)
	}
	for (q in 1:dim(quant)[2]){
		if (q==dim(quant)[2]%/%2+1){
			lty <- 1
		} else {
			lty <- 2
		}
		if (q==1 && missing(real_fun)) {
			plot(c(0,times_samp),c(quant[1,q],quant[,q]),type='S',xlab='t', ylab=label, ylim=ylim,main=label, lty=lty, ...)
		} else {
			lines(c(0,times_samp),c(quant[1,q],quant[,q]),type='S',lty=lty)
		}
	}
	dev.off()
}

# Create multiple PDF files with Regression Function plots for different variables
plot_to_pdf_many <- function(n, prefix='', filename, real_fun_many, quant_list, label, ylim, common_scale=TRUE, ...) {
	if (missing(ylim)){
		ylim <- range(unlist(quant_list),na.rm=TRUE, finite=TRUE)
		if (!missing(real_fun_many)) {
			ylim <- range(c(sapply(1:n,function(i)real_fun_many(times,i)),ylim),na.rm=TRUE, finite=TRUE)
		}
	}
	for (i in 1:n){
		if (!missing(real_fun_many)) {
			if (common_scale) {
				plot_to_pdf(prefix=prefix,filename=paste0(filename,i), real_fun=function(x)real_fun_many(x,i), quant=quant_list[[i]], label=paste(label,i), ylim=ylim, ...)
			} else {
				plot_to_pdf(prefix=prefix,filename=paste0(filename,i,'_sc'), real_fun=function(x)real_fun_many(x,i), quant=quant_list[[i]], label=paste(label,i), ...)
			}
		} else {
			if (common_scale) {
				plot_to_pdf(filename=paste0(filename,i), quant=quant_list[[i]], label=paste(label,i), ylim=ylim,...)
			} else {
				plot_to_pdf(filename=paste0(filename,i,'_sc'), quant=quant_list[[i]], label=paste(label,i),...)
			}
		}
	}
}

# Calculate the 2.5%, 50% and 97.5% quantiles for the MCMC sample (the input parameter "samp" is the result of MCMC algorithm)				    
quant_for_samp <- function(samp) {
	quant_val <- c(0.025,0.5,0.975)
	if (is.array(samp)) {
		quant <- t(apply(samp, 2, quantile, quant_val, na.rm=TRUE))
	} else {
		quant <- matrix(quantile(samp, quant_val, na.rm=TRUE),1,length(quant_val))
	}
	return(quant) 
}

# Run and analyze the model, and produce the plots of the estimated Regression Functions					    
run_model <- function(do.samp=F, do.analysis=F, do.plot=F, do.plot.real=F, common_scale=F, do.cpo=F, do.dic=F,do.loglik=do.cpo||do.dic, do.survival=F) {
	dyn.load(dyn_lib_name)
	# MCMC Sampling
	# This step will only be performed if do.samp=T
	# If do.samp=F then it's assumed that MCMC sampling have already been performed 
	# and the function is called only to analyze the results
	if (do.samp) {
		cat('Sampling...\n')
		times_samp <<- seq(0,tau,tau/ntimes_samp)[-1] # exclude 0
		covar_fun <<- function(t){
			n_int <<- sapply(t, function(t) min(which(times>=t)))
			return(array(covar[,,n_int],c(dim(covar)[1:2],length(t))))
		}
		covar_samp <<- gen_data_stepwise_approx(times_samp, covar_fun)
		samp <<- fit_car_model(niter, dat, times_samp, covar_samp, covlim, matr)
	}
					 
	# Analysis of MCMC sample (calculating quantile estimates based on the sample)
	# This step will only be performed if do.analysis=T
	# if do.analysis=F, this step will be assumed to have already been performed earlier and will be skipped
	if (do.analysis) {
		cat('Analysing\n')
		if (nburnin==0) {
			samp1 <<- samp
		} else {
			samp1 <<- samp[-(1:nburnin),1:(ntimes_samp*(1+ncov+nfrail)),drop=F]
		}
		samp_bl <<- samp1[,1:ntimes_samp,drop=F]
		samp_rf <<- vector("list", ncov)
		for (i in 1:ncov){
			samp_rf[[i]] <<- samp1[,(1:ntimes_samp-1)*ncov+i+ntimes_samp,drop=F]
		}
		samp_fr <<- vector("list", nfrail)
		for (i in 1:nfrail){
			samp_fr[[i]] <<- samp1[,(1:ntimes_samp-1)*nfrail+i+ntimes_samp*(1+ncov),drop=F]
		}
		quant_bl <<- quant_for_samp(samp_bl)
		quant_rf_list <<- vector("list", ncov)
		for (i in 1:ncov){
			quant_rf_list[[i]] <<- quant_for_samp(samp_rf[[i]])
		}
		quant_fr_list <<- vector("list", nfrail)
		for (i in 1:nfrail){
			quant_fr_list[[i]] <<- quant_for_samp(samp_fr[[i]])
		}
		
	}
	# Plotting the Regression Function estimates
	# This step will be performed only if do.plot=T
	# The plots will be saved to PDF files				 
	# if do.plot=F then PDF files with Regression Function estimates will not be produced
	if (do.plot) {
		cat('Producing plots\n')
		if (do.plot.real) {
			plot_to_pdf(filename='baseline', real_fun=baseline_fun, quant=quant_bl, label='Baseline')
			reg_fun_many <<- function(x,i) {
				if (i==1) {
					return(reg_fun_1(x))
				} else if (i==2) {
					return(reg_fun_2(x))
				} else {
					return(0*x)
				}
			}
			plot_to_pdf_many(n=ncov, filename='regfun', real_fun_many = reg_fun_many, quant_list = quant_rf_list, label='Regression function',common_scale=common_scale)
			
			frail_fun_many <<- function(x,i) 0*x+fr[i]
			plot_to_pdf_many(n=nfrail, filename='frail', real_fun_many = frail_fun_many, quant_list = quant_fr_list, label='Frailty',common_scale=common_scale)
		} else {
			plot_to_pdf(filename='baseline', quant=quant_bl, label='Baseline')
			plot_to_pdf_many(n=ncov, filename='regfun', quant_list = quant_rf_list, label='Regression function',common_scale=common_scale)
			plot_to_pdf_many(n=nfrail, filename='frail', quant_list = quant_fr_list, label='Frailty',common_scale=common_scale)
		}
		
		
	}
	# Calculating Log-Likelihood (required for calculating Bayesian Model Comparison Criteria CPO and DIC				 
	# This step will be performed only if do.loglik=T
	# Use do.loglik=F if log-likelihood has already been calculated or calculation of CPO and DIC is not required
	# NOTE: This step is computationally heavy, so it's recommended you set do.loglik=F if you don't need it calculated				 
	if (do.loglik) {
		cat('Computing likelihood\n')
		ev_int <<- sapply(dat$time[dat$event==1], function(t) min(which(times_samp>=t)))
		event_matr <<- matrix(0,nobs,ntimes_samp)
		event_matr[which(dat$event==1),] <<- t(sapply(ev_int, function(int) times_samp == times_samp[int]))
		atrisk_matr <<- sapply(times_samp, function(t) dat$enter<=t & dat$time>=t)
		times_samp_prev <<- c(0,times_samp)[-length(times_samp)-1]
		enter_matr <<- matrix(dat$enter, nobs, ntimes_samp, byrow=F)
		prev_matr <<- matrix(times_samp_prev, nobs, ntimes_samp, byrow=T)
		times_matr <<- matrix(times, nobs, ntimes_samp, byrow=T)
		evtime_matr <<- matrix(dat$time, nobs, ntimes_samp, byrow=F)
		atrisk_times <<- pmin(times_matr,evtime_matr)-max(enter_matr,prev_matr)
		loglik_fun <<- function(i) {
			bl <- samp_bl[i,]
			rf <- t(sapply(samp_rf, function(samp)samp[i,]))
			fr <- t(sapply(samp_fr, function(samp)samp[i,]))
			hz <- compute_hazard(bl, covar_samp, rf, fr, dat$region)
			return(loglik_marginal(hz,event_matr,atrisk_times))
		}			
		loglik <<- sapply(1:nrow(samp1),loglik_fun)
		bl_mean <<- colMeans(samp_bl)
		rf_mean <<- t(sapply(samp_rf, colMeans))
		fr_mean <<- t(sapply(samp_fr, colMeans))
		haz_of_mean <<- compute_hazard(bl_mean, covar_samp, rf_mean, fr_mean, dat$region)
		loglik_of_mean <<- loglik_marginal(haz_of_mean, event_matr, atrisk_times)
	}
	# Calculate CPO criterion value			       
	if (do.cpo) {
		cat('Computing CPO\n')
		cpo_terms <<- 1/exp(loglik)
		cpo <<- 1/apply(cpo_terms,1,mean)
		lcpo <<- sum(log(cpo))	
		write(lcpo, file = paste0('out/lcpo_',as.integer(ntimes_samp),'_',as.integer(niter),'.txt'))
	}
	# Calculate DIC cretetion value			       
	if (do.dic) {
		cat('Computing DIC\n')
		dev <<- -2*apply(loglik, 1, mean)
		dev_sum <<- sum(dev)
		# computing deviance at mean parameters
		npar_eff <<- dev_sum + 2*sum(loglik_of_mean)
		dic <<- dev_sum + npar_eff
		write(c(dic,npar_eff), file = paste0('out/dic_',as.integer(ntimes_samp),'_',as.integer(niter),'.txt'))
	}
	# Plot a "median" Survival Function, i.e. Survival Function for an individual with median values of all covariates			       
	if (do.survival) {
		cat('Computing survival\n')
		bl_median <<- quant_bl[,'50%']
		rf_median <<- matrix(t(sapply(quant_rf_list,function(quant)quant[,'50%'])),ncov,ntimes)
		fr_median <<- matrix(t(sapply(quant_fr_list,function(quant)quant[,'50%'])),nfrail,ntimes)
		median_of_fr <<- colMedians(fr_median)
		fr_to_use <<- median_of_fr
		cv_median <<- apply(covar_samp,3,colMedians)
		times_samp_prev <<- c(0,times_samp)[-length(times_samp)-1]
		library(survival)
		fit <<- survfit(formula = Surv(time,event)~stage,data=dat)
		pdf(paste0('img/median_survival','_',as.integer(ntimes_samp),'-',as.integer(niter),'.pdf'))
		plot(fit,conf.int=F,mark.time=F,xlab='t',ylab='Survival',main='Survival',lty=c(1,2))
		for (i in 0:1) {
			cv_to_use <<- cv_median
			cv_to_use[3,] <<- i # stage
			rfcv_median <<- colSums(rf_median*cv_to_use)
			hz_median <<- fr_to_use + rfcv_median + bl_median
			hzdt_median <<- hz_median * (times_samp-times_samp_prev)
			surv_median <<- exp(-cumsum(hzdt_median))
			nsteps <<- max(as.integer(1000/ntimes_samp),1)
			tm <<- seq(0,tau, tau/(nsteps*ntimes_samp))[-1]
			hzt <<- rep(hz_median,rep(nsteps,ntimes_samp))
			tp <<- c(0,tm)[-(length(tm)+1)]
			hzt_dt <<- hzt*(tm-tp)
			surv_to_plot <<- exp(-cumsum(hzt_dt))
			lines(c(0,tm),c(1,surv_to_plot),type='l',lty=(i+1))
		}
		dev.off()
	}
	cat('Done\n')
	dyn.unload(dyn_lib_name)
}

# Calculate medians of all columns in a matrix					      
colMedians <- function(matr) {
	apply(matr,2,median)
}

# Calculate marginal log-likelihood 
loglik_marginal <- function(haz,ev,ar_time) {
	loghaz <- apply(log(haz)*ev,1,sum)
	loglik <- loghaz-apply(haz*ar_time,1,sum)
}
