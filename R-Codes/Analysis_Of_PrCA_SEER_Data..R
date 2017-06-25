source('run_mod.R')

dat <- read.table('dat/dat.txt',header=T,sep='\t')
dat$region <- (dat$county-1)/2
# dat$time[dat$time==0] <- 0.00001
# dat$time[!dat$event] <- runif(length(dat$time),dat$time, dat$time+1)[!dat$event]
dat$enter <- 0
dat$time <- dat$time+1
tau <- max(dat$time)

model_matr <- model.matrix(~ age + race + stage + marriage, data=dat)[,-1] # exclude intercept
nobs <- nrow(model_matr)
ncov <- ncol(model_matr)
cv <- model_matr

ntimes_samp <- 72
ntimes <- ntimes_samp

covar <- array(cv, dim=c(nobs,ncov,ntimes))
covlim <- t(apply(cv, 2, function(x)range(c(x,0))))

times <- seq(0,tau,tau/ntimes)[-1]

niter <- 100
nburnin <- as.integer(niter/4)

regions <- read.table('dat/counties.txt', sep='\t')
regions <- regions[order(regions[,1]),]
nfrail <- nrow(regions)
matr <- as.matrix(apply(regions, 1, function(a) {r <- rep(0,nfrail); r[a[!is.na(a)]+1] <- 1; return(r)}))

dyn_lib_name <- 'car_model.so'
                  
run_model(do.samp=T,do.analysis=T,do.plot=T,do.plot.real=F,do.cpo=F)
