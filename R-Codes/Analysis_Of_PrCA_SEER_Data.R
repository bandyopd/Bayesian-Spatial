source('run_mod.R')

# reaing in data table
dat <- read.table('dat.txt',header=T,sep='\t')
dat$region <- (dat$county-1)/2
# dat$time[dat$time==0] <- 0.00001
# dat$time[!dat$event] <- runif(length(dat$time),dat$time, dat$time+1)[!dat$event]
dat$enter <- 0
dat$time <- dat$time+1
tau <- max(dat$time)

# transforming data into a matrix format with dummy variables
model_matr <- model.matrix(~ age + race + stage + marriage, data=dat)[,-1] # exclude intercept
nobs <- nrow(model_matr)
ncov <- ncol(model_matr)
cv <- model_matr

# number of time intervals considered in the model
ntimes_samp <- 72
ntimes <- ntimes_samp

# creating a 3-dimentional array of covariate function values 
# (a value for each variable, for each individual at each time point)
# since the covariates in data are constant (not time-varying), 
# the same value of the covariate is repeated for each time point
covar <- array(cv, dim=c(nobs,ncov,ntimes))
covlim <- t(apply(cv, 2, function(x)range(c(x,0))))

# equidistant time points                  
times <- seq(0,tau,tau/ntimes)[-1]

# number of iterations and burn-in iteration for the MCMC algorithm                  
niter <- 100
nburnin <- as.integer(niter/4)

# reading in counties adjustment data and transforming it into an adjustment matrix format                  
regions <- read.table('counties.txt', sep='\t')
regions <- regions[order(regions[,1]),]
nfrail <- nrow(regions)
matr <- as.matrix(apply(regions, 1, function(a) {r <- rep(0,nfrail); r[a[!is.na(a)]+1] <- 1; return(r)}))

# the name of the dynamic library with the compiled C functions                  
dyn_lib_name <- 'car_model.so'

# running the model                  
run_model(do.samp=T,do.analysis=T,do.plot=T,do.cpo=T,do.dic=T)
