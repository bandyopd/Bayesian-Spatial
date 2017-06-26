# Bayesian Spatial Additive Hazard Model
## General Description

This model was developed by me and my supervisors during my Research work at University of Windsor.

The full text of the Thesis with details about the model, theoretical background, and results of Simulations and application of the model to real data, can be downloaded here: <a href="http://scholar.uwindsor.ca/etd/4965" target="_blank">Bayesian Spatial Additive Hazard Model</a>.

The computationally heavy task of fitting the model was written in *C* and compiled into a shared library that can be used by *R* for data analysis.

The source codes are available under the *./C-Codes* and *./R-Codes* directories.

The *C* code was designed for maximum reusability, hence generic methods like Gibbs Sampler and Metropolis-Hastings (Markov Chain Monte Carlo) were written in a form allowing to use them for any probability distribution function (see *./C-Codes/src/mcmc.h* and *./C-Codes/src/mcmc.c*).

The *C* implementation uses multiple-precision data types from open-source libraries *MPFR* and *GMP*. These libraries allow us to operate with numbers much closer to zero than it is allowed by the standard floating-point types of *R* and *C* languages, which is required by the model.

Because of using multiple-precision data types from the *MPRF* and *GMP* libraries, we had to implement some standard functions like CDF of Gamma distribution and Householder method for solving systems of linear equations, which are readily available only for standard floating-point data types and not for multiple-precision types. While implementing these functions, we used the codes from the open-source *GSL* library as a point of reference.

Data used in the Thesis for analysis is not included in the repository because of the copyright restrictions, however we include the dataset representing the adjustment structure of the counties in Louisiana. The first value in each row identifies the county, and the rest of the values list all counties adjucent to this county. This dataset has been created manually based on the map of the state Louisiana.

## Using the Codes

The simulations and data analysis for the research have been done in *R*, while the computational part of the model have been coded in *C*. This has been achieved through compiling *C* codes into a shared library that can be loaded from *R*.
 
The *C* code uses functions from three open-source libraries that have to be downloaded and linked to the *C* program:
- GNU Scientific Library (GSL) <a href="https://www.gnu.org/software/gsl" target="_blank">https://www.gnu.org/software/gsl</a>
- GNU Multiple Precision Arithmetic Library (GMP) <a href="https://gmplib.org/" target="_blank"> https://gmplib.org </a>
- GNU Multiple Precision Floating-Point Computations Library with Correct Rounding (MPFR) <a href="www.mpfr.org" target="_blank"> www.mpfr.org </a>

Once the above listed libraries have been downloaded, you should compile *./C-Codes/car_model.c* into a shared library (you will need to link the open source libraries and all the codes from *./C-Codes/src* to the *car_model.c* during compilation). You might also need to specify the location of the header file *R.h* on your system (usually within the *R* home directory under *include*, and might need to be installed if it didn't come with your distribution of *R*).

Once the shared library is created, the C functions can be accessed from R (see *./R-Codes/run_mod.R*). The main function in *run_mod.R* is *run_model* that has various parameters specifying the actions to be performed (the meaning of the function parameters can be found in the commentaries within the code). 

The function expects several input variables to be defined:
- *dyn_lib_name* - the name of (and/or path to) the dynamic (shared) library containing the compiled *C* functions
- *tau* - the length of the study period
- *nobs* - number of observations in the data
- *ncov* - number of covariates
- *nfrail* - number of frailty term
- *times* - a vector oftime points at which the covariate values are specified
- *covar* - an array of the values of the covariates at the above defined time points (the covariate functions are assumed to be constant between the points)
- *covlim* - an array of ranges of admissible values for covariats
- *ntimes_samp* - the number of time intervals considered in the model (for MCMC sampling); the time intervals are assumed to be equidistant; note that for the purpose of fitting a model, the covariate functions defined above will be approximated by piecewise constant functions on the time intervals defined by the constant *ntimes_samp* which in general are different from the intervals defined by the vector *times*
- *matr* - the adjustment matrix of spacial regions

The program *./R-Codes/Analysis_Of_PrCA_SEER_Data.R* demonstrates how the model can be applied to real data. The data used in this program is the SEER Prostate Cancer data described in the Thesis. As was already mentioned before, the data itself is not included in this repository due to copyright restrictions, however the counties adjustment data used in the program can be found in *./data/counties.txt*

After the model is run, the sample produced by the MCMC algorithm will be saved in the following variables:
- *samp_bl* - sample for baseline hazard: rows represent different sample iterations, columns represent different time points
- *samp_rf* - a list of samples for regression functions (each element of the list is a matrix of values, representing a sample for one variable: rows represent different sample iterations, columns represent different time points)
- *samp_fr* - sample for frailties (each element of the list is a matrix of values, representing a sample for one frailty: rows represent different sample iterations, columns represent different time points)

Estimated quantiles are stored in the following variables:
- *quant_bl* - quantiles for baseline hazard: rows represent 2.5%, 50% and 97.5% quantiles, columns represent different time points 
- *quant_rf_list* - a list of quantiles for regressions functions (each element of the list is a matrix of values: rows represent 2.5%, 50% and 97.5% quantiles, columns represent different time points)
- *quant_fr_list* - a list of quantiles for frailties (each element of the list is a matrix of values: rows represent 2.5%, 50% and 97.5% quantiles, columns represent different time points)
