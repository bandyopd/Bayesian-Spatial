# Bayesian-Spatial
## C and R implementation of Bayesian Spatial Additive Hazard Model

This model was developed by me during my Research work at University of Windsor.

The full text of the Thesis with details about the model, theoretical background, and results of Simulations and application of the model to real data, can be found [here] (http://scholar.uwindsor.ca/etd/4965/).

The computationally heavy task of fitting the model was written in *C* and compiled into a shared library that can be used by *R* for data analysis.

The source codes are available under the *./C-Codes* and *./R-Codes* directories.

The *C* code was designed for maximum reusability, hence generic methods like Gibbs Sampler and Metropolis-Hastings (Markov Chain Monte Carlo) were written in a form allowing to use them for any probability distribution function (see *./C-Codes/src/mcmc.h* and *./C-Codes/src/mcmc.c*).

The *C* implementation uses multiple-precision data types from open-source libraries *MPFR* and *GMP*. These libraries allow to operate with numbers much closer to zero than it is allowed by the standard floating-point types of *R* and *C* languages, which is required by the model.

Because of using *MPRF* and *GMP* libraries I had to implement some standard functions like CDF of Gamma distribution and Householder method of solving systems of linear equations, which are readily available only for standard floating-point data types and not for multiple-precision types.

Data used in the Thesis for analysis is not included in the repository because of the copyright restrictions.
