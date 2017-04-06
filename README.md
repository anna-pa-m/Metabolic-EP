# Metabolic-EP
Expectation Propagation algorithm for metabolic networks
=======
# Authors: 
Alfredo Braunstein, Anna Paola Muntoni and Andrea Pagnani
# Date:
6th April 2017


# Description
This is an implementation of the Expectation Propagation algorithm for studying the space of solution of constrained metabolic fluxes. 
The main outputs of the m-function are the means and the variances of truncated Gaussian distributions that approximate the marginal probability density of observing a flux, given a stoichiometric matrix and a measure of the intakes/uptakes. This is part of the work:

"An analytic approximation of the feasible space of metabolic networks" - A. Braunstein, A. Muntoni, A. Pagnani - doi:10.1038/ncomms14915

Input
- S: stoichiometric matrix of "Nm metabolites" x "Nr reactions"
- b: vector of Nm intakes/uptakes
- nuinf, nusup: lower and upper bounds for each metabolic flux
- Beta: inverse variance of the noise, if any. Otherwise a "large" number (ex. 1e9)
- damp: damping coefficient (from 0 to 1) applied to the update of means "a" and variances "d" of approximating Gaussians
        Ex. "new a" = damp * "new a" + (1 - damp) * "old a"
- max_iter: maximum number of iterations (ex. 1e3)
- minvar, maxvar: lower and upper bounds for the variances "d" of the approximation. (ex. 1e-50, 1e50)
- precision:  precision required to stop the algorithm (ex. 1e-9)

 
Input (optional) to fix an experimental profile
- av_exp: mean of the experimental profile
- var_exp: variance of the experimental profile
- exp_i: index of the measured flux
If no experimental evidence is available, set av_exp = 0, var_exp = 0 and exp_i = 0.	


Output
- mu: vector parametrizing the mean of the posterior distribution
- s: vector parametrizing the variance of the posterior distribution
- a: vector containing the means of the approximated priors
- d: vector containing the variances of the approximated priors
- av: averages of the truncated Gaussians of the approximation
- va: variances of the truncated Gaussians of the approximation
- t: running time



Julia Version

Input
- S: stoichiometric matrix of "Nm metabolites" x "Nr reactions" (either sparse or dense)
- b: vector of Nm intakes/uptakes
- nuinf, nusup: lower and upper bounds for each metabolic flux

Optional Arguments
- beta::Real: inverse variance of the noise [default = 1e-7]
- verbose::Bool verbosity of output [default = true]
- damp::Real damping coefficient (from 0 to 1) applied to the update of means "a" and variances "d" of approximating Gaussians Ex. "new a" = damp * "new a" + (1 - damp) * "old a"; [default = 0.9]
- espsilonconv:  precision required to stop the algorithm [default 1e-6]
- maxiter::Int: maximum number of iterations [default = 2000]
- minvar::Real, maxvar::Real: lower and upper bounds for the variances "d" of the approximation. [default: 1e-50, 1e50]


Output
(μ,s, av, var) where:

- μ: vector parametrizing the mean of the posterior distribution
- s: vector parametrizing the variance of the posterior distribution
- av: averages of the truncated Gaussians of the approximation
- va: variances of the truncated Gaussians of the approximation


Test model: iJR904 model for Escherichia Coli. https://doi.org/10.1093/nar/gkv1049 
