# Metabolic-EP
Expectation Propagation algorithm for metabolic networks
=======
# Authors:
Alfredo Braunstein, Anna Paola Muntoni and Andrea Pagnani


# Description

This is an implementation of the Expectation Propagation algorithm for
studying the space of solution of constrained metabolic fluxes.  The
main outputs of the m-function are the means and the variances of
truncated Gaussian distributions that approximate the marginal
probability density of observing a flux, given a stoichiometric matrix
and a measure of the intakes/uptakes. This is part of the work:

["An analytic approximation of the feasible space of metabolic
networks"](http://www.nature.com/articles/ncomms14915) -
A. Braunstein, A. Muntoni, A. Pagnani - Nature Communications 8,
Article number: 14915 (2017) - doi:10.1038/ncomms14915

There are two implementations: one in matlab (under folder `matlab`), the second in [Julia](http://julialang.org).


Matlab Version
==============

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
=============

Installing the package

``julia> Pkg.clone("https://github.com/anna-pa-m/Metabolic-EP/","MetabolicEP.jl")``.

Otherwise, if you do not want to use the package manager, from a local copy of  the directory ``src`` in this repository, you can
``julia> include("dirtosource/src/MetabolicEP.jl"); using MetabolicEP``

It works with version 0.5, and 0.4 (with some warnings).

Typical usage is

``julia> res=metabolicEP(S,b,numin,numax)``

The output in res is of type ``EPout`: there are several fields:
-   ``μ::Vector``: A parameter linked to the mean of the posterior probability
-   ``σ::Vector``: A parameter linked to the std  of the posterior probability
-   ``av::Vector``: The mean posterior probability
-   ``va::Vector``: The variance of the posterior probability
-   ``sol::EPFields``: The internal field status. From this value we can
restart the sampling from a specific state.
-   ``status::Symbol``: either ``:converged`` or ``:unconverged``.


Input (required)
----
- `S`: MxN matrix (either sparse or dense) please note that if you input a dense version, the algorithm is slighlty more efficient. Dense matrices can be create from sparse ones with ``full(S)``.
- `b`: a vector of M intakes/uptakes
- `nuinf`: a vector of lengh N of upper bounds.
- `nusup`: a vector of lengh N of lower bounds.


Input (optional argument).
----
- `beta` (inverse temperature::``Real``): default 10^7
- `verbose` (``true`` or ``false``): default ``true``
- `damp` (∈ (0,1) newfield = damp * oldfield + (1-damp)* newfield): default 0.9  
- `epsconv` (convergence criterion): default 1e-6
- `maxiter` (maximum number of iterations): default 2000
- `maxvar`  (threshold on maximum variance): default 1e50
- `minvar`  (threshold on minimum variance): default 1e-50
- `solution` (start from solution. Is of type ``EPout``): default: ``nothing``
- `expval` (fix to posterior probability of mean and/or variance to
values): default ``nothing``. expval can be either at
``Tuple{Float64,Float64,Int}`` or a
``Vector{Tuple{Float64,Float64,Int}}``. Values can be fixed as
``expval=(0.2,0.4,4)`` meaning that for flux index 4 the mean is set to 0.2
and the variance to 0.4. Fixing more values ``expval=[(0.2, 0.3, 4),
(0.4, nothing, 5)]``: in this case, we fix the posterior of flux 4 to
0.2 (mean) and 0.3 (variance), while for flux 5 we fix the mean to 0.4
and we keep the variance free.

[COBRA](https://github.com/opencobra/COBRA.jl) compatibility
---

We developed a COBRA compatibility so that now models can be loaded
with the `COBRA.loadModel()`
[utility](https://opencobra.github.io/COBRA.jl/stable/functions.html#loadModel)
metabolicEP can be also run passing a `LPproblem`
[type](https://opencobra.github.io/COBRA.jl/stable/functions.html#LPproblem)
as returned by `loadModel`.

COBRA has not yet been updated to `julia v1.0`. For this reason the compatibility has been temporarily removed.

Reading matlab metabolic reconstruction (.mat files)
---

There is a small convenience reader for metabolic reconstructions in
matlab format (.mat).  It can be invoked as:

``julia> met=ReadMatrix("nomefile.mat")``

The output `met` is of type ``MetNet`` whose fields are:
- ``N::Int`` number of fluxes
- ``M::Int`` number of metabolites
- ``S::SparseMatrixCSC{Float64,Int}`` Stoichiometric matrix M x N sparse
- ``b::Array{Float64,1}``  right hand side of equation  S ν = b (vector of size M)
- ``c::Array{Float64,1}`` reaction index of biomass (vector of size N)
- ``lb::Array{Float64,1}`` fluxes lower bound N elements vector
- ``ub::Array{Float64,1}``  fluxes upper bound N elements vector
- ``genes::Array{String,1}``  gene names N elements vector
- ``rxnGeneMat::SparseMatrixCSC{Float64,Int}``  
- ``grRules::Array{String,1}``  gene-reaction rule N elements vector of strings (and / or allowed)
- ``mets::Array{String,1}``  metabolites short-name M elements
- ``rxns::Array{String,1}``  reactions short-name N elements
- ``metNames::Array{String,1}``  metabolites long-names M elements
- ``metFormulas::Array{String,1}`` metabolites formula M elements
- ``rxnNames::Array{String,1}``  reactions long-names N elements
- ``rev::Array{Bool,1}``  reversibility of reactions N elements
-  ``subSystems::Array{String,1}``  cellular component of fluxes N elements


Test model (in folder data): iJR904 model for Escherichia
Coli. https://doi.org/10.1093/nar/gkv1049
