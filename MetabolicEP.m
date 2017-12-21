function [mu, s, a,d, av, va, Cov, t] = MetabolicEP(S, b, nuinf, nusup, Beta, damp, max_iter, minvar, maxvar, precision, av_exp, var_exp, exp_i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%																			
%%% Expectation Propagation algorithm for metabolic networks								
%%% 																				
%%% Authors: Alfredo Braunstein, Andrea Pagnani and Anna Paola Muntoni
%%% Date: April 2017
%%%
%%% Input
%%% - S: stoichiometric matrix of "Nm metabolites" x "Nr reactions".
%%% - b: vector of Nm intakes/uptakes
%%% - nuinf, nusup: lower and upper bounds for each metabolic flux
%%% - Beta: inverse variance of the noise, if any. Otherwise a "large" number (ex. 1e9)
%%% - damp: damping coefficient (from 0 to 1) applied to the update of means "a" and variances "d" of approximating Gaussians
%%%	    Ex. "new a" = damp * "new a" + (1 - damp) * "old a"
%%% - max_iter: maximum number of iterations
%%% - minvar, maxvar: lower and upper bounds for the variances "d" of the approximation. (ex. 1e-50, 1e50)
%%% - precision:  precision required to stop the algorithm (ex. 1e-9)
%%%
%%% Input (optional) to fix an experimental profile
%%% - av_exp: mean of the experimental profile
%%% - var_exp: variance of the experimental profile
%%% - exp_i: index of the measured flux
%%% If no experimental evidence is available, set av_exp = 0, var_exp = 0 and exp_i = 0.	
%%%
%%%
%%% Output
%%% - mu: vector parametrizing the mean of the posterior distribution
%%% - s: vector parametrizing the variance of the posterior distribution
%%% - a: vector containing the means of the approximated priors
%%% - d: vector containing the variances of the approximated priors
%%% - av: averages of the truncated Gaussians of the approximation
%%% - va: variances of the truncated Gaussians of the approximation
%%% - Cov: convariance matrix of the fluxes
%%% - t: running time
%%%
%%% The marginal probability density of the n-th flux is a truncated Gaussian N(mu(n),s(n)) in the interval [nuinf, nusup]. It has average av(n) and variance va(n)
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nr = size(S,2);
Nm = size(S,1);
a = zeros(Nr,1,'double');
d = ones(Nr,1,'double');
av = zeros(Nr,1,'double');
va = ones(Nr,1,'double');
mu = zeros(Nr,1,'double');
s = ones(Nr,1,'double');
factor = max(max(abs(nuinf)), max(abs(nusup)));
nusup = nusup / factor;
nuinf = nuinf / factor;
b = b / factor;
av_exp = av_exp /factor;
var_exp = var_exp / factor^2;
D = sparse(1:Nr, 1:Nr, 1./d);
KK = Beta * (S' * S);
KB = Beta * S' * b;
err = 100;
iter = 0;
n = size(nuinf,1);


tic
while (err > precision && iter < max_iter)
    iter = iter + 1;

    % fast computation of the means and variances of the truncated Gaussians 
    I1 = inv(KK + D);
    v = (I1) * (KB + D*a);
    I = diag(I1);
    I = min(I,d);
    s1 = min(maxvar, max(minvar, 1./I - 1./d));
    %s = max(minvar,1./s1);
    s = 1./s1;
    mu = (v-(a.*I)./d)./(1-I./d);
    mu(I == d) = 0.5*(nuinf(I == d) + nusup(I == d));
  
    % compute means and variances of the tilted distributions
    s05 = s.^0.5;
    x0 = (nuinf-mu)./s05; 
    x1 = (nusup-mu)./s05;     
    [z,eps] = compute_mom_vec(x0, x1);    
    oldvar = va;
    oldav = av;
    av = mu + z .* s05;
    va = max(0, s .* (1 + eps));
    % insert experimental evidences
    if(exp_i)
        av(exp_i) = av_exp;
    	va(exp_i) = var_exp;
    end
    % moments matching (update "a" and "d")
    err = max(max(abs(av-oldav)), max(abs(va-oldvar)));
    new_d = 1./(1./va - 1./s);
    new_d = min(maxvar, max(minvar, new_d));
    new_a = av + (av - mu) .* new_d.* s1;    
    a = damp * a + (1-damp) * new_a;
    d = damp * d + (1-damp) * new_d;
    %D = diag(1./d);
    D = sparse(1:Nr, 1:Nr, 1./d);

    if(exp_i)
        err = max(err, abs(a(exp_i) -new_a(exp_i)) + abs(d(exp_i) - new_d(exp_i)));
    end
    fprintf('it:%i err:%e beta:%e\n', iter, err, Beta);
end

toc
t=toc;
mu = mu * factor;
s = s * factor^2;
a = a*factor;
d = d*factor^2;
av = av*factor;
va = va*factor^2;
D = sparse(1:Nr, 1:Nr, 1./d);
Cov = inv(KK + D);
fprintf('Time %f\n', t);

end

function y = Phi(x)
    y = 0.5 * (1 + erf(x/sqrt(2)));
end

function y = phi(x)
    y = 1/sqrt(2*pi) * exp(-x.^2/2);
end

function [z, z1] = compute_mom_vec(xinf, xsup)
    n = size(xinf, 1);
    z = zeros(n,1,'double');
    z1 = zeros(n,1,'double');
    for i = 1:n
        [z(i),z1(i)]=compute_mom5d(xinf(i),xsup(i));
    end
end

%{
function [scra1, scra12] = compute_mom3d(xinf, xsup)
    minval = min(abs(xinf), abs(xsup));
    sgn = sign(xinf*xsup);
    if minval <= 8. || sgn <= 0
        Phisup   = Phi(xsup);
        phisup   = phi(xsup);
        Phiinf   = Phi(xinf);
        phiinf   = phi(xinf);
        scra1 = (phiinf - phisup)/(Phisup - Phiinf);
        scra2 = (xinf * phiinf - xsup*phisup)/(Phisup -Phiinf);
    else
        delta2 = (xsup^2 - xinf^2)*0.5;
        if delta2 > 40.
            scra1 = xinf^3/(xinf^2-1);
            scra2 = xinf^4/(xinf^2-1);
        else
            scra1 = (xinf*xsup)^3 * (1. - exp(delta2)) / ( exp(delta2)*(1-xinf^2)*xsup^3 - xinf^3*(1-xsup^2));
            scra2 = (xinf*xsup)^3 * (xsup - xinf*exp(delta2)) / ( exp(delta2)*(1-xinf^2)*xsup^3 - xinf^3*(1-xsup^2));
        end
    end
    scra12=scra2-scra1^2;
end
%}

function [scra1, scra12] = compute_mom5d(xinf, xsup)
    if xsup - xinf < 1e-10
        scra1 = 0.5*(xsup + xinf);
        scra12 = -1;
        return;
    end
    
    if min(abs(xinf), abs(xsup)) <= 6. || xinf*xsup <= 0
        Phisup   = Phi(xsup);
        phisup   = phi(xsup);
        Phiinf   = Phi(xinf);
        phiinf   = phi(xinf);
        scra1 = (phiinf - phisup)/(Phisup - Phiinf);
        scra2 = (xinf * phiinf - xsup*phisup)/(Phisup - Phiinf);        
    else
        delta2 = (xsup^2 - xinf^2)*0.5;
        if delta2 > 40.
            scra1 = xinf^5/(3 - xinf^2 + xinf^4);
            scra2 = xinf^6/(3 - xinf^2 + xinf^4);
        else
            scra1 = (xinf*xsup)^5 * (1. - exp(delta2)) / (-exp(delta2)*(3.0-xinf^2 + xinf^4)*xsup^5 + xinf^5*(3-xsup^2 + xsup^4));
            scra2 = (xinf*xsup)^5 * (xsup - xinf*exp(delta2)) / (-exp(delta2)*(3.0-xinf^2 + xinf^4)*xsup^5 + xinf^5*(3-xsup^2 + xsup^4));
        end
    end
    scra12=scra2-scra1^2;
end
