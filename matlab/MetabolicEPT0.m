function [mu, s, a,d, av, va, Cov, t] = MetabolicEPT0(S, b, nuinf, nusup, damp, max_iter, minvar, maxvar, precision, precision_lin)

%%% Expectation Propagation algorithm for metabolic networks (implementation for \beta -> +\infty)
%%%
%%% Input
%%% - S: stoichiometric matrix of "Nm metabolites" x "Nr reactions".
%%% - b: vector of Nm intakes/uptakes
%%% - nuinf, nusup: lower and upper bounds for each metabolic flux
%%% - damp: damping coefficient (from 0 to 1) applied to the update of means "a" and variances "d" of approximating Gaussians
%%%	    Ex. "new a" = damp * "new a" + (1 - damp) * "old a"
%%% - max_iter: maximum number of iterations
%%% - minvar, maxvar: lower and upper bounds for the variances "d" of the approximation. (ex. 1e-50, 1e50)
%%% - precision:  precision required to stop the algorithm (ex. 1e-9)
%%% - precision_lin: precision required on the linear constraints: (S*nu - b)*(S*nu - b)/Nr <= precision_lin
%%%	
%%% Output
%%% - mu: vector parametrizing the mean of the posterior distribution
%%% - s: vector parametrizing the variance of the posterior distribution
%%% - a: vector containing the means of the approximated priors
%%% - d: vector containing the variances of the approximated priors
%%% - av: averages of the truncated Gaussians of the approximation
%%% - va: variances of the truncated Gaussians of the approximation
%%% - Cov: covariance matrix of the fluxes
%%% - t: running time
%%%
%%% The marginal probability density of the n-th flux is a truncated Gaussian N(mu(n),s(n)) in the interval [nuinf, nusup]. It has average av(n) and variance va(n)

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
err = 100;
iter = 0;
err_avva = 100;
err_lin = 100;

% row echelon form
[A, idxd] = rref(cat(2, S, b));
Nd = length(idxd);
A = A(1:Nd,:);
idxf = setdiff(1:Nr,idxd);
Nf = length(idxf);
C = A(:,idxf);
Cp = C';
y = A(:,end);
basis = zeros(Nr,Nf);
basis(idxd,:) = -A(:,idxf);
basis(idxf,:) = eye(Nf, Nf);
Df = sparse(1:Nf, 1:Nf, 1./d(idxf));
Dd = sparse(1:Nd, 1:Nd, 1./d(idxd));

for i = 1:Nr
    if(ismember(i,idxd))
        posd(i) = find(idxd == i);
        posf(i) = 0;
    else
        posd(i) = 0;
        posf(i) = find(idxf == i);
    end
end

tic
while ((err_avva > precision || err_lin > precision_lin) && iter < max_iter)
    iter = iter + 1;

    I1 = inv(Df + Cp*Dd*C);
    v = (I1) * (Df*a(idxf) + Cp*Dd*(y - a(idxd)));
    I = diag(I1);
    oldvar = va;
    oldav = av;
    for i = 1:Nr
        if(ismember(i,idxf))
            idx = posf(i);
            ss = min(maxvar, max(minvar, I(idx)));
            vv = v(idx);
        else	
            idx = posd(i);
            x = -C(idx,:)';
            ss = min(maxvar, max(minvar, x'*I1*x));
            vv = x'*v + y(idx);
        end
        s(i) = min(maxvar, max(minvar, 1/(1/ss - 1/d(i))));
        mu(i) = s(i) * (vv/ss - a(i)/d(i));
        mu(i) = (vv-(a(i)*ss)/d(i))/(1-ss/d(i));
        if(ss == d(i))
    		mu(i) = 0.5*(nuinf(i) + nusup(i));
        end
        
        s05(i) = s(i)^0.5;
        x0 = (nuinf(i)-mu(i))/s05(i); 
        x1 = (nusup(i)-mu(i))/s05(i);     
        [z,eps] = compute_mom_vec(x0, x1);    
        av(i) = mu(i) + z * s05(i);
        va(i) = max(0, s(i) * (1 + eps));
    end

    err_avva = max(max(abs(av-oldav)), max(abs(va-oldvar)));
    err_lin = 1/Nr * (S*av - b)'*(S*av - b);

    % moments matching (update "a" and "d")
    new_d = 1./(1./va - 1./s);
    new_d = min(maxvar, max(minvar, new_d));
    new_a = av + (av - mu) .* new_d.* (1./s);    
    olda = a;
    oldd = d;
    a = damp * a + (1-damp) * new_a;
    d = damp * d + (1-damp) * new_d;
    Df = sparse(1:Nf, 1:Nf, 1./d(idxf));
    Dd = sparse(1:Nd, 1:Nd, 1./d(idxd));
    %if(~(mod(iter,100)))
    fprintf('it:%i err_lin:%e conv_EP:%e\n', iter, err_lin, err_avva);
    %end
end
t = toc;
mu = mu * factor;
s = s * factor^2;
a = a*factor;
d = d*factor^2;
av = av*factor;
va = va*factor^2;
nusup = nusup * factor;
nuinf = nuinf * factor;
b = b * factor;
Df = sparse(1:Nf, 1:Nf, 1./d(idxf));
Dd = sparse(1:Nd, 1:Nd, 1./d(idxd));
M = Df + Cp*Dd*C;
Cov = inv(M);
Cov = basis * Cov * basis';

if(iter < max_iter && ~(isnan(err_lin)))
    fprintf('EP converged\nit:%i err:%e conv_EP:%e\n', iter, err_lin, err_avva);
else
    fprintf('EP has not converged\nit:%i err:%e conv_EP:%e\n', iter, err_lin, err_avva);
end


end

function R = funcg(x)

    R = 1/x - 1/x^3 + 3/x^5 - 15/x^7 + 105/x^9;
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

















