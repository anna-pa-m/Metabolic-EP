%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Run Expectation Propagation algorithm for COBRA models of metabolic networks
%%%
%%% 1) Load and put your data in "model" variable (MATLAB format is required) 
%%% 2) Choose input parameters for EP algorithm (see MetabolicEP m-file for a detailed description)
%%% 3) Run this script
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	



model = Ec_iJR904;
Beta=1e7;
damping=0.9;
precision=1e-6;
maxit=2000;
minvar=1e-50;
maxvar=1e50;
av_exp = 0;
va_exp = 0;
exp_i = 0;


[mu, s, a, d, av, va, t_EP]  = MetabolicEP(full(model.S),model.b,model.lb,model.ub,Beta, damping, maxit, minvar, maxvar, precision,  av_exp, va_exp, exp_i);

%% fix biomass flux of E.Coli 

% uncontrained run
model = Ec_iJR904;
index_glc = strmatch('D Glucose exchange', model.rxnNames);
model.lb(index_glc) = -43;
exp_i = 0;
av_exp = 0;
va_exp = 0;

Beta=1e7;
damping=0.9;
precision=1e-6;
maxit=2000;
minvar=1e-50;
maxvar=1e50;

[mu_free, s_free, a_free,d_free, av_free, va_free, t_EP_free]  = MetabolicEP(full(model.S),model.b,model.lb,model.ub,Beta, damping, maxit, minvar, maxvar, precision, av_exp, va_exp, exp_i);

% constrained run
model = Ec_iJR904;
index_glc = strmatch('D Glucose exchange', model.rxnNames);
model.lb(index_glc) = -43;
exp_i = strmatch('BiomassEcoli', model.rxnNames);
av_exp = 0.92;
va_exp = 0.0324;

Beta=1e7;
damping=0.9;
precision=1e-5;
maxit=2000;
minvar=1e-50;
maxvar=1e50;
[mu, s, a, d, av, va, t_EP]  = MetabolicEP(full(model.S),model.b,model.lb,model.ub,Beta, damping, maxit, minvar, maxvar, precision, av_exp, va_exp, exp_i);

%% check 
va_aux = (s(exp_i) * d(exp_i))/(s(exp_i) + d(exp_i));
av_aux = (mu(exp_i) * d(exp_i) + a(exp_i)* s(exp_i))/(s(exp_i) + d(exp_i));

[av_exp av_aux va_exp va_aux]







