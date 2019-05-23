%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% Run Expectation Propagation algorithm for COBRA models of metabolic networks
%%%
%%% 1) Load and put your data in "model" variable (MATLAB format is required) 
%%% 2) Pre-process your model 
%%% 3) Choose input parameters for EP algorithm (see MetabolicEP m-file for a detailed description)
%%% 4) Run this script
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	


load('../test/ecoli_core_model.mat');
pmodel = pre_processing(model);
Beta=1e10;
damping=0.9;
precision=1e-6;
maxit=2000;
minvar=1e-50;
maxvar=1e50;
av_exp = 0;
va_exp = 0;
exp_i = 0;

[mu, s, a, d, av, va, Cov, t_EP]  = MetabolicEP(full(pmodel.S),pmodel.b,pmodel.lb,pmodel.ub,Beta, damping, maxit, minvar, maxvar, precision,  av_exp, va_exp, exp_i);

% plot Biomass marginal
idx_bm = strmatch('Biomass_Ecoli_core_w_GAM', pmodel.rxns);
color = [0, 0 ,1];
plot_fluxmarginal(-1e3, 1e3, mu(idx_bm), s(idx_bm), pmodel.lb(idx_bm),pmodel.ub(idx_bm), ...
                av(idx_bm) , sqrt(va(idx_bm)), color );
%% Beta -> \infty implementation

precision_lin = 1e-7;
[muT0, sT0, aT0, dT0, avT0, vaT0, CovT0, t_EPT0] = MetabolicEPT0(full(pmodel.S), pmodel.b, pmodel.lb, pmodel.ub, damping, maxit, minvar, maxvar, precision, precision_lin);

plot_fluxmarginal(-1e3, 1e3, muT0(idx_bm), sT0(idx_bm), pmodel.lb(idx_bm),pmodel.ub(idx_bm), ...
                avT0(idx_bm) , sqrt(vaT0(idx_bm)), color );

%% Ex: fix biomass flux of E.Coli 

clear all
% uncontrained run
load('../data/Ec_iJR904.mat')
index_glc = strmatch('D Glucose exchange', Ec_iJR904.rxnNames);
Ec_iJR904.lb(index_glc) = -43;
model = pre_processing(Ec_iJR904);
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







