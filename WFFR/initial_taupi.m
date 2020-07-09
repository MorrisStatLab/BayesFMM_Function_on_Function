function [tau,PI,alpha,a_tau,b_tau,a_pi,b_pi]=initial_taupi(betans,Vbetans,model,wavespecs,MCMCspecs,meanop)
% This function compute the initial values of tau and pi. we do not let T
% depend on ij, but let T
%           to be T_ijk. so that tau_ijk can be computed.
%
%  Input: 
%         betans-- the betans estimate.
%         Vbetans-- the variance of betans.
%         p -- matrix containing starting values for pi, pi indexed by ij.
%
%
%  Functions needed: repvec.m : same as Splus function rep applied to vectors
%
%%#zhu#% This code is written by assuming that \pi_{ij} depend on
%%only i and j, not k. \tau_{ijk}.
% see my detailed derivation. T_{ijk}=tau_{ijk}/V_{ijk}.

% (1) Find initial values of tau, PI.
J_temp=1:wavespecs.J;
pstart0=rep(repmat(0.8*2.^(-(J_temp-1)),model.p,1),wavespecs.Kj);
Zeta=betans./sqrt(Vbetans);
tau0=Vbetans.*max((Zeta).^2-1,1e-6);
BigTau=tau0./Vbetans; 
O=min(pstart0./(1-pstart0).*(1+BigTau).^(-.5).*exp(Zeta.^2/2.*(BigTau./(BigTau+1))),MCMCspecs.maxO);
alpha=O./(O+1);
%meanop=uneqkron(wavespecs.Kj)*diag(1./wavespecs.Kj);
PI=alpha*meanop; 

% n_nosmooth=MCMCspecs.nj_nosmooth;
% if n_nosmooth>0    
%    PI(:,1:n_nosmooth)=repmat(1-MCMCspecs.minp,model.p,n_nosmooth); %%% Set P===1
%    PI=max(min(PI,1-MCMCspecs.minp),MCMCspecs.minp);   
% end

% (2) Now determine Hypo-parameters for Inv-Gamma prior for tau_{ij}, and Beta
% prior for PI_{ij}
tau=tau0*meanop;
tau(tau<1e-6)=1e-6; % try to avoid 0 values of tau.
[a_tau,b_tau]=InvGamma_param(tau,MCMCspecs.tau_prior_var,0,linspace(0.01,30,100));       

% (3) Now determine Hypo-parameters for beta priors for PI.
%[a_pi,b_pi]=Beta_para(PI,MCMCspecs.PI_prior_var,0,0,linspace(0.001,1,100));
a_pi=PI;b_pi=1-PI;