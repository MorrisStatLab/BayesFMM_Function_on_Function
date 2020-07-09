%% Applied the wavfmm4_v2_zhu model to colon cancer data.
%%%%run1: IG prior for tau depend on ij.
addpath('/Users/markmeyer/Dropbox/Research/Dissertation/WFMM4');

%%
load('/Users/markmeyer/Dropbox/Research/Dissertation/WFMM4/coloncancer_data2.mat');

%%
MCMCspecs.B=100;
MCMCspecs.burnin=200;
MCMCspecs.thin=1;
MCMCspecs.propsdTheta=1.3;

%% MCMCspecs.nj_nosmooth=1;  % can be 0, if is 0, don't do any constraint for pi_{ij}, bigTau_{ij}.
MCMCspecs.minp=1e-14;
MCMCspecs.maxO=1e20;
MCMCspecs.minVC=1e-20;

%%
MCMCspecs.VC0_thresh=1e-4;
MCMCspecs.delta_theta=1e-4;
MCMCspecs.thetaMLE_maxiter=10^6;
MCMCspecs.EmpBayes_maxiter=10^6;
MCMCspecs.time_update=100;

%%
wavespecs.nlevels=8;
wavespecs.wavelet='db8';
wavespecs.wtmode='sym';

%%
MCMCspecs.tau_prior_var =1e3; % the variance of tau_{ijk} when finding prior parameters for tau_{ijk}.
MCMCspecs.tau_prior_idx =1; % 1 indicate that a_tau and b_tau depend on ij, 0 indicate that they depend on jk. 
MCMCspecs.PI_prior_var  =0.06; % this range should be in [0.02 0.09].
MCMCspecs.sampleU       = 1;

%%
model.C=ones(size(Y,1),1);
model.Hstar=0;

%% Function Test set up
sampleU     = MCMCspecs.sampleU;
get_sigma   = 1;
blocksize   = NaN;
paramroute  = false; %'/Users/markmeyer/Dropbox/Research/Dissertation/Paper 1/MCMC Samples'

%%
tic;
res=wavfmm4_v5_zhu(Y,model,wavespecs,MCMCspecs,sampleU,get_sigma,blocksize,paramroute); 
time_spend=toc;

%%
save('C:\Documents and Settings\hzhu1\My Documents\MATLAB\RobustWFMM\wavfmm4_v2_zhu\Colon\colon_wav4v2_result1.mat');

%% Applied the wavfmm4_v2_zhu model to colon cancer data.
%%%%run2: IG prior for tau depend on jk.
clear;
addpath('/Users/markmeyer/Dropbox/Research/Dissertation/WFMM4');
load('/Users/markmeyer/Dropbox/Research/Dissertation/WFMM4/coloncancer_data2.mat');

%%
MCMCspecs.B=100;
MCMCspecs.burnin=200;
MCMCspecs.thin=5;
MCMCspecs.propsdTheta=1.3;

%MCMCspecs.nj_nosmooth=1;  % can be 0, if is 0, don't do any constraint for pi_{ij}, bigTau_{ij}.
MCMCspecs.minp=1e-14;
MCMCspecs.maxO=1e20;
MCMCspecs.minVC=1e-20;

MCMCspecs.VC0_thresh=1e-4;
MCMCspecs.delta_theta=1e-4;
MCMCspecs.thetaMLE_maxiter=10^6;
MCMCspecs.EmpBayes_maxiter=10^6;
MCMCspecs.time_update=100;

wavespecs.nlevels=8;
wavespecs.wavelet='db8';

MCMCspecs.tau_prior_var=1e3; % the variance of tau_{ijk} when finding prior parameters for tau_{ijk}.
MCMCspecs.tau_prior_idx=0; % 1 indicate that a_tau and b_tau depend on ij, 0 indicate that they depend on jk. 
MCMCspecs.PI_prior_var=0.06; % this range should be in [0.02 0.09].


model.C=ones(size(Y,1),1);
model.Hstar=0;

tic;
res=wavfmm4_v5_zhu(Y,model,wavespecs,MCMCspecs,1); 
time_spend=toc;

save('C:\Documents and Settings\hzhu1\My Documents\MATLAB\RobustWFMM\wavfmm4_v2_zhu\Colon\colon_wav4v2_result2.mat');

%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Analysis of the above 2 results.
% result1: priors for tau_{ijk} depend on ij. 
clear;
load('C:\Documents and Settings\hzhu1\My Documents\MATLAB\RobustWFMM\wavfmm4_v2_zhu\Colon\colon_wav4v2_result1.mat');
addpath('C:\Documents and Settings\hzhu1\My Documents\MATLAB\RobustWFMM\wavfmm8_v2_zhu\Colon');
addpath('C:\Documents and Settings\hzhu1\My Documents\MATLAB\RobustWFMM\wavfmm4_v2_zhu');
position=linspace(0,1,256);
coloncancer_result_ana(res,position,1); % this function need to be fixed for the new results. add in betans.
cd('C:\Documents and Settings\hzhu1\My Documents\MATLAB\RobustWFMM\wavfmm4_v2_zhu\Colon\Result1');

% result2: priors for tau_{ijk} depend on jk.
clear;
load('Z:\Matlabwork\wavfmm_zhu\wavfmm4_v2_zhu\Colon\colon_wav4v2_result2.mat');
addpath('C:\Documents and Settings\hzhu1\My Documents\MATLAB\RobustWFMM\wavfmm8_v2_zhu\Colon');
addpath('C:\Documents and Settings\hzhu1\My Documents\MATLAB\RobustWFMM\wavfmm4_v2_zhu');
position=linspace(0,1,256);
coloncancer_result_ana(res,position,1); % this function need to be fixed for the new results. add in betans.
cd('C:\Documents and Settings\hzhu1\My Documents\MATLAB\RobustWFMM\wavfmm4_v2_zhu\Colon\Result2');
