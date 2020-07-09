function res = opwavfm_old(Y,model,wavespecsy,wavespecsc,pcaspecsx,MCMCspecs) 
%function res = opwavfm(Y,model,wavespecsy,wavespecsc,pcaspecsx,MCMCspecs,get_sigma,saveraw) 
%==========================================================================
% Note: wavfmm4_gbf.m is based on wavfmm4_update.m.
% This version samples tau's and pi's and allows for a more general basis
% function for the outcome variable. Here we implement a two phase basis 
% expansion using first wavelets and then PCA decomposition. 
%==========================================================================
%==========================================================================

%% get parameters %%
[model.n,model.p]   = size(model.X);   %%% n=# functions, p=# covariates
model.c             = size(model.C,2);               %%% assumed to be iid mean zero with a
n                   = model.n;
w                   = size(model.W,2);
V                   = model.V;
T                   = size(Y,2);
fsize               = model.p - (1 + w);

%% Step 1: set up levels and determine Zstar^(0)
levels          = sort(unique(Y)); % find levels of Y
model.levels    = levels;
L               = length(levels);
Yvec            = reshape(Y,1,n*T);
[~,I]           = sort(Yvec);
Ysort           = I/(size(Yvec,2)+1);
Z               = norminv(Ysort,0,1);
Z               = reshape(Z,n,T);
initial.Z0      = Z;

%% set priors for cutpoints
muk     = zeros(1,L-1);
sigk    = 100*ones(1,L-1);

%% Step 2: Apply DWT to project Zstar into wavelet space (do this outside function?)
[D,wavespecsy]  = dwt_rows(Z,wavespecsy);
fprintf('\n Finished projecting data to wavelet-PCA space.\n \n');

wavespecsy.K    = sum(wavespecsy.Kj);
meanop          = uneqkron(wavespecsy.Kj)*diag(1./wavespecsy.Kj);
expandop        = uneqkron(wavespecsy.Kj)';

tic

%% set initial values for wavfm %%
model.H     = 0;
W           = GetW(model,D);
[theta_mle,theta_sd,betans,Vbetans,Wv,] = regression_mle(W,model,D,wavespecsy,MCMCspecs);

%% save the initial values %%
initial.theta_mle   = theta_mle;
initial.beta_mle    = betans;

%% check whether X matrix is orthogonal %%
%  if it is, then some calculation shortcuts are available later
if (sum(sum(abs(diag(diag(W.XtX))-W.XtX)>1e-15))>0)
    orthogonal = 0;
else
    orthogonal = 1;
end
model.orthogonal_X  = orthogonal;

theta               = theta_mle;
theta_flag          = (theta_mle>=MCMCspecs.VC0_thresh); % If less than VC0_thresh, treat it as zero and don't update in MH.(stepsize=0)
propsd_Theta        = MH_stepsize(theta_sd,theta_mle,MCMCspecs.propsdTheta,0).*theta_flag; % if theta_flag=1, then step size=0.

%% set initial values for shrinkage hyperparameters pi and tau %%
[tau,PI,alpha,a_tau,b_tau,a_pi,b_pi] = initial_taupi(betans,Vbetans,model,wavespecsy,MCMCspecs,meanop);
initial.PI      = PI;
initial.tau     = tau;
initial.a_tau   = a_tau;
initial.b_tau   = b_tau;
initial.a_pi    = a_pi;
initial.b_pi    = b_pi;

%% set regularization matrices %%
PiMat           = PI*expandop; %#zhu#%  expand from p by J to p by K by repeating within same level.
TauMat          = tau*expandop./Vbetans; %#zhu#%  In this code, do we always assume covT=1?

bstar           = betans.*alpha.*(TauMat./(TauMat+1));  %%% Posterior expected values of beta.
Wv.L2           = Get_L2(bstar,Wv,wavespecsy); %% Update L2 for starting values of beta
fprintf('\n Now have starting values for regularization parameters. \n \n');

%% compute vague proper priors for theta based on empirical Bayes %%
[prior_Theta_a,prior_Theta_b] = EmpBayes_Theta(theta_mle,model,MCMCspecs);
fprintf('\n Starting Values Initialized. \n \n'),toc;

p           = model.p;
K           = wavespecsy.K;
B           = MCMCspecs.B;
burnin      = MCMCspecs.burnin;
thin        = MCMCspecs.thin;
blocksize   = B;

%% set matrices for sampler %%
MCMC_alpha          = NaN(blocksize,p*K); 
MCMC_bstar          = NaN(blocksize,p*K);
MCMC_theta          = NaN(blocksize,size(theta,1)*size(theta,2));
MCMC_flag_theta     = NaN(blocksize,K); %% indicator function of acceptance of new set of VCs
MCMC_beta           = NaN(blocksize,V*T);

MCMC_tau    = NaN(blocksize,p*wavespecsy.J); % in the memory, we only save a block of the samples. 
MCMC_pi     = NaN(blocksize,p*wavespecsy.J);

MCMC_cuts   = NaN(blocksize,L-1);
% MCMC_Z      = NaN(blocksize,n*K);

ii              = 0; % counter for MCMC samples to output
tic

%% MCMC loop: B is desired # samples
for i = 1:(B*thin+burnin)
    
    %% (1) update cutpoints
    cuts    = UpdateCuts(Y, Z, muk, sigk, L);

    %% (2) update beta
    [bstar,gamma,alpha]     = UpdateBetaNoOrthog(bstar, Vbetans, PiMat, TauMat, Wv, model, wavespecsy, MCMCspecs);
    Wv.L2                   = Get_L2(bstar, Wv, wavespecsy); %#zhu#% beta is updated here, so does L2
    
    %% (3) update theta. q_jk, and s_jk, D should be Zi* also might need to update W
    [theta,flag_theta,Vbetans,Wv]   = UpdateTheta(bstar, theta, Vbetans, Wv, D, W, model, prior_Theta_a, prior_Theta_b, propsd_Theta, wavespecsy, MCMCspecs);
    
    %% (4) update tau_ijk(as well as TauMat) and PI_ij(as well as PiMat)
    tau     = Update_tau(bstar, gamma, a_tau, b_tau, meanop);
    TauMat  = tau*expandop./Vbetans;
    PI      = Update_pi(gamma, a_pi, b_pi, wavespecsy);
    PiMat   = PI*expandop;
    
    %% (5) project beta and theta back into data-space to get Z from Z* (D)
    % SAVE projected beta and theta so we don't have to do it later
    bstem           = bstar(1:fsize,:);
    [beta, sigma]   = BTPost(bstem, theta, wavespecsy, wavespecsx, pcaspecsx);
    
    %% (6) update latent variable Z, need s_jk from UpdateTheta
    [Z, W]  = UpdateZ(Y, model, beta, sigma, cuts);
    
    %% (7) perform DWT on Z
    [D, ~]  = dwt_rows(Z, wavespecsy);
    
    %% update W, may or may not need to do this
    % W       = GetW(model,D);
    
    %%  record MCMC samples %%
    if   (i > burnin) && (mod(i-burnin,thin) == 0)   %% Save MCMC samples of beta in single matrix
            ii                      = ii+1; % this is the real row number among the B samples
            is                      = mod(ii-1,blocksize)+1; % this is the row number in the block
            MCMC_bstar(is,:)        = reshape(bstar',1,p*K); %#zhu#% each row=[beta(1,j,k),j,k|beta(2,j,k),j,k|...|beta(p,j,k),j,k]
            MCMC_beta(is,:)         = reshape(beta',1,V*T);
            MCMC_alpha(is,:)        = reshape(alpha',1,p*K);
            MCMC_theta(is,:)        = reshape(theta,1,size(theta,1)*size(theta,2));
            MCMC_flag_theta(is,:)   = flag_theta; 
            MCMC_tau(is,:)          = reshape(tau,1,p*wavespecsy.J);   % each block is one j, contains p values
            MCMC_pi(is,:)           = reshape(PI,1,p*wavespecsy.J); % [pi_{11},...pi(1J1);pi_{21},...,pi{2J2};..],J blocks, each block has p values
            MCMC_cuts(is,:)         = cuts;
%             MCMC_Z(is,:)            = reshape(Z,1,n*K);
    end
    
    %% print sampler updates %%
    if mod(i, MCMCspecs.time_update) == 0
       fprintf('\n %d \n',i),toc;
    end
end

fprintf('\n Done with MCMC \n');
fprintf('\n  Now projecting results back to data space. \n')

%% project MCMC samples back to data space %%
%  Will need to edit post processor         %
% postout = PostProcess_FFR(MCMC_beta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wavespecsc,wavespecsy,pcaspecsx,twosample,get_sigma);


%% save results %%
res.MCMC_beta   = MCMC_beta;
res.MCMC_theta  = MCMC_theta;
res.MCMC_tau    = MCMC_tau;
res.MCMC_pi     = MCMC_pi;

%% put results into "res" %%
res.model       = model;
res.wavespecsy  = wavespecsy;
res.wavespecsx  = wavespecsc;
res.pcaspecsx   = pcaspecsx;
res.MCMCspecs   = MCMCspecs;
res.theta_flag  = theta_flag; % theta_flag=0 indicates treating the theta component as 0 all the time
res.betans      = betans; %#zhu#% note that this betans has not been updated during MCMC
res.Vbetans     = Vbetans; %#zhu#% note that this Vbetans has been updated during MCMC

res.bstarhat    = postout.bstarhat;
res.bstar025CI  = postout.bstar_025CI;
res.bstar975CI  = postout.bstar_975CI;
res.alphahat    = postout.alphahat;
res.psi         = postout.psi;
res.blocksize   = blocksize; % the block size: how many interations of MCMC samples we save to text each time

res.bhat        = postout.bhat;
res.Q025_bhat   = postout.Q025_bhat;
res.Q975_bhat   = postout.Q975_bhat;
res.ahat        = postout.ahat;
res.Q025_ahat   = postout.Q025_ahat;
res.Q975_ahat   = postout.Q975_ahat;
res.pwCI        = postout.pwCI;
res.MAPs        = postout.MAPs;
res.UMAPs       = postout.UMAPs;
res.LMAPs       = postout.LMAPs;

res.thetahat    = postout.thetahat;
res.theta_025CI = postout.theta_025CI;
res.theta_975CI = postout.theta_975CI;
res.acpt_rate   = postout.accept_rate_theta;
res.theta_MLE   = theta_mle;
res.theta_MLEsd = theta_sd;

res.tauhat      = postout.tauhat;
res.pihat       = postout.pihat;

%==========================================================================
%res.Y=Y; % we actually don't have to save Y.
%==========================================================================
res.Wv          = Wv; %#zhu#% this is the Wv of very last iteration
res.D           = D;
res.initial     = initial;
