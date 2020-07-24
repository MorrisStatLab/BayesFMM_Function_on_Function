function res = hwfmm_inits(Y,model,wpspecs,MCMCspecs,initial,sampleU) 
%function res = hwfmm(Y,model,wavespecsy,wavespecsc,pcaspecsx,MCMCspecs,sampleU,get_sigma,twosample) 
%==========================================================================
% Note: hwfmm.m is based on wavfmm4_gbf.m.
% This version samples tau's and pi's and allows for a more general basis
% function for the outcome variable. Here we implement a two phase basis 
% expansion using first wavelets and then PCA decomposition. 
%==========================================================================
%==========================================================================

%% Step 1: Apply DWT to project observed functions into wavelet-pca space
% [D,wavespecsy]  = dwt_rows(Y,wavespecsy);
% fprintf('\n Finished projecting data to wavelet-PCA space.\n \n');

% unlist wpspecs components
wpspecsy        = wpspecs.wpspecsy;
wpspecsx        = wpspecs.wpspecsx;

% extract y-space Kj
% default of length two since each step of hackett decomposition is DWT level 1
Kj      = zeros(2,wpspecsy.dwtAtLevel(wpspecs.nlevels));
k       = wpspecsy.dwtAtLevel(wpspecs.nlevels);
for i = 1:wpspecsy.dwtAtLevel(wpspecs.nlevels)
    Kj(:,i) = wpspecsy.trackerJacker{k}.Kj;
    k       = k+1;
end
Kj      = reshape(Kj,1,2*wpspecsy.dwtAtLevel(wpspecs.nlevels));

%% set up elements for initial values
[model.n,model.p]   = size(model.X);    %%% n=# functions, p=# columns of X
model.c             = size(model.C,2);  %%% assumed to be iid mean zero with a
wpspecs.K   = sum(Kj);
wpspecs.J   = length(Kj);
wpspecs.Kj  = Kj;
wpspecs.T   = wpspecsy.trackerJacker{1}.T;
wpspecs.V   = wpspecsx.trackerJacker{1}.T;
meanop      = uneqkron(Kj)*diag(1./Kj);
expandop    = uneqkron(Kj)';
fsize       = wpspecs.V; %% use V for fsize, may need to change if we explore data reduction in X
wsize       = model.p - fsize;

tic
%% ========================================================================
if isempty(model.Z{1})||(sum(sum(model.Z{1})) == 0)
    model.H     = 0;
    W           = GetW(model,Y);
    [theta_mle,theta_sd,beta_mle,Vbetans,Wv,] = regression_mle(W,model,Y,wpspecs,MCMCspecs);
else
    model.H     = length(model.Z);           %%% H=# of sets of random effects    
    model.m     = zeros(1,model.H);       %%% # of random effects per set %#zhu#% When z{1}=0, model.m=[];
    for h = 1:model.H
        model.m(h)  = size(model.Z{h},2);  %#zhu#% When z{1}=0, this loop will not be excuted since model.H=0
    end
    %%% Compute the initial value of theta using MLE and standard
    %%% deviation.
    model.M     = sum(model.m);
    W           = GetW(model,Y);
    theta0      = chi2rnd(1,model.H+model.c,wpspecs.K);
    [theta_mle,~,theta_sd,~,beta_mle,Vbetans,u_mle,Wv,converge] = mixed_mle(W,theta0,model,Y,wpspecs,MCMCspecs); %#ok<NASGU>   
    initial.u_mle = u_mle;
end

%% save the initial values
betans              = [initial.beta; beta_mle((fsize+1):end,:)];
initial.theta_mle   = theta_mle;
initial.beta_mle    = beta_mle;

%% ========================================================================
%%%% Check whether X matrix is orthogonal
%%%% If it is, then some calculation shortcuts are available later
if (sum(sum(abs(diag(diag(W.XtX))-W.XtX)>1e-15))>0)
    orthogonal = 0;
else
    orthogonal = 1;
end
model.orthogonal_X  = orthogonal;

theta               = theta_mle;
theta_flag          = (theta_mle>=MCMCspecs.VC0_thresh); % If less than VC0_thresh, treat it as zero and don't update in MH.(stepsize=0)
propsd_Theta        = MH_stepsize(theta_sd,theta_mle,MCMCspecs.propsdTheta,0).*theta_flag; % if theta_flag=1, then step size=0.

%% Set initial values for shrinkage hyperparameters pi and tau. 
[tau,PI,alpha,a_tau,b_tau,a_pi,b_pi] = initial_taupi(betans,Vbetans,model,wpspecs,MCMCspecs,meanop);
initial.PI      = PI;
initial.tau     = tau;
initial.a_tau   = a_tau;
initial.b_tau   = b_tau;
initial.a_pi    = a_pi;
initial.b_pi    = b_pi;

%% get PiMat, TauMat, L2, initial value of beta %%
PiMat           = PI*expandop; %#zhu#%  expand from p by J to p by K by repeating within same level.
TauMat          = tau*expandop./Vbetans; %#zhu#%  In this code, do we always assume covT=1?

beta            = betans.*alpha.*(TauMat./(TauMat+1));  %%% Posterior expected values of beta.
Wv.L2           = Get_L2(beta,Wv,wpspecs); %% Update L2 for starting values of beta

%% separate functional from scalar covariates and constrain surface %% DO AS LAST STEP BEFORE SAMPLER??
betah               = wavphc(beta(1:fsize,:), wpspecsx, wpspecsy);
histMark            = ones(size(betah));
histCoef            = wavphc(histMark, wpspecsx, wpspecsy);
scalCoef            = ones(wsize,wpspecs.K);
model.sampCoef      = [histCoef; scalCoef];

%% update betans with constrained initial values %%
beta(1:fsize,:)     = betah;

fprintf('\n Now have starting values for regularization parameters. \n \n');

%% Compute vague proper priors for theta based on empirical Bayes.
[prior_Theta_a,prior_Theta_b] = EmpBayes_Theta(theta_mle,model,MCMCspecs);
fprintf('\n Starting Values Initialized. \n \n'),toc;

%%
p           = model.p;
K           = wpspecs.K;
B           = MCMCspecs.B;
burnin      = MCMCspecs.burnin;
thin        = MCMCspecs.thin;
blocksize   = B;

%%
MCMC_alpha          = NaN(blocksize,p*K); 
MCMC_beta           = NaN(blocksize,fsize*K); %% functional coefficients
MCMC_zeta           = NaN(blocksize,wsize*K); %% scalar coefficients
MCMC_theta          = NaN(blocksize,size(theta,1)*size(theta,2));
MCMC_flag_theta     = NaN(blocksize,K); %% indicator function of acceptance of new set of VCs
if (model.H>0)&&(MCMCspecs.sampleU==1) 
    MCMC_U = NaN(blocksize,sum(model.m)*K);
end

%%
MCMC_tau    = NaN(blocksize,p*wpspecs.J); % in the memory, we only save a block of the samples. 
MCMC_pi     = NaN(blocksize,p*wpspecs.J);

ii              = 0; % counter for MCMC samples to output
tic
%% MCMC Loop: B is desired # samples
for i = 1:(B*thin+burnin)
        
    %% (1) update beta %%
    [beta, gamma, alpha]    = UpdateBetaNoOrthog_hist(beta, Vbetans, PiMat, TauMat, Wv, model, wpspecs, MCMCspecs);
    Wv.L2                   = Get_L2(beta,Wv,wpspecs);
    
    %% (2) Update theta. q_jk, and s_jk %%
    [theta,flag_theta,Vbetans,Wv] = UpdateTheta(beta,theta,Vbetans,Wv,Y,W,model,prior_Theta_a,prior_Theta_b,propsd_Theta,wpspecs,MCMCspecs);   
    
    %% (3) Update U when needed %%
    if (model.H > 0) && (sampleU == 1)
       U = UpdateU(beta,theta,model,Y,wpspecs);
    end
    
    %% (4) Update tau_ijk(as well as TauMat) and PI_ij(as well as PiMat) %%
    tau     = Update_tau(beta,gamma,a_tau,b_tau,meanop);
    TauMat  = tau*expandop./Vbetans;
    PI      = Update_pi(gamma,a_pi,b_pi,wpspecs);
    PiMat   = PI*expandop;
    
    %%
    %%%%%%  Record MCMC samples.%%%%%%    
    if   (i > burnin) && (mod(i-burnin,thin) == 0)   %% Save MCMC samples of beta in single matrix.
            ii                      = ii+1; % this is the real row number among the B samples.
            is                      = mod(ii-1,blocksize)+1; % this is the row number in the block.
            MCMC_beta(is,:)         = reshape(beta(1:fsize,:),1,fsize*K); % record functional covariates
            MCMC_zeta(is,:)         = reshape(beta((fsize+1):end,:),1,wsize*K); % record scalar covariates
            MCMC_alpha(is,:)        = reshape(alpha',1,p*K); % P(gamma == 1)
            MCMC_theta(is,:)        = reshape(theta,1,size(theta,1)*size(theta,2));
            MCMC_flag_theta(is,:)   = flag_theta; 
            MCMC_tau(is,:)          = reshape(tau,1,p*wpspecs.J);   % each block is one j, contains p values.
            MCMC_pi(is,:)           = reshape(PI,1,p*wpspecs.J); % [pi_{11},...pi(1J1);pi_{21},...,pi{2J2};..],J blocks, each block has p values.            
            if (model.H > 0) && (sampleU == 1) 
               MCMC_U(is,:) = reshape(U',1,numel(U));    
            end            
    end
    if mod(i, MCMCspecs.time_update) == 0
       fprintf('\n %d \n',i),toc;
    end
    
end

fprintf('\n Done with MCMC \n');

res.MCMC_beta           = MCMC_beta;   %% functional covariates
res.MCMC_zeta           = MCMC_zeta;   %% intercept plus scalar covariates (if any)
res.MCMC_alpha          = MCMC_alpha;
res.MCMC_theta          = MCMC_theta;
res.MCMC_flag_theta     = MCMC_flag_theta;
res.MCMC_tau            = MCMC_tau;
res.MCMC_pi             = MCMC_pi;
res.theta               = theta;
res.model               = model;
res.initial             = initial;
res.wpspecs             = wpspecs;
%% here: 1/16/14 %%


% fprintf('\n  Now projecting results back to data space. \n')
% 
% %%%%%  Project MCMC samples back to data space
% 
% if sampleU==0        
%     postout = PostProcess_FFR(MCMC_beta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wavespecsc,wavespecsy,pcaspecsx,twosample,get_sigma,sampleU);
% else
%     postout = PostProcess_FFR(MCMC_beta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wavespecsc,wavespecsy,pcaspecsx,twosample,get_sigma,sampleU,MCMC_U);
% end
% res.MCMC_beta   = MCMC_beta;
% res.MCMC_theta  = MCMC_theta;
% res.MCMC_tau    = MCMC_tau;
% res.MCMC_pi     = MCMC_pi;
% if (model.H > 0) && (sampleU == 1)
%      res.MCMC_U = MCMC_U;
% end
% 
% %%%%% put results into "res" %%%%
% res.model       = model;
% res.wavespecsy  = wavespecsy;
% res.wavespecsx  = wavespecsc;
% res.pcaspecsx   = pcaspecsx;
% res.MCMCspecs   = MCMCspecs;
% res.theta_flag  = theta_flag; % theta_flag=0 indicates treating the theta component as 0 all the time. 
% res.betans      = betans; %#zhu#% note that this betans has not been updated during MCMC.
% res.Vbetans     = Vbetans; %#zhu#% note that this Vbetans has been updated during MCMC.
% 
% res.bstarhat    = postout.bstarhat;
% res.bstar025CI  = postout.bstar_025CI;
% res.bstar975CI  = postout.bstar_975CI;
% res.alphahat    = postout.alphahat;
% res.psi         = postout.psi;
% res.blocksize   = blocksize; % the block size: how many interations of MCMC samples we save to text each time.
% res.paramroute  = paramroute;
% 
%     res.bhat        = postout.bhat;
%     res.Q025_bhat   = postout.Q025_bhat;
%     res.Q975_bhat   = postout.Q975_bhat;
%     res.ahat        = postout.ahat;
%     res.Q025_ahat   = postout.Q025_ahat;
%     res.Q975_ahat   = postout.Q975_ahat;
% 
% res.thetahat    = postout.thetahat;
% res.theta_025CI = postout.theta_025CI;
% res.theta_975CI = postout.theta_975CI;
% res.acpt_rate   = postout.accept_rate_theta;
% res.theta_MLE   = theta_mle;
% res.theta_MLEsd = theta_sd;
% 
% res.tauhat      = postout.tauhat;
% res.pihat       = postout.pihat;
% 
% if (model.H > 0) && (sampleU == 1)
%     % res.MCMC_U=MCMC_U;
%     res.Uhat       = postout.uhat;
%     res.U_025CI    = postout.u_025CI;
%     res.U_975CI    = postout.u_975CI;
%     res.g_Ut       = postout.g_Ut;
%     res.g_Ut025    = postout.g_Ut025;
%     res.g_Ut975    = postout.g_Ut975;
% end
% if get_sigma == 1
%     res.Sigma = postout.Sigma;
% end
% %==========================================================================
% %res.Y=Y; % we actually don't have to save Y.
% %==========================================================================
% res.Wv          = Wv; %#zhu#% this is the Wv of very last iteration.
% res.D           = D;
% res.initial     = initial;
