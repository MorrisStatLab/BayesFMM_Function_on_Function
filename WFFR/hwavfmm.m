function res = hwavfmm(Y,model,wavespecsy,wavespecsc,pcaspecsx,MCMCspecs,sampleU,get_sigma,twosample,saveraw,blocksize,paramroute) 
%function res = wavfmm4_gbf(Y,model,wavespecsy,wavespecsc,pcaspecsx,MCMCspecs,sampleU,get_sigma,twosample,blocksize,paramroute) 
%==========================================================================
% Note: wavfmm4_gbf.m is based on wavfmm4_update.m.
% This version samples tau's and pi's and allows for a more general basis
% function for the outcome variable. Here we implement a two phase basis 
% expansion using first wavelets and then PCA decomposition. 
%==========================================================================
% Input: Y     = (n x T) data matrix; each row is one observed function on 
%                 an equally spaced grid of size T 
%
%        model = structure containing model information, with elements:
%                X = (n x p) design matrix for fixed effects functions;
%                Z = cell array of H matrices, each containing (n x m_h) design
%                matrix for set of random effects functions;
%                H = # of levels of random effects
%                Hstar = # of levels for which columns are conditionally
%                independent? (so indicate conditionally independent blocks)
%                Q: Does conditional independence also condition on
%                random effects at higher-up levels?
%                V = n x v matrix? = Functions with 1's in given column of v
%                share same residual covariance matrix S, so v = number
%                of unique residual covariance matrices estimated.
%   wavespecsy = structure containing wavelet settings, including:
%                wavelet = wavelet basis to use
%                nlevels = number of levels of decomposition to do
%   MCMCspecs = parameters for MCMC iterations.  
%                minVC = minimum size of VC (tolerance limit) (in EmpBayes) default 0.001
%                epsilon = threshhold on MOM est. of VC below which VC is assumed to be 0 (in EmpBayes) default 0.001
%                delta_theta = multiple for prior on theta: "number of datasets of information" in prior (in
%                             wavefmm4) default 0.0001
%                xi_theta = multiple of MOM est of VC to use as proposal variance
%                                   (in wavefmm4) -- default 1.5
%      sampleU = if 1, sample random effect coefficient U. If 0, do not
%                sample U.
%    get_sigma = if 1, we compute the covariance estimation for U_(t), and
%                E(t) by idwt during the post-process. 
%   export2txt =
%   exportsize =
%   paramroute =
% Output:   res = data structure with elements: 
%           wavelet = wavelet basis used
%           nlevels = number of wavelet levels
%           betans = nonshrunken mean wavelet coefficient estimates (p x T matrix)
%           Vbetans = variance of nonshrunken wavelet coefficient estimates (p x T)
%           pi = empirical bayes estimates of proportion of nonzero wavelet coefficients per level (p x nlevels+1)
%           Tau = empirical bayes estimates of variance of nonzero wavelet coefficients per level (p x nlevels+1)
%           Gamma = matrix of posterior probabilities of nonzero mean-level coefficients (p x T)
%           betahat = matrix of shrinkage estimates of mean functions' wavelet coefficient (p x T)
%           ghat = matrix of smoothed mean function estimates (p x T)
%           ghatns = matrix of unsmoothed mean function estimates (p x T)
%           theta_q = matrix (H x *) of variance components
%           theta_s = vector (1 x *) of residual variance components
%                           *=1,J, or K if q(s)form=0, 1, or 2,
%                           respectively
%==========================================================================

%%% Check inputs %%%
if nargin < 12;
    paramroute  = false;
end;

if nargin < 11;
    blocksize   = NaN;
end;

if nargin <10;
    saveraw     = 0;
end;

%%% Step 1: Apply DWT to project observed functions into wavelet-pca space
[D,wavespecsy]  = dwt_rows(Y,wavespecsy);
fprintf('\n Finished projecting data to wavelet-PCA space.\n \n');

wavespecsy.K    = sum(wavespecsy.Kj);
meanop          = uneqkron(wavespecsy.Kj)*diag(1./wavespecsy.Kj);
expandop        = uneqkron(wavespecsy.Kj)';

tic
%==========================================================================
[model.n,model.p]       = size(model.X);   %%% n=# functions, p=# covariates
model.c                 = size(model.C,2);               %%% assumed to be iid mean zero with a
if isempty(model.Z{1})||(sum(sum(model.Z{1})) == 0)
    model.H     = 0; 
    W           = GetW(model,D);   
    [theta_mle,theta_sd,betans,Vbetans,Wv,] = regression_mle(W,model,D,wavespecsy,MCMCspecs);    
else
    model.H     = length(model.Z);           %%% H=# of sets of random effects    
    model.m     = zeros(1,model.H);       %%% # of random effects per set %#zhu#% When z{1}=0, model.m=[];
    for h = 1:model.H
        model.m(h)  = size(model.Z{h},2);  %#zhu#% When z{1}=0, this loop will not be excuted since model.H=0
    end
    %%% Compute the initial value of theta using MLE and standard
    %%% deviation.
    model.M     = sum(model.m);
    W           = GetW(model,D);
    theta0      = chi2rnd(1,model.H+model.c,wavespecsy.K);
    [theta_mle,theta_var,theta_sd,beta_mle,betans,Vbetans,u_mle,Wv,converge] = mixed_mle(W,theta0,model,D,wavespecsy,MCMCspecs); %#ok<NASGU>   
    initial.u_mle = u_mle;
end

% restrict surface of initial values incorporating constraint %
% set up mods for constraint %
% modx    = wpspecsx.trackerJacker{wpspecsx.totalDWT}.Kj(1);
% mody    = wpspecsy.trackerJacker{wpspecsy.totalDWT}.Kj(1);

% zero out lower triangle for each node in Y and X %
% for i = 1:size(eW,1)
%     % mod i
%     mi  = mod(i,modx);
%     for j = 1:size(eW,2)
%         % mod j
%         mj  = mod(j,mody);
%         if mj < mi
%             eW(i,j) = 0;
%         end;
%     end;
% end;

% save the initial values
initial.theta_mle   = theta_mle;
initial.betans      = betans;
initial.beta_mle    = betans; % these two are the same.
%==========================================================================
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

%%%% Set initial values for shrinkage hyperparameters pi and tau. 
[tau,PI,alpha,a_tau,b_tau,a_pi,b_pi] = initial_taupi(betans,Vbetans,model,wavespecsy,MCMCspecs,meanop);
initial.PI      = PI;
initial.tau     = tau;
initial.a_tau   = a_tau;
initial.b_tau   = b_tau;
initial.a_pi    = a_pi;
initial.b_pi    = b_pi;
%==========================================================================
PiMat           = PI*expandop; %#zhu#%  expand from p by J to p by K by repeating within same level.
TauMat          = tau*expandop./Vbetans; %#zhu#%  In this code, do we always assume covT=1?

beta            = betans.*alpha.*(TauMat./(TauMat+1));  %%% Posterior expected values of beta.
Wv.L2           = Get_L2(beta,Wv,wavespecsy); %% Update L2 for starting values of beta
fprintf('\n Now have starting values for regularization parameters. \n \n');

%%%% Compute vague proper priors for theta based on empirical Bayes.
[prior_Theta_a,prior_Theta_b] = EmpBayes_Theta(theta_mle,model,MCMCspecs);
fprintf('\n Starting Values Initialized. \n \n'),toc;

p       = model.p;
K       = wavespecsy.K;
B       = MCMCspecs.B;
burnin  = MCMCspecs.burnin;
thin    = MCMCspecs.thin;

if isnan(blocksize)
    blocksize = B;
end

MCMC_alpha          = NaN(blocksize,p*K); 
MCMC_beta           = NaN(blocksize,p*K);
MCMC_theta          = NaN(blocksize,size(theta,1)*size(theta,2));
MCMC_flag_theta     = NaN(blocksize,K); %% indicator function of acceptance of new set of VCs
if (model.H>0)&&(MCMCspecs.sampleU==1) 
    MCMC_U = NaN(blocksize,sum(model.m)*K);
end

MCMC_tau    = NaN(blocksize,p*wavespecsy.J); % in the memory, we only save a block of the samples. 
MCMC_pi     = NaN(blocksize,p*wavespecsy.J);

ii              = 0; % counter for MCMC samples to output
currentblock    = 0;
tic
%%% MCMC Loop: B is desired # samples
for i = 1:(B*thin+burnin)
        
    % (1) Update beta
    [beta,gamma,alpha]  = UpdateBetaNoOrthog_HFLM(beta,Vbetans,PiMat,TauMat,Wv,model,wavespecsy,MCMCspecs); 
    Wv.L2               = Get_L2(beta,Wv,wavespecsy); %#zhu#% beta is updated here, so does L2. 
    
    % (2) Update theta. q_jk, and s_jk.    
    [theta,flag_theta,Vbetans,Wv] = UpdateTheta(beta,theta,Vbetans,Wv,D,W,model,prior_Theta_a,prior_Theta_b,propsd_Theta,wavespecsy,MCMCspecs);   
    
    % (3) Update U when needed.
    if (model.H > 0) && (sampleU == 1)     
       U = UpdateU(beta,theta,model,D,wavespecsy);       
    end
    
    % (4) Update tau_ijk(as well as TauMat) and PI_ij(as well as PiMat).    
    tau     = Update_tau(beta,gamma,a_tau,b_tau,meanop);
    TauMat  = tau*expandop./Vbetans;    
    PI      = Update_pi(gamma,a_pi,b_pi,wavespecsy);
    PiMat   = PI*expandop;
    
    %%%%%%  Record MCMC samples.%%%%%%    
    if   (i > burnin) && (mod(i-burnin,thin) == 0)   %% Save MCMC samples of beta in single matrix.
            ii                      = ii+1; % this is the real row number among the B samples.
            is                      = mod(ii-1,blocksize)+1; % this is the row number in the block.
            MCMC_beta(is,:)         = reshape(beta',1,p*K); %#zhu#% each row=[beta(1,j,k),j,k|beta(2,j,k),j,k|...|beta(p,j,k),j,k]
            MCMC_alpha(is,:)        = reshape(alpha',1,p*K);
            MCMC_theta(is,:)        = reshape(theta,1,size(theta,1)*size(theta,2));
            MCMC_flag_theta(is,:)   = flag_theta; 
            MCMC_tau(is,:)          = reshape(tau,1,p*wavespecsy.J);   % each block is one j, contains p values.
            MCMC_pi(is,:)           = reshape(PI,1,p*wavespecsy.J); % [pi_{11},...pi(1J1);pi_{21},...,pi{2J2};..],J blocks, each block has p values.            
            if (model.H > 0) && (sampleU == 1) 
               MCMC_U(is,:) = reshape(U',1,numel(U));    
            end
            
        %    if any(paramroute)&&(is==blocksize)  % when the block is full, export it to text.
        %       currentblock=currentblock+1;
        %       dlmwrite(strcat(paramroute,'/MCMC_beta_rec.txt'),MCMC_beta,'-append','precision','%12.6e');
        %       dlmwrite(strcat(paramroute,'/MCMC_alpha_rec.txt'),MCMC_alpha,'-append','precision','%12.6e');
        %       dlmwrite(strcat(paramroute,'/MCMC_theta_rec.txt'),MCMC_theta,'-append','precision','%12.6e');
        %       dlmwrite(strcat(paramroute,'/MCMC_flag_theta_rec.txt'),MCMC_flag_theta,'-append','precision','%12.6e');
        %       dlmwrite(strcat(paramroute,'/MCMC_tau_rec.txt'),MCMC_tau,'-append','precision','%12.6e');
        %       dlmwrite(strcat(paramroute,'/MCMC_pi_rec.txt'),MCMC_pi,'-append','precision','%12.6e');
        %       dlmwrite(strcat(paramroute,'/MCMC_U_rec.txt'),MCMC_U,'-append','precision','%12.6e');  
        %       % I think there is some problem with dlmwrite, the when read MCMC_U, it always gives one vector.               
        %    end         
    end
    if mod(i, MCMCspecs.time_update) == 0
       fprintf('\n %d \n',i),toc;
    end
    
end

if saveraw;
    save FPC_results_raw_noshrinkage.mat;
end;

% if any(paramroute)&&(is<blocksize)&&(currentblock==floor(B/blocksize)) % if B is not multiple of blocksize, we still need to save the remaider.
%    dlmwrite(strcat(paramroute,'/MCMC_beta_rec.txt'),MCMC_beta(1:is,:),'-append','precision','%12.6e');
%    dlmwrite(strcat(paramroute,'/MCMC_alpha_rec.txt'),MCMC_alpha(1:is,:),'-append','precision','%12.6e');
%    dlmwrite(strcat(paramroute,'/MCMC_theta_rec.txt'),MCMC_theta(1:is,:),'-append','precision','%12.6e');
%    dlmwrite(strcat(paramroute,'/MCMC_flag_theta_rec.txt'),MCMC_flag_theta(1:is,:),'-append','precision','%12.6e');
%    dlmwrite(strcat(paramroute,'/MCMC_tau_rec.txt'),MCMC_tau(1:is,:),'-append','precision','%12.6e');
%    dlmwrite(strcat(paramroute,'/MCMC_pi_rec.txt'),MCMC_pi(1:is,:),'-append','precision','%12.6e');
%    dlmwrite(strcat(paramroute,'/MCMC_U_rec.txt'),MCMC_U(1:is,:),'-append','precision','%12.6e');                
% end

fprintf('\n Done with MCMC \n');
fprintf('\n  Now projecting results back to data space. \n')

%%%%%  Project MCMC samples back to data space

if ~any(paramroute) % no output of samples
    if sampleU==0        
%         [betahat,beta_025CI,beta_975CI,alphahat,accept_rate_theta,ghat,ghatns,Q025_ghat,Q975_ghat,thetahat,theta_025CI,theta_975CI,tauhat,pihat,...
%         Sigma,]=...
%         PostProcess_inMem(MCMC_beta,MCMC_alpha,MCMC_flag_theta,MCMC_theta,MCMC_tau,MCMC_pi,betans,theta,model,wavespecsy,get_sigma,sampleU);
        postout = PostProcess_FFR(MCMC_beta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wavespecsc,wavespecsy,pcaspecsx,twosample,get_sigma,sampleU);
    else
%         [betahat,beta_025CI,beta_975CI,alphahat,accept_rate_theta,ghat,ghatns,Q025_ghat,Q975_ghat,thetahat,theta_025CI,theta_975CI,tauhat,pihat,...
%         Sigma,uhat,u_025CI,u_975CI,g_Ut,g_Ut025,g_Ut975]=...
%         PostProcess_inMem(MCMC_beta,MCMC_alpha,MCMC_flag_theta,MCMC_theta,MCMC_tau,MCMC_pi,betans,theta,model,wavespecsy,get_sigma,sampleU,MCMC_U);
        postout = PostProcess_FFR(MCMC_beta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wavespecsc,wavespecsy,pcaspecsx,twosample,get_sigma,sampleU,MCMC_U);
    end
    res.MCMC_beta   = MCMC_beta;
    res.MCMC_theta  = MCMC_theta;
    res.MCMC_tau    = MCMC_tau;
    res.MCMC_pi     = MCMC_pi;
    if (model.H > 0) && (sampleU == 1)
         res.MCMC_U = MCMC_U;
    end

else %  with output of samples
        MCMC_beta   = [];
        MCMC_alpha  = [];
        MCMC_theta  = [];
        MCMC_pi     = [];
        MCMC_flag_theta = [];
        MCMC_tau    = [];
        MCMC_pi     = [];
        clear MCMC_beta MCMC_alpha MCMC_theta MCMC_flag_theta MCMC_tau MCMC_pi;
        if (model.H>0)&&(sampleU==1)     
            MCMC_U=[];
            clear  MCMC_U;
        end        
        [betahat,beta_025CI,beta_975CI,alphahat,accept_rate_theta,ghat,ghatns,Q025_ghat,Q975_ghat,thetahat,theta_025CI,theta_975CI,tauhat,pihat,...
        Sigma,uhat,u_025CI,u_975CI,g_Ut,g_Ut025,g_Ut975]=...
        PostProcess(paramroute,betans,theta,model,wavespecsy,B,get_sigma,sampleU);    
end

%%%%% put results into "res" %%%%
res.model       = model;
res.wavespecsy  = wavespecsy;
res.wavespecsx  = wavespecsc;
res.pcaspecsx   = pcaspecsx;
res.MCMCspecs   = MCMCspecs;
res.theta_flag  = theta_flag; % theta_flag=0 indicates treating the theta component as 0 all the time. 
res.betans      = betans; %#zhu#% note that this betans has not been updated during MCMC.
res.Vbetans     = Vbetans; %#zhu#% note that this Vbetans has been updated during MCMC.

res.bstarhat    = postout.bstarhat;
res.bstar025CI  = postout.bstar_025CI;
res.bstar975CI  = postout.bstar_975CI;
res.alphahat    = postout.alphahat;
res.psi         = postout.psi;
res.blocksize   = blocksize; % the block size: how many interations of MCMC samples we save to text each time.
res.paramroute  = paramroute;

if twosample == 1;
    res.Q025_bhat0  = postout.Q025_bhat0;
    res.Q975_bhat0  = postout.Q975_bhat0;
    res.bhat0       = postout.bhat0;
    res.Q025_bhat1  = postout.Q025_bhat1;
    res.Q975_bhat1  = postout.Q975_bhat1;
    res.bhat1       = postout.bhat1;
    res.Q025_diff   = postout.Q025_diff;
    res.Q975_diff   = postout.Q975_diff;
    res.diff        = postout.diff;
    res.Q025_ahat0  = postout.Q025_ahat0;
    res.Q975_ahat0  = postout.Q975_ahat0;
    res.ahat0       = postout.ahat0;
    res.Q025_ahat1  = postout.Q025_ahat1;
    res.Q975_ahat1  = postout.Q975_ahat1;
    res.ahat1       = postout.ahat1;
    res.psi0        = postout.psi0;
    res.psi1        = postout.psi1;
else
    res.bhat        = postout.bhat;
    res.Q025_bhat   = postout.Q025_bhat;
    res.Q975_bhat   = postout.Q975_bhat;
    res.ahat        = postout.ahat;
    res.Q025_ahat   = postout.Q025_ahat;
    res.Q975_ahat   = postout.Q975_ahat;
end;

res.thetahat    = postout.thetahat;
res.theta_025CI = postout.theta_025CI;
res.theta_975CI = postout.theta_975CI;
res.acpt_rate   = postout.accept_rate_theta;
res.theta_MLE   = theta_mle;
res.theta_MLEsd = theta_sd;

res.tauhat      = postout.tauhat;
res.pihat       = postout.pihat;

if (model.H > 0) && (sampleU == 1)
    % res.MCMC_U=MCMC_U;
    res.Uhat       = postout.uhat;
    res.U_025CI    = postout.u_025CI;
    res.U_975CI    = postout.u_975CI;
    res.g_Ut       = postout.g_Ut;
    res.g_Ut025    = postout.g_Ut025;
    res.g_Ut975    = postout.g_Ut975;
end
if get_sigma == 1
    res.Sigma = postout.Sigma;
end
%==========================================================================
%res.Y=Y; % we actually don't have to save Y.
%==========================================================================
res.Wv          = Wv; %#zhu#% this is the Wv of very last iteration.
res.D           = D;
res.initial     = initial;
