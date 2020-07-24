function res=wavfmm4_v3_zhu(Y,model,wavespecs,MCMCspecs,sampleU) 
%function res=wavfmm4_v3_zhu(D,model,wavespecs,MCMCspecs,sampleU) 
%==========================================================================
% the above function will be reused later on, when dwt is used.
%==========================================================================
%function res=wavfmm4_zhu(True,model,wavespecs,MCMCspecs) 
% Wavefmm4_v3_zhu: Based on wavfmm4_zhu, add IG prior to \tau_{ij}(in stead of tau_{ijk} in v2), and add prior to 
%                  pi_{ij}.
%  
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
%        wavespecs=structure containing wavelet settings, including:
%                wavelet = wavelet basis to use
%                nlevels = number of levels of decomposition to do
%        MCMCspecs=parameters for MCMC iterations.  
%                minVC = minimum size of VC (tolerance limit) (in EmpBayes) default 0.001
%                epsilon = threshhold on MOM est. of VC below which VC is assumed to be 0 (in EmpBayes) default 0.001
%                delta_theta = multiple for prior on theta: "number of datasets of information" in prior (in
%                             wavefmm4) default 0.0001
%                xi_theta = multiple of MOM est of VC to use as proposal variance
%                                   (in wavefmm4) -- default 1.5
%       sampleU: if 1, sample random effect coefficient U. If 0, do not
%       sample U.
%
% Output:   res= data structure with elements: 
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
%%% Step 1: Apply DWT to project observed functions into wavelet space
[D,wavespecs]=dwt_rows(Y,wavespecs);   
fprintf('\n Finished projecting data to wavelet space.\n \n');
%#zhu#% the above two lines will be uncommented later on when we introduce
%dwt and idwt to the function. The following line is temporary for
%simulation purpose. 
%D=True.D;
meanop=uneqkron(wavespecs.Kj)*diag(1./wavespecs.Kj);
expandop=uneqkron(wavespecs.Kj)';
%==========================================================================
[model.n,model.p]=size(model.X);   %%% n=# functions, p=# covariates
model.c=size(model.C,2);               %%% assumed to be iid mean zero with a
if model.Z{1}==0
    model.H=0; 
    W=GetW(model,D);   
    [theta_mle,theta_sd,betans,Vbetans,Wv,]=regression_mle(W,model,D,wavespecs,MCMCspecs);    
else
    model.H=length(model.Z);           %%% H=# of sets of random effects    
    model.m=repmat(0,1,model.H);       %%% # of random effects per set %#zhu#% When z{1}=0, model.m=[];
    for h=1:model.H
        model.m(h)=size(model.Z{h},2);  %#zhu#% When z{1}=0, this loop will not be excuted since model.H=0
    end
    %%% Compute the initial value of theta using MLE and standard
    %%% deviation.
    model.M=sum(model.m);
    W=GetW(model,D);
    theta0=chi2rnd(1,model.H+model.c,wavespecs.K);
    [theta_mle,theta_var,theta_sd,beta_mle,betans,Vbetans,u_mle,Wv,converge]=mixed_mle(W,theta0,model,D,wavespecs,MCMCspecs); %#ok<NASGU>   
    initial.u_mle=u_mle;
end
% save the initial values
initial.theta_mle=theta_mle;
initial.betans=betans;
initial.beta_mle=betans; % these two are the same.
%==========================================================================
%%%% Check whether X matrix is orthogonal
%%%% If it is, then some calculation shortcuts are available later
if (sum(sum(abs(diag(diag(W.XtX))-W.XtX)>1e-15))>0)
    orthogonal=0;
else
    orthogonal=1;
end
model.orthogonal_X=orthogonal;

theta=theta_mle;
theta_flag=theta_mle>=MCMCspecs.VC0_thresh; % If less than VC0_thresh, treat it as zero and don't update in MH.(stepsize=0)
propsd_Theta=MH_stepsize(theta_sd,theta_mle,MCMCspecs.propsdTheta,0).*theta_flag; % if theta_flag=1, then step size=0.

%%%% Set initial values for shrinkage hyperparameters pi and tau. 
%[tau,PI,alpha,a_tau,b_tau,a_pi,b_pi]=initial_taupi(betans,Vbetans,model,wavespecs,MCMCspecs);
[tau,PI,alpha,a_tau,b_tau,a_pi,b_pi]=initial_taupi(betans,Vbetans,model,wavespecs,MCMCspecs,meanop);
initial.PI=PI;
initial.tau=tau;
initial.a_tau=a_tau;
initial.b_tau=b_tau;
initial.a_pi=a_pi;
initial.b_pi=b_pi;
%==========================================================================
PiMat=PI*expandop; %#zhu#%  expand from p by J to p by K by repeating within same level.
TauMat=tau*expandop./Vbetans; %#zhu#%  In this code, do we always assume covT=1?

beta=betans.*alpha.*(TauMat./(TauMat+1));  %%% Posterior expected values of beta.
Wv.L2=Get_L2(beta,Wv,wavespecs); %% Update L2 for starting values of beta
fprintf('\n Now have starting values for regularization parameters. \n \n');

%%%% Compute vague proper priors for theta based on empirical Bayes.
[prior_Theta_a,prior_Theta_b]=EmpBayes_Theta(theta_mle,model,MCMCspecs);
fprintf('\n Starting Values Initialized. \n \n');

p=model.p;
K=wavespecs.K;
B=MCMCspecs.B;
burnin=MCMCspecs.burnin;
thin=MCMCspecs.thin;

MCMC_alpha=NaN(B,p*K); 
MCMC_beta=NaN(B,p*K);
MCMC_theta=NaN(B,size(theta,1)*size(theta,2));
MCMC_flag_theta=NaN(B,K); %% indicator function of acceptance of new set of VCs
if (model.H>0)&&(sampleU==1) 
    MCMC_U=NaN(B,sum(model.m)*K);
end

MCMC_tau=NaN(B,p*wavespecs.J);
MCMC_pi=NaN(B,p*wavespecs.J);

ii=0; % counter for MCMC samples to output
%%% MCMC Loop: B is desired # samples
for i=1:(B*thin+burnin)
        
    % (1) Update beta
    [beta,gamma,alpha]=UpdateBetaNoOrthog(beta,Vbetans,PiMat,TauMat,Wv,model,wavespecs,MCMCspecs); 
    Wv.L2=Get_L2(beta,Wv,wavespecs); %#zhu#% beta is updated here, so does L2. 
    
    % (2) Update theta. q_jk, and s_jk.    
    [theta,flag_theta,Vbetans,Wv]=UpdateTheta(beta,theta,Vbetans,Wv,D,W,model,prior_Theta_a,prior_Theta_b,propsd_Theta,wavespecs,MCMCspecs);   
    %theta./theta_start
    
    % (3) Update U when needed.
    if (model.H>0)&&(sampleU==1)     
       U=UpdateU(beta,theta,model,D,wavespecs);       
    end
    
    % (4) Update tau_ijk(as well as TauMat) and PI_ij(as well as PiMat).    
    %tau=Update_tau(tau,beta,gamma,a_tau,b_tau,MCMCspecs,wavespecs);
    tau=Update_tau(beta,gamma,a_tau,b_tau,meanop);
    TauMat=tau*expandop./Vbetans;    
    PI=Update_pi(gamma,a_pi,b_pi,wavespecs);
    PiMat=PI*expandop;
    
    %%%%%%  Record MCMC samples.       
    if (i>burnin)&&(mod(i-burnin,thin)==0)      %% Save MCMC samples of beta in single matrix -- 1 column per iteration
            ii=ii+1;
            MCMC_beta(ii,:)=reshape(beta',1,p*K); %#zhu#% each row=[beta(1,j,k),j,k|beta(2,j,k),j,k|...|beta(p,j,k),j,k]
            MCMC_alpha(ii,:)=reshape(alpha',1,p*K);
            MCMC_theta(ii,:)=reshape(theta,1,size(theta,1)*size(theta,2));
            MCMC_flag_theta(ii,:)=flag_theta; 
            MCMC_tau(ii,:)=reshape(tau,1,p*wavespecs.J);   % each block is one j, contains p values.
            MCMC_pi(ii,:)=reshape(PI,1,p*wavespecs.J); % [pi_{11},...pi(1J1);pi_{21},...,pi{2J2};..],J blocks, each block has p values.            
            if (model.H>0)&&(sampleU==1) 
            MCMC_U(ii,:)=reshape(U',1,numel(U));    
            end
            
    end
    if mod(i, MCMCspecs.time_update)==0
            fprintf('\n %d \n',i);          
    end
    
end

fprintf('\n Done with MCMC \n');
fprintf('\n  Now projecting results back to data space. \n')
%%%%%  Project MCMC samples back to data space
get_sigma=0;
[betahat,beta_025CI,beta_975CI,alphahat,accept_rate_theta,ghat,ghatns,Q025_ghat,Q975_ghat,MCMC_g,thetahat,theta_025CI,theta_975CI,tauhat,pihat,Sigmat,Sigma,Rho]=...
    PostProcess(MCMC_beta,MCMC_alpha,MCMC_flag_theta,MCMC_theta,MCMC_tau,MCMC_pi,betans,theta,model,wavespecs,MCMCspecs,get_sigma);
%%%%% put results into "res"
res.model=model;
res.wavespecs=wavespecs;
res.MCMCspecs=MCMCspecs;
res.theta_flag=theta_flag; % theta_flag=0 indicates treating the theta component as 0 all the time. 
res.betans=betans; %#zhu#% not that this betans has not been updated during MCMC.
res.Vbetans=Vbetans; %#zhu#% note that this Vbetans has been updated during MCMC.

res.betahat=betahat;
res.beta025CI=beta_025CI;
res.beta975CI=beta_975CI;
res.alphahat=alphahat;
res.MCMC_beta=MCMC_beta;

res.MCMC_g=MCMC_g; %output MCMC samples of data-space betas
res.ghatns=ghatns;
res.ghat=ghat;
res.Q025_ghat=Q025_ghat;
res.Q975_ghat=Q975_ghat;

res.thetahat=thetahat;
res.MCMC_theta=MCMC_theta;
res.theta_025CI=theta_025CI;
res.theta_975CI=theta_975CI;
res.acpt_rate=accept_rate_theta;
res.theta_MLE=theta_mle;
res.theta_MLEsd=theta_sd;

res.MCMC_tau=MCMC_tau;
res.MCMC_pi=MCMC_pi;
res.tauhat=tauhat;
res.pihat=pihat;

if (model.H>0)&&(sampleU==1)
     res.MCMC_U=MCMC_U;
end
res.Rho=Rho;
res.Sigma=Sigma;
res.Wv=Wv; %#zhu#% this is the Wv of very last iteration.
%==========================================================================
res.Y=Y; %#zhu#% this line will be uncommented later on when we introduce
%dwt and idwt to the function
%==========================================================================
res.D=D;
res.initial=initial;
