function [ initial ] = WavInitVals( Y, wavespecsy, model, MCMCspecs )
%WAVINITVALS Find inital values for coefficient sampler in OPWAVFM
%
%   Created: 3/7/2014
%   By: Mark John Meyer

%% set model initial values (put in a function)
[D,wavespecsy]  = dwt_rows(Y,wavespecsy);

wavespecsy.K        = sum(wavespecsy.Kj);
meanop              = uneqkron(wavespecsy.Kj)*diag(1./wavespecsy.Kj);
expandop            = uneqkron(wavespecsy.Kj)';
initial.expandop    = expandop;
initial.meanop      = meanop;
initial.wavespecsy  = wavespecsy;
initial.D           = D;

tic

%% find ML Estimates
% model.H     = 0;
W           = GetW(model,D);
[theta_mle,theta_sd,betans,Vbetans,Wv,] = regression_mle(W,model,D,wavespecsy,MCMCspecs);    

%% save the initial values
initial.theta_mle   = theta_mle;
initial.Vbetans     = Vbetans;
initial.beta_mle    = betans; % these two are the same.
initial.W           = W;

%% Check whether X matrix is orthogonal
%  If it is, then some calculation shortcuts are available later
if (sum(sum(abs(diag(diag(W.XtX))-W.XtX)>1e-15))>0)
    orthogonal = 0;
else
    orthogonal = 1;
end
model.orthogonal_X  = orthogonal;

theta               = theta_mle;
theta_flag          = (theta_mle>=MCMCspecs.VC0_thresh); % If less than VC0_thresh, treat it as zero and don't update in MH.(stepsize=0)
propsd_Theta        = MH_stepsize(theta_sd,theta_mle,MCMCspecs.propsdTheta,0).*theta_flag; % if theta_flag=1, then step size=0.

%% save initial values
initial.theta           = theta;
initial.theta_flag      = theta_flag;
initial.propsd_Theta    = propsd_Theta;

%% Set initial values for shrinkage hyperparameters pi and tau. 
[tau,PI,alpha,a_tau,b_tau,a_pi,b_pi] = initial_taupi(betans,Vbetans,model,wavespecsy,MCMCspecs,meanop);
initial.PI      = PI;
initial.tau     = tau;
initial.a_tau   = a_tau;
initial.b_tau   = b_tau;
initial.a_pi    = a_pi;
initial.b_pi    = b_pi;
initial.PiMat   = PI*expandop; %#zhu#%  expand from p by J to p by K by repeating within same level.

%% get TauMat
TauMat          = tau*expandop./Vbetans; %#zhu#%  In this code, do we always assume covT=1?
initial.TauMat  = TauMat;

%%
bstar           = betans.*alpha.*(TauMat./(TauMat+1));  %%% Posterior expected values of beta.
Wv.L2           = Get_L2(bstar,Wv,wavespecsy); %% Update L2 for starting values of beta
initial.bstar   = bstar;
initial.Wv      = Wv;
fprintf('\n Now have starting values for regularization parameters. \n \n');

%% Compute vague proper priors for theta based on empirical Bayes.
[prior_Theta_a,prior_Theta_b]   = EmpBayes_Theta(theta_mle,model,MCMCspecs);
initial.prior_Theta_a           = prior_Theta_a;
initial.prior_Theta_b           = prior_Theta_b;
fprintf('\n Starting Values Initialized. \n \n'),toc;

end

