function [MCMC_contrast,cbeta,LB,UB,se]=get_contrast(C,MCMC_beta,p,alpha)

%%% [mean,LB,UB]=get_contrast(C,MCMC_beta,alpha) : Compute contrast and 1-alpha credible intervals on specified contrast
%%%
%%% Inputs: MCMC_beta:  B x p*T matrix of samples of beta from MCMC
%%%         C:          c x p matrix, each row specifying c contrasts
%%%         p:          number of fixed effects functions
%%%         alpha:      specify pointwise one-sided type I error (0.05 by default)
%%%
%%% Outputs:MCMC_contrast B x c*T matrix of samples of C*beta from MCMC
%%%         cbeta:        c x T matrix of estimates of contrast
%%%         LB:           c x T matrix of pointwise alpha quantiles
%%%         UB:           c x T matrix of pointwise 1-alpha quantiles
%%%

if (nargin<4)
    alpha=0.05;
end;

B=size(MCMC_beta,1);
T=size(MCMC_beta,2)/p;
c=size(C,1);
MCMC_contrast=repmat(0,B,c*T);
for (i=1:B)
    MCMC_contrast(i,:)=reshape((C*reshape(MCMC_beta(i,:),T,p)')',1,c*T);
end;
MCMC_contrast=sort(MCMC_contrast);
LB=reshape(MCMC_contrast(floor(alpha*size(MCMC_contrast,1)),:)',T,c)';
UB=reshape(MCMC_contrast(floor((1-alpha)*size(MCMC_contrast,1)),:)',T,c)';
cbeta=reshape(mean(MCMC_contrast)',T,c)';
se=reshape(sqrt(var(MCMC_contrast))',T,c)';
