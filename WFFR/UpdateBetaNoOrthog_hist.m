function [beta,gamma,alpha]=UpdateBetaNoOrthog_hist(beta, Vbeta, PI, T, Wv, model, wpspecs, MCMCspecs)

%%%%    UpdateBetaNoOrthog_hist(Vbetans,PI,T,XvX,Xvd,model,wlevels): Updates Fixed effects Beta from f(Beta|theta,pi,T,D) 
%%%%    for a Historical Wavelet-based Functional Mixed Model.
%%%%        * Updates Betas one-at-a-time
%%%%        * Assumes non-orthogonal design matrix X
%%%%
%%%%    Input:  (p x K) matrix beta: matrix of sampled betas from previous
%%%%                                MCMC iteration.
%%%%            (p x K) matrix Vbeta: matrix of variances of betans
%%%%            (p x K) matrix pi: prior probabilities of nonzero
%%%%                                coefficients
%%%%            (p x K) matrix T: prior variances for nonzero coefficients
%%%%            cell array of J or K (p x p) matrices XvX: X'(Sig_jk)^(-1)X
%%%%            cell array of J or K (p x 1) vectors Xvd: X'(Sig_jk)d_jk
%%%%            model: contains model information
%%%%            wpspecs: contains details of wavelet-packet decomposition
%%%%            sampCoef: indicator matrix of which coefficients to sample
%%%%                from beta for both functional (and thereby constrained)
%%%%                and scalar covariates
%%%%
%%%%    Output: (p x K) matrix beta, containing samples values of fixed
%%%%                                    effects (functional and scalar)
%%%%            (p x K) matrix gamma, containing indicators for whether
%%%%                                    coefficient is nonzero
%%%%            (p x K) matrix alpha, containing the posterior probability
%%%%                                of nonzero wavelet coefficients
%%%%
%%%%    Functions needed: UpdateBeta(betans,Vbetans,PI,T,MCMCspecs);

%% function parameters %%
sampCoef    = model.sampCoef;
p           = model.p;
K           = wpspecs.K;

%% 
XvX         = Wv.XvX;
Xvd         = Wv.Xvd;

%% declare matrices %%
Btilde      = Vbeta.*Xvd;
gamma       = zeros(p,K);
alpha       = zeros(p,K);
Betansi     = zeros(p,K);

%% sample desired coefficents by row %%
for i=1:p  %#zhu#% compute beta_i conditional on beta_{-i} here, it is the correct as we derived, for each j,k, beta_{ijk} updated alternatively.
    ind     = sampCoef(i,:);
    indi    = find(ind == 1);
    Si      = sum(ind);
    Bi      = zeros(1,Si);

    %% build Bi using only desired coefs %%
    for j = 1:size(indi,2)
        k       = indi(j);
        Bi(j)   = XvX(i, :, k)*beta(:, k); %#zhu#% these are the stuff to conditional on.
    end
    %%
    Betansi(i,indi)     = Btilde(i,indi) + beta(i,indi) - Vbeta(i,indi).*Bi;

    %% update beta, gamma, and alpha (alpha is P(gamma == 1))
    [beta(i,indi), gamma(i,indi), alpha(i,indi)]   = UpdateBeta(Betansi(i,indi), Vbeta(i,indi), PI(i,indi), T(i,indi), MCMCspecs);
end
