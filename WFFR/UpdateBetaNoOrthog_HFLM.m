function [ beta, gamma, alpha ] = UpdateBetaNoOrthog_HFLM( beta, Vbeta, PI, T, Wv, wpspecs, MCMCspecs )

%%%%    UpdateBetaNoOrthog(Vbetans,PI,T,XvX,Xvd,model,wlevels): Updates Fixed effects Beta from f(Beta|theta,pi,T,D) 
%%%%    for a Wavelet-based Functional Mixed Model.
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
%%%%            wlevels: contains number of wavelet coefficients per level
%%%%
%%%%    Output: (p x K) matrix beta, containing samples values of fixed
%%%%                                    effects
%%%%            (p x K) matrix gamma, containing indicators for whether
%%%%                                    coefficient is nonzero
%%%%            (p x K) matrix alpha, containing the posterior probability
%%%%                                of nonzero wavelet coefficients
%%%%
%%%%    Functions needed: UpdateBeta(betans,Vbetans,PI,T);


% Vbeta = Vbetans;
% PI    = PiMat;
% T     = TauMat;

% p=model.p;
% K=wavespecs.K;

p       = size(beta,1);
K       = size(beta,2);
pv      = 1:p;

%% set up mods %%
modx            = wpspecs.wpspecsx.trackerJacker{wpspecs.wpspecsx.totalDWT}.Kj(1);
mody            = wpspecs.wpspecsy.trackerJacker{wpspecs.wpspecsy.totalDWT}.Kj(1);
pvm             = mod(pv,modx);
pvm(pvm == 0)   = modx;

%%
XvX     = Wv.XvX;
Xvd     = Wv.Xvd;

%% update beta
Btilde  = Vbeta.*Xvd;
gamma   = zeros(p,K);
alpha   = zeros(p,K);
Betansi = zeros(p,K);

%%
for i = 1:p  %#zhu#% compute beta_i conditional on beta_{-i} here, it is the correct as we derived, for each j,k, beta_{ijk} updated alternatively.
%     Bi=zeros(1,K);

%%
    mi  = mod(i,modx);
    if mi == 0; mi = modx; end;
%%
    for k = 1:K
        mk  = mod(k,mody);
        if mk == 0; mk = mody; end;
%         Bi(k)=XvX(i,:,k)*beta(:,k); %#zhu#% these are the stuff to conditional on.
        % NEED TO MOD Bi AND beta %
%%
        if mk >= mi
            %% find index ends for constraint %%
            from    = pv(pvm == 1);
            to      = pv(pvm == mk);
            %% build index vector %%
            ind = zeros(to(1),length(from));
            for r = 1:length(from)
                ind(:,r) = from(r):to(r);
            end
            ind     = reshape(ind,1,to(1)*length(from));
            %% implement Jeff's sampler sampling only hist coefs %%
            Bi              = XvX(i,ind,k)*beta(ind,k); % MODIFY :, must only depend on desired elements of beta %
            Betansi(i,k)    = Btilde(i,k) + beta(i,k) - Vbeta(i,k).*Bi;
            [ beta(i,k), gamma(i,k), alpha(i,k) ]   = UpdateBeta( Betansi(i,k), Vbeta(i,k), PI(i,k), T(i,k), MCMCspecs );%#zhu#% I suspect a problem in this code.
        end
    end
    % Can we just pull this into the loop above? %
%     Betansi(i,:)=Btilde(i,:)+beta(i,:)-Vbeta(i,:).*Bi;
%     [beta(i,:),gamma(i,:),alpha(i,:)]=UpdateBeta(Betansi(i,:),Vbeta(i,:),PI(i,:),T(i,:),MCMCspecs);%#zhu#% I suspect a problem in this code.
end % i loop

%[Btilde(:,176),Betansi(:,176),beta(:,176)]

%beta=reshape(MCMC_beta(1,:),p,306);
%beta(:,176)



%% NEED MATRIX OF ZEROS AND ONES THAT CORRESPOND TO WHICH COEFFICIENTS NEED TO BE SAMPLED %%
%% CAN ALSO STREAMLINE SAMPLER BY INCORPORATING INTERCEPT AND SCALAR COVARIATES HERE TOO %%


