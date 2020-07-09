function [beta,gamma,alpha]=UpdateBetaNoOrthog(beta,Vbeta,PI,T,Wv,model,wavespecs,MCMCspecs)

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


p=model.p;
K=wavespecs.K;

XvX=Wv.XvX;
Xvd=Wv.Xvd;

Btilde=Vbeta.*Xvd;
gamma=zeros(p,K);
alpha=zeros(p,K);
Betansi=zeros(p,K);
for i=1:p  %#zhu#% compute beta_i conditional on beta_{-i} here, it is the correct as we derived, for each j,k, beta_{ijk} updated alternatively.
    Bi=zeros(1,K);
    for k=1:K
        Bi(k)=XvX(i,:,k)*beta(:,k); %#zhu#% these are the stuff to conditional on.
    end
    
    Betansi(i,:)=Btilde(i,:)+beta(i,:)-Vbeta(i,:).*Bi;
    [beta(i,:),gamma(i,:),alpha(i,:)]=UpdateBeta(Betansi(i,:),Vbeta(i,:),PI(i,:),T(i,:),MCMCspecs);%#zhu#% I suspect a problem in this code.
end % i loop

%[Btilde(:,176),Betansi(:,176),beta(:,176)]

%beta=reshape(MCMC_beta(1,:),p,306);
%beta(:,176)


