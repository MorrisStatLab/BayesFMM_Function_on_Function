function [theta,flag,Vbetans,Wv]=UpdateTheta(beta,theta,Vbetans,Wv,D,W,model,a,b,propsd_Theta,wavespecs,MCMCspecs)
%%%%    UpdateTheta(model,beta,theta,L1,L2,W,a,b,propsd_theta): Updates Variance components conditional on fixed effects beta
%%%%         using Metropolis step for each column of theta in Wavelet-based Functional Mixed Model.
%%%%
%%%%    Input:  model = structure containing model information.           
%                   X = (n x p) design matrix for fixed effects functions;
%                   Z = cell array of H matrices, each containing (n x m_h) design
%                       matrix for set of random effects functions; 
%               beta: p x K matrix containing values of fixed effects from
%                       current MCMC iteration
%               theta: H+c x * matrix containing current values of variance
%                       components
%                           (*=1 if model.covQ=0, J if model.covQ=1, K if model.covQ=2)
%               betans: p x K matrix of non-shrinkage estimates of fixed
%                           effect coefficients
%               Vbetans: p x K matrix containing estimated variances of each of the
%                       non-shrinkage fixed effects estimates
%               L1: 1 x * matrix containing values of log|Sigma_jk|
%               L2: 1 x * matrix containing values of 
%                   (d_jk-X beta_jk)'Sigma(Theta_jk)^{-1} (d_jk-X beta_jk)
%                   (summed over k if covQ=1)
%
%               W = structure containing:
%                   XtX=X'X
%                   XtZ=X'Z
%                   XtD=X'D
%                   ZtD=Z'D
%               a = H+1 x * matrix containing first inverse gamma parameter
%                   for prior on theta
%               b = H+1 x * matrix containing second inverse gamma
%                   parameter for prior on theta
%               propsd_theta = H+c x * matrix of proposal standard
%                   deviations for random walk Metropolis for theta
%       Output:
%               theta = new value of theta
%               newtheta = 1 x * matrix containing indicator if "new"
%                           theta accepted by Metropolis
%               XvX, Xvd, dvd, XvZ, ZvZ, Zvd = updated values of
%                           X'(V^(-1))X, etc
%               betans, Vbetans = updated values of betans, Vbetans
%               L1, L2: new values of L1 and L2
%
%%%%    Functions needed: 
%%%%      GetGCP.m, Get_L2.m, rep.m                  

p=model.p;
K=wavespecs.K;
theta_old=theta;
Vbetans_old=Vbetans;
Wv_old=Wv;

%minVar=MCMCspecs.minVC;
%theta_new=max(normrnd(theta,propsd_Theta),minVar);%#zhu#% problem here.
%we should use log proposal. 
theta_new=max(exp(normrnd(log(theta_old),propsd_Theta)),MCMCspecs.minVC);

[betans_new,Vbetans_new,Wv_new]=GetGCP_byblocks(theta_new,D,W,model,wavespecs,1,MCMCspecs);

Wv_new.L2=Get_L2(beta,Wv_new,wavespecs);

A1=.5*(Wv_old.L1-Wv_new.L1+Wv_old.L2-Wv_new.L2); % log likelihood ratio
A2=sum((a+1).*(log(theta_old)-log(theta_new))+b.*(theta_old.^(-1)-theta_new.^(-1))); % log prior ratio.
Prop_ratio=sum(log(theta_new)-log(theta_old)); % for each j, f(sig_new|sig_old)/f(sig_old|sig_new)             

log_ratio=A1+A2+Prop_ratio;
if any(isnan(log_ratio)==1)||any(isinf(log_ratio)) 
    error('Update Theta M-H logratio contains Inf or NaN');
end
flag=log(rand(1,K))<log_ratio; % with probability ratio, accept the new theta.
flag1=repmat(flag,size(theta,1),1);
flag2=repmat(flag,p,1);

theta=flag1.*theta_new+(1-flag1).*theta_old; %#zhu#% Theta is updated
%betans=flag2.*betans_new+(1-flag2).*betans_old;  %#zhu#% betans will not
%be used in the future MCMCs hence need not to be updated.
Vbetans=flag2.*Vbetans_new+(1-flag2).*Vbetans_old;

Wv.L1=flag.*Wv_new.L1+(1-flag).*Wv_old.L1; %#zhu#% Start to update Wv based on new theta.
Wv.L2=flag.*Wv_new.L2+(1-flag).*Wv_old.L2; 

for j=1:K  % update Wv.   
    Wv.XvX(:,:,j)=flag(j)*Wv_new.XvX(:,:,j)+(1-flag(j))*Wv_old.XvX(:,:,j);
    if model.H~=0
        Wv.XvZ(:,:,j)=flag(j)*Wv_new.XvZ(:,:,j)+(1-flag(j))*Wv_old.XvZ(:,:,j);
        Wv.ZvZ(:,:,j)=flag(j)*Wv_new.ZvZ(:,:,j)+(1-flag(j))*Wv_old.ZvZ(:,:,j);       
        Wv.Zvd(:,j)=flag(j)*Wv_new.Zvd(:,j)+(1-flag(j))*Wv_old.Zvd(:,j);
    end
    Wv.Xvd(:,j)=flag(j)*Wv_new.Xvd(:,j)+(1-flag(j))*Wv_old.Xvd(:,j);    
    Wv.dvd(j)=flag(j)*Wv_new.dvd(j)+(1-flag(j))*Wv_old.dvd(j);        
end
