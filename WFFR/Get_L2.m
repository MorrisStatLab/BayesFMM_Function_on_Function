function L2=Get_L2(beta,Wv,wavespecs)

%%%%%   [L2]=Get_L2(model,wlevels,XvX,Xvd,dvd,beta)
%%%%%   Compute exponential part of normal likelihood (conditioning on
%%%%%   beta and theta)
%%%%%       L2(beta_jk,theta_jk) =  (d_jk-X beta_jk)'Sigma(theta_jk)^(-1)(d_jk-X beta_jk)
%%%%%   for (j,k)= (1,1), ..., (J,K_J)
%%%%%   
%%%%%   
%%%%%
%
% Input:    
%           wlevels = (J x 1) vector -- specifies how many wavelet coefficients
%                                   per wavelet level.
%           model = structure with integer elements:
%               X = (n x p) design matrix for fixed effects functions;
%               Z = cell array containg H matrices (each n x m_h); design
%                       matrices for random effects functions; 
%               covQ = assumptions on within-curve random effect variance components q
%               covS = assumptions on within-curve residual variance components s
%                           (0=homoscesastic, 1=homoscedastic within
%                           levels, 2=heteroscedastic)
%           beta = (p x K) matrix of coefficients.
%           XvX = cell array with K elements, each containing X' Sigma_{jk}^(-1) X
%           Xvd, dvd:  similarly defined

K=wavespecs.K;
L2=NaN(1,K);

XvX=Wv.XvX;
Xvd=Wv.Xvd;
dvd=Wv.dvd;

for j=1:K
      XvXk=XvX(:,:,j);
      L2(j)=dvd(:,j)+beta(:,j)'*XvXk*beta(:,j)-2*beta(:,j)'*Xvd(:,j);
end

        