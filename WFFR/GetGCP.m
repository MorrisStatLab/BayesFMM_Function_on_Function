function [betans,Vbetans,Wv]=GetGCP(theta,W,model,wavespecs,Update_betans,D,MCMCspecs)
%function [XvX,Xvd,dvd,XvZ,ZvZ,Zvd,betans,Vbetans,L1]=GetGCP(model,theta,wlevels,W)
%%%%% [XvX,Xvd,dvd,XvZ,ZvZ,Zvd,L1]=GetGCP(model,theta,wlevels,W)
%%%%%   Compute generalized cross-products matrices
%%%%%           X'Siginv_jk X       X'Siginv_jk Z       X'Siginv_jk d_jk
%%%%%           Z'Siginv_jk X       Z'Siginv_jk Z       Z'Signinv_jk d_jk
%%%%%           d_jk'Siginv_jk X    d_jk'Siginv_jk Z    d_jk'Siginv_jk d_jk
%%%%%   as well as L1=log|Sig_jk|
%%%%%   for (j,k)= (1,1), ..., (J,K_J)
%%%%%   
%%%%%   using Sweep method suggested in Wolfinger, Tobias, and Sall, 1994
%%%%%
%
% Input:    
%           wlevels = (J x 1) vector -- specifies how many wavelet coefficients
%                                   per wavelet level.
%           model = structure with integer elements:
%               X = (n x p) design matrix for fixed effects functions;
%               Z = cell array containg H matrices (each n x m_h); design
%                       matrices for random effects functions; 
%           theta = (H+c x K) matrix of within-curve variance component
%                           estimates.
%                           
%           W = structure containing:
%
%                 XtX=X'X
%                 XtZ=X'Z
%                 XtD=X'D
%                 ZtD=Z'D
%
%           Wv = structure containing:
% Output:   XvX = cell array with K elements, each containing X' Sigma_{jk}^(-1) X
%           Xvd, dvd, XvZ, ZvZ, Zvd: similarly defined (except make
%               Xvd,Zvd,dvd be matrices of size p x K, m x K, and 1 x K,
%               respectively)
%           L1 = log|Sigma|
%
%  
%
%   Functions needed: repvec.m -- same as Splus function "rep"
%                     sweep.m -- perform sweep operation
%
%
%   Takes 0.09 seconds on test data with covQ=covS=1 (0.14 if write out
%                   redundant XtX for all k)
%         1.04 seconds on tests data with covQ=covS=2
%
%   Changed 1/13/05 to accommodate no random effects (Z{1}=0), and
%   to allow residual covariance matrix S to vary by some covariate
%   C, with c discrete levels (only changed covQ=covS=2 case)
%   Note: %#zhu#%: In the current version, I removed covQ=covS=0 and 1 cases. 
%   And XvX, XvZ, ZvZ are put in 3-d arrays(rather than cell arrays).
%   Reference: (Sweep operator) Wolfinger, R. 1994, Computing Gaussian
%   Liklihoods and their derivatives for general linear mixed models.
epsilon=MCMCspecs.VC0_thresh; %%% threshhold below which VC is assumed 0.
K=wavespecs.K;
n=model.n;
p=model.p;
c=model.c;
C=model.C;
Xrange=1:p;

if model.H~=0
    m=model.m;
    M=sum(m);
    theta_q=theta(1:model.H,:);
    theta_s=theta(model.H+1:end,:);    
    
    
    XvZ=NaN(p,M,K);
    ZvZ=NaN(M,M,K);
    Zvd=NaN(M,K);
    Zrange=(Xrange(end)+1):(Xrange(end)+M);
    Drange=(Zrange(end)+1):(Zrange(end)+1); %#zhu#% This is because we compute Wv for each D_{j,k}, so the no. of index is one.
else 
    theta_s=theta(model.H+1:end,:);  
    Drange=(Xrange(end)+1):(Xrange(end)+1);
end

L1=NaN(1,K);
XvX=NaN(p,p,K);
Xvd=NaN(p,K);
dvd=NaN(1,K);

temp_W=cell(c,1);
for l=1:c        %#zhu#% When c!=1, extract that group of data and compute them separately, this may help remove unnecessary columns of X or Z, save time.
    temp_model=model;
    temp_model.X=model.X(C(:,l)==1,:);     
    if model.H~=0
        for h=1:model.H
            temp_model.Z{h}=model.Z{h}(C(:,l)==1,:);
        end
    end
    temp_D=D(C(:,l)==1,:);
    temp_W{l}=GetW(temp_model,temp_D);
end

 for j=1:K        
      
      temp_indicator=0;
      if (model.H>0)&&(max(theta_q(:,j))>epsilon)
         temp_indicator=1;
      end
      
      if temp_indicator==1
            if c==1 
                W0=[W.XtX,W.XtZ,W.XtD(:,j);W.XtZ',W.ZtZ,W.ZtD(:,j);W.XtD(:,j)',W.ZtD(:,j)',W.DtD(j,j)]/theta_s(j);%#zhu#% DtD(j,j) is that for one column of D.
                %temp1=eye(M)+diag(repvec(theta_q(:,j),m))*(W.ZtZ/(theta_s(j)));           
            else
                W0=0;
                for l=1:c
                   W0=W0+[temp_W{l}.XtX,temp_W{l}.XtZ,temp_W{l}.XtD(:,j);temp_W{l}.XtZ',temp_W{l}.ZtZ,temp_W{l}.ZtD(:,j);temp_W{l}.XtD(:,j)',temp_W{l}.ZtD(:,j)',temp_W{l}.DtD(j,j)]/theta_s(l,j);
                end %#zhu#% theta_s(l,j), when c!=1, theta_s should be c by K, because different block share different error.
            end %#zhu#% 
            L_chol=diag(repvec(sqrt(theta_q(:,j)),m)); % L is supposed to be a lower triangular matrix, chol decomp of G(random effect cov). 
            temp1=eye(M)+L_chol*W0(Zrange,Zrange)*L_chol; %#zhu#% For this part, refer Wolfinger,1994. 
            temp2=L_chol*W0(Zrange,:);
            A=[temp1,temp2;temp2',W0];
            L1(j)=0;
            [temp i_sort]=sort(diag(A(1:M,1:M))); %%%% sweep largest diagonals first
            piv=NaN(M,1);
            for ii=1:M   %#zhu#% default sort is to sort in ascending order. 
                A=sweep(A,i_sort(M-ii+1));
                piv(i_sort(M-ii+1),1)=1/A(i_sort(M-ii+1),i_sort(M-ii+1));
            end
            if c==1
                L1(j)=sum(log(piv))+n*log(theta_s(j));
            else
                L1(j)=sum(log(piv))+sum(log(model.C*theta_s(:,j))); %#zhu#% When c>1, the R matrix in sweep operator is not homogeneous in diagonal.
            end
            W1=A((M+1):end,(M+1):end);
            XvX(:,:,j)=W1(Xrange,Xrange);
            XvZ(:,:,j)=W1(Xrange,Zrange);
            ZvZ(:,:,j)=W1(Zrange,Zrange);    
            Xvd(:,j)=W1(Xrange,Drange);
            Zvd(:,j)=W1(Zrange,Drange);
            dvd(:,j)=W1(Drange,Drange);    
      else  %%% must change here too 1/13/05  %#zhu#% When model.H=0 or max(theta_q(:,j))<espsilon, there is no random effect. Sigma_{jk}=s_{j,k}I
            XvX(:,:,j)=temp_W{1}.XtX/theta_s(1,j);
            %XvZ(:,:,j)=temp_W{1}.XtZ/theta_s(1,j);
            %ZvZ(:,:,j)=temp_W{1}.ZtZ/theta_s(1,j);
            Xvd(:,j)=temp_W{1}.XtD(:,j)/theta_s(1,j);
            %Zvd(:,j)=temp_W{1}.ZtD(:,j)/theta_s(1,j);
            dvd(:,j)=temp_W{1}.DtD(j,j)/theta_s(1,j);
            L1(j)=sum(C(:,1))*log(theta_s(1,j));
            if c>1
                for l=2:c
                    XvX(:,:,j)=XvX(:,:,j)+temp_W{l}.XtX/theta_s(l,j);
                    %XvZ(:,:,j)=XvZ(:,:,j)+temp_W{l}.XtZ/theta_s(l,j);
                    %ZvZ(:,:,j)=ZvZ(:,:,j)+temp_W{l}.ZtZ/theta_s(l,j);
                    Xvd(:,j)=Xvd(:,j)+temp_W{l}.XtD(:,j)/theta_s(l,j);
                    %Zvd(:,j)=Zvd(:,j)+temp_W{l}.ZtD(:,j)/theta_s(l,j);
                    dvd(:,j)=dvd(:,j)+temp_W{l}.DtD(j,j)/theta_s(l,j);
                    L1(j)=L1(j)+sum(C(:,l))*log(theta_s(l,j));
                end
            end
      end
 end % j for loop

%%% Now add part to update betans and Vbetans (check if necessary).

if Update_betans==1
    betans=NaN(p,K);
    Vbetans=NaN(p,K);
    for j=1:K
          XvXk=XvX(:,:,j);
          %L2(kctr)=dvd{kctr}+beta(:,kctr)'*XvXk*beta(:,kctr)-2*beta(:,kctr)'*Xvd{kctr};
          temp=[XvXk,Xvd(:,j);Xvd(:,j)',dvd(:,j)]; %#zhu#% for every j or k, either d will change or both XvX and d will change.
          for ii=1:p
             temp=sweep(temp,ii);
          end
          betans(:,j)=temp(1:p,p+1); %#zhu#% the upper right corner is the estimated MLE. This is not computed as the iterative method stated on Page 188, Morris2006.
            %Vbetans(:,kctr)=diag(temp(1:p,1:p));
          Vbetans(:,j)=(diag(XvXk)).^(-1);  %%% Vbeta should be computed marginally %#zhu#% Note this is diag(X_i'\Sigma_{j,k}^{-1}X_i).^(-1).
                                             %%% i.e. (diag(XvX))^(-1), not
                                            %%% diag( (XvX)^(-1) ) 12/2/03
        % Vbetans(:,j)=diag(temp(1:p,1:p)); % This gives higher variance.
    end
else
    Vbetans=repmat(0,p,K);
    for k=1:K
        Vbetans(:,k)=diag(XvX(:,:,k)).^(-1);
    end
    betans=[];
end

Wv.XvX=XvX;
Wv.Xvd=Xvd;
if model.H~=0
    Wv.XvZ=XvZ;
    Wv.ZvZ=ZvZ;
    Wv.Zvd=Zvd;
end
Wv.dvd=dvd;
Wv.L1=L1;
