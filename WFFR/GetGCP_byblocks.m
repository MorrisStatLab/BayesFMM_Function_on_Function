function [betans,Vbetans,Wv]=GetGCP_byblocks(theta,D,W,model,wavespecs,Update_betans,MCMCspecs)
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
%               covQ = assumptions on within-curve random effect variance components q
%               covS = assumptions on within-curve residual variance components s
%                           (0=homoscesastic, 1=homoscedastic within levels, 2=heteroscedastic)
%               Hstar = number of Z matrices whose columns are independent
%                           subjects. (Assumed to be first Hstar<=H Z
%                           matrices).
%   
%           theta = (H+1 x *) matrix of within-curve variance component
%                           estimates.                                          
%                           (*=1 if model.covQ=0, J if model.covQ=1, K if model.covQ=2)
%
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
if model.H~=0
    model.M=sum(model.m);
    Wv.XvZ=NaN(model.p,model.M,wavespecs.K);
    Wv.ZvZ=NaN(model.M,model.M,wavespecs.K);  %#zhu#% save in cell form is not necessary, can use 3-arrays.    
    Wv.Zvd=NaN(model.M,wavespecs.K);
end
Wv.XvX=NaN(model.p,model.p,wavespecs.K);
Wv.Xvd=NaN(model.p,wavespecs.K);  %#zhu#% why Xvd,Zvd are not index by j=1,..K?
Wv.dvd=NaN(1,wavespecs.K);
Wv.L1=NaN(1,wavespecs.K);

Hstar=model.Hstar;
if (Hstar==0)
    [betans,Vbetans,Wv]=GetGCP(theta,W,model,wavespecs,Update_betans,D,MCMCspecs);
else %#zhu#% The case of Hstar~=0 and H>1 are Not tested.
    error('The case of Hstar>0 will be added in here');
%     m_ctr=0;  %%%% Ctr for which subject-level random effect is updated
%     m_star=sum(model.m(1:Hstar));   %%% total number of subject random effects.
%     for hstar=1:Hstar  %%% Have not tested for Hstar>1.
%         for m=1:model.m(hstar)
%             selected=(model.Z{hstar}(:,m)==1);    %#zhu#% Is it possible that select=[]? Assume not.    
%             temp_model.n=sum(selected);   %%%% pick off rows corresponding to this cluster
%             temp_X=model.X(selected,:);
%             %%%%% begin added 1/13/05
%             %%%%% end added 1/13/05
%             if sum(selected)>1
%                 X_cols=sum(abs(temp_X))>0;    %%% breaks when n=1. %#zhu#% A indicator that selects the columns that are not all zeros for the sub-X matrix.
%             else
%                 X_cols=abs(temp_X)>0;   %%%% discard columns of X containing all zeros (temporarily)
%             end
%             temp_model.p=sum(X_cols);
%             temp_model.X=temp_X(:,X_cols);
%             temp_model.H=model.H-Hstar+1;
%             temp_model.m=repmat(0,model.H-Hstar+1,1);
%             temp_model.VC0_thresh=model.VC0_thresh;
%             temp_model.Z{1}=model.Z{hstar}(selected,m);  %%% Single random effect for subject.
%             temp_model.m(1)=1;
%             temp_model.C=model.C(selected,:);
%             temp_model.c=model.c;
%             Z_cols=m_ctr+m;
%             m_ctr2=m_star; %%%% counter for individual random effects.
%             
%                 
%                for h=(Hstar+1):model.H %%%% Random effects nested w/in subject, if necessary
%                     temp_Z=model.Z{h}(selected,:);
%                     Zh_cols{h}=sum(abs(temp_Z))>0;
%                     temp_model.Z{h-Hstar+1}=temp_Z(:,Zh_cols{h});%#zhu#%  I don't understand Line 105 to Line 113
%                     temp_model.m(h-Hstar+1)=sum(Zh_cols{h});
%                     temp=(1:model.m(h));
%                     Z_cols=[Z_cols,m_ctr2+temp(Zh_cols{h})];
%                     m_ctr2=m_ctr2+model.m(h);
%                 end
%            
%             temp_W=GetW(temp_model,D(selected,:));
%             [betans,Vbetans,temp_Wv]=GetGCP(theta([hstar,(Hstar+1):end],:),temp_W,temp_model,wavespecs,0,D(selected,:));
%             Wv.Xvd(X_cols,:)=Wv.Xvd(X_cols,:)+temp_Wv.Xvd;
%             Wv.Zvd(Z_cols,:)=Wv.Zvd(Z_cols,:)+temp_Wv.Zvd;
%             Wv.dvd=Wv.dvd+temp_Wv.dvd;
%             Wv.L1=Wv.L1+temp_Wv.L1;
%             for k=1:wavespecs.K
%                 Wv.XvX{k}(X_cols,X_cols)=Wv.XvX{k}(X_cols,X_cols)+temp_Wv.XvX{k};
%                 Wv.XvZ{k}(X_cols,Z_cols)=Wv.XvZ{k}(X_cols,Z_cols)+temp_Wv.XvZ{k};
%                 Wv.ZvZ{k}(Z_cols,Z_cols)=Wv.ZvZ{k}(Z_cols,Z_cols)+temp_Wv.ZvZ{k};
%             end; %% k loop
%         end; %% m loop
%         m_ctr=m_ctr+model.m(hstar);
%     end %% hstar loop
end %%% Hstar>0 else

if (Update_betans==1)
    betans=repmat(0,model.p,wavespecs.K);
    Vbetans=betans;
    for k=1:wavespecs.K   
         Vbetans(:,k)=diag(Wv.XvX(:,:,k)).^(-1);
         betans(:,k)=pinv(Wv.XvX(:,:,k))*Wv.Xvd(:,k);
    end
else
    Vbetans=NaN(model.p,wavespecs.K);
    for k=1:wavespecs.K
        Vbetans(:,k)=diag(Wv.XvX(:,:,k)).^(-1);
    end
    betans=0;
end