function [Theta_mle,Theta_var,Theta_sd,beta_mle,betans,Vbetans,u_mle,Wv,converge]=mixed_mle(W,theta0,model,D,wavespecs,MCMCspecs)

%MIXED   Computes MLE by Henderson's Mixed Model Equations Algorithm.
%======================================================================
% Syntax:  
%   [Theta_mle,Theta_var,Theta_sd,beta_mle,betans,u_mle,loglik,converge]=mixed_mle(W,theta0
%   ,model,wavespecs,MCMCspecs);
%
%======================================================================
% Model: Y=X*b+Z*u+e,
%        u=(u_1',...,u_H')',
%        E(u)=0, Var(u)=diag(sigma^2_h*I_{m_h}), i=1,...,H.
%        E(e)=0, Var(e)=sigma^2_{H+1}*I_n,
%        Var(y)=Sig=sum_{h=1}^{H+1} sigma^2_h*Sig_h.
%        We assume normality and independence of u and e.
%
% Inputs:  
%   W       - The W matrix contains [X'X,X'Z,X'D;...] obtained from GetW()
%   theta0  - Initial values of variance components, H+c by K.
%   model   - model parameters. A structure contains X, Z, n, p ...
%   wavespecs
%           - wavelet transform parameters.     
%   MCMCspecs-contain the maximum iterations allowed.
%======================================================================
% Outputs: 
%   Theta_mle - Estimated vector of variance components, H+c by K 
%               Each column contains:(sigma^2_1,..., sigma^2_{H+c})'.              
%   b_mle     - p by K matrix of estimated fixed effects beta.           
%   u_mle     - M by K matrix of EBLUP's of random effects U,
%               each column contains: u=[u_1;...;u_r]. 
%   Theta_var - H+1 by H+1 by K array, Fisher information matrix for variance components;
%   loglik    - K by 1 vecoter, Log-likelihood evaluated at the estimated
%               parameters;
%   converge  - K by 1 vector, 1 indicates converged, 0 not converged.
%======================================================================
% REFERENCES
%
% Searle, S.R., Cassela, G., McCulloch, C.E.: Variance Components. 
% John Wiley & Sons, INC., New York, 1992. (pp. 275-286).
% 
% Witkovsky, V.: MATLAB Algorithm mixed.m for solving 
% Henderson's Mixed Model Equations.
% Technical Report, Institute of Measurement Science, 
% Slovak Academy of Sciences, Bratislava, Dec. 2001.
% See http://www.mathpreprints.com.
%  
% The algorithm mixed.m is a simplied version of Viktor Witkovsky's code 
% available at http://www.mathworks.com/matlabcentral/fileexchange
% see the Statistics Category.
%======================================================================
% Note by Hongxiao Zhu: We use the result that c==1, i.e. treat error
% variance to be sig_e^2I_n. In case of c>1, we just let sig_c^2==sig_e^2.
% To extend this code to c>1 case, we need to rederive all estimation 
% equations. Refer Searle 1992 and Witkovsky's notes for derivation. 
% I haven't found an easy way to solve it yet.
%======================================================================
K=wavespecs.K;
[n,p]=size(model.X);
m=model.m;
c=model.c;
H=model.H;
M=sum(m);
Im=eye(M);

Theta_mle=NaN(H+c,K);
beta_mle=NaN(p,K);
u_mle=NaN(M,K);
Theta_var=NaN(H+1,H+1,K); %#zhu#% here we assume R=sig0^2I_n. hence c==1
Theta_sd=NaN(H+c,K);
converge=NaN(K,1);

for j=1:K 
    
    a=[W.XtD(:,j);W.ZtD(:,j)];     
    s2=theta0(:,j);   
    %======================================================================
    %	START OF THE MAIN LOOP
    %======================================================================
    epss=1e-6; % Given precission for stopping rule
    crit=1;
    loops=0;
    while (crit>epss)&&(loops<MCMCspecs.thetaMLE_maxiter)
          loops=loops+1;
          sigaux=s2;
          s0=mean(s2(H+1:H+c));
          Dd=diag(repvec(s2(1:H),m));
          V=s0*Im+W.ZtZ*Dd;         
          WW=s0*pinv(V);
          A=[W.XtX W.XtZ*Dd;W.XtZ' V];
          bb=pinv(A)*a;
          b=bb(1:p);
          v=bb(p+1:p+M);
          u=Dd*v;
    %======================================================================
    % 	ESTIMATION OF ML OF VARIANCE COMPONENTS 
    %======================================================================
          hupp=0;
          for h=1:H
              hlow=hupp+1;
              hupp=hupp+m(h);
              Whh=WW(hlow:hupp,hlow:hupp);      
              uh=u(hlow:hupp);
              utu=uh'*uh;             
              s2(h)=max(min(utu/(m(h)-trace(Whh)),1e20),1e-20);              
          end                
          s2(H+1:H+c)=max(min((W.DtD(j,j)-b'*W.XtD(:,j)-u'*W.ZtD(:,j))/n,1e20),1e-20); %#zhu#% assume that R is sig_e^2*In. If not, need new derivations.
          crit=norm(sigaux-s2);         
    end
    
    if loops==MCMCspecs.thetaMLE_maxiter
        converge(j)=0;
    elseif crit<=epss
        converge(j)=1;
    end
    Theta_mle(:,j)=s2;
    %======================================================================
    %	END OF THE MAIN LOOP
    %======================================================================    
    s0=mean(s2(H+1:H+c));
    Dd=diag(repvec(s2(1:H),m));
    V=s0*Im+W.ZtZ*Dd;    
    WW=s0*pinv(V);   
    A=[W.XtX W.XtZ*Dd;W.XtZ' V];    
    bb=pinv(A)*a;
    b=bb(1:p);
    v=bb(p+1:p+M);
    u=Dd*v;
    beta_mle(:,j)=b;
    u_mle(:,j)=u;
    %======================================================================
    %	FISHER INFORMATION MATRIX FOR VARIANCE COMPONENTS
    %======================================================================
    Is2=NaN(H+1);  %#zhu#% again, we assume c==1.
    Is2(H+1:H+1,H+1:H+1)=max((n-M+trace(WW*WW))/(s2(H+1)^2),1e-20);
   
    hupp=0;
    for h=1:H
       hlow=hupp+1;
       hupp=hupp+m(h);
       trhh=trace(WW(hlow:hupp,hlow:hupp));
       trsum=0;
       lupp=0;
       for l=1:H,
           llow=lupp+1;
           lupp=lupp+m(l);
           tr=trace(WW(hlow:hupp,llow:lupp)*WW(llow:lupp,hlow:hupp));
           trsum=trsum+tr;
           Is2(h,l)=((h==l)*(m(h)-2*trhh)+tr)/(s2(h)*s2(l));
       end
       Is2(H+1,h)=(trhh-trsum)/(s2(H+1)*s2(h));
       Is2(h,H+1)=Is2(H+1,h);
    end
    asymp_var=max(2*pinv(Is2),MCMCspecs.minVC); % the asymptotic variance is the inverse of fisher Information matrix. 
    Theta_var(:,:,j)=asymp_var;          % For this part see Page 6 of Witkovsky's notes and Page 239 of Searle 1992.
    Theta_sd(1:H+1,j)=sqrt(diag(asymp_var));
    if c>1
     Theta_sd(H+2:H+c,j)=repmat(sqrt(asymp_var(H+1,H+1)),c-1,1);
    end
end % end of for j=1:K.

[betans,Vbetans,Wv]=GetGCP_byblocks(Theta_mle,D,W,model,wavespecs,1,MCMCspecs);

fprintf('\n Convergence rate for all (j,k) is %4.2f \n \n',1-sum(converge==0)/length(converge));
%======================================================================
%	EOF MIXED.M
%======================================================================