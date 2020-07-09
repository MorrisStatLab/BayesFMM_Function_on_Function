function [theta_mle,theta_sd,betans,Vbetans,Wv,beta_mle,Var_beta]=regression_mle(W,model,D,wavespecs,MCMCspecs)
% This function compute the maximum likelihood estimation of the model:
% d_{jk}=XB_{jk}+E_{jk}
% where E(E_{jk})=0, Var(E_{jk})=\sigma_{jk}^2.
% This is the special case of mixed effect model d_{jk}=XB_{jk}+ZU_{jk}+E_{jk}
% when Z=0.
[beta_mle,theta_mle0,XtX_inv]=sweep_simple_regress(W,model.n);
if model.c>1
    theta_mle=repmat(theta_mle0',model.c,1);
else
    theta_mle=theta_mle0';
end
Var_beta=diag(XtX_inv)*theta_mle0'; % The variance of beta_mle
%betans=beta_mle;  % This is just an approximation, this is not the real betans defined in Morris 2006.
%Vbetans=(diag(W.XtX)).^(-1)*theta_mle0'; % Vbetans_i(jk)=(X(i)'s_{jk}^(-1)X(i))^(-1)=s_{jk}*(X(i)'*X(i))^(-1).
theta_sd=sqrt(2/model.n)*theta_mle;

[betans,Vbetans,Wv]=GetGCP_byblocks(theta_mle,D,W,model,wavespecs,1,MCMCspecs);




    