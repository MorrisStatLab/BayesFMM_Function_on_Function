function delta=MH_stepsize(theta_sd,theta_mle,multiple,j)
% This function compute the step size for Metropolis-Hastings update(proposal distribution) 
% based on MLE estimations of Theta and their standard deviation estimates.
% Input:
%       theta_sd:  standard error.
%       theta_mle: MLE estimate.
%       multiple: the multiple of the standard deviation to control step
%       size.
%       j: j=0, compute for all j,k. 
%          j=j, compute for the jth column only. in this case, theta_sd and
%          theta_mle is the jth column of the whole matrix.
% Output:
%       delta: the step size used for log normal proposal. 
%       log sigma^2_new=log sigma^2_old+eps, where eps~N(0, delta^2), for
%       all j,k. 

if j==0
    delta=sqrt(log(0.5*(1+sqrt(1+4*theta_sd.^2*multiple^2./theta_mle.^2))));
elseif j>0
    delta=sqrt(log(0.5*(1+sqrt(1+4*theta_sd.^2*multiple^2./theta_mle.^2))));
end