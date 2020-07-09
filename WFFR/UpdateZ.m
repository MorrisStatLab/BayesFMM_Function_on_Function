function [ Zstar, W ] = UpdateZ( Y, model, beta, sigma, cuts )
%UPDATEZSTAR Update values of latent variable Z for OPWAVFM using
%               truncated gaussians based on data-space X, beta, and
%               variance estimates of Zstar
%
%   Created: 2/22/2014
%   By: Mark John Meyer

%% set parameters
n       = model.n;
K       = wavespecsy.K;
L       = model.L;
XB      = model.X*beta;
Zstar   = NaN(n,K);
% bound   = [-Inf cuts Inf];

%% update Zstar by columns of D, assume independence
for k = 1:K
    Yk      = Y();
    Zk      = XB(:,k);
    a       = [-Inf cuts];
    b       = [cuts Inf];
    
    
    
%     for l = 0:L
%         j       = l+1;
%         a       = bound(j);
%         b       = bound(j+1);
%         
%     end
    
end


end

