function [ yStar ] = UpdateLV( beta, Y, model, cuts )
%UPDATELV Update Latent Variable for OPWAVFM using truncated normals.
%           Normals truncated based on values of cuts. Currently only
%           designed for ordinal Y = {0, 1, 2}.
%
%   Created: 3/7/2014
%   By: Mark John Meyer

%% get parameters
X   = model.orgX;
N   = model.n;
O   = model.O;

%% find E(Y*)
EY  = X*beta;

%% vectorize Y and EY, set up yStar
Yvec    = reshape(Y,1,N*O);
EYvec   = reshape(EY,1,N*O);
yStar   = NaN(size(Yvec));
Nvec    = size(Yvec,2);

%% sample from truncated normals

% Y == 0 %
rand0   = 0 + (normcdf(cuts(1),EYvec,1)-0).*rand(1,Nvec);
norm0   = norminv(rand0,EYvec,1);
yStar(Yvec == 0)    = norm0(Yvec == 0);

% Y == 1 %
rand1   = normcdf(cuts(1),EYvec,1) + (normcdf(cuts(2),EYvec,1)-normcdf(cuts(1),EYvec,1)).*rand(1,Nvec);
norm1   = norminv(rand1,EYvec,1);
yStar(Yvec == 1)    = norm1(Yvec == 1);

% Y == 2 %
rand2   = normcdf(cuts(2),EYvec,1) + (1-normcdf(cuts(2),EYvec,1)).*rand(1,Nvec);
norm2   = norminv(rand2,EYvec,1);
yStar(Yvec == 2)    = norm2(Yvec == 2);

%% reshape yStar
% yStar   = reshape(yStar,N,O);

end

