function [ b ] = UpdateCuts( yStar, Y, a )
%UPDATECUTS Update cutpoints (a, b) for OPWAVFM first fixing a
%           for identifiability then sampling b from Unif(c1,c2)
%
%   Created 3/7/2014
%   By: Mark John Meyer

%% vectorize Y
N       = size(Y,1);
O       = size(Y,2);
Yvec    = reshape(Y,1,N*O);

%% find c1, c2 using yStar, Yvec, and a (usually set a = 0)
c1      = max(max(yStar(Yvec == 1)),a);
c2      = min(yStar(Yvec == 2));

%% sample b|y*,y,a ~ U(c1,c2)
b       = c1 + (c2-c1)*rand;

end

