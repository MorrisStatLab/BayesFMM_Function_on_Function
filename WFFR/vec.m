function [ vA ] = vec( A )
%VEC Vectorizes the matrix A
%   returning a 1 x nrow(A)*ncol(A)
%   vector
%
%   Created: 5/2/2014
%   By: Mark John Meyer

nrow    = size(A,1);
ncol    = size(A,2);
vA      = reshape(A,1,nrow*ncol);

end

