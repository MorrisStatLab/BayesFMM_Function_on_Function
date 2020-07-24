function A=rep(a,ntimes)

% rep(a,ntimes): constructs a matrix repeating each column a(i) ntimes(i)
%               times
%     Needs function uneqkron.m
%

A=a*uneqkron(ntimes)';


    