function [ tri ] = triFun( seq, mu, sig )
%TRIFUN Generates a non-central triangle function.
%       MU is the centrality parameter, SIG is the
%       scale parameter, and SEQ is the support of
%       the function.
%   
%   Created: 3/31/2014
%   By: Mark John Meyer
    
    tri     = max(1 - abs((seq-mu)/sig),0);
end

