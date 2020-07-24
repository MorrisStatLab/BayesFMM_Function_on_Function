function [ samp ] = rsample( N, K, replacement )
%RSAMPLE Sample N elements of K with or without replacement.
%          N and K are integer values. Sampling is based on
%          RANDPERM(K), type 'help RANDPERM' for details.
%
%   SAMP = RSAMPLE(N, K) gives N samples from K without
%          replacement, NOTE N must be <= K.
%
%   SAMP = RSAMPLE(N, K, replacement) with replacement = true,
%          gives N samples from K with replacement.
%
%   Created: 10/10/2012
%   By: Mark John Meyer

if nargin < 3;
    replacement = false;
end;

if replacement;
    samp   = zeros(1,N);
    for i = 1:N;
        temp_samp   = randperm(K);
        samp(i)     = temp_samp(1);
    end;
else
    if N > K;
        error('Requested more objects sampled than allowed without replacement, N > K');
    end;
    temp_samp   = randperm(K);
    samp        = temp_samp(1:N);
end;

end

