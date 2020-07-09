function [ beta ] = wavphc( beta, wpspecsx, wpspecsy )
%WAVPHC Enforce historical constraint in wavelet-packet space
%   Detailed explanation goes here

% set up mods for constraint %
modx    = wpspecsx.trackerJacker{wpspecsx.totalDWT}.Kj(1);
mody    = wpspecsy.trackerJacker{wpspecsy.totalDWT}.Kj(1);

% zero out lower triangle for each node in Y and X %
for i = 1:size(beta,1)
    % mod i
    mi  = mod(i,modx);
    if mi == 0; mi = modx; end;
    for j = 1:size(beta,2)
        % mod j
        mj  = mod(j,mody);
        if mj == 0; mj = mody; end;
        if mj < mi
            beta(i,j) = 0;
        end;
    end;
end;


end

