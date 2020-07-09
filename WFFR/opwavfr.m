function [ res ] = opwavfr( Y,model,wavespecsy,MCMCspecs,type )
%OPWAVFR Performs Ordinal Probit Wavelet-based Functional Regression for
%           either a set of scalar predictors or functional predictors.
%           Calls OPWAVFM for function-on-scalar regression and OPWAVFFR 
%           for function-on-function regression. Default is function-on-
%           scalar. Use type = 'ffr' for function-on-function.
%
%   Created: 3/29/2014
%   By: Mark John Meyer

if nargin < 5
    type    = 'fm';
end

switch type
    % OPWAVFM
    case 'fm'
        res         = opwavfm(Y,model,wavespecsy,MCMCspecs);
    % OPWAVFFR
    case 'ffr'
        wavespecsx  = model.wavespecsx;
        pcaspecsx   = model.pcaspecsx;
        res         = opwavffr(Y,model,wavespecsy,wavespecsx,pcaspecsx,MCMCspecs);
end

end

