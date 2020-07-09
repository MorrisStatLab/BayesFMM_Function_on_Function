function [ beta, sigma ] = BTPost( bstar, theta, wavespecsy, wavespecsx, pcaspecsx )
%BTPOST Performs post processing on beta* and theta
%       returning beta and sigma in the data space
%
%   Created: 2/24/2014
%   By: Mark John Meyer

    %% Unlist x-space PCA specs
    pcalevelx   = pcaspecsx.pca_level;            
    scorex      = pcaspecsx.score;
    coefx       = pcaspecsx.coef;
    meansx      = pcaspecsx.mean_mat;
    Tx          = wavespecsx.T;
    
    %% Update wavelet specs
    if wavespecsx.compress == 1;
        wavespecsx.Kj = wavespecsx.Kj_all;
    end;

    %% Project beta and theta back into Y space
    yidwt       = idwt_rows(bstar,wavespecsy);
    sigma       = idwt_rows(theta,wavespecsy);

    %% PCA and Wavelet inverse
    pcaCol              = pcaspecsx.output(pcaspecsx.output(:,1) == pcalevelx,2);
    temp                = zeros(size(scorex,2),Tx);
    temp(1:pcaCol,:)    = yidwt;
    Betai               = (temp'/coefx + meansx)';                
    if wavespecsx.compress == 1;
        temp                        = zeros(size(Betai',1),length(wavespecsx.keep));
        temp(:,wavespecsx.keep==1)  = Betai';
        Betai                       = temp;
    end;
    beta       = idwt_rows(Betai,wavespecsx)';

end

