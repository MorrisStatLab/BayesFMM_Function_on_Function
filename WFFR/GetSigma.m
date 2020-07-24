function [Sigmat,Sigma,Rho]=GetSigma(theta,wavespecs)

%%%% Compute data-space correlation matrix (and variance function) corresponding to diagonal
%%%%        wavelet-space covariance matrix.
%%%%
%%%%    Inputs: 
%              cov = assumption on variance components
%                               (0=homoscesastic, 1=homoscedastic within levels, 2=heteroscedastic)
%%%%           theta: vector of length * (1, J, or K) containing
%%%%                    wavelet-space variance components
%%%%            wavespecs: structure with elements:
%%%%                nlevels = number of levels of decomposition
%%%%                wavelet = name of wavelet filter used
%%%%                k = number of wavelet coefficients per level.
%%%%
%%%%            T: number of points sampled per function (resolution of
%%%%                grid).
%%%%
%%%%    Outputs:    Sigma = 1 x n vector of variances
%%%%                Rho = n x n correlation matrix.
%%%%
%%%%    Functions called: Form_2d_dwt, waverec2

A=form_2d_dwt(theta,wavespecs.Kj);
Sigmat=waverec2(A,[[wavespecs.Kj;wavespecs.T],[wavespecs.Kj;wavespecs.T]],wavespecs.wavelet);

Sigma=diag(Sigmat);
Rho=diag(Sigma)^(-.5)*Sigmat*diag(Sigma)^(-.5);

