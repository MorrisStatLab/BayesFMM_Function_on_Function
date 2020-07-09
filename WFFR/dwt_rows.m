function [D,wavespecs]=dwt_rows(infile,wavespecs)

% dwt_rows(infile,wavespecs): Compute DWT on rows of infile, using wavelet and number of levels specified.
%
% Input: infile : matrix, with each row containing a function on which to perform DWT (n x k)
%        wavespecs: structure containing wavelet settings, including
%           nlevels: Number of levels of wavelet transform
%           wavelet: Name of wavelet basis to use
%        wtmode: a string, indicates which type of Discrete wavelet transform
%        extension mode to use. Can be: 'sym' (default),'per','zpd'.. See
%        help dwtmode.
% Output: D : matrix with each row containing the wavelet coefficients for the function (n x k)
%         K : vector of number of coefficients per level (1 x nlevels)
%         T : number of sampled points in function
% hello
dwtmode(wavespecs.wtmode);
[d,K]=wavedec(infile(1,:),wavespecs.nlevels,wavespecs.wavelet);
D=NaN(size(infile,1),sum(K)-size(infile,2));
T=K(end);
K=K(1:(end-1))';
wavespecs.Kj=K;
wavespecs.T=T;
wavespecs.J=length(K);
wavespecs.K=sum(K);

D(1,:)=d;
for i=2:size(infile,1)
    [D(i,:),]=wavedec(infile(i,:),wavespecs.nlevels,wavespecs.wavelet); 
end


