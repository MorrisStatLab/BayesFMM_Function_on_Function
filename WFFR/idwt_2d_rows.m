function Y=idwt_2d_rows(infile,wavespecs)

% dwt_rows(infile,nlevels,wavelet): Compute IDWT on rows of infile, using wavelet and 
%                                   number of levels specified.
%
% Input: infile : matrix, with each row containing wavelet coefficients for function on which to perform IDWT (n x k)
%        nlevels: Number of levels of wavelet transform
%       wavespecs: structure with information on wavelet settings,
%                   including:
%           wavelet: Name of wavelet basis to use
%
% Output: D : matrix with each row containing a function (n x T)
%         K : vector of number of coefficients per level (1 x nlevels)
%         T : number of sampled points per function 
%
% Once can change the extension mode  to periodict by using dwtmode('per').

Y=NaN(size(infile,1),wavespecs.t1*wavespecs.t2);
for i=1:size(infile,1)
    y=waverec2(infile(i,:),wavespecs.S,wavespecs.wavelet);
    Y(i,:)=reshape(y,1,wavespecs.t1*wavespecs.t2);
end



