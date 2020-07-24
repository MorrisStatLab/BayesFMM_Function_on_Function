function [W,Kj]=Get_DWT(wavelet,n,boundary,extend,nlevels);

%%%% [W,Kj]=Get_DWT(wavelet,n,nlevels,boundary): Compute W matrix
%%%%
%%%% Input: wavelet = name of wavelet filter
%%%%        n = length of signal
%%%%        nlevels = number of levels of decomposition (log_2((n-1)/(L-1)+1))) by
%%%%                    default, following Percival and Walden)
%%%%        boundary = boundary condition (reflection by default)
%%%%        settings = 1 if reflecting at endpoint, 0 if past
%%%%        extend = 1 means keep all floor((n_j-1)/2+N) coefficients, like
%%%%                    Matlab, 
%%%%                = 0 keep only floor(n_j/2) wavelet coeffs at level j
%%%%                = 10 keep all coeffs, then at end only keep middle
%%%%                    floor(n_j/2)
%%%%  Based on method described on page 116 of Vidakovic's 1999 textbook.
%%%%    Here, we split H_j = H_{jA} + H_{jB}, where H_{jA} contains the
%%%%    part of H_j that is based on "real" data, i.e. is independent of
%%%%    the boundary conditions.  The matrix H_{jB} then contains the 
%%%%    part of the filter that differs depending on the boundary
%%%%    conditions used.  This makes it easy to adapt the method to
%%%%    different boundary conditions.  Here we allow three: periodic,
%%%%    reflection, and pad with zeros
%%%%


if nargin<3
    boundary='reflection';
end;

if nargin<4
    extend=1;
end;

h=dbwavf(wavelet)*sqrt(2);
if nargin<5
    nlevels=floor(log((n-1)/(length(h)-1)+1)/log(2));
end;

N=length(h)/2;
J=nlevels;

nn=n;
Kj=repmat(0,J+1,1);
for (j=1:J);
    H1A=get_H1A(nn,h,extend); %%% Scaling filter
    G1A=get_H1A(nn,reverse(h),extend);  %%% Wavelet Filter
    H1B=get_H1B(nn,h,boundary,extend);
    G1B=get_H1B(nn,reverse(h),boundary,extend);
    H{j}=H1A+H1B; 
    G{j}=G1A+G1B;
    nn=size(H{j},1);
    Kj(J+2-j)=size(G{j},1);
    
    W=1;
    for (k=j:-1:1)
        W=[W*H{k};G{k}];
    end;

end;
Kj(1)=nn;

if extend==10
    [W,Kj]=junk(W,Kj,n);
end;



    figure(1)
    colormap(gray)
    imagesc(W'*W);
    colorbar
