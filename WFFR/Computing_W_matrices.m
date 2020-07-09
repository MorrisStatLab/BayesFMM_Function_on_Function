
wavelet='db3'
n=660
boundary='periodic  '
%boundary='reflection'
nlevels=8
extend=1
[W,Kj]=Get_DWT(wavelet,n,boundary,1,nlevels);
[W_0,Kj_0]=Get_DWT(wavelet,n,boundary,0,nlevels);
[W_10,Kj_10]=Get_DWT(wavelet,n,boundary,10,nlevels); %%%% do extended version, then keep only middle floor(n/2) coefficients/level

colormap(gray)
Sig=W_10*W_10';
sig=diag(Sig);
Rho=diag(sig.^(-1/2))*Sig*diag(sig.^(-1/2));
imagesc(Rho)
colorbar
eig(Rho)

y=normrnd(repmat(0,n,1),repmat(1,n,1));

[a,b]=dwt(y,wavelet);
[d,Kj]=wavedec(y,nlevels,wavelet);
[yy]=waverec(d,Kj,wavelet);
temp=[d,W*y];

figure(1)
imagesc(W_1'*W_1)
colorbar

figure(2)
imagesc(W_0'*W_0)
colorbar

