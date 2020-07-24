function Rho=get_Rho(Sig)

%%% compute correlation matrix given covariance matrix

sig=diag(Sig);
Rho=diag(sig.^(-1/2))*Sig*diag(sig.^(-1/2));
