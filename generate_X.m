%%%%%%%%%%%%%%%%%%%%%%%%
% Michele Zemplenyi
% 6/13/2020
% File to generate simulated exposure matrix
% using an autoregressive(1) covariance structure
%%%%%%%%%%%%%%%%%%%%%%%%

cd('C:/Users/Michele/Dropbox/Brent Coull/Code for Github/'); % path to save output
N = 400; % number of subjects
T = 90; % number of days of exposure
sigma       = 4;
rho         = 0.98;

% generate ar(1) covariance pattern %%
ar1Corr     = eye(T); % eye = identity matrix
for i = 1:T,
    for j = (i+1):T,
        ar1Corr(i,j) = rho^(j-i);
        ar1Corr(j,i) = rho^(j-i);
    end;
end;

ar1Cov      = sigma*ar1Corr;

% generate x(v) based on AR(1) estimated cov %%
simX        = NaN(N,T);
muX         = zeros(T,1);
for i = 1:N;
    simX(i,:)   = mvnrnd(muX,ar1Cov);
end;
csvwrite(sprintf('simX_N%d_T%d.csv', N, T), simX);