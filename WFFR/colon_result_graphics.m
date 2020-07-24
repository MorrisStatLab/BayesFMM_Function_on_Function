%% Created:     8/20/2012
%% Modified:    8/20/2012
%% Graphic Script for Morris and Caroll 2006 Reproduction %%
% Add path and load data
addpath('/Users/markmeyer/Dropbox/Research/Dissertation/WFMM4');
load('/Users/markmeyer/Dropbox/Research/Dissertation/WFMM4/colon_wav4v2_result1.mat');

%% Plots for each variable of X with 95% credible interval %%
position    = linspace(0,1,256);

%% Graph of each fucntion with 95% credible interval
for i = 1:9,
    subplot(3,3,i)
    j = i+1;
    plot(position,res.Q025_ghat(j,:),'k--',position,res.ghat(j,:),'k-',position,res.Q975_ghat(j,:),'k--');
end;

%% Histogram of acceptance rate
hist(res.acpt_rate)