function [ res ] = flmm_cluster( Y, X, Z, FDR, TSspecs, MCMCspecs, wavespecsx, wavespecsy, gbf, pcaspecsx)
%function [ res ] = flmm(  )
%==========================================================================
% Note: flmm is based on wavfmm4_gbf.m.
% flmm performs function on function regression as described in Meyer,
% Coull, and Morris 2013. GBF stands for Generalized Basis Function. The
% user can specify a DWT only or PCA + DWT basis expansion of the predictor
% function. The outcome function is modeled using DWT. Model assumes T = V.
%==========================================================================
%   Detailed explanation goes here
    
    %% Internal model clock %%
    tic;

    %% Check inputs %%
    if nargin < 9;
        %% Default to PCA on X
        gbf                         = 1;
    end;
    
    % Default wavelet specs for Y-space
    if nargin < 8;
        %% wavespecsy
        wavespecsy.nlevels          = 3;
        wavespecsy.wavelet          = 'db8';
        wavespecsy.wtmode           = 'sym';
    end;
    
    % Default wavelet specs for X-space
    if nargin < 7;
        %% wavespecsx
        wavespecsx.comp_level       = 0.99;
        wavespecsx.alph             = [0.9,0.95,0.99,0.995,0.999,0.9999];
        wavespecsx.xticks           = [5 10 15 20 25 30 35 40 45 50 55];
        wavespecsx.nlevels          = 3;
        wavespecsx.wavelet          = 'db8';
        wavespecsx.wtmode           = 'sym';
        wavespecsx.ndim             = 1;
    end;
        
    % Default sampler specs
    if nargin < 6;
        %% MCMC specs
        MCMCspecs.B                 = 2000;
        MCMCspecs.burnin            = 5000;     % 10000
        MCMCspecs.thin              = 1;
        MCMCspecs.propsdTheta       = 1;
        MCMCspecs.nj_nosmooth       = 1;        % can be 0, if is 0, don't do any constraint for pi_{ij}, bigTau_{ij}.
        MCMCspecs.minp              = 1e-14;
        MCMCspecs.maxO              = 1e20;
        MCMCspecs.minVC             = 1e-20;
        MCMCspecs.VC0_thresh        = 1e-4;
        MCMCspecs.delta_theta       = 1e-4;
        MCMCspecs.thetaMLE_maxiter  = 10^6;
        MCMCspecs.EmpBayes_maxiter  = 10^6;
        MCMCspecs.time_update       = 100;
        MCMCspecs.tau_prior_var     = 1e3;      % the variance of tau_{ijk} when finding prior parameters for tau_{ijk}.
        MCMCspecs.tau_prior_idx     = 1;        % 1 indicate that a_tau and b_tau depend on ij, 0 indicate that they depend on jk. 
        MCMCspecs.PI_prior_var      = 0.06;     % this range should be in [0.02 0.09].
     end;
    
    % Assume one sample if not otherwise specified
    if nargin < 5;
        TSspecs.twosample           = 0;
    end;

    % Error for Missing FDR settings
    if nargin < 4;
        error('FDR delta and alpha required');
    end;

    %% Center and Scale X by time
    N           = size(X,1);
    t           = size(X,2);
    scaleX      = zeros(N,t);
    for i = 1:t,
        time_i          = X(:,i);
        tempStd         = std(time_i);
        tempMean        = mean(time_i);
        scaleX(:,i)     = (time_i-tempMean)/tempStd;
    end;
%     fprintf('\n Done scaling X.\n \n');

    %% Center and Scale Y by time
    scaleY      = zeros(N,t);
    for i = 1:t,
        time_i          = Y(:,i);
        tempStd         = std(time_i);
        tempMean        = mean(time_i);
        scaleY(:,i)     = (time_i-tempMean)/tempStd;
    end;
%     fprintf('\n Done scaling Y.\n \n');

    %% Check if random effects exist set sampleu accordingly %%
    if isempty(Z);
        MCMCspecs.sampleU           = 0;
    else
        MCMCspecs.sampleU           = 1;
    end;
    
    %% Check and assign Two Sample Specs
    twosample   = TSspecs.twosample;
    fprintf('\n Twosample = %d \n',twosample);
    if twosample == 1 && length(fieldnames(TSspecs)) < 2;
        error('Grouping variable required for two sample model');
    end;
    if twosample == 1;
        group   = TSspecs.group;
    end;
    
    %%
    alph                    = wavespecsx.alph;
    xticks                  = wavespecsx.xticks;
    
    %% Compression
    [compRes, compSet, D_all, C, keep, wavespecsc] = Wavelet_compress_v1208_rectangular(scaleX,wavespecsx,1,alph,xticks);

    %% Select Compression Level
    comp_level              = wavespecsx.comp_level;
    
    %% Additional Compression Output
    compOut.compRes         = compRes;
    compOut.compSet         = compSet;
    compOut.C               = C;

    %% Update wavespecs
    % D is compressed, wavelet transformed data
    wavespecsc.compOut      = compOut;
    wavespecsc.keep         = keep(alph == comp_level,:);
    wavespecsc.compress     = 1;
    D                       = D_all(:,wavespecsc.keep == 1);

    %% get new Kj
    [wavespecsc.Kj_comp, wavespecsc.Kj_all] = get_Kj_compress(wavespecsc);
    wavespecsc.Kj_all       = wavespecsc.Kj_all';
    wavespecsc.J_all        = length(wavespecsc.Kj_all);
    wavespecsc.J            = length(wavespecsc.Kj_comp);
    wavespecsc.K            = sum(wavespecsc.Kj_comp);
    
    %% Perform PCA on X after DWT %%
    if gbf == 1;
        % Default PCA specs for X-space
        if nargin < 10;
            pcaspecsx.pca_cutoff    = [0.9, 0.95, 0.975, 0.99, 0.995, 0.999, 0.9999];
            pcaspecsx.pca_level     = 0.99;
        end;
        pca_cutoff              = pcaspecsx.pca_cutoff;
        [coef, score, eigen]    = princomp(D);
        comp_var                = cumsum(eigen)./sum(eigen);
        pca_keep                = zeros(size(pca_cutoff,2),1);
        for i = 1:size(pca_cutoff,2),
            pca_keep(i)         = find(comp_var > pca_cutoff(i), 1 );
        end;
        pca_out                 = [pca_cutoff' pca_keep];

        %% Select PCA level %%
        pca_level               = pcaspecsx.pca_level;

        %% Select PCs to use %%
        X                       = score(:,1:pca_out(pca_out(:,1) == pca_level,2));

        %% Generate column mean matrix %%
        col_means               = mean(D);
        col_mean_mat            = zeros(t,size(D,2));
        for i = 1:t,
            col_mean_mat(i,:)   = col_means;
        end;

        %% Save PCA specs %%
        pcaX                    = X;
        pcaspecsx.pca           = 1;
        pcaspecsx.score         = score;
        pcaspecsx.coef          = coef;
        pcaspecsx.eigen         = eigen;
        pcaspecsx.mean_mat      = col_mean_mat;
        pcaspecsx.output        = pca_out;
    end;
    
    %% Twosample Design
    if twosample == 1;
        if gbf == 1;
            %% Set up big design matrix %%
            tsX         = zeros(N, 2*size(pcaX,2));
            group0      = find(group == 0);
            group1      = find(group == 1);

            %% Fill in big design matrix %%
            tsX(group1, 1:size(pcaX,2))                     = pcaX(group0, :);
            tsX(group0, (size(pcaX,2)+1):(2*size(pcaX,2)))  = pcaX(group1, :);
            
            %% Assign tsX to X
            X           = tsX;
        else
            %% Set up big design matrix %%
            tsX         = zeros(N, 2*size(D,2));
            group0      = find(group == 0);
            group1      = find(group == 1);

            %% Fill in big design matrix %%
            tsX(group0, 1:size(D,2))                    = D(group0, :);
            tsX(group1, (size(D,2)+1):(2*size(D,2)))    = D(group1, :);
            
            %% Assign tsX to D
            D           = tsX;
        end;
    end;

    %% Model Set up %%
    if gbf == 1;
        model.X     = X;
    else
        model.X     = D;
    end;
    model.Z{1}  = Z;
    model.C     = ones(size(scaleY,1),1);
    model.H     = 1;
    model.Hstar = 0;
    if length(FDR) == 1;
        model.delt  = FDR(1);
        model.alf   = 0.05;
    end;
    if length(FDR) == 2;
        model.delt  = FDR(1);
        model.alf   = FDR(2);
    end;
    if length(FDR) > 2;
        error('Too many FDR inputs');
    end;

    sampleU     = MCMCspecs.sampleU;
    get_sigma   = sampleU;

    %% Run Model %%
    %%% consider blocksize and paramroute or just remove entirely %%%
    res         = wavfmm4_cluster(scaleY,model,wavespecsy,wavespecsc,pcaspecsx,MCMCspecs,sampleU,get_sigma,twosample);
    res.time    = toc;

end

