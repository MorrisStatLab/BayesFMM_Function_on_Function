function [ res ] = flmm_sim( Y, X, Z, FDR, SIMspecs, TSspecs, MCMCspecs, wavespecsx, wavespecsy, gbf, saveraw, pcaspecsx)
%function [ res ] = flmm(  )
%==========================================================================
% Note: flmm is based on wavfmm4_gbf.m.
% flmm performs function on function regression as described in Meyer,
% Coull, and Morris 2013. GBF stands for Generalized Basis Function. The
% user can specify a DWT only or PCA + DWT basis expansion of the predictor
% function. The outcome function is modeled using DWT. Model assumes T = V.
%==========================================================================
%   Detailed explanation goes here
    
    %% Save Unprocessed Output %%
    if nargin < 11;
        saveraw                     = 0;
    end;

    %% Check inputs %%
    if nargin < 10;
        %% Default to PCA on X
        gbf                         = 1;
    end;
    
    % Default wavelet specs for Y-space
    if nargin < 9;
        %% wavespecsy
        wavespecsy.nlevels          = 3;
        wavespecsy.wavelet          = 'db8';
        wavespecsy.wtmode           = 'sym';
    end;
    
    % Default wavelet specs for X-space
    if nargin < 8;
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
    if nargin < 7;
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
    if nargin < 6;
        TSspecs.twosample           = 0;
    end;

    % Error for missing Simulation Specs
    if nargin < 5;
        error('Simulation specs required for flmm_sim');
    end;
    
    % Error for Missing FDR settings
    if nargin < 4;
        error('FDR delta and alpha required');
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
    
    %% Center and Scale X and Y
    N           = size(X,1);
    t           = size(X,2);
    scaleX      = zeros(N,t);
    for i = 1:t,
        time_i          = X(:,i);
        tempStd         = std(time_i);
        tempMean        = mean(time_i);
        scaleX(:,i)     = (time_i-tempMean)/tempStd;
    end;
    fprintf('\n Done scaling X.\n \n');

    %% Center and Scale Y by time
%     scaleY      = zeros(N,t);
%     for i = 1:t,
%         time_i          = Y(:,i);
%         tempStd         = std(time_i);
%         tempMean        = mean(time_i);
%         scaleY(:,i)     = (time_i-tempMean)/tempStd;
%     end;
%     fprintf('\n Done scaling Y.\n \n');

    %% Check if random effects exist set sampleu accordingly %%
    if isempty(Z);
        MCMCspecs.sampleU           = 0;
    else
        MCMCspecs.sampleU           = 1;
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
        if nargin < 12;
            pcaspecsx.pca_cutoff    = [0.9, 0.95, 0.975, 0.99, 0.995, 0.999, 0.9999];
            pcaspecsx.pca_level     = 0.99;
        end;
        pca_cutoff              = pcaspecsx.pca_cutoff;
        [coef, score, eigen]    = princomp(D);
        comp_var                = cumsum(eigen)./sum(eigen);
        pca_keep                = zeros(size(pca_cutoff,2),1);
        for i = 1:size(pca_cutoff,2),
            pca_keep(i)         = find(comp_var > pca_cutoff(i), 1);
        end;
        pca_out                 = [pca_cutoff' pca_keep];

        %% Select PCA level %%
        pca_level               = pcaspecsx.pca_level;

        %% Select PCs to use %%
        X                       = score(:,1:pca_out(pca_out(:,1) == pca_level,2));

        %% Generate column mean matrix %%
        col_means               = mean(D);
        col_mean_mat            = zeros(t,size(D,2));
        for i = 1:t;
            col_mean_mat(i,:)   = col_means;
        end;

        %% Save PCA specs %%
        pcaspecsx.X             = X;
        pcaspecsx.pca           = 1;
        pcaspecsx.score         = score;
        pcaspecsx.coef          = coef;
        pcaspecsx.eigen         = eigen;
        pcaspecsx.mean_mat      = col_mean_mat;
        pcaspecsx.output        = pca_out;
    else
        pcaspecsx.pca            = 0;
    end;
    
    %% Twosample Design Matrix %%
    if twosample == 1;
        if gbf == 1;
            %% Set up big design matrix %%
            X0          = X(group == 0,:);
            X1          = X(group == 1,:);
            twoSize     = size(X0,2) + size(X1,2);
            X           = zeros(N, twoSize);

            %% Fill in big design matrix %%
            X(group == 0, 1:size(X0,2))                 = X0;
            X(group == 1, (size(X0,2)+1):(twoSize))     = X1;            
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
    if twosample == 1;
        %% Allows for group-specific intercepts %%
        a0  = zeros(N,1);
        a1  = zeros(N,1);
        a0(group == 0)  = 1;
        a1(group == 1)  = 1;
        if gbf == 1;
            model.X     = [a0 a1 X];
        else
            model.X     = [a0 a1 D];
        end;
    else
        %% Include intercept %%
        a   = ones(N,1);
        if gbf == 1;
            model.X     = [a X];
        else
            model.X     = [a D];
        end;
    end;
    model.Z{1}  = Z;
    model.C     = ones(size(Y,1),1);
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
    if length(FDR) == 3;
        model.delt  = FDR(1:2);
        model.alf   = FDR(3);
    end;
    if length(FDR) > 3;
        error('Too many FDR inputs');
    end;

    sampleU     = MCMCspecs.sampleU;
    get_sigma   = sampleU;

    %% Run Model %%
    %%% consider blocksize and paramroute or just remove entirely %%%
    res         = wavfmm4_sim(Y,model,wavespecsy,wavespecsc,pcaspecsx,MCMCspecs,sampleU,get_sigma,twosample,saveraw);

    %% RMSE %%
    %  extract data estimate and use true surface to calculate RMSE %
    trueB       = SIMspecs.trueB;
    bhat        = res.bhat;
    t           = size(bhat,1);
    one         = ones(t,1);
    mse         = one'*(bhat-trueB).^2*one;
    
    res.mse     = mse;

end

