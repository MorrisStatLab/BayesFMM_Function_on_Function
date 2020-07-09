function [ pcaspecsx, wavespecsc ] = wPC( X, wavespecsx, pcaspecsx )
%WPC Performs wavelet-PC decomposition on X
%
%   Created: 2/24/2014
%   By: Mark John Meyer

%%
    alph                    = wavespecsx.alph;
    xticks                  = wavespecsx.xticks;
    
    %% Compression
    [compRes, compSet, D_all, C, keep, wavespecsc] = Wavelet_compress_v1208_rectangular(X,wavespecsx,1,alph,xticks);

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

    %% get new Kj, save wavespecs
    [wavespecsc.Kj_comp, wavespecsc.Kj_all] = get_Kj_compress(wavespecsc);
    wavespecsc.Kj_all       = wavespecsc.Kj_all';
    wavespecsc.J_all        = length(wavespecsc.Kj_all);
    wavespecsc.J            = length(wavespecsc.Kj_comp);
    wavespecsc.K            = sum(wavespecsc.Kj_comp);

    %% PCA
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
    Xwpc                    = score(:,1:pca_out(pca_out(:,1) == pca_level,2));

    %% Generate column mean matrix %%
    col_means               = mean(D);
    col_mean_mat            = zeros(size(X,2),size(D,2));
    for i = 1:size(X,2);
        col_mean_mat(i,:)   = col_means;
    end;

    %% Save PCA specs %%
    pcaspecsx.X             = Xwpc;
    pcaspecsx.pca           = 1;
    pcaspecsx.score         = score;
    pcaspecsx.coef          = coef;
    pcaspecsx.eigen         = eigen;
    pcaspecsx.mean_mat      = col_mean_mat;
    pcaspecsx.output        = pca_out;
end

