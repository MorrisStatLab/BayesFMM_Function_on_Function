function postout = ...
            PostProcess_HFLM(MCMC_beta,MCMC_zeta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wpspecs)
    %% Outputs:
    %  thetahat,theta_025CI,theta_975CI,tauhat,pihat,Sigma,uhat,u_025CI,u_975CI,g_Ut,g_Ut025,g_Ut975,ghatns
    %% Eventual Inputs
    %  MCMC_beta,MCMC_zeta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wpspecs
    
    %% Test:
%     MCMC_beta           = res.MCMC_beta;
%     MCMC_zeta           = res.MCMC_zeta;
%     model               = res.model;
%     wpspecs             = res.wpspecs;
%     MCMC_tau            = res.MCMC_tau;
%     MCMC_pi             = res.MCMC_pi;
%     MCMC_flag_theta     = res.MCMC_flag_theta;
%     MCMC_theta          = res.MCMC_theta;
%     wpspecs.V           = res.wpspecs.wpspecsy.trackerJacker{1}.T;
%     wpspecs.T           = res.wpspecs.wpspecsx.trackerJacker{1}.T;

    %% function parameters
%     sampCoef    = model.sampCoef;
    p           = model.p;
    K           = wpspecs.K;
    V           = wpspecs.V; %% use V for fsize, may need to change if we explore data reduction in X
    T           = wpspecs.T;
    fsize       = size(model.Dx,2);
    w           = model.p - fsize; %% wsize
    wpspecsy    = wpspecs.wpspecsy;
    wpspecsx    = wpspecs.wpspecsx;
    keep        = model.keep; %% threshold size
%     H           = model.H;

    %%
    delt        = model.delt;
    alf         = model.alf;

    %% obtain wavelet-packet space estimates %%
    % functional covariates %
%     postout.bstar               = reshape(mean(MCMC_beta)',K,V)';
%     postout.bstar_025CI         = reshape(quantile(MCMC_beta,0.025)',K,V)';
%     postout.bstar_975CI         = reshape(quantile(MCMC_beta,0.975)',K,V)';
    
    % scalar covariates %
%     postout.zstar               = reshape(mean(MCMC_zeta)',K,w)';
%     postout.zstar_025CI         = reshape(quantile(MCMC_zeta,0.025)',K,w)';
%     postout.zstar_975CI         = reshape(quantile(MCMC_zeta,0.975)',K,w)';

    postout.accept_rate_theta   = mean(MCMC_flag_theta);
    postout.tauhat              = reshape(mean(MCMC_tau),p,wpspecs.J);
    postout.pihat               = reshape(mean(MCMC_pi)',p,wpspecs.J);
    postout.alpha               = reshape(mean(MCMC_alpha)',p,K)';
    
    %% perform IDWPT in both directions %%
    B       = size(MCMC_beta,1);
    Beta    = NaN(B,V*T);
    Zeta    = NaN(B,w*T);
%%
    for i = 1:B;
        %% extract coefficients of interest %%
        Betam       = reshape(MCMC_beta(i,:),fsize,K);
        Zm          = reshape(MCMC_zeta(i,:)',K,w)';
        
        %% add back thresholded coefs %%
        Betat           = zeros(K,K);
        Betat(keep,:)   = Betam;
        
        %% project Beta and Zeta back into Y space %%
        BZ          = [Betat; Zm];
        BZidwpt     = idwpt_rows(BZ,wpspecsy);
        
        yidwpt      = BZidwpt(1:K,:);
        Zidwpt      = BZidwpt((K+1):end,:);
                                    
        %% project Beta back into X space %%
        Bidwpt      = idwpt_rows(yidwpt',wpspecsx)';
        
        %% update Beta and Zeta %%
        Beta(i,:)   = reshape(Bidwpt,1,V*T);
        Zeta(i,:)   = reshape(Zidwpt,1,w*T);
        
        if mod(i,100) == 0;
            fprintf('\n Done projecting MCMC samples M = %d \n',i);
        end;
    end;

    %% average and find credible intervals %%
    postout.beta        = Beta;
    postout.zeta        = Zeta;
    postout.Q025_bhat   = reshape(quantile(Beta,0.025),V,T);
    postout.Q975_bhat   = reshape(quantile(Beta,0.975),V,T);
    postout.bhat        = reshape(mean(Beta),V,T);
    postout.Q025_zhat   = reshape(quantile(Zeta,0.025),T,w)';
    postout.Q975_zhat   = reshape(quantile(Zeta,0.975),T,w)';
    postout.zhat        = reshape(mean(Zeta),T,w)';
    fprintf('\n Done averaging and finding CIs.\n \n');
    
    %% coefficients for inference %%
    bh      = zeros(V,T);
    for i = 1:T;
        for j = i:V
            bh(i,j) = 1;
        end
    end
    bFlag   = reshape(bh,1,V*T);
    bINF    = Beta(:,bFlag == 1);
    
    %% implement FDR %%
    R       = sum(bFlag);
    pst     = NaN(R,1);
    lFDR    = NaN(R,1);
    for r = 1:R;
        Betar       = abs(bINF(:,r));
        pst(r)      = (1/B)*size(find(Betar > delt),1);
        if pst(r) == 1;
            pst(r)  = 1 - (2*B)^(-1);
        end;
        lFDR(r)     = 1 - pst(r);
    end;
    if sum(pst) == 0;
        error('No coef > delta');
    end;
    pr      = sort(pst,'descend');
    rstar   = cumsum(1-pr)./linspace(1,R,R)';
    gam     = find(rstar <= alf, 1, 'last' );
    if isempty(gam);
        phi = 1;
    else
        phi = pr(gam);
    end;
    psi     = pst >= phi;
    
    %% insert psi back in to full vector %%
    psif                = zeros(size(bFlag'));
    pstf                = zeros(size(bFlag'));
    psif(bFlag == 1)    = psi;
    pstf(bFlag == 1)    = pst;

    postout.psi     = reshape(psif,V,T);
    postout.pst     = reshape(pstf,V,T);
    fprintf('\n Done calculating FDR.\n \n');

    %% implement MAPs %%
    [MAPSb, upper_CIb, lower_CIb]   = jointband_maps(bINF,alf);
    MAPS                    = zeros(size(bFlag));
    MAPS(bFlag == 1)        = MAPSb;
    upper_CI                = zeros(size(bFlag));
    upper_CI(bFlag == 1)    = upper_CIb;
    lower_CI                = zeros(size(bFlag));
    lower_CI(bFlag == 1)    = lower_CIb;

    postout.MAPs                = reshape(MAPS,V,T);
    postout.UMAPs               = reshape(upper_CI,V,T);
    postout.LMAPs               = reshape(lower_CI,V,T);
    fprintf('\n Done calculating MAPs.\n \n');

    %% additional output of interest %%
%     if sampleU == 1
%         M                   = model.M;
%         postout.uhat        = reshape(mean(MCMC_U)',K,M);
%         postout.u_025CI     = reshape(quantile(MCMC_U,0.025)',K,M)';
%         postout.u_975CI     = reshape(quantile(MCMC_U,0.975)',K,M)';
%         g_Ut                = NaN(M,T);
%         g_Ut025             = NaN(M,T);
%         g_Ut975             = NaN(M,T);
%         for j = 1:M
%             gUtsample       = idwt_rows(MCMC_U(:,((j-1)*K+1):(j*K)),wavespecsy);
%             g_Ut(j,:)       = mean(gUtsample);
%             g_Ut025(j,:)    = quantile(gUtsample,0.025);
%             g_Ut975(j,:)    = quantile(gUtsample,0.975);
%         end
%         postout.g_Ut        = g_Ut;
%         postout.g_Ut025     = g_Ut025;
%         postout.g_Ut975     = g_Ut975;
%         fprintf('\n Done with mixed effects output.\n \n');
%     end
    
    %% output thetas %%
    thetahat                = reshape(mean(MCMC_theta),size(theta,1),size(theta,2));
    postout.thetahat        = thetahat;
    postout.theta_025CI     = reshape(quantile(MCMC_theta,0.05),size(theta,1),size(theta,2));
    postout.theta_975CI     = reshape(quantile(MCMC_theta,0.95),size(theta,1),size(theta,2));
    fprintf('\n Done with regularization parameters.\n \n');
    
    %% output sigmas if desired %%
%     if (get_sigma==1)
%         Sigma.td        = NaN(T,T,H+model.c); 
%         Sigma.diagtd    = NaN(T,H+model.c);
%         Sigma.Rhotd     = NaN(T,T,H+model.c);    
%         for h = 1:model.H
%             [Sigma.td(:,:,h),Sigma.diagtd(:,h),Sigma.Rhotd(:,:,h)]  = GetSigma(thetahat(h,:),wavespecsy);
%         end
%         for cc = 1:model.c
%             [Sigma.td(:,:,model.H+cc),Sigma.diagtd(:,model.H+cc),Sigma.Rhotd(:,:,model.H+cc)]   = GetSigma(thetahat(model.H+cc,:),wavespecsy);
%         end
%         postout.Sigma   = Sigma;
%     else
%         postout.Sigma=[];
%     end

