function postout = ...
            PostProcess_cluster(MCMC_beta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wavespecsc,wavespecsy,pcaspecsx,twosample,get_sigma,sampleU,MCMC_U)
    %% Outputs:
    %  thetahat,theta_025CI,theta_975CI,tauhat,pihat,Sigma,uhat,u_025CI,u_975CI,g_Ut,g_Ut025,g_Ut975,ghatns
    %% Inputs
    %  MCMC_theta,betans,theta,get_sigma,sampleU,MCMC_U,MCMC_tau,MCMC_pi
    
    %% Test:
%     MCMC_beta     = res.MCMC_beta;
%     model         = res.model;
%     wavespecsy    = res.wavespecs;
%     MCMC_tau      = res.MCMC_tau;
%     MCMC_pi       = res.MCMC_pi;

    %%
    p           = model.p;
    T           = wavespecsy.T; % post PCA+DWT value of T, total Y time
    K           = wavespecsy.K;
    Tx          = wavespecsc.T;
    H           = model.H;
    delt        = model.delt;
    alf         = model.alf;

    %% Obtain wavelet/PCA space estimates
    postout.bstarhat            = reshape(mean(MCMC_beta)',K,p)';
    postout.bstar_025CI         = reshape(quantile(MCMC_beta,0.025)',K,p)';
    postout.bstar_975CI         = reshape(quantile(MCMC_beta,0.975)',K,p)';
    postout.alphahat            = reshape(mean(MCMC_alpha)',K,p)';
    postout.accept_rate_theta   = mean(MCMC_flag_theta);

    postout.tauhat              = reshape(mean(MCMC_tau),p,wavespecsy.J);
    postout.pihat               = reshape(mean(MCMC_pi)',p,wavespecsy.J);
    
    %% Unlist x-space PCA specs
    if pcaspecsx.pca == 1;
        scorex      = pcaspecsx.score;
        coefx       = pcaspecsx.coef;
        meansx      = pcaspecsx.mean_mat;
        pcalevelx   = pcaspecsx.pca_level;
    end;
    
    %% Update wavelet specs
    if wavespecsc.compress == 1;
        wavespecsc.Kj = wavespecsc.Kj_all;
    end;

    %% Perform IDWT in both directions
    %  Consider both one and two sample case
    %  For twosample, perform inference
    %  twosample defaults to 0
    B       = size(MCMC_beta,1);
    if twosample == 0;
        Beta    = NaN(B,Tx*Tx);
        for i = 1:B;
            Betai   = reshape(MCMC_beta(i,:)',K,p)';
            
            %% Project back into Y space
            yidwt   = idwt_rows(Betai,wavespecsy);
                        
            %% PCA and Wavelet inverse
            if pcaspecsx.pca == 1;
                pcaCol              = pcaspecsx.output(pcaspecsx.output(:,1) == pcalevelx,2);
                temp                = zeros(size(scorex,2),Tx);
                temp(1:pcaCol,:)    = yidwt;
                Betai               = (temp'/coefx + meansx)';                
            end;
            if wavespecsc.compress == 1;
                temp                        = zeros(size(Betai',1),length(wavespecsc.keep));
                temp(:,wavespecsc.keep==1)  = Betai';
                Betai                       = temp;
            end;
            Betai       = idwt_rows(Betai,wavespecsc)';
            Beta(i,:)   = reshape(Betai,1,Tx*Tx);

%             if mod(i,100) == 0;
%                 fprintf('\n Done projecting MCMC samples M = %d \n',i);
%             end;       
        end;
    else
        Beta1   = NaN(B,Tx*Tx);
        Beta2   = NaN(B,Tx*Tx);
        Diff    = NaN(B,Tx*Tx);        
        for i = 1:B;
            Betai = reshape(MCMC_beta(i,:)',K,p)';
            
            %% Split samples
            group_p     = p/2; % assumes equal number of time points for each group
            beta_temp1  = Betai(1:group_p,:);
            beta_temp2  = Betai((group_p+1):end,:);

            %% Project back into Y space
            yidwt1      = idwt_rows(beta_temp1,wavespecsy);
            yidwt2      = idwt_rows(beta_temp2,wavespecsy);

            %%
            if pcaspecsx.pca ==1;
                %%  get # of PCA columns
                pcaCol              = pcaspecsx.output(pcaspecsx.output(:,1) == pcalevelx,2);
                
                %% Error check: make sure pcaCol == group_p
                if pcaCol ~= group_p;
                    error('PCA dimensions do not match group dimensions')
                end;

                %% insert zeros for unused components
                temp1               = zeros(size(scorex,2),Tx);
                temp2               = zeros(size(scorex,2),Tx);
                temp1(1:pcaCol,:)   = yidwt1;
                temp2(1:pcaCol,:)   = yidwt2;

                %% project back into data space
                beta_temp1          = (temp1'/coefx + meansx)';
                beta_temp2          = (temp2'/coefx + meansx)';
            else
                beta_temp1          = yidwt1;
                beta_temp2          = yidwt2;
            end;
            if wavespecsc.compress == 1;
                % first surface
                temp1                           = zeros(size(beta_temp1',1),length(wavespecsc.keep));
                temp1(:,wavespecsc.keep==1)     = beta_temp1';
                beta_temp1                      = temp1;
                % second surface
                temp2                           = zeros(size(beta_temp2',1),length(wavespecsc.keep));
                temp2(:,wavespecsc.keep==1)     = beta_temp2';
                beta_temp2                      = temp2;
            end;
            xidwt1      = idwt_rows(beta_temp1,wavespecsc)';
            xidwt2      = idwt_rows(beta_temp2,wavespecsc)';
            
            %% IDWT on X scale
            Beta1(i,:)  = reshape(xidwt1,1,Tx*Tx);
            Beta2(i,:)  = reshape(xidwt2,1,Tx*Tx);
            Diff(i,:)   = Beta1(i,:) - Beta2(i,:);
                       
%             if mod(i,100) == 0;
%                 fprintf('\n Done projecting MCMC samples M = %d \n',i);
%             end;
        end;
    end;
    
    %% Implement FDR
    if twosample == 0;
        R       = size(Beta,2);
        pst     = NaN(R,1);
        lFDR    = NaN(R,1);
        for r = 1:R;
            Betar       = abs(Beta(:,r));
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
    else
        R   = size(Diff,2);
        pst     = NaN(R,1);
        lFDR    = NaN(R,1);
        for r = 1:R;
            Diffr       = abs(Diff(:,r));
            pst(r)      = (1/B)*size(find(Diffr > delt),1);
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
    end;
    
    postout.psi     = reshape(psi,Tx,Tx);
%     fprintf('\n Done calculating FDR.\n \n');
    
    %% Then average and find credible interval, update postout
    
    if twosample == 0;
        postout.Q025_bhat   = reshape(quantile(Beta,0.025),Tx,Tx);
        postout.Q975_bhat   = reshape(quantile(Beta,0.975),Tx,Tx);
        postout.bhat        = reshape(mean(Beta),Tx,Tx);
    else
        postout.Q025_bhat1  = reshape(quantile(Beta,0.025),Tx,Tx);
        postout.Q975_bhat1  = reshape(quantile(Beta,0.975),Tx,Tx);
        postout.bhat1       = reshape(mean(Beta),Tx,Tx);
        postout.Q025_bhat2  = reshape(quantile(Beta,0.025),Tx,Tx);
        postout.Q975_bhat2  = reshape(quantile(Beta,0.975),Tx,Tx);
        postout.bhat2       = reshape(mean(Beta),Tx,Tx);
        postout.Q025_diff   = reshape(quantile(Diff,0.025),Tx,Tx);
        postout.Q975_diff   = reshape(quantile(Diff,0.975),Tx,Tx);
        postout.diff        = reshape(mean(Diff),Tx,Tx);
    end;
%     fprintf('\n Done averaging and finding CIs.\n \n');

    %% Additional output of interest
    if sampleU == 1
        M                   = model.M;
        postout.uhat        = reshape(mean(MCMC_U)',K,M);
        postout.u_025CI     = reshape(quantile(MCMC_U,0.025)',K,M)';
        postout.u_975CI     = reshape(quantile(MCMC_U,0.975)',K,M)';
        g_Ut                = NaN(M,T);
        g_Ut025             = NaN(M,T);
        g_Ut975             = NaN(M,T);
        for j = 1:M
            gUtsample       = idwt_rows(MCMC_U(:,((j-1)*K+1):(j*K)),wavespecsy);
            g_Ut(j,:)       = mean(gUtsample);
            g_Ut025(j,:)    = quantile(gUtsample,0.025);
            g_Ut975(j,:)    = quantile(gUtsample,0.975);
        end
        postout.g_Ut        = g_Ut;
        postout.g_Ut025     = g_Ut025;
        postout.g_Ut975     = g_Ut975;
%         fprintf('\n Done with mixed effects output.\n \n');
    end
    
    %%
    thetahat                = reshape(mean(MCMC_theta),size(theta,1),size(theta,2));
    postout.thetahat        = thetahat;
    postout.theta_025CI     = reshape(quantile(MCMC_theta,0.05),size(theta,1),size(theta,2));
    postout.theta_975CI     = reshape(quantile(MCMC_theta,0.95),size(theta,1),size(theta,2));
%     fprintf('\n Done with regularization parameters.\n \n');
    
    %%
    if (get_sigma==1)
        Sigma.td        = NaN(T,T,H+model.c); 
        Sigma.diagtd    = NaN(T,H+model.c);
        Sigma.Rhotd     = NaN(T,T,H+model.c);    
        for h = 1:model.H
            [Sigma.td(:,:,h),Sigma.diagtd(:,h),Sigma.Rhotd(:,:,h)]  = GetSigma(thetahat(h,:),wavespecsy);
        end
        for cc = 1:model.c
            [Sigma.td(:,:,model.H+cc),Sigma.diagtd(:,model.H+cc),Sigma.Rhotd(:,:,model.H+cc)]   = GetSigma(thetahat(model.H+cc,:),wavespecsy);
        end
        postout.Sigma   = Sigma;
    else
        postout.Sigma=[];
    end

