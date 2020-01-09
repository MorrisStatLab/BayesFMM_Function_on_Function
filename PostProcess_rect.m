function postout = ...
            PostProcess_rect(MCMC_beta,MCMC_alpha,MCMC_flag_theta,MCMC_tau,MCMC_pi,MCMC_theta,theta,model,wavespecsc,wavespecsy,pcaspecsx,twosample,get_sigma,sampleU,MCMC_U)
    
    %%
    p           = model.p;
    nScalarCov  = model.nScalarCov; 
    T           = wavespecsy.T; % number of (compressed) sites in Y space
    K           = wavespecsy.K;
    Tx          = wavespecsc.T; % pertains to number of "time" points (compressed) in x space,
    H           = model.H;
    delt        = model.delt;
    alf         = model.alf;
    stdXt        = model.stdXt; 
    stdY         = model.stdY;  
    if twosample == 1;
        delt0   = model.delt0;
        delt1   = model.delt1;
    end;

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
        pcalevelx   = pcaspecsx.pca_level;            
        scorex      = pcaspecsx.score;
        coefx       = pcaspecsx.coef;
        meansx      = pcaspecsx.mean_mat;
    end;
    
    %% Update wavelet specs
    if wavespecsc.compress == 1;
        wavespecsc.Kj = wavespecsc.Kj_all;
    end;
    
    %% Prepare matrix of scaling factors (stdY/stdXt) to reverse the 
    % scaling of Y and Xt once mcmc samples are returned to data space 
    sf = NaN(Tx,T); % matrix to hold scaling factors 
    for t = 1:T 
        for v = 1:Tx
            sf(v,t) = stdY(t)/stdXt(v);
        end
    end
    sfv = reshape(sf,1, Tx* T); % vectorize matrix

    %% Perform IDWT in both directions
    %  Consider both one and two sample case
    %  For twosample, perform inference
    %  twosample defaults to 0
    B       = size(MCMC_beta,1);  %1000 
    if twosample == 0;
        Beta    = NaN(B,Tx*T); 
        A       = NaN(B,T); 
        BetaScalar = NaN(B, T*nScalarCov); 
        for i = 1:B;
            Betai   = reshape(MCMC_beta(i,:)',K,p)'; % p x K (because tranpose)
            
            %% Project back into Y space
            As          = Betai(1,:); % pertains to intercept
            A1          = idwt_rows(As,wavespecsy);
            if nScalarCov > 0;
                betaSc   = Betai(2:(1+nScalarCov), :); % add 1 to account for intercept 
                betaSc_idwt = idwt_rows(betaSc, wavespecsy);               
                beta_temp = Betai((2+nScalarCov):end, :); % corresponds to functional exposure part of surface
                yidwt       = idwt_rows(beta_temp,wavespecsy);
            else 
                beta_temp   = Betai(2:end,:); % no scalar covariates to account for
                yidwt       = idwt_rows(beta_temp,wavespecsy);
            end;
              
            
            
            %% Remove and update intercept and scalar covariates
            A(i,:)      = A1;
            if nScalarCov > 0;
                BetaScalar(i, :) = reshape(betaSc_idwt, 1, T*nScalarCov);               
            end;     

            %% PCA and Wavelet inverse - just performed on time-varying exposure stored in yidwt
            if pcaspecsx.pca == 1;
                pcaCol              = pcaspecsx.output(pcaspecsx.output(:,1) == pcalevelx,2);
                temp                = zeros(size(scorex,2),Tx); 
                temp(1:pcaCol,:)    = yidwt;
                Betai               = (temp'/coefx + meansx)';                
            end;
            if pcaspecsx.pca == 0; 
                Betai = yidwt;
            end;
            if wavespecsc.compress == 1;
                temp                        = zeros(size(Betai',1),length(wavespecsc.keep));
                temp(:,wavespecsc.keep==1)  = Betai';
                Betai                       = temp;
            end;
            Betai       = idwt_rows(Betai,wavespecsc)';
            Beta(i,:)   = reshape(Betai,1,Tx*T); 
            
            if mod(i,100) == 0;
                fprintf('\n Done projecting MCMC samples M = %d \n',i);
            end;       
        end;% end loop over the B MCMC samples
       
    % Below pertains to two-sample case which MZ did not edit
    else
        Beta0   = NaN(B,Tx*Tx);
        Beta1   = NaN(B,Tx*Tx);
        Diff    = NaN(B,Tx*Tx);
        A0      = NaN(B,Tx);
        A1      = NaN(B,Tx);
        for i = 1:B;
            Betai = reshape(MCMC_beta(i,:)',K,p)';
            
            %% Split off surfaces
            As          = Betai(1:2,:);
            cut_p       = size(pcaspecsx.X,2);
            beta_temp0  = Betai(3:(cut_p+2),:);
            beta_temp1  = Betai((cut_p+3):end,:);

            %% Project back into Y space
            As          = idwt_rows(As,wavespecsy);
            yidwt0      = idwt_rows(beta_temp0,wavespecsy);
            yidwt1      = idwt_rows(beta_temp1,wavespecsy);

            %% Remove and update intercepts
            A0(i,:)     = As(1,:);
            A1(i,:)     = As(2,:);

            %%
            if pcaspecsx.pca ==1;
                %%  get # of PCA columns
                pcaCol              = pcaspecsx.output(pcaspecsx.output(:,1) == pcalevelx,2);
                
                %% insert zeros for unused components
                temp0               = zeros(size(scorex,2),Tx);
                temp1               = zeros(size(scorex,2),Tx);
                temp0(1:pcaCol,:)   = yidwt0;
                temp1(1:pcaCol,:)   = yidwt1;

                %% project back into data space
                beta_temp0          = (temp0'/coefx + meansx)';
                beta_temp1          = (temp1'/coefx + meansx)';
            else
                beta_temp0          = yidwt0;
                beta_temp1          = yidwt1;
            end;
            if wavespecsc.compress == 1;
                %% first surface
                temp0                           = zeros(size(beta_temp0',1),length(wavespecsc.keep));
                temp0(:,wavespecsc.keep==1)     = beta_temp0';
                beta_temp0                      = temp0;
                %% second surface
                temp1                           = zeros(size(beta_temp1',1),length(wavespecsc.keep));
                temp1(:,wavespecsc.keep==1)     = beta_temp1';
                beta_temp1                      = temp1;
            end;
            %% X-space IDWT
            xidwt0      = idwt_rows(beta_temp0,wavespecsc)';
            xidwt1      = idwt_rows(beta_temp1,wavespecsc)';
            
            %% IDWT on X scale
            Beta0(i,:)  = reshape(xidwt0,1,Tx*Tx);
            Beta1(i,:)  = reshape(xidwt1,1,Tx*Tx);
            Diff(i,:)   = Beta1(i,:) - Beta0(i,:);
                       
            if mod(i,100) == 0;
                fprintf('\n Done projecting MCMC samples M = %d \n',i);
            end;
        end;
    end;
    
    %% Implement FDR
    if twosample == 0;
        if length(delt) == 1;
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
                fprintf('\n No coef > delta found for BFDR.\n \n');
                psi = zeros(Tx*T, 1); 
                %error('No coef > delta'); % uncomment if you want code to error out if no coefficients pass BFDR threshold             
            else
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
        end;
        if length(delt) == 2;
            dStrict = delt(1);
            dLax    = delt(2);
            
            %% Strict FDR %%
            R       = size(Beta,2);
            pst     = NaN(R,1);
            lFDR    = NaN(R,1);
            for r = 1:R;
                Betar       = abs(Beta(:,r));
                pst(r)      = (1/B)*size(find(Betar > dStrict),1);
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
            psiStrict   = pst >= phi;

            %% Lax FDR %%
            R       = size(Beta,2);
            pst     = NaN(R,1);
            lFDR    = NaN(R,1);
            for r = 1:R;
                Betar       = abs(Beta(:,r));
                pst(r)      = (1/B)*size(find(Betar > dLax),1);
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
            psiLax   = pst >= phi;
        end;
%  Below pertains to two sample case
    else
        %% FDR on Diff %%
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
        
        %% FDR on Beta0 %%
        R0      = size(Beta0,2);
        pst0    = NaN(R0,1);
        lFDR0   = NaN(R0,1);
        for r = 1:R0;
            Betar0      = abs(Beta0(:,r));
            pst0(r)     = (1/B)*size(find(Betar0 > delt0),1);
            if pst0(r) == 1;
                pst0(r)  = 1 - (2*B)^(-1);
            end;
            lFDR0(r)     = 1 - pst0(r);
        end;
        if sum(pst0) == 0;
            error('No coef > delta');
        end;
        pr0     = sort(pst0,'descend');
        rstar0  = cumsum(1-pr0)./linspace(1,R0,R0)';
        gam0    = find(rstar0 <= alf/2, 1, 'last' );
        if isempty(gam0);
            phi0    = 1;
        else
            phi0    = pr0(gam0);
        end;
        psi0    = pst0 >= phi0;

        %% FDR on Beta1 %%
        R1      = size(Beta1,2);
        pst1    = NaN(R1,1);
        lFDR1   = NaN(R1,1);
        for r = 1:R1;
            Betar1      = abs(Beta1(:,r));
            pst1(r)     = (1/B)*size(find(Betar1 > delt1),1);
            if pst1(r) == 1;
                pst1(r)  = 1 - (2*B)^(-1);
            end;
            lFDR1(r)     = 1 - pst1(r);
        end;
        if sum(pst1) == 0;
            error('No coef > delta');
        end;
        pr1     = sort(pst1,'descend');
        rstar1  = cumsum(1-pr1)./linspace(1,R1,R1)';
        gam1    = find(rstar1 <= alf/2, 1, 'last' );
        if isempty(gam1);
            phi1    = 1;
        else
            phi1    = pr1(gam1);
        end;
        psi1    = pst1 >= phi1;
    end;
    
    if length(delt) == 1;
        postout.psi     = reshape(psi,Tx,T); 
    end;
    if length(delt) == 2;
        postout.psiStrict   = reshape(psiStrict,Tx,T);
        postout.psiLax      = reshape(psiLax,Tx,T);
    end;
    fprintf('\n Done calculating FDR.\n \n');
    
    %% Then average and find credible interval, update postout
    if twosample == 0;
        Q025_bhat           = quantile(Beta,0.025);
        Q975_bhat           = quantile(Beta,0.975);
        pwCI                = ones(size(Q975_bhat));
        for i = 1:length(pwCI)
           if Q025_bhat(i) < 0 && Q975_bhat(i) > 0
               pwCI(i)      = 0;
           end
        end
        
        postout.pwCI        = reshape(pwCI,Tx,T);
        postout.Q025_bhat   = reshape(Q025_bhat,Tx,T);
        postout.Q975_bhat   = reshape(Q975_bhat,Tx,T);
        postout.bhat        = reshape(mean(Beta),Tx,T);
        postout.Beta        = Beta; % can be useful to retain for running FDR sensitivities, but is a large object      

        postout.Q025_ahat   = quantile(A,0.025);
        postout.Q975_ahat   = quantile(A,0.975);
        postout.ahat        = mean(A);
        
        if nScalarCov > 0 ;
            postout.betahatScalar = reshape(mean(BetaScalar),nScalarCov, T); 
        end;
    else
        %% group 0 %%
        Q025_bhat0          = quantile(Beta0,0.025);
        Q975_bhat0          = quantile(Beta0,0.975);
        pwCI0               = ones(size(Q975_bhat0));
        for i = 1:length(pwCI0)
           if Q025_bhat0(i) < 0 && Q975_bhat0(i) > 0
               pwCI0(i)      = 0;
           end
        end
        
        postout.pwCI0       = reshape(pwCI0,Tx,T);
        postout.Q025_bhat0  = reshape(Q025_bhat0,Tx,T);
        postout.Q975_bhat0  = reshape(Q975_bhat0,Tx,T);
        postout.bhat0       = reshape(mean(Beta0),Tx,T);
        
        %% group 1 %%
        Q025_bhat1          = quantile(Beta1,0.025);
        Q975_bhat1          = quantile(Beta1,0.975);
        pwCI1               = ones(size(Q975_bhat1));
        for i = 1:length(pwCI1)
           if Q025_bhat1(i) < 0 && Q975_bhat1(i) > 0
               pwCI1(i)      = 0;
           end
        end
        
        postout.pwCI1       = reshape(pwCI1,Tx,T);
        postout.Q025_bhat1  = reshape(Q025_bhat1,Tx,T);
        postout.Q975_bhat1  = reshape(Q975_bhat1,Tx,T);
        postout.bhat1       = reshape(mean(Beta1),Tx,T);
        
        %% diff %%
        Q025_diff           = quantile(Diff,0.025);
        Q975_diff           = quantile(Diff,0.975);
        pwCId               = ones(size(Q975_bhatd));
        for i = 1:length(pwCId)
           if Q025_diff(i) < 0 && Q975_diff(i) > 0
               pwCId(i)      = 0;
           end
        end
        
        postout.pwCId       = reshape(pwCId,Tx,T);
        postout.Q025_diff   = reshape(Q025_diff,Tx,T);
        postout.Q975_diff   = reshape(Q975_diff,Tx,T);
        postout.diff        = reshape(mean(Diff),Tx,T);

        %% intercepts %%
        postout.Q025_ahat0  = quantile(A0,0.025);
        postout.Q975_ahat0  = quantile(A0,0.975);
        postout.ahat0       = mean(A0);

        postout.Q025_ahat1  = quantile(A1,0.025);
        postout.Q975_ahat1  = quantile(A1,0.975);
        postout.ahat1       = mean(A1);

        %% group specific psi %%
        postout.psi0        = reshape(psi0,Tx,T);
        postout.psi1        = reshape(psi1,Tx,T);
    end;
    fprintf('\n Done averaging and finding CIs.\n \n');
    

    %% calculate MAPs %%
    if twosample == 0
        [MAPS, upper_CI, lower_CI]  = jointband_maps(Beta,alf);
        
        postout.MAPs                = reshape(MAPS,Tx,T);
        postout.UMAPs               = reshape(upper_CI,Tx,T);
        postout.LMAPs               = reshape(lower_CI,Tx,T);
    else
        %% group 0 %%
        [MAPS0, upper_CI0, lower_CI0]   = jointband_maps(Beta0,alf/2);
        
        postout.MAPs0                   = reshape(MAPS0,Tx,T);
        postout.UMAPs0                  = reshape(upper_CI0,Tx,T);
        postout.LMAPs0                  = reshape(lower_CI0,Tx,T);
        
        %% group 1 %%
        [MAPS1, upper_CI1, lower_CI1]   = jointband_maps(Beta1,alf/2);
        
        postout.MAPs1                   = reshape(MAPS1,Tx,T);
        postout.UMAPs1                  = reshape(upper_CI1,Tx,T);
        postout.LMAPs1                  = reshape(lower_CI1,Tx,T);

        %% diff %%
        [MAPSd, upper_CId, lower_CId]   = jointband_maps(Diff,alf);
        
        postout.MAPsD                   = reshape(MAPSd,Tx,T);
        postout.UMAPsD                  = reshape(upper_CId,Tx,T);
        postout.LMAPsD                  = reshape(lower_CId,Tx,T);

    end
    fprintf('\n Done calculating MAPs.\n \n');

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
        fprintf('\n Done with mixed effects output.\n \n');
    end
    
    %%
    thetahat                = reshape(mean(MCMC_theta),size(theta,1),size(theta,2));
    postout.thetahat        = thetahat;
    postout.theta_025CI     = reshape(quantile(MCMC_theta,0.05),size(theta,1),size(theta,2));
    postout.theta_975CI     = reshape(quantile(MCMC_theta,0.95),size(theta,1),size(theta,2));
    fprintf('\n Done with regularization parameters.\n \n');
    
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

