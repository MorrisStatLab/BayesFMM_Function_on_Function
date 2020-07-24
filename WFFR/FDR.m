function postout = ...
            FDR(MCMC_beta,model,wavespecsc,wavespecsy,pcaspecsx,twosample)
    %% Outputs:
    %  FDR regions and posterior probabilities for one to two surfaces
    %% Inputs
    %  MCMC_beta,model,wavespecsc,wavespecsy,pcaspecsx,twosample from a run
    %  of flmm
    
    %%
    p           = model.p;
    K           = wavespecsy.K;
    Tx          = wavespecsc.T;
    delt        = model.delt;
    alf         = model.alf;
    if twosample == 1;
        delt0   = model.delt0;
        delt1   = model.delt1;
    end;
    
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

    %% Perform IDWT in both directions
    %  Consider both one and two sample case
    %  For twosample, perform inference
    %  twosample defaults to 0
    B       = size(MCMC_beta,1);
    if twosample == 0;
        Beta    = NaN(B,Tx*Tx);
        A       = NaN(B,Tx);
        for i = 1:B;
            Betai   = reshape(MCMC_beta(i,:)',K,p)';
            
            %% Project back into Y space
            As          = Betai(1,:);
            beta_temp   = Betai(2:end,:);
            A1          = idwt_rows(As,wavespecsy);
            yidwt       = idwt_rows(beta_temp,wavespecsy);
            
            %% Remove and update intercept
            A(i,:)      = A1;
                        
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

            if mod(i,100) == 0;
                fprintf('\n Done projecting MCMC samples M = %d \n',i);
            end;       
        end;
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
    
    postout.psi     = reshape(psi,Tx,Tx);
    postout.pst     = reshape(pst,Tx,Tx);
    if twosample == 1;
        postout.psi0        = reshape(psi0,Tx,Tx);
        postout.pst0        = reshape(pst0,Tx,Tx);
        postout.psi1        = reshape(psi1,Tx,Tx);
        postout.pst1        = reshape(pst1,Tx,Tx);
    end;
    fprintf('\n Done calculating FDR.\n \n');
    
%     %% Calculate 2-sided p-value %%
%     if twosample == 0;
%         R       = size(Beta,2);
%         plt0    = zeros(R,1);
%         pgt0    = zeros(R,1);
%         p0      = zeros(R,1);
%         for r = 1:R;
%             Betar       = Beta(:,r);
%             plt0(r)     = (1/B)*size(find(Betar < 0),1);
%             pgt0(r)     = (1/B)*size(find(Betar > 0),1);
%             p0(r)       = 2*min(plt0(r),pgt0(r));
%         end;
%     else
%         %% Difference Bayesian p-value %%
%         R       = size(Diff,2);
%         plt0d   = zeros(R,1);
%         pgt0d   = zeros(R,1);
%         p0d     = zeros(R,1);
%         for r = 1:R;
%             Diffr       = Diff(:,r);
%             plt0d(r)    = (1/B)*size(find(Diffr < 0),1);
%             pgt0d(r)    = (1/B)*size(find(Diffr > 0),1);
%             p0d(r)      = 2*min(plt0d(r),pgt0d(r));
%         end;
%         
%         %% Beta0 Bayesian p-value %%
%         R0      = size(Beta0,2);
%         plt00   = zeros(R0,1);
%         pgt00   = zeros(R0,1);
%         p00     = zeros(R0,1);
%         for r = 1:R0;
%             Betar0      = Beta0(:,r);
%             plt00(r)    = (1/B)*size(find(Betar0 < 0),1);
%             pgt00(r)    = (1/B)*size(find(Betar0 > 0),1);
%             p00(r)      = 2*min(plt00(r),pgt00(r));
%         end;
% 
%         %% Beta1 Bayesian p-value %%
%         R1      = size(Beta1,2);
%         plt01   = zeros(R1,1);
%         pgt01   = zeros(R1,1);
%         p01     = zeros(R1,1);
%         for r = 1:R1;
%             Betar1      = Beta1(:,r);
%             plt01(r)    = (1/B)*size(find(Betar1 < 0),1);
%             pgt01(r)    = (1/B)*size(find(Betar1 > 0),1);
%             p01(r)      = 2*min(plt01(r),pgt01(r));
%         end;
%     end;
%     
%     %% update output for bayesian p-value %%
%     if twosample == 0;
%         postout.p0      = p0;
%     else
%         postout.p0d     = reshape(p0d,Tx,Tx);
%         postout.p00     = reshape(p00,Tx,Tx);
%         postout.p01     = reshape(p01,Tx,Tx);
%     end;
    
    %% compute pointwise CI %%
    
    if twosample == 0;
        Q025_bhat           = quantile(Beta,0.025);
        Q975_bhat           = quantile(Beta,0.975);
        pwCI                = ones(size(Q975_bhat));
        for i = 1:length(pwCI)
           if Q025_bhat(i) < 0 && Q975_bhat(i) > 0
               pwCI(i)      = 0;
           end
        end
        
        postout.pwCI        = reshape(pwCI,Tx,Tx);
        postout.Q025_bhat   = reshape(Q025_bhat,Tx,Tx);
        postout.Q975_bhat   = reshape(Q975_bhat,Tx,Tx);
        postout.bhat        = reshape(mean(Beta),Tx,Tx);
        

        postout.Q025_ahat   = quantile(A,0.025);
        postout.Q975_ahat   = quantile(A,0.975);
        postout.ahat        = mean(A);
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
        
        postout.pwCI0       = reshape(pwCI0,Tx,Tx);
        postout.Q025_bhat0  = reshape(Q025_bhat0,Tx,Tx);
        postout.Q975_bhat0  = reshape(Q975_bhat0,Tx,Tx);
        postout.bhat0       = reshape(mean(Beta0),Tx,Tx);
        
        %% group 1 %%
        Q025_bhat1          = quantile(Beta1,0.025);
        Q975_bhat1          = quantile(Beta1,0.975);
        pwCI1               = ones(size(Q975_bhat1));
        for i = 1:length(pwCI1)
           if Q025_bhat1(i) < 0 && Q975_bhat1(i) > 0
               pwCI1(i)      = 0;
           end
        end
        
        postout.pwCI1       = reshape(pwCI1,Tx,Tx);
        postout.Q025_bhat1  = reshape(Q025_bhat1,Tx,Tx);
        postout.Q975_bhat1  = reshape(Q975_bhat1,Tx,Tx);
        postout.bhat1       = reshape(mean(Beta1),Tx,Tx);
        
        %% diff %%
        Q025_diff           = quantile(Diff,0.025);
        Q975_diff           = quantile(Diff,0.975);
        pwCId               = ones(size(Q975_diff));
        for i = 1:length(pwCId)
           if Q025_diff(i) < 0 && Q975_diff(i) > 0
               pwCId(i)      = 0;
           end
        end
        
        postout.pwCId       = reshape(pwCId,Tx,Tx);
        postout.Q025_diff   = reshape(Q025_diff,Tx,Tx);
        postout.Q975_diff   = reshape(Q975_diff,Tx,Tx);
        postout.diff        = reshape(mean(Diff),Tx,Tx);

        %% intercepts %%
        postout.Q025_ahat0  = quantile(A0,0.025);
        postout.Q975_ahat0  = quantile(A0,0.975);
        postout.ahat0       = mean(A0);

        postout.Q025_ahat1  = quantile(A1,0.025);
        postout.Q975_ahat1  = quantile(A1,0.975);
        postout.ahat1       = mean(A1);

    end;
    fprintf('\n Done averaging and finding CIs.\n \n');
    
    
    %% calculate MAPs %%
    if twosample == 0
        [MAPS, upper_CI, lower_CI]  = jointband_maps(Beta,alf);
        
        postout.MAPs                = reshape(MAPS,Tx,Tx);
        postout.UMAPs               = reshape(upper_CI,Tx,Tx);
        postout.LMAPs               = reshape(lower_CI,Tx,Tx);
    else
        %% group 0 %%
        [MAPS0, upper_CI0, lower_CI0]   = jointband_maps(Beta0,alf/2);
        
        postout.MAPs0                   = reshape(MAPS0,Tx,Tx);
        postout.UMAPs0                  = reshape(upper_CI0,Tx,Tx);
        postout.LMAPs0                  = reshape(lower_CI0,Tx,Tx);
        
        %% group 1 %%
        [MAPS1, upper_CI1, lower_CI1]   = jointband_maps(Beta1,alf/2);
        
        postout.MAPs1                   = reshape(MAPS1,Tx,Tx);
        postout.UMAPs1                  = reshape(upper_CI1,Tx,Tx);
        postout.LMAPs1                  = reshape(lower_CI1,Tx,Tx);

        %% diff %%
        [MAPSd, upper_CId, lower_CId]   = jointband_maps(Diff,alf);
        
        postout.MAPsD                   = reshape(MAPSd,Tx,Tx);
        postout.UMAPsD                  = reshape(upper_CId,Tx,Tx);
        postout.LMAPsD                  = reshape(lower_CId,Tx,Tx);

    end

    

