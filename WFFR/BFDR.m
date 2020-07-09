function [ reshape_psi ] = BFDR( Beta, betamax, percSignal, alf, method, nProbes, nDays, outPath )
%BFDR For a specified delta (delt) and alpha (alf) computes Bayesian False
%       Discovery Rate using MCMC samples Beta
%
%   Created: 3/11/2014
%   By: Mark John Meyer
    makeFigures = 0;
    
    delt = betamax * (percSignal/100);
    B       = size(Beta,1); % number of samples
    R       = size(Beta,2); % number of coefficients
    pst     = NaN(R,1);
    lFDR    = NaN(R,1);
    for r = 1:R
        Betar       = abs(Beta(:,r)); % take absolute value of all samples for single coefficient 
        pst(r)      = (1/B)*size(find(Betar > delt),1); % find fraction of samples that abs(coefficient) exceeds delta  
        if pst(r) == 1;  % if all the samples exceed the delta then record a fraction close to 1
            pst(r)  = 1 - (2*B)^(-1); % upper limit of the fraction (will be close to 1) 
        end;
        lFDR(r)     = 1 - pst(r); % inversely related to the fraction of samples that exceed delta
    end
    if sum(pst) == 0  % if no coefficients have samples that exceed the delta
        psi     = zeros(size(pst)); % no coefficients marked as significant
        fprintf('\n No p(s,t) > delta, psi set to zero \n');
    else
        pr      = sort(pst,'descend'); 
        rstar   = cumsum(1-pr)./linspace(1,R,R)'; %linspace generates R points between 1 and R
        gam     = find(rstar <= alf, 1, 'last' );
        if isempty(gam)
            phi = 1;
        else
            phi = pr(gam);
        end
        psi     = pst > phi; 
        reshape_psi     = reshape(psi,nDays,nProbes);
       % reshape_pst    = reshape(pst,nDays,nProbes);
    end    
    if(makeFigures == 1)
        %% make figures (move this section inside the last "end" if you want it
        figure
        colormap(hot)
        x = [1 nDays];
        y = [1 nProbes];
        imagesc(x,y,reshape_psi')
        colorbar
        set(gca, 'Ydir', 'normal')
        set(gca,'FontSize',14)
        xlabel('Gestational Age (days)','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
        ylabel('Genomic Position','FontSize',14) % note this is really "t" but for my pres. I change to site "s"
    %     ax = gca;
    %     ax.XDisplayLabels = num2cell(100:200);
         title(sprintf('BFDR for %s Est. Surface: %d%% of Max. Signal', method, percSignal),'FontSize',14)
        saveas(gcf, sprintf('%s%s_BFDR_%d.png',outPath, method,percSignal))
        
        figure
        colormap(hot)
        x = [1 nDays];
        y = [1 nProbes];
        imagesc(x,y,reshape_pst')
        colorbar
        set(gca, 'Ydir', 'normal')
        set(gca,'FontSize',14)
        xlabel('Gestational Age (days)','FontSize',14) %note this is really "v" but for my presentation I change to time "t"
        ylabel('Genomic Position','FontSize',14) % note this is really "t" but for my pres. I change to site "s"
    %     ax = gca;
    %     ax.XDisplayLabels = num2cell(100:200);
        title(sprintf('Pst for %s Est. Surface: %d%% of Max. Signal', method, percSignal),'FontSize',14)
        saveas(gcf, sprintf('%s%s_Pst_%d.png',outPath, method,percSignal))
    end
end
    

    
    


