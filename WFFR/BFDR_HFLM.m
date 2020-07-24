function [ psi, pst ] = BFDR_HFLM( Beta, V, T, delt, alf, bFlag )
%BFDR_HFLM Compute BFDR in the historical functional setting
%
%   Created: 4/9/2014
%   By: Mark John Meyer

    %% parameters %%
    B       = size(Beta,1); % number of samples

    %% coefficients for inference %%
    if nargin < 6
        bh      = zeros(V,T);
        for i = 1:T;
            for j = i:V
                bh(i,j) = 1;
            end
        end
        bFlag   = reshape(bh,1,V*T);
    end
    bINF    = Beta(:,bFlag == 1);
    
    %% implement FDR %%
    R       = sum(bFlag);
    pst     = NaN(R,1);
    lFDR    = NaN(R,1);
    %%
    for r = 1:R;
        Betar       = abs(bINF(:,r));
        pst(r)      = (1/B)*size(find(Betar > delt),1);
        if pst(r) == 1;
            pst(r)  = 1 - (2*B)^(-1);
        end;
        lFDR(r)     = 1 - pst(r);
    end;
    %%
    if sum(pst) == 0
        psi     = zeros(size(pst));
        fprintf('\n No p(s,t) > delta, psi set to zero \n');
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

        %% insert psi back in to full vector %%
        psif                = zeros(size(bFlag'));
        pstf                = zeros(size(bFlag'));
        psif(bFlag == 1)    = psi;
        pstf(bFlag == 1)    = pst;

        psi                 = reshape(psif,V,T);
        pst                 = reshape(pstf,V,T);
    end

end

