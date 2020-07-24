function [ res ] = idwpt_rows( infile, wpspecs )
%IDWPT_ROWS Compute IDWPT on rows of infile using IDWT_ROWS
%   Detailed explanation goes here

    %% unlist elements of wpspecs
    nlevels         = wpspecs.nlevels;
    dwtAtLevel      = wpspecs.dwtAtLevel;
    trackerJacker   = wpspecs.trackerJacker;
    nrows           = size(infile,1);

    %% reconstruct signal by rows
    for i = 1:nlevels
        %%
        f           = (nlevels-i)+1;
        w           = dwtAtLevel(f);
        k           = w;
        workinf     = zeros(nrows,w*trackerJacker{w}.T);
        %%
        for j = 1:w
            %% extract appropriate wavelet specs and wavelet coefs
            wavespecs   = trackerJacker{k};
            K           = wavespecs.K;
            T           = wavespecs.T;
            
            %% set indices for IDWT
            fromK       = K*(j-1)+1;
            toK         = K*j;
            indK        = fromK:toK;
            
            %% set indices for update
            fromT       = T*(j-1)+1;
            toT         = T*j;
            indT        = fromT:toT;
            
            %% perform IDWT on correct subset
            infw        = infile(:,indK);
            Y           = idwt_rows(infw,wavespecs);

            %% update coef matrix and counters
            k                   = k + 1;
            workinf(:,indT)     = Y;
        end
            %% update infile to resize
            %  resizing will occur if non-Haar wavelets are used
            %  if Haar are used, and the signal is dyadic
            %  dimmensions will remain the same
            infile  = workinf;
    end
    
    %%
    res     = infile;
end

