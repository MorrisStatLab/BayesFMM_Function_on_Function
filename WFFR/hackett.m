function [ MH, decompCount, waveTracker ] = hackett( M, wavespecs, currentLevel, desiredLevel, decompCount, waveTracker )
%HACKETT Generates Wavelet Packets by recursively applying DWT_ROWS to the
%        matrix M untile desired level of decompositon is reached
%   Detailed explanation goes here
    % check level
    if currentLevel == desiredLevel;
        MH              = M;
    else
        %% first decomposition
        [ M, wavespecs1 ]    = dwt_rows(M,wavespecs);
        A   = M(:,1:wavespecs1.Kj(1));
        D   = M(:,(wavespecs1.Kj(1)+1):(sum(wavespecs1.Kj)));
        
        %% update current level of decomposition
        currentLevel        = currentLevel + 1;
%         sprintf('Current Level: %d\n',currentLevel)
        
        decompCount                             = decompCount + 1;
        waveTracker{decompCount, currentLevel}  = wavespecs1;
%         sprintf('Decomp Count: %d\n',decompCount)        
            
        %% call hackett function
        [ A, decompCount, waveTracker ]   = hackett(A, wavespecs1, currentLevel, desiredLevel, decompCount, waveTracker);
        [ D, decompCount, waveTracker ]   = hackett(D, wavespecs1, currentLevel, desiredLevel, decompCount, waveTracker);
        
        % return M
        MH                  = [ A D ];
    end;
end

