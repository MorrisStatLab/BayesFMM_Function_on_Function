function [ D, wpspecs ] = dwpt_rows( infile, wpspecs )
%DWPT_ROWS Performs Wavelet Packet transformation on the rows of infile
%          using specifications from wpsecs using hacked Wavelet Packets
%          or hacketts.
%   Detailed explanation goes here
    %% wavelet set up %%
    wavespecs.wavelet   = wpspecs.wavelet;
    wavespecs.wtmode    = wpspecs.wtmode;
    desiredLevel        = wpspecs.nlevels;

    %% set starting values for hacketts %%
    currentLevel        = 0;
    decompCount         = 0;
    wavespecs.nlevels   = 1;
    
    %% set up wavelet tracker %%
    dcAtLevel           = zeros(desiredLevel,1);
    for l = 1:desiredLevel;
        dcAtLevel(l) = 2^(l-1);
    end;
    numbDecomp          = sum(dcAtLevel);
    waveTracker         = cell(numbDecomp,desiredLevel);

    %% run the hacketts %%
    [ D, dc, waveTracker ] = hackett( infile, wavespecs, currentLevel, desiredLevel, decompCount, waveTracker );

    %% re-arrange output from wavelet tracker %%
    %    in a more meaningful way of course    %
    trackerJacker       = cell(numbDecomp,1);
    trackerCount        = 1;
    for d = 1:desiredLevel;
        for n = 1:numbDecomp;
            if isempty(waveTracker{n,d}) == 0;
                trackerJacker{trackerCount}     = waveTracker{n,d};
                trackerCount                    = trackerCount + 1;
            end;
        end;
    end;
    
    %% return results %%
    wpspecs.trackerJacker   = trackerJacker;
    wpspecs.dwtAtLevel      = dcAtLevel;
    wpspecs.totalDWT        = dc;

end

