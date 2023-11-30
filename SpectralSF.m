function  [PitchMaxPeaksRatio, ChangeBtw_MaxAndPitch_Power, MaxPeak_Above_Third_Harmony, PitchPeak_Lower_than_Next2Peaks, PitchPeak_20PreLower_than_Next3Peaks]=SpectralSF(frame,fs,frame_pitch)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ***   Special Frames Conditions Per Frame  ***
% Input:  Current frame of voiced segment 
%         Fs sampling frequency
%         Pitch value of the Current frame 
% Output: 1. Ratio between Pitch Peak and Max Peak
%         the following conditions return 1 for true condition and 0 for false condition:
%         2. 50% Change Between Max and Pitch Power
%         3. Max Peak above Third harmony
%         4. Pitch Peak Lower than at least Next 2 Peaks
%         5. Pitch Peak Lower in 20% than at least Next 3 Peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pxx, freq]=periodogram(frame,[],length(frame),fs,'onesided','power');
Space=round(0.6*frame_pitch/max(diff(freq)));
[pks,loc]=findpeaks(pxx,'MinPeakDistance',Space,'MinPeakHeight',max(pxx)*0.025);
[~, index]=min(abs(freq(loc)-3000));
if index<2 % At least two peaks are required for analysis, Therefore this frame does not count 
    ChangeBtw_MaxAndPitch_Power = nan;
    MaxPeak_Above_Third_Harmony = nan;
    PitchPeak_Lower_than_Next2Peaks = nan;
    PitchPeak_20PreLower_than_Next3Peaks = nan;
    PitchMaxPeaksRatio = nan;
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 50% Change Between Max and Pitch Power
    [~,ind]=min(abs(frame_pitch-freq(loc)));
    [maxPEAK,maxLOC]=max(pxx(loc));
    if (maxPEAK~=pxx(loc(ind)) && pxx(loc(ind))<0.50*maxPEAK)
        ChangeBtw_MaxAndPitch_Power=1;
    else
        ChangeBtw_MaxAndPitch_Power=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Max Peak above Third harmony
    if (ceil(freq(loc(max_ind))/freq(loc(ind)))>3)
        MaxPeak_Above_Third_Harmony=1;
    else
        MaxPeak_Above_Third_Harmony=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pitch Peak Lower than at least Next 2 Peaks
    count=0;
    if max_ind-ind>1
        for i=1:max_ind-ind-1
            if pxx(loc(ind))<pxx(loc(ind+i));
                count=count+1;
            end
        end
        if count>=2
            PitchPeak_Lower_than_Next2Peaks=1;
        else
            PitchPeak_Lower_than_Next2Peaks=0;
        end
    else
        PitchPeak_Lower_than_Next2Peaks=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pitch Peak Lower in 20% than at least Next 3 Peaks
    count=0;
    if max_ind-ind>1
        for i=1:max_ind-ind-1
            if 0.8*pxx(loc(ind))<pxx(loc(ind+i))
                count=count+1;
            end
        end
        if count>=3
            PitchPeak_20PreLower_than_Next3Peaks=1;
        else
            PitchPeak_20PreLower_than_Next3Peaks=0;
        end
    else
        PitchPeak_20PreLower_than_Next3Peaks=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Ratio between Pitch Peak and Max Peak
    PitchMaxPeaksRatio=pxx(ind)/maxPEAK;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot of the spectrum for validation
    %         plot(freq(1:100),pxx(1:100),freq(loc),pks,'*')
    %         pause
end