function  [AllPksSlope, PitchLastPksSlope, MaxMinPksSlope, maxPEAK, PitchPeak]=SpectralSlope(frame,fs,frame_pitch)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ***     Spectral Slope Per Frame    ****
% Inputs: Current frame of voiced segment 
%         Fs sampling frequency
%         Pitch value of the Current frame 
% Outputs: The slope between all peaks in the frame
%          The slope between the peak that correlate to the pitch to the last peak in the frame.
%          The slope between the maximum peak to the minimum peak 
%          The power of the maximum peak 
%          The power of the peak that correlate to the pitch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pxx, freq]=periodogram(frame,[],length(frame),fs,'onesided','power');
Space=round(0.6*frame_pitch/max(diff(freq))); 
[pks, loc]=findpeaks(pxx,'MinPeakDistance',Space,'MinPeakHeight',max(pxx)*0.025);
[~, index]=min(abs(freq(loc)-3000));

if index<2 % At least two peaks are required for analysis, Therefore this frame does not count
    AllPksSlope = nan;
    PitchLastPksSlope = nan;
    MaxMinPksSlope = nan;
    maxPEAK = nan;
    PitchPeak = nan;
else 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % slope of max to min peaks
    [maxPEAK,max_ind]=max(pxx(loc));
    [minPEAK,~]=min(pxx(loc(max_ind:end)));
    min_ind=find(pxx(loc)==minPEAK);
    p_max=polyfit([freq(loc(max_ind)),freq(loc(min_ind))],[maxPEAK,minPEAK],1);
    MaxMinPksSlope=p_max(1); % slope
    % Plot of the spectrum for validation
%     plot(freq(1:100),pxx(1:100),freq(loc),pks,'*',freq(loc(max_ind)),maxPEAK,'o',freq(loc(min_ind)),minPEAK,'o')
%     pause 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % slope of all peaks
    p_all=polyfit(freq(loc),pxx(loc),1);
    AllPksSlope=p_all(1); % slope
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pitch peaks
    [~,ind]=min(abs(frame_pitch-freq(loc))); 
    PitchPeak=pxx(loc(ind));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % slope of pitch to last peaks 
    p_last=polyfit([freq(loc(ind)),freq(end)],[pxx(loc(ind)),pxx(end)],1);
    PitchLastPksSlope=p_last(1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
end