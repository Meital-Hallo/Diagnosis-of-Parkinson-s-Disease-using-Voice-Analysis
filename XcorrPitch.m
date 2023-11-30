function [pitch, AutoAmp] = XcorrPitch (frame, fs)
% finding fundamental frequency of signal

MinPk_dis_samples=0.006*fs;
[c, lags] = xcorr(frame, 'coeff');

% Take only non-negative lag values and Convert lag to seconds.
c = c(lags>=0);  lags = lags(lags>=0);  lags_s = lags/fs;
c(find(lags_s>0.015))=0; %delete any freq low 50 Hz

[pks,lcs] = findpeaks(c, 'MinPeakDistance', MinPk_dis_samples, 'MinPeakHeight',0.1*max(c)); 
%plot(lags_s ,c, lcs/fs , pks,'r*')

if ~isempty(pks)
    [~, ind_max]=max(pks); % Finding location of maximum amplitude frequency
    MaxPeakLag=lags_s(lcs(ind_max))/ind_max; %taking max amplitude location divided by its index
    FirstPeakLag=lags_s(lcs(1)); %taking the first max freq
    
    if (MaxPeakLag>0.95*FirstPeakLag && MaxPeakLag<1.05*FirstPeakLag) %if point is in range, it's considered mulitple of lower freq
        pitch=1/FirstPeakLag; AutoAmp=c(lcs(1));
    else
        pitch=1/(MaxPeakLag*ind_max); AutoAmp=c(lcs(ind_max));
    end
    
else %if pks is empty - insert freq of 0
    pitch=0; AutoAmp=0; 
end

%figure(2); plot(lags_s ,c, lcs/fs , pks,'r*'); title(['F0= ' num2str(pitch)]);

end