function [pitch, AutoPeakAmp] = PitchPerFrame(A ,fs)

[m,n]=size(A); i_run=1:n;

for i=i_run
    frame=A(:,i);
    [pitch(i), AutoPeakAmp(i)] = XcorrPitchVoice (frame, fs);
     %plot((0:m-1)/fs, A(:,i)); xlabel('Time [sec]') %title(['Frame=' num2str(i), '  F0=' num2str(pitch(i))]);
end
%pitch=medfilt1(pitch);

% to correct wrong pitch
% for i=i_run 
%     if (pitch(i)>0)
%         frame=A(:,i);
%         [~, ~ ,val]= FindPeaks_Func2(frame, fs, pitch(i), 0.8 ,0.4, i);
%         if (val==1 && pitch(i)>145 && AutoPeakAmp(i)<0.75)
%             pitch(i)=(pitch(i)/2);
%         end
%     end
% end


%figure(2); plot(pitch, 'r'); xlabel('Number Of Frame'); ylabel('Frequency (Hz)'); xlim([0, length(pitch)]); title(['Pitch Variation Through Sustain Vowel - Median Pitch: ',  num2str(median(pitch))]);
end


function [pitch, AutoAmp] = XcorrPitchVoice (frame, fs)
% finding fundamental frequency of signal, A could be matrix or vector

MinPk_dis_samples=0.007*fs;
[c, lags] = xcorr(frame, 'coeff');

% Take only non-negative lag values and Convert lag to seconds.
c = c(lags>=0);  lags = lags(lags>=0);  lags_s = lags/fs;
c(find(lags_s>0.015))=0; %delete any freq low 50 Hz
c(find(lags_s<0.0025))=0; %delete any freq above 400 Hz

[pks,lcs] = findpeaks(c, 'MinPeakDistance', MinPk_dis_samples, 'MinPeakHeight',0.4*max(c)); 

if ~isempty(pks)
    [~, ind_max]=max(pks); % Finding location of maximum amplitude frequency
    MaxPeakLag=lags_s(lcs(ind_max))/ind_max; %taking max amplitude location divided by its index
    FirstPeakLag=lags_s(lcs(1)); %taking the first max freq
    
    if (MaxPeakLag>0.95*FirstPeakLag && MaxPeakLag<1.05*FirstPeakLag) %if point is in range, it's considered mulitple of lower freq
        pitch=1/FirstPeakLag; AutoAmp=c(lcs(1));
        if ind_max==2
            keyboard;
        end
    else
        pitch=1/(MaxPeakLag*ind_max); AutoAmp=c(lcs(ind_max));
    end
    
else %if pks is empty - insert freq of 0
    pitch=0; AutoAmp=0; 
end

figure(2); 
subplot(211); plot((0:1:length(frame)-1)/fs, frame)
subplot(212); plot(lags_s ,c, lcs/fs , pks,'r*');

end



