function [Abnorm2T0, Abnorm2Ti, LowAutoCorr, LowNumPeaks, TotPeaksRatio, TotNumSpecFrames, RatioNumSpecFrames] = SpecialSeg(A, fs, window_size_sec, F0, AutoPeakAmp, PitchForFrame, PATHO, SPKR_NUM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ***  Counter of Special Frames   ***
% Input:  A record after VAD and division to frames
%         Fs sampling frequency
%         Size of frame in sec
%         F0 median Pitch of the entire record
%         AutoPeakAmp is the amplitude of autocorrelation
%         PATHO is the label of patology (3=HC, 1=PD)
%         SPKR_NUM id the speaker number
% Output: Number of frames with abnormal periods comparing to Pitch
%         Number of frames with abnormal periods comparing to Pitch per frame
%         Number of frames with low autocorrelation amplitude
%         Number of frames with low number peaks
%         Detect abnormal number of peaks
%         Total Number of abnormal frames
%         Ratio between Total Number of abnormal frames to the total number of frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,n_frames]=size(A); 
Abnorm2T0=0; Abnorm2Ti=0; LowAutoCorr=0; LowNumPeaks=0;
%SpecFrameIndex=zeros(1,n);
%ErrName=[];

for i=1:n_frames
    if PitchForFrame(i)>0
        Data= A(:,i);
        [~, ind] = FindPeaks_Func2(Data, fs, PitchForFrame(i), 0.75 ,0.35, i);
        
        if length(ind)>1
            periods=(diff(ind));
            
            %abnormal T(i) comparing to F0
            AbnormalPeriodsF0=length(find(periods>1.3/F0 | periods<0.7/F0 )); 
            if (AbnormalPeriodsF0~=0)
                Abnorm2T0=Abnorm2T0+1;
                %SpecFrameIndex(i)=1;
                %ErrName=[ErrName, ' Abnormal T0'];
            end
            
            %abnormal 2*T(i) comparing to F(i)
            AbnormalPeriodsFi=length(find(periods>1.3/PitchForFrame(i) | periods<0.7/PitchForFrame(i) ));
            if (AbnormalPeriodsFi~=0)
                Abnorm2Ti=Abnorm2Ti+1;
                %SpecFrameIndex(i)=1;
                %ErrName=[ErrName, ' Abnormal T(frame)'];
            end
            
            %low amplitude in autocorrelation
            if (AutoPeakAmp(i)<0.5) 
                LowAutoCorr=LowAutoCorr+1;
                %SpecFrameIndex(i)=1;
                %ErrName=[ErrName, ' Low AutoCorr ' num2str(AutoPeakAmp(i))];
            end
            
            %low number of peaks Ttot/T0(i)
            if  (PitchForFrame(i)>80 && length(periods)<floor((ind(end)-ind(1))*PitchForFrame(i)) |  sum(periods)<0.4*window_size_sec)
                LowNumPeaks=LowNumPeaks+1;
                %SpecFrameIndex(i)=1;
                %ErrName=[ErrName, ' Abnormal Number of Peaks'];
            end
            
            % Ratio between the number of periods that is suppose to be
            PeaksRatio(i)=length(periods)/round((ind(end)-ind(1))*PitchForFrame(i));
        else
            LowNumPeaks=LowNumPeaks+1;
            %SpecFrameIndex(i)=1;
            %ErrName=[ErrName, ' 0-1 Peak'];
        end
        
        
    else 
    Abnorm2Ti=nan; 
    TotPeaksRatio=nan;
    Abnorm2T0=nan;
    LowAutoCorr=nan;
    LowNumPeaks=nan;
 
    end
end

%Abnorm2Ti=Abnorm2Ti/n_frames;
TotNumSpecFrames=Abnorm2T0+Abnorm2Ti+LowAutoCorr+LowNumPeaks;
RatioNumSpecFrames=TotNumSpecFrames/n_frames;   
TotPeaksRatio=entropy(PeaksRatio);

% Option To Save figures of abnormal frames automatically.    
    %m=length(Data);
    %t=(0:1:m-1)/fs;
    %if  SpecFrameIndex(i)==1 & PATHO=='HC     '
    %    filename=[PATHO, ' Speaker ', num2str(SPKR_NUM),'-' f_name, ' Frame-', num2str(i),' ~~' ErrName]
    %    saveas(gcf,['C:\Users\Meital Hallo\Desktop\fig\HC\' ,filename, '.fig'])
    %end 
    %if  SpecFrameIndex(i)==1 & PATHO=='PD     '
    %    filename=[PATHO, ' Speaker ', num2str(SPKR_NUM),'-' f_name, ' Frame-', num2str(i),' ~~' ErrName]
    %    saveas(gcf,['C:\Users\Meital Hallo\Desktop\fig\PD\' ,filename, '.fig'])
    %end
    %ErrName=[];   
end

