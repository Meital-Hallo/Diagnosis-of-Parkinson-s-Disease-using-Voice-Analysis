function [HNR, NHR, SNR]=SpectrumFeatures(frame,frame_pitch,fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ***     HNR, NHR and SNR Per Frame     ****
% Inputs: Current frame of voiced segment
%         Pitch value of the Current frame
%         Fs sampling frequency
% Outputs: HNR
%          NHR
%          SNR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prm=0;
[pxx, freq]=periodogram(frame,[],length(frame),fs,'onesided','power');
Space=round(0.6*frame_pitch/max(diff(freq)));
[pks,loc]=findpeaks(pxx,'MinPeakDistance',Space,'MinPeakHeight',max(pxx)*0.025);
[~, index]=min(abs(freq(loc)-3000)); % Finds indexes of peaks which are placed up to 3000 Hz
Enoise=0;
if index<2 % At least two peaks are required for analysis, Therefore this frame does not count
    HNR = nan;
    NHR = nan;
    SNR = nan;
else
    for j=1:index-1
        ind_noise=find(pxx(loc(j)+1:loc(j+1)-1)<=(mean(pxx(loc(j)+1:loc(j+1)-1)))*1.4);
        if isempty(ind_noise)
            ind_noise=find(pxx(loc(j)+1:loc(j+1)-1)<=(min(pxx(loc(j)+1:loc(j+1)-1)))*10); %Aggravated condition
        end
        if sum(pxx(loc(j)+ind_noise(1):loc(j)+ind_noise(end)))> min(pxx(loc))*8 % In case when part of the harmony is identified as noise
            Enoise=Enoise;
            
        else
            Enoise=Enoise+sum(pxx(loc(j)+ind_noise(1):loc(j)+ind_noise(end)));
        end
        % This plot visualize the separation capability between harmony and noise regions
        %         z=linspace(0,max(pxx),length(freq));
        %         f=ones(1,length(freq));
        %         figure(1);plot(freq(1:100),pxx(1:100),freq(loc),pks,'*',f*freq(loc(j)+ind_noise(1)),z,'g',f*freq(loc(j)+ind_noise(end)),z,'g');
        %         legend('Power Spectrum','Energy peaks','start of noise region',' end of noise region');
        %         xlabel('frequency (Hz)');
        %         ylabel('Power');
        %         pause
    end
    
    ind_start_of_sig=find(pxx(1:loc(1)) <= pxx(loc(1))*0.05); % finds the index of the beginning of the first harmony
    if isempty(ind_start_of_sig)
        ind_start_of_sig=find(pxx(1:loc(1)) <= pxx(loc(1))*0.5);
        if isempty(ind_start_of_sig) % Happens when the first peak is distorted, Therefore this frame does not count
            HNR = nan;
            NHR = nan;
            SNR = nan;
            prm=1;
        end
    end
    ind_end_of_sig=find(diff(pxx(loc(j+1):end)) > 0); % Finds the index of the end of the last harmony
    
    if  Enoise==0
        HNR = nan;
        NHR = nan;
        SNR = nan;
        prm=1;
    end
    
    if prm==0
        % Calculates the power of the harmony parts and the power of the overall spectrum
        Eall=sum(pxx((ind_start_of_sig(end)):((loc(j+1)+ind_end_of_sig(1)-1))));
        Esig=Eall-Enoise;
        
        HNR = 10*log10(Esig/Enoise);
        NHR = 10*log10(Enoise/Eall);
        SNR = 10*log10(Esig/Eall);
        
    end
end
