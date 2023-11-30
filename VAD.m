function [Speech, Result, Voiced, UnVoiced, AllSilence] = VAD(Record ,fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This VAD detects Voiced Unvoiced and Silence segments in speech
% and removes silence at the start and end.
% Outputs: Speech = vector without the slience at the start and end.
%          Result - a vector with the classification of the frame - UV=0.5 V=1 S=0
%          Voiced - a vector with only voiced segments.
%          UnVoiced - a vector with only UnVoiced segments.
%          AllSilence - a vector with only Silence segments of the entire record.
%
%    Note: The size of the frame most remain 0.02msec, else need to chage
%    the thresholds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

window_size_sampels=0.02*fs;
Record=vec2frame(Record, window_size_sampels, 0, 0);
[m,n_frames]=size(Record); speech=Record(:); Record=Record./max(speech); 

zcd = dsp.ZeroCrossingDetector;  
for i=1:n_frames
    frame=Record(:,i);
    
    %Error of LPC
    [~,LPC_ERR(i)]=lpc(frame,16); 
    
    %ZeroCrossing
    ZC(i)=double(zcd(frame)); 
    
    %R1 of autocorrelation
    [acf,lags] =xcorr(frame, 'none'); acf = acf(lags>=0);
    R1(i)=abs(acf(2));
    
    %STFT - spectrogram
    [s,~, ~] =spectrogram(frame, 'yaxis');
    MAX_STFT(i)=max(abs(s(:)));
    
    %Pitch
    [pitch(i), ~] = XcorrPitch(frame, fs);
end

%normalization
RMS=rms(Record); 
R1=R1/max(R1); R1=medfilt1(R1);
LPC_ERR=LPC_ERR/max(LPC_ERR); LPC_ERR=medfilt1(LPC_ERR);
MAX_STFT=MAX_STFT/max(MAX_STFT) ;MAX_STFT=medfilt1(MAX_STFT,5);
ZC=medfilt1(ZC);
AZC=ZC;

%Threshold values
ZCV=(max(ZC)-min(ZC))*0.5+min(ZC);
R1V=(max(R1)-min(R1))*0.01+min(R1);
eV=(max(LPC_ERR)-min(LPC_ERR))*0.005+min(LPC_ERR);
MAXV=(max(MAX_STFT)-min(MAX_STFT))*0.01+min(MAX_STFT);
RMSV=(max(RMS)-min(RMS))*0.1+min(RMS);

if (max(AZC)<100 && mean(RMS)>0.1) %if it has only ha without silence
    ZCV=ZCV/max(ZC); ZC=ZC/max(ZC); RMS=RMS/max(RMS); RMSV=RMSV/max(RMS);
    for i=1:n_frames
        if (MAX_STFT(i)>0.9*MAXV && R1(i)>0.9*R1V && RMS(i)>0.3*RMSV && (pitch(i)>60 && pitch(i)<700))
            Result(i)=1; %Voiced
        elseif (ZC(i)>0.5*ZCV && MAX_STFT(i)>0.9*MAXV && LPC_ERR(i)>0.5*eV && R1(i)>0.1*R1V && pitch(i)<700)
            Result(i)=0.5; %Unvoiced
        else
            Result(i)=0; %Silence
        end
    end
else
    ZCV=ZCV/max(ZC); ZC=ZC/max(ZC); RMS=RMS/max(RMS); RMSV=RMSV/max(RMS);
    for i=1:n_frames
        if  (ZC(i)<1.1*ZCV) && MAX_STFT(i)>MAXV && R1(i)>0.9*R1V && RMS(i)>1.2*RMSV && (pitch(i)>60 && pitch(i)<700)
            Result(i)=1; %Voiced
        elseif (ZC(i)>0.2*ZCV && MAX_STFT(i)>0.8*MAXV && RMS(i)>1.2*RMSV && LPC_ERR(i)>0.3*eV && R1(i)>0.02*R1V && (pitch(i)<60 || pitch(i)>700)) %0.4*ZCV
            Result(i)=0.5; %Unvoiced
        else
            Result(i)=0; %Silence
        end
    end  
end
Result=medfilt1(Result, 5);

% figure;
% window_size_sec= window_size_sampels/fs;
% t=(0:(m*n_frames-1))/fs;
% tt=0.5*window_size_sec:window_size_sec:(m*n_frames-0.5*window_size_sec)/fs;
% plot(t, Record(:), tt,LPC_ERR/max(LPC_ERR), tt, ZC/max(ZC),  tt, MAX_STFT/max(MAX_STFT),tt, R1/max(R1), tt, RMS ,tt, pitch/max(pitch) ,tt, Result, 'r*-'); 
% legend('Speech', 'LPC ERROR', 'Zero Crossing', 'STFT', 'R1', 'RMS', 'pitch' ,'VAD Result'); ylim([ min(Record(:)), max(Record(:))]); xlim([0 tt(end)]);
% xlabel('Time (sec)')



V_UV_ind=find(Result==0.5 | Result==1);
START_POINT=V_UV_ind(1);
END_POINT=V_UV_ind(end);
Result=Result(START_POINT:END_POINT);
Record=Record(:,START_POINT:END_POINT);
Voiced=Record(:, Result==1); Voiced=(Voiced(:));
UnVoiced=Record(:, Result==0.5); UnVoiced=(UnVoiced(:));
AllSilence=(Record(:, Result==0)); AllSilence=(AllSilence(:));
Speech=Record(:);


end

% res=A(:,(Result>0));
% plot( Speech(:)) 
% soundsc(Speech(:), fs)
