function [FtrVec, FtrName]=CalcFeatures(FtrIndex,SpeechDat, Fs, Pathology, SPKR_ID)
% Calcualte features from speech file

global T_frame;
global t_overlap;

% General Parameters
Tw=T_frame;
Ts=t_overlap;

% Normalize the signal
if (FtrIndex~=3 && FtrIndex~=9) %VAD needs the record as it is.
    SpeechDat=SpeechDat-mean(SpeechDat);
    SpeechDat=SpeechDat/rms(SpeechDat);
end

% features selection
switch  (FtrIndex)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NHR , HNR , SNR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        disp(' ***  NHR, HNR, SNR ***')
        clear Ceoff;
        [  frames ] =shir_vec2frames(SpeechDat,Fs,Tw,Ts);
        [ ~, n_frames]=size(frames);
        for i_frame=1:n_frames
            curr_frame=frames(:,i_frame);
            frame_after_threshold=threshold(curr_frame);
            FramePitch(i_frame)=pitch_period(frame_after_threshold,Fs);
            [HNR(i_frame), NHR(i_frame), SNR(i_frame)]=SpectrumFeatures(curr_frame,FramePitch(i_frame),Fs);
        end
        %deleting nan values
        HNR(isnan(HNR))=[];
        NHR(isnan(NHR))=[];
        SNR(isnan(SNR))=[];
        
        FtrVec(1,:)=HNR;
        FtrName(1,:)=['HNR' num2str(1)];
        
        FtrVec(2,:)=NHR;
        FtrName(2,:)=['NHR' num2str(2)];
        
        FtrVec(3,:)=SNR;
        FtrName(3,:)=['SNR' num2str(3)];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Spectral Slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
    case 2
        disp(' ***  Spectral Slope  ***')
        clear Ceoff;
        [  frames ] =shir_vec2frames(SpeechDat,Fs,Tw,Ts);
        [ ~, n_frames]=size(frames);
        for i_frame=1:n_frames
            curr_frame=frames(:,i_frame);
            frame_after_threshold=threshold(curr_frame);
            FramePitch(i_frame)=pitch_period(frame_after_threshold,Fs);
            [AllPksSlope(i_frame), PitchLastPksSlope(i_frame), MaxMinPksSlope(i_frame),MaxPeak(i_frame),PitchPeak(i_frame)]=SpectralSlope(curr_frame,Fs,FramePitch(i_frame));
        end
        % deleting nan values
        AllPksSlope(isnan(AllPksSlope))=[];
        PitchLastPksSlope(isnan(PitchLastPksSlope))=[];
        MaxMinPksSlope(isnan(MaxMinPksSlope))=[];
        MaxPeak(isnan(MaxPeak))=[];
        PitchPeak(isnan(PitchPeak))=[];
        
        FtrVec(1,:)=AllPksSlope;
        FtrName(1,:)=['AllPksSlope      ' num2str(1)];
         
        FtrVec(2,:)=PitchLastPksSlope;
        FtrName(2,:)=['PitchLastPksSlope' num2str(2)]; 
        
        FtrVec(3,:)=MaxMinPksSlope;
        FtrName(3,:)=['MaxMinPksSlope   ' num2str(3)];
        
        % slope of Max Pks
        p_max = polyfit((1:length(MaxPeak))*Tw,MaxPeak,1);
        FtrVec(4,1)=p_max(1); %slope
        FtrName(4,:)=['MaxPksSlope      ' num2str(4)];

        % slope of Pitch Pks
        p_pitch = polyfit((1:length(PitchPeak))*Tw,PitchPeak,1);
        FtrVec(5,1)=p_pitch(1); %slope
        FtrName(5,:)=['PitchPksSlope    ' num2str(5)];                            
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Speech Rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
   
    case 3 
        disp(' *** Speech Rate ***');
        clear Ceoff;
        [SpeechRate_V_UV,SpeechRate_V_ALL,SpeechRate_UV_ALL]=SpeechRate(SpeechDat,Fs);
        
        FtrVec(1,1)=SpeechRate_V_UV;
        FtrName(1,:)=['SpeechRate_V_UV   ' num2str(1)];
        
        FtrVec(2,1)=SpeechRate_V_ALL;
        FtrName(2,:)=['SpeechRate_V_ALL  ' num2str(2)];
       
        FtrVec(3,1)=SpeechRate_UV_ALL;
        FtrName(3,:)=['SpeechRate_UV_ALL ' num2str(3)];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Jitter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    case 4  
        disp(' ***   Jitter ***');
        clear Ceoff;
        frames=vec2frame(SpeechDat, Tw*Fs, Ts*Tw*Fs, 0);
        [PitchForFrame, AutoPeakAmp]= PitchPerFrame(frames, Fs);
        Pitch=median(PitchForFrame); [~, n_frames]=size(frames);
        [jitter_absolute,jitter_relative]=Jitter(frames(:,ceil(0.05*n_frames):n_frames*0.85), Fs, Pitch, AutoPeakAmp, PitchForFrame);
        
        FtrVec(1,:)=jitter_absolute;
        FtrName(1,:)=['jitter_absolute' num2str(1)];
        
        FtrVec(2,:)=jitter_relative;
        FtrName(2,:)=['jitter_relative' num2str(2)];
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Shimmer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    case 5
        disp(' ***   Shimmer ***');
        clear Ceoff;
        frames=vec2frame(SpeechDat, Tw*Fs, Ts*Tw*Fs, 0);
        [PitchForFrame, AutoPeakAmp]= PitchPerFrame(frames, Fs); Pitch=median(PitchForFrame); [~, n_frames]=size(frames);
        [shimmer_dB, shimmer_abs]=Shimmer(frames(:,ceil(0.05*n_frames):n_frames*0.85), Fs, Pitch, AutoPeakAmp, PitchForFrame);
        
        FtrVec(1,:)=shimmer_dB;
        FtrName(1,:)=['shimmer_db ' num2str(1)];
        
        FtrVec(2,:)=shimmer_abs;
        FtrName(2,:)=['shimmer_abs' num2str(2)];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Special Frames - Frequency Domain 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
     case 6
        disp(' ***  Special frames Freq  ***')
        clear Ceoff;
        [  frames ] =shir_vec2frames(SpeechDat,Fs,Tw,Ts);
        [ ~, n_frames]=size(frames);
        for i_frame=1:n_frames
            curr_frame=frames(:,i_frame);
            frame_after_threshold=threshold(curr_frame);
            FramePitch(i_frame)=pitch_period(frame_after_threshold,Fs);
            [ChangeBtw_MaxAndPitch_Power(i_frame), MaxPeak_Above_Third_Harmony(i_frame), PitchPeak_Lower_than_Next2Peaks(i_frame), PitchPeak_20PreLower_than_Next3Peaks(i_frame), PitchMaxPeaksRatio(i_frame)]=SpectralSF(curr_frame,Fs,FramePitch);
        end
        % deleting nan values
        ChangeBtw_MaxAndPitch_Power(isnan(ChangeBtw_MaxAndPitch_Power))=[];
        MaxPeak_Above_Third_Harmony(isnan(MaxPeak_Above_Third_Harmony))=[];
        PitchPeak_Lower_than_Next2Peaks(isnan(PitchPeak_Lower_than_Next2Peaks))=[];
        PitchPeak_20PreLower_than_Next3Peaks(isnan(PitchPeak_20PreLower_than_Next3Peaks))=[];
        PitchMaxPeaksRatio(isnan(PitchMaxPeaksRatio))=[];
        
        FtrVec(1,1)=mean(ChangeBtw_MaxAndPitch_Power);
        FtrName(1,:)=['ChangeBtw_MaxAndPitch_Power         ' num2str(1)];
        
        FtrVec(2,1)=mean(MaxPeak_Above_Third_Harmony);
        FtrName(2,:)=['MaxPeak_Above_Third_Harmony         ' num2str(2)];
        
        FtrVec(3,1)=mean(PitchPeak_Lower_than_Next2Peaks);
        FtrName(3,:)=['PitchPeak_Lower_than_Next2Peaks     ' num2str(3)];
        
        FtrVec(4,1)=mean(PitchPeak_20PreLower_than_Next3Peaks);
        FtrName(4,:)=['PitchPeak_20PreLower_than_Next3Peaks' num2str(4)];
        
        FtrVec(5,1)=mean(PitchMaxPeaksRatio);
        FtrName(5,:)=['PitchMaxPeaksRatio(mean)            ' num2str(5)];
        
        FtrVec(6,1)=std(PitchMaxPeaksRatio);
        FtrName(6,:)=['PitchMaxPeaksRatio(std)             ' num2str(6)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Special frames - Time domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 7
        disp(' ***  Special Frames ***');
        clear Ceoff;
        frames=vec2frame(SpeechDat, Tw*Fs, Ts*Tw*Fs, 0);
        [PitchForFrame, AutoPeakAmp]= PitchPerFrame(frames, Fs);
        Pitch=median(PitchForFrame);[~, n_frames]=size(frames);
        [Abnorm2T0, Abnorm2Ti, LowAutoCorr, LowNumPeaks, TotPeaksRatio, TotNumSpecFrames, RatioNumSpecFrames] = SpecialSeg(frames(:,ceil(0.1*n_frames):n_frames*0.8), Fs, Tw, Pitch, AutoPeakAmp, PitchForFrame, Pathology, SPKR_ID);
        
        FtrVec(1,1)=Abnorm2T0;
        FtrName(1,:)=['Abnorm2T0         ' num2str(1)];
        
        FtrVec(2,1)=Abnorm2Ti;
        FtrName(2,:)=['Abnorm2Ti         ' num2str(2)];
        
        FtrVec(3,1)=LowAutoCorr;
        FtrName(3,:)=['LowAutoCorr       ' num2str(3)];
        
        FtrVec(4,1)=LowNumPeaks;
        FtrName(4,:)=['LowNumPeaks       ' num2str(4)];
        
        FtrVec(5,1)=TotPeaksRatio;
        FtrName(5,:)=['TotPeaksRatio     ' num2str(5)];
        
        FtrVec(6,1)=TotNumSpecFrames;
        FtrName(6,:)=['TotNumSpecFrames  ' num2str(6)];
        
        FtrVec(7,1)=RatioNumSpecFrames;
        FtrName(7,:)=['RatioNumSpecFrames' num2str(7)];
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   'S' Features
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    case 8
        disp(' ***  S Feature:  PeriodSumDiffSum,  DTW,  SAmpDiffSum, SAmpSTD, SIndDiffSum ***');
        clear Ceoff;
        frames=vec2frame(SpeechDat, Tw*Fs, Ts*Tw*Fs, 0);
        [PitchForFrame, AutoPeakAmp]= PitchPerFrame(frames, Fs);
        Pitch=median(PitchForFrame); [~, n_frames]=size(frames);
        [PeriodEnergyDiffEnt, DTW, SAmpEnt, SAmpMean, SIndDiff]=S_features(frames(:,ceil(0.05*n_frames):n_frames*0.85), Fs, Pitch, AutoPeakAmp, PitchForFrame);

        FtrVec(1,:)=PeriodEnergyDiffEnt; 
        FtrName(1,:)=['PeriodEnergyEnt' num2str(1)];
        
        FtrVec(2,:)=DTW;
        FtrName(2,:)=['DTW            ' num2str(2)];
        
        FtrVec(3,:)=SAmpEnt;
        FtrName(3,:)=['SAmpDiffEnt    ' num2str(3)];
        
        FtrVec(4,:)=SAmpMean;
        FtrName(4,:)=['SAmpMean       ' num2str(4)];
        
        FtrVec(5,:)=SIndDiff;
        FtrName(5,:)=['SIndDiff       ' num2str(5)];
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Lengths Statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    case 9
        disp(' ***   Speech lengths, entropy and SNR ***');
        clear Ceoff;      
         [V_len_sec, UV_len_sec, S_len_sec, V_entropy, UV_entropy, S_entropy]=SpeechLengths(SpeechDat,Tw, Fs);

        FtrVec(1,1)=V_entropy;
        FtrName(1,:)=['V_entropy ' num2str(1)];
         
        FtrVec(2,1)=UV_entropy; 
        FtrName(2,:)=['UV_entropy' num2str(2)];
 
        FtrVec(3,1)=S_entropy;         
        FtrName(3,:)=['S_entropy ' num2str(3)];

        FtrVec(4,1)=std(V_len_sec);
        FtrName(4,:)=['V_len_sec ' num2str(5)];

        FtrVec(5,1)=std(UV_len_sec);
        FtrName(5,:)=['UV_len_sec' num2str(5)];
 
        FtrVec(6,1)=std(S_len_sec);
        FtrName(6,:)=['S_len_sec ' num2str(6)];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Pitch Gradient and std
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 10
        disp(' *** Pitch Gradient and std ***');
        clear Ceoff;      
        [  frames ] =shir_vec2frames(SpeechDat,Fs,Tw,Ts);
        [ ~, n_frames]=size(frames);
        for i_frame=1:n_frames
            curr_frame=frames(:,i_frame);
            frame_after_threshold=threshold(curr_frame);
            pitch(1,i_frame)=pitch_period(frame_after_threshold,Fs);
        end
        median_pitch=median(pitch);
        pitch(find(pitch>400))=median_pitch;
        p = polyfit(1:length(pitch),pitch,1);
        
        FtrVec(1,1)=p(1); %gradient
        FtrName(1,:)=['Pitch Gradient' num2str(1)];
        
        FtrVec(2,1)=std(pitch);
        FtrName(2,:)=['Pitch STD     ' num2str(2)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  RMS Gradient and std
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 11
        disp(' *** RMS Gradient and std ***');
        clear Ceoff; 
        [  frames ] =shir_vec2frames(SpeechDat,Fs,Tw,Ts);
        RMS_changes=rms(frames);
        p = polyfit(1:length(RMS_changes),RMS_changes,1);
        
        FtrVec(1,1)=p(1); % gradient
        FtrName(1,:)=['RMS Gradient' num2str(1)];
        
        FtrVec(2,1)=std(RMS_changes);
        FtrName(2,:)=['RMS STD     ' num2str(2)];
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Autocorrelation Amplitude Per Frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 12
        disp(' ***   Autocorrelation Amplitude Per Frame ***');
        clear Ceoff;
        frames=vec2frame(SpeechDat,Tw*Fs, Ts*Tw*Fs, 0);
        [~, n_frames]=size(frames);
        [~, AutoPeakAmp]= PitchPerFrame(frames(:,ceil(0.05*n_frames):n_frames*0.85), Fs); 
        AutoPeakAmp=medfilt1(AutoPeakAmp);
        FtrVec(1,:)=AutoPeakAmp;
        FtrName(1,:)=['AutoCorrAmp' num2str(1)];
        
       plot(1:length(AutoPeakAmp), AutoPeakAmp, 1:length(AutoPeakAmp),ones(1,length(AutoPeakAmp))*mean(AutoPeakAmp), '-r'); xlabel('Number of Frame'); title(['HC Mean Value: ', num2str(mean(AutoPeakAmp))]); %xlim([1, length(AutoPeakAmp)]);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Autocorrelation Amplitude Per Record
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 13        
        disp(' ***   Autocorrelation Amplitude Per Record ***');
        clear Ceoff;
        [c, lags] = xcorr(SpeechDat, 'coeff'); MinPk_dis_samples=60;
        c = c(lags>=0);  lags = lags(lags>=0);  lags_s = lags/Fs;
        c(find(lags_s<0.0025))=0; %delete any freq above 400 Hz
        [pks,lcs] = findpeaks(c, 'MinPeakDistance', MinPk_dis_samples, 'MinPeakHeight',0.2*max(c), 'Npeaks' ,3);
        MinPk_dis_samples=lcs(find(pks==max(pks)));
        [pks,lcs] = findpeaks(c, 'MinPeakDistance', 0.9*MinPk_dis_samples, 'MinPeakHeight',0.1*max(c), 'Npeaks' ,7);
        
        FtrVec(1,1)=pks(1);
        FtrName(1,:)=['FirstAmp   ' num2str(1)];        
        
        FtrVec(2,1)=pks(2);
        FtrName(2,:)=['SecAmp     ' num2str(2)]; 
        
        FtrVec(3,1)=pks(1)-pks(4);
        FtrName(3,:)=['Diff1-4    ' num2str(3)];

        FtrVec(4,1)=sum(diff(pks(1:6)));
        FtrName(4,:)=['SixPksDiff ' num2str(4)];
        
        FtrVec(5,1)=sum(abs(diff(diff(lcs(1:6)/Fs))));
        FtrName(5,:)=['SixLcsDiff ' num2str(5)];
        
        P=polyfit(lags_s ,c',1);              
        FtrVec(6,1)=P(1);
        FtrName(6,:)=['SlopeaLL   ' num2str(6)];
        
        plot(lags_s ,c); hold on; plot( lcs/Fs , pks,'mo', 'MarkerSize', 7.5); hold on; plot(lags_s, polyval(P,lags_s), 'g',  'LineWidth',2.1); 
%         title('AutoCorrelation Per Record'); xlabel('Time (sec)'); xlim([0, 0.1]); legend('Autocorr Result', 'Slope', 'Peaks'); %xlim([0, length(lags_s)/Fs]);
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Entropy of speech
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 14
        disp(' ***   Entropy  ***');
        clear Ceoff;
        frames=vec2frame(SpeechDat,Tw*Fs, Ts*Tw*Fs, 0);
        [~, n_frames]=size(frames); frames=frames(:,ceil(0.05*n_frames):n_frames*0.85); [~, n_frames]=size(frames); 
        for i=1:n_frames
            ent(i)=entropy(frames(:, i));
        end
        FtrVec(1,:)=(ent); 
        FtrName(1,:)=['Entropy' num2str(1)];
        
        plot(ent); xlabel('Number of Frame'); title('Entropy'); xlim([1, length(ent)]);

end
