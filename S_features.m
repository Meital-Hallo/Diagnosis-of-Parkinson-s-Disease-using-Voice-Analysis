function [PeriodEnergyDiffEnt, DTW, SAmpEnt, SAmpMean, SIndDiff] = S_features(A, fs, F0, AutoPeakAmp, PitchForFrame)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ***     Function for features of S from findpeaks per frame     ****
% Inputs: A is a record after VAD
%         Fs sampaling rate
%         The median pitch of the record F0
%         The autocorrelation amplitude
%         Pitch of each frame
% Outputs: Entropy of the difference between the sum of the periods in each frame
%          DTW between the first and last periods of the frame
%          Entropy of the changes of the peaks in each frame
%          Mean of the peaks in each frame
%          The deviation of the periods in S
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[~,n]=size(A); A=A/max(A(:));
for i=1:n
    Frame=A(:,i);
    [S_func, S_pks ,S_ind] = S_PeaksAutocorr(Frame, fs, PitchForFrame(i), 0.8 ,0.4, i);
    if PitchForFrame(i)>0 && AutoPeakAmp(i)>0.4 && length(S_ind)>2
        %figure; plot(((1:length(S_func)))/fs, S_func, 'k-', S_ind/fs, S_pks ,'r*'); title(' S Function'); legend('S', 'Peaks');
        for j=1:length(S_pks)-1
            PeriodSum(j)=sum(abs(S_func(S_ind(j):S_ind(j+1))));
            plot((0:length(S_func)-1)/fs, S_func, (S_ind(j):S_ind(j+1))/fs, S_func(S_ind(j):S_ind(j+1)))
        end
        DTW1=S_func(S_ind(1):S_ind(2));
        DTW2=S_func(S_ind(end-1):S_ind(end));
        DTW(i)=dtw(DTW1, DTW2);
        
        SAmpEnt(i)=entropy(abs(diff(S_pks))); 
        
        S_MaxAmp=S_pks; S_MaxAmp(find(max(S_MaxAmp)))=[];
        SAmpMean(i)=mean(S_MaxAmp); 
        
        SIndDiff(i)=sum(abs(diff(diff(S_ind))))/fs;
        
        PeriodEnergyDiffEnt(i)=entropy(abs(diff(abs(PeriodSum)))); 
        
    else
        PeriodEnergyDiffEnt(i)=nan;
        SAmpMean(i)=nan;
        DTW(i)=nan;
        SAmpEnt(i)=nan;
        SIndDiff(i)=nan;
    end
end

PeriodEnergyDiffEnt(isnan(PeriodEnergyDiffEnt))=[];
DTW(isnan(DTW))=[];
SAmpMean(isnan(SAmpMean))=[];
SAmpEnt(isnan(SAmpEnt))=[];
SIndDiff(isnan(SIndDiff))=[];

end

function [S, S_pks ,S_ind] = S_PeaksAutocorr(Frame, fs, F0, MinPeakDistance, MinPeakHeight, NumberOfFrame)

m=length(Frame);
t=(0:1:m-1)/fs;

%finding Area around maximal point of the frame
Area_size_samples=round(0.5*fs/F0); %choosing entire period around max peak
[~, MaxInd_sample]=max(Frame);

if (MaxInd_sample-Area_size_samples>0 && MaxInd_sample+Area_size_samples<=m)%if max point is located at center
    Area=Frame((MaxInd_sample-Area_size_samples):(MaxInd_sample+Area_size_samples));
elseif (MaxInd_sample+Area_size_samples>m) %if max point is located at the end of the frame
    Area=Frame((m-2*Area_size_samples):m);
else %if max point is located at the end of the frame
    Area=Frame(1:(1+2*Area_size_samples));
end

%Runing with window on the frame
for i=1:m-length(Area)
    S(i)=sum(Frame(i:(i+length(Area)-1)).*Area);
end

S=S/max(S);
[S_pks,S_ind]=findpeaks(S,'MinPeakDistance',MinPeakDistance*fs/F0, 'minpeakheight',MinPeakHeight*max(S));

figure; 
subplot(2,1,1); plot(t, Frame); title('Frame'); xlabel('Time (sec)'); ylabel('Amplitude')
subplot(2,1,2); plot(((1:length(S)))/fs, S, 'k-', S_ind/fs, S_pks ,'r*'); title(' S Function'); legend('S', 'Peaks');
end
