function [pks,lcs, val] = FindPeaks_Func2(Frame, fs, F0, MinPeakDistance, MinPeakHeight, NumberOfFrame)

m=length(Frame);
t=(0:1:m-1)/fs;

%finding Area around maximal point of the frame
Area_size_samples=round(0.5*fs/F0); %choosing entire period around max peak
[~, MaxInd_sample]=max(Frame);

if (MaxInd_sample-Area_size_samples>0 && MaxInd_sample+Area_size_samples<=m)%if max point is located at center
    Area=Frame((MaxInd_sample-Area_size_samples):(MaxInd_sample+Area_size_samples));
    %figure(1); plot(t,Frame,((MaxInd_sample-Area_size_samples):(MaxInd_sample-Area_size_samples+length(Area)-1))/fs,Area,'r-'); title(['The Chosen Window ,', 'F0=' num2str(F0)]); xlabel('Time (sec)'); legend('Frame', 'The Chosen Window')
    [~, MaxInd_sample]=max(Area);%need this index for correction
elseif (MaxInd_sample+Area_size_samples>m) %if max point is located at the end of the frame
    Area=Frame((m-2*Area_size_samples):m);
    %figure(1); plot(t,Frame,((m-2*Area_size_samples):m)/fs,Area,'r-');title(['The Chosen Window ,', 'F0=' num2str(F0)]);xlabel('Time (sec)'); legend('Frame', 'The Chosen Window')
    [~, MaxInd_sample]=max(Area); %need this index for correction
else %if max point is located at the end of the frame
    Area=Frame(1:(1+2*Area_size_samples));
    %figure(1); plot(t,Frame,(1:(1+2*Area_size_samples))/fs,Area,'r-');title(['The Chosen Window ,', 'F0=' num2str(F0)]);xlabel('Time (sec)'); legend('Frame', 'The Chosen Window')
    [~, MaxInd_sample]=max(Area); %need this index for correction
    if (MaxInd_sample==1)
        [~,locc]=findpeaks(Area,'SortStr','descend');
        MaxInd_sample=locc(1);
    end
end

k=0; %Runing with window on the frame
for i=1:m-length(Area)
    S(i)=sum(Frame(i:(i+length(Area)-1)).*Area); 
    %figure(1); plot(t, Frame,(i:length(Area)+k)/fs, Area,'r--'); k=k+1;
end

[S_MaxAmp,S_MaxIndex]=findpeaks(S,'MinPeakDistance',MinPeakDistance*fs/F0, 'minpeakheight',MinPeakHeight*max(S)); 
%[RealPeaks, RealIndex]=findpeaks(Frame,fs,'MinPeakDistance',MinPeakDistance/F0, 'minpeakheight',MinPeakHeight*max(Frame));

%shifting and fixing the index of the peaks by searching for the maximal peak around it
if ~(isempty(S_MaxAmp))
    FixInd=S_MaxIndex+MaxInd_sample; k=1;
    for i=FixInd
        if FixInd(1)-38>1 && FixInd(end)+38<m
            [pks(k), lcss]=max(Frame(FixInd(k)-38:FixInd(k)+38));
            lcs(k)=FixInd(k)-39+lcss;
            k=k+1;
        else
            val= min(FixInd(1)-1, m-FixInd(end)-1);
            [pks(k), lcss]=max(Frame(FixInd(k)-val:FixInd(k)+val));
            lcs(k)=FixInd(k)-val-2+lcss;
            k=k+1;
        end
    end
    lcs=lcs/fs;
    
    %fixing peaks in the start and end of the frame incase of malfunction.
    if (length(lcs)>=2 && S_MaxAmp(1)<0.7*S_MaxAmp(2) && lcs(1)<0.006) %start
        [pks(1), lcss] =max(Frame(1:lcs(1)*fs));
        lcs(1)=lcss/fs;
        %disp('correction')
    end

    if (length(lcs)>=2 && lcs(end)>0.045 && abs(S_MaxAmp(end-1)-S_MaxAmp(end))<0.65*S_MaxAmp(end-1))
        [pks(end), lcss] =max(Frame(lcs(end)*fs-35:end-35));
        lcs(end)=lcs(end)-36/fs+lcss/fs;
        %disp('correction MAX')
    end
    
    S_MaxAmp=S_MaxAmp/max(S_MaxAmp);
    %return value that helps correct the pitch
    GR1_S_MaxAmp=mean(S_MaxAmp(1:2:end));
    GR2_S_MaxAmp=mean(S_MaxAmp(2:2:end));
    if abs(GR1_S_MaxAmp-GR2_S_MaxAmp)>0.5
        val=1;
    else
        val=0;
    end
        
    subplot(2,1,1); plot(t, Frame, ((1:length(S))+MaxInd_sample-1)/fs, S, '-', (S_MaxIndex+MaxInd_sample-1)/fs, S(S_MaxIndex), 'r*');
    title(['Running Window Results of Frame '  num2str(NumberOfFrame)]); xlabel('Time (sec)'); ylabel('Amplitude')
    subplot(2,1,2); plot(t, Frame,lcs, pks,'y*'); % plot(t, Frame, RealIndex, RealPeaks, 'r*',lcs, pks,'y*');
    title(['Peaks of Frame '  num2str(NumberOfFrame)]);  xlabel('Time (sec)'); ylabel('Amplitude');
    
else
    val=1; pks=[]; lcs=[];
end


