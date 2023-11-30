function [jitter_absolute,jitter_relative]  = Jitter(A, fs, F0 ,AutoPeakAmp, PitchForFrame)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ***     Jitter per frame     ****
% Inputs: A is a record after VAD
%         Fs sampaling rate
%         The median pitch of the record F0
%         The autocorrelation amplitudes
%         Pitch of each frame
% Outputs: jitter of the periods
%          jitter relative to the number of periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,n_frames]=size(A); jitter_absolute=zeros(1,n_frames); jitter_relative=zeros(1,n_frames);

for i=1:n_frames
    frame= A(:,i);
    if (AutoPeakAmp(i)>0.5)  
        [~, ind ,~]= FindPeaks_Func2(frame, fs, PitchForFrame(i), 0.75 ,0.3, i);
        if (length(ind)>2) %number of peaks must exceed 2 cause we cant divide in 0 (num_of_periods-1)
            periods=diff(ind);
            num_of_periods=length(periods);
            jitta=sum(abs(diff(periods)));
            sum_of_Ti=sum(periods);
            jitter_absolute(i)=jitta/(num_of_periods-1);
            jitter_relative(i)=100*(jitter_absolute(i)/(sum_of_Ti/num_of_periods));
        else
            jitter_absolute(i)=NaN;
            jitter_relative(i)=NaN;
        end       
    else
        jitter_absolute(i)=NaN;
        jitter_relative(i)=NaN;
        
    end
end

jitter_absolute(isnan(jitter_absolute))=[];
jitter_relative(isnan(jitter_relative))=[];