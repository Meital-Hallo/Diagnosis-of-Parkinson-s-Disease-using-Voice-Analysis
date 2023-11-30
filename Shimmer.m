function [shimmer_dB, shimmer_abs]= Shimmer(A, fs, F0, AutoPeakAmp ,PitchPerFrame)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       ***     Shimmer per frame     ****
% Inputs: A is a record after VAD
%         Fs sampaling rate
%         The median pitch of the record F0
%         The autocorrelation amplitude
%         Pitch of each frame
% Outputs: Shimmer of the amplitudes in dB
%          Shimmer relative to the number of periods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[~,n]=size(A); shimmer_dB=zeros(1,n); shimmer_abs=zeros(1,n);

for i=1:n
    frame= A(:,i);
    if (AutoPeakAmp(i)>0.5)
        [pks, ~, ~]= FindPeaks_Func2(-frame, fs, PitchPerFrame(i), 0.75, 0.3, i);
        if length(pks)>2 %atleast three peaks
            num_of_periods=length(pks)-1;
            shimma= log_division(pks); 
            relative_amp=sum(abs(diff(pks)));
            shimmer_dB(i)=shimma/(num_of_periods-1);
            shimmer_abs(i)=relative_amp/(num_of_periods-1);
            
            if (isinf(shimmer_dB(i)))
                keyboard
            end
            
        else
            shimmer_dB(i)=NaN;
            shimmer_abs(i)=NaN;
        end 
    else
        shimmer_dB(i)=NaN;
        shimmer_abs(i)=NaN;
    end
end

shimmer_dB(isnan(shimmer_dB))=[];
shimmer_abs(isnan(shimmer_abs))=[];  

end



function [ result ] = log_division(array)
%function that returns the division of adjected numbers in array in dB

m=length(array);
A1=array(1:m-1);
A2=array(2:m);
result=sum(abs(20*log10(A2./A1)));
end

