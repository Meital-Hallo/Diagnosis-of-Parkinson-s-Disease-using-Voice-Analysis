function pitch_frame=pitch_period(frame_after_threshold,fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ***  Pitch Per Frame   ***
% Input:  current frame after threshold function 
%         Fs sampling frequency
% Output: Pitch Per Frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R=autocorrelation_func(frame_after_threshold);
R_length=length(R);
max_val=max(R);
R_axis=(1:length(R))/(fs);
MinPk_dis_samples=0.007*(fs);
[pks,loc] = findpeaks((1/max_val)*(R),'SortStr','descend','MinPeakDistance',MinPk_dis_samples);
pitch_frame=(fs)/loc(1);
loc1=loc(1);
% in a case where the highest peak is located at harmonic period instead of pitch period
% a threshold value is defined, which is equal to the location of the maximum multiplied by a constant, 
% in a range of 0.8 to 1.2.
% If the sampling number of the second highest maximum multiplied by a period number (2 or 3) is within that range, 
% it is assigned as the pitch period.
if length(loc)>=2
   if loc(1)>loc(2)
      if ((loc(2)*2)<(1.2*loc(1))) && ((loc(2)*2)>(0.8*loc(1))) 
         pitch_frame=(fs)/loc(2);
         loc1=loc(2);
      else
         pitch_frame=(fs)/loc(1);
         loc1=loc(1);
      end
  end
   if loc(1)>loc(2)
      if ((loc(2)*3)<(1.2*loc(1))) && ((loc(2)*3)>(0.8*loc(1)))
         pitch_frame=(fs)/loc(2);
         loc1=loc(2);
      end
  end
end


