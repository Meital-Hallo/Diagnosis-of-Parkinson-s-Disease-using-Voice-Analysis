function [SpeechRate_V_UV,SpeechRate_V_ALL,SpeechRate_UV_ALL]=SpeechRate(speechIN,fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: a record, with silence in the begning and end
%        Fs sampling frequency
% Output: Duration of unvoiced segments divided by the duration of all speech signal
%         Duration of voiced segments divided by the duration of all speech signal
%         Duration of voiced segments divided by the duration of unvoiced segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tw=0.02; % window duration
[Speech, ~, Voiced, UnVoiced, ~] = VAD(speechIN ,fs);
   
VOICED_region_in_sec=length(Voiced)*Tw;
UNVOICED_region_in_sec=length(UnVoiced)*Tw;
ALL_region_in_sec=length(Speech)*Tw;
    
SpeechRate_V_UV=VOICED_region_in_sec/UNVOICED_region_in_sec;
SpeechRate_V_ALL=VOICED_region_in_sec/ALL_region_in_sec;
SpeechRate_UV_ALL=UNVOICED_region_in_sec/ALL_region_in_sec;

