 function [V_len_sec, UV_len_sec, S_len_sec, V_entropy, UV_entropy, S_entropy ] = SpeechLengths (Speech,window_size_sec, fs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: a record, with silence in the begning and end.
% Output: V_length_sec is the lengths of voiced segments and V_length_entropy is its entropy
%         UV_length_sec is the lengths of Unvoiced segments and UV_length_entropy is its entropy
%         S_length_sec is the lengths of silence segments and S_length_entropy is its entropy
%         SNR between median RMS of voiced to min RMS of silence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

window_size_sampels=window_size_sec*fs;

% VAD
[~, Result, Voiced, ~, AllSilence] = VAD(Speech, fs);
Result=medfilt1(Result);

%Evaluating Speech Rate and lengths entropy
index=find(diff(Result)~=0); %finding the index of change between the types
index=[index, length(Result)]; %adding the index of the end
lengths=diff([0,index]);
v=1; uv=1; s=1; L=1; unvoice_length=0;voice_length=0; silence_length=0;
%plot(1:length(Result), Result, 'r*-');

for i=index(1:end)
    if Result(i)==1
        voice_length(v)=lengths(L);
        v=v+1;
    end
    
    if Result(i)==0.5
        unvoice_length(uv)=lengths(L);
        uv=uv+1;
    end
    if Result(i)==0 && i~=index(1) && i~=index(end)
        silence_length(s)=lengths(L);
        s=s+1;
    end
    
    L=L+1;
end

V_len_sec=voice_length*0.02;
UV_len_sec=unvoice_length*0.02;
S_len_sec=silence_length*0.02;
V_entropy=entropy(V_len_sec);
UV_entropy=entropy(UV_len_sec);
S_entropy=entropy(S_len_sec);


end


