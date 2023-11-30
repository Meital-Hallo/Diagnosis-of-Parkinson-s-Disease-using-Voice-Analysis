function frames= shir_vec2frames(sample,fs,Tw,Ts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: a record, without silence in the begning and end
%        Fs sampling frequency
%        Tw window dutation
%        Ts shifting dutation
% Output: frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nw=round(fs*Tw);
Ns=round(fs*(Tw-Ts));
n_frames=floor((length(sample)-Nw)/Ns+1);
frames=zeros(Nw,n_frames);
j=1;
for i=1:n_frames
    frames(:,i)=sample(j:Nw+j-1);
    j=Ns*i;
end

