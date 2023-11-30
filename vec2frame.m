function [frameData] = vec2frame(fname,window_size_sampels, window_overlap,window)
% function that turns audio vector into frames without zero padding
% if window=1 => hanning, window=0 => no window
% if fname is array already, it should be array of rows!!

window_size_sampels=round(window_size_sampels);
if window==1
    window=hamming(window_size_sampels);
else
    window=1;
end

if ischar(fname)
    vec=audioread(fname); vec=zscore(vec);
else
    vec=fname; vec=zscore(vec);
end

m=length(vec);
if window_overlap==0 %determine the number of frames
    numFrames = floor(m/window_size_sampels)+1;
else
    numFrames = floor(m/(window_size_sampels-window_overlap*fs));
end

frameData = zeros(window_size_sampels, numFrames-1);  %allocate memory leaving out the end frame which is incomplete
startAtIdx = 1;

if (numFrames<=1)
    frameData=vec(:,1);
else
    for k=1:numFrames-1
        if k~=numFrames
            frameData(:,k) = window(:,1).*vec(startAtIdx:(startAtIdx+window_size_sampels-1),1);
        end
        startAtIdx = k*(1-window_overlap)*window_size_sampels; %overlap between the frames
    end
end

