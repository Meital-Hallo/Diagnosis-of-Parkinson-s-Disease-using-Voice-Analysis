function curr_frame=threshold(frame) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ***  Clipping Threshold Algorithm Per Frame  ***
% Input:  current frame
% Output: current frame after clipping threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frame_length=length(frame);
for i=1:frame_length
    threshold_val=max(frame)*0.3;
    if frame(i) >= threshold_val
        frame(i)= frame(i)- threshold_val;
    else if frame(i)<= (threshold_val)*(-1)
            frame(i)= frame(i)+ threshold_val;
        else if abs(frame(i))< threshold_val
                frame(i)=0;
            end
        end
    end
end

curr_frame=frame;