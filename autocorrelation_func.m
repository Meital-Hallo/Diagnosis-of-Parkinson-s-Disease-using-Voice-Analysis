function R=autocorrelation_func(frame_mat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: frame_mat is a matrix that contains VOICE segments 
% Output: Autocorrelation Result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[m_samples,n_frames]=size(frame_mat);
mat=zeros(m_samples,n_frames);
s=0;d=0;
for i=1:n_frames
    for j=1:m_samples
        for k=1:m_samples-d
           mat(j,i)= s+frame_mat(k,i).*frame_mat(k+d,i);
           s=mat(j,i);
        end
        d=d+1;s=0;
    end
    d=0;
end
R=mat;
