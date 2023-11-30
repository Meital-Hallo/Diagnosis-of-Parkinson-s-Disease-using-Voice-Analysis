function [StatVec, Stat_str] =CalcStati(FtrVec,i_stat)

%---------------------------------
% calcualte statistics of feature vector
% 1=Median       6=MaxMinRatio   
% 2=Mean         7=MaxMeanRatio 
% 3=Std          8=MaxMedianRatio  
% 4=Min          9=Entropy 
% 5=Max         10=Nothing
%---------------------------------

[FtrDim, n_frames]=size(FtrVec);
switch i_stat
    
    case 1 %median
        for i_c=1:FtrDim
            StatVec(i_c)=median(FtrVec(i_c,:));
        end
        Stat_str='median';       
    
    case 2 %mean
        for i_c=1:FtrDim
            StatVec(i_c)=mean(FtrVec(i_c,:));
        end
        Stat_str='mean';
        
    case 3 %std
        for i_c=1:FtrDim
            StatVec(i_c)=std(FtrVec(i_c,:));
        end
        Stat_str='std';  
   
    case 4 %min
        for i_c=1:FtrDim
            StatVec(i_c)=min(FtrVec(i_c,:));
        end
        Stat_str='min';
        
    case 5 %max
        for i_c=1:FtrDim
            StatVec(i_c)=max(FtrVec(i_c,:));
        end
        Stat_str='max';
        
    case 6 %MaxMinRatio
        for i_c=1:FtrDim
            StatVec(i_c)=max(FtrVec(i_c,:))-min(FtrVec(i_c,:));
        end
        Stat_str='MaxMinRatio';
        
     case 7 %MaxMeanRatio
        for i_c=1:FtrDim
            StatVec(i_c)=max(FtrVec(i_c,:))-mean(FtrVec(i_c,:));
        end
        Stat_str='MaxMeanRatio';        
        
     case 8 %MaxMedianRatio
        for i_c=1:FtrDim
            StatVec(i_c)=max(FtrVec(i_c,:))-median(FtrVec(i_c,:));
        end
        Stat_str='MaxMedianRatio';   
                
    case 9 %Entropy
        for i_c=1:FtrDim
            StatVec(i_c)=entropy(real(FtrVec(i_c,:)));
        end
        Stat_str='Entropy';
    
    case 10 %Nothing - In case he feature returns one value
            StatVec=FtrVec;
        Stat_str='X';
         
end

end
