function [FtrName]=FtrNameSTR(ftr_indx)
  
switch (ftr_indx)
    case 1
        FtrName='NHR_HNR_SNR';
    case 2
        FtrName='Spectral_Slope';
    case 3
        FtrName='SPEECH_RATE';
    case 4
        FtrName='Jitter';
    case 5 
        FtrName='Shimmer';
    case 6 
        FtrName='Abnormal_Frames_Freq';
    case 7
        FtrName='Abnormal_Frames_Time';
    case 8
        FtrName='S_Features';
    case 9
        FtrName='Lengths_Statistics';
    case 10
        FtrName='Pitch_Gradient_STD';
    case 11
        FtrName='RMS_Gradient_STD';
    case 12
        FtrName='Frame_AutoCor_Amp';
    case 13
        FtrName='Word_AutoCor_Amp';
    case 14
        FtrName='Entropy';

end

