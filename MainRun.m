clc; clear; close all

% global parameters
global T_frame; global t_overlap; global NormFtrMatrix;
global SaySTR; global START_FILE;

% Location defenitions
PARKINSON_DATABASE_CSV='C:\Users\Meital Hallo\Desktop\Project\DATA\PARKINSON_HEALTHY_SICK_UPDATE.CSV';
ResultsDirName='C:\Users\Meital Hallo\Desktop\Project\Results';
TextDataFileName='MachineLearningResults.txt';
outFTRsFileName=['C:\Users\Meital Hallo\Desktop\Project\Results\\TEST_Ftr_' SaySTR '_' date '_' ];

% Parameters
T_frame=0.05; %sec
t_overlap=0; %sec
START_FILE=1;
PatList=[1,3];
%--------------------
NormFtrMatrix=0;
SaySTR='MPT' ; %MPT CNT PATAKA SPONT T_ILND PIC_DESC
FtrListRUN=[13]; % Ftr index
Ftr_statList=[1 2 3 4 5 6 7 8 9]; 
        %---------------------------------
        % calcualte statistics of feature vector
        % 1=Median       6=MaxMinRatio
        % 2=Mean         7=MaxMeanRatio
        % 3=Std          8=MaxMedianRatio
        % 4=Min          9=Entropy
        % 5=Max         10=Nothing
        %---------------------------------

if (length(FtrListRUN)>1)
    FtrNameTxt='Combined';
else
    FtrNameTxt=FtrNameSTR(FtrListRUN);
end
FTRsFileNameRUN=[outFTRsFileName  FtrNameTxt  '.txt' ];

% Features Extraction & Classification
[FtrName]=CreateFtrMatrix(FtrListRUN,Ftr_statList,FTRsFileNameRUN,PARKINSON_DATABASE_CSV,PatList,SaySTR);

[Accuracy, Spec, Sens, NumPDRec, NumHCRec, NspkrPD, NspkrHC] = Machine_Learning(PatList,[1,3],FTRsFileNameRUN, FtrNameTxt);
disp(['Max Accuracy: ' num2str(Accuracy)]); 

% Add results to text file
f=fopen([ResultsDirName '\' TextDataFileName ],'at');
if (f>0)
    fprintf(f,'%s\n',[' Date ' date]);
    List=['Num of Frames Of PD ' num2str(NumPDRec) ', Of HC ' num2str(NumHCRec), ' FtrIndex=' num2str(FtrListRUN)];
    fprintf(f,'%s\n',List);
    fprintf(f,'%s\n',['FtrName=' FtrNameTxt, '  Statistics Index=' num2str(Ftr_statList)]);
    List=[' Number of Speakers in Each Pathology: ' num2str([NspkrPD NspkrHC])];
    fprintf(f,'%s\n',List);
    List=['  SVM', '         KNN'];
    fprintf(f,'%s\n',List);
    List=[ '% Correct Results=' num2str(Accuracy)];
    fprintf(f,'%s\n',List);
    List=[' Sensitivity: ' num2str(Sens)]; % The number of predicted PD from real PD
    fprintf(f,'%s\n',List);
    List=[' Specificity: ' num2str(Spec)];% The number of predicted HC from real HC
    fprintf(f,'%s\n',List);
    fprintf(f,'%s\n',' *************************************** \n');
    fclose(f);
end
