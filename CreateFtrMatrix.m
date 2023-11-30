function [n_files]=CreateFtrMatrix(FtrList,Ftr_statList,OutdBFileName,CsvFileName,PatList,SaySTR)
% This Function Create EXCEL MATRIX OF ACOUSTIC FEATURES

global START_FILE;
MIN_SPEECH_LEN_SEC=0.1;
HD_TEXT_flg=1; % create header when =1 only (first loop)
First_Line_flg=1; %only when runs for the first time enters the loops of headers
 
%Reading Excel with speaker data
[Age, Gender, Pathology, PathoCode, WavFileNames, n_files, SPKR_ID]=ProcessEXCEL(CsvFileName,PatList,SaySTR);
f_dB=fopen(OutdBFileName,'wt');

for i_file=START_FILE:n_files
    
    if (START_FILE>1)
        disp('    WARNING!! --- NOT STARTING FROM FIRST FILE ---')
    end
    
    % Reading the address of record in line i_file and removing spaces from end
    FileName=char(WavFileNames(i_file,:));
    c=FileName(end);
    while (c~='v') & (c~='V')    
        FileName=FileName(1:end-1);
        c=FileName(end);
    end
       
    disp(['=  iFile=' num2str(i_file) ' out of ' num2str(n_files) ' progress =' num2str(i_file/n_files*100) '%']);
    disp(['=  Reading file - ' FileName ' pathology is ' Pathology(i_file,1:7) ' **** ']);
    
    % Reading the audio and VADing
    [SpeechIN, Fs]=audioread(FileName);
    [~, ~, Voiced, ~, ~]=VAD(SpeechIN ,Fs);
    
    % Checking the length of the record
    SpeechLEN=length(Voiced)/Fs;
    if (SpeechLEN<MIN_SPEECH_LEN_SEC)
        disp('Seg too Short');
        keyboard;
    end
    
 %-------------------------------------------------------
 % loop on all features
 %-------------------------------------------------------
    for i_Ftr=1:length(FtrList)
        clear FtrVec; clear FtrName
        FtrIndex=FtrList(i_Ftr);
        if (FtrIndex==3 || FtrIndex==9) % These features need full data without VAD
            SpeechDat=SpeechIN;
        else
            SpeechDat=Voiced;
        end
        
        % calculate features
        [FtrVec, FtrName]=CalcFeatures(FtrIndex,SpeechDat, Fs, Pathology, SPKR_ID);
        [FtrDim, n_frames]=size(FtrVec); 
        FtrDimm(i_Ftr)=FtrDim;          
         
        %------ calculating statistics-----------  
        for n_stat=1:length(Ftr_statList)  %Looping on the list of wanted statistics
            i_stat=Ftr_statList(n_stat);
            [StatVec, StatNameSTR]=CalcStati(FtrVec, i_stat); %calculate statistics
            STATCeoff(i_Ftr,n_stat,1:FtrDim)=StatVec;
            
            % Create header with Statistics names
            for i_c=1:FtrDim
                if (HD_TEXT_flg==1)
                    HD_TEXT_flg=0;
                    HD_TEXT=[StatNameSTR '_' FtrName(i_c,:)];
                else
                    HD_TEXT=[HD_TEXT ',' StatNameSTR '_' FtrName(i_c,:)];
                end
            end
        end
    end % iftr
    
    %first line headers
    if (First_Line_flg==1)
        First_Line_flg=0;
        Text1=sprintf('%s','FileName, SpkrID, Saying, PATHOLOGY_NAME, PATHOLOGY_ID, FrameNum ,');
        Headers=[Text1 HD_TEXT];
        e=fprintf(f_dB,'%s\n',Headers);
    end
    
    %write speaker data into text file
    coll1=char(WavFileNames(i_file,:));
    coll2=num2str(SPKR_ID(i_file)); 
    coll3=SaySTR;
    coll4= Pathology(i_file,:);
    coll5=num2str(PathoCode(i_file));
    FileLineText=[coll1 ',' coll2 ',' coll3 ',' coll4 ',' coll5];
    
    CeoffTxt='';
    %write numeric data to text file
        for i_Ftr=1:length(FtrList)
            for n_stat=1:length(Ftr_statList)
                for i_c=1:FtrDimm(i_Ftr)
                    CeoffTxt=[CeoffTxt num2str(STATCeoff(i_Ftr,n_stat,i_c)) ',' ];
                end
            end
        end
        
        LineText=[FileLineText ', -1 ,' CeoffTxt];
        e=fprintf(f_dB,'%s\n',LineText);
        
end %end loop on all files
fclose(f_dB);
end