function  [Age, Gender, Pathology, PathoCode, WavFileNames, number_files, SPKR_ID]=ProcessEXCEL(CsvFileName,PatList,SaySTR)

%---  Load CSV file Data ---
ExcelData=csvimport(CsvFileName);
[n_rows, ~]=size(ExcelData);
number_files=0; i=1;

% parameters of excel
Dir_Collon=1;
FileName_collon=2;
Text_collon=3;
spkr_col=4;
patho_collon=5; %in old Excel 6
patho_code_collon=6;
gender_collon=7; %in old Excel 11
age_collon=8; %in old Excel 12

for i_row=1:n_rows-1
    
    % Reading type of Saying
    ToSayIt=cell2mat(ExcelData(i_row+1,Text_collon)); 
    c=ToSayIt(end);
    while (c==' ')
        ToSayIt=ToSayIt(1:end-1);
        c=ToSayIt(end);
    end
     
    SPKR_IDx=cell2mat(ExcelData(i_row+1,spkr_col)); 
    
    % find only those records with SaySRT
    if (strcmp(ToSayIt,SaySTR) | strcmp('ALL',SaySTR))
        Dirname=cell2mat(ExcelData(i_row+1,Dir_Collon));
        % fix space at end
        c=Dirname(end);
        while (c==' ') % (c==26)
            Dirname=Dirname(1:end-1);
            c=Dirname(end);
        end
        
        %construct file name
        Fname=[Dirname '\' num2str(cell2mat(ExcelData(i_row+1,FileName_collon))) '.wav'];
        if  (strfind(Fname,'.wav.wav') > 0)
            Fname=Fname(1:end-4);
        end
        disp(Fname)
        FileName=char(Fname);
        c=FileName(end);
        while (c~='v') & (c~='V')
            FileName=FileName(1:end-1);
            c=FileName(end);
        end
        WavFileNames(i,1:length(Fname))=FileName;
        
        % pathology number and name
        PathoCode(i)=cell2mat(ExcelData(i_row+1,patho_code_collon));
        disp(['Running file ' num2str(i_row) '  Name=' Fname ' PathoCode=' num2str(PathoCode(i)) ' Spkr= ' num2str(SPKR_IDx)]);
        if PathoCode(i)==1
            Pathology(i,1:7) ='PD     ';
        elseif  PathoCode(i)==3
            Pathology(i,1:7) ='HC     ';
        else
            Pathology(i,1:7) ='Unknown';
        end

        % Reading the Age of the speaker
        tmp=cell2mat(ExcelData(i_row+1,age_collon)); 
        if isnumeric(tmp)
            Age(i)=tmp;
        else
            Age(i)=str2num(tmp);
        end
        
        % Reading the gender of the speaker
        Gender(i)=cell2mat(ExcelData(i_row+1,gender_collon));
        
        % reading speaker number and cheaking if it is string
        if (isstr(SPKR_IDx))
            SPKR_ID(i)=str2num(SPKR_IDx);
        else
            SPKR_ID(i)=SPKR_IDx;
        end
        
        number_files=number_files+1;
        i=i+1;
    end
end

disp(['******************** Found ' num2str(number_files) ' Of Speech ' SaySTR ' ******************']);

end


