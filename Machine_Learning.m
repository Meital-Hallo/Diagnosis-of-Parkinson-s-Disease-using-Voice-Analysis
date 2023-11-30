function [Accuracy, Spec, Sens, NumPDRec, NumHCRec, NspkrPD, NspkrHC] = Machine_Learning(Patho_class, ML_TYPE, INdBFileName, FtrNameTxt)
%ML_TYPE -> 1=KNN, 3=SVM

global NormFtrMatrix;
global SaySTR;
PathoCol=5;SpkrCol=2;
DATA_START_COL=7; 

if string(SaySTR)=="MPT"
    KNN=11;
else
    KNN=5;
end

% initialize values
NumPDRec=0; NumHCRec=0; 
Accuracy_SVM=-1; Sense_SVM=-1; Spec_SVM=-1;
Accuracy_knn=-1; Spec_knn=-1; Sense_knn=-1;

% read CSV features file
disp(['  Features File Name: ' INdBFileName]);
disp('     >>>>  Reading Excel File  <<<');
ExcelData=csvimport(INdBFileName);
[n_rows, n_cols]=size(ExcelData);
while isempty(cell2mat(ExcelData(3,n_cols))) & (n_cols>DATA_START_COL)
    n_cols=n_cols-1;
end
if (n_cols<DATA_START_COL)
    disp('No data'); pause;
end

disp('Checking Start Col -->> This Should Be The Data:')
ExcelData(2:10,DATA_START_COL)

% reading features and labels
disp([' Preparing data. Number of Rows=' num2str(n_rows)]);
for i=1:n_rows-1
    SpkrID(i)=cell2mat(ExcelData(i+1,SpkrCol));
    label(i)=cell2mat(ExcelData(i+1,PathoCol));
    if (label(i)==Patho_class(1)) % SICK=3
        NumPDRec=NumPDRec+1;
    else  % HEALTHY=1
        NumHCRec=NumHCRec+1;
    end
    %reading the featues into ftrMat
    for i_col=DATA_START_COL:n_cols
        if ischar(cell2mat(ExcelData(i+1,i_col)))
            ftrMat(i,(i_col-DATA_START_COL+1))=str2num(cell2mat(ExcelData(i+1,i_col)));
        else
            ftrMat(i,(i_col-DATA_START_COL+1))=cell2mat(ExcelData(i+1,i_col));
        end
    end
end
[DataRow, DataCol]=size(ftrMat);

if size(ftrMat, 2)==1 
    SpkrID(find(isnan(ftrMat)))=[];
    label(isnan(ftrMat))=[];
    ftrMat(isnan(ftrMat))=[];
else
    SpkrID(find(isnan(ftrMat(: ,1))))=[];
    label(isnan(ftrMat(: ,1)))=[];
    ftrMat(isnan(ftrMat(: ,1)), :)=[];   
end   
    

%Normalizing features
if (NormFtrMatrix==1)
    for i_col=1:DataCol
        Col=ftrMat(:,i_col);  v1=std(Col);  Col=Col-mean(Col);
        if (v1~=0)
            Col=Col/v1;
        end
        ftrMat(:,i_col)=Col;
    end
end

ftrMat=ftrMat(:, ([2]));


% counts how much speakers of each pathology
SpkrLIST=unique(SpkrID); %Returns all speakers number with no repitition sorted
for i_spkr=1:length(SpkrLIST)
    SpkrLoc=find(SpkrLIST(i_spkr)==SpkrID);
    SpkrLIST_Patho(i_spkr)=label(SpkrLoc(1));
end
NspkrPD=sum(SpkrLIST_Patho==Patho_class(1));
NspkrHC=sum(SpkrLIST_Patho==Patho_class(2));

disp([' Pathology(' num2str(Patho_class(1)) ') Samples = ' num2str(NumPDRec) ' Pathology(' num2str(Patho_class(2)) ') Samples = ' num2str(NumHCRec) ]);
disp([' Number of PD Speakers (Patho is ' num2str(Patho_class(1)) ') is ' num2str(NspkrPD) ]);
disp([' Number of HC Speakers (Patho is ' num2str(Patho_class(2)) ') is ' num2str(NspkrHC) ]);

k_fold=length(ftrMat); 
ind=crossvalind('kfold',k_fold, k_fold); %leave 1-out
disp(' *** Starting Classifiers *** ')

for i=1:k_fold
    test_ind=(ind==i); % return ones where ind is equal to i.
    train_ind=~test_ind; % return ones where ind is not-equal to i.
    
    % building test and train groups with the indexs abtained above
    x_train=ftrMat(train_ind,:);
    x_test=ftrMat(test_ind,:);
    y_train=label(train_ind);
    y_test=label(test_ind);

    if (sum(ML_TYPE==1)>0)
        disp([' <<<  KNN Iteration Number ' num2str(i), '   >>> '])
        KNN_Mdl = fitcknn(x_train, y_train,'NumNeighbors', KNN);
        KNN_Lab = predict(KNN_Mdl, x_test);
        CM_KNN(:,:,i) = confusionmat(y_test,KNN_Lab, 'order', [3 ;1] )/size(x_test, 1);
    end
    
    if (sum(ML_TYPE==3)>0)
        disp([' <<<  SVM Iteration Number ' num2str(i), '   >>> '])
        SVM_mdl = fitcsvm(x_train, y_train, 'KernelFunction','polynomial','boxconstraint',1);
        SVM_Lab = predict(SVM_mdl, x_test);
        CM_SVM(:,:,i) = confusionmat(y_test, SVM_Lab, 'order', [3 ;1] )/size(x_test, 1);
    end
end

if (sum(ML_TYPE==1)>0)
    MeanCM_knn= sum(CM_KNN, 3)/k_fold;
    Accuracy_knn=trace(MeanCM_knn)*100;
    Spec_knn=MeanCM_knn(2,2)/(MeanCM_knn(2,1)+MeanCM_knn(2,2));
    Sense_knn=MeanCM_knn(1,1)/(MeanCM_knn(1,1)+MeanCM_knn(1,2));
end

if (sum(ML_TYPE==3)>0)
    MeanCM_SVM= sum(CM_SVM, 3)/k_fold;
    Accuracy_SVM=trace(MeanCM_SVM)*100;
    Sense_SVM=MeanCM_SVM(2,2)/(MeanCM_SVM(2,1)+MeanCM_SVM(2,2));
    Spec_SVM=MeanCM_SVM(1,1)/(MeanCM_SVM(1,1)+MeanCM_SVM(1,2));
end

Accuracy=[Accuracy_SVM, Accuracy_knn];
Sens=[Sense_SVM, Sense_knn];
Spec=[Spec_SVM, Spec_knn];

if size(ftrMat, 2)==1 % 1D feature
    histogram(ftrMat(find(label==1)), 20, 'BinWidth', 0.03); % sick
    hold on
    histogram(ftrMat(find(label==3)), 20, 'BinWidth',0.03);
    title(['Plot of ' ,FtrNameTxt, ' ' num2str(max(Accuracy)), '%' ]) 
    legend('PD' , 'HC')
    xlabel('Feature Values');
    ylabel('Number of Records');
end

if size(ftrMat,2)==2 % 2D feature
    figure; gscatter(ftrMat(:,1), ftrMat(:,2), label, ['b', 'r']); grid on;
    title(['Features Case: ', FtrNameTxt, ' ' num2str(max(Accuracy)), '%' ]) 
    legend('PD' , 'HC')
    %Name=FtrNameSTR(str2num(FtrNameTxt(1))); xlabel(Name);
    %Name=FtrNameSTR(str2num(FtrNameTxt(end-1:end))); ylabel(Name);
end
 
if size(ftrMat,2)==3 % 3D feature
      figure; plot3(ftrMat(find(label==1), 1), ftrMat(find(label==1), 2), ftrMat(find(label==1), 3), 'b.', 'MarkerSize', 15); grid on; hold on;
      plot3(ftrMat(find(label==3), 1), ftrMat(find(label==3), 2), ftrMat(find(label==3), 3), 'r.', 'MarkerSize', 15); 
      title(['3D Plot of ' ,FtrNameTxt, ' ' num2str(max(Accuracy)), '%' ]) 
      legend('PD' , 'HC')
end
 
end

