function [accuracy, Predict_test_Label,SpecificitySensitivity,weightall,FeatureAll,SelectFeatureIndex,Specif_Sensit_acc_fold,MaxDim]=SVMMVPA(Cfg)
% using support vector machines for classification
% accuracy The correct rate of classification
% P is specific and sensitive
Fold=Cfg.Fold;
strSelectFeature=Cfg.FeatureSelectionStr;
Data=Cfg.Data;
nDim=Cfg.ContrPCA;

Label=Cfg.Label;

modelParemter=Cfg.modelParemter;
ID=findstr(modelParemter,'-t');
FishZStr=Cfg.FishZStr;
strFeatureReduce=Cfg.FeatureReduceStr;
SelectFeatureNum=Cfg.SelectFeatureNum;
SelectFeatureIndex=[];
if ~strcmp(strFeatureReduce,'PCA (Principal Component Analysis)')
    MaxDim='NODimensionalityReduction';
else
    for i=1:max(Fold)    
        train_data=Data;
        train_data(find(Fold==i),:)=[];
        if strcmp(FishZStr,'Normalization for Each Feature')
            train_data=zscore(train_data);
        end
        if strcmp(strFeatureReduce,'PCA (Principal Component Analysis)')
            [coeff,DataScore,latent,tsquared,explained] = pca(train_data);
            RContributionLatent=num2cell([explained, cumsum(explained)]);
            RContributionLatent=cell2mat(RContributionLatent);
            MaxDimPCA(i)=min(find(RContributionLatent(:,2)>nDim));
        end
    end
    MaxDim=max(MaxDimPCA);
end
%%
%20200113 edit 
if ~strcmp(strFeatureReduce,'PCA (Principal Component Analysis)')
    if ~strcmp(strSelectFeature,'None') && ~strcmp(strSelectFeature,'Lasso')
        %Feature selection
        weightall=zeros(size(Data,2),length(SelectFeatureNum),max(Fold));
        FeatureAll=zeros(size(Data,2),length(SelectFeatureNum),max(Fold));
        SelectFeatureIndex=zeros(size(Data,2),max(Fold));
    else
        %No feature selection
        weightall=zeros(size(Data,2),max(Fold));
        SelectFeatureIndex=[];
        FeatureAll=zeros(size(Data,2),max(Fold));
    end
else
    weightall=[];
    FeatureAll=[];
end
%%
ReLabel=[];
train_data=[];
Real_Label1=[];
predict_Label1=[];
Real_Label2=[];
predict_Label2=[];

ReLabel=unique(Label);%true label 0 1; 1 -1
h=waitbar(0,'Runing...');
for i=1:max(Fold)
    str=['Runing...',num2str(i/max(Fold)*100),'%'];
    waitbar(i/max(Fold),h,str);
    IX=[];
    weight=[];
    test_data=Data(find(Fold==i),:);
    test_label=Label(find(Fold==i));
    train_label=Label;
    train_data=Data;
    train_label(find(Fold==i),:)=[];%i fold as test data, other as train data
    train_data(find(Fold==i),:)=[];
    if strcmp(FishZStr,'Normalization for Each Feature')
        meantrain_data=mean(train_data);
        stdtrain_data=std(train_data);
        train_data=zscore(train_data);
        for izscore=1:size(train_data,2)
            test_data(:,izscore)=(test_data(:,izscore)-meantrain_data(izscore))./stdtrain_data(izscore);
        end
    end
    if strcmp(strFeatureReduce,'PCA (Principal Component Analysis)')
        [coeff,DataScore,latent,tsquared,explained] = pca(train_data);
        test_data1=bsxfun(@minus,test_data,mean(train_data,1))*coeff;
        train_data=DataScore(:,1:MaxDim);
        test_data=test_data1(:,1:MaxDim);
    end
    model = svmtrain(train_label,train_data,modelParemter);
    if length(ReLabel)==2 && modelParemter(ID+3)=='0'
        weight=model.SVs'*model.sv_coef;%compute weight
        weightall(:,i)=weight; %all weight
    end
    %%
    switch strSelectFeature
        case 'None'%

            [predict_label, acc, dec_values] =svmpredict(test_label,test_data, model);
            
            accuracy(i)=acc(1);
            for IDRelabel=1:length(ReLabel)
                Specif_Sensit_acc_fold(i,IDRelabel)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
            end
            Specif_Sensit_acc_fold(i,length(ReLabel)+1)=length(find(predict_label==test_label))/length(test_label)*100;
            PredictLabel(find(Fold==i),1)=predict_label;
            RealLabel(find(Fold==i),1)=test_label;
            DecValues(find(Fold==i),1)=dec_values(:,1);         
        case 'Lasso'
            %%
            %20200113 edit
            [B,FitInfo]=lasso(train_data,train_label,'CV',10);
            IndexMinMSE=FitInfo.IndexMinMSE;
            train_data=train_data(:,find(B(:,IndexMinMSE)~=0));
            model = svmtrain(train_label,train_data,modelParemter);
            FeatureAll(find(B(:,IndexMinMSE)~=0),i)=1;
            if length(ReLabel)==2 && modelParemter(ID+3)
                weight=[];
                weight=model.SVs'*model.sv_coef;%compute weight
                weightall(find(B(:,IndexMinMSE)~=0),i)=weight; %all weight
                
            end
            test_data=test_data(:,find(B(:,IndexMinMSE)~=0));
            [predict_label, acc, dec_values] =svmpredict(test_label,test_data, model);
            accuracy(i)=acc(1);
          
            for IDRelabel=1:length(ReLabel)
                Specif_Sensit_acc_fold(i,IDRelabel)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
            end
            Specif_Sensit_acc_fold(i,length(ReLabel)+1)=length(find(predict_label==test_label))/length(test_label)*100;
            PredictLabel(find(Fold==i),1)=predict_label;
            RealLabel(find(Fold==i),1)=test_label;
            DecValues(find(Fold==i),1)=dec_values(:,1);
        case 'Weight'
            %%
            weightsvm=abs(weight);
            [B,IX]=sort(weightsvm,'descend');
            size(SelectFeatureIndex)
            size(IX)
            SelectFeatureIndex(:,i)=IX;
            for j=1:length(SelectFeatureNum)
                 weight=[];
                kWeight=SelectFeatureNum(j);
                FFData=train_data(:,IX(1:kWeight));
                Data_tes=test_data(:,IX(1:kWeight));
                model = svmtrain(train_label,FFData,modelParemter);
                FeatureAll(1:kWeight,j,i)=1;
                if length(ReLabel)==2 && modelParemter(ID+3)
                    weight=model.SVs'*model.sv_coef;
                    weightall(1:kWeight,j,i)=weight;
                end
                [predict_label, acc, dec_values] =svmpredict(test_label,Data_tes, model);                
                
                for IDRelabel=1:length(ReLabel)
                    Specif_Sensit_acc_fold(i,IDRelabel,j)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
                end
                Specif_Sensit_acc_fold(i,length(ReLabel)+1,j)=length(find(predict_label==test_label))/length(test_label)*100;
                
                accuracy(i,j)=acc(1);
                
                PredictLabel(find(Fold==i),j)=predict_label;
                RealLabel(find(Fold==i),j)=test_label;
                DecValues(find(Fold==i),:,j)=dec_values;

            end
        case 'F-Score'           
            [IX,ind] = fget(train_label,train_data);
            SelectFeatureIndex(:,i)=IX;
            for j=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(j);
                Data_testrank=test_data(:,IX(1:kWeight));
                model = svmtrain(train_label,train_data(:,IX(1:kWeight)),modelParemter);
                FeatureAll(1:kWeight,j,i)=1;
                if length(ReLabel)==2 && modelParemter(ID+3)
                    weight=[];
                    weight=model.SVs'*model.sv_coef;%compute weight
                    weightall(1:kWeight,j,i)=weight;%all weight
%                     FeatureAll(IX(1:kWeight),j,i)=1;
                end
                
                [predict_label, acc, dec_values] =svmpredict(test_label,Data_testrank, model);
                
                accuracy(i,j)=acc(1);
                for IDRelabel=1:length(ReLabel)
                    Specif_Sensit_acc_fold(i,IDRelabel,j)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
                end
                %acc of each fold
                Specif_Sensit_acc_fold(i,length(ReLabel)+1,j)=length(find(predict_label==test_label))/length(test_label)*100;
                
                PredictLabel(find(Fold==i),j)=predict_label;
                RealLabel(find(Fold==i),j)=test_label;
                DecValues(find(Fold==i),:,j)=dec_values;
            end
        case 'Relieff'
            [IX, w] = reliefF( train_data,train_label,size(train_data,2));
            SelectFeatureIndex(:,i)=IX;
            for j=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(j);
                Data_testrank=test_data(:,IX(1:kWeight));
                FeatureAll(1:kWeight,j,i)=1;
                model = svmtrain(train_label,train_data(:,IX(1:kWeight)),modelParemter);
                if length(ReLabel)==2 && modelParemter(ID+3)
                    weight=[];
                    weight=model.SVs'*model.sv_coef;
                    weightall(1:kWeight,j,i)=weight;
                end
                [predict_label, acc, dec_values] =svmpredict(test_label,Data_testrank, model);
                
                accuracy(i,j)=acc(1);
                for IDRelabel=1:length(ReLabel)
                    Specif_Sensit_acc_fold(i,IDRelabel,j)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
                end
                %acc of each fold
                Specif_Sensit_acc_fold(i,length(ReLabel)+1,j)=length(find(predict_label==test_label))/length(test_label)*100;                
                PredictLabel(find(Fold==i),j)=predict_label;
                RealLabel(find(Fold==i),j)=test_label;
                DecValues(find(Fold==i),:,j)=dec_values;
            end
        case 'Mrmr'          
            IX = mRMR(train_data,train_label,size(train_data,2));
            SelectFeatureIndex(:,i)=IX;
            for j=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(j);
                Data_testrank=test_data(:,IX(1:kWeight));
                model = svmtrain(train_label,train_data(:,IX(1:kWeight)),modelParemter);
                FeatureAll(1:kWeight,j,i)=1;
                if length(ReLabel)==2 && modelParemter(ID+3)
                    weight=[];
                    weight=model.SVs'*model.sv_coef;
                    weightall(1:kWeight,j,i)=weight;                    
                end
                [predict_label, acc, dec_values] =svmpredict(test_label,Data_testrank, model);
                
                accuracy(i,j)=acc(1);
                for IDRelabel=1:length(ReLabel)
                    Specif_Sensit_acc_fold(i,IDRelabel,j)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
                end
                %acc of each fold
                Specif_Sensit_acc_fold(i,length(ReLabel)+1,j)=length(find(predict_label==test_label))/length(test_label)*100;
                
                PredictLabel(find(Fold==i),j)=predict_label;
                RealLabel(find(Fold==i),j)=test_label;
                DecValues(find(Fold==i),:,j)=dec_values;
            end
    end
end
close(h)

if strcmp(strSelectFeature,'None')
    for IDRelabel=1:length(ReLabel)
        SpecificitySensitivity(IDRelabel)=length(find(RealLabel(find(RealLabel==ReLabel(IDRelabel)))==PredictLabel(find(RealLabel==ReLabel(IDRelabel)))))/length(find(RealLabel==ReLabel(IDRelabel)))*100;
    end
    Predict_test_Label=[RealLabel PredictLabel DecValues];
else
    for j=1:size(RealLabel,2) %j means the times of feature selecting
        for IDRelabel=1:length(ReLabel)
            SpecificitySensitivity(IDRelabel,j)=length(find(RealLabel(find(RealLabel(:,j)==ReLabel(IDRelabel)),j)==PredictLabel(find(RealLabel(:,j)==ReLabel(IDRelabel)),j)))/length(find(RealLabel(:,j)==ReLabel(IDRelabel)))*100;
        end
        Predict_test_Label(:,:,j)=[RealLabel(:,j) PredictLabel(:,j) DecValues(:,:,j)];
    end
end