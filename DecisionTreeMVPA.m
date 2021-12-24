function [accuracy, Predict_test_Label,SpecificitySensitivity,weightall,FeatureAll,SelectFeatureIndex,Specif_Sensit_acc_fold,MaxDim]=DecisionTreeMVPA(Cfg)
%Use Decision Tree do classify
Fold=Cfg.Fold;
strSelectFeature=Cfg.FeatureSelectionStr;
Data=Cfg.Data;
SelectFeatureNum=Cfg.SelectFeatureNum;
% if ~strcmp(strSelectFeature,'None')
%     weightall=zeros(SelectFeatureNum(length(SelectFeatureNum)),length(SelectFeatureNum),max(Fold));
% else
%     weightall=zeros(size(Data,2),max(Fold));
% end


Label=Cfg.Label;
FishZStr=Cfg.FishZStr;
nDim=Cfg.ContrPCA;
strFeatureReduce=Cfg.FeatureReduceStr;


%%
%20200113 edit
if ~strcmp(strFeatureReduce,'PCA (Principal Component Analysis)')
    if ~strcmp(strSelectFeature,'None') && ~strcmp(strSelectFeature,'Lasso')
        %Feature selection
        weightall=zeros(size(Data,2),length(SelectFeatureNum),max(Fold));
        FeatureAll=zeros(size(Data,2),length(SelectFeatureNum),max(Fold));
        SelectFeatureIndex=zeros(size(Data,2),max(Fold));
    else
        %Feature selection
        weightall=zeros(size(Data,2),max(Fold));
        FeatureAll=zeros(size(Data,2),max(Fold));
        SelectFeatureIndex=[];
    end
else
    weightall=[];
    FeatureAll=[];
    SelectFeatureIndex=[];
end
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
ReLabel=[];
FData=[];
ReLabel=unique(Label);
h=waitbar(0,'Runing...');
for i=1:max(Fold)
    str=['Runing...',num2str(i/max(Fold)*100),'%'];
    waitbar(i/max(Fold),h,str);
    IX=[];
    test_data=Data(find(Fold==i),:);
    test_label=Label(find(Fold==i));
    train_label=Label;
    train_data=Data;
    train_label(find(Fold==i),:)=[];
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
    %%
    switch strSelectFeature
        case 'None'%
            %%
            model = fitctree(train_data, train_label);
            [predict_label,dec_values] =predict(model, test_data);
            for IDRelabel=1:length(ReLabel)
                Specif_Sensit_acc_fold(i,IDRelabel,1)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
            end
            Specif_Sensit_acc_fold(i,length(ReLabel)+1,1)=length(find(predict_label==test_label))/length(test_label)*100;
            acc= length(find(predict_label== test_label))/length(test_label)*100;
            accuracy(i)=acc;
            
            PredictLabel(find(Fold==i),1)=predict_label;
            RealLabel(find(Fold==i),1)=test_label;
            DecValues(find(Fold==i),:,1)=dec_values;
        case 'Lasso'
            %%
            %20200113 edit
            [B,FitInfo]=lasso(train_data,train_label,'CV',10);
            IndexMinMSE=FitInfo.IndexMinMSE;
            train_data=train_data(:,find(B(:,IndexMinMSE)~=0));
            model = fitctree(train_data, train_label);
            test_data=test_data(:,find(B(:,IndexMinMSE)~=0));
            [predict_label,dec_values] =predict(model, test_data);
            
            acc= length(find(predict_label== test_label))/length(test_label)*100;
            accuracy(i)=acc;
            
            for IDRelabel=1:length(ReLabel)
                Specif_Sensit_acc_fold(i,IDRelabel,1)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
            end
            Specif_Sensit_acc_fold(i,length(ReLabel)+1,1)=length(find(predict_label==test_label))/length(test_label)*100;
            
            PredictLabel(find(Fold==i),1)=predict_label;
            RealLabel(find(Fold==i),1)=test_label;
            DecValues(find(Fold==i),:,1)=dec_values;
            
        case 'F-Score'
            [IX,ind] = fget(train_label,train_data);
            SelectFeatureIndex(:,i)=IX;
            for j=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(j);
                Data_testrank=test_data(:,IX(1:kWeight));
                FeatureAll(1:kWeight,j,i)=1;
                model = fitctree(train_data(:,IX(1:kWeight)),train_label);
                [predict_label,dec_values] =predict(model,Data_testrank);
                for IDRelabel=1:length(ReLabel)
                    Specif_Sensit_acc_fold(i,IDRelabel,j)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
                end
                Specif_Sensit_acc_fold(i,length(ReLabel)+1,j)=length(find(predict_label==test_label))/length(test_label)*100;
                acc= length(find(predict_label== test_label))/length(test_label)*100;
                accuracy(i,j)=acc;
                
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
                model = fitctree(train_data(:,IX(1:kWeight)),train_label);
                [predict_label,dec_values] =predict(model,Data_testrank);
                acc= length(find(predict_label==test_label))/length(test_label)*100;
                accuracy(i,j)=acc;

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
                FeatureAll(1:kWeight,j,i)=1;
                model = fitctree(train_data(:,IX(1:kWeight)),train_label);
                [predict_label,dec_values] =predict(model,Data_testrank);
                acc= length(find(predict_label==test_label))/length(test_label)*100;
                accuracy(i,j)=acc;

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
%%
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
