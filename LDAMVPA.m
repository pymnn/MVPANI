function [accuracy, Predict_test_Label,SpecificitySensitivity,weightall,FeatureAll,SelectFeatureIndex,Specif_Sensit_acc_fold,MaxDim]=LDAMVPA(Cfg)
%ldatrain function in the 'decoding_toolbox_v3'
%by Martin Hebart & Kai Goergen was called in this function
%ldapredict function in the 'decoding_toolbox_v3'
%by Martin Hebart & Kai Goergen was called in this function
Fold=Cfg.Fold;
strSelectFeature=Cfg.FeatureSelectionStr;
Data=Cfg.Data;
FishZStr=Cfg.FishZStr;
nDim=Cfg.ContrPCA;
strFeatureReduce=Cfg.FeatureReduceStr;
SelectFeatureNum=Cfg.SelectFeatureNum;
if ~strcmp(strSelectFeature,'None')
    weightall=zeros(SelectFeatureNum(length(SelectFeatureNum)),length(SelectFeatureNum),max(Fold));
else
    weightall=zeros(size(Data,2),max(Fold));
end
IX=[];
Label=Cfg.Label;

param.shrinkage='none';   %set param of lda
ReLabel=[];
train_data=[];

ReLabel=unique(Label);%true label
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
FeatureAll=[];
if ~strcmp(strFeatureReduce,'PCA (Principal Component Analysis)')
    if ~strcmp(strSelectFeature,'None') && ~strcmp(strSelectFeature,'Lasso')
        %Feature selection
        weightall=zeros(size(Data,2),length(SelectFeatureNum),max(Fold));
        FeatureAll=zeros(size(Data,2),length(SelectFeatureNum),max(Fold));
        SelectFeatureIndex=zeros(size(Data,2),max(Fold));
    else
        %Feature selection
        weightall=zeros(size(Data,2),max(Fold));
        SelectFeatureIndex=[];
    end
else
    weightall=[];
    SelectFeatureIndex=[];
end
h=waitbar(0,'Runing...');
for i=1:max(Fold)
    str=['Runing...',num2str(i/max(Fold)*100),'%'];
    waitbar(i/max(Fold),h,str);
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
            model = ldatrain(train_label,train_data,param);
            [predict_label dec_values acc] = ldapredict(test_label,test_data,model);
            Specif_Sensit_acc_fold(i,1)=length(find(predict_label(find(predict_label==test_label))==ReLabel(1)))/length(find(test_label==ReLabel(1)))*100;
            Specif_Sensit_acc_fold(i,2)=length(find(predict_label(find(predict_label==test_label))==ReLabel(2)))/length(find(test_label==ReLabel(2)))*100;
            Specif_Sensit_acc_fold(i,3)=length(find(predict_label==test_label))/length(test_label)*100;
            acc=acc*100;
            accuracy(i)=acc;
            
            PredictLabel(find(Fold==i),1)=predict_label;
            RealLabel(find(Fold==i),1)=test_label;
            DecValues(find(Fold==i),1)=dec_values(:,1);
        case 'Lasso'
            %%
            %20200113 edit
            [B,FitInfo]=lasso(train_data,train_label,'CV',10);
            IndexMinMSE=FitInfo.IndexMinMSE;
            train_data=train_data(:,find(B(:,IndexMinMSE)~=0));
            test_data=test_data(:,find(B(:,IndexMinMSE)~=0));
            
            model = ldatrain(train_label,train_data,param);
            [predict_label dec_values acc] = ldapredict(test_label,test_data,model);
            
            acc=acc*100;
            accuracy(i)=acc;
            Specif_Sensit_acc_fold(i,1)=length(find(predict_label(find(predict_label==test_label))==ReLabel(1)))/length(find(test_label==ReLabel(1)))*100;
            Specif_Sensit_acc_fold(i,2)=length(find(predict_label(find(predict_label==test_label))==ReLabel(2)))/length(find(test_label==ReLabel(2)))*100;
            Specif_Sensit_acc_fold(i,3)=length(find(predict_label==test_label))/length(test_label)*100;
            
            PredictLabel(find(Fold==i),1)=predict_label;
            RealLabel(find(Fold==i),1)=test_label;
            DecValues(find(Fold==i),1)=dec_values(:,1);
        case 'F-Score'
            [IX,ind] = fget(train_label,train_data);
            SelectFeatureIndex(:,i)=IX;
            for j=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(j);
                Data_testrank=test_data(:,IX(1:kWeight));
                model = ldatrain(train_label,train_data(:,IX(1:kWeight)),param);
                [predict_label dec_values acc] = ldapredict(test_label,Data_testrank,model);
                
                for IDRelabel=1:length(ReLabel)
                    Specif_Sensit_acc_fold(i,IDRelabel,j)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
                end
                Specif_Sensit_acc_fold(i,length(ReLabel)+1,j)=length(find(predict_label==test_label))/length(test_label)*100;
                
                acc=acc*100;
                accuracy(i,j)=acc;
                
                PredictLabel(find(Fold==i),j)=predict_label;
                RealLabel(find(Fold==i),j)=test_label;
                DecValues(find(Fold==i),j)=dec_values(:,1);
            end
            
        case 'Relieff'
            [IX, w] = reliefF( train_data,train_label,size(train_data,2));
            SelectFeatureIndex(:,i)=IX;
            for j=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(j);
                Data_testrank=test_data(:,IX(1:kWeight));
                model = ldatrain(train_label,train_data(:,IX(1:kWeight)),param);
                [predict_label dec_values acc] = ldapredict(test_label,Data_testrank,model);
                
                for IDRelabel=1:length(ReLabel)
                    Specif_Sensit_acc_fold(i,IDRelabel,j)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
                end
                Specif_Sensit_acc_fold(i,length(ReLabel)+1,j)=length(find(predict_label==test_label))/length(test_label)*100;
                
                acc=acc*100;
                accuracy(i,j)=acc;
                
                PredictLabel(find(Fold==i),j)=predict_label;
                RealLabel(find(Fold==i),j)=test_label;
                DecValues(find(Fold==i),j)=dec_values(:,1);
            end
        case 'Mrmr'
            IX = mRMR(train_data,train_label,size(train_data,2));
            SelectFeatureIndex(:,i)=IX;
            for j=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(j);
                Data_testrank=test_data(:,IX(1:kWeight));
                model = ldatrain(train_label,train_data(:,IX(1:kWeight)),param);
                [predict_label dec_values acc] = ldapredict(test_label,Data_testrank,model);
                
                for IDRelabel=1:length(ReLabel)
                    Specif_Sensit_acc_fold(i,IDRelabel,j)=length(find(predict_label(find(predict_label==test_label))==ReLabel(IDRelabel)))/length(find(test_label==ReLabel(IDRelabel)))*100;
                end
                Specif_Sensit_acc_fold(i,length(ReLabel)+1,j)=length(find(predict_label==test_label))/length(test_label)*100;
                
                acc=acc*100;
                accuracy(i,j)=acc;
                
                PredictLabel(find(Fold==i),j)=predict_label;
                RealLabel(find(Fold==i),j)=test_label;
                DecValues(find(Fold==i),j)=dec_values(:,1);
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
