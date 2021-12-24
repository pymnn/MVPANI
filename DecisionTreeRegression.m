function [weightall,FeatureAll,IX,rr,pp,Predict_test_Label,MaxDim]=DecisionTreeRegression(Cfg)
%Use Decision Tree do regression
Fold=Cfg.Fold;
Label=Cfg.Label;
Data=Cfg.Data;
strSelectFeature=Cfg.FeatureSelectionStr;
SelectFeatureNum=Cfg.SelectFeatureNum;
FishZStr=Cfg.FishZStr;
nDim=Cfg.ContrPCA;
strFeatureReduce=Cfg.FeatureReduceStr;
SelectFeatureIndex=[];
weightall=[];
FeatureAll=[];
if ~strcmp(strSelectFeature,'None') && ~strcmp(strFeatureReduce,'PCA (Principal Component Analysis)')
    %Feature selection
    weightall=zeros(size(Data,2),length(SelectFeatureNum),max(Fold));
    FeatureAll=zeros(size(Data,2),length(SelectFeatureNum),max(Fold));
    SelectFeatureIndex=zeros(size(Data,2),max(Fold));
    pptesty=zeros(size(Fold,1),length(SelectFeatureNum));
    testyy=zeros(size(Fold,1),length(SelectFeatureNum));
else
    pptesty=zeros(size(Fold,1),1);
    testyy=zeros(size(Fold,1),1);
end
IX=[];
if strcmp(strFeatureReduce,'PCA (Principal Component Analysis)')
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
else
    MaxDim='NODimensionalityReduction';
end
h=waitbar(0,'Runing...');
for i=1:max(Fold)
    %%
    str=['Runing...',num2str(i/max(Fold)*100),'%'];
    waitbar(i/max(Fold),h,str);
    weightt=[];
    test_data=Data(find(Fold==i),:);
    test_label=Label(find(Fold==i),1);
    train_label=Label(:,1);
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
            model = fitrtree(train_data, train_label);
            [predict_label,dec_values] =predict(model, test_data);
            pptesty(find(Fold==i),1)=predict_label;
            testyy(find(Fold==i),1)=test_label;
        case 'Lasso'
            %%
            %20200113 edit
            [B,FitInfo]=lasso(train_data,train_label,'CV',10);
            IndexMinMSE=FitInfo.IndexMinMSE;
            train_data=train_data(:,find(B(:,IndexMinMSE)~=0));
            model = fitctree(train_data, train_label);
            test_data=test_data(:,find(B(:,IndexMinMSE)~=0));
            [predict_label,dec_values] =predict(model, test_data);
            
            pptesty(find(Fold==i))=predict_label;
            testyy(find(Fold==i))=test_label;
        case 'F-Score'
            weightt=[];
            [IX,ind] = fget(train_label,train_data);
            SelectFeatureIndex(:,i)=IX;
            for jNum=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(jNum);
                test_datarank=test_data(:,IX(1:kWeight));
                model = fitrtree(train_data(:,IX(1:kWeight)), train_label);
                [predict_label,dec_values] =predict(model, test_datarank);
                pptesty(find(Fold==i),jNum)=predict_label;
                testyy(find(Fold==i),jNum)=test_label;
            end
        case 'Mrmr'
            IX = mRMR(train_data,train_label,size(train_data,2));
            SelectFeatureIndex(:,i)=IX;
            for jNum=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(jNum);
                FeatureAll(1:kWeight,jNum,i)=1;
                test_datarank=test_data(:,IX(1:kWeight));
                model = fitrtree(train_data(:,IX(1:kWeight)), train_label);
                [predict_label,dec_values] =predict(model, test_datarank);
                pptesty(find(Fold==i),jNum)=predict_label;
                testyy(find(Fold==i),jNum)=test_label;
            end
        case 'Relieff'
            weightt=[];
            [IX, w] = reliefF( train_data,train_label,size(train_data,2));
            SelectFeatureIndex(:,i)=IX;
            for jNum=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(jNum);
                FeatureAll(1:kWeight,jNum,i)=1;
                test_datarank=test_data(:,IX(1:kWeight));
                model = fitrtree(train_data(:,IX(1:kWeight)), train_label);
                [predict_label,dec_values] =predict(model, test_datarank);
                pptesty(find(Fold==i),jNum)=predict_label;
                testyy(find(Fold==i),jNum)=test_label;
            end
    end
end
close(h)
%%
if strcmp( strSelectFeature,'None')
    [r,p]=corrcoef(testyy(:,1),pptesty(:,1));
    r1=r(1,2);
    p1=p(1,2);
    Predict_test_Label=[pptesty testyy];
else
    for k=1:size(pptesty,2)
        [r,p]=corrcoef(testyy(:,k),pptesty(:,k));
        r1(k)=r(1,2);
        p1(k)=p(1,2);
        Predict_test_Label(:,:,k)=[pptesty(:,k) testyy(:,k)];
    end
end
rr=r1';
pp=p1';
