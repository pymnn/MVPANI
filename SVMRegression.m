function [weightall,FeatureAll,SelectFeatureIndex,rr,pp,Predict_test_Label,MaxDim]=SVMRegression(Cfg)
%Tt is use data (feature) to predict the label using svm method.
%MaxDim is remainded number of the data after PCA.
Fold=Cfg.Fold;
Label=Cfg.Label;
Data=Cfg.Data;
modelParemter=Cfg.modelParemter;
SelectFeatureNum=Cfg.SelectFeatureNum;
strSelectFeature=Cfg.FeatureSelectionStr;
FishZStr=Cfg.FishZStr;
strFeatureReduce=Cfg.FeatureReduceStr;
nDim=Cfg.ContrPCA;
SelectFeatureIndex=[];
weightall=[];
FeatureAll=[];
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
        pptesty=zeros(size(Fold,1),length(SelectFeatureNum));
        testyy=zeros(size(Fold,1),length(SelectFeatureNum));
    else
        %Feature selection
        weightall=zeros(size(Data,2),max(Fold));
        SelectFeatureIndex=[];
        pptesty=zeros(size(Fold,1),1);
        testyy=zeros(size(Fold,1),1);
    end
else
    weightall=[];
end
h=waitbar(0,'Runing...');
for i=1:max(Fold)
    %%
    str=['Runing...',num2str(i/max(Fold)*100),'%'];
    waitbar(i/max(Fold),h,str);
    IX=[];
    weight=[];
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
    model = svmtrain(train_label,train_data,modelParemter);
    weight=model.SVs'*model.sv_coef;
    %%
    switch strSelectFeature
        case 'None'%
            weightall(:,i)=weight;
            [predict_label,acc, dec_values] = svmpredict(test_label,test_data,model);
            pptesty(find(Fold==i))=predict_label;
            testyy(find(Fold==i))=test_label;
        case 'Lasso'
            %%
            %20200113 edit
            [B,FitInfo]=lasso(train_data,train_label,'CV',10);
            IndexMinMSE=FitInfo.IndexMinMSE;
            train_data=train_data(:,find(B(:,IndexMinMSE)~=0));
            model = svmtrain(train_label,train_data,modelParemter);
            weight=[];
            weight=model.SVs'*model.sv_coef;%compute weight
            weightall(find(B(:,IndexMinMSE)~=0),i)=weight; %all weight
            test_data=test_data(:,find(B(:,IndexMinMSE)~=0));
            [predict_label, acc, dec_values] =svmpredict(test_label,test_data, model);
            pptesty(find(Fold==i))=predict_label;
            testyy(find(Fold==i))=test_label;
        case 'Weight'
            %%
            weightsvm=abs(weight);
            [B,IX]=sort(weightsvm,'descend');
            SelectFeatureIndex(:,i)=IX;
            weight=[];
            for jNum=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(jNum);
                Ftrain_data=train_data(:,IX(1:kWeight));
                Data_tes=test_data(:,IX(1:kWeight));
                model = svmtrain(train_label,Ftrain_data,modelParemter);
                weight=model.SVs'*model.sv_coef;
                weightall(1:kWeight,jNum,i)=weight;
                [predict_label, acc, dec_values] =svmpredict(test_label,Data_tes, model);
                pptesty(find(Fold==i),jNum)=predict_label;
                testyy(find(Fold==i),jNum)=test_label;
            end
            
        case 'F-Score'
            weight=[];
            [IX,ind] = fget(train_label,train_data);
            Frank=train_data(:,IX);
            SelectFeatureIndex(:,i)=IX;
            for jNum=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(jNum);
                test_datarank=test_data(:,IX);
                model = svmtrain(train_label,Frank(:,1:kWeight),modelParemter);
                weight=model.SVs'*model.sv_coef;
                weightall(1:kWeight,jNum,i)=weight;
                [predict_label, acc, dec_values] =svmpredict(test_label,test_datarank(:,1:kWeight), model);
                pptesty(find(Fold==i),jNum)=predict_label;
                testyy(find(Fold==i),jNum)=test_label;
            end
        case 'Relieff'
            [IX, w] = reliefF( train_data,train_label,size(train_data,2));
            Frank=train_data(:,IX);
            SelectFeatureIndex(:,i)=IX;
            for jNum=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(jNum);
                FeatureAll(1:kWeight,jNum,i)=1;
                test_datarank=test_data(:,IX);
                model = svmtrain(train_label,Frank(:,1:kWeight),modelParemter);
                weight=model.SVs'*model.sv_coef;
                weightall(1:kWeight,jNum,i)=weight;
                [predict_label, acc, dec_values] =svmpredict(test_label,test_datarank(:,1:kWeight), model);
                pptesty(find(Fold==i),jNum)=predict_label;
                testyy(find(Fold==i),jNum)=test_label;
            end
        case 'Mrmr'
            IX = mRMR(train_data,train_label,size(train_data,2));
            Frank=train_data(:,IX);
            SelectFeatureIndex(:,i)=IX;
            for jNum=1:length(SelectFeatureNum)
                kWeight=SelectFeatureNum(jNum);
                FeatureAll(1:kWeight,jNum,i)=1;
                test_datarank=test_data(:,IX);
                model = svmtrain(train_label,Frank(:,1:kWeight),modelParemter);
                weight=model.SVs'*model.sv_coef;
                weightall(1:kWeight,jNum,i)=weight;
                [predict_label, acc, dec_values] =svmpredict(test_label,test_datarank(:,1:kWeight), model);
                pptesty(find(Fold==i),jNum)=predict_label;
                testyy(find(Fold==i),jNum)=test_label;
            end
    end
end
close(h)
%%
if strcmp( strSelectFeature,'None')
    [r,p]=corrcoef(testyy,pptesty);
    r1=r(1,2);
    p1=p(1,2);
    Predict_test_Label=[testyy pptesty];
else
    for k=1:size(pptesty,2)
        [r,p]=corrcoef(testyy(:,k),pptesty(:,k));
        r1(k)=r(1,2);
        p1(k)=p(1,2);
        Predict_test_Label(:,:,k)=[testyy(:,k) pptesty(:,k) ];
    end
end
rr=r1';
pp=p1';

