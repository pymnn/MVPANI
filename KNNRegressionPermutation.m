function [rr,pp]=KNNRegressionPermutation(Cfg)
strSelectFeature=Cfg.FeatureSelectionStr;
SelectFeatureNum=Cfg.SelectFeatureNum;
Fold=Cfg.Fold;
label=Cfg.Label;
Data=Cfg.Data;
PermutationNum=Cfg.PermutationNum;
NumNeighbor=Cfg.NumNeighbor;
randLabel=Cfg.RandLabel;
FishZStr=Cfg.FishZStr;
nDim=Cfg.ContrPCA;
strFeatureReduce=Cfg.FeatureReduceStr;
if isempty(gcp('nocreate'))
    NumCore=inputdlg('the number of parallel core£º','Core',[1 40],{'4'});
    if length(NumCore)
        NumCore=str2num(cell2mat(NumCore));
        parpool('local',NumCore);
    else
        error('error:please click < ok > command');
    end
end

if ~strcmp(strSelectFeature,'None')
    pptesty=zeros(size(Fold,1),length(SelectFeatureNum));
    testyy=zeros(size(Fold,1),length(SelectFeatureNum));
else
    pptesty=zeros(size(Fold,1));
    testyy=zeros(size(Fold,1));
end
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
parfor iper=1:PermutationNum
    Tempptesty=[];
    Temptestyy=[];
    rtem=[];
    ptem=[];
    for i=1:max(Fold)
        test_data=Data(find(Fold==i),:);
        test_label=label(find(Fold==i));
        train_label=randLabel(:,iper);
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
        
        switch strSelectFeature
            case 'None'%
                Factor = ClassificationKNN.fit(train_data, train_label, 'NumNeighbors', NumNeighbor);
                [predict_label, dec_values] = predict(Factor, test_data);
                Tempptesty(find(Fold==i),1)=predict_label;
                Temptestyy(find(Fold==i),1)=test_label;
            case 'Lasso'
                %%
                %20200113 edit
                [B,FitInfo]=lasso(train_data,train_label,'CV',10);
                IndexMinMSE=FitInfo.IndexMinMSE;
                train_data=train_data(:,find(B(:,IndexMinMSE)~=0));
                test_data=test_data(:,find(B(:,IndexMinMSE)~=0));
                model = ClassificationKNN.fit(train_data, train_label, 'NumNeighbors', NumNeighbor);
                [predict_label, dec_values] = predict(model, test_data);
                
                Tempptesty(find(Fold==i),1)=predict_label;
                Temptestyy(find(Fold==i),1)=test_label;
            case 'F-Score'
                [IX,ind] = fget(train_label,train_data);
                for jNum=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(jNum);
                    Frank=train_data(:,IX(1:kWeight));
                    test_datarank=test_data(:,IX(1:kWeight));
                    Factor = ClassificationKNN.fit(Frank, train_label, 'NumNeighbors', NumNeighbor);
                    [predict_label, dec_values] = predict(Factor, test_datarank);
                    Tempptesty(find(Fold==i),jNum)=predict_label;
                    Temptestyy(find(Fold==i),jNum)=test_label;
                end
            case 'Relieff'
                [IX, w] = reliefF( train_data,train_label,size(train_data,2));
                for jNum=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(jNum);
                    Frank=train_data(:,IX(1:kWeight));
                    test_datarank=test_data(:,IX(1:kWeight));
                    Factor = ClassificationKNN.fit(Frank, train_label, 'NumNeighbors', NumNeighbor);
                    [predict_label, dec_values] = predict(Factor, test_datarank);
                    Tempptesty(find(Fold==i),jNum)=predict_label;
                    Temptestyy(find(Fold==i),jNum)=test_label;
                end
            case 'Mrmr'
                IX = mRMR(train_data,train_label,size(train_data,2));
                for jNum=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(jNum);
                    Frank=train_data(:,IX(1:kWeight));
                    test_datarank=test_data(:,IX(1:kWeight));
                    Factor = ClassificationKNN.fit(Frank, train_label, 'NumNeighbors', NumNeighbor);
                    [predict_label, dec_values] = predict(Factor, test_datarank);
                    Tempptesty(find(Fold==i),jNum)=predict_label;
                    Temptestyy(find(Fold==i),jNum)=test_label;
                end
        end
    end
    for k=1:size(Tempptesty,2)
        [r,p]=corrcoef(Temptestyy(:,k),Tempptesty(:,k));
        rtem(k)=r(1,2);
        ptem(k)=p(1,2);
    end
    rr(iper,:)=rtem;
    pp(iper,:)=ptem;
end
if ~isempty(gcp('nocreate'))
    delete(gcp);
end