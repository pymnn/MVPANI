function [accuracy,DecValues]=logisticPermutation(Cfg)
strSelectFeature=Cfg.FeatureSelectionStr;
FishZStr=Cfg.FishZStr;
Fold=Cfg.Fold;
Label=Cfg.Label;

ReLabel=unique(Label);%the real label

Data=Cfg.Data;
PermutationNum=Cfg.PermutationNum;

nDim=Cfg.ContrPCA;
strFeatureReduce=Cfg.FeatureReduceStr;
SelectFeatureNum=Cfg.SelectFeatureNum;

rand_label=Cfg.RandLabel;
RandLabel=rand_label;
RandLabel(find(rand_label==max(ReLabel)))=1;
RandLabel(find(rand_label==min(ReLabel)))=0;

if isempty(gcp('nocreate'))
    NumCore=inputdlg('the number of parallel core£º','Core',[1 40],{'4'});
    if length(NumCore)
        NumCore=str2num(cell2mat(NumCore));
        parpool('local',NumCore);
    else
        error('error:please click < ok > command');
    end
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
    TemAcc=[];
    TemDecValues=[];
    for i=1:max(Fold)
        predict=[];
        test_data=Data(find(Fold==i),:);
        test_label=Label(find(Fold==i));
        train_label=RandLabel(:,iper);
        
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
                theta = glmfit(train_data, [train_label ones(size(train_label,1),1)], 'binomial', 'link', 'logit');
                dec_values = glmval(theta,test_data, 'logit');
                predict(dec_values>0.5)=1;
                predict(dec_values<0.5)=0;
                predict_label=predict';
                
                acc= length(find(predict_label== test_label))/length(test_label)*100;
                TemAcc(i)=acc;
                TemDecValues(find(Fold==i),1)=dec_values(:,1);
            case 'Lasso'
                %%
                %20200113 edit
                [B,FitInfo]=lasso(train_data,train_label,'CV',10);
                IndexMinMSE=FitInfo.IndexMinMSE;
                train_data=train_data(:,find(B(:,IndexMinMSE)~=0));
                test_data=test_data(:,find(B(:,IndexMinMSE)~=0));
                
                model = glmfit(train_data, [train_label ones(size(train_label,1),1)], 'binomial', 'link', 'logit');
                dec_values = glmval(model,test_data, 'logit');
                predict(dec_values>0.5)=1;
                predict(dec_values<0.5)=0;
                predict_label=predict';
                TemAcc(i) = size(find(predict_label==(test_label)),1)*100/size(test_data,1);
                TemDecValues(find(Fold==i),1)=dec_values(:,1);
            case 'F-Score'
                [IX,ind] = fget(train_label,train_data);
                for j=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(j);
                    Frank=train_data(:,IX(1:kWeight));
                    Data_testrank=test_data(:,IX(1:kWeight));
                    theta = glmfit(Frank, [train_label ones(size(train_label,1),1)], 'binomial', 'link', 'logit');
                    dec_values = glmval(theta,Data_testrank, 'logit');
                    predict(dec_values>0.5)=1;
                    predict(dec_values<0.5)=0;
                    predict_label=predict';
                    acc= length(find(predict_label== test_label))/length(test_label)*100;
                    TemAcc(j,i)=acc;
                    TemDecValues(find(Fold==i),1)=dec_values(:,1);
                end
                
            case 'Relieff'
                [IX, ind] = reliefF( train_data,train_label,size(train_data,2));
                for j=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(j);
                    Frank=train_data(:,IX(1:kWeight));
                    Data_testrank=test_data(:,IX(1:kWeight));
                    theta = glmfit(Frank, [train_label ones(size(train_label,1),1)], 'binomial', 'link', 'logit');
                    dec_values = glmval(theta,Data_testrank, 'logit');
                    predict(dec_values>0.5)=1;
                    predict(dec_values<0.5)=0;
                    predict_label=predict';
                    acc= length(find(predict_label== test_label))/length(test_label)*100;
                    TemAcc(j,i)=acc;
                    TemDecValues(find(Fold==i),1)=dec_values(:,1);
                end
            case 'Mrmr'
                [IX,ind] = fget(train_label,train_data);
                IX = mRMR(train_data,train_label,size(train_data,2));
                for j=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(j);
                    Frank=train_data(:,IX(1:kWeight));
                    Data_testrank=test_data(:,IX(1:kWeight));
                    theta = glmfit(Frank, [train_label ones(size(train_label,1),1)], 'binomial', 'link', 'logit');
                    dec_values = glmval(theta,Data_testrank, 'logit');
                    predict(dec_values>0.5)=1;
                    predict(dec_values<0.5)=0;
                    predict_label=predict';
                    acc= length(find(predict_label== test_label))/length(test_label)*100;
                    TemAcc(j,i)=acc;
                    TemDecValues(find(Fold==i),1)=dec_values(:,1);
                end
        end
    end
    
    switch strSelectFeature
        case {'None','Lasso'}%
            class(iper)=mean(TemAcc);
        case {'F-Score','Relieff','Mrmr'}
            class1(:,iper)=mean(TemAcc,2);
    end
    DecValues(:,:,iper)=TemDecValues;
end
switch strSelectFeature
    case {'None','Lasso'}%
        accuracy=class;
    case {'F-Score','Relieff','Mrmr'}
        accuracy=class1;
end
if ~isempty(gcp('nocreate'))
    delete(gcp);
end