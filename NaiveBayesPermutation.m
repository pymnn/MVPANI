function [accuracy,DecValues]=NaiveBayesPermutation(Cfg)
strSelectFeature=Cfg.FeatureSelectionStr;
FeatureNum=Cfg.FeatureNum;
FeatureNum=eval(['[',FeatureNum,']']);
Fold=Cfg.Fold;
label=Cfg.Label;
Data=Cfg.Data;
FishZStr=Cfg.FishZStr;
strFeatureReduce=Cfg.FeatureReduceStr;
PermutationNum=Cfg.PermutationNum;
randLabel=Cfg.RandLabel;
DecValues=[];
nDim=Cfg.ContrPCA;

ReLabel=unique(label);%the real label

if isempty(gcp('nocreate'))
    NumCore=inputdlg('the number of parallel core£º','Core',[1 40],{'4'});
    if length(NumCore)
        NumCore=str2num(cell2mat(NumCore));
        parpool('local',NumCore);
    else
        error('error:please click < ok > command');
    end
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
class=[];
parfor iper=1:PermutationNum
    TemAcc=[];
    TemDecValues=[];
    permlabel=randLabel(:,iper);
    for i=1:max(Fold)    
        test_data=Data(find(Fold==i),:);
        test_label=label(find(Fold==i));
        train_label=permlabel;
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
        index1=[];
        index2=[];
        for jNB=size(train_data,2):-1:1
            index1=find(train_label==ReLabel(1));
            index2=find(train_label==ReLabel(2));
            if var(train_data(index1,jNB))==0 || var(train_data(index2,jNB))==0
                train_data(:,jNB)=[];
                test_data(:,jNB)=[];
            end
        end
        switch strSelectFeature
            case 'None'%
                Factor = NaiveBayes.fit(train_data, train_label);
                [dec_values,predict_label] = posterior(Factor, test_data);
                acc= length(find(predict_label== test_label))/length(test_label)*100;
                TemAcc(i)=acc;
                TemDecValues(1,find(Fold==i))=dec_values(:,1);
            case 'Lasso'
                %%
                %20200113 edit
                [B,FitInfo]=lasso(train_data,train_label,'CV',10);
                IndexMinMSE=FitInfo.IndexMinMSE;
                train_data=train_data(:,find(B(:,IndexMinMSE)~=0));
                test_data=test_data(:,find(B(:,IndexMinMSE)~=0));
                
                model = NaiveBayes.fit(train_data, train_label);
                [dec_values,predict_label] = posterior(model, test_data); 
                
                acc= length(find(predict_label== test_label))/length(test_label)*100;
                TemAcc(i)=acc;
                TemDecValues(1,find(Fold==i))=dec_values(:,1);
            case 'F-Score'       
                [IX,ind] = fget(train_label,train_data);
                SelectFeatureNum=floor(FeatureNum*size(train_data,2));
                for j=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(j);
                    Frank=train_data(:,IX(1:kWeight));
                    FC_tes=test_data(:,IX(1:kWeight));
                    model = NaiveBayes.fit(Frank, train_label);
                    [dec_values,predict_label] = posterior(model, FC_tes); 
                    acc= length(find(predict_label== test_label))/length(test_label)*100;
                    TemAcc(j,i)=acc(1);
                    TemDecValues(j,find(Fold==i))=dec_values(:,1);
                end
            case 'Relieff'
                [IX, ind] = reliefF( train_data,train_label,size(train_data,2));
                 SelectFeatureNum=floor(FeatureNum*size(train_data,2));
                for j=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(j);
                    Frank=train_data(:,IX(1:kWeight));
                    FC_tes=test_data(:,IX(1:kWeight));
                    
                    Factor = NaiveBayes.fit(train_data(:,IX(1:kWeight)), train_label);
                    [dec_values,predict_label] = posterior(Factor, FC_tes);
                    
                    acc= length(find(predict_label== test_label))/length(test_label)*100;
                    TemAcc(j,i)=acc(1);
                    %acc of each fold
                    TemDecValues(j,find(Fold==i),:)=dec_values;
                end
            case 'Mrmr'
                IX = mRMR(train_data,train_label,size(train_data,2));
                 SelectFeatureNum=floor(FeatureNum*size(train_data,2));
                for j=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(j);
                    Frank=train_data(:,IX(1:kWeight));
                    FC_tes=test_data(:,IX(1:kWeight));
                    
                    Factor = NaiveBayes.fit(train_data(:,IX(1:kWeight)), train_label);
                    [dec_values,predict_label] = posterior(Factor, FC_tes);
                    
                    acc= length(find(predict_label== test_label))/length(test_label)*100;
                    TemAcc(j,i)=acc(1);
                    %acc of each fold
                    TemDecValues(j,find(Fold==i),:)=dec_values;
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