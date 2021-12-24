function [accuracy,DecValues]=SVMPermutation(Cfg)
%using svm method to prediction, do permutation
%    Data: the data for prediction
%    label: the value is 1,-1 or 1,0
strSelectFeature=Cfg.FeatureSelectionStr;
SelectFeatureNum=Cfg.SelectFeatureNum;
Fold=Cfg.Fold;
label=Cfg.Label;
Data=Cfg.Data;
modelParemter=Cfg.modelParemter;
PermutationNum=Cfg.PermutationNum;
FishZStr=Cfg.FishZStr;
strFeatureReduce=Cfg.FeatureReduceStr;
nDim=Cfg.ContrPCA;
randLabel=Cfg.RandLabel;
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
            MaxDimPCA(i)=max(find(RContributionLatent(:,2)<nDim));
        end
    end
    MaxDim=max(MaxDimPCA);
else
    MaxDim='NODimensionalityReduction';
end

parfor iper=1:PermutationNum
    TemAcc=[];
    TemDecValues=[];
    permlabel=randLabel(:,iper);
    %     h(iper)=waitbar(0,'Runing...')
    for i=1:max(Fold)
        %         str=['Runing...',num2str(i/max(Fold)*100),'%'];
        %         waitbar(i/max(Fold),h(iper),str);
        test_data=Data(find(Fold==i),:);
        test_label=label(find(Fold==i));
        train_data=Data;
        train_label=permlabel;
        train_label(find(Fold==i),:)=[];%delete test data,residues will be used as train data
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
        weightt=model.SVs'*model.sv_coef;%compute weight
        switch strSelectFeature
            case 'None'%
                model = svmtrain(train_label,train_data,modelParemter);
                [predict_label, accuracy, dec_values] =svmpredict(test_label,test_data, model);
                TemAcc(i)=accuracy(1);
                TemDecValues(1,find(Fold==i))=dec_values(:,1);
            case 'Lasso'
                %%
                %20200113 edit
                [B,FitInfo]=lasso(train_data,train_label,'CV',10);
                IndexMinMSE=FitInfo.IndexMinMSE;
                train_data=train_data(:,find(B(:,IndexMinMSE)~=0));
                model = svmtrain(train_label,train_data,modelParemter);
                weight=model.SVs'*model.sv_coef;%compute weight
                test_data=test_data(:,find(B(:,IndexMinMSE)~=0));
                [predict_label, acc, dec_values] =svmpredict(test_label,test_data, model);
                TemAcc(i)=acc(1);
                TemDecValues(1,find(Fold==i))=dec_values(:,1);
            case 'Weight'
                weightsvm=abs(weightt);
                [B,IX]=sort(weightsvm,'descend');
                for j=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(j);
                    FF=train_data(:,IX(1:kWeight));
                    FC_tes=test_data(:,IX(1:kWeight));
                    model = svmtrain(train_label,FF,modelParemter);
                    [predict_label, accuracy, dec_values] =svmpredict(test_label,FC_tes, model);
                    TemAcc(j,i)=accuracy(1);
                    TemDecValues(j,find(Fold==i),:)=dec_values;
                    %j means the number of feature,iper means the number of permutation
                    %i means fold
                end
            case 'F-Score'
                [IX,ind] = fget(train_label,train_data);
                for j=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(j);
                    Frank=train_data(:,IX(1:kWeight));
                    FC_tes=test_data(:,IX(1:kWeight));
                    model = svmtrain(train_label,Frank,modelParemter);
                    [predict_label, accuracy, dec_values] =svmpredict(test_label,FC_tes, model);
                    TemAcc(j,i)=accuracy(1);
                    %j means the number of feature,iper means the number of permutation
                    %i means fold
                    TemDecValues(j,find(Fold==i),:)=dec_values;
                end
            case 'Relieff'
                [IX, ind] = reliefF( train_data,train_label,size(train_data,2));
                for j=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(j);
                    Data_testrank=test_data(:,IX(1:kWeight));
                    model = svmtrain(train_label,train_data(:,IX(1:kWeight)),modelParemter);
                    [predict_label, accuracy, dec_values] =svmpredict(test_label,Data_testrank, model);
                    TemAcc(j,i)=accuracy(1);
                    %acc of each fold
                    TemDecValues(j,find(Fold==i),:)=dec_values;
                end
            case 'Mrmr'
                IX = mRMR(train_data,train_label,size(train_data,2));
                for j=1:length(SelectFeatureNum)
                    kWeight=SelectFeatureNum(j);
                    Data_testrank=test_data(:,IX(1:kWeight));
                    model = svmtrain(train_label,train_data(:,IX(1:kWeight)),modelParemter);
                    [predict_label, accuracy, dec_values] =svmpredict(test_label,Data_testrank, model);
                    TemAcc(j,i)=accuracy(1);
                    %acc of each fold
                    TemDecValues(j,find(Fold==i),:)=dec_values;
                end
        end
    end
    % close(h(iper))
    switch strSelectFeature
        case {'None','Lasso'}%
            class(iper)=mean(TemAcc);
        case {'F-Score','Weight','Relieff','Mrmr'}
            class1(:,iper)=mean(TemAcc,2);
    end
    DecValues(:,:,iper)=TemDecValues;
end
% close(h)
switch strSelectFeature
    case {'None','Lasso'}%
        accuracy=class;
    case {'F-Score','Weight','Relieff','Mrmr'}
        accuracy=class1;
end
if ~isempty(gcp('nocreate'))
    delete(gcp);
end