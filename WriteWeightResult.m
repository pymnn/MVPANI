function WriteWeightResult(WeightAll,IX,Cfg)
%Output WeightAll, WeightImgMean and WeightMean
%WeightAll means the weight of each feature in each fold
%WeightMean means the mean weight of each feature in all folds
%WeightImgMean means the mean weight of each feature in all folds in the
%mask(the mask may be image or matrix)
MaskPath=Cfg.MaskPath;
sel_CheckWeightAll=Cfg.sel_CheckWeightAll;
sel_CheckWeightI=Cfg.sel_CheckWeightI;
OutputDir=Cfg.OutputDir;
nDim1=Cfg.Mask.nDim1;
nDim2=Cfg.Mask.nDim2;
nDim3=Cfg.Mask.nDim3;
IndexColum=Cfg.IndexColum;
SelectFeatureNum=Cfg.SelectFeatureNum;
FeatureNum=Cfg.DataNum;
strFeatureReduce=Cfg.FeatureReduceStr;

if strcmp(strFeatureReduce,'PCA (Principal Component Analysis)')
    if sel_CheckWeightAll
        %MeanWeight
        MeanWeight=mean(WeightAll,2);
        save([OutputDir,filesep,'WeightFeature.mat'],'MeanWeight')
    end
    if sel_CheckWeightI
        %output weight of each fold
%         save([OutputDir,filesep,'WeightAll.mat'],'WeightAll');
        for i=1:size(WeightAll,2)
            WeightImg=WeightAll(:,i);
            save([OutputDir,filesep,'WeightFeatureFold_',num2str(i),'.mat'],'WeightImg')
        end
    end
else
    if isempty(MaskPath)
%         save([OutputDir,filesep,'WeightAll.mat'],'WeightAll');
        if sel_CheckWeightAll
            %MeanWeight
            if length(size(WeightAll))>2
                %had done feature selection
                WeightImg=zeros(FeatureNum,size(WeightAll,2),size(WeightAll,3));
                for i=1:size(WeightAll,2)
                    for j=1:size(WeightAll,3)
                        WeightImg(IndexColum(IX(1:SelectFeatureNum(i),j)),i,j)=WeightAll(1:SelectFeatureNum(i),i,j);
                    end
                end
                MeanWeight=mean(WeightImg,3);
                for iWeight=1:size(MeanWeight,2)
                    MeanWeightFeature=MeanWeight(:,iWeight);
                    save([OutputDir,filesep,'WeightSelectFeature_',num2str(iWeight),'.mat'],'MeanWeightFeature')
                end
            else
                %not feature selection
                WeightImg=zeros(FeatureNum,1);
                MeanWeight=mean(WeightAll,2);
                WeightImg(IndexColum)=MeanWeight;
                save([OutputDir,filesep,'WeightFeature.mat'],'WeightImg')
            end
        end
        if sel_CheckWeightI %output all weight
            if length(size(WeightAll))>2
                %had done feature selection
                WeightImg=zeros(FeatureNum,1);
                for i=1:size(WeightAll,2)
                    for j=1:size(WeightAll,3)
                        WeightImg(IndexColum(IX(1:SelectFeatureNum(i),j)))=WeightAll(1:SelectFeatureNum(i),i,j);
                        save([OutputDir,filesep,'WeightSelectFeature_',num2str(i),'_Fold_',num2str(j),'.mat'],'WeightImg')
                    end
                end
            else
                %not feature selection
                WeightImg=zeros(FeatureNum,1);
                for i=1:size(WeightAll,2)
                    WeightImg(IndexColum)=WeightAll(:,i);
                    save([OutputDir,filesep,'WeightFeatureFold_',num2str(i),'.mat'],'WeightImg')
                end
            end
        end
    else
        %MaskPath not empty
%         save([OutputDir,filesep,'WeightAll.mat'],'WeightAll');
        [Path,Name,fileType]=fileparts(MaskPath);
        if strcmp(fileType,'.img')|| strcmp(fileType,'.nii')
            %if the mask is image, the weight will be mapped in the brain
            if sel_CheckWeightAll
                %output the mean weight
                if length(size(WeightAll))>2
                    %had done feature selection
                    WeightImg=zeros(nDim1*nDim2*nDim3,size(WeightAll,2),size(WeightAll,3));
                    for i=1:size(WeightAll,2)
                        %i denote select feature,j denote fold
                        for j=1:size(WeightAll,3)
                            WeightImg(IndexColum(IX(1:SelectFeatureNum(i),j)),i,j)=WeightAll(1:SelectFeatureNum(i),i,j);
                        end
                    end
                    MeanWeight=mean(WeightImg,3);
                    for iWeight=1:size(MeanWeight,2)
                        MeanWeightFeature=MeanWeight(:,iWeight);
                        MeanWeightFeature=reshape(MeanWeightFeature,nDim1,nDim2,nDim3);
                        VMask=spm_vol(MaskPath);
                        v=spm_create_vol(VMask);
                        v.dt=[16 0];
                        v.fname=[OutputDir,filesep,'WeightImgSelectFeature_',num2str(iWeight),'.nii'];
                        spm_write_vol(v,MeanWeightFeature);
                        save([OutputDir,filesep,'WeightImgSelectFeature_',num2str(iWeight),'.mat'],'MeanWeightFeature');
                    end
                else
                    %No feature selection
                    WeightImg=zeros(1,nDim1*nDim2*nDim3);
                    MeanWeight=mean(WeightAll,2);
                    WeightImg(IndexColum)=MeanWeight;
                    WeightImg=reshape(WeightImg,nDim1,nDim2,nDim3);
                    VMask=spm_vol(MaskPath);
                    v=spm_create_vol(VMask);
                    v.dt=[16 0];
                    v.fname=[OutputDir,filesep,'WeightImgFeature.nii'];
                    spm_write_vol(v,WeightImg);
                    save([OutputDir,filesep,'WeightImgFeature.mat'],'WeightImg');
                end
            end
            if sel_CheckWeightI
                %output each fold weight
                %                 save([OutputDir,filesep,'WeightAll.mat'],'WeightAll');
                if length(size(WeightAll))>2
                    for i=1:size(WeightAll,2)
                        %i denote select feature,j denote fold
                        for j=1:size(WeightAll,3)
                            WeightImg=zeros(1,nDim1*nDim2*nDim3);
                            WeightImg(IndexColum(IX(1:SelectFeatureNum(i),j)))=WeightAll(1:SelectFeatureNum(i),i,j);
                            WeightImg=reshape(WeightImg,nDim1,nDim2,nDim3);
                            VMask=spm_vol(MaskPath);
                            v=spm_create_vol(VMask);
                            v.dt=[16 0];
                            v.fname=[OutputDir,filesep,'WeightImgSelectFeature_',num2str(i),'_Fold_',num2str(j),'.nii'];
                            save([OutputDir,filesep,'WeightImgSelectFeature_',num2str(i),'_Fold_',num2str(j),'.mat'],'WeightImg');
                            spm_write_vol(v,WeightImg);
                            clear WeightImg VMask v
                        end
                    end
                else
                    for i=1:size(WeightAll,2)
                        WeightImg=zeros(1,nDim1*nDim2*nDim3);
                        WeightImg(IndexColum)=WeightAll(:,i);
                        WeightImg=reshape(WeightImg,nDim1,nDim2,nDim3);
                        VMask=spm_vol(MaskPath);
                        v=spm_create_vol(VMask);
                        v.dt=[16 0];
                        v.fname=[OutputDir,filesep,'WeightImgFeatureFold_',num2str(i),'.nii'];
                        save([OutputDir,filesep,'WeightImgFeatureFold_',num2str(i),'.mat'],'WeightImg')
                        spm_write_vol(v,WeightImg);
                    end
                end
            end
        end
        if strcmp(fileType,'.txt')|| strcmp(fileType,'.mat')
            %if the mask is txt, the weight will be mapped in the txt
%             save([OutputDir,filesep,'WeighttAll.mat'],'WeightAll');
            if sel_CheckWeightAll
                if length(size(WeightAll))>2
                    %had done feature selection
                    WeightImg=zeros(nDim1*nDim2,size(WeightAll,2),size(WeightAll,3));
                    for i=1:size(WeightAll,2)
                        %i denote select feature,j denote fold
                        for j=1:size(WeightAll,3)
                            WeightImg(IndexColum(IX(1:SelectFeatureNum(i),j)),i,j)=WeightAll(1:SelectFeatureNum(i),i,j);
                        end
                    end
                    MeanWeight=mean(WeightImg,3);
                    for iWeight=1:size(MeanWeight,2)
                        MeanWeightFeature=MeanWeight(:,iWeight);
                        MeanWeightFeature=reshape(MeanWeightFeature,nDim1,nDim2);
                        save([OutputDir,filesep,'WeightSelectFeature_',num2str(iWeight),'.mat'],'MeanWeightFeature');
                    end
                else
                    WeightImg=zeros(1,nDim1*nDim2);
                    MeanWeight=mean(WeightAll,2);
                    WeightImg(IndexColum)=MeanWeight;
                    WeightImg=reshape(WeightImg,nDim1,nDim2);
                    save([OutputDir,filesep,'WeightFeature.mat'],'WeightImg')
                end
            end
            if sel_CheckWeightI
                
                if length(size(WeightAll))>2
                    for i=1:size(WeightAll,2)
                        for j=1:size(WeightAll,3)
                            WeightImg=zeros(1,nDim1*nDim2);
                            WeightImg(IndexColum(IX(1:SelectFeatureNum(i),j)))=WeightAll(1:SelectFeatureNum(i),i,j);
                            WeightImg=reshape(WeightImg,nDim1,nDim2);
                            save([OutputDir,filesep,'WeightSelectFeature_',num2str(i),'_Fold_',num2str(j),'.mat'],'WeightImg')
                        end
                    end
                else
                    for i=1:size(WeightAll,2)
                        WeightImg=zeros(1,nDim1*nDim2);
                        WeightImg(IndexColum)=WeightAll(:,i);
                        WeightImg=reshape(WeightImg,nDim1,nDim2);
                        save([OutputDir,filesep,'WeightFeatureFold_',num2str(i),'.mat'],'WeightImg')
                    end
                end
            end
        end
    end
end
