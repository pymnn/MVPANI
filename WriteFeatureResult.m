function WriteFeatureResult(FeatureAll,IX,Cfg)
%Output which features are selected
%FeatureAll means the selected feature in each fold
%FeatureMean means the pecentage of the selected feature in all folds
%FeatureImgMean means the pecentage of the selected feature in all folds in the
%mask(the mask may be image or matrix)
MaskPath=Cfg.MaskPath;
sel_CheckWeightAll=Cfg.sel_CheckWeightAll;
sel_CheckWeightI=Cfg.sel_CheckWeightI;
strSelectFeature=Cfg.FeatureSelectionStr;
OutputDir=Cfg.OutputDir;
nDim1=Cfg.Mask.nDim1;
nDim2=Cfg.Mask.nDim2;
nDim3=Cfg.Mask.nDim3;
IndexColum=Cfg.IndexColum;
SelectFeatureNum=Cfg.SelectFeatureNum;
FeatureNum=Cfg.DataNum;
strFeatureReduce=Cfg.FeatureReduceStr;
if ~strcmp(strSelectFeature,'None')
% only for feature selection    
    if isempty(MaskPath)
%         save([OutputDir,filesep,'FeatureAll.mat'],'FeatureAll');
        if sel_CheckWeightAll
            %the pecentage of the selected feature
            %had done feature selection
            FeatureImg=zeros(FeatureNum,size(FeatureAll,2),size(FeatureAll,3));
            for i=1:size(FeatureAll,2)
                for j=1:size(FeatureAll,3)
                    FeatureImg(IndexColum(IX(1:SelectFeatureNum(i),j)),i,j)=FeatureAll(1:SelectFeatureNum(i),i,j);
                end
            end
            MeanFeature=mean(FeatureImg,3);
            for iFeature=1:size(MeanFeature,2)
                SelectFeatureMean=MeanFeature(:,iFeature);
                save([OutputDir,filesep,'SelectFeaturePercentage_',num2str(iFeature),'.mat'],'SelectFeatureMean')
            end
        end
        if sel_CheckWeightI 
            %output which features are selected in each fold 
            %i denote select feature,j denote fold
            SelectFeature=zeros(FeatureNum,1);
            for i=1:size(FeatureAll,2)
                for j=1:size(FeatureAll,3)
                    SelectFeature(IndexColum(IX(1:SelectFeatureNum(i),j)))=FeatureAll(1:SelectFeatureNum(i),i,j);
                    save([OutputDir,filesep,'SelectFeature_',num2str(i),'_Fold_',num2str(j),'.mat'],'SelectFeature')
                end
            end
        end
    else
        [Path,Name,fileType]=fileparts(MaskPath);
        if strcmp(fileType,'.img')|| strcmp(fileType,'.nii')
            %if the mask is image, the Feature will be mapped in the brain
            if sel_CheckWeightAll
                %output the percentage of the selected features
                %had done feature selection;size(FeatureAll,2):feature
                %number;size(FeatureAll,3)): fold number
                SelectedFeatureImg=zeros(nDim1*nDim2*nDim3,size(FeatureAll,2),size(FeatureAll,3));
                for i=1:size(FeatureAll,2)
                    %i denote select feature,j denote fold
                    for j=1:size(FeatureAll,3)
                        SelectedFeatureImg(IndexColum(IX(1:SelectFeatureNum(i),j)),i,j)=FeatureAll(1:SelectFeatureNum(i),i,j);
                    end
                end
                MeanFeature=mean(SelectedFeatureImg,3);
                for iFeature=1:size(MeanFeature,2)
                    SelectFeature=MeanFeature(:,iFeature);
                    SelectFeature=reshape(SelectFeature,nDim1,nDim2,nDim3);
                    VMask=spm_vol(MaskPath);
                    v=spm_create_vol(VMask);
                    v.dt=[16 0];
                    v.fname=[OutputDir,filesep,'SelectFeaturePercentageImg_',num2str(iFeature),'.nii'];
                    spm_write_vol(v,SelectFeature);
                    save([OutputDir,filesep,'SelectFeaturePercentageImg_',num2str(iFeature),'.mat'],'SelectFeature');
                end
            end
            if sel_CheckWeightI
                %output which features are selected in each fold 
%                 save([OutputDir,filesep,'FeatureAll.mat'],'FeatureAll');
                for i=1:size(FeatureAll,2)
                    %i denote select feature,j denote fold
                    for j=1:size(FeatureAll,3)
                        SelectedFeatureImg=zeros(1,nDim1*nDim2*nDim3);
                        SelectedFeatureImg(IndexColum(IX(1:SelectFeatureNum(i),j)))=FeatureAll(1:SelectFeatureNum(i),i,j);
                        SelectedFeatureImg=reshape(SelectedFeatureImg,nDim1,nDim2,nDim3);
                        VMask=spm_vol(MaskPath);
                        v=spm_create_vol(VMask);
                        v.dt=[16 0];
                        v.fname=[OutputDir,filesep,'SelectFeatureImg_',num2str(i),'_Fold_',num2str(j),'.nii'];
                        save([OutputDir,filesep,'SelectFeatureImg_',num2str(i),'_Fold_',num2str(j),'.mat'],'SelectedFeatureImg');
                        spm_write_vol(v,SelectedFeatureImg);
                        clear SelectedFeatureImg VMask v
                    end
                end
            end
        end
        if strcmp(fileType,'.txt')|| strcmp(fileType,'.mat')
            %if the mask is txt, the weight will be mapped in the txt
            if sel_CheckWeightAll
                %had done feature selection
                FeatureImg=zeros(nDim1*nDim2,size(FeatureAll,2),size(FeatureAll,3));
                for i=1:size(FeatureAll,2)
                    %i denote select feature,j denote fold
                    for j=1:size(FeatureAll,3)
                        FeatureImg(IndexColum(IX(1:SelectFeatureNum(i),j)),i,j)=FeatureAll(1:SelectFeatureNum(i),i,j);
                    end
                end
                MeanFeature=mean(FeatureImg,3);
                for iFeature=1:size(MeanFeature,2)
                    SelectFeature=MeanFeature(:,iFeature);
                    SelectFeature=reshape(SelectFeature,nDim1,nDim2);
                    save([OutputDir,filesep,'SelectFeaturePercentage_',num2str(iFeature),'.mat'],'SelectFeature');
                end
            end
            if sel_CheckWeightI
%                 save([OutputDir,filesep,'FeaturetAll.mat'],'FeatureAll');
                for i=1:size(FeatureAll,2)
                    for j=1:size(FeatureAll,3)
                        FeatureImg=zeros(1,nDim1*nDim2);
                        FeatureImg(IndexColum(IX(1:SelectFeatureNum(i),j)))=FeatureAll(1:SelectFeatureNum(i),i,j);
                        FeatureImg=reshape(FeatureImg,nDim1,nDim2);
                        save([OutputDir,filesep,'SelectFeature_',num2str(i),'_Fold_',num2str(j),'.mat'],'FeatureImg')
                    end
                end
            end
        end
    end
end

