function [Data,IndexColum]=read_data(Cfg)
%This function can be used to read data (such as .nii, .img, .mat, .txt, 
% .gii).
%After NAN removed or normalization according to row for the input feature,
%the Data and its index will be returned.
GroupFile=Cfg.GroupFile;
OutputDir=Cfg.OutputDir;

strSelectFeature=Cfg.FeatureSelectionStr;
FileNum=numel(Cfg.GroupFile);
TempData=[];
DataTxt=[];
Data=[];
MaskPath=Cfg.MaskPath;
[Path,Name,MaskType]=fileparts( MaskPath);
Mask=[];
if ~isempty(MaskPath)
    if strcmp(MaskType,'.nii') || strcmp(MaskType,'.img')
        VMask=spm_vol(MaskPath);
        Mask=spm_read_vols(VMask);
    end
    if strcmp(MaskType,'.txt')
        Mask=load(MaskPath);
    end
    if strcmp(MaskType,'.mat')
        Mask=cell2mat(struct2cell(load(MaskPath)));
    end
end
%%
%read data
disp(sprintf('>>>>>>>>>>>>Start reading data>>>>>>>>>>>>>'))
IndexColum=[];
[Path,Name,fileType]=fileparts(GroupFile{1});

switch fileType
    case '.txt'
        for iD=1:FileNum
            DataTxt=load(GroupFile{iD});
            if ~isempty(MaskPath)
                %If there are mask, the data will be transformed based on the mask.
                nDim1=size(Mask,1);
                nDim2=size(Mask,2);
                Mask=reshape(Mask,1,nDim1*nDim2);
                DataTxt=reshape(DataTxt,1,nDim1*nDim2);
                TempData(iD,:)=DataTxt.*Mask;
            else
                %the read data is a row vector, and there are no mask
                TempData(iD,:)=DataTxt;
            end
        end
        Data=[Data;TempData];
    case '.mat'
        for iD=1:FileNum
            DataTxt=cell2mat(struct2cell(load(GroupFile{iD})));
            if ~isempty(MaskPath)
                %If there are mask, the data will be transformed based on the mask.
                nDim1=size(Mask,1);
                nDim2=size(Mask,2);
                Mask=reshape(Mask,1,nDim1*nDim2);
                DataTxt=reshape(DataTxt,1,nDim1*nDim2);
                TempData(iD,:)=DataTxt.*Mask;
            else
                TempData(iD,:)=DataTxt;
            end
        end
        Data=[Data;TempData];
    case '.gii'
        for iD=1:FileNum
            X=gifti(GroupFile{iD});
            DataTxt=X.cdata;
            if ~isempty(MaskPath)
                %If there are mask, the data will be transformed based on the mask.
                nDim1=size(Mask,1);
                nDim2=size(Mask,2);
                DataTxt=reshape(DataTxt,1,nDim1*nDim2);
                TempData(iD,:)=DataTxt.*Mask;
            else
                TempData(iD,:)=DataTxt;
            end
        end
        Data=[Data;TempData];
    case {'.img','.nii'}
        nDim1=size(Mask,1);
        nDim2=size(Mask,2);
        nDim3=size(Mask,3);
        Mask=reshape(Mask,1,nDim1*nDim2*nDim3);
        for i=1:FileNum
            V=spm_vol(GroupFile{i});
            image=spm_read_vols(V);
            image=reshape(image,1,nDim1*nDim2*nDim3);
            image=image.*Mask;
            img(i,:)=image;
        end
        Data=[Data;img];
end
for i=1:size(Data,1)
    tempData=Data(i,:);
    tempData(isnan(tempData)) = 0;%NuN will set as 0
    tempData(find(tempData==Inf))=0;
    Data(i,:)=tempData;
end

%delete the colum with 0
ZeroColum=all(Data==0,1); % such as 0 1 0 0
NotZeroColum=~ZeroColum;
IndexColum=find(NotZeroColum);
Data(:,find(all(Data==0,1)))=[];

FishZStr=Cfg.FishZStr;
%feature change
switch FishZStr
    case 'Normalization for Each Sample'
        Data=zscore(Data,0,2);
end
save([OutputDir,filesep,'RawData.mat'], 'Data');%the result of adding Mask
disp(sprintf('>>>>>>>>>>>>Reading data finished>>>>>>>>>>>>>'))
end
