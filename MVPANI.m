function varargout = MVPANI(varargin)
% MVPANI MATLAB code for MVPANI.fig
%      MVPANI, by itself, creates a new MVPANI or raises the existing
%      singleton*.
%
%      H = MVPANI returns the handle to a new MVPANI or the handle to
%      the existing singleton*.
%
%      MVPANI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MVPANI.M with the given input arguments.
%
%      MVPANI('Property','Value',...) creates a new MVPANI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MVPANI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MVPANI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MVPANI

% Last Modified by GUIDE v2.5 04-Aug-2020 22:09:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MVPANI_OpeningFcn, ...
    'gui_OutputFcn',  @MVPANI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before MVPANI is made visible.
function MVPANI_OpeningFcn(hObject, eventdata, handles, varargin)
Release='V1.0_20190604';
if ispc
    UserName =getenv('USERNAME');
else
    UserName =getenv('USER');
end
Datetime=fix(clock);
fprintf('Welcome: %s, %.4d-%.2d-%.2d %.2d:%.2d \n', UserName,Datetime(1),Datetime(2),Datetime(3),Datetime(4),Datetime(5));
fprintf('MVPANI Release=%s\n',Release);
fprintf('Copyright(c) 2019~2025 \nTianjin Key Lab of Functional Imaging, Tianjin Medical University, China\n');
fprintf('1 Guangdong Road, Hexi District, Tianjin, China 300203\n');
fprintf('http://funi.tmu.edu.cn\n');
fprintf('-----------------------------------------------------------\n');
fprintf('Citing Information:\nIf you think MVPANI is useful for your work, please cite the following paper in your work:\n')
fprintf('Reference: Peng Y, Zhang X, Li Y, Su Q, Wang S, Liu F, Yu C, Liang M. ')
fprintf('MVPANI: A Toolkit With Friendly Graphical User Interface for Multivariate Pattern Analysis of Neuroimaging Data. Front.Neurosci. 2020,14:545\n');


axes(handles.axes_logo);
% axis image;
MVPANIPath=fileparts(which('MVPANI.m'));
[A, map, alpha] = imread(fullfile(MVPANIPath,'LOGO.png'));
h = imshow(A, map);
set(h, 'AlphaData', alpha);

handles.Cfg.MethodSelectionStr=' Classification';
handles.Cfg.ClassAlgorithmStr='SVM (Support Vector Machine)';
handles.Cfg.modelParemter='-s 0 -t 0';
handles.Cfg.FeatureSelectionStr='None';
handles.Cfg.FeatureReduceStr='None';
handles.Cfg.FishZStr='None';
handles.Cfg.IsSearchLight=0;
handles.Cfg.IscheckboxROC=0;
handles.Cfg.IsWeightI=0;
handles.Cfg.IsWeightAll=0;
handles.Cfg.IscheckboxSenSpeFold=0;
handles.Cfg.ValueClassAlgorithm=1;
handles.Cfg.ValueFishZ=1;
handles.Cfg.selValueFeature=1;
handles.Cfg.stlFeatureReduce=1;
handles.Cfg.selValueMethod=1;
handles.Cfg.MaskEntry=[];
handles.Cfg.InputFoldEntry=[];
handles.Cfg.strSelectFold='radFold';
handles.Cfg.SelectFeatureNum=[];
handles.Cfg.MvpaResults={};
handles.Cfg.DataNum=[];
handles.Cfg.ContrPCA=95;
handles.Cfg.MethodSelectionStr=' Classification';
handles.Cfg.OutputDir=pwd;
handles.Cfg.PermutationNum=1000;
handles.Cfg.edFold=10;
handles.Cfg.GroupFeatureFusionList={};
handles.Cfg.FeatureNum='';
handles.Cfg.SphereRadius=6;
handles.Cfg.IscheckboxSenSpeFold=0;
handles.Cfg.MaskPath={};
handles.Cfg.LabelPath={};
handles.Cfg.GroupListboxStr={};

handles.output = hObject;
handles.GroupCells={};
handles.GroupLabel=[];
handles.GroupNum=0;
handles.DataNum=[];
handles.Fold=[];
handles.CreateLabel=[];
handles.ImageCells={};
handles.ImageLabel={};
handles.LabelFold={};
handles.FileData2=[];
handles.FileData=[];
handles.edit12=[];
handles.Data=[];
handles.Mask={};
handles.GroupFile={};
handles.SelectFeatureNum=[];
handles.MeanAccuracy=[]; %true accuracy
handles.MvpaResults={};
handles.CurDir=pwd;
handles.modelParemter='-s 0 -t 0';
handles.SearchLightRealResults=[];
handles.FeatureSelectionStr='None';
handles.FeatureReduceStr='None';
handles.FishZStr='None';
handles.ContrPCA=95;
handles.MethodSelectionStr=' Classification';
handles.OutputDir=pwd;
handles.IsSearchLight=0;
handles.IscheckboxROC=0;
handles.IsWeightI=0;
handles.IsWeightAll=0;
handles.IscheckboxSenSpeFold=0;
%GUI1 and GUI2 %20190111
% handles.output=hObject;
guidata(hObject, handles);%old
%
% UIWAIT makes MVPANI wait for user response (see UIRESUME)
% uiwait(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = MVPANI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in GroupListbox.


% --- Executes on button press in GroupAdd.
function GroupAdd_Callback(hObject, eventdata, handles)
GroupFile=handles.GroupFile;
FileNum=numel(handles.GroupFile);

[FileName,Path] = uigetfile({'*.nii',' *.nii';...
    '*.img',' *.img';...
    '*.mat',' *.mat';...
    '*.txt',' *.txt';...
    },'Select data files', 'MultiSelect', 'on');
if isnumeric(Path)
    return
end
[Path,Name,fileType]=fileparts( fullfile(Path,FileName{1}));

for i=FileNum+1:length(FileName)+FileNum
    GroupFile{i}=[Path,filesep,FileName{i-FileNum}];
end
handles.GroupFile=GroupFile;
handles.Cfg.GroupFile=GroupFile;

for i=1:size(FileName,2)
    StringOne={sprintf('(%s) %s',  FileName{i},['ID:',sprintf('%03d',i),filesep,'Tol:',num2str(size(FileName,2))],['Path:', Path,filesep,FileName{i}])};
    AddString(handles.GroupListbox, StringOne);
    guidata(hObject, handles);
end
handles.Cfg.GroupListboxStr=get(handles.GroupListbox,'String');
guidata(hObject, handles);


% --- Executes on button press in GroupRemove.
function GroupRemove_Callback(hObject, eventdata, handles)
% hObject    handle to GroupRemove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Value1=get(handles.GroupListbox, 'Value');
if Value1==0
    return
end
handles.Cfg.GroupFile(Value1)=[];
RemoveString(handles.GroupListbox, Value1);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ExpressionEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExpressionEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes during object creation, after setting all properties.
function OutputEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to OutputEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OutputButton.
function OutputButton_Callback(hObject, eventdata, handles)
% hObject    handle to OutputButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FilePath=uigetdir(handles.CurDir, 'Pick Output Directory');
if isnumeric(FilePath)
    return
end
handles.CurDir=fileparts(FilePath);
set(handles.OutputEntry, 'String', FilePath);
handles.OutputDir=FilePath;
handles.Cfg.OutputDir=handles.OutputDir;

guidata(hObject, handles);



% --- Executes on button press in ComputeButton.
function ComputeButton_Callback(hObject, eventdata, handles)
% hObject    handle to ComputeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Label=[];
sel_CheckWeightI=get(handles.WeightI, 'value');%output weight of each fold
sel_CheckWeightAll=get(handles.WeightAll, 'value'); %output mean weight

handles.Cfg.sel_CheckWeightAll=sel_CheckWeightAll;
handles.sel_CheckWeightAll=sel_CheckWeightAll;
handles.Cfg.sel_CheckWeightI=sel_CheckWeightI;
handles.sel_CheckWeightI=sel_CheckWeightI;

sel_checkboxROC=get(handles.checkboxROC, 'value');
sel_CheckboxSearchLight=get(handles.checkboxSearchLight, 'value');
handles.Cfg.sel_checkboxROC=sel_checkboxROC;
handles.sel_checkboxROC=sel_checkboxROC;
handles.Cfg.sel_CheckboxSearchLight=sel_CheckboxSearchLight;
handles.sel_CheckboxSearchLight=sel_CheckboxSearchLight;
sel_checkboxSenSpeFold=handles.Cfg.IscheckboxSenSpeFold;%output Specificity,Sensitivity and accuracy of each fold
%%
%read label
LabelPath=get(handles.LabelEntry, 'String');
[LabPath,LabName,LabelFileType]=fileparts(LabelPath);
if strcmp(LabelFileType,'.txt')
    Label=load(LabelPath);
end
if strcmp(LabelFileType,'.mat')
    Label=cell2mat(struct2cell(load(LabelPath)));
end
if strcmp(LabelFileType,'.xlsx')
    Label=xlsread(LabelPath);
end
if strcmp(LabelFileType,'.xls')
    Label=xlsread(LabelPath);
end
handles.Label=Label;
handles.LabelPath=LabelPath;
handles.Cfg.LabelPath=handles.LabelPath;
handles.Cfg.Label=handles.Label;
%%
%read group
GroupListboxStr=get(handles.GroupListbox,'String');
for i=1:size(GroupListboxStr,1)
    PathPos=strfind(GroupListboxStr{i,:},'Path:');
    GroupListbox=GroupListboxStr{i,:};
    GroupFile{i}=GroupListbox(PathPos+5:length(GroupListbox)-2);
end
handles.GroupFile=GroupFile;
handles.Cfg.GroupFile=GroupFile;
handles.GroupListboxStr=GroupListboxStr;
handles.Cfg.GroupListboxStr=GroupListboxStr;

%%
OutputDir=get(handles.OutputEntry, 'String');
handles.Cfg.OutputDir=OutputDir;
%%
%read mask
MaskPath=get(handles.MaskEntry, 'String');
handles.MaskPath=MaskPath;
handles.Cfg.MaskPath=handles.MaskPath;
if ~isempty(MaskPath)
    [FilePath,Name,MaskType]=fileparts(MaskPath);
    if strcmp(MaskType,'.nii') || strcmp(MaskType,'.img')
        VMask=spm_vol(MaskPath);
        Mask=spm_read_vols(VMask);
    end
    if strcmp(MaskType,'.txt')
        Mask=load(MaskPath);
    end
    if strcmp(MaskType,'.txt') || strcmp(MaskType,'.mat')
        Mask=importdata(MaskPath);
    end
    handles.Mask.img=Mask;
    handles.Cfg.Mask.img=Mask;
    handles.Cfg.Mask.nDim1=size(Mask,1);
    handles.Cfg.Mask.nDim2=size(Mask,2);
    handles.Cfg.Mask.nDim3=size(Mask,3);
    handles.Mask.nDim1=size(Mask,1);
    handles.Mask.nDim2=size(Mask,2);
    handles.Mask.nDim3=size(Mask,3);
else
    handles.Mask.img=[];
    handles.MaskPath='';
    handles.Cfg.MaskPath='';
    handles.Cfg.Mask.img=[];
    handles.Cfg.Mask.nDim1=0;
    handles.Cfg.Mask.nDim2=0;
    handles.Cfg.Mask.nDim3=0;
    handles.Mask.nDim1=0;
    handles.Mask.nDim2=0;
    handles.Mask.nDim3=0;
end
%%
[Data,IndexColum]=read_data(handles.Cfg);
handles.Cfg.IndexColum=IndexColum;
handles.Cfg.Data=Data;
handles.Cfg.DataNum=size(Data,2);
handles.DataNum=size(Data,2);
DataNum=size(Data,2);
MvpaResults.FeatureNum=DataNum;
MvpaResults.IndexColum=IndexColum;
%%
%read Fold
SelectFold=get(handles.BttSelectFold, 'SelectedObject');
strSelectFold=get(SelectFold,'tag');

handles.Cfg.strSelectFold =strSelectFold;
handles.strSelectFold =strSelectFold;
if strcmp(strSelectFold,'radFold')
    KFold=str2double(get(handles.edFold, 'String'));
    handles.Cfg.edFold=KFold;
    Fold=crossvalind('kfold',size(Data,1),KFold);
    handles.Cfg.Fold=Fold;
    handles.Fold=Fold;
    MvpaResults.Fold.strSelectFold='Random fold';
    MvpaResults.Fold.KFold=KFold;
else
    FoldPath=get(handles.InputFoldEntry,'String');
    [FPath,FoldName,FoldFileType]=fileparts(FoldPath);
    if strcmp(FoldFileType,'.txt')
        Fold=load(FoldPath);
    end
    if strcmp(FoldFileType,'.mat')
        Fold=cell2mat(struct2cell(load(FoldPath)));
    end
    if strcmp(FoldFileType,'.xlsx')
        Fold=xlsread(FoldPath);
    end
    if strcmp(FoldFileType,'.xls')
        Fold=xlsread(FoldPath);
    end
    handles.Fold=Fold;
    handles.Cfg.Fold=handles.Fold;
    handles.FoldPath=FoldPath;
    handles.Cfg.FoldPath=handles.FoldPath;
    MvpaResults.Fold.FoldPath=FoldPath;
    MvpaResults.Fold.FoldValue=Fold;
end
MvpaResults.MaskPath=MaskPath;
MvpaResults.Label.LabelPath=LabelPath;
MvpaResults.Label.LabelValue=Label;
%%
%feature change

ValueFishZ=get(handles.FishZ, 'value');
str=get(handles.FishZ, 'String');
strFishZ=str{ValueFishZ};

selFeatureReduce=get(handles.FeatureReduce, 'value');
str=get(handles.FeatureReduce, 'String');
FeatureReduceStr=str{selFeatureReduce};

%%
%feature selection
%none no feature selection
%weight：select features with max weight；%F-Score;
%there are two methods to select feature.(1)based on number; (2)based on %
selValueFeature=get(handles.FeatureSelection, 'value');
str=get(handles.FeatureSelection, 'String');
strSelectFeature=str{selValueFeature};
handles.Cfg.FeatureSelectionStr=strSelectFeature;

selValueMethod=get(handles.popmMethod, 'value');
str=get(handles.popmMethod, 'String');
strSelectMethod=str{selValueMethod};

if ~strcmp(strSelectFeature,'None') && ~strcmp(strSelectFeature,'Lasso')
    SelectFeatureNum=get(handles.editFeatureNum,'String');
    handles.FeatureNum=SelectFeatureNum;
    handles.Cfg.FeatureNum=SelectFeatureNum;
    SelectFeatureNum=eval(['[',SelectFeatureNum,']']);
    
    % input 50:50:1000,means we first select the max features, then step is 50, until select the max 1000 features
    NumType=SelectFeatureNum(1);
    
    if NumType>1% based on number to select feature
    else% based on % to select feature,input the values such as 0.2：0.1：0.8, the max value is 1
        SelectFeatureNum=floor(SelectFeatureNum*DataNum);
    end
    handles.SelectFeatureNum=SelectFeatureNum;
    handles.Cfg.SelectFeatureNum=handles.SelectFeatureNum;
    MvpaResults.SelectFeatureNum=SelectFeatureNum;
    
end
MvpaResults.Method=strSelectMethod;
MvpaResults.SelectFeature=strSelectFeature;
MvpaResults.FeatureTansformation=strFishZ;

%%
tic

if ~sel_CheckboxSearchLight % no select SearchLight
    %%
    %the methods of selectiong fold, there are two methods: input fold
    %file ; or input the numbler of k fold
    if  strcmp(strSelectMethod,' Regression')
        
        strClassAlgorithm=handles.Cfg.ClassAlgorithmStr;%select different  classification algorithms
        switch strClassAlgorithm
            case 'SVM (Support Vector Machine)'
                disp(sprintf('>>>>>>>>>>>>Start SVM regression>>>>>>>>>>>>>'));
                [WeightAll,FeatureAll,IX,r,p,predLabel,MaxDim]=SVMRegression(handles.Cfg);
                WriteWeightResult(WeightAll,IX,handles.Cfg);
                WriteFeatureResult(WeightAll,IX,handles.Cfg);
                WriteCoeffResult(r,handles.Cfg);
                MvpaResults.ModelParemter=handles.Cfg.modelParemter;
                MvpaResults.r=r;
                MvpaResults.p=p;
                MvpaResults.predLabel=predLabel;
                MvpaResults.ClassAlgorithm=strClassAlgorithm;
                if strcmp(FeatureReduceStr,'PCA (Principal Component Analysis)')
                    nDimContrPCA= handles.Cfg.ContrPCA;
                    DimensionalityReduction=sprintf('%s, nDimContrPCA=%d %%, NumPCA=%d',FeatureReduceStr,nDimContrPCA,MaxDim);
                    MvpaResults.DimensionalityReduction=DimensionalityReduction;
                else
                    MvpaResults.DimensionalityReduction=FeatureReduceStr;
                end
                disp( sprintf('>>>>>>>>>>>>SVM regression finished>>>>>>>>>>>>>'));
            case 'RF (Random Forest)'
                disp( sprintf('>>>>>>>>>>>>Start RandomForest regression>>>>>>>>>>>>>'));
                [WeightAll,FeatureAll,IX,r,p,predLabel,MaxDim]=RandomForestRegression(handles.Cfg);
                WriteCoeffResult(r,handles.Cfg);
                WriteFeatureResult(WeightAll,IX,handles.Cfg);
                MvpaResults.r=r;
                MvpaResults.p=p;
                MvpaResults.predLabel=predLabel;
                NumTree=handles.Cfg.NumTree;
                ClassAlgorithm=sprintf('%s, NumTree=%d',strClassAlgorithm,NumTree);
                MvpaResults.ClassAlgorithm=ClassAlgorithm;
                if strcmp(FeatureReduceStr,'PCA (Principal Component Analysis)')
                    nDimContrPCA= handles.Cfg.ContrPCA;
                    DimensionalityReduction=sprintf('%s, nDimContrPCA=%d %%, NumPCA=%d',FeatureReduceStr,nDimContrPCA,MaxDim);
                    MvpaResults.DimensionalityReduction=DimensionalityReduction;
                else
                    MvpaResults.DimensionalityReduction=FeatureReduceStr;
                end
                disp(sprintf('>>>>>>>>>>>>RandomForest regression finished>>>>>>>>>>>>>'));
            case 'KNN (K-Nearest Neighbor)'
                disp(sprintf('>>>>>>>>>>>>Start KNN regression>>>>>>>>>>>>>'));
                [WeightAll,FeatureAll,IX,r,p,predLabel,MaxDim]=KNNRegression(handles.Cfg);
                WriteCoeffResult(r,handles.Cfg);
                WriteFeatureResult(WeightAll,IX,handles.Cfg);
                MvpaResults.r=r;
                MvpaResults.p=p;
                MvpaResults.predLabel=predLabel;
                NumNeighbor=handles.Cfg.NumNeighbor;
                ClassAlgorithm=sprintf('%s, NumNeighbor=%d',strClassAlgorithm,NumNeighbor);
                MvpaResults.ClassAlgorithm=ClassAlgorithm;
                if strcmp(FeatureReduceStr,'PCA (Principal Component Analysis)')
                    nDimContrPCA= handles.Cfg.ContrPCA;
                    DimensionalityReduction=sprintf('%s, nDimContrPCA=%d %%, NumPCA=%d',FeatureReduceStr,nDimContrPCA,MaxDim);
                    MvpaResults.DimensionalityReduction=DimensionalityReduction;
                else
                    MvpaResults.DimensionalityReduction=FeatureReduceStr;
                end
                disp(sprintf('>>>>>>>>>>>>KNN regression finished>>>>>>>>>>>>>'));
            case 'DT (Decision Tree)'
                disp(sprintf('>>>>>>>>>>>>Start DecisionTree regression>>>>>>>>>>>>>'));
                [WeightAll,FeatureAll,IX,r,p,predLabel,MaxDim]=DecisionTreeRegression(handles.Cfg);
                WriteCoeffResult(r,handles.Cfg);
                WriteFeatureResult(WeightAll,IX,handles.Cfg);
                MvpaResults.r=r;
                MvpaResults.p=p;
                MvpaResults.predLabel=predLabel;
                MvpaResults.ClassAlgorithm=strClassAlgorithm;
                if strcmp(FeatureReduceStr,'PCA (Principal Component Analysis)')
                    nDimContrPCA= handles.Cfg.ContrPCA;
                    DimensionalityReduction=sprintf('%s, nDimContrPCA=%d %%, NumPCA=%d',FeatureReduceStr,nDimContrPCA,MaxDim);
                    MvpaResults.DimensionalityReduction=DimensionalityReduction;
                else
                    MvpaResults.DimensionalityReduction=FeatureReduceStr;
                end
                disp(sprintf('>>>>>>>>>>>>DecisionTree regression finished>>>>>>>>>>>>>'));
        end
        disp(sprintf('Correlation coefficient: %.4f',r))
    end
    if strcmp(strSelectMethod,' Classification')
        strClassAlgorithm=handles.Cfg.ClassAlgorithmStr;%select different  classification algorithms

        switch strClassAlgorithm
            case 'SVM (Support Vector Machine)'
                disp(sprintf('>>>>>>>>>>>>Start SVM classification>>>>>>>>>>>>>'));
                [Accuracy,Predict_test_Label,SpecificitySensitivity,WeightAll,FeatureAll,IX,Specif_Sensit_acc_fold,MaxDim]=SVMMVPA(handles.Cfg);
                MeanAccuracy=mean(Accuracy);
                
                WriteWeightResult(WeightAll,IX,handles.Cfg)
                WriteFeatureResult(FeatureAll,IX,handles.Cfg)

                MvpaResults.ModelParemter=handles.Cfg.modelParemter;

                disp('>>>>>>>>>>>>SVM classification finished>>>>>>>>>>>>>');
            case 'KNN (K-Nearest Neighbor)'
                disp(sprintf('>>>>>>>>>>>>Start KNN classification>>>>>>>>>>>>>'));
                [Accuracy, Predict_test_Label,SpecificitySensitivity,WeightAll,FeatureAll,IX,Specif_Sensit_acc_fold,MaxDim]=KNNMVPA(handles.Cfg);
                MeanAccuracy=mean(Accuracy);
                WriteFeatureResult(FeatureAll,IX,handles.Cfg)
            
                NumNeighbor=handles.Cfg.NumNeighbor;
                ClassAlgorithm=sprintf('%s, NumNeighbor=%d',strClassAlgorithm,NumNeighbor);
                                                       
                disp(sprintf('>>>>>>>>>>>>KNN classification finished>>>>>>>>>>>>>'));
            case 'LR (Logistic Regression)'
                fprintf('The LR classification is processing,-----\n')
                [Accuracy, Predict_test_Label,SpecificitySensitivity,WeightAll,FeatureAll,IX,Specif_Sensit_acc_fold,MaxDim]=logistic(handles.Cfg);
                MeanAccuracy=mean(Accuracy);
                WriteFeatureResult(FeatureAll,IX,handles.Cfg);
                
                disp('>>>>>>>>>>>>LR classification finished>>>>>>>>>>>>>');
                guidata(hObject, handles);
            case 'LDA (Linear Discriminant Analysis)'
                fprintf('The LDA classification is processing,-----')
                [Accuracy, Predict_test_Label,SpecificitySensitivity,WeightAll,FeatureAll,IX,Specif_Sensit_acc_fold,MaxDim]=LDAMVPA(handles.Cfg);
                MeanAccuracy=mean(Accuracy);
                WriteFeatureResult(WeightAll,IX,handles.Cfg);
                                
                disp('>>>>>>>>>>>>LDA classification finished>>>>>>>>>>>>>');
            case 'NB (Naive Bayes)'
                disp(sprintf('>>>>>>>>>>>>Start NaiveBayes classification>>>>>>>>>>>>>'));
                [Accuracy, Predict_test_Label,SpecificitySensitivity,WeightAll,FeatureAll,IX,Specif_Sensit_acc_fold,MaxDim]=NaiveBayesMVPA(handles.Cfg);
                MeanAccuracy=mean(Accuracy);
                WriteFeatureResult(WeightAll,IX,handles.Cfg);
                                            
                disp(sprintf('>>>>>>>>>>>>NaiveBayes classification finished>>>>>>>>>>>>>'));
            case 'RF (Random Forest)'
                disp(sprintf('>>>>>>>>>>>>Start RandomForest classification>>>>>>>>>>>>>'));
                [Accuracy, Predict_test_Label,SpecificitySensitivity,WeightAll,FeatureAll,IX,Specif_Sensit_acc_fold,MaxDim]=RandomForestMVPA(handles.Cfg);
                MeanAccuracy=mean(Accuracy);
                WriteFeatureResult(WeightAll,IX,handles.Cfg);
                
                NumTree=handles.Cfg.NumTree;
                ClassAlgorithm=sprintf('%s, NumTree=%d',strClassAlgorithm,NumTree);
                               
                disp(sprintf('>>>>>>>>>>>>RandomForest classification finished>>>>>>>>>>>>>'));
            case 'DT (Decision Tree)'
                disp(sprintf('>>>>>>>>>>>>Start DecisionTree classification>>>>>>>>>>>>>'));
                [Accuracy, Predict_test_Label,SpecificitySensitivity,WeightAll,FeatureAll,IX,Specif_Sensit_acc_fold,MaxDim]=DecisionTreeMVPA(handles.Cfg);
                MeanAccuracy=mean(Accuracy);
                WriteFeatureResult(WeightAll,IX,handles.Cfg);
                
                disp(sprintf('>>>>>>>>>>>>DecisionTree classification finished>>>>>>>>>>>>>'));
        end
        %%%%%%%%%%
        
        MvpaResults.MeanAccuracy=MeanAccuracy;
        MvpaResults.Predict_test_Label=Predict_test_Label;
        MvpaResults.ClassAlgorithm=strClassAlgorithm;
        if strcmp(FeatureReduceStr,'PCA (Principal Component Analysis)')
            nDimContrPCA= handles.Cfg.ContrPCA;
            DimensionalityReduction=sprintf('%s, nDimContrPCA=%d %%, NumPCA=%d',FeatureReduceStr,nDimContrPCA,MaxDim);
            MvpaResults.DimensionalityReduction=DimensionalityReduction;
        else
            MvpaResults.DimensionalityReduction=FeatureReduceStr;
        end
        %%%%%%%%%%%%%%%
        Label=handles.Cfg.Label;
        ReLabel=unique(Label);%true label 0 1; 1 -1
        if length(ReLabel)<3
            
            AUC=WriteROCACCResult(MeanAccuracy,Predict_test_Label,handles.Cfg,strClassAlgorithm);
            MvpaResults.AUC=AUC;
            if sel_checkboxSenSpeFold
                MvpaResults.SpecificitySensitivityAccEachFold=Specif_Sensit_acc_fold;
            end
            SpecificitySensitivity=mean(Specif_Sensit_acc_fold,1);
            MvpaResults.SpecificitySensitivity=SpecificitySensitivity;
        else
            
            if sel_checkboxSenSpeFold
                MvpaResults.AccuracyEachClassEachFold=Specif_Sensit_acc_fold;
            end
            SpecificitySensitivity=mean(Specif_Sensit_acc_fold,1);
            MvpaResults.MeanAccuracyEachClass=SpecificitySensitivity;
            %%
            %confusion_matrix
            for i=1:size(Predict_test_Label,3)
                confusion_matrix=[];
                for X=1:length(ReLabel)
                    for Y=1:length(ReLabel)
                        confusion_matrix(X,Y)=length(find(Predict_test_Label(find(Predict_test_Label(:,1,1)==ReLabel(X)),2,i)==ReLabel(Y)));
                    end
                end
                ConfusionMatrix(:,:,i)=confusion_matrix;
  
                textStrings = num2str(confusion_matrix(:),'%d');  %# Create strings from the matrix values
                textStrings = strtrim(cellstr(textStrings));
                
                h=figure;
                imagesc(confusion_matrix)
                xlabel('Predicted Label','FontSize',15)
                ylabel('True Label','FontSize',15)
                set(gca,'FontSize',15)
                set(gca, 'LineWidth',2)
                
                [x,y] = meshgrid(1:length(ReLabel));   %# Create x and y coordinates for the strings
                hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                    'HorizontalAlignment','center');
                midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
                textColors = repmat(confusion_matrix(:) > midValue,1,3);  %# Choose white or black for the
                
                set(gca,'XTick',1:length(ReLabel),...                         %# Change the axes tick marks
                    'XTickLabel',num2cell(ReLabel),...  %#   and tick labels
                    'YTick',1:length(ReLabel),...
                    'YTickLabel',num2cell(ReLabel),...
                    'TickLength',[0 0]);
                colorbar
                saveas(h,[OutputDir,filesep,'ConfusionMatrix_',num2str(i),'.bmp'])
                saveas(h,[OutputDir,filesep,'ConfusionMatrix_',num2str(i),'.fig'])
                
            end
            MvpaResults.confusion_matrix=ConfusionMatrix;        
        end
                
        handles.Cfg.MvpaResults=MvpaResults;
        handles.MeanAccuracy=MeanAccuracy;
        handles.Cfg.MeanAccuracy=MeanAccuracy;
        handles.MvpaResults=MvpaResults;
        disp(sprintf('Classification Accuracy: %.2f %%',MeanAccuracy))
%         disp(sprintf('Max Classification Accuracy: %.2f %%',max(MeanAccuracy)))
    end
    MvpaResults.GroupFile=handles.Cfg.GroupFile;
    save([OutputDir,filesep,'MvpaResults.mat'],'MvpaResults');
end

if sel_CheckboxSearchLight % select SearchLight
    %%
    %the methods of selectiong fold, there are two methods: input fold
    %file ; or input the numbler of k fold
    SelectFold=get(handles.BttSelectFold, 'SelectedObject');
    strSelectFold=get(SelectFold,'tag');
    if strcmp(strSelectFold,'radFold')
        KFold=str2double(get(handles.edFold, 'String'));
        Fold=crossvalind('kfold',size(Data,1),KFold);
    end
    handles.Cfg.SphereRadius=get(handles.editSphere,'String');
    handles.SphereRadius=get(handles.editSphere,'String');
    SphereRadius=str2num(handles.Cfg.SphereRadius);
    
    %
    
    cfg = decoding_defaults;
    cfg.results.overwrite = 1;
    for i=1:length(GroupFile)
        cfg.files.name{i}=(GroupFile{i});
    end
    
    cfg.analysis = 'searchlight';
    if ~exist([OutputDir,filesep,'results',filesep,'searchlight',filesep])
        mkdir([OutputDir,filesep,'results',filesep,'searchlight',filesep]);
    end
    cfg.results.dir = [OutputDir,filesep,'results',filesep,'searchlight',filesep];
    cfg.files.mask = MaskPath;
    cfg.files.chunk = Fold;
    cfg.files.label = Label;
    cfg.software=spm('ver');
    cfg.design = make_design_cv(cfg);
    cfg.design.unbalanced_data = 'ok';
    cfg.searchlight.unit = 'voxels';
    cfg.searchlight.radius = SphereRadius; % this will yield a searchlight radius of 12mm.
    cfg.feature_selection.check_datatrans_ok = true;
    if  strcmp(strSelectMethod,' Classification')
        cfg.decoding.train.classification.model_parameters = [handles.Cfg.modelParemter '  -q'];%'-s 0 -t 0 -c 1 -b 0 -q';
        cfg.decoding.method = 'classification'; % because feature selection doesn't work with kernel method
        save([OutputDir,filesep,'results',filesep,'searchlight',filesep,'cfg.mat'],'cfg');
        [results,cfg] = decoding(cfg);
        load([OutputDir,filesep,'results',filesep,'searchlight',filesep,'res_accuracy_minus_chance.mat']);
        handles.SearchLightRealResults=results;
    end
    if  strcmp(strSelectMethod,' Regression')
        cfg.decoding.train.regression.model_parameters = [handles.Cfg.modelParemter '  -q'];%'-s 0 -t 0 -c 1 -b 0 -q';
        cfg.decoding.method = 'regression';
        cfg.results.output = {'corr'};
        save([OutputDir,filesep,'results',filesep,'searchlight',filesep,'cfg.mat'],'cfg');
        [results,cfg] = decoding(cfg);
        load([OutputDir,filesep,'results',filesep,'searchlight',filesep,'res_corr.mat']);
        handles.SearchLightRealResults=results;
    end
    
end
toc
guidata(hObject, handles);

function AddString(ListboxHandle, NewCell)
StringCell=get(ListboxHandle, 'String');
StringCell=[StringCell; NewCell];
set(ListboxHandle, 'String', StringCell, 'Value', numel(StringCell));

function RemoveString(ListboxHandle, Value)
StringCell=get(ListboxHandle, 'String');
StringCell(Value)=[];
if isempty(StringCell)
    Value=0;
end
if Value > numel(StringCell)
    Value=Value-1;
end
if ~isempty(StringCell)   
    for i=1:numel(StringCell)
        S1=strfind(StringCell{i},'ID:');
        S2=strfind(StringCell{i},'Tol:');
        S3=strfind(StringCell{i},'(Path:');
        str=StringCell{i};
        str1=regexprep(StringCell{i}, str(S1:S2-2),sprintf('ID:%03d', i));
        StringCell{i}=regexprep(str1, str(S2:S3-1), sprintf('Tol:%d', numel(StringCell)));
    end
end

set(ListboxHandle, 'String', StringCell, 'Value', Value);

function Volume3D=w_MEAN(Volume4D)
if ndims(Volume4D)==3
    Volume3D=Volume4D;
elseif ndims(Volume4D)==4
    Volume3D=mean(Volume4D, 4);
end

function Volume3D=w_STD(Volume4D)
if ndims(Volume4D)==3
    Volume3D=zeros(size(Volume4D));
elseif ndims(Volume4D)==4
    Volume3D=std(Volume4D,0,4);
end

function Volume4D=w_REPMAT(Volume3D, T)
if ndims(Volume3D)==3
    Volume4D=repmat(Volume3D, [1, 1, 1, T]);
elseif ndims(Volume3D)==4
    Volume4D=Volume3D;
end

function V=w_CORR(V1, V2, Flag)
if strcmpi(Flag, 'temporal')
    [n1, n2, n3, n4]=size(V1);
    V1=reshape(V1, [], n4);
    V2=reshape(V2, [], n4);
    V=zeros(n1*n2*n3, 1);
    for i=1:n1*n2*n3
        V(i, 1)=corr(V1(i,:)', V2(i,:)');
    end
    V=reshape(V, [n1, n2, n3]);
elseif strcmpi(Flag, 'spatial')
    if ndims(V1)==4 && ndims(V2)==4
        n4=size(V1, 4);
        V1=reshape(V1, [], n4);
        V2=reshape(V2, [], n4);
        V=zeros(n4, 1);
        for i=1:n4
            V(i, 1)=corr(V1(:, i), V2(:, i));
        end
    elseif ndims(V1)==4 && ndims(V2)==3
        n4=size(V1, 4);
        V1=reshape(V1, [], n4);
        V=zeros(n4, 1);
        for i=1:n4
            V(i, 1)=corr(V1(:, i), V2(:));
        end
    elseif ndims(V1)==3 && ndims(V2)==4
        n4=size(V2, 4);
        V1=reshape(V2, [], n4);
        V=zeros(n4, 1);
        for i=1:n4
            V(i, 1)=corr(V1(:), V2(:, i));
        end
    else
        V=corr(V1(:), V2(:));
    end
end



% --- Executes during object creation, after setting all properties.
function edFold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edFold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Accuracy.
function Accuracy_Callback(hObject, eventdata, handles)
% hObject    handle to Accuracy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Accuracy


% --- Executes on button press in Specificity.
function Specificity_Callback(hObject, eventdata, handles)
% hObject    handle to Specificity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Specificity


% --- Executes on button press in Sensibility.
function Sensibility_Callback(hObject, eventdata, handles)
% hObject    handle to Sensibility (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Sensibility


% --- Executes on button press in Feature.
function Feature_Callback(hObject, eventdata, handles)
% hObject    handle to Feature (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Feature


% --- Executes on button press in radFold.
function radFold_Callback(hObject, eventdata, handles)
set(handles.edFold,'Enable','on')
set(handles.InputFoldEntry,'Enable','off')


function MaskEntry_Callback(hObject, eventdata, handles)
% hObject    handle to MaskEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaskEntry as text
%        str2double(get(hObject,'String')) returns contents of MaskEntry as a double


% --- Executes during object creation, after setting all properties.
function MaskEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaskEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btMask.
function btMask_Callback(hObject, eventdata, handles)
% hObject    handle to btMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,FilePath] = uigetfile({'*.nii',' *.nii';...
    '*.img',' *.img';...
    '*.mat',' *.mat';...
    '*.txt',' *.txt';...
    },'Select Mask file');
if isnumeric(FilePath)
    return
end
if ~isempty(FilePath)
    MaskPath=fullfile(FilePath,filesep, FileName);
    handles.MaskPath=MaskPath;
    handles.Cfg.MaskPath=handles.MaskPath;
    set(handles.MaskEntry, 'String',MaskPath);
    guidata(hObject, handles);
end


% --- Executes on button press in bttFoldMy.
function bttFoldMy_Callback(hObject, eventdata, handles)
[FileName,FilePath] =uigetfile({'*.xlsx',' *.xlsx';...
    '*.xls',' *.xls';...
    '*.mat',' *.mat';...
    '*.txt',' *.txt';...
    },'Select fold file');
if isnumeric(FilePath)
    return
end
handles.FoldPath=[FilePath FileName];
handles.Cfg.FoldPath=handles.FoldPath;
handles.Cfg.strSelectFold='bttFoldMy';
set(handles.InputFoldEntry, 'String',[FilePath, FileName]);
handles.Cfg.InputFoldEntry=[FilePath, FileName];
set(handles.InputFoldEntry,'Enable','on')
set(handles.edFold,'Enable','off')
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of bttFoldMy


% --- Executes during object creation, after setting all properties.
function FileType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on selection change in ClassAlgorithm.
function ClassAlgorithm_Callback(hObject, eventdata, handles)
% hObject    handle to ClassAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global modelParemter;
modelParemter=handles.modelParemter;
ClassAlgorithm=get(hObject, 'value');
str=get(handles.ClassAlgorithm, 'String');
strClassAlgorithm=str{ClassAlgorithm};
selValueMethod=get(handles.popmMethod, 'value');
str=get(handles.popmMethod, 'String');
strSelectMethod=str{selValueMethod};
FeatureSelectionValue=get(handles.FeatureSelection,'value');

if  strcmp(strSelectMethod,' Regression')
    %     set(handles.OutAcc,'Enable','off');
    switch strClassAlgorithm
        case 'SVM (Support Vector Machine)'%
            SVMParemter('Visible','on')
            handles.Cfg.modelParemter=modelParemter;
            handles.modelParemter=modelParemter;
            ID=findstr(modelParemter,'-t');

            if modelParemter(ID+3)=='0'
                set(handles.FeatureSelection,'string',{'None','Lasso','Relieff','Weight'});
            else
                if FeatureSelectionValue==4
                    set(handles.FeatureSelection,'value',1);
                    set(handles.editFeatureNum,'Enable','off');
                end
                set(handles.FeatureSelection,'string',{'None','Lasso','Relieff'});
            end
            set(handles.WeightI, 'Enable','on');
            set(handles.WeightAll, 'Enable','on');
            
            if handles.Cfg.IscheckboxROC
                set(handles.checkboxROC, 'value',1);
            end
            if handles.Cfg.IsWeightAll
                set(handles.WeightAll, 'value',1);
            end
            if handles.Cfg.IsWeightI
                set(handles.WeightI, 'value',1);
            end
        case 'KNN (K-Nearest Neighbor)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value');
            if FeatureSelectionValue==4
                set(handles.FeatureSelection,'value',1)
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','Relieff'});
            NumNeighbor=inputdlg('Input the number of neighbors：','KNN',[1 40],{'11'});
            if length(NumNeighbor)
                NumNeighbor=str2num(cell2mat(NumNeighbor));
                handles.Cfg.NumNeighbor=NumNeighbor;
            else
                error('error: please click < ok > command')
            end
            set(handles.WeightI, 'Enable','off');
            set(handles.WeightAll, 'Enable','off');
            set(handles.checkboxROC, 'Enable','off');
            set(handles.checkboxSenSpeFold,'Enable','off');
            set(handles.WeightI,'value',0);
            set(handles.WeightAll,'value',0);
            set(handles.checkboxROC,'value',0);
            set(handles.checkboxSenSpeFold,'value',0);
            handles.Cfg.IsWeightI=0;
            handles.Cfg.IsWeightAll=0;
            handles.Cfg.IscheckboxSenSpeFold=0;
            handles.Cfg.IscheckboxROC=0;
        case 'RF (Random Forest)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value');
            if FeatureSelectionValue==4
                set(handles.FeatureSelection,'value',1)
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','Relieff'});
            NumTree=inputdlg('Input the number of Trees：','RandomForest',[1 40],{'500'});
            if length(NumTree)
                NumTree=str2num(cell2mat(NumTree));
                handles.Cfg.NumTree=NumTree;
            else
                error('error: please click < ok > command')
            end
            set(handles.WeightI, 'Enable','off');
            set(handles.WeightAll, 'Enable','off');
            set(handles.checkboxROC, 'Enable','off');
            set(handles.checkboxSenSpeFold,'Enable','off');
            set(handles.WeightI,'value',0);
            set(handles.WeightAll,'value',0);
            set(handles.checkboxROC,'value',0);
            set(handles.checkboxSenSpeFold,'value',0);
            handles.Cfg.IsWeightI=0;
            handles.Cfg.IsWeightAll=0;
            handles.Cfg.IscheckboxSenSpeFold=0;
            handles.Cfg.IscheckboxROC=0;
        case 'LR (Logistic Regression)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value');
            if FeatureSelectionValue==4
                set(handles.FeatureSelection,'value',1)
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','Relieff'});
            set(handles.WeightI, 'Enable','off');
            set(handles.WeightAll, 'Enable','off');
            set(handles.checkboxROC, 'Enable','off');
            set(handles.checkboxSenSpeFold,'Enable','off');
            set(handles.WeightI,'value',0);
            set(handles.WeightAll,'value',0);
            set(handles.checkboxROC,'value',0);
            set(handles.checkboxSenSpeFold,'value',0);
            handles.Cfg.IsWeightI=0;
            handles.Cfg.IsWeightAll=0;
            handles.Cfg.IscheckboxSenSpeFold=0;
            handles.Cfg.IscheckboxROC=0;
        case 'NB (Naive Bayes)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value');
            if FeatureSelectionValue==4
                set(handles.FeatureSelection,'value',1)
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','Relieff'});
            set(handles.WeightI, 'Enable','off');
            set(handles.WeightAll, 'Enable','off');
            set(handles.checkboxROC, 'Enable','off');
            set(handles.checkboxSenSpeFold,'Enable','off');
            set(handles.WeightI,'value',0);
            set(handles.WeightAll,'value',0);
            set(handles.checkboxROC,'value',0);
            set(handles.checkboxSenSpeFold,'value',0);
            handles.Cfg.IsWeightI=0;
            handles.Cfg.IsWeightAll=0;
            handles.Cfg.IscheckboxSenSpeFold=0;
            handles.Cfg.IscheckboxROC=0;
        case 'DT (Decision Tree)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value');
            if FeatureSelectionValue==4
                set(handles.FeatureSelection,'value',1)
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','Relieff'});
            set(handles.WeightI, 'Enable','off');
            set(handles.WeightAll, 'Enable','off');
            set(handles.checkboxROC, 'Enable','off');
            set(handles.checkboxSenSpeFold,'Enable','off');
            set(handles.WeightI,'value',0);
            set(handles.WeightAll,'value',0);
            set(handles.checkboxROC,'value',0);
            set(handles.checkboxSenSpeFold,'value',0);
            handles.Cfg.IsWeightI=0;
            handles.Cfg.IsWeightAll=0;
            handles.Cfg.IscheckboxSenSpeFold=0;
            handles.Cfg.IscheckboxROC=0;
        case 'LDA (Linear Discriminant Analysis)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value');
            if FeatureSelectionValue==4
                set(handles.FeatureSelection,'value',1)
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','Relieff'});
            set(handles.WeightI, 'Enable','off');
            set(handles.WeightAll, 'Enable','off');
            set(handles.checkboxROC, 'Enable','off');
            set(handles.checkboxSenSpeFold,'Enable','off');
            set(handles.WeightI,'value',0);
            set(handles.WeightAll,'value',0);
            set(handles.checkboxROC,'value',0);
            set(handles.checkboxSenSpeFold,'value',0);
            handles.Cfg.IsWeightI=0;
            handles.Cfg.IsWeightAll=0;
            handles.Cfg.IscheckboxSenSpeFold=0;
            handles.Cfg.IscheckboxROC=0;
    end
else
    switch strClassAlgorithm
        case 'SVM (Support Vector Machine)'%
%             setappdata(SVMParemter,'editP',value) 
            SVMParemter('Visible','on')
            handles.Cfg.modelParemter=modelParemter;
            handles.modelParemter=modelParemter;
            ID=findstr(modelParemter,'-t');
            LabelPath=handles.Cfg.LabelPath;
            if isempty(LabelPath)
                return
            end
            [FilePath,Name,Extension]=fileparts(LabelPath);
            if strcmp(Extension,'.xls') || strcmp(Extension,'.xlsx')
                Label=xlsread(LabelPath);
            end
            if strcmp(Extension,'.mat')
                Label=cell2mat(struct2cell(load(LabelPath)));
            end
            if strcmp(Extension,'.txt')
                Label=load(LabelPath);
            end
            ReLabel=unique(Label);%true label 0 1; 1 -1
            if modelParemter(ID+3)=='0' && length(ReLabel)==2
                set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff','Weight'});
            else               
                if FeatureSelectionValue==5
                    set(handles.FeatureSelection,'value',1);
                    set(handles.editFeatureNum,'Enable','off');
                end
                set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff'});
            end
            set(handles.WeightI, 'Enable','on');
            set(handles.WeightAll, 'Enable','on');
            set(handles.checkboxROC, 'Enable','on');
            set(handles.checkboxSenSpeFold,'Enable','on');
            %             if handles.Cfg.IscheckboxSenSpeFold
            %                 set(handles.checkboxSenSpeFold,'value',1);
            %             end
            %             if handles.Cfg.IscheckboxROC
            %                 set(handles.checkboxROC, 'value',1);
            %             end
            %             if handles.Cfg.IsWeightAll
            %                 set(handles.WeightAll, 'value',1);
            %             end
            %             if handles.Cfg.IsWeightI
            %                 set(handles.WeightI, 'value',1);
            %             end
        case 'KNN (K-Nearest Neighbor)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value')
            if FeatureSelectionValue==5
                set(handles.FeatureSelection,'value',1);
                set(handles.editFeatureNum,'Enable','off');
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff'});
            NumNeighbor=inputdlg('Input the number of neighbors：','KNN',[1 40],{'11'});
            if length(NumNeighbor)
                NumNeighbor=str2num(cell2mat(NumNeighbor));
                handles.Cfg.NumNeighbor=NumNeighbor;
            else
                error('error: please click < ok > command');
            end
            %%
%             set(handles.WeightI, 'Enable','off');
%             set(handles.WeightAll, 'Enable','off');
%             set(handles.WeightI,'value',0);
%             set(handles.WeightAll,'value',0);
%%
            set(handles.checkboxROC, 'Enable','on');
            set(handles.checkboxSenSpeFold,'Enable','on');
            handles.Cfg.IsWeightAll=0;
            handles.Cfg.IsWeightI=0;
            
        case 'RF (Random Forest)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value');
            if FeatureSelectionValue==5
                set(handles.FeatureSelection,'value',1);
                set(handles.editFeatureNum,'Enable','off');
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff'});
            NumTree=inputdlg('Input the number of Trees：','RandomForest',[1 40],{'500'});
            if length(NumTree)
                NumTree=str2num(cell2mat(NumTree));
                handles.Cfg.NumTree=NumTree;
            else
                error('error: please click < ok > command')
            end
            %%
%             set(handles.WeightI, 'Enable','off');
%             set(handles.WeightAll, 'Enable','off');
%             set(handles.WeightI,'value',0);
%             set(handles.WeightAll,'value',0);
%             handles.Cfg.IsWeightI=0;
%             handles.Cfg.IsWeightAll=0;
            %%
            set(handles.checkboxROC, 'Enable','on');
            set(handles.checkboxSenSpeFold,'Enable','on');       
            
        case 'LR (Logistic Regression)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value');
            if FeatureSelectionValue==5
                set(handles.FeatureSelection,'value',1);
                set(handles.editFeatureNum,'Enable','off');
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff'});
            %%
%             set(handles.WeightI, 'Enable','off');
%             set(handles.WeightAll, 'Enable','off');
%             set(handles.WeightI,'value',0);
%             set(handles.WeightAll,'value',0);
            %%
            set(handles.checkboxROC, 'Enable','on');
            set(handles.checkboxSenSpeFold,'Enable','on');
            %             handles.Cfg.IsWeightI=0;
            %             handles.Cfg.IsWeightAll=0;
            %             handles.Cfg.IscheckboxROC = 1;
            %             handles.Cfg.IscheckboxSenSpeFold = 1;
            
            if handles.Cfg.IscheckboxSenSpeFold
                set(handles.checkboxSenSpeFold,'value',1);
            end
            if handles.Cfg.IscheckboxROC
                set(handles.checkboxROC, 'value',1);
            end
            if handles.Cfg.IsWeightAll
                set(handles.WeightAll, 'value',1);
            end
            if handles.Cfg.IsWeightI
                set(handles.WeightI, 'value',1);
            end
        case 'NB (Naive Bayes)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value');
            if FeatureSelectionValue==5
                set(handles.FeatureSelection,'value',1);
                set(handles.editFeatureNum,'Enable','off');
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff'});
            %%
%             set(handles.WeightI, 'Enable','off');
%             set(handles.WeightAll, 'Enable','off');
%             set(handles.WeightI,'value',0);
%             set(handles.WeightAll,'value',0);
            %%
            set(handles.checkboxROC, 'Enable','on');
            set(handles.checkboxSenSpeFold,'Enable','on');
            %             handles.Cfg.IsWeightI=0;
            %             handles.Cfg.IsWeightAll=0;
            %             handles.Cfg.IscheckboxROC = 1;
            %             handles.Cfg.IscheckboxSenSpeFold = 1;
            
            if handles.Cfg.IscheckboxSenSpeFold
                set(handles.checkboxSenSpeFold,'value',1);
            end
            if handles.Cfg.IscheckboxROC
                set(handles.checkboxROC, 'value',1);
            end
            if handles.Cfg.IsWeightAll
                set(handles.WeightAll, 'value',1);
            end
            if handles.Cfg.IsWeightI
                set(handles.WeightI, 'value',1);
            end
        case 'DT (Decision Tree)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value');
            if FeatureSelectionValue==5
                set(handles.FeatureSelection,'value',1);
                set(handles.editFeatureNum,'Enable','off');
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff'});
%%
%             set(handles.WeightI, 'Enable','off');
%             set(handles.WeightAll, 'Enable','off');
%             set(handles.WeightI,'value',0);
%             set(handles.WeightAll,'value',0);
%             handles.Cfg.IsWeightI=0;
%             handles.Cfg.IsWeightAll=0;
         %%   
            set(handles.checkboxROC,'value',0);
            handles.Cfg.IscheckboxROC=0;
            handles.Cfg.IscheckboxSenSpeFold = 1;
            set(handles.checkboxROC, 'Enable','on');
            set(handles.checkboxSenSpeFold,'Enable','on');
            
            if handles.Cfg.IscheckboxSenSpeFold
                set(handles.checkboxSenSpeFold,'value',1);
            end
            if handles.Cfg.IscheckboxROC
                set(handles.checkboxROC, 'value',1);
            end
            if handles.Cfg.IsWeightAll
                set(handles.WeightAll, 'value',1);
            end
            if handles.Cfg.IsWeightI
                set(handles.WeightI, 'value',1);
            end
            
        case 'LDA (Linear Discriminant Analysis)'
            FeatureSelectionValue=get(handles.FeatureSelection,'value');
            if FeatureSelectionValue==5
                set(handles.FeatureSelection,'value',1);
                set(handles.editFeatureNum,'Enable','off');
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff'});
            %%
%             set(handles.WeightI, 'Enable','off');
%             set(handles.WeightAll, 'Enable','off');            
%             set(handles.WeightI,'value',0);
%             set(handles.WeightAll,'value',0);           
%             handles.Cfg.IsWeightI=0;
%             handles.Cfg.IsWeightAll=0;
            %%
            handles.Cfg.IscheckboxROC = 1;
            handles.Cfg.IscheckboxSenSpeFold = 1;
            set(handles.checkboxROC, 'Enable','on');
            set(handles.checkboxSenSpeFold,'Enable','on');
            if handles.Cfg.IscheckboxSenSpeFold
                set(handles.checkboxSenSpeFold,'value',1);
            end
            if handles.Cfg.IscheckboxROC
                set(handles.checkboxROC, 'value',1);
            end
            if handles.Cfg.IsWeightAll
                set(handles.WeightAll, 'value',1);
            end
            if handles.Cfg.IsWeightI
                set(handles.WeightI, 'value',1);
            end
    end
end
handles.Cfg.ClassAlgorithmStr =strClassAlgorithm;
handles.Cfg.ValueClassAlgorithm=ClassAlgorithm;

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function FeatureSelection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FeatureSelection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function PermutationNum_Callback(hObject, eventdata, handles)
handles.Cfg.PermutationNum =str2double(get(hObject,'String'));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function PermutationNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PermutationNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in FileType.
function FileType_Callback(hObject, eventdata, handles)
% hObject    handle to FileType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns FileType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileType


% --- Executes during object creation, after setting all properties.
function editFeatureNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFeatureNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function pushbuT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbuT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in FeatureSelection.
function FeatureSelection_Callback(hObject, eventdata, handles)
selValueFeature=get(hObject, 'value');
str=get(hObject, 'String');
strSelectFeature=str{selValueFeature};
handles.Cfg.FeatureSelectionStr =strSelectFeature;
handles.Cfg.selValueFeature=selValueFeature;
if ~strcmp(strSelectFeature,'None')
    set(handles.editFeatureNum,'Enable','on')
end
if strcmp(strSelectFeature,'None') 
    set(handles.editFeatureNum,'Enable','off')
end
if strcmp(strSelectFeature,'Lasso')
    set(handles.editFeatureNum,'Enable','off')
end
guidata(hObject, handles);


% --- Executes on button press in radLabel.
function radLabel_Callback(hObject, eventdata, handles)
% hObject    handle to radLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radLabel
[FileName,FilePath] =uigetfile({'*.xls',' *.xls';...
    '*.xlsx',' *.xlsx';...
    '*.mat',' *.mat';...
    '*.txt',' *.txt';...
    },'Select the file');
if isnumeric(FilePath)
    return
end
% [handles.LabelXlsDir, Name]=fileparts(FilePath);

D=fullfile(FilePath, FileName);
[Path,Name,FileType]=fileparts(D);

if strcmp(FileType,'.txt')
    LabelData=load(D);
end
if strcmp(FileType,'.mat')
    LabelData=cell2mat(struct2cell(load(D)));
end
if strcmp(FileType,'.xlsx')
    LabelData=xlsread(D);
end
if strcmp(FileType,'.xls')
    LabelData=xlsread(D);
end

handles.CreateLabel=LabelData;
handles.LabelPath=FilePath;

StringOne={sprintf('{(%s) %s', FilePath, FileName)};
guidata(hObject, handles);


% --- Executes during object deletion, before destroying properties.
function uipanel8_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to uipanel8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in FishZ.
function FishZ_Callback(hObject, eventdata, handles)
ValueFishZ=get(hObject, 'value');
str=get(hObject, 'String');
strFishZ=str{ValueFishZ};
handles.Cfg.FishZStr =strFishZ;
handles.Cfg.ValueFishZ=ValueFishZ;

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function FishZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FishZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ROC.
function ROC_Callback(hObject, eventdata, handles)
% hObject    handle to ROC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Permute.
function Permute_Callback(hObject, eventdata, handles)
% hObject    handle to Permute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Permute

sel_Permute=get(handles.Permute, 'value')

if sel_Permute
    set(handles.PermutationNum,'Enable','on')
end
% --- Executes during object creation, after setting all properties.
function Permute_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Permute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in DeleteAllfile.
function DeleteAllfile_Callback(hObject, eventdata, handles)
% hObject    handle to DeleteAllfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FileNum=numel(handles.Cfg.GroupFile);
handles.GroupFile={};
while FileNum>0
    RemoveString(handles.GroupListbox, FileNum);
    guidata(hObject, handles);
    FileNum=FileNum-1;
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function LabelEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LabelEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btLabel.
function btLabel_Callback(hObject, eventdata, handles)
% hObject    handle to btLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,FilePath] =uigetfile({'*.xlsx',' *.xlsx';...
    '*.xls',' *.xls';...
    '*.mat',' *.mat';...
    '*.txt',' *.txt';...
    },'Select label file');
if isnumeric(FilePath)
    return
end

LabelPath=[FilePath FileName];
handles.LabelPath=LabelPath;
handles.Cfg.LabelPath=handles.LabelPath;
set(handles.LabelEntry, 'String',[FilePath, FileName]);
guidata(hObject, handles);


% --- Executes on button press in btLoad.
function handles=btLoad_Callback(hObject, eventdata, handles)
global selType;
global selKernelFunction;
global strSelectMethod;
global strGama;
global strDegree;
global strCost;
global strNu;
global strCoef;
global strP;
global modelParemter;


[filename, pathname] = uigetfile({'*.mat'}, 'Load Parameters From');
if isnumeric(pathname)
    return
end
load([pathname,filename]);
handles.Cfg=Cfg;
modelParemter=Cfg.modelParemter;
modelParemter1=strrep(modelParemter,' ','');
IDStart=findstr(modelParemter1,'-s');
IDEnd=findstr(modelParemter1,'-t');
selType=modelParemter1(IDStart+2:IDEnd-1);
IDStart=findstr(modelParemter1,'-t');
IDEnd=findstr(modelParemter1,'-c');
selKernelFunction=modelParemter1(IDStart+2:IDEnd-1);
IDStart=findstr(modelParemter1,'-c');
IDEnd=findstr(modelParemter1,'-d');
strCost=modelParemter1(IDStart+2:IDEnd-1);
IDStart=findstr(modelParemter1,'-d');
IDEnd=findstr(modelParemter1,'-g');
strDegree=modelParemter1(IDStart+2:IDEnd-1);
IDStart=findstr(modelParemter1,'-g');
IDEnd=findstr(modelParemter1,'-n');
strGama=modelParemter1(IDStart+2:IDEnd-1);
IDStart=findstr(modelParemter1,'-n');
IDEnd=findstr(modelParemter1,'-r');
strNu=modelParemter1(IDStart+2:IDEnd-1);
IDStart=findstr(modelParemter1,'-r');
IDEnd=findstr(modelParemter1,'-p');
strCoef=modelParemter1(IDStart+2:IDEnd-1);
IDStart=findstr(modelParemter1,'-p');
strP=modelParemter1(IDStart+2:length(modelParemter1));

%%
%select method Classification/Regression
set(handles.popmMethod, 'value',handles.Cfg.selValueMethod);
MethodSelectionStr=handles.Cfg.MethodSelectionStr;
strSelectMethod=MethodSelectionStr;
%%
set(handles.MaskEntry,'String',Cfg.MaskPath);
set(handles.LabelEntry,'String',Cfg.LabelPath);
set(handles.GroupListbox,'String',handles.Cfg.GroupListboxStr, 'Value', 1);%GroupListbox如何load进去
set(handles.OutputEntry,'String',Cfg.OutputDir);
set(handles.PermutationNum,'String',num2str(handles.Cfg.PermutationNum));
FeatureSelectionStr=handles.Cfg.FeatureSelectionStr;
if ~strcmp(handles.Cfg.FeatureSelectionStr,'None')
    set(handles.editFeatureNum,'Enable', 'on','Value',1);
    set(handles.editFeatureNum,'String',Cfg.FeatureNum);
end

%%
%feature change
ValueFishZ=Cfg.ValueFishZ;
set(handles.FishZ, 'value',ValueFishZ);
%%
%reduce feature
set(handles.FeatureReduce,'value',handles.Cfg.stlFeatureReduce);
%feature selection
selValueFeature=handles.Cfg.selValueFeature;
%%
% algorithm
ValueClassAlgorithm=handles.Cfg.ValueClassAlgorithm;
set(handles.ClassAlgorithm,'value',ValueClassAlgorithm);
ClassAlgorithm=get(handles.ClassAlgorithm, 'value');
str=get(handles.ClassAlgorithm, 'String');
strClassAlgorithm=str{ClassAlgorithm};
if strcmp(MethodSelectionStr,' Classification')
    if strcmp(strClassAlgorithm,'SVM (Support Vector Machine)')
        
        modelParemter=handles.Cfg.modelParemter;
        ID=findstr(modelParemter,'-t');
        LabelPath=handles.Cfg.LabelPath;
        [FilePath,Name,Extension]=fileparts(LabelPath);
        if strcmp(Extension,'.xls') || strcmp(Extension,'.xlsx')
            Label=xlsread(LabelPath);
        end
        if strcmp(Extension,'.mat')
            Label=cell2mat(struct2cell(load(LabelPath)));
        end
        if strcmp(Extension,'.txt')
            Label=load(LabelPath);
        end
        ReLabel=unique(Label);%true label 0 1; 1 -1
        if modelParemter(ID+3)=='0' && length(ReLabel)==2
            set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff','Weight'});
        else
            set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff'});
        end
        
        set(handles.checkboxROC, 'Enable','on');
        set(handles.WeightI, 'Enable','on');
        set(handles.WeightAll,'Enable','on');
        get(handles.WeightAll,'Enable')
        set(handles.checkboxSenSpeFold,'Enable','on');
    else
        set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff'});
        if strcmp(strClassAlgorithm,'DT (Decision Tree)')
            set(handles.checkboxROC, 'Enable','on');
            set(handles.WeightI, 'Enable','off');
            set(handles.WeightAll,'Enable','off');
            set(handles.checkboxSenSpeFold,'Enable','off');
        else
            set(handles.checkboxROC, 'Enable','on');
            set(handles.WeightI, 'Enable','off');
            set(handles.WeightAll,'Enable','off');
            set(handles.checkboxSenSpeFold,'Enable','on');
        end
    end
    set(handles.FeatureSelection,'value',selValueFeature);
else

    set(handles.ClassAlgorithm,'string',{'SVM (Support Vector Machine)',...
        'KNN (K-Nearest Neighbor)','RF (Random Forest)','DT (Decision Tree)'});
    
    set(handles.checkboxROC, 'Enable','off');
    set(handles.checkboxSenSpeFold,'Enable','off');
    if strcmp(strClassAlgorithm,'SVM (Support Vector Machine)')
        set(handles.WeightI, 'Enable','on');
        set(handles.WeightAll,'Enable','on');
        modelParemter=handles.Cfg.modelParemter;
        ID=findstr(modelParemter,'-t');      
        if modelParemter(ID+3)=='0'
            set(handles.FeatureSelection,'string',{'None','Lasso','Relieff','Weight'});
        else
            set(handles.FeatureSelection,'string',{'None','Lasso','Relieff'});
        end
    else
        set(handles.WeightI, 'Enable','off');
        set(handles.WeightAll,'Enable','off');
        set(handles.FeatureSelection,'string',{'None','Lasso','Relieff'});
    end
    set(handles.FeatureSelection,'value',selValueFeature);
end


%%

%%
strSelectFold=handles.Cfg.strSelectFold;
if strcmp(strSelectFold, 'bttFoldMy')
    set(handles.bttFoldMy,'Enable', 'on');
    set(handles.bttFoldMy,'value', 1);
    set(handles.InputFoldEntry,'Enable','on')
    set(handles.InputFoldEntry,'String',Cfg.InputFoldEntry);
    set(handles.edFold,'Enable','off');
else
    set(handles.radFold,'Enable', 'on');
    set(handles.edFold,'Enable', 'on');
    set(handles.radFold,'value', 1);
    set(handles.edFold,'String',num2str(handles.Cfg.edFold));
    set(handles.InputFoldEntry,'Enable','off');
end
if handles.Cfg.IsSearchLight
    set(handles.editSphere,'Enable','on')
    set(handles.checkboxSearchLight,'value',1)
    set(handles.editSphere,'String',Cfg.SphereRadius);
else
    set(handles.checkboxSearchLight,'value',0)
end

if handles.Cfg.IscheckboxROC
    set(handles.checkboxROC,'value',1)
else
    set(handles.checkboxROC,'value',0)
end
if handles.Cfg.IsWeightI
    set(handles.WeightI,'value',1)
else
    set(handles.WeightI,'value',0)
end
handles.Cfg.IsWeightI
if handles.Cfg.IsWeightAll
    set(handles.WeightAll,'value',1)
else
    set(handles.WeightAll,'value',0)
end
if handles.Cfg.IscheckboxSenSpeFold
    set(handles.checkboxSenSpeFold,'value',1)
else
    set(handles.checkboxSenSpeFold,'value',0)
end
guidata(hObject, handles);

% UpdateDisplay(handles);

% --- Executes on button press in btSave.
function btSave_Callback(hObject, eventdata, handles)
Cfg=handles.Cfg;
Cfg.Data=[];
[filename, FilePath] = uiputfile({'*.mat'}, 'Save Parameters As');
if isnumeric(FilePath)
    return
end
save(['',FilePath,filename,''], 'Cfg');


% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton8


% --- Executes on button press in radiobutton9.
function radiobutton9_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton9


% --- Executes on button press in checkboxSearchLight.
function checkboxSearchLight_Callback(hObject, eventdata, handles)
global strSelectMethod;

warndlg('Searchlight only used in NITIF files')
ClassAlgorithmValue=get(handles.ClassAlgorithm,'value');
if ClassAlgorithmValue>1
    set(handles.ClassAlgorithm,'value',1)
end
if get(hObject,'Value')
    handles.Cfg.IsSearchLight = 1;
    set(handles.editSphere,'Enable','on')
    handles.Cfg.SphereRadius=get(handles.editSphere,'String');
    handles.SphereRadius=get(handles.editSphere,'String');
    set(handles.FeatureReduce,'Enable','off');
    set(handles.FishZ,'Enable','off');
    set(handles.FeatureSelection,'Enable','off');
    set(handles.WeightI, 'Enable','off');
    set(handles.WeightAll, 'Enable','off');
    set(handles.checkboxSenSpeFold,'Enable','off');
    set(handles.checkboxSenSpeFold,'value',0);
    set(handles.checkboxROC,'value',0);
    set(handles.checkboxROC, 'Enable','off');
    set(handles.editFeatureNum,'Enable','off');
    set(handles.ClassAlgorithm,'string',{'SVM (Support Vector Machine)'});
else
    handles.Cfg.IsSearchLight = 0;
    set(handles.editSphere,'Enable','off')
    set(handles.FeatureReduce,'Enable','on');
    set(handles.FishZ,'Enable','on');
    set(handles.FeatureSelection,'Enable','on');
    set(handles.WeightI, 'Enable','on');
    set(handles.WeightAll, 'Enable','on');
    set(handles.checkboxSenSpeFold,'Enable','on');
    set(handles.checkboxROC, 'Enable','on');
    if  strcmp(strSelectMethod,' Regression')
        set(handles.ClassAlgorithm,'string',{'SVM (Support Vector Machine)',...
            'KNN (K-Nearest Neighbor)','RF (Random Forest)','DT (Decision Tree)'});
    else
        set(handles.ClassAlgorithm,'string',{'SVM (Support Vector Machine)',...
            'KNN (K-Nearest Neighbor)','RF (Random Forest)','DT (Decision Tree)',...
            'LDA (Linear Discriminant Analysis)','LR (Logistic Regression)','NB (Naive Bayes)'});
    end
end
guidata(hObject, handles);



function editSphere_Callback(hObject, eventdata, handles)
handles.Cfg.SphereRadius =str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Sphere_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Sphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function editSphere_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushbuttonInputLabel.
function pushbuttonInputLabel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonInputLabel (see GCBO)

[FileName,FilePath] =uigetfile({'*.xls',' *.xls';...
    '*.xlsx',' *.xlsx';...
    '*.mat',' *.mat';...
    '*.txt',' *.txt';...
    }, 'Select random label files','MultiSelect', 'on');

if isnumeric(FilePath)
    return
end
%%
%20171202edit
[FilePath,Name,fileType]=fileparts(fullfile(FilePath,FileName{1}));
for i=1:length(FileName)
    InputLabelFile{i}=[FilePath,filesep,FileName{i}];
end
handles.InputLabelFile=InputLabelFile;
handles.Cfg.InputLabelFile=handles.InputLabelFile;

guidata(hObject, handles);

% --- Executes on button press in PermutationRun.
function PermutationRun_Callback(hObject, eventdata, handles)
global modelParemter; %

OutputDir=get(handles.OutputEntry, 'String');
handles.Cfg.OutputDir=OutputDir;
PermutationNum=get(handles.PermutationNum, 'String');
PermutationNum=str2double(PermutationNum);
handles.Cfg.PermutationNum=PermutationNum;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read data,mask,label and fold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%read label
LabelPath=get(handles.LabelEntry, 'String');
[LabPath,LabName,LabelFileType]=fileparts(LabelPath);
if strcmp(LabelFileType,'.txt')
    Label=load(LabelPath);
end
if strcmp(LabelFileType,'.mat')
    Label=cell2mat(struct2cell(load(LabelPath)));
end
if strcmp(LabelFileType,'.xlsx')
    Label=xlsread(LabelPath);
end
if strcmp(LabelFileType,'.xls')
    Label=xlsread(LabelPath);
end
handles.Label=Label;
handles.LabelPath=LabelPath;
handles.Cfg.LabelPath=handles.LabelPath;
handles.Cfg.Label=handles.Label;


%read group
GroupListboxStr=get(handles.GroupListbox,'String');
for i=1:size(GroupListboxStr,1)
    PathPos=strfind(GroupListboxStr{i,:},'Path:');
    GroupListbox=GroupListboxStr{i,:};
    GroupFile{i}=GroupListbox(PathPos+5:length(GroupListbox)-2);
end
handles.GroupFile=GroupFile;
handles.Cfg.GroupFile=GroupFile;
handles.GroupListboxStr=GroupListboxStr;
handles.Cfg.GroupListboxStr=GroupListboxStr;

%%
%read mask
MaskPath=get(handles.MaskEntry, 'String');
handles.MaskPath=MaskPath;
handles.Cfg.MaskPath=handles.MaskPath;
if ~isempty(MaskPath)
    [FilePath,Name,MaskType]=fileparts(MaskPath);
    if strcmp(MaskType,'.nii') || strcmp(MaskType,'.img')
        VMask=spm_vol(MaskPath);
        Mask=spm_read_vols(VMask);
    end
    if strcmp(MaskType,'.txt')
        Mask=load(MaskPath);
    end
    if strcmp(MaskType,'.txt') || strcmp(MaskType,'.mat')
        Mask=importdata(MaskPath);
    end
    handles.Mask.img=Mask;
    handles.Cfg.Mask.img=Mask;
    handles.Cfg.Mask.nDim1=size(Mask,1);
    handles.Cfg.Mask.nDim2=size(Mask,2);
    handles.Cfg.Mask.nDim3=size(Mask,3);
    handles.Mask.nDim1=size(Mask,1);
    handles.Mask.nDim2=size(Mask,2);
    handles.Mask.nDim3=size(Mask,3);
else
    handles.Mask.img=[];
    handles.MaskPath='';
    handles.Cfg.MaskPath='';
    handles.Cfg.Mask.img=[];
    handles.Cfg.Mask.nDim1=0;
    handles.Cfg.Mask.nDim2=0;
    handles.Cfg.Mask.nDim3=0;
    handles.Mask.nDim1=0;
    handles.Mask.nDim2=0;
    handles.Mask.nDim3=0;
end
%read data
[Data,IndexColum]=read_data(handles.Cfg);
handles.Cfg.IndexColum=IndexColum;
handles.Cfg.Data=Data;
handles.Cfg.DataNum=size(Data,2);
handles.DataNum=size(Data,2);
DataNum=size(Data,2);

%read Fold
SelectFold=get(handles.BttSelectFold, 'SelectedObject');
strSelectFold=get(SelectFold,'tag');

handles.Cfg.strSelectFold =strSelectFold;
handles.strSelectFold =strSelectFold;
if strcmp(strSelectFold,'radFold')
    KFold=str2double(get(handles.edFold, 'String'));
    handles.Cfg.edFold=KFold;
    Fold=crossvalind('kfold',size(Data,1),KFold);
    handles.Cfg.Fold=Fold;
    handles.Fold=Fold;
else
    FoldPath=get(handles.InputFoldEntry,'String');
    [FPath,FoldName,FoldFileType]=fileparts(FoldPath);
    if strcmp(FoldFileType,'.txt')
        Fold=load(FoldPath);
    end
    if strcmp(FoldFileType,'.mat')
        Fold=cell2mat(struct2cell(load(FoldPath)));
    end
    if strcmp(FoldFileType,'.xlsx')
        Fold=xlsread(FoldPath);
    end
    if strcmp(FoldFileType,'.xls')
        Fold=xlsread(FoldPath);
    end
    handles.Fold=Fold;
    handles.Cfg.Fold=handles.Fold;
    handles.FoldPath=FoldPath;
    handles.Cfg.FoldPath=handles.FoldPath;
end
%%
%20161202edit
sel_buttonPermLabel=get(handles.BttSelectPermLabel, 'SelectedObject');
strSelectPermLabel=get(sel_buttonPermLabel,'tag');
if strcmp(strSelectPermLabel,'InputLabel')
    InputLabelFile=handles.Cfg.InputLabelFile;
    RandLabel=importdata(InputLabelFile);
    handles.Cfg.RandLabel=RandLabel;
end
if strcmp(strSelectPermLabel,'CreateLabel')%create random label
    label=handles.Label;
    for iD=1:PermutationNum
        randn=randperm(length(label));
        RandLabel(:,iD)=label(randn);
    end
    if exist([OutputDir,filesep,'RandLabel.mat'],'file')
        delete([OutputDir,filesep,'RandLabel.mat']);
    end
    save([OutputDir,filesep,'RandLabel.mat'],'RandLabel');
    handles.Cfg.RandLabel=RandLabel;
end

%%
selValueMethod=get(handles.popmMethod, 'value');
str=get(handles.popmMethod, 'String');
strSelectMethod=str{selValueMethod};
if handles.Cfg.IsSearchLight
    %% if select searchlight, permutation will do
    if isempty(gcp('nocreate'))
        NumCore=inputdlg('the number of parallel core：','Core',[1 40],{'4'});
        if length(NumCore)
            NumCore=str2num(cell2mat(NumCore));
            parpool('local',NumCore);
        else
            error('error: please click ok command');
        end
    end
    tic
    disp(sprintf('>>>>>>>Permutation of SearchLight will be begin>>>>>>>>'));
    SphereRadius=str2num(handles.Cfg.SphereRadius);
    if  strcmp(strSelectMethod,' Classification')
        parfor i=1:PermutationNum
            permpath=[OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_' num2str(i)];
            if ~exist([permpath,filesep,'res_accuracy_minus_chance.mat'],'file')
                cfg = decoding_defaults;
                cfg.results.overwrite = 1;
                for j=1:length(GroupFile)
                    cfg.files.name{j}=(GroupFile{j});
                end
                cfg.files.mask =MaskPath;
                cfg.analysis = 'searchlight';
                cfg.decoding.software='libsvm';
                cfg.software=spm('ver');
                cfg.files.chunk = Fold;
                cfg.files.label =RandLabel(:,i);
                new_folder = [OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_' num2str(i)]; %each permutation will be saved in a new dirctory
                mkdir(new_folder);
                cfg.results.dir=[OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_' num2str(i)];
                cfg.design = make_design_cv(cfg);
                cfg.design.unbalanced_data = 'ok';
                cfg.searchlight.radius = SphereRadius % this will yield a searchlight radius of 12mm.
                cfg.feature_selection.check_datatrans_ok = true;
                cfg.decoding.method = 'classification';
                cfg.decoding.train.classification.model_parameters =[handles.Cfg.modelParemter ' -q'];% '-s 0 -t 0 -c 1 -b 0 -q';
                decoding(cfg);
            end
        end
        %%
        %corrected based on voxel
        load([OutputDir,filesep,'results',filesep,'searchlight',filesep,'res_accuracy_minus_chance.mat']);%the true accuracy of searchlight
        SearchLightRealResultsData=results.accuracy_minus_chance.output+50;
        for i=1:PermutationNum
            load([OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_' num2str(i),filesep,'res_accuracy_minus_chance.mat']);
            out=results.accuracy_minus_chance.output+50;
            randnout(1,i)=max(out);
        end
        %FWE corrected
        for i=1:length(SearchLightRealResultsData)
            acc=SearchLightRealResultsData(i,1);
            GreaterNumbers=0;
            for ii=1:PermutationNum
                if randnout(1,ii)>acc
                    GreaterNumbers=GreaterNumbers+1;
                end
            end
            P_value(i,1)=GreaterNumbers/PermutationNum;
        end
        save([OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_FWE_Pvalues.mat'],'P_value');
        
        PP=P_value;
        SearchLightRealResultsDataFWE=SearchLightRealResultsData;
        SearchLightRealResultsDataFWE(find(PP>0.05))=NaN;
        
        mask_index=results.mask_index;
        v=spm_vol(MaskPath);
        img=spm_read_vols(v);
        img(find(img==0))=-100;
        img(mask_index)=P_value;
        
        v1=spm_create_vol(v);
        v1.dt=[64,0];
        v1.fname=[OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_FWE_Pvalues.nii'];
        spm_write_vol(v1,img);
        
        v2=spm_create_vol(v);
        v2.fname=[OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_FWE_corrected_P_0.05_accuracy.nii'];
        v2.dt=[64,0];
        v=spm_vol(MaskPath);
        img=spm_read_vols(v);
        img(find(img==0))=-100;
        img(mask_index)=SearchLightRealResultsDataFWE;
        spm_write_vol(v2,img);
        
        SearchLightRealResultsDataFWE=SearchLightRealResultsData;       
        uncor_P_value=str2num(cell2mat(inputdlg('FWE corrected P value：','p',[1 40],{'0.001'})));

        SearchLightRealResultsDataFWE(find(P_value>uncor_P_value))=NaN;
        
        v=spm_vol(MaskPath);
        img=spm_read_vols(v);
        img(find(img==0))=-100;
        img(mask_index)=SearchLightRealResultsDataFWE;
        v1=spm_create_vol(v);
        v1.dt=[64,0];
        v1.fname=[OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_FWE_corrected_P_',num2str(uncor_P_value),'_accuracy.nii'];
        spm_write_vol(v1,img);
        
        delete(gcp)
        toc
        disp('>>>>>>>Permutation of searchlight for classification  finished>>>>>>>>');
        
    else
        parfor i=1:PermutationNum
            permpath=[OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_' num2str(i)];
            if ~exist([permpath,filesep,'res_corr.mat'],'file')
                cfg = decoding_defaults;
                cfg.results.overwrite = 1;
                for j=1:length(GroupFile)
                    cfg.files.name{j}=(GroupFile{j});
                end
                cfg.files.mask =MaskPath;
                cfg.analysis = 'searchlight';
                cfg.decoding.software='libsvm';
                cfg.software=spm('ver');
                cfg.files.chunk = Fold;
                cfg.files.label =RandLabel(:,i);
                new_folder = [OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_' num2str(i)]; %each permutation will be saved in a new dirctory
                mkdir(new_folder);
                cfg.results.dir=[OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_' num2str(i)];
                cfg.design = make_design_cv(cfg);
                cfg.design.unbalanced_data = 'ok';
                cfg.searchlight.radius = SphereRadius % this will yield a searchlight radius of 12mm.
                cfg.feature_selection.check_datatrans_ok = true;
                cfg.decoding.method = 'regression';
                cfg.decoding.train.classification.model_parameters =[handles.Cfg.modelParemter ' -q'];% '-s 0 -t 0 -c 1 -b 0 -q';
                cfg.results.output = {'corr'};
                decoding(cfg);
            end
        end
        %%
        %corrected based on voxel
        %%
        %corrected based on corr
        load([OutputDir,filesep,'results',filesep,'searchlight',filesep,'res_corr.mat']);%the true accuracy of searchlight
        SearchLightRealResultsData=results.corr.output;
        for i=1:PermutationNum
            load([OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_' num2str(i),filesep,'res_corr.mat']);
            out=results.corr.output;
            randnout(1,i)=max(out);
        end
        for i=1:length(SearchLightRealResultsData)
            acc=SearchLightRealResultsData(i,1);
            GreaterNumbers=0;
            for ii=1:PermutationNum
                if randnout(1,ii)>acc
                    GreaterNumbers=GreaterNumbers+1;
                end
            end
            P_value(i,1)=GreaterNumbers/PermutationNum;
        end
        %%
        save([OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_corr_FWE_Pvalues.mat'],'P_value');
        
        mask_index=results.mask_index;
        v=spm_vol(MaskPath);
        img=spm_read_vols(v);
        img(find(img==0))=-100;
        img(mask_index)=P_value;
     
        v1=spm_create_vol(v);
        v1.dt=[64,0];
        v1.fname=[OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_corr_FWE_Pvalues.nii'];
        spm_write_vol(v1,img);
        
        SearchLightRealResultsDataFWE=SearchLightRealResultsData;   
        SearchLightRealResultsDataFWE(find(P_value>0.05))=NaN;
        v=spm_vol(MaskPath);
        img=spm_read_vols(v);
        img(find(img==0))=-100;
        img(mask_index)=SearchLightRealResultsDataFWE;

        v1=spm_create_vol(v);
        v1.dt=[64,0];
        v1.fname=[OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_corr_FWE_corrected_P_0.05_coefficient.nii'];
        spm_write_vol(v1,img);
        
        uncor_P_value=str2num(cell2mat(inputdlg('FWE corrected P value for corr：','p',[1 40],{'0.001'})));
        SearchLightRealResultsDataFWE=SearchLightRealResultsData;  
        SearchLightRealResultsDataFWE(find(P_value>uncor_P_value))=NaN;
        v=spm_vol(MaskPath);
        img=spm_read_vols(v);
        img(find(img==0))=-100;
        img(mask_index)=SearchLightRealResultsDataFWE;
        
        v1=spm_create_vol(v);
        v1.dt=[64,0];
        v1.fname=[OutputDir,filesep,'results',filesep,'permutation',filesep,'permutation_corr_FWE_corrected_P_',num2str(uncor_P_value),'_coefficient.nii'];
        spm_write_vol(v1,img);
                
        delete(gcp)
        toc
        disp('>>>>>>>Permutation of searchlight for regression finished>>>>>>>>');
        
    end
else
    %%
    load([OutputDir,filesep,'MvpaResults.mat']);
    strSelectFeature=handles.Cfg.FeatureSelectionStr;
    handles.Cfg.FeatureSelectionStr=MvpaResults.SelectFeature;
    if ~strcmp(strSelectFeature,'None') & ~strcmp(strSelectFeature,'Lasso')
        handles.Cfg.SelectFeatureNum=MvpaResults.SelectFeatureNum;
    end
    if isempty(gcp('nocreate'))
        NumCore=inputdlg('the number of parallel core：','Core',[1 40],{'4'});
        if length(NumCore)
            NumCore=str2num(cell2mat(NumCore));
            parpool('local',NumCore);
        else
            error('error:please click < ok > command');
        end
    end
    %%
    
    %select difference algorithm
    ClassAlgorithm=get(handles.ClassAlgorithm, 'value');
    str=get(handles.ClassAlgorithm, 'String');
    strClassAlgorithm=str{ClassAlgorithm};
    if  strcmp(strSelectMethod,' Classification')
        MeanAccuracy=MvpaResults.MeanAccuracy;  %the true accurcay
        switch strClassAlgorithm
            case 'SVM (Support Vector Machine)'%
                disp(sprintf('>>>>>>>Start permutation of svm>>>>>>>>'));
                [PermutationClass,DecValues]=SVMPermutation(handles.Cfg);
                disp(sprintf('>>>>>>>Permutation of svm finished>>>>>>'));
            case 'KNN (K-Nearest Neighbor)'
                disp(sprintf('>>>>>>>Start permutation of KNN>>>>>>>>'));
                [PermutationClass,DecValues]=KNNPermutation(handles.Cfg);
                disp(sprintf('>>>>>>>Permutation of KNN finished>>>>>>>>'));
            case 'DT (Decision Tree)'
                disp(sprintf('>>>>>>>Start permutation of DecisionTree>>>>>>>>'));
                [PermutationClass,DecValues]=DecisionTreePermutation(handles.Cfg);
                disp(sprintf('>>>>>>>Permutation of DecisionTree finished>>>>>>>>'));
            case 'LR (Logistic Regression)'
                disp(sprintf('>>>>>>>Start permutation of LR>>>>>>>>'));
                [PermutationClass,DecValues]=logisticPermutation(handles.Cfg);
                disp(sprintf('>>>>>>>Permutation of LR finished>>>>>>>>'));
            case 'LDA (Linear Discriminant Analysis)'
                disp(sprintf('>>>>>>>Start permutation of LDA>>>>>>>>'));
                [PermutationClass,DecValues]=LDAPermutation(handles.Cfg);
                disp(sprintf('>>>>>>>Permutation of LDA finished>>>>>>>>'));
            case 'NB (Naive Bayes)'
                disp(sprintf('>>>>>>>Start permutation of NaiveBayes>>>>>>>>'));
                [PermutationClass,DecValues]=NaiveBayesPermutation(handles.Cfg);
                disp(sprintf('>>>>>>>Permutation of NaiveBayes finished>>>>>>>>'));
            case 'RF (Random Forest)'
                disp(sprintf('>>>>>>>Start permutation of RandomForest>>>>>>>>'));
                [PermutationClass,DecValues]=RandomForestPermutation(handles.Cfg);
                disp(sprintf('>>>>>>>Permutation of RandomForest finished>>>>>>>>'));
        end
        %%
        if ~strcmp(strSelectFeature,'None')%feature selection,PermutationClass = feature* permutation
            PermutationMeanClass=PermutationClass;
            for i=1:size(PermutationMeanClass,1)
                Permutation_FWE_voxel_p(i)=length(PermutationMeanClass(i,(find(PermutationMeanClass(i,:)>=MeanAccuracy(i)))))/PermutationNum;
                h=figure(i);
                histogram(PermutationMeanClass(i,:),[0:10:100])
                hold on
                plot([MeanAccuracy(i) MeanAccuracy(i)],[max(hist(PermutationMeanClass(i,:))) 0],'r', 'LineWidth',3);
                xlabel('Classification Accuracy','LineWidth',3);
                ylabel('Frequency','LineWidth',3);
                set(gca, 'LineWidth',2)
                if isempty(PermutationMeanClass(i,(find(PermutationMeanClass(i,:)>=MeanAccuracy(i)))))
                    title(['p<',num2str(1/PermutationNum),'(1/',num2str(PermutationNum),')'],'fontsize',15);
                else
                    title(['p=',num2str(Permutation_FWE_voxel_p(i)),'(',...
                        num2str(length(PermutationMeanClass(i,(find(PermutationMeanClass(i,:)>=MeanAccuracy(i)))))),...
                        '/',num2str(PermutationNum),')'],'fontsize',15);
                end
                saveas(h,[OutputDir,filesep,'PermutationClassAccuracy_feature_',num2str(i),'.bmp']);
                saveas(h,[OutputDir,filesep,'PermutationClassAccuracy_feature_',num2str(i),'.fig']);
            end
        else  %no feature selection,PermutationClass = permutation
            PermutationMeanClass=PermutationClass;
            Permutation_FWE_voxel_p=length(PermutationMeanClass(find(PermutationMeanClass>=MeanAccuracy)))/PermutationNum;
            h=figure;
            histogram(PermutationMeanClass,[0:10:100])
            hold on
            plot([MeanAccuracy MeanAccuracy],[max(hist(PermutationMeanClass)) 0],'r', 'LineWidth',3)
            xlabel('Classification Accuracy','LineWidth',3);
            ylabel('Frequency','LineWidth',3);
            set(gca, 'LineWidth',2)
            if isempty(PermutationMeanClass(find(PermutationMeanClass>=MeanAccuracy)))
                title(['p<',num2str(1/PermutationNum),...
                    '(1/',num2str(PermutationNum),')'],'fontsize',15);       
            else  
                title(['p=',num2str(Permutation_FWE_voxel_p),...
                    '(',num2str(length(PermutationMeanClass(find(PermutationMeanClass>=MeanAccuracy)))),...
                    '/',num2str(PermutationNum),')'],'fontsize',15);
            end
            saveas(h,[OutputDir,filesep,'PermutationClassAccuracy.bmp'])
            saveas(h,[OutputDir,filesep,'PermutationClassAccuracy.fig'])
        end
        MvpaPermResults.PermutationMeanClass=PermutationMeanClass';
        MvpaPermResults.DecValues=DecValues;
        MvpaPermResults.IndexColum=IndexColum;
        MvpaPermResults.Permutation_p=Permutation_FWE_voxel_p;
        handles.MvpaPermResults=MvpaPermResults;
        handles.Cfg.MvpaPermResults=MvpaPermResults;
        save([OutputDir,filesep,'MvpaPermResults.mat'],'MvpaPermResults');%
    else
        % regression
        MeanAccuracy=MvpaResults.r;  %the true correlation coefficient
        switch strClassAlgorithm
            case 'SVM (Support Vector Machine)'%
                disp(sprintf('>>>>>>>Start permutation of svm>>>>>>>>'));
                [rlab,plab]=SVMRegressionPermutation(handles.Cfg);
                disp(sprintf('>>>>>>>Permutation of svm finished>>>>>>>>'));
            case 'KNN (K-Nearest Neighbor)'
                disp(sprintf('>>>>>>>Start permutation of KNN>>>>>>>>'));
                [rlab,plab]=KNNRegressionPermutation(handles.Cfg);
                disp(sprintf('>>>>>>>Permutation of KNN finished>>>>>>>>'));
            case 'DT (Decision Tree)'
                disp(sprintf('>>>>>>>Start permutation of DecisionTree>>>>>>>>'));
                [rlab,plab]=DecisionTreeRegressionPermutation(handles.Cfg);
                disp(sprintf('>>>>>>>Permutation of DecisionTree finished>>>>>>>>'));
            case 'RF (Random Forest)'
                disp(sprintf('>>>>>>>Start permutation of RandomForest>>>>>>>>'));
                [rlab,plab]=RandomForestRegressionPermutation(handles.Cfg);
                disp(sprintf('>>>>>>>Permutation of RandomForest finished>>>>>>>>'));
        end
        
        %%
        switch strClassAlgorithm
            case {'SVM (Support Vector Machine)','KNN (K-Nearest Neighbor)', 'DT (Decision Tree)','RF (Random Forest)'}
                
                if ~strcmp(strSelectFeature,'None')%feature selection,PermutationClass = feature* permutation
                     PermutationMeanClass=rlab;
                    for ilab=1:size(rlab,2)
                        Permutation_FWE_voxel_p(ilab)=length(PermutationMeanClass((find(PermutationMeanClass(:,ilab)>=MeanAccuracy(ilab))),ilab))/length(PermutationMeanClass(:,ilab));
                        h=figure(ilab);
                        histogram(PermutationMeanClass(:,ilab),[-1:0.1:1])
                        hold on
                        plot([MeanAccuracy(ilab) MeanAccuracy(ilab)],[max(hist(PermutationMeanClass(:,ilab))) 0],'r', 'LineWidth',3)
                        xlabel('Correlation Coefficient','LineWidth',3);
                        ylabel('Frequency','LineWidth',3);
                              
                        if isempty(PermutationMeanClass((find(PermutationMeanClass(:,ilab)>=MeanAccuracy(ilab))),ilab))
                            title(['p<',num2str(1/PermutationNum),'(1/',num2str(PermutationNum),')'],'fontsize',15);%(length(PermutationMeanClass(i,(find(PermutationMeanClass(i,:)>MeanAccuracy(i)))))+1)
                        else
                            title(['p=',num2str(Permutation_FWE_voxel_p(ilab)),'(',...
                                num2str(length(PermutationMeanClass((find(PermutationMeanClass(:,ilab)>=MeanAccuracy(ilab))),ilab))),...
                                '/',num2str(PermutationNum),')'],'fontsize',15);
                        end
                          
                        set(gca, 'LineWidth',2)
                        saveas(h,[OutputDir,filesep,'PermutationCorrelationCoefficient_feature_',num2str(ilab),'.bmp'])
                        saveas(h,[OutputDir,filesep,'PermutationCorrelationCoefficient_feature_',num2str(ilab),'.fig'])
                    end
                else %no feature selection,PermutationClass = permutation
                    PermutationMeanClass=rlab;
                    Permutation_FWE_voxel_p=length(PermutationMeanClass((find(PermutationMeanClass>=MeanAccuracy))))/length(PermutationMeanClass);
                    h=figure;
                    histogram(PermutationMeanClass,[-1:0.1:1])
                    hold on
                    plot([MeanAccuracy MeanAccuracy],[max(hist(PermutationMeanClass)) 0],'r', 'LineWidth',3)
                    xlabel('Correlation Coefficient','LineWidth',3);
                    ylabel('Frequency','LineWidth',3);
                    if isempty(PermutationMeanClass(find(PermutationMeanClass>=MeanAccuracy)))
                        title(['p<',num2str(1/PermutationNum),...
                            '(1/',num2str(PermutationNum),')'],'fontsize',15);    
                    else 
                        title(['p=',num2str(Permutation_FWE_voxel_p),...
                            '(',num2str(length(PermutationMeanClass(find(PermutationMeanClass>=MeanAccuracy)))),...
                            '/',num2str(PermutationNum),')'],'fontsize',15);
                    end
                    
                    set(gca, 'LineWidth',2)
                    saveas(h,[OutputDir,filesep,'PermutationCorrelationCoefficient.bmp'])
                    saveas(h,[OutputDir,filesep,'PermutationCorrelationCoefficient.fig'])
                end
                MvpaPermResults.r=rlab;
                MvpaPermResults.p=plab;
                MvpaPermResults.Permutation_p=Permutation_FWE_voxel_p;
                handles.MvpaPermResults=MvpaPermResults;
                handles.Cfg.MvpaPermResults=MvpaPermResults;
                save([OutputDir,filesep,'MvpaPermResults.mat'],'MvpaPermResults');%the result of permtuation
        end
    end
end



% --- Executes on button press in WeightAll.
function WeightAll_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IsWeightAll = 1;
else
    handles.Cfg.IsWeightAll = 0;
end
guidata(hObject, handles);


% --- Executes on button press in checkboxROC.
function checkboxROC_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IscheckboxROC = 1;
else
    handles.Cfg.IscheckboxROC = 0;
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ClassAlgorithm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClassAlgorithm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CreateLabel.
function radiobuttCreateLabel_Callback(hObject, eventdata, handles)
label=handles.CreateLabel;
randn=randperm(length(label));
randLabel=label(randn);
OutputDir=get(handles.OutputEntry, 'String');
save([OutputDir,filesep,'randLabel.mat'],'randLabel');
handles.permutationLabel=randLabel;
guidata(hObject, handles);


% --- Executes on button press in InputLabel.
function InputLabel_Callback(hObject, eventdata, handles)
[FileName,FilePath] =uigetfile({'*.xls',' *.xls';...
    '*.xlsx',' *.xlsx';...
    '*.mat',' *.mat';...
    '*.txt',' *.txt';...
    },'Select the file');
%     }, 'MultiSelect', 'on');

if isnumeric(FilePath)
    return
end
InputLabelFile=[FilePath,FileName];

handles.InputLabelFile=InputLabelFile;
handles.Cfg.InputLabelFile=InputLabelFile;

guidata(hObject, handles);



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton28.
function pushbutton28_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes during object creation, after setting all properties.
function RedueDim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RedueDim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OpenFeature.
function OpenFeature_Callback(hObject, eventdata, handles)
global n;
global GroupFeatureFusionList;
global FeatureFusionMask;
if isempty(n)
    n=0;
    GroupFeatureFusionList={};
    FeatureFusionMask={};
end
n=n+1;
[GroupFeatureFusionList(n,:),FeatureFusionMask{n}]=OpenData('visible','on');
guidata(hObject, handles);


% --- Executes on button press in btVote.
function btVote_Callback(hObject, eventdata, handles)
[FileName,FilePath] = uigetfile({'*.mat',' *.mat' }, 'MultiSelect', 'on');

if isnumeric(FilePath)
    return
end
[FilePath,Name,fileType]=fileparts( fullfile(FilePath,FileName{1}));
for i=1:length(FileName)
    GroupFile{i}=[FilePath,filesep,FileName{i}];
end
disp('----------------the vote beginning----------------')
selValueMethod=get(handles.popmMethod, 'value');
str=get(handles.popmMethod, 'String');
strSelectMethod=str{selValueMethod};
if  strcmp(strSelectMethod,' Regression')
    for iD=1:length(FileName)
        load(GroupFile{iD});
        predLabel=MvpaResults.predLabel;
        if length(size(predLabel))>2 %had done feature selection
            true_label=predLabel(:,2,1);
            for i=1:size(predLabel,3)
                pred_label(:,i,iD)=predLabel(:,1,i);
            end
        else
            true_label=predLabel(:,2);
            pred_label(:,iD)=predLabel(:,1);
        end
    end
    if length(size(predLabel))>2 %had done feature selection
        PredictLabelVote=mean(pred_label,2);
        for i=1:size(PredictLabelVote,3)
            [rr,pp]=corrcoef(PredictLabelVote(:,:,i),true_label);
            r(i)=rr(1,2);
            p(i)=pp(1,2);
        end
    else
        PredictLabelVote=mean(pred_label,2);
        [rr,pp]=corrcoef(PredictLabelVote,true_label);
        r(i)=rr(1,2);
        p(i)=pp(1,2);
    end
    
    VoteResults.pred_label_vote=PredictLabelVote;
    VoteResults.pred_label=pred_label;
    VoteResults.true_label=true_label;
    VoteResults.r_vote=r;
    VoteResults.p_vote=p;
else
    
    for iD=1:length(FileName)
        load(GroupFile{iD});
        Predict_test_Label=MvpaResults.Predict_test_Label;
        if length(size(Predict_test_Label))>2 %had done feature selection
            true_label=Predict_test_Label(:,1,1);
            for i=1:size(Predict_test_Label,3)
                pred_label(:,iD,i)=Predict_test_Label(:,2,i);
            end
        else
            true_label=Predict_test_Label(:,1);
            pred_label(:,iD)=Predict_test_Label(:,2);
        end
    end
    ReLabel=unique(true_label);%true label 0 1; 1 2
    if length(size(Predict_test_Label))>2 %had done feature selection
        for iDlab=1:size(pred_label,3)
            for ilab=1:size(pred_label,1)
                if length(find(pred_label(ilab,:,iDlab)==ReLabel(1)))>=length(find(pred_label(ilab,:,iDlab)==ReLabel(2)))
                    pred_label_vote(ilab,iDlab)=ReLabel(1);
                else
                    pred_label_vote(ilab,iDlab)=ReLabel(2);
                end
            end
        end
    else
        
        for ilab=1:size(pred_label,1)
            if length(find(pred_label(ilab,:)==ReLabel(1)))>=length(find(pred_label(ilab,:)==ReLabel(2)))
                pred_label_vote(ilab,1)=ReLabel(1);
            else
                pred_label_vote(ilab,1)=ReLabel(2);
            end
        end
    end
    for i=1:size(pred_label_vote,2)
        Accuracy_vote(i)=length(find(pred_label_vote(:,i)==true_label))/length(true_label);
        sensitivity_specificity(i)=length(find(pred_label_vote(find(true_label==ReLabel(1)),i)==true_label(find(true_label==ReLabel(1)))))/length(true_label(find(true_label==ReLabel(1))));
        Sensitivity_vote(i)=length(find(pred_label_vote(find(true_label==ReLabel(2)),i)==true_label(find(true_label==ReLabel(2)))))/length(true_label(find(true_label==ReLabel(2))));
    end
    VoteResults.pred_label_vote=pred_label_vote;
    VoteResults.pred_label=pred_label;
    VoteResults.true_label=true_label;
    VoteResults.Accuracy_vote=Accuracy_vote;
    VoteResults.sensitivity_specificity=sensitivity_specificity;
    VoteResults.Sensitivity_vote=Sensitivity_vote;
end
OutputDir=get(handles.OutputEntry, 'String');
save([OutputDir,filesep,'VoteResults.mat'],'VoteResults');
disp('----------------the vote finished----------------')

% --- Executes on selection change in popmMethod.
function popmMethod_Callback(hObject, eventdata, handles)
global strSelectMethod;
selValueMethod=get(hObject, 'value');
str=get(hObject, 'String');
strSelectMethod=str{selValueMethod};
handles.Cfg.MethodSelectionStr =strSelectMethod;
handles.Cfg.selValueMethod=selValueMethod;
handles.MethodSelectionStr =strSelectMethod;
handles.selValueMethod=selValueMethod;

if  strcmp(strSelectMethod,' Regression')
    %     set(handles.ClassAlgorithm,'string','Select a method');
    ClassAlgorithmValue=get(handles.ClassAlgorithm,'value');
    if ClassAlgorithmValue>4
        set(handles.ClassAlgorithm,'value',1)
    end
    set(handles.ClassAlgorithm,'string',{'SVM (Support Vector Machine)',...
        'KNN (K-Nearest Neighbor)','RF (Random Forest)','DT (Decision Tree)'});
    %     set(handles.OutAcc, 'Enable','off');
    set(handles.checkboxSenSpeFold,'Enable','off');
    set(handles.checkboxSenSpeFold,'value',0);
    set(handles.checkboxROC,'value',0);
    set(handles.checkboxROC, 'Enable','off');
    set(handles.WeightI, 'Enable','on');
    set(handles.WeightAll, 'Enable','on');
    handles.Cfg.IscheckboxSenSpeFold=0;
    handles.Cfg.IscheckboxROC=0;
    
    
    %
    set(handles.DecisionVote,'value',1);  
    set(handles.DecisionVote,'string',{'None','Weighted Vote'});
    
else
    set(handles.ClassAlgorithm,'string',{'SVM (Support Vector Machine)',...
        'KNN (K-Nearest Neighbor)','RF (Random Forest)','DT (Decision Tree)',...
        'LDA (Linear Discriminant Analysis)','LR (Logistic Regression)','NB (Naive Bayes)'});
    set(handles.checkboxROC, 'Enable','on');
    set(handles.checkboxSenSpeFold,'Enable','on');
    %     set(handles.OutAcc, 'Enable','on');
    set(handles.WeightI, 'Enable','on');
    set(handles.WeightAll, 'Enable','on');
    
    set(handles.DecisionVote,'value',1);
    set(handles.DecisionVote,'string',{'None','Unweighted Vote','Weighted Vote'});
    
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popmMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popmMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in CreateLabel.
function CreateLabel_Callback(hObject, eventdata, handles)
% hObject    handle to CreateLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of CreateLabel


% --- Executes on key press with focus on MVPANIFigure or any of its controls.
function MVPANIFigure_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to MVPANIFigure (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in FeatureReduce.
function FeatureReduce_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureReduce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns FeatureReduce contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FeatureReduce
selFeatureReduce=get(hObject, 'value');
str=get(hObject, 'String');
strFeatureReduce=str{selFeatureReduce};
handles.Cfg.FeatureReduceStr =strFeatureReduce;
handles.Cfg.stlFeatureReduce=selFeatureReduce;
strClassAlgorithm=handles.Cfg.ClassAlgorithmStr;%select different  classification algorithms
FeatureSelectionValue=get(handles.FeatureSelection,'value');

if strcmp(strFeatureReduce,'PCA (Principal Component Analysis)')
    nDim=inputdlg('Input the % of variance kept after dimension reduction:：','nDim',[1 40],{'95'});
    if ~isempty(nDim)
        nDim=str2num(cell2mat(nDim));
        handles.Cfg.ContrPCA=nDim;
        handles.ContrPCA=nDim;
    else
        error('error: please click < ok > command');
    end
    FeatureSelectionValue=get(handles.FeatureSelection,'value');
    if FeatureSelectionValue>1
        set(handles.FeatureSelection,'value',1)
    end
    set(handles.FeatureSelection,'string',{'None'});
    set(handles.editFeatureNum,'Enable','off');
else
    if strcmp(strClassAlgorithm,'SVM (Support Vector Machine)')
        modelParemter=handles.Cfg.modelParemter;
        ID=findstr(modelParemter,'-t');
        
        LabelPath=handles.Cfg.LabelPath;
        [FilePath,Name,Extension]=fileparts(LabelPath);
        if strcmp(Extension,'.xls') || strcmp(Extension,'.xlsx')
            Label=xlsread(LabelPath);
        end
        if strcmp(Extension,'.mat')
            Label=cell2mat(struct2cell(load(LabelPath)));
        end
        if strcmp(Extension,'.txt')
            Label=load(LabelPath);
        end
        ReLabel=unique(Label);%true label 0 1; 1 -1
        if modelParemter(ID+3)=='0' && length(ReLabel)==2
            set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff','Weight'});
        else
            if FeatureSelectionValue==5
                set(handles.FeatureSelection,'value',1);
                set(handles.editFeatureNum,'Enable','off');
            end
            set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff'});
        end
        set(handles.WeightI, 'Enable','on');
        set(handles.WeightAll, 'Enable','on');
        set(handles.checkboxROC, 'Enable','on');
        set(handles.checkboxSenSpeFold,'Enable','on');
        
        
    else
        set(handles.FeatureSelection,'string',{'None','Lasso','F-Score','Relieff',});
    end
    
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function FeatureReduce_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FeatureReduce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in WeightI.
function WeightI_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IsWeightI = 1;
else
    handles.Cfg.IsWeightI = 0;
end
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function PermutationRun_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PermutationRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkboxSenSpeFold.
function checkboxSenSpeFold_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.Cfg.IscheckboxSenSpeFold = 1;
else
    handles.Cfg.IscheckboxSenSpeFold = 0;
end
guidata(hObject, handles);


% --- Executes on button press in FeatureFusion.
function FeatureFusion_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureFusion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global GroupFeatureFusionList;
global FeatureFusionMask;
global n;
disp('----------------Start feature fusion----------------')
OutputDir=get(handles.OutputEntry, 'String');
FileNum=size(GroupFeatureFusionList,2);
FeatureFusionData=[];
for iGroup=1:size(GroupFeatureFusionList,1)
    TempData1=[];
    DataTxt1=[];
    Data1=[];
    MaskPath=FeatureFusionMask{iGroup};
    GroupFile=GroupFeatureFusionList(iGroup,:);
    Mask=[];
    
    if ~isempty(MaskPath)
        [FilePath,Name,MaskType]=fileparts(MaskPath);
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
    %
    %% input group1
    [FilePath,Name,fileType]=fileparts( GroupFile{1});
    switch fileType
        case '.txt'
            for iD=1:FileNum
                DataTxt1=load(GroupFile{iD});
                nDim1=size(DataTxt1,1);
                nDim2=size(DataTxt1,2);
                if ~isempty(MaskPath)
                    %If there are mask, the data will be transformed based on the mask.
                    Mask=reshape(Mask,1,nDim1*nDim2);
                    DataTxt1=reshape(DataTxt1,1,nDim1*nDim2);
                    TempData1(iD,:)=DataTxt1.*Mask;
                else
                    TempData1(iD,:)=DataTxt1;
                end
                clear DataTxt1
            end
            Data1=[Data1;TempData1];
            clear TempData1
        case '.mat'
            for iD=1:FileNum
                DataTxt1=cell2mat(struct2cell(load(GroupFile{iD})));
                nDim1=size(DataTxt1,1);
                nDim2=size(DataTxt1,2);
                if ~isempty(MaskPath)
                    %If there are mask, the data will be transformed based on the mask.
                    DataTxt1=reshape(DataTxt1,1,nDim1*nDim2);
                    TempData1(iD,:)=DataTxt1.*Mask;
                else
                    TempData1(iD,:)=DataTxt1;
                end
                clear DataTxt1
            end
            Data1=[Data1;TempData1];
            clear TempData1
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
                [FilePath,Name,fileType]=fileparts(GroupFile{i});
                img(i,:)=image;
                clear image
            end
            Data1=[Data1;img];
            clear img
    end
    FeatureFusionData=[FeatureFusionData Data1];
end

for i=1:size(FeatureFusionData,1)
    tempData1=FeatureFusionData(i,:);
    tempData1(isnan(tempData1)) = 0;%NuN will set as 0
    FeatureFusionData(i,:)=tempData1;
end
%delete the colum with 0
FeatureFusionData(:,find(all(FeatureFusionData==0,1)))=[];
for isub=1:size(FeatureFusionData,1)
    FeatureFusion=FeatureFusionData(isub,:);
    save([OutputDir,filesep,'FeatureFusion_',sprintf('%03d',isub),'.mat'],'FeatureFusion');
end
fid = fopen([OutputDir,filesep,'FeatureFusionFile.txt'],'wt');
for iGroup=1:size(GroupFeatureFusionList,1)
    GroupFile=GroupFeatureFusionList(iGroup,:);
    fprintf(fid,'%s\n',['Group',num2str(iGroup)]);
    for iFile=1:size(GroupFile,2)
        fprintf(fid,'%s\n',GroupFile{iFile});
    end
    fprintf(fid,'\n');
end
fclose(fid);
clear global GroupFeatureFusionList;
clear global FeatureFusionMask
clear global n;
disp('----------------Feature fusion finished----------------')
% FusionData=[FusionData; Data1]
% handles.Cfg.FeatureFusionData=FusionData;



function InputFoldEntry_Callback(hObject, eventdata, handles)
% hObject    handle to InputFoldEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InputFoldEntry as text
%        str2double(get(hObject,'String')) returns contents of InputFoldEntry as a double


% --- Executes during object creation, after setting all properties.
function InputFoldEntry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InputFoldEntry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edFold_Callback(hObject, eventdata, handles)
handles.Cfg.edFold =str2double(get(hObject,'String'));



function editFeatureNum_Callback(hObject, eventdata, handles)
SelectFeatureNum=get(handles.editFeatureNum,'String');
handles.FeatureNum=SelectFeatureNum;
handles.Cfg.FeatureNum=SelectFeatureNum;
guidata(hObject, handles);


% --- Executes when user attempts to close MVPANIFigure.
function MVPANIFigure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to MVPANIFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes during object creation, after setting all properties.
function axes_logo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_logo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes_logo


% --- Executes on selection change in DecisionVote.
function DecisionVote_Callback(hObject, eventdata, handles)
% hObject    handle to DecisionVote (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns DecisionVote contents as cell array
%        contents{get(hObject,'Value')} returns selected item from DecisionVote
ValueDecisionVote=get(hObject, 'value');
str=get(hObject, 'String');
strDecisionVote=str{ValueDecisionVote};

%%
[FileName,FilePath] = uigetfile({'*.mat',' *.mat' }, 'MultiSelect', 'on');

if isnumeric(FilePath)
    return
end
[FilePath,Name,fileType]=fileparts( fullfile(FilePath,FileName{1}));
for i=1:length(FileName)
    GroupFile{i}=[FilePath,filesep,FileName{i}];
end

selValueMethod=get(handles.popmMethod, 'value');
str=get(handles.popmMethod, 'String');
strSelectMethod=str{selValueMethod};
if  strcmp(strSelectMethod,' Regression')
    disp('----------------the weighted vote beginning----------------')
    %%
    for iD=1:length(FileName)
        load(GroupFile{iD});
        predLabel=MvpaResults.predLabel;
        if length(size(predLabel))>2 %had done feature selection
            true_label=predLabel(:,2,1);
            for i=1:size(predLabel,3)
                pred_label(:,i,iD)=predLabel(:,1,i);
            end
        else
            true_label=predLabel(:,2);
            pred_label(:,iD)=predLabel(:,1);
        end
    end
    if length(size(predLabel))>2 %had done feature selection
        PredictLabelVote=mean(pred_label,3);
        for i=1:size(PredictLabelVote,2)
            [rr,pp]=corrcoef(PredictLabelVote(:,i),true_label);
            r(i)=rr(1,2);
            p(i)=pp(1,2);
        end
    else
        PredictLabelVote=mean(pred_label,2);
        [rr,pp]=corrcoef(PredictLabelVote,true_label);
        r=rr(1,2);
        p=pp(1,2);
    end
    
    WeightedVoteResults.pred_label_vote=PredictLabelVote;
    WeightedVoteResults.pred_label=pred_label;
    WeightedVoteResults.true_label=true_label;
    WeightedVoteResults.r_vote=r;
    WeightedVoteResults.p_vote=p;
    OutputDir=get(handles.OutputEntry, 'String');
    save([OutputDir,filesep,'WeightedVoteResults.mat'],'WeightedVoteResults');
    disp('----------------the weighted vote finished----------------')
else
    %for classification
    switch strDecisionVote
        %vote or decision
        case 'Unweighted Vote'
            disp('----------------the unweighted vote beginning----------------')
            for iD=1:length(FileName)
                load(GroupFile{iD});
                Predict_test_Label=MvpaResults.Predict_test_Label;
                if length(size(Predict_test_Label))>2 %had done feature selection
                    true_label=Predict_test_Label(:,1,1);
                    for i=1:size(Predict_test_Label,3)
                        pred_label(:,i,iD)=Predict_test_Label(:,2,i);
                    end
                else
                    true_label=Predict_test_Label(:,1);
                    pred_label(:,iD)=Predict_test_Label(:,2);
                end
            end

            ReLabel=unique(true_label);%true label 0 1; 1 2
            if length(size(Predict_test_Label))>2 %had done feature selection
                for iDlab=1:size(pred_label,2)%dim of feature selection
                    for ilab=1:size(pred_label,1) %dim of subject
                        nums=zeros(1,size(ReLabel,1));
                        pred_label_one=pred_label(ilab,iDlab,:);
                        for i=1:length(ReLabel)
                            nums(i)=length(find(pred_label_one==ReLabel(i)));
                        end
                        [value,Id]=max(nums);
                        pred_label_vote(ilab,iDlab)=ReLabel(Id);
                    end
                end
            else
                for ilab=1:size(pred_label,1)
                    nums=zeros(1,size(ReLabel,1));
                    pred_label_one=pred_label(ilab,:);
                    for i=1:length(ReLabel)
                        nums(i)=length(find(pred_label_one==ReLabel(i)));
                    end
                    [value,Id]=max(nums);
                    pred_label_vote(ilab,1)=ReLabel(Id);
                end
            end
            for i=1:size(pred_label_vote,2)
                Accuracy_vote(i)=length(find(pred_label_vote(:,i)==true_label))/length(true_label);
                for IDRelabel=1:length(ReLabel)
                    sensitivity_specificity(i,IDRelabel)=length(find(pred_label_vote(find(true_label==ReLabel(IDRelabel)),i)==true_label(find(true_label==ReLabel(IDRelabel)))))/length(true_label(find(true_label==ReLabel(IDRelabel))));
                end
            end
            OutputDir=get(handles.OutputEntry, 'String');

            if length(ReLabel)>2 %muti-class
                for i=1:size(pred_label_vote,2)
                    confusion_matrix=[];
                    for X=1:length(ReLabel)
                        for Y=1:length(ReLabel)
                            confusion_matrix(X,Y)=length(find(pred_label_vote(find(true_label==ReLabel(X)),i)==ReLabel(Y)));
                        end
                    end
                    ConfusionMatrix(:,:,i)=confusion_matrix;
                    %%%%%%%%%%%%
                    %%%%%%%%%%%%
                    textStrings = num2str(confusion_matrix(:),'%d');  %# Create strings from the matrix values
                    textStrings = strtrim(cellstr(textStrings));
                    
                    h=figure;
                    imagesc(confusion_matrix)
                    xlabel('Predicted Label','FontSize',15)
                    ylabel('True Label','FontSize',15)
                    set(gca,'FontSize',15)
                    set(gca, 'LineWidth',2)

                    
                    [x,y] = meshgrid(1:length(ReLabel));   %# Create x and y coordinates for the strings
                    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                        'HorizontalAlignment','center');
                    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
                    textColors = repmat(confusion_matrix(:) > midValue,1,3);  %# Choose white or black for the
                    
                    set(gca,'XTick',1:length(ReLabel),...                         %# Change the axes tick marks
                        'XTickLabel',num2cell(ReLabel),...  %#   and tick labels
                        'YTick',1:length(ReLabel),...
                        'YTickLabel',num2cell(ReLabel),...
                        'TickLength',[0 0]);
                    colorbar
                    saveas(h,[OutputDir,filesep,'ConfusionMatrix_',num2str(i),'.bmp'])
                    saveas(h,[OutputDir,filesep,'ConfusionMatrix_',num2str(i),'.fig'])

                end
                UnweightedVoteResults.confusion_matrix=ConfusionMatrix;
                
                UnweightedVoteResults.accuracy_each_class=sensitivity_specificity;
            else
                UnweightedVoteResults.specificity_sensitivity=sensitivity_specificity;
            end
            
            UnweightedVoteResults.pred_label_vote=pred_label_vote;
            UnweightedVoteResults.pred_label=pred_label;
            UnweightedVoteResults.true_label=true_label;
            UnweightedVoteResults.Accuracy_vote=Accuracy_vote;
            
            
            save([OutputDir,filesep,'UnweightedVoteResults.mat'],'UnweightedVoteResults');
            disp('----------------the unweighted vote finished----------------')
        case 'Weighted Vote'
            
%             [FileName,FilePath] = uigetfile({'*.mat',' *.mat' }, 'MultiSelect', 'on');
%             
%             if isnumeric(FilePath)
%                 return
%             end
%             [FilePath,Name,fileType]=fileparts( fullfile(FilePath,FileName{1}));
%             for i=1:length(FileName)
%                 GroupFile{i}=[FilePath,filesep,FileName{i}];
%             end
            disp('----------------the weighted vote beginning----------------')
            selValueMethod=get(handles.popmMethod, 'value');
            str=get(handles.popmMethod, 'String');
            strSelectMethod=str{selValueMethod};
            
            
            load(GroupFile{1})
            strClassAlgorithm1=MvpaResults.ClassAlgorithm;
            strClassAlgorithm2=strsplit(strClassAlgorithm1,',');
            strClassAlgorithm=strClassAlgorithm2{1};
            
            switch strClassAlgorithm

%                 case 'DT (Decision Tree)' 
%                     %|'NB (Naive Bayes)'
%                     %read value
%                     disp('----------------There is not function for the DT (Decision Tree)----------------')
%                     return
%                 case 'NB (Naive Bayes)'
%                     %|'NB (Naive Bayes)'
%                     %read value
%                     disp('----------------There is not function for the NB (Naive Bayes)----------------')
%                     return
                case 'SVM (Support Vector Machine)'
                    %%
                    for iD=1:length(FileName)
                        load(GroupFile{iD});
                        Predict_test_Label=MvpaResults.Predict_test_Label;
                        if length(size(Predict_test_Label))>2 %had done feature selection
                            true_label=Predict_test_Label(:,1,1);
                            for i=1:size(Predict_test_Label,3)                                
                                pred_label(:,:,i,iD)=Predict_test_Label(:,3:size(Predict_test_Label,2),i);
                            end
                        else
                            true_label=Predict_test_Label(:,1);
                            pred_label(:,:,iD)=Predict_test_Label(:,3:size(Predict_test_Label,2));
                        end
                    end
                    
                    ReLabel=unique(true_label,'stable');%true label 0 1; 1 2
                    if length(size(Predict_test_Label))>2 %had done feature selection
                        mean_pred_label=mean(pred_label,4);
                        for iDlab=1:size(pred_label,3)%the numbers of feature selection
                            for ilab=1:size(pred_label,1)%the numbers of subjects
                                nums=zeros(length(ReLabel),1);
                                IdMatrix=zeros(length(ReLabel));
                                mean_pred_label_value=mean_pred_label(ilab,:,iDlab);
                                ind=find(triu(ones(length(ReLabel)))-eye(length(ReLabel)));
                                IdMatrix(ind)=mean_pred_label_value;
                                for i=1:size(IdMatrix,1)-1
                                    for j=i+1:size(IdMatrix,2)
                                        if IdMatrix(i,j)>0
                                            nums(i)=nums(i)+1;
                                        else
                                            nums(j)=nums(j)+1;
                                        end
                                    end
                                end
                                [value,id_nums]=max(nums);
                                pred_label_vote(ilab,iDlab)=ReLabel(id_nums);
                            end
                        end
                    else
                        
                        mean_pred_label=mean(pred_label,3);                       
                        for ilab=1:size(pred_label,1)
                            nums=zeros(size(ReLabel,1),1);
                            IdMatrix=zeros(size(ReLabel,1));
                            mean_pred_label_value=mean_pred_label(ilab,:);
                            ind=find(triu(ones(size(ReLabel,1))-eye(size(ReLabel,1))))
                            IdMatrix(ind)=mean_pred_label_value;
                            for i=1:size(IdMatrix,1)-1
                                for j=i+1:size(IdMatrix,2)
                                    if IdMatrix(i,j)>0
                                        nums(i)=nums(i)+1;
                                    else
                                        nums(j)=nums(j)+1;
                                    end
                                end
                            end
                            [value,id_nums]=max(nums);
                            pred_label_vote(ilab,1)=ReLabel(id_nums);
                        end
                        
                    end
                    otherwise 
                        %read value
                        for iD=1:length(FileName)
                            load(GroupFile{iD});
                            Predict_test_Label=MvpaResults.Predict_test_Label;
                            if length(size(Predict_test_Label))>2 %had done feature selection
                                true_label=Predict_test_Label(:,1,1);
                                for i=1:size(Predict_test_Label,3)
                                    
                                    pred_label(:,:,i,iD)=Predict_test_Label(:,3:size(Predict_test_Label,2),i);
                                end
                            else
                                true_label=Predict_test_Label(:,1);
                                pred_label(:,:,iD)=Predict_test_Label(:,3:size(Predict_test_Label,2));
                            end
                        end
                        ReLabel=unique(true_label,'stable');%true label 0 1; 1 2
                        if length(size(Predict_test_Label))>2 %had done feature selection
                            mean_pred_label=mean(pred_label,4);
                            for iDlab=1:size(pred_label,3)
                                for ilab=1:size(pred_label,1)
                                    [value,id_mean_pred_label]=sort(mean_pred_label(ilab,:,iDlab),'descend');
                                    pred_label_vote(ilab,iDlab)=ReLabel(id_mean_pred_label(1));
                                end
                            end
                        else
                            mean_pred_label=mean(pred_label,3);
                            for ilab=1:size(pred_label,1)
                                [value,id_mean_pred_label]=sort(mean_pred_label(ilab,:),'descend');
                                pred_label_vote(ilab,1)=ReLabel(id_mean_pred_label(1));
                            end
                        end
                        
            end
            OutputDir=get(handles.OutputEntry, 'String');
            for i=1:size(pred_label_vote,2)
                Accuracy_vote(i)=length(find(pred_label_vote(:,i)==true_label))/length(true_label);
                for IDRelabel=1:length(ReLabel)
                    sensitivity_specificity(i,IDRelabel)=length(find(pred_label_vote(find(true_label==ReLabel(IDRelabel)),i)==true_label(find(true_label==ReLabel(IDRelabel)))))/length(true_label(find(true_label==ReLabel(IDRelabel))));
                end
            end
            if length(ReLabel)>2 %muti-class
                for i=1:size(pred_label_vote,2)
                    confusion_matrix=[];
                    for X=1:length(ReLabel)
                        for Y=1:length(ReLabel)
                            confusion_matrix(X,Y)=length(find(pred_label_vote(find(true_label==ReLabel(X)),i)==ReLabel(Y)));
                        end
                    end
                    ConfusionMatrix(:,:,i)=confusion_matrix;
                    %%%%%%%%%%%%
                    %%%%%%%%%%%%
                    textStrings = num2str(confusion_matrix(:),'%d');  %# Create strings from the matrix values
                    textStrings = strtrim(cellstr(textStrings));
                    
                    h=figure;
                    imagesc(confusion_matrix)
                    xlabel('Predicted Label','FontSize',15)
                    ylabel('True Label','FontSize',15)
                    set(gca,'FontSize',15)
                    set(gca, 'LineWidth',2)

                    
                    [x,y] = meshgrid(1:length(ReLabel));   %# Create x and y coordinates for the strings
                    hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                        'HorizontalAlignment','center');
                    midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
                    textColors = repmat(confusion_matrix(:) > midValue,1,3);  %# Choose white or black for the
                    
                    set(gca,'XTick',1:length(ReLabel),...                         %# Change the axes tick marks
                        'XTickLabel',num2cell(ReLabel),...  %#   and tick labels
                        'YTick',1:length(ReLabel),...
                        'YTickLabel',num2cell(ReLabel),...
                        'TickLength',[0 0]);
                    colorbar
                    saveas(h,[OutputDir,filesep,'ConfusionMatrix_',num2str(i),'.bmp'])
                    saveas(h,[OutputDir,filesep,'ConfusionMatrix_',num2str(i),'.fig'])

                end
                WeightedVoteResults.confusion_matrix=ConfusionMatrix;
                
                WeightedVoteResults.accuracy_each_class=sensitivity_specificity;
            else
                WeightedVoteResults.specificity_sensitivity=sensitivity_specificity;
            end
            WeightedVoteResults.pred_label_vote=pred_label_vote;
            WeightedVoteResults.pred_label=pred_label;
            WeightedVoteResults.true_label=true_label;
            WeightedVoteResults.Accuracy_vote=Accuracy_vote;
            
            save([OutputDir,filesep,'WeightedVoteResults.mat'],'WeightedVoteResults');
            disp('----------------the weighted vote finished----------------')
    end
end


% --- Executes during object creation, after setting all properties.
function DecisionVote_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DecisionVote (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
