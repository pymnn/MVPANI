function varargout = SVMParemter(varargin)
% SVMPAREMTER M-file for SVMParemter.fig
%      SVMPAREMTER, by itself, creates a new SVMPAREMTER or raises the existing
%      singleton*.
%
%      H = SVMPAREMTER returns the handle to a new SVMPAREMTER or the handle to
%      the existing singleton*.
%
%      SVMPAREMTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SVMPAREMTER.M with the given input arguments.
%
%      SVMPAREMTER('Property','Value',...) creates a new SVMPAREMTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SVMParemter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SVMParemter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SVMParemter

% Last Modified by GUIDE v2.5 22-Apr-2019 19:19:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SVMParemter_OpeningFcn, ...
                   'gui_OutputFcn',  @SVMParemter_OutputFcn, ...
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


% --- Executes just before SVMParemter is made visible.
function SVMParemter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SVMParemter (see VARARGIN)

% Choose default command line output for SVMParemter
% handles.
handles.output = hObject;
global selType;
global selTypeNum;
global selKernelFunction;
global strSelectMethod;
global strGama;
global strDegree;
global strCost;
global strNu;
global strCoef;
global strP;
global modelParemter;
%%
modelParemter1=strrep(modelParemter,' ','');
IDStart=findstr(modelParemter1,'-s');
IDEnd=findstr(modelParemter1,'-t');
selType=str2num(modelParemter1(IDStart+2:IDEnd-1));
IDStart=findstr(modelParemter1,'-t');
IDEnd=findstr(modelParemter1,'-c');
selKernelFunction=str2num(modelParemter1(IDStart+2:IDEnd-1));
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


if isempty(selKernelFunction)
    set(handles.KernelFunction, 'value',1);
    selKernelFunction=1;
else
    set(handles.KernelFunction, 'value',selKernelFunction+1);
end

if isempty(strGama)
    set(handles.editGama,'string',0.1);
else
    set(handles.editGama,'string',strGama);
end

if isempty(strDegree)
    set(handles.editDegree,'string',3);
else
    set(handles.editDegree,'string',strDegree);
end

if isempty(strCost)
    set(handles.editCost,'string',1);
else
    set(handles.editCost,'string',strCost);
end

if isempty(strNu)
    set(handles.editNu,'string',0.5);
else
    set(handles.editNu,'string',strNu);
end

if isempty(strCoef)
    set(handles.editCoef,'string',0);
else
    set(handles.editCoef,'string',strCoef);
end

if isempty(strP)
    set(handles.editP,'string',0.1);
else
    set(handles.editP,'string',strP);
end

if isempty(selType)
    set(handles.Type,'value',1);
    selTypeNum=1;
else
    selTypeNum=selType+1;
end
if  strcmp(strSelectMethod,' Regression')
    if selTypeNum<=4
        set(handles.Type,'value',1);
        selTypeNum=1;
    else
        set(handles.Type,'value',selTypeNum-3);
        selTypeNum=get(handles.Type,'value');
    end
    set(handles.Type,'String',{'3--e-SVR','4--v-SVR'});
    
else
    if selTypeNum>3
        set(handles.Type,'value',1);
        selTypeNum=1;
    else
        set(handles.Type,'value',selTypeNum);
        selTypeNum=get(handles.Type,'value');
    end
    set(handles.Type,'String',{'0--C-SVC','1--v-SVC','2--one-class-SVM'});
end
set(handles.KernelFunction, 'String',{'0--linear', '1--polynomial','2--RBF','3--sigmoid'});
%%
%Enable 
% selType
% TypeString=get(handles.Type,'string');
% selTypeString=TypeString(selType);
selTypeString=[];
% 

if  strcmp(strSelectMethod,' Regression')
    switch selTypeNum
        case 1
            selTypeString='3--e-SVR';
        case 2
            selTypeString='4--v-SVR';
    end
else
    switch selTypeNum
        case 1
            selTypeString='0--C-SVC';
        case 2
            selTypeString='1--v-SVC';
        case 3
            selTypeString='2--one-class-SVM';
    end
end

if  strcmp(selTypeString,'0--C-SVC')|| strcmp(selTypeString,'3--e-SVR')  || strcmp(selTypeString,'4--v-SVR')   
    set(handles.editCost,'Enable','on');    
end

if  strcmp(selTypeString,'1--v-SVC') || strcmp(selTypeString,'2--one-class-SVM')  || strcmp(selTypeString,'4--v-SVR')   
    set(handles.editNu,'Enable','on');    
    set(handles.editCost,'Enable','off'); 
end

if  strcmp(selTypeString,'3--e-SVR')    
    set(handles.editP,'Enable','on');  
    set(handles.editCost,'Enable','off'); 
end

switch selKernelFunction
    case 1
        set(handles.editDegree,'Enable','on');
        set(handles.editCoef,'Enable','on');
        set(handles.editGama,'Enable','on');
    case 2
        set(handles.editGama,'Enable','on');
    case 3
        set(handles.editGama,'Enable','on');
        set(handles.editCoef,'Enable','on');
end


try
    uiwait(handles.SVMfigure);
catch
    uiresume(handles.SVMfigure);
end
guidata(hObject, handles);
% UIWAIT makes SVMParemter wait for user response (see UIRESUME)
% uiwait(handles.SVMfigure);


% --- Outputs from this function are returned to the command line.
function varargout = SVMParemter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.SVMfigure);



% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Type.
function Type_Callback(hObject, eventdata, handles)
% hObject    handle to Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Type
%the type of svm
global selType;
global selKernelFunction;
selType=get(handles.Type, 'value');
str=get(handles.Type, 'String');
strType=str{selType};
switch strType
    case '0--C-SVC'%
        set(handles.editCost,'Enable','on','Value',1);
        set(handles.editNu,'Enable','off','Value',0);
        set(handles.editP,'Enable','off','Value',0);
    case '1--v-SVC'
         set(handles.editNu,'Enable','on','Value',1);
         set(handles.editCost,'Enable','off','Value',0);
         set(handles.editP,'Enable','off','Value',0);
    case '2--one-class-SVM'
        set(handles.editNu,'Enable','on','Value',1);
        set(handles.editCost,'Enable','off','Value',0);
        set(handles.editP,'Enable','off','Value',0);
    case '3--e-SVR'
        set(handles.editCost,'Enable','on','Value',1);
        set(handles.editP,'Enable','on','Value',1);
        set(handles.editNu,'Enable','off','Value',0);
    case '4--v-SVR'    
        set(handles.editCost,'Enable','on','Value',1);
        set(handles.editNu,'Enable','on','Value',1);
        set(handles.editP,'Enable','off','Value',0);
end


% --- Executes during object creation, after setting all properties.
function Type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in KernelFunction.
function KernelFunction_Callback(hObject, eventdata, handles)
global selKernelFunction;
selKernelFunction=get(handles.KernelFunction, 'value');
str=get(handles.KernelFunction, 'String');
strKernelFunction=str{selKernelFunction};
switch selKernelFunction
    case 1
        set(handles.editDegree, 'Enable','off');
        set(handles.editGama, 'Enable','off');
        set(handles.editCoef, 'Enable','off');
    case 2
        set(handles.editDegree, 'Enable','on');
        set(handles.editGama, 'Enable','on');
        set(handles.editCoef, 'Enable','on');
    case 3
        set(handles.editDegree, 'Enable','off');
        set(handles.editGama, 'Enable','on');
    case 4
        set(handles.editDegree, 'Enable','off');
        set(handles.editGama, 'Enable','on');
        set(handles.editCoef, 'Enable','on');
end


% --- Executes on button press in BtnParameter.
function BtnParameter_Callback(hObject, eventdata, handles)
% hObject    handle to BtnParameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global modelParemter;
global strGama;
global strDegree;
global strCost;
global strNu;
global strCoef;
global strP;
global strType;
%select the tpye of svm

selType=get(handles.Type, 'value');
strSVM_type=get(handles.Type, 'String');
strType=strSVM_type{selType};
%select kernal function
selKernelFunction=get(handles.KernelFunction, 'value');
str=get(handles.KernelFunction, 'String');
strKernelFunction=str{selKernelFunction};

%select other paremeter
strGama=get(handles.editGama,'String');
strDegree=get(handles.editDegree,'String');
strCost=get(handles.editCost,'String');

strNu=get(handles.editNu,'String');
strCoef=get(handles.editCoef,'String');
strP=get(handles.editP,'String');
if length(strSVM_type)==3
    modelParemter=['-s ', num2str(selType-1), ' -t ', num2str(selKernelFunction-1),...
        ' -c ',strCost,'-d ',strDegree,' -g ',strGama, ' -n ',strNu,' -r ',strCoef,' -p ',strP];
else
    modelParemter=['-s ', num2str(selType-1)+3, ' -t ', num2str(selKernelFunction-1),...
        ' -c ',strCost,'-d ',strDegree,' -g ',strGama, ' -n ',strNu,' -r ',strCoef,' -p ',strP];
end

uiresume(handles.SVMfigure);


% --- Executes on selection change in Degree.
function Degree_Callback(hObject, eventdata, handles)
% hObject    handle to Degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Degree contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Degree


% --- Executes during object creation, after setting all properties.
function Degree_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Degree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Gama.
function Gama_Callback(hObject, eventdata, handles)
% hObject    handle to Gama (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Gama contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Gama


% --- Executes during object creation, after setting all properties.
function Gama_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Gama (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function SVMfigure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SVMfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function editDegree_Callback(hObject, eventdata, handles)
% hObject    handle to editDegree (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDegree as text
%        str2double(get(hObject,'String')) returns contents of editDegree as a double


% --- Executes during object creation, after setting all properties.
function editCost_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCost (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function editCoef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCoef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function editNu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCoef_Callback(hObject, eventdata, handles)
selKernelFunction=get(handles.KernelFunction, 'value');
str=get(handles.KernelFunction, 'String');
strKernelFunction=str{selKernelFunction};
switch selKernelFunction
    case 2
        set(handles.editDegree, 'Enable','on');
        set(handles.editGama, 'Enable','on');
        set(handles.editCoef, 'Enable','on');
end
