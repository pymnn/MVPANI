function varargout = OpenData(varargin)
%OPENDATA M-file for OpenData.fig
%      OPENDATA, by itself, creates a new OPENDATA or raises the existing
%      singleton*.
%
%      H = OPENDATA returns the handle to a new OPENDATA or the handle to
%      the existing singleton*.
%
%      OPENDATA('Property','Value',...) creates a new OPENDATA using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to OpenData_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      OPENDATA('CALLBACK') and OPENDATA('CALLBACK',hObject,...) call the
%      local function named CALLBACK in OPENDATA.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OpenData

% Last Modified by GUIDE v2.5 10-Jan-2019 23:00:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OpenData_OpeningFcn, ...
                   'gui_OutputFcn',  @OpenData_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before OpenData is made visible.
function OpenData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for OpenData
handles.GroupFile={};
handles.MaskPath=[];
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes OpenData wait for user response (see UIRESUME)
% uiwait(handles.figOpenFeature);
    % UIWAIT makes rest_ROIList_gui wait for user response (see UIRESUME)
    try
        uiwait(handles.figOpenFeature);
    catch
        uiresume(handles.figOpenFeature);
    end

% --- Outputs from this function are returned to the command line.
function varargout = OpenData_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.GroupFile;
varargout{2} =handles.MaskPath;
% varargout{3}=handles.FileNum;

delete(handles.figOpenFeature);


% --- Executes on button press in Group1Remove.
function Group1Remove_Callback(hObject, eventdata, handles)
Value=get(handles.GroupListbox1, 'Value');
if Value==0
    return
end
% handles.GroupCells(Value)=[];
% handles.GroupLabel(Value)=[];
handles.GroupFile(Value)=[];
handles.FileNum=numel(handles.GroupFile);
RemoveString(handles.GroupListbox1, Value);
guidata(hObject, handles);


% --- Executes on button press in buttAddGrop1.
function buttAddGrop1_Callback(hObject, eventdata, handles)
%%
GroupFile=handles.GroupFile;
FileNum=numel(handles.GroupFile);

[FileName,Path] = uigetfile({'*.nii','*.nii';...
    '*.img','*.img';...
    '*.mat','*.mat';...
    '*.txt','*.txt';...
    }, 'MultiSelect', 'on');

[Path,Name,fileType]=fileparts( fullfile(Path,FileName{1}));

for i=FileNum+1:length(FileName)+FileNum
    GroupFile{i}=[Path,filesep,FileName{i-FileNum}];
end
handles.GroupFile=GroupFile;
for i=1:size(FileName,2)
    pathNameRead=fullfile(Path,FileName{i});
    StringOne={sprintf('(%s) %s',  FileName{i},['ID:',sprintf('%03d',i),filesep,'Tol:',num2str(size(FileName,2))],['Path:', Path,filesep,FileName{i}])};

    AddString(handles.GroupListbox1, StringOne);
    guidata(hObject, handles);
end
handles.FileNum=numel(handles.GroupFile);
% handles.Cfg.GroupListboxStr=get(handles.GroupListbox,'String');
guidata(hObject, handles);



% --- Executes on selection change in GroupListbox1.
function GroupListbox1_Callback(hObject, eventdata, handles)
% hObject    handle to GroupListbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns GroupListbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GroupListbox1


% --- Executes during object creation, after setting all properties.
function GroupListbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GroupListbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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


% --- Executes on button press in DeleteAllfile1.
function DeleteAllfile1_Callback(hObject, eventdata, handles)
FileNum=numel(handles.GroupFile);
handles.GroupFile={};
while FileNum>0
    RemoveString(handles.GroupListbox1, FileNum);
    guidata(hObject, handles);
    FileNum=FileNum-1;
end
handles.FileNum=numel(handles.GroupFile);
guidata(hObject, handles);



function MaskEntry1_Callback(hObject, eventdata, handles)
% hObject    handle to MaskEntry1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaskEntry1 as text
%        str2double(get(hObject,'String')) returns contents of MaskEntry1 as a double


% --- Executes during object creation, after setting all properties.
function MaskEntry1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaskEntry1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btMask1.
function btMask1_Callback(hObject, eventdata, handles)
[FileName,Path] = uigetfile({'*.nii','*.nii';...
    '*.img','*.img';...
    '*.mat','*.mat';...
    '*.txt','*.txt';...
    },'Select the mask file');
if isnumeric(Path)
    return
end
if ~isempty(Path)    
    handles.MaskPath=[Path, FileName];
    set(handles.MaskEntry1, 'String',[Path, FileName]);
    StringOne={sprintf('{(%s) %s', Path, FileName)};
    guidata(hObject, handles);
end




% --- Executes on button press in btnDone.
function btnDone_Callback(hObject, eventdata, handles)

if ~numel(handles.GroupFile)
    error('Error: Please input data')
end
uiresume(handles.figOpenFeature);
