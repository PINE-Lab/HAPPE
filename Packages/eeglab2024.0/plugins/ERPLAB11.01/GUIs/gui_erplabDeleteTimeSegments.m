function varargout = gui_erplabDeleteTimeSegments(varargin)
% gui_erplabDeleteTimeSegments MATLAB code for gui_erplabDeleteTimeSegments.fig
%      gui_erplabDeleteTimeSegments, by itself, creates a new gui_erplabDeleteTimeSegments or raises the existing
%      singleton*.
%
%      H = gui_erplabDeleteTimeSegments returns the handle to a new gui_erplabDeleteTimeSegments or the handle to
%      the existing singleton*.
%
%      gui_erplabDeleteTimeSegments('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in gui_erplabDeleteTimeSegments.M with the given input arguments.
%
%      gui_erplabDeleteTimeSegments('Property','Value',...) creates a new gui_erplabDeleteTimeSegments or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_erplabDeleteTimeSegments_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_erplabDeleteTimeSegments_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui_erplabDeleteTimeSegments

% Last Modified by GUIDE v2.5 16-Jun-2021 12:34:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_erplabDeleteTimeSegments_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_erplabDeleteTimeSegments_OutputFcn, ...
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

% --- Executes just before gui_erplabDeleteTimeSegments is made visible.
function gui_erplabDeleteTimeSegments_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui_erplabDeleteTimeSegments (see VARARGIN)

% Choose default command line output for gui_erplabDeleteTimeSegments
handles.output = []; % hObject;


% Handle input parameters from ERPLABWORKINGMEMORY
try
    handles.timeThresholdMS             = varargin{1}{1};
    handles.startEventcodeBufferMS      = varargin{1}{2};
    handles.endEventcodeBufferMS        = varargin{1}{3};
    handles.ignoreUseEventcodes         = varargin{1}{4};
    handles.ignoreUseType               = varargin{1}{5};
    handles.displayEEG                  = varargin{1}{6};
    handles.ignoreBoundary              = varargin{1}{7}; 
    
catch
    % Default values for GUI
    handles.timeThresholdMS             = 0;
    handles.startEventcodeBufferMS      = 0;
    handles.endEventcodeBufferMS        = 0;
    handles.ignoreUseEventcodes         = [];
    handles.ignoreUseType               = 'ignore';
    handles.displayEEG                  = false;
    handles.ignoreBoundary              = false;
end


set(handles.editbox_timeThresholdMS...
    , 'String', num2str(handles.timeThresholdMS));
set(handles.editbox_startEventcodeBufferMS ...
    , 'String', num2str(handles.startEventcodeBufferMS));
set(handles.editbox_endEventcodeBufferMS ...
    , 'String', handles.endEventcodeBufferMS);
set(handles.editbox_ignoreEventCodes, ...
    'String', num2str(handles.ignoreUseEventcodes));
set(handles.checkbox_displayEEG, ...
    'Value', handles.displayEEG);
set(handles.checkbox_ignoreBE, ...
    'Value', handles.ignoreBoundary); 



% Set window title
windowTitle = ['ERPLAB ' geterplabversion() '   -   Delete Time Segments GUI'];
set(handles.gui_chassis, 'Name', windowTitle);      

handles = painterplab(handles);                     % Set color GUI
handles = setfonterplab(handles);                   % Set font size



% Run intialization procedures
initialize_gui(hObject, handles, false);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui_erplabDeleteTimeSegments wait for user response (see UIRESUME)
uiwait(handles.gui_chassis);


% --- Outputs from this function are returned to the command line.
function varargout = gui_erplabDeleteTimeSegments_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% The figure can be deleted now
delete(handles.gui_chassis);
pause(0.5)




% --- Executes on button press in pushbutton_deleteTimeSegment.
function pushbutton_deleteTimeSegment_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to pushbutton_deleteTimeSegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% mass = handles.metricdata.density * handles.metricdata.volume;
% set(handles.mass, 'String', mass);

% Command-line feedback to user
display('Deleting time segments...');


% Save the input variables to output
handles.output = {                      ...
    handles.timeThresholdMS,            ...
    handles.startEventcodeBufferMS,     ...
    handles.endEventcodeBufferMS,       ...
    handles.ignoreUseEventcodes,        ...
    handles.ignoreUseType,              ...
    handles.displayEEG,                 ...           
    handles.ignoreBoundary}; 

% Update handles structure
guidata(hObject, handles);
uiresume(handles.gui_chassis);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Command-line feedback to user
disp('User selected Cancel')

% Clear all input variables
handles.output = []; 

% Update handles structure
guidata(hObject, handles);
uiresume(handles.gui_chassis);



% --------------------------------------------------------------------
function initialize_gui(fig_handle, handles, isreset) %#ok<*INUSD>
% If the metricdata field is present and the pushbutton_cancel flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to pushbutton_cancel the data.



% Update handles structure
guidata(handles.gui_chassis, handles);




function editbox_timeThresholdMS_Callback(hObject, eventdata, handles)
% hObject    handle to editbox_timeThresholdMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editbox_timeThresholdMS as text
%        str2double(get(hObject,'String')) returns contents of editbox_timeThresholdMS as a double

% Use `str2num` (vs `str2double`) to handle both string arrray input and
% single string/character input

% returns contents of editbox_timeThresholdMS as a double
handles.timeThresholdMS = str2double(get(hObject,'String')); 

% Save the new replace channels value
guidata(hObject,handles)


function editbox_startEventcodeBufferMS_Callback(hObject, eventdata, handles)
% hObject    handle to editbox_startEventcodeBufferMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editbox_startEventcodeBufferMS as text
%        str2double(get(hObject,'String')) returns contents of editbox_startEventcodeBufferMS as a double

% returns contents of editbox_timeThresholdMS as a double
handles.startEventcodeBufferMS = str2double(get(hObject,'String')); 

% Save the new value
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function editbox_startEventcodeBufferMS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editbox_startEventcodeBufferMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editbox_endEventcodeBufferMS_Callback(hObject, eventdata, handles)
% hObject    handle to editbox_endEventcodeBufferMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editbox_endEventcodeBufferMS as text
%        str2double(get(hObject,'String')) returns contents of editbox_endEventcodeBufferMS as a double

% returns contents of editbox_endeventcodebufferms as a double
handles.endEventcodeBufferMS = str2double(get(hObject,'String')); 

% Save the new value
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function editbox_endEventcodeBufferMS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editbox_endEventcodeBufferMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editbox_ignoreEventCodes_Callback(hObject, eventdata, handles)
% hObject    handle to editbox_ignoreEventCodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ec = get(hObject,'String');

ec = colonrange2num(ec);

% Strip any non-numeric token and replace w/ whitespace (' ')
editString               = regexprep(ec, '[\D]', ' ');


handles.ignoreUseEventcodes = str2num(editString)  %#ok<ST2NM>

% Display corrected eventcode string back to GUI
set(handles.editbox_ignoreEventCodes, 'String', editString);

% Save the new replace channels value
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function editbox_ignoreEventCodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editbox_ignoreEventCodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_displayEEG.
function checkbox_displayEEG_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_displayEEG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_displayEEG
% returns contents of editbox_endeventcodebufferms as a double
handles.displayEEG = get(hObject,'Value'); 

% Save the new value
guidata(hObject,handles);


% --- Executes when user attempts to close gui_chassis.
function gui_chassis_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to gui_chassis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isequal(get(handles.gui_chassis, 'waitstatus'), 'waiting')
      %The GUI is still in UIWAIT, us UIRESUME
      handles.output = '';
      %Update handles structure
      guidata(hObject, handles);
      uiresume(handles.gui_chassis);
else
    % Hint: delete(hObject) closes the figure
    delete(hObject);
end


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%web('https://github.com/lucklab/erplab/wiki/Continuous-EEG-Preprocessing#delete-time-segments',...
%    '-browser');
web('https://github.com/lucklab/erplab/wiki/Preprocessing-Continuous-EEG-Data#delete-time-segments',...
    '-browser');



% --- Executes on selection change in popupIgnoreUseOption.
function popupIgnoreUseOption_Callback(hObject, eventdata, handles)
% hObject    handle to popupIgnoreUseOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents              = cellstr(get(hObject,'String'));  % Returns popupIgnoreUseOption contents as cell array
handles.ignoreUseType = contents{get(hObject,'Value')};  % Returns selected item from popupIgnoreUseOption

% Save the new value
guidata(hObject,handles);

% Hints: contents = cellstr(get(hObject,'String')) returns popupIgnoreUseOption contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupIgnoreUseOption


% --- Executes during object creation, after setting all properties.
function popupIgnoreUseOption_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupIgnoreUseOption (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_ignoreBE.
function checkbox_ignoreBE_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ignoreBE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ignoreBoundary = get(hObject,'Value'); 
guidata(hObject,handles);

% Hint: get(hObject,'Value') returns toggle state of checkbox_ignoreBE
