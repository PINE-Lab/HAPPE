function varargout = pop_REST_reref(varargin)
% POP_REST_REREF MATLAB code for pop_REST_reref.fig
%      POP_REST_REREF, by itself, creates a new POP_REST_REREF or raises the existing
%      singleton*.
%
%      H = POP_REST_REREF returns the handle to a new POP_REST_REREF or the handle to
%      the existing singleton*.
%
%      POP_REST_REREF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POP_REST_REREF.M with the given input arguments.
%
%      POP_REST_REREF('Property','Value',...) creates a new POP_REST_REREF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pop_REST_reref_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pop_REST_reref_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pop_REST_reref

% Last Modified by GUIDE v2.5 24-Oct-2016 14:54:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @pop_REST_reref_OpeningFcn, ...
                   'gui_OutputFcn',  @pop_REST_reref_OutputFcn, ...
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


% --- Executes just before pop_REST_reref is made visible.
function pop_REST_reref_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pop_REST_reref (see VARARGIN)

% Choose default command line output for pop_REST_reref
handles.output = hObject;

handles.cfg.InputLeadfieldDir = [];          % Lead field 
handles.cfg.SelectChanns.chanlist = [];      % Channel list
handles.cfg.SelectChanns.cellchannames = []; % Channel names
handles.cfg.OrigReferFlag = 1;      % Original reference flag
handles.cfg.RetainChannsFlag = 1;   % Retain remaining channels?
handles.cfg.REST_EEG = [];          % EEG data refered to REST
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes pop_REST_reref wait for user response (see UIRESUME)
% uiwait(handles.pop_REST_reref);


% --- Outputs from this function are returned to the command line.
function varargout = pop_REST_reref_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function LeadField_edit1_Callback(hObject, eventdata, handles)
% hObject    handle to LeadField_edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LeadField_edit1 as text
%        str2double(get(hObject,'String')) returns contents of LeadField_edit1 as a double

handles.cfg.InputLeadfieldDir = get(hObject,'String');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function LeadField_edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LeadField_edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectLeadFieldFile_pushbutton1.
function SelectLeadFieldFile_pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to SelectLeadFieldFile_pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile({'*.txt;*.dat;*.xls;*.xlsx','Leadfield Files (*.txt;*.dat;*.xls;*.xlsx)';'*.*', 'All Files (*.*)';}, ...
    'Select a leadfield file');
if ~(filename==0)
    handles.cfg.InputLeadfieldDir =[pathname filename];
    set(handles.LeadField_edit1 ,'String',handles.cfg.InputLeadfieldDir);
    guidata(hObject, handles);
end

% --- Executes on button press in Run_pushbutton3.
function Run_pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to Run_pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try EEG = evalin('base','EEG');
    if ~isempty(EEG.data)
        % ----------------------------------------------
        % Load lead field --> G
        disp('----------------------------');
        disp('Loading Lead Field...');
        disp(handles.cfg.InputLeadfieldDir);
        if ~isempty(handles.cfg.InputLeadfieldDir)
            fpath = handles.cfg.InputLeadfieldDir;
            [~, ~, ext] = fileparts(fpath);
            switch ext
                case '.xlsx'
                    G = xlsread(fpath);
                case '.xls'
                    G = xlsread(fpath);
                case '.dat'
                    G = load(fpath);
                case '.txt'
                    G = load(fpath,'-ascii');
            end
            if sum(isnan(G(:)))>0
                errordlg('NaN is contained in lead field matrix!!!!','Data Error');
                return
            end
        else
            % calculating leadfield at once
            disp('Calculating leadfield based on 3-concentric spheres headmodel at once...');
            % -------------
            % use xyz coordinates in the EEG.chanlocs
            if isfield(EEG.chanlocs,'X') && isfield(EEG.chanlocs,'Y') && isfield(EEG.chanlocs,'Z')
                if ~isempty(EEG.chanlocs(1).X) && ~isempty(EEG.chanlocs(1).Y) &&~isempty(EEG.chanlocs(1).Z)
                    channs = handles.cfg.SelectChanns.chanlist; % selected channs
                    xyz_elec = zeros(length(channs),3);
                    for nc = 1:length(channs)
                        xyz_elec(nc,1) = EEG.chanlocs(channs(nc)).X;
                        xyz_elec(nc,2) = EEG.chanlocs(channs(nc)).Y;
                        xyz_elec(nc,3) = EEG.chanlocs(channs(nc)).Z;
                    end
                else
                   errordlg('EEG coordinates (EEG.chanlocs.X/Y/Z) are empty, please select lead field file OR load channel locations in EEGLAB first!!!!','Data Error');
                   return
                end
            else
                   errordlg('EEG coordinates (EEG.chanlocs.X/Y/Z) are empty, please select lead field file OR load channel locations in EEGLAB first!!!!','Data Error');
                return
            end
            % -------------------
            % load fixed dipoles and define their oritations
            % it can be defined by a file with dipole coordinates
            [ProgramPath, ~, ~] = fileparts(which('pop_REST_reref.m'));
            xyz_dipoles = load([ProgramPath,filesep,'corti869-3000dipoles.dat']);
            
            % Calculate the dipole orientations.
            xyz_dipOri           = bsxfun ( @rdivide, xyz_dipoles, sqrt ( sum ( xyz_dipoles .^ 2, 2 ) ) );
            xyz_dipOri ( 2601: 3000, 1 ) = 0;
            xyz_dipOri ( 2601: 3000, 2 ) = 0;
            xyz_dipOri ( 2601: 3000, 3 ) = 1;
            % ------------------
            % define headmodel
            headmodel        = [];
            headmodel.type   = 'concentricspheres';
            headmodel.o      = [ 0.0000 0.0000 0.0000 ];
            headmodel.r      = [ 0.8700,0.9200,1];
            headmodel.cond   = [ 1.0000,0.0125,1];
            headmodel.tissue = { 'brain' 'skull' 'scalp' };
            % -------------------
            % calculate leadfield
            [G,~] = dong_calc_leadfield3(xyz_elec,xyz_dipoles,xyz_dipOri,headmodel);
            G = G';
        end
        disp(['Lead Field Matrix: ',num2str(size(G,1)),' sources X ',num2str(size(G,2)),' channels']);
        % ----------------------------------------------
        % Load EEG data
        disp('Loading EEG data...');
        try disp(['Current data set: ',EEG.setname]);catch;end;
        if ~isempty(handles.cfg.SelectChanns.chanlist)
            channs = handles.cfg.SelectChanns.chanlist;
            if length(size(EEG.data)) == 3
                OrigData = EEG.data(channs,:);
                disp('********EEG.data is 3D epoched data!!!! Default of data demension is channels X timepoints X epochs!!!');
                disp('********Reshape to channels X timepoints');
            else
                OrigData = EEG.data(channs,:);
            end
        else
            errordlg('Please select re-referencing channels first!!!!','Data Error');
            return
        end
        disp(['EEG data: ',num2str(size(OrigData,1)),' channels X ',num2str(size(OrigData,2)),' time points'])
        % ---------------------------------------------- 
        if size(OrigData,1) == size(G,2)        
            disp('Checking original reference...');
            OrigReferFlag = handles.cfg.OrigReferFlag;   % Original reference flag
            switch OrigReferFlag
                case 1 % average
                    OrigData = OrigData - repmat(mean(OrigData),size(OrigData,1),1);
                case 2 % re-refer to average
                    OrigData = OrigData - repmat(mean(OrigData),size(OrigData,1),1);
            end
            % refer to REST
            disp('Re-referencing to REST...');
            REST_EEG =  rest_refer(OrigData,G);
            if length(size(EEG.data)) == 3
               REST_EEG = reshape(REST_EEG,size(REST_EEG,1),size(EEG.data,2),size(EEG.data,3));
               disp('********Reshape to channels X timepoints X epochs!!!!');
            end
            handles.cfg.REST_EEG = REST_EEG;       % EEG data refered to REST
            % Update handles structure
            guidata(hObject, handles);
            disp('Completed...');
        else
            errordlg('No. of Channels of lead field matrix and data are NOT equal!!!','Data Error');
            return;
        end
    else
        errordlg('EEG is empty!!!!','Data Error');
        return;
    end
catch
    errordlg('Failed to find EEG data OR other mistakes (no channel selected; channel with no location, etc...)!!!','Data Error');
    return;
end;

% --- Executes on button press in Help_pushbutton4.
function Help_pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to Help_pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
msgbox({'Re-referencing to REST steps: ';...
        '   [1] Select original reference;';...
        '       Average: average reference.';...
        '       A Fixed Electrode: reference is a fixed a electrode, e.g. Cz.';...
        '   [2] Select Channels: Select channels you want to re-reference;';...
        '       Retain un-selected channels?: if you want to keep un-selected channels in the data, check the box;';...
        '   [3] If channel location are available, Press button ''Run'' directly;';... 
        '       OR you can selcet a lead field file, and then Press button ''Run'';';...
        '    (a) Load channel location in EEGLAB';...
        '        Edit--> channel locations --> OK';...
        '        Default head model is a 3-concentric sphere head model,more details can be seen in Dong et al., 2017';...
        '        headmodel.cond = [1,0.0125,1]';...
        '        headmodel.tissue = { ''brain'',''skull'', ''scalp'' }';...
        '        headmodel.r = [0.87,0.92,1]';...
        '    (b) Select a lead field file (has been calculated and saved as *.txt/*.xls/*.xlsx/*.dat):';...
        '        The size of lead field matrix should be No. of sources X No. channels';...
        '   [4] Press button ''OK'' to save the re-referencing data to workspace (ALLEEG).In EEGLAB, click ''Datasets''-->''*_REST''';...
        },'Help');

% --------------------------------------------------------------------
function Help_1_Callback(hObject, eventdata, handles)
% hObject    handle to Help_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in Channel_listbox1.
function Channel_listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to Channel_listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Channel_listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Channel_listbox1


% --- Executes during object creation, after setting all properties.
function Channel_listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Channel_listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SelectChann_pushbutton5.
function SelectChann_pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to SelectChann_pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try EEG = evalin('base','EEG');
    if ~isempty(EEG) && length(size(EEG.data)) == 2
        if isempty(EEG.chanlocs)
            Nchanns = size(EEG.data,1);
            [chanlist,~,cellchannames] = pop_chansel(cellstr(num2str((1:Nchanns)')),'withindex','off');
        else
            [chanlist,~,cellchannames] = pop_chansel({EEG.chanlocs.labels},'withindex','on');
            Nchanns = length(chanlist);
        end
        handles.cfg.SelectChanns.chanlist = chanlist;
        handles.cfg.SelectChanns.cellchannames = cellchannames;
        set(handles.Channel_listbox1,'string',cellchannames);
        set(handles.Num_Of_Channs_text5,'string',num2str(Nchanns));
        % Update handles structure
        guidata(hObject, handles);
    elseif ~isempty(EEG) && length(size(EEG.data)) == 3
        warndlg('EEG.data is 3D epoched data!!!! Default of data demension is channel X timepoints X epochs!!!','Data');
        if isempty(EEG.chanlocs)
            Nchanns = size(EEG.data,1);
            [chanlist,~,cellchannames] = pop_chansel(cellstr(num2str((1:Nchanns)')),'withindex','off');
        else
            [chanlist,~,cellchannames] = pop_chansel({EEG.chanlocs.labels},'withindex','on');
            Nchanns = length(chanlist);
        end
        handles.cfg.SelectChanns.chanlist = chanlist;
        handles.cfg.SelectChanns.cellchannames = cellchannames;
        set(handles.Channel_listbox1,'string',cellchannames);
        set(handles.Num_Of_Channs_text5,'string',num2str(Nchanns));
        % Update handles structure
        guidata(hObject, handles);
        
    else
        errordlg('EEG is empty!!!!','Data Error');
        return;
    end
catch
    errordlg('Failed to find EEG data, Please load data first!!!','Data Error');
    return;
end;

% --- Executes on button press in OK_pushbutton6.
function OK_pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to OK_pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.cfg.REST_EEG)
    try ALLEEG = evalin('base','ALLEEG');
    catch
        ALLEEG = [];
    end
    EEG = evalin('base','EEG');
    
    flag1 = handles.cfg.RetainChannsFlag; % retain remaining channles?
    if flag1 == 0
        % -------------------------
        % ONLY save channels which re-reference to REST (Delete unselected channels)
        EEG.data = handles.cfg.REST_EEG;
        EEG.ref = 'REST';
        try
        [~, Origdatfile1, ext] = fileparts(EEG.datfile);
        EEG.datfile = [Origdatfile1,'_REST',ext];
        catch
            disp('''EEG.datfile'' may be empty!');
        end
        try
        [~, Origdatfile2, ext] = fileparts(EEG.filename);
        EEG.filename = [Origdatfile2,'_REST',ext];
        catch
            disp('''EEG.filename'' may be empty!');
        end
        try
        [~, Origdatfile3, ext] = fileparts(EEG.setname);
        EEG.setname = [Origdatfile3,'_REST',ext];
        catch
            disp('''EEG.setname'' may be empty!');
        end
        EEG.chanlocs = EEG.chanlocs(handles.cfg.SelectChanns.chanlist);
        EEG.nbchan = length(handles.cfg.SelectChanns.chanlist);
        EEG.comments = 'Re-referencing to REST';
        % ----------------
    else
        % Retain un-selected channels (e.g. EMG, EOG, bad channels etc)
        if length(size(EEG.data)) == 3
            EEG.data(handles.cfg.SelectChanns.chanlist,:,:) = handles.cfg.REST_EEG;
        else
            EEG.data(handles.cfg.SelectChanns.chanlist,:) = handles.cfg.REST_EEG;
        end
        EEG.ref = 'REST';
        try
        [~, Origdatfile1, ext] = fileparts(EEG.datfile);
        EEG.datfile = [Origdatfile1,'_REST',ext];
        catch
            disp('''EEG.datfile'' may be empty!');
        end
        try
        [~, Origdatfile2, ext] = fileparts(EEG.filename);
        EEG.filename = [Origdatfile2,'_REST',ext];
        catch
            disp('''EEG.filename'' may be empty!');
        end
        try
        [~, Origdatfile3, ext] = fileparts(EEG.setname);
        EEG.setname = [Origdatfile3,'_REST',ext];
        catch
            disp('''EEG.setname'' may be empty!');
        end
        EEG.comments = 'Re-referencing to REST';
        % --------------------
    end
    ALLEEG = [ALLEEG,EEG];
    CURRENTSET = evalin('base','CURRENTSET');
    assignin('base','CURRENTSET',length(ALLEEG));
    assignin('base','ALLEEG',ALLEEG);
    assignin('base','EEG',EEG);
    disp('');
    disp('A new dataset was created with the REST reference');
    disp('');
    eeglab('redraw');% redraw EEGLAB
    close(pop_REST_reref);
    
else
    errordlg('Re-referencing to REST is failed,please check the inputs!!!','Data Error');
    return;
end

% --- Executes during object creation, after setting all properties.
function uipanel_OrigReference_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_OrigReference (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in uipanel_OrigReference.
function uipanel_OrigReference_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_OrigReference 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
TagName = get(hObject,'Tag');
switch TagName
    case 'radiobutton1_Average'
        handles.cfg.OrigReferFlag = 1;
    case 'radiobutton2_OneFixed'
        handles.cfg.OrigReferFlag = 2;
end
% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function RESTwebsite_1_Callback(hObject, eventdata, handles)
% hObject    handle to RESTwebsite_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('http://www.neuro.uestc.edu.cn/rest/');


% --- Executes on button press in checkbox_RetainChanns.
function checkbox_RetainChanns_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_RetainChanns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_RetainChanns
value1 = get(hObject,'Value');
handles.cfg.RetainChannsFlag = value1;
% Update handles structure
guidata(hObject, handles);
