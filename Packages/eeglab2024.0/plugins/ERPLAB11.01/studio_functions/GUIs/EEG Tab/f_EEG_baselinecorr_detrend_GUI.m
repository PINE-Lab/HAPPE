%Author: Guanghui ZHANG--zhang.guanghui@foxmail.com
%Center for Mind and Brain
%University of California, Davis
%Davis, CA, USA
% 2024

% EEGLAB Studio

function varargout = f_EEG_baselinecorr_detrend_GUI(varargin)

global observe_EEGDAT;
addlistener(observe_EEGDAT,'eeg_two_panels_change',@eeg_two_panels_change);
addlistener(observe_EEGDAT,'count_current_eeg_change',@count_current_eeg_change);
addlistener(observe_EEGDAT,'Reset_eeg_panel_change',@Reset_eeg_panel_change);


%%---------------------------gui-------------------------------------------
[version reldate,ColorB_def,ColorF_def,errorColorF_def] = geterplabstudiodef;
if nargin == 0
    fig = figure(); % Parent figure
    EEG_basecorr_detrend_box = uiextras.BoxPanel('Parent', fig, 'Title', 'Baseline Correction & Linear Detrend (Epoched EEG)', 'Padding', 5,'BackgroundColor',ColorB_def); % Create boxpanel
elseif nargin == 1
    EEG_basecorr_detrend_box = uiextras.BoxPanel('Parent', varargin{1}, 'Title', 'Baseline Correction & Linear Detrend (Epoched EEG)', 'Padding', 5,'BackgroundColor',ColorB_def);
else
    EEG_basecorr_detrend_box = uiextras.BoxPanel('Parent', varargin{1}, 'Title', 'Baseline Correction & Linear Detrend (Epoched EEG)', 'Padding', 5, 'FontSize', varargin{2},'BackgroundColor',ColorB_def);
end

gui_eeg_blc_dt = struct();
try
    FonsizeDefault = varargin{2};
catch
    FonsizeDefault = [];
end
if isempty(FonsizeDefault)
    FonsizeDefault = f_get_default_fontsize();
end
EEG_blc_dt_gui(FonsizeDefault);
varargout{1} = EEG_basecorr_detrend_box;


    function EEG_blc_dt_gui(FonsizeDefault)
        [version reldate,ColorB_def,ColorF_def,errorColorF_def] = geterplabstudiodef;
        
        Enable_label = 'off';
        gui_eeg_blc_dt.blc_dt = uiextras.VBox('Parent',EEG_basecorr_detrend_box,'Spacing',1,'BackgroundColor',ColorB_def);
        
        %%Measurement type
        gui_eeg_blc_dt.blc_dt_type_title = uiextras.HBox('Parent',  gui_eeg_blc_dt.blc_dt,'Spacing',1,'BackgroundColor',ColorB_def);
        uicontrol('Style', 'text','Parent', gui_eeg_blc_dt.blc_dt_type_title,...
            'String','Type:','FontWeight','bold','FontSize',FonsizeDefault ,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.blc_dt_option = uiextras.HBox('Parent',  gui_eeg_blc_dt.blc_dt,'Spacing',1,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.blc = uicontrol('Style', 'radiobutton','Parent', gui_eeg_blc_dt.blc_dt_option,...
            'String','Baseline Correction','callback',@baseline_correction_EEG,'Value',1,'Enable',Enable_label,'FontSize',FonsizeDefault,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.blcp.KeyPressFcn= @EEG_blcorrdetrend_presskey;
        gui_eeg_blc_dt.dt = uicontrol('Style', 'radiobutton','Parent', gui_eeg_blc_dt.blc_dt_option,...
            'String','Linear detrend','callback',@detrend_EEG,'Value',0,'Enable',Enable_label,'FontSize',FonsizeDefault,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.dt.KeyPressFcn= @EEG_blcorrdetrend_presskey;
        gui_eeg_blc_dt.EEGTab_baseline_detrend{1} = gui_eeg_blc_dt.blc.Value;
        
        %%Baseline period: Pre, post whole custom
        gui_eeg_blc_dt.blc_dt_baseline_period_title = uiextras.HBox('Parent',  gui_eeg_blc_dt.blc_dt,'Spacing',1,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.blc_dt_title = uicontrol('Style', 'text','Parent', gui_eeg_blc_dt.blc_dt_baseline_period_title,...
            'String','Baseline Period:','FontWeight','bold','FontSize',FonsizeDefault,'BackgroundColor',ColorB_def);
        
        gui_eeg_blc_dt.blc_dt_bp_option = uiextras.HBox('Parent',  gui_eeg_blc_dt.blc_dt,'Spacing',1,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.pre = uicontrol('Style', 'radiobutton','Parent', gui_eeg_blc_dt.blc_dt_bp_option,...
            'String','Pre','callback',@pre_EEG,'Value',1,'Enable',Enable_label,'FontSize',FonsizeDefault,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.pre.KeyPressFcn= @EEG_blcorrdetrend_presskey;
        gui_eeg_blc_dt.post = uicontrol('Style', 'radiobutton','Parent', gui_eeg_blc_dt.blc_dt_bp_option,...
            'String','Post','callback',@post_EEG,'Value',0,'Enable',Enable_label,'FontSize',FonsizeDefault,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.post.KeyPressFcn= @EEG_blcorrdetrend_presskey;
        gui_eeg_blc_dt.whole = uicontrol('Style', 'radiobutton','Parent', gui_eeg_blc_dt.blc_dt_bp_option,...
            'String','Whole','callback',@whole_EEG,'Value',0,'Enable',Enable_label,'FontSize',FonsizeDefault,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.whole.KeyPressFcn= @EEG_blcorrdetrend_presskey;
        gui_eeg_blc_dt.blc_dt_bp_option_cust = uiextras.HBox('Parent',  gui_eeg_blc_dt.blc_dt,'Spacing',1,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.custom = uicontrol('Style', 'radiobutton','Parent', gui_eeg_blc_dt.blc_dt_bp_option_cust,...
            'String','Custom (ms) [start stop]','callback',@custom_EEG,'Value',0,'Enable',Enable_label,'FontSize',FonsizeDefault,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.custom.KeyPressFcn= @EEG_blcorrdetrend_presskey;
        gui_eeg_blc_dt.custom_edit = uicontrol('Style', 'edit','Parent', gui_eeg_blc_dt.blc_dt_bp_option_cust,...
            'String','','callback',@custom_edit,'Enable',Enable_label,'FontSize',FonsizeDefault);
        gui_eeg_blc_dt.custom_edit.KeyPressFcn= @EEG_blcorrdetrend_presskey;
        set(gui_eeg_blc_dt.blc_dt_bp_option_cust, 'Sizes',[160  100]);
        if gui_eeg_blc_dt.pre.Value==1
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2} = 1;
        elseif gui_eeg_blc_dt.post.Value==1
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2} = 2;
        elseif gui_eeg_blc_dt.whole.Value==1
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2} = 3;
        elseif gui_eeg_blc_dt.custom.Value==1
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2} = str2num(gui_eeg_blc_dt.custom_edit.String);
        else
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2} =1;
        end
        %%Bin and channels selection
        gui_eeg_blc_dt.blc_dt_bin_chan_title = uiextras.HBox('Parent',  gui_eeg_blc_dt.blc_dt,'Spacing',1,'BackgroundColor',ColorB_def);
        uicontrol('Style', 'text','Parent', gui_eeg_blc_dt.blc_dt_bin_chan_title,...
            'String','Chan Selection:','FontWeight','bold','FontSize',FonsizeDefault,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.blc_bin_chan_option = uiextras.HBox('Parent',  gui_eeg_blc_dt.blc_dt,'Spacing',1,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.all_bin_chan = uicontrol('Style', 'radiobutton','Parent', gui_eeg_blc_dt.blc_bin_chan_option,...
            'String','All (Recommended)','callback',@All_bin_chan,'Value',1,'Enable',Enable_label,'FontSize',FonsizeDefault,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.all_bin_chan.KeyPressFcn= @EEG_blcorrdetrend_presskey;
        gui_eeg_blc_dt.Selected_bin_chan = uicontrol('Style', 'radiobutton','Parent', gui_eeg_blc_dt.blc_bin_chan_option,...
            'String','Selected chan','callback',@Selected_bin_chan,'Value',0,'Enable',Enable_label,'FontSize',FonsizeDefault,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.Selected_bin_chan.KeyPressFcn= @EEG_blcorrdetrend_presskey;
        set(gui_eeg_blc_dt.blc_bin_chan_option, 'Sizes',[135  175]);
        gui_eeg_blc_dt.EEGTab_baseline_detrend{3} = gui_eeg_blc_dt.all_bin_chan.Value;
        %%Cancel and advanced
        gui_eeg_blc_dt.other_option = uiextras.HBox('Parent',gui_eeg_blc_dt.blc_dt,'Spacing',1,'BackgroundColor',ColorB_def);
        uiextras.Empty('Parent', gui_eeg_blc_dt.other_option,'BackgroundColor',ColorB_def);
        gui_eeg_blc_dt.Cancel = uicontrol('Parent',gui_eeg_blc_dt.other_option,'Style','pushbutton',...
            'String','Cancel','callback',@Cancel_blc_dt,'Enable',Enable_label,'FontSize',FonsizeDefault,'BackgroundColor',[1 1 1]);
        uiextras.Empty('Parent', gui_eeg_blc_dt.other_option);
        gui_eeg_blc_dt.apply = uicontrol('Style','pushbutton','Parent',gui_eeg_blc_dt.other_option,...
            'String','Apply','callback',@apply_blc_dt,'Enable',Enable_label,'FontSize',FonsizeDefault,'BackgroundColor',[1 1 1]);
        uiextras.Empty('Parent', gui_eeg_blc_dt.other_option);
        set(gui_eeg_blc_dt.other_option, 'Sizes',[15 105  30 105 15]);
        
        set(gui_eeg_blc_dt.blc_dt,'Sizes',[18 25 15 25 25 15 25 30]);
        
        estudioworkingmemory('EEGTab_baseline_detrend',0);
    end
%%*************************************************************************
%%*******************   Subfunctions   ************************************
%%*************************************************************************

%%--------------------------------setting for amplitude--------------------
    function  baseline_correction_EEG(~,~)
        if isempty(observe_EEGDAT.EEG)
            observe_EEGDAT.count_current_eeg=1;
            return;
        end
        %%first checking if the changes on the other panels have been applied
        [messgStr,eegpanelIndex] = f_check_eegtab_panelchanges();
        if ~isempty(messgStr) && eegpanelIndex~=15
            observe_EEGDAT.EEG_two_panels = observe_EEGDAT.EEG_two_panels+1;%%call the functions from the other panel
        end
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.apply.ForegroundColor = [1 1 1];
        EEG_basecorr_detrend_box.TitleColor= [ 0.5137    0.7569    0.9176];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [1 1 1];
        estudioworkingmemory('EEGTab_baseline_detrend',1);
        gui_eeg_blc_dt.blc.Value =1;
        gui_eeg_blc_dt.dt.Value = 0;
        gui_eeg_blc_dt.blc_dt_title.String = 'Baseline Period:';
    end

%%--------------------------Setting for phase------------------------------
    function detrend_EEG(~,~)
        if isempty(observe_EEGDAT.EEG)
            observe_EEGDAT.count_current_eeg=1;
            return;
        end
        %%first checking if the changes on the other panels have been applied
        [messgStr,eegpanelIndex] = f_check_eegtab_panelchanges();
        if ~isempty(messgStr) && eegpanelIndex~=15
            observe_EEGDAT.EEG_two_panels = observe_EEGDAT.EEG_two_panels+1;%%call the functions from the other panel
        end
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.apply.ForegroundColor = [1 1 1];
        EEG_basecorr_detrend_box.TitleColor= [ 0.5137    0.7569    0.9176];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [1 1 1];
        estudioworkingmemory('EEGTab_baseline_detrend',1);
        gui_eeg_blc_dt.dt.Value = 1;
        gui_eeg_blc_dt.blc.Value =0;
        gui_eeg_blc_dt.blc_dt_title.String = 'Calculate Trend During:';
    end

%%----------------Setting for "pre"----------------------------------------
    function pre_EEG(~,~)
        if isempty(observe_EEGDAT.EEG)
            observe_EEGDAT.count_current_eeg=1;
            return;
        end
        %%first checking if the changes on the other panels have been applied
        [messgStr,eegpanelIndex] = f_check_eegtab_panelchanges();
        if ~isempty(messgStr) && eegpanelIndex~=15
            observe_EEGDAT.EEG_two_panels = observe_EEGDAT.EEG_two_panels+1;%%call the functions from the other panel
        end
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.apply.ForegroundColor = [1 1 1];
        EEG_basecorr_detrend_box.TitleColor= [ 0.5137    0.7569    0.9176];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [1 1 1];
        estudioworkingmemory('EEGTab_baseline_detrend',1);
        gui_eeg_blc_dt.pre.Value=1;
        gui_eeg_blc_dt.post.Value=0;
        gui_eeg_blc_dt.whole.Value=0;
        gui_eeg_blc_dt.custom.Value=0;
        gui_eeg_blc_dt.custom_edit.Enable = 'off';
        if observe_EEGDAT.EEG.times(1)>=0
            CUstom_String = '';
        else
            CUstom_String = num2str([observe_EEGDAT.EEG.times(1),0]);
        end
        gui_eeg_blc_dt.custom_edit.String = CUstom_String;
    end


%%----------------Setting for "post"---------------------------------------
    function post_EEG(~,~)
        if isempty(observe_EEGDAT.EEG)
            observe_EEGDAT.count_current_eeg=1;
            return;
        end
        %%first checking if the changes on the other panels have been applied
        [messgStr,eegpanelIndex] = f_check_eegtab_panelchanges();
        if ~isempty(messgStr) && eegpanelIndex~=15
            observe_EEGDAT.EEG_two_panels = observe_EEGDAT.EEG_two_panels+1;%%call the functions from the other panel
        end
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.apply.ForegroundColor = [1 1 1];
        EEG_basecorr_detrend_box.TitleColor= [ 0.5137    0.7569    0.9176];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [1 1 1];
        estudioworkingmemory('EEGTab_baseline_detrend',1);
        gui_eeg_blc_dt.pre.Value=0;
        gui_eeg_blc_dt.post.Value=1;
        gui_eeg_blc_dt.whole.Value=0;
        gui_eeg_blc_dt.custom.Value=0;
        gui_eeg_blc_dt.custom_edit.Enable = 'off';
        if observe_EEGDAT.EEG.times(end)<=0
            CUstom_String = '';
        else
            CUstom_String = num2str([0 observe_EEGDAT.EEG.times(end)]);
        end
        gui_eeg_blc_dt.custom_edit.String = CUstom_String;
    end

%%----------------Setting for "whole"--------------------------------------
    function whole_EEG(~,~)
        if isempty(observe_EEGDAT.EEG)
            observe_EEGDAT.count_current_eeg=1;
            return;
        end
        %%first checking if the changes on the other panels have been applied
        [messgStr,eegpanelIndex] = f_check_eegtab_panelchanges();
        if ~isempty(messgStr) && eegpanelIndex~=15
            observe_EEGDAT.EEG_two_panels = observe_EEGDAT.EEG_two_panels+1;%%call the functions from the other panel
        end
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.apply.ForegroundColor = [1 1 1];
        EEG_basecorr_detrend_box.TitleColor= [ 0.5137    0.7569    0.9176];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [1 1 1];
        estudioworkingmemory('EEGTab_baseline_detrend',1);
        gui_eeg_blc_dt.pre.Value=0;
        gui_eeg_blc_dt.post.Value=0;
        gui_eeg_blc_dt.whole.Value=1;
        gui_eeg_blc_dt.custom.Value=0;
        gui_eeg_blc_dt.custom_edit.Enable = 'off';
        CUstom_String = num2str([observe_EEGDAT.EEG.times(1) observe_EEGDAT.EEG.times(end)]);
        gui_eeg_blc_dt.custom_edit.String = CUstom_String;
    end

%%----------------Setting for "custom"-------------------------------------
    function custom_EEG(~,~)
        if isempty(observe_EEGDAT.EEG)
            observe_EEGDAT.count_current_eeg=1;
            return;
        end
        %%first checking if the changes on the other panels have been applied
        [messgStr,eegpanelIndex] = f_check_eegtab_panelchanges();
        if ~isempty(messgStr) && eegpanelIndex~=15
            observe_EEGDAT.EEG_two_panels = observe_EEGDAT.EEG_two_panels+1;%%call the functions from the other panel
        end
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.apply.ForegroundColor = [1 1 1];
        EEG_basecorr_detrend_box.TitleColor= [ 0.5137    0.7569    0.9176];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [1 1 1];
        estudioworkingmemory('EEGTab_baseline_detrend',1);
        gui_eeg_blc_dt.pre.Value=0;
        gui_eeg_blc_dt.post.Value=0;
        gui_eeg_blc_dt.whole.Value=0;
        gui_eeg_blc_dt.custom.Value=1;
        gui_eeg_blc_dt.custom_edit.Enable = 'on';
    end

%%----------------input baseline period defined by user--------------------
    function custom_edit(Source,~)
        if isempty(observe_EEGDAT.EEG)
            observe_EEGDAT.count_current_eeg=1;
            return;
        end
        %%first checking if the changes on the other panels have been applied
        [messgStr,eegpanelIndex] = f_check_eegtab_panelchanges();
        if ~isempty(messgStr) && eegpanelIndex~=15
            observe_EEGDAT.EEG_two_panels = observe_EEGDAT.EEG_two_panels+1;%%call the functions from the other panel
        end
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.apply.ForegroundColor = [1 1 1];
        EEG_basecorr_detrend_box.TitleColor= [ 0.5137    0.7569    0.9176];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [1 1 1];
        estudioworkingmemory('EEGTab_baseline_detrend',1);
        
        lat_osci = str2num(Source.String);
        if isempty(lat_osci)
            msgboxText =  ['Baseline Correction & Linear Detrend (Epoched EEG) - Invalid input for "baseline range"'];
            titlNamerro = 'Warning for EEG Tab';
            estudio_warning(msgboxText,titlNamerro);
            return;
        end
        if numel(lat_osci) ==1
            msgboxText =  ['Baseline Correction & Linear Detrend (Epoched EEG) - Wrong baseline range. Please, enter two values'];
            titlNamerro = 'Warning for EEG Tab';
            estudio_warning(msgboxText,titlNamerro);
            return;
        end
        if lat_osci(1)>= lat_osci(2)
            msgboxText =  ['Baseline Correction & Linear Detrend (Epoched EEG) - The first value must be smaller than the second one'];
            titlNamerro = 'Warning for EEG Tab';
            estudio_warning(msgboxText,titlNamerro);
            return;
        end
        if lat_osci(2) > observe_EEGDAT.EEG.times(end)
            msgboxText =  ['Baseline Correction & Linear Detrend (Epoched EEG) - Second value must be smaller than',32,num2str(observe_EEGDAT.EEG.times(end))];
            titlNamerro = 'Warning for EEG Tab';
            estudio_warning(msgboxText,titlNamerro);
            return;
        end
        if lat_osci(1) < observe_EEGDAT.EEG.times(1)
            msgboxText =  ['Baseline Correction & Linear Detrend (Epoched EEG) - First value must be larger than',32,num2str(observe_EEGDAT.EEG.times(1))];
            titlNamerro = 'Warning for EEG Tab';
            estudio_warning(msgboxText,titlNamerro);
            return;
        end
        
    end

%%---------------------Setting for all chan and bin------------------------
    function All_bin_chan(~,~)
        if isempty(observe_EEGDAT.EEG)
            observe_EEGDAT.count_current_eeg=1;
            return;
        end
        %%first checking if the changes on the other panels have been applied
        [messgStr,eegpanelIndex] = f_check_eegtab_panelchanges();
        if ~isempty(messgStr) && eegpanelIndex~=15
            observe_EEGDAT.EEG_two_panels = observe_EEGDAT.EEG_two_panels+1;%%call the functions from the other panel
        end
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.apply.ForegroundColor = [1 1 1];
        EEG_basecorr_detrend_box.TitleColor= [ 0.5137    0.7569    0.9176];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [1 1 1];
        estudioworkingmemory('EEGTab_baseline_detrend',1);
        gui_eeg_blc_dt.all_bin_chan.Value = 1;
        gui_eeg_blc_dt.Selected_bin_chan.Value = 0;
    end

%%----------------Setting for selected bin and chan------------------------
    function Selected_bin_chan(~,~)
        if isempty(observe_EEGDAT.EEG)
            observe_EEGDAT.count_current_eeg=1;
            return;
        end
        %%first checking if the changes on the other panels have been applied
        [messgStr,eegpanelIndex] = f_check_eegtab_panelchanges();
        if ~isempty(messgStr) && eegpanelIndex~=15
            observe_EEGDAT.EEG_two_panels = observe_EEGDAT.EEG_two_panels+1;%%call the functions from the other panel
        end
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.apply.ForegroundColor = [1 1 1];
        EEG_basecorr_detrend_box.TitleColor= [ 0.5137    0.7569    0.9176];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [0.5137    0.7569    0.9176];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [1 1 1];
        estudioworkingmemory('EEGTab_baseline_detrend',1);
        gui_eeg_blc_dt.all_bin_chan.Value = 0;
        gui_eeg_blc_dt.Selected_bin_chan.Value = 1;
    end
%%--------------------------Setting for plot-------------------------------
    function apply_blc_dt(~,~)
        if isempty(observe_EEGDAT.EEG)
            observe_EEGDAT.count_current_eeg=1;
            return;
        end
        %%first checking if the changes on the other panels have been applied
        [messgStr,eegpanelIndex] = f_check_eegtab_panelchanges();
        if ~isempty(messgStr) && eegpanelIndex~=15
            observe_EEGDAT.EEG_two_panels = observe_EEGDAT.EEG_two_panels+1;%%call the functions from the other panel
        end
        EEGArray =  estudioworkingmemory('EEGArray');
        if isempty(EEGArray)
            EEGArray =  length(observe_EEGDAT.ALLEEG);
            observe_EEGDAT.EEG = observe_EEGDAT.ALLEEG(end);
            observe_EEGDAT.CURRENTSET = EEGArray;
            estudioworkingmemory('EEGArray',EEGArray);
        end
        try
            if gui_eeg_blc_dt.pre.Value==1
                BaselineMethod = 'pre';
            elseif  gui_eeg_blc_dt.post.Value==1
                BaselineMethod = 'post';
            elseif  gui_eeg_blc_dt.whole.Value==1
                BaselineMethod = 'all';
            elseif  gui_eeg_blc_dt.custom.Value ==1
                BaselineMethod = str2num(gui_eeg_blc_dt.custom_edit.String);
            end
        catch
            BaselineMethod = 'pre';
        end
        %%Check the baseline period defined by the custom.
        if gui_eeg_blc_dt.custom.Value ==1
            if isempty(BaselineMethod)
                msgboxText =  ['Baseline Correction & Linear Detrend (Epoched EEG) - Invalid input for baseline range; Please Cancel two values'];
                titlNamerro = 'Warning for EEG Tab';
                estudio_warning(msgboxText,titlNamerro);
                observe_EEGDAT.eeg_panel_message =2;
                return;
            end
            if numel(BaselineMethod) ==1
                msgboxText =  ['Baseline Correction & Linear Detrend (Epoched EEG) - Wrong baseline range. Please, enter two values'];
                titlNamerro = 'Warning for EEG Tab';
                estudio_warning(msgboxText,titlNamerro);
                observe_EEGDAT.eeg_panel_message =2;
                return;
            end
            if BaselineMethod(1)>= BaselineMethod(2)
                msgboxText =  ['Baseline Correction & Linear Detrend (Epoched EEG) - The first value must be smaller than the second one'];
                titlNamerro = 'Warning for EEG Tab';
                estudio_warning(msgboxText,titlNamerro);
                observe_EEGDAT.eeg_panel_message =2;
                return;
            end
            if roundn(BaselineMethod(2),-3) > roundn(observe_EEGDAT.EEG.times(end),-3)
                msgboxText =  ['Baseline Correction & Linear Detrend (Epoched EEG) - Second value must be smaller than',32,num2str(observe_EEGDAT.EEG.times(end))];
                titlNamerro = 'Warning for EEG Tab';
                estudio_warning(msgboxText,titlNamerro);
                observe_EEGDAT.eeg_panel_message =2;
                return;
            end
            if roundn(BaselineMethod(1),-3) < roundn(observe_EEGDAT.EEG.times(1),-3)
                msgboxText =  ['Baseline Correction & Linear Detrend (Epoched EEG) - First value must be larger than',32,num2str(observe_EEGDAT.EEG.times(1))];
                titlNamerro = 'Warning for EEG Tab';
                estudio_warning(msgboxText,titlNamerro);
                observe_EEGDAT.eeg_panel_message =2;
                return;
            end
        end
        %%Run the function based on the defined parameters
        %%--------------Loop start for removeing baseline for the selected EEGsets------------
        if gui_eeg_blc_dt.dt.Value ==1
            Suffix_str = '_detrend';
        else
            Suffix_str = '_baselinecorr';
        end
        
        %%%%-------------------Loop fpor baseline correction---------------
        estudioworkingmemory('f_EEG_proces_messg','Baseline Correction & Linear Detrend (Epoched EEG)');
        observe_EEGDAT.eeg_panel_message =1; %%Marking for the procedure has been started.
        
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 1 1 1];
        gui_eeg_blc_dt.apply.ForegroundColor = [0 0 0];
        EEG_basecorr_detrend_box.TitleColor= [0.05,0.25,0.50];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [1 1 1];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [0 0 0];
        estudioworkingmemory('EEGTab_baseline_detrend',0);
        
        gui_eeg_blc_dt.EEGTab_baseline_detrend{1} = gui_eeg_blc_dt.blc.Value;
        if gui_eeg_blc_dt.pre.Value==1
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2} = 1;
        elseif gui_eeg_blc_dt.post.Value==1
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2} = 2;
        elseif gui_eeg_blc_dt.whole.Value==1
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2} = 3;
        elseif gui_eeg_blc_dt.custom.Value==1
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2} = str2num(gui_eeg_blc_dt.custom_edit.String);
        else
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2} =1;
        end
        gui_eeg_blc_dt.EEGTab_baseline_detrend{3} = gui_eeg_blc_dt.all_bin_chan.Value;
        
        ALLEEG = observe_EEGDAT.ALLEEG;
        ChanArray = estudioworkingmemory('EEG_ChanArray');
        if isempty(ChanArray) || any(ChanArray(:)<=0) || any(ChanArray(:)>observe_EEGDAT.EEG.nbchan)
            ChanArray = [1:observe_EEGDAT.EEG.nbchan];
        end
        
        ALLEEG_out = [];
        for Numofeeg = 1:numel(EEGArray)
            EEG = ALLEEG(EEGArray(Numofeeg));
            if EEG.trials==1
                msgboxText =  'Baseline Correction & Linear Detrend (Epoched EEG) cannot work for continous EEG';
                titlNamerro = 'Warning for EEG Tab';
                estudio_warning(msgboxText,titlNamerro);
                observe_EEGDAT.eeg_panel_message =2;
                return;
            end
            if gui_eeg_blc_dt.all_bin_chan.Value == 1
                ChanArray = [1:EEG.nbchan];
            else
                if any(ChanArray(:)>observe_EEGDAT.EEG.nbchan)
                    ChanArray = [1:EEG.nbchan];
                end
            end
            
            if gui_eeg_blc_dt.dt.Value ==1
                [EEG LASTCOM] = pop_eeglindetrend( EEG, 'Baseline', BaselineMethod,'ChanArray',ChanArray, 'History','gui' );
            else
                [EEG LASTCOM]= pop_blceeg( EEG , 'Baseline', BaselineMethod,'ChanArray',ChanArray,...
                    'Saveas', 'off','History','gui');
            end
            if isempty(LASTCOM)
                observe_EEGDAT.eeg_panel_message =2;
                return;
            end
            EEG = eegh(LASTCOM, EEG);
            if Numofeeg==1
                eegh(LASTCOM);
            end
            [ALLEEG_out,~,~] = pop_newset(ALLEEG_out, EEG, length(ALLEEG_out), 'gui', 'off');
        end%%Loop end for the selected ERset
        
        Answer = f_EEG_save_multi_file(ALLEEG_out,1:numel(EEGArray),Suffix_str);
        if isempty(Answer)
            observe_EEGDAT.eeg_panel_message =2;
            return;
        end
        
        if ~isempty(Answer{1})
            ALLEEG_out = Answer{1};
            Save_file_label = Answer{2};
        end
        for Numofeeg = 1:numel(EEGArray)
            EEG = ALLEEG_out(Numofeeg);
            if Save_file_label
                [pathstr, file_name, ext] = fileparts(EEG.filename);
                EEG.filename = [file_name,'.set'];
                [EEG, LASTCOM] = pop_saveset(EEG,'filename', EEG.filename, 'filepath',EEG.filepath,'check','on');
                EEG = eegh(LASTCOM, EEG);
                if Numofeeg==1
                    eegh(LASTCOM);
                end
            else
                EEG.filename = '';
                EEG.saved = 'no';
                EEG.filepath = '';
            end
            [ALLEEG,~,~] = pop_newset(ALLEEG, EEG, length(ALLEEG), 'gui', 'off');
        end
        estudioworkingmemory('f_EEG_BLS_Detrend',{BaselineMethod,0,1});
        observe_EEGDAT.ALLEEG = ALLEEG;
        try
            Selected_EEG_afd =  [length(observe_EEGDAT.ALLEEG)-numel(EEGArray)+1:length(observe_EEGDAT.ALLEEG)];
            observe_EEGDAT.CURRENTSET = length(observe_EEGDAT.ALLEEG)-numel(EEGArray)+1;
        catch
            Selected_EEG_afd = length(observe_EEGDAT.ALLEEG);
            observe_EEGDAT.CURRENTSET = length(observe_EEGDAT.ALLEEG);
        end
        observe_EEGDAT.EEG = observe_EEGDAT.ALLEEG(observe_EEGDAT.CURRENTSET);
        estudioworkingmemory('EEGArray',Selected_EEG_afd);
        assignin('base','EEG',observe_EEGDAT.EEG);
        assignin('base','CURRENTSET',observe_EEGDAT.CURRENTSET);
        assignin('base','ALLEEG',observe_EEGDAT.ALLEEG);
        
        observe_EEGDAT.count_current_eeg=1;
        observe_EEGDAT.eeg_panel_message =2;
    end


%%-----------------Setting for save option---------------------------------
    function Cancel_blc_dt(~,~)
        if isempty(observe_EEGDAT.EEG)
            observe_EEGDAT.count_current_eeg=1;
            return;
        end
        %%first checking if the changes on the other panels have been applied
        [messgStr,eegpanelIndex] = f_check_eegtab_panelchanges();
        if ~isempty(messgStr) && eegpanelIndex~=15
            observe_EEGDAT.EEG_two_panels = observe_EEGDAT.EEG_two_panels+1;%%call the functions from the other panel
        end
        try
            methodtype =   gui_eeg_blc_dt.EEGTab_baseline_detrend{1};
        catch
            methodtype=1;
            gui_eeg_blc_dt.EEGTab_baseline_detrend{1}=1;
        end
        if isempty(methodtype) || numel(methodtype)~=1 || (methodtype~=0 && methodtype~=1)
            methodtype=1;
            gui_eeg_blc_dt.EEGTab_baseline_detrend{1}=1;
        end
        gui_eeg_blc_dt.blc.Value =methodtype;
        gui_eeg_blc_dt.dt.Value = ~methodtype;
        if gui_eeg_blc_dt.blc.Value==1
            gui_eeg_blc_dt.blc_dt_title.String = 'Baseline Period:';
        else
            gui_eeg_blc_dt.blc_dt_title.String = 'Calculate Trend During:';
        end
        %%baseline period
        try
            bsperiod =   gui_eeg_blc_dt.EEGTab_baseline_detrend{2};
        catch
            bsperiod=1;
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2}=1;
        end
        if isempty(bsperiod) || (numel(bsperiod)~=1 && numel(bsperiod)~=2)
            bsperiod=1;
            gui_eeg_blc_dt.EEGTab_baseline_detrend{2}=1;
        end
        
        if numel(bsperiod)==1
            if bsperiod~=1 && bsperiod~=2 && bsperiod~=3
                bsperiod=1;
                gui_eeg_blc_dt.EEGTab_baseline_detrend{2}=1;
            end
            if bsperiod==2
                gui_eeg_blc_dt.pre.Value=0;
                gui_eeg_blc_dt.post.Value=1;
                gui_eeg_blc_dt.whole.Value=0;
            elseif   bsperiod==3
                gui_eeg_blc_dt.pre.Value=0;
                gui_eeg_blc_dt.post.Value=0;
                gui_eeg_blc_dt.whole.Value=1;
            else
                gui_eeg_blc_dt.pre.Value=1;
                gui_eeg_blc_dt.post.Value=0;
                gui_eeg_blc_dt.whole.Value=0;
            end
            gui_eeg_blc_dt.custom.Value=0;
            gui_eeg_blc_dt.custom_edit.Enable = 'off';
            gui_eeg_blc_dt.custom_edit.String = '';
        elseif numel(bsperiod)==2
            gui_eeg_blc_dt.pre.Value=0;
            gui_eeg_blc_dt.post.Value=0;
            gui_eeg_blc_dt.whole.Value=0;
            gui_eeg_blc_dt.custom.Value=1;
            gui_eeg_blc_dt.custom_edit.Enable = 'on';
            if any(bsperiod> observe_EEGDAT.EEG.times(end)) || any(bsperiod< observe_EEGDAT.EEG.times(1))
                bsperiod = [];
                gui_eeg_blc_dt.EEGTab_baseline_detrend{2}=[];
            end
            gui_eeg_blc_dt.custom_edit.String = num2str(bsperiod);
        end
        
        %%bin & chan selection
        try
            all_bin_chan = gui_eeg_blc_dt.EEGTab_baseline_detrend{3};
        catch
            all_bin_chan=1;
        end
        if isempty(all_bin_chan) || numel(all_bin_chan)~=1 || (all_bin_chan~=0&& all_bin_chan~=1)
            gui_eeg_blc_dt.EEGTab_baseline_detrend{3}=1;
            all_bin_chan=1;
        end
        gui_eeg_blc_dt.all_bin_chan.Value = all_bin_chan;
        gui_eeg_blc_dt.Selected_bin_chan.Value = ~all_bin_chan;
        
        estudioworkingmemory('EEGTab_baseline_detrend',0);
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 1 1 1];
        gui_eeg_blc_dt.apply.ForegroundColor = [0 0 0];
        EEG_basecorr_detrend_box.TitleColor= [0.05,0.25,0.50];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [1 1 1];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [0 0 0];
    end


%%-------------------Setting for the whole panel of fitering based on ALLEEG and CURRENTEEG--------------
    function count_current_eeg_change(~,~)
        if observe_EEGDAT.count_current_eeg~=21
            return;
        end
        EEGUpdate = estudioworkingmemory('EEGUpdate');
        if isempty(EEGUpdate) || numel(EEGUpdate)~=1 || (EEGUpdate~=0 && EEGUpdate~=1)
            EEGUpdate = 0;  estudioworkingmemory('EEGUpdate',0);
        end
        if  isempty(observe_EEGDAT.EEG) || (~isempty(observe_EEGDAT.EEG)&& observe_EEGDAT.EEG.trials==1)|| EEGUpdate==1
            Enable_Label = 'off';
        else
            Enable_Label = 'on';
        end
        gui_eeg_blc_dt.blc.Enable = Enable_Label;
        gui_eeg_blc_dt.dt.Enable = Enable_Label;
        gui_eeg_blc_dt.apply.Enable = Enable_Label;
        gui_eeg_blc_dt.Cancel.Enable = Enable_Label;
        gui_eeg_blc_dt.pre.Enable= Enable_Label;
        gui_eeg_blc_dt.post.Enable= Enable_Label;
        gui_eeg_blc_dt.whole.Enable= Enable_Label;
        gui_eeg_blc_dt.custom.Enable= Enable_Label;
        gui_eeg_blc_dt.custom_edit.Enable = Enable_Label;
        gui_eeg_blc_dt.apply.Enable = Enable_Label;
        gui_eeg_blc_dt.Cancel.Enable = Enable_Label;
        gui_eeg_blc_dt.all_bin_chan.Enable = Enable_Label;
        gui_eeg_blc_dt.Selected_bin_chan.Enable = Enable_Label;
        if ~isempty(observe_EEGDAT.EEG) && observe_EEGDAT.EEG.trials ==1
            EEG_basecorr_detrend_box.TitleColor= [0.7500    0.7500    0.75000];
        else
            EEG_basecorr_detrend_box.TitleColor= [0.0500    0.2500    0.5000];
        end
        
        if gui_eeg_blc_dt.custom.Value==1
            gui_eeg_blc_dt.custom_edit.Enable = 'on';
        else
            gui_eeg_blc_dt.custom_edit.Enable = 'off';
        end
        if  isempty(observe_EEGDAT.EEG) || (~isempty(observe_EEGDAT.EEG)&& observe_EEGDAT.EEG.trials==1)|| EEGUpdate==1
            observe_EEGDAT.count_current_eeg=22;
            return;
        end
        EEGArray =  estudioworkingmemory('EEGArray');
        if isempty(EEGArray) || any(EEGArray> length(observe_EEGDAT.ALLEEG))
            EEGArray =  length(observe_EEGDAT.ALLEEG);
            observe_EEGDAT.EEG = observe_EEGDAT.ALLEEG(end);
            estudioworkingmemory('EEGArray',EEGArray);
            observe_EEGDAT.CURRENTSET = EEGArray;
        end
        
        if gui_eeg_blc_dt.custom.Value==1
            baseline = str2num(gui_eeg_blc_dt.custom_edit.String);
            if ~isempty(baseline)
                if any(baseline>observe_EEGDAT.EEG.times(end)) || any(baseline<observe_EEGDAT.EEG.times(1))
                    gui_eeg_blc_dt.custom_edit.String = '';
                end
            end
        end
        observe_EEGDAT.count_current_eeg=22;
    end

%%--------------press return to execute "Apply"----------------------------
    function EEG_blcorrdetrend_presskey(~,eventdata)
        keypress = eventdata.Key;
        ChangeFlag =  estudioworkingmemory('EEGTab_baseline_detrend');
        if ChangeFlag~=1
            return;
        end
        if strcmp (keypress, 'return') || strcmp (keypress , 'enter')
            apply_blc_dt();
            estudioworkingmemory('EEGTab_baseline_detrend',0);
            gui_eeg_blc_dt.apply.BackgroundColor =  [ 1 1 1];
            gui_eeg_blc_dt.apply.ForegroundColor = [0 0 0];
            EEG_basecorr_detrend_box.TitleColor= [0.05,0.25,0.50];%% the default is [0.0500    0.2500    0.5000]
            gui_eeg_blc_dt.Cancel.BackgroundColor =  [1 1 1];
            gui_eeg_blc_dt.Cancel.ForegroundColor = [0 0 0];
        else
            return;
        end
    end

%%--------------------------Reset with default parameters------------------
    function Reset_eeg_panel_change(~,~)
        if observe_EEGDAT.Reset_EEG_paras_panel~=17
            return;
        end
        estudioworkingmemory('EEGTab_baseline_detrend',0);
        gui_eeg_blc_dt.apply.BackgroundColor =  [ 1 1 1];
        gui_eeg_blc_dt.apply.ForegroundColor = [0 0 0];
        EEG_basecorr_detrend_box.TitleColor= [0.05,0.25,0.50];%% the default is [0.0500    0.2500    0.5000]
        gui_eeg_blc_dt.Cancel.BackgroundColor =  [1 1 1];
        gui_eeg_blc_dt.Cancel.ForegroundColor = [0 0 0];
        gui_eeg_blc_dt.blc.Value =1;
        gui_eeg_blc_dt.dt.Value = 0;
        gui_eeg_blc_dt.pre.Value=1;
        gui_eeg_blc_dt.post.Value=0;
        gui_eeg_blc_dt.whole.Value=0;
        gui_eeg_blc_dt.custom.Value=0;
        gui_eeg_blc_dt.custom_edit.Enable = 'off';
        gui_eeg_blc_dt.custom_edit.String = '';
        gui_eeg_blc_dt.all_bin_chan.Value = 1;
        gui_eeg_blc_dt.Selected_bin_chan.Value = 0;
        gui_eeg_blc_dt.blc_dt_title.String = 'Baseline Period:';
        observe_EEGDAT.Reset_EEG_paras_panel=18;
    end
end
%Progem end: