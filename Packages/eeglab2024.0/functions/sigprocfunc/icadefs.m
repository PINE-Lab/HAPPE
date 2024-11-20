% ICADEFS - function to read in a set of EEGLAB system-wide (i.e. lab-wide)
%             or working directory-wide constants and preferences. Change the 
%             way these are defined in the master icadefs.m file (usually
%             in dir eeglab/functions/sigprocfunc) or make a custom copy of 
%             the icadefs.m file in a project directory. Then, calling functions 
%             that call icadefs from an EEGLAB session in that working directory 
%             will read the local copy, which may set preferences different from 
%             the system-wide copy.
%
% Author: Arnaud Delorme, Scott Makeig, SCCN/INC/UCSD, La Jolla, 05-20-97 

% Copyright (C) 05-20-97 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This file is part of EEGLAB, see http://www.eeglab.org
% for the documentation and details.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

% ----------------------------------------------------------------------
% ------ EEGLAB DEFINITION - YOU MAY CHANGE THE TEXT BELOW -------------
% ----------------------------------------------------------------------



EEGOPTION_PATH = ''; % if empty, the home folder of the current user is used
                     % Note that this may create problems under Windows
                     % when unicode characters are part of the user name
                     % In this case, enter the path name manually here.

YDIR  = 1;                  % positive potential up = 1; negative up = -1 
                            % for most ERP plots

HZDIR = 'up';               % ascending freqs = 'up'; descending = 'down' 
                            % (e.g., timef/newtimef frequency direction)
                            
% Checking MATLAB version
clear version
tmpvers = version;
indp = find(tmpvers == '.');
if str2num(tmpvers(indp(1)+1)) >= 1, tmpvers = [ tmpvers(1:indp(1)) '0' tmpvers(indp(1)+1:end) ]; end
indp = find(tmpvers == '.');
VERS = str2num(tmpvers(1:indp(2)-1));                            

% font size
tmpComputer   = computer;
tmpScreenSize = get(0, 'ScreenSize');

% Graph Definitions
DEFAULT_COLORMAP = 'jet';

if VERS < 8.04
    PLOT_LINEWIDTH   = 2;
    PLOT_LINEWIDTH_S = 1;
    
    % AXES FONTSIZE
    AXES_FONTSIZE   = 10;                % Axis labels and legend font size
    AXES_FONTSIZE_S = AXES_FONTSIZE - 2; % Axis labels and legend font size Small
    AXES_FONTSIZE_L = 16;                % Axis labels and legend font size Large
    
    % GUI FONTSIZE
    GUI_FONTSIZE    = 10;               % graphic interface font size
    GUI_FONTSIZE_S  = GUI_FONTSIZE - 2; % graphic interface font size Small
    GUI_FONTSIZE_L  = GUI_FONTSIZE + 2; % graphic interface font size Large
   
    % TEXT FONTSIZE
    TEXT_FONTSIZE = 10;                  % Miscellaneous font sizes
    TEXT_FONTSIZE_S = TEXT_FONTSIZE - 2; % Miscellaneous font sizes Small
    TEXT_FONTSIZE_L = TEXT_FONTSIZE + 2; % Miscellaneous font sizes Large
    
elseif VERS >= 8.04
  
    if strcmpi(tmpComputer(1:3), 'MAC')
      
        PLOT_LINEWIDTH   = 1;
        PLOT_LINEWIDTH_S = 0.5;
      
        %scale up fontsizes on higher resolution mac screens
        retinaDisplay = false;
        if tmpScreenSize(3) >= 1920 % bump fontsize only for the highest retina res settings
            retinaDisplay = true; %comment this out if you don't want fontsizes increased at high display resolutions
            %disp('Mac OSX retina display detected. If this is not desired comment out line 83 of icadefs.m');
        end
        
        % AXES FONTSIZE
        if retinaDisplay
          AXES_FONTSIZE   = 12;                 % Axis labels and legend font size
        else
          AXES_FONTSIZE   = 9;                 % Axis labels and legend font size
        end
        AXES_FONTSIZE_S = AXES_FONTSIZE - 2; % Axis labels and legend font size Small
        AXES_FONTSIZE_L = 12.5;              % Axis labels and legend font size Large
        
        % GUI FONTSIZE
        if retinaDisplay
          GUI_FONTSIZE    = 14;                % graphic interface font size
        else
          GUI_FONTSIZE    = 12;                % graphic interface font size
        end
        GUI_FONTSIZE_S  = GUI_FONTSIZE - 2; % graphic interface font size Small
        GUI_FONTSIZE_L  = GUI_FONTSIZE + 2; % graphic interface font size Large
        
        % TEXT FONTSIZE
        if retinaDisplay
          TEXT_FONTSIZE   = 14;                 % Miscellaneous font sizes
        else
          TEXT_FONTSIZE   = 12;                 % Miscellaneous font sizes
        end
        TEXT_FONTSIZE_S = TEXT_FONTSIZE - 2; % Miscellaneous font sizes Small
        TEXT_FONTSIZE_L = TEXT_FONTSIZE + 4; % Miscellaneous font sizes Large
    else
        PLOT_LINEWIDTH   = 1;
        PLOT_LINEWIDTH_S = 0.5;
        
        % AXES FONTSIZE
        AXES_FONTSIZE   = 9;                 % Axis labels and legend font size
        AXES_FONTSIZE_S = AXES_FONTSIZE - 2; % Axis labels and legend font size Small
        AXES_FONTSIZE_L = 12.5;              % Axis labels and legend font size Large
        
        % GUI FONTSIZE
        GUI_FONTSIZE    = 10;                % graphic interface font size
        GUI_FONTSIZE_S  = GUI_FONTSIZE - 2; % graphic interface font size Small
        GUI_FONTSIZE_L  = GUI_FONTSIZE + 2; % graphic interface font size Large
        
        % TEXT FONTSIZE
        TEXT_FONTSIZE   = 10;                 % Miscellaneous font sizes
        TEXT_FONTSIZE_S = TEXT_FONTSIZE - 2; % Miscellaneous font sizes Small
        TEXT_FONTSIZE_L = TEXT_FONTSIZE + 4; % Miscellaneous font sizes Large
    end
end

clear retinaDisplay tmpScreenSize tmpComputer tmpvers indp;

% the eeg_options.m file also contains additional options

% ----------------------------------------------------------------------
% ------------------------ END OF DEFINITIONS --------------------------
% ----------------------------------------------------------------------

% INSERT location of ica executable (UNIX ONLY) for binica.m below
if ~isdeployed
    % ICA binary file in functions/supportfiles
    ICABINARY = 'ica_linux'; 
    tmpComputer = computer;
    if strcmpi(tmpComputer(1:3), 'MAC')
        ICABINARY = 'ica_osx';
    elseif strcmpi(tmpComputer(1:2), 'PC')
        ICABINARY = 'binica.exe';
    end
    clear tmpComputer
else
    ICABINARY = 'ica_linux';
end

try
    set(0,'defaultaxesfontsize',AXES_FONTSIZE);
    set(0,'defaulttextfontsize',TEXT_FONTSIZE);
    set(0,'DefaultUicontrolFontSize',GUI_FONTSIZE);
catch
    % most likely Octave here
end

TUTORIAL_URL = 'http://sccn.ucsd.edu/wiki/EEGLAB'; % online version
DEFAULT_SRATE = 256.0175;      % default local sampling rate (rarely used)
DEFAULT_TIMLIM = [-1000 2000]; % default local epoch limits (ms)

% Set EEGLAB figure and GUI colors
% --------------------------------
lowscreendepth = 0;
if ~exist('OCTAVE_VERSION', 'builtin')
    if get(0, 'screendepth') <=8 % if mono or 8-bit color
	lowscreendepth = 1; 
    end
end
if lowscreendepth    
    %fprintf('icadefs(): Setting display parameters for mono or 8-bit color\n');
    BACKCOLOR           = [1 1 1];    % Background figure color 
    BACKEEGLABCOLOR     = [1 1 1];    % EEGLAB main window background
    GUIBUTTONCOLOR      = [1 1 1];    % Buttons colors in figures
    GUIPOPBUTTONCOLOR   = [1 1 1];    % Buttons colors in GUI windows
    GUIBACKCOLOR        = [1 1 1];    % GUI background color
    GUITEXTCOLOR        = [0 0 0];      % GUI foreground color for text    
    PLUGINMENUCOLOR     = [.5 0 .5];  % plugin menu color

else % if full color screen
    BACKCOLOR           = [.93 .96 1];    % EEGLAB Background figure color 
    BACKEEGLABCOLOR     = [.66 .76 1];    % EEGLAB main window background
    GUIBUTTONCOLOR      = BACKEEGLABCOLOR;% Buttons colors in figures
    GUIPOPBUTTONCOLOR   = BACKCOLOR;      % Buttons colors in GUI windows
    GUIBACKCOLOR        = BACKEEGLABCOLOR;% EEGLAB GUI background color <---------
    GUITEXTCOLOR        = [0 0 0.4];      % GUI foreground color for text
    PLUGINMENUCOLOR     = [.5 0 .5];      % plugin menu color
end


% THE FOLLOWING PARAMETERS WILL BE DEPRECATED IN LATER VERSIONS
% -------------------------------------------------------------

SHRINKWARNING = 1;          % Warn user about the shrink factor in topoplot() (1/0)

MAXENVPLOTCHANS   = 264;  % maximum number of channels to plot in envproj.m
MAXPLOTDATACHANS  = 264;  % maximum number of channels to plot in dataplot.m
MAXPLOTDATAEPOCHS = 264;  % maximum number of epochs to plot in dataplot.m
MAXEEGPLOTCHANS   = 264;  % maximum number of channels to plot in eegplot.m
MAXTOPOPLOTCHANS  = 264;  % maximum number of channels to plot in topoplot.m

DEFAULT_ELOC  = 'chan.locs'; % default electrode location file for topoplot.m
DEFAULT_EPOCH = 10;       % default epoch width to plot in eegplot(s) (in sec)

SC  =  ['binica.sc'];           % Master .sc script file for binica.m
                                % MATLAB will use first such file found
                                % in its path of script directories.
                                % Copy to pwd to alter ICA defaults
