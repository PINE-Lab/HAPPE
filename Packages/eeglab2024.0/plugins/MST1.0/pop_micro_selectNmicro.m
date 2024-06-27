% [] = pop_micro_selectNmicro() - Select active number of microstates.
%
%  Select active number of microstates. By default this function will plot
%  a new figure with four measures of fit of the best segmentation for
%  each number of microstates.
%  The user can also directly define the active number of microstates
%  without plotting measures.
%  Note: Use the function MicroPlotFitmeas to plot measures of fit in an
%  existing figure.
%
% Usage:
%   >> EEG = pop_micro_selectNmicro ( EEG ); % pop up window
%   >> EEG = pop_micro_selectNmicro( EEG, 'key1', 'val1', 'key2', 'val2', ....)
%
%  Please cite this toolbox as:
%  Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (2018).
%  Microstate EEGlab toolbox: An introductionary guide. bioRxiv.
%
%  Inputs
%  EEG      - EEG-lab EEG structure (channels x samples (x epochs)) with
%             EEG.microstate field.
%
%  Optional input:
%  'plot_range'  - The range of numbers of microstates to be displayed (as 
%                  vector). If empty (default), measures for all microstate
%                  numbers will be plotted.
%  'Measures'    - Cell array of strings defining, which measures of fit to
%                  plot. Can also be a single string. Default is 'ALL',
%                  which plots all measures. Possible strings: 'CV', 'GEV', 
%                  'W', 'KL', 'KL_nrm' and 'ALL'.
%  'do_subplots' - If set, the created figure will contain a subplot for
%                  each measure of fit, instead of plotting all in one
%                  plot. 1 for subplots, 0 for all in on plot (default).
%  'Nmicro'      - If defined, this is made the active number of microstates
%                  used. If this is set, no figure will be created. Is
%                  empty by default. Can be an interger or empty.
%
% Outputs:
%   EEG    - Output dataset. This function saves output in the substruct
%            OUTEEG.microstate, which contains info general and algorithm-
%            specific settings as well as the results of segemntation.
%
% Authors:
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% April 2017.
%
% Copyright (C) 2017  Andreas Trier Poulsen, atpo@dtu.dk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function [EEG, com] = pop_micro_selectNmicro(EEG,varargin)
%% Error check and initialisation
com = '';

if nargin < 1
    help pop_micro_selectNmicro;
    return;
end

% check whether necessary microstate substructures exist.
if ~isfield(EEG,'microstate')
   error('No microstate data present. Run microstate segmentation first.') 
end
    
if nargin < 2
    % pop-up window in case no further input is given
    settings = input_popup();
    if strcmp(settings,'cancel')
        return
    end
else
    settings = check_settings(varargin);
end;

Nmicro = settings.Nmicro;


%% Check and ordering of selected measures of fit
if ischar(settings.Measures)
    if strcmp(settings.Measures,'ALL')
        Measures = {'CV', 'GEV', 'W', 'KL_nrm', 'KL'};
    else
        Measures = {settings.Measures};
    end
else
    Measures = settings.Measures;
end


% check if measures of fit are present
for m = 1:length(Measures)
    if ~isfield(EEG.microstate.Res, Measures{m})
        error('Could not find the measure %s in EEG.microstate.Res.',Measures{m})
    end
end


%% Pop-up for plotting measures of fit and selecting active number of microstates
if isempty(Nmicro)
    Nmicro = fitmeas_popup(EEG, Measures, settings);
    if strcmp(Nmicro,'cancel')
        return
    end
end


%% Selecting active number of microstates
Nmicro_ind = find(EEG.microstate.algorithm_settings.Nmicrostates == Nmicro);

if isempty(Nmicro_ind)
    error('The selected number of microstates, "%d" could not be found.', Nmicro)
else
    EEG.microstate.prototypes = EEG.microstate.Res.A_all{Nmicro_ind};
    EEG.microstate.labels = EEG.microstate.Res.L_all{Nmicro_ind};
    EEG.microstate.Res.K_act = Nmicro;
end


%% Define command string (only for Nmicro)
com = sprintf('%s = pop_micro_selectNmicro( %s, ''Nmicro'', %d );', ...
    inputname(1), inputname(1), Nmicro);

% com = settings_to_string(com,settings);
% com = [com ' );'];
end

% -------------- Pop-ups-------------- %
function settings = input_popup()
% Function for creating popup window to input algorithm settings
%

%% Create Inputs for popup
% plot_range
style.plot_range = 'edit';
line.plot_range = { {'Style' 'text' 'string' 'Range of microstates to plot:'}, ...
    {'Style' style.plot_range 'string' '' 'tag' 'plot_range'},... %end of first line
    {'Style' 'text' 'string' '(If empty; entire range will be plotted).'},...
    {} }; %end of second line
geo.plot_range = {[1 .3] [1 .3]};

% Measures title string
meas_str = 'Measures of fit to be plotted:';
line.meastitle = { {} {'Style' 'text' 'string' meas_str} };
geo.meastitle = {1 .3};

% CV
style.CV = 'checkbox';
CV_tipstr = 'Cross validation criterion.';
line.CV = { {'Style' style.CV 'value' 1 'string' 'CV' ...
    'tooltipstring' CV_tipstr 'tag' 'CV'} {} };
geo.CV = {[1 1]};

% GEV
style.GEV = 'checkbox';
GEV_tipstr = 'Global explained variance.';
line.GEV = { {'Style' style.GEV 'value' 1 'string' 'GEV' ...
    'tooltipstring' GEV_tipstr 'tag' 'GEV'} {} };
geo.GEV = {[1 1]};

% W
style.W = 'checkbox';
W_tipstr = 'Dispersion.';
line.W = { {'Style' style.W 'value' 1 'string' 'W' ...
    'tooltipstring' W_tipstr 'tag' 'W'} {} };
geo.W = {[1 1]};

% KL
style.KL = 'checkbox';
KL_tipstr = 'Krzanowski-Lai criterion.';
line.KL = { {'Style' style.KL 'value' 1 'string' 'KL' ...
    'tooltipstring' KL_tipstr 'tag' 'KL'} {} };
geo.KL = {[1 1]};

% KL_nrm
style.KL_nrm = 'checkbox';
KLnrm_tipstr = 'Normalised Krzanowski-Lai criterion.';
line.KL_nrm = { {'Style' style.KL_nrm 'value' 0 'string' 'KL_nrm' ...
    'tooltipstring' KLnrm_tipstr 'tag' 'KL_nrm'} {} };
geo.KL_nrm = {[1 1]};

% Do subplot?
style.do_subplots = 'checkbox';
% sub_tipstr = 'Plot measures on seperate subplots?';
line.do_subplots = { {'Style' style.do_subplots 'value' 1 'string'...
    'Plot measures on seperate subplots?' ...
    'tag' 'do_subplots'} {} };
%     'tooltipstring' sub_tipstr 'tag' 'do_subplots'} {} };
geo.do_subplots = {[1 1]};

% Nmicro
style.Nmicro = 'edit';
line.Nmicro = { {'Style' 'text' 'string' 'Select active number of microstates:'}, ...
    {'Style' style.Nmicro 'string' '' 'tag' 'Nmicro'},... %end of first line
    {'Style' 'text' 'string' '(If set; no figure will be created).'},...
    {} }; %end of second line
geo.Nmicro = {[1 .3] [1 .3]};

%% Order inputs for GUI
geometry = [geo.plot_range geo.meastitle geo.CV geo.GEV geo.W ...
    geo.KL geo.KL_nrm {1} geo.do_subplots {1} geo.Nmicro];
uilist = [line.plot_range line.meastitle line.CV line.GEV line.W ...
    line.KL line.KL_nrm {{}} line.do_subplots {{}} line.Nmicro];
% geometry = [geo.info geo.plot_range ]%geo.meastitle geo.CV geo.GEV geo.W];% ...
% %     geo.KL {1} geo.do_subplots geo.Nmicro];
% uilist = [line.info line.plot_range ]%line.meastitle line.CV line.GEV line.W];% ...
% %     line.KL {{}} line.do_subplots line.Nmicro];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_selectNmicro'');', 'Select active number of microstates -- pop_micro_selectNmicro()');


%% Interpret output from popup
if isstruct(pop_out)
    settings = struct;
    settings = interpret_popup(pop_out, settings, style);
    
    % Save selected measures of fit in Measures
    all_measures = {'CV', 'GEV', 'W', 'KL_nrm' 'KL'};
    settings.Measures = {};
    for m = 1:length(all_measures)
        if settings.(all_measures{m})
            settings.Measures{end+1} = all_measures{m};
        end
        settings = rmfield(settings, all_measures{m}); % removing field from settings
    end
else
    settings = 'cancel';
end

end

function Nmicro = fitmeas_popup(EEG, Measures, settings)
% Creates a pop-up with plot(s) for measures of fit, and a button to
% select number of microstates.
Nmeas = length(Measures);
Nmicro = 'cancel'; % if user presses cancel button.
if isempty(settings.plot_range)
    plot_range = EEG.microstate.algorithm_settings.Nmicrostates;
else
    plot_range = settings.plot_range;
end

% Find plot_idx and check that all elements in plot_range are in Nmicrostates
[~, plot_idx] = ismember(plot_range, EEG.microstate.algorithm_settings.Nmicrostates);
if sum(plot_idx==0)
   warning(['The following requested microstates numbers are present in '...
       'EEG.microstate.algorithm_settings.Nmicrostates. These will be skipped:'])
   disp(plot_range(plot_idx==0))
end


%% Create figure
h = figure('Units', 'normalized','position',[.2 .2 .6 .6], 'Visible','off');
if settings.do_subplots
    for m = 1:Nmeas
        subplot(Nmeas+1,1,m)
        MicroPlotFitmeas(EEG.microstate.Res, {Measures{m}}, plot_range, plot_idx)
    end
else
    subplot('position',[0.1 .3 .85 .65])
    MicroPlotFitmeas(EEG.microstate.Res, Measures, plot_range, plot_idx)
end


%% Nmicro input
Nmicro_txt = uicontrol('Style','text','Units', 'normalized',...
    'Position',[.2 .15 .3 .04], 'HorizontalAlignment', 'right',...
    'String','No. of microstates: ');
Nmicro_edit = uicontrol('Style','edit', 'Units', 'normalized',...
    'Position',[.51 .15 .05 .04], 'HorizontalAlignment', 'left',...
    'BackgroundColor', [1 1 1], 'String', ['  ' num2str(EEG.microstate.Res.K_act)]);


%% Buttons
callback_cancel = 'close gcbf';
btn_cancel = uicontrol('Style', 'pushbutton', 'String', 'Cancel',...
    'Units', 'normalized', 'Position', [.69 .1 .1 .04],...
    'Callback', callback_cancel);

callback_ok = 'set(gcbo, ''userdata'', 1);';
btn_ok = uicontrol('Style', 'pushbutton', 'String', 'OK',...
    'Units', 'normalized', 'Position', [.8 .1 .1 .04],...
    'tag', 'ok', 'Callback', callback_ok);

set(h,'Visible','on');
waitfor( findobj('parent', h, 'tag', 'ok'), 'userdata');

if ~(ishandle(h)), return; end % Check if figure still exist

Nmicro = eval(get(Nmicro_edit, 'string'));
close(h);


end

% -------------- helper functions -------------- %
function settings = interpret_popup(pop_out, settings, style, popmenu)
% Interpret output from pop_up window, "pop_out", and arrange it in
% "settings" struct. The fields in "style" should be the same as in "pop_out"
% (defined as tags in inputgui.m). The struct popmenu is optional and only
% needed if popmenus are used in the pop_up window.

names = fieldnames(style);
for i = 1:length(names)
    switch style.(names{i})
        case 'edit'
            if isempty(pop_out.(names{i})) % empty?
                settings.(names{i}) = [];
            else
                settings.(names{i}) = eval(pop_out.(names{i}));
            end
        case 'checkbox'
            settings.(names{i}) = pop_out.(names{i});
        case 'popupmenu'
            settings.(names{i}) = popmenu.(names{i}){pop_out.(names{i})};
    end
end

end

function settings = check_settings(vargs)
%% check settings
% Checks settings given as optional inputs for MicroPlot.
% Undefined inputs is set to default values.
varg_check = {   'Measures'  {'string' 'cell'}  []  'ALL' ;
    'plot_range' 'real' [] [];
    'do_subplots' 'integer' [0 1] 1;
    'Nmicro' 'integer' [] []};
settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end; % check for error
end

function com = settings_to_string(com,settings)
% Adds settings struct to existing com string in the form 'key1', 'val1',
% 'key2', 'val2' ... .
% Can handle structs, strings, vectors and scalars. I.e. not matrices.

names = fieldnames(settings);

for i = 1:length(names)
    if isstruct(settings.(names{i})) % struct?
        com = settings_to_string(com,settings.(names{i}));
    elseif isempty(settings.(names{i})) % empty?
        com = [ com sprintf(', ''%s'', []', names{i}) ];
    elseif ischar(settings.(names{i})) % string?
        com = [ com sprintf(', ''%s'', ''%s''', names{i}, settings.(names{i})) ];
    elseif length(settings.(names{i})) > 1 % vector?
        N_elements = length(settings.(names{i}));
        range = max(settings.(names{i})) - min(settings.(names{i})) + 1;
        if  N_elements == range % write vetor as 'min_value:max_value'
            com = [ com sprintf(', ''%s'', %g:%g', names{i}, ...
                min(settings.(names{i})), max(settings.(names{i}))) ];
        else % write vector with individual elements
            com = [ com sprintf(', ''%s'', [%g', names{i}, ...
                settings.(names{i})(1))];
            for n = 2:N_elements
                com = [ com sprintf(',%g', settings.(names{i})(n))];
            end
            com = [ com ']'];
        end
    else % scalar
        com = [ com sprintf(', ''%s'', %g', names{i}, settings.(names{i})) ];
    end
end

end
