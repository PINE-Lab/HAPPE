% pop_micro_stats() - calculates microstate statistics.
%
% Usage:
%   >> EEG = pop_micro_stats ( EEG ); % pop up window
%   >> EEG = pop_micro_stats ( EEG, 'key1', 'val1', 'key2', 'val2' ... )
%
% Please cite this toolbox as:
% Poulsen, A. T., Pedroni, A., Langer, N., &  Hansen, L. K. (2018).
% Microstate EEGlab toolbox: An introductionary guide. bioRxiv.
%
% Inputs:
%   EEG     - EEG-lab EEG structure (channels x samples (x epochs)) with
%             .microstate.fit.bestLabel (created by MicroFit.m)
%
% Optional inputs:
%  'label_type' - Use microstate labels obtained from: 'segmentation'
%                 (default) or 'backfit'.
%  'epoch'      - Timewindow of analysis (vector of timeframe indices). If
%                 empty, all time samples will be used (default).
%  'polarity'   - Account for polarity when calculating the global map
%                 dissimilarity. Typically off for spontaneous EEG and on
%                 for ERP data (default = 0).
%
% Outputs:
%  EEG.microstate.stats      - Structure of microstate parameters per trial.
%  EEG.microstate.stats.avgs - Structure of microstate parameters mean /
%                              standard deviation over trials (only for
%                              epoched data).
%
% Authors:
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, DTU Compute, Cognitive systems.
%
% Andreas Pedroni, andreas.pedroni@uzh.ch
% University of Z�rich, Psychologisches Institut, Methoden der
% Plastizit�tsforschung. 
%
% September 2017.
%
% See also: eeglab

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

function [EEG, com] = pop_micro_stats(EEG, varargin)
%% Error check and initialisation
if nargin < 1
    help pop_micro_stats;
    return;
end

% check whether necessary microstate substructures exist.
if ~isfield(EEG,'microstate')
   error('No microstate data present. Run microstate segmentation first.') 
end

com = '';


%% pop-up window in case no further input is given
if nargin < 2
    settings = stats_popup();
    if strcmp(settings,'cancel')
        return
    end
else
    settings = check_settings(varargin,EEG);
end


%% Select data
if strcmp(settings.label_type,'segmentation')
    if ischar(EEG.microstate.data)
        data = EEG.data;
    else
        data = EEG.microstate.data;
    end
    labels = EEG.microstate.labels;
else
    data = EEG.data;
    labels = EEG.microstate.fit.labels;
end
    
% Prototypes
A = EEG.microstate.prototypes;


%% Subtract data of interest (epoch)
if isempty(settings.epoch)
    epoch = 1:size(EEG.data,2); % all samples
else
    epoch = settings.epoch;
end
data = data(:,epoch,:);
labels = labels(:,epoch);


%% Run stats
Mstats = MicroStats(data,A,labels,settings.polarity,EEG.srate);


%% Save stats
if isfield(EEG.microstate,'stats')
   warning('Overwriting existing stats struct in EEG.microstate') 
end
EEG.microstate.stats = Mstats;
EEG.microstate.stats.label_type = settings.label_type;


%% Define command string
if isempty(settings.epoch)
    % removing epoch field if empty
    settings = rmfield(settings,'epoch');
end
com = sprintf('%s = pop_micro_stats( %s', inputname(1), inputname(1));
com = settings_to_string(com,settings);
com = [com ' );'];

end

% ------------------------------ Pop-up ---------------------------------%
function settings = stats_popup()
% Function for creating popup window to input settings
%

%% Create Inputs for popup
% Select labels to smooth
style.label_type = 'popupmenu';
dropdown_label = {'Microstate segmentation' 'Backfitting prototypes to EEG'}; % For dropdown menu
popmenu.label_type = {'segmentation' 'backfit'}; % Corresponding calls for pop-function
label_str = dropdown_label{1}; %string for popupmenu
for l = 2:length(dropdown_label); label_str = [label_str '|' dropdown_label{l}]; end
line.label_type = { {'Style' 'text' 'string' 'Use labels obtained from:'}, ...
    {'Style' style.label_type 'string' label_str 'tag' 'label_type' 'value' 1} };
geo.label_type = {[1 1]};

% Epoch
style.epoch = 'edit';
epoch_tipstr = 'Vector of timeframes. Leave empty to use entire range.';
line.epoch = { {'Style' 'text' 'string' 'Timewindow of analysis (vector of timeframe indices)', ...
    'tooltipstring' epoch_tipstr}, ...
    {'Style' style.epoch 'string' '' 'tag' 'epoch'} };
geo.epoch = {[1 .2]};

% Polarity
style.polarity = 'checkbox';
pol_tipstr = 'Spontaneous EEG typically ignore polarity (off). Typically on for ERP data.';
line.polarity = { {'Style' 'text' 'string' 'Account for polarity when fitting.' ...
    'tooltipstring' pol_tipstr}, ...
    {'Style' style.polarity 'value' 0 'tag' 'polarity'} };
geo.polarity = {[1 .2]};


%% Order inputs for GUI
geometry = [geo.label_type  geo.epoch geo.polarity];
uilist = [line.label_type line.epoch line.polarity];


%% Create Popup
[~,~,~,pop_out] = inputgui( geometry, uilist, ...
    'pophelp(''pop_micro_stats'');', 'calculates microstate stats -- pop_micro_stats()');
 

%% Interpret output from popup
if isstruct(pop_out)
    settings = struct;
    settings = interpret_popup(pop_out, settings, style, popmenu);
else
    settings = 'cancel';
end
end
% ----------------------------------------------------------------------- %

% -------------------------- Helper functions --------------------------- %
function settings = check_settings(vargs, EEG)
%% Check settings
% Checks settings given as optional inputs for MicroStats.
% Undefined inputs is set to default values.
varg_check = { 'label_type'  'string'    []         'segmentation';
    'epoch'  'integer'    []         1:size(EEG.data,2);
    'polarity'  'integer'    []         0};
settings = finputcheck( vargs, varg_check);
if ischar(settings), error(settings); end % check for error
end

function settings = interpret_popup(pop_out, settings, style, popmenu)
%%
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

function com = settings_to_string(com,settings)
%%
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