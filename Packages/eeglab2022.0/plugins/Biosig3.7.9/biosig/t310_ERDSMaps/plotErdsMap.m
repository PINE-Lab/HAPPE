function plotErdsMap(r, varargin)
% Displays time-frequency (ERDS) maps.
%
% This function plots ERDS maps as calculated by calcErdsMap.m.
%
% Usage:
%   plotErdsMap(r);
%
% Input parameters:
%   r ... Input structure calculated with calcErdsMap.
%
% Optional input parameters (variable argument list):
%   't_range' ... Time range to plot <1x2>. Specify start and end points within a 
%                 trial (in s) to plot only a specific time range.
%                 Default: The whole time range is plotted.

% Copyright by Clemens Brunner
% $Revision: 0.8 $ $Date: 03/10/2009 14:46:00 $
% E-Mail: clemens.brunner@tugraz.at

% Revision history:
%   0.8: Add 'range' option to plot only a specific segment of time.
%   0.7: Another complete rewrite to adapt to new versions of the map function.
%   0.6: Complete rewrite of the function using axes commands that makes the
%        whole thing much more customizable.

% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 59 Temple
% Place - Suite 330, Boston, MA  02111-1307, USA.

if (nargin < 1)
    error('No input ERDS map specified.');
end;

% Default parameters, can be overwritten by varargin
t_range = [];  % Plot the whole time range

% Overwriting default values with optional input parameters
if ~isempty(varargin)  % Are there optional parameters available?
   k = 1;
   while k <= length(varargin)
       if strcmp(varargin{k}, 't_range')
           t_range = varargin{k + 1};
           k = k + 2;
       else  % Ignore unknown parameters
            k = k + 2;
       end;
   end;
end;

if ~isfield(r, 'refmethod')
    r.refmethod = 'classic';
end;

% Does the range interval lie inside the calculated time segment?
if ~isempty(t_range)
    if numel(t_range) ~= 2
        error('Argument t_range must be a <1x2> vector containing the start and end point (in s).');
    end;
    if t_range(1) >= t_range(2)
        error('First element of t_range must be less than the second element.');
    end;
    if t_range(1) < r.t_plot(1) || t_range(2) > r.t_plot(end)
        error('Argument t_range must lie inside calculated time segment.');
    end;
    
    % Cut out time segment to plot
    [temp, pos] = find(r.t_plot >= t_range(1) & r.t_plot <= t_range(2));
    r.t_plot = r.t_plot(pos);
    for k = 1:length(r.ERDS)
        r.ERDS{k}.erds = r.ERDS{k}.erds(pos,:);
    end;
    if ~strcmp(r.sig, 'none')
        for k = 1:length(r.ERDS)
            r.ERDS{k}.cl = r.ERDS{k}.cl(pos,:);
            r.ERDS{k}.cu = r.ERDS{k}.cu(pos,:);
        end;
    end;
    
end;

border = 0.1;  % Border around figure
border_plots = 0.01;  % Border around each plot
plot_area = 1 - 2 * border;

% Topographic layout
if isfield(r, 'montage')
    plot_index = find(r.montage' == 1);
    n_rows = size(r.montage, 1);
    n_cols = size(r.montage, 2);
else  % create default layout
    plot_index = 1:length(r.ERDS);
    n_cols = ceil(sqrt(length(r.ERDS)));
    if (length(r.ERDS) > 2)
        n_rows = n_cols;
    else
        n_rows = 1;
    end;
end;

i_width = plot_area / n_cols;  % Width of one subplot
i_height = plot_area / n_rows;  % Height of one subplot
font_size = 1/32;  % Default normalized axes font size

f = figure;
set(f, 'PaperOrientation', 'landscape');
if ~exist('OCTAVE_VERSION','builtin')
	set(f, 'PaperType', 'A4');
end; 
set(f, 'PaperUnits', 'centimeters');
set(f, 'PaperPosition', [1, 1, 27.7, 19]);
set(f, 'Color', [1 1 1]);
set(f, 'DefaultAxesFontUnits', 'normalized');
set(f, 'DefaultAxesFontSize', font_size);

% If significance information exists, plot only significant data
if ~strcmp(r.sig, 'none')
    for chn = 1:length(r.ERDS)
        sig_matrix = (r.ERDS{chn}.cl > 0 & r.ERDS{chn}.cu > 0) | ...
                     (r.ERDS{chn}.cl < 0 & r.ERDS{chn}.cu < 0);
        r.ERDS{chn}.erds = sig_matrix .* r.ERDS{chn}.erds;
    end;
end;

% Invert color map so that ERS is blue and ERD is red
load erdscolormap;
colormap(erdcolormap);

counter_total = 1;  % Iterates through all rows and columns
counter_plots = 1;  % Iterates through all subplots

a = cell(1, length(plot_index));  % Contains the axes of the ERDS subplots

for i_rows = 1:n_rows
    for i_cols = 1:n_cols
        if sum(counter_total == plot_index) == 1
            a{counter_plots} = axes('position', [border + border_plots + i_width * (i_cols - 1), plot_area + border + border_plots - i_height * i_rows, i_width - border_plots, i_height - border_plots]);
            set(f, 'CurrentAxes', a{counter_plots});

            if strcmp(r.refmethod, 'absolute')
                imagesc(r.t_plot, r.f_plot, r.ERDS{counter_plots}.erds');
            else
                imagesc(r.t_plot, r.f_plot, r.ERDS{counter_plots}.erds', [-1, 1.5]);
            end;

            %set(gca, 'Tag', num2str(counter_plots));

            set(a{counter_plots}, 'ydir', 'normal');
            if sum(counter_total + n_cols == plot_index) == 1  % Draw x-labels only if there is no subplot below
                set(a{counter_plots}, 'XTickLabel', '');
            end;
            if sum(counter_total - 1 == plot_index) == 1 && i_cols ~= 1  % Draw y-labels only if there is no subplot to the left
                set(a{counter_plots}, 'YTickLabel', '');
            end;

            % Draw lines for reference interval and cue
            v = axis;
            line([r.ref(1), r.ref(1)], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k');
            line([r.ref(2), r.ref(2)], [v(3), v(4)], 'LineStyle', ':', 'Color', 'k');
            if isfield(r, 'cue')
                line([r.cue, r.cue], [v(3), v(4)], 'Color', 'k');
            end;
            
            counter_plots = counter_plots + 1;
        end;
        counter_total = counter_total + 1;
    end;
end;

% Heading
if isfield(r, 'heading')  % Recording date
    axes('position', [border + border_plots, 1 - 3/4* border, 1 - 2 * (border + border_plots), border], 'visible', 'off');
    text(0.5, 0, r.heading, 'FontUnits', 'normalized', 'FontSize', 1/4, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'Bottom', 'Interpreter', 'none', 'FontWeight', 'bold');
end;
%line([0, 1], [0, 0], 'Color', 'k');

% Additional text
axes('position', [border + border_plots, 3/4*border, 1 - 2 * (border + border_plots), 3/4*border], 'visible', 'off');
temp_str{1} = ['{\bfERDS maps 0.7 (', upper(r.ERDS_method), ').}{\rm Calculated on ', r.date_calc, '.}'];

if length(r.classes) > 1
    classes_str = '[';
    for k = 1:length(r.classes)
        classes_str = [classes_str, num2str(r.classes(k))];
        if k < length(r.classes)
            classes_str = [classes_str, ', '];
        end;
    end;
    classes_str = [classes_str, ']'];
else
    classes_str = num2str(r.classes);
end;
t_str = ['[', num2str(r.t(1)), ', ', num2str(r.t(2)), ', ', num2str(r.t(3)), ']s'];

f_borders_str = '[';
for k = 1:length(r.f_borders)
    f_borders_str = [f_borders_str, num2str(r.f_borders(k))];
    if k < length(r.f_borders)
        f_borders_str = [f_borders_str, ', '];
    end;
end;
f_borders_str = [f_borders_str, ']Hz'];

if length(r.f_bandwidths) > 1
    f_bandwidths_str = '[';
    for k = 1:length(r.f_bandwidths)
        f_bandwidths_str = [f_bandwidths_str, num2str(r.f_bandwidths(k))];
        if k < length(r.f_bandwidths)
            f_bandwidths_str = [f_bandwidths_str, ', '];
        end;
    end;
    f_bandwidths_str = [f_bandwidths_str, ']'];
else
    f_bandwidths_str = num2str(r.f_bandwidths);    
end;
f_bandwidths_str = [f_bandwidths_str, 'Hz'];

if length(r.f_steps) > 1
    f_steps_str = '[';
    for k = 1:length(r.f_steps)
        f_steps_str = [f_steps_str, num2str(r.f_steps(k))];
        if k < length(r.f_steps)
            f_steps_str = [f_steps_str, ', '];
        end;
    end;
    f_steps_str = [f_steps_str, ']'];
else
    f_steps_str = num2str(r.f_steps);    
end;
f_steps_str = [f_steps_str, 'Hz'];

ref_str = ['[', num2str(r.ref(1)), ', ', num2str(r.ref(2)), ']s'];

temp_str{2} = ['Trials: ', num2str(r.n_trials), ', classes: ', classes_str, ', fs: ', num2str(r.fs), 'Hz, time: ', t_str, ', ref: ', ref_str];
temp_str{3} = ['f borders: ', f_borders_str, ', f bandwidths: ', f_bandwidths_str, ', f steps: ', f_steps_str, ', '];
if strcmp(r.sig, 'boot')
    temp_str{3} = [temp_str{3}, 'Bootstrap significance test (\alpha = ', num2str(r.alpha), ')'];
elseif strcmp(r.sig, 'boxcox')
    temp_str{3} = [temp_str{3}, 'Box-Cox significance test (\alpha = ', num2str(r.alpha), ', \lambda = ', num2str(r.lambda), ')'];
else
    temp_str{3} = [temp_str{3}, 'no significance test.'];
end;
if exist('OCTAVE_VERSION','builtin')
	text([0,0,0], [0,1,2], temp_str, 'FontUnits', 'normalized', 'FontSize', 1/6, 'VerticalAlignment', 'Top', 'Interpreter', 'Tex');
else
	text(0, 0, temp_str, 'FontUnits', 'normalized', 'FontSize', 1/6, 'VerticalAlignment', 'Top', 'Interpreter', 'Tex');
end
%line([0, 1], [0, 0], 'Color', 'k');

%hndl = findobj('parent', gcf, 'type', 'axes');
%for a = 1:length(hndl)
%    set(findobj('parent', hndl(a)), 'ButtonDownFcn', 'zoomMap(r)');
%end;