function plotAveVar(r, varargin)
% Displays the mean and standard deviation of each channel.
%
% This function displays the mean and standard deviation of each channel in a
% topographical layout as calculated by calcAveVar.m.
%
% Usage:
%   plotAveVar(r);
%
% Input parameters:
%   r ... Input structure calculated with calcAveVar.
%
% Optional input parameters (variable argument list):
%   't_range' ... Time range to plot <1x2>. Specify start and end points within a 
%                 trial (in s) to plot only a specific time range.
%                 Default: The whole time range is plotted.


% Copyright by Clemens Brunner
% $Revision: 0.5 $ $Date: 03/12/2009 14:43:05 $
% E-Mail: clemens.brunner@tugraz.at

% Revision history:
%   0.5: Add 'range' option to plot only a specific segment of time.
%   0.4: Another complete rewrite for the new toolbox.
%   0.3: Complete rewrite of the function using axes commands that makes the
%        whole thing much more customizable

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
    error('No input average/variance map specified.');
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
    for k = 1:length(r.AveVar)
       r.AveVar{k}.average = r.AveVar{k}.average(pos);
       r.AveVar{k}.stdev = r.AveVar{k}.stdev(pos);
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
    plot_index = 1:length(r.AveVar);
    n_cols = ceil(sqrt(length(r.AveVar)));
    if (length(r.AveVar) > 2)
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
set(f, 'PaperType', 'A4');
set(f, 'PaperUnits', 'centimeters');
set(f, 'PaperPosition', [1, 1, 27.7, 19]);
set(f, 'Color', [1 1 1]);
set(f, 'DefaultAxesFontUnits', 'normalized');
set(f, 'DefaultAxesFontSize', font_size);

counter_total = 1;  % Iterates through all rows and columns
counter_plots = 1;  % Iterates through all subplots

peaks1 = [0 0];  % Contains global minimum and maximum of averages over all subplots
peaks2 = [0 0];  % Contains global minimum and maximum of stds over all subplots

a1 = cell(1, length(plot_index));  % Contains the axes of the average subplots
a2 = cell(1, length(plot_index));  % Contains the axes of the variance subplots

for i_rows = 1:n_rows
    for i_cols = 1:n_cols
        if sum(counter_total == plot_index) == 1
            a1{counter_plots} = axes('position', [border + border_plots + i_width * (i_cols - 1), plot_area + border + border_plots - i_height * i_rows, i_width - border_plots, 2 * (i_height - border_plots) / 3]);
            a2{counter_plots} = axes('position', [border + border_plots + i_width * (i_cols - 1), plot_area + border + border_plots - i_height * i_rows + 2 * (i_height - border_plots) / 3, i_width - border_plots, (i_height - border_plots) / 3]);
            set(f, 'CurrentAxes', a1{counter_plots});
            set(gca, 'FontSize', font_size/2*3);
            plot(r.t_plot, r.AveVar{counter_plots}.average);
            set(f, 'CurrentAxes', a2{counter_plots});
            set(gca, 'FontSize', font_size/1*3);
            plot(r.t_plot, r.AveVar{counter_plots}.stdev.^2, 'r');

            % Identical scaling for each plot
            if min(r.AveVar{counter_plots}.average) < peaks1(1)
                peaks1(1) = min(r.AveVar{counter_plots}.average);
            end;
            if max(r.AveVar{counter_plots}.average) > peaks1(2)
                peaks1(2) = max(r.AveVar{counter_plots}.average);
            end;
            if min(r.AveVar{counter_plots}.stdev.^2) < peaks2(1)
                peaks2(1) = min(r.AveVar{counter_plots}.stdev.^2);
            end;
            if max(r.AveVar{counter_plots}.stdev.^2) > peaks2(2)
                peaks2(2) = max(r.AveVar{counter_plots}.stdev.^2);
            end;

            if sum(counter_total + n_cols == plot_index) == 1  % Draw x-labels only if there is no subplot below
                set(a1{counter_plots}, 'XTickLabel', '');
            end;
            if sum(counter_total - 1 == plot_index) == 1 && i_cols ~= 1  % Draw y-labels only if there is no subplot to the left
                set(a1{counter_plots}, 'YTickLabel', '');
                set(a2{counter_plots}, 'YTickLabel', '');
            end;
            
            set(a2{counter_plots}, 'XTickLabel', '');
            
            counter_plots = counter_plots + 1;
        end;
        counter_total = counter_total + 1;
    end;
end;

for k = 1:length(a1)
    axis(a1{k}, [r.t_plot(1), r.t_plot(end), peaks1(1), peaks1(2)]);
    % Draw line for cue
    if isfield(r, 'cue')
        set(f, 'CurrentAxes', a1{k});
        v = axis;
        line([r.cue, r.cue], [v(3), v(4)], 'Color', 'k');
    end;
end;
for k = 1:length(a2)
    axis(a2{k}, [r.t_plot(1), r.t_plot(end), peaks2(1), peaks2(2)]);
    % Draw line for cue
    if isfield(r, 'cue')
        set(f, 'CurrentAxes', a2{k});
        v = axis;
        line([r.cue, r.cue], [v(3), v(4)], 'Color', 'k');
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
temp_str{1} = ['{\bfAverage/variance maps 0.7.}{\rm Calculated on ', r.date_calc, '.}'];

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

temp_str{2} = ['Trials: ', num2str(r.n_trials), ', classes: ', classes_str, ', fs: ', num2str(r.fs), 'Hz, time: ', t_str];
text(0, 0, temp_str, 'FontUnits', 'normalized', 'FontSize', 1/6, 'VerticalAlignment', 'Top', 'Interpreter', 'Tex');
%line([0, 1], [0, 0], 'Color', 'k');

%hndl = findobj('parent', gcf, 'type', 'axes');
%for a = 1:length(hndl)
%    set(findobj('parent', hndl(a)), 'ButtonDownFcn', 'zoomMap(r)');
%end;