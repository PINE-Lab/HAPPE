% Template for using the EEG analysis toolbox.
%
% This batch file serves as a template for the EEG analysis toolbox. It shows
% examples of how to call specific functions, such as plotting ERDS maps. Just
% copy and paste the commands you need into your own batch file and adapt the
% parameters. For help on specific commands, see the corresponding help texts.

% Copyright by Clemens Brunner
% $Revision: 0.3 $ $Date: 10/09/2008 15:10:29 $
% E-Mail: clemens.brunner@tugraz.at

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

%% Load data
name = 'sample';
path_name = './';
data = '*.gdf';

channels = 0;  % 0 ... All channels
classes = [1];  % [] ... All classes
method = 'bp';  % 'bp' or 'fft'
refmethod = 'classic';  % 'classic' or 'trial'

[s, h] = sload([path_name, data], channels, 'OVERFLOWDETECTION:OFF');
s(isnan(s)) = 0;

%% ERDS maps
r1 = calcErdsMap(s, h, [0, 0.05, 5], [5, 40], 'method', method, 'class', classes, 'ref', [0.25, 0.75], 'f_bandwidths', [2], 'f_steps', [1], 'sig', 'boxcox', 'lambda', 1, 'alpha', 0.05, 'heading', name, 'montage', [1 1 1 1], 'cue', 3, 'refmethod', refmethod);
plotErdsMap(r1);

%% Average/variance
r2 = calcAveVar(s, h, [0, 0.05, 8], 'class', classes, 'heading', name, 'montage', [1 1 1 1], 'cue', 3);
plotAveVar(r2);

%% Combination of ERDS maps and average/variance
r3 = calcCombiMap(s, h, [0, 0.05, 8], [5, 18, 40], 'method', method, 'class', classes, 'ref', [0.5, 1.5], 'f_bandwidths', [2, 4], 'f_steps', [1, 1], 'sig', 'boxcox', 'lambda', 1, 'alpha', 0.05, 'heading', name, 'montage', [1 1 1 1], 'cue', 3);
plotCombiMap(r3);