function r = calcAveVar(s, h, t, varargin)
% Calculates the mean and variance of each channel.
%
% This function calculates the mean and variance of each channel. The results
% can be displayed with the function plotAveVar.
%
% Usage:
%   r = calcAveVar(s, h, t);
%
% Input parameters:
%   s         ... Input signal (as obtained by sload) <TxC>.
%   h         ... Header structure (as obtained by sload) <1x1 struct>. Only the
%                 following fields are required: h.SampleRate, h.TRIG,
%                 h.Classlabel.
%   t         ... Start point, time resolution and end point within a trial (in 
%                 s) <1x3>. If the second value is 0, the time resolution
%                 corresponds to 1/fs.
%                 Example:
%                   t = [0, 0.25, 8];
%                   This corresponds to a trial that starts 0s after the trial
%                   start event and lasts until second 8. 4 values per second
%                   are calculated.
%
% Optional input parameters (variable argument list):
%   'class'    ... List of classes used in the calculation <1xM>.
%                  Default: all available classes are used.
%   'bandpass' ... Bandpass filter cutoff frequencies (in Hz) <1 x 2>.
%                  Default: No bandpass filter.
%   'heading'  ... Heading of the plot <string>.
%                  Default: No heading is used.
%   'montage'  ... Topographic layout of channels <NxM>. This matrix consists of
%                  zeros and ones. The channels are arranged in N rows and M
%                  columns on the plot, and they are located where the values of
%                  the matrix are equal to 1.
%                  Default: A rectangular layout is used.
%   'cue'      ... Draws a vertical line at the location of the cue (in s)
%                  <1x1>.
%                  Default: No cue is drawn.
%
% Output parameter:
%   r ... Structure containing the results <1x1 struct>.

% Copyright by Clemens Brunner
% $Revision: 0.61 $ $Date: 02/23/2009 15:54:05 $
% E-Mail: clemens.brunner@tugraz.at

% Revision history:
%   0.61: Rounded all conversions from time to samples.
%   0.60: Complete rewrite of the code to adapt it to the new toolbox.
%   0.50: Triggering of artifact-selected data always works now
%   0.40: Implement all optional arguments as a variable argument list

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
    error('No input signal specified.');
end;
if (nargin < 2)
    error('No header structure specified.');
end;
if (nargin < 3)
    error('Trial timing definition not specified (start, time resolution, end).');
end;

if length(t) == 2  % If only start and end time are provided, use full resolution
    t = [t(1), 0, t(2)];
end;

% Default parameters, can be overwritten by varargin
class = [];  % All classes are used in the calculation
bandpass = [];  % No bandpass filter
heading = [];  % Default heading
montage = [];  % Default montage
cue = [];  % Default do not draw cue

% Overwriting default values with optional input parameters
if ~isempty(varargin)  % Are there optional parameters available?
    k = 1;
    while k<=length(varargin)
        if strcmp(varargin{k}, 'class')
            class = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'montage')
            montage = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'bandpass')
            bandpass = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'cue')
            cue = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'heading')
            heading = varargin{k + 1};
            k = k + 2;
        else  % Ignore unknown parameters
            k = k + 2;
        end;
    end;
end;

if isempty(class)
    class = unique(h.Classlabel);
end;

fs = h.SampleRate;
triallen = round((t(3) - t(1)) * fs) + 1;  % Trial length (in samples)
n_trials = length(h.TRIG(ismember(h.Classlabel, class)));

if t(2) ~= 0
    t_vec_r = t(1):t(2):t(3);  % Time vector, reduced resolution
end;

if ~isempty(bandpass)
    b = fir1(fs, bandpass/(fs/2));
    s = filtfilt(b, 1, s);
end;

% Trigger data
s_t = trigg(s, h.TRIG(ismember(h.Classlabel, class)), round(t(1)*fs), round(t(3)*fs));

for chn = 1:size(s, 2)
    temp = reshape(s_t(chn,:)', triallen, n_trials)';
    
    if t(2) ~= 0
        idx = round(t_vec_r*fs) + 1;  % Indices of reduced resolution in original time vector
        average = mean(temp);
        stdev = std(temp);
        
        r.AveVar{chn}.average = average(idx);
        r.AveVar{chn}.stdev = stdev(idx);
    else
        r.AveVar{chn}.average = mean(temp);
        r.AveVar{chn}.stdev = std(temp);
    end;
    
    
end;

% Create time vector
if t(2) ~= 0
    r.t_plot = t(1):t(2):t(3);
else  % Full time resolution
    r.t_plot = t(1):1/fs:t(3);
end;
r.t = t;
r.fs = fs;
if isfield(h, 'FILE')
    if isfield(h.FILE, 'Name') && isfield(h.FILE, 'Ext')
        if length(h.FILE) > 1
            for k = 1:length(h.FILE)
                r.fname{k} = [h.FILE(k).Name, '.', h.FILE(k).Ext];
            end;
        else
            r.fname = [h.FILE.Name, '.', h.FILE.Ext];
        end;
    end;
end;
r.date_calc = datestr(clock);
if isfield(h, 'T0')  % Recording date
    r.date_rec = datestr(h.T0);
end;
r.n_trials = n_trials;
if ~isempty(heading)
    r.heading = heading;
end;
if ~isempty(class)
    r.classes = class;
else
    r.classes = unique(h.Classlabel);
end;
if ~isempty(montage)
    r.montage = montage;
end;
if ~isempty(cue)
    r.cue = cue;
end;