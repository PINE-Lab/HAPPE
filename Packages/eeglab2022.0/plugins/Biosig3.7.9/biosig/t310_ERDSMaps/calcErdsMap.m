function r = calcErdsMap(s, h, t, f_borders, varargin)
% Calculates time-frequency (ERDS) maps.
%
% This function calculates time-frequency (ERDS) maps by using either bandpower,
% FFT or wavelets to estimate the power in specific frequency bands. Maps can be
% calculated for more than one channel at once.
%
% Usage:
%   r = calcErdsMap(s, h, t, f_borders);
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
%   f_borders ... Frequency borders (in Hz) <1xF>. Contains the borders of the
%                 frequency bands and can be used with the optional parameters
%                 'f_bandwidths' and 'f_steps' (see below).
%                 Examples:
%                   f_borders = [6, 30];
%                   The maps are calculated from 6Hz to 30Hz with the default
%                   bandwidth and in the default frequency step size.
%
%                   f_borders = [4, 12, 20, 40];
%                   The maps are calculated for the segments 4-12Hz, 12-20Hz,
%                   and 20-40Hz with the default bandwidths and step sizes.
%
% Optional input parameters (variable argument list):
%   'method'       ... Calculation method <string>. User one of the following 
%                      methods: 'bp', 'fft', 'wavelet'.
%                      Default: 'bp'.
%   'f_bandwidths' ... Bandwidths for the segments specified in f_borders (in
%                      Hz) <1xF-1>.
%                      Default: 2Hz in all segments.
%   'f_steps'      ... Frequency step sizes for the segments specified in
%                      f_borders (in Hz) <1xF-1>.
%                      Default: 1Hz in all segments.
%   'class'        ... List of classes used in the calculation <1xM>.
%                      Default: all available classes are used.
%   'ref'          ... Reference interval (in s) <1x2>.
%                      Default: the whole trial is used as reference.
%   'submean'      ... Subtract the mean signal to suppress evoked components 
%                      <1x1 logical>.
%                      Default: true.
%   'sig'          ... Method to calculate the significance <string>. User one
%                      of the following methods: 'boot', 'boxcox' or 'none'.
%                      Default: 'none'.
%   'lambda'       ... Parameter of the Box-Cox transform <1x1>. If lambda is 0,
%                      the transform is a log-transform.
%                      Default: 0.
%   'alpha'        ... Significance level <1x1>. If 'sig' is set to 'none', this 
%                      value is ignored.
%                      Default: 0.01.
%   'heading'      ... Heading of the plot <string>.
%                      Default: No heading is used.
%   'montage'      ... Topographic layout of channels <NxM>. This matrix
%                      consists of zeros and ones. The channels are arranged in
%                      N rows and M columns on the plot, and they are located 
%                      where the values of the matrix are equal to 1.
%                      Default: A rectangular layout is used.
%   'cue'          ... Draws a vertical line at the location of the cue (in s)
%                      <1x1>.
%                      Default: No cue is drawn.
%   'wide_trials'  ... Use samples before and after the trial definition to get 
%                      rid of border effects when using the FFT method 
%                      <1x1 logical>. Note that samples before the beginning of 
%                      the first and after the end of the last trial are
%                      required.
%                      Default: true.
%   'refmethod'    ... Calculation mode <string>. 'classic' uses the classical
%                      approach with an averaged reference interval. 'trial'
%                      uses an individual reference for each trial. 'absolute'
%                      does not use references, but calculates bandpower maps
%                      instead.
%                      Default: 'classic'.
%
% Output parameter:
%   r ... Structure containing the results <1x1 struct>.

% Copyright by Clemens Brunner
% $Revision: 0.99 $ $Date: 10/29/2009 13:13:00 $
% E-Mail: clemens.brunner@tugraz.at

% Revision history:
%   0.99: Added new calculation mode that uses an individual reference for each
%         trial (right now only working with BP method).
%   0.91: Add 'wide_trials' option for FFT maps.
%   0.90: Clean up code, move unnecessary parameters to plot function. Enhance
%         functionality such as the support for different spacings in multiple
%         frequency bands.
%   0.80: Complete rewrite (make this function a wrapper for bandpower, FFT, and
%         wavelet maps).
%   0.70: Implement all optional arguments as a variable argument list.

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
if (nargin < 4)
    error('Frequency definition not specified.');
end;

if length(t) == 2  % If only start and end time are provided, use full resolution
    t = [t(1), 0, t(2)];
end;

% Default parameters, can be overwritten by varargin
method = 'bp';  % Use bandpower method to calculate the maps
f_bandwidths = 2 * ones(1, length(f_borders) - 1);  % Use 2Hz bands for all segments
f_steps = ones(1, length(f_borders) - 1);  % Use step size of 1Hz for all segments
class = [];  % All classes are used in the calculation
ref = [t(1), t(3)];  % Take the whole trial as reference
submean = true;  % Subtract mean before calculation
sig = 'none';  % No significance test
alpha = 0.01;  % Default significance level of 1%
lambda = 0;  % Default Box-Cox transform: log transform
heading = [];  % Default heading
montage = [];  % Default montage
cue = [];  % Default do not draw cue
wide_trials = true;  % Default use additional samples outside of trial for FFT method
refmethod = 'classic';
warning('The parameter "refmethod" is currently only working with the BP method.');

% Overwriting default values with optional input parameters
if ~isempty(varargin)  % Are there optional parameters available?
    k = 1;
    while k <= length(varargin)
        if strcmp(varargin{k}, 'method')
            method = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'f_bandwidths')
            f_bandwidths = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'f_steps')
            f_steps = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'class')
            class = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'ref')
            ref = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'submean')
            submean = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'sig')
            sig = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'lambda')
            lambda = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'alpha')
            alpha = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'heading')
            heading = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'montage')
            montage = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'cue')
            cue = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'wide_trials')
            wide_trials = varargin{k + 1};
            k = k + 2;
        elseif strcmp(varargin{k}, 'refmethod')
            refmethod = varargin{k + 1};
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

% Test if the input frequency range makes sense
if f_borders(1) < 0
    error('Starting frequency (%d) cannot be negative.', f_borders(1));
end;
if f_borders(end) > fs/2
    error('End frequency (%d) must be smaller than or equal to the Nyquist frequency (%d).', f_borders(end), fs/2);
end;
if ~issorted(f_borders)
    error('Frequencies must be monotonically increasing.');
end;
if length(f_bandwidths) ~= length(f_borders) - 1
    error('The number of bandwidths is incorrect.');
end;
if length(f_steps) ~= length(f_borders) - 1
    error('The number of frequency steps is incorrect.');
end;

% Test if the input start and end values make sense
if h.TRIG(1) + t(1)*fs < 0
    error('Not enough data before first trial available.');
end;
if h.TRIG(end)+t(3)*fs > size(s,1)
    error('Not enough data after the last trial available.');
end;
if t(3) <= t(1)
    error('Trial end point must be greater than starting point.');
end;

% Test if the reference interval is within the trial time window
if ref(1) < t(1) || ref(2) > t(3)
    error('Reference interval must be within the trial time window.')
end;

switch lower(method)
    case 'bp'
        r = calcErdsMapBP(s, h, t, f_borders, f_bandwidths, f_steps, class, ref, submean, sig, lambda, alpha, refmethod);
    case 'fft'
        r = calcErdsMapFFT(s, h, t, f_borders, f_bandwidths, f_steps, class, ref, submean, sig, lambda, alpha, wide_trials);
    case 'wavelet'
        error('The wavelet method has not been implemented yet.');
    otherwise
        error('Unknown calculation method %s.', method);
end;

r.f_borders = f_borders;
r.f_bandwidths = f_bandwidths;
r.f_steps = f_steps;

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
r.class = class;
r.ref = ref;
r.submean = submean;
r.sig = sig;
r.alpha = alpha;
r.lambda = lambda;
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
r.ERDS_method = method;