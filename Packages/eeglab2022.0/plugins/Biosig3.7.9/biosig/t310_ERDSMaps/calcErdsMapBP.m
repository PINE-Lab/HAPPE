function r = calcErdsMapBP(s, h, t, f_borders, f_bandwidths, f_steps, class, ref, submean, sig, lambda, alpha, refmethod)
% Calculates time-frequency (ERDS) maps based on the bandpower (BP) method.
%
% This function calculates time-frequency (ERDS) maps by using a bandpass filter
% and subsequent squaring of the signals to estimate the power in specific
% frequency bands. Maps can be calculated for more than one channel at once.
% Note that this function should not be used directly, use calcErdsMap instead
% and specify 'bp' as the calculation method.
%
% Usage:
%   r = calcErdsMapBP(s, h, t, f_borders, f_bandwidths, f_steps, class, ref, 
%                     submean, sig, lambda, alpha);
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
%   f_bandwidths ... Bandwidths for the segments specified in f_borders (in Hz)
%                    <1xF-1>.
%   f_steps      ... Frequency step sizes for the segments specified in
%                    f_borders (in Hz) <1xF-1>.
%   class        ... List of classes used in the calculation <1xM>.
%   ref          ... Reference interval (in s) <1x2>.
%   submean      ... Subtract the mean signal to suppress evoked components 
%                    <1x1 logical>.
%   sig          ... Method to calculate the significance <string>. User one of
%                    the following methods: 'boot', 'boxcox' or 'none'.
%   lambda       ... Parameter of the Box-Cox transform <1x1>. If lambda is 0,
%                    the transform is a log-transform.
%   alpha        ... Significance level <1x1>. If 'sig' is set to 'none', this 
%                    value is ignored.
%   refmethod    ... Calculation mode <string>. 'classic' uses the classical
%                    approach with an averaged reference interval. 'trial' uses
%                    an individual reference for each trial.
%
% Output parameter:
%   r ... Structure containing the results <1x1 struct>.

% Copyright by Clemens Brunner
% $Revision: 0.30 $ $Date: 10/29/2009 13:11:00 $
% E-Mail: clemens.brunner@tugraz.at

% Revision history:
%   0.30: Added new calculation mode that uses an individual reference for each
%         trial.
%   0.21: Rounded all conversions from time to samples
%   0.20: Modifications concerning the interface and computations.
%   0.10: Initial version.

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

fs = h.SampleRate;

n_butter = 4;  % Filter order of Butterworth filter

triallen = round((t(3) - t(1)) * fs) + 1;  % Trial length (in samples)
if t(2) ~= 0
    t_vec_r = t(1):t(2):t(3);  % Time vector, reduced resolution
end;

% Determine indices of reference interval
ref(1) = round(ref(1) * fs) + 1 - round(t(1) * fs);
ref(2) = round(ref(2) * fs) + 1 - round(t(1) * fs);

f_plot = [];  % Center frequencies in the plot
f_low = [];  % Lower frequency border for each center frequency
f_up = [];  % Upper frequency border for each center frequency
for k = 1:length(f_borders) - 1
    f_plot = [f_plot, f_borders(k):f_steps(k):f_borders(k + 1)-f_steps(k)];
    f_low = [f_low, f_borders(k) - f_bandwidths(k)/2:f_steps(k):f_borders(k + 1) - f_steps(k) - f_bandwidths(k)/2];
    f_up = [f_up, f_borders(k) + f_bandwidths(k)/2:f_steps(k):f_borders(k + 1) - f_steps(k) + f_bandwidths(k)/2];
end;
f_plot = [f_plot, f_borders(end)];
f_low = [f_low, f_borders(end) - f_bandwidths(end)/2];
f_up = [f_up, f_borders(end) + f_bandwidths(end)/2];

fn  = fs/2;

for chn = 1:size(s,2)  % Loop over all channels
       
    if submean  % Subtract evoked components?
        s_t = trigg(s(:, chn), h.TRIG(ismember(h.Classlabel, class)), round(t(1)*fs), round(t(3)*fs));
        temp = reshape(s_t, triallen, length(s_t)/triallen);  % Reshape to samples x trials
        average = mean(temp, 2);
        temp = zeros(length(s), 1);
        for k = 1:length(h.TRIG)
            temp(h.TRIG(k) + round(t(1)*fs) : h.TRIG(k) + round(t(1)*fs) + length(average)-1) = average;
        end;
        s(:, chn) = s(:, chn) - temp;
    end;
    
    if strcmp(refmethod, 'classic')  % 1 reference per frequency band
        refp = zeros(1, length(f_plot));
    else  % References for each trial per frequency band
        if submean  % s_t was already defined
            refp = zeros(length(s_t)/triallen, length(f_plot));
        else
            s_t = trigg(s(:, chn), h.TRIG(ismember(h.Classlabel, class)), round(t(1)*fs), round(t(3)*fs));
            refp = zeros(length(s_t)/triallen, length(f_plot));
        end;
    end;
    erds = zeros(triallen, length(f_plot));
    
    for k = 1:length(f_plot)  % Loop over frequency bands
        [b, a] = butter(n_butter, [f_low(k), f_up(k)]./fn);
        smooth_length = ceil(2*fs/f_low(k));
        s_f = filtfilt(ones(1, smooth_length)/smooth_length, 1, filtfilt(b, a, s(:,chn)).^2);


        % Trigger data
        s_t = trigg(s_f, h.TRIG(ismember(h.Classlabel, class)), round(t(1)*fs), round(t(3)*fs));

        if strcmp(refmethod, 'classic')  % Use classical calculation scheme with trial-averaged reference
            % The variable pre_erds is needed for the bootstrap statistics, it contains the single-trial ERDS values
            % pre_erds{frequency band} = <time x trials>
            pre_erds{k} = reshape(s_t, triallen, length(s_t)/triallen);
            activity = mean(pre_erds{k}, 2);  % Average activity power over all trials
            refp(k) = mean(mean(pre_erds{k}(ref(1):ref(2),:)));  % Average reference power over all trials and the reference time segment
            erds(:,k) = activity./refp(k) - 1;  % Calculate ERDS
            pre_erds{k} = pre_erds{k}./refp(k) - 1;
        elseif strcmp(refmethod, 'trial')  % Use alternative calculation scheme with trial-individual references
            pre_erds{k} = reshape(s_t, triallen, length(s_t)/triallen);
            refp(:,k) = mean(pre_erds{k}(ref(1):ref(2),:));  % Average reference power for each trial
            pre_erds{k} = pre_erds{k}./repmat(refp(:,k)',size(pre_erds{k},1),1) - 1;
            erds(:,k) = mean(pre_erds{k},2);
        elseif strcmp(refmethod, 'absolute')  % Calculate BP maps
            pre_erds{k} = reshape(s_t, triallen, length(s_t)/triallen);
            erds(:,k) = mean(pre_erds{k},2);
        end;
    end;
    
    if t(2) ~= 0
        idx = round(t_vec_r*fs) + 1;  % Indices of reduced resolution in original time vector
        r.ERDS{chn}.erds = erds(idx,:);
        for k = 1:length(f_plot)
            pre_erds{k} = pre_erds{k}(idx,:);
        end;
    else
        r.ERDS{chn}.erds = erds;
    end;
    if ~strcmp(refmethod, 'absolute')  % We don't need reference power values for BP maps
        r.ERDS{chn}.refp = refp;
    end;
    
    switch lower(sig)

        case 'boot'
            for k = 1:length(f_plot)
                [temp, cl(:,k), cu(:,k)] = bootts(pre_erds{k}', 300, alpha);
                if t(2) ~= 0  % Confidence intervals should be smoothed
                    smooth_length = ceil(2*fs/f_low(k) / (triallen/length(idx)));
                else
                    smooth_length = ceil(2*fs/f_low(k));
                end;
                r.ERDS{chn}.cl(:,k) = filtfilt(ones(1, smooth_length)/smooth_length, 1, cl(:,k));
                r.ERDS{chn}.cu(:,k) = filtfilt(ones(1, smooth_length)/smooth_length, 1, cu(:,k));
            end;
            
        case 'boxcox'
            % FIXME: A normal distribution with known sigma is assumed, which is not
            % quite correct. It should be replaced with the Student's t-distribution.
            % This is implemented in boxcox2, but it has to be tested
            % first. In the meantime, it is safe to use boxcox.
            z = erfinv(1 - alpha) * sqrt(2);  % Estimate z for confidence interval (mu - z * sigma)

            for k = 1:length(f_plot)
                % All data points must be >= 0
                offset = floor(min(min(pre_erds{k})));
                if offset < 0
                    pre_erds{k} = pre_erds{k} + abs(offset);  % Add offset
                end;
                if lambda == 0
                    erds_t = log(pre_erds{k});
                    erds_t_m = mean(erds_t, 2);  % Mean over trials
                    erds_t_s = std(erds_t, 0, 2);  % Standard deviation over trials
                    ser = erds_t_s ./ sqrt(size(erds_t, 2));
                    cl_t = erds_t_m - z * ser;
                    cu_t = erds_t_m + z * ser;
                    r.ERDS{chn}.cl(:,k) = exp(cl_t) - abs(offset);
                    r.ERDS{chn}.cu(:,k) = exp(cu_t) - abs(offset);
                else
                   erds_t = (pre_erds{k}.^lambda - 1)./lambda;
                   erds_t_m = mean(erds_t, 2);  % Mean over trials
                   erds_t_s = std(erds_t, 0, 2);  % Standard deviation over trials
                   ser = erds_t_s ./ sqrt(size(erds_t, 2));
                   cl_t = erds_t_m - z * ser;
                   cu_t = erds_t_m + z * ser;
                   r.ERDS{chn}.cl(:,k) = (cl_t * lambda + 1).^(1/lambda) - abs(offset);
                   r.ERDS{chn}.cu(:,k) = (cu_t * lambda + 1).^(1/lambda) - abs(offset);
                end;
            end;
        case 'boxcox2'  % This version should be the correct one (using unknown variance and a t-test)
            
            for k = 1:length(f_plot)
                % All data points must be >= 0
                offset = floor(min(min(pre_erds{k})));
                if offset < 0
                    pre_erds{k} = pre_erds{k} + abs(offset);  % Add offset
                end;
                if lambda == 0
                    erds_t = log(pre_erds{k});
                    erds_t_m = mean(erds_t, 2);  % Mean over trials
                    erds_t_s = std(erds_t, 0, 2);  % Standard deviation over trials
                    ser = erds_t_s ./ sqrt(size(erds_t, 2));
                    z = tinv((1 - alpha / 2), size(erds_t, 2) - 1);
                    cl_t = erds_t_m - z * ser;
                    cu_t = erds_t_m + z * ser;
                    r.ERDS{chn}.cl(:,k) = exp(cl_t) - abs(offset);
                    r.ERDS{chn}.cu(:,k) = exp(cu_t) - abs(offset);
                else
                    erds_t = (pre_erds{k}.^lambda - 1)./lambda;
                    erds_t_m = mean(erds_t, 2);  % Mean over trials
                    erds_t_s = std(erds_t, 0, 2);  % Standard deviation over trials
                    ser = erds_t_s ./ sqrt(size(erds_t, 2));
                    z = tinv((1 - alpha / 2), size(erds_t, 2) - 1);
                    cl_t = erds_t_m - z * ser;
                    cu_t = erds_t_m + z * ser;
                    r.ERDS{chn}.cl(:,k) = (cl_t * lambda + 1).^(1/lambda) - abs(offset);
                    r.ERDS{chn}.cu(:,k) = (cu_t * lambda + 1).^(1/lambda) - abs(offset);
                end;
            end;
    end;
end;

r.f_plot = f_plot;
r.f_low = f_low;
r.f_up = f_up;
r.n_trials = length(h.TRIG(ismember(h.Classlabel, class)));
r.refmethod = refmethod;