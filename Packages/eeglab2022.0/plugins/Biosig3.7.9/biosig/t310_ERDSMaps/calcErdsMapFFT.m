function r = calcErdsMapFFT(s, h, t, f_borders, f_bandwidths, f_steps, class, ref, submean, sig, lambda, alpha, wide_trials)
% Calculates time-frequency (ERDS) maps based on the FFT.
%
% This function calculates time-frequency (ERDS) maps by using the FFT to 
% estimate the power in specific frequency bands. Maps can be calculated for
% more than one channel at once. Note that this function should not be used
% directly, use calcErdsMap instead and specify 'fft' as the calculation method.
%
% Usage:
%   r = calcErdsMapFFT(s, h, t, f_borders, f_bandwidths, f_steps, class, ref, 
%                      submean, sig, lambda, alpha, wide_trials);
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
%                    Default: 0.
%   alpha        ... Significance level <1x1>. If 'sig' is set to 'none', this
%                    value is ignored.
%   wide_trials  ... Use samples before and after the trial definition to get
%                    rid of border effects when using the FFT method 
%                    <1x1 logical>. Note that samples before the beginning of 
%                    the first and after the end of the last trial are required.
%
% Output parameter:
%   r ... Structure containing the results <1x1 struct>.

% Copyright by Clemens Brunner
% $Revision: 0.21 $ $Date: 02/23/2009 15:54:39 $
% E-Mail: clemens.brunner@tugraz.at

% Revision history:
%   0.21: Rounded all conversions from time to samples.
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

triallen = round((t(3) - t(1)) * fs) + 1;  % Trial length (in samples)
if t(2) ~= 0
    t_step=round(t(2)*fs)/fs;
    t_vec_r = t(1):t(2):t(3);  % Time vector, reduced resolution
end;

% Determine indices of reference interval
if t(2) ~= 0
    ref(1) = round(ref(1) / t(2)) + 1 - round(t(1) / t(2));
    ref(2) = round(ref(2) / t(2)) + 1 - round(t(1) / t(2));
else
    ref(1) = round(ref(1) * fs) + 1 - round(t(1) * fs);
    ref(2) = round(ref(2) * fs) + 1 - round(t(1) * fs);
end;

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

for chn = 1:size(s,2)  % Loop over all channels

    % Trigger data
    if wide_trials
        % Number of samples before the beginning and after the end of each trial
        n_add = round(1/f_low(1)*fs);
        s_t = trigg(s(:, chn), h.TRIG(ismember(h.Classlabel, class)), round(t(1)*fs) - n_add, round(t(3)*fs) + n_add);
        triallen = triallen + 2 * n_add;
    else
        s_t = trigg(s(:, chn), h.TRIG(ismember(h.Classlabel, class)), round(t(1)*fs), round(t(3)*fs));
    end;
            
    if submean
        temp = reshape(s_t, triallen, length(s_t)/triallen);  % Reshape to samples x trials
        average = mean(temp, 2);
        temp = temp - repmat(average, 1, size(temp,2));
        s_t = reshape(temp, 1, size(temp,1)*size(temp,2));
    end;
    
    if t(2) ~= 0
        erds = zeros(length(t_vec_r), length(f_plot));  % One frequency band per column
    else
        erds = zeros(triallen, length(f_plot));  % One frequency band per column
    end;
    refp = zeros(1, length(f_plot));  % One reference value per column
    
    s_reshape = reshape(s_t, triallen, length(s_t)/triallen);
    
    % The variable pre_erds is needed for the bootstrap statistics, it contains the single-trial ERDS values
    % pre_erds{frequency band} = <time x trials>
    pre_erds = cell(1, length(f_plot));  % Create empty cells, one for each frequency band
    for k = 1:length(f_plot)  % Fill the cells with zeros
        if t(2) ~= 0
            pre_erds{k} = zeros(length(t_vec_r), size(s_reshape,2));
        else
            pre_erds{k} = zeros(triallen, size(s_reshape,2));
        end;
    end;
    
    if wide_trials
        triallen = triallen - 2 * n_add;
    end;
    
    % Loop over different frequency ranges
    for k = 1:length(f_borders) - 1
        % Calculation of the required frequency resolution and the number of points for the FFT
        f_res = min(f_steps(k),f_bandwidths(k)/2);
        
        % Required number of points for the FFT
        nfft = 2^nextpow2(fs/f_res);  % Take the next higher power of 2 to avoid using DFT instead of FFT
        f = fs/2*linspace(0,1,nfft/2);
        
        % Create Hanning windows with suitable length
        % (depending on the lower frequency bound of the current frequency range)
        % Alternatively, the windows could be of constant length
        size_window = round(1/f_low(f_plot == f_borders(k))*fs) * 2;
        window_hann = repmat(hanning(size_window), 1, size(s_reshape,2));
        
        % Zero-padding at the beginning and end of the trials
        if wide_trials
            s_padded = s_reshape(n_add - ceil(size_window/2) + 1:triallen + n_add + ceil(size_window/2), :);  % Signal has already been padded, but with actual data instead of zeros.
        else
            s_padded = [zeros(ceil(size_window/2), size(s_reshape,2))', s_reshape', zeros(ceil(size_window/2), size(s_reshape,2))']';
        end;
            
        % Calculate the FFT for each time step
        idx = 1;  % Number of sample in trial (full time resolution)
        counter = 1;  % Number of sample in the ERDS of the trial (reduced time resolution)
        while idx <= triallen
            s_padded_fft = fft(s_padded(idx:idx+size_window-1,:).*window_hann, nfft)/nfft;
            s_padded_fft = s_padded_fft.*conj(s_padded_fft);
            
            start_f = find(f_plot == f_borders(k));
            if k < length(f_borders) - 1
                end_f = find(f_plot == f_borders(k+1)) - 1;
            else
                end_f = find(f_plot == f_borders(k+1));
            end;
            
            for l = start_f:end_f
                pre_erds{l}(counter,:) = sum(s_padded_fft(find(f>=f_low(l) & f<=f_up(l)),:));
            end;

            if t(2) ~= 0
                idx = idx + t_step*fs;
            else
                idx = idx + 1;
            end;
            counter = counter + 1;
        end;
    end;
    
    for k = 1:length(f_plot)
        activity = mean(pre_erds{k}, 2);  % Activity power
        refp(k) = mean(mean(pre_erds{k}(ref(1):ref(2),:)));  % Reference power
        erds(:,k) = activity./refp(k) - 1;
        pre_erds{k} = pre_erds{k}./refp(k) - 1;
    end;
    r.ERDS{chn}.erds = erds;
    r.ERDS{chn}.refp = refp;

    switch lower(sig)

        case 'boot'
            for k = 1:length(f_plot)
                [temp, cl(:,k), cu(:,k)] = bootts(pre_erds{k}', 300, alpha);
                smooth_length = ceil(2*fs/f_low(k));
                r.ERDS{chn}.cl(:,k) = filtfilt(ones(1, smooth_length)/smooth_length, 1, cl(:,k));
                r.ERDS{chn}.cu(:,k) = filtfilt(ones(1, smooth_length)/smooth_length, 1, cu(:,k));
            end;

        case 'boxcox'
            % t-test formula from  Hans-Jochen Bartsch 20th edition
            % page 729
            % degrees of freedom is number of samples - 1 
            df = size(pre_erds{1},2) - 1; % degrees of freedom

            for k = 1:length(f_plot)
                % All data points must be >= 0
                offset = floor(min(min(pre_erds{k})));
                if offset < 0
                    pre_erds{k} = pre_erds{k} + abs(offset);  % Add offset
                end;
                % x_mean ... meanof the test samples
                % x0 ... predicted mean
                % ser ... standartderifation/sqrt(samples)
                % following statment must be true if x0 is a possible mean
                % abs((x_mean - x0)/ser) <= t_crit
                % calculation of the coeffidenz interval
                % abs(x_mean - x0) <= t_crit.*ser
                % c_t = x_mean +- t_crit.*ser
                if lambda == 0
                    erds_t = log(pre_erds{k});
                    erds_t_m = mean(erds_t, 2);  % Mean over trials
                    erds_t_s = std(erds_t, 0, 2);  % Standard deviation over trials
                    ser = erds_t_s ./ sqrt(df + 1);
                    t_crit = tinv((1 - alpha / 2), df);
                    cl_t = erds_t_m - t_crit .* ser;
                    cu_t = erds_t_m + t_crit .* ser;
                    r.ERDS{chn}.cl(:,k) = exp(cl_t) - abs(offset);
                    r.ERDS{chn}.cu(:,k) = exp(cu_t) - abs(offset);
                else
                    erds_t = (pre_erds{k}.^lambda - 1)./lambda;
                    erds_t_m = mean(erds_t, 2);  % Mean over trials
                    erds_t_s = std(erds_t, 0, 2);  % Standard deviation over trials
                    ser = erds_t_s ./ sqrt(df + 1);
                    t_crit = tinv((1 - alpha / 2), df);
                    cl_t = erds_t_m - t_crit .* ser;
                    cu_t = erds_t_m + t_crit .* ser;
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