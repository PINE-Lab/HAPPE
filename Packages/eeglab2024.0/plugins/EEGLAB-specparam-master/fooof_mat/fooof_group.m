% fooof_group() - Run the fooof model on a group of neural power spectra.
%
% Usage:
%   fooof_results = fooof_group(freqs, psds, f_range, settings);
%
% Inputs:
%   freqs           = row vector of frequency values
%   psds            = matrix of power values, which each row representing a spectrum
%   f_range         = fitting range (Hz)
%   settings        = fooof model settings, in a struct, including:
%       settings.peak_width_limts
%       settings.max_n_peaks
%       settings.min_peak_height
%       settings.peak_threshold
%       settings.aperiodic_mode
%       settings.verbose
%
% Outputs:
%   fooof_results   = fooof model ouputs, in a struct, including:
%       fooof_results.aperiodic_params
%       fooof_results.peak_params
%       fooof_results.gaussian_params
%       fooof_results.error
%       fooof_results.r_squared
%
% Notes
%   Not all settings need to be set. Any settings that are not
%     provided as set to default values. To run with all defaults,
%     input settings as an empty struct.

function fooof_results = fooof_group(freqs, psds, f_range, settings, return_model)
    if ~exist('return_model', 'var')
        return_model = false; 
    end

    % Check settings - get defaults for those not provided
    settings = fooof_check_settings(settings);

    % Initialize object to collect FOOOF results
    fooof_results = [];

    % Run FOOOF across the group of power spectra
    for psd = psds
        cur_results = fooof(freqs, psd', f_range, settings, return_model);
        fooof_results = [fooof_results, cur_results];
    end

end