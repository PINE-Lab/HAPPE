% NOTE: THIS IS A MODIFICATION OF fooof_group.m and fooof.m from 
% EEGLAB's specparam plugin

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

function fooof_results = calcFooof(freqs, psds, f_range, settings, return_model, return_pwe)
    if ~exist('return_model', 'var')
        return_model = false; 
    end
    if ~exist('return_pwe', 'var')
        return_pwe = false; 
    end
    % Check settings - get defaults for those not provided
    settings = fooof_check_settings(settings);

    % Initialize object to collect FOOOF results
    fooof_results = [];

    % Run FOOOF across the group of power spectra
    for psd = psds
        cur_results = fooofHelper(freqs, psd', f_range, settings, return_model, return_pwe);
        fooof_results = [fooof_results, cur_results];
    end

end

% fooof() - Fit the FOOOF model on a neural power spectrum.
%
% Usage:
%   >> fooof_results = fooof(freqs, power_spectrum, f_range, settings, return_model);
%
% Inputs:
%   freqs           = row vector of frequency values
%   power_spectrum  = row vector of power values
%   f_range         = fitting range (Hz)
%   settings        = fooof model settings, in a struct, including:
%       settings.peak_width_limts
%       settings.max_n_peaks
%       settings.min_peak_height
%       settings.peak_threshold
%       settings.aperiodic_mode
%       settings.verbose
%   return_model    = boolean of whether to return the FOOOF model fit, optional
%
% Outputs:
%   fooof_results   = fooof model ouputs, in a struct, including:
%       fooof_results.aperiodic_params
%       fooof_results.peak_params
%       fooof_results.gaussian_params
%       fooof_results.error
%       fooof_results.r_squared
%       if return_model is true, it also includes:
%            fooof_results.freqs
%            fooof_results.power_spectrum
%            fooof_results.fooofed_spectrum
%            fooof_results.ap_fit
%
% Notes
%   Not all settings need to be defined by the user.
%     Any settings that are not provided are set to default values.
%     To run with all defaults, input settings as an empty struct.

function fooof_results = fooofHelper(freqs, power_spectrum, f_range, settings, return_model, return_pwe)

    % Check settings - get defaults for those not provided
    settings = fooof_check_settings(settings);

    % Convert inputs
    freqs = py.numpy.array(freqs);
    power_spectrum = py.numpy.array(power_spectrum);
    f_range = py.list(f_range);

    % Initialize FOOOF object
    fm = py.fooof.FOOOF(settings.peak_width_limits, ...
                        settings.max_n_peaks, ...
                        settings.min_peak_height, ...
                        settings.peak_threshold, ...
                        settings.aperiodic_mode, ...
                        settings.verbose);

    % Run FOOOF fit
    fm.fit(freqs, power_spectrum, f_range)

    % Extract outputs
    fooof_results = fm.get_results();
    fooof_results = fooof_unpack_results(fooof_results);
    
    % Re-calculating r-squared
    %   r_squared doesn't seem to get computed properly (in NaN).
    %   It is unclear why this happens, other than the error can be traced
    %   back to the internal call to `np.cov`, and fails when this function
    %   gets two arrays as input.
    %   Therefore, we can simply recalculate r-squared
    coefs = corrcoef(double(py.array.array('d', fm.power_spectrum)), ...
                     double(py.array.array('d', fm.fooofed_spectrum_))); % ERROR PRONE
    fooof_results.r_squared = coefs(2);

    % Also return the actual model fit, if requested
    %   This will default to not return model, if variable not set
    if exist('return_model', 'var') && return_model

        % Get the model, and add outputs to fooof_results
        model_out = fooof_get_model(fm);
        for field = fieldnames(model_out)'
            fooof_results.(field{1}) = model_out.(field{1});
        end

    end

    % Also return pointwise error, if requested
    %   This will default to not return pwe, if variable not set
    if exist('return_pwe', 'var') && return_pwe

        % Get the pointwise error, and add outputs to fooof_results
        pwe_out = py.fooof.analysis.compute_pointwise_error_fm(fm, false, true);
        fooof_results.pointwise_error = double(pwe_out.tolist());

    end

end
