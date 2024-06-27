function EEG = eeg_fooof(EEG, datatype, ids, epoch_range, percent,  f_range, settings)
    % Author: The Voytek Lab and Brian Barry 
    % Calls FOOOF wrapper on spectral data from EEGLAB
    
    % calls pop_spectopo.m to retrieve desired channel or component spectra 
    % calls fooof_group even if only one component queried.
    
    % Inputs:
    %   EEG
    %   epoch_range = [min_ms, max_ms]
    %   percent  = [float 0 to 100] percent of the data to sample for computing the spectra. Values < 100 speed up the computation. {default: 100}
    %   f_range         = spectral computation range AND fitting range (Hz)
    %   ids = ICA components to include
    %   settings        = fooof model settings, in a struct, including:
    %       settings.peak_width_limits
    %       settings.max_n_peaks
    %       settings.min_peak_height
    %       settings.peak_threshold
    %       settings.aperiodic_mode
    %       settings.verbose
    
    % Current outputs (default from fooof_mat)
    %   fooof_results   = fooof model ouputs, in a struct, including:
    %       fooof_results.aperiodic_params
    %       fooof_results.peak_params
    %       fooof_results.gaussian_params
    %       fooof_results.error
    %       fooof_results.r_squared
    %       and with return_model hard-coded as true, it also includes:
    %            fooof_results.freqs
    %            fooof_results.power_spectrum
    %            fooof_results.fooofed_spectrum
    %            fooof_results.ap_fit
    
    if ~exist('epoch_range', 'var')
        epoch_range = [EEG.xmin*1000, EEG.xmax*1000];
    end
    
    if ~exist('percent', 'var')
        percent = 100;
    end

    if ~exist('settings', 'var')
        settings = struct();
    end
    
    % compute spectrum and fit
    % extracting power spectrum w/ pop_spectopo
    % psds as ___specdB = matrix of power values, with each row representing a spectrum)
    if datatype == "component"
        [eegspecdB, specfreqs, compeegspecdB] = pop_spectopo(EEG, 0, ...
            epoch_range, 'EEG', 'percent', percent , 'freq', [10], ...
            'freqrange',f_range, 'plotchan', 0, 'icacomps', ids, ...
            'nicamaps', 0,'electrodes','off', 'plot', 'off'); % icacomps actually useless
        compeegspecdB = compeegspecdB(ids,:); % only selecting desired rows (or should it be cols?) since pop_spectopo returns all components
        specdata = arrayfun(@(y) 10^(y/10), compeegspecdB'); % reshaping + undoing the 10*log10(power) transformation
        specfreqs = specfreqs';  %reshaping frequencies
        fooof_results = cell([size(EEG.icaweights,1), 1]); % indexed by component
        
        % computing FOOOF
        fooof_results_temp = fooof_group(specfreqs, specdata, f_range, settings, true); 
        for i = 1:numel(ids)
            comp_i = ids(i);
            fooof_results{comp_i} = fooof_results_temp(i);
        end

    else % channel case
        [eegspecdB, specfreqs] = pop_spectopo(EEG, 1, epoch_range, 'EEG', ...
            'percent', percent, 'freq', [10], 'freqrange', f_range, ...
            'electrodes','off', 'plot', 'off');
        eegspecdB = eegspecdB(ids,:);
        specdata = arrayfun(@(y) 10^(y/10), eegspecdB'); % reshaping + undoing the 10*log10(power) transformation
        specfreqs = specfreqs';  % reshaping frequencies
        fooof_results = cell([size(EEG.chanlocs,2), 1]);
        
        % computing FOOOF
        fooof_results_temp = fooof_group(specfreqs, specdata, f_range, settings, true); 
        for i = 1:numel(ids)
            comp_i = ids(i);
            fooof_results{comp_i} = fooof_results_temp(i);
        end 
    end

    EEG.etc.FOOOF_results = fooof_results; 