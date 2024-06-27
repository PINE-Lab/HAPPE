function STUDY = std_fooof(STUDY, ALLEEG, clusters, fit_mode, f_range, settings)
    % Author: The Voytek Lab and Brian Barry 
    % Calls FOOOF wrapper on spectral data from EEGLAB 
    
    % TODO: extract study design to accommodate fooof_group()
    % - incorporate statistics 
    
    % For fooof related docs see: https://github.com/fooof-tools/fooof_mat/blob/master/fooof_mat/
    
    % For relevant EEGLAB related docs see:
    % - For study: https://github.com/sccn/eeglab/blob/develop/functions/studyfunc/std_specplot.m
    
    % Inputs:
    % STUDY, ALLEEG
    % clusters is a range or array of clusters, e.g. 3:14
    % fit_mode: if 'group' (default), performs fooof_group for each design variable spectrum
    %    if 'across design' performs fooof_group on spectrum of averaged design variables
    %    if 'individual' performs fooof on each averaged design variable

    % FOOOF specific inputs:  
    %   f_range         = fitting range (Hz)
    %   psds = matrix of power values, which each row representing a spectrum (in single case power_spectrum  = row vector of power values)
    %   settings        = fooof model settings, in a struct, including:
    %       settings.peak_width_limts
    %       settings.max_n_peaks
    %       settings.min_peak_height
    %       settings.peak_threshold
    %       settings.aperiodic_mode
    %       settings.verbose
    %    return_model    = boolean of whether to return the FOOOF model fit, optional
    
    % Current outputs (default from fooof_mat)
    %   fooof_results   = fooof model ouputs, in a struct, including:
    %       fooof_results.aperiodic_params
    %       fooof_results.peak_params
    %       fooof_results.gaussian_params
    %       fooof_results.error
    %       fooof_results.r_squared
    %       for single case only: if return_model is true, it also includes:
    %            fooof_results.freqs
    %            fooof_results.power_spectrum
    %            fooof_results.fooofed_spectrum
    %            fooof_results.ap_fit
    
    if ~exist('fit_mode', 'var')
        fit_mode = 'group';
    end

    if ~exist('settings', 'var')
        settings = struct();
    end
    
    if ~exist('plot_model', 'var')
        plot_model = false; % shape needs to be reexamined
    end
    
    design_var = STUDY.design(STUDY.currentdesign).variable.value; %cell array of design variables
    std_fooof_results = cell([numel(STUDY.cluster), 1]); % indexed by cluster
    
    for c = clusters
        [STUDY, specdata, specfreqs] = std_specplot(STUDY,ALLEEG, 'clusters', c, 'noplot', 'on'); 
        
        if strcmpi(fit_mode, 'group')
            fooof_results_c = cell([numel(design_var), 1]); %fooof results for a particular cluster, arranged by design variable
            for v = 1:numel(design_var) 
                results_v = fooof_group(specfreqs, specdata{v}, f_range, settings, true); %specdata at design variable v, shape {powers x trials}
                fooof_results_c{v} = results_v;
            end
            std_fooof_results{c} = fooof_results_c;

        elseif strcmpi(fit_mode, 'across design')
            design_spec = cell([numel(design_var), 1]); %this is what you see when you plot spectrum for a cluster in EEGLAB
            for v = 1:numel(design_var)
                spec_mean = mean(specdata{v}, 2); 
                design_spec{v} = spec_mean; 
            end
            fooof_results_c = fooof_group(specfreqs, horzcat(design_spec{:}), f_range, settings, true); %horzcat makes designspec dims |psds| x #design variables
            std_fooof_results{c} = fooof_results_c;
        
        % elseif strcmpi(fit_mode, 'individual')
        %     fooof_results_c = cell([numel(design_var), 1]);
        %     for v = 1:numel(design_var) 
        %         spec_mean = mean(specdata{v}, 2);
        %         results_v = fooof(specfreqs, spec_mean, f_range, settings, return_model);
        %         fooof_results_c{v} = results_v;
        %     end
        %     std_fooof_results{c} = fooof_results_c;

        end
    end

    STUDY.etc.FOOOF_results =  std_fooof_results;