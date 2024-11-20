function std_fooofplot(STUDY, ALLEEG, cluster_, f_range, log_freqs, settings)
    % plots FOOOF spectral fit for all conditions (averaged over components)
    
    if ~isfield(STUDY.etc, 'FOOOF_results')
        disp('FOOOF data has not yet been saved to STUDY structure.')
    end

    design_var = STUDY.design(STUDY.currentdesign).variable.value; %cell array of design variables
    [STUDY, specdata, specfreqs] = std_specplot(STUDY,ALLEEG, 'clusters', cluster_, 'noplot', 'on'); 
    
    % UNCOMMENT FOR TESTING
    % load('~/Desktop/IversenLab/external/data/spectra/dip_only/brian_diponly_3_spectra.mat');
    % design_var = {'a','b','c','d'};
    % settings.verbose = false;
    %%%%%%%%
    
    design_spec = cell([numel(design_var), 1]); %averaged spectra – this is what you see when you call Study -> Edit/plot component clusters
    for v = 1:numel(design_var)
        spec_mean = mean(specdata{v}, 2);
        design_spec{v} = spec_mean; 
    end
    fooof_results_c = fooof_group(specfreqs, horzcat(design_spec{:}), f_range, settings, true);

    % finding dimensions of subplot
    dim = sqrt(numel(design_var));
    if floor(dim) == dim
        a = dim; b = dim;
    else
        if dim - floor(dim) < 1
            a = floor(dim) + 1; b = floor(dim) + 1;
        else
            a = floor(dim); b = floor(dim) + 1;
        end
    end
    
    %plotting
    for plt_i = 1:numel(design_var)
        subplot(a, b, plt_i);
        fooof_plot(fooof_results_c(plt_i), log_freqs, false);
        title(design_var{plt_i});
    end
    hold off