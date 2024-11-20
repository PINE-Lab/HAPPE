function [STUDY LASTCOM] = pop_std_fooof(STUDY, ALLEEG)
    % Author: The Voytek Lab and Brian Barry 
    % GUI for FOOOF wrapper on spectral data from EEGLAB

    % For fooof related docs see: https://github.com/fooof-tools/fooof_mat/blob/master/fooof_mat/
    
    % For relevant EEGLAB related docs see:
    % https://github.com/sccn/eeglab/blob/develop/functions/popfunc/pop_spectopo.m
    
    % Should have
        % settings.peak_width_limits sets the possible lower- and upper-bounds for the fitted peak widths.
        % settings.max_n_peaks sets the maximum number of peaks to fit.
        % settings.min_peak_height sets an absolute limit on the minimum height (above aperiodic) for any extracted peak.
        % settings.peak_threshold sets a relative threshold above which a peak height must cross to be included in the model.
        % settings.aperiodic_mode defines the approach to use to parameterize the aperiodic component.

    uilist = { { 'style' 'text' 'string' 'Clusters:' } ...
            { 'style' 'edit' 'string' '' } ...
            { 'style' 'text' 'string' 'Fit mode ("group" or "across design"):' } ...
            { 'style' 'edit' 'string' '"group"' } ...
            { 'style' 'text' 'string' 'Frequency range to fit:' } ...
            { 'style' 'edit' 'string' '' } ...
            ... % Now FOOOF settings
            { 'style' 'text' 'string' '                     FOOOF settings (optional)' 'fontweight' 'bold' }...       
            { 'style' 'text' 'string' 'peak_width_limits' } ...
            { 'style' 'edit' 'string' '' } ...
            { 'style' 'text' 'string' 'max_n_peaks' } ...
            { 'style' 'edit' 'string' '' } ...
            { 'style' 'text' 'string' 'min_peak_height' } ...
            { 'style' 'edit' 'string' '' } ...
            { 'style' 'text' 'string' 'peak_threshold' } ...
            { 'style' 'edit' 'string' '' } ...
            { 'style' 'text' 'string' 'aperiodic_mode' } ...
            { 'style' 'edit' 'string' "'fixed'" } ... %want to make a checkmark later
            { 'style' 'text' 'string' 'verbose (boolean)' } ...
            { 'style' 'edit' 'string' 'false' } };
    uigeom = { [9 4] [9 3] [9 3] [1] [9 3] [9 3] [9 3] [9 3] [9 3] [9 3]};
    [result, usrdat, sres2, sres] = inputgui( 'uilist', uilist, 'geometry', uigeom, 'title', 'FOOOF EEG - pop_std_fooof_eeg()', 'helpcom', 'pophelp(''pop_fooof_eeg'');', 'userdata', 0); %currently ignoring usrdat, sres2, sres
    params = {}; %parameters for std_fooof_eeg w/o FOOOF settings
    settings_keys = {'peak_width_limits','max_n_peaks','min_peak_height','peak_threshold','aperiodic_mode','verbose'};
    settings = struct(); %can be empty
    for i = 1:length(result)
        if i < 4
            param_curr = eval( [ '[' result{i} ']' ] );
            params{end+1} = param_curr;
        else
            if ~isempty(eval( [ '[' result{i} ']' ] ))
                settings.(settings_keys{i-3}) = eval( [ '[' result{i} ']' ] );
            end
        end
    end
    
    if ~isempty(params)
        STUDY = std_fooof(STUDY, ALLEEG, params{1}, params{2}, params{3}, settings);
        LASTCOM = sprintf('STUDY = std_fooof(STUDY, ALLEEG, [ %s ], "%s", [%d %d],  %s)', sprintf('%d ', params{1}), params{2}, params{3}(1), params{3}(2), struct2str(settings)); 
    else                                                                                     % clusters           %fit mode      % f_range1        f_range2    settings  
        LASTCOM = '';
    end