function eeg_fooofplot(EEG, epoch_range, percent, datatype, id, f_range, log_freqs, settings)
    % Plot FOOOF results for individual ICs channels
    % id: component or channel id
    %   datatype = "channel" or "component"
    
    if ~exist('epoch_range', 'var')
        epoch_range = [EEG.xmin*1000, EEG.xmax*1000];
    end
    
    if ~exist('percent', 'var')
        percent = 100;
    end

    if ~exist('settings', 'var')
        settings = struct();
    end

    if datatype == "component"
        [eegspecdB, specfreqs, compeegspecdB] = pop_spectopo(EEG, 0, epoch_range, 'EEG', 'percent', percent , 'freq', [10], 'freqrange',f_range, 'plotchan', 0, 'icacomps', id, 'nicamaps', 0,'electrodes','off', 'plot', 'off'); %icacomps doesn't do anything
        specfreqs = specfreqs';  %reshaping frequencies
        compeegspecdB = compeegspecdB(id,:); % only selecting desired rows since pop_spectopo returns all components
        specdata = arrayfun(@(y) 10^(y/10), compeegspecdB); %undoing the 10*log10(power) transformation
        fooof_results = fooof(specfreqs, specdata, f_range, settings, true);
    else %channel case
        [eegspecdB, specfreqs] = pop_spectopo(EEG, 1, epoch_range, 'EEG' , 'percent', percent, 'freq', [10], 'freqrange',f_range,'electrodes','off', 'plot', 'off');
        specfreqs = specfreqs';  %reshaping frequencies
        eegspecdB = eegspecdB(id,:); 
        specdata = arrayfun(@(y) 10^(y/10), eegspecdB); % reshaping + undoing the 10*log10(power) transformation
        
        fooof_results = fooof(specfreqs, specdata, f_range, settings, true);
    end
    %plotting
    fooof_plot(fooof_results, log_freqs, true);
    title(datatype + " " + int2str(id));
    hold off
    