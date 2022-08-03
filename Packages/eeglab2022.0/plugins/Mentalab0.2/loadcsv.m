function [EEG, ORN, com] = loadcsv(filepath)
    com = '';
    EEG = [];
    EEG = eeg_emptyset;

    ORN = [];
    ORN = eeg_emptyset;

    orn_srate = 20; % Sampling rate of ORN data

    [directory, filename, ~] = fileparts(filepath);
    idx_final_underscore = find(filename == '_', 1, 'last');
    name = extractBefore(filename, idx_final_underscore);
    
    fullfile_EXG = fullfile(directory, append(name, '_ExG.csv'));
    EXG = importdata(fullfile_EXG, ',', 1);
    eeg_data = sortrows(EXG.data, 1);
    eeg_data = eeg_data(:, 2:end)';
    eeg_ch_names = EXG.textdata(2:end);
    eeg_timestamps = EXG.data(:, 1);
    
    fullfile_ORN = fullfile(directory, append(name, '_ORN.csv'));
    ORN = importdata(fullfile_ORN, ',', 1);
    orn_data = sortrows(ORN.data, 1);
    orn_data = orn_data(:, 2:end)';
    orn_ch_names = ORN.textdata(2:end);
    orn_timestamps = ORN.data(:, 1);
    
    fullfile_marker = fullfile(directory, append(name, '_Marker.csv'));
    marker = readtable(fullfile_marker); % works better when file contains strings and numerics than importdata
    eeg_marker = [];
    orn_marker = [];
    if numel(marker) > 0
        [eeg_marker, orn_marker] = getMarkerIdxs(eeg_timestamps, orn_timestamps, marker);
    end

    sample_rate = getSamplingRate(eeg_timestamps);
    
    % Convert to EEGLAB structure
    eeg_chanlocs = struct('labels', eeg_ch_names);
    EEG = pop_importdata('dataformat', 'array', 'nbchan', ...
        size(eeg_data, 1), 'data', eeg_data, 'setname', 'raw_eeg', ...
        'srate', sample_rate, 'xmin', 0, 'chanlocs', eeg_chanlocs);
    EEG = eeg_checkset(EEG);

    orn_chanlocs = struct('labels', orn_ch_names);
    ORN = pop_importdata('dataformat', 'array', 'nbchan', 9, 'data', ...
        orn_data, 'setname', 'raw_orn', 'srate', orn_srate, 'xmin', 0, ...
        'chanlocs', orn_chanlocs);
    ORN = eeg_checkset(ORN);
    
    if (isempty(eeg_marker))
        return;
    end

    EEG = pop_importevent(EEG, 'event', eeg_marker, 'fields', ...
        {'latency', 'type'}, 'timeunit', NaN);
    EEG = eeg_checkset(EEG);

    ORN = pop_importevent(ORN, 'event', orn_marker, 'fields', ...
        {'latency', 'type'}, 'timeunit', NaN);
    ORN = eeg_checkset(ORN);
end    