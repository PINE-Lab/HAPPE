function [EEG, ORN, com] = pop_loadfile(filepath, varargin)  
    com = '';

    if nargin > 0  % File provided in call
        runMain(filepath)
        return;
    end

    [filename, path] = getFileFromUser();
    [EEG, ORN, com] = runMain([path, filename]);
end


function [EEG, ORN, com] = runMain(filepath)
    [directory, filename, ext] = fileparts(filepath);
    if (contains(ext, 'csv'))
        checkFolderContents(filename, directory);
    end
    [EEG, ORN, com] = loadBINCSV(filepath, ext);

    channelNameList = requestChannelLabels(EEG);
    if length(channelNameList) > 3
        for n = 1:length(channelNameList)
            EEG.chanlocs(n).labels = channelNameList{n};
        end
    end

    EEG.filename = filename;
    ORN.filename = filename;
end


function [filename, path] = getFileFromUser() 
    [filename, path] = uigetfile({'*.BIN;*.CSV;' 'All BIN and CSV files';}, ...
        'Select a BIN or CSV file'); 
    if filename == 0
        error('---> File selection cancelled.')
    end

    if (~contains(filename, '_ExG.csv', 'IgnoreCase', true) ...
            && ~contains(filename, '_ORN.csv', 'IgnoreCase', true) ...
            && ~contains(filename, '_Marker.csv', 'IgnoreCase', true) ...
            && ~contains(filename, '.bin', 'IgnoreCase', true)) 
        error(['---> Error on: "' filename '". Unsuitable file type.' ...
            ' Please select a BIN file or CSV file with suffix "_ExG", "_ORN" or "_Marker".'])
    end
end


function [EEG, ORN, com] = loadBINCSV(filepath, ext)
    if (contains(ext, "bin", 'IgnoreCase', true))
        [EEG, ORN, com] = loadbin(filepath);
    else
        [EEG, ORN, com] = loadcsv(filepath);
    end
end


function checkFolderContents(filename, directory) % Will be CSV file
    idx_final_underscore = find(filename == '_', 1, 'last');
    name = extractBefore(filename, idx_final_underscore);

    exg_orn_marker = false(1, 3);
    files = dir(directory);
    for i = 1:length(files)
        file_i = files(i).name;
        if (contains(file_i, [name '_ExG.csv'], 'IgnoreCase', true)...
                || contains(file_i, [name '_ExG'], 'IgnoreCase', true))
            exg_orn_marker(1) = true;
        elseif (contains(file_i, [name '_ORN.csv'], 'IgnoreCase', true)...
                || contains(file_i, [name '_ORN'], 'IgnoreCase', true))
            exg_orn_marker(2) = true;
        elseif (contains(file_i, [name '_Marker.csv'], 'IgnoreCase', true)...
                || contains(file_i, [name '_Marker'], 'IgnoreCase', true))
            exg_orn_marker(3) = true;
        end
    end

    if (~exg_orn_marker(1))
        error(['---> Selected directory does not contain: ' name '_ExG.csv.'])
    end
    if (~exg_orn_marker(2))
         error(['---> Selected directory does not contain: ' name '_ORN.csv.'])
    end
    if (~exg_orn_marker(3))
         error(['---> Selected directory does not contain: ' name '_Marker.csv.'])
    end
end


function channelNameList = requestChannelLabels(EEG)
    prompt = cell(1, EEG.nbchan);
    definput = cell(1, EEG.nbchan);
    for i = 1:EEG.nbchan
        prompt(i) = {['Channel ' num2str(i) ':']};
        definput(i) = {EEG.chanlocs(1,i).labels};
    end

    dlgtitle = 'Channel Labels';
    dims = [1 12];
    channelNameList = inputdlg(prompt, dlgtitle, dims, definput, 'on');
end
