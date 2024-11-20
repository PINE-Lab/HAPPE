function all_data = read_data(data_folder, subj_divider, ignore_contains, conditions, value_columns)

all_data = table();
files = dir(fullfile(data_folder,'**/*.csv'));
n_read = 0;
n_skip = 0;
parfor i = 1:length(files)
    baseFileName = lower(files(i).name);
    % fullFileName = fullfile(data_folder, baseFileName);
    fullFileName = fullfile(files(i).folder, files(i).name);
    split_arr = split(baseFileName, lower(subj_divider)); 
    subj = convertCharsToStrings(split_arr{1});

    if any(arrayfun(@(x) contains(lower(baseFileName), lower(x)), ignore_contains))
        % fprintf(1, 'Skipping %s\n', baseFileName);
        n_skip = n_skip + 1;
        continue
    end
    
    % fprintf(1, 'Now reading file: %s\n', baseFileName);
    n_read = n_read + 1;
    
    data = readtable(fullFileName,'PreserveVariableNames', true, 'Delimiter', ',');
    trial_mask = arrayfun(@(x) ~contains(lower(x{:}), "average"), data{:,'Row'});
    data = data(trial_mask,:);

    cond = "other";
    for j = 1:length(conditions)
        if contains(lower(baseFileName), lower(conditions(j)))
            cond = conditions(j);
            break
        end
    end

    cond = convertCharsToStrings(cond);

    n_trials = length(data{:,'Row'});

    sub_data = horzcat( ...
        array2table(repelem(subj, n_trials, 1), 'VariableNames', "id"), ...
        array2table(repelem(cond, n_trials, 1), 'VariableNames', "cond"), ...
        cell2table(data{:,'Row'}, 'VariableNames', "trial"));

    for j = 1:length(value_columns)
        vc = value_columns(j);
        sub_data = horzcat(sub_data, table(data{:,vc}, 'VariableNames', vc))
    end
    
    all_data = vertcat(all_data, sub_data);
end

fprintf(1, 'Finished reading data; loaded %d files, skipped %d files\n', n_read, n_skip);

end