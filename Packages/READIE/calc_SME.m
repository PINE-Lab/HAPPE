function result_SME = calc_SME(all_data, value_column, num_iterations, save_root)

fprintf(1, 'Now Calculating trial level standard measurement error (SME) ...\n');

result_SME = table();
uq_ids = unique(all_data{:,'id'});

for i = 1:length(uq_ids)
    tmp_data = all_data(all_data{:,'id'} == uq_ids(i),:);
    uq_conds = unique(tmp_data{:,'cond'});
    for j = 1:length(uq_conds)
        sub_data = tmp_data(tmp_data{:,'cond'}==uq_conds(j),:);

        mean_amplitudes = sub_data{:,value_column};
        n_trials = length(mean_amplitudes);
        aSME = std(mean_amplitudes) / sqrt(n_trials);

        bootstrap_sample_means = zeros(1, num_iterations);
        for k = 1:num_iterations
            bootstrap_indices = randi(length(mean_amplitudes), [1, length(mean_amplitudes)]);
            bootstrap_sample_means(k) = mean(mean_amplitudes(bootstrap_indices));
        end
        bSME = std(bootstrap_sample_means);

        row = cell2table({uq_ids(i), uq_conds(j), n_trials, mean(mean_amplitudes), std(mean_amplitudes), aSME, bSME});
        result_SME = vertcat(result_SME, row);

    end
end
result_SME.Properties.VariableNames = {'id', 'cond', 'n_trials', 'mean', 'sd', 'aSEM', 'bSEM'};

save_folder = fullfile(save_root, value_column);
mkdir(save_folder);
full_filename = fullfile(save_folder, "SME_triallevel.csv");
writetable(result_SME, full_filename, 'Delimiter', ',')

fprintf(1, 'Trial level SME results saved in file: %s\n', full_filename);

end