function [all_eff_new, summary_eff] = calc_effectsize_overall( ...
    all_data, value_column, num_iterations, save_root)

fprintf(1, 'Now Calculating overall effect size ...\n');

G = findgroups(all_data(:,'id'));
% all_result_eff = table();

conds = all_data{:, 'cond'};
vals = all_data{:, value_column};
uq_conds = unique(conds);
uq_combs = nchoosek(uq_conds, 2);
size_comb = size(uq_combs);

% colnames = ["iteration", "sample", "seed", uq_conds'];
all_conds = uq_conds';
if size_comb(1) > 0
    for k=1:size_comb(1)
        all_conds = [all_conds, sprintf("%s_%s", uq_combs(k,1), uq_combs(k,2))];
    end
end


colnames = ["seed", all_conds];

all_eff_mtx = zeros(num_iterations, length(colnames));

parfor i = 1:num_iterations
    seed = randi(10000);
    eff_row = zeros(1, length(colnames));
    eff_row(find(colnames=="seed",1)) = seed;
    
    op_effect = @(vals, conds) get_vals_overall(vals, conds, seed, false);
    
    % Apply the operation on the data using splitapply
    comp_res = splitapply(op_effect, vals, conds, G);
    res = cell2mat(comp_res(:,1));
    eff_names = comp_res{1,2};
    
    for k=1:length(uq_conds)
        eff_idx = find(eff_names==uq_conds(k), 1);
        eff = meanEffectSize(res(:, eff_idx));
        col = string(uq_conds(k));
        idx = find(colnames==col,1);
        eff_row(idx) = eff.Effect;
        % tmp = table(eff.Effect, 'VariableNames', uq_conds(k));
        % eff_row = horzcat(eff_row, tmp);
    end
    size_comb = size(uq_combs);
    for k=1:size_comb(1)
        eff_idx1 = find(eff_names==uq_combs(k,1), 1);
        eff_idx2 = find(eff_names==uq_combs(k,2), 1);
        eff = meanEffectSize(res(:, eff_idx1), res(:, eff_idx2));
        col = sprintf("%s_%s", uq_combs(k,1), uq_combs(k,2));
        idx = find(colnames==col,1);
        eff_row(idx) = eff.Effect;
    end
        % all_result_eff = vertcat(all_result_eff, eff_row)
    
    all_eff_mtx(i,:) = eff_row;

end

all_eff_new = array2table(all_eff_mtx, 'VariableNames', colnames);


summ_cols = {'cond_s', 'mean', 'sd', 'upper_pct', 'lower_pct'};
summary_eff = table('Size',[length(colnames),length(summ_cols)],... 
	    'VariableNames', summ_cols, ...
        'VariableTypes', ["string", repelem("double", length(summ_cols)-1)]);
% 
% 
summary_eff.cond_s = colnames';
summary_eff.mean = mean(all_eff_new{:,:}, 1)';
summary_eff.sd = std(all_eff_new{:,:}, 1)';
summary_eff.upper_pct = prctile(all_eff_new{:,:}, 97.5, 1)';
summary_eff.lower_pct = prctile(all_eff_new{:,:}, 2.5, 1)';
summary_eff = summary_eff(summary_eff.cond_s~="seed",:);

save_folder = fullfile(save_root, value_column);
mkdir(save_folder);

full_filename = fullfile(save_folder, "eff_overall.csv");
summary_filename = fullfile(save_folder, "eff_summary_overall.csv");

writetable(all_eff_new, full_filename, 'Delimiter', ',')
writetable(summary_eff, summary_filename, 'Delimiter', ',')
fprintf(1, 'Overall effect size results saved in file: %s\n', full_filename);
fprintf(1, 'Overall effect size results summary saved in file: %s\n', summary_filename);
end


