function [all_rel_new, summary_reliability] = calc_reliability_triallevel(all_data, value_column, values_sequence, num_iterations, save_root)

fprintf(1, 'Now Calculating trial level reliability ...\n');

G = findgroups(all_data(:,'id'));

conds = all_data{:, 'cond'};
vals = all_data{:, value_column};

uq_conds = unique(conds)';
colnames = get_colnames(uq_conds, values_sequence, false);
colnames = ["seed", colnames];

all_rel_mtx = zeros(num_iterations,length(colnames));

all_n = zeros(1,length(colnames));

% run the first iteration, get n values
corrs_row = zeros(1, length(colnames));
seed = randi(10000);
corrs_row(find(colnames=="seed",1)) = seed;
i=1;

for j = 1:length(values_sequence)
    s = values_sequence(j);
    op_reliability = @(vals, conds) get_vals_triallevel(vals, conds, uq_conds, s, seed, true);

    % Apply the operation on the data using splitapply
    comp_res = splitapply(op_reliability, vals, conds, G); % size(res) = n_id*(n_cond*2)
    % uq_conds = unique(conds);
    res = cell2mat(comp_res(:,1));
    corr_names = comp_res{1,2};

    corrs = corrcoef(res, 'Rows', 'pairwise'); % Running the correlations
    sb_corrs = (2*corrs)./(1+abs(corrs)); % Conducting the Spearman-Brown Formula
    
    for g = 1:length(uq_conds)
        col = string(uq_conds{g}) + '_' + string(s);
        idx = find(colnames==col,1);
        corr_idx = find(contains(corr_names, uq_conds{g}), 2);
        corrs_row(idx) = sb_corrs(corr_idx(1),corr_idx(2));
        all_n(idx) = sum(~isnan(res(:,corr_idx(1))));
    end
end
all_rel_mtx(i,:) = corrs_row;


parfor i = 2:num_iterations
    seed = randi(10000);
    corrs_row = zeros(1, length(colnames));
    corrs_row(find(colnames=="seed",1)) = seed;

    

    for j = 1:length(values_sequence)
        s = values_sequence(j);
        op_reliability = @(vals, conds) get_vals_triallevel(vals, conds, uq_conds, s, seed, true);

        % Apply the operation on the data using splitapply
        comp_res = splitapply(op_reliability, vals, conds, G); % size(res) = n_id*(n_cond*2)
        % uq_conds = unique(conds);
        res = cell2mat(comp_res(:,1));
        corr_names = comp_res{1,2};
    
        corrs = corrcoef(res, 'Rows', 'pairwise'); % Running the correlations
        sb_corrs = (2*corrs)./(1+abs(corrs)); % Conducting the Spearman-Brown Formula
        
        for g = 1:length(uq_conds)
            col = string(uq_conds{g}) + '_' + string(s);
            idx = find(colnames==col,1);
            corr_idx = find(contains(corr_names, uq_conds{g}), 2);
            corrs_row(idx) = sb_corrs(corr_idx(1),corr_idx(2));
        end
    end
    all_rel_mtx(i,:) = corrs_row;
end

all_rel_new = array2table(all_rel_mtx, 'VariableNames', colnames);

summ_cols = {'cond_s', 'mean', 'sd', 'upper_pct', 'lower_pct', 'n'};
summary_reliability = table('Size',[length(colnames),length(summ_cols)],... 
	    'VariableNames', summ_cols, ...
        'VariableTypes', ["string", repelem("double", length(summ_cols)-1)]);
% 
% 
summary_reliability.cond_s = colnames';
summary_reliability.mean = mean(all_rel_new{:,:}, 1)';
summary_reliability.sd = std(all_rel_new{:,:}, 1)';
summary_reliability.upper_pct = prctile(all_rel_new{:,:}, 97.5, 1)';
summary_reliability.lower_pct = prctile(all_rel_new{:,:}, 2.5, 1)';
summary_reliability.n = all_n';
summary_reliability = summary_reliability(summary_reliability.cond_s~="seed",:);

save_folder = fullfile(save_root, value_column);
mkdir(save_folder);

all_res_filename = fullfile(save_folder, "rel_triallevel.csv");
summary_filename = fullfile(save_folder, "rel_summary_triallevel.csv");

writetable(all_rel_new, all_res_filename, 'Delimiter', ',')
writetable(summary_reliability, summary_filename, 'Delimiter', ',')

fprintf(1, 'Trial level reliability results saved in file: %s\n', all_res_filename);
fprintf(1, 'Trial level reliability results summary saved in file: %s\n', summary_filename);

end
