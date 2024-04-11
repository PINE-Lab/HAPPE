function plot_from_summary(summary_table, mode, values_sequence, value_column, conditions, colors, save_root)
% all_conds = uq_conds';
mode = string(mode);
fprintf(1, 'Now plotting for %s ...\n', mode);

comb_conds = [];
if length(conditions) > 0
    uq_combs = nchoosek(conditions, 2);
    size_comb = size(uq_combs);
    for k=1:size_comb(1)
        comb_conds = [comb_conds, sprintf("%s_%s", uq_combs(k,1), uq_combs(k,2))];
    end
end

c = 1;
errors_upper = summary_table.upper_pct - summary_table.mean;
errors_lower = summary_table.mean - summary_table.lower_pct;
figure('visible','off'); 
for i=1:length(conditions)
    cond_mask = contains(summary_table.cond_s, conditions(i));
    for j=1:length(comb_conds)
         cond_mask = cond_mask & ~contains(summary_table.cond_s, comb_conds(j));
    end
    % cond_mask = arrayfun(@(x) contains(lower(x{:}), lower(conditions(i))), summary_table.cond_s);
    sub_summ = summary_table(cond_mask,:);
    errorbar(values_sequence, sub_summ.mean, errors_lower(cond_mask), errors_upper(cond_mask), ...
        'o-', 'LineWidth', 1, 'DisplayName', conditions(i), 'Color', colors(c, :));
    c = c + 1;
    hold on;
end
if lower(mode) ~= "reliability" && lower(mode) ~= "rel"
    for i=1:length(comb_conds)
        cond_mask = contains(summary_table.cond_s, comb_conds(i));
        % cond_mask = arrayfun(@(x) contains(lower(x{:}), lower(conditions(i))), summary_table.cond_s);
        sub_summ = summary_table(cond_mask,:);
        errorbar(values_sequence, sub_summ.mean, errors_lower(cond_mask), errors_upper(cond_mask), ...
            'o-', 'LineWidth', 1, 'DisplayName', comb_conds(i), 'Color', colors(c, :));
        c = c + 1;
        hold on;
    end
end

legend('Location','best')
hold off;

if ~isnan(save_root)
    img_path = fullfile(save_root, value_column, sprintf("rel_by_cond %s.tiff", mode));
    saveas(gcf, img_path);
    
    fprintf(1, 'Trial level %s plot saved as file: %s\n', mode, img_path);
else
    shg;
end



end