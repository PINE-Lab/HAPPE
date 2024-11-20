
function output = get_vals_triallevel(vals, conds, uq_conds, s, seed, split)

groups = unique(conds);
G = findgroups(conds);

res = splitapply(@(vals)select_average_triallevel(vals, s, seed, split), vals, G);

if split
    bins = [1,2];
    colnames = get_colnames([1,2], uq_conds, true);
    vals_row = NaN(1, length(colnames));
    for i = 1:length(groups)
        for j = 1: length(bins)
            col = string(groups{i}) + '_' + string(bins(j));
            idx = find(colnames==col,1);
            vals_row(idx) = res{i}{bins(j)};
        end
    end
else
    colnames = uq_conds;
    vals_row = NaN(1, length(colnames));
    for i = 1:length(groups)
        idx = find(colnames==groups(i),1);
        vals_row(idx) = res{i}{1};
    end
end

output = {vals_row, colnames};

end 
