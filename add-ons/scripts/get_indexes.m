%% get_indexes(values, list)
% returns a list of indexes of requested values in array
% gives index of closest element if value not in array

function indeces = get_indexes(values, list, range)
    if exist('range','var')
        if range
            start_end = get_indeces_from_list(values,list);
            indeces   = start_end(1):start_end(2);
        else
            indeces = get_indeces_from_list(values,list);
        end
    else
        indeces = get_indeces_from_list(values,list);
    end
end

function idxs = get_indeces_from_list(vals,lst)
    idxs = zeros(1, length(vals));
    for i = 1:length(vals)
        [~, idxs(i)] = min(abs(vals(1,i)-lst));
    end
end