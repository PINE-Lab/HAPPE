function output = select_average_overall(vals, seed, split)

    rng(seed);

    if ~split
        % add a simple bootstrap too the mean function, otherwise values
        % will all be the same
        output = {{mean(randsample(vals, length(vals), true))}};
        return
    end

    % Subsample data
    n = numel(vals); % Calculate the total number of elements in 'vals'
    half_size = ceil(n / 2); % Compute the size of half the 'vals' array, rounding up if needed

    % Generate a random sample (subset) of half_size elements from 'vals'
    subtrials_1 = randsample(vals, half_size, false);
    % get the other half 
    subtrials_2 = setdiff(vals, subtrials_1);

    output = {{mean(subtrials_1), mean(subtrials_2)}};
%     else
%         output = {{nan, nan}};
%     end
end

