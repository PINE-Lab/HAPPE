function output = select_average_triallevel(vals, s, seed, split)

if length(vals) >= s

    % Set the seed for reproducibility
    rng(seed);

    % Subsample data
    data_sub = randsample(vals,s,false);

    if ~split
        output = {{mean(data_sub)}};
        return
    end

    % Divide them into half 
    subtrials_1 = data_sub(1:ceil(length(data_sub)/2));
    subtrials_2 = data_sub(ceil(length(data_sub)/2) + 1: end);

    output = {{mean(subtrials_1), mean(subtrials_2)}};

else
    output = {{nan, nan}};
end

end
