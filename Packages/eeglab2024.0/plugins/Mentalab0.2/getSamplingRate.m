function sample_rate = getSamplingRate(exg_timestamp)

    i = 0;
    current_sec = 0;
    for n = 1:length(exg_timestamp)
        ts = exg_timestamp(n);
        if floor(ts) ~= current_sec
            if i > 249
                break;
            end
            current_sec = floor(ts);
            i = 0;
        end
        i = i + 1;
    end

    if i < 300 % arbitrary buffer against oversampling
        sample_rate = 250;
    elseif i < 600
        sample_rate = 500;
    else 
        sample_rate = 1000;
    end
end

