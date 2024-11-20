function s = struct2str(struct_)
    s = 'struct(';
    fn = fieldnames(struct_);
    for k =1:numel(fn)
        el = struct_.(fn{k});
        s = strcat(s, sprintf("'%s', ", fn{k} )); %cat the field
        if isscalar(el)
            s = strcat(s, sprintf('%d, ', el));
        elseif isstring(el) | ischar(el)
            s = strcat(s, sprintf("'%s', ", el));
        elseif islogical(el)
            s = strcat(s, sprintf("%s, ", string(el)));
        else
            a2s = sprintf('%d ', el); %converting array to str
            s = strcat(s, sprintf('[%s], ', a2s));
        end

        if k == numel(fn)
            s = strip(s, 'right', ' '); %doing this because string slicing not working
            s = strip(s, 'right', ',');
            s = strcat(s, ')');
        end
    end
