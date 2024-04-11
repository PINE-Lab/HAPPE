function colnames = get_colnames(arr1, arr2, reverse)

colnames = strings(1, length(arr1)*length(arr2));

for i = 1:length(arr1)
    for j = 1:length(arr2)
        if reverse
            colnames(length(arr1)*(j-1)+i) = string(arr2(j)) + "_" + string(arr1(i));
        else
            colnames(length(arr1)*(j-1)+i) = string(arr1(i)) + "_" + string(arr2(j));
        end
    end
end

end