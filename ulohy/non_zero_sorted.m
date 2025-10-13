function [non_zero_sorted_values] = non_zero_sorted(values, size)
    non_zero_sorted_values = [];
    i = 1;
    for row = 1:length(values)
        if values(row, 2) ~= 0
            non_zero_sorted_values(i, :) = [values(row, 1), values(row,2)/size];
            i = i+1;
        end
    end
    non_zero_sorted_values = sortrows(non_zero_sorted_values, -2, "descend");
end

