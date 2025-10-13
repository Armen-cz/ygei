function [image] = inverse_zigzag(array)
    % transforms array into square matrix 
    % (array needs to have a whole number square root)
    % input: array [m^2, 1]
    % output: matrix [m, m]
    i = 1;
    j = 1;
    m = round(sqrt(length(array)));
    n = round(sqrt(length(array)));
    image = zeros(m, n);
    array_i = 1;
    veta = false; % rule veta, if one of the conditions for 'out_of_bounds' is
    % true, the next iteration has to be valid

    up = true;
    while array_i <= m*n
        if up % goes up right
            image(i, j) = array(array_i);
            array_i = array_i + 1;

            if (j >= n) && ~veta % if column index would be out of bounds
                up = false;
                i = i + 1;
                veta = true;

            elseif (i <= 1) && ~veta % if row index would be out of bounds
                up = false;
                j = j + 1;
                veta = true;

            else
                i = i - 1; % if correct or rule veta
                j = j + 1;
                veta = false;
            end

        else % goes left down
            image(i, j) = array(array_i);
            array_i = array_i + 1;

            if (i >= m) && ~veta % if row index would be out of bounds
                up = true;
                j = j + 1;
                veta = true;

            elseif (j <= 1) && ~veta % if column index would be out of bounds
                up = true;
                i = i + 1;
                veta = true;

            else
                j = j - 1; % if correct or rule veta
                i = i + 1;
                veta = false;
            end
        end
    end
end
