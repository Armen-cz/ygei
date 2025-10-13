function [array] = zigzag(image)
    % transforms square matrix to an array in zigzag pattern
    % input: square matrix [m, m]
    % output: array [m^2, 1]
    [m, n] = size(image);
    i = 1;
    j = 1;
    array = zeros(m*n, 1);
    array_i = 1;
    veta = false; % rule veta, if one of the conditions for 'out_of_bounds' is
    % true, the next iteration has to be valid
    
    up = true;
    while array_i <= m*n
        if up % goes up right
            array(array_i) = image(i, j);
            array_i = array_i + 1;
    
            if (j>=n) & ~veta % if column index would be out of bounds
                up = false;
                i = i + 1;
                veta = true;
    
            elseif (i<=1) & ~veta % if row index would be out of bounds
                up = false;
                j = j + 1;
                veta = true;
    
            else
                i = i - 1; % if correct or rule veta
                j = j + 1;
                veta = false;
            end
    
        else % goes left down
            array(array_i) = image(i, j);
            array_i = array_i + 1;
    
            if (i>= m) & ~veta % if row index would be out of bounds
                up = true;
                j = j + 1;
                veta = true;
            
            elseif (j <= 1) & ~veta % if column index would be out of bounds
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

