function [image] = inverse_zigzag(array)

    i = 1;
    j = 1;
    m = round(sqrt(length(array)));
    n = round(sqrt(length(array)));
    image = zeros(m, n);
    array_i = 1;
    veta = false; % same rule: skip one movement after border collision
    up = true;

    while array_i <= m*n
        if up % goes up right
            image(i, j) = array(array_i);
            array_i = array_i + 1;

            if (j >= n) && ~veta % right border
                up = false;
                i = i + 1;
                veta = true;

            elseif (i <= 1) && ~veta % top border
                up = false;
                j = j + 1;
                veta = true;

            else
                i = i - 1;
                j = j + 1;
                veta = false;
            end

        else % goes left down
            image(i, j) = array(array_i);
            array_i = array_i + 1;

            if (i >= m) && ~veta % bottom border
                up = true;
                j = j + 1;
                veta = true;

            elseif (j <= 1) && ~veta % left border
                up = true;
                i = i + 1;
                veta = true;

            else
                j = j - 1;
                i = i + 1;
                veta = false;
            end
        end
    end
end
