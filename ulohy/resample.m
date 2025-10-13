function [new_img] = resample(image, step)
% resampling function
% input: raster, step size (level of resampling)
% output: resampled raster

    [m, n] = size(image);
    new_m = ceil(m / step); % rounded for indexing
    new_n = ceil(n / step); % rounded for indexing
    new_img = zeros(new_m, new_n);

    i = 1;
    for x = 1:step:m
        j = 1;
        for y = 1:step:n
            x_end = min(x + step - 1, m); % to not go over max index value
            y_end = min(y + step - 1, n);
            
            % calculating new value of pixel
            resampled_value = sum(sum(image(x:x_end, y:y_end))) / (step^2);
            new_img(i, j) = resampled_value;

            j = j + 1;
        end
        i = i + 1;
    end
end