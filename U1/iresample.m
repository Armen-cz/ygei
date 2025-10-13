function [new_img] = iresample(image, step, max_size_x, max_size_y)
% inverse resampling function (increases size of a raster)
% input: raster, step size (level of resampling),
%        max size to ensure the maximum size of a created raster
% output: resampled raster
    
    if nargin < 3
        max_size_x = Inf;
    end
    if nargin < 4
        max_size_y = Inf;
    end

    [m, n] = size(image);
    new_m = min(m * step, max_size_x);
    new_n = min(n * step, max_size_y);
    new_img = zeros(new_m, new_n);

    i = 1;
    for x = 1:m
        j = 1;
        for y = 1:n
            % assinging new values from existing raster

            new_img(i:min(i+(step-1), max_size_x), j:min(j+(step-1), max_size_y)) = image(x, y);

            j = j + step;
        end
        i = i + step;
    end
end
