workspace; 
format long g
clc; clear
fontSize = 20;

% changes if the raster is filled/compact
fill =  true;

% reads an image using 3 RGB bands
I = imread("TM25_sk2.jpg");

% pixel coordinates of map corners (used for geotiff)
y_top = 579;
y_bottom = 4954;
x_right = 4570;
x_left = 395;

% crops the image to fit the map only
I = I(y_top:y_bottom, x_left:x_right, :);

% creates a new raster of zeros
only_green = zeros([size(I, 1:2)]);

% fills the raster with values
for y = 1:size(I,1)
    for x = 1:size(I,2)

        % contours
        if (I(y,x,1) > 182 && I(y,x,1) < 227) && ...
           (I(y,x,2) > 155 && I(y,x,2) < 203) && ...
           (I(y,x,3) > 107 && I(y,x,3) < 153)
            only_green(y, x) = 255;
        
        % dark forest
        elseif (I(y,x,1) > 222 && I(y,x,1) < 238) && ...
           (I(y,x,2) > 230 && I(y,x,2) < 250) && ...
           (I(y,x,3) > 160 && I(y,x,3) < 220)
            only_green(y, x) = 255;  
        
        % light forest
        elseif (I(y,x,1) > 240 && I(y,x,1) < 250) && ...
           (I(y,x,2) > 240 && I(y,x,2) < 253) && ...
           (I(y,x,3) > 200 && I(y,x,3) < 220)
            only_green(y, x) = 255;
        end
    end
end

%imshow(only_green, [])

if fill
    faded_gauss = imgaussfilt(only_green, 5); % used filter for fill
else
    faded_gauss = imgaussfilt(only_green, 2); % used filter for no fill
end
% imshow(faded_gauss, [])

% figure
faded_std = stdfilt(only_green); % did not really work for me
% imshow(faded_std, [])

average=fspecial('average',[9,9]); % exact 9x9 average filter, 
% similar but slightly worse results than gauss filter
faded_average=imfilter(only_green, average);

% imshow(round(faded_average), [])

% transfers double to single
rgb_fade = im2single(faded_gauss);
% imshow(rgb_fade, [])

% segmentation to 2 categories
[L,C] = imsegkmeans(rgb_fade, 2, NumAttempts=10);
% imshow(L, [])

% fills thin gaps and smooths smaller holes
if fill
    mask = imbinarize(L);
    se = strel('disk', 10);
    closed = imclose(mask, se);
else
    closed = imbinarize(L);
end

% fills all holes
filled = imfill(closed, 'holes');
% imshow(filled, [])

% removes small particles under 50 pixels
clean = bwareaopen(filled, 50);

figure
imshow(clean, [])
title('filtered image')


[row, col] = find(clean);
coords = [row, col];

% uložení do souboru
save('lesy.mat', 'coords');

out_file = "lesy_fill.tif";

imwrite(clean, out_file); % saves the file in tif -> this is then
                                                   % used in python script

