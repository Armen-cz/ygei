workspace; 
format long g
clc; clear
fontSize = 20;

I = imread("TM25_sk2.jpg");

lab_I = rgb2lab(I);
ab = lab_I(:, :, 2:3); % bere kanály a a b (2:3), zelená-červená, žlutá-modrá

only_green = zeros([size(I, 1:2)]);


for y = 1:size(I,1)
    for x = 1:size(I,2)
        % vrstevnice na zelenou
        if (I(y,x,1) > 182 && I(y,x,1) < 227) && ...
           (I(y,x,2) > 155 && I(y,x,2) < 203) && ...
           (I(y,x,3) > 107 && I(y,x,3) < 153)
            only_green(y, x) = 255;
        elseif (I(y,x,1) > 222 && I(y,x,1) < 238) && ...
           (I(y,x,2) > 230 && I(y,x,2) < 250) && ...
           (I(y,x,3) > 160 && I(y,x,3) < 220)
            only_green(y, x) = 255;   
        elseif (I(y,x,1) > 240 && I(y,x,1) < 250) && ...
           (I(y,x,2) > 240 && I(y,x,2) < 253) && ...
           (I(y,x,3) > 200 && I(y,x,3) < 220)
            only_green(y, x) = 255;
        end
    end
end

ab = im2single(ab);
% ab_3(:,:,1:2) = ab;
% ab_3(:,:,3) = zeros(size(ab(:,:,1)));
% figure;
% subplot(1,3,1)
% imshow(I(1100:1500,1100:1500,1), [])
% title('R channel (green–red)')
% 
% subplot(1,3,2)
% imshow(I(1100:1500,1100:1500,2), [])
% title('G channel (blue–yellow)')
% 
% subplot(1,3,3)
% imshow(I(1100:1500,1100:1500,3), [])
% title('B channel (blue–yellow)')

figure
% imshow(imresize(only_green, 0.5))

rgb_fade = imgaussfilt(only_green, 5); % nějak otestovat nejlepší hodnotu (,x) 
% imshow(rgb_fade, [])
% ab = stdfilt(ab, ones(3));
% rgb_fade = stdfilt(only_green);
% ab = fftshift(log(1+abs(ff2(ab)))); imshow(ab, [])

%%% gabor příprava
% [numRows,numCols,~] = size(I);
% 
% wavelengthMin = 4/sqrt(2);
% wavelengthMax = hypot(numRows,numCols);
% n = floor(log2(wavelengthMax/wavelengthMin));
% wavelength = 2.^(0:(n-2)) * wavelengthMin;
% 
% deltaTheta = 45;
% orientation = 0:deltaTheta:(180-deltaTheta);
% 
% g = gabor([2, 6, 9, 12],orientation);
% [g, ~] = imgaborfilt(ab,g);
% g = max(g, [], 3);

rgb_fade = im2single(rgb_fade);
[L,C] = imsegkmeans(rgb_fade, 2, NumAttempts=10);
% [L,C] = imsegkmeans([ab, g], 5, NumAttempts=10); % s gabor

mask = imbinarize(L);
se = strel('disk', 10);

closed = imclose(mask, se);
filled = imfill(closed, 'holes');
clean = bwareaopen(filled, 50);

imshow(clean)
title('Filled and cleaned mask')


% imshow(L, [])

%imgaussfilt()

figure;
% subplot(1,3,1)
% imshow(I(1100:1500,1100:1500,1), [])
% title('R channel (green–red)')
% 
% subplot(1,3,2)
% imshow(I(1100:1500,1100:1500,2), [])
% title('G channel (blue–yellow)')