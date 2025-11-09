workspace; 
format long g
clc; clear
fontSize = 20;

I = imread("TM25_sk2.jpg");

lab_I = rgb2lab(I);
ab = lab_I(:, :, 2); % bere kanály a a b (2:3), zelená-červená, žlutá-modrá

ab = im2single(ab);
% ab_3(:,:,1:2) = ab;
% ab_3(:,:,3) = zeros(size(ab(:,:,1)));
% imshow([ab_3])

% ab = imgaussfilt(ab, 4); % nějak otestovat nejlepší hodnotu (,x) 

% ab = stdfilt(ab, ones(3));

% ab = fftshift(log(1+abs(ff2(ab)))); imshow(ab, [])

%%% gabor příprava
[numRows,numCols,~] = size(I);

wavelengthMin = 4/sqrt(2);
wavelengthMax = hypot(numRows,numCols);
n = floor(log2(wavelengthMax/wavelengthMin));
wavelength = 2.^(0:(n-2)) * wavelengthMin;

deltaTheta = 45;
orientation = 0:deltaTheta:(180-deltaTheta);

g = gabor([2, 6, 9, 12],orientation);
[g, ~] = imgaborfilt(ab,g);
g = max(g, [], 3);


%[L,C] = imsegkmeans(ab, 6, NumAttempts=10);
[L,C] = imsegkmeans([ab, g], 5, NumAttempts=10); % s gabor
imshow(L, [])

%imgaussfilt()

