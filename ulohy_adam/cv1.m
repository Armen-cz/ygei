clc; clear; format long g

%Load

fig1= imread("C:\Users\adamk\Documents\.CVUT\7_semestr\geoinformatika\ygei\cviceni\cv1\Image1.bmp");
fig2= imread('C:\Users\adamk\Documents\.CVUT\7_semestr\geoinformatika\ygei\cviceni\cv1\Image2.bmp');
%imshow(fig1)


R = double(fig2(:,:,1)); %Double, aby se s tím dalo počítat
G = double(fig2(:,:,2));
B = double(fig2(:,:,3));

Y = 0.2990*R + 0.5870*G + 0.1140*B;
Cb = -0.1687*R - 0.3313*G + 0.5000*B + 128;
Cr = 0.5000*R - 0.4187*G - 0.0813*B + 128;

%Quantisation matrix
Qy = [16 11 10 16 24 40 51 61;
12 12 14 19 26 58 60 55;
14 13 16 24 40 87 69 56;
14 17 22 29 51 87 80 62;
18 22 37 26 68 109 103 77;
24 35 55 64 81 104 113 92;
49 64 78 87 103 121 120 101;
72 92 95 98 112 100 103 99];

% Chrominance matrix
Qc = [17 18 24 47 66 99 99 99
18 21 26 66 99 99 99 99
24 26 56 99 99 99 99 99
47 69 99 99 99 99 99 99
99 99 99 99 99 99 99 99
99 99 99 99 99 99 99 99
99 99 99 99 99 99 99 99
99 99 99 99 99 99 99 99];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JPEG COMPRESSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Store components
Y_old = Y;
Cb_old = Cb;
Cr_old = Cr;

[m_old, n_old] = size(Y_old);

% interval transfer
Y = 2*Y - 255;
Cb = 2*Cb - 255;
Cr = 2*Cr - 255;

% raster resampling
resample_size = 2; % used for resample and inverse resample

Y_res = resample(Y, resample_size);
Cb_res = resample(Cb, resample_size);
Cr_res = resample(Cr, resample_size);

% Division to submatrices
[m, n] = size(Y_res);

% compression factor - 1. krok vytvářející kompresi
q = 50;
Qyf = 50 * Qy / q;
Qcf = 50 * Qc / q;

% Process rows
for i = 1:8:(m-7)
    % Process columns
    for j = 1:8:(n-7)
        % Get submatrices
        Ys = Y_res(i:i+7, j:j+7);
        Cbs = Cb_res(i:i+7, j:j+7); 
        Crs = Cr_res(i:i+7, j:j+7);   

        % DCT vstup: raster(img-1b), výstup: raster(imgT-1b)
        % apply dct
        Ys_dct = dct(Ys);
        Cbs_dct = dct(Cbs);
        Crs_dct = dct(Crs);

        % quantization
        Ys_q = round(Ys_dct ./ Qyf);
        Cb_q = round(Cbs_dct ./ Qcf);
        Cr_q = round(Crs_dct ./ Qcf);

        % update transformed matrix
        Y_res(i:i+7, j:j+7) = Ys_q;
        Cb_res(i:i+7, j:j+7) = Cb_q;
        Cr_res(i:i+7, j:j+7) = Cr_q;


    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JPEG DECOMPRESSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Process rows
for i = 1:8:(m-7)
    % Process columns
    for j = 1:8:(n-7)
        % Get submatrices
        Ys = Y_res(i:i+7, j:j+7);
        Cbs = Cb_res(i:i+7, j:j+7); 
        Crs = Cr_res(i:i+7, j:j+7);  

        % dequantization
        Ys_dq = Ys .* Qyf;
        Cbs_dq = Cbs .* Qcf;
        Crs_dq = Crs .* Qcf;

        % IDCT vstup: raster(imgT-1b), výstup: raster(img-1b)
        % apply idct
        Ys_idct = idct(Ys_dq);
        Cbs_idct = idct(Cbs_dq);
        Crs_idct = idct(Crs_dq);

        % update transformed matrix
        Y_res(i:i+7, j:j+7) = Ys_idct;
        Cb_res(i:i+7, j:j+7) = Cbs_idct;
        Cr_res(i:i+7, j:j+7) = Crs_idct;

    end
end

% Resampling back to the original size 
% (to have same size for standard deviation calculations)
Y = iresample(Y_res, resample_size, m_old, n_old);
Cb = iresample(Cb_res, resample_size, m_old, n_old);
Cr = iresample(Cr_res, resample_size, m_old, n_old);

% Transform interval
Y = 1/2 * (Y + 255);
Cb = 1/2 * (Cb + 255);
Cr = 1/2 * (Cr + 255);

% YCbCr to RGB
R_new = Y + 1.4020*(Cr-128);
G_new = Y - 0.3441*(Cb-128) - 0.7141*(Cr-128);
B_new = Y + 1.7720*(Cb-128) - 0.0001*(Cr-128);

% double to uint8
Ri=uint8(R_new);
Gi=uint8(G_new);
Bi=uint8(B_new);

% create new raster
new_raster(:,:,1) = Ri;
new_raster(:,:,2) = Gi;
new_raster(:,:,3) = Bi;

imshow(new_raster)

% compute standard deviations for rgb components
dR = R - R_new;
dG = G - G_new;
dB = B - B_new;

sigmaR = sqrt(sum(sum(dR.^2))/(m*n));
sigmaG = sqrt(sum(sum(dG.^2))/(m*n));
sigmaB = sqrt(sum(sum(dB.^2))/(m*n));

function [img_t] = dct(img)
% discrete cosine transform

img_t = img;

% Process lines
for u = 0:7

    % compute cu
    if (u == 0)
        cu = sqrt(2)/2;
    else
        cu = 1;
    end
% Process columns
    for v = 0:7

        %compute cv
        if (v == 0)
            cv = sqrt(2)/2;
        else
            cv = 1;
        end
        
        % compute sum
        fuv = 0;
        % process lines
        for x = 0:7

            %process columns
            for y = 0:7

                fuv = fuv + 1/4 * cu *cv * img(x+1, y+1) * ...
                    cos((2*x+1)*u*pi/16)*cos((2*y+1)*v*pi/16);

            end
        end
        % update raster
        img_t(u+1, v+1) = fuv;

    end

end
% end of function
end


function [img] = idct(img_t)
% inverse discrete cosine transform

img = img_t;
% Process lines

% process lines
for x = 0:7

    %process columns
    for y = 0:7

        % compute sum
        fxy = 0;

        for u = 0:7
        
            % compute cu
            if (u == 0)
                cu = sqrt(2)/2;
            else
                cu = 1;
            end

            for v = 0:7
        
                %compute cv
                if (v == 0)
                    cv = sqrt(2)/2;
                else
                    cv = 1;
                end

                fxy=fxy+1/4*cu*cv*(img_t(u+1, v+1)*cos((2*x+1)*u*pi/16)*cos((2*y+1)*v*pi/16));

            end
        end

        % update raster
        img(x+1, y+1) = fxy;

    end

end
% end of function
end


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