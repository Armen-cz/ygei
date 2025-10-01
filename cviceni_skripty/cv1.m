clc 
clear
format long g

%Load

fig1= imread("C:\Users\adamk\Documents\.CVUT\7_semestr\geoinformatika\ygei\cviceni\cv1\Image1.bmp");
fig2= imread('C:\Users\adamk\Documents\.CVUT\7_semestr\geoinformatika\ygei\cviceni\cv1\Image2.bmp');
%imshow(fig1)


R = double(fig1(:,:,1)); %Double, aby se s tím dalo počítat
G = double(fig1(:,:,2));
B = double(fig1(:,:,3));

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

% Division to submatrices
[m, n] = size(Y);

% Process rows
for i = 1:8:(m-7)
    % Process columns
    for j = 1:8:(n-7)
        % Get submatrices
        Ys = Y(i:i+7, j:j+7);
        Cbs = Cb(i:i+7, j:j+7); 
        Crs = Cr(i:i+7, j:j+7);   

        % DCT vstup: raster(img-1b), výstup: raster(imgT-1b)
        % apply dct
        Ys_dct = dct(Ys);
        Cbs_dct = dct(Cbs);
        Crs_dct = dct(Crs);

        % compression factor - 1. krok vytvářející kompresi
        q = 50;
        Qyf = 50 * Qy / q;
        Qcf = 50 * Qc / q;

        % quantization
        Ys_q = round(Ys_dct ./ Qyf);
        Cb_q = round(Cbs_dct ./ Qcf);
        Cr_q = round(Crs_dct ./ Qcf);

        % update transformed matrix
        Y(i:i+7, j:j+7) = Ys_q;
        Cb(i:i+7, j:j+7) = Cb_q;
        Cr(i:i+7, j:j+7) = Cr_q;


    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JPEG DECOMPRESSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dequantize

% InvDCT

% YCbCr -> RGB

% Standard deviations for RGB components



function [img_t] = dct(img)

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
                    cos((2*(x+1)+1)*u*pi/16)*cos((2*(y+1)+1)*v*pi/16);

            end
        end
        % update raster
        img_t(u+1, v+1) = fuv;

    end

end
% end of function
end