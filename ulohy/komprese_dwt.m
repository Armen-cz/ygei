clc; clear; format long g

%Load

fig1= imread('Image1.bmp');
fig2= imread('Image2.bmp');
figShrek= imread('shrek.png');
%imshow(fig1)


R = double(figShrek(:,:,1)); %Double, aby se s tím dalo počítat
G = double(figShrek(:,:,2));
B = double(figShrek(:,:,3));

[LL, LH, HL, HH] = dwt2d(R);
R_mega_new = idwt2d(LL, LH, HL, HH);    

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
resample_size = 1; % used for resample and inverse resample

Y_res = resample(Y, resample_size);
Cb_res = resample(Cb, resample_size);
Cr_res = resample(Cr, resample_size);

% Resampling back to the original size 
% (to have same size for standard deviation calculations) and for better
% results from zig-zag sequence and huffman coding
Y = iresample(Y_res, resample_size, m_old, n_old);
Cb = iresample(Cb_res, resample_size, m_old, n_old);
Cr = iresample(Cr_res, resample_size, m_old, n_old);

% Direct wavelet transform
[Y_LL, Y_LH, Y_HL, Y_HH] = dwt2d(Y);
[Cb_LL, Cb_LH, Cb_HL, Cb_HH] = dwt2d(Cb);
[Cr_LL, Cr_LH, Cr_HL, Cr_HH] = dwt2d(Cr);

% Division to submatrices
[m, n] = size(Y_LL);

% compression factor - 1. krok vytvářející kompresi
q = 10;
Qyf = 50 * Qy / q;
Qcf = 50 * Qc / q;

% Process rows
for i = 1:8:(m-7)
    % Process columns
    for j = 1:8:(n-7)
        % Get submatrices
        Ys = Y_LL(i:i+7, j:j+7);
        Cbs = Cb_LL(i:i+7, j:j+7); 
        Crs = Cr_LL(i:i+7, j:j+7);   

        % quantization
        Ys_q = round(Ys ./ Qyf);
        Cb_q = round(Cbs ./ Qcf);
        Cr_q = round(Crs ./ Qcf);

        
        % update transformed matrix
        Y_LL(i:i+7, j:j+7) = Ys_q;
        Cb_LL(i:i+7, j:j+7) = Cb_q;
        Cr_LL(i:i+7, j:j+7) = Cr_q;


    end
end

if m == n % only works with square pictures

    %%%% Tady bude cik-cak (zig-zag) sekvence
    Y_zig = zigzag(Y_LL);
    Cb_zig = zigzag(Cb_LL);
    Cr_zig = zigzag(Cr_LL);

    %Huffman coding, z cik-cak rovnou do huffmana
    [Y_huffman, Y_codes] = my_huffman(Y_zig);
    [Cb_huffman, Cb_codes] = my_huffman(Cb_zig);
    [Cr_huffman, Cr_codes] = my_huffman(Cr_zig);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JPEG DECOMPRESSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Huffman decoding
if m == n
Y_LL = my_ihuffman(Y_huffman, Y_codes);
Cb_LL = my_ihuffman(Cb_huffman, Cb_codes);
Cr_LL = my_ihuffman(Cr_huffman,Cr_codes);

% inverse zig-zag back to m x m matrix
Y_zig = inverse_zigzag(Y_LL);
Cb_zig = inverse_zigzag(Cb_LL);
Cr_zig = inverse_zigzag(Cr_LL);
end

clear Y_LL; clear Cb_LL; clear Cr_LL

% Process rows
for i = 1:8:(m-7)
    % Process columns
    for j = 1:8:(n-7)
        % Get submatrices
        Ys = Y_zig(i:i+7, j:j+7);
        Cbs = Cb_zig(i:i+7, j:j+7); 
        Crs = Cr_zig(i:i+7, j:j+7);  

        % dequantization
        Ys_dq = Ys .* Qyf;
        Cbs_dq = Cbs .* Qcf;
        Crs_dq = Crs .* Qcf;

        % update transformed matrix
        Y_LL(i:i+7, j:j+7) = Ys_dq;
        Cb_LL(i:i+7, j:j+7) = Cbs_dq;
        Cr_LL(i:i+7, j:j+7) = Crs_dq;

    end
end

% inverse direct wavelet transform
Y = idwt2d(Y_LL, Y_LH, Y_HL, Y_HH); 
Cb = idwt2d(Cb_LL, Cb_LH, Cb_HL, Cb_HH); 
Cr = idwt2d(Cr_LL, Cr_LH, Cr_HL, Cr_HH); 

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


function [LL, LH, HL, HH] = dwt2d(img)
    % processes columns (n columns -> n/2)
    [m, n] = size(img);
    L = zeros(m, n/2);
    H = zeros(m, n/2);
    
    % creates low pass and high pass filter
    for i = 1:m
        k = 1;
        for j = 1:2:n
            L(i, k) = (img(i,j) + img(i,j+1)) / sqrt(2);   % lowpass
            H(i, k) = (img(i,j) - img(i,j+1)) / sqrt(2);   % highpass
            k = k + 1;
        end
    end

    % processes rows creates LL, LH, HL and HH filters
    LL = zeros(m/2, n/2);
    LH = zeros(m/2, n/2);
    HL = zeros(m/2, n/2);
    HH = zeros(m/2, n/2);
    
    for j = 1:n/2
        k = 1;
        for i = 1:2:m
            % lowpass rows
            LL(k,j) = (L(i,j) + L(i+1,j)) / sqrt(2);  
            LH(k,j) = (L(i,j) - L(i+1,j)) / sqrt(2); 
    
            % highpass rows
            HL(k,j) = (H(i,j) + H(i+1,j)) / sqrt(2); 
            HH(k,j) = (H(i,j) - H(i+1,j)) / sqrt(2); 
    
            k = k + 1;
        end
    end
end

function [img] = idwt2d(LL, LH, HL, HH)
    % Inverse 2D Haar DWT
    [m_LL, n_LL] = size(LL);
    m = m_LL * 2;
    n = n_LL * 2;

    L = zeros(m, n/2);
    H = zeros(m, n/2);
    img = zeros(m, n);
    
    for j = 1:n/2
        k = 1;
        for i = 1:m/2
            % Recombine columns
            L(2*i-1,j) = (LL(i,j) + LH(i,j)) / sqrt(2);
            L(2*i,  j) = (LL(i,j) - LH(i,j)) / sqrt(2);
    
            H(2*i-1,j) = (HL(i,j) + HH(i,j)) / sqrt(2);
            H(2*i,  j) = (HL(i,j) - HH(i,j)) / sqrt(2);
            k = k + 1;
        end
    end
    
    for i = 1:m
        k = 1;
        for j = 1:n/2
            img(i,2*j-1) = (L(i,j) + H(i,j)) / sqrt(2);
            img(i,2*j)   = (L(i,j) - H(i,j)) / sqrt(2);
            k = k + 1;
        end
    end
end
