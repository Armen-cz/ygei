clc; clear; format long g

%Load
fig1= imread('stromy_barevne.png');
fig2= imread('fotografie.png');
fig3 = imread("vektor.png");

%imshow(fig1)

transform_type = "DFT"; % Choose DCT or DFT

R = double(fig3(:,:,1)); %Double, aby se s tím dalo počítat
G = double(fig3(:,:,2));
B = double(fig3(:,:,3));

% RGB to YCbCr
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

% Division to submatrices
[m, n] = size(Y);

% compression factor - 1. krok vytvářející kompresi
q = 70;
Qyf = 50 * Qy / q;
Qcf = 50 * Qc / q;

% Process rows
for i = 1:8:(m-7)
    % Process columns
    for j = 1:8:(n-7)
        % Get submatrices
        Ys = Y(i:i+7, j:j+7);
        Cbs = Cb(i:i+7, j:j+7); 
        Crs = Cr(i:i+7, j:j+7);  

        switch upper(transform_type) % Do only chosen transformation

            case "DCT"
            % DCT vstup: raster(img-1b), výstup: raster(imgT-1b)
            % apply dct
            Ys_t = dct(Ys);
            Cbs_t = dct(Cbs);
            Crs_t = dct(Crs);

            case "DFT"
            %DFT vstup: raster(img-1b, výstup: raster(imgT-1b)
            Ys_t = real(dft(Ys));
            Cbs_t = real(dft(Cbs));
            Crs_t = real(dft(Crs));

        end

        % quantization
        Ys_q = round(Ys_t ./ Qyf);
        Cb_q = round(Cbs_t ./ Qcf);
        Cr_q = round(Crs_t ./ Qcf);

        
        % update transformed matrix
        Y(i:i+7, j:j+7) = Ys_q;
        Cb(i:i+7, j:j+7) = Cb_q;
        Cr(i:i+7, j:j+7) = Cr_q;


    end
end

if m == n % only works with square pictures

    %%%% Tady bude cik-cak (zig-zag) sekvence
    Y_zig = zigzag(Y);
    Cb_zig = zigzag(Cb);
    Cr_zig = zigzag(Cr);
    
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
Y = my_ihuffman(Y_huffman, Y_codes);
Cb = my_ihuffman(Cb_huffman, Cb_codes);
Cr = my_ihuffman(Cr_huffman,Cr_codes);

% cik-cak inverse sekvence - udělat stále v tom if statement
Y_izig = inverse_zigzag(Y);
Cb_izig = inverse_zigzag(Cb);
Cr_izig = inverse_zigzag(Cr);
end

clear Y; clear Cb; clear Cr

% Process rows
for i = 1:8:(m-7)
    % Process columns
    for j = 1:8:(n-7)
        % Get submatrices
        Ys = Y_izig(i:i+7, j:j+7);
        Cbs = Cb_izig(i:i+7, j:j+7); 
        Crs = Cr_izig(i:i+7, j:j+7);  

        % dequantization
        Ys_dq = Ys .* Qyf;
        Cbs_dq = Cbs .* Qcf;
        Crs_dq = Crs .* Qcf;

        switch upper(transform_type) % does only chosen transformation

            case "DCT"
            % IDCT vstup: raster(imgT-1b), výstup: raster(img-1b)
            % apply idct
            Ys_it = idct(Ys_dq);
            Cbs_it = idct(Cbs_dq);
            Crs_it = idct(Crs_dq);
            
            case "DFT"
            % IDFT vstup: raster(imgT-1b), výstup: raster(img-1b)
            Ys_it = real(idft(Ys_dq));
            Cbs_it = real(idft(Cbs_dq));
            Crs_it = real(idft(Crs_dq));
        end

        % update transformed matrix
        Y(i:i+7, j:j+7) = Ys_it;
        Cb(i:i+7, j:j+7) = Cbs_it;
        Cr(i:i+7, j:j+7) = Crs_it;

    end
end

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
%imwrite(new_raster,"novy.png")

% compute standard deviations for rgb components
dR = R - R_new;
dG = G - G_new;
dB = B - B_new;

sigmaR = round(sqrt(sum(sum(dR.^2))/(m*n)),3);
sigmaG = round(sqrt(sum(sum(dG.^2))/(m*n)),3);
sigmaB = round(sqrt(sum(sum(dB.^2))/(m*n)),3);

%fprintf("$ %1d \\times %1d $ & %2d & %6.3f & %6.3f & %6.3f \\\\ \\hline\n", resample_size, resample_size, q, sigmaR, sigmaG, sigmaB)

function [img_t] = dct(img)
% discrete cosine transformation

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


function [img_t] = dft(img)
% discrete cosine transformation

[m, n] = size(img);
clear j; % to make sure 'j' is a complex unit
img_t = img;

% Process lines
for u = 0:7
% Process columns
    for v = 0:7
        % compute sum
        fuv = 0;
        % process lines
        for x = 0:7

            %process columns
            for y = 0:7
                fuv = fuv + img(x+1, y+1)*exp(-j*2*pi*((u*x/m) + (v*y/n)));
            end
        end
        % update raster
        img_t(u+1, v+1) = fuv;

    end

end
% end of function
end


function [img] = idft(img_t)
% inverse discrete cosine transform

[m, n] = size(img_t);
img = img_t;
% Process lines

% process lines
for x = 0:7

    %process columns
    for y = 0:7

        % compute sum
        fxy = 0;

        for u = 0:7

            for v = 0:7

                fxy = fxy + 1/(m*n) * img_t(u+1, v+1)*exp(j*2*pi*((u*x/m) + (v*y/n)));

            end
        end

        % update raster
        img(x+1, y+1) = fxy;

    end

end
% end of function
end