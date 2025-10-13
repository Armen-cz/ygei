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

% Division to submatrices
[m, n] = size(Y);

% compression factor - 1. krok vytvářející kompresi
q = 50;
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

        % % DCT vstup: raster(img-1b), výstup: raster(imgT-1b)
        % % apply dct
        % Ys_t = dct(Ys);
        % Cbs_t = dct(Cbs);
        % Crs_t = dct(Crs);

        %DFT vstup: raster(img-1b, výstup: raster(imgT-1b)
        Ys_t = real(dft(Ys));
        Cbs_t = real(dft(Cbs));
        Crs_t = real(dft(Crs));

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
    % Nejlépe nějaká fce
    
    %Huffman coding, z cik-cak rovnou do huffmana
    [Y_huffman, Y_codes] = my_huffman(Y);
    [Cb_huffman, Cb_codes] = my_huffman(Cb);
    [Cr_huffman, Cr_codes] = my_huffman(Cr);

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

end

% Process rows
for i = 1:8:(m-7)
    % Process columns
    for j = 1:8:(n-7)
        % Get submatrices
        Ys = Y(i:i+7, j:j+7);
        Cbs = Cb(i:i+7, j:j+7); 
        Crs = Cr(i:i+7, j:j+7);  

        % dequantization
        Ys_dq = Ys .* Qyf;
        Cbs_dq = Cbs .* Qcf;
        Crs_dq = Crs .* Qcf;

        % % IDCT vstup: raster(imgT-1b), výstup: raster(img-1b)
        % % apply idct
        % Ys_it = idct(Ys_dq);
        % Cbs_it = idct(Cbs_dq);
        % Crs_it = idct(Crs_dq);

        % IDFT vstup: raster(imgT-1b), výstup: raster(img-1b)
        Ys_it = real(idft(Ys_dq));
        Cbs_it = real(idft(Cbs_dq));
        Crs_it = real(idft(Crs_dq));

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

% compute standard deviations for rgb components
dR = R - R_new;
dG = G - G_new;
dB = B - B_new;

sigmaR = sqrt(sum(sum(dR.^2))/(m*n));
sigmaG = sqrt(sum(sum(dG.^2))/(m*n));
sigmaB = sqrt(sum(sum(dB.^2))/(m*n));

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

function [non_zero_sorted_values] = non_zero_sorted(values, size)
    non_zero_sorted_values = [];
    i = 1;
    for row = 1:length(values)
        if values(row, 2) ~= 0
            non_zero_sorted_values(i, :) = [values(row, 1), values(row,2)/size];
            i = i+1;
        end
    end
    non_zero_sorted_values = sortrows(non_zero_sorted_values, -2, "descend");
end


function [huffman_values, codes] = my_huffman(values)
    % input: matrix of values (for example 8x8)
    % output: matrix of cells of the same size as input with huffman coding
    [m, n] = size(values);
    
    % gets the count of all the values
    all_value_counts = tabulate(reshape(values,1,[]));

    % removes values with 0 occurence and sorts the list
    huffman_nodes = non_zero_sorted(all_value_counts, length(values)^2);
    huffman_nodes = num2cell(huffman_nodes); % better for nodes
    
    % builds the huffman tree
    while size(huffman_nodes,1) > 1
        % sorting to have the two with the lowest propability at the top
        huffman_nodes = sortrows(huffman_nodes, 2);
    
        % takes the smallest nodes
        left_node = huffman_nodes(1, :);
        right_node = huffman_nodes(2, :);
    
        % creates a new node
        new_symbol = {left_node{1}, right_node{1}}; % store children
        new_prob = left_node{2} + right_node{2};
        new_node = {new_symbol, new_prob};
    
        % removes old and creating a new one
        huffman_nodes(1:2, :) = [];
        huffman_nodes(end+1, :) = new_node;
    end
    
    % the final Huffman tree
    huffman_tree = huffman_nodes{1,1};
    
    % defining map for better variable storing
    queue = {huffman_tree, ''};
    codes = containers.Map('KeyType','double','ValueType','char');
    
    % creating huffman codes
    while length(queue) > 0
        % Pops the first element
        node = queue{1,1};
        prefix = queue{1,2};
        queue(1,:) = [];  % dequeue
    
        if iscell(node)
            % Left branch appends '0'
            queue(end+1,:) = {node{1}, strcat(prefix,'0')};
            % Right branch appends '1'
            queue(end+1,:) = {node{2}, strcat(prefix,'1')};
        else
            if iscell(node)
                node = node{1};
            end
            % Stores code in Map
            codes(node) = prefix;
        end
    end
    
    huffman_values = cell(m, n);
    
    for i = 1:m
        for j = 1:n
            codeStr = codes(values(i,j));
            huffman_values{i,j} = codeStr - '0'; % convert ascii to numbers
        end
    end
end

function [values] = my_ihuffman(huffman_values, codes)
    % input: matrix of huffman values (for example 8x8) and code map
    % output: matrix of cells of the same size as input with original
    % values
    [m, n] = size(huffman_values);
    values = zeros([m, n]);

    % builds inverse mapping: code string to symbol
    inv_codes = containers.Map('KeyType','char','ValueType','double');
    code_keys = codes.keys;
    for k = 1:length(code_keys)
        sym = code_keys{k};
        code = codes(sym);
        inv_codes(code) = sym;
    end

    for i = 1:m
        for j = 1:n
            code_bits = huffman_values{i,j};     
            code_str = char(code_bits + '0');  % converts numeric to ascii
            values(i,j) = inv_codes(code_str);
        end
    end
end
