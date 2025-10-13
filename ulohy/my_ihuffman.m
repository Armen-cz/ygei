
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
    values = reshape(values, [m^2, 1]);
end