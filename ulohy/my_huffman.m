function [huffman_values, codes] = my_huffman(values)
    % input: matrix of values (for example 8x8)
    % output: matrix of cells of the same size as input with huffman coding
    %         codes map

    % reshapes values from zigzag function (array)
    values = reshape(values, [round(sqrt(length(values))), round(sqrt(length(values)))]);
    
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
            % left branch appends '0'
            queue(end+1,:) = {node{1}, strcat(prefix,'0')};
            % right branch appends '1'
            queue(end+1,:) = {node{2}, strcat(prefix,'1')};
        else
            if iscell(node)
                node = node{1};
            end
            % stores code in Map
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