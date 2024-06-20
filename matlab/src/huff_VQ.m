% function to compute the Huffman code, returns the length of the huffman
% coded signal

% test the function 
%{
clear all, close all;
% made to take in compress code vectors

N = 1024;

% Example usage:
data_train = [1,3,4,6,6,6,3]; %randi([1, N], 1, 40000); % random normal num
data_test =  [3,3,4,6,6,10]; %randi([1, N], 1, 40000);
symbols = 1:N;
p = build_frequency_table(data_train, N);
dict = huffmandict(symbols, p);
code = huffmanenco(data_test, dict);

decodedSig = huffmandeco(code, dict);
isequal(data_test, decodedSig) % Should return true

lenOriginal = log2(N)*size(data_test,2);
lenHuff = size(code,2);
RCR = (1 - lenHuff/lenOriginal)*100;
CR = lenOriginal/lenHuff;

disp(['CR = ' num2str(CR)])
%}

function len_huff = huff_VQ(j_test, j_train, N)
    j_train_flat = round(reshape(j_train,[1,size(j_train,1)*size(j_train,2)])); 
    j_test_flat = round(reshape(j_test,[1,size(j_test,1)*size(j_test,2)])); 

    symbols = 1:N;
    p = build_frequency_table(j_train_flat, N);
    dict = huffmandict(symbols, p);
    code = huffmanenco(j_test_flat, dict);

    decoded = huffmandeco(code, dict);

    disp('is reconstruct equal: ');
    disp(isequal(j_test_flat, decoded));

    %lenOriginal = log2(N)*size(j_test_flat,2);
    len_huff = size(code,2);
    %cr = lenOriginal/lenHuff;
end


function probs = build_frequency_table(data, N)
    frequency_table = containers.Map('KeyType', 'double', 'ValueType', 'double');
    probs = zeros(N,1);
    for num = data
        if isKey(frequency_table, num)
            frequency_table(num) = frequency_table(num) + 1;
        else
            frequency_table(num) = 1;
        end
        probs(num) = frequency_table(num)/size(data,2);
    end
end

