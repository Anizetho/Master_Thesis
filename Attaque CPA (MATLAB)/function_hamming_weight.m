function [H] = function_hamming_weight(V)
    bin_hypothetic_inter2 = dec2bin(V);
    hamming_weight = sum(bin_hypothetic_inter2.' == '1');
    H = reshape(hamming_weight,size(V));
end

