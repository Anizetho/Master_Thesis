function [V] = function_hypothetic_intermediate(state_byte)
    %% Calculate hypothetical intermediate values
    n_traces = size(state_byte,1);
    keyz = 0:255;
    len_keyz = length(keyz);
    
    hypothetic_inter1 = zeros(n_traces,len_keyz);
    for k = 1:len_keyz
        hypothetic_inter1(:,k) = bitxor(state_byte,keyz(k));
    end
    
    V = sbox(hypothetic_inter1);
    
end