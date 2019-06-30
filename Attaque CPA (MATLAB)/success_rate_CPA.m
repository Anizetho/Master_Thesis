%clear
%clc
%close all
mode = 32;
NUMBER_TRACES = 20000;
byte_to_attack = 14;

%% Load Traces
Nth_measurement = strcat(int2str(NUMBER_TRACES),'_traces');
formatSpec = 'In progress %i%% ... \n';
pourcent = 2;
fprintf(formatSpec,pourcent)
path = '..\Data\AES_256\Measurement_';
totalPath = strcat(path, Nth_measurement,'\');
traceFile = strcat(totalPath,'traces.txt');
plaintextFile = strcat(totalPath,'plaintexts.mat');
keyFile = strcat(totalPath,'keys_unique.mat');
ciphertextsFile = strcat(totalPath,'ciphertext.txt');
ciphertextsData = importdata(ciphertextsFile);
load(traceFile)
load(plaintextFile);
load(keyFile);

%%
% Fs = 1e9;
% Fc = 48e6;
% Wn = pi*Fc/Fs;
% N_order = 100;
% B = fir1(N_order, Wn, 'low', hamming(N_order+1));
% traces = filter(B,1,traces,[],2);
%%

n_total_traces  = size(traces,1);
nb_samples      = size(traces,2);
nb_sub_samples  = 400;

%% Correct data
row = 1;
ciphertexts = zeros(n_total_traces,16);
for i=1:n_total_traces
    for j=1:16
        if size(ciphertextsData{i},2) == 32
            TwoLetters = ciphertextsData{i}(((2*j)-1):(2*j));
            ciphertexts(i,j) = hex2dec(TwoLetters);
        else 
            error(row,j) = i;
            if j == 16
                row = row + 1;
            end
        end
    end            
end

sizeCal = n_total_traces - size(error,1);
newciphertexts = zeros(sizeCal, 16);
newplaintexts = zeros(sizeCal, 16);
newtraces = zeros(sizeCal, size(traces,2));
line=1;
for i=1:n_traces
    if ciphertexts(i,1) ~= 0 
        newciphertexts(line,:) = ciphertexts(i,:);
        newplaintexts(line,:) = plaintexts(i,:);
        newtraces(line,:) = traces(i,:);
        line = line+1;
    end
end

clear ciphertexts;
clear plaintexts;
clear traces;

ciphertexts = newciphertexts;
plaintexts = newplaintexts;
traces = newtraces;
n_total_traces  = size(traces,1);

%% Success Rate
n_traces = 8000;
success_rate = zeros(1, n_traces);


n_exp = 100;
for exp=1:n_exp
    fprintf('Iter : %d \n', exp)
    %% Select randomly subset of traces
    index1  = randperm(n_total_traces-1, n_traces)+1; 
    T       = traces(index1,1:nb_sub_samples);
    %% Select Key Byte to Attack
    %byte_to_attack  = 1;
    key_byte        = keys_unique(byte_to_attack);
    pt_byte         = plaintexts(index1, byte_to_attack);
    %% Calculate hypothetical intermediate values
    sbox_cipher = sbox(ciphertexts(index1 - 1, byte_to_attack));
    V = function_hypothetic_intermediate(pt_byte);
    dist =  bitxor(sbox_cipher,V);    
    bin_dist = dec2bin(dist); 
    Dist_hamming_desordonnate = sum(bin_dist.' == '1');
    Dist_hamming = reshape(Dist_hamming_desordonnate,[n_traces,256]);
    H = Dist_hamming;
    %% Correlation versus Traces
    R = corr(H,T);
    [~, ind] = max(abs(R(:)));
    [~, time_leakage] = ind2sub(size(R),ind);
    
%     fprintf('Real  Key : 0x%s  \n', dec2hex(key(byte_to_attack), 2));
%     fprintf('Found Key : 0x%s  \n', dec2hex(key_found, 2));
%     pause(1);
    tic
    for n = 3:n_traces  
        R = corr(H(1:n,:),T(1:n, time_leakage));
        [~, ind] = max(abs(R(:)));
        [key_found, ~] = ind2sub(size(R),ind);
        key_found = key_found - 1;
        if key_found == key_byte
           success_rate(n) = success_rate(n) + 1;
        end
    end
    toc
end

plot(1:n_traces, success_rate)
title('Success Rate of CPA vs Number of Traces')
xlabel('Number of Traces')
ylabel('Success Rate [%]')
grid on
hold on