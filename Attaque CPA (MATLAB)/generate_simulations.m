clear
clc

%% The user defines the number of traces
NUMBER_TRACES = 100;
generate_mask = 1;

%% To load the data
Nth_measurement = strcat(int2str(NUMBER_TRACES),'_traces');
path = '..\Data\AES_256\Simulations_';
mkdir(strcat(path, Nth_measurement));
file = fopen('..\Data\AES_256\FileForName.txt','w');
fprintf(file,'%s\n',Nth_measurement);
fclose(file);


%% Inputs
n_keys          = 1;
n_ppk           = NUMBER_TRACES;
n_plaintexts    = n_keys*n_ppk;

size_key        = 32;

%% Generate key
keyHexReal      = ["00";"01";"02";"03";"04";"05";"06";"07";"08";"09";"0a";"0b";"0c";"0d";"0e";"0f";"10";"11";"12";"13";"14";"15";"16";"17";"18";"19";"1a";"1b";"1c";"1d";"1e";"1f"];
keys_unique_Dec = hex2dec(keyHexReal);
keys_others     = randi([0 255], n_plaintexts-1, 32);
keys            = [keys_unique_Dec' ; keys_others];

% Keys in decimal
save(strcat(path, Nth_measurement,'\keysDec.mat'),'keys')

% Keys in hexadecimal
keys_rshp = reshape(keys.', n_plaintexts*size_key,1);
keys_hex1 = dec2hex(keys_rshp,2);
keys_hex  = reshape(keys_hex1.', n_plaintexts*2*size_key,1);
keys_hex  = reshape(keys_hex, 2*size_key, n_plaintexts).';
save(strcat(path, Nth_measurement,'\keysHex.mat'),'keys_hex')
fid = fopen(strcat(path, Nth_measurement,'\keysHex.txt'),'w');
for i =1:n_plaintexts
    fprintf(fid,keys_hex(i,:));
    fprintf(fid, '\n');
end
fclose(fid);


%% Generate mask
if generate_mask == 1
    mask_unique     = ["10";"20";"30";"40";"50";"60";"70";"80";"90";"a0";"b0";"c0";"d0";"e0";"f1";"01";"11";"21";"31";"41";"51";"61";"71";"81";"91";"a1";"b1";"c1";"d1";"e1";"f1";"01"];
    mask_unique_dec = hex2dec(mask_unique');
    masks           = repmat(mask_unique_dec, n_ppk, 1);
    idx             = randperm(n_plaintexts)';
    masks           = masks(idx,:);
    save(strcat(path, Nth_measurement,'\masks_unique.mat'),'mask_unique_dec');
        
    masks_rshp = reshape(masks.', n_plaintexts*size_key,1);
    masks_hex1 = dec2hex(masks_rshp,2);
    masks_hex  = reshape(masks_hex1.', n_plaintexts*2*size_key,1);
    masks_hex  = reshape(masks_hex, 2*size_key, n_plaintexts).';

    fid = fopen(strcat(path, Nth_measurement,'\masks.txt'),'w');
    for i =1:n_plaintexts
        fprintf(fid,masks_hex(i,:));
        fprintf(fid, '\n');
    end
    fclose(fid);

end

%% Generate plaintexts
plaintextHex = ["00";"11";"22";"33";"44";"55";"66";"77";"88";"99";"aa";"bb";"cc";"dd";"ee";"ff"];
plaintextDec = hex2dec(plaintextHex');
plaintexts = randi([0 255], n_plaintexts-1, 16);
plaintexts = [plaintextDec ; plaintexts];
save(strcat(path, Nth_measurement,'\plaintexts.mat'),'plaintexts')

plaintexts_rshp = reshape(plaintexts.', n_plaintexts*16,1);
plaintexts_hex1 = dec2hex(plaintexts_rshp,2);
plaintexts_hex  = reshape(plaintexts_hex1.', n_plaintexts*32,1);
plaintexts_hex  = reshape(plaintexts_hex, 32,n_plaintexts).';

fid = fopen(strcat(path, Nth_measurement,'\pt.txt'),'w');
for i =1:n_plaintexts 
    fprintf(fid,plaintexts_hex(i,:));
    fprintf(fid, '\n');
end
fclose(fid);



%% Computation results
fid = fopen(strcat(path, Nth_measurement,'\ciphertexts.txt'),'a');
for i=1:n_plaintexts
    % Plaintexts
    plaintextDec = plaintexts(i,:);
    [nb_traces, nb_bytes]=size(plaintextDec);
    
    % Keys
    keyDecReal = keys(i,:);
    keyDecMask = mask_unique_dec;
    keyDecFake = bitxor(keyDecReal, keyDecMask);
    % Key scheduling
    
    round_key_256 = zeros(8,32);
    round_key_256(1,:) = keyDecReal;
    for i = 1:7
        round_key_256(i+1,:) = key_schedule(round_key_256(i,:),i);
    end


    % Init Round
    % For test
    XOR = bitxor(plaintextDec, keyDecReal(1:16));

    % 13 rounds for AES-256
    for r=1:7
        if r==1
            SboxResult = sbox(XOR);
            ShiftRowsResult = shiftrows(SboxResult');
            %test = conversion(ShiftRowsResult)
            MixColumnResult = mixcolumns(ShiftRowsResult);
            NewKey = round_key_256(r,(17:32));
            XOR = bitxor(MixColumnResult,NewKey);
        else
            for k=1:2
                SboxResult = sbox(XOR);
                ShiftRowsResult = shiftrows(SboxResult');
                MixColumnResult = mixcolumns(ShiftRowsResult);
                if k==1
                    NewKey = round_key_256(r,(1:16));
                else
                    NewKey = round_key_256(r,(17:32));
                end
                XOR = bitxor(MixColumnResult,NewKey);
            end
        end
    end


    % Finish Round
    SboxResult = sbox(XOR);
    ShiftRowsResult = shiftrows(SboxResult');
    NewKey = round_key_256(8,(1:16));
    XOR = bitxor(ShiftRowsResult,NewKey);
    Ciphertext = XOR;
    CiphertextHex = dec2hex(Ciphertext,2);
    
    CiphertextHex  = reshape(CiphertextHex.', 32,1);
    CiphertextHex  = CiphertextHex.';
    fprintf(fid,CiphertextHex);
    fprintf(fid, '\n');

end
