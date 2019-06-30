clear
clc

%% The user defines the number of traces
NUMBER_TRACES = 52000;
generate_mask = 1;

%% To load the data
Nth_measurement = strcat(int2str(NUMBER_TRACES),'_traces');
path = '..\Data\AES_256\Measurement_';
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
%keys_unique = randi([0 255], n_keys, 32);
keyHexReal  = ["00";"01";"02";"03";"04";"05";"06";"07";"08";"09";"0a";"0b";"0c";"0d";"0e";"0f";"10";"11";"12";"13";"14";"15";"16";"17";"18";"19";"1a";"1b";"1c";"1d";"1e";"1f"];
keys_unique = hex2dec(keyHexReal);
keys        = repmat(keys_unique', n_ppk, 1);
labels      = repmat((1:n_keys)',n_ppk,1);
idx         = randperm(n_plaintexts)';
keys        = keys(idx,:);
labels      = labels(idx);
save(strcat(path, Nth_measurement,'\labels.mat'),'labels');

keys_rshp = reshape(keys.', n_plaintexts*size_key,1);
keys_hex1 = dec2hex(keys_rshp,2);
keys_hex  = reshape(keys_hex1.', n_plaintexts*2*size_key,1);
keys_hex  = reshape(keys_hex, 2*size_key, n_plaintexts).';

fid = fopen(strcat(path, Nth_measurement,'\keys.txt'),'w');
for i =1:n_plaintexts
    fprintf(fid,keys_hex(i,:));
    fprintf(fid, '\n');
end
fclose(fid);

keys_hex_fpga = "00" + keys_hex(:,1:32) + "01" + keys_hex(:,33:64);

fid = fopen(strcat(path, Nth_measurement,'\keys_fpga.txt'),'w');
for i =1:n_plaintexts
    fprintf(fid,keys_hex_fpga(i,:));
    fprintf(fid, '\n');
end
fclose(fid);

save(strcat(path, Nth_measurement,'\keys_unique.mat'),'keys_unique');

keys_unique_rshp = reshape(keys_unique.', n_keys*size_key,1);
key_unique_hex = dec2hex(keys_unique_rshp,2);
key_unique_hex  = reshape(key_unique_hex.', n_keys*2*size_key,1);
key_unique_hex  = reshape(key_unique_hex, 2*size_key, n_keys).';

fid = fopen(strcat(path, Nth_measurement,'\keys_unique.txt'),'w');
for i =1:n_keys
    fprintf(fid,key_unique_hex(i,:));
    fprintf(fid, '\n');
end
fclose(fid);

%% Generate mask
if generate_mask == 1
    %mask_unique = randi([0 255], n_keys, 32);
    mask_unique     = ["10";"20";"30";"40";"50";"60";"70";"80";"90";"a0";"b0";"c0";"d0";"e0";"f1";"01";"11";"21";"31";"41";"51";"61";"71";"81";"91";"a1";"b1";"c1";"d1";"e1";"f1";"01"];
    mask_unique_dec = hex2dec(mask_unique');
    masks           = repmat(mask_unique_dec, n_ppk, 1);
    idx             = randperm(n_plaintexts)';
    masks           = masks(idx,:);

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

    masks_hex_fpga = "02" + masks_hex(:,1:32) + "03" + masks_hex(:,33:64);

    fid = fopen(strcat(path, Nth_measurement,'\masks_fpga.txt'),'w');
    for i =1:n_plaintexts
        fprintf(fid,masks_hex_fpga(i,:));
        fprintf(fid, '\n');
    end
    fclose(fid);

    save(strcat(path, Nth_measurement,'\masks_unique.mat'),'mask_unique_dec');

    masks_unique_rshp = reshape(mask_unique_dec.', n_keys*size_key,1);
    mask_unique_hex = dec2hex(masks_unique_rshp,2);
    mask_unique_hex  = reshape(mask_unique_hex.', n_keys*2*size_key,1);
    mask_unique_hex  = reshape(mask_unique_hex, 2*size_key, n_keys).';

    fid = fopen(strcat(path, Nth_measurement,'\masks_unique.txt'),'w');
    for i =1:n_keys
        fprintf(fid,mask_unique_hex(i,:));
        fprintf(fid, '\n');
    end
    fclose(fid);
    header_plaintexts_fpga = "04";
    
elseif generate_mask == 0
    header_plaintexts_fpga = "03";

end

%% Generate plaintexts
plaintextHex = ["00";"11";"22";"33";"44";"55";"66";"77";"88";"99";"aa";"bb";"cc";"dd";"ee";"ff"];
plaintextDec = hex2dec(plaintextHex');
plaintexts = randi([0 255], n_plaintexts-1, 16);
plaintexts = [plaintextDec ; plaintexts];
save(strcat(path, Nth_measurement,'\plaintexts.mat'),'plaintexts')

plaintexts = reshape(plaintexts.', n_plaintexts*16,1);
plaintexts_hex1 = dec2hex(plaintexts,2);
plaintexts_hex = reshape(plaintexts_hex1.', n_plaintexts*32,1);
plaintexts_hex = reshape(plaintexts_hex, 32,n_plaintexts).';

plaintexts_fpga = header_plaintexts_fpga + plaintexts_hex;

fid = fopen(strcat(path, Nth_measurement,'\pt.txt'),'w');
for i =1:n_plaintexts 
    fprintf(fid,plaintexts_hex(i,:));
    fprintf(fid, '\n');
end
fclose(fid);

fid = fopen(strcat(path, Nth_measurement,'\pt_fpga.txt'),'w');
for i =1:n_plaintexts 
    fprintf(fid,plaintexts_fpga(i,:));
    fprintf(fid, '\n');
end
fclose(fid);