%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         CPA attack (AES-256)   HW
%                               with Filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

%% Load Traces
%% Configuration :
% - The user defines the AES's mode (128=16 ; 256=32)
% - The user defines the number of traces that he would like test
mode = 32;
NUMBER_TRACES = 41000;
Enable_Faking_Implementation = 0;
sampling_interval_ns = 2;

%% Load the data
Nth_measurement = strcat(int2str(NUMBER_TRACES),'_traces');
formatSpec = 'In progress %i%% ... \n';
pourcent = 2;
fprintf(formatSpec,pourcent)
path = '..\Data\AES_256\Measurement_';
totalPath = strcat(path, Nth_measurement,'\');
traceFile = strcat(totalPath,'traces.txt');
traces = load(traceFile);
plaintextFile = strcat(totalPath,'plaintexts.mat');
load(plaintextFile);
keyFile = strcat(totalPath,'keys_unique.mat');
load(keyFile);

plaintexts = plaintexts([7000:17000],:);
traces = traces([7000:17000],:);

nb_traces  = size(traces,1);
nb_samples = size(traces,2);
times = (0:1:nb_samples-1)*(1/(1e9/2)); % en seconde
times = 1e9*times;                      % en ns



%% Faking key
if Enable_Faking_Implementation == 1
    maskFile = strcat(totalPath,'masks_unique.mat');
    load(maskFile);
    Fake_Key = bitxor(keys_unique',mask_unique_dec);
    Fake_Key_hex = conversion(Fake_Key);
end

%% Digital Filter
Fs = 1e9/2;
Fc = 48e6;
Wn = pi*Fc/Fs;
N_order = 100; % impacte : délai en temporel
B = fir1(N_order, Wn, 'low', hamming(N_order+1));
raw_traces_filt = filter(B,1,traces,[],2);

%% Discrete Fourier Transform
N  = nb_samples;
fd  = -Fs/2:Fs/N:Fs/2-Fs/N; % en Hz
fd = fd/1e6;                % en MHz

fft_trace = fft(traces(1000,:),[],2);
fft_trace_filt = fft(raw_traces_filt(1000,:),[],2);

%% Parameters plot
fontSizeLabel  = 16;
fontSizeAxis   = 14;
sampling_delete = 65;
TimeToPlot = times(sampling_delete:size(times,2));
TimeToPlot = TimeToPlot - TimeToPlot(1,1);

%% Plot
figure
subplot(2,1,1)
plot(TimeToPlot, traces(1000,1:size(times,2)-sampling_delete+1));
title('Trace brute', 'fontsize', fontSizeLabel)
xlabel('Temps [ns]', 'fontsize', fontSizeLabel)
ylabel('Tension [mV]', 'fontsize', fontSizeLabel)
grid on;
ylim([-110 70]);
set(gca,'FontSize', fontSizeAxis)

subplot(2,1,2)
plot(TimeToPlot, raw_traces_filt(1000,sampling_delete:size(times,2)));
title('Trace filtrée', 'fontsize', fontSizeLabel)
xlabel('Temps [ns]', 'fontsize', fontSizeLabel)
ylabel('Tension [mV]', 'fontsize', fontSizeLabel)
grid on;
ylim([-110 70]);
set(gca,'FontSize', fontSizeAxis)


figure
plot(fd, abs(fftshift(fft_trace)))
title('FFT - Trace brute')
xlabel('Fréquence [MHz]')
ylabel('Amplitude [/]')
%ylim([0 5500]);
grid on
set(gca)

figure
plot(fd, abs(fftshift(fft_trace_filt)))
title('FFT - Trace filtrée')
xlabel('Fréquence [MHz]')
ylabel('Amplitude [/]')
ylim([0 5500]);
grid on
set(gca,'FontSize', fontSizeAxis)


%% Simulate the algorithm
keys = [0:255];
len_keys = length(keys');
[nb_traces,nb_sampling] = size(traces);
nb_columns = (mode/2)*len_keys;
hypothetic_inter1 = zeros(nb_traces,nb_columns);
for byte=1:16
    hypothetic_inter1(:,[((256*mode)-(256*((mode+1)-byte)))+1:((256*mode)-(256*(mode-byte)))]) = bitxor(plaintexts(:,byte),keys);
end
hypothetic_inter2 = sbox(hypothetic_inter1);
bin_hypothetic_inter2 = dec2bin(hypothetic_inter2);
Weight_Hamm_desordonnate = sum(bin_hypothetic_inter2.' == '1');
Weight_Hamm = reshape(Weight_Hamm_desordonnate,[nb_traces,size(hypothetic_inter2,2)]);
clc;
pourcent = 10;
fprintf(formatSpec,pourcent)
%toc;

%% Calculate the correlation coefficient
coefficient = corr(Weight_Hamm, raw_traces_filt);
Result = coefficient';
clc;
pourcent = 30;
fprintf(formatSpec,pourcent)
%toc;
        
%% Key values
Key_values = zeros(1,mode/2);
for i=1:16
    n_1 = (i*256)-255;
    n=i*256;
    [maximum,indice] = max(max(abs(Result(:,[n_1:n]))));
    Secret_Key = keys(indice);
    Key_values(1,i)=Secret_Key;
end
clc;
pourcent = 50;
fprintf(formatSpec,pourcent)
%toc;

%% Second round
if mode == 32
    % True Result
    XOR = bitxor(plaintexts,Key_values(1:16));
    SboxResult = sbox(XOR);
    ShiftRowsResult = shiftrows(SboxResult);
    MixColumnResult = mixcolumns(ShiftRowsResult);
    hypothetic_inter1_secondRd = zeros(nb_traces,nb_columns);
    for byte=1:16
        hypothetic_inter1_secondRd(:,[((256*mode)-(256*((mode+1)-byte)))+1:((256*mode)-(256*(mode-byte)))]) = bitxor(MixColumnResult(:,byte),keys);
    end
    hypothetic_inter2_secondRd = sbox(hypothetic_inter1_secondRd);
    clc;
    pourcent = 75;
    fprintf(formatSpec,pourcent)

    for byte=1:16
        Col = [((256*mode)-(256*((mode+1)-byte)))+1:((256*mode)-(256*(mode-byte)))];
        dist(:,Col) =  bitxor(SboxResult(:,byte),hypothetic_inter2_secondRd(:,Col));
    end
    bin_dist = dec2bin(dist); % Jusqu'ici Okay !
    Dist_hamming_desordonnate = sum(bin_dist.' == '1');
    Dist_hamming = reshape(Dist_hamming_desordonnate,[nb_traces,nb_columns]);
    clc;
    pourcent = 95;
    fprintf(formatSpec,pourcent)
    %toc;

    %% Calculate the correlation coefficient
    coefficient_secondRd = corr(Dist_hamming, raw_traces_filt);
    Result_secondRd = coefficient_secondRd';
    %toc;

    %% Key values
    Key_values_secondRd = zeros(1,16);
    for i=1:16
        n_1 = (i*256)-255;
        n=i*256;
        [maximum,indice] = max(max(abs(Result_secondRd(:,[n_1:n]))));
        Secret_Key = keys(indice);
        Key_values_secondRd(1,i)=Secret_Key;
    end
end
clc;
pourcent = 100;
fprintf(formatSpec,pourcent)
%toc;

%% To compare the keys
if mode == 16
    ResponseObtained = Key_values';
    ResponseExcepted = keys_unique';
elseif mode == 32
    Key_values = [Key_values, Key_values_secondRd];
    ResponseObtained = Key_values';
    if Enable_Faking_Implementation == 0
        ResponseExcepted = keys_unique';    
    elseif Enable_Faking_Implementation == 1
        ResponseExcepted = Fake_Key';   
    end
end


%% Print result
clc;
error = 0;
fprintf('=> Results : \n');  

% Round 1
fprintf('-> Attack Round 1 \n');  
fprintf('N°  Byte  : ');
for i=1:mode/2
    if i<10
        fprintf('%i     ', i);  
    else    
        fprintf('%i    ', i); 
    end
end
fprintf('\n');
fprintf('Real  Key : ');
for i=1:mode/2
    fprintf('0x%s  ', dec2hex(ResponseExcepted(i), 2));  
end
fprintf('\n');
fprintf('Found Key : ');
for i=1:mode/2
    if ResponseExcepted(i) ==  ResponseObtained(i)
        fprintf('0x%s  ', dec2hex(ResponseObtained(i), 2));
    else
        fprintf(2,'0x%s  ', dec2hex(ResponseObtained(i), 2));
        error = error + 1;
    end
end

% Round 2
fprintf('\n \n');
fprintf('-> Attack Round 2 \n'); 
fprintf('N°  Byte  : ');
for i=1:mode/2
    if i<10
        fprintf('%i     ', i);  
    else    
        fprintf('%i    ', i); 
    end
end
fprintf('\n');
fprintf('Real  Key : ');  
for i=17:32
    fprintf('0x%s  ', dec2hex(ResponseExcepted(i), 2));  
end
fprintf('\n');
fprintf('Found Key : ');
for i=17:32
    if ResponseExcepted(i) ==  ResponseObtained(i)
        fprintf('0x%s  ', dec2hex(ResponseObtained(i), 2));
    else
        fprintf(2,'0x%s  ', dec2hex(ResponseObtained(i), 2));
        error = error + 1;
    end
end

fprintf('\n \n');
formatSpec = '==> There are %i byte(s) of error (on %i bytes) for %i traces. \n';
fprintf(formatSpec,error, mode, nb_traces)


