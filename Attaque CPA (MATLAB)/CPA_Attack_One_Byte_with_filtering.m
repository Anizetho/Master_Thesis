%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         CPA attack (AES-256)   HW                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc ;
tic;

%% Configuration :
% - The user defines the AES's mode (128=16 ; 256=32)
% - The user defines the number of traces that he would like test
mode = 32;
NUMBER_TRACES = 40000;
sampling_interval_ns = 2;
Enable_Attack = 1;
Enable_Plot_traces = 0;
Enable_Faking_Implementation = 1;
Enable_plot_correlation = 0;
Enable_plot_corr_vs_keys = 1;
byte = 14 ;

%% Load the data
Nth_measurement = strcat(int2str(NUMBER_TRACES),'_traces');
formatSpec = 'In progress %i%% ... \n';
pourcent = 2;
fprintf(formatSpec,pourcent)
path = '..\Data\AES_256\Measurement_';
totalPath = strcat(path, Nth_measurement,'\');
traceFile = strcat(totalPath,'traces.txt');
load(traceFile)
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
%freqz(B, 1, nb_samples);
%legend(result_hamming,'Hamming')
%hfvt = fvtool(B,1);
%legend(hfvt,'Hamming')

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

%% Plot time domain and frequency domain
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
hypothetic_inter1 = bitxor(plaintexts(:,byte),keys);
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

figure
subplot(2,1,1)
for i=1:nb_traces
    plot(raw_traces_filt(i,:));
    hold on;
end
grid on;
title('Traces')
xlabel('Time sample')
ylabel('Voltage [V]')

subplot(2,1,2)
plot(coefficient);
hold on;
grid on;
title('Correlation')
xlabel('Time sample')
ylabel('SNR')


        
%% Key values
Key_values = zeros(1,mode/2);
[maximum,indice] = max(max(abs(Result)));
Secret_Key = keys(indice);
Key_values =Secret_Key;
clc;
pourcent = 50;
fprintf(formatSpec,pourcent)
%toc;

%% Plot correlation vs keys
if Enable_plot_corr_vs_keys == 1
    figure()
    [ValuesMaxCorr,indiceValuesMaxCorr] = max(abs(Result));
    X = linspace(0,255,256);
    plot(X, Result(indiceValuesMaxCorr(indice),:))
    hold on;
    CorrMax = Result(indiceValuesMaxCorr(indice),Secret_Key+1);
    plot(Secret_Key,CorrMax, "r*");
    hold on
    text(Secret_Key+4,CorrMax,'Corrélation maximum');
    title('Coefficient de corrélation vs Clés possibles')
    xlabel('Valeurs de clés possibles')
    ylabel('Coefficient de corrélation')
    %ylim([-0.5 0.5])
    grid on
    grid minor

    %% Plot correlation vs time
    figure()
    subplot(2,2,1)
    X = linspace(0,nb_sampling*sampling_interval_ns-2,(nb_sampling*sampling_interval_ns)/2);
    plot(X,coefficient(38,:))
    title('Clé 38')
    xlabel('Temps [ns]')
    ylabel('Coefficient de corrélation')
    ylim([-0.2 0.2])
    grid on;
    hold on;
    subplot(2,2,2)
    plot(X, coefficient(38+1,:))
    ylim([-0.2 0.2])
    title('Clé 39')
    xlabel('Temps [ns]')
    ylabel('Coefficient de corrélation')
    grid on;
    hold on;
    subplot(2,2,3)
    plot(X, coefficient(237+1,:))
    ylim([-0.2 0.2])
    title('Clé 237')
    xlabel('Temps [ns]')
    ylabel('Coefficient de corrélation')
    grid on;
    hold on;
    subplot(2,2,4)
    plot(X, coefficient(237+2,:))
    ylim([-0.2 0.2])
    title('Clé 238')
    xlabel('Temps [ns]')
    ylabel('Coefficient de corrélation')
    grid on;
    hold on;
end

%% To compare the keys
if mode == 16
    ResponseObtained = Key_values';
    ResponseExcepted = keys_unique';
elseif mode == 32
    ResponseObtained = Key_values';
    ResponseExcepted = keys_unique';    
end


%% Print result
clc;
error = 0;
fprintf('=> Results : \n');  

% Round 1
fprintf('-> Attack Round 1 \n');  
formatSpec = 'N°  Byte  : %i';
fprintf(formatSpec,byte)
fprintf('\n');
fprintf('Real  Key : ');
fprintf('0x%s  ', dec2hex(ResponseExcepted(byte), 2));  
fprintf('\n');
fprintf('Found Key : ');
if ResponseExcepted(byte) ==  ResponseObtained
    fprintf('0x%s  ', dec2hex(ResponseObtained, 2));
else
    fprintf(2,'0x%s  ', dec2hex(ResponseObtained, 2));
    error = error + 1;
end
fprintf('\n');


%% Plot Correlation vs Number of Traces
R_part = zeros(nb_traces-5, len_keys);
D = zeros(1,256);
% pour trouver l'échantillon à annalyser : 
All_SAMPLING = max(abs(coefficient));
[maximum_sampling, ind_sampling_max] = max(All_SAMPLING);
for i=6:nb_traces
    T_part = raw_traces_filt(1:i,:);
    H_part = Weight_Hamm(1:i,:);
    A = corr(H_part, T_part);
    B = A(:,ind_sampling_max);
    
%     B = max(abs(A(:,150)))
%     B = max(A');
%     C = min(A');
%     for j = 1:256
%        if B(1,j) > abs(C(1,:))
%            D(1,j) = B(1,j);
%        else
%            D(1,j) = C(1,j);
%        end
%     end
    R_part(i-5,:) = B;
end

nb_traces = size(R_part, 1);
Key_values = 13;
p1 = plot(1:(nb_traces), R_part(:,Key_values+1),'g', 'DisplayName','Clé secrète');
hold on
plot(1:(nb_traces), R_part(:,1:Key_values), 'k', 1:(nb_traces), R_part(:,Key_values+2:end), 'k')
p2 = plot(1:(nb_traces), R_part(:,Key_values-10),'k', 'DisplayName','Autres clés');
title('Coefficient de corrélation vs nombre de traces (pour des traces filtrées)')
xlabel('Nombre de traces')
ylabel('Coefficient de corrélation')
legend([p1 p2])
grid on
grid minor
