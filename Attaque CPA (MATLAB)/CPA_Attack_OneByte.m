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
NUMBER_TRACES = 51000;
sampling_interval_ns = 2;
Enable_Attack = 1;
Enable_Plot_traces = 0;
Enable_Faking_Implementation = 1;
Enable_plot_correlation = 0;
Enable_plot_corr_vs_keys = 1;
Enable_plot_One_trace = 1;
byte = 5 ;

%% Load the data
Nth_measurement = strcat(int2str(NUMBER_TRACES),'_traces');
formatSpec = 'In progress %i%% ... \n';
pourcent = 2;
fprintf(formatSpec,pourcent)
path = '..\Data\AES_256\Measurement_';
totalPath = strcat(path, Nth_measurement,'\');
traceFile = strcat(totalPath,'traces.txt');
load(traceFile)
%load('traces.txt')
%traces_brute = csvread(traceFile);
%traces=31.2500E-6*traces_brute;
plaintextFile = strcat(totalPath,'plaintexts.mat');
load(plaintextFile);
keyFile = strcat(totalPath,'keys_unique.mat');
load(keyFile);

plaintexts = plaintexts([5000:40000],:);
traces = traces([5000:40000],:);

if Enable_Faking_Implementation == 1
    maskFile = strcat(totalPath,'masks_unique.mat');
    load(maskFile);
    Fake_Key = bitxor(keys_unique',mask_unique_dec);
    Fake_Key_hex = conversion(Fake_Key);
end

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
coefficient = corr(Weight_Hamm, traces);
Result = coefficient';
clc;
pourcent = 30;
fprintf(formatSpec,pourcent)
%toc;

%% Plot data's
if Enable_plot_One_trace == 1
    figure
    nb_samples = size(traces,2);
    times = (0:1:nb_samples-1)*(1/(1e9/2)); % en seconde
    times = 1e9*times;                      % en ns
    plot(times, traces(1000,:));
    grid on;
    title('Trace pour implémentation séquentielle')
    xlabel('Temps [ns]')
    ylabel('Tension [mV]')
end

figure
subplot(2,1,1)
for i=1:nb_traces
    plot(traces(i,:));
    hold on;
end
grid on;
title('Trace pour implémentation séquentielle')
xlabel('Temps [ns]')
ylabel('Tension [mV]')

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
    Fake_Key_Attacked = Fake_Key(byte);
    True_Key_Attacked = keys_unique(byte);
    [ValuesMaxCorr,indiceValuesMaxCorr] = max(abs(Result));
    X = linspace(0,255,256);
    plot(X, Result(indiceValuesMaxCorr(indice),:))
    hold on;
    CorrMax = Result(indiceValuesMaxCorr(indice),Secret_Key+1);
    p1 = plot(Secret_Key,CorrMax, "m*", 'DisplayName','Clé révélée');
    hold on
    text(Secret_Key+4,CorrMax,'Corrélation maximum');
    title('Coefficient de corrélation vs Clés possibles')
    xlabel('Valeurs de clés possibles')
    ylabel('Coefficient de corrélation')
    minus = min(min(Result));
    maxus = max(max(Result));
    vertical = [minus, maxus + maxus*0.2];
    p2 = plot(Fake_Key_Attacked*ones(1,2), vertical ,'r--', 'DisplayName','Fausse clé');
    hold on
    p3 = plot(True_Key_Attacked*ones(1,2), vertical ,'g--', 'DisplayName','Vraie clé');
    hold on;
    grid on
    grid minor
    ylim([minus (maxus + maxus*0.2)])
    legend([p1 p2 p3])

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
R_part = zeros(nb_traces-3, len_keys);
% pour trouver l'échantillon à annalyser : 
All_SAMPLING = max(abs(coefficient));
[maximum_sampling, ind_sampling_max] = max(All_SAMPLING);
for i=4:nb_traces
    T_part = traces(1:i,:);
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
    R_part(i-3,:) = B;
end
figure()
plot(1:(nb_traces-3), R_part(:,1:Key_values), 'k', 1:(nb_traces-3), R_part(:,Key_values+2:end), 'k')
hold on
p1 = plot(1:(nb_traces-3), R_part(:,Key_values-10),'k', 'DisplayName','Autres clés');
hold on
p2 = plot(1:(nb_traces-3), R_part(:,Key_values+1),'m', 'DisplayName','Clé révélée');
hold on
p3 = plot(1:(nb_traces-3), R_part(:,Fake_Key_Attacked),'r', 'DisplayName','Fausse clé');
hold on
p4 = plot(1:(nb_traces-3), R_part(:,True_Key_Attacked),'g', 'DisplayName','Vraie clé');
hold on
title('Coefficient de corrélation vs nombre de traces')
xlabel('Nombre de traces')
ylabel('Coefficient de corrélation')
legend([p1 p2 p3 p4])
ylim([-0.4 0.4])
grid on
grid minor
