clear;
clc ;
tic;

%% Configuration :
% - The user defines the AES's mode (128=16 ; 256=32)
% - The user defines the number of traces that he would like test
mode = 32;
NUMBER_TRACES = 52000;
plot_power = 0;

% %% Plot data from oscilloscope (normal)
% traces_brute = load('testTraces2.csv');
% traces=31.2500E-6*traces_brute;
% plaintextFile = strcat(totalPath,'plaintexts.mat');
% load(plaintextFile);
% keyFile = strcat(totalPath,'keys_unique.mat');
% load(keyFile);
% 
% [y,x]=size(traces);
% OneTrace = traces(129,:);
% X = x;
% Y = OneTrace;
% figure()
% plot(A)

%% Load the data
Nth_measurement = strcat(int2str(NUMBER_TRACES),'_traces');
path = '..\Data\AES_256\Measurement_';
totalPath = strcat(path, Nth_measurement,'\');
traceFile = strcat(totalPath,'traces.txt');
load(traceFile)
timeFile = strcat(totalPath,'time.txt');
load(timeFile)

tracesX = time(1,:);
%% Plot 2D 4 traces
figure
subplot(2,2,1)
plot(tracesX, traces(1,:))
title('Première trace')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on
subplot(2,2,2)
plot(tracesX, traces(2,:))
title('Deuxième trace')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on
subplot(2,2,3)
plot(tracesX, traces(3,:))
title('Troisième trace')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on
subplot(2,2,4)
plot(tracesX, traces(4,:))
title('Quatrième trace')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on


%% Plot 2D 4 traces
figure
subplot(2,2,1)
plot(tracesX, traces(20,:))
title('Trace 20')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on
subplot(2,2,2)
plot(tracesX, traces(21,:))
title('Trace 21')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on
subplot(2,2,3)
plot(tracesX, traces(22,:))
title('Trace 22')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on
subplot(2,2,4)
plot(tracesX, traces(23,:))
title('Trace 23')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on


%% Plot 2D 4 traces
figure
subplot(2,2,1)
plot(tracesX, traces(300,:))
title('Trace 300')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on
subplot(2,2,2)
plot(tracesX, traces(301,:))
title('Trace 301')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on
subplot(2,2,3)
plot(tracesX, traces(302,:))
title('Trace 302')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on
subplot(2,2,4)
plot(tracesX, traces(303,:))
title('Trace 303')
xlabel("Nombre d'échantillons")
ylabel('Tension')
grid on

if plot_power == 1
    figure()
    plot(traces(1,[250:1250]))
    title('Puissance statique')
    xlabel('Temps [ns]')
    ylabel('Tension [mv]')
    ylim([-20 20]);
    grid on;
    hold on ;
    plot(traces(4000,[250:1250]))
    title('Puissance statique & dynamique')
    xlabel('Temps [ns]')
    ylabel('Tension [mv]')
    ylim([-20 20]);
    grid on;
    hold on ;
    legend("Puissance statique", "Puissance statique & dynamique")
end
