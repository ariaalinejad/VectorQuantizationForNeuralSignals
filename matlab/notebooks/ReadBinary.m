% This file reads the data from the HC-1 dataset and saves it as .bin. The
% hc-1 dataset can be downloaded at: https://crcns.org/data-sets/hc/hc-1. 
% The LoadBinary is neccesarry for runnin this file. Go to:
% https://portal.nersc.gov/project/crcns/download/hc-2
% make a profile, and download: crcns-hc2-scripts.zip


%clc, close all, clear all;

[Data OrigIndex] = LoadBinary('C:\Users\ariaa\OneDrive - NTNU\Masteroppgave\Code & Data\matlab\src\read_hc1\d533101.dat', [2,3,4,5,6,7,8]);

offset = 2048;
fs = 10000;
Data_extracellular = Data(1:4, :) - offset;
Data_intracellular = Data(5:7, :);

%sprintf('Size extracellular data: %d, %d \n', size(Data_extracellular)')
%sprintf('Size intracellular data: %d, %d \n', size(Data_intracellular)')

figure(1);
for i = 1:4
    subplot(2,2,i);
    plot(Data_extracellular(i,1:480000));
end

t_max = size(Data_intracellular, 2)/fs;
stepsize = t_max/( size(Data_intracellular, 2) - 1 );
t = 1:stepsize:t_max;

M = 4800;
figure(2);
subplot(3,1,1);
plot(t(1:M), Data_intracellular(1,1:M));
title('Intracellular signal')
subplot(3,1,2);
plot(t(1:M), Data_intracellular(2,1:M));
title('Membrane potential')
subplot(3,1,3);
plot(t(1:M), Data_intracellular(3,1:M));
title('Injected current')
xlabel('seconds [s]');

%% Plot intra against reconstructed:
n_1 = 30690;
n_2 = 30710;
n_3 = 600001;


% f1 = 300;  % Cutoff frequency in Hz
% f2 = 4000;
% fs = 10000; % Sampling rate in Hz
% [b, a] = butter(2, [f1,f2]/(fs/2), "bandpass"); % Design the filter
% %freqz(b, a, [], fs); % Plot magnitude and phase responses
% intra_filt = filter(b, a, Data_intracellular(1,:));
% plot(intra_filt);
% intra_filt = intra_filt(n_3+n_1:n_3+n_2);


intra_normal = -Data_intracellular(1,n_3+n_1:n_3+n_2) - max(-Data_intracellular(1,n_3+n_1:n_3+n_2));
extra_normal = temp(n_1:n_2,1) - max(temp(n_1:n_2,1));
re_normal = temp_hat(n_1:n_2,1) - max(temp_hat(n_1:n_2,1));


% intra_normal = -(Data_intracellular(1,n_3+n_1:n_3+n_2) - Data_intracellular(1,n_3+n_1));%mean(Data_intracellular(1,n_3+n_1:n_3+n_2+1000)));
[V_intra, I_intra] = max(abs(intra_normal));
[V_extra, I_extra] = max(abs(extra_normal));
[V_re, I_re] =       max(abs(re_normal));


delta = I_intra - I_extra;
extra_scaled = extra_normal./V_extra .* V_intra;

delta = I_intra - I_re;
re_scaled = re_normal./V_re .* V_intra;


figure(3);
% plot(Data_extracellular(1,n_1:n_2));
plot(intra_normal);
hold on
plot(extra_scaled)
plot(re_scaled)
hold off
legend('intra', 'extra', 're')


n_1 = 1;
n_2 = 50000;
intra_normal = -Data_intracellular(1,n_3+n_1:n_3+n_2) - mean(-Data_intracellular(1,n_3+n_1:n_3+n_2));
extra_normal = temp(n_1:n_2,1) - mean(temp(n_1:n_2,1));
re_normal = temp_hat(n_1:n_2,1) - mean(temp_hat(n_1:n_2,1));


[V_intra, I_intra] = max(abs(intra_normal));
[V_extra, I_extra] = max(abs(extra_normal));
[V_re, I_re] =       max(abs(re_normal));

delta = I_intra - I_extra;
extra_scaled = extra_normal./V_extra .* V_intra;

delta = I_intra - I_re;
re_scaled = re_normal./V_re .* V_intra;

figure(4);
% plot(Data_extracellular(1,n_1:n_2));
plot(intra_normal);
hold on
plot(extra_scaled)
plot(re_scaled)
hold off
legend('intra', 'extra', 're')

P_d = norm(intra_normal);
P_error = norm(intra_normal - extra_scaled);
SNDR_intra = 20*log10(P_d/P_error);
disp(['SNDR (extra): ', num2str(SNDR_intra, '%.2f'), ' dB']);

P_d = norm(intra_normal);
P_error = norm(intra_normal - re_scaled);
SNDR_intra = 20*log10(P_d/P_error);
disp(['SNDR (reconstruct): ', num2str(SNDR_intra, '%.2f'), ' dB']);


%% Save to bin

% Write the data to a binary file
fileID = fopen('hc1_d533101_extracellular.bin', 'wb');
fwrite(fileID, Data_extracellular, 'double');
fclose(fileID);

%% Save to bin 240000 samples
fileID = fopen('hc1_d533101_intracellular_240000.bin', 'wb');
fwrite(fileID, Data_intracellular(1,240000:480000), 'double');
fclose(fileID);

%% Save to bin 600000 samples
fileID = fopen('hc1_d533101_intracellular_600000.bin', 'wb');
fwrite(fileID, Data_intracellular(1,600000+1:1200000), 'double');
fclose(fileID);

%% Save to bin 240000 samples for sq-dict
fileID = fopen('hc1_d533101_intracellular_sq-dict_240000.bin', 'wb');
fwrite(fileID, Data_intracellular(1,240000*3+1:240000*4), 'double');
fclose(fileID);

%% Save extracellular as .mat
x = Data_extracellular;
save('C:\Users\ariaa\OneDrive - NTNU\Masteroppgave\Code & Data\data\mat\hc1_d533101_extracellular.mat', 'x');

%% Save intracellular as .mat
x = Data_intracellular;
save('C:\Users\ariaa\OneDrive - NTNU\Masteroppgave\Code & Data\data\mat\hc1_d533101_intracellular.mat', 'x');