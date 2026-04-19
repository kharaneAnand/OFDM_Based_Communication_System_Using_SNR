%% ================================================
%% SETUP.m — Run ONCE on TX PC, copy files to RX
%% ================================================
clear; clc; j = 1i;

fs  = 1e6;
fc  = 866e6;

%% Long preamble
N_FFT = 64;
L_k   = [1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,...
         0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1];
vsc   = zeros(1, N_FFT-length(L_k));
LP_freq = [vsc(1:6), L_k, vsc(7:11)];

%% Payload
rng(42); data_Payload_1 = randi([0 3], 1, 48);
rng(43); data_Payload_2 = randi([0 3], 1, 48);

%% PN sequence for sync — Gold code style, 127 chips
rng(99);
pn_seq = 2*randi([0 1], 1, 127) - 1;   % ±1 BPSK

save('OFDM_config.mat', 'fs', 'fc', 'LP_freq', ...
     'data_Payload_1', 'data_Payload_2', 'pn_seq');
fprintf('OFDM_config.mat saved — copy to RX PC\n');