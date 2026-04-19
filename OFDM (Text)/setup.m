%% 
%% ================================================
%% OFDM_Setup.m — Run ONCE on TX PC, copy files to RX
%% ================================================
clear; clc; j = 1i;

fs  = 1e6;
fc  = 866e6;

%% Long preamble
N_FFT = 64;
L_k   = [1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,...
         0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1];
vsc     = zeros(1, N_FFT-length(L_k));
LP_freq = [vsc(1:6), L_k, vsc(7:11)];

%% ── TEXT MESSAGE (max 24 characters) ────────────────────────────────
msg = 'HELLO OFDM IITJ GRP10';   % edit this — max 24 chars

%% Pad to exactly 24 characters with spaces
msg_padded = char(32 * ones(1, 24));
msg_padded(1:min(length(msg),24)) = msg(1:min(length(msg),24));
fprintf('Message: "%s"\n', msg_padded);

%% ASCII → bits (8 bits per char = 192 bits total)
bits = reshape(dec2bin(msg_padded, 8)' - '0', 1, []);   % [1x192]

%% Bit pairs → QPSK symbols (0–3)
%% bi2de returns a column vector — transpose to row for TX compatibility
symbols_all    = bi2de(reshape(bits, 2, [])', 'left-msb').';  % [1x96] row
data_Payload_1 = symbols_all(1:48);    % [1x48] row
data_Payload_2 = symbols_all(49:96);   % [1x48] row
%% PN sequence
rng(99);
pn_seq = 2*randi([0 1], 1, 127) - 1;

save('OFDM_config.mat', 'fs', 'fc', 'LP_freq', ...
     'data_Payload_1', 'data_Payload_2', 'pn_seq', ...
     'msg_padded', 'bits');
fprintf('OFDM_config.mat saved — copy to RX PC\n');