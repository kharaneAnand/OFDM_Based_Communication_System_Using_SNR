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

%% ── TEXT MESSAGE ─────────────────────────────────────────────────────
%% Repetition coding: each bit sent 3 times
%% 192 bits available → 64 original bits → 8 characters max
msg = 'ABDe ugh';   % max 8 characters — edit freely

%% Pad to exactly 8 characters
msg_padded = char(32 * ones(1, 8));
msg_padded(1:min(length(msg),8)) = msg(1:min(length(msg),8));
fprintf('Message : "%s"\n', msg_padded);

%% ASCII → bits (8 bits x 8 chars = 64 bits)
bits_orig = reshape(dec2bin(msg_padded, 8)' - '0', 1, []);   % [1x64]

%% Repetition encode — each bit repeated 3 times
%% 64 bits → 192 bits
bits = repelem(bits_orig, 3);   % [1x192]
fprintf('Bits original : %d | After rep-3 coding : %d\n', ...
        length(bits_orig), length(bits));

%% Bits → QPSK symbols (bit pairs → 0,1,2,3)
symbols_all    = bi2de(reshape(bits, 2, [])', 'left-msb').';   % [1x96]
data_Payload_1 = symbols_all(1:48);
data_Payload_2 = symbols_all(49:96);

%% PN sequence
rng(99);
pn_seq = 2*randi([0 1], 1, 127) - 1;

save('OFDM_config.mat', 'fs', 'fc', 'LP_freq', ...
     'data_Payload_1', 'data_Payload_2', 'pn_seq', ...
     'msg_padded', 'bits', 'bits_orig');
fprintf('OFDM_config.mat saved — copy to RX PC\n');