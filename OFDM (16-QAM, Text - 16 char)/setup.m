%% ================================================
%% OFDM_Setup.m — 16-QAM 16 CHARS SOFT-COMBINING
%% ================================================
clear; clc; j = 1i;

fs  = 1e6;
fc  = 866e6;

N_FFT = 64;
L_k   = [1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,...
         0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1];
vsc     = zeros(1, N_FFT-length(L_k));
LP_freq = [vsc(1:6), L_k, vsc(7:11)];

%% ── TEXT MESSAGE — 16 characters ────────────────────────────────────
%% 16 chars × 8 bits = 128 bits / 4 bits per sym = 32 symbols per copy
%% 384 total symbols / 32 = 12 copies for soft combining
%% 12x averaging gives ~10.8 dB SNR improvement
msg = 'IITJ SDR 16QAM!!';   % exactly 16 chars

msg_padded = char(32 * ones(1, 16));
msg_padded(1:min(length(msg),16)) = msg(1:min(length(msg),16));
fprintf('Message : "%s"\n', msg_padded);

%% ASCII → bits: 16 × 8 = 128 bits
bits_orig = reshape(dec2bin(msg_padded, 8)' - '0', 1, []);   % [1x128]

%% 12 copies interleaved: 128 × 12 = 1536 bits / 4 = 384 symbols — exact fit
bits = repmat(bits_orig, 1, 12);   % [1x1536]
fprintf('bits_orig=%d | 12 copies=%d | symbols=%d\n', ...
        length(bits_orig), length(bits), length(bits)/4);

symbols_all = bi2de(reshape(bits, 4, [])', 'left-msb').';   % [1x384]

data_Payload_1 = symbols_all(1:48);
data_Payload_2 = symbols_all(49:96);
data_Payload_3 = symbols_all(97:144);
data_Payload_4 = symbols_all(145:192);
data_Payload_5 = symbols_all(193:240);
data_Payload_6 = symbols_all(241:288);
data_Payload_7 = symbols_all(289:336);
data_Payload_8 = symbols_all(337:384);

rng(99);
pn_seq = 2*randi([0 1], 1, 127) - 1;

save('OFDM_config.mat', 'fs', 'fc', 'LP_freq', ...
     'data_Payload_1', 'data_Payload_2', 'data_Payload_3', 'data_Payload_4', ...
     'data_Payload_5', 'data_Payload_6', 'data_Payload_7', 'data_Payload_8', ...
     'pn_seq', 'msg_padded', 'bits_orig');
fprintf('OFDM_config.mat saved\n');