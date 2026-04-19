%% ================================================
%% OFDM_Setup.m — Run ONCE on TX PC
%% ================================================
clear; clc; j = 1i;

fs  = 1e6;
fc  = 866e6;

N_FFT = 64;
L_k   = [1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,...
         0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1];
vsc     = zeros(1, N_FFT-length(L_k));
LP_freq = [vsc(1:6), L_k, vsc(7:11)];

%% ── TEXT MESSAGE — 32 characters ────────────────────────────────────
msg = 'HELLO FROM IITJ GROUP10 SDR OTA';

msg_padded = char(32 * ones(1, 32));
msg_padded(1:min(length(msg),32)) = msg(1:min(length(msg),32));
fprintf('Message : "%s"\n', msg_padded);

%% ASCII → bits: 32 × 8 = 256 bits
bits_orig = reshape(dec2bin(msg_padded, 8)' - '0', 1, []);   % [1x256]

%% ── INTERLEAVED rep-3 encoding ───────────────────────────────────────
%% Instead of: [b1 b1 b1 b2 b2 b2 ... b256 b256 b256]  (copies adjacent)
%% We use:     [b1 b2 ... b256 | b1 b2 ... b256 | b1 b2 ... b256]
%% This spreads the 3 copies of each bit 256 positions apart
%% Bad subcarrier at position X only corrupts copy-1 of bit X
%% Copies 2 and 3 of bit X are at X+256 and X+512 — safe
bits = [bits_orig, bits_orig, bits_orig];   % [1x768] interleaved rep-3
fprintf('bits_orig=%d | after interleaved rep-3=%d\n', ...
        length(bits_orig), length(bits));

%% Bits → QPSK symbols
symbols_all = bi2de(reshape(bits, 2, [])', 'left-msb').';   % [1x384]

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
fprintf('OFDM_config.mat saved — copy to RX PC\n');
