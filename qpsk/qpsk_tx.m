%% =========================================
%% FINAL QPSK TX (BURST MODE + PREAMBLE)
%% =========================================

clear; clc;

%% PARAMETERS
fs = 1e6;
fc = 866e6;
M = 4;
sps = 4;
rolloff = 0.35;
span = 10;

fprintf("TX running @ %.0f MHz\n", fc/1e6);

%% DATA
numSymbols = 2000;

rng(1);
preambleBits = randi([0 1], 200, 1);   % KNOWN PREAMBLE

rng(0);
payloadBits = randi([0 1], numSymbols*2, 1);

txBits = [preambleBits; payloadBits];

symbols = qammod(txBits, M, ...
    'InputType','bit', ...
    'UnitAveragePower', true);

%% RRC FILTER
rrc = rcosdesign(rolloff, span, sps, 'sqrt');
txSig = upfirdn(symbols, rrc, sps);

%% NORMALIZE
txSig = txSig / max(abs(txSig));

%% BURST STRUCTURE (IMPORTANT)
guard = zeros(2000,1);
txFrame = [guard; txSig; guard];

%% USRP TX
tx = comm.SDRuTransmitter( ...
    'Platform','B200', ...
    'SerialNum','34D8D87', ...
    'CenterFrequency', fc, ...
    'MasterClockRate', 20e6, ...
    'InterpolationFactor', 20, ...
    'Gain', 45);

fprintf("Transmitting bursts...\n");

while true
    tx(txFrame);
end

release(tx);