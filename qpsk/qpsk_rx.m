%% =========================================
%% FINAL QPSK RX (BURST + SYNC + CORRECT BER)
%% =========================================
clear; clc;
%% PARAMETERS
fs = 1e6;
fc = 866e6;
M = 4;
sps = 4;
rolloff = 0.35;
span = 10;
numSymbols = 2000;
%% REFERENCE DATA (MUST MATCH TX)
rng(1);
preambleBits = randi([0 1], 200, 1);
rng(0);
refBits = randi([0 1], numSymbols*2, 1);
preambleSym = qammod(preambleBits, M, ...
    'InputType','bit', ...
    'UnitAveragePower', true);
%% RTL RECEIVER
rx = comm.SDRRTLReceiver( ...
    'CenterFrequency', fc, ...
    'SampleRate', fs, ...
    'SamplesPerFrame', 16384, ...
    'OutputDataType', 'double', ...
    'EnableTunerAGC', false, ...
    'TunerGain', 35);
%% FILTER + SYNCHRONIZATION
rrc = rcosdesign(rolloff, span, sps, 'sqrt');
carrierSync = comm.CarrierSynchronizer( ...
    'Modulation','QPSK', ...
    'SamplesPerSymbol', sps, ...
    'NormalizedLoopBandwidth', 0.02);
timingSync = comm.SymbolSynchronizer( ...
    'Modulation','PAM/PSK/QAM', ...
    'SamplesPerSymbol', sps, ...
    'NormalizedLoopBandwidth', 0.02);
fprintf("RX running...\n");
while true

    rxRaw = rx();

    %% SIGNAL DETECTION
    power = mean(abs(rxRaw).^2);
    if power < 0.005
        fprintf("No signal\n");
        continue;
    end

    %% SYNC CHAIN
    rxC = carrierSync(rxRaw);
    rxFilt = filter(rrc, 1, rxC);
    rxSync = timingSync(rxFilt);

    if length(rxSync) < 500
        continue;
    end

    %% REMOVE FILTER TRANSIENT
    rxSync = rxSync(span:end-span);

    %% COARSE PHASE CORRECTION
    phaseEst = angle(mean(rxSync.^4)) / 4;
    rxSync = rxSync * exp(-1j*phaseEst);

    %% FRAME SYNC USING PREAMBLE
    corr = abs(conv(rxSync, flip(conj(preambleSym))));
    [~, idx] = max(corr);

    startIdx = idx - length(preambleSym) + 1;

    if startIdx < 1 || startIdx + length(preambleSym) + numSymbols - 1 > length(rxSync)
        continue;
    end

    %% EXTRACT FRAME
    rxFrame = rxSync(startIdx : startIdx + length(preambleSym) + numSymbols - 1);

    %% REMOVE PREAMBLE
    rxPayload = rxFrame(length(preambleSym)+1:end);

    %% NORMALIZE
    rxPayload = rxPayload / rms(rxPayload);

    %% TRY ALL 4 QPSK ROTATIONS (FIX MAPPING ISSUE)
    bestBER = inf;

    for k = 0:3

        rotated = rxPayload * exp(1j * k * pi/2);

        rxBits = qamdemod(rotated, M, ...
            'OutputType','bit', ...
            'UnitAveragePower', true);

        L = min(length(rxBits), length(refBits));
        errors = sum(rxBits(1:L) ~= refBits(1:L));
        ber = errors / L;

        if ber < bestBER
            bestBER = ber;
            bestSymbols = rotated;
        end
    end

    %% DISPLAY RESULTS
    fprintf("BER = %.6f\n", bestBER);

    figure(1); clf;
    scatter(real(bestSymbols), imag(bestSymbols), 10, 'filled');
    grid on;
    axis([-2 2 -2 2]); axis square;
    title(sprintf('QPSK FINAL | BER = %.6f', bestBER));
    drawnow;

end
release(rx);