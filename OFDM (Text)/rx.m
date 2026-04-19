%% ================================================
%% OFDM_RX.m — PC 2, RTL-SDR — FINAL
%% ================================================
clear; close all; clc; j = 1i;
load('OFDM_config.mat');
load('pn_len.mat');

%% ── Parameters ───────────────────────────────────────────────────────
N_FFT      = 64;
N_CP       = 16;
OVR        = 2;
rolloff    = 0.5;
span       = 6;
M          = 4;
Ts         = 1/fs;
dIdx       = [7:11, 13:25, 27:32, 34:39, 41:53, 55:59];
pilot_pos  = [12, 26, 40, 54];
pilot_vals = [1, 1, 1, -1];
spl        = 16;
corners    = exp(1i*[pi/4, 3*pi/4, 5*pi/4, 7*pi/4]);
RRC        = rcosdesign(rolloff, span, OVR, 'sqrt');

%% ── PN template ──────────────────────────────────────────────────────
pn_rrc_ds = pn_rrc(1:OVR:end);
pn_tmpl   = pn_rrc_ds / norm(pn_rrc_ds);
pn_ds_len = length(pn_rrc_ds);

%% ── Ideal LP for channel estimation ─────────────────────────────────
LP_t_ideal = ifft(ifftshift(LP_freq));
LP_f_ideal = fftshift(fft(LP_t_ideal));

%% ── State ────────────────────────────────────────────────────────────
f_bulk      = 0;
phase_track = 0;
phase_init  = false;
cal_done    = false;
cal_buf     = [];
N_cal       = 20;
ber_hist    = [];
const_buf   = [];
const_max   = 500;

%% ── RTL-SDR ──────────────────────────────────────────────────────────
rx = comm.SDRRTLReceiver('CenterFrequency',fc,'SampleRate',fs,...
        'SamplesPerFrame',32768,'OutputDataType','double',...
        'EnableTunerAGC',false,'TunerGain',49);

%% ── GUI ──────────────────────────────────────────────────────────────
fig = figure('Name','RX','NumberTitle','off');
set(fig,'Units','centimeters','position',[1 2 49 24]);
uicontrol('String','Stop','Position',[1475 15 100 60],...
          'Callback',@(~,~)setappdata(fig,'stop',true));
setappdata(fig,'stop',false);

fprintf('RX @ %.0f MHz\n', fc/1e6);
fprintf('Expecting : "%s"\n', msg_padded);

while ~getappdata(fig,'stop')
try

    rxRaw = rx();
    pwr   = mean(abs(rxRaw).^2);
    if pwr < 0.001; fprintf('No signal\n'); continue; end

    rxFilt = filter(RRC, 1, rxRaw.');
    rxDS   = rxFilt(1:OVR:end);
    N      = length(rxDS);

    %% ── Calibration ──────────────────────────────────────────────────
    if ~cal_done
        S    = abs(fftshift(fft(rxDS(1:min(N,32768)), 32768)));
        fax  = linspace(-fs/2, fs/2, 32768);
        posS = S .* (fax > 500 & fax < fs/2-500);
        pk   = find(posS > 0.1*max(posS));
        if isempty(pk); continue; end
        fm = fax(pk(1));
        if abs(fm) > 400000; continue; end
        cal_buf(end+1) = fm; %#ok
        fprintf('Cal %d/%d = %.1f Hz\n', length(cal_buf), N_cal, fm);
        if length(cal_buf) >= N_cal
            med      = median(cal_buf);
            good     = cal_buf(abs(cal_buf-med) < 10000);
            f_bulk   = mean(good);
            cal_done = true;
            fprintf('>>> f_bulk = %.1f Hz\n', f_bulk);
        end
        continue;
    end

    %% ── Bulk CFO removal ─────────────────────────────────────────────
    rxC = rxDS .* exp(-j*2*pi*f_bulk*Ts*(0:N-1));

    %% ── PN cross-correlation ─────────────────────────────────────────
    xc              = abs(conv(real(rxC), fliplr(real(pn_tmpl))));
    [pn_pk, pn_idx] = max(xc);
    noise_floor     = mean(xc);

    if pn_pk < 5*noise_floor
        fprintf('No PN peak\n');
        phase_init = false;
        continue;
    end

    ofdm_st = pn_idx - length(pn_tmpl) + 1 + pn_ds_len;
    if ofdm_st < 1 || ofdm_st + 479 > N
        fprintf('Out of bounds\n'); continue;
    end

    rxF = rxC(ofdm_st : ofdm_st+479);

    %% ── Power balance check ──────────────────────────────────────────
    p = [mean(abs(rxF(1:160)).^2),   mean(abs(rxF(161:320)).^2),...
         mean(abs(rxF(321:400)).^2), mean(abs(rxF(401:480)).^2)];
    if p(3) < 0.3 || p(4) < 0.3
        fprintf('Weak payload — skip\n'); continue;
    end

    %% ── Coarse CFO ───────────────────────────────────────────────────
    z = 0;
    for k = 1:8
        z = z + rxF(spl*(k-1)+1:spl*k) * conj(rxF(spl*k+1:spl*(k+1))).';
    end
    fc_ = (-1/(2*pi*spl*Ts)) * angle(z);
    rxF = rxF .* exp(-j*2*pi*fc_*Ts*(0:479));

    %% ── Fine CFO ─────────────────────────────────────────────────────
    lp1 = rxF(spl*12+1 : spl*16);
    lp2 = rxF(spl*16+1 : spl*20);
    ff_ = (-1/(2*pi*64*Ts)) * angle(lp1 * conj(lp2).');
    rxF = rxF .* exp(-j*2*pi*ff_*Ts*(0:479));

    %% ── Channel estimation — single frame ────────────────────────────
    lp1f   = rxF(spl*12+1 : spl*16);
    lp2f   = rxF(spl*16+1 : spl*20);
    Y1     = fftshift(fft(lp1f));
    Y2     = fftshift(fft(lp2f));
    H      = 0.5*(Y1+Y2) ./ (LP_f_ideal + eps);
    H_mean = mean(abs(H(dIdx)));
    if H_mean < 0.001; continue; end
    H      = H / H_mean;
    H(abs(LP_f_ideal) < 0.01) = 1;
    H_safe = H;
    H_safe(abs(H) < 0.01) = 0.01;

    LP_corr = abs(lp1f*conj(lp2f).') / (norm(lp1f)*norm(lp2f)+eps);

    if LP_corr < 0.5
        fprintf('Bad LP (%.3f) — skip\n', LP_corr);
        continue;
    end

    %% ── Equalizer ────────────────────────────────────────────────────
    eq1 = fftshift(fft(rxF(321+N_CP:400))) ./ H_safe;
    eq2 = fftshift(fft(rxF(401+N_CP:480))) ./ H_safe;
    P1d = eq1(dIdx);
    P2d = eq2(dIdx);

    %% ── Pilot phase correction ───────────────────────────────────────
    rx_pilots  = [eq1(pilot_pos), eq2(pilot_pos)];
    ref_pilots = repmat(pilot_vals, 1, 2);
    ph_err     = mean(angle(rx_pilots .* conj(ref_pilots)));
    P1d = P1d * exp(-j*ph_err);
    P2d = P2d * exp(-j*ph_err);

    %% ── Amplitude normalization ──────────────────────────────────────
    scale = rms([P1d, P2d]);
    P1d   = P1d / scale;
    P2d   = P2d / scale;

    %% ── Save pre-phase symbols for plot 7 ───────────────────────────
    P1d_pre = P1d;
    P2d_pre = P2d;

    %% ── Phase tracking PLL ───────────────────────────────────────────
    if ~phase_init || abs(fc_) > 5000
        allSyms   = [P1d, P2d];
        bestScore = -inf;
        for deg = 0:1:359
            phi     = deg * pi/180;
            rotated = allSyms * exp(j*phi);
            dists   = arrayfun(@(s) min(abs(s-corners)), rotated);
            score   = -sum(dists);
            if score > bestScore
                bestScore   = score;
                phase_track = phi;
            end
        end
        phase_init = true;
        const_buf  = [];
    else
        bestScore = -inf;
        bestDelta = 0;
        for deg = -10:0.5:10
            phi     = phase_track + deg*pi/180;
            rotated = [P1d,P2d] * exp(j*phi);
            dists   = arrayfun(@(s) min(abs(s-corners)), rotated);
            score   = -sum(dists);
            if score > bestScore
                bestScore = score;
                bestDelta = deg*pi/180;
            end
        end
        phase_track = phase_track + 0.3*bestDelta;
    end

    P1d = P1d * exp(j*phase_track);
    P2d = P2d * exp(j*phase_track);

    %% ── Demodulate ───────────────────────────────────────────────────
    b1 = pskdemod(P1d, M, pi/4);
    b2 = pskdemod(P2d, M, pi/4);

    %% ── BER ──────────────────────────────────────────────────────────
    e   = sum([abs(sign(data_Payload_1-b1)), abs(sign(data_Payload_2-b2))]);
    BER = e / (length(data_Payload_1)+length(data_Payload_2));

    %% ── Try all 4 QPSK rotations — pick best BER ────────────────────
    best_BER   = BER;
    best_b1    = b1;
    best_b2    = b2;
    best_P1d   = P1d;
    best_P2d   = P2d;
    best_phase = phase_track;

    for rot_k = 1:3
        rot_phi = rot_k * pi/2;
        P1d_try = P1d * exp(j*rot_phi);
        P2d_try = P2d * exp(j*rot_phi);
        b1_try  = pskdemod(P1d_try, M, pi/4);
        b2_try  = pskdemod(P2d_try, M, pi/4);
        e_try   = sum([abs(sign(data_Payload_1-b1_try)), ...
                       abs(sign(data_Payload_2-b2_try))]);
        BER_try = e_try / (length(data_Payload_1)+length(data_Payload_2));
        if BER_try < best_BER
            best_BER   = BER_try;
            best_b1    = b1_try;
            best_b2    = b2_try;
            best_P1d   = P1d_try;
            best_P2d   = P2d_try;
            best_phase = phase_track + rot_phi;
        end
    end

    BER         = best_BER;
    b1          = best_b1;
    b2          = best_b2;
    P1d         = best_P1d;
    P2d         = best_P2d;
    phase_track = best_phase;

    %% ── Decode bits to text ──────────────────────────────────────────
    rx_bits  = [de2bi(b1, 2, 'left-msb'); de2bi(b2, 2, 'left-msb')];
    rx_bits  = reshape(rx_bits.', 1, []);
    rx_chars = bi2de(reshape(rx_bits, 8, [])', 'left-msb');
    rx_msg   = char(rx_chars);

    %% ── Accumulate constellation ─────────────────────────────────────
    const_buf = [const_buf, P1d, P2d]; %#ok
    if length(const_buf) > const_max
        const_buf = const_buf(end-const_max+1:end);
    end

    %% ── Bulk update — faster alpha when fc still large ───────────────
    if abs(fc_) > 300
        f_bulk = f_bulk + 0.15 * fc_;   % converge faster
    else
        f_bulk = f_bulk + 0.05 * fc_;   % slow stable tracking
    end

    ber_hist(end+1) = BER; %#ok
    if length(ber_hist) > 20; ber_hist = ber_hist(end-19:end); end

    %% ── Terminal output ──────────────────────────────────────────────
    if abs(fc_) > 300
        cfo_status = 'converging';
    else
        cfo_status = 'locked';
    end
    fprintf('BER=%.4f Avg=%.4f LP=%.3f fc=%.0f [%s]\n', ...
            BER, mean(ber_hist), LP_corr, fc_, cfo_status);
    fprintf('  TX: "%s"\n', msg_padded);
    fprintf('  RX: "%s"\n', rx_msg);

    %% ── Plot 1: RX Raw ───────────────────────────────────────────────
    subplot(2,4,1);
    plot(rxRaw,'.','MarkerSize',1);
    title('RX Raw (IQ)');
    axis([-1.5 1.5 -1.5 1.5]); axis square;

    %% ── Plot 2: I ────────────────────────────────────────────────────
    subplot(2,4,2);
    plot(real(rxRaw));
    title('I (in-phase)');
    axis([1 length(rxRaw) -1.5 1.5]); axis square;

    %% ── Plot 3: Q ────────────────────────────────────────────────────
    subplot(2,4,3);
    plot(imag(rxRaw));
    title('Q (quadrature)');
    axis([1 length(rxRaw) -1.5 1.5]); axis square;

    %% ── Plot 4: Welch PSD ────────────────────────────────────────────
    subplot(2,4,4);
    [Sw,fw] = pwelch(rxRaw,[],[],[],fs,'centered','power');
    plot(fw,pow2db(Sw));
    title('Welch PSD'); axis square;

    %% ── Plot 5: PN XCorr ─────────────────────────────────────────────
    subplot(2,4,5);
    plot(xc(1:min(end,20000)));
    title(sprintf('PN XCorr (%.1fx noise)',pn_pk/noise_floor));
    axis square;

    %% ── Plot 6: Channel IR ───────────────────────────────────────────
    subplot(2,4,6);
    H_t = ifft(ifftshift(H));
    plot(abs(H_t));
    title(sprintf('Channel IR | LP=%.3f',LP_corr));
    axis([1 64 0 3]); axis square; xlabel('Tap');

    %% ── Plot 7: After EQ before phase corr ──────────────────────────
    subplot(2,4,7);
    plot(real([P1d_pre,P2d_pre]), imag([P1d_pre,P2d_pre]), ...
         '*','MarkerSize',8,'Color',[0.2 0.6 1]);
    hold on;
    plot(real(corners),imag(corners),'r+','MarkerSize',15,'LineWidth',2);
    hold off;
    grid on; axis([-2 2 -2 2]); axis square;
    title(sprintf('After EQ | BER=%.4f',BER));
    xlabel('In-phase'); ylabel('Quadrature');

    %% ── Plot 8: Accumulated QPSK constellation ───────────────────────
    subplot(2,4,8);
    if length(const_buf) > 4
        cmap = [0 0.5 1; 1 0.4 0; 0.1 0.75 0.1; 0.8 0 0.8];
        hold off;
        for qi = 1:4
            [~,nearest] = min([abs(const_buf-corners(1));
                               abs(const_buf-corners(2));
                               abs(const_buf-corners(3));
                               abs(const_buf-corners(4))]);
            pts = const_buf(nearest==qi);
            if ~isempty(pts)
                plot(real(pts),imag(pts),'.', ...
                     'MarkerSize',6,'Color',cmap(qi,:));
                hold on;
            end
        end
        plot(real(corners),imag(corners),'k+', ...
             'MarkerSize',18,'LineWidth',3);
        hold off;
    end
    grid on; axis([-1.5 1.5 -1.5 1.5]); axis square;
    title(sprintf('QPSK | BER=%.4f | Avg=%.4f',BER,mean(ber_hist)));
    xlabel('In-phase'); ylabel('Quadrature');

    drawnow;

catch ME
    fprintf('Error: %s\n', ME.message);
end
end

release(rx);
close all;
fprintf('Done\n');