%% ================================================
%% OFDM_RX.m — 16-QAM 16 CHARS SOFT COMBINING x12
%% ================================================
clear; close all; clc; j = 1i;
load('OFDM_config.mat');
load('pn_len.mat');

N_FFT      = 64;
N_CP       = 16;
OVR        = 2;
rolloff    = 0.5;
span       = 6;
M          = 16;
Ts         = 1/fs;
dIdx       = [7:11, 13:25, 27:32, 34:39, 41:53, 55:59];
pilot_pos  = [12, 26, 40, 54];
pilot_vals = [1, 1, 1, -1];
spl        = 16;
RRC        = rcosdesign(rolloff, span, OVR, 'sqrt');
frameSamps = 960;
N_copies   = 12;   % 12 copies for soft combining
N_syms     = 32;   % symbols per copy (128 bits / 4 bits per sym)

%% ── 16-QAM constellation ─────────────────────────────────────────────
r           = [-3,-1,1,3] / sqrt(10);
[RI,RQ]     = meshgrid(r,r);
qam16_const = RI(:).' + 1j*RQ(:).';

%% ── PN template ──────────────────────────────────────────────────────
pn_rrc_ds = pn_rrc(1:OVR:end);
pn_tmpl   = pn_rrc_ds / norm(pn_rrc_ds);
pn_ds_len = length(pn_rrc_ds);

%% ── Ideal LP ─────────────────────────────────────────────────────────
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
const_max   = 1000;

data_all = [data_Payload_1, data_Payload_2, data_Payload_3, data_Payload_4,...
            data_Payload_5, data_Payload_6, data_Payload_7, data_Payload_8];

rx = comm.SDRRTLReceiver('CenterFrequency',fc,'SampleRate',fs,...
        'SamplesPerFrame',32768,'OutputDataType','double',...
        'EnableTunerAGC',false,'TunerGain',49);

fig = figure('Name','RX 16-QAM','NumberTitle','off');
set(fig,'Units','centimeters','position',[1 2 49 24]);
uicontrol('String','Stop','Position',[1475 15 100 60],...
          'Callback',@(~,~)setappdata(fig,'stop',true));
setappdata(fig,'stop',false);

fprintf('RX @ %.0f MHz  [16-QAM 16-char 12x soft combining]\n', fc/1e6);
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
        fprintf('No PN peak\n'); phase_init = false; continue;
    end

    ofdm_st = pn_idx - length(pn_tmpl) + 1 + pn_ds_len;
    if ofdm_st < 1 || ofdm_st + frameSamps - 1 > N
        fprintf('Out of bounds\n'); continue;
    end

    rxF = rxC(ofdm_st : ofdm_st + frameSamps - 1);

    %% ── Power balance ────────────────────────────────────────────────
    bad = false;
    for pi_ = 1:8
        seg_s = 320 + (pi_-1)*80 + 1;
        seg_e = 320 + pi_*80;
        if mean(abs(rxF(seg_s:seg_e)).^2) < 0.3
            bad = true; break;
        end
    end
    if bad; fprintf('Weak payload — skip\n'); continue; end

    %% ── Coarse CFO ───────────────────────────────────────────────────
    z = 0;
    for k = 1:8
        z = z + rxF(spl*(k-1)+1:spl*k) * conj(rxF(spl*k+1:spl*(k+1))).';
    end
    fc_ = (-1/(2*pi*spl*Ts)) * angle(z);
    rxF = rxF .* exp(-j*2*pi*fc_*Ts*(0:frameSamps-1));

    %% ── Fine CFO ─────────────────────────────────────────────────────
    lp1 = rxF(spl*12+1 : spl*16);
    lp2 = rxF(spl*16+1 : spl*20);
    ff_ = (-1/(2*pi*64*Ts)) * angle(lp1 * conj(lp2).');
    rxF = rxF .* exp(-j*2*pi*ff_*Ts*(0:frameSamps-1));

    %% ── Channel estimation ───────────────────────────────────────────
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
        fprintf('Bad LP (%.3f) — skip\n', LP_corr); continue;
    end

    %% ── Equalize all 8 symbols ───────────────────────────────────────
    eq_all = zeros(8, 64);
    Pd     = zeros(8, 48);
    for pi_ = 1:8
        seg_s         = 320 + (pi_-1)*80 + N_CP + 1;
        seg_e         = 320 + pi_*80;
        eq_sym        = fftshift(fft(rxF(seg_s:seg_e))) ./ H_safe;
        eq_all(pi_,:) = eq_sym;
        Pd(pi_,:)     = eq_sym(dIdx);
    end

    %% ── Per-symbol pilot phase correction ────────────────────────────
    for pi_ = 1:8
        sym_pilots    = eq_all(pi_, pilot_pos);
        ph_sym        = mean(angle(sym_pilots .* conj(pilot_vals)));
        Pd(pi_,:)     = Pd(pi_,:) * exp(-j*ph_sym);
        eq_all(pi_,:) = eq_all(pi_,:) * exp(-j*ph_sym);
    end

    %% ── Amplitude normalization ──────────────────────────────────────
    scale = rms(Pd(:));
    Pd    = Pd / scale;

    Pd_pre = Pd;

    %% ── Phase tracking PLL ───────────────────────────────────────────
    allSyms = reshape(Pd.', 1, []);

    if ~phase_init || abs(fc_) > 5000
        bestScore = -inf;
        for deg = 0:2:358
            phi      = deg * pi/180;
            rotated  = allSyms * exp(j*phi);
            dist_mat = abs(repmat(rotated.', 1, 16) - ...
                           repmat(qam16_const, length(rotated), 1));
            score    = -sum(min(dist_mat, [], 2));
            if score > bestScore
                bestScore   = score;
                phase_track = phi;
            end
        end
        phase_init = true;
        const_buf  = [];
    else
        bestScore = -inf; bestDelta = 0;
        for deg = -15:0.5:15
            phi      = phase_track + deg*pi/180;
            rotated  = allSyms * exp(j*phi);
            dist_mat = abs(repmat(rotated.', 1, 16) - ...
                           repmat(qam16_const, length(rotated), 1));
            score    = -sum(min(dist_mat, [], 2));
            if score > bestScore
                bestScore = score; bestDelta = deg*pi/180;
            end
        end
        phase_track = phase_track + 0.3*bestDelta;
    end
    Pd = Pd * exp(j*phase_track);

    %% ── Hard demod for BER and 4-rotation fix ────────────────────────
    B = zeros(8, 48);
    for pi_ = 1:8
        B(pi_,:) = qamdemod(Pd(pi_,:), M, 'UnitAveragePower', true);
    end
    b_all = reshape(B.', 1, []);

    ref_bits    = de2bi(data_all, 4, 'left-msb');
    ref_bits    = reshape(ref_bits.', 1, []);
    rx_bits_all = de2bi(b_all, 4, 'left-msb');
    rx_bits_all = reshape(rx_bits_all.', 1, []);
    BER_sym     = sum(abs(rx_bits_all - ref_bits)) / length(ref_bits);

    %% ── 4-rotation fix ───────────────────────────────────────────────
    best_BER = BER_sym;
    best_Pd  = Pd;
    best_ph  = phase_track;
    for rk = 0:3
        rp   = rk*pi/2;
        Pd_t = Pd * exp(j*rp);
        B_t  = zeros(8,48);
        for pi_=1:8
            B_t(pi_,:) = qamdemod(Pd_t(pi_,:), M, 'UnitAveragePower', true);
        end
        b_t  = reshape(B_t.', 1, []);
        rb_t = de2bi(b_t, 4, 'left-msb');
        rb_t = reshape(rb_t.', 1, []);
        bt   = sum(abs(rb_t - ref_bits)) / length(ref_bits);
        if bt < best_BER
            best_BER=bt; best_Pd=Pd_t; best_ph=phase_track+rp;
        end
    end
    BER_sym     = best_BER;
    Pd          = best_Pd;
    phase_track = best_ph;

    %% ── Decision-directed residual phase correction ──────────────────
    B2 = zeros(8,48);
    for pi_=1:8
        B2(pi_,:) = qamdemod(Pd(pi_,:), M, 'UnitAveragePower', true);
    end
    b2          = reshape(B2.', 1, []);
    ideal_syms  = qammod(b2, M, 'UnitAveragePower', true);
    recv_syms   = reshape(Pd.', 1, []);
    resid_ph    = mean(angle(recv_syms .* conj(ideal_syms)));
    Pd          = Pd * exp(-j*resid_ph);
    phase_track = phase_track - resid_ph;

    %% ── SOFT COMBINING: average 12 copies of each symbol ────────────
    %% This is the key improvement over hard majority vote
    %% Averaging N complex symbols reduces noise by sqrt(N) = sqrt(12) = 10.8 dB
    %% Each copy is 32 symbols — interleaved 256 bits apart in the bit stream
    Pd_flat  = reshape(Pd.', 1, []);   % [1x384] row-major
    Pd_soft  = zeros(1, N_syms);       % [1x32] averaged symbols
    for ci = 1:N_copies
        Pd_soft = Pd_soft + Pd_flat((ci-1)*N_syms+1 : ci*N_syms);
    end
    Pd_soft = Pd_soft / N_copies;      % complex average — maximally combines SNR

    %% ── Demodulate soft-combined symbols ─────────────────────────────
    b_soft   = qamdemod(Pd_soft, M, 'UnitAveragePower', true);   % [1x32]

    %% ── Decode to text ───────────────────────────────────────────────
    rx_bits  = de2bi(b_soft, 4, 'left-msb');   % [32x4]
    rx_bits  = reshape(rx_bits.', 1, []);       % [1x128]

    e_rep    = sum(abs(rx_bits - bits_orig));
    BER_rep  = e_rep / length(bits_orig);

    rx_chars = bi2de(reshape(rx_bits, 8, [])', 'left-msb');
    rx_msg   = char(rx_chars(:).');
    safe_msg = regexprep(rx_msg, '[^\x20-\x7E]', '?');

    %% ── Accumulate constellation ─────────────────────────────────────
    if BER_rep == 0 && abs(fc_) < 100
        const_buf = [const_buf, Pd_soft]; %#ok
        if length(const_buf) > const_max
            const_buf = const_buf(end-const_max+1:end);
        end
    end

    %% ── Adaptive bulk update ─────────────────────────────────────────
    if abs(fc_) > 500;      f_bulk = f_bulk + 0.3*fc_;
    elseif abs(fc_) > 100;  f_bulk = f_bulk + 0.15*fc_;
    else;                   f_bulk = f_bulk + 0.05*fc_;
    end

    ber_hist(end+1) = BER_rep; %#ok
    if length(ber_hist) > 20; ber_hist = ber_hist(end-19:end); end

    if abs(fc_) > 100; cfo_st='converging'; else; cfo_st='locked'; end

    fprintf('BER-sym=%.4f | BER-rep=%.4f | Avg=%.4f | LP=%.3f | fc=%.0f [%s]\n',...
            BER_sym, BER_rep, mean(ber_hist), LP_corr, fc_, cfo_st);
    fprintf('  TX: "%s"\n', msg_padded);
    fprintf('  RX: "%s"\n', safe_msg);

    %% ── Plots ────────────────────────────────────────────────────────
    subplot(2,4,1);
    plot(rxRaw,'.','MarkerSize',1);
    title('RX Raw (IQ)'); axis([-1.5 1.5 -1.5 1.5]); axis square;

    subplot(2,4,2);
    plot(real(rxRaw)); title('I (in-phase)');
    axis([1 length(rxRaw) -1.5 1.5]); axis square;

    subplot(2,4,3);
    plot(imag(rxRaw)); title('Q (quadrature)');
    axis([1 length(rxRaw) -1.5 1.5]); axis square;

    subplot(2,4,4);
    [Sw,fw] = pwelch(rxRaw,[],[],[],fs,'centered','power');
    plot(fw,pow2db(Sw)); title('Welch PSD'); axis square;

    subplot(2,4,5);
    plot(xc(1:min(end,20000)));
    title(sprintf('PN XCorr (%.1fx noise)',pn_pk/noise_floor)); axis square;

    subplot(2,4,6);
    H_t = ifft(ifftshift(H)); plot(abs(H_t));
    title(sprintf('Channel IR | LP=%.3f',LP_corr));
    axis([1 64 0 3]); axis square; xlabel('Tap');

    subplot(2,4,7);
    pre_pts = reshape(Pd_pre.', 1, []);
    plot(real(pre_pts),imag(pre_pts),'.','MarkerSize',4,'Color',[0.2 0.6 1]);
    hold on;
    plot(real(qam16_const),imag(qam16_const),'r+','MarkerSize',10,'LineWidth',2);
    hold off;
    grid on; axis([-2 2 -2 2]); axis square;
    title(sprintf('After EQ (all 384 syms) | BER-sym=%.4f',BER_sym));
    xlabel('In-phase'); ylabel('Quadrature');

    subplot(2,4,8);
    if length(const_buf) > 4
        plot(real(const_buf),imag(const_buf),'.','MarkerSize',6,...
             'Color',[0.2 0.7 1]);
        hold on;
        plot(real(qam16_const),imag(qam16_const),'k+','MarkerSize',15,'LineWidth',2);
        hold off;
    else
        text(0,0,'Accumulating...','HorizontalAlignment','center',...
             'Color','white','FontSize',12);
    end
    grid on; axis([-1.5 1.5 -1.5 1.5]); axis square;
    title(sprintf('16-QAM soft x12 | BER-rep=%.4f | Avg=%.4f',...
                  BER_rep, mean(ber_hist)),'FontSize',8);
    xlabel(sprintf('TX: "%s"', msg_padded),'FontSize',8);
    ylabel(sprintf('RX: "%s" | %d pts', safe_msg, ...
                   length(const_buf)),'FontSize',8);

    drawnow;

catch ME
    fprintf('Error: %s\n', ME.message);
end
end

release(rx);
close all;
fprintf('Done\n');