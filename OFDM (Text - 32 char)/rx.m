%% ================================================
%% OFDM_RX.m — PC 2, RTL-SDR — 32 CHARS INTERLEAVED REP-3
%% ================================================
clear; close all; clc; j = 1i;
load('OFDM_config.mat');
load('pn_len.mat');

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
frameSamps = 960;

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

    %% ── Save pre-phase for plot 7 ────────────────────────────────────
    Pd_pre = Pd;

    %% ── Phase tracking PLL ───────────────────────────────────────────
    allSyms = reshape(Pd.', 1, []);
    if ~phase_init || abs(fc_) > 5000
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
        bestScore = -inf; bestDelta = 0;
        for deg = -20:0.5:20
            phi     = phase_track + deg*pi/180;
            rotated = allSyms * exp(j*phi);
            dists   = arrayfun(@(s) min(abs(s-corners)), rotated);
            score   = -sum(dists);
            if score > bestScore
                bestScore = score; bestDelta = deg*pi/180;
            end
        end
        phase_track = phase_track + 0.3*bestDelta;
    end
    Pd = Pd * exp(j*phase_track);

    %% ── Demodulate ───────────────────────────────────────────────────
    B = zeros(8, 48);
    for pi_ = 1:8
        B(pi_,:) = pskdemod(Pd(pi_,:), M, pi/4);
    end
    b_all = reshape(B.', 1, []);   % [1x384] row-major

    %% ── BER at symbol level ──────────────────────────────────────────
    e_sym   = sum(abs(sign(data_all - b_all)));
    BER_sym = e_sym / length(data_all);

    %% ── 4-rotation fix ───────────────────────────────────────────────
    best_BER = BER_sym;
    best_B   = B;
    best_Pd  = Pd;
    best_ph  = phase_track;
    for rk = 0:3
        rp   = rk*pi/2;
        Pd_t = Pd * exp(j*rp);
        B_t  = zeros(8,48);
        for pi_=1:8; B_t(pi_,:) = pskdemod(Pd_t(pi_,:),M,pi/4); end
        b_t  = reshape(B_t.', 1, []);
        et   = sum(abs(sign(data_all - b_t)));
        bt   = et/length(data_all);
        if bt < best_BER
            best_BER=bt; best_B=B_t; best_Pd=Pd_t;
            best_ph=phase_track+rp;
        end
    end
    BER_sym     = best_BER;
    B           = best_B;
    Pd          = best_Pd;
    b_all       = reshape(B.', 1, []);
    phase_track = best_ph;

    %% ── Decision-directed residual phase correction ──────────────────
    ideal_syms  = pskmod(b_all, M, pi/4);
    recv_syms   = reshape(Pd.', 1, []);
    resid_ph    = mean(angle(recv_syms .* conj(ideal_syms)));
    Pd          = Pd * exp(-j*resid_ph);
    phase_track = phase_track - resid_ph;

    %% ── Get received bits ────────────────────────────────────────────
    rx_bits = de2bi(b_all, 2, 'left-msb');   % [384x2]
    rx_bits = reshape(rx_bits.', 1, []);      % [1x768]

    %% ── INTERLEAVED majority vote decode ─────────────────────────────
    %% TX sent: [copy1_b1..b256 | copy2_b1..b256 | copy3_b1..b256]
    %% So the 3 copies of bit i are at positions i, i+256, i+512
    %% Bad subcarrier only corrupts 1 copy — majority vote always wins
    n_orig   = 256;
    bits_dec = zeros(1, n_orig);
    for bi = 1:n_orig
        c1 = rx_bits(bi);           % copy 1 of bit bi
        c2 = rx_bits(bi + 256);     % copy 2 of bit bi — 256 apart
        c3 = rx_bits(bi + 512);     % copy 3 of bit bi — 512 apart
        bits_dec(bi) = round((c1 + c2 + c3) / 3);
    end

    e_rep   = sum(abs(bits_dec - bits_orig));
    BER_rep = e_rep / length(bits_orig);

    rx_chars = bi2de(reshape(bits_dec, 8, [])', 'left-msb');
    rx_msg   = char(rx_chars);

    %% ── Accumulate constellation ─────────────────────────────────────
    if BER_rep == 0 && abs(fc_) < 100
        pts_new   = reshape(Pd.', 1, []);
        const_buf = [const_buf, pts_new]; %#ok
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
    fprintf('  RX: "%s"\n', rx_msg);

    %% ── Plot 1 ───────────────────────────────────────────────────────
    subplot(2,4,1);
    plot(rxRaw,'.','MarkerSize',1);
    title('RX Raw (IQ)'); axis([-1.5 1.5 -1.5 1.5]); axis square;

    %% ── Plot 2 ───────────────────────────────────────────────────────
    subplot(2,4,2);
    plot(real(rxRaw)); title('I (in-phase)');
    axis([1 length(rxRaw) -1.5 1.5]); axis square;

    %% ── Plot 3 ───────────────────────────────────────────────────────
    subplot(2,4,3);
    plot(imag(rxRaw)); title('Q (quadrature)');
    axis([1 length(rxRaw) -1.5 1.5]); axis square;

    %% ── Plot 4 ───────────────────────────────────────────────────────
    subplot(2,4,4);
    [Sw,fw] = pwelch(rxRaw,[],[],[],fs,'centered','power');
    plot(fw,pow2db(Sw)); title('Welch PSD'); axis square;

    %% ── Plot 5 ───────────────────────────────────────────────────────
    subplot(2,4,5);
    plot(xc(1:min(end,20000)));
    title(sprintf('PN XCorr (%.1fx noise)',pn_pk/noise_floor)); axis square;

    %% ── Plot 6 ───────────────────────────────────────────────────────
    subplot(2,4,6);
    H_t = ifft(ifftshift(H)); plot(abs(H_t));
    title(sprintf('Channel IR | LP=%.3f',LP_corr));
    axis([1 64 0 3]); axis square; xlabel('Tap');

    %% ── Plot 7 ───────────────────────────────────────────────────────
    subplot(2,4,7);
    pre_pts = reshape(Pd_pre.', 1, []);
    plot(real(pre_pts),imag(pre_pts),'*','MarkerSize',5,'Color',[0.2 0.6 1]);
    hold on;
    plot(real(corners),imag(corners),'r+','MarkerSize',15,'LineWidth',2);
    hold off;
    grid on; axis([-2 2 -2 2]); axis square;
    title(sprintf('After EQ | BER-sym=%.4f',BER_sym));
    xlabel('In-phase'); ylabel('Quadrature');

    %% ── Plot 8 ───────────────────────────────────────────────────────
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
                plot(real(pts),imag(pts),'.','MarkerSize',4,'Color',cmap(qi,:));
                hold on;
            end
        end
        plot(real(corners),imag(corners),'k+','MarkerSize',20,'LineWidth',3);
        theta = linspace(0,2*pi,100);
        for qi=1:4
            plot(real(corners(qi))+0.4*cos(theta),...
                 imag(corners(qi))+0.4*sin(theta),...
                 '--','Color',[0.5 0.5 0.5],'LineWidth',0.5);
        end
        hold off;
    else
        text(0,0,'Accumulating...','HorizontalAlignment','center',...
             'Color','white','FontSize',12);
    end
    grid on; axis([-1.5 1.5 -1.5 1.5]); axis square;
    title(sprintf('QPSK | BER-sym=%.4f | BER-rep=%.4f',...
                  BER_sym,BER_rep),'FontSize',8);
    xlabel(sprintf('TX: "%s"',msg_padded),'FontSize',7);
    ylabel(sprintf('RX: "%s" | %d pts',...
                   strrep(rx_msg,'\','\\'),length(const_buf)),'FontSize',7);

    drawnow;

catch ME
    fprintf('Error: %s\n', ME.message);
end
end

release(rx);
close all;
fprintf('Done\n');
