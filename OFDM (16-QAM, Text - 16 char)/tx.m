%% ================================================
%% OFDM_TX.m — PC 1, USRP B200 — 16-QAM
%% ================================================
clear; close all; clc; j = 1i;
load('OFDM_config.mat');

N_FFT   = 64;
N_CP    = 16;
OVR     = 2;
rolloff = 0.5;
span    = 6;
M       = 16;   % 16-QAM
pilot   = [1,1,1,-1];
vsc11   = zeros(1,11);
RRC     = rcosdesign(rolloff, span, OVR, 'sqrt');

%% ── Short preamble ───────────────────────────────────────────────────
S_k = sqrt(13/6)*[0,0,1+j,0,0,0,-1-j,0,0,0,1+j,0,0,0,-1-j,0,0,0,-1-j,...
      0,0,0,1+j,0,0,0,0,0,0,0,-1-j,0,0,0,-1-j,0,0,0,1+j,0,0,0,1+j,...
      0,0,0,1+j,0,0,0,1+j,0,0];
vsc  = zeros(1, N_FFT-length(S_k));
SP_f = [vsc(1:6), S_k, vsc(7:11)];
SP_t = ifft(ifftshift(SP_f));
SP   = repmat(SP_t(1:N_CP), 1, 10);

%% ── Long preamble ────────────────────────────────────────────────────
LP_t = ifft(ifftshift(LP_freq));
LP   = [LP_t(33:64), LP_t, LP_t];

%% ── Payload builder — 16-QAM ─────────────────────────────────────────
function sym = make_ofdm_sym(data, pilot, vsc11)
    %% qammod with UnitAveragePower — normalizes constellation to power=1
    d = qammod(data, 16, 'UnitAveragePower', true);
    f = [vsc11(1:6),d(1:5),pilot(1),d(6:18),pilot(2),...
         d(19:24),0,d(25:30),pilot(3),d(31:43),pilot(4),...
         d(44:48),vsc11(7:11)];
    t   = ifft(ifftshift(f));
    sym = [t(49:64), t];
end

P1 = make_ofdm_sym(data_Payload_1, pilot, vsc11);
P2 = make_ofdm_sym(data_Payload_2, pilot, vsc11);
P3 = make_ofdm_sym(data_Payload_3, pilot, vsc11);
P4 = make_ofdm_sym(data_Payload_4, pilot, vsc11);
P5 = make_ofdm_sym(data_Payload_5, pilot, vsc11);
P6 = make_ofdm_sym(data_Payload_6, pilot, vsc11);
P7 = make_ofdm_sym(data_Payload_7, pilot, vsc11);
P8 = make_ofdm_sym(data_Payload_8, pilot, vsc11);

ofdm_frame = [SP, LP, P1, P2, P3, P4, P5, P6, P7, P8];
fprintf('Frame = %d samples\n', length(ofdm_frame));

pn_up  = upsample(pn_seq, OVR);
pn_rrc = filter(RRC, 1, pn_up);
pn_rrc = pn_rrc / max(abs(pn_rrc));

ofdm_up  = upsample(ofdm_frame, OVR);
ofdm_rrc = filter(RRC, 1, ofdm_up);
ofdm_rrc = ofdm_rrc / max(abs(ofdm_rrc));

burst = [pn_rrc, ofdm_rrc];
burst = burst / max(abs(burst));

pn_rrc_len = length(pn_rrc);
save('pn_len.mat', 'pn_rrc', 'pn_rrc_len');
fprintf('pn_len.mat saved\n');
fprintf('Transmitting : "%s"\n', msg_padded);

gap   = zeros(1, 500);
txOut = [zeros(1,5000), repmat([burst,gap],1,20), zeros(1,5000)].';

tx = comm.SDRuTransmitter('Platform','B200','SerialNum','34D8D87',...
        'CenterFrequency',fc,'MasterClockRate',20e6,...
        'InterpolationFactor',20,'Gain',89);

fig = figure('Name','TX');
uicontrol('String','Stop','Position',[80 50 100 60],...
          'Callback',@(~,~)setappdata(fig,'stop',true));
setappdata(fig,'stop',false);
set(fig,'Units','centimeters','position',[3 3 7 6]);

fprintf('TX @ %.0f MHz\n', fc/1e6);
while ~getappdata(fig,'stop')
    tx(txOut); drawnow;
end
release(tx);
fprintf('TX stopped\n');