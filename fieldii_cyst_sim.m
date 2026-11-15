%% FIELD II SIMULATION: focused-transmit anechoic cyst phantom
% Output: fieldII_cyst_RF.mat (read by das_beamform_cyst.m)
% Requires Field II and the Signal Processing Toolbox (resample, hilbert).
%
% Timing convention used by both scripts:
%   Transmit delays are set explicitly with xdc_focus_times so that the
%   earliest element fires at t = 0 in every transmit event:
%       delay_i = (r_max - r_i)/c,   r_i = |element_i - focal point|
%   A single point scatterer at the focus measures the fixed pulse lag.
%   All transmit events are placed on one common time axis at the
%   high simulation rate (before decimation), so no sub-sample rounding
%   error is introduced between transmits.

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
QUICK = false;            % true: much faster, slightly lower quality (use to test the pipeline)

f0     = 2.717391e6;      % centre frequency (Hz)
fs     = 10.85914e6;      % stored sampling frequency (Hz)
c      = 1540;            % speed of sound (m/s)
lambda = c/f0;

N_elements     = 64;
element_width  = 0.25e-3;
element_height = 3e-3;
pitch          = 0.30e-3;             % placeholder, replace with the real probe pitch
kerf = pitch - element_width;
assert(kerf >= 0, 'Pitch must be >= element width.');

z_focus  = 60e-3;         % transmit focal depth
tx_x_max = 11e-3;         % transmit foci span -tx_x_max ... +tx_x_max

SNR_dB        = -15;      % default saved noisy RF; DAS script can override for SNR sweeps
tail_duration = 20e-6;    % signal-free tail appended after the last echo

if QUICK
    OS = 4;   N_tx = 17;  N_scatterers = 8000;
else
    OS = 10;  N_tx = 31;  N_scatterers = 16000;
end
fs_sim     = fs*OS;
tx_focus_x = linspace(-tx_x_max, tx_x_max, N_tx);

fprintf('Centre frequency : %.3f MHz\n', f0/1e6);
fprintf('Simulation rate  : %.2f MHz (stored at %.3f MHz)\n', fs_sim/1e6, fs/1e6);
fprintf('Wavelength       : %.3f mm\n', lambda*1e3);
fprintf('Transmit events  : %d\n', N_tx);

%% ---------------- PHANTOM ----------------
rng(10,'twister');

x_min = -16e-3;  x_max = 16e-3;
z_min =  28e-3;  z_max = 90e-3;

x = x_min + (x_max-x_min)*rand(N_scatterers,1);
y = zeros(N_scatterers,1);
z = z_min + (z_max-z_min)*rand(N_scatterers,1);
amplitudes = randn(N_scatterers,1);

cyst_x = 0;  cyst_z = 55e-3;  cyst_radius = 5e-3;
inside_cyst = (x-cyst_x).^2 + (z-cyst_z).^2 <= cyst_radius^2;
amplitudes(inside_cyst) = 0;

% Three point targets for position and resolution checks: [x z]
point_targets = [0 38e-3; -7e-3 72e-3; 7e-3 72e-3];
point_amp     = 15;

phantom_positions = [x y z; point_targets(:,1), zeros(3,1), point_targets(:,2)];
amplitudes        = [amplitudes; point_amp*ones(3,1)];

area_mm2 = (x_max-x_min)*(z_max-z_min)*1e6;
fprintf('\nScatterers: %d (%.1f per mm^2), inside cyst removed: %d\n', ...
    N_scatterers, N_scatterers/area_mm2, sum(inside_cyst));

assert(all(isfinite(phantom_positions(:))) && all(isfinite(amplitudes)));

%% ---------------- FIELD II INITIALISATION ----------------
try, field_end; catch, end
field_init(0);
set_field('c', c);
set_sampling(fs_sim);

% Attenuation 0.5 dB/(MHz cm). Field II units: dB/m at att_f0, slope in dB/(m Hz).
attenuation_dB_per_MHz_cm = 0.5;
freq_att = attenuation_dB_per_MHz_cm*100/1e6;
att      = freq_att*f0;
set_field('att', att);
set_field('freq_att', freq_att);
set_field('att_f0', f0);
set_field('use_att', 1);

%% ---------------- APERTURES, PULSES ----------------
emit_aperture = xdc_linear_array(N_elements, element_width, element_height, ...
    kerf, 1, 1, [0 0 0]);
receive_aperture = xdc_linear_array(N_elements, element_width, element_height, ...
    kerf, 1, 1, [0 0 0]);

t_ir = 0:1/fs_sim:2/f0;
L_ir = numel(t_ir);
win  = 0.5*(1 - cos(2*pi*(0:L_ir-1)/(L_ir-1)));
impulse_response = sin(2*pi*f0*t_ir).*win;
excitation       = impulse_response;

xdc_impulse(emit_aperture, impulse_response);
xdc_impulse(receive_aperture, impulse_response);
xdc_excitation(emit_aperture, excitation);
xdc_apodization(emit_aperture, 0, ones(1,N_elements));
xdc_apodization(receive_aperture, 0, ones(1,N_elements));

%% ---------------- TRANSMIT GEOMETRY AND DELAYS ----------------
element_x = ((1:N_elements) - (N_elements+1)/2)*pitch;     % 1 x N

tx_rmax   = zeros(N_tx,1);
tx_delays = zeros(N_tx,N_elements);
for tx = 1:N_tx
    r = sqrt((element_x - tx_focus_x(tx)).^2 + z_focus^2);
    tx_rmax(tx)      = max(r);
    tx_delays(tx,:)  = (tx_rmax(tx) - r)/c;     % smallest delay is exactly 0
end

%% ---------------- PULSE LAG CALIBRATION ----------------
% Point scatterer at the focus of an on-axis transmit. The envelope peak on
% the two central channels gives the fixed lag of the pulse-echo response.
fprintf('\nCalibrating pulse lag...\n');
r_cal     = sqrt(element_x.^2 + z_focus^2);
r_cal_max = max(r_cal);
xdc_focus_times(emit_aperture, 0, (r_cal_max - r_cal)/c);

[v_cal, t_cal] = calc_scat_multi(emit_aperture, receive_aperture, [0 0 z_focus], 1);
v_cal = to_samples_by_channels(double(v_cal), N_elements);

ch        = [N_elements/2, N_elements/2+1];
env_cal   = abs(hilbert(v_cal(:,ch)));
[~, pk]   = max(env_cal);                         % 1 x 2
t_peak    = t_cal + (pk-1)/fs_sim;
t_expect  = (r_cal_max + r_cal(ch))/c;
pulse_lag = mean(t_peak - t_expect);

fprintf('Calibrated pulse lag: %.4f us\n', pulse_lag*1e6);
if pulse_lag < 0.3e-6 || pulse_lag > 2.5e-6
    warning('Pulse lag is outside the expected 0.3 to 2.5 us range. Inspect v_cal.');
end

%% ---------------- SIMULATE ALL TRANSMIT EVENTS ----------------
raw    = cell(N_tx,1);
t0_all = zeros(N_tx,1);
t_loop = tic;

for tx = 1:N_tx
    xdc_focus_times(emit_aperture, 0, tx_delays(tx,:));

    [v, t0] = calc_scat_multi(emit_aperture, receive_aperture, ...
        phantom_positions, amplitudes);

    raw{tx}    = single(to_samples_by_channels(double(v), N_elements));
    t0_all(tx) = t0;

    el = toc(t_loop);
    fprintf('Transmit %2d/%d done | elapsed %.1f min | ETA %.1f min\n', ...
        tx, N_tx, el/60, el/tx*(N_tx-tx)/60);
end

xdc_free(emit_aperture);
xdc_free(receive_aperture);
field_end;

%% ---------------- COMMON TIME AXIS, DECIMATION ----------------
t_min    = min(t0_all);
pad_sim  = round((t0_all - t_min)*fs_sim);          % front padding at the simulation rate
len_sim  = cellfun(@(v) size(v,1), raw);
L_sim    = ceil(max(pad_sim + len_sim)/OS)*OS;
N_sig    = L_sim/OS;                                 % samples that can contain echoes
tail_samples = ceil(tail_duration*fs);
N_samples    = N_sig + tail_samples;

rf_clean = zeros(N_samples, N_elements, N_tx, 'single');
for tx = 1:N_tx
    buf = zeros(L_sim, N_elements);
    buf(pad_sim(tx) + (1:len_sim(tx)), :) = double(raw{tx});
    dec = resample(buf, 1, OS);                      % anti-alias filter + decimation
    rf_clean(1:N_sig,:,tx) = single(dec(1:N_sig,:));
end
clear raw buf dec

time_axis = t_min + (0:N_samples-1)'/fs;

% Scale to unit peak (Field II output is ~1e-26, which is hostile to float maths)
rf_scale = 1/max(abs(rf_clean(:)));
rf_clean = rf_clean*rf_scale;

fprintf('\nRF dimensions: %d x %d x %d (samples x channels x transmits)\n', ...
    size(rf_clean,1), size(rf_clean,2), size(rf_clean,3));

%% ---------------- NOISE ----------------
% Signal power over the echo part of the record only (excludes the tail)
tmp = double(rf_clean(1:N_sig,:,:));
signal_power = mean(tmp(:).^2);
clear tmp

% Use the same echo-record power definition as DAS_adaptive_CF.m.
noise_power = signal_power/10^(SNR_dB/10);

rng(20,'twister');
noise_unit = randn(size(rf_clean),'single');         % unit-variance noise realisation
noise      = sqrt(noise_power)*noise_unit;
rf_noisy   = rf_clean + noise;

fprintf('Target channel SNR  : %.2f dB\n', SNR_dB);
fprintf('Measured channel SNR: %.2f dB\n', ...
    10*log10(signal_power/mean(double(noise(:)).^2)));
fprintf('Noise model: zero-mean white Gaussian channel noise, common variance across channels/transmits.\n');

%% ---------------- SAVE ----------------
save('fieldII_cyst_RF.mat', ...
    'rf_clean','rf_noisy','noise','noise_unit','signal_power','noise_power', ...
    'fs','fs_sim','OS','f0','c','lambda', ...
    'time_axis','t_min','t0_all','N_sig','tail_samples','rf_scale', ...
    'pulse_lag','tx_delays','tx_rmax','element_x', ...
    'phantom_positions','amplitudes','point_targets', ...
    'tx_focus_x','pitch','z_focus', ...
    'element_width','element_height', ...
    'N_elements','N_tx','N_scatterers','SNR_dB', ...
    'cyst_x','cyst_z','cyst_radius','inside_cyst', ...
    '-v7.3');
fprintf('\nSaved fieldII_cyst_RF.mat\n');

%% ---------------- QUICK LOOK ----------------
figure('Color','w','Name','Phantom');
keep = ~inside_cyst;
plot(x(keep)*1e3, z(keep)*1e3, '.', 'MarkerSize', 2, 'Color', [0.45 0.45 0.45]); hold on;
plot(point_targets(:,1)*1e3, point_targets(:,2)*1e3, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
th = linspace(0, 2*pi, 200);
plot((cyst_x + cyst_radius*cos(th))*1e3, (cyst_z + cyst_radius*sin(th))*1e3, 'r', 'LineWidth', 1.5);
set(gca,'YDir','reverse'); axis equal; grid on;
xlabel('Lateral position (mm)'); ylabel('Depth (mm)');
title('Phantom: scatterers, cyst boundary, point targets');

mid = round(N_tx/2);
env_ch = abs(hilbert(double(rf_clean(:,:,mid))));
figure('Color','w','Name','Channel data');
imagesc(1:N_elements, time_axis*c/2*1e3, 20*log10(env_ch/max(env_ch(:)) + 1e-6));
colormap gray; clim([-60 0]); colorbar;
xlabel('Receive element'); ylabel('Approximate depth (mm)');
title(sprintf('Clean channel envelope, transmit %d', mid));

figure('Color','w','Name','Clean vs noisy channel');
rx_show = N_elements/2;
s_clean = double(rf_clean(:,rx_show,mid));
s_noisy = double(rf_noisy(:,rx_show,mid));
plot(time_axis*1e6, s_noisy, 'Color', [0.7 0.7 0.7]); hold on;
plot(time_axis*1e6, s_clean, 'k');
xlabel('Time (\mus)'); ylabel('Amplitude (unit peak)');
legend('Noisy','Clean'); grid on;
title(sprintf('Channel %d, transmit %d, SNR %.0f dB', rx_show, mid, SNR_dB));

%% ---------------- LOCAL FUNCTIONS ----------------
function v = to_samples_by_channels(v, N)
    % Field II returns samples x channels; transpose if needed.
    if size(v,2) ~= N && size(v,1) == N
        v = v.';
    end
    if size(v,2) ~= N
        error('Unexpected RF dimensions: %d x %d', size(v,1), size(v,2));
    end
end
