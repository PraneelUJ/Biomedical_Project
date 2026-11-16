
clear;
clc;
close all;

%% PARAMETERS

f0 = 2.717391e6;

fs = 10.85914e6;
OS = 10;
fs_sim = fs*OS;

c = 1540;
lambda = c/f0;

N_elements = 64;
N_tx = 32;
N_scatterers = 10000;

element_width = 0.25e-3;
element_height = 3e-3;
pitch = 0.30e-3;                 % Assumed pitch
kerf = pitch-element_width;

z_focus = 60e-3;
tx_focus_x = linspace(-8e-3,8e-3,N_tx);

SNR_dB = 0;

%% PHANTOM

rng(10);

x_min = -15e-3;
x_max =  15e-3;

z_min = 25e-3;
z_max = 90e-3;

x = x_min+(x_max-x_min)*rand(N_scatterers,1);
y = zeros(N_scatterers,1);
z = z_min+(z_max-z_min)*rand(N_scatterers,1);

cyst_x = 0;
cyst_z = 55e-3;
cyst_radius = 5e-3;

inside_cyst = ...
    (x-cyst_x).^2+(z-cyst_z).^2 <= cyst_radius^2;

amplitudes = randn(N_scatterers,1);
amplitudes(inside_cyst) = 0;

phantom_positions = [x y z];

%% FIELD II INITIALIZATION

field_init(0);

set_sampling(fs_sim);
set_field('c',c);

% Attenuation: 0.5 dB/(MHz cm)
attenuation_dB_per_MHz_cm = 0.5;

freq_att = attenuation_dB_per_MHz_cm*100/1e6;
att = freq_att*f0;

set_field('att',att);
set_field('freq_att',freq_att);
set_field('att_f0',f0);
set_field('use_att',1);

%% TRANSDUCER APERTURES

focus = [0 0 0];

emit_aperture = xdc_linear_array( ...
    N_elements,element_width,element_height,kerf, ...
    1,1,focus);

receive_aperture = xdc_linear_array( ...
    N_elements,element_width,element_height,kerf, ...
    1,1,focus);

%% EXCITATION AND IMPULSE RESPONSE

t_ir = (0:1/fs_sim:2/f0).';

pulse_window = hann(numel(t_ir));

impulse_response = ...
    sin(2*pi*f0*t_ir).*pulse_window;

excitation = ...
    sin(2*pi*f0*t_ir).*pulse_window;

% Set all aperture pulses after setting the sampling frequency.
xdc_impulse(emit_aperture,impulse_response.');
xdc_impulse(receive_aperture,impulse_response.');
xdc_excitation(emit_aperture,excitation.');

%% CALIBRATE PULSE DELAY USING A POINT SCATTERER

fprintf('Calibrating pulse timing...\n');

calibration_point = [0 0 z_focus];

xdc_focus(emit_aperture,0,calibration_point);

[v_cal,t_cal] = calc_scat_multi( ...
    emit_aperture,receive_aperture, ...
    calibration_point,1);

v_cal = double(v_cal);

% Ensure samples are in rows and receive channels in columns.
if size(v_cal,2) ~= N_elements && size(v_cal,1)==N_elements
    v_cal = v_cal.';
end

if size(v_cal,2) ~= N_elements
    error('Unexpected calibration RF dimensions.');
end

% Average envelope of the two central receive channels.
central_channels = [N_elements/2,N_elements/2+1];

calibration_envelope = mean( ...
    abs(hilbert(v_cal(:,central_channels))),2);

[~,peak_index] = max(calibration_envelope);

measured_peak_time = t_cal+(peak_index-1)/fs_sim;

% At the focal point:
% transmit propagation = z_focus/c
% receive propagation  = z_focus/c

expected_propagation_time = 2*z_focus/c;

pulse_delay = measured_peak_time-expected_propagation_time;

fprintf('Measured RF peak time : %.6f us\n',measured_peak_time*1e6);
fprintf('Expected propagation  : %.6f us\n', ...
    expected_propagation_time*1e6);
fprintf('Calibrated pulse delay: %.6f us\n',pulse_delay*1e6);

if pulse_delay < 0
    warning('Measured pulse delay is negative. Inspect calibration RF.');
end

%% SIMULATE ALL TRANSMISSIONS

rf_cell = cell(N_tx,1);
time_start = zeros(N_tx,1);

fprintf('\nStarting Field II simulation at %.3f MHz\n', ...
    fs_sim/1e6);

for tx = 1:N_tx

    focus_point = [tx_focus_x(tx),0,z_focus];

    xdc_focus(emit_aperture,0,focus_point);

    [v,t0] = calc_scat_multi( ...
        emit_aperture,receive_aperture, ...
        phantom_positions,amplitudes);

    v = double(v);

    if size(v,2) ~= N_elements && size(v,1)==N_elements
        v = v.';
    end

    if size(v,2) ~= N_elements
        error('Unexpected RF dimensions in transmission %d.',tx);
    end

    % Anti-alias filtering and decimation.
    % resample compensates for its own FIR filter delay.
    v_dec = resample(v,1,OS);

    rf_cell{tx} = single(v_dec);
    time_start(tx) = t0;

    if mod(tx,8)==0 || tx==N_tx
        fprintf('Completed %d/%d transmissions\n',tx,N_tx);
    end

end

%% ALIGN TRANSMISSIONS

t_global_start = min(time_start);

offsets = round((time_start-t_global_start)*fs);

signal_lengths = cellfun(@(v) size(v,1),rf_cell);

tail_duration = 20e-6;
tail_samples = ceil(tail_duration*fs);

N_samples = max(offsets+signal_lengths)+tail_samples;

rf_clean = zeros(N_samples,N_elements,N_tx,'single');

for tx = 1:N_tx

    v = rf_cell{tx};

    first_sample = offsets(tx)+1;
    last_sample = first_sample+size(v,1)-1;

    rf_clean(first_sample:last_sample,:,tx) = v;

end

time_axis = t_global_start+(0:N_samples-1)'/fs;

%% CALCULATE CHANNEL NOISE POWER

% Calculate signal power only over actual RF samples.
% Do not include zero padding or the added noise-only tail.

total_signal_energy = 0;
total_signal_samples = 0;

for tx = 1:N_tx

    v = double(rf_cell{tx});

    total_signal_energy = total_signal_energy+sum(v(:).^2);
    total_signal_samples = total_signal_samples+numel(v);

end

signal_power = total_signal_energy/total_signal_samples;

noise_power = signal_power/10^(SNR_dB/10);

rng(20);

noise = single(sqrt(noise_power)* ...
    randn(size(rf_clean),'single'));

rf_noisy = rf_clean+noise;

measured_channel_SNR = 10*log10( ...
    signal_power/mean(double(noise(:)).^2));

fprintf('\n========== SIMULATION RESULTS ==========\n');

fprintf('Output sampling frequency : %.5f MHz\n',fs/1e6);
fprintf('RF dimensions              : %d x %d x %d\n', ...
    size(rf_clean,1),size(rf_clean,2),size(rf_clean,3));

fprintf('Requested channel SNR      : %.2f dB\n',SNR_dB);
fprintf('Measured channel SNR       : %.2f dB\n', ...
    measured_channel_SNR);

%% SAVE DATA

save('fieldII_cyst_RF.mat', ...
    'rf_clean','rf_noisy','noise', ...
    'fs','fs_sim','OS','f0','c','lambda', ...
    'time_axis','time_start','t_global_start','offsets', ...
    'pulse_delay', ...
    'phantom_positions','amplitudes', ...
    'tx_focus_x','pitch','z_focus', ...
    'element_width','element_height', ...
    'N_elements','N_tx','N_scatterers','SNR_dB', ...
    'cyst_x','cyst_z','cyst_radius','inside_cyst', ...
    '-v7.3');

%% RELEASE FIELD II

xdc_free(emit_aperture);
xdc_free(receive_aperture);

field_end;

fprintf('\nSaved fieldII_cyst_RF.mat\n');
