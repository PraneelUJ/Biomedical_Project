
function [rf_out, details] = adaptive_frequency_bank( ...
    rf, fs, f_low, f_high, noise_idx, win_length, min_gain)

% ADAPTIVE FREQUENCY BANK ULTRASOUND RF DENOISING
%
% Inputs:
% rf         : Samples x receive channels
% fs         : RF sampling frequency (Hz)
% f_low      : Lower useful frequency (Hz)
% f_high     : Upper useful frequency (Hz)
% noise_idx  : Indices of noise-only samples
% win_length : Local power estimation window (odd integer)
% min_gain   : Minimum Wiener gain
%
% Outputs:
% rf_out     : Denoised RF data
% details    : Intermediate processing information

    rf = double(rf);

    [Ns, Nc] = size(rf);

    if f_low <= 0 || f_high >= fs/2 || f_low >= f_high
        error('Invalid frequency limits.');
    end

    if isempty(noise_idx) || any(noise_idx < 1) || ...
            any(noise_idx > Ns)
        error('Provide valid noise-only sample indices.');
    end

    if mod(win_length,2) == 0
        error('win_length must be odd.');
    end

    % Remove channel-wise DC offset
    rf = rf - mean(rf,1);

    % Frequency vector
    NFFT = 2^nextpow2(Ns);

    f = (0:NFFT-1)' * fs/NFFT;

    % Fold frequencies into positive absolute frequencies
    f_abs = min(f, fs-f);

    % Four band centers
    centers = linspace(f_low, f_high, 4);

    % Gaussian bandwidth
    sigma_f = (f_high-f_low)/4;

    % Global passband taper
    passband = zeros(NFFT,1);

    idx = f_abs >= f_low & f_abs <= f_high;

    passband(idx) = 1;

    % Smooth transition near passband edges
    edge_width = 0.05*(f_high-f_low);

    left = f_abs >= f_low-edge_width & ...
           f_abs < f_low;

    right = f_abs > f_high & ...
            f_abs <= f_high+edge_width;

    passband(left) = 0.5 - 0.5*cos( ...
        pi*(f_abs(left)-f_low+edge_width)/edge_width);

    passband(right) = 0.5 + 0.5*cos( ...
        pi*(f_abs(right)-f_high)/edge_width);

    % Construct normalized overlapping band windows
    W = zeros(NFFT,4);

    for k = 1:4
        W(:,k) = exp(-0.5*((f_abs-centers(k))/sigma_f).^2);
    end

    W = W ./ (sum(W,2)+eps);

    W = W .* passband;

    % FFT of input
    X = fft(rf,NFFT,1);

    rf_out = zeros(Ns,Nc);

    gains = zeros(Ns,Nc,4);
    noise_power = zeros(Nc,4);

    for k = 1:4

        % Frequency decomposition
        Xk = X .* W(:,k);

        xk_full = ifft(Xk,[],1,'symmetric');

        xk = xk_full(1:Ns,:);

        % Noise estimation per channel
        pn = mean(xk(noise_idx,:).^2,1);

        noise_power(:,k) = pn(:);

        % Local power estimation
        local_power = movmean(xk.^2, ...
            win_length,1);

        % Estimate local signal power
        signal_power = max(local_power-pn,0);

        % Adaptive Wiener gain
        gain = signal_power ./ ...
            (local_power+eps);

        % Prevent excessive suppression
        gain = max(gain,min_gain);

        % Smooth gain to reduce abrupt changes
        gain = movmean(gain,5,1);

        % Apply adaptive gain
        xk_filtered = gain .* xk;

        % Reconstruct
        rf_out = rf_out + xk_filtered;

        gains(:,:,k) = gain;

    end

    details.gains = gains;
    details.noise_power = noise_power;
    details.centers = centers;
    details.frequency = f;
    details.windows = W;

end
