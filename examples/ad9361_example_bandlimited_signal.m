% Convert the ad9361_example_bandlimited_signal.py to Matlab
% which extend the example: 
% https://github.com/analogdevicesinc/pyadi-iio/blob/master/examples/ad9361_example.py
% to use a bandlimited signal.
%
% Author: Germain PHAM, @cyber-g, Aug. 2023


% Create radio
tx = sdrtx('AD936x','IPAddress','analog.local')
rx = sdrrx('AD936x','IPAddress','analog.local')

% Configure properties
BasebandSampleRate          = 61.44e6;
Fc                          = 2e9;
enabled_channels            = 1;
tx.BasebandSampleRate       = BasebandSampleRate
rx.BasebandSampleRate       = BasebandSampleRate
tx.CenterFrequency          = Fc
rx.CenterFrequency          = Fc
tx.ChannelMapping           = enabled_channels
rx.ChannelMapping           = enabled_channels
tx.Gain                     = -30
rx.GainSource               = "AGC Slow Attack"
rx.OutputDataType           = "double" 
% 'double' — Double-precision floating point values are scaled to the range of
% [-1, 1]. The object derives this value from the sign-extended 16 bits received
% from the board.

% Create a complex-valued bandlimited signal
fs    = BasebandSampleRate;
N     = 1024;
bw    = 5e6;
ts    = 1/fs;
t     = 0:ts:(N-1)*ts;

DAC_full_scale_int = 2^14;

i = generate_filtered_noise([N,1], bw, fs) * DAC_full_scale_int;
q = generate_filtered_noise([N,1], bw, fs) * DAC_full_scale_int;
iq = i + 1j * q;

% Send data
% Double-precision floating point — Complex values in the range of [-1, 1].
% Since the AD9361/AD9364 RF chip has a 12-bit DAC, numbers of magnitude less
% than 0.0625 are lost.
transmitRepeat(iq)

% Wait for 1 second to allow for TX to settle
pause(1)

% Collect data
rx.SamplesPerFrame = N;
fshift = (-N/2:N/2-1)*(fs/N); % zero-centered frequency range
for r = 1:20
    x = rx()
    Pxx_den = abs(fft(x)).^2/N;
    clf()
    semilogy(fftshift(f), fftshift(Pxx_den))
    ylim([1e-7, 1e2])
    xlabel("frequency [Hz]")
    ylabel("PSD [''V**2/Hz'']")
    grid on
    drawnow()
    pause(0.1)
title("Received signal")


% Plot the signal sent to DACs
figure()
Pxx_den = abs(fft(iq)).^2/N;
semilogy(fftshift(f), fftshift(Pxx_den))
ylim([1e-7, 1e3])
xlabel("frequency [Hz]")
ylabel("PSD [''V**2/Hz'']")
grid on
title("Signal sent to DACs")



function taps = compute_fir_kaiser(fpass,fstop,fs,atten_level)
%compute_fir_kaiser - Compute the coefficients of a Kaiser FIR filter
% 
% taps = compute_fir_kaiser(fpass,fstop,fs,atten_level)
% 
% Parameters
% ----------
% fpass : float
%     Passband frequency [Hz] 
% fstop : float
%     Stopband frequency [Hz]
% fs : float
%     Sampling frequency [Hz]
% atten_level : float, optional
%     Attenuation level [dB]. 
%     The default is 60.
% 
% Returns
% -------
% taps : numpy array
%     Filter coefficients
% 
% Author: Germain Pham, @cyber-g, Aug. 2023

    % Check if the attenuation level is provided, otherwise set it to 60 dB
    if nargin < 4
        atten_level = 60; % dB
    end

    % Compute the filter order
    % [N,Wn,BTA,FILTYPE] = KAISERORD(F,A,DEV,Fs)
    [N,Wn,beta,FILTYPE] = kaiserord([fpass fstop],[1 0],[1 10^(-atten_level/20)],fs);

    % Design the filter
    taps = fir1(N, Wn, FILTYPE, kaiser(N+1,beta), 'noscale');

end

function noise_sig = generate_filtered_noise(size_array, BW, FS)
%generate_filtered_noise - Generate a filtered noise signal
% 
% noise_sig = generate_filtered_noise(size_array, BW, FS)
% 
% Parameters
% ----------
% size_array : row vector of int
%     matrix size
% BW : float
%     Bandwidth of the signal [Hz]
% FS : float
%     Sampling frequency      [Hz]
% 
% Returns
% -------
% noise_sig : numpy array
%     Filtered noise signal
% 
% Author: Germain Pham, @cyber-g, Aug. 2023

    % Generate the noise
    noise = randn(size_array);

    % Instanciate the filter
    taps  = compute_fir_kaiser(0, BW/2, FS);

    % Filter the noise, use convolve in order to get transients fade in/out
    % making the signal more circular (useful for spectral analysis and cyclic
    % transmission)
    noise_sig = conv(taps, noise);

    % Normalize the signal to have unitary magnitude
    noise_sig = noise_sig/max(abs(noise_sig));

end