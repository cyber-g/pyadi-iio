# Extend the example: 
# https://github.com/analogdevicesinc/pyadi-iio/blob/master/examples/ad9361_example.py
# to use a bandlimited signal.
#
# Author: Germain PHAM, @cyber-g, Aug. 2023

import time

import adi
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.fftpack import fftshift

def compute_fir_kaiser(fpass,fstop,fs,atten_level=60):
    """
    Compute the coefficients of a Kaiser FIR filter

    Parameters
    ----------
    fpass : float
        Passband frequency [Hz] 
    fstop : float
        Stopband frequency [Hz]
    fs : float
        Sampling frequency [Hz]
    atten_level : float, optional
        Attenuation level [dB]. 
        The default is 60.

    Returns
    -------
    taps : numpy array
        Filter coefficients

    Author: Germain Pham, @cyber-g, Aug. 2023
    Inspired by: 
    https://github.com/unpingco/Python-for-Signal-Processing/blob/master/Filtering_Part2.ipynb
    """
    
    # Compute the filter order
    M,beta= signal.kaiserord(atten_level, (fstop-fpass)/(fs/2.))
    
    # Design the filter
    taps = signal.firwin(M,(fstop+fpass)/2.,window=('kaiser',beta),fs=fs)
    
    return taps

def generate_filtered_noise(size_array, BW, FS):
    """
    Generate a filtered noise signal

    Parameters
    ----------
    size_array : tuple of int
        matrix size
    BW : float
        RF Bandwidth of the signal (bilateral)  [Hz]
    FS : float
        Sampling frequency                      [Hz]

    Returns
    -------
    noise_sig : numpy array
        Filtered noise signal

    Author: Germain Pham, @cyber-g, Aug. 2023
    """
    # Generate the noise
    noise = np.random.normal(size=size_array)
    
    # Instanciate the filter
    taps = compute_fir_kaiser(BW/2, (BW/2)*1.1, FS)

    # Filter the noise, use convolve in order to get transients fade in/out
    # making the signal more circular (useful for spectral analysis and cyclic
    # transmission)
    noise_sig = signal.convolve(taps, noise)

    # Normalize the signal to have unitary magnitude
    noise_sig = noise_sig/np.max(np.abs(noise_sig))
    
    return noise_sig

def print_all_properties(sdr):
    """
    Print all properties of the radio

    Parameters
    ----------
    sdr : object
        Radio object

    Returns
    -------
    None.

    Properties are found in pyadi-iio/adi/ad936x.py

    Author: Germain Pham, @cyber-g, Aug. 2023
    """
    
    print('The SDR properties are:')    
    print('\t loopback: ______________ ',(sdr.loopback))
    print('\t gain_control_mode_chan0: ',(sdr.gain_control_mode_chan0))
    print('\t rx_hardwaregain_chan0: _ ',(sdr.rx_hardwaregain_chan0))
    print('\t tx_hardwaregain_chan0: _ ',(sdr.tx_hardwaregain_chan0))
    print('\t rx_rf_bandwidth: _______ ',(sdr.rx_rf_bandwidth))
    print('\t tx_rf_bandwidth: _______ ',(sdr.tx_rf_bandwidth))
    print('\t sample_rate: ___________ ',(sdr.sample_rate))
    print('\t rx_lo: _________________ ',(sdr.rx_lo))
    print('\t tx_lo: _________________ ',(sdr.tx_lo))
    print('\t gain_control_mode_chan1: ',(sdr.gain_control_mode_chan1))
    print('\t rx_hardwaregain_chan1: _ ',(sdr.rx_hardwaregain_chan1))
    print('\t tx_hardwaregain_chan1: _ ',(sdr.tx_hardwaregain_chan1))

# Create radio
sdr = adi.ad9361(uri="ip:analog.local")

# Configure properties
sdr.filter                  = 'phaser/LTE20_MHz.ftr'

# Print all properties
print_all_properties(sdr)

# Configure some parameters of the radio
sdr.rx_lo                   = int(2e9)
sdr.tx_lo                   = int(2e9)
sdr.tx_cyclic_buffer        = True
sdr.tx_hardwaregain_chan0   = -30
sdr.gain_control_mode_chan0 = "slow_attack"

# Configuration data channels
sdr.rx_enabled_channels = [0]
sdr.tx_enabled_channels = [0]

# Create a complex-valued bandlimited signal
fs    = int(sdr.sample_rate)
N     = 1024
bw    = 5e6 # anything below 19.36 MHz (tx_rf_bandwidth & rx_rf_bandwidth)
ts    = 1 / float(fs)
t     = np.arange(0, N * ts, ts)

DAC_full_scale_int = 2 ** 14

i = generate_filtered_noise(N, bw, fs) * DAC_full_scale_int
q = generate_filtered_noise(N, bw, fs) * DAC_full_scale_int
iq = i + 1j * q

# Send data
sdr.tx(iq)

# Wait for 1 second to allow for TX to settle
time.sleep(1)

# Collect data
for r in range(20):
    x = sdr.rx()
    f, Pxx_den = signal.welch(x, fs, return_onesided=False)
    plt.clf()
    plt.semilogy(fftshift(f), fftshift(Pxx_den))
    plt.ylim([1e-9, 1e3])
    plt.xlabel("frequency [Hz]")
    plt.ylabel("PSD [V**2/Hz]")
    plt.grid()
    plt.draw()
    plt.pause(0.05)
    time.sleep(0.1)
plt.title("Received signal")

# Plot the signal sent to DACs
plt.figure()
f, Pxx_den = signal.welch(iq, fs, return_onesided=False)
plt.semilogy(fftshift(f), fftshift(Pxx_den))
plt.ylim([1e-9, 1e3])
plt.xlabel("frequency [Hz]")
plt.ylabel("PSD [V**2/Hz]")
plt.grid()
plt.title("Signal sent to DACs")

plt.show()
