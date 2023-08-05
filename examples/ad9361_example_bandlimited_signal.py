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
    M,beta= signal.fir_filter_design.kaiserord(atten_level, (fstop-fpass)/(fs/2.))
    
    # Design the filter
    taps = signal.firwin(M,(fstop+fpass)/2.,window=('kaiser',beta),nyq=fs/2.)
    
    return taps

def generate_filtered_noise(size_array, BW, FS):
    """
    Generate a filtered noise signal

    Parameters
    ----------
    size_array : tuple of int
        matrix size
    BW : float
        Bandwidth of the signal [Hz]
    FS : float
        Sampling frequency      [Hz]

    Returns
    -------
    noise_sig : numpy array
        Filtered noise signal

    Author: Germain Pham, @cyber-g, Aug. 2023
    """
    # Generate the noise
    noise = np.random.normal(size=size_array)
    
    # Instanciate the filter
    taps = compute_fir_kaiser(0, BW, FS)

    # Filter the noise, use convolve in order to get transients fade in/out
    # making the signal more circular (useful for spectral analysis and cyclic
    # transmission)
    noise_sig = signal.convolve(taps, noise)

    # Normalize the signal to have unitary magnitude
    noise_sig = noise_sig/np.max(np.abs(noise_sig))
    
    return noise_sig



# Create radio
sdr = adi.ad9361(uri="ip:analog.local")

# Configure properties
sdr.rx_rf_bandwidth         = int(4e6)
sdr.sample_rate             = int(6e6)
sdr.rx_lo                   = int(2e9)
sdr.tx_lo                   = int(2e9)
sdr.tx_cyclic_buffer        = True
sdr.tx_hardwaregain_chan0   = -30
sdr.gain_control_mode_chan0 = "slow_attack"

# Configuration data channels
sdr.rx_enabled_channels = [0]
sdr.tx_enabled_channels = [0]

# Read properties
print("RX LO %s" % (sdr.rx_lo))

# Create a sinewave waveform
fs    = int(sdr.sample_rate)
N     = 1024
fc    = int(1e6 / (fs / N)) * (fs / N)
ts    = 1 / float(fs)
t     = np.arange(0, N * ts, ts)

DAC_full_scale_int = 2 ** 14

i = np.cos(2 * np.pi * t * fc) * DAC_full_scale_int
q = np.sin(2 * np.pi * t * fc) * DAC_full_scale_int
iq = i + 1j * q

# Send data
sdr.tx(iq)

# Collect data
for r in range(20):
    x = sdr.rx()
    f, Pxx_den = signal.periodogram(x, fs)
    plt.clf()
    plt.semilogy(f, Pxx_den)
    plt.ylim([1e-7, 1e2])
    plt.xlabel("frequency [Hz]")
    plt.ylabel("PSD [V**2/Hz]")
    plt.draw()
    plt.pause(0.05)
    time.sleep(0.1)

plt.show()
