import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import freqz

def plot_initial_extrema(ctx):
    plt.figure()
    plt.plot(ctx.grid, np.zeros(ctx.gridSize))
    plt.scatter([ctx.grid[i] for i in ctx.ext],np.zeros(len(ctx.ext)), label="Initial extrema")
    plt.xlabel("Normalized frequency")
    plt.ylabel("Amplitude")
    plt.title("Initial guess of extremal frequencies")
    plt.legend()
    plt.grid()

def plot_error(ctx, iteration=None):
    plt.figure()
    plt.plot(ctx.grid, ctx.E, label="Weight error")
    plt.scatter([ctx.grid[i] for i in ctx.ext], [ctx.E[i] for i in ctx.ext])
    plt.xlabel("Normalized frequency")
    plt.ylabel("Weighted error E(f)")
    if iteration is not None:
        plt.title(f"Remez error — iteration {iteration}")
    else:
        plt.title("Remez weighted error")
    plt.legend()
    plt.grid()


def plot_time_signal(signal,filtered,num_samples=65536):
    N = min(num_samples,len(signal))
    plt.figure()
    plt.plot(signal[:N],label="Input")
    plt.plot(filtered[:N],label="Filtered")
    plt.xlabel("Sample")
    plt.ylabel("Amplitude")
    plt.title("Signal before and after FIR")
    plt.legend()
    plt.grid()

def calculate_spectrum(signal,fs,num_samples=65536):
    N = min(num_samples,len(signal))
    signal = signal[:N]

    window = np.hanning(N)
    spectrum = np.fft.rfft(signal * window)
    frequencies = np.fft.rfftfreq(N,1.0 / fs)

    magnitude = np.abs(spectrum)

    magnitude_db = 20 * np.log10(
        np.maximum(
            magnitude,
            1e-12
        )
    )

    return frequencies,magnitude_db

def plot_spectrum(signal,filtered,fs):
    f_input,spectrum_input = calculate_spectrum(signal,fs)
    f_filtered,spectrum_filtered = calculate_spectrum(filtered,fs)
    plt.figure()
    plt.plot(f_input,spectrum_input,label="Input")
    plt.plot(f_filtered,spectrum_filtered,label="Filtered")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude [dB]")
    plt.title("Spectrum before and after FIR")
    plt.legend()
    plt.grid()

def plot_remez_comparison(h_my,h_ref):
    w,H_my = freqz(h_my,worN=4096)
    _,H_ref = freqz(h_ref,worN=4096)
    f = w / np.pi
    plt.figure()
    plt.plot(f,np.abs(H_my),label="Our Remez")
    plt.plot(f,np.abs(H_ref),label="scipy.signal.remez")
    plt.xlabel("Normalized frequency")
    plt.ylabel("|H(f)|")
    plt.title("Comparison of frequency responses")
    plt.legend()
    plt.grid()

def plot_response_error(h_my,h_ref):
    w,H_my = freqz(h_my,worN=4096)
    _,H_ref = freqz(h_ref,worN=4096)
    f = w / np.pi
    error = np.abs(H_my) - np.abs(H_ref)
    plt.figure()
    plt.plot(f,error)
    plt.axhline(0,linewidth=0.8)
    plt.xlabel("Normalized frequency")
    plt.ylabel("H_our - H_reference")
    plt.title("Frequency response difference")
    plt.grid()