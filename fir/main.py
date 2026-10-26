import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import remez as scipy_remez

from remez import remez
from fir import apply_fir
from plots import (
    plot_initial_extrema,
    plot_time_signal,
    plot_spectrum,plot_remez_comparison
)

FS = 62314000
NUM_TAPS = 11

BANDS = [
    0.0,
    0.97,
    0.98,
    1.0
]

DESIRED = [
    1.0,
    1.0,
    0.0,
    0.0
]
SCIPY_DESIRED = [
    1.0,
    0.0
]
WEIGHT = [
    1.0,
    1.0
]

print("Designing filter using our Remez...")

h_my, ctx = remez(
    NUM_TAPS,
    BANDS,
    DESIRED,
    WEIGHT,
    grid_density=16
)

print("\nOur coefficients:")
print(h_my)

plot_initial_extrema(ctx)

print("\nDesigning reference filter using scipy.signal.remez...")

h_ref = scipy_remez(NUM_TAPS,BANDS,SCIPY_DESIRED,weight=WEIGHT,fs=2.0)

print("\nReference coefficients:")
print(h_ref)

coefficient_error = h_my - h_ref

max_coefficient_error = np.max(np.abs(coefficient_error))



print("\nCOEFFICIENT COMPARISON")
print("Maximum coefficient error:",max_coefficient_error)


plot_remez_comparison(h_my,h_ref)


symmetry_error = np.max(np.abs(h_my - h_my[::-1]))

print("\nFILTER SYMMETRY")
print("Maximum symmetry error:",symmetry_error)

# signal = np.fromfile("Glonass_L1_CA_SIM_Fs_62Msps_4ms.dat",dtype=np.complex64)

IF_HZ = 9540000

signal = np.fromfile(
    "Glonass_L1_CA_SIM_Fs_62Msps_4ms.dat",
    dtype=np.complex64
)

n = np.arange(len(signal), dtype=np.float64)

# Перенос частоты IF в окрестность нуля
mixer = np.exp(-1j * 2.0 * np.pi * IF_HZ * n / FS)

baseband = signal * mixer

# Фильтрация уже перенесённого сигнала
filtered = apply_fir(baseband, h_my)

filtered.astype(np.complex64).tofile("signal_filtered.dat")

print("\nSIGNAL")
print("Number of samples:",len(signal))
print("Minimum:",signal.min())
print("Maximum:",signal.max())
print("Mean:",signal.mean())

filtered = apply_fir(signal,h_my)

filtered.astype(np.complex64).tofile("signal_filtered.dat")

print("\nFiltered signal saved to signal_filtered.dat")

plot_time_signal(signal,filtered,num_samples=65536)

plot_spectrum(signal,filtered,FS)

plt.show()