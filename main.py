import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import remez as scipy_remez

from remez import remez
from fir import apply_fir
from plots import (
    plot_initial_extrema,
    plot_error,
    plot_response_error,
    plot_time_signal,
    plot_spectrum,plot_remez_comparison
)

FS = 4000000
NUM_TAPS = 5

BANDS = [
    0.0,
    0.45,
    0.55,
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
plot_error(ctx)

print("\nDesigning reference filter using scipy.signal.remez...")

h_ref = scipy_remez(NUM_TAPS,BANDS,SCIPY_DESIRED,weight=WEIGHT,fs=2.0)

print("\nReference coefficients:")
print(h_ref)

coefficient_error = h_my - h_ref

max_coefficient_error = np.max(np.abs(coefficient_error))



print("\nCOEFFICIENT COMPARISON")
print("Maximum coefficient error:",max_coefficient_error)


plot_remez_comparison(h_my,h_ref)
plot_response_error(h_my,h_ref)



symmetry_error = np.max(np.abs(h_my - h_my[::-1]))

print("\nFILTER SYMMETRY")
print("Maximum symmetry error:",symmetry_error)

signal = np.fromfile("signal_source.dat",dtype=np.complex64)


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