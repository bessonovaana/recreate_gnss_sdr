import numpy as np


def apply_fir(signal, h):
    signal = np.asarray(signal, dtype=np.complex64)
    h = np.asarray(h, dtype=np.float32)
    filtered = np.zeros_like(signal)
    for n in range(len(signal)):
        s = 0.0 + 0.0j
        for k in range(len(h)):
            if n - k >= 0:
                s += signal[n - k] * h[k]
        filtered[n] = s
    return filtered