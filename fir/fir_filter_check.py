import numpy as np
import matplotlib.pyplot as plt

fs = 4000000   

x = np.fromfile("signal_source.dat", dtype=np.float32)
iq = x[::2] + 1j*x[1::2]

N = 65536

spectrum = np.fft.fftshift(
    np.fft.fft(iq[:N])
)

freq = np.linspace(-fs/2, fs/2, N)
plt.subplot(2,1,1)
plt.plot(freq/1e6, 20*np.log10(abs(spectrum)))
plt.xlabel("MHz")
plt.ylabel("dB")
plt.title("signal_source.dat")
plt.grid()

x2 = np.fromfile("signal_filtered.dat", dtype=np.float32)
iq2 = x2[::2] + 1j*x2[1::2]

N2 = 65536

spectrum2 = np.fft.fftshift(
    np.fft.fft(iq2[:N2])
)

freq2 = np.linspace(-fs/2, fs/2, N2)

plt.subplot(2,1,2)
plt.plot(freq2/1e6, 20*np.log10(abs(spectrum2)))
plt.xlabel("MHz")
plt.ylabel("dB")
plt.title("our_fir.dat")
plt.grid()
plt.show()