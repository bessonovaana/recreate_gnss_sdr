import numpy as np
import matplotlib.pyplot as plt


f1 = 20
f2 =10
t = np.linspace(0,1, 1000)
phi1 = 2*np.pi
phi2 = np.pi/2

X = np.sin(2*np.pi*t*f1+phi1)
Y = np.sin(2*np.pi*t*f2+phi2)

shift = X * Y

plt.subplot(2,2,1); plt.plot(t, X);    plt.title('X, 20 Гц')
plt.subplot(2,2,2); plt.plot(t, Y);    plt.title('Y, 10 Гц')
plt.subplot(2,1,2); plt.plot(t, shift); plt.title('XY — пики на 10 и 30 Гц')
plt.tight_layout(); plt.show()