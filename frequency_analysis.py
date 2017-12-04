from scipy.fftpack import fft, fftfreq, fftshift
import matplotlib.pyplot as plt
import numpy as np
def linestyles():
    for i in ['dotted', 'solid', 'dashed']:
        yield i
linestyle = linestyles()
def plot_frequency(DATA, timelist, label='', QOBJ=True):
    if QOBJ:
        N = len(DATA.expect[1])
        y = DATA.expect[1]
        yf = fft(y)
        T = abs(timelist[1]-timelist[0])

        xf = fftfreq(N, T)
        xf = fftshift(xf)
        yplot = fftshift(yf)
    else:
        N = len(DATA)
        y = DATA
        yf = fft(y)
        T = abs(timelist[1]-timelist[0])

        xf = fftfreq(N, T)
        xf = -fftshift(xf)
        yplot = fftshift(yf)
        print "Sample spacing {}, sample rate {}".format(T, 1/(N*T))
    # sample spacing

    plt.plot(2*np.pi*xf, 1.0/N * np.abs(yplot), label=label, linestyle=linestyle.next())
    plt.grid()
    plt.ylabel("Coherence frequency spectrum")
    plt.xlabel("frequency")
