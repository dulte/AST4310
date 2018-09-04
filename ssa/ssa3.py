import numpy as np
import matplotlib.pyplot as plt
from scipy import constants


def planck(temp,wavelength):
    c = 2.99792e10
    h = 6.62607e-27 #Planck constant (erg s)
    k = 1.38e-16
    return (2*h*c**2/wavelength**5)*1./(np.exp(h*c/(wavelength*k*temp))-1)

def plot_planck():
    wav = np.arange(1000,20801,200)
    b = np.zeros(wav.shape)
    plt.xlabel(r'wavelength $\lambda / \AA$', size=14)
    plt.ylabel(r'Planck function', size=14)
    plt.xlim(0,20800)
    for T in range(8000,5000-1,-200):
        b[:] = planck(T, wav[:]*1e-8)
        plt.plot(wav,b,'-')
    plt.show()
if __name__ == '__main__':
    plot_planck()
