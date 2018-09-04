import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
try:
    import seaborn as sns
    sns.set_style("darkgrid")
except ImportError as e:
    print("To get the full joy out of this program install the following modules")
    print(e)

rc('font',**{'family':'serif'})

"""
Schadeenium:
"""
def partition_function_E(temp):
    chiion = np.array([7, 16, 31, 51]) # Schadee ionization energies into numpy array
    k = 8.61734e-5 # Boltzmann constant in eV/deg
    u = np.zeros(4) # declare a 4 zero-element array
    for r in range(4):
        for s in range(chiion[r]):
            u[r] = u[r] + np.exp(-s/k/temp)
    return u # returns all the values of u array


def boltz_E(temp, r, s):
    u = partition_function_E(temp)
    KeV = 8.61734e-5 # This constant does need to be defined here again if it was before
    relnrs = 1. / u[int(r) - 1] * np.exp(-(int(s) - 1) / (KeV * temp))
    return relnrs


def saha_E(temp, elpress, ionstage):
    kerg = 1.380658e-16
    kev = 8.61734e-5 #Boltzmann constant (eV/deg)
    h = 6.62607e-27 #Planck constant (erg s)
    elmass = 9.109390e-28 #electron mass (g)
    keVT = kev * temp
    kergT = kerg * temp
    eldens = elpress / kergT
    chiion = np.array([7, 16, 31, 51 ])
    u = partition_function_E(temp)
    u = np.append(u, 2) # With this command we are adding a new element to the array
    sahaconst = (2. * np.pi * elmass * kergT / (h**2))**1.5 * 2. / eldens
    nstage = np.zeros(5)
    nstage[0] = 1.
    # We set the first element of the array to a value 1
    for r in range(4):
        nstage[r + 1] = nstage[r] * sahaconst * u[r + 1] / u[r] * np.exp(-chiion[r] / keVT)
        ntotal = np.sum(nstage)
        nstagerel = nstage / ntotal
    return nstagerel[ionstage-1]


def sahabolt_E(temp, elpress, ion, level):
    return saha_E(temp, elpress, ion) * boltz_E(temp, ion, level)


def plot_payne():
    temp = np.arange(0,30001,1000)
    #print temp
    pop = np.zeros((5,31))
    for T in np.arange(1,31):
        for r in np.arange(1,5):
            pop[r,T] = sahabolt_E(temp[T],131.,r,1)
    labellst = ['ground stage', 'first ion stage', 'second ion stage', 'third ion stage']
    #print pop
    plt.figure(0)
    # ground-state plot
    for i in range(1,5):
        plt.plot(temp,pop[i,:], label=labellst[i-1])
    plt.xlabel('temperature', size=14)
    plt.ylabel('population', size=14)
    plt.yscale('log')
    plt.ylim([1e-3, 1.1])
    plt.legend(loc='best')
    plt.show()
"""
Hydrogen:
"""
def sahabolt_H(temp,elpress,level):
    kerg = 1.380658e-16
    keV = 8.61734e-5 #Boltzmann constant (eV/deg)
    h = 6.62607e-27 #Planck constant (erg s)
    elmass = 9.109390e-28 #electron mass (g)
    keVT = keV*temp
    kergT = kerg*temp
    eldens = elpress/kergT
    # energy levels and weights for hydrogen
    nrlevels = 100
    # reasonable partition function cut-off value
    g = np.zeros((2,nrlevels))
    # declarations weights (too many for proton)

    chiexc = np.zeros((2,nrlevels))
    # declaration excitation energies (idem)
    for s in range(nrlevels):
        g[0,s] = 2.*(s+1.)**2.
        # statistical weights
        chiexc[0,s] = 13.598*(1.-1./(s+1.)**2.)
    # excitation weights
    g[1,0] = 1.
    # statistical weights free proton
    chiexc[1,0] = 0.
    # partition functions
    u = np.zeros([2])
    for s in range(nrlevels):
        u[0] = u[0] + g[0,s]*np.exp(-chiexc[0,s]/keVT)
    u[1] = g[1,0]
    # Saha
    sahaconst = (2*np.pi*elmass*kergT /(h*h))**(1.5)*2./eldens
    nstage = np.zeros(2)
    nstage[0] = 1.
    nstage[1] = nstage[0] * sahaconst * u[1]/u[0] * np.exp(-13.598/keVT)
    ntotal = np.sum(nstage)
    # sum both stages = total hydrogen density
    # Boltzmann
    nlevel = nstage[0]*g[0,level-1]/u[0]*np.exp(-chiexc[0,level-1]/keVT)
    nlevelrel = nlevel/ntotal
    # fraction of total hydrogen density
    for s in range(6):
        print(s+1, g[0,s], chiexc[0,s], g[0,s]*np.exp(-chiexc[0,s]/keVT))
    #print
    for s in range(0,nrlevels,10):
        print (s+1, g[0,s], chiexc[0,s], g[0,s]*np.exp(-chiexc[0,s]/keVT))


    return nlevelrel

if __name__ == '__main__':
    #print(sahabolt_H(5000,1e2,1))
