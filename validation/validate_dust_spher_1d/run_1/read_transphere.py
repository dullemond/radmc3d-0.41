import numpy as np

def read_envelope(file='envstruct.dat'):
    with open(file,'r') as f:
        str  = f.readline()
        nr   = int(str)
        data = np.loadtxt(f)
        r, rho, temp = data.T
    return r, rho, temp

def read_spectrum(file='spectrum.dat'):
    cc  = 2.9979245800000e10      # Light speed             [cm/s]
    with open(file,'r') as f:
        str  = f.readline()
        nrfr = int(str)
        data = np.loadtxt(f)
        freq, spectrum = data.T
    lam = 1e4*cc/freq
    return lam, spectrum
