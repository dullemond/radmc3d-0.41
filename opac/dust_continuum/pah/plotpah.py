from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math
from makepah import *

mp   = 1.6726e-24     # Mass of proton          [g]

nf   = 1000
lam0 = 0.03
lam1 = 1000.
#lam0 = 2
#lam1 = 40.
lam  = lam0 * (lam1/lam0)**(np.linspace(0.,1.,nf))

nc   = 40.
nh   = 16.

kappa_neutral = pah(nc,nh,0,lam)
kappa_charged = pah(nc,nh,1,lam)

sig_neutral = kappa_neutral * (12*mp)
sig_charged = kappa_charged * (12*mp)

#
# Reproducing Figure 2 of Li & Draine 2001
#
plt.figure(1)
plt.plot(lam,sig_neutral,label='neutral')
plt.plot(lam,sig_charged,label='charged')
plt.xscale('log')
plt.yscale('log')
plt.axis([2,40,1e-22,1e-19])   # Range of Fig. 2 top of LD2001
plt.show()

plt.figure(2)
plt.plot(lam,sig_neutral,label='neutral')
plt.plot(lam,sig_charged,label='charged')
plt.xscale('log')
plt.yscale('log')
plt.axis([0.03,300,1e-23,1e-16])   # Range of Fig. 2 bottom of LD2001
plt.show()


