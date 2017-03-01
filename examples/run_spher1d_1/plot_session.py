from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import math

# Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy.image import *    
from radmc3dPy.analyze import *  
from radmc3dPy.natconst import * 

#
# Make sure to have run
#
#   radmc3d mctherm
#
# to compute the dust temperature before you run this plotting session.
#
# Now plot the temperature profile
#
a    = readData(dtemp=True,binary=False)
r    = a.grid.x[:]
temp = a.dusttemp[:,0,0,0]
plt.figure(1)
plt.plot(r/au,temp)
plt.xlabel('r [au]')
plt.ylabel('T [K]')
plt.show()

#
# Now make sure to have run 
#
#   radmc3d sed
#
# to get the spectral energy distribution
#
s    = readSpectrum()
plt.figure(2)
lammic = s[:,0]
flux   = s[:,1]
nu     = 1e4*cc/lammic
nufnu  = nu*flux
nulnu  = nufnu*4*math.pi*pc*pc
plt.plot(lammic,nulnu/ls)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\lambda$ [$\mu$m]')
plt.ylabel(r'$\nu L_{\nu}$ [$L_{\odot}$]')
plt.axis([1e-1,1e4, 1e-6,1e1])
plt.show()

#
# Now make sure to have run 
#
#   radmc3d image lambda 10
#
# to get the spectral energy distribution
#
im   = readImage()
plt.figure(3)
plotImage(im,log=True,maxlog=6,au=True,bunit='inu')
plt.show()
