import problem_setup as p
import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
from matplotlib import cm
from radmc3dPy.image import *    # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy.analyze import *  # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from readquant import *

icell = 0                        # This is a one-cell example, so set icell=0
q     = readquant(2)
temp  = q["temp"]
tdist = q["tdist"][icell,:]

plt.plot(temp,temp*tdist)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('T [K]')
plt.ylabel('T P(T)')
plt.axis([1e0, 1e4, 1e-8, 1e1])

plt.show()
