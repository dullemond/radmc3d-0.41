import numpy as np
from matplotlib import pyplot as plt
from radmc3dPy.image import *    # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory
from radmc3dPy.analyze import *  # Make sure that the shell variable PYTHONPATH points to the RADMC-3D python directory

#
# Make sure to have done the following beforhand:
#
#  First compile RADMC-3D
#  Then run:
#   python problem_setup.py
#   radmc3d mctherm
#

d=readData(dtemp=True)

plt.imshow(d.dusttemp[:,:,16,0].T)
#plt.imshow(d.dusttemp[:,16,:,0].T)

plt.show()
