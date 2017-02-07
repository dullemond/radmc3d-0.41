#
# First do, for example:
#   radmc3d mctherm
#   radmc3d image lambda 0.4e3 theta 45 phi 45
#
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import math
import numpy as np
from radmc3dPy.image import *
from plotpoldir import *
au  = 1.49598e13     # Astronomical Unit       [cm]

a = readImage()
plotImage(a,cmap=cm.hot,au=True)
plotpoldir(a.x/au,a.y/au,a.image[:,:,0,0],a.image[:,:,1,0],a.image[:,:,2,0])



