#
# Import NumPy for array handling
#
import numpy as np
import scipy
import math
from matplotlib import cm
import matplotlib.pyplot as plt
import os
import read_transphere as tp
from radmc3dPy import *
from radmc3dPy.natconst import *
#
# Import plotting libraries (start Python with ipython --matplotlib)
#
#from mpl_toolkits.mplot3d import axes3d
#from matplotlib import pyplot as plt
#
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
cc  = 2.9979245800000e10      # Light speed             [cm/s]
#
# Monte Carlo parameters
#
nphot    = 100000
#
# Grid parameters
#
nr       = 100
#
# Model parameters
#
rin      = 5*au
rout     = 100*au
rho0     = 3e-15
prho     = -2.e0
#
# Star parameters
#
mstar    = ms
rstar    = rs
tstar    = ts
#
# Parameters for transphere
#
nriter   = 30             # Maximum nr of iterations
convcrit = 0.001          # Convergence criterion
ncst     = 10             # Nr of rays for star
ncex     = 30             # Nr of rays between star and Rin
ncnr     = 1              # Nr of rays per radial grid point
itypemw  = 1              # Type of mu weighting
idump    = 1              # Dump convergence history
#
# Make the coordinates
#
# Note: The way the xi grid is made is slightly non-standard, but is
#       done this way to be consistent with problem_setup.pro (the IDL version)
#
ri       = rin * (rout/rin)**(np.linspace(0.,nr-1,nr+1)/(nr-1.0))
rc       = 0.5e0 * ( ri[0:nr] + ri[1:nr+1] )
#
# Make the dust density model
#
rr       = rc
rhod     = rho0 * (rr/au)**prho
#
# Read the dust opacity file
#
with open('dustkappa_silicate.inp','r') as f:
    str = f.readline()
    str = f.readline()
    opac_nf = int(str)
    data = np.loadtxt(f)
    opac_lam, opac_kabs, opac_ksca = data.T
#
# Write the wavelength_micron.inp file
#
lam1     = 0.1e0
lam2     = 7.0e0
lam3     = 25.e0
lam4     = 1.0e4
n12      = 20
n23      = 100
n34      = 30
lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
lam      = np.concatenate([lam12,lam23,lam34])
freq     = 1e4*cc/lam
nlam     = lam.size
#
#-----------------------------------------------------------------
#                 FIRST THE TRANSPHERE MODEL
#-----------------------------------------------------------------
#
# Write the frequency file (for TRANSPHERE)
# NOTE: File name slightly changed, also in ../src_transphere/
#
with open('frequency_hz.inp','w+') as f:
    f.write('%d\n'%(nlam))
    np.savetxt(f,freq.T,fmt=['%13.6e'])
#
# Write the dustopac file (for TRANSPHERE)
#
kabs = scipy.interp(lam,opac_lam,opac_kabs)
with open('dustopac_1.inp','w+') as f:
    f.write('%d %d\n\n'%(nlam,1))
    for i in range(nlam):
        f.write('%13.6e\n'%(kabs[i]))
    f.write('\n')
    for i in range(nlam):
        f.write('%13.6e\n'%(0.0))
#
# Write the stars.inp file (for TRANSPHERE)
#
with open('starinfo.inp','w+') as f:
    f.write('1\n')
    f.write('%13.6e\n'%(rstar))
    f.write('%13.6e\n'%(mstar))
    f.write('%13.6e\n'%(tstar))
#
# Write the transphere.inp file (for TRANSPHERE)
#
with open('transphere.inp','w+') as f:
    f.write('2\n')
    f.write('%d\n'%(nriter))
    f.write('%13.6e\n'%(convcrit))
    f.write('%d\n'%(ncst))
    f.write('%d\n'%(ncex))
    f.write('%d\n'%(ncnr))
    f.write('%d\n'%(itypemw))
    f.write('%d\n'%(idump))
#
# Write the envelope structure (for TRANSPHERE)
#
with open('envstruct.inp','w+') as f:
    f.write('%d\n\n'%(nr))
    for ir in range(nr):
        f.write('%13.6e %13.6e %13.6e\n'%(rc[ir],rhod[ir],0.0))
#
# Dust opacity control file (for TRANSPHERE)
#
with open('dustopac.inp','w+') as f:
    f.write('1               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('-1              Way in which this dust species is read\n')
    f.write('1               Extension of name of dustopac_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
#
# Now run TRANSPHERE
#
os.system('../src_transphere/transphere')
#
# Read the TRANSPHERE result 
#
tp_r, tp_rhod, tp_temp = tp.read_envelope()
tp_lam, tp_spectrum    = tp.read_spectrum()
tp_nu                  = 1e4*cc/tp_lam
tp_nufnu               = tp_nu*tp_spectrum
tp_nulnu               = tp_nufnu*4*math.pi*pc*pc

#
#-----------------------------------------------------------------
#              NOW THE EQUIVALENT RADMC-3D MODEL
#-----------------------------------------------------------------
#
nx       = nr
ny       = 1
nz       = 1
xi       = ri
yi       = np.array([0.,math.pi])
zi       = np.array([0.,math.pi*2])
xc       = 0.5e0 * ( xi[0:nx] + xi[1:nx+1] )
pstar    = [0.,0.,0.]
#
# Write the wavelength file (for RADMC-3D)
#
with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
#
# Write the stars.inp file
#
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
    f.write('\n%13.6e\n'%(-tstar))
#
# Write the grid file
#
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('100\n')                     # Coordinate system
    f.write('0\n')                       # gridinfo
    f.write('1 0 0\n')                   # Include x,y,z coordinate
    f.write('%d %d %d\n'%(nx,ny,nz))     # Size of grid
    np.savetxt(f,xi.T,fmt=['%13.6e'])    # X coordinates (cell walls)
    np.savetxt(f,yi.T,fmt=['%13.6e'])    # Y coordinates (cell walls)
    np.savetxt(f,zi.T,fmt=['%13.6e'])    # Z coordinates (cell walls)
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = rhod.ravel(order='F')         # Create a 1-D view, fortran-style indexing
    np.savetxt(f,data.T,fmt=['%13.6e'])  # The data
#
# Dust opacity control file
#
with open('dustopac.inp','w+') as f:
    f.write('2               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('silicate        Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
#
# Now run RADMC-3D to determine the dust temperature
#
os.system('radmc3d mctherm')
#
# Read dust temperature
#
rmc_temp = analyze.readData(dtemp=True,binary=False)
#
# Now run RADMC-3D to determine the SED
#
os.system('radmc3d sed')
#
# Read the SED
#
rmc_sed   = analyze.readSpectrum()
rmc_lam   = rmc_sed.T[0,:]
rmc_nu    = 1e4*cc/rmc_lam
rmc_nufnu = rmc_nu*rmc_sed.T[1,:]
rmc_nulnu = rmc_nufnu*4*math.pi*pc*pc


#-----------------------------------------------------------------
#                         COMPARE THE TWO
#-----------------------------------------------------------------
plt.figure(1)
plt.plot(tp_r/au,tp_temp,label='TRANSPHERE')
plt.plot(rmc_temp.grid.x/au,rmc_temp.dusttemp[:,0,0,0],label='RADMC-3D')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('r [au]')
plt.ylabel('T [K]')
plt.autoscale(enable=True,axis='x',tight=True)
plt.autoscale(enable=True,axis='y',tight=True)
plt.legend(loc='upper right')
#plt.savefig('dust_temperature_comparison.png')
plt.show()

plt.figure(2)
plt.plot(tp_lam,tp_nulnu/ls,label='TRANSPHERE')
plt.plot(rmc_lam,rmc_nulnu/ls,label='RADMC-3D')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\lambda$ [$\mu$m]')
plt.ylabel(r'$\nu L_{\nu}$ [$L_{\odot}$]')
plt.axis([1e-1,1e3, 1e-6,1e1])
plt.legend(loc='upper left')
#plt.savefig('sed_comparison.png')
plt.show()
