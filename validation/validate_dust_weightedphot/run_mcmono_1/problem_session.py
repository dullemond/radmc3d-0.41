import numpy as np
import math
from matplotlib import cm
import matplotlib.pyplot as plt
import os
from radmc3dPy import *
from radmc3dPy.natconst import *
from solvetdust import *
#
# Monte Carlo parameters
#
nphot    = 1000000
#
# Grid parameters
#
nx       = 128
ny       = 32
nz       = 32
sizex    = 40*au
sizey    = 10*au
sizez    = 10*au
#
# Model parameters
#
radius   = 5*au
rho0     = 1e-23   # Optically thin
#rho0     = 1e-16   # Optically thick
temp0    = 10.     # A dummy temperature (since we test mcmono)
sx1      = -30*au
sx2      = +30*au
#
# Star parameters
#
mstar1    = ms
rstar1    = rs
tstar1    = ts
pstar1    = np.array([sx1,0.,0.])
mstar2    = ms
rstar2    = rs
tstar2    = 0.5*ts
pstar2    = np.array([sx2,0.,0.])
#
# Wavelength for mcmono
#
lammono   = 1.
#
# Make the coordinates
#
xi       = np.linspace(-sizex,sizex,nx+1)
yi       = np.linspace(-sizey,sizey,ny+1)
zi       = np.linspace(-sizez,sizez,nz+1)
xc       = 0.5 * ( xi[0:nx] + xi[1:nx+1] )
yc       = 0.5 * ( yi[0:ny] + yi[1:ny+1] )
zc       = 0.5 * ( zi[0:nz] + zi[1:nz+1] )
#
# Make the dust density model
#
qq       = np.meshgrid(xc,yc,zc,indexing='ij')
xx       = qq[0]
yy       = qq[1]
zz       = qq[2]
rr       = np.sqrt(xx**2+yy**2+zz**2)
rhod     = rho0 + rr*0.0
tempd    = temp0 + rr*0.0
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
nlam     = lam.size
#
# Write the wavelength file
#
with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
#
# Write the mcmono wavelength file
#
with open('mcmono_wavelength_micron.inp','w+') as f:
    f.write('1\n')
    f.write('%13.6e\n'%(lammono))
#
# Write the stars.inp file
#
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('2 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar1,mstar1,pstar1[0],pstar1[1],pstar1[2]))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar2,mstar2,pstar2[0],pstar2[1],pstar2[2]))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
    f.write('\n%13.6e\n'%(-tstar1))
    f.write('\n%13.6e\n'%(-tstar2))
#
# Write the grid file
#
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('0\n')                       # Coordinate system
    f.write('0\n')                       # gridinfo
    f.write('1 1 1\n')                   # Include x,y,z coordinate
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
# Write the temperature file
#
with open('dust_temperature.dat','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nx*ny*nz))           # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = tempd.ravel(order='F')        # Create a 1-D view, fortran-style indexing
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

#----------------------------------------------------------------------

#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot_mono = %d\n'%(nphot))
    f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
    f.write('mc_weighted_photons = 0')
#
# Now run RADMC-3D to determine the mean intensity
#
#os.system('radmc3d mcmono')
os.system('radmc3d mcmono setthreads 4')   # Run with OpenMP parallellization
#
# Read mean intensity
#
with open('mean_intensity.out','r') as f:
    iformat = int(f.readline())
    assert iformat==2, "Cannot read this format"
    nrcells = int(f.readline())
    nf = int(f.readline())
    freq = np.zeros(nf)
    for inu in range(nf):
        freq[inu] = float(f.readline())
    data = np.loadtxt(f)
    meanint_nowgt = np.reshape(data,(nx,ny,nz,nf),order='F')
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot_mono = %d\n'%(nphot))
    f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
    f.write('mc_weighted_photons = 1')
#
# Now run RADMC-3D to determine the mean intensity
#
#os.system('radmc3d mcmono')
os.system('radmc3d mcmono setthreads 4')   # Run with OpenMP parallellization
#
# Read mean intensity
#
with open('mean_intensity.out','r') as f:
    iformat = int(f.readline())
    assert iformat==2, "Cannot read this format"
    nrcells = int(f.readline())
    nf = int(f.readline())
    freq = np.zeros(nf)
    for inu in range(nf):
        freq[inu] = float(f.readline())
    data = np.loadtxt(f)
    meanint_wgt = np.reshape(data,(nx,ny,nz,nf),order='F')
#
#-----------------------------------------------------------------
#      COMPUTE MEAN INTENSITY FOR BOTH SOURCES
#-----------------------------------------------------------------
x        = xc
nx       = x.size
bnu1     = bplanck(freq,tstar1)
bnu2     = bplanck(freq,tstar2)
mi1      = (bnu1/4)*(rstar1/(x-sx1))**2
mi2      = (bnu2/4)*(rstar2/(x-sx2))**2
mi       = mi1+mi2
# 
#-----------------------------------------------------------------
#    COMPARE THE TWO, AND OVERPLOT OPTICALLY THIN TEMPERATURES
#-----------------------------------------------------------------
# 
plt.figure(1)
plt.plot(x/au,meanint_wgt[:,16,16,0],label='Weighted')
plt.plot(x/au,meanint_nowgt[:,16,16,0],label='Unweighted')
plt.plot(x/au,mi1,label='MI1')
plt.plot(x/au,mi2,label='MI2')
plt.plot(x/au,mi,label='MItot')
plt.yscale('log')
plt.xlabel('x [au]')
plt.ylabel(r'$J_\nu [\mathrm{erg/s/cm}^2\mathrm{/Hz/ster}]$')
plt.legend(loc='lower left')
#plt.savefig('mean_intensity_comparison.png')
plt.show()


