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
nphot    = 10000000
#
# Grid parameters
#
nx       = 32
ny       = 32
nz       = 32
sizex    = 10*au
sizey    = 10*au
sizez    = 10*au
#
# Model parameters
#
radius   = 0.7*au
#rho0     = 1e-23   # Optically thin
rho0     = 1e-13   # Optically thick
#
# Background temperature
#
tbg       = 600.0
fdil      = 1.0
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
rhod     = rho0 * np.exp(-(rr**2/radius**2)/2.0)
#
# Make the wavelength grid
#
# lam1     = 0.1e0
# lam2     = 7.0e0
# lam3     = 25.e0
# lam4     = 1.0e4
# n12      = 20
# n23      = 100
# n34      = 30
# lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
# lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
# lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
# lam      = np.concatenate([lam12,lam23,lam34])
# nlam     = lam.size
#
# Or read the wavelength grid from the opacity file
#
cc       = 2.9979245800000e10      # Light speed             [cm/s]
with open('dustkappa_silicate.inp','r') as f:
    iformat = int(f.readline())
    nlam    = int(f.readline())
    data    = np.loadtxt(f)
    lam, kabs, ksca = data.T
#
# Write the wavelength file
#
with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
#
# Write the external radiation field
#
cc  = 2.9979245800000e10      # Light speed             [cm/s]
nu  = 1e4*cc/lam
jnu = fdil*bplanck(nu,tbg)
with open('external_source.inp','w+') as f:
    f.write('2\n')
    f.write('%d\n\n'%(nlam))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
    f.write('\n')
    np.savetxt(f,jnu.T,fmt=['%13.6e'])
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
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
    f.write('mc_weighted_photons = 1\n')
    f.write('iseed = -6521\n')
#
# Now run RADMC-3D to determine the dust temperature
#
#os.system('radmc3d mctherm')
os.system('radmc3d mctherm setthreads 4')   # Run with OpenMP parallellization
#
# Read dust temperature
#
temp_wgt = analyze.readData(dtemp=True,binary=False)
#
#-----------------------------------------------------------------
#      COMPUTE OPTICALLY THIN TEMPERATURE FOR BOTH SOURCES
#-----------------------------------------------------------------
x        = temp_wgt.grid.x
nx       = x.size
t2       = np.zeros(nx)
cc       = 2.9979245800000e10      # Light speed             [cm/s]
with open('dustkappa_silicate.inp','r') as f:
    iformat = int(f.readline())
    nf      = int(f.readline())
    data    = np.loadtxt(f)
    lam, kabs, ksca = data.T
nu = 1e4*cc/lam
jnu = fdil*bplanck(nu,tbg)
t2 = np.zeros(nx) + solvetdust_jnu(nu,kabs,jnu)

#-----------------------------------------------------------------
#    COMPARE THE TWO, AND OVERPLOT OPTICALLY THIN TEMPERATURES
#-----------------------------------------------------------------

plt.figure(1)
plt.plot(temp_wgt.grid.x/au,temp_wgt.dusttemp[:,16,16],label='Weighted')
plt.plot(x/au,t2,label='T2')
plt.xlabel('x [au]')
plt.ylabel('T [K]')
plt.legend(loc='lower right')
plt.axis([-10, 10, 0, 800])
plt.savefig('dust_temperature_comparison.png')
plt.show()


