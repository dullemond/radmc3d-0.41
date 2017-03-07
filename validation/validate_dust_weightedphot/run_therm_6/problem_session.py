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
rho0     = 1e-19   # Optically thin
#rho0     = 1e-16   # Optically thick
sx1      = -30*au
sx2      = +30*au
#
# Star parameters
#
mstar    = ms
rstar    = rs
tstar    = ts
pstar    = np.array([sx1,0.,0.])
#
# Heat source: model this as a blob with a luminosity similar
# to a star 
#
ttempl   = 0.5*ts
#rrtempl  = 0.1*au    # Narrow Gaussian ---> Behaves like a point source (like run_1)
rrtempl  = 2*au      # Wide Gaussian ---> Behaves like a smoothed-out stellar population
rtempl   = rs
ltempl   = ls * (rtempl/rs)**2 * (ttempl/ts)**4
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
rr       = np.sqrt((xx-sx2)**2+yy**2+zz**2)
rhod     = rho0 + rr*0.0
#
# Make the internal heat source distribution model
#
heatsrc  = np.exp(-rr**2/(2*rrtempl**2))
volcell  = (xi[1]-xi[0])**3
dummy    = heatsrc*volcell
tot      = dummy.sum()
heatsrc  = heatsrc * ltempl / tot
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
# Write the stars.inp file
#
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
    f.write('\n%13.6e\n'%(-tstar))
#
# Make the internal heat source
#
with open('heatsource.inp','w') as f:
    f.write('1\n')
    f.write('%d\n'%(nx*ny*nz))
    data = heatsrc.ravel(order='F')      # Create a 1-D view, fortran-style indexing
    np.savetxt(f,data.T,fmt=['%13.6e'])  # The data
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
    f.write('mc_weighted_photons = 0\n')
#
# Now run RADMC-3D to determine the dust temperature
#
#os.system('radmc3d mctherm')
os.system('radmc3d mctherm setthreads 4')   # Run with OpenMP parallellization
#
# Read dust temperature
#
temp_nowgt = analyze.readData(dtemp=True,binary=False)
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
    f.write('mc_weighted_photons = 1\n')
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
#      WHERE WE ASSUME THE SECOND SOURCE TO ALSO BE A POINT
#     SOURCE IN SPITE OF THE FACT THAT IT IS A SMOOTH SOURCE
#-----------------------------------------------------------------
x        = temp_wgt.grid.x
nx       = x.size
t1       = np.zeros(nx)
t2       = np.zeros(nx)
cc       = 2.9979245800000e10      # Light speed             [cm/s]
with open('dustkappa_silicate.inp','r') as f:
    iformat = int(f.readline())
    nf      = int(f.readline())
    data    = np.loadtxt(f)
    lam, kabs, ksca = data.T
nu = 1e4*cc/lam
qheat = heatsrc[:,16,16]/rhod[:,16,16] + 1.
for i in range(nx):
    dist1 = x[i]-sx1
    t1[i] = solvetdust_bbstar(nu,kabs,rstar,tstar,dist1)
    t2[i] = solvetdust_qheat(nu,kabs,qheat[i])

#-----------------------------------------------------------------
#                         COMPARE THE TWO
#-----------------------------------------------------------------

plt.figure(1)
plt.plot(temp_wgt.grid.x/au,temp_wgt.dusttemp[:,16,16],label='Weighted')
plt.plot(temp_nowgt.grid.x/au,temp_nowgt.dusttemp[:,16,16],label='Unweighted')
plt.plot(x/au,t1,label='T1')
plt.plot(x/au,t2,label='T2')
plt.xlabel('x [au]')
plt.ylabel('T [K]')
plt.legend(loc='upper right')
#plt.savefig('dust_temperature_comparison.png')
plt.show()


