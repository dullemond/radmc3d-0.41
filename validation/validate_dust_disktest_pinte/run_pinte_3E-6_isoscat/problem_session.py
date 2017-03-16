import numpy as np
import math
from matplotlib import cm
import matplotlib.pyplot as plt
import scipy as sp
from scipy import interpolate
import os
from radmc3dPy import *
from radmc3dPy.natconst import *


rstar  = 2.e0 * rs      # radius of the star
mstar  = 1.0e0 * ms     # mass of the star
tstar  = 4000e0         # temperature of the star

nr     = 120		# number of points in radius
nt     = 60		# number of points in theta in regular grid
ntvac  = 20             # number of extra points in theta for vacuum
hrgrid = 0.4		# Grid span in Theta (measured negatively from equator)

nlevr  = 4		# higher resolution at level #
nspanr = 3		# number of steps to resolve
nstepr = 2		# number of sub-steps to create

r0     = 100*au         # 
rin    = 0.1e0 * au	# inner radius = 1 AU	
rout   = 4.e2 * au	# outer radius = 1000 AU
mdust  = 3e-6*ms        # Dust mass
plsig  = -1.5e0         # Powerlaw for Sigma_dust 
h0     = 10*au          # Disk height at r0
plh    = 1.125          # Powerlaw for 

nphot  = 100000000        # Nr of photon packages

#
# Find the opacity and wavelength lambda 
#
with open('dustkappa_1.inp') as f:
    iformat = int(f.readline())
    nlam    = int(f.readline())
    data    = np.loadtxt(f)
    lam, kabs, ksca, gsca = data.T
#
# Get the opacity at 0.81 micron
#
f     = interpolate.interp1d(lam,kabs+ksca)
kappa = f([0.81])[0]
#
# Make the R-grid
#
nextra = nspanr*nlevr*(2**nstepr-1)
assert nextra <= nr-5, "Sorry not enough nr for this refinement level"
rc     = rin*(rout/rin)**(np.linspace(0.0,1.0,nr-nextra))
for ilev in range(1,nlevr+1):
    for ispan in range(nspanr,0,-1):
        rins = np.zeros(2**nstepr-1)
        fr   = (rc[ispan]/rc[ispan-1])**(0.5**nstepr)
        for i in range(1,2**nstepr):
            rins[i-1]=rc[ispan-1]*fr**i
        r  = rc
        rc = np.append(r[0:ispan],rins)
        rc = np.append(rc,r[ispan:])
#
# Make the Theta-grid:
#
thmax  = math.pi/2.e0
thmin  = math.pi/2.e0 - hrgrid
theta  = (thmax-thmin)*np.linspace(0.,nt-1.,nt)/(nt-1.0+0.50)+thmin
the    = (thmin)*(np.linspace(0.,ntvac-1.,ntvac)+0.5)/(1.0*ntvac)
theta  = np.append(the,theta)
ntextr = 0
thetac = np.append(theta,math.pi-theta[::-1])
#
# Make the cell walls
#
nnr    = rc.size
nnt    = thetac.size
ri     = np.append(rc[0],np.sqrt(rc[1:nnr]*rc[0:nnr-1]))
ri     = np.append(ri,rc[nnr-1])
ti     = np.append(0.,0.5e0*(thetac[1:nnt]+thetac[0:nnt-1]))
ti     = np.append(ti,math.pi)
phii   = np.array([0.,2*math.pi])
#
# Compute the cell volumes
#
vr     = ri[1:]**3-ri[0:-1]**3
vt     = np.abs(np.cos(ti[0:-1])-np.cos(ti[1:]))
vvr, vvt = np.meshgrid(vr,vt,indexing='ij')
vol    = (2.*math.pi/3.)*vvr*vvt
#
# Make 2-D grid
#
qq     = np.meshgrid(rc,thetac,indexing='ij')
rr     = qq[0]
tt     = qq[1]
#
# Set up the density model
#
rcyl   = np.abs(rr*np.sin(tt))
zz     = np.abs(rr*np.cos(tt))
hh     = h0 * (rcyl/r0)**plh
sig0   = 1.0
rho    = (sig0/(math.sqrt(2*math.pi)*hh))*np.exp(-zz**2/(2*hh**2))*(rcyl/au)**plsig
#
# Renormalize to get the right mass
#
dum    = rho*vol
mas    = dum.sum()
rho    = rho * mdust / mas
#
# Write the grid file
#
with open('amr_grid.inp','w+') as f:
    f.write('1\n')                       # iformat
    f.write('0\n')                       # AMR grid style  (0=regular grid, no AMR)
    f.write('100\n')                     # Coordinate system
    f.write('0\n')                       # gridinfo
    f.write('1 1 0\n')                   # Include x,y,z coordinate
    f.write('%d %d %d\n'%(nnr,nnt,1))    # Size of grid
    np.savetxt(f,ri.T,fmt=['%13.6e'])    # R coordinates (cell walls)
    np.savetxt(f,ti.T,fmt=['%13.6e'])    # Theta coordinates (cell walls)
    np.savetxt(f,phii.T,fmt=['%13.6e'])  # Phi coordinates (cell walls)
#
# Write the density file
#
with open('dust_density.inp','w+') as f:
    f.write('1\n')                       # Format number
    f.write('%d\n'%(nnr*nnt))            # Nr of cells
    f.write('1\n')                       # Nr of dust species
    data = rho.ravel(order='F')          # Create a 1-D view, fortran-style indexing
    np.savetxt(f,data.T,fmt=['%13.6e'])  # The data
#
# Write the wavelength file
#
with open('wavelength_micron.inp','w+') as f:
    f.write('%d\n'%(nlam))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
#
# Make the stars.inp
#
with open('stars.inp','w+') as f:
    f.write('2\n')
    f.write('1 %d\n\n'%(nlam))
    f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,0.,0.,0.))
    np.savetxt(f,lam.T,fmt=['%13.6e'])
    f.write('\n%13.6e\n'%(-tstar))
#
# Read Pinte density
#
#with open('density_radial_from_pinte','r') as f:
#    npr   = int(f.readline())
#    datar = np.loadtxt(f)
#with open('density_vertical_from_pinte','r') as f:
#    npz   = int(f.readline())
#    dataz = np.loadtxt(f)
#
# Dust opacity control file
#
with open('dustopac.inp','w+') as f:
    f.write('2               Format number of this file\n')
    f.write('1               Nr of dust species\n')
    f.write('============================================================================\n')
    f.write('1               Way in which this dust species is read\n')
    f.write('0               0=Thermal grain\n')
    f.write('1               Extension of name of dustkappa_***.inp file\n')
    f.write('----------------------------------------------------------------------------\n')
#
# Write the radmc3d.inp control file
#
with open('radmc3d.inp','w+') as f:
    f.write('nphot = %d\n'%(nphot))
    f.write('modified_random_walk = 1\n')
    f.write('scattering_mode_max = 1\n')

print "WARNING: Isotropic scattering!"

#----------------------------------------------------------------------

imid = nnt/2

#
# Compute the optical depth
#
dr    = ri[1:]-ri[0:-1]
tau   = (rho[:,imid]*dr*kappa).sum()
print 'Optical depth along midplane = ',tau

#
# Plot comparison of density with Pinte
#
#plt.figure(1)
#plt.plot(rc/au,rho[:,imid],label='RADMC-3D')
#plt.plot(datar.T[0],datar.T[1]*0.01*(mdust/(3e-8*ms)),label='Pinte Paper')
#plt.xscale('log')
#plt.yscale('log')
#plt.xlabel('r [au]')
#plt.ylabel(r'$\rho_d$ [g/cm$^3$]')
#plt.legend(loc='lower right')
#
#f        = interpolate.interp2d(thetac,rc/au,rho[:,:])
#rcyl     = 200
#zcyl     = (math.pi/2-thetac)*rcyl
#rsph     = np.sqrt(rcyl**2+zcyl**2)
#tsph     = math.pi/2-np.arctan(zcyl/rcyl)
#rhodust200  = np.zeros(nnt)
#for iz in range(nnt): rhodust200[iz] = f(tsph[iz],rsph[iz])
#
#plt.figure(2)
#plt.plot(math.pi/2-thetac,rhodust200,label='RADMC-3D')
#plt.plot(dataz.T[0]/200.,dataz.T[1]*0.01*(mdust/(3e-8*ms)),label='Pinte Paper')
#plt.yscale('log')
#plt.axis([-1., 1., 1e-30, 1e-19])
#plt.xlabel('z/r')
#plt.ylabel(r'$\rho_d$ [g/cm$^3$]')
#plt.legend(loc='upper right')
#plt.show()

#----------------------------------------------------------------------

#
# Now run RADMC-3D to determine the dust temperature
#
os.system('radmc3d mctherm')
#os.system('radmc3d mctherm setthreads 4')   # Run with OpenMP parallellization

#----------------------------------------------------------------------

with open('midplane_temp_tau1e5_isoscat_mcfost','r') as f:
    data = np.loadtxt(f)
pinte_r, pinte_mid_temp = data.T

with open('vert0.2_temp_tau1e5_isoscat_mcfost','r') as f:
    data = np.loadtxt(f)
pinte_vert_02_z, pinte_vert_02_temp = data.T

with open('vert200_temp_tau1e5_isoscat_mcfost','r') as f:
    data = np.loadtxt(f)
pinte_vert_200_z, pinte_vert_200_temp = data.T


a=analyze.readData(dtemp=True,ddens=True,binary=False)

imid = a.grid.y.size/2
rin=0.1*au
plt.figure(3)
plt.plot((a.grid.x-rin)/au,a.dusttemp[:,imid,0,0],label='RADMC-3D')
plt.plot(pinte_r-rin/au,pinte_mid_temp,label='MCFOST')
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$r-r_{\mathrm{in}}$ [au]')
plt.ylabel(r'$T_d$ [K]')
plt.axis([0.00001, 400., 10, 2000.])
plt.legend(loc='lower left')
#plt.savefig('compare_tdust_midplane.png')
plt.show()

f     = interpolate.interp1d(rc/au,np.linspace(0,nnr-1,nnr))
ir02  = int(f([0.2])[0])+1
ir200 = int(f([200])[0])

f        = interpolate.interp2d(thetac,rc/au,a.dusttemp[:,:,0,0])
rcyl02   = 0.2
zcyl02   = (math.pi/2-thetac)*rcyl02
rsph     = np.sqrt(rcyl02**2+zcyl02**2)
tsph     = math.pi/2-np.arctan(zcyl02/rcyl02)
tdust02  = np.zeros(nnt)
for iz in range(nnt): tdust02[iz] = f(tsph[iz],rsph[iz])
rcyl200  = 200
zcyl200  = (math.pi/2-thetac)*rcyl200
rsph     = np.sqrt(rcyl200**2+zcyl200**2)
tsph     = math.pi/2-np.arctan(zcyl200/rcyl200)
tdust200 = np.zeros(nnt)
for iz in range(nnt): tdust200[iz] = f(tsph[iz],rsph[iz])

with open('vert0.2_temp_tau1e5_isoscat_radmc3d','w') as f:
    for iz in range(0,nnt/2):
        f.write("%13.6e %13.6e\n"%(zcyl02[iz],tdust02[iz]))

with open('vert200_temp_tau1e5_isoscat_radmc3d','w') as f:
    for iz in range(0,nnt/2):
        f.write("%13.6e %13.6e\n"%(zcyl200[iz],tdust200[iz]))

plt.figure(4)
plt.plot(zcyl02/0.2,tdust02,label='RADMC-3D')
plt.plot(pinte_vert_02_z/(0.2),pinte_vert_02_temp,label='MCFOST')
plt.yscale('log')
plt.xlabel(r'$z/r$')
plt.ylabel(r'$T_d$ [K]')
plt.axis([-0.5,0.5, 100., 1000.])
plt.legend(loc='lower left')
#plt.savefig('compare_tdust_vert_0.2au.png')
plt.show()

plt.figure(5)
plt.plot(zcyl200/200,tdust200,label='RADMC-3D')
plt.plot(pinte_vert_200_z/(200.),pinte_vert_200_temp,label='MCFOST')
plt.yscale('log')
plt.xlabel(r'$z/r$')
plt.ylabel(r'$T_d$ [K]')
plt.axis([-0.5,0.5, 10., 100.])
plt.legend(loc='lower left')
#plt.savefig('compare_tdust_vert_200au.png')
plt.show()

