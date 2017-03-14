from makedustopac import *
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy import interpolate

agraincm  = 1 * 1e-4     # Grain size in cm
optconst  = "Draine_Si_sUV" # The optical constants name
matdens   = 3.5           # The material density in gram / cm^3
extrapol  = True          # Extrapolate optical constants beyond its wavelength grid, if necessary
verbose   = False         # If True, then write out status information
ntheta    = 181           # Number of scattering angle sampling points
chop      = 1             # Chop off forward scattering and renormalize

optconstfile= optconst+'.lnk'

#
# Set up a wavelength grid upon which we want to compute the opacities
#
with open(optconstfile,'r') as f:
    data = np.loadtxt(f)
lamcm = data.T[0] * 1e-4
#
# Set up an angular grid for which we want to compute the scattering matrix Z
#
theta     = np.linspace(0.,180.,ntheta)

#-----------------------------------------------------------------------------
#
# Now make the opacity with the bhmie code
#
print "Running the code. Please wait..."
opac       = compute_opac_mie(optconstfile,matdens,agraincm,lamcm,theta=theta,
                              extrapolate=extrapol,chopforward=chop,verbose=verbose)
#
# Now write it out to a RADMC-3D opacity file:
#
# ...The full scattering matrix file
#
print "Writing the opacity to scatmat file"
write_radmc3d_scatmat_file(opac,optconst)
#
# ...Only the opacity file with simple scattering info
#
print "Writing the opacity to kappa file"
write_radmc3d_kappa_file(opac,optconst)
#
# Now that RADMC-3D does not like it when both files are there, so
# you must choose whether you want to do the full scattering or not
# (advice: yes, do the full scattering, which is safer as it is more
# accurate).
#

#-----------------------------------------------------------------------------

plt.figure(1)
plt.plot(opac["lamcm"]*1e4,opac["kabs"],label='abs')
plt.plot(opac["lamcm"]*1e4,opac["kscat"],label='scat')
plt.plot(opac["lamcm"]*1e4,opac["kabs"]+opac["kscat"],label='tot')
plt.xlabel(r'$\lambda [\mu\mathrm{m}]$')
plt.ylabel(r'$\kappa [\mathrm{cm}^2/\mathrm{g}]$')
plt.xscale('log')
plt.yscale('log')
plt.axis([1e-1, 1e3, 1e-2, 2e4])
plt.legend(loc='lower left')
#plt.savefig('opacity_kabs_kscat.png')
plt.show()

plt.figure(2)
plt.plot(opac["lamcm"]*1e4,opac["kscat"]/(opac["kabs"]+opac["kscat"]),label='albedo')
plt.plot(opac["lamcm"]*1e4,opac["gscat"],label='g')
plt.xlabel(r'$\lambda [\mu\mathrm{m}]$')
plt.xscale('log')
plt.axis([1e-1, 1e2, 0, 1])
plt.legend(loc='lower left')
#plt.savefig('opacity_albedo_and_g.png')
plt.show()

nlam   = opac["lamcm"].size
f      = interpolate.interp1d(opac["lamcm"],np.linspace(0,nlam-1,nlam),fill_value='extrapolate')
ilam   = f([1e-4])[0]
iilam  = int(ilam)
epslam = ilam-iilam
zscat  = (1.-epslam)*opac["zscat"][iilam] + epslam*opac["zscat"][iilam+1]
lam0cm = ( (1.-epslam)*opac["lamcm"][iilam] + epslam*opac["lamcm"][iilam+1] ) *1e4

plt.figure(3)
plt.plot(opac["theta"],zscat.T[0],label='Z11')
plt.xlabel(r'$\lambda [\mu\mathrm{m}]$')
plt.ylabel(r'$Z_{11}$')
plt.yscale('log')
#plt.savefig('opacity_phasefunction.png')
plt.show()

