from radmc3dPy.image import *
from plot import *
from natconst import *
import os
# First:
# Run python problem_setup.py
# Run radmc3d mctherm
#
lammic = 0.5     # Wavelength in micron
incl   = 55      # Inclination 

lammicstr = '{}'.format(lammic)
inclstr   = '{}'.format(incl)

with open('radmc3d.inp','w') as f:
    f.write('nphot = 100000\n')
    f.write('scattering_mode_max = 1\n')
os.system('radmc3d image lambda '+lammicstr+' nostar incl '+inclstr+' simplescat')
im_iso = readImage()
with open('radmc3d.inp','w') as f:
    f.write('nphot = 100000\n')
    f.write('scattering_mode_max = 5\n')
os.system('radmc3d image lambda '+lammicstr+' nostar incl '+inclstr+' simplescat')
im_aniso = readImage()
maxlog   = 3.
floorint = im_aniso.image.max()/10.**maxlog
#image_aniso = np.log10(im_aniso.image[:,:,0]+floorint)
#image_iso   = np.log10(im_iso.image[:,:,0]+floorint)
image_aniso = im_aniso.image[:,:,0]
image_iso   = im_iso.image[:,:,0]
vmax = image_aniso.max()
vmin = image_aniso.min()
plt.figure()
plt.plot(im_aniso.x,image_aniso[:,im_aniso.ny/2],label='Aniso')
plt.plot(im_aniso.x,image_iso[:,im_iso.ny/2],label='Iso')
plt.legend()
plt.figure()
plt.plot(im_aniso.x,image_aniso[im_iso.nx/2,:],label='Aniso')
plt.plot(im_aniso.x,image_iso[im_iso.nx/2,:],label='Iso')
plt.legend()
extent=np.array([-1,1,-1,1])*im_aniso.x.max()/au
plt.figure()
plt.imshow(image_aniso.T,extent=extent,vmin=vmin,vmax=vmax,origin='lower')
plt.colorbar()
plt.xlabel('x [au]')
plt.ylabel('y [au]')
plt.title('Anisotopic scattering')
plt.figure()
plt.imshow(image_iso.T,extent=extent,vmin=vmin,vmax=vmax,origin='lower')
plt.colorbar()
plt.xlabel('x [au]')
plt.ylabel('y [au]')
plt.title('Isotropic scattering')

from radmc3dPy.analyze import *
o    = readOpac(ext=['mix'],scatmat=[True])
nlam = len(o.wav[0])
ilam = int(np.interp(lammic,o.wav[0],np.linspace(0,nlam-1,nlam)))
theta_scat_forward   = 90. - incl
theta_scat_backward  = 90. + incl
nang = o.nang[0]
itheta_scat_forward   = int(np.interp(theta_scat_forward,o.scatang[0],np.linspace(0,nang-1,nang)))
itheta_scat_backward  = int(np.interp(theta_scat_backward,o.scatang[0],np.linspace(0,nang-1,nang)))
ratio_from_opacitylaw = o.z11[0][ilam,itheta_scat_forward]/o.z11[0][ilam,itheta_scat_backward]
nxhalf = image_aniso.shape[0]/2
nyhalf = image_aniso.shape[1]/2
ratio_from_image      = image_aniso[nxhalf,:nyhalf].max()/image_aniso[nxhalf,nyhalf:].max()

print('Ratio from scattering phase function = {}'.format(ratio_from_opacitylaw))
print('Ratio from image                     = {}'.format(ratio_from_image))
print('  (they should be similar, but do not need to be identical, given geometry, interpolation etc.)')

plt.show()
