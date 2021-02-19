from radmc3dPy.image import *
from plot import *
from natconst import *
import os
# Run python problem_setup.py
# Run radmc3d mctherm
with open('radmc3d.inp','w') as f:
    f.write('nphot = 100000\n')
    f.write('scattering_mode_max = 1\n')
os.system('radmc3d image lambda 0.5 nostar incl 55 simplescat')
im_iso = readImage()
with open('radmc3d.inp','w') as f:
    f.write('nphot = 100000\n')
    f.write('scattering_mode_max = 5\n')
os.system('radmc3d image lambda 0.5 nostar incl 55 simplescat')
im_aniso = readImage()
maxlog   = 6.
floorint = im_aniso.image.max()/10.**maxlog
plt.figure()
plt.plot(im_aniso.x,np.log10(im_aniso.image[:,im_aniso.ny/2,0]+floorint),label='Aniso')
plt.plot(im_aniso.x,np.log10(im_iso.image[:,im_iso.ny/2,0]+floorint),label='Iso')
plt.legend()
plt.figure()
plt.plot(im_aniso.x,np.log10(im_aniso.image[im_iso.nx/2,:,0]+floorint),label='Aniso')
plt.plot(im_aniso.x,np.log10(im_iso.image[im_iso.nx/2,:,0]+floorint),label='Iso')
plt.legend()
extent=np.array([-1,1,-1,1])*im_aniso.x.max()/au
plt.figure()
plt.imshow(np.log10(im_aniso.image[:,:,0].T+floorint),extent=extent,origin='lower')
plt.colorbar()
plt.xlabel('x [au]')
plt.ylabel('y [au]')
plt.title('Anisotopic scattering')
plt.figure()
plt.imshow(np.log10(im_iso.image[:,:,0].T+floorint),extent=extent,origin='lower')
plt.colorbar()
plt.xlabel('x [au]')
plt.ylabel('y [au]')
plt.title('Isotropic scattering')
plt.show()
