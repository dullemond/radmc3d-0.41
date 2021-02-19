from radmc3dPy.image import *
from plot import *
from natconst import *
import os
# Run python problem_setup.py
# Run radmc3d mctherm
os.system('radmc3d image lambda 0.5 nostar incl 60')
imfull = readImage()
os.system('radmc3d image lambda 0.5 nostar incl 60 maxnrscat 1')
imsing = readImage()
os.system('radmc3d image lambda 0.5 nostar incl 60 simplescat')
imsmpl = readImage()
maxlog   = 6.
floorint = imfull.image.max()/10.**maxlog
plt.figure()
plt.plot(imfull.x,np.log10(imfull.image[:,imfull.ny/2,0]+floorint),label='Full')
plt.plot(imfull.x,np.log10(imsing.image[:,imsing.ny/2,0]+floorint),label='Sing')
plt.plot(imfull.x,np.log10(imsmpl.image[:,imsmpl.ny/2,0]+floorint),label='Smpl')
plt.legend()
plt.figure()
plt.plot(imfull.x,np.log10(imfull.image[imfull.nx/2,:,0]+floorint),label='Full')
plt.plot(imfull.x,np.log10(imsing.image[imsing.nx/2,:,0]+floorint),label='Sing')
plt.plot(imfull.x,np.log10(imsmpl.image[imsmpl.nx/2,:,0]+floorint),label='Smpl')
plt.legend()
extent=np.array([-1,1,-1,1])*imfull.x.max()/au
plt.figure()
plt.imshow(np.log10(imfull.image[:,:,0].T+floorint),extent=extent,origin='lower')
plt.colorbar()
plt.xlabel('x [au]')
plt.ylabel('y [au]')
plt.title('Full Monte Carlo scattering')
plt.figure()
plt.imshow(np.log10(imsing.image[:,:,0].T+floorint),extent=extent,origin='lower')
plt.colorbar()
plt.xlabel('x [au]')
plt.ylabel('y [au]')
plt.title('Monte Carlo single scattering')
plt.figure()
plt.imshow(np.log10(imsmpl.image[:,:,0].T+floorint),extent=extent,origin='lower')
plt.colorbar()
plt.xlabel('x [au]')
plt.ylabel('y [au]')
plt.title('Simple single scattering')
plt.show()
