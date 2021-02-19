from radmc3dPy.image import *
from plot import *
from bplanck import *
from raytrace_spher import *
from natconst import *
# Run python problem_setup.py
# Run radmc3d mctherm
#### Run radmc3d image lambda 0.5 nostar
# Run radmc3d image lambda 0.5 nostar simplescat
im = readImage()
with open('stars.inp','r') as f:
    f.readline()
    nlam = int(f.readline().split()[1])
    f.readline()
    rstar = float(f.readline().split()[0])
    f.readline()
    for i in range(nlam): f.readline()
    f.readline()
    tstar = -float(f.readline())
with open('amr_grid.inp','r') as f:
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    nr = int(f.readline().split()[0])
    ri = np.zeros(nr+1)
    for ir in range(nr+1):
        ri[ir] = float(f.readline())
with open('dustkappa_silicate.inp','r') as f:
    f.readline()
    nlamkap = int(f.readline().split()[0])
    data = np.loadtxt(f)
    lamkap = data[:,0]
    kapsca = data[:,2]
lammic     = 1e4*cc/im.freq[0]
kappa_scat = np.interp(lammic,lamkap,kapsca)         # Kappa_scat at lambda of image
with open('dust_density.inp','r') as f:
    f.readline()
    nrr = int(f.readline().split()[0])
    f.readline()
    assert nrr==nr,'Grid sizes do not match'
    rhodust = np.loadtxt(f)
tauscati = np.zeros_like(ri)
for ir in range(nr):
    tauscati[ir+1] = tauscati[ir] + (ri[ir+1]-ri[ir])*rhodust[ir]*kappa_scat
assert tauscati[-1]<1e-2, 'Sorry, only very optically thin test allowed (tauscat<1e-2)'
r       = (ri[1:]*ri[:-1])**0.5                   # Cell centers
bnustar = bplanck(im.freq[0],tstar)
lnustar = 4*np.pi*rstar**2 * np.pi * bnustar
fnustar = lnustar / (4*np.pi*r**2)
Jnustar = fnustar / (4*np.pi)
anuscat = rhodust*kappa_scat
snuscat = Jnustar
jnuscat = anuscat*snuscat
x       = im.x
image   = np.zeros(len(x))
for ix in range(len(x)):
    image[ix] = raytrace_spher(ri,anuscat,jnuscat,np.abs(x[ix]),0.0)

maxlog = 6.
floorint = im.image.max()/10.**maxlog

plt.figure()
plt.plot(np.log10(im.image[:,im.ny/2,0]+floorint),label='RADMC-3D')
plt.plot(np.log10(image+floorint),'o',label='Python')
plt.legend()
plt.show()

