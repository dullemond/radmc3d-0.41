import math
import scipy.optimize as op
import numpy as np

def bplfunc(nu,t):
    cc  = 2.9979245800000e10
    hh  = 6.6262e-27
    kk  = 1.3807e-16
    bpl = 0.e0
    if t > 0.0:
        x = hh*nu/(kk*t)
        if x < 1.e2:
            if x > 1.e-3:
                bpl = 2.0*hh*nu*nu*nu/(cc*cc*(math.exp(x)-1.e0))
            else:
                bpl = 2.0*nu*nu*kk*t/(cc*cc)
    return bpl

def bplanck(nu,t):
    if type(nu)==np.ndarray:
        n = nu.size
        bnu = np.zeros(n)
        for i in range(n):
            bnu[i] = bplfunc(nu[i],t)
    else:
        bnu = bplfunc(nu,t)
    return bnu

def heatcool(temp,nu,kabs,jnu):
    bnu  = bplanck(nu,temp)
    f    = kabs*(jnu-bnu)
    dnu  = np.abs(nu[1:]-nu[0:-1])
    fav  = 0.5*(f[1:]+f[0:-1])
    tot  = (fav*dnu).sum()
    return tot

def solvetdust_bbstar(nu,kabs,rstar,tstar,dist):
    bpl  = bplanck(nu,tstar)
    lnu  = math.pi*bpl * 4*math.pi * rstar**2
    fnu  = lnu / ( 4*math.pi * dist**2 )
    jnu  = fnu / ( 4*math.pi )
    temp = op.brentq(heatcool,1e-2,1e4,args=(nu,kabs,jnu))
    return temp

def solvetdust_jnu(nu,kabs,jnu):
    temp = op.brentq(heatcool,1e-2,1e4,args=(nu,kabs,jnu))
    return temp

# def test():
#     cc   = 2.9979245800000e10      # Light speed             [cm/s]
#     nf   = 100
#     lam  = 10.0**np.linspace(-1,3,nf)
#     nu   = 1e4*cc/lam
#     kabs = 1e3/lam       # Powerlaw opacity kappa_abs propto 1/lambda
#     #kabs = 1e3 + 0*lam   # Grey opacity
#     rstar= 6.96e10
#     tstar= 5.78e3
#     dist = 1.49598e13
#     temp = solvetdust_bbstar(nu,kabs,rstar,tstar,dist)
#     print("%13.6e\n"%(temp))
