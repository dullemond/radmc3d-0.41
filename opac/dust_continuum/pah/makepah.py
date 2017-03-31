import numpy as np
import math

def pah(nc,nh,z,lam,lab=False,xgra=0.01):
    """
    PAH Model of Li & Draine (2001), ApJ, 554, 778

    nc   = Nr of C atoms of this PAH
    nh   = Nr of H atoms of this PAH
    z    = 0 or 1 is charge of this PAH
    lam  = wavelength grid in micron
    lab  = False = Observed fitted values, = True = lab values
    xgra = Mixing ratio of pure PAH and graphite
    """
    #
    # A constant
    #
    mp   = 1.6726e-24     # Mass of proton          [g]
    #
    # Ratio of N_H / N_C 
    #
    hc   = (1.0*nh)/(1.0*nc)
    #
    # The inverse wavelength
    #
    x = 1./lam
    #
    # Initialize the arrays
    #
    stot = x*0.
    spah = x*0.
    #
    # Which PAH feature values to use
    #
    if lab:
        e62 = 1.   # Lab measurement
        e77 = 1.   # Lab measurement
        e86 = 1.   # Lab measurement
    else:
        e62 = 3.   # Observed value
        e77 = 2.   # Observed value
        e86 = 2.   # Observed value
    #
    # The parameters of the Drude profiles (D&L2001 Table 1)
    #
    pah_lam    = np.array([0.0722, 0.2175,   3.3,   6.2,   7.7,    8.6, 
                           11.3,     11.9,  12.7,   16.4, 18.3,   21.2,
                           23.1,     26.0])
    pah_gamma  = np.array([0.195,  0.217,  0.012,  0.032, 0.091, 0.047,
                           0.018,  0.025,  0.024,  0.010, 0.036, 0.038,
                           0.046,   0.69])
    pah_gamlam = np.array([0.0141, 0.0473,  0.04,  0.2,   0.7,   0.4,
                           0.2,    0.3,     0.3,   0.16,  0.66,  0.81, 
                           1.07,   18.]) * 1e-4
    pah_signeu = np.array([7.97e7, 1.23e7, 197*hc, 19.6*e62, 60.9*e77, 34.7*e86*hc,
                           427*hc/3, 72.7*hc/3, 167*hc/3, 5.52, 6.04, 10.8,
                           2.78, 15.2]) * 1e-20
    pah_sigcha = np.array([7.97e7, 1.23e7, 44.7*hc, 157*e62, 548*e77, 242* e86*hc,
                           400*hc/3, 61.4*hc/3, 149*hc/3, 5.52, 6.04, 10.8,
                           2.78, 15.2]) * 1e-20
    #
    # Nr of benzene rings and cutoff wavelength (L&D2001 Appendix A)
    #
    if nc<40:
        Mbr = 0.3*nc
    else:
        Mbr = 0.4*nc
    if z>0:
        lam_c  = 1. / (2.282 * Mbr**(-0.5) + 0.889 )  # Ionized PAH
    else:
        lam_c  = 1. / (3.804 * Mbr**(-0.5) + 1.052 )  # Neutral PAH
    #
    # Cutoff formula
    #
    if z>0:
        y = lam_c / lam   # Ionized PAH
    else:      
        y = lam_c / lam   # Neutral PAH
    cutoff = (1./math.pi) * np.arctan(1000.*(y-1.)**3/y) + 0.5
    #
    # Which PAH feature parameters to take
    #
    if z>0:
        pah_sig = pah_sigcha
    else:
        pah_sig = pah_signeu
    #
    # The PAH features in the Infrared
    #
    for imode in range(2,pah_lam.size):
        spah += (2/math.pi)*pah_gamlam[imode]*pah_sig[imode]/((lam/pah_lam[imode]-pah_lam[imode]/lam)**2+pah_gamma[imode]**2)
    #
    # The Optical/UV features
    #
    imode = 0
    s0 = (2/math.pi)*pah_gamlam[imode]*pah_sig[imode]/((lam/pah_lam[imode]-pah_lam[imode]/lam)**2+pah_gamma[imode]**2)
    imode = 1
    s1 = (2/math.pi)*pah_gamlam[imode]*pah_sig[imode]/((lam/pah_lam[imode]-pah_lam[imode]/lam)**2+pah_gamma[imode]**2)
    #
    # We ignore carbon opacity for ultrashort wavelength
    #
    # Adding the various continuum components and the PAH features
    #
    i = (x>17.25)
    stot[i] = (126.-6.4943*17.25)*(x[i]/17.25)**(-1)*1e-18 # CPD's replacement of graphite
    i = (x>15.) & (x<=17.25)
    stot[i] = (126.-6.4943*x[i])*1e-18
    i = (x>10.) & (x<=15)
    stot[i] = (-3.+1.35*x[i])*1e-18 + s0[i]
    i = (x>7.7) & (x<=10)
    stot[i] = (66.302-24.367*x[i]+2.950*x[i]**2-0.1057*x[i]**3)*1e-18
    i = (x>5.9) & (x<=7.7)
    stot[i] = (1.8687+0.1905*x[i]+0.4175*(x[i]-5.9)**2-0.0437*(x[i]-5.9)**3)*1e-18 + s1[i]
    i = (x>3.3) & (x<=5.9)
    stot[i] = (1.8687+0.1905*x[i])*1e-18 + s1[i]
    i = (x<=3.3)
    stot[i] = cutoff[i]*34.58*10.**(-18.-3.431/x[i]) + spah[i]
    #
    # Now convert to opacity
    #
    kappa_pah = stot / ( 12 * mp )
    #
    # The model of L&D2001, Fig 2, shows a continuum in the NIR that
    # is not, so far, in the stot here. This is a 1% contribution 
    # by graphite. 
    #
    # Make a dummy graphite model (very simple!)
    #
    kappa_gra = 1e4*lam**(-1.6)+5e2/((33./lam)**2.1+(lam/33.)**2.1)
    kappa_gra[kappa_gra>1e5] = 1e5
    #
    # Mix these two opacities
    #
    kappa = (1.-xgra) * kappa_pah + xgra * kappa_gra
    #
    # Return this result
    #
    return kappa



def makepah(nc,nh,z,lam,lab=False,xgra=0.01):
    """
    Call the PAH routine and write the opacity file.

    nc   = Nr of C atoms of this PAH
    nh   = Nr of H atoms of this PAH
    z    = 0 or 1 is charge of this PAH
    lam  = wavelength grid in micron
    lab  = False = Observed fitted values, = True = lab values
    xgra = Mixing ratio of pure PAH and graphite
    """
    kappa = pah(nc,nh,z,lam,lab=lab,xgra=xgra)
    with open('dustkappa_pah.inp','w') as f:
        f.write('2\n')
        f.write('%4d\n'%(lam.size))
        for inu in range(lam.size):
            f.write('%13.6e %13.6e %13.6e\n'%(lam[inu],kappa[inu],0.))



def example():
    """
    Simple example of how to set up a wavelength grid and compute the
    PAH opacity for that grid.
    """
    mp   = 1.6726e-24     # Mass of proton          [g]
    nf   = 1000
    lam0 = 0.03
    lam1 = 1000.
    lam  = lam0 * (lam1/lam0)**(np.linspace(0.,1.,nf))
    nc   = 40.
    nh   = 16.
    z    = 0
    makepah(nc,nh,z,lam,lab=False,xgra=0.01)
