import numpy as np

def raytrace_spher(ri,alphanu,jnu,b,intens0):
    """
    Supersimple raytracer for 1-D spherical coordinates. Simple first-order integration.

    ARGUMENTS:
      ri          The radii of the interfaces of the cells (array is 1 longer than nr of cells)
      alphanu     The alpha_nu extinction coefficient at the cell centers
      jnu         The j_nu source function at the cell centers
      b           The impact parameter of the ray
      itens0      The background intensity

    RETURNS:
      intens      The intensity emerging from the cloud

    EXAMPLE:
    from raytrace_spher import *
    import matplotlib.pyplot as plt
    rin     = 1.0
    rout    = 5.0
    nr      = 30
    ri      = rin * (rout/rin)**np.linspace(0.,1.,nr) # Cell walls
    r       = (ri[1:]*ri[:-1])**0.5                   # Cell centers
    alphanu = np.ones_like(r)*3e-1
    srcnu   = r**(-0.5)
    jnu     = srcnu*alphanu
    intens0 = 0.0
    nx      = 1000
    x       = np.linspace(-1.,1.,nx)*ri.max()
    image   = np.zeros(nx)
    for ix in range(nx):
        image[ix] = raytrace_spher(ri,alphanu,jnu,np.abs(x[ix]),intens0)
    plt.figure()
    plt.plot(x,image)
    plt.show()
    """
    assert b>=0, "Negative b is not allowed"
    assert ri[1]>ri[0], "ri must be monotonically increasing"
    scross = np.zeros(len(ri))
    ir0 = 0
    for ir in range(len(ri)):
        if ri[ir]>b:
            scross[ir] = np.sqrt(ri[ir]**2-b**2)
        else:
            ir0 = ir+1
    sprev  = -scross[-1]
    intens = intens0
    if b<ri[-1]:
        for ir in range(len(ri)-2,ir0-1,-1):
            scurr  = -scross[ir]
            ds     = scurr-sprev
            dtau   = alphanu[ir]*ds
            src    = jnu[ir]/alphanu[ir]
            intens = intens*np.exp(-dtau) + (1.0-np.exp(-dtau))*src
            sprev  = scurr
        if ir0>0:
            ir     = ir0-1
            scurr  = scross[ir0]
            ds     = scurr-sprev
            dtau   = alphanu[ir]*ds
            src    = jnu[ir]/alphanu[ir]
            intens = intens*np.exp(-dtau) + (1.0-np.exp(-dtau))*src
            sprev  = scurr
        else:
            sprev  = scross[0]
        for ir in range(ir0,len(ri)-1):
            scurr  = scross[ir+1]
            ds     = scurr-sprev
            dtau   = alphanu[ir]*ds
            src    = jnu[ir]/alphanu[ir]
            intens = intens*np.exp(-dtau) + (1.0-np.exp(-dtau))*src
            sprev  = scurr
    return intens
