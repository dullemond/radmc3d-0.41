import numpy as np
def readopac(ext='1'):
    filename = 'dustkapscatmat_'+ext+'.inp'
    with open(filename,'r') as f:
        str = f.readline()
        str = f.readline()
        str = f.readline()
        str = f.readline()
        str = f.readline()
        str = f.readline()
        iformat = int(f.readline())
        nlam    = int(f.readline())
        nang    = int(f.readline())
        lam     = np.zeros(nlam)
        kabs    = np.zeros(nlam)
        ksca    = np.zeros(nlam)
        gsca    = np.zeros(nlam)
        ang     = np.zeros(nang)
        zmatrix = np.zeros((nlam,nang,6))
        str = f.readline()
        for inu in range(nlam):
            dum = f.readline().split()
            lam[inu]  = float(dum[0])
            kabs[inu] = float(dum[1])
            ksca[inu] = float(dum[2])
            gsca[inu] = float(dum[3])
        str = f.readline()
        for iang in range(nang):
            ang[iang] = float(f.readline())
        str = f.readline()
        for inu in range(nlam):
            for iang in range(nang):
                dum = f.readline().split()
                for iz in range(6):
                    zmatrix[inu,iang,iz] = float(dum[iz])
    q = {}
    q['lam'] = lam
    q['kabs'] = kabs
    q['ksca'] = ksca
    q['gsca'] = gsca
    q['ang'] = ang
    q['zmatrix'] = zmatrix
    return q

