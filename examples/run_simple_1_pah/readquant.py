import numpy as np

def readquant(ispec):
    filename = "tempdistrib_%d"%(ispec)+".dat"
    with open(filename) as f:
        str=f.readline()
        str=f.readline()
        ncells=int(str)
        str=f.readline()
        ntemp=int(str)
        temp=np.zeros(ntemp)
        str=f.readline()
        for it in range(ntemp):
            str=f.readline()
            temp[it]=float(str)
        str=f.readline()
        tdist=np.zeros((ncells,ntemp))
        for icell in range(ncells):
            str=f.readline().split()
            for it in range(ntemp):
                tdist[icell,it] = float(str[it])
    return {"temp":temp,"tdist":tdist}
