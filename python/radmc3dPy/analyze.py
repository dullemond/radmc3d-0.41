"""This module contains classes and functions to read and write input/output data for RADMC-3D and
to do some simple analysis/diagnostics of the model.

"""
try:
    import numpy as np
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'


try:
    import matplotlib.pylab as plb
except:
    print ' WARNING'
    print ' matploblib.pylab cannot be imported ' 
    print ' To used the visualization functionality of the python module of RADMC-3D you need to install matplotlib'
    print ' Without matplotlib you can use the python module to set up a model but you will not be able to plot things or'
    print ' display images'

import subprocess as sp
import sys, os, copy
from radmc3dPy.natconst import *
import radmc3dPy.crd_trans  as crd_trans 
from staratm import StellarAtm


class radmc3dGrid(object):
    """ Class for spatial and frequency grid used by RADMC-3D.

    Attributes
    ----------

    act_dim    : ndarray 
                A three element vector the i-th element is 1 if the i-th dimension is active, otherwize the i-th element is zero
    
    crd_sys    : {'sph', 'cyl', 'car'} 
                coordinate system of the spatial grid
    
    nx         : int
                Number of grid points in the x (cartesian) / r (cylindrical) / r (spherical) dimension
    
    ny         : int
                Number of grid points in the y (cartesian) / theta (cylindrical) / theta (spherical) dimension
    
    nz         : int
                Number of grid points in the z (cartesian) / z (cylindrical) / phi (spherical) dimension
    
    nxi        : int
                Number of cell interfaces in the x (cartesian) / r (cylindrical) / r (spherical) dimension
    
    nyi        : int
                Number of cell interfaces in the y (cartesian) / theta (cylindrical) / theta (spherical) dimension
    
    nzi        : int
                Number of cell interfaces in the z (cartesian) / z (cylindrical) / phi (spherical) dimension
    
    nwav       : int
                Number of wavelengths in the wavelength grid
    
    nfreq      : int
                Number of frequencies in the grid (equal to nwav)
    
    x          : ndarray
                Cell centered x (cartesian) / r (cylindrical) / r (spherical)  grid points
    
    y          : ndarray
                Cell centered y (cartesian) / theta (cylindrical) / theta (spherical)  grid points
    
    z          : ndarray
                Cell centered z (cartesian) / z (cylindrical) / phi (spherical)  grid points
    
    xi         : ndarray
                Cell interfaces in the x (cartesian) / r (cylindrical) / r (spherical)  dimension
    
    yi         : ndarray
                Cell interfaces in the y (cartesian) / theta (cylindrical) / theta (spherical)  dimension
    
    zi         : ndarray
                Cell interfaces in the z (cartesian) / z (cylindrical) / phi (spherical)  dimension
    
    wav        : ndarray
                Wavelengh  grid
    
    freq       : ndarray
                Frequency  grid


    """
# --------------------------------------------------------------------------------------------------

    def __init__(self):

        self.crd_sys = 'sph'
        self.act_dim = [1,1,1]
        self.nx    = -1
        self.ny    = -1
        self.nz    = -1
        self.nxi   = -1
        self.nyi   = -1
        self.nzi   = -1
        self.nwav  = -1
        self.nfreq = -1
        self.x     = -1
        self.y     = -1
        self.z     = -1
        self.xi    = -1
        self.yi    = -1
        self.zi    = -1
        self.xx    = -1
        self.yy    = -1
        self.zz    = -1
        self.wav   = -1
        self.freq  = -1

# --------------------------------------------------------------------------------------------------
    def makeWavelengthGrid(self, wbound=None, nw=None, ppar=None):
        """Creates the wavelength/frequency grid.

        Parameters
        ----------

        wbound : list 
                 Contains the wavelength boundaries of the wavelength grid (should contain at least two elements)

        nw     : list 
                 Contains len(wbound)-1 elements containing the number of wavelengths between the bounds
                 set by wbound
        
        ppar   : dictionary, optional
                 Contains all input parameters with the parameter names as keys  
        """
        
        if ppar:
            if not wbound: wbound = ppar['wbound']
            if not nw: nw = ppar['nw']

        if (wbound==None)|(nw==None):
            if (ppar==None): 
                print 'ERROR!'
                print 'Either the boundaries or the number of gridpoints has not be specified in the wavelength grid'
                return
            
        self.nwav = nw[0]
        self.wav  = wbound[0] * (wbound[1]/wbound[0])**(np.arange(nw[0], dtype=np.float64) / nw[0])

        for ipart in range(1,len(nw)-1): 
            dum      = wbound[ipart] * (wbound[ipart+1]/wbound[ipart])**(np.arange(nw[ipart], dtype=np.float64) / nw[ipart])
            self.wav = np.append(self.wav, dum)

        ipart      = len(nw)-1
        dum        = wbound[ipart] * (wbound[ipart+1]/wbound[ipart])**(np.arange(nw[ipart], dtype=np.float64) / (nw[ipart]-1.))
        self.wav   = np.append(self.wav, dum)
        self.nwav  = self.wav.shape[0]
        cc = 29979245800.
        self.freq  = cc / self.wav * 1e4
        self.nfreq = self.nwav

# --------------------------------------------------------------------------------------------------
    def writeWavelengthGrid(self, fname='', old=False):
        """Wriites the wavelength grid to a file (e.g. wavelength_micron.inp).

        Parameters
        ----------
        
        fname  : str, optional
                 File name into which the wavelength grid should be written. If omitted 'wavelength_micron.inp' will be used
        
        old    : bool, optional
                 If set to True the file format of the previous, 2D version of radmc will be used
        """
       
        if not old:
            if fname=='':
                fname = 'wavelength_micron.inp'

            print 'Writing '+fname
            wfile = open(fname, 'w')
            wfile.write('%d\n'%self.nwav)
            for ilam in range(self.nwav):
                wfile.write('%.9e\n'%self.wav[ilam])
            wfile.close()
        else:
            if fname=='':
                fname='frequency.inp'
            try :
                wfile = open(fname, 'w')
            except:
                print 'Error!' 
                print fname+' cannot be opened!'
                return 
           
            print 'Writing '+fname
            wfile.write("%d\n"%self.nfreq)
            wfile.write(" \n")
            #
            # Reverse the order of the frequency grid as it is ordered in frequency in radmc
            #
            freq = self.freq[::-1]
            for i in range(self.nfreq):
                wfile.write("%.7e\n"%freq[i])
            
            wfile.close()             
# --------------------------------------------------------------------------------------------------
    def makeSpatialGrid(self,crd_sys=None,xbound=None,ybound=None,zbound=None,nxi=None,nyi=None,nzi=None,ppar=None):
        """Calculates the spatial grid.

        Parameters
        ----------
        
        crd_sys : {'sph','car'}  
                    Coordinate system of the spatial grid
        
        xbound  : list 
                    (with at least two elements) of boundaries for the grid along the first dimension
        
        ybound  : list 
                    (with at least two elements) of boundaries for the grid along the second dimension
        
        zbound  : list 
                    (with at least two elements) of boundaries for the grid along the third dimension
        
        nxi     : int
                    Number of grid points along the first dimension. List with len(xbound)-1 elements with 
                    nxi[i] being the number of grid points between xbound[i] and xbound[i+1]
        
        nyi     : int
                    Same as nxi but for the second dimension
        
        nzi     : int
                    Same as nxi but for the third dimension
        
        ppar    : Dictionary containing all input parameters of the model (from the problem_params.inp file)
                   if ppar is set all keyword arguments that are not set will be taken from this dictionary
        """

        self.act_dim = [1,1,1]
        if ppar:
            if not crd_sys : crd_sys = ppar['crd_sys']
            self.crd_sys =  crd_sys
           
            if not xbound : 
                if ppar.has_key('xbound'):
                    xbound = ppar['xbound']
                else:
                    print ' No boundary for the first dimension is given, first dimension is deactivated.'
                    self.act_dim[0] = 0
            if not nxi:
                if ppar.has_key('nx'):
                    if (type(ppar['nx']).__name__!='list'): 
                        ppar['nx'] = [ppar['nx']]
                    nxi = [i+1 for i in ppar['nx']] #nxi = ppar['nx']+1
                    if ppar['nx'][0]==0:
                        self.act_dim[0] = 0
                else:
                    self.act_dim[0] = 0


            if not ybound : 
                if ppar.has_key('ybound'):
                    ybound = ppar['ybound']
                else:
                    print ' No boundary for the second dimension is given, second dimension is deactivated.'
                    self.act_dim[1] = 0
            if not nyi:
                if ppar.has_key('ny'):
                    if (type(ppar['ny']).__name__!='list'): 
                        nyi = [ppar['ny']+1]
                        ppar['ny'] = [ppar['ny']]

                    else:
                        ppar['ny'] = ppar['ny']
                        nyi = [i+1 for i in ppar['ny']] #ppar['ny']+1
                    
                    if ppar['ny'][0]==0:
                        self.act_dim[1] = 0
                else:
                    self.act_dim[1] = 0

            if not zbound : 
                if ppar.has_key('zbound'):
                    zbound = ppar['zbound']
                else:
                    print ' No boundary for the third dimension is given, third dimension is deactivated.'
                    self.act_dim[2] = 0
            if not nzi:
                if (ppar.has_key('nz'))&(ppar['nz']>0.):
                    if (type(ppar['nz']).__name__!='list'): 
                        ppar['nz'] = [ppar['nz']]
                    nzi = [i+1 for i in ppar['nz']] #nzi = ppar['nz']+1
                    if ppar['nz'][0]==0:
                        self.act_dim[2] = 0
                else:
                    self.act_dim[2] = 0
                    nzi = [0]


        if (crd_sys=='car'):
#
# First check whether the grid boundaries are specified
#
            if (xbound==None): 
                print 'ERROR'
                print 'Boundaries on the cartesian x-axis is not specified'
                print 'Without the boundaries no grid can be created'
                return
            
            if (ybound==None): 
                print 'ERROR'
                print 'Boundaries on the cartesian y-axis is not specified'
                print 'Without the boundaries no grid can be created'
                return
            if (zbound==None): 
                print 'ERROR'
                print 'Boundaries on the cartesian z-axis is not specified'
                print 'Without the boundaries no grid can be created'
                return
            
            if ((nxi==None)|(nyi==None)|(nzi==None)):
                print 'ERROR'
                print 'Number of grid points is not specified'
                return

#
# Type checking 
#

            if (type(nxi).__name__=='int'):  nxi = [nxi]
            if (type(nyi).__name__=='int'):  nyi = [nyi]
            if (type(nzi).__name__=='int'):  nzi = [nzi]

#
# Create the x-axis
#
            if (len(nxi)>1): 
                self.nxi = sum(nxi)
                self.nx  = self.nxi-1
                self.xi  = xbound[0] + (xbound[1] - xbound[0])*(np.arange(nxi[0], dtype=np.float64)/float(nxi[0]))
                for ipart in range(1,len(nxi)-1):
                    dum = xbound[ipart] + (xbound[ipart+1] - xbound[ipart])*(np.arange(nxi[ipart], dtype=np.float64)/float(nxi[ipart]))
                    self.xi = append(self.xi, dum)

                ipart = len(nxi)-1 
                dum = xbound[ipart] + (xbound[ipart+1] - xbound[ipart])*(np.arange(nxi[ipart], dtype=np.float64)/float(nxi[ipart]-1))
                self.xi = append(self.xi, dum)
                self.x  = 0.5*(self.xi[0:self.nx] + self.xi[1:self.nx+1])
            else:
                if self.act_dim[0]==1:
                    self.nxi = nxi[0]
                    self.xi = xbound[0] + (xbound[1] - xbound[0])*(np.arange(self.nxi, dtype=np.float64)/float(self.nxi-1.))
                    self.nx = self.nxi-1
                    self.x  = 0.5*(self.xi[0:self.nx] + self.xi[1:self.nx+1])
                else:
                    self.x = [0.]
                    self.xi = [0., 0.,]
                    self.nx = 1
                    self.nxi = 2

#
# Create the y-ayis
#
            if (len(nyi)>1): 
                self.nyi = sum(nyi)
                self.ny  = self.nyi-1
                self.yi  = ybound[0] + (ybound[1] - ybound[0])*(np.arange(nyi[0], dtype=np.float64)/float(nyi[0]))
                for ipart in range(1,len(nyi)-1):
                    dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*(np.arange(nyi[ipart], dtype=np.float64)/float(nyi[ipart]))
                    self.yi = append(self.yi, dum)

                ipart = len(nyi)-1 
                dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*(np.arange(nyi[ipart], dtype=np.float64)/float(nyi[ipart]-1))
                self.yi = append(self.yi, dum)
                self.y  = 0.5*(self.yi[0:self.ny] + self.yi[1:self.ny+1])
            else:
                if self.act_dim[0]==1:
                    self.nyi = nyi[0]
                    self.yi = ybound[0] + (ybound[1] - ybound[0])*(np.arange(self.nyi, dtype=np.float64)/float(self.nyi-1.))
                    self.ny = self.nyi-1
                    self.y  = 0.5*(self.yi[0:self.ny] + self.yi[1:self.ny+1])
                else:
                    self.y = [0.]
                    self.yi = [0., 0.,]
                    self.ny = 1
                    self.nyi = 2


#
# Create the z-azis
#
            if (len(nzi)>1): 
                self.nzi = sum(nzi)
                self.nz  = self.nzi-1
                self.zi  = zbound[0] + (zbound[1] - zbound[0])*(np.arange(nzi[0], dtype=np.float64)/float(nzi[0]))
                for ipart in range(1,len(nzi)-1):
                    dum = zbound[ipart] + (zbound[ipart+1] - zbound[ipart])*(np.arange(nzi[ipart], dtype=np.float64)/float(nzi[ipart]))
                    self.zi = append(self.zi, dum)

                ipart = len(nzi)-1 
                dum = zbound[ipart] + (zbound[ipart+1] - zbound[ipart])*(np.arange(nzi[ipart], dtype=np.float64)/float(nzi[ipart]-1))
                self.zi = append(self.zi, dum)
                self.z  = 0.5*(self.zi[0:self.nz] + self.zi[1:self.nz+1])
            else:
                if self.act_dim[0]==1:
                    self.nzi = nzi[0]
                    self.zi = zbound[0] + (zbound[1] - zbound[0])*(np.arange(self.nzi, dtype=np.float64)/float(self.nzi-1.))
                    self.nz = self.nzi-1
                    self.z  = 0.5*(self.zi[0:self.nz] + self.zi[1:self.nz+1])
                else:
                    self.z = [0.]
                    self.zi = [0., 0.0]
                    self.nz = 1
                    self.nzi = 2



        if (crd_sys=='sph'): 
#
# r->x, theta->y, phi-z            
#
            if (xbound==None): 
                print 'ERROR'
                print 'Boundaries on the radius is not specified'
                print 'Without the boundaries no grid can be created'
                return

            if (ybound==None): ybound = [0.0, np.pi]
            if (zbound==None): zbound = [0.0, 2.0*np.pi]

            if ((nxi==None)|(nyi==None)|(nzi==None)):
                print 'ERROR'
                print 'Number of grid points is not specified'
                return

#
# Type checking (what is in the dimension numbers)
#

            if (type(nxi).__name__=='int'):  nxi = [nxi]
            if (type(nyi).__name__=='int'):  nyi = [nyi]
            if (type(nzi).__name__=='int'):  nzi = [nzi]
#
# Create the x axis
#
            if (len(nxi)>1): 
                self.nxi = sum(nxi)
                self.nx  = self.nxi-1
                self.xi  = xbound[0] * (xbound[1] / xbound[0])**(np.arange(nxi[0], dtype=np.float64)/float(nxi[0]))
                for ipart in range(1,len(nxi)-1):
                    dum = xbound[ipart] * (xbound[ipart+1] / xbound[ipart])**(np.arange(nxi[ipart], dtype=np.float64)/float(nxi[ipart]))
                    self.xi = np.append(self.xi, dum)

                ipart = len(nxi)-1 
                dum = xbound[ipart] * (xbound[ipart+1] / xbound[ipart])**(np.arange(nxi[ipart], dtype=np.float64)/float(nxi[ipart]-1))
                self.xi = np.append(self.xi, dum)
                self.x  = np.sqrt(self.xi[0:self.nx] * self.xi[1:self.nx+1])
            else:
                if self.act_dim[0]==1:
                    self.nxi = nxi[0]
                    self.xi = xbound[0] * (xbound[1] / xbound[0])**(np.arange(self.nxi, dtype=np.float64)/float(self.nxi-1.))
                    self.nx = self.nxi-1
                    self.x  = np.sqrt(self.xi[0:self.nx] * self.xi[1:self.nx+1])
                else:
                    self.x = [0.]
                    self.xi = [0., 0.,]
                    self.nx = 1
                    self.nxi = 2
                
            ## This has to be done properly
            #if ppar.has_key('xres_nlev'):
                #ri_ext = array([self.xi[0], self.xi[ppar['xres_nspan']]])
                #for i in range(ppar['xres_nlev']):
                    #dum_ri = ri_ext[0] + (ri_ext[1]-ri_ext[0]) * np.arange(ppar['xres_nstep']+1, dtype=np.float64) / float(ppar['xres_nstep'])
                    #print ri_ext[0:2]/au
                    #print dum_ri/au
                    #ri_ext_old = array(ri_ext)
                    #ri_ext = array(dum_ri)
                    #ri_ext = append(ri_ext,ri_ext_old[2:])
                    #print ri_ext/au
                    #print '----------'
                    
                #r_ext = (ri_ext[1:] + ri_ext[:-1]) * 0.5

                #self.xi = append(ri_ext, self.xi[ppar['xres_nspan']+1:])
                #self.x = append(r_ext, self.x[ppar['xres_nspan']:])
                #self.nx = self.x.shape[0]
                #self.nxi = self.xi.shape[0]

            # Refinement of the inner edge of the grid
            # This has to be done properly
            if ppar.has_key('xres_nlev'):
                if ppar['xres_nlev']>0:
                    ri_ext = np.array([self.xi[0], self.xi[ppar['xres_nspan']]])
                    for i in range(ppar['xres_nlev']):
                        dum_ri = ri_ext[0] + (ri_ext[1]-ri_ext[0]) * np.arange(ppar['xres_nstep']+1, dtype=np.float64) / float(ppar['xres_nstep'])
                        #print ri_ext[0:2]/au
                        #print dum_ri/au
                        ri_ext_old = np.array(ri_ext)
                        ri_ext = np.array(dum_ri)
                        ri_ext = np.append(ri_ext,ri_ext_old[2:])
                        
                    r_ext = (ri_ext[1:] + ri_ext[:-1]) * 0.5

                    self.xi = np.append(ri_ext, self.xi[ppar['xres_nspan']+1:])
                    self.x = np.append(r_ext, self.x[ppar['xres_nspan']:])
                    self.nx = self.x.shape[0]
                    self.nxi = self.xi.shape[0]
                
#
# Create the y axis
#
            if (len(nyi)>1):
                
                # Check if we go to the full [0,pi] interval or only use the upper half-plane [0, pi/2]
                
                if ybound[len(ybound)-1]!=np.pi/2.:
                    self.nyi = sum(nyi)+1
                    self.ny  = self.nyi-1
                    self.yi  = ybound[0] + (ybound[1] - ybound[0])*(np.arange(nyi[0], dtype=np.float64)/float(nyi[0])) 

                    for ipart in range(1,len(nyi)-1):
                        # Now make sure that pi/2 will be a cell interface
                        # 
                        # BUGFIX! 16-05-2012
                        # The grid was not symmetric to pi/2 when the grid contained multiple sections (i.e. len(nyi)>1)
                        # This is now fixed
                        if (ybound[ipart]<np.pi/2.):
                            dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*(np.arange(nyi[ipart], dtype=np.float64)/float(nyi[ipart]))
                        else:
                            if (ybound[ipart]==np.pi/2.):
                                dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((np.arange(nyi[ipart]+1, dtype=np.float64))/(float(nyi[ipart])))
                            else:
                                dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((np.arange(nyi[ipart], dtype=np.float64)+1.)/float(nyi[ipart]))

                        self.yi = np.append(self.yi, dum)

                    ipart   = len(nyi)-1 
                    if len(nyi)==2:
                        dum     = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((np.arange(nyi[ipart]+1, dtype=np.float64))/(float(nyi[ipart])))
                    else:
                        dum     = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((np.arange(nyi[ipart], dtype=np.float64)+1.)/float(nyi[ipart]))

                else:
                    self.nyi = sum(nyi)+1
                    self.ny  = self.nyi-1
                    self.yi  = ybound[0] + (ybound[1] - ybound[0])*(np.arange(nyi[0], dtype=np.float64)/float(nyi[0]))                
                    for ipart in range(1,len(nyi)-1):
                        # Now make sure that pi/2 will be a cell interface
                        # 
                        # BUGFIX! 16-05-2012
                        # The grid was not symmetric to pi/2 when the grid contained multiple sections (i.e. len(nyi)>1)
                        # This is now fixed
                        if (ybound[ipart]<np.pi/2.):
                            dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*(np.arange(nyi[ipart], dtype=np.float64)/float(nyi[ipart]))
                        else:
                            dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((np.arange(nyi[ipart]+1, dtype=np.float64))/(float(nyi[ipart])))
                        
                        self.yi = np.append(self.yi, dum)

                    ipart   = len(nyi)-1 

                    if len(nyi)==2:
                        dum     = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((np.arange(nyi[ipart]+1, dtype=np.float64))/(float(nyi[ipart])))
                    else:
                        dum     = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((np.arange(nyi[ipart], dtype=np.float64)+1.)/float(nyi[ipart]))



        
                
                self.yi = np.append(self.yi, dum)
                self.y  = 0.5*(self.yi[0:self.ny] + self.yi[1:self.ny+1])

            else:
                if self.act_dim[1]==1:
                    self.nyi = nyi[0]
                    self.yi = ybound[0] + (ybound[1] - ybound[0])*(np.arange(self.nyi, dtype=np.float64)/float(self.nyi-1.))
                    self.ny = self.nyi-1
                    self.y  = 0.5*(self.yi[0:self.ny] + self.yi[1:self.ny+1])
                else:
                    self.y = [0.]
                    self.yi = [0., 0.,]
                    self.ny = 1
                    self.nyi = 2
#
# Create the z axis

            if (len(nzi)>1):
                self.nzi = sum(nzi)
                self.nz  = self.nzi-1

                self.zi  = zbound[0] + (zbound[1] - zbound[0])*(np.arange(nzi[0], dtype=np.float64)/float(nzi[0]))                
                for ipart in range(1,len(nzi)-1):
                    dum = zbound[ipart] + (zbound[ipart+1] - zbound[ipart])*(np.arange(nzi[ipart], dtype=np.float64)/float(nzi[ipart]))
                    self.zi = np.append(self.zi, dum)
                ipart   = len(nzi)-1 
                dum     = zbound[ipart] + (zbound[ipart+1] - zbound[ipart])*(np.arange(nzi[ipart], dtype=np.float64)/float(nzi[ipart]-1))
                self.zi = np.append(self.zi, dum)
                self.z  = 0.5*(self.zi[0:self.nz] + self.zi[1:self.nz+1])
            else:
                if self.act_dim[2]==1:
                    self.nzi = nzi[0]
                    self.zi = zbound[0] + (zbound[1] - zbound[0])*(np.arange(self.nzi, dtype=np.float64)/float(self.nzi-1))
                    self.nz = self.nzi-1
                    self.z  = 0.5*(self.zi[0:self.nz] + self.zi[1:self.nz+1])
                else:
                    self.z = np.array([0.])
                    self.zi = np.array([0., np.pi*2.])
                    self.nz = 1
                    self.nzi = 2
            
        #if (crd_sys!='sph'):
            #print 'WARNING:'
            #print 'Currently only spherical coordinate system is supported!'
            #return

#       Create Multi-D versions of the coordinate grids (useful for
#       model making or plotting)

        qq       = np.meshgrid(self.x,self.y,self.z,indexing='ij')
        self.xx  = qq[0]
        self.yy  = qq[1]
        self.zz  = qq[2]




# --------------------------------------------------------------------------------------------------
    def writeSpatialGrid(self, fname='', old=False):
        """Writes the wavelength grid to a file (e.g. amr_grid.inp).

        Parameters
        ----------
        
        fname : str, optional
                File name into which the spatial grid should be written. If omitted 'amr_grid.inp' will be used. 
        
        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
        
        #
        # Write the spatial grid for radmc3d
        #
        if not old:
            if fname=='':
                fname = 'amr_grid.inp'

            print 'Writing '+fname
            wfile = open(fname, 'w')
            wfile.write('%d\n'%1)                    # Format number
            wfile.write('%d\n'%0)                    # AMR self.style (0=regular self. NO AMR)
            if self.crd_sys=='car':
                wfile.write('%d\n'%0)                  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if self.crd_sys=='sph':
                wfile.write('%d\n'%100)                  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            if self.crd_sys=='cyl':
                wfile.write('%d\n'%200)                  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
            wfile.write('%d\n'%0)                    # Gridinfo
            
            wfile.write('%d %d %d \n'%(self.act_dim[0], self.act_dim[1], self.act_dim[2]))       # Which dimension is active
            wfile.write('%d %d %d \n'%(self.nx,self.ny,self.nz))    # Grid size (x,y,z or r,phi,theta, or r,phi,z)
            for i in range(self.nxi): wfile.write('%.9e\n'%self.xi[i])
            for i in range(self.nyi): wfile.write('%.9e\n'%self.yi[i])
            for i in range(self.nzi): wfile.write('%.9e\n'%self.zi[i])
            wfile.close()
        #
        # Write the spatial grid for radmc
        #
        else:

            fname='radius.inp'
            try :
                wfile = open(fname, 'w')
            except:
                print 'Error!' 
                print fname+' cannot be opened!'
                return

            print 'Writing '+fname
            x = np.sqrt(self.xi[1:] * self.xi[:-1])
            wfile.write("%d\n"%self.nx)
            wfile.write(" \n")
            for i in range(self.nx):
                wfile.write("%.7e\n"%x[i])
            wfile.close()
            
            fname='theta.inp'
            try :
                wfile = open(fname, 'w')
            except:
                print 'Error!' 
                print fname+' cannot be opened!'
                return

            print 'Writing '+fname
            wfile.write("%d 1\n"%(self.ny/2))
            wfile.write(" \n")
            for i in range(self.ny/2):
                wfile.write("%.7e\n"%self.y[i])
            wfile.close()
            

# --------------------------------------------------------------------------------------------------
    def readGrid(self, fname='', old=False):
        """Reads the spatial (amr_grid.inp) and frequency grid (wavelength_micron.inp).
        
        Parameters
        ----------

        fname : str, optional
                File name from which the spatial grid should be read. If omitted 'amr_grid.inp' will be used. 

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """


#
# Natural constants
#

        cc = 29979245800.

# 
# Read the radmc3d format
#
        if not old:
            if fname=='':
                fname = 'amr_grid.inp'
    # 
    # Read the spatial grid 
    #
            try :
                rfile = open(fname, 'r')
            except:
                print 'Error!' 
                print 'amr_grid.inp was not found!'
                return 
        
            form        = float(rfile.readline())
            grid_style  = float(rfile.readline())
            crd_system  = int(rfile.readline())
            if crd_system<100:
                self.crd_sys = 'car'
            elif ((crd_system>=100)&(crd_system<200)):
                self.crd_sys = 'sph'
            elif ((crd_system>=200)&(crd_system<300)):
                self.crd_sys = 'cyl'
            else:
                rfile.close()
                print 'ERROR'
                print ' unsupported coordinate system in the amr_grid.inp file'
                print crd_system
                return

            grid_info   = float(rfile.readline())
            dum         = rfile.readline().split()
            self.act_dim = [int(dum[i]) for i in range(len(dum))]
            dum         = rfile.readline().split()
            self.nx,self.ny,self.nz    = int(dum[0]), int(dum[1]), int(dum[2])
            self.nxi,self.nyi,self.nzi = self.nx+1, self.ny+1, self.nz+1

            self.xi           = np.zeros(self.nx+1, dtype=np.float64)
            self.yi           = np.zeros(self.ny+1, dtype=np.float64)
            self.zi           = np.zeros(self.nz+1, dtype=np.float64)
           
            for i in range(self.nxi): self.xi[i] = float(rfile.readline())
            for i in range(self.nyi): self.yi[i] = float(rfile.readline())
            for i in range(self.nzi): self.zi[i] = float(rfile.readline())

            if self.crd_sys=='car':
                self.x = (self.xi[0:self.nx] +  self.xi[1:self.nx+1]) * 0.5
                self.y = (self.yi[0:self.ny] +  self.yi[1:self.ny+1]) * 0.5
                self.z = (self.zi[0:self.nz] +  self.zi[1:self.nz+1]) * 0.5
            else: 
                self.x = np.sqrt(self.xi[0:self.nx] * self.xi[1:self.nx+1])
                self.y = (self.yi[0:self.ny] +  self.yi[1:self.ny+1]) * 0.5
                self.z = (self.zi[0:self.nz] +  self.zi[1:self.nz+1]) * 0.5

            rfile.close()

            # Create Multi-D versions of the coordinate grids (useful for
            # model making or plotting)

            qq       = np.meshgrid(self.x,self.y,self.z,indexing='ij')
            self.xx  = qq[0]
            self.yy  = qq[1]
            self.zz  = qq[2]

    # 
    # Read the frequency grid 
    #

            try :
                rfile = open('wavelength_micron.inp', 'r')
            except:
                print 'Error!' 
                print 'wavelength_micron.inp was not found!'
                return 

            self.nwav = int(rfile.readline())
            self.nfreq = self.nwav
            self.wav  = np.zeros(self.nwav, dtype=np.float64)

            for i in range(self.nwav): self.wav[i] = float(rfile.readline())

            self.freq = cc / self.wav * 1e4

            rfile.close()
#
# Read the old radmc format
#
        else:
            self.crd_sys = 'sph'
            self.act_dim = [1,1,0]
        
            #
            # Read the radial grid
            #
            try: 
                rfile = open('radius.inp')
            except:
                print 'Error!' 
                print 'radius.inp was not found!'
                return 

            self.nx  = int(rfile.readline())
            self.nxi = self.nx + 1
            dum = rfile.readline()
            self.x  = np.zeros(self.nx, dtype=float)
            self.xi = np.zeros(self.nxi, dtype=float)
            for i in range(self.nx):
                self.x[i] = float(rfile.readline())
            self.xi[1:-1] = 0.5 * (self.x[1:] + self.x[:-1])
            self.xi[0]    = self.x[0] - (self.xi[1] - self.x[0])
            self.xi[-1]   = self.x[-1] + (self.x[-1] - self.xi[-2])
            rfile.close()

            
            #
            # Read the poloidal angular grid
            #
            try: 
                rfile = open('theta.inp')
            except:
                print 'Error!' 
                print 'theta.inp was not found!'
                return 

    
            ny = int(rfile.readline().split()[0])
            self.ny = ny*2
            self.nyi = self.ny+1
            dum = rfile.readline()

            self.y  = np.zeros(self.ny, dtype=float)
            self.yi = np.zeros(self.nyi, dtype=float)
            self.yi[0] = 0.
            self.yi[-1] = np.pi
            self.yi[ny] = np.pi*0.5

            for i in range(self.ny/2):
                 self.y[i]           = float(rfile.readline())
                 self.y[self.ny-1-i] = np.pi-self.y[i]

            self.yi[1:-1] = 0.5 * (self.y[1:] + self.y[:-1])
            self.yi[ny] = np.pi*0.5
            rfile.close()

            #
            # Create the azimuthal grid
            #

            self.nz = 1
            self.zi = np.array([0., 2.*np.pi], dtype=float)

            # 
            # Read the frequency grid 
            #
            try: 
                rfile = open('frequency.inp')
            except:
                print 'Error!' 
                print 'frequency.inp was not found!'
                return 

            self.nfreq = int(rfile.readline())
            self.nwav  = self.nfreq
            dum = rfile.readline()
            self.freq = np.zeros(self.nfreq, dtype=float)
            self.wav  = np.zeros(self.nfreq, dtype=float)
            for i in range(self.nfreq):
                self.freq[i] = float(rfile.readline())
                self.wav[i]  = cc/self.freq[i]*1e4
            rfile.close()
# --------------------------------------------------------------------------------------------------
    def getCellVolume(self):
        """Calculates the volume of grid cells.

        """

        if self.crd_sys=='sph':

            if self.act_dim[0]==0:
                print '----------------------------------------------------------'
                print 'ERROR'
                print 'The r-dimension of a spherical grid is switched off'
                print '----------------------------------------------------------'
            elif self.act_dim[1]==0:
                #print '----------------------------------------------------------'
                #print 'ERROR'
                #print 'The theta-dimension of a spherical grid is switched off'
                #print 'This model (ppdisk) is not perpared for such grid style'
                #print '----------------------------------------------------------'

                if self.act_dim[2]==0:
                    vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                    diff_r3   = self.xi[1:]**3 - self.xi[:-1]**3
                    diff_cost = 2.0
                    diff_phi  = 2.*np.pi
                    for ix in range(self.nx):
                        vol[ix,0,0] = 1./3. * diff_r3[ix] * diff_cost * diff_phi

                else: 
                    vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                    diff_r3   = self.xi[1:]**3 - self.xi[:-1]**3
                    diff_cost = 2.0
                    diff_phi  = self.zi[1:] - self.zi[:-1] 
                    for ix in range(self.nx):
                        for iz in range(self.nz):
                            vol[ix,0,iz] = 1./3. * diff_r3[ix] * diff_cost * diff_phi[iz]

            elif self.act_dim[2]==0:
                vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                diff_r3   = self.xi[1:]**3 - self.xi[:-1]**3
                diff_cost = np.cos(self.yi[:-1]) - np.cos(self.yi[1:])
                diff_phi  = 2.*np.pi
                for ix in range(self.nx):
                    for iy in range(self.ny):
                        vol[ix,iy,:] = 1./3. * diff_r3[ix] * diff_cost[iy] * diff_phi

            else:
                vol = np.zeros([self.nx, self.ny, self.nz], dtype=np.float64)
                diff_r3   = self.xi[1:]**3 - self.xi[:-1]**3
                diff_cost = np.cos(self.yi[:-1]) - np.cos(self.yi[1:])
                diff_phi  = self.zi[1:] - self.zi[:-1] 
                for ix in range(self.nx):
                    for iy in range(self.ny):
                        vol[ix,iy,:] = 1./3. * diff_r3[ix] * diff_cost[iy] * diff_phi
        else:
            print 'ERROR!'
            print "coordinate system '" + self.crd_sys+ "' is not yet supported"
            return 0

        return vol
# --------------------------------------------------------------------------------------------------
class radmc3dData(object):
    """RADMC-3D data class.
        Reading and writing dust density/temperature, gas density/temperature/velocity,
        generating a legacy vtk file for visualization.

    
    Attributes
    ----------
    
    grid      : radmc3dGrid 
                Instance of the radmc3dGrid class, contains the spatial and frequency grids
    
    rhodust   : ndarray
                Dust density in g/cm^3 
    
    dusttemp  : ndarray
                Dust temperature in K 
    
    rhogas    : ndarray
                Gas density in g/cm^3
    
    ndens_mol : ndarray
                Number density of the molecule [molecule/cm^3]
    
    ndens_cp  : ndarray
                Number density of the collisional partner [molecule/cm^3]
    
    gasvel    : ndarray
                Gas velocity in cm/s 
    
    gastemp   : ndarray
                Gas temperature in K
    
    vturb     : ndarray
                Mictroturbulence in cm/s
    
    taux      : ndarray
                Optical depth along the x (cartesian) / r (cylindrical) / r (spherical) dimension
    
    tauy      : ndarray
                Optical depth along the y (cartesian) / theta (cylindrical) / theta (spherical) dimension
    
    tauz      : ndarray
                Optical depth along the z (cartesian) / z (cylindrical) / phi (spherical) dimension
    
    sigmadust : ndarray
                Dust surface density in g/cm^2
    
    sigmagas  : ndarray
                Gas surface density in molecule/cm^2 (or g/cm^2 depending on the dimension of rhogas)
    """
    
    def __init__(self, grid=None):

        if grid:
            self.grid = copy.deepcopy(grid)
        else:
            self.grid = radmc3dGrid()

        self.rhodust   = -1
        self.dusttemp  = -1
        self.rhogas    = -1
        self.ndens_mol = -1
        self.ndens_cp  = -1
        self.gasvel    = -1
        self.gastemp   = -1
        self.vturb     = -1
        self.taux      = -1
        self.tauy      = -1
        self.tauz      = -1
        self.sigmadust = -1
        self.sigmagas  = -1
# --------------------------------------------------------------------------------------------------
    def _scalarfieldWriter(self, data=None, fname='', binary=True):
        """Writes a scalar field to a file.

        Parameters
        ----------
        
        data   : ndarray
                Scalar variable to be written
        
        fname  : str
                Name of the file containing a scalar variable
        
        binary : bool
                If True the file will be in binary format, if False the file format is formatted ASCII text
        
        """

        wfile = open(fname, 'w')
        if binary:
            if len(data.shape)==3:
                hdr = np.array([1, 8, self.grid.nx*self.grid.ny*self.grid.nz], dtype=int)
            elif len(data.shape)==4:
                hdr = np.array([1, 8, self.grid.nx*self.grid.ny*self.grid.nz,  data.shape[3]], dtype=int)
            hdr.tofile(wfile)
            # Now we need to flatten the dust density array since the Ndarray.tofile function writes the 
            # array always in C-order while we need Fortran-order to be written
            if len(data.shape)==4:
                data = np.swapaxes(data,0,3)
                data = np.swapaxes(data,1,2)
                data.tofile(wfile)
            elif len(data.shape)==3:
                data = np.swapaxes(data,0,2)
                data.tofile(wfile)
            else:
                print 'ERROR'
                print 'Unknown array shape  : '
                print data.shape
                return
        else:
            
            if len(data.shape)==3:
                hdr = np.array([1, self.grid.nx*self.grid.ny*self.grid.nz], dtype=int)
                hdr.tofile(wfile, sep=" ", format="%d\n")
                # Now we need to flatten the dust density array since the Ndarray.tofile function writes the 
                # array always in C-order while we need Fortran-order to be written
                data = np.swapaxes(data,0,2)
                data.tofile(wfile, sep=" ", format="%.9e\n")


            elif len(data.shape)==4:
                hdr = np.array([1, self.grid.nx*self.grid.ny*self.grid.nz,  self.rhodust.shape[3]], dtype=int)
                hdr.tofile(wfile, sep=" ", format="%d\n")
                # Now we need to flatten the dust density array since the Ndarray.tofile function writes the 
                # array always in C-order while we need Fortran-order to be written
                data = np.swapaxes(data,0,3)
                data = np.swapaxes(data,1,2)
                data.tofile(wfile, sep=" ", format="%.9e\n")
            else:
                print 'ERROR'
                print 'Unknown array shape  : '
                print data.shape
                return
            
        wfile.close()

# --------------------------------------------------------------------------------------------------
    def _scalarfieldReader(self, fname='', binary=True):
        """Reads a scalar field from file.

        Parameters
        ----------
        
        fname  : str 
                Name of the file containing a scalar variable
        
        binary : bool
                If True the file is in binary format, if False the file format is formatted ASCII text
        
        Returns
        -------
        
        Returns a numpy Ndarray with the scalar field
        """

        if binary:
            # hdr[0] = format number
            # hdr[1] = data precision (4=single, 8=double)
            # hdr[2] = nr of cells
            # hdr[3] = nr of dust species
            hdr = np.fromfile(fname, count=4, dtype=int)
            if hdr[2]!=(self.grid.nx*self.grid.ny*self.grid.nz):
                print ' ERROR'
                print ' Number of grid points in '+fname+' is different from that in amr_grid.inp'
                print npoints
                print hdr[2]
                return

            if hdr[1]==8:
                data = np.fromfile(fname, count=-1, dtype=np.float64)
            elif hdr[1]==4:
                data = np.fromfile(fname, count=-1, dtype=float)
            else:
                print 'ERROR'
                print 'Unknown datatype in '+fname
                return
            

            if data.shape[0]==(hdr[2]+3):
                data = np.reshape(data[3:], [1, self.grid.nz,self.grid.ny,self.grid.nx])
            elif data.shape[0]==(hdr[2]*hdr[3]+4):
                data = np.reshape(data[4:], [hdr[3],self.grid.nz,self.grid.ny,self.grid.nx])

            
            #data = reshape(data, [hdr[3],self.grid.nz,self.grid.ny,self.grid.nx])
            # We need to change the axis orders as Numpy always writes binaries in C-order while RADMC-3D
            # uses Fortran-order
            data = np.swapaxes(data,0,3)
            data = np.swapaxes(data,1,2)

        else:
            rfile = -1
            try :
                rfile = open(fname, 'r')
            except:
                print 'Error!' 
                print fname+' was not found!'
                
             
            if (rfile!=(-1)):

                hdr = np.fromfile(fname, count=3, sep="\n", dtype=int)
                
                if ((self.grid.nx * self.grid.ny * self.grid.nz)!=hdr[1]):
                    print 'Error!'
                    print 'Number of grid points in amr_grid.inp is not equal to that in '+fname
                else:

                    data = np.fromfile(fname, count=-1, sep="\n", dtype=np.float64)
                    if data.shape[0]==hdr[1]+2:
                        data = np.reshape(data[2:], [1, self.grid.nz,self.grid.ny,self.grid.nx])
                    elif data.shape[0]==hdr[1]*hdr[2]+3:
                        data = np.reshape(data[3:], [hdr[2],self.grid.nz,self.grid.ny,self.grid.nx])
                    # We need to change the axis orders as Numpy always reads  in C-order while RADMC-3D
                    # uses Fortran-order
                    data = np.swapaxes(data,0,3)
                    data = np.swapaxes(data,1,2)
            
            else:
                data = -1

            if rfile!=(-1):
                rfile.close()
        return data

# --------------------------------------------------------------------------------------------------
    def  getTauOneDust(self, idust=0, axis='', kappa=0.):
        """Calculates the optical depth of a single dust species along any given combination of the axes.

        Parameters
        ----------
        
        idust : int
                Index of the dust species whose optical depth should be calculated
        
        axis  : str
                Name of the axis/axes along which the optical depth should be calculated 
                (e.g. 'x' for the first dimension or 'xyz' for all three dimensions)
        
        kappa : float
                Mass extinction coefficients of the dust species at the desired wavelength
        
        Returns
        -------
        
        Returns a dictionary with the following keys
            
            taux  : ndarray
                    optical depth along the first dimension
            tauy  : ndarray
                    optical depth along the second dimension
            
            (tauz is not yet implemented)
        """
    
        # Check along which axis should the optical depth be calculated
        do_taux = False
        do_tauy = False

        if axis.find('x')>=0 : do_taux = True
        if axis.find('y')>=0 : do_tauy = True
      
        # Calculate the optical depth along the x-axis (r in spherical coordinates)
        if do_taux:
            taux = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64)
            diff_x    = self.grid.xi[1:] - self.grid.xi[:-1]
            taux[0,:,:] = self.rhodust[0,:,:,idust] * kappa * diff_x[0] 
            for ix in range(1,self.grid.nx):
                taux[ix,:,:] = taux[ix-1,:,:] + self.rhodust[ix,:,:,idust] * kappa * diff_x[ix] 
        else:
            taux = [-1.]
        
        # Calculate the optical depth along the theta in spherical coordinates
        # Warning the formulation below is valid only in spherical coordinate sytem

        dum_x = np.zeros([self.grid.nx, self.grid.nz], dtype=np.float64)
        for iz in range(self.grid.nz):
            dum_x[:,iz] = self.grid.x

        if do_tauy:
            tauy = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64)
            diff_y    = self.grid.yi[1:] - self.grid.yi[:-1]
            tauy[:,0,:] = self.rhodust[:,0,:,idust] * kappa * diff_y[0] * dum_x
            for iy in range(1,self.grid.ny):
                tauy[:,iy,:] = tauy[:,iy-1,:] + self.rhodust[:,iy,:,idust] * kappa * diff_y[iy] * dum_x
        else:
            tauy = [-1.]
        return {'taux':taux, 'tauy':tauy}

# --------------------------------------------------------------------------------------------------
    def  getTau(self, idust=[], axis='xy', wav=0., kappa=None, old=False):
        """Calculates the optical depth along any given combination of the axes.

        Parameters
        ----------
        
        idust : list
                List of dust component indices whose optical depth should be calculated
                If multiple indices are set the total optical depth is calculated summing 
                over all dust species in idust
        
        axis  : str
                Name of the axis/axes along which the optical depth should be calculated 
                (e.g. 'x' for the first dimension or 'xyz' for all three dimensions)
        
        wav   : float
                Wavelength at which the optical depth should be calculated
        
        kappa : bool
                If set it should be a list of mass extinction coefficients at the desired wavelength
                The number of elements in the list should be equal to that in the idust keyword

        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
        # Check if the input idust indices can be found in rhoudust 
        if len(self.rhodust.shape)==3: 
            ndust = 1
        else:
            ndust = self.rhodust.shape[3]

        if len(idust)==0:
            idust = np.arange(ndust)

        if max(idust)>ndust:
            print 'ERROR'
            print ' There are less number of dust species than some of the indices in idust'
            return -1


        scatmat = None
        # If the kappa keyword is set it should be used during the optical depth calculation
        if kappa:
            # Safety check
            if len(kappa)!=len(idust):
                print 'ERROR'
                print ' The number of kappa values should be identical to the number of specified dust species '
                return -1
        else: 
            if not old:
                # Read the master opacity file to get the dustkappa file name extensions
                dum = radmc3dDustOpac()
                mo  = dum.readMasterOpac()
                dummy_ext = mo['ext']
                scatmat = mo['scatmat']
                if len(dummy_ext)<=max(idust):
                    print 'ERROR'
                    print 'There are less dust species specified in dustopac.inp than some of the specified idust indices'
                    return -1
                else:
                    ext = [dummy_ext[i] for i in idust]

        if axis.find('x')>=0:
            self.taux = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64) 
        if axis.find('y')>=0:
            self.tauy = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=np.float64) 

        for i in idust:

            if kappa==None:
                if old:
                    opac = readOpac(ext=[("%d"%(i+1))], old=True)
                else:
                    opac = readOpac(ext=ext[i], scatmat=scatmat)
                
                if opac.ext==[]:
                    return -1
                else:
                    kabs = 10.**np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0]), np.log10(opac.kabs[0]))
                if opac.ksca[0][0]>0:
                    ksca = 10.**np.interp(np.log10(np.array(wav)), np.log10(opac.wav[0]), np.log10(opac.ksca[0]))
                else:
                    ksca = np.array(kabs)*0.
            
                print ' Opacity at '+("%.2f"%wav)+'um : ', kabs+ksca
                dum  = self.getTauOneDust(i, axis=axis, kappa=kabs + ksca)
            else:
                dum  = self.getTauOneDust(i, axis=axis, kappa=kappa[i])

            if axis.find('x')>=0:
                self.taux = self.taux + dum['taux']
            if axis.find('y')>=0:
                self.tauy = self.tauy + dum['tauy']

# --------------------------------------------------------------------------------------------------
    def readDustDens(self, fname='', binary=True, old=False):
        """Reads the dust density.

        Parameters
        ----------
        
        fname : str, optional
                Name of the file that contains the dust density. If omitted 'dust_density.inp' is used
                (or if binary=True the 'dust_density.binp' is used).
        
        binary : bool, optional 
                If true the data will be read in binary format, otherwise the file format is ascii
        
        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
   
        if (self.grid.nx==-1):
            self.grid.readGrid(old=old)
            
        print 'Reading dust density'

        # 
        # Read radmc3d output 
        #
        if not old:

            if binary:
                if fname=='':
                    fname = 'dust_density.binp'
            else:
                if fname=='':
                    fname = 'dust_density.inp'
                
            self.rhodust = self._scalarfieldReader(fname=fname, binary=binary)
        # 
        # Read the output of the previous 2d version of the code
        #
        
        else:
            try :
                rfile = open('dustdens.inp', 'r')
            except:
                print 'Error!' 
                print 'dustdens.inp was not found!'
                return 
           
            dum = rfile.readline().split()
            ndust = int(dum[0])
            nr    = int(dum[1])
            nt    = int(dum[2])
            imirt = int(dum[3])

            rfile.readline()
            self.rhodust = np.zeros([nr,nt*2,1,ndust], dtype=float)

            for idust in range(ndust):
                for ix in range(nr):
                    for iy in range(nt):
                        self.rhodust[ix,iy,0,idust] = float(rfile.readline())
                        self.rhodust[ix,self.grid.ny-1-iy,0,idust] = self.rhodust[ix,iy,0,idust]

            rfile.close()
            
# --------------------------------------------------------------------------------------------------
    def readDustTemp(self, fname='', binary=True, old=False):
        """Reads the dust temperature.

        Parameters
        ----------
        
        fname : str, optional
                Name of the file that contains the dust temperature. 
        
        binary : bool, optional 
                If true the data will be read in binary format, otherwise the file format is ascii
        """
       

        if (self.grid.nx==-1):
            self.grid.readGrid(old=old)
            
        print 'Reading dust temperature'

        if not old:
            if binary:
                if fname=='':
                    fname = 'dust_temperature.bdat'
            else:
                if fname=='':
                    fname = 'dust_temperature.dat'
                
            self.dusttemp = self._scalarfieldReader(fname=fname, binary=binary)
        else:
            try :
                rfile = open('dusttemp_final.dat', 'r')
            except:
                print 'Error!' 
                print 'dusttemp_final.dat was not found!'
                return 
           
            dum = rfile.readline().split()
            ndust = int(dum[0])
            nr    = int(dum[1])
            nt    = int(dum[2])
            imirt = int(dum[3])

            rfile.readline()
            self.dusttemp = np.zeros([nr,nt*2,1,ndust], dtype=float)

            for idust in range(ndust):
                rfile.readline()
                for ix in range(nr):
                    for iy in range(nt):
                        self.dusttemp[ix,iy,0,idust] = float(rfile.readline())
                        self.dusttemp[ix,self.grid.ny-1-iy,0,idust] = self.dusttemp[ix,iy,0,idust]

            rfile.close()
# --------------------------------------------------------------------------------------------------
    def readGasVel(self, fname='', binary=True):
        """Reads the gas velocity.  
        
        Parameters
        -----------
        
        fname : str, optional
                Name of the file that contains the gas velocity
                If omitted 'gas_velocity.inp' (if binary=True 'gas_velocity.binp')is used.
        
        binary : bool
                If true the data will be read in binary format, otherwise the file format is ascii

        """

        if binary:
            if fname=='':
                fname = 'gas_velocity.binp'
            if (self.grid.nx==-1):
                self.grid.readGrid()

            print 'Reading gas velocity'
            
            try :
                rfile = open(fname, 'r')
            except:
                print 'Error!' 
                print fname+' was not found!'

            if (rfile!=(-1)):            
                hdr = np.fromfile(fname, count=3, dtype=int)
                if (hdr[2]!=self.grid.nx*self.grid.ny*self.grid.nz):
                    print 'ERROR'
                    print 'Number of grid points in '+fname+' is different from that in amr_grid.inp'
                    print self.grid.nx, self.grid.ny, self.grid.nz
                    print hdr[1]
                    return

                if hdr[1]==8:
                    self.gasvel = np.fromfile(fname, count=-1, dtype=np.float64)
                elif hdr[1]==4:
                    self.gasvel = np.fromfile(fname, count=-1, dtype=float)
                else:
                    print 'ERROR'
                    print 'Unknown datatype in '+fname
                    return
                self.gasvel = np.reshape(self.gasvel[3:], [self.grid.nz,self.grid.ny,self.grid.nx,3])
                self.gasvel = np.swapaxes(self.gasvel, 0, 2)

            else:
                self.gasvel=-1
            

        else:
            if fname=='':
                fname = 'gas_velocity.inp'

            if (self.grid.nx==-1):
                self.grid.readGrid()

            print 'Reading gas velocity'

            rfile = -1

            try :
                rfile = open(fname, 'r')
            except:
                print 'Error!' 
                print fname+' was not found!'
                
            if (rfile!=(-1)):            
                dum = rfile.readline()
                dum = int(rfile.readline())
                
                if ((self.grid.nx * self.grid.ny * self.grid.nz)!=dum):
                    print 'Error!'
                    print 'Number of self.grid.points in amr_grid.inp is not equal to that in gas_velocity.inp'
                else:
                    
                    self.gasvel = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz, 3], dtype=np.float64)
                    
                    for k in range(self.grid.nz):
                        for j in range(self.grid.ny):
                            for i in range(self.grid.nx):
                                dum = rfile.readline().split()
                                self.gasvel[i,j,k,0] = float(dum[0])
                                self.gasvel[i,j,k,1] = float(dum[1])
                                self.gasvel[i,j,k,2] = float(dum[2])
    #                            self.gasvel[i,j,k,:] = [float(dum[i]) for i in range(3)]

            else:
                self.gasvel = -1                            

            rfile.close()
# --------------------------------------------------------------------------------------------------
    def readVTurb(self, fname='', binary=True):
        """Reads the turbulent velocity field. 
        
        Parameters
        ----------
        
        fname : str, optional 
                Name of the file that contains the turbulent velocity field
                If omitted 'microturbulence.inp' (if binary=True 'microturbulence.binp') is used.
        
        binary : bool 
                If true the data will be read in binary format, otherwise the file format is ascii
        """
        
        if (self.grid.nx==-1):
            self.grid.readGrid()
            
        print 'Reading microturbulence'

        if binary:
            if fname=='':
                fname = 'microturbulence.binp'
        else:
            if fname=='':
                fname = 'microturbulence.inp'
            
        self.vturb = self._scalarfieldReader(fname=fname, binary=binary)
       
# --------------------------------------------------------------------------------------------------
    def readGasDens(self,ispec='',binary=True):
        """Reads the gas density.

        Parameters
        ----------
        
        ispec : str 
                File name extension of the 'numberdens_ispec.inp' (or if binary=True 'numberdens_ispec.binp') file.


        binary : bool 
                If true the data will be read in binary format, otherwise the file format is ascii

        """
        
        if (self.grid.nx==-1):
            self.grid.readGrid()
           

        if binary:
            fname = 'numberdens_'+ispec+'.binp'
        else:
            fname = 'numberdens_'+ispec+'.inp'
            
        print 'Reading gas density ('+fname+')'
        self.ndens_mol = self._scalarfieldReader(fname=fname, binary=binary)
       
# --------------------------------------------------------------------------------------------------

    def readGasTemp(self, fname='', binary=True):
        """Reads the gas temperature.

        Parameters
        ----------
        
        fname : str,optional
                Name of the file that contains the gas temperature. If omitted 'gas_temperature.inp' 
                (or if binary=True 'gas_tempearture.binp') is used.
        
        binary : bool
                If true the data will be read in binary format, otherwise the file format is ascii
        """
      
        if (self.grid.nx==-1):
            self.grid.readGrid()
            
        print 'Reading gas temperature'

        if binary:
            if fname=='':
                fname = 'gas_temperature.binp'
        else:
            if fname=='':
                fname = 'gas_temperature.inp'
            
        self.gastemp = self._scalarfieldReader(fname=fname, binary=binary)

# --------------------------------------------------------------------------------------------------
    def writeDustDens(self, fname='', binary=True, old=False):
        """Writes the dust density.

        Parameters
        ----------
        
        fname : str, optional
                Name of the file into which the dust density should be written. If omitted 'dust_density.inp' is used.
        
        binary : bool
                If true the data will be written in binary format, otherwise the file format is ascii
    
        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
      
        # 
        # Write dust density for radmc3d
        #
        if not old:
            if fname=='':
                if binary:
                    fname = 'dust_density.binp'
                else:
                    fname = 'dust_density.inp'

            print 'Writing '+fname

            self._scalarfieldWriter(data=self.rhodust, fname=fname, binary=binary)

        # 
        # Write dust density for the previous 2D version of the code
        #
        else:
            if self.rhodust.shape[2]>1:
                print 'ERROR'
                print 'You are trying to write a 3D dust density structure for a 2D model'
                print 'The "old" keyword is set meaning that the input is meant for the previous 2D version of radmc'
                return

            if fname=='':
                fname = 'dustdens.inp'
            try :
                wfile = open(fname, 'w')
            except:
                print 'Error!' 
                print fname+' cannot be opened!'
                return 
          
            wfile.write("%d %d %d 1\n"%(self.rhodust.shape[3], self.grid.nx, self.grid.ny/2))# self.rhodust.shape[0], self.rhodust.shape[1]))
            wfile.write(" \n")
            for idust in range(self.rhodust.shape[3]):
                for ix in range(self.grid.nx):
                    for iy in range(self.grid.ny/2):
                        wfile.write("%.7e\n"%self.rhodust[ix,iy,0,idust])

            wfile.close()
# --------------------------------------------------------------------------------------------------
    def writeDustTemp(self, fname='', binary=True):
        """Writes the dust density.

        Parameters
        ----------
        
        fname : str, optional
                Name of the file into which the dust density should be written. If omitted 'dust_density.inp' is used.
        
        binary : bool
                If true the data will be written in binary format, otherwise the file format is ascii
        """
        if fname=='':
            if binary:
                fname = 'dust_temperature.bdat'
            else:
                fname = 'dust_temperature.dat'

        print 'Writing '+fname
        self._scalarfieldWriter(data=self.dusttemp, fname=fname, binary=binary)
    
# --------------------------------------------------------------------------------------------------
    def writeGasDens(self, fname='', ispec='',binary=True):
        """Writes the gas density.

        Parameters
        ----------
        
        fname  : str, optional
                 Name of the file into which the data will be written. If omitted "numberdens_xxx.inp" and
                 "numberdens_xxx.binp" will be used for ascii and binary format, respectively (xxx is the name of the molecule).
        
        ispec  : str
                 File name extension of the 'numberdens_ispec.inp' (if binary=True 'numberdens_ispec.binp') 
                 file into which the gas density should be written
        
        binary : bool
                 If true the data will be written in binary format, otherwise the file format is ascii
        """
        if ispec=='':
            print 'ERROR'
            print 'ispec keyword was not specified. This keyword is required to generate the '
            print "output file name 'numberdens_ispec.dat'" 
            return -1
        else:
            if fname=='':
                if binary:
                    fname = 'numberdens_'+ispec+'.binp'
                else:
                    fname = 'numberdens_'+ispec+'.inp'

            print 'Writing '+fname
            self._scalarfieldWriter(data=self.ndens_mol, fname=fname, binary=binary)
        
       
# --------------------------------------------------------------------------------------------------
    def writeGasTemp(self, fname='', binary=True):
        """Writes the gas temperature.

        Parameters
        ----------
        
        fname : str, optional 
                Name of the file into which the gas temperature should be written. If omitted 
                'gas_temperature.inp' (if binary=True 'gas_tempearture.binp') is used.
        
        binary : bool
                If true the data will be written in binary format, otherwise the file format is ascii
        """
        if fname=='':
            if binary:
                fname = 'gas_temperature.binp'
            else:
                fname = 'gas_temperature.inp'

        print 'Writing '+fname
        self._scalarfieldWriter(data=self.gastemp, fname=fname, binary=binary)
   
# --------------------------------------------------------------------------------------------------
    def writeGasVel(self, fname='', binary=True):
        """Writes the gas velocity.

        Parameters
        ----------
        
        fname  : str, optional
                Name of the file into which the gas temperature should be written. 
                If omitted 'gas_velocity.inp' (if binary=True 'gas_velocity.binp') is used.
        
        binary : bool
                If true the data will be written in binary format, otherwise the file format is ascii
        """
   
        if binary:
            if fname=='':
                fname = 'gas_velocity.binp'

            wfile = open(fname, 'w')
            hdr = np.array([1, 8, self.grid.nx*self.grid.ny*self.grid.nz], dtype=int)
            hdr.tofile(wfile)
            # Now we need to change the axis orders since the Ndarray.tofile function writes the 
            # array always in C-order while we need Fortran-order to be written
            self.gasvel = np.swapaxes(self.gasvel,0,2)
            self.gasvel.tofile(wfile)

            # Switch back to the original axis order
            self.gasvel = np.swapaxes(self.gasvel,0,2)
            wfile.close()
        else:
            if fname=='':
                fname = 'gas_velocity.inp'

            wfile = open(fname, 'w')
            
            wfile.write('%d\n'%1)
            wfile.write('%d\n'%(self.grid.nx*self.grid.ny*self.grid.nz))

            for iz in range(self.grid.nz):
                for iy in range(self.grid.ny):
                    for ix in range(self.grid.nx):
                        wfile.write("%9e %9e %9e\n"%(self.gasvel[ix,iy,iz,0], self.gasvel[ix,iy,iz,1], self.gasvel[ix,iy,iz,2]))
                    
            wfile.close()
        print 'Writing '+fname
# --------------------------------------------------------------------------------------------------
    def writeVTurb(self, fname='', binary=True):
        """Writes the microturbulence file.

        Parameters
        ----------
        
        fname : str, optional
                Name of the file into which the turubulent velocity field should be written. 
                If omitted 'microturbulence.inp' (if binary=True 'microturbuulence.binp') is used.
        
        binary : bool
                If true the data will be written in binary format, otherwise the file format is ascii
        """
   
        if fname=='':
            if binary:
                fname = 'microturbulence.binp'
            else:
                fname = 'microturbulence.inp'

        print 'Writing '+fname
        self._scalarfieldWriter(data=self.vturb, fname=fname, binary=binary)


# --------------------------------------------------------------------------------------------------
    def writeVTK(self, vtk_fname='', ddens=False, dtemp=False, idust=[0], \
                          gdens=False, gvel=False, gtemp=False):
        """Writes physical variables to a legacy vtk file.

        Parameters
        ----------
        
        vtk_fname : str
                    Name of the file to be written, if not specified 'radmc3d_data.vtk' will be used
        
        ddens     : bool
                    If set to True the dust density will be written to the vtk file
        
        dtemp     : bool
                    If set to True the dust temperature will be written to the vtk file
        
        idust     : list
                    List of indices that specifies which dust component should be written 
                    if not set then the first dust species (zero index) will be used
        
        gdens     : bool
                    If set to True the gas density will be written to the vtk file
        
        gtemp     : bool
                    If set to True the gas temperature will be written to the vtk file
        
        gvel      : bool
                    If set to True the gas velocity will be written to the vtk file
        """

        if (vtk_fname==''):
            vtk_fname = 'radmc3d_data.vtk'
        else:
            vtk_fname = str(vtk_fname)

#
# Get the grid 
#
        
        x  = self.grid.xi
        # For the theta axis I leave out the poles
        #  The current cell type is hexahedron and the cells near the pole are
        #    rather tetrahedra than hexahedra and this is not yet implemented
        y  = np.array(self.grid.yi[1:self.grid.nyi-1])
        z  = self.grid.zi
        nxi = x.shape[0]
        nyi = y.shape[0]
        nzi = z.shape[0]


#
# Gas velocity field (Should be corner centered)
# TODO
#  The lines below should be double checked and re-implemented as
#  the re-mapping of the cell centered velocity field to the cell corners
#  is physically not correct and very messy... 
#
        if gvel:
            vgas = np.zeros([nxi,nyi,nzi,3], dtype=np.float64)
            vgas[0:nxi-1,0:nyi,0:nzi-1,:]  = self.gasvel[:,1:nyi+1,:,:]
            vgas[nxi-1,:,:,:] = vgas[nxi-2,:,:,:]
            vgas[:,nyi-1,:,:] = vgas[:,nyi-2,:,:]
            vgas[:,:,nzi-1,:] = vgas[:,:,nzi-2,:]

# 
# Header 
# 

        
        wfile = open(vtk_fname, 'w')
        wfile.write('%s\n'%'# vtk DataFile Version 3.0')
        wfile.write('%s\n'%'RADMC-3D Data')
        wfile.write('%s\n'%'ASCII')
        wfile.write('%s\n'%'DATASET UNSTRUCTURED_GRID')

#
# Write out the coordinates of the cell corners 
# 
        wfile.write('%s\n'%('POINTS '+str(nxi*nyi*nzi).strip()+' double'))
        print 'Writing POINTS: '
        for ix in range(nxi):
            print ix, nxi
            for iy in range(nyi):
                for iz in range(nzi):
                    crd = crd_trans.ctrans_sph2cart([x[ix],z[iz],y[iy]])
                    wfile.write('%.9e %9e %9e\n'%(crd[0], crd[1], crd[2]))
            
# ---------------------------------------------------------------------------------------------
# Write out the indices of the cell interface mesh that define a
# hexahedron (VTK cell type #12)
# 
# The indexing of a hexahedron is as follows
#
#                  7________6
#                 /|      / |               
#                / |     /  |               
#               4_------5   |             z ^   ^ y
#               |  3____|___2               |  /
#               | /     |  /                | /
#               |/      | /                 |/
#               0-------1                   0-----> x
#
# ---------------------------------------------------------------------------------------------
 
        wfile.write('%s %d %d\n'%('CELLS ', ((nxi-1)*(nyi-1)*(nzi-1)), ((nxi-1)*(nyi-1)*(nzi-1))*9))


        for ix in range(nxi-1):
            print 'Writing CELL COORDINATES: ', ix, self.grid.nxi-2
            for iy in range(nyi-1):
                for iz in range(nzi-1):                
                
                    id1 = nzi*nyi*ix     + nzi*iy     + iz
                    id2 = nzi*nyi*ix     + nzi*(iy+1) + iz
                    id4 = nzi*nyi*ix     + nzi*iy     + ((iz+1) % (nzi-1))
                    id3 = nzi*nyi*ix     + nzi*(iy+1) + ((iz+1) % (nzi-1))
                    id5 = nzi*nyi*(ix+1) + nzi*iy     + iz
                    id6 = nzi*nyi*(ix+1) + nzi*(iy+1) + iz
                    id7 = nzi*nyi*(ix+1) + nzi*(iy+1) + ((iz+1) % (nzi-1))
                    id8 = nzi*nyi*(ix+1) + nzi*iy     + ((iz+1) % (nzi-1))
                
                
                    line = np.array([8,id1,id2,id3,id4,id5,id6,id7,id8])
                    line.tofile(wfile, sep=' ', format='%d')
                    wfile.write('\n')
#
# Now write out the type of each cell (#12)
#
        wfile.write('%s %d\n'%('CELL_TYPES', ((nxi-1)*(nyi-1)*(nzi-1))))

        for ix in range(nxi-1):
            for iy in range(nyi-1):
                for iz in range(nzi-1):      
                    wfile.write('%d\n'%12)
# 
# Now write out the corner centered velocities
#
                
        if gvel:
            wfile.write('%s %d\n'%('POINT_DATA', (nxi*nyi*nzi)))
            wfile.write('%s\n'%'VECTORS gas_velocity double')
            for ix in range(nxi):
                print 'Writing velocity : ', ix, nxi-1
                for iy in range(nyi):
                    for iz in range(nzi):      
                        vsph = np.array([vgas[ix,iy,iz,0],vgas[ix,iy,iz,2],vgas[ix,iy,iz,1]])
                        vxyz = crd_trans.vtrans_sph2cart([x[ix],z[iz],y[iy]], vsph)
                
                        wfile.write('%.9e %.9e %.9e\n'%(vxyz[0], vxyz[1], vxyz[2]))
                      
# 
# Write out the cell centered scalars
# 
        wfile.write('%s %d\n'%('CELL_DATA', ((nxi-1)*(nyi-1)*(nzi-1))))

    
        if ddens:
            for ids in idust:
                wfile.write('%s\n'%('SCALARS dust_density_'+str(int(ids))+' double'))
                wfile.write('%s\n'%'LOOKUP_TABLE default')

                for ix in range(nxi-1):
                    print 'Writing dust density : ', ix, nxi-2
                    for iy in range(nyi-1):
                        for iz in range(nzi-1):
                            wfile.write('%.9e\n'%self.rhodust[ix,iy,iz,ids])

        if dtemp:
            for ids in idust:
                wfile.write('%s\n'%('SCALARS dust_temperature_'+str(int(ids))+' double'))
                wfile.write('%s\n'%'LOOKUP_TABLE default')

                for ix in range(nxi-1):
                    print 'writing dust temperature : ', ix, nxi-2
                    for iy in range(nyi-1):
                        for iz in range(nzi-1):
                            wfile.write('%.9e\n'%self.dusttemp[ix,iy,iz,ids])


        if gdens:
            wfile.write('%s\n'%'SCALARS gas_numberdensity double')
            wfile.write('%s\n'%'LOOKUP_TABLE default')

            for ix in range(nxi-1):
                print 'writing gas density : ', ix, nxi-2
                for iy in range(nyi-1):
                    for iz in range(nzi-1):
                        wfile.write('%.9e\n'%self.ndens_mol[ix,iy,iz])

        if gtemp:
            for ids in idust:
                wfile.write('%s\n'%('SCALARS gas_temperature double'))
                wfile.write('%s\n'%'LOOKUP_TABLE default')

                for ix in range(nxi-1):
                    print 'writing dust temperature : ', ix, nxi-2
                    for iy in range(nyi-1):
                        for iz in range(nzi-1):
                            wfile.write('%.9e\n'%self.gastemp[ix,iy,iz])
                            
# --------------------------------------------------------------------------------------------------
# Close the file
# --------------------------------------------------------------------------------------------------
        wfile.close()

# --------------------------------------------------------------------------------------------------
    def getSigmaDust(self, idust=-1):
        """Calculates the dust surface density.
        
        Parameters
        ----------
        
        idust : int, optional
                Index of the dust species for which the surface density should be calculated 
                if omitted the calculated surface density will be the sum over all dust species
        """

        # Calculate the volume of each grid cell
        vol  = self.grid.getCellVolume()
        # Dustmass in each grid cell
        if len(self.rhodust)>3:
            if idust>=0:
                mass = vol * self.rhodust[:,:,:,idust]
            else:
                mass = vol * self.rhodust.sum(3)
        else:
            mass = vol * self.rhodust

        # Calculate the surface of each grid facet in the midplane
        surf     = np.zeros([self.grid.nx, self.grid.nz], dtype=np.float64)
        diff_r2  = (self.grid.xi[1:]**2 - self.grid.xi[:-1]**2) * 0.5
        diff_phi = self.grid.zi[1:] - self.grid.zi[:-1]
        for ix in range(self.grid.nx):
            surf[ix,:] = diff_r2[ix] * diff_phi

        
        # Now get the surface density 
        dum = np.squeeze(mass.sum(1))
        self.sigmadust = dum / np.squeeze(surf)

# --------------------------------------------------------------------------------------------------
    def getSigmaGas(self):
        """Calculates the gas surface density.
        This method uses radmc3dData.rhogas to calculate the surface density, thus the 
        unit of surface density depends on the unit of radmc3dData.rhogas (g/cm^2 or molecule/cm^2)
        """

        # Calculate the volume of each grid cell
        vol  = self.grid.getCellVolume()
        # Total number of molecules / gas mass in each grid cell
        mass = vol * self.rhogas
        # Calculate the surface are of each grid facet in the midplane
        surf     = np.zeros([self.grid.nx, self.grid.nz], dtype=np.float64)
        diff_r2  = (self.grid.xi[1:]**2 - self.grid.xi[:-1]**2) * 0.5
        diff_phi = self.grid.zi[1:] - self.grid.zi[:-1]  
        for ix in range(self.grid.nx):
            surf[ix,:] = diff_r2[ix] * diff_phi


        # Now get the surface density 
        dum = np.squeeze(mass.sum(1))
        self.sigmagas = dum / np.squeeze(surf)

# --------------------------------------------------------------------------------------------------
class radmc3dRadSources(object):
    """Class of the radiation sources.
    Currently discrete stars and continuous starlike source, the latter only in spherical coordinates.


    Attributes
    ----------

    wav          : ndarray
                    Wavelength for the stellar spectrum
    
    freq         : ndarray
                    Frequency for the stellar spectrum
    
    nwav         : int
                    Number of wavelenghts in the stellar spectrum
    
    nfreq        : int
                    Number of frequencies in the stellar spectrum

    mstar        : list
                    List of stellar masses
    
    tstar        : list
                    List of stellar effective temperatures
    
    rstar        : list
                    List of stellar radii
    
    lstar        : list 
                    List of stellar luminosities 
    
    nstar        : int
                    Number of stars
    
    pstar        : list
                    Each element of the list contains a three element list, the cartesian coordinates of the stars
    
    fnustar      : ndarray
                    Stellar spectrum (flux@1pc in erg/s/cm/cm/Hz)

    csdens       : ndarray
                    Stellar density for continuous starlike source
    
    csntemplate  : int
                    Number of stellar templates
    
    cstemp       : ndarray
                    Stellar template
    
    cstemptype   : int
                    Stellar template type 1 - Blackbody given by the effective temperature 2 - Frequency dependent spectrum
    
    cststar      : ndarray
                    Stellar effective temperature
    
    csmstar      : ndarray
                    Stellar mass
    
    csrstar      : ndarray
                    Stellar radius

    tacc         : ndarray
                    Effective temperature of a viscous accretion disk as a function of cylindrical radius
    
    accrate      : float
                    Accretion rate of the viscous accretion disk [g/s]
    
    fnuaccdisk   : ndarray
                    Spatially integrated frequency-dependent flux density of the accretion disk @ 1pc distance
    
    tspot        : float
                    Temperature of the hot spot / boundary layer on the stellar surface
    
    starsurffrac : float
                    Fraction of the stellar surface covered by the hot spot / boundary layer
    
    fnustpot     : ndarray
                        Frequency-dependent flux density of the hot spot / boundary layer @ 1pc distance

    """
    def __init__(self, ppar=None, grid=None):


        # Spatial and frequency grid
        self.grid     = grid

        # Discrete stars

        self.mstar     = []  
        self.tstar     = []  
        self.rstar     = []  
        self.lstar     = []  
        self.nstar     = 0
        self.pstar     = []
        self.fnustar   = []

        # Continuous starlike source 
        self.csdens    = []
        self.csntemplate = 0
        self.cstemp    = []
        self.cstemptype = 1
        self.cststar   = []
        self.csmstar   = []
        self.csrstar   = []

        # Viscous accretion disk
        self.incl_accretion = False
        self.tacc     = []
        self.accrate  = 0.
        self.fnuaccdisk = []

        # Hot spot / boundary layer - to model accretion in YSOs
        self.tspot    = 0.
        self.starsurffrac = 0.
        self.fnuspot  = []

        if ppar:
            if type(ppar['mstar']).__name__=='list':
                self.mstar = ppar['mstar']
            else:
                self.mstar = [ppar['mstar']]
                
            if type(ppar['tstar']).__name__=='list':
                self.tstar = ppar['tstar']
            else:
                self.tstar = [ppar['tstar']]

            if type(ppar['rstar']).__name__=='list':
                self.rstar = ppar['rstar']
            else:
                self.rstar = [ppar['rstar']]
            
            for istar in range(self.nstar):
                self.lstar.append(4.*np.pi*self.rstar[istar]**2. * ss* self.tstar[istar]**4.)
            self.pstar = ppar['pstar']

            if ppar.has_key('incl_cont_stellarsrc'):
                self.incl_accretion = ppar['incl_cont_stellarsrc']
            else:
                self.incl_accretion = False

            if ppar.has_key('accrate'):
                self.accrate = ppar['accrate']
            else:
                self.accrate = 0.
# --------------------------------------------------------------------------------------------------

    def findPeakStarspec(self):

        """Calculates the peak wavelength of the stellar spectrum.
       
        Returns
        -------
        
        The peak wavelength of the stellar spectrum in nu*Fnu for all 
            stars as a list
        """
   
        pwav = []
   
        for istar in range(self.nstar):
            ii = (self.fnustar[:,istar]*self.grid.freq).argmax()
            pwav.append(self.grid.wav[ii])

        return pwav

# --------------------------------------------------------------------------------------------------
    def readStarsinp(self, fname=''):
        """Reads the data of discrete stellar sources from the stars.inp file.

        Parameters
        ----------
        
        fname : str, optional
                File name of the file that should be read (if omitted stars.inp will be used)
        """
        
        if fname=='':
            fname = 'stars.inp'

        try:
            rfile = open(fname, 'r')
        except:
            print ' ERROR '
            print fname+' cannot be opened '
            return

        dum = rfile.readline()
        iformat = int(dum)
        if iformat!=2:
            print ' ERROR '
            print ' Unknown file format '
            print ' Format number : ', iformat
            rfile.close()
            return

        dum = rfile.readline().split()
        self.nstar = int(dum[0])
        self.grid  = radmc3dGrid()
        self.grid.nwav  = int(dum[1])
        self.grid.nfreq = self.grid.nwav
        self.rstar = []
        self.mstar = []
        self.tstar = []
        for istar in range(self.nstar):
            dum = rfile.readline().split()
            self.rstar.append(float(dum[0]))
            self.mstar.append(float(dum[1]))
            self.pstar.append([float(dum[2]), float(dum[3]), float(dum[4])])
        
        dum = rfile.readline()
        wav = []
        for ilam in range(self.grid.nwav):
            dum = rfile.readline()
            wav.append(float(dum))

        self.grid.wav = np.array(wav, dtype=float)
        self.grid.freq = cc/self.grid.wav*1e4
        dum = rfile.readline()
        dum = rfile.readline()
       
         
        # If we have only the stellar temperature
        if float(dum)<0:
            self.tstar.append(-float(dum))
            for istar in range(1,self.nstar):
                self.tstar.append(-float(dum))
            # 
            # Now calculates the stellar spectrum
            #
            self.getStarSpectrum()
        
        # If we have the full frequency-dependent spectrum 
        else:
            self.fnustar = np.zeros([self.grid.nfreq, self.nstar], dtype=float)
            self.fnustar[0,0] = float(dum)
            for ifreq in range(1,self.grid.nfreq):
                self.fnustar[ifreq, istar] = float(rfile.readline())
            
            for istar in range(1,self.nstar):
                for ifreq in range(self.grid.nfreq):
                    self.fnustar[ifreq, istar] = float(rfile.readline())

        rfile.close()



# --------------------------------------------------------------------------------------------------
    def writeStarsinp(self, ppar=None, wav=None, freq=None, old=False):
        """Writes the input file for discrete stellar sources (stars.inp).

        Parameters
        ----------
        
        ppar  : dictionary
                Dictionary containing all parameters of the model (only mandatory if accretion is switched on)

        wav   : ndarray, optional
                Wavelength grid for the stellar spectrum
        
        freq  : ndarray, optional
                Frequency grid for the stellar spectrum (either freq or wav should be set)
        
        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """

        if not old:
            if freq!=None:
                self.grid.wav  = cc/np.array(freq)*1e4
                self.grid.freq = np.array(freq)
                self.grid.nwav = self.grid.wav.shape[0]
                self.grid.nfreq = self.grid.nwav

            if wav!=None:
                self.grid.wav = np.array(wav)
                self.grid.freq = cc/self.grid.wav*1e4
                self.grid.nwav = self.grid.wav.shape[0]
                self.grid.nfreq = self.grid.nwav
            
            self.nstar = len(self.rstar)
            #self.pstar = ppar['pstar']
           

            print 'Writing stars.inp'
            wfile = open('stars.inp', 'w')
            wfile.write('%d\n'%2)
            wfile.write('%d %d\n'%(self.nstar,self.grid.nwav))
            
            if (self.nstar>1):
                for istar in range(self.nstar):
                    wfile.write('%.9e %.9e %.9e %.9e %.9e\n'%(self.rstar[istar], self.mstar[istar],
                        self.pstar[istar][0],self.pstar[istar][1],self.pstar[istar][2]))
            else:
                wfile.write('%.9e %.9e %.9e %.9e %.9e\n'%(self.rstar[0], self.mstar[0],
                    self.pstar[0],self.pstar[1],self.pstar[2]))

            wfile.write('%s\n'%' ')
            for ilam in range(self.grid.nwav):
                wfile.write('%.9e\n'%self.grid.wav[ilam])
            wfile.write('%s\n'%' ')

            # If accretion is active write the sum of the spot and the star emission
            # NOTE, for now the spot emission is only added to the first star in the list
            if self.incl_accretion:
                # Get the stellar spectrum
                self.getStarSpectrum(tstar=self.tstar, rstar=self.rstar)
                # Get the spot spectrum
                self.getSpotSpectrum(ppar=ppar)
                
                # Write out the spectrum of the first star with the additional spot luminosity  
                for ilam in range(self.grid.nwav):
                    wfile.write('%.9e\n'%(self.fnustar[ilam,0]+self.fnuspot[ilam]))
                
                # Now write the spectrum of all other discrete stars without the spot emission
                for istar in range(1,self.nstar):
                    for ilam in range(self.grid.nwav):
                        wfile.write('%.9e\n'%(self.fnustar[ilam,istar]))


            else:
                self.getStarSpectrum(ppar=ppar)
                for istar in range(self.nstar):
                    if self.staremis_type[istar].strip().lower()=="blackbody":
                        wfile.write('%.9e\n'%(-self.tstar[istar]))
                    else:
                        for ilam in range(self.grid.nwav):
                            wfile.write('%.9e\n'%(self.fnustar[ilam,istar]))
            wfile.close()
        else:
            if freq!=None:
                self.grid.wav  = cc/np.array(freq)*1e4
                self.grid.freq = np.array(freq)
                self.grid.nwav = self.grid.wav.shape[0]
                self.grid.nfreq = self.grid.nwav

            if wav!=None:
                self.grid.wav = np.array(wav)
                self.grid.freq = cc/self.grid.wav*1e4
                self.grid.nwav = self.grid.wav.shape[0]
                self.grid.nfreq = self.grid.nwav
           

            self.nstar = len(self.rstar)
            self.getStarSpectrum(ppar=ppar)
            
            if self.grid.freq[-1]<self.grid.freq[0]:
                self.grid.freq = self.grid.freq[::-1]
                self.grid.wav  = self.grid.wav[::-1]
                for istar in range(self.nstar):
                    self.fnustar[:,istar] = self.fnustar[::-1,istar]
                

            print 'Writing starinfo.inp'
            fname = 'starinfo.inp'
            try :
                wfile = open(fname, 'w')
            except:
                print 'Error!' 
                print fname+' cannot be opened!'
                return 
            wfile.write("1\n")
            wfile.write("%.7e\n"%ppar['rstar'][0])
            wfile.write("%.7e\n"%ppar['mstar'][0])
            wfile.write("%.7e\n"%ppar['tstar'][0])
            wfile.close()

            print 'Writing starspectrum.inp'
            fname = 'starspectrum.inp'
            try :
                wfile = open(fname, 'w')
            except:
                print 'Error!' 
                print fname+' cannot be opened!'
                return 


            wfile.write("%d\n"%self.grid.nfreq)
            wfile.write(" \n")
            for i in range(self.grid.nfreq):
                wfile.write("%.7e %.7e\n"%(self.grid.freq[i], self.fnustar[i,0]))
            wfile.close()

# --------------------------------------------------------------------------------------------------
    def getStarSpectrum(self, tstar=None, rstar=None, lstar=None, mstar=None, ppar=None, grid=None):
        """Calculates a blackbody stellar spectrum.

        Parameters
        ----------
        
        tstar : list
                Effective temperature of the stars in [K]
        
        rstar : list
                Radius of the stars in [cm]
        
        lstar : list
                Bolometric luminosity of the star [erg/s] (either rstar or lstar should be given)
        
        mstar : list
                Stellar mass in [g] (only required if an atmosphere model is used to calculate logg)

        ppar  : dictionary
                Dictionary containing all input parameters
        
        grid  : radmc3dGrid, optional
                An instance of a radmc3dGrid class containing the spatial and wavelength grid
        """
#
# Check the input which parameters are set and which should be calculated
#

        if grid!=None:
            self.grid = grid

        if ppar!=None:
            if tstar==None:
                tstar = ppar['tstar']
            if rstar==None:
                rstar = ppar['rstar']
            if mstar==None:
                self.mstar = ppar['mstar']

        if mstar:
            self.mstar = mstar

        if tstar:
            if type(tstar).__name__!='list':
                tstar = [tstar]
            dum1 = len(tstar)
            if lstar and rstar: 
                print 'ERROR'
                print ' Only two of the input variables tstar, rstar, lstar should be set not all three'
                return 0
            elif lstar:
                if len(lstar)!=dum1:
                    print 'ERROR'
                    print 'lstar and tstar have different number of elements'
                    return 0
                else:
                    self.tstar = np.array(tstar)
                    self.lstar = np.array(lstar)
                    self.nstar = self.lstar.shape[0]
                    self.rstar = np.sqrt(self.lstar / (4.*np.pi*ss*self.tstar**4.))
        else:
            if lstar and rstar:
                if len(lstar)!=len(rstar):
                    print 'ERROR'
                    print 'lstar and rstar have different number of elements'
                    return 0
                else:
                    self.lstar = np.array(lstar)
                    self.rstar = np.array(rstar)
                    self.nstar = self.rstar.shape[0]
                    self.tstar = (self.lstar / (4.*np.pi*ss*self.rstar**2.))**0.25

            
        #
        # If we take blackbody spectrum
        #
       
        if ppar:
            if ppar.has_key('staremis_type'):
                self.staremis_type = ppar['staremis_type']
            else:
                self.staremis_type = []
                for i in range(self.nstar):
                    self.staremis_type.append("blackbody")
        else:
            self.staremis_type = []
            for i in range(self.nstar):
                self.staremis_type.append("blackbody")


        self.fnustar   = np.zeros([self.grid.nwav, len(self.tstar)], dtype=np.float64)
        for istar in range(self.nstar):
            if self.staremis_type[istar].strip().lower()=="blackbody":
                self.fnustar[:,istar]   = 2.*hh*self.grid.freq**3./cc**2/(np.exp(hh*self.grid.freq/kk/self.tstar[istar])-1.0) * \
                        np.pi * self.rstar[istar]**2. / pc**2.
            elif self.staremis_type[istar].strip().lower()=="kurucz":
                sta = StellarAtm()
                dum = sta.getAtmModel(teff=self.tstar[istar], mstar=self.mstar[istar], rstar=self.rstar[istar], iwav=self.grid.wav, model="kurucz")
                self.fnustar[:,istar]   = dum['lnu'] / (4. * np.pi * pc**2)
            
            elif self.staremis_type[istar].strip().lower()=="nextgen":
                sta = StellarAtm()
                dum = sta.getAtmModel(teff=self.tstar[istar], mstar=self.mstar[istar], rstar=self.rstar[istar], iwav=self.grid.wav, model="nextgen")
                self.fnustar[:,istar]   = dum['lnu'] / (4. * np.pi * pc**2)
            
            else:
                print 'ERROR'
                print 'Unknown stellar atmosphere model : ', self.staremis_type[istar]
                return
        
# --------------------------------------------------------------------------------------------------
    def getAccdiskTemperature(self, ppar=None, grid=None):
        """Calculates the effective temperature of a viscous accretion disk.

        Parameters
        ----------
        
        ppar : dictionary
            Dictionary containing all input parameters keys should include
                * mstar   : stellar mass
                * rstar   : stellar radius
                * accrate : accretion rate

            NOTE, that for the calculation of the effective disk temperature only the first
                star is used if more than one values are given in mstar and rstar. 
    
        grid : radmc3dGrid, optional 
            An instance of a radmc3dGrid class containing the spatial and wavelength grid

        """
      
        if grid!=None:
            self.grid = grid

        self.tacc = ((3.0 * gg * ppar['mstar'][0] * ppar['accrate']) / (8.0 * np.pi * self.grid.x**3 * ss) * \
                    (1.0 - (ppar['rstar'][0]/self.grid.x)**0.5))**0.25

# --------------------------------------------------------------------------------------------------
    def getSpotSpectrum(self, ppar=None, grid=None):
        """Calculates the spectrum of a hot spot / boundary layer on the stellar surface

        Parameters
        ----------

        ppar : dictionary
            Dictionary containing all input parameters keys should include
                * mstar   : stellar mass
                * rstar   : stellar radius
                * accrate : accretion rate

            NOTE, that for the calculation of the effective disk temperature only the first
                star is used if more than one values are given in mstar and rstar. 
        
        grid : radmc3dGrid, optional 
            An instance of a radmc3dGrid class containing the spatial and wavelength grid
        """
#
# Check the input which parameters are set and which should be calculated
#

        if grid!=None:
            self.grid = grid

        # Calculate the total spot luminosity assuming boundary layer accretion, i.e.
        # that half of the total accretion luminosity is released from the boundary layer

        if ppar.has_key('accrate'):
            tot_acclum = 0.5 * gg * ppar['mstar'][0] * ppar['accrate'] / ppar['rstar'][0]
            spotsurf   = 4.* np.pi * ppar['rstar'][0]**2 * ppar['starsurffrac']
            self.starsurffrac = ppar['starsurffrac']
            if spotsurf==0.:
                self.tspot = 0.
                
                # Now calculate the spot spectrum (i.e. the flux density @ 1pc distance)
                self.fnuspot = np.zeros(self.grid.nfreq, dtype=float)
            else:
                self.tspot = (tot_acclum / spotsurf / ss)**0.25

                # Now calculate the spot spectrum (i.e. the flux density @ 1pc distance)
                self.fnuspot = np.pi * ppar['rstar'][0]**2 * ppar['starsurffrac'] / pc**2 * \
                                2. * hh * self.grid.freq**3/cc**2 / (np.exp(hh*self.grid.freq/kk/self.tspot)-1.0)
        else:
            self.fnuspot = np.zeros(self.grid.nfreq, dtype=float)

# --------------------------------------------------------------------------------------------------
    def getTotalLuminosities(self, readInput=True):
        """Calcultes the frequency integrated luminosities of all radiation sources.

        
        Parameters
        ----------
        
        readInput - bool, optional
                    If true the input files of the radiation sources are read and the the total luminosities
                    are calculated from them. If readInput is set to False, the luminosities are calculated
                    by semi-analytic spectra.
        
        Returns
        -------
        
        Returns a dictionary with the following keys
            * lnu_star    : Luminosity of the discrete stars
            * lnu_accdisk : Luminosity of the accretion disk
            * lnu_spot    : Luminosity of the hot spot / boundary layer on the stellar surface
        """
        

        # If readIpnut
        if readInput:
            self.readStarsinp()

            # Note the negative sign in dnu is there because the frequency array is ordered in wavelength not in frequency
            dnu = -(self.grid.freq[1:] - self.grid.freq[:-1])

            res = {}
            # Calculate the stellar luminosity
            res['lnu_star'] = np.zeros(self.nstar, dtype=float)

            for istar in range(self.nstar):
                res['lnu_star'][istar] = 0.5 * ((self.fnustar[1:,istar] + self.fnustar[:-1,istar]) * dnu).sum() * 4.*np.pi*pc**2

            # Calculate the luminosity in the continuous stellar sources
            csDensFound = False
            csTempFound = False
            if os.path.exists('stellarsrc_density.inp'):
                self.readStellarsrcDensity(fname='stellarsrc_density.inp', binary=False)
                csDensFound = True
            if os.path.exists('stellarsrc_density.binp'):
                self.readStellarsrcDensity(fname='stellarsrc_density.binp', binary=True)
                csDensFound = True
            if os.path.exists('stellarsrc_templates.inp'):
                self.readStellarsrcTemplates()
                csTempFound = True
            if (csDensFound)&(csTempFound):
                vol = self.grid.getCellVolume()
                
                factor = 4. * np.pi * np.pi / 4./ np.pi 
                dnu = abs(self.grid.freq[1:] - self.grid.freq[:-1])
                lum = 0. 
                for itemp in range(self.csntemplate):
                    for ix in range(self.grid.nx):
                        for iy in range(self.grid.ny):
                            for iz in range(self.grid.nz):
                                if self.csdens[ix,iy,iz,itemp]>0.:
                                    expterm = (hh*self.grid.freq/kk/(-self.cststar[ix])).clip(-600,600)
                                    bb = 2.*hh*self.grid.freq**3/cc**2/(np.exp(expterm)-1.)
                                    dum = bb * np.pi * vol[ix,iy,iz] * 4. * np.pi * self.csdens[ix,iy,iz,itemp]
                                    lum = lum + ( (dum[1:] + dum[:-1])*0.5*dnu ).sum() 


                res['lnu_accdisk'] = lum 

            else:
                res['lnu_spot'] = 0.
                res['lnu_accdisk'] = 0.
        else:
            
            # Note the negative sign in dnu is there because the frequency array is ordered in wavelength not in frequency
            dnu = -(self.grid.freq[1:] - self.grid.freq[:-1])

            res = {}
            # Calculate the stellar luminosity
            res['lnu_star'] = np.zeros(self.nstar, dtype=float)

            for istar in range(self.nstar):
                res['lnu_star'][istar] = 4.*np.pi*self.rstar[istar]**2*ss*self.tstar[istar]**4

            if self.accrate==0.:
                print 'Viscsous accretion is switched off'
                res['lnu_spot'] = 0.
                res['lnu_accdisk'] = 0.
            if not self.incl_accretion:
                print 'Viscsous accretion is switched off'
                res['lnu_spot'] = 0.
                res['lnu_accdisk'] = 0.
            else:
                # Calculate the spot luminosity
                res['lnu_spot'] = 0.5 * gg * self.mstar[0] * self.accrate / self.rstar[0]#4.*np.pi*self.rstar[istar]**2*starsurffac*ss*self.tspot**4
                # Calculate the accretion disk luminosity
                res['lnu_accdisk'] = 0.5 * gg * self.mstar[0] * self.accrate / self.rstar[0]
                #0.5 * ((self.fnuaccdisk[1:] + self.fnuaccdisk[:-1]) * dnu).sum() * 2. * np.pi *pc**2 


        return res
# --------------------------------------------------------------------------------------------------
    def getAccdiskSpectra(self, ppar=None, grid=None, incl=0.):
        """Calculates the emergent spectra of an optically thick accretion disk at face-on orientation (incl=0deg).
        
        Parameters
        ----------

        ppar : dictionary
            Dictionary containing all input parameters keys should include
                * mstar   : stellar mass
                * rstar   : stellar radius
                * accrate : accretion rate

                NOTE, that for the calculation of the effective disk temperature only the first
                    star is used if more than one values are given in mstar and rstar. 
        
        incl : float, optional
                Inclination angle in degrees at which the spectrum be calculated  (default - 0deg)
        
        grid : radmc3dGrid, optional 
            An instance of a radmc3dGrid class containing the spatial and wavelength grid

        """

        if ppar.has_key('accrate'):
            if ppar['accrate']>0.:
                self.getAccdiskTemperature(ppar=ppar, grid=grid)
                fnuaccdisk = np.zeros([self.grid.nx, self.grid.nwav], dtype=float)
                for ix in range(self.grid.nx):
                    dum = hh*self.grid.freq/kk/self.tacc[ix]
                    dum = dum.clip(-600., 600.)
                    bb = 2.*hh*self.grid.freq**3/cc**2 / (np.exp(np.float64(dum))-1.0)
                    fnuaccdisk[ix,:] =  bb * np.pi*(self.grid.xi[ix+1]**2 - self.grid.xi[ix]**2)  / pc**2 
                self.fnuaccdisk = fnuaccdisk.sum(0) 
            else:
                self.fnuaccdisk = np.zeros([self.grid.nwav], dtype=float)
        else:
            self.fnuaccdisk = np.zeros([self.grid.nwav], dtype=float)

# --------------------------------------------------------------------------------------------------
    def getAccdiskStellarTemplates(self, ppar=None, grid=None):
        """Calculates the stellar template for continuous starlike sources for modeling a viscous accretion disk.
       

        Parameters
        ----------

        ppar : dictionary
            Dictionary containing all input parameters keys should include:
                * mstar   : stellar mass
                * rstar   : stellar radius
                * accrate : accretion rate

                NOTE, that for the calculation of the effective disk temperature only the first
                    star is used if more than one values are given in mstar and rstar. 
        
        incl : float, optional
                Inclination angle in degrees at which the spectrum be calculated  (default - 0deg)
        
        grid : radmc3dGrid, optional 
            An instance of a radmc3dGrid class containing the spatial and wavelength grid

        """

       
        if self.incl_accretion:
            self.getAccdiskTemperature(ppar=ppar, grid=grid)
            self.cstemptype  = 1
            self.cststar     = -self.tacc
            self.csmstar     = self.cststar * 0. + 1.
            self.csrstar     = self.cststar * 0. + 1.
            self.csntemplate = self.grid.nx
        else:

            self.cstemptype  = 1
            self.cststar     = np.zeros(self.grid.nx)
            self.csmstar     = self.cststar * 0. 
            self.csrstar     = self.cststar * 0. 
            self.csntemplate = self.grid.nx
# --------------------------------------------------------------------------------------------------
    def getAccdiskStellarDensity(self, ppar=None, grid=None):
        """Calculates the stellar density for continuous starlike sources for modeling a viscous accretion disk.
        
        Parameters
        ----------

        ppar : dictionary
            Dictionary containing all input parameters keys should include:
                * mstar   : stellar mass
                * rstar   : stellar radius
                * accrate : accretion rate

                NOTE, that for the calculation of the effective disk temperature only the first
                    star is used if more than one values are given in mstar and rstar. 
        
        incl : float, optional
                Inclination angle in degrees at which the spectrum be calculated  (default - 0deg)
        
        grid : radmc3dGrid, optional 
            An instance of a radmc3dGrid class containing the spatial and wavelength grid
        
        """

        if grid!=None:
            self.grid = grid

        self.csdens = np.zeros([self.grid.nx, self.grid.ny, self.grid.nz, self.grid.nx], dtype=float)
        vol         = self.grid.getCellVolume()

        if self.incl_accretion:
            if self.grid.crd_sys != 'sph':
                print 'ERROR'
                print ' Viscous accretion is currently available only in spherical coordinate system'
                return False
            else:
                if abs(self.grid.yi[self.grid.ny]-np.pi)<1e-8:

                    for ix in range(self.grid.nx):

                        dA = 2.0 * (self.grid.xi[ix+1]**2 - self.grid.xi[ix]**2) * np.pi * (self.grid.zi[1:] - self.grid.zi[:-1]) / (2.*np.pi)
                        dV = vol[ix,self.grid.ny/2-1,:] +  vol[ix,self.grid.ny/2,:] 
                        self.csdens[ix,self.grid.ny/2-1,:,ix] = dA/dV/(4.*np.pi)
                        self.csdens[ix,self.grid.ny/2,:,ix]   = dA/dV/(4.*np.pi)

                elif abs(self.grid.yi[self.grid.ny]-np.pi/2.)<1e-8:
                    for ix in range(self.grid.nx):
                        dA = 2.0 * (self.grid.xi[ix+1]**2 - self.grid.xi[ix]**2) * np.pi * (self.grid.zi[1:] - self.grid.zi[:-1]) / (2.*np.pi)
                        dV = vol[ix,self.grid.ny-1,:] * 2.
                        self.csdens[ix,self.grid.ny-1,:,ix] = dA/dV/(4.*np.pi)
    
        return True

# --------------------------------------------------------------------------------------------------
    def readStellarsrcTemplates(self,fname='stellarsrc_templates.inp'):
        """Reads the stellar template of a continuous starlike source.

        Parameters
        ----------

        fname : str, optional 
                Name of the file from which the stellar templates will be read. If omitted the default
                'stellarsrc_templates.inp' will be used.
        """

        rfile = -1
        try :
            rfile = open(fname, 'r')
        except:
            print 'Error!' 
            print fname+' was not found!'
        
        if (rfile!=(-1)):
            
            self.grid = readGrid()
            
            hdr = np.fromfile(fname, count=3, sep="\n", dtype=int)
            
            if (self.grid.nwav!=hdr[2]):
                print 'Error!'
                print 'Number of grid points in wavelength_micron.inp is not equal to that in '+fname
            else:
                self.csntemplate = hdr[1] 
                dum = ''

                # Read the header
                for i in range(3):
                    dum = rfile.readline()

                # Read the frequency grid
                for i in range(hdr[2]):
                    dummy = float(rfile.readline())
                    if abs((cc/dummy*1e4 - self.grid.wav[i]) / self.grid.wav[i])>1e-4:
                        print 'ERROR'
                        print 'The wavelength grid in wavelength_micron.inp is different than that in '+fname
                        print cc/dummy*1e4, self.grid.wav[i]

                dum = rfile.readline()
                if float(dum)>0:
                    self.cstemp  = np.zeros([self.csntemplate,self.grid.nwav], dtype=float)
                    self.cststar = []
                    self.csrstar = []
                    self.csmstar = []
                    
                    self.cstemp[0,0] = float(dum)
                    for inu in range(1,self.grid.nwav):
                        dum = rfile.readline()
                        self.cstemp[0,inu] = float(dum)

                    for itemp in range(1, self.csntemplate):
                        for inu in range(self.grid.nwav):
                            dum = rfile.readline()
                            self.cstemp[itemp,inu] = float(dum)
                    
                else:
                    self.cstemp  = []
                    self.cststar = np.zeros(self.csntemplate, dtype=float)
                    self.csrstar = np.zeros(self.csntemplate, dtype=float)
                    self.csmstar = np.zeros(self.csntemplate, dtype=float)
                    
                    self.cststar[0] = float(dum)
                    dum = rfile.readline()
                    self.csrstar[0] = float(dum)
                    dum = rfile.readline()
                    self.csmstar[0] = float(dum)

                    for i in range(1,self.csntemplate):
                        dum = rfile.readline()
                        self.cststar[i] = float(dum)
                        dum = rfile.readline()
                        self.csrstar[i] = float(dum)
                        dum = rfile.readline()
                        self.csmstar[i] = float(dum)



                #data = np.fromfile(fname, count=-1, sep="\n", dtype=np.float64)
                #print data[3:].shape
                #print hdr[2]*self.grid.nz*self.grid.ny*self.grid.nx
                #print self.grid.nx, self.grid.ny, self.grid.nz, hdr[2]
                #exit()
                #data = reshape(data[3:], [hdr[2],self.grid.nz,self.grid.ny,self.grid.nx])
                ## We need to change the axis orders as Numpy always reads  in C-order while RADMC-3D
                ## uses Fortran-order
                #data = swapaxes(data,0,3)
                #data = swapaxes(data,1,2)
        
        else:
            data = -1

        if rfile!=(-1):
            rfile.close()
    


# --------------------------------------------------------------------------------------------------
    def writeStellarsrcTemplates(self,fname='stellarsrc_templates.inp'):
        """Writes the stellar template of a continuous starlike source.

        Parameters
        ----------
        
        fname : str, optional
                Name of the file into which the stellar templates will be written. If omitted the default
                'stellarsrc_templates.inp' will be used.
        """
        # First check if we'd need to write anything at al

        if len(self.cststar)==0:
            if len(self.cstemp)==0:
                if os.path.exists('stellarsrc_templates.inp'):
                    print 'The continuous starlike source seems to be inactive (zero input luminosity)'
                    print ' still stellarsrc_templates.inp file is present in the current working directory.'
                    dum = raw_input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0]=='y':
                        os.system('rm stellarsrc_templates.inp')
                    return
                return
            else:
                if abs(self.cstemp).max()==0.:
                    if os.path.exists('stellarsrc_templates.inp'):
                        print 'The continuous starlike source seems to be inactive (zero input luminosity)'
                        print ' still stellarsrc_templates.inp file is present in the current working directory.'
                        dum = raw_input('Can it be deleted (yes/no)')
                        if dum.strip().lower()[0]=='y':
                            os.system('rm stellarsrc_templates.inp')
                        return
                    return

        else:
            if abs(self.cststar).max()==0.:
                if os.path.exists('stellarsrc_templates.inp'):
                    print 'The continuous starlike source seems to be inactive (zero input luminosity)'
                    print ' still stellarsrc_templates.inp file is present in the current working directory.'
                    dum = raw_input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0]=='y':
                        os.system('rm stellarsrc_templates.inp')
                    return
                return

        
        print 'Writing '+fname
        wfile = open(fname, 'w')
        # Format number
        wfile.write("%d\n"%1)
        # Nr of templates
        wfile.write("%d\n"%self.csntemplate)
        # Nr of wavelengths
        wfile.write("%d\n"%self.grid.nwav)
        # Write the wavelength grid (in micron!)
        for ilam in range(self.grid.nwav):
            wfile.write("%.9e\n"%self.grid.freq[ilam])

        # Now write the templates
        # Similar to the discrete stellar imput if the first number is negative it means that 
        #  instead of a full frequency-dependent spectrum only the blackbody temperature is given.
        #  Thus I'd only need to give the temperatures as negative numbers and radmc-3d will take care
        #  of calculating the Planck-function. This will save some harddisk space.

        if self.cstemptype==1:
            for itemp in range(self.csntemplate):
                # Effective temperature
                if self.cststar[itemp]>0:
                    wfile.write("%.9e\n"%(-self.cststar[itemp]))
                else:
                    wfile.write("%.9e\n"%self.cststar[itemp])
                # "Stellar radius"
                wfile.write("%.9e\n"%self.csrstar[itemp])
                # "Stellar mass"
                wfile.write("%.9e\n"%self.csmstar[itemp])
        elif self.cstemptype==2:
            for itemp in range(self.csntemplate):
                for inu in range(self.grid.nwav):
                    wfile.write("%.9e\n"%self.cstemp[itemp,inu])
        else:
            print 'ERROR'
            print 'Unknown cstemptype for the continuous starlike source'
            print self.cstemptype
            wfile.close()
            return False

        wfile.close()

# --------------------------------------------------------------------------------------------------
    def readStellarsrcDensity(self, fname='', binary=False):
        """Reads the stellar density of a continuous starlike source.

        Parameters
        ----------
        
        fname : str, optional
                Name of the file from which the stellar templates will be read. If omitted the default
                'stellarsrc_templates.inp' will be used.

        binary : bool, optional
                If True the file should contain a C-style binary stream, if False the file should be 
                written as formatted ASCII
        """
        self.grid = readGrid()

        if binary:
            hdr = np.fromfile(fname, count=4, dtype=int)
            
            if hdr[2]!=(self.grid.nx*self.grid.ny*self.grid.nz):
                print ' ERROR'
                print ' Number of grid points in '+fname+' is different from that in amr_grid.inp'
                print npoints
                print hdr[2]
                return

            if hdr[1]==8:
                data = np.fromfile(fname, count=-1, dtype=np.float64)
            elif hdr[1]==4:
                data = np.fromfile(fname, count=-1, dtype=float)
            else:
                print 'ERROR'
                print 'Unknown datatype in '+fname
                return
            
            data = np.reshape(data[4:], [hdr[3],self.grid.nz,self.grid.ny,self.grid.nx])
            data = np.swapaxes(data,0,3)
            data = np.swapaxes(data,1,2)
        else:
            rfile = -1
            try :
                rfile = open(fname, 'r')
            except:
                print 'Error!' 
                print fname+' was not found!'
            
            if (rfile!=(-1)):

                hdr = np.fromfile(fname, count=3, sep="\n", dtype=int)
                
                if ((self.grid.nx * self.grid.ny * self.grid.nz)!=hdr[1]):
                    print 'Error!'
                    print 'Number of grid points in amr_grid.inp is not equal to that in '+fname
                else:

                    data = np.fromfile(fname, count=-1, sep="\n", dtype=np.float64)
                    data = np.reshape(data[3:], [hdr[2],self.grid.nz,self.grid.ny,self.grid.nx])
                    # We need to change the axis orders as Numpy always reads  in C-order while RADMC-3D
                    # uses Fortran-order
                    data = np.swapaxes(data,0,3)
                    data = np.swapaxes(data,1,2)
            
            else:
                data = -1

            if rfile!=(-1):
                rfile.close()
        
        self.csdens = data
# --------------------------------------------------------------------------------------------------
    def writeStellarsrcDensity(self, fname='', binary=False):
        """Writes the stellar density of a continuous starlike source.

        Parameters 
        ----------
        
        fname : str, optional
                Name of the file into which the stellar templates will be written. If omitted the default
                'stellarsrc_templates.inp' will be used.

        binary : bool, optional
                If True the output will be written in a C-style binary stream, if False the output will be 
                formatted ASCII
        """
        
        # First check if we'd need to write anything at al 
        # Both the stellar temperature and the stellar template spectrum arrays are empty 
        # So the continuous starlike source needs to be deactivated 
        # Let's check if there are any files, and if so ask them to be deleted
        if len(self.cststar)==0:
            if len(self.cstemp)==0:
                if (os.path.exists('stellarsrc_density.inp'))|(os.path.exists('stellarsrc_density.binp')):
                    print 'The continuous starlike source seems to be inactive (zero input luminosity)'
                    print ' still stellarsrc_density.inp/stellarsrc_density.binp file is present in the current '
                    print ' working directory.'
                    dum = raw_input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0]=='y':
                        os.system('rm stellarsrc_density.inp')
                        os.system('rm stellarsrc_density.binp')
                    return
                return
            else:
                if abs(self.cstemp).max()==0.:
                    if (os.path.exists('stellarsrc_density.inp'))|(os.path.exists('stellarsrc_density.binp')):
                        print 'The continuous starlike source seems to be inactive (zero input luminosity)'
                        print ' still stellarsrc_density.inp/stellarsrc_density.binp file is present in the current '
                        print ' working directory.'
                        dum = raw_input('Can it be deleted (yes/no)')
                        if dum.strip().lower()[0]=='y':
                            os.system('rm stellarsrc_density.inp')
                            os.system('rm stellarsrc_density.binp')
                        return
                    return

        else:
            if abs(self.cststar).max()==0.:
                if (os.path.exists('stellarsrc_density.inp'))|(os.path.exists('stellarsrc_density.binp')):
                    print 'The continuous starlike source seems to be inactive (zero input luminosity)'
                    print ' still stellarsrc_density.inp/stellarsrc_density.binp file is present in the current '
                    print ' working directory.'
                    dum = raw_input('Can it be deleted (yes/no)')
                    if dum.strip().lower()[0]=='y':
                        os.system('rm stellarsrc_density.inp')
                        os.system('rm stellarsrc_density.binp')
                    return
                return


        if binary:
            if fname.strip()=='':
                fname = 'stellarsrc_density.binp'
            print 'Writing '+fname
            wfile = open(fname, 'w')
            hdr = np.array([1, 8, self.grid.nx*self.grid.ny*self.grid.nz, self.csntemplate], dtype=int)
            hdr.tofile(wfile)
            data = np.array(self.csdens)
            data = np.swapaxes(data,0,3)
            data = np.swapaxes(data,1,2)
            data.tofile(wfile)

        else:
            if fname.strip()=='':
                wfile = open('stellarsrc_density.inp', 'w')
            else:
                wfile = open(fname, 'w')
            hdr = np.array([1, self.grid.nx*self.grid.ny*self.grid.nz, self.csntemplate], dtype=int)
            hdr.tofile(wfile, sep=" ", format="%d\n")
            data = np.array(self.csdens)
            data = np.swapaxes(data,0,3)
            data = np.swapaxes(data,1,2)
            data.tofile(wfile, sep=" ", format="%.9e\n")
        
        wfile.close()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class radmc3dDustOpac(object):
    """
    Class to handle dust opacities.


    Attributes
    ----------

    wav     : list
                Each element of the list contains an ndarray with the wavelength grid
    
    freq    : list
                Each element of the list contains an ndarray with the frequency grid
    
    nwav    : list
                Each element of the list contains an integer with the number of wavelengths
    
    kabs    : list
                Each element of the list contains an ndarray with the absorption coefficient per unit mass 
    
    ksca    : list
                Each element of the list contains an ndarray with the scattering coefficient per unit mass
    
    phase_g : list
                Each element of the list contains an ndarray with the hase function
    
    ext     : list
                Each element of the list contains a string wht the file name extension of the duskappa_ext.Kappa file
    
    therm   : list
                Each element of the list contains a bool, if it is set to False the dust grains are quantum-heated (default: True)
    
    idust   : lisintt
                Each element of the list contains an integer with the index of the dust species in the dust density distribution array

    scatmat : list
                Each element is a boolean indicating whether the dust opacity table includes (True) the full scattering matrix or not (False)

    nang    : list
                Each element is a string, containing the number of scattering angles in the scattering matrix if its given

    scatang : list
                Each element is a numpy ndarray containing the scattering angles in the scattering matrix if its given

    z11     : list
                Each element is a numpy ndarray containing the (1,1) element of the scattering angles in the scattering matrix if its given

    z12     : list
                Each element is a numpy ndarray containing the (1,2) element of the scattering angles in the scattering matrix if its given

    z22     : list
                Each element is a numpy ndarray containing the (2,2) element of the scattering angles in the scattering matrix if its given
    
    z33     : list
                Each element is a numpy ndarray containing the (3,3) element of the scattering angles in the scattering matrix if its given
    
    z34     : list
                Each element is a numpy ndarray containing the (3,4) element of the scattering angles in the scattering matrix if its given
    
    z44     : list
                Each element is a numpy ndarray containing the (4,4) element of the scattering angles in the scattering matrix if its given

    """
# --------------------------------------------------------------------------------------------------
    def __init__(self):

        self.wav      = []
        self.freq     = []
        self.nwav     = []
        self.nfreq    = []
        self.kabs     = []
        self.ksca     = []
        self.phase_g  = []
        self.ext      = []
        self.idust    = []
        self.therm    = []
        self.scatmat  = []
        self.z11      = [] 
        self.z12      = [] 
        self.z22      = [] 
        self.z33      = [] 
        self.z34      = [] 
        self.z44      = [] 
        self.scatang  = []
        self.nang     = []

# --------------------------------------------------------------------------------------------------
    def readOpac(self, ext=[''], idust=None, scatmat=None, old=False):
        """Reads the dust opacity files.

        Parameters
        ----------
        
        ext  : list
                File name extension (file names should look like 'dustkappa_ext.inp')
        
        idust: list
                Indices of the dust species in the master opacity file (dustopac.inp') - starts at 0 
        
        scatmat: list
                If specified, its elements should be booleans indicating whether the opacity file 
                contains also the full scattering matrix (True) or only dust opacities (False)
        
        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """
       
        # Check the input keywords and if single strings are given convert them to lists
        # This assumes, though, that there is a single dust opacity file or dust species, though!!
        if (type(ext).__name__=='str'):  ext = [ext]
        if scatmat!=None:
            if (type(scatmat).__name__=='str'):  scatmat = [scatmat]
        else:
            # If the scatmat keyword is not given (i.e. if it is None) then assume that 
            # it is False for all dust species
            scatmat = []
            if idust!=None:
                for i in range(len(idust)):
                    scatmat.append(False)
            else:
                for i in range(len(ext)):
                    scatmat.append(False)
        
        if idust!=None:
            if (type(idust).__name__=='int'):  idust = [idust]

        if (len(ext)==1)&(ext[0]!=''):
            if idust!=None:
                print 'ERROR'
                print 'Either idust or ext should be specified, but not both'
                print idust
                print ext
                return [-1]
        
        # Read the master dust opacity file to get the dust indices and dustkappa file name extensions
        mopac = self.readMasterOpac()

        # Find the file name extensions in the master opacity file if idust is specified instead of ext
        if idust:
            ext = []
            for ispec in idust:
                if (ispec+1)>len(mopac['ext']):    
                    print 'ERROR'
                    print 'No dust species found at index ', ispec
                    return [-1]
                else:
                    ext.append(mopac['ext'][ispec])

        # If only the extension is specified look for the master opacity file and find the index of this dust species
        #  or set the index to -1 if no such dust species is present in the master opacity file
        else:
            idust = []
            for iext in ext:
                try:
                    dum2 = mopac['ext'].index(iext)
                except:
                    dum2 = -1
                idust.append(dum2)
        
        # Now read all dust opacities
        for i in range(len(ext)):
            if scatmat[i]:
                try:
                    rfile = open('dustkapscatmat_'+ext[i]+'.inp', 'r')
                except:
                    print 'ERROR'
                    print ' No dustkapscatmat_'+ext[i]+'.inp file was found'
                    return -1
                
                print 'Reading dustkapscatmat_'+ext[i]+'.inp ....'

                self.ext.append(ext[i])
                
                # Read the header/comment field
                dum = rfile.readline()
                while dum.strip()[0]=='#':
                    dum = rfile.readline()


                #for j in range(6):
                    #dum = rfile.readline()

                # Read the file format
                iformat = int(dum)
                #iformat = int(rfile.readline())
                if iformat!=1:
                    print 'ERROR'
                    print 'Format number of the file dustkapscatmat_'+ext[i]+'.inp (iformat='+("%d"%iformat)+') is unkown'
                    return [-1]

                # Read the number of wavelengths in the file
                dum = int(rfile.readline())
                self.nwav.append(dum)
                self.nfreq.append(dum)
                self.idust.append(idust[i])
                idu = len(self.nwav)-1
                
                # Read the scattering angular grid
                self.nang.append(int(rfile.readline()))
                wav     = np.zeros(self.nwav[idu], dtype = np.float64)
                kabs    = np.zeros(self.nwav[idu], dtype = np.float64)
                ksca    = np.zeros(self.nwav[idu], dtype = np.float64)
                phase_g = np.zeros(self.nwav[idu], dtype = np.float64)
                scatang = np.zeros(self.nang[idu], dtype=np.float64)
                z11     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
                z12     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
                z22     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
                z33     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
                z34     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
                z44     = np.zeros([self.nwav[idu], self.nang[idu]], dtype=np.float64) 
            
                print 'Reading the opacities..'
                dum = rfile.readline()
                for ilam in range(self.nwav[idu]):
                    dum      = rfile.readline().split()
                    wav[ilam]  = float(dum[0])
                    kabs[ilam] = float(dum[1])
                    ksca[ilam] = float(dum[2])
                    phase_g[ilam] = float(dum[3])

                print 'Reading the angular grid..'
                dum = rfile.readline()
                for iang in range(self.nang[idu]):
                    dum        = rfile.readline()
                    scatang[iang] = float(dum)

                print 'Reading the scattering matrix..'
                for ilam in range(self.nwav[idu]):
                    dum = rfile.readline()
                    for iang in range(self.nang[idu]):
                        dum        = rfile.readline().split()
                        z11[ilam,iang] = float(dum[0])
                        z12[ilam,iang] = float(dum[1])
                        z22[ilam,iang] = float(dum[2])
                        z33[ilam,iang] = float(dum[3])
                        z34[ilam,iang] = float(dum[4])
                        z44[ilam,iang] = float(dum[5])
                
                self.wav.append(wav)
                self.freq.append(cc/wav*1e4)
                self.kabs.append(kabs)
                self.ksca.append(ksca)
                self.phase_g.append(phase_g)
                self.scatang.append(scatang)
                self.z11.append(z11)
                self.z12.append(z12)
                self.z22.append(z22)
                self.z33.append(z33)
                self.z34.append(z34)
                self.z44.append(z44)
               
                rfile.close()
            else:
                if not old:
                    try:
                        rfile = open('dustkappa_'+ext[i]+'.inp', 'r')
                    except:
                        print 'ERROR'
                        print ' No dustkappa_'+ext[i]+'.inp file was found'
                        return -1

                    self.ext.append(ext[i])

                    # Read the file format
                    iformat = int(rfile.readline())
                    if (iformat<1)|(iformat>3):
                        print 'ERROR'
                        print 'Unknown file format in the dust opacity file'
                        rfile.close()
                        return -1


                    # Read the number of wavelengths in the file
                    dum = rfile.readline()
                    self.nwav.append(int(dum))
                    self.nfreq.append(int(dum))
                    self.idust.append(idust[i])
                    idu = len(self.nwav)-1

                    # If only the absorption coefficients are specified
                    if iformat==1:
                        wav = np.zeros(self.nwav[idu], dtype=np.float64)
                        kabs = np.zeros(self.nwav[idu], dtype=np.float64)
                        for ilam in range(self.nwav[idu]):
                            dum = rfile.readline().split()
                            wav[ilam] = float(dum[0])
                            kabs[ilam] = float(dum[1])
                        self.wav.append(wav)
                        self.freq.append(cc/wav*1e4)
                        self.kabs.append(kabs)
                        self.ksca.append([-1])
                        self.phase_g.append([-1])
                    # If the absorption and scattering coefficients are specified
                    elif iformat==2:
                        wav = np.zeros(self.nwav[idu], dtype=np.float64)
                        kabs = np.zeros(self.nwav[idu], dtype=np.float64)
                        ksca = np.zeros(self.nwav[idu], dtype=np.float64)
                        for ilam in range(self.nwav[idu]):
                            dum = rfile.readline().split()
                            wav[ilam] = float(dum[0])
                            kabs[ilam] = float(dum[1])
                            ksca[ilam] = float(dum[2]) 
                        self.wav.append(wav)
                        self.freq.append(cc/wav*1e4)
                        self.kabs.append(kabs)
                        self.ksca.append(ksca)
                        self.phase_g.append([-1])
                    
                    # If the absorption and scattering coefficients and also the scattering phase function are specified
                    elif iformat==3:
                        wav = np.zeros(self.nwav[idu], dtype=np.float64)
                        kabs = np.zeros(self.nwav[idu], dtype=np.float64)
                        ksca = np.zeros(self.nwav[idu], dtype=np.float64)
                        phase_g = np.zeros(self.nwav[idu], dtype=np.float64)
                        for ilam in range(self.nwav[idu]):
                            dum = rfile.readline().split()
                            wav[ilam] = float(dum[0])
                            kabs[ilam] = float(dum[1])
                            ksca[ilam] = float(dum[2])
                            phase_g[ilam] = float(dum[3]) 
                        
                        self.wav.append(wav)
                        self.freq.append(cc/wav*1e4)
                        self.kabs.append(kabs)
                        self.ksca.append(ksca)
                        self.phase_g.append(phase_g)
               
                    rfile.close()
                else:
                    try:
                        rfile = open('dustopac_'+ext[i]+'.inp', 'r')
                    except:
                        print 'ERROR'
                        print ' No dustopac_'+ext[i]+'.inp file was found'
                        return -1
                 
                    freq = np.fromfile('frequency.inp', count=-1, sep="\n", dtype=float)
                    nfreq = int(freq[0])
                    freq = freq[1:]

                    self.ext.append(ext[i])
                    dum   = rfile.readline().split()
                    if int(dum[0])!=nfreq:
                        print 'ERROR'
                        print 'dustopac_'+ext[i]+'.inp contains a different number of frequencies than frequency.inp'
                        return

                    wav     = cc/freq*1e4
                    kabs    = np.zeros(nfreq, dtype = float)
                    ksca    = np.zeros(nfreq, dtype = float)

                    dum     = rfile.readline()
                    for ilam in range(nfreq):
                        kabs[ilam] = float(rfile.readline())
                    dum     = rfile.readline()
                    for ilam in range(nfreq):
                        ksca[ilam] = float(rfile.readline())
                    
                    rfile.close()

                    self.wav.append(wav[::-1])
                    self.freq.append(freq[::-1])
                    self.kabs.append(kabs[::-1])
                    self.ksca.append(ksca[::-1])
                    self.phase_g.append([-1])


                    
        return 0 
#--------------------------------------------------------------------------------------------------------------------
    def makeOpac(self, ppar=None, wav=None, old=False):
        """Createst the dust opacities using a Mie code distributed with RADMC-3D.
        
        Parameters
        ----------
        
        ppar  : dictionary
                Parameters of the simulations
        
        wav   : ndarray, optional
                Wavelength grid on which the mass absorption coefficients should be calculated
    
        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """

    #
    # Create the wavelength grid if it is not specified
    #
        if wav==None:
            grid = radmc3dGrid()
            grid.makeWavelengthGrid(ppar=ppar)
            wav = grid.wav

    #
    # Do we need to mix the opacities?
    #
        if type(ppar['lnk_fname']).__name__=='str':
            ppar['lnk_fname'] = [ppar['lnk_fname']]

        if len(ppar['lnk_fname'])>1:
            #if old:
                #print 'ERROR'
                #print 'Making old (RADMC) style opacities is not finished for more than one dust species'
                #exit() 

            ext = []
            for idust in range(len(ppar['lnk_fname'])):
                
                # makedust needs the lnk file to be sorted in wavelength so create a dummy file 
                # which contains the sorted optical constants 
                try:
                    rfile = open(ppar['lnk_fname'][idust], 'r')
                except:
                    print 'ERROR'
                    print ppar['lnk_fname'][idust] + ' could not be opened'
                    return

                try:
                    w = []
                    n = []
                    k = []
                    dum = rfile.readline()
                    while len(dum)>0:
                        dum = dum.split()
                        w.append(dum[0])
                        n.append(dum[1])
                        k.append(dum[2])
                        dum = rfile.readline()

                    rfile.close()
                except:
                    print 'ERROR'
                    print ppar['lnk_fname'][idust] + ' could not be read'
                    return

                w = np.array(w, dtype=float)
                n = np.array(n, dtype=float)
                k = np.array(k, dtype=float)

                if float(w[0])>float(w[w.shape[0]-1]):
                    w = w[::-1]
                    n = n[::-1]
                    k = k[::-1]

                #Write out the dummy file containing the sorted optical constants
                wfile = open('opt_const.dat', 'w')
                for iwav in range(w.shape[0]):
                    wfile.write("%s %s %s \n"%(w[iwav], n[iwav], k[iwav]))
                wfile.close()

                # Run makedust
                self.runMakedust(freq=cc/wav*1e4, gmin=ppar['gsmin'], gmax=ppar['gsmax'], ngs=ppar['ngs'], \
                        lnk_fname='opt_const.dat', gdens=ppar['gdens'][idust])

                # Change the name of makedust's output
                for igs in range(ppar['ngs']):
                    dum = sp.Popen('mv dustkappa_'+str(igs+1)+'.inp dustkappa_idust_'+str(idust+1)+'_igsize_'+str(igs+1)+'.inp', shell=True).wait()
                    ext.append('idust_'+str(idust+1)+'_igsize_'+str(igs+1))

                os.remove('opt_const.dat')

            # Mix the opacity of different dust species for a given grain size if mixing is requested
            if ppar.has_key('mixabun'):
                if len(ppar['mixabun'])==len(ppar['lnk_fname']):
                    ext = []
                    for igs in range(ppar['ngs']):
                        mixnames = ['dustkappa_igsize_'+str(igs+1)+'.inp']
                        mixspecs = [['dustkappa_idust_'+str(idust+1)+'_igsize_'+str(igs+1)+'.inp' for idust in range(len(ppar['lnk_fname']))]]
                        self.mixOpac(mixnames=mixnames, mixspecs=mixspecs, mixabun=[ppar['mixabun']])
                    
                        ext.append('igsize_'+str(igs+1))
                else:
                    print 'ERROR'
                    print ' mixabun and lnk_fname should have the same number of elements.'
                    print ' To disable mixing either set mixabun to an empty list ([]) or comment it out in the problem_params.inp file'
                    return
            
            therm = [True for i in range(len(ext))]
            self.writeMasterOpac(ext=ext, therm=therm, scattering_mode_max=ppar['scattering_mode_max'], old=old)
           
            if old:
                self.makeopacRadmc2D(ext=ext)

        else:
            # makedust needs the lnk file to be sorted in wavelength so create a dummy file 
            # which contains the sorted optical constants 
            try:
                rfile = open(ppar['lnk_fname'][0], 'r')
            except:
                print 'ERROR'
                print ppar['lnk_fname'][0] + ' could not be opened'
                return
            try:
                w = []
                n = []
                k = []
                dum = rfile.readline()
                while len(dum)>0:
                    dum = dum.split()
                    w.append(dum[0])
                    n.append(dum[1])
                    k.append(dum[2])
                    dum = rfile.readline()

                rfile.close()
            except:
                print 'ERROR'
                print ppar['lnk_fname'][0] + ' could not be read'
                return

            w = np.array(w, dtype=float)
            n = np.array(n, dtype=float)
            k = np.array(k, dtype=float)

            if float(w[0])>float(w[w.shape[0]-1]):
                w = w[::-1]
                n = n[::-1]
                k = k[::-1]
            
            # Write out the dummy file containing the sorted optical constants
            wfile = open('opt_const.dat', 'w')
            for iwav in range(w.shape[0]):
                wfile.write("%s %s %s \n"%(w[iwav], n[iwav], k[iwav]))
            wfile.close()

            # Run makedust
            self.runMakedust(freq=cc/wav*1e4, gmin=ppar['gsmin'], gmax=ppar['gsmax'], ngs=ppar['ngs'], \
                    lnk_fname='opt_const.dat', gdens=ppar['gdens'][0])

            # Change the name of makedust's output
            ext = []
            therm = []
            for igs in range(ppar['ngs']):
                dum = sp.Popen('mv dustkappa_'+str(igs+1)+'.inp dustkappa_idust_1_igsize_'+str(igs+1)+'.inp', shell=True).wait()
                ext.append('idust_1_igsize_'+str(igs+1))
                therm.append(True)


#            # Change the name of makedust's output 
#            dum = sp.Popen('mv dustkappa_1.inp dustkappa_idust_1_igsize_1.inp', shell=True).wait()
#            os.remove('opt_const.dat')


            self.writeMasterOpac(ext=ext, therm=therm, scattering_mode_max=ppar['scattering_mode_max'], old=old)
            if old:
                self.makeopacRadmc2D(ext=ext)
       
        # Clean up and remove dust.inp and frequency.inp
        os.remove('dust.inp')
        if not old:
            os.remove('frequency.inp')
# --------------------------------------------------------------------------------------------------
    def mixOpac(self, ppar=None, mixnames=[], mixspecs=[], mixabun=[], writefile=True):
        """Mixes dust opacities.


        Parameters
        -----------
        ppar      : dictionary, optional
                    All parameters of the actual model setup.
        
        mixnames  : list, optional
                    Names of the files into which the mixed dust opacities will be written (not needed if writefile=False)
        
        mixspecs  : list, optional
                    Names of the files from which the dust opacities are read (not needed if readfile=False)
        
        mixabun   : list, optional
                    Abundances of different dust species
        
        writefile : bool
                    If False the mixed opacities will not be written out to files given in mixnames.  
           
        NOTE, either ppar or  mixname, mixspecs, and mixabun should be set. 

        """

        if writefile:
            if len(mixnames)==0:
                if ppar!=None:
                    mixnames = ppar['mixnames']
                else:
                    print 'ERROR'
                    print ' Neither ppar nor mixnames are set in mixOpac '
                    return

        if len(mixspecs)==0:
            if ppar!=None:
                mixspecs = ppar['mixspecs']
            else:
                print 'ERROR'
                print ' Neither ppar nor mixspecs are set in mixOpac '
                return
            
        if len(mixabun)==0:
            if ppar!=None:
                mixabun = ppar['mixabun']
            else:
                print 'ERROR'
                print ' Neither ppar nor mixabun are set in mixOpac '
                return

        mwav  = []
        mcabs = []
        mcsca = []
        for i in range(len(mixnames)):
            #
            # Read the dust opacities to be mixed for composite dust species #1
            #
            ocabs = []
            ocsca = []
            ogsym = []
            oform = 0
            for j in range(len(mixspecs[i])):
                try:
                    rfile=open(mixspecs[i][j], 'r')
                    form   = int(rfile.readline())
                    nwav   = int(rfile.readline())
                    dw     = np.zeros(nwav, dtype=float)
                    dcabs  = np.zeros(nwav, dtype=float)
                    dcsca  = np.zeros(nwav, dtype=float)
                    gsym   = np.zeros(nwav, dtype=float)
                    if form==1:
                        if ((oform==0)|(oform==1)):
                            oform = 1
                        else:
                            print ' '
                            print 'WARNING'
                            print ' You are trying to mix opacity tables with different formats'
                            print ' Some of the tables contain scattering coefficients while (format>=2) while other do not (format=1)'
                            print ' If you wish to continue mixing will only be done for the absorption and the output opacity table'
                            print ' will have a format number of 1.'
                            dum = raw_input('Do you wish to continue (1-yes, 0-no) ?')
                            if dum.strip()!='1':
                                return

                        for iwav in range(nwav):
                            dum = rfile.readline().split()
                            dw[iwav], dcabs[iwav] = float(dum[0]), float(dum[1])
                    if form==2:
                        if ((oform==0)|(oform==2)):
                            oform=2
                        else:
                            print ' '
                            print 'WARNING'
                            print ' You are trying to mix opacity tables with different formats'
                            print ' Some of the tables contain scattering coefficients while (format>=2) while other do not (format=1)'
                            print ' If you wish to continue mixing will only be done for the absorption and the output opacity table'
                            print ' will have a format number of 1.'
                            dum = raw_input('Do you wish to continue (1-yes, 0-no) ?')
                            if dum.strip()!='1':
                                return
                        for iwav in range(nwav):
                            dum = rfile.readline().split()
                            dw[iwav], dcabs[iwav], dcsca[iwav] = float(dum[0]), float(dum[1]), float(dum[2])
                    if form==3:
                        if ((oform==0)|(oform==3)):
                            oform=3
                        else:
                            print ' '
                            print 'WARNING'
                            print ' You are trying to mix opacity tables with different formats'
                            print ' Some of the tables contain scattering coefficients while (format>=2) while other do not (format=1)'
                            print ' If you wish to continue mixing will only be done for the absorption and the output opacity table'
                            print ' will have a format number of 1.'
                            dum = raw_input('Do you wish to continue (1-yes, 0-no) ?')
                            if dum.strip()!='1':
                                return
                        for iwav in range(nwav):
                            dum = rfile.readline().split()
                            dw[iwav], dcabs[iwav], dcsca[iwav], gsym[iwav] = float(dum[0]), float(dum[1]), float(dum[2]), float(dum[3])
                    if form>3:
                        print ' '
                        print ' ERROR'
                        print ' Unsupported dust opacity table format (format number: '+form+')'
                        print ' Currently only format number 1 and 2 are supported'
                        return
                    rfile.close()

                    if dw[1]<dw[0]:
                        print ' Dust opacity table seems to be sorted in frequency instead of wavelength'
                        print ' Reversing the arrays'
                        dw = dw[::-1]
                        dcabs = dcabs[::-1]
                        dcsca = dcsca[::-1]
                except:
                    print 'ERROR'
                    print mixspecs[i][j]+ ' could not be read'
                    return

                if j==0:
                    ocabs = np.array(dcabs) * mixabun[i][j]
                    ocsca = np.array(dcsca) * mixabun[i][j]
                    ogsym = np.array(gsym) * mixabun[i][j]
                    nwav0 = dw.shape[0]
                    owav  = np.array(dw)
                else:
                    #
                    # Interpolate dust opacities to the wavelength grid of the first dust species
                    #
                    ii = ( (owav>=dw[0])&(owav<=dw[nwav-1]) )
                    il = (owav<dw[0]) 
                    ih = (owav>dw[nwav-1])
                    dum = np.zeros(nwav0, dtype=float)
                    dum[ii] = 10.**np.interp(np.log10(owav[ii]), np.log10(dw), np.log10(dcabs))
                   
                    # Edwtrapolate the absorption coefficients using linear fit in log-log space (i.e. fitting a polinomial) for short wavelengths
                    der = np.log10(dcabs[1]/dcabs[0]) / np.log10(dw[1]/dw[0])
                    dum[il] = 10.**(np.log10(dcabs[0]) + np.log10(dw[0]/owav[il]))
                    
                    # Edwtrapolate the absorption coefficients using linear fit in log-log space (i.e. fitting a polinomial) for long wavelengths
                    der = np.log10(dcabs[nwav-1]/dcabs[nwav-2]) / np.log10(dw[nwav-1]/dw[nwav-2])
                    dum[ih] = 10.**(np.log10(dcabs[nwav-1]) + np.log10(owav[il]/dw[nwav-1]))
                 
                    ocabs = ocabs + np.array(dum) * mixabun[i][j]
                    
                    if oform==2:
                        # Do the inter-/extrapolation of for the scattering coefficients
                        dum = np.zeros(nwav0, dtype=float)
                        dum[ii] = 10.**np.interp(np.log10(owav[ii]), np.log10(dw), np.log10(dcsca))
                       
                        der = np.log10(dcsca[1]/dcsca[0]) / np.log10(dw[1]/dw[0])
                        dum[il] = 10.**(np.log10(dcsca[0]) + np.log10(dw[0]/owav[il]))
                        
                        der = np.log10(dcsca[nwav-1]/dcsca[nwav-2]) / np.log10(dw[nwav-1]/dw[nwav-2])
                        dum[ih] = 10.**(np.log10(dcsca[nwav-1]) + np.log10(owav[il]/dw[nwav-1]))
                       
                        ocsca = ocsca + array(dum) * mixabun[i][j]

                    if oform==3:
                        # Do the inter-/extrapolation of for the scattering phase function
                        dum = np.zeros(nwav0, dtype=float)
                        dum[ii] = 10.**np.interp(np.log10(owav[ii]), np.log10(dw), np.log10(gsym))
                       
                        der = np.log10(gsym[1]/gsym[0]) / np.log10(dw[1]/dw[0])
                        dum[il] = 10.**(np.log10(gsym[0]) + np.log10(dw[0]/owav[il]))
                        
                        der = np.log10(gsym[nwav-1]/gsym[nwav-2]) / np.log10(dw[nwav-1]/dw[nwav-2])
                        dum[ih] = 10.**(np.log10(gsym[nwav-1]) + np.log10(owav[il]/dw[nwav-1]))
                       
                        ogsym = ogsym + np.array(dum) * mixabun[i][j]


       
            #
            # Write out the mixed dust opacities
            #
            wfile = open(mixnames[i], 'w')
            wfile.write("%d\n"%oform) 
            wfile.write("%d\n"%owav.shape[0])
            if oform==1:
                for iwav in range(owav.shape[0]):
                    wfile.write("%.9e %.9e\n"%(owav[iwav], ocabs[iwav]))
            if oform==2:
                for iwav in range(owav.shape[0]):
                    wfile.write("%.9e %.9e %.9e\n"%(owav[iwav], ocabs[iwav], ocsca[iwav]))
            if oform==3:
                for iwav in range(owav.shape[0]):
                    wfile.write("%.9e %.9e %.9e %.9e\n"%(owav[iwav], ocabs[iwav], ocsca[iwav], ogsym[iwav]))

        return 
# --------------------------------------------------------------------------------------------------
    def  readMasterOpac(self):
        """Reads the master opacity file 'dustopac.inp'. 
        It reads the dustkappa filename extensions (dustkappa_ext.inp) corresponding to dust species indices

        Returns
        -------
        
        Returns a dictionary with the following keys:
            
            *ext   : list of dustkappa file name extensions
            
            *therm : a list of integers specifying whether the dust grain is thermal or quantum heated 
            (0 - thermal, 1 - quantum heated)
        """
        
        try: 
            rfile = open('dustopac.inp', 'r')
        except:
            print 'Error'
            print ' No dustopac.inp file was found'
            return -1

       
        # file format
        dum = rfile.readline()
        # nr of dust species
        ndust = int(rfile.readline().split()[0])
        # Comment line
        dum = rfile.readline()

        ext = []
        therm= []
        scatmat = []
        for idust in range(ndust):
            # Check if we have dust opacities also for the full scattering matrix
            dum = rfile.readline().split()
            if int(dum[0])==1:
                scatmat.append(False)
            elif int(dum[0])==10:
                scatmat.append(True)

            # Check if the dust grain is thermal or quantum heated
            dum = int(rfile.readline().split()[0])
            if dum==0:
                therm.append(True)
            else:
                therm.append(False)
            # Dustkappa filename extension
            dum = rfile.readline().split()[0]
            ext.append(dum)
            #Comment line
            dum = rfile.readline()
        rfile.close()

        return {'ext':ext, 'therm':therm, 'scatmat':scatmat}
# --------------------------------------------------------------------------------------------------
    def  writeMasterOpac(self, ext=None, therm=None, scattering_mode_max=1, old=False):
        """Writes the master opacity file 'dustopac.inp'. 

        Parameters
        ----------
        
        ext   : list
                List of dustkappa file name extensions
        
        therm : list
                List of integers specifying whether the dust grain is thermal or quantum heated
                (0-thermal, 1-quantum)
        
        old   : bool, optional
                If set to True the file format of the previous, 2D version of radmc will be used
        """

        print 'Writing dustopac.inp'
       
        if not ext:
            print 'ERROR'
            print 'No file name extension is specified. Without it dustopac.inp cannot be written'
            return -1
        else:
            if (type(ext).__name__=='str'):  ext = [ext]

        if therm:
            if (type(therm).__name__=='int'): therm = [therm]
            if (len(ext)!=len(therm)):
                print 'ERROR'
                print ' The number of dust species in ext and in therm are different'
                return -1
        else:
            # If therm is not specified it is assumed that all grains are thermal, no quantum heating

            therm = [True for i in range(len(ext))]

        wfile = open('dustopac.inp', 'w')

        # File format
        wfile.write('%-15s %s\n'%('2', 'Format number of this file'))
        # Number of dust species
        wfile.write('%-15s %s\n'%(str(len(ext)), 'Nr of dust species'))
        # Separator
        wfile.write('%s\n'%'============================================================================')

        if not old:
            for idust in range(len(ext)):
                # Dust opacity will be read from a file
                if scattering_mode_max<5:
                    wfile.write('%-15s %s\n'%('1', 'Way in which this dust species is read'))
                else:
                    wfile.write('%-15s %s\n'%('10', 'Way in which this dust species is read'))

                # Check if the dust grain is thermal or quantum heated
                if therm:
                    if therm[idust]:
                        wfile.write('%-15s %s\n'%('0', '0=Thermal grain, 1=Quantum heated'))
                else:
                    wfile.write('%-15s %s\n'%('1', '0=Thermal grain, 1=Quantum heated'))

                # Dustkappa filename extension
                wfile.write('%s %s %s\n'%(ext[idust], '    ', 'Extension of name of dustkappa_***.inp file'))
                # Separator
                wfile.write('%s\n'%'----------------------------------------------------------------------------')
        else:
            for idust in range(len(ext)):
                # Dust opacity will be read from a file
                wfile.write('%-15s %s\n'%('-1', 'Way in which this dust species is read (-1=File)'))

                # Check if the dust grain is thermal or quantum heated
                wfile.write('%-15s %s\n'%('0', '0=Thermal grain, 1=Quantum heated'))
                # Dustkappa filename extension
                wfile.write('%d %s %s\n'%((idust+1), '    ', 'Extension of name of dustopac_***.inp file'))
                # Separator
                wfile.write('%s\n'%'----------------------------------------------------------------------------')

        wfile.close()

# --------------------------------------------------------------------------------------------------
    def  makeopacRadmc2D(self, ext=None):
        """
        Creates dust opacities (dustopac_*.inp files) for the previous 2D version of radmc
        It takes the input dust opacity files and interpolates them onto the used frequency grid

        Parameters
        ----------

            ext : list
                  List of dustkappa file name extensions, i.e. the input file name has to be named
                  as dustkappa_ext[i].inp

        """
        cc = 29979245800.
        
        if type(ext).__name__=='str':
            ext = [ext]
      
        # 
        # Read the frequency.inp file
        #
        freq = np.fromfile('frequency.inp', count=-1, sep="\n", dtype=float)
        nfreq = int(freq[0])
        freq = freq[1:]
        freq = freq[::-1]
        wav  = cc/freq*1e4 
        o    = self.readOpac(ext=ext)

        #
        # Check if the frequency grid is ordered in frequency or in wavelength
        #
        worder = False
        if freq[-1]<freq[0]:
            worder = True

        for i in range(len(ext)):
            kabs = np.zeros(nfreq, dtype=float)
            ksca = np.zeros(nfreq, dtype=float)
            ish  = (wav<self.wav[i][0])
            ilo  = (wav>self.wav[i][-1])
            ii   = ((wav>=self.wav[i][0])&(wav<=self.wav[i][-1])) 
           
            # 
            # Do logarithmic interpolation for the overlapping wavelenght domain
            #
            kabs[ii] = 10.**np.interp(np.log10(wav[ii]), np.log10(self.wav[i]), np.log10(self.kabs[i]))
            if len(self.ksca[i])>1:
                ksca[ii] = 10.**np.interp(np.log10(wav[ii]), np.log10(self.wav[i]), np.log10(self.ksca[i]))

            # 
            # Do the long wavelength part
            #
            if ilo.__contains__(True):
                x1  = np.log10(self.wav[i][-1])
                x0  = np.log10(self.wav[i][-2])
            
                y1  = np.log10(self.kabs[i][-1])
                y0  = np.log10(self.kabs[i][-2])
                der = (y1-y0) / (x1-x0)
                kabs[ilo] = 10.**(y1 + der*(np.log10(wav[ilo])-x1))

                y1  = np.log10(self.ksca[i][-1])
                y0  = np.log10(self.ksca[i][-2])
                der = (y1-y0) / (x1-x0)
                ksca[ilo] = 10.**(y1 + der*(np.log10(wav[ilo])-x1))

            #
            # Do the shorter wavelength
            # 
            if ish.__contains__(True):
                kabs[ish] = self.kabs[0][0]
                ksca[ish] = self.ksca[0][0]

            #
            # Now write the results to file
            # 
            fname = 'dustopac_'+("%d"%(i+1))+'.inp'
            try :
                wfile = open(fname, 'w')
            except:
                print 'Error!' 
                print fname+' cannot be opened!'
                return 
            
            print 'Writing '+fname

            wfile = open(fname, 'w')
            wfile.write("%d 1\n"%nfreq)
            
            wfile.write(" \n")
            #
            # Reverse the order of kabs,ksca as they are ordered in frequency in radmc 
            #
            if worder:
                x = kabs[::-1]
            else:
                x = kabs
            for ilam in range(nfreq):
                wfile.write("%.7e\n"%x[ilam])
            
            wfile.write(" \n")
            if worder:
                x = ksca[::-1]
            else:
                x = ksca
            for ilam in range(nfreq):
                wfile.write("%.7e\n"%x[ilam])
            wfile.close()
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------------------------------------------------------------------------------------------------
    def runMakedust(self, freq=None, gmin=None, gmax=None, ngs=None, lnk_fname=None, gdens=None):
        """Interface function to the F77 code makedust to calculate mass absorption coefficients. 

        Parameters 
        ----------
        freq       : ndarray
                    Contains the frequency grid on which the opacities should be calculated
        
        gmin       : float
                    Minimum grain size
        
        gmax       : float
                    Maximum grain size
        
        ngs        : int
                    Number of grain sizes
        
        gdens      : float
                    Density of the dust grain in g/cm^3
        
        lnk_faname : str
                    Name of the file in which the optical constants are stored
    
        Returns
        -------

        Returns an ndarray with [nfreq,ngs] dimensions containing the resulting opacities
        """

#
# Calculate the grain sizes
#
        if ngs>1:
            gsize = gmin * (gmax/gmin)**(np.arange(ngs, dtype=np.float64)/(float(ngs)-1.))
        else:
            gsize = [gmin]

#
# Write the frequency.inp file
#
        wfile = open('frequency.inp', 'w')
        wfile.write("%d\n"%freq.shape[0])
        wfile.write("  \n")
        for i in range(freq.shape[0]):
            wfile.write("%.10e\n"%freq[i])
        wfile.close()

#
# Write the dust.inp file (makedust main control file)
#
        wfile = open('dust.inp', 'w')
        for igs in range(ngs):
            wfile.write("%s\n"%lnk_fname)
            wfile.write("%s\n"%"MIE")
            wfile.write("%d %f %f %f %d %f %f %f\n"%(1,0.0,np.log10(gsize[igs]), np.log10(gsize[igs]),1.,-3.5,gdens,-2.0))
        wfile.close()

#
# Run the Mie-code
#
        dum = sp.Popen('makedust', shell=True).wait()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class radmc3dPar(object):
    """Parameter class for a RADMC-3D model.

    
    Attributes
    ----------

    ppar   : dictionary
            Contains parameter values with parameter names as keys 
    
    pdesc  : dictionary
            Contains parameter description (comments in the parameter file) with parameter names as keys
    
    pblock : dictionary
            Contains the block names in the parameter file and parameter names as values 
    
    pvalstr: dictionary
            Contains parameter values as strings with parameter names as keys
    
    """

    def __init__(self):

        self.ppar = {}
        self.pdesc = {}
        self.pblock = {}
        self.pvalstr = {}
# --------------------------------------------------------------------------------------------------
    def readPar(self, fname=''):
        """Reads a parameter file.
        The parameters in the files should follow python syntax


        Parameters
        ----------
        
        fname  : str, optional
                File name to be read (if omitted problem_params.inp is used)

        Returns
        -------
            Returns a dictionary with the parameter names as keys
            
        """

        if fname=='':
            fname = 'problem_params.inp'

        try:
            rfile = open(fname, 'r')
        except:
            return 

        cchar  = '#'
        lbchar = ""

        # Add pi to the local variables
        pi = np.pi
    # ------------------------------------------------------------------------------------------------------------------------
    # First read every line that is not commented (i.e. does not begin with a comment character)
    # ------------------------------------------------------------------------------------------------------------------------
        dumlist = []
        dumline = '-'

        dumline = rfile.readline()
        while dumline!='':
            # First check if the line is commented out, in which case ignore it
            comment = False
            if dumline[0]==cchar:
                if dumline.find('Block')<0:
                    comment = True

            # OK, the line is not commented out, now check if it contains a '=' sign (if not ignore the line)
            if not comment:
                # Check if we have an empty line in which case also ignore it
                if dumline.strip()!='':
                    dumlist.append(dumline)

            # Read the next line
            dumline = rfile.readline()
        
        rfile.close()

    # ------------------------------------------------------------------------------------------------------------------------
    # After every line in the file was read try to decode the lines to 
    #  [variable name] = [variable value] # [comment]
    # also try to catch if an expression has been broken into multiple lines
    # ------------------------------------------------------------------------------------------------------------------------

        varlist = []
        iline = 0
        while iline<len(dumlist):
            # First check if the line contains an '=' sign if not we have a problem 
            #  expression broken into multiple lines are should already be dealt with
            ind = dumlist[iline].find('=')
            if ind<=0:
                if dumlist[iline].find('Block')<=0:
                    print 'ERROR'
                    print ' Invalid expression in line ', iline
                    print dumlist[iline]
                    print dumlist[iline+1]
                    return
                else:
                    if dumlist[iline].find(':')<=0:
                        print 'ERROR'
                        print 'Invalid block identified'
                        print 'The syntax of the block name field is :'
                        print ' # Block : Blockname '
                        return
                    else:
                        blockname = dumlist[iline].split(':')[1].strip()

            else:
                # The line contains a '=' sign and a variable name, so let's check if the
                #  value expression is broken into multiple lines
                vlist = dumlist[iline].split('=')
                lbind = vlist[1].find('\\')
                cind  = vlist[1].find('#')

                # The line is full not broken
                if lbind==-1:
                    # Check if there is a comment field
                    if cind>=0:
                        vlist = [vlist[0], vlist[1][:cind], vlist[1][cind+1:], blockname]
                    else:
                        vlist = [vlist[0], vlist[1][:cind], ' ', blockname]
                    
                    varlist.append(vlist)
                # The value expression is broken into multiple lines; take all lines and join the pieces
                else:
                    # Check if there is any comment in the line 
                    inBrokenLine = False
                    if cind>=0:
                        # Part of the line is commented, now check if the line break is before or after the comment character
                        if lbind>cind:
                            # The line break is in the comment field so there is no real line break
                            vlist = [vlist[0], vlist[1][:cind], vlist[1][cind+1:], blockname]
                        else:
                            # The line break is before the comment character 
                            inBrokenLine = True
                            expr = vlist[1][:lbind]
                            com  = vlist[1][cind+1:]
                    else: 
                        inBrokenLine = True
                        expr  = vlist[1][:lbind]
                        com   = ' '

                    if inBrokenLine:
                        # Now gather all other pieces of this line
                        iline2 = 0
                        while inBrokenLine:
                            iline2 = iline2 + 1
                            dummy = dumlist[iline + iline2]
                            # Search for comments
                            cind2 = dummy.find('#')
                            # Search for another line break
                            lbind2 = dummy.find('\\')

    # TODO:
    # At the moment I neglect the possiblity that the second line in a broken long line begins
    # with a linebreak or commented out

                            # There is comment
                            if cind2>0:

                                # There is line break
                                if lbind2>0:
                                    # The line break is commented out
                                    if lbind2>cind:
                                        expr = expr + dummy[:cind2].strip()
                                        com  = com  + dummy[cind2+1:]
                                        inBrokenLine = False
                                    else:
                                        # The line break is not commented out
                                        expr = expr + dummy[:lbind2].strip()
                                        com  = com + dummy[cind2+1:]
                                else:
                                    #There is no line break
                                    expr = expr + dummy[:cind2].strip()
                                    com  = com  + dummy[cind2+1:]
                                    inBrokenLine = False

                            # There is no comment
                            else:
                                # There is a line break
                                if lbind2>0:
                                    expr = expr + dummy[:lbind2].strip()
                                    com  = com + dummy[cind2+1:]
                                    
                                #There is no line break
                                else:
                                    expr = expr + dummy[:cind2].strip()
                                    com  = com  + ' '
                                    inBrokenLine = False
                        iline = iline + iline2 
                        vlist = [vlist[0], expr, com, blockname]
                        varlist.append(vlist)

            iline = iline + 1
    # ------------------------------------------------------------------------------------------------------------------------
    # Now evaluate the expressions in the value field and make the final dictionary
    # ------------------------------------------------------------------------------------------------------------------------
        self.ppar = {}
        glob = globals()
        glob['pi'] = np.pi
        loc  = locals()
        loc['pi'] = np.pi
        for i in range(len(varlist)):
            try:
                val= eval(varlist[i][1], glob)
                glob[varlist[i][0].strip()] = val
            except:
                try:
                    val= eval(varlist[i][1], loc)
                    loc[varlist[i][0].strip()] = val
                except:
                    print 'Unknown expression "'+varlist[i][1]+'"'
            self.ppar[varlist[i][0].strip()] = val
            self.pvalstr[varlist[i][0].strip()] = varlist[i][1].strip()
            self.pdesc[varlist[i][0].strip()] = varlist[i][2].strip()
            self.pblock[varlist[i][0].strip()] = varlist[i][3].strip()
        return

# --------------------------------------------------------------------------------------------------
    def setPar(self,parlist=[]):
        """Sets a parameter in the radmc3DPar class.
        If the paramter is already defined its value will be modified

        Parameters
        ----------
        
        parlist : list
                  If the parameter is already defined parlist should be a two element
                  list 1) parameter name, 2) parameter expression/value as a string

                  If the parameter is not yet defined parlist should be a four element
                  list 1) parameter name, 2) parameter expression/value as a string
                  3) Parameter description (= comment field in the parameter file)
        """



        parname = parlist[0].strip()

        # 
        # Check whether or not the parameter is already defined
        #
        new_par = False
        if len(parlist)==2:
            if not self.ppar.keys().__contains__(parname):
                print ' ERROR'
                print ' The argument of radmc3dPar.setPar() should be a four element list if a new'
                print ' parameter is defined 1) parameter name, 2) parameter expression/value as a string'
                print ' 3) Parameter description (= comment field in the parameter file)'
                print ' 4) The name of the block in which the parameter must be placed in the problem_params.inp file'
                return
        else:
            new_par = True

        # Add pi to the local variables
        pi = np.pi
        # 
        # Add the parameter to the dictionaries /change its value
        #
        glob = globals()
        glob['pi'] = np.pi
        loc = locals()
        loc['pi'] = np.pi
        

        try:
            self.ppar[parname] = eval(parlist[1].strip(), glob)
            glob[parname] = self.ppar[parname]
        except Exception, e:
            print e
            try:
                self.ppar[parname] = eval(parlist[1].strip(), loc)
                loc[parname] = self.ppar[parname]
            except Exception, e:
                print 'Unknown expression '+parlist[1].strip()
                print e
                return

        self.pvalstr[parname] = parlist[1].strip()
        
        if new_par:
            if not self.pdesc.has_key(parname):
                self.pdesc[parname] = parlist[2].strip()
            if len(parlist)==4:
                if not self.pblock.has_key(parname):
                    self.pblock[parname] = parlist[3].strip()


# --------------------------------------------------------------------------------------------------
    def loadDefaults(self, model='', ppar={}, reset=True):
        """Sets all parameters to default values.

        Parameters
        ----------
        model : str
                Model name whose paraemters should also be loaded
        
        ppar  : dictionary
                Contains parameter values as string and parameter names as keys
                Default values will be re-set to the values in this dictionary

        reset : bool
                If True the all class attributes will be re-initialized before
                the default values would be loaded. I.e. it will remove all entries
                from the dictionary that does not conain default values either in this
                function or in the optional ppar keyword argument
        """

        if reset:
            self.ppar = {}
            self.pvarstr = {}
            self.pdesc = {}
            self.pblock = {}

        #
        # Radiation sources
        #
        self.setPar(['incl_disc_stellarsrc', 'True', '# Switches on (True) or off (False) discrete stellar sources)', 'Radiation sources'])
        self.setPar(['mstar', '[1.0*ms]', '# Mass of the star(s)', 'Radiation sources'])
        self.setPar(['rstar','[2.0*rs]', '# Radius of the star(s)', 'Radiation sources'])
        self.setPar(['tstar','[4000.0]', '# Effective temperature of the star(s) [K]', 'Radiation sources'])
        self.setPar(['pstar','[0.0, 0.0, 0.0]', '# Position of the star(s) (cartesian coordinates)', 'Radiation sources'])
        self.setPar(['staremis_type','["blackbody"]', '# Stellar emission type ("blackbody", "kurucz", "nextgen")', 'Radiation sources'])
        self.setPar(['incl_cont_stellarsrc', 'False', '# Switches on (True) or off (False) continuous stellar sources )', 'Radiation sources'])
        #
        # Grid parameters
        #
        self.setPar(['crd_sys', "'sph'", '  Coordinate system used (car/cyl)', 'Grid parameters']) 
        self.setPar(['nx', '50', '  Number of grid points in the first dimension', 'Grid parameters']) 
        self.setPar(['ny', '30', '  Number of grid points in the second dimension', 'Grid parameters'])
        self.setPar(['nz', '36', '  Number of grid points in the third dimension', 'Grid parameters'])
        self.setPar(['xbound', '[1.0*au, 100.*au]', '  Boundaries for the x grid', 'Grid parameters'])
        self.setPar(['ybound', '[0.0, pi]', '  Boundaries for the y grid', 'Grid parameters'])
        self.setPar(['zbound', '[0.0, 2.0*pi]', '  Boundraries for the z grid', 'Grid parameters'])
        self.setPar(['xres_nlev', '3', 'Number of refinement levels (spherical coordinates only', 'Grid parameters'])
        self.setPar(['xres_nspan', '3', 'Number of the original grid cells to refine (spherical coordinates only)', 'Grid parameters'])
        self.setPar(['xres_nstep', '3', 'Number of grid cells to create in a refinement level (spherical coordinates only)', 'Grid parameters'])
        self.setPar(['wbound', '[0.1, 7.0, 25., 1e4]', '  Boundraries for the wavelength grid', 'Grid parameters'])
        self.setPar(['nw', '[19, 50, 30]', '  Number of points in the wavelength grid', 'Grid parameters'])

        #
        # Dust opacity
        #
        self.setPar(['lnk_fname', "['/disk2/juhasz/Data/JPDOC/astrosil/astrosil_WD2001_new.lnk', '/disk2/juhasz/Data/JPDOC/carbon/A/cel600.lnk']", ' ', 'Dust opacity'])
        self.setPar(['gdens', '[3.6, 1.8]', ' Bulk density of the materials in g/cm^3', 'Dust opacity'])
        self.setPar(['gsmin', '0.1', ' Minimum grain size', 'Dust opacity'])
        self.setPar(['gsmax', '10.0', ' Maximum grain size', 'Dust opacity'])
        self.setPar(['ngs', '1', ' Number of grain sizes', 'Dust opacity'])
        self.setPar(['gsdist_powex', '-3.5', ' Grain size distribution power exponent', 'Dust opacity'])
        self.setPar(['mixabun',       '[0.75, 0.25]', ' Mass fractions of the dust componetns to be mixed', 'Dust opacity'])
        self.setPar(['dustkappa_ext',"['silicate']", ' ', 'Dust opacity'])
        
        #
        # Gas line RT 
        #
        self.setPar(['gasspec_mol_name', "['co']", '  Name of the gas species - the extension of the molecule_EXT.inp file', 'Gas line RT'])
        self.setPar(['gasspec_mol_abun', '[1e-4]', '  Abundance of the molecule', 'Gas line RT']) 
        self.setPar(['gasspec_mol_dbase_type',"['leiden']", '  leiden or linelist', 'Gas line RT'])
        self.setPar(['gasspec_colpart_name', "['h2']", '  Name of the gas species - the extension of the molecule_EXT.inp file', 'Gas line RT'])
        self.setPar(['gasspec_colpart_abun', '[1e0]', '  Abundance of the molecule', 'Gas line RT']) 
        #self.setPar(['gasspec_vturb', '0.1e5', '  Microturbulence', 'Gas line RT'])
        #self.setPar(['writeGasTemp', 'False', '  Whether or not to write a separate gas temperature file (gas_temperature.inp) if such function exists in the model', 'Gas line RT'])
        #
        # Code parameters
        #
        self.setPar(['nphot', 'long(3e5)', '  Nr of photons for the thermal Monte Carlo', 'Code parameters'])
        self.setPar(['nphot_scat','long(3e5)', '  Nr of photons for the scattering Monte Carlo (for images)', 'Code parameters'])
        self.setPar(['nphot_spec','long(1e5)', '  Nr of photons for the scattering Monte Carlo (for spectra)', 'Code parameters'])
        self.setPar(['scattering_mode_max','1', '  0 - no scattering, 1 - isotropic scattering, 2 - anizotropic scattering', 'Code parameters'])
        self.setPar(['lines_mode', '-1', '  Line raytracing mode', 'Code parameters'])
        self.setPar(['istar_sphere', '0', '  1 - take into account the finite size of the star, 0 - take the star to be point-like', 'Code parameters'])
        self.setPar(['itempdecoup', '1', '  Enable for different dust components to have different temperatures', 'Code parameters'])
        self.setPar(['tgas_eq_tdust', '1', '  Take the dust temperature to identical to the gas temperature', 'Code parameters'])
        self.setPar(['rto_style', '1', '  Format of outpuf files (1-ascii, 2-unformatted f77, 3-binary', 'Code parameters'])
        self.setPar(['modified_random_walk', '0', '  Switched on (1) and off (0) modified random walk', 'Code parameters'])
        #
        # Model parameters
        #
        if model!='':
            try:
                mdl = __import__(model)
            except:
                try:
                    mdl  = __import__('radmc3dPy.models.'+model, fromlist=['']) 
                except:
                    print 'ERROR'
                    print ' '+model+'.py could not be imported'
                    print ' The model files should either be in the current working directory or'
                    print ' in the radmc3d python module directory'
                    return

            modpar = mdl.getDefaultParams()
            for i in range(len(modpar)):
                dum = modpar[i]
                if len(dum)==3:
                    dum.append('Model '+model)
                self.setPar(dum)
        
# --------------------------------------------------------------------------------------------------
    def printPar(self):
        """Prints the parameters of the current model.
        
        """
        
        #
        # First get the unique block names 
        #

        blocknames = ['Radiation sources', 'Grid parameters', 'Dust opacity', 'Gas line RT', 'Code parameters']
        for key in self.pblock.keys():
            dum = self.pblock[key]
            if not blocknames.__contains__(dum):
                blocknames.append(dum)

       
        #
        # Get the parameter block names and distionary keys
        #
        par_keys = self.pblock.keys()
        par_block = self.pblock.values()

        #
        # Print the parameters by blocks 
        #
        for iblock in blocknames:
            print ('%s'%'# -------------------------------------------------------------------------------------------------------------------------')
            txt = '# Block: '+iblock
            print ('%s'%txt)
            print ('%s'%'# -------------------------------------------------------------------------------------------------------------------------')
           

            keys = []
            for i in range(len(par_block)):
                if par_block[i]==iblock:
                    keys.append(par_keys[i])

            keys.sort()
            for key in keys:
                print (key.ljust(25) + ' = ' + self.pvalstr[key].strip() + '  # ' + self.pdesc[key].strip())
# --------------------------------------------------------------------------------------------------
    def writeParfile(self, fname=''):
        """Writes a parameter file.

        Parameters
        ----------
        
        fname  : str, optional
                File name to be read (if omitted problem_params.inp is used)

        """
        
        if fname=='':
            fname = 'problem_params.inp'

        print 'Writing '+fname
    
        #
        # First get the uniq block names 
        #

        blocknames = ['Radiation sources', 'Grid parameters', 'Dust opacity', 'Gas line RT', 'Code parameters']
        for key in self.pblock.keys():
            dum = self.pblock[key]
            if not blocknames.__contains__(dum):
                blocknames.append(dum)

        
        try :
            wfile = open(fname, 'w')
        except:
            print ' ERROR '
            print ' Cannot create '+fname 
            return
        #
        # Write header
        #

        wfile.write('%s\n'%'###########################################################################################################################')
        wfile.write('%s\n'%'# RADMC-3D PARAMETER SETUP')
        wfile.write('%s\n'%'# Created by the python module of RADMC-3D')
        wfile.write('%s\n'%'###########################################################################################################################')
       
        #
        # Get the parameter block names and distionary keys
        #
        par_keys = self.pblock.keys()
        par_block = self.pblock.values()

        #
        # Write the parameterfile
        #
        for iblock in blocknames:
            wfile.write('%s\n'%'# -------------------------------------------------------------------------------------------------------------------------')
            txt = '# Block: '+iblock
            wfile.write('%s\n'%txt)
            wfile.write('%s\n'%'# -------------------------------------------------------------------------------------------------------------------------')
           

            keys = []
            for i in range(len(par_block)):
                if par_block[i]==iblock:
                    keys.append(par_keys[i])

            keys.sort()
            for key in keys:
                wfile.write(key.ljust(25) + ' = ' + self.pvalstr[key].strip() + '  # ' + self.pdesc[key].strip() + '\n')

# --------------------------------------------------------------------------------------------------
# Functions for an easy compatibility with the IDL routines
# --------------------------------------------------------------------------------------------------
def readOpac(ext=[''], idust=None, scatmat=None, old=False):
    """Reads the dust opacity files.
    This function is an interface to radmc3dDustOpac.readOpac()

    Parameters
    ----------
    ext   : list
            Each element of the list is be a string, the file name extension (file names should look like 'dustkappa_ext.inp')
    
    idust : list
            Each element of the list is an integer, the index of the dust species in the master opacity file (dustopac.inp')

    scatmat: list
            If specified, its elements should be booleans indicating whether the opacity file 
            contains also the full scattering matrix (True) or only dust opacities (False)
        
    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used
    
    Returns
    -------
        Returns an instance of the radmc3dDustOpac class 
    """


    res = radmc3dDustOpac()
    res.readOpac(ext=ext, idust=idust, scatmat=scatmat, old=old)
    
    return res
# --------------------------------------------------------------------------------------------------
# Functions for an easy compatibility with the IDL routines
# --------------------------------------------------------------------------------------------------
def readData(ddens=False, dtemp=False, gdens=False, gtemp=False, gvel=False, ispec=None, vturb=False, binary=True, old=False):
    """Reads the physical variables of the model (e.g. density, velocity, temperature).

    Parameters
    ----------

    ddens : bool
            If True dust density will be read (all dust species and grain sizes)
    
    dtemp : bool
            If True dust temperature will be read (all dust species and grain sizes)
    
    gdens : bool
            If True gas density will be read (NOTE: the gas density will be number density in 1/cm^3)
    
    gtemp : bool
            If True gas temperature will be read (all dust species and grain sizes)
    
    gvel  : bool
            If True the velocity field will be read
    
    ispec : str
            Name of the molecule in the 'molecule_ispec.inp' filename
        
    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used

    Returns
    -------
    Returns an instance of the radmc3dData class 
    """

    res = radmc3dData()
    if ddens: res.readDustDens(binary=binary, old=old)
    if dtemp: res.readDustTemp(binary=binary, old=old)
    if gvel: res.readGasVel(binary=binary)
    if gtemp: res.readGasTemp(binary=binary)
    if vturb: res.readVTurb(binary=binary)
    if gdens:
        if not ispec:
            print 'ERROR'
            print 'No gas species is specified!'
            print 'The ispec input keyword should be set to the name of the gas species as it appears in '
            print ' numberdens_gasspecname.inp'
            return 0
        else:
            res.readGasDens(ispec=ispec,binary=binary)

    return res

# --------------------------------------------------------------------------------------------------
def readGrid():
    """Reads the spatial and frequency grid.
    This function is an interface to radmc3dGrid.readGrid().

    Returns
    -------

    Returns an instance of the radmc3dGrid class 
    """

    grid = radmc3dGrid()
    grid.readGrid()

    return grid

# --------------------------------------------------------------------------------------------------
def readParams():
    """Reads the problem_params.inp file.
    This function is an interface to radmc3dPar.readPar().

    Returns
    -------
    Returns an instance of the radmc3dPar class 

    """

    dum = radmc3dPar()
    dum.readPar()
    return dum
# --------------------------------------------------------------------------------------------------
def writeDefaultParfile(model='', fname=''):
    """Writes a parameter file (problem_params.inp) with default parameters for a given model.

    Parameters
    ----------
    
    model : str
            Name of the model whose parameter should be written to the file

    fname : str, optional
            Name of the parameter file to be written (if omitted problem_params.inp will be used)
    """
    
    if model=='':
        print ' ERROR '
        print ' No model name is given '
        return

    dum  = radmc3dPar()
    dum.loadDefaults(model=model)
    dum.writeParfile()
# --------------------------------------------------------------------------------------------------
def readSpectrum(fname='', old=False):
    """Reads the spectrum / SED


    Parameters
    -----------
    fname : str, optional
            Name of the file to be read
        
    old   : bool, optional
            If set to True the file format of the previous, 2D version of radmc will be used


    Returns
    -------

        Returns an ndarray with [Nwavelength, 2] dimensions 
        [Nwavelength,0] is the wavelength / velocity and
        [Nwavelength,1] is the flux density
        
    """
  
    if not old:
        if fname.strip()=='':
            fname = 'spectrum.out'

        rfile = open(fname, 'r')
        # Read the format number
        dum = rfile.readline()
        # Read the number of wavelengths 
        nwav = int(rfile.readline())
        # Read a blank line
        dum = rfile.readline()
        
        res = np.zeros([nwav, 2], dtype=np.float64)
        for iwav in range(nwav):
            dum = rfile.readline().split()
            res[iwav,0] = float(dum[0])
            res[iwav,1] = float(dum[1])
        rfile.close()

    else:
        if fname.strip()=='':
            fname = 'spectrum.dat'

        try :
            rfile = open(fname, 'r')
        except:
            print 'Error!' 
            print fname+' could not be read!'
            return 

        # Read the number of wavelengths 
        nwav = int(rfile.readline())
        rfile.readline()
        cc = 29979245800.
       
        res = np.zeros([nwav,2], dtype=float)
        for iwav in range(nwav):
            dum = rfile.readline().split()
            res[iwav,0] = cc/float(dum[0])*1e4
            res[iwav,1] = float(dum[1])

        rfile.close()

    return res

# --------------------------------------------------------------------------------------------------
def getDensVstruct(data=None, vmean_temp=False, ispec_tgas=0, gsize=None, idust=None, mstar=0.):
    """Calculates the vertical hydrostatic equilibrium

    Parameters
    ----------
    data        : radmc3dData
                  An instance of the radmc3DData class containing the density structure of the model

    vmean_temp  : bool
                  If True (T(z) = T(-z) = 0.5*(T(z) + T(-z))) if False (T(z)!=T(-z)) 
    
    idust       : list
                  List of dust indices whose structure must be calculated
    
    mstar       : float
                  Stellar mass
    
    ispec_tgas  : int
                  Index of dust species whose temperature is taken to be the gas temperature
    
    gsize       : ndarray
                  Dust grain sizes - If specified, the gas temperature is calculated as the average temperature
                  of all dust grains in the grid cell weighted by the total surface area of dust grains with given
                  size - NOTE: this approach assumes that all dust grains of a given size have the same bulk density

    Returns
    -------
    Returns an ndarray with the dust density
    """
        
    # Fix the mean molecular weight to 2.3
    mu = 2.3

    # Pre-calculate some constants
    A  = mu*mp*gg*mstar / kk
    cost  = np.cos(data.grid.y)
    costi = np.cos(data.grid.yi)

    if not mstar:
        print 'ERROR'
        print ' You should specify the stellar mass (mstar = ??)'
        return

    if idust==None:
        print ' No dust index was given for which the vertical structure should be calculated'
        print ' So we do for all dust species'
        idust = range(data.rhodust.shape[3])
    else:
        if (type(idust).__name__=='int')| (type(idust).__name__=='float'):
            idust = [idust]
    # To improve the smoothness of the temperature structure, if the density structure is
    #  symmetric to the disk midplane we use T_new(theta) = T_new(pi-theta) = 0.5 * (T(theta) + T(pi-theta))
    if vmean_temp:

        if abs(data.grid.yi[data.grid.nyi-1]-np.pi/2.)<1e-8:
            print 'ERROR'
            print "Cannot average temperature in the vertical direction if theta mirroring is active"
            return None
        else:
            print ' Smoothing the vertical temperature structure by averaging the temperature of the '
            print " two half planes above and below the disk midplane"
            dusttemp = np.zeros(data.dusttemp.shape, dtype=np.float64)
            for iy in range(data.grid.ny/2):
                print iy
                dusttemp[:,iy,:,:] = 0.5 * (data.dusttemp[:,iy,:,:] + data.dusttemp[:,data.grid.ny-1-iy,:,:])
                dusttemp[:,data.grid.ny-1-iy,:,:] = dusttemp[:,iy,:,:]
    # Calculate the vertical hydrostatic equilibrium for the two half space (z<0, z>0) separately
    else:
        dusttemp = data.dusttemp
      
    #rho_new = np.zeros(data.rhodust.shape, dtype=np.float64)
    rho_new = np.array(data.rhodust)
    if len(gsize)!=0:
        mean_dusttemp = np.zeros([data.grid.nx, data.grid.ny, data.grid.nz], dtype=np.float64)
        w             = np.zeros(data.rhodust.shape, dtype=np.float64)
        for ispec in idust:
            w[:,:,:,ispec] = gsize[ispec]**2 * (data.rhodust[:,:,:,ispec] / gsize[ispec]**3) 
        
        wnorm = w.sum(3)
        for ispec in idust:
            w[:,:,:,ispec] = w[:,:,:,ispec]/wnorm

        for ispec in idust:
            mean_dusttemp = mean_dusttemp + data.dusttemp[:,:,:,ispec] * w[:,:,:,ispec]

    # Loop over all dust species where we should calculate the vertical structure
    for ispec in idust:
        rho_new[:,:,:,ispec] = 0.
        for ir in range(data.grid.nx):
            print ir, data.grid.nx-1
            r     = data.grid.x[ir]
            z     = r * cost
            zi    = r * costi
            dz    = z[:-1] - z[1:]
            const = A / r**3

            # Do we have theta mirroring active?
            if abs(data.grid.yi[data.grid.nyi-1]-np.pi/2.)<1e-8:
                for ip in range(data.grid.nz):
                    dlgrho  = np.log(data.rhodust[ir,1:,ip,ispec]) - np.log(data.rhodust[ir,:-1,ip,ispec])
                    if len(gsize)!=0:
                        temp    = mean_dusttemp[ir,:,ip] 
                    else:
                        temp    = dusttemp[ir,:,ip,ispec]
                    
                    it = data.grid.ny-1
                    temp[it] = 0.5 * (temp[it] + temp[it-1])    

                    dlgtemp = np.log(temp[1:]) - np.log(temp[:-1])
                    zpt     = z/temp
                    zpt     = 0.5*(zpt[1:] + zpt[:-1])


                    # Calculate the normalized (rho[z=0] = 1.0) density
                    rho_new[ir,data.grid.ny-1,ip,ispec] = 1.0

                    for it in range(data.grid.ny-1)[::-1]:
                        rho_new[ir,it,ip,ispec] = rho_new[ir,it+1,ip,ispec] * np.exp(-(const*zpt[it] + dlgtemp[it]/dz[it])*dz[it])
                
                    rho_new = rho_new.clip(1e-90, 1e90)
                  
                    # Now re-normalize the surface density to the input value
                    sigma = (data.rhodust[ir,:,ip,ispec] * (zi[1:] - zi[:-1])).sum()
                    sigma_new = (rho_new[ir,:,ip,ispec] * (zi[1:] - zi[:-1])).sum()

                    rho_new[ir,:,ip,ispec] = rho_new[ir,:,ip,ispec] * sigma / sigma_new

            else:
                for ip in range(data.grid.nz):
                    dlgrho  = np.log(data.rhodust[ir,1:,ip,ispec]) - np.log(data.rhodust[ir,:-1,ip,ispec])
                    if len(ispec_weights)!=0:
                        temp    = (dusttemp[ir,:,ip,ispec]*ispec_weights).sum()
                    else:
                        temp    = dusttemp[ir,:,ip,ispec]
                    dlgtemp = np.log(temp[1:]) - np.log(temp[:-1])
                    zpt     = z/temp
                    zpt     = 0.5*(zpt[1:] + zpt[:-1])

                    # Calculate the normalized (rho[z=0] = 1.0) density
                    rho_new[ir,data.grid.ny/2-1,ip,ispec] = 1.0
                    rho_new[ir,data.grid.ny/2,ip,ispec] = 1.0

                    for it in range(data.grid.ny/2)[::-1]:
                        rho_new[ir,it-1,ip,ispec] = rho_new[ir,it,ip,ispec] * exp(-(const*zpt[it] + dlgtemp[it]/dz[it])*dz[it])
                    for it in range(data.grid.ny/2, data.grid.ny-1)[::1]:
                        rho_new[ir,it,ip,ispec] = rho_new[ir,it-1,ip,ispec] * exp((const*zpt[it-1] + dlgtemp[it-1]/dz[it-1])*dz[it-1])
                   
                    rho_new = rho_new.clip(1e-90, 1e90)

                    # Now re-normalize the surface density to the input value
                    sigma = (data.rhodust[ir,:,ip,ispec] * (zi[1:] - zi[:-1])).sum()
                    sigma_new = (rho_new[ir,:,ip,ispec] * (zi[1:] - zi[:-1])).sum()

                    rho_new[ir,:,ip,ispec] = rho_new[ir,:,ip,ispec] * sigma / sigma_new
            
            print rho_new[ir,data.grid.ny/2-1,ip,ispec]

    return rho_new

# --------------------------------------------------------------------------------------------------
def plotSpectrum(a,ev=False,kev=False,hz=False,micron=False,jy=False,lsun=False,
                 lnu=False,nulnu=False,fnu=False,nufnu=False,dpc=1.e0,
                 oplot=False,xlg=False,ylg=False,obs=False,
                 mol=None,iline=None):
    """Plot the spectrum / SED 

    Parameters
    ----------
    ev              = True --> frequency in electronvolt (default=Hz)

    kev             = True --> frequency in kiloelectronvolt (default=Hz)

    micron          = True --> wavelength in micron (default=Hz)

    jy              = True --> Flux in Jansky

    lnu             = True --> L_nu (default L_nu)

    nulnu           = True --> nu*L_nu (default F_nu)

    lsun            = True --> nu*L_nu in units of solar luminosity

    dpc             = Distance of observer in units of parsec
                      Default: 1 pc

    oplot           = True --> Plot without refreshing subplot

    xlg             = True --> logarithmic x-axis

    ylg             = True --> logarithmic y-axis

    obs             = True --> Treat the spectrum as an observation
                               (i.e. do not scale with dpc^(-2))

    mol             = (optional) Molecule data (see radmc3dMolecule class)
                      This is required if you want to plot a line spectrum
                      with on the x-axis the radial velocity in km/s

    iline           = (if set) the index of the line (of mol; starting,
                      as in RADMC-3D, with the index 1) which shall act
                      as the 0 km/s wavelength reference. If iline is set
                      the x axis will be in km/s (overriding other settings)

    """
    #
    # Make sure to have the plotting library loaded
    #
    from matplotlib import pyplot as plt
    #
    # Basic
    #
    lam    = a[:,0]
    fluxnu = a[:,1]
    #
    # Calculate frequency in Hz
    #
    cc    = 2.9979245800000e10     # Light speed             [cm/s]
    freq  = 1e4*cc/lam
    #
    # Default: frequency in Hz
    #
    xcoord = freq
    xtitle = '$\\nu [\mathrm{Hz}]$'
    #
    # If ev: electronvolt
    #
    if ev:
        xcoord = 4.13568842841e-15 * freq
        xtitle = '$h\\nu [\mathrm{eV}]$'
    #
    # If kev: kiloelectronvolt
    #
    if kev:
        xcoord = 4.13568842841e-18 * freq
        xtitle = '$h\\nu [\mathrm{KeV}]$'
    #
    # If micron
    #
    if micron:
        xcoord = lam
        xtitle = '$\lambda [\mu\mathrm{m}]$'
    #
    # Plot nuFnu or Fnu (same with Lnu)? And what about Fnu vs Lnu?
    #
    # Default:
    sed=True
    ylum=False
    # The flags:
    if jy:
        sed=False
    if fnu:
        sed=False
        ylum=False
    if lnu:
        sed=False
        ylum=True
    if nulnu:
        sed=True
        ylum=True
    if fnu:
        sed=False
        ylum=False
    if nufnu:
        sed=True
        ylum=False
    if jy:
        ylum=False
    if lsun:
        ylum=True
        sed=True
    #
    # If iline is set, then override the above and use instead the line
    # as a reference and use km/s as x-axis
    #
    if iline is not None:
        if mol is None:
            print "Error in plotSpectrum(): if you specify iline, you must give a molecule with mol=..."
            return
        else:
            freq0  = mol.freq[iline-1]
            xcoord = 2.9979245800000e+10*(freq0-freq)/freq0/1.e5
            xtitle = '$\Delta v [\mathrm{km/h}]$'
    #
    # Which plot to make? Lum or flux?
    #
    if not ylum:
        #
        # Plot spectrum as flux at a certain distance
        #
        if not obs:
            distfact = 1.0 / (dpc**2)
        else:
            distfact = 1.0
        #
        # Set the vertical axis name
        #
        if not jy:
            if not sed:
                lumfact=1.0
                ytitle='$F_{\\nu}\; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{Hz}^{-1}\, \mathrm{s}^{-1}]$'
            else:
                lumfact=1.0*freq
                ytitle='$\\nu F_{\\nu}\; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$'
        else:
            if not sed:
                lumfact=1e+23
                ytitle='$F_{\\nu} [Jy]$'
            else:
                lumfact=1e+23*freq
                ytitle='$\\nu F_{\\nu} [JyHz]$'
    else:
        #
        # Plot spectrum as luminosity
        #
        if not obs:
            distfact = 1.1965280793e38   # = 4*pi*(1 parsec)^2 = 1.19d38 cm^2
        else:
            distfact = dpc**2 * 1.1965280793e38

        if not sed:
            lumfact=1.e0
            ytitle='$L_{\\nu}\; [\mathrm{erg}\,\mathrm{Hz}^{-1}\, \mathrm{s}^{-1}]$'
        else:
            if not lsun:
                lumfact = 1.0*freq
                ytitle  = '$\\nu L_{\\nu}\; [\mathrm{erg}\, \mathrm{s}^{-1}]$'
            else:
                lumfact = freq * 2.5956986e-34
                ytitle  = '$\\nu L_{\\nu}\; [L_{\odot}]$'

    #
    # The data on the y axis
    #
    ycoord = distfact*lumfact*fluxnu
    #
    # If not oplot, then reset the subplot and set the axes
    #
    if not oplot:
        plt.cla()
        if xlg:
            plt.xscale('log')
        if ylg:
            plt.yscale('log')
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
    #
    # Now plot
    #
    plt.plot(xcoord,ycoord)


# --------------------------------------------------------------------------------------------------
class radmc3dMolecule:
    """
    RADMC-3D molecule class
    Based on the Leiden LAMDA database, but is in principle generic

    NOTE: For now only the levels and lines are included, not the 
          collision rates. 

    Attributes
    ----------
    name            = The name as listed in the molecule file
    molweight       = Molecular weight in units of proton mass
    nlevels         = Nr of energy levels
    nlines          = Nr of lines
    energycminv     = Energy[ilevel] of level ilevel in 1/cm
    energy          = Energy[ilevel] of level ilevel in erg
    wgt             = Statistical weight[ilevel] of level ilevel
    jrot            = Quantum rotational J[ilevel] of level ilevel
    iup             = ilevel of upper level of line iline (starting with 0)
    ilow            = ilevel of lower level of line iline (starting with 0)
    aud             = Einstein A up low of line iline in 1/second
    freq            = Frequency of line iline in Hz
    lam             = Wavelength of line iline in micron

    """

    def __init__(self):
        self.name        = ""
        self.molweight   = 0.0
        self.nlevels     = 0
        self.nlines      = 0
        self.energycminv = 0.0
        self.energy      = 0.0
        self.wgt         = 0.0
        self.jrot        = 0.0
        self.iup         = 0
        self.ilow        = 0
        self.aud         = 0.0
        self.freq        = 0.0
        self.lam         = 0.0

    # --------------------------------------------------------------------------------------------------
    def read(self,mol):
        """Read the molecule_<mol>.inp file

        The file format is the format of the Leiden LAMDA molecular database

        Parameters
        ----------
        mol             = molecule name (e.g. 'co')

        """

        filename = 'molecule_'+mol+'.inp'
        with open(filename,'r') as f:
            dum             = f.readline()
            dum             = f.readline().split()
            self.name       = dum[0]
            dum             = f.readline()
            self.molweight  = float(f.readline())
            dum             = f.readline()
            self.nlevels    = int(f.readline())
            dum             = f.readline()
            self.energycminv= np.zeros(self.nlevels)
            self.energy     = np.zeros(self.nlevels)
            self.wgt        = np.zeros(self.nlevels)
            self.jrot       = np.zeros(self.nlevels)
            for i in range(self.nlevels):
                dum                 = f.readline().split()
                self.energycminv[i] = float(dum[1])
                self.energy[i]      = float(dum[1])*1.9864847851996e-16  # const=h*c
                self.wgt[i]         = float(dum[2])
                self.jrot[i]        = float(dum[3])
            dum             = f.readline()
            self.nlines     = int(f.readline())
            dum             = f.readline()
            self.iup        = np.zeros(self.nlines)
            self.ilow       = np.zeros(self.nlines)
            self.aud        = np.zeros(self.nlines)
            self.freq       = np.zeros(self.nlines)
            self.lam        = np.zeros(self.nlines)
            for i in range(self.nlines):
                dum            = f.readline().split()
                self.iup[i]    = int(dum[1])   # Use as index: [iup-1]
                self.ilow[i]   = int(dum[2])   # Use as index: [ilow-1]
                self.aud[i]    = float(dum[3])
                self.freq[i]   = float(dum[4])*1e9
                self.lam[i]    = 2.9979245800000e+14/self.freq[i]

# --------------------------------------------------------------------------------------------------
def readMol(mol):
    """ Wrapper around the radmc3dMolecule.read() method

        Parameters
        ----------
        mol             = molecule name (e.g. 'co')

    """

    m = radmc3dMolecule()
    m.read(mol)
    return m

