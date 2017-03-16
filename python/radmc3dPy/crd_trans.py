"""
This module contains functions for coordinate transformations (e.g. rotation).
"""
try:
    import numpy as np
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'


#import scipy.linalg as la

def ctrans_sph2cyl(crd=None, reverse=False):
    """Transform coordinates between spherical to cylindrical systems

    Parameters
    ----------
    crd      : ndarray
               Three element array containing the input coordinates
               [r,phi,theta] or [r, phi, z] if reverse is True

    reverse  : bool
               If True it calculates the inverse trasnformation
               (cylindrical -> spherical). In this case crd should be [r,phi,z]

    Returns
    --------
    Returns a numpy array of [Nr,Nphi,Ntheta,3] dimensions containing the cylindrical 
       coordinates [rcyl, z, phi]
    """
    
    nr     = r.shape[0]
    nt     = theta.shape[0]
    result = np.zeros([nr, nt], dtype=np.float64)
    for ix in range(grid.nr):
        result[ix,:,0] = r[ix] * np.sin(theta)
        result[ix,:,1] = r[ix] * np.cos(theta)

    return result


def ctrans_sph2cart(crd=[0,0,0], reverse=False):
    """Transform coordinates between spherical to cartesian systems

    Parameters
    ----------
    crd      : ndarray
               Three element array containing the input
               coordinates [x,y,z] or [r,theta,phi] by default
               the coordinates assumed to be in the cartesian system
    
    reverse  : bool
               If True calculates the inverse trasnformation
               (cartesian -> spherical). In this case crd should be [r,theta,phi]

    Returns
    -------
    Returns a three element array containg the output coordinates [r,theta,phi] or [x,y,z]
    """
    if (reverse==False):
        r     = crd[0]
        theta = crd[1] + 1e-50
        phi   = crd[2]

        x = np.sin(theta) * np.cos(phi) * r
        y = np.sin(theta) * np.sin(phi) * r
        z = np.cos(theta) * r

        crdout = [x, y, z]

    else:

        x = crd[0]
        y = crd[1]
        z = crd[2]

        r     = np.sqrt(x**2 + y**2 + z**2)
        phi   = np.arccos(x / np.sqrt(x**2 + y**2) + 1e-90)
        theta = np.arccos(z / r)

        if (y<0.): phi = 2.0 * np.pi - phi
        
        crdout = [r, theta, phi]
        

    return crdout

def vtrans_sph2cart(crd=[0,0,0], v=[0,0,0], reverse=False):
    """Transform velocities between spherical to cartesian systems

    Parameters
    ----------
    crd     : ndarray
              Three element array containing the input
              coordinates [x,y,z] or [r,theta,phi] by default
              the coordinates assumed to be in the cartesian system

    v       : ndarray
              Three element array containing the input
              velocities in the same coordinate system as crd


    reverse : bool
              If True it calculates the inverse trasnformation (cartesian -> spherical)

    Returns
    -------

    Returns a three element array containg the output velocities [vr,vphi,vtheta] or [vx,vy,vz]


    """
    
    #NOTE!!!!! The velocities in the spherical system are not angular velocities!!!!
    #v[1] = dphi/dt * r
    #v[2] = dtheta/dt * r

    if (reverse==False):
        r      = crd[0]
        theta  = crd[1]
        phi    = crd[2]
        
        vr     = v[0]
        vtheta = v[1]
        vphi   = v[2]
        
        vx     = vr*np.sin(theta)*np.cos(phi) - vphi*np.sin(phi) + vtheta*np.cos(theta)*np.cos(phi)
        vy     = vr*np.sin(theta)*np.sin(phi) + vphi*np.cos(phi) + vtheta*np.cos(theta)*np.sin(phi)
        vz     = vr*np.cos(theta) - vtheta*np.sin(theta)

        vout   = [vx,vy,vz]

    else:
        
       # crd_sph = ctrans_sph2cart(crd, reverse=True)
        #r       = crd_sph[0]
        #theta   = crd_sph[1]
        #phi     = crd_sph[2]
        
        #a       = [[np.sin(theta)*np.cos(phi), -np.sin(phi), np.cos(theta)*np.cos(phi)],\
                   #[np.sin(theta)*np.sin(phi), np.cos(phi), np.cos(theta)*np.sin(phi)],\
                   #[np.cos(theta), 0., -np.sin(theta)]]
        
        #a       = [[np.sin(theta)*np.cos(phi), np.cos(theta)*np.cos(phi), -np.sin(phi)],\
                   #[np.sin(theta)*np.sin(phi), np.cos(theta)*np.sin(phi), np.cos(phi)],\
                   #[np.cos(theta), -np.sin(theta),0.]]


        #a       = np.array(a, dtype=np.float64)

        #vout = np.linalg.solve(a,v)


        #
        # New stuff
        #
        vout = np.zeros(3, dtype=float)
        r    = np.sqrt((crd**2).sum())
        rc   = np.sqrt(crd[0]**2 + crd[1]**2)

        # Vr
        vout[0] = (crd * v).sum()/r             
        # Vtheta
        vout[1] = (crd[2]*(crd[0]*v[0] + crd[1]*v[1]) - v[2]*rc**2) / (r*rc)
        # Vphi
        vout[2] = (crd[0]*v[1] - crd[1]*v[0])/rc


    return vout

def csrot(crd=None, ang=None, xang=0.0, yang=0.0, zang=0.0, deg=False):
    """ Performs coordinate system rotation.

    Parameters
    ----------

    crd : numpy ndarray
          three element vector containing the coordinates of a
          given point in a cartesian system

    ang : list, numpy ndarray
          three element array, angles of rotation around the x,y,z axes

    deg : float, optional 
          If True angles should be given in degree instead of radians (as by default)
 
   
    Returns
    -------
    list
        Returns a three element list with the rotated coordinates


    Notes
    -----

    Rotation matrices
    
    Around the x-axis:
    
    .. math:: 
         \\left(\\begin{matrix} 
                 1 & 0 & 0 \\\\
                 0 & cos(\\alpha) & -sin(\\alpha)\\\\
                 0 & sin(\\alpha) & cos(\\alpha)
                 \\end{matrix}\\right)

    Around the y-axis:

    .. math:: 
         \\left(\\begin{matrix} 
                 cos(\\beta) & 0 & -sin(\\beta) \\\\
                 0 & 1 & 0\\\\
                 sin(\\beta)& 0 & cos(\\beta)
                 \\end{matrix}\\right)

    Around the z-axis

    .. math:: 
         \\left(\\begin{matrix} 
                 cos(\\gamma) & -sin\\gamma) & 0 \\\\
                 sin(\\gamma) &  cos(\\gamma) & 0 \\\\
                 0  & 0 & 1
                 \\end{matrix}\\right)
 

    """
    crd_new = np.zeros(len(crd), dtype=np.float64)

    if ang==None:
        if (xang==0.)&(yang==0.)&(zang==0.):
            return crd

    if (ang!=None):
        xang = ang[0]
        yang = ang[1]
        zang = ang[2]

#
# Convert degree into radian if the angles are given in degree
#

    if (deg==True):
        xang = xang / 180.0 * np.pi
        yang = yang / 180.0 * np.pi
        zang = zang / 180.0 * np.pi

#
# Rotation around the x axis
#

    if (xang!=0.0):
        dumx = crd[0]
        dumy = np.cos(xang)*crd[1]  - np.sin(xang)*crd[2]
        dumz = np.sin(xang)*crd[1] + np.cos(xang)*crd[2] 
    
        crd_new = [dumx, dumy, dumz]

#
# Rotation around the y axis
#

    if (yang!=0.0):
        dumx = np.cos(yang)*crd[0] + np.sin(yang)*crd[2]
        dumy = crd[1]
        dumz = -np.sin(yang)*crd[0] + np.cos(yang)*crd[2]
        
        crd_new = [dumx, dumy, dumz] 

#
# Rotation around the z axis
#
  
    if (zang!=0.0):
        dumx = np.cos(zang)*crd[0] - np.sin(zang)*crd[1] + 0.0
        dumy = np.sin(zang)*crd[0] + np.cos(zang)*crd[1] + 0.0
        dumz = crd[2]
        
        crd_new = [dumx, dumy, dumz]



    return crd_new

def vrot(crd=None, v=None, ang=None):
    """Rotates a vector in spherical coordinate system.
    First transforms the vector to cartesian coordinate system, then does the rotation then 
    makes the inverse transformation

    Parameters 
    ----------
    crd  : ndarray
           Three element array containing the coordinates of a
           given point in the cartesian system

    v    : ndarray
           Three element array, angles of rotation around the x,y,z axes

    ang  : ndarray 
           Three element arrray containing the angles to rotate around the x, y, z, axes, respectively

    """
# Convert the position vector to cartesian coordinate system
    crd_xyz   = ctrans_sph2cart(crd=crd)
# Convert the velocity vector to cartesian coordinate system
    v_xyz     = vtrans_sph2cart(crd=crd, v=v)
# Rotate the vector
    v_xyz_rot = csrot(crd=v_xyz, ang=ang)
# Transform the rotated vector back to the spherical coordinate system
    v_rot     = vtrans_sph2cart(crd=crd_xyz, v=v_xyz_rot, reverse=True) 
    
    return v_rot
