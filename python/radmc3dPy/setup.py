"""This module contains functions to set up a RADMC-3D model for dust and/or line simulations.
For help on the syntax or functionality of each function see the help of the individual functions
"""

try:
    import numpy as np
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'

import subprocess as sp
import os, sys, copy

from radmc3dPy.natconst import *
import radmc3dPy.analyze as analyze
import inspect
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def problemSetupDust(model='', binary=True, writeDustTemp=False, old=False, dfunc=None, dfpar=None, **kwargs):
    """
    Function to set up a dust model for RADMC-3D 
    
    Parameters
    ----------
    model           : str
                      Name of the model that should be used to create the density structure.
                      The file should be in a directory from where it can directly be imported 
                      (i.e. the directory should be in the PYTHON_PATH environment variable or
                      it should be in the current working directory)
                      and the file name should be 'model_xxx.py', where xxx stands for the string
                      that should be specified in this variable
    
    binary          : bool, optional
                      If True input files will be written in binary format, if False input files are
                      written as formatted ascii text. 

    writeDustTemp   : bool, optional
                      If True a separate dust_temperature.inp/dust_tempearture.binp file will be
                      written under the condition that the model contains a function getDustTemperature() 
        
    old             : bool, optional
                      If set to True the input files for the old 2D version of radmc will be created

    dfunc           : function, optional
                      Decision function for octree-like amr tree building. It should take linear arrays of 
                      cell centre coordinates (x,y,z) and cell half-widhts (dx,dy,dz) in all three dimensions,
                      a radmc3d model, a dictionary with all parameters from problem_params.inp and an other 
                      keyword argument (**kwargs). It should return a boolean ndarray of the same length as 
                      the input coordinates containing True if the cell should be resolved and False if not. 
                      An example for the implementation of such decision function can be found in radmc3dPy.analyze
                      module (radmc3dPy.analyze.gdensMinMax()).

    dfpar           : dictionary
                      Dicionary of keyword arguments to be passed on to dfunc. These parameters will not be written
                      to problem_params.inp. Parameters can also be passed to dfunc via normal keyword arguments 
                      gathered in **kwargs, however all keyword arguments in **kwargs will be written to problem_params.inp

    **kwargs        : Any varible name in problem_params.inp can be used as a keyword argument.
                      At first all variables are read from problem_params.in to a dictionary called ppar. Then 
                      if there is any keyword argument set in the call of problem_setup_dust the ppar dictionary 
                      is searched for this key. If found the value belonging to that key in the ppar dictionary 
                      is changed to the value of the keyword argument. If no such key is found then the dictionary 
                      is simply extended by the keyword argument. Finally the problem_params.inp file is updated
                      with the new parameter values.
 
      
    Notes
    -----

    Files written by problemSetupDust() for RADMC-3D
        
        * dustopac.inp             : Dust opacity master file.
        
        * wavelength_micron.inp    : Wavelength grid.
        
        * amr_grid.inp             : Spatial grid.
        
        * stars.inp                : Input radiation field (discrete stellar sources).
        
        * stellarsrc_density.inp   : Input radiation field (continuous stellar sources).
        
        * stellarsrc_templates.inp : Input radiation field (continuous stellar sources).
        
        * dust_density.inp         : Dust density distribution.
        
        * radmc3d.inp              : Parameters for RADMC-3D (e.g. Nr of photons to be used, scattering type, etc).

    """
  
    # Read the parameters from the problem_params.inp file 
    modpar = analyze.readParams()

    # Make a local copy of the ppar dictionary
    ppar = modpar.ppar


    if model=='':
        print 'ERROR'
        print 'No model name is given'
        return

    if not ppar:
        print 'problem_params.inp was not found'
        return 

    if ppar['grid_style'] != 0:
        if old:
            print 'ERROR'
            print 'problemSetupDust was called with the old switch, meaning to create a model setup'
            print 'for the predecessor radmc code, and with the AMR activated'
            print 'radmc does not support mesh refinement'
            return

# --------------------------------------------------------------------------------------------
# If there is any additional keyword argument (**kwargs) then check
#   if there is such key in the ppar dictionary and if is change its value that of
#   the keyword argument. If there is no such key in the ppar dictionary then add the keyword
#   to the dictionary
# --------------------------------------------------------------------------------------------
    if binary:
        modpar.setPar(['rto_style', '3', '', ''])

    if kwargs:
        for ikey in kwargs.keys():
            modpar.ppar[ikey] = kwargs[ikey]
            
            if type(kwargs[ikey]) is float:
                modpar.setPar([ikey, ("%.7e"%kwargs[ikey]), '', ''])
            elif type(kwargs[ikey]) is int:
                modpar.setPar([ikey, ("%d"%kwargs[ikey]), '', ''])
            elif type(kwargs[ikey]) is str:
                modpar.setPar([ikey, kwargs[ikey], '', ''])
            elif type(kwargs[ikey]) is list:
                dum = '['
                for i in range(len(kwargs[ikey])):
                    if type(kwargs[ikey][i]) is float:
                        dum = dum + ("%.7e"%kwargs[ikey][i])
                    elif type(kwargs[ikey][i]) is int:
                        dum = dum + ("%d"%kwargs[ikey][i])
                    elif type(kwargs[ikey][i]) is str:
                        dum = dum + (kwargs[ikey][i])
                    else:
                        print ' ERROR '
                        print ' Unknown data type in '+ikey
                        print kwargs[ikey][i]

                    if i<len(kwargs[ikey])-1:
                        dum = dum + ', '
                dum = dum + (']') 
                modpar.setPar([ikey, dum, '', ''])

        modpar.writeParfile()
        ppar = modpar.ppar

# --------------------------------------------------------------------------------------------
# Create the grid
# --------------------------------------------------------------------------------------------
    #
    # Check if AMR is activated or not
    #
    if ppar['grid_style'] == 1:
        grid = analyze.radmc3dOctree()
        
        # Pass all parameters from dfpar to ppar
        if dfpar is not None:
            for ikey in dfpar.keys():
                ppar[ikey] = dfpar[ikey]

        # Spatial grid
        grid.makeSpatialGrid(ppar=ppar, dfunc=dfunc, model=model, **kwargs)
    else:
        grid = analyze.radmc3dGrid()
        # Spatial grid
        grid.makeSpatialGrid(ppar=ppar)
    # Wavelength grid
    grid.makeWavelengthGrid(ppar=ppar)

# --------------------------------------------------------------------------------------------
# Dust opacity
# --------------------------------------------------------------------------------------------
    if ppar.has_key('dustkappa_ext'):
        opac=analyze.radmc3dDustOpac()
        #Master dust opacity file
        opac.writeMasterOpac(ext=ppar['dustkappa_ext'], scattering_mode_max=ppar['scattering_mode_max'], old=old)
        if old:
            #Frequency grid
            grid.writeWavelengthGrid(old=old)
            # Create the dust opacity files
            opac.makeopacRadmc2D(ext=ppar['dustkappa_ext'])
    else:
        #if old:
            #print 'ERROR'
            #print 'Calculating dust opacities for radmc (2D version) is not yet implemented)'
        opac=analyze.radmc3dDustOpac()
        # Calculate the opacities and write the master opacity file
        opac.makeOpac(ppar=ppar,old=old)

# --------------------------------------------------------------------------------------------
# Try to get the specified model
# --------------------------------------------------------------------------------------------
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
# --------------------------------------------------------------------------------------------
# Create the input radiation field (stars at this point) 
# --------------------------------------------------------------------------------------------

    radSources = analyze.radmc3dRadSources(ppar=ppar, grid=grid)
    radSources.getStarSpectrum(tstar=ppar['tstar'], rstar=ppar['rstar'])


    # Check if the model has functions to set up continuous starlike sources 
    if ppar.has_key('incl_cont_stellarsrc'):
        stellarsrcEnabled = ppar['incl_cont_stellarsrc']
       
        if stellarsrcEnabled:
            if dir(mdl).__contains__('getStellarsrcDensity'):
                if callable(getattr(mdl, 'getStellarsrcDensity')):
                    if dir(mdl).__contains__('getStellarsrcTemplates'):
                        if callable(getattr(mdl, 'getStellarsrcTemplates')):
                            stellarsrcEnabled = True
                        else: 
                            stellarsrcEnabled = False
                    else: 
                        stellarsrcEnabled = False
                else: 
                    stellarsrcEnabled = False
            else: 
                stellarsrcEnabled = False
   
    else:
        stellarsrcEnabled = False

    if stellarsrcEnabled:
        dum = mdl.getStellarsrcTemplates(grid=grid, ppar=ppar)
        if dum[0,0]<=0.:
            radSources.csntemplate = dum.shape[0]
            radSources.cstemp = []
            radSources.cstemptype = 1 
            radSources.cststar = dum[:,0]
            radSources.csrstar = dum[:,1]
            radSources.csmstar = dum[:,2]
        else:
            radSources.csntemplate = dum.shape[0]
            radSources.cstemp = dum
            radSources.cstemptype = 2 
            radSources.cststar = []
            radSources.csrstar = []
            radSources.csmstar = []


        radSources.csdens = mdl.getStellarsrcDensity(grid=grid, ppar=ppar)
# --------------------------------------------------------------------------------------------
# Create the dust density distribution 
# --------------------------------------------------------------------------------------------
    data = analyze.radmc3dData(grid)
    if dir(mdl).__contains__('getDustDensity'):
        if callable(getattr(mdl, 'getDustDensity')):
            data.rhodust = mdl.getDustDensity(grid=grid, ppar=ppar)
        else:
            print 'WARNING'
            print ' '+model+'.py does not contain a getDustDensity() function, therefore, '
            print ' dust_density.inp cannot be written'
            return 
    else:
        print 'WARNING'
        print ' '+model+'.py does not contain a getDustDensity() function, therefore, '
        print ' dust_density.inp cannot be written'
        return 
# --------------------------------------------------------------------------------------------
# Create the dust temperature distribution if the model has such function
# --------------------------------------------------------------------------------------------
    if writeDustTemp:
        if dir(mdl).__contains__('getDustTemperature'):
            if callable(getattr(mdl, 'getDustTemperature')):
                data.dusttemp = mdl.getDustTemperature(grid=grid, ppar=ppar)
            else:
                print 'WARNING'
                print ' '+model+'.py does not contain a getDustTemperature() function, therefore, '
                print ' dust_temperature.dat cannot be written'
                return 
        else:
            print 'WARNING'
            print ' '+model+'.py does not contain a getDustTemperature() function, therefore, '
            print ' dust_temperature.dat cannot be written'
            return 
    #data.rhodust = mdl.get_temperature(grid=grid, ppar=ppar) * ppar['dusttogas']
# --------------------------------------------------------------------------------------------
# Now write out everything 
# --------------------------------------------------------------------------------------------

    if ppar['grid_style'] == 1:
        #Frequency grid
        grid.writeWavelengthGrid()
        #Spatial grid
        grid.writeSpatialGrid()

    else:
        #Frequency grid
        grid.writeWavelengthGrid(old=old)
        #Spatial grid
        grid.writeSpatialGrid(old=old)
    #Input radiation field
    radSources.writeStarsinp(ppar=ppar, old=old)
    
    # Continuous starlike sources
    if stellarsrcEnabled:
        radSources.writeStellarsrcTemplates()
        radSources.writeStellarsrcDensity(binary=binary)

    #totlum = radSources.getTotalLuminosities()
    print '-------------------------------------------------------------'
    print 'Luminosities of radiation sources in the model :'
    
    totlum = radSources.getTotalLuminosities(readInput=True)
    print 'As calculated from the input files :'
    print 'Stars : '
    print ("  Star #%d + hotspot        : %.6e"%(0, totlum['lnu_star'][0]))
    for istar in range(1,radSources.nstar):
        print ("  Star #%d               : %.6e"%(istar, totlum['lnu_star'][istar]))
    print ("Continuous starlike source : %.6e"%totlum['lnu_accdisk'])
    print ' '
    print '-------------------------------------------------------------'

    #Dust density distribution
    if ppar['grid_style'] == 1:
        data.writeDustDens(binary=binary, octree=True)
    else:
        data.writeDustDens(binary=binary, old=old)

    #Dust temperature distribution
    if writeDustTemp:
        if ppar['grid_style'] == 1:
            data.writeDustTemp(binary=binary, octree=True)
        else:
            data.writeDustTemp(binary=binary)
    #radmc3d.inp
    if not old:
        writeRadmc3dInp(modpar=modpar)
    else:
        writeRadmcInp(modpar=modpar)


# --------------------------------------------------------------------------------------------
# Calculate optical depth for diagnostics purposes
# --------------------------------------------------------------------------------------------
    #radSources = radmc3dRadSources(ppar=ppar, grid=grid)
    #radSources.readStarsinp()
    #pwav = radSources.findPeakStarspec()[0]
    #data.getTau(wav=pwav, usedkappa=False)
    #print 'Radial optical depth at '+("%.2f"%pwav)+'um : ', data.taux.max()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def problemSetupGas(model='', fullsetup=False, binary=True,  writeGasTemp=False, dfunc=None, dfpar=None, **kwargs):
    """
    Function to set up a gas model for RADMC-3D 
    
    Parameters
    ----------

    model           : str
                      Name of the model that should be used to create the density structure
                      the file should be in a directory from where it can directly be imported 
                      (i.e. the directory should be in the PYTHON_PATH environment variable, or
                      it should be the current working directory)
                      and the file name should be 'MODELNAME.py', where MODELNAME stands for the string
                      that should be specified in this variable

    fullsetup       : bool, optional
                      If False only the files related to the gas simulation is written out
                      (i.e. no grid, stellar parameter file and radmc3d master command file is written)
                      assuming that these files have already been created for a previous continuum simulation.
                      If True the spatial and wavelength grid as well as the input radiation field
                      and the radmc3d master command file will be (over)written. 

    binary          : bool, optional
                      If True input files will be written in binary format, if False input files are
                      written as formatted ascii text. 

    writeGasTemp    : bool, optional
                      If True a separate gas_temperature.inp/gas_tempearture.binp file will be
                      written under the condition that the model contains a function get_gas_temperature() 
    
    dfunc           : function, optional
                      Decision function for octree-like amr tree building. It should take linear arrays of 
                      cell centre coordinates (x,y,z) and cell half-widhts (dx,dy,dz) in all three dimensions,
                      a radmc3d model, a dictionary with all parameters from problem_params.inp and an other 
                      keyword argument (**kwargs). It should return a boolean ndarray of the same length as 
                      the input coordinates containing True if the cell should be resolved and False if not. 
                      An example for the implementation of such decision function can be found in radmc3dPy.analyze
                      module (radmc3dPy.analyze.gdensMinMax()). 

    dfpar           : dictionary
                      Dicionary of keyword arguments to be passed on to dfunc. These parameters will not be written
                      to problem_params.inp. Parameters can also be passed to dfunc via normal keyword arguments 
                      gathered in **kwargs, however all keyword arguments in **kwargs will be written to problem_params.inp

    **kwargs        : Any varible name in problem_params.inp can be used as a keyword argument.
                      At first all variables are read from problem_params.in to a dictionary called ppar. Then 
                      if there is any keyword argument set in the call of problem_setup_gas the ppar dictionary 
                      is searched for such key. If found the value belonging to that key in the ppar dictionary 
                      is changed to the value of the keyword argument. If no such key is found then the dictionary 
                      is simply extended by the keyword argument. Finally the problem_params.inp file is updated
                      with the new parameter values.
                      Any additional keyword argument for the octree AMR mesh generation should also be passed here.

       
    Notes
    -----
       
    Files written by problemSetupGas()
        
        
        * lines.inp             : Line mode master command file.
        
        * numberdens_xxx.inp    : Number density of molecule/atomic species 'xxx'
        
        * gas_velocity.inp      : Gas velocity
        
        * microturbulence.inp   : The standard deviation of the Gaussian line profile caused by turbulent 
                                broadening.
        
        * gas_temperature.inp   : Gas temperature (which may be different from the dust temperature). If
                                tgas_eq_tdust is set to zero in radmc3d.inp the gas temperature in this
                                file will be used instead of the dust temperature. 

        If fullsetup is set to True the following additional files will be created

        * amr_grid.inp          : Spatial grid.
        
        * wavelength_micron.inp : Wavelength grid.
        
        * stars.inp             : Input radiation field.
        
        * radmc3d.inp           : Parameters for RADMC-3D (e.g. Nr of photons to be used, scattering type, etc).
        
    """
    # Read the parameters from the problem_params.inp file 
    modpar = analyze.readParams()

    # Make a local copy of the ppar dictionary
    ppar = modpar.ppar
    

    if not ppar:
        print 'ERROR'
        print 'problem_params.inp was not found'
        return
    
# --------------------------------------------------------------------------------------------
# If there is any additional keyword argument (**kwargs) then check
#   if there is such key in the ppar dictionary and if is change its value that of
#   the keyword argument. If there is no such key in the ppar dictionary then add the keyword
#   to the dictionary
# --------------------------------------------------------------------------------------------
    if binary:
        modpar.setPar(['rto_style', '3', '', ''])

    if kwargs:
        for ikey in kwargs.keys():
            modpar.ppar[ikey] = kwargs[ikey]
            
            if type(kwargs[ikey]) is float:
                modpar.setPar([ikey, ("%.7e"%kwargs[ikey]), '', ''])
            elif type(kwargs[ikey]) is int:
                modpar.setPar([ikey, ("%d"%kwargs[ikey]), '', ''])
            elif type(kwargs[ikey]) is str:
                modpar.setPar([ikey, kwargs[ikey], '', ''])
            elif type(kwargs[ikey]) is list:
                dum = '['
                for i in range(len(kwargs[ikey])):
                    if type(kwargs[ikey][i]) is float:
                        dum = dum + ("%.7e"%kwargs[ikey][i])
                    elif type(kwargs[ikey][i]) is int:
                        dum = dum + ("%d"%kwargs[ikey][i])
                    elif type(kwargs[ikey][i]) is str:
                        dum = dum + (kwargs[ikey][i])
                    else:
                        print ' ERROR '
                        print ' Unknown data type in '+ikey
                        print kwargs[ikey][i]

                    if i<len(kwargs[ikey])-1:
                        dum = dum + ', '
                dum = dum + (']') 
                modpar.setPar([ikey, dum, '', ''])

        modpar.writeParfile()
        ppar = modpar.ppar
            
            
# --------------------------------------------------------------------------------------------
# Try to get the specified model
# --------------------------------------------------------------------------------------------
    try:
        import os
        imp_path = os.getcwd()
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

# --------------------------------------------------------------------------------------------
# If the current working directory is empty (i.e. no dust setup is present) then
#   we must make a complete setup and dump the spatial and wavelength grids as well
#   as the parameters in the radmc3d.inp file
# --------------------------------------------------------------------------------------------
    if fullsetup:

# --------------------------------------------------------------------------------------------
# Create the grid
# --------------------------------------------------------------------------------------------
        #
        # Check if AMR is activated or not
        #
        if ppar['grid_style'] == 1:
            grid = analyze.radmc3dOctree()
            # Pass all parameters from dfpar to ppar
            if dfpar is not None:
                for ikey in dfpar.keys():
                    ppar[ikey] = dfpar[ikey]
            # Spatial grid
            grid.makeSpatialGrid(ppar=ppar, dfunc=dfunc, model=model, **kwargs)
        else:
            grid = analyze.radmc3dGrid()
            # Spatial grid
            grid.makeSpatialGrid(ppar=ppar)
    
        # Wavelength grid
        grid.makeWavelengthGrid(ppar=ppar)


# --------------------------------------------------------------------------------------------
# Create the input radiation field (stars at this point) 
# --------------------------------------------------------------------------------------------

        radSources = analyze.radmc3dRadSources(ppar=ppar, grid=grid)
        radSources.getStarSpectrum(tstar=ppar['tstar'], rstar=ppar['rstar'])
        
        # Check if the model has functions to set up continuous starlike sources
        if ppar.has_key('incl_accretion'):
            stellarsrcEnabled = ppar['incl_accretion']
        else:
            stellarsrcEnabled = False
        if dir(mdl).__contains__('getStellarsrcDensity'):
            if callable(getattr(mdl, 'getStellarsrcDensity')):
                if dir(mdl).__contains__('getStellarsrcTemplates'):
                    if callable(getattr(mdl, 'getStellarsrcTemplates')):
                        stellarsrcEnabled = True
                    else: 
                        stellarsrcEnabled = False
                else: 
                    stellarsrcEnabled = False
            else: 
                stellarsrcEnabled = False
        else: 
            stellarsrcEnabled = False

        if stellarsrcEnabled:
            dum = mdl.getStellarsrcTemplates(grid=grid, ppar=ppar)
            if dum[0,0]<0.:
                radSources.csntemplate = dum.shape[0]
                radSources.cstemp = []
                radSources.cstemptype = 1 
                radSources.cststar = dum[:,0]
                radSources.csrstar = dum[:,1]
                radSources.csmstar = dum[:,2]
            else:
                radSources.csntemplate = dum.shape[0]
                radSources.cstemp = dum
                radSources.cstemptype = 2 
                radSources.cststar = []
                radSources.csrstar = []
                radSources.csmstar = []

            radSources.csdens = mdl.getStellarsrcDensity(grid=grid, ppar=ppar)


# --------------------------------------------------------------------------------------------
# Now write out everything 
# --------------------------------------------------------------------------------------------

        #Frequency grid
        grid.writeWavelengthGrid()
        #Spatial grid
        grid.writeSpatialGrid()
        #Input radiation field
        radSources.writeStarsinp(ppar=ppar)
        # Continuous starlike sources
        if stellarsrcEnabled:
            radSources.writeStellarsrcTemplates()
            radSources.writeStellarsrcDensity(binary=binary)
    
        print '-------------------------------------------------------------'
        print 'Luminosities of radiation sources in the model :'
        
        totlum = radSources.getTotalLuminosities(readInput=True)
        print 'As calculated from the input files :'
        print 'Stars : '
        print ("  Star #%d + hotspot        : %.6e"%(0, totlum['lnu_star'][0]))
        for istar in range(1,radSources.nstar):
            print ("  Star #%d               : %.6e"%(istar, totlum['lnu_star'][istar]))
        print ("Continuous starlike source : %.6e"%totlum['lnu_accdisk'])
        print ' '
        print '-------------------------------------------------------------'

        
        #radmc3d.inp
        writeRadmc3dInp(ppar=ppar)
# --------------------------------------------------------------------------------------------
# If the current working directory contains already a dust setup then we can use the
#   already existing grid files 
# --------------------------------------------------------------------------------------------
    else:
        grid=analyze.readGrid()
    
# --------------------------------------------------------------------------------------------
# Create the gas density distribution 
# --------------------------------------------------------------------------------------------
    # Create the data structure
    data = analyze.radmc3dData(grid)
    # Calculate the gas density and velocity
    # NOTE: the density function in the model sub-modules should provide the gas volume density
    #       in g/cm^3 but RADMC-3D needs the number density in 1/cm^3 so we should convert the
    #       output of the get_density() function to number density using ppar['gasspecMolAbun']
    #       which is the abundance of the gas species with respect to hydrogen divided by the
    #       mean molecular weight
    if dir(mdl).__contains__('getGasDensity'):
        if callable(getattr(mdl, 'getGasDensity')):
            if ppar['grid_style'] == 1:
                data.rhogas = mdl.getGasDensity(grid=grid, ppar=ppar)
    else:
        print 'WARNING'
        print ' '+model+'.py does not contain a getGasDensity() function, therefore, '
        print ' numberdens_***.inp cannot be written'
        return 
       
# --------------------------------------------------------------------------------------------
# Create the molecular abundance
# --------------------------------------------------------------------------------------------
    #if ppar.has_key('gasspec_mol_abun'):
        #data.rhogas = data.rhogas * ppar['gasspec_mol_abun']
    #else:
    if dir(mdl).__contains__('getGasAbundance'):
        if callable(getattr(mdl, 'getGasAbundance')):
            for imol in range(len(ppar['gasspec_mol_name'])):
                gasabun = mdl.getGasAbundance(grid=grid, ppar=ppar, ispec=ppar['gasspec_mol_name'][imol])
                data.ndens_mol = data.rhogas / (2.4 * mp) * gasabun 

                # Write the gas density
                if ppar['grid_style'] == 1:
                    data.writeGasDens(ispec=ppar['gasspec_mol_name'][imol], binary=binary, octree=True)
                else:
                    data.writeGasDens(ispec=ppar['gasspec_mol_name'][imol], binary=binary)

            if abs(ppar['lines_mode'])>2:
                for icp in range(len(ppar['gasspec_colpart_name'])):
                    gasabun = mdl.getGasAbundance(grid=grid, ppar=ppar, ispec=ppar['gasspec_colpart_name'][icp])
                    data.ndens_mol = data.rhogas / (2.4*mp) * gasabun 
                    # Write the gas density
                    data.writeGasDens(ispec=ppar['gasspec_colpart_name'][icp], binary=binary)

    else:
        print 'WARNING'
        print ' '+model+'.py does not contain a getGasAbundance() function, and no "gasspec_mol_abun" '
        print ' parameter is found in the problem_setup.inp file. numberdens_***.inp cannot be written'
        return

# --------------------------------------------------------------------------------------------
# Get the gas velocity field
# --------------------------------------------------------------------------------------------
    if dir(mdl).__contains__('getVelocity'):
        if callable(getattr(mdl, 'getVelocity')):
            data.gasvel = mdl.getVelocity(grid=grid, ppar=ppar)
            # Write the gas velocity
            if ppar['grid_style'] == 1:
                data.writeGasVel(binary=binary, octree=True) 
            else:
                data.writeGasVel(binary=binary) 
    else:
        print 'WARNING'
        print ' '+model+'.py does not contain a getVelocity() function, therefore, '
        print ' gas_velocity.inp cannot be written'
        return
    
# --------------------------------------------------------------------------------------------
# Get the kinetik gas temperature
# --------------------------------------------------------------------------------------------
    # Write the gas temperature if specified 
    if writeGasTemp:
        if dir(mdl).__contains__('getGasTemperature'):
            if callable(getattr(mdl, 'getGasTemperature')):
                data.gastemp = mdl.getGasTemperature(grid=grid, ppar=ppar)
                # Write the gas temperature
                if ppar['grid_style'] == 1:
                    data.writeGasTemp(binary=binary, octree=True) 
                else:
                    data.writeGasTemp(binary=binary) 
        else:
            print 'WARNING'
            print ' '+model+'.py does not contain a getGasTemperature() function, therefore, '
            print ' gas_temperature.inp cannot be written'
            return

# --------------------------------------------------------------------------------------------
# Get the turbulent velocity field
# --------------------------------------------------------------------------------------------
    if dir(mdl).__contains__('getVTurb'):
        if callable(getattr(mdl, 'getVTurb')):
            data.vturb = mdl.getVTurb(grid=grid, ppar=ppar)
            # Write the turbulent velocity field
            if ppar['grid_style'] == 1:
                data.writeVTurb(binary=binary, octree=True) 
            else:
                data.writeVTurb(binary=binary) 
    else:
        data.vturb = np.zeros([grid.nx, grid.ny, grid.nz], dtype=float64)
        data.vturb[:,:,:] = 0.
        data.writeVTurb(binary=binary)
# --------------------------------------------------------------------------------------------
# Write the line RT control file
# --------------------------------------------------------------------------------------------
    # Write the lines.inp the main control file for the line RT
    writeLinesInp(ppar=ppar)
# --------------------------------------------------------------------------------------------------
def writeRadmcInp(modpar=None, nphot=None):
    """Writes the radmc.inp master command file for the 2D version of radmc

    Parameters
    ----------
    ppar   : dictionary
             Contains all parameters of a radmc run.

    nphot  : int
             Number of photons used for the MC simulation
    """

    if nphot==None:
        ppar = modpar.ppar
        nphot = ppar['nphot']

    fname = 'radmc.inp'
    try :
        wfile = open(fname, 'w')
    except:
        print 'Error!' 
        print fname+' cannot be opened!'
        return 

    wfile.write("nphot       =    %d\n"%(nphot))
    wfile.write("iseed       =    -17933201\n")
    wfile.write("imethod     =    2\n")
    wfile.write("ifast       =    0\n")
    wfile.write("enthres     =      0.010000000\n")
    wfile.write("cntdump     =    100000000\n")
    wfile.write("irestart    =        0\n")
    wfile.write("itempdecoup =        1\n")
    wfile.write("iquantum    =        0\n")
    wfile.write("istarsurf   =        1\n")
    wfile.write("nphotdiff   =       15\n")
    wfile.write("errtol      =    1.0000000e-10\n")
    wfile.write("nvstr       =        0\n")
    wfile.write("vserrtol    =     0.0100000\n")
    wfile.write("ivstrt      =        1\n")
    wfile.write("ntemp       =     3000\n")
    wfile.write("temp0       =      0.010000000\n")
    wfile.write("temp1       =        1000000.0\n")
    wfile.close()

    #fname = 'raytrace.inp'
    #try :
        #wfile = open(fname, 'w')
    #except:
        #print 'Error!' 
        #print fname+' cannot be opened!'
        #return 

    #wfile.write("nrphiinf    =       32\n")
    #wfile.write("nrrayextra  =      -20\n")
    #wfile.write("imethod     =        1\n")
    #wfile.write("nrref       =       10\n")
    #wfile.write("dbdr        =        1\n")
    #wfile.write("inclination =       45.0000\n")
    #wfile.close()
# --------------------------------------------------------------------------------------------------
def writeRadmc3dInp(modpar=None):
    """Writes the radmc3d.inp master command file for RADMC-3D

    Parameters
    ----------
    ppar   : dictionary
             Contains all parameters of a RADMC-3D run.

    """

    # Select those keywords, whose block name is 'Code parameters'
    ppar = {}

    for ikey in modpar.pblock.keys():
        if modpar.pblock[ikey]=='Code parameters':
            ppar[ikey] = modpar.ppar[ikey]

    print 'Writing radmc3d.inp'

    wfile = open('radmc3d.inp', 'w')
    keys = ppar.keys()
    keys.sort()
    for key in keys:
        #if modpar.ppar.has_key(key):
        wfile.write('%s %d\n'%(key+' =',ppar[key]))

    wfile.close()

# --------------------------------------------------------------------------------------------------
def writeLinesInp(ppar=None):
    """Writes the lines.inp master command file for line simulation in RADMC-3D

    Parameters
    ----------
    
    ppar   : dictionary,
            Contains all parameters of a RADMC-3D run 

    """

    # Do a consistency check
    n1 = len(ppar['gasspec_mol_name'])
    n2 = len(ppar['gasspec_mol_abun'])
    n3 = len(ppar['gasspec_mol_dbase_type'])

    if ((n1!=n2)|(n2!=n3)):
        print ' ERROR '
        print ' gasspec_mol_name, gasspec_mol_abun and gasspec_mol_dbase_type have different number of elements'
        return 

    if ppar.has_key('gasspec_colpart_name') & ppar.has_key('gasspec_colpart_abun'):
        n4 = len(ppar['gasspec_colpart_name'])
        n5 = len(ppar['gasspec_colpart_abun'])
    else:
        n4 = 0
        n5 = 0

    if (n4!=n5):
        print ' ERROR '
        print ' gasspec_colpart_name and gasspec_colpart_abun have different number of elements'
        return 


    print 'Writing lines.inp'
    wfile = open('lines.inp', 'w')
    # File format
    wfile.write("%d\n"%ppar['lines_mode'])
    # Nr of gas species
    wfile.write("%d\n"%n1)
    # Gas species name and database type

    if abs(ppar['lines_mode'])<=2:
        for imol in range(n1):
            wfile.write("%s %s %d %d %d\n"%(ppar['gasspec_mol_name'][imol], ppar['gasspec_mol_dbase_type'][imol], 0, 0, 0))
    else:
        for imol in range(n1):
            wfile.write("%s %s %d %d %d\n"%(ppar['gasspec_mol_name'][imol], ppar['gasspec_mol_dbase_type'][imol], 0, 0, n4))

        if n4>0:
            for icp in range(n4):
                wfile.write("%s\n"%ppar['gasspec_colpart_name'][icp])
        else:
            print ' ERROR'
            print ' An NLTE line excitation method is selected (lines_mode='+("%d"%ppar["lines_mode"])+'), but no collisional'
            print ' partner is given in the parameter file. '
            wfile.close()
            return

    wfile.close()

# --------------------------------------------------------------------------------------------------
def validateModel(model='', dustModel=False, gasModel=False, writeDustTemp=False, octree=False):
    """
    Function to validate a model. It checks three things: 1) whether or not the model can be imported,
    2) whether the model has all the function to be used as dust and/or gas model, 3) if it has the right
    number of arguments. The function names tested are getDefaultParams, getDustDensity, getGasDensity, 
    getGasAbundance, getVTurb, getVelocity, getDustTempearture (optional).

    Parameters
    ----------

    model       : str
                  Name of the model to be tested

    dustModel   : bool
                  If True the existence of functions getDustDensity() and getDustTemperature() will be checked.
                  The latter is only checked if writeDustTemp is set to True.

    gasModel    : bool
                  If True the existence of functions getGasDensity(), getGasAbundance(), getVTurb(), getVelocity()
                  will be checked.

    writeDustTemp: bool
                   If True the existence of the function getDustTemperature() will be checked.

    octree      : bool
                  If True the number of argument of the model functions will be checked. For regular grids only two 
                  arguments should be present for the grid instance and for the parameter dictionary (grid, ppar). 
                  For a model to be used with octree AMR three additional arguments for the three spatial coordiantes
                  (x,y,z) should be present. The argument sequence should then be x, y, z, grid, ppar.

    Returns
    -------
    A boolean True if the model is valid and False if it is not. 

    """

    #
    # First check if the model can be imported
    #
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


    isValid = True
    #
    # Now check the function names in the model
    #
    fnamelist = [f[0] for f in inspect.getmembers(mdl) if inspect.isfunction(f[1])]
    
    if 'getDefaultParams' not in fnamelist:
        print 'ERROR'
        print model + ' does not contain a function to provide default parameters (getDefaultParams)'
        isValid = False

    if dustModel:
        if 'getDustDensity' not in fnamelist:
            print 'ERROR'
            print model + ' does not contain a function to provide dust density (getDustDensity)'
            isValid = False
           
        if writeDustTemp:
            if 'getDustTemperature' not in fnamelist:
                print 'ERROR'
                print model + ' does not contain a function to provide dust temperature (getDustTemperature)'
                print 'yet the setup function has been called with the option to write the dust temperature.'
                isValid = False

    if gasModel:
        if 'getGasDensity' not in fnamelist:
            print 'ERROR'
            print model + ' does not contain a function to provide gas density (getGasDensity)'
            isValid = False

        if 'getGasAbundance' not in fnamelist:
            print 'ERROR'
            print model + ' does not contain a function to provide molecular abundance (getGasAbundance)'
            isValid = False
        
        if 'getVTurb' not in fnamelist:
            print 'ERROR'
            print model + ' does not contain a function to provide turbulent velocity (getVTurb)'
            isValid = False

        if 'getVelocity' not in fnamelist:
            print 'ERROR'
            print model + ' does not contain a function to provide gas velocity (getVelocity)'
            isValid = False

    #
    # Check the number of arguments
    #
    if dustModel:
        arglist = inspect.getargspec(mdl.getDustDensity).args
        argnames = ''
        if len(arglist)>0:
            argnames = arglist[0]
            for iarg in arglist[1:]:
                argnames += ', '+iarg
        if octree:
            if len(arglist)<5:
                print 'ERROR'
                print model+'.getDustDensity() has only '+("%d"%len(arglist))+' arguments : ', argnames
                print 'To use octree the argument list should be :'
                print 'x=None, y=None, z=None, grid=None, ppar=None)'
                isValid = False
        else:
            if len(arglist)<2:
                print 'ERROR'
                print model+'.getDustDensity() has only '+("%d"%len(arglist))+' arguments : ', argnames
                print 'The minimal argument list of a model function should be :'
                print 'grid=None, ppar=None)'
                isValid = False


        if writeDustTemp:
            arglist = inspect.getargspec(mdl.getDustTemperature).args
            argnames = ''
            if len(arglist)>0:
                argnames = arglist[0]
                for iarg in arglist[1:]:
                    argnames += ', '+iarg
            if octree:
                if len(arglist)<5:
                    print 'ERROR'
                    print model+'.getDustTemperature() has only '+("%d"%len(arglist))+' arguments : ', argnames
                    print 'To use octree the argument list should be :'
                    print 'x=None, y=None, z=None, grid=None, ppar=None)'
                    isValid = False
            else:
                if len(arglist)<2:
                    print 'ERROR'
                    print model+'.getDustTemperature() has only '+("%d"%len(arglist))+' arguments : ', argnames
                    print 'The minimal argument list of a model function should be :'
                    print 'grid=None, ppar=None)'
                    isValid = False

    if gasModel:
        arglist = inspect.getargspec(mdl.getGasDensity).args
        argnames = ''
        if len(arglist)>0:
            argnames = arglist[0]
            for iarg in arglist[1:]:
                argnames += ', '+iarg
        if octree:
            if len(arglist)<5:
                print 'ERROR'
                print model+'.getGasDensity() has only '+("%d"%len(arglist))+' arguments : ', argnames
                print 'To use octree the argument list should be :'
                print 'x=None, y=None, z=None, grid=None, ppar=None)'
                isValid = False
        else:
            if len(arglist)<2:
                print 'ERROR'
                print model+'.getGasDensity() has only '+("%d"%len(arglist))+' arguments : ', argnames
                print 'The minimal argument list of a model function should be :'
                print 'grid=None, ppar=None)'
                isValid = False


        arglist = inspect.getargspec(mdl.getGasAbundance).args
        argnames = ''
        if len(arglist)>0:
            argnames = arglist[0]
            for iarg in arglist[1:]:
                argnames += ', '+iarg
        if octree:
            if len(arglist)<5:
                print 'ERROR'
                print model+'.getGasAbundance() has only '+("%d"%len(arglist))+' arguments : ', argnames
                print 'To use octree the argument list should be :'
                print 'x=None, y=None, z=None, grid=None, ppar=None)'
                isValid = False
        else:
            if len(arglist)<2:
                print 'ERROR'
                print model+'.getGasAbundance() has only '+("%d"%len(arglist))+' arguments : ', argnames
                print 'The minimal argument list of a model function should be :'
                print 'grid=None, ppar=None)'
                isValid = False

        arglist = inspect.getargspec(mdl.getVTurb).args
        argnames = ''
        if len(arglist)>0:
            argnames = arglist[0]
            for iarg in arglist[1:]:
                argnames += ', '+iarg
        if octree:
            if len(arglist)<5:
                print 'ERROR'
                print model+'.getVTurb() has only '+("%d"%len(arglist))+' arguments : ', argnames
                print 'To use octree the argument list should be :'
                print 'x=None, y=None, z=None, grid=None, ppar=None)'
                isValid = False
        else:
            if len(arglist)<2:
                print 'ERROR'
                print model+'.getVTurb() has only '+("%d"%len(arglist))+' arguments : ', argnames
                print 'The minimal argument list of a model function should be :'
                print 'grid=None, ppar=None)'
                isValid = False


        arglist = inspect.getargspec(mdl.getVelocity).args
        argnames = ''
        if len(arglist)>0:
            argnames = arglist[0]
            for iarg in arglist[1:]:
                argnames += ', '+iarg
        if octree:
            if len(arglist)<5:
                print 'ERROR'
                print model+'.getVelocity() has only '+("%d"%len(arglist))+' arguments : ', argnames
                print 'To use octree the argument list should be :'
                print 'x=None, y=None, z=None, grid=None, ppar=None)'
                isValid = False
        else:
            if len(arglist)<2:
                print 'ERROR'
                print model+'.getVelocity() has only '+("%d"%len(arglist))+' arguments : ', argnames
                print 'The minimal argument list of a model function should be :'
                print 'grid=None, ppar=None)'
                isValid = False

    
    #
    # If it passed all tests so far then formally the model should be OK. There is no guarantee, though
    # that it will work properly. 
    #
    return isValid
