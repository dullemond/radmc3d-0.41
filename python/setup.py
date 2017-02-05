#!/usr/bin/env python

from distutils.core import setup
import os, sys
from subprocess import Popen, PIPE

def findFiles(src_dir, *wildcards):

    src_dir = src_dir.strip()
    while src_dir[-1]=='/':
        src_dir = src_dir[:-1]


    # Find all directory names
    dirList = Popen(['find '+src_dir+' -name "*"'], shell=True, \
            stdout=PIPE).communicate()[0].split()
   

    foundList = []
    for i in range(len(dirList)):
        #if (os.path.isdir(dirList[i])&(dirList[i].strip()!=src_dir)):
        if os.path.isdir(dirList[i]):
            # Find the appropriate files within each directory
            fileList = []
            for wc in wildcards:
                dum = Popen(['ls -1 '+dirList[i]+'/'+wc], shell=True, \
                        stdout=PIPE, stderr=PIPE).communicate()[0].split()
                #dum = Popen(['find '+dirList[i]+' -name "'+wc+'"'], shell=True, \
                        #stdout=PIPE).communicate()[0].split()

                if len(dum)>0:
                    fileList.extend(dum)

            if len(fileList)>0:
                foundList.append((dirList[i], fileList))

    return foundList    
   
fileList = findFiles('./radmc3dPy', '*.py')

python_files = []
for i in range(len(fileList)):
    for j in range(len(fileList[i][1])):
        python_files.append(fileList[i][1][j])

moduleNames = []
packageNames = []
for i in range(len(python_files)):

    ind1 = python_files[i].strip()[::-1].find('/')
    dum  = python_files[i].strip()[-ind1:-3]

    if dum.strip()!='__init__':
        moduleNames.append('radmc3dPy.'+dum)
    else:
        sdum = python_files[i].split('/')[1:-1]

        txt = sdum[0]
        if len(sdum)>1:
            for imod in range(1,len(sdum)):
                txt += '.'+sdum[imod]

        packageNames.append(txt)


setup(name='radmc3dPy',
        version='0.28.2',
      description='Python module for RADMC3D',
      author='Attila Juhasz',
      author_email='juhasz@ast.cam.ac.uk',
      packages=packageNames)

