#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Concatenates cubes and then plots files that show the time evolution of variables. Not parallelised, slow and memory-intensive
Code developed by Hamish Gordon,Jesus Vergara Temprado and Kirsty Pringle
h.gordon@leeds.ac.uk
eejvt@leeds.ac.uk
K.Pringle@leeds.ac.uk

Aerosol modellers group
Institute for climate and atmospheric science (ICAS)
University of Leeds 2016

"""

import sys
import numpy as np
#import cartopy.crs as ccrs
import time
import iris
from glob import glob
#import matplotlib.pyplot as plt
#import matplotlib.ticker as mticker
#import iris.plot as iplt
#from iris.coord_systems import RotatedGeogCS, GeogCS
#from iris.analysis.cartography import unrotate_pole
#from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
#import iris.quickplot as qplt
#from iris.experimental.equalise_cubes import equalise_attributes
import datetime
#from matplotlib.backends.backend_pdf import PdfPages
from scipy.io import netcdf
import os
#import pdb
#import getpass
#import multiprocessing

rose = 'u-cs093'
stashcode='m01s34i119'
def make_directories(nameofdir):
    newdir = os.path.join(nameofdir)
    print(newdir)
    try:
        os.mkdir(newdir)
    except OSError as e:
        if e.errno == 17:
            print('Folder {} already exists'.format(newdir))
            pass
        else:
            raise OSError
    return('Created Folder: {}'.format(newdir))


def concatenate_nc_files(smallncfolder, bigncfolder, stashcode):
    print('Begin Concatenation')
    print('Save Location: {}'.format(bigncfolder))
    cubes = iris.load(smallncfolder+'*'+stashcode+'*.nc')
    print(cubes)
    cube = cubes[0]
    print(cube)
    time_array = cube.coord('time').points
    print('Concatenated {} files'.format(len(time_array)))
    saving_name = bigncfolder+'Rgn_'+cube.name()+'_'+stashcode+'.nc'
    iris.save(cube,saving_name, netcdf_format="NETCDF4")
    print('File {} Saved'.format(saving_name))
    print('\n')

def delete_smallnc(smallncfolder, stashcode):
    stashfiles = sorted(glob(smallncfolder+'*'+stashcode+'*')) 
    for stashfile in stashfiles:
        os.remove(stashfile)
    print('Deleted {} small nc files for stash {}'.format(len(stashfiles), stashcode))

#rosefolder = '/jet/home/ding0928/python_analysis/Han_connect/nc_flie/'+rose+'/'
rosefolder = '/ocean/projects/atm200005p/ding0928/nc_file_full/'+rose+'/'

bigncfolder = rosefolder+'full_nc_files/'
smallncfolder = rosefolder+'small_nc_files/'


make_directories(bigncfolder)
concatenate_nc_files(smallncfolder, bigncfolder, stashcode)
delete_smallnc(smallncfolder, stashcode)
sys.exit()



