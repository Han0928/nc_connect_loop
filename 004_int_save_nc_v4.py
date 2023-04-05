#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Concatenates cubes and then plots files that show the time evolution of variables.
"""

import sys
import numpy as np
import time
import iris
from glob import glob
import datetime
from scipy.io import netcdf
import os

stashcodes = ['m01s34i101', 'm01s34i103', 'm01s34i107', 'm01s34i113', 'm01s34i119']
file_path = "/jet/home/ding0928/python_analysis/Han_connect/"

latbottom = 35.0
lattop = 45.0  

rose = 'u-ct706'
files_directory_UKCA='/jet/home/ding0928/cylc-run/'+rose+'/share/cycle/'
days=[ str('0720'), str('0722'), str('0724'), str('0726'), str('0728'), str('0730'), str('0801'), str('0803'), str('0805'), str('0807'), str('0809')]
filechunks = ['pb','pc','pd','pe']
run = '20140720T0000Z'


def lat_range(cell):
    return latbottom < cell < lattop

def height_level_range(cell):
    return 0 <= cell <= 40

def make_directories(nameofdir):
    newdir = os.path.join(files_directory_UKCA, nameofdir)
    try:
        os.mkdir(newdir)
    except FileExistsError:
        print('Folder {} already exists'.format(newdir))
    except OSError as e:
        raise OSError
    return('Created Folder: {}'.format(newdir))

def save_small_nc_files(bigarray, ncfolder, stashcodes, timepointslist):
    print('Begin Saving')
    print('Save Location: {}'.format(ncfolder))
    for timepoint in timepointslist:
        for stashcode in stashcodes:
            cubes_to_save = []
            for cubes_list in bigarray:
                for cube in cubes_list:
                    if np.all((cube.attributes['STASH'] == stashcode) & (cube.coord('time').points[0] == timepoint)):
                        cubes_to_save.append(cube)

            if len(cubes_to_save) > 0:
                saving_name = ncfolder+'Rgn_'+str(stashcode)+'_'+str(timepoint)+'.nc'
                print('saving', saving_name)
                iris.save(cubes_to_save, saving_name, netcdf_format="NETCDF4")
    return 'Saving Complete'

for iday in days:
    files_directory=files_directory_UKCA+'2014'+iday+'T0000Z/Regn1/resn_1/RA2M/um/' 
    pp_files1=[sorted(glob(files_directory+'*saa*'+chunk+'*')) for chunk in filechunks]
    pp_files = pp_files1
    date = run[0:8]
    year=date[0:4]

    rosefolder = '/jet/home/ding0928/python_analysis/Han_connect/nc_flie/'+rose+'/'
    ncfolder = rosefolder+'small_nc_files/'

    tacc3hr=0   
    stashconstrs = []
    for stashcode in stashcodes:
        stashconstrs.append(iris.AttributeConstraint(STASH=stashcode))
    print(stashconstrs)

    # Load the cubes
    cubes = []
    for stashconstr in stashconstrs:
        cubes.append(iris.load((pp_files[0])[0], stashconstr))
        if len(pp_files[0]) > 1:
            fileindex = 1
            for step_file in (pp_files[0])[1:]:
                morecubes = [iris.load(step_file, constr) for constr in stashconstrs]

                print('loading cubes '+str(step_file))
                if len(pp_files) > 1:
                    for filelist in pp_files[1:]:
                        print('loading cubes '+str(filelist[fileindex]))
                        morecubel = iris.load(filelist[fileindex], stashconstr)
                        for morecube in morecubel:
                            morecubes.append(morecube)
                i = 0
                for cube in morecubes:
                    if cube and cube[0].attributes['STASH'] == stashcode:
                        cubes[-1].append(cube)
                    i += 1


    # Process the cubes
    bigarray = []
    timepointslist = []
    print('Begin Cube Data Processing')
    for i, stashcode in enumerate(stashcodes):
        stashcode_cubes = cubes[i]
        cl = iris.cube.CubeList()
        for cube in stashcode_cubes:
            cube.remove_coord('forecast_reference_time')
            try:
                cube.remove_coord('surface_altitude')
            except Exception:
                pass
            if tacc3hr == 1:
                print(cube.coord('time').points)
                print(cube.coord('time').bounds)
                cube.coord('time').bounds = None
                iris.util.promote_aux_coord_to_dim_coord(cube, 'time')
            timepointslist.append(cube.coord('time').points)
            time = cube.coord('time')
            dates = time.units.num2date(time.points)
            if len(dates) > 1:
                for sub_cube in cube.slices_over('time'):
                    print(sub_cube)
                    cl.append(sub_cube)
            else:
                cl.append(cube)
        bigarray.append(cl)

    print(bigarray) 
    print(ncfolder)
    make_directories(rosefolder)
    make_directories(ncfolder)
    save_small_nc_files(bigarray, ncfolder, stashcodes, timepointslist)

sys.exit()
