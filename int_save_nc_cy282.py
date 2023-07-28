# new version of the code!!!
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Concatenates cubes and then plots files that show the time evolution of variables.
"""

import numpy as np
import iris
from glob import glob
import os

stashcode = 'm01s34i103' #'m01s34i103' 
rose = 'u-cy282'
header='cy282a.'
file_chunk='pm'
months = ['2014jan', '2014feb', '2014mar', '2014apr', '2014may', '2014jun',
          '2014jul', '2014aug', '2014sep', '2014oct', '2014nov', '2014dec']

files_directory_UKCA = '/jet/home/ding0928/cylc-run/' + rose + '/share/data/History_Data/'
rosefolder = '/ocean/projects/atm200005p/ding0928/nc_file_full/' + rose + '/'
ncfolder = rosefolder + 'small_nc_files/'

def make_directories(nameofdir):
    newdir = os.path.join(files_directory_UKCA, nameofdir)
    try:
        os.mkdir(newdir)
    except OSError as e:
        if e.errno == 17:
            print('Folder {} already exists'.format(newdir))
        else:
            raise OSError
    return 'Created Folder: {}'.format(newdir)

def save_small_nc_files(bigarray, ncfolder, rose, stashcode, month):
    print('Begin Saving')
    print('Save Location: {}'.format(ncfolder))
    i = 0
    for cube in bigarray[0]:
        saving_name = ncfolder + rose + '_' + stashcode + '_' + month + '.nc'
        print('saving', saving_name)
        iris.save(cube, saving_name, netcdf_format="NETCDF4")
        i += 1
    return 'Saving Complete'

for month in months:
    # Update the file path and filename generation with 'month' instead of 'iday'
    files_directory = files_directory_UKCA
    # print('test',files_directory + header+ file_chunk+ month[-3:])
    pp_files1 = [sorted(glob(files_directory + header+ file_chunk+'*' + month[-3:]))]
    print('pp_files1', pp_files1)
    print('pp_files1[0]_test', pp_files1[0])

    rosefolder = '/ocean/projects/atm200005p/ding0928/nc_file_full/' + rose + '/'
    ncfolder = rosefolder + 'small_nc_files/'

    tacc3hr = 0

    stashconstr = iris.AttributeConstraint(STASH=stashcode)
    cubes = iris.load((pp_files1[0])[0], stashconstr)
    print('loading cubes ' + str((pp_files1[0])[0]))

    if len(pp_files1) > 1:
        for filelist in pp_files1[1:]:
            print('loading cubes ' + str(filelist[0]))
            cubel = iris.load(filelist[0], stashconstr)
            for cube1 in cubel:
                cubes.append(cube1)
    print('cubes', cubes)

    coord_names = [coord.name() for coord in cubes[0].coords()]
    print(coord_names)

    bigarray = []
    timepointslist = []
    print('Begin Cube Data Processing')
    for cube in cubes:
        cl = iris.cube.CubeList()
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

    if len(pp_files1[0]) > 1:
        fileindex = 1
        for step_file in (pp_files1[0])[1:]:
            morecubes = iris.load(step_file, stashconstr)
            print('loading cubes ' + str(step_file))
            if len(pp_files1) > 1:
                for filelist in pp_files1[1:]:
                    print('loading cubes ' + str(filelist[fileindex]))
                    morecubel = iris.load(filelist[fileindex], stashconstr)
                    for morecube in morecubel:
                        morecubes.append(morecube)
            i = 0
            for cube in morecubes:
                print(cube.coord('time').points)
                for j in timepointslist[i]:
                    it = 0
                    for k in cube.coord('time').points:
                        if k == j:
                            print('removing', j)
                            if it == 0:
                                cube = cube[1:, ...]
                            elif it == 1 and len(cube.coord('time').points) == 2:
                                cube = cube[0:1, ...]
                            else:
                                raise ValueError('Overlapping fudged removal failure')
                            print(cube)
                        it += 1

                cube.remove_coord('forecast_reference_time')
                try:
                    cube.remove_coord('surface_altitude')
                except Exception:
                    pass
                if tacc3hr == 1:
                    cube.coord('time').bounds = None
                    iris.util.promote_aux_coord_to_dim_coord(cube, 'time')

                if cube.name() == (bigarray[i])[0].name() and cube.attributes['STASH'] == (bigarray[i])[0].attributes['STASH']:
                    if len(cube.coord('time').points) > 1:
                        for sub_cube in cube.slices_over('time'):
                            bigarray[i].append(sub_cube)
                    else:
                        bigarray[i].append(cube)
                else:
                    for cubelist in bigarray:
                        if cube.name() == cubelist[0].name() and cube.attributes['STASH'] == cubelist[0].attributes['STASH']:
                            if len(cube.coord('time').points) > 1:
                                for sub_cube in cube.slices_over('time'):
                                    bigarray[i].append(sub_cube)
                            else:
                                cubelist.append(cube)
                            break
                for j in cube.coord('time').points:
                    array = timepointslist[i]
                    print(j)
                    timepointslist[i] = np.append(array, j)
                i += 1
            print(timepointslist)
            fileindex += 1

    print(bigarray)
    print(timepointslist)
    print(np.shape(bigarray))
    print(np.shape(timepointslist))

    make_directories(rosefolder)
    make_directories(ncfolder)
    save_small_nc_files(bigarray, ncfolder, rose, stashcode, month)
sys.exit()

#old version of the code, works for u-ct706 version of the model, save a backup copy.
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

stashcode='m01s38i405'  #P m01s00i408
# regional model location+ stashcode
# nuc_num='m01s34i101', ait_num='m01s34i103', acc_num='m01s34i107',cor_num_='m01s34i113',aitin_num='m01s34i119'
# nuc_diam ='m01s38i401', ait_diam ='m01s38i402',acc_diam ='m01s38i403',coar_diam ='m01s38i404',ait_diam_inso = 'm01s38i405'
# theta = 'm01s00i004'
# air_potential_temperature = 'm01s00i408'

# stashcodes = ['m01s34i101', 'm01s34i103', 'm01s34i107', 'm01s34i113', 'm01s34i119']
rose = 'u-cy282'
file_path = "/jet/home/ding0928/python_analysis/Han_connect/"
months = ['pm2014jan', 'pm2014feb', 'pm2014mar', 'pm2014apr', 'pm2014may', 'pm2014jun', 'pm2014jul', 'pm2014aug', 'pm2014sep', 'pm2014oct', 'pm2014nov', 'pm2014dec']
files_directory_UKCA = '/jet/home/ding0928/cylc-run'+rose+'/share/data/History_Data/'


rosefolder = '/ocean/projects/atm200005p/ding0928/nc_file_full/'+rose+'/'
ncfolder = rosefolder + 'small_nc_files/'

# def make_directories(nameofdir):
#     newdir = os.path.join(files_directory_UKCA, nameofdir)
#     try:
#         os.mkdir(newdir)
#     except OSError as e:
#         if e.errno == 17:
#             print('Folder {} already exists'.format(newdir))
#             pass
#         else:
#             raise OSError
#     return('Created Folder: {}'.format(newdir))

# def save_small_nc_files(bigarray, ncfolder, rose, stashcode, month):
#     print('Begin Saving')
#     print('Save Location: {}'.format(ncfolder))
#     i = 0
#     for cube in bigarray[0]:
#         saving_name = ncfolder + rose + '_' + stashcode + '_' + month + '.nc'
#         print('saving', saving_name)
#         iris.save(cube, saving_name, netcdf_format="NETCDF4")
#         i += 1
#     return 'Saving Complete'

# for month in months:
#     # Update the file path and filename generation with 'month' instead of 'iday'
#     files_directory = files_directory_UKCA + month + '/'
#     pp_files1 = [sorted(glob(files_directory + chunk + month[-3:] + '*')) for chunk in ['']]

#     print('pp_files1[0]_test', pp_files1[0])
#     pp_files = pp_files1

# for month in months:
#     # Update the file path and filename generation with 'month' instead of 'iday'
#     files_directory = files_directory_UKCA + month + '/'
#     pp_files1 = [sorted(glob(files_directory + chunk + month[-3:] + '*')) for chunk in ['']]

#     print('pp_files1[0]_test', pp_files1[0])
#     pp_files = pp_files1

#     rosefolder = '/ocean/projects/atm200005p/ding0928/nc_file_full/' + rose + '/'
#     ncfolder = rosefolder + 'small_nc_files/'

#     tacc3hr = 0

#     stashconstr = iris.AttributeConstraint(STASH=stashcode)
#     cubes=iris.load((pp_files[0])[0],stashconstr)
#     print('loading cubes '+str(pp_files[0][0]))


#     if len(pp_files) > 1:
#         for filelist in pp_files[1:]:
#             print('loading cubes '+str(filelist[0]))
#             cubel = iris.load(filelist[0],stashconstr)
#             for cube1 in cubel:
#                 cubes.append(cube1)
#     print('cubes',cubes)

#     coord_names = [coord.name() for coord in cubes[0].coords()]
#     print(coord_names)

#     bigarray=[]
#     timepointslist=[]
#     print('Begin Cube Data Processing')
#     for cube in cubes:
#         cl = iris.cube.CubeList()
#         cube.remove_coord('forecast_reference_time')
#         try:    
#             cube.remove_coord('surface_altitude')
#         except Exception:
#             pass
#         if tacc3hr==1:
#             print(cube.coord('time').points)
#             print(cube.coord('time').bounds)
#             cube.coord('time').bounds = None
#             iris.util.promote_aux_coord_to_dim_coord(cube,'time')

#         timepointslist.append(cube.coord('time').points)
#         time = cube.coord('time')
#         dates = time.units.num2date(time.points)
#         if len(dates) > 1:
#             for sub_cube in cube.slices_over('time'):
#                 print(sub_cube)
#                 cl.append(sub_cube)
#         else:
#             cl.append(cube)
#         bigarray.append(cl)

#     if len(pp_files[0]) > 1:
#         fileindex=1
#         for step_file in (pp_files[0])[1:]:
#             morecubes=iris.load(step_file,stashconstr)        
#             print('loading cubes '+str(step_file))
#             if len(pp_files) > 1:
#                 for filelist in pp_files[1:]:
#                     print('loading cubes '+str(filelist[fileindex]))
#                     morecubel =iris.load(filelist[fileindex],stashconstr)
#                     for morecube in morecubel:
#                         morecubes.append(morecube)
#             i=0
#             for cube in morecubes:
#                 print(cube.coord('time').points)
#                 for j in timepointslist[i]:
#                     it=0
#                     for k in cube.coord('time').points:
#                         if k ==j:
#                             print('removing',j)
#                             if it==0:
#                                 cube = cube[1:,...]
#                             elif it==1 and len(cube.coord('time').points)==2:
#                                 cube = cube[0:1,...]
#                             else:
#                                 raise ValueError('Overlapping fudged removal failure')
#                             print(cube)
#                         it=it+1

#                 cube.remove_coord('forecast_reference_time')
#                 try:
#                     cube.remove_coord('surface_altitude')
#                 except Exception:
#                     pass
#                 if tacc3hr==1:
#                     cube.coord('time').bounds=None
#                     iris.util.promote_aux_coord_to_dim_coord(cube,'time')

#                 if cube.name() ==(bigarray[i])[0].name() and cube.attributes['STASH']==(bigarray[i])[0].attributes['STASH']:
#                     if len(cube.coord('time').points) >1:
#                         for sub_cube in cube.slices_over('time'):
#                             bigarray[i].append(sub_cube)
#                     else:
#                         bigarray[i].append(cube)
#                 else:
#                     for cubelist in bigarray:
#                         if cube.name() ==cubelist[0].name() and cube.attributes['STASH']==cubelist[0].attributes['STASH']:
#                             if len(cube.coord('time').points) >1:
#                                 for sub_cube in cube.slices_over('time'):
#                                     bigarray[i].append(sub_cube)
#                             else:
#                                 cubelist.append(cube)
#                             break
#                 for j in cube.coord('time').points:
#                     array=timepointslist[i]
#                     print(j)
#                     timepointslist[i] = np.append(array,j)
#                 i=i+1
#             print(timepointslist)
#             fileindex=fileindex+1
    
#     print(bigarray)
#     print(timepointslist)
#     print(np.shape(bigarray))
#     print(np.shape(timepointslist))

#     make_directories(rosefolder)
#     make_directories(ncfolder)
#     save_small_nc_files(bigarray, ncfolder, stashcode, timepointslist[0])

