''' This script is used to save the data in the small nc files, but in a simplified way.
    The data is saved in a way that is easier to read and understand.
'''
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
rose = 'u-ct706'

# Define latitudinal range of interest
latbottom = 35.0
lattop = 45.0

# Define the days and filechunks of interest
days = ['0720', '0722', '0724', '0726', '0728', '0730', '0801', '0803', '0805', '0807', '0809']
filechunks = ['pb','pc','pd','pe']
run = '20140720T0000Z'

# Define the function to check if a cell is within the latitudinal range of interest
def lat_range(cell):
    return latbottom < cell < lattop

for iday in days:
    # Define the file directory for each day
    files_directory = '/jet/home/ding0928/cylc-run/'+rose+'/share/cycle/2014'+iday+'T0000Z/Regn1/resn_1/RA2M/um/'

    # Load the cubes for each stashcode
    cubes = []
    # for chunk in filechunks:
    #     for stashcode in stashcodes:
    #         file_name = f"{files_directory}umnsaa_{chunk}006_Rgn_{stashcode}.nc"
    #         try:
    #             cubes.append(iris.load_cube(file_name))
    #         except:
    #             print(f"Unable to load {file_name}")
    for chunk in filechunks:
        for stashcode in stashcodes:
            file_name = f"{files_directory}umnsaa_{chunk}006_Rgn_{stashcode}.nc"
            try:
                cube = iris.load_cube(file_name)
                print(f"Loaded {file_name}")
                cubes.append(cube)
            except:
                print(f"Unable to load {file_name}")


    # Process the cubes and save them in .nc files for each time slot
    for time_point in cubes[0].coord('time').points:
        time_cube_list = []
        for cube in cubes:
            time_cube = cube.extract(iris.Constraint(time=time_point, latitude=lat_range))
            time_cube_list.append(time_cube)
        combined_cubes = iris.cube.CubeList(time_cube_list).merge_cube()
        combined_cubes.attributes['STASH'] = ','.join(stashcodes)
        combined_cubes.attributes['time'] = str(combined_cubes.coord('time').cell(time_point).point)
        saving_name = f"{file_path}{rose}/small_nc_files/Rgn_{time_point}.nc"
        iris.save(combined_cubes, saving_name, netcdf_format="NETCDF4")
