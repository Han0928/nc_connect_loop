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
files_directory_UKCA='/jet/home/ding0928/cylc-run/'+rose+'/share/cycle/'
# Define latitudinal range of interest
latbottom = 35.0
lattop = 45.0

# Define the days and filechunks of interest
days = ['0720', '0722', '0724', '0726', '0728', '0730', '0801', '0803', '0805', '0807', '0809']
filechunks = ['pb','pc','pd','pe']
run = '20140720T0000Z'

rosefolder = '/jet/home/ding0928/python_analysis/Han_connect/nc_flie/'+rose+'/'
ncfolder = rosefolder+'small_nc_files/'


# Define the function to check if a cell is within the latitudinal range of interest
def lat_range(cell):
    return latbottom < cell < lattop

# Define the function to make directories if they don't already exist
def make_directories(nameofdir):
    newdir = os.path.join(file_path, rose, nameofdir)
    try:
        os.makedirs(newdir, exist_ok=True)
        print('Created Folder: {}'.format(newdir))
    except OSError as e:
        print(f"Error creating folder {newdir}: {e}")

for iday in days:
    # Define the file directory for each day
    files_directory = f'/jet/home/ding0928/cylc-run/{rose}/share/cycle/2014{iday}T0000Z/Regn1/resn_1/RA2M/um/'

    # Load the cubes for each stashcode
    cubes = []
    for chunk in filechunks:
        for stashcode in stashcodes:
            file_name = f"{files_directory}umnsaa_{chunk}006_Rgn_{stashcode}.nc"
            try:
                cubes.append(iris.load_cube(file_name))
            except:
                print(f"Unable to load {file_name}")

    # Process the cubes and save them in a single .nc file for each time slot
    for cube in cubes:
        cube = cube.extract(iris.Constraint(coord_values={'latitude': lat_range}))
        time = cube.coord('time').points[0]
        saving_name = f"{file_path}{rose}/small_nc_files/Rgn_{time}.nc"
        iris.save(cube, saving_name, netcdf_format="NETCDF4") 

    # Make the necessary directories
    make_directories("small_nc_files")
