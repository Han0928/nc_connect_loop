'''
I figured out why 1970 is there, it's jsut a reference time,really doesn't matter. 
this script is for analyzing connected_nc file, trying to make 1.time_series 2. spatial gif map.

--Han 02.28
'''
import numpy as np
import iris,sys,glob
import pandas as pd
from netCDF4 import Dataset
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy.interpolate import interp1d, RegularGridInterpolator
import matplotlib.pyplot as plt
import numpy.ma as ma
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import iris.plot as iplt
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
import imageio
import iris.plot as iplt
import datetime
import os, sys

import numpy as np
import netCDF4 as nc

rose = 'u-cs093'
stashcode='m01s34i101' 
rosefolder = '/jet/home/ding0928/python_analysis/Han_connect/nc_flie/'+rose+'/'
bigncfolder = rosefolder+'full_nc_files/'
smallncfolder = rosefolder+'small_nc_files/'

# Open the netCDF file and read the "time" variable
nc_file = nc.Dataset(bigncfolder+ "Rgn_mass_fraction_of_ammonia_in_air_m01s34i076.nc")
time_var = nc_file.variables["time"]

# Extract the time values from the "time" variable and convert them to a datetime array
time_array = nc.num2date(time_var[:], time_var.units)

# Print the resulting datetime array
print(time_array)


# so the first functino should be load the nc file and convert the 1970 time into 2014 time.
def ion_topo(cube, level, ylabel, map_min, map_max, saving_name, i):
    fig = plt.figure(figsize=(16, 12))
    ax = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
    data = cube.extract(iris.Constraint(model_level_number=level)).data
    # data = cube.extract(iris.Constraint(time=i)).data
    print('data',data)
    pl = ax.pcolormesh(data[i], vmin=map_min, vmax=map_max)
    pl.cmap.set_under('k')
    ax.add_feature(cfeature.LAND, facecolor='0.9')
    ax.add_feature(cfeature.OCEAN, facecolor='w')
    ax.add_feature(cfeature.COASTLINE, linewidth=2)
    ax.add_feature(cfeature.BORDERS, linestyle='--')
    # ax.set_xlim(lonmin,lonmax)
    # ax.set_ylim(latmin,latmax)
    ax.text(0.3, 1.10, 'SO$_{2}$ (pptv)', transform=ax.transAxes, fontsize=22, fontweight='bold', va='top')
    ax.text(0.1, 1.0, 'SO$_{2}$ (pptv)', transform=ax.transAxes, fontsize=22, fontweight='bold', va='top')
    plt.text(-104.99, 39.48, '*Denver', fontsize=22, color='red')
    cmap = plt.get_cmap('jet')
    cmap.set_under('w')
    cmap.set_over('k')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlabel_style = {'size': 22}
    gl.ylabel_style = {'size': 22}
    cb = plt.colorbar(pl, label=ylabel, extend="max", shrink=0.7)
    cb.set_label(label=ylabel, fontsize=22)
    cb.ax.tick_params(labelsize=22)
     # add timestamp on top left of plot
    timestamp = cube.coord('time').points[i]
    timestamp_str = datetime.datetime.utcfromtimestamp(timestamp).strftime('%Y-%m-%d %H:%M:%S')
    plt.text(0.05, 0.95, timestamp_str, transform=ax.transAxes, fontsize=18, fontweight='bold', va='top') 
    plt.savefig('/jet/home/ding0928/python_analysis/nc_long_series/' + saving_name + '_timestep_' + str(i) + '.png', dpi=100)

path='/jet/home/ding0928/python_analysis/Han_connect/nc_flie/u-cs093/full_nc_files/Rgn_mass_fraction_of_sulfur_dioxide_in_air_m01s34i072.nc'
so2_cube = iris.load_cube(path,iris.AttributeConstraint(STASH='m01s34i072'))
print('so2_cube',so2_cube)
so2_conc = so2_cube*0.029/0.064*1e12 #ppt
print(so2_conc.coord('time')) 

for j in [2]: #vertical level
    for i in range(2): #time slice
        ion_topo(so2_conc, j, 'SO$_{2}$ (pptv)', 0, 8, 'so2_level' + str(j), i)

plt.show()



# images = []
# for j in [2, 10]:
#     for i in range(12):
#         filename = '/jet/home/ding0928/python_analysis/ion_analysis/ion_lev' + str(j) + '_timestep_' + str(i) + '.png'
#         images.append(imageio.imread(filename))

#     imageio.mimsave('/jet/home/ding0928/python_analysis/ion_analysis/ion_lev'+ str(j)+'.gif', images, duration=0.5)

