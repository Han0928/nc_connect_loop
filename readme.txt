int_save_nc_34101-ocean.py this is the final version of extract stashcode from rose-id(such as u-ct706)into separate .nc files,
next step is connect them together: concatenate_101.py, running this you can get a full  concatenate cubes of data and then save the concatenated data as a new NetCDF file, the script deletes the smaller NetCDF files that were used for the concatenation

Now, for the global model id(u-cy282),I also want to do the same thing, I name it as: int_save_nc_cy282.py

besides: 004_int_save_nc_v6/5/4/3/2.py looks like failed version of loop stashcode using a array,but I didn't delete them, so you can just ignore them(004***.py)
