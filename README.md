The purpose of this repository is to 
1. extract stashcode into .nc files and save them in a small_nc folder. 
2. these small files are then connected into a big full file for the entire simulation time series, and the small ones are deleted.

Working Code
In this folder, you will find the final version of the working code:
int_save_nc_34101-ocean.py - for purpose 1
concatenate_101.py - for purpose 2

Other Files
The following files are iterations of the for loop and a jupyter version, but they are not currently functioning correctly and require further debugging:
004_int_save_nc_vXXX.py
small_nc_XXX.ipynb
These files were created because the loop takes too much time over pb/c/d/e and 5 stashcodes, 
so the code was transferred to a Jupyter notebook. However, further work is required to make these files functional.
