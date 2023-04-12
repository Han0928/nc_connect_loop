import os
import datetime

# loop through all files in the current directory
for file_name in os.listdir('.'):
    # only process files with extension '.nc'
    if file_name.endswith('.nc'):
        # get the file's modification time in seconds since the epoch
        mtime = os.path.getmtime(file_name)
        # convert the epoch time to human-readable format
        mtime_str = datetime.datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')
        # print the file name and its modified time in human-readable format
        print(file_name, mtime_str)

