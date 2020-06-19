#!/usr/bin/env python

import os
from lsl import astro
day = range(58950,58974)   #58887) #x = 58876
"""
k= astro.get_julian_from_sys()
mjd=k-2400000.5
x=int(mjd)-1
#print x
"""
for x in day:
    path_dir = '/leo/savin/wide_lwasv/oims_wide/images/'+str(x) 
    if os.path.isdir(path_dir):
       os.chdir("/leo/savin/wide_lwasv/oims_wide/standard_sub/")
       os.system("python  oims2hdf5.py "+str(x))
       os.system("python  wide_imager_lwasv.py "+str(x))
       os.system("python filter_list.py "+str(x)) 
       os.system("python response_source.py "+str(x))


       os.system("python  wide_imager_lwasv_stokesv.py "+str(x))
       os.system("python filter_list_stokesv.py "+str(x)) 
    
    else:
        print "No .oims file for the day"
