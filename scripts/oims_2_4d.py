"""
Converting an oims image file to hdf5 file format
"""

import glob
import os
import numpy as np
import h5py
import time
import sys
from OrvilleImageDB import OrvilleImageDB

mjd = sys.argv[1]

path3 = '/leo/savin/wide_lwasv/oims_wide/standard_sub/images_hdf5_dif/'+str(mjd)+'/'
try: 
    os.mkdir(path3)
except OSError:
               print 'directry already exists'

dir = '/leo/savin/wide_lwasv/oims_wide/images/'
   
path1 = dir + str(mjd)+'/*210000_30.137MHz_49.837MHz.oims'
files1 = sorted(glob.glob(path1))

for filename in files1:
    
    print filename
    db = OrvilleImageDB(filename,'r') #define the image data base

    #getting parameters from the input file

    ints = db.nint #number of integration
    station =  db.header.station
    stokes = db.header.stokes_params
    inp_flag = db.header.flags
    file_start = db.header.start_time
    file_end = db.header.stop_time

    ngrid = db.header.ngrid #size of images (x-axis)
    psize = db.header.pixel_size #angular size of a pixel (at zenith)
    nchan = db.header.nchan # number of channels
   
    file_ext = os.path.basename(filename)
    h5_ext = os.path.splitext(file_ext)[0]
  
    outfile = h5py.File(path3+h5_ext+'.hdf5','a') 
    dset1 = outfile.create_dataset('image',(10,nchan,4,ngrid,ngrid),dtype='float32')
    dset2 = outfile.create_dataset('time',(10,4),dtype='float32')
    dset3 = outfile.create_dataset('header',(10,11),dtype='float32')
   

    for i in range(10):
        hdr, dat = db.read_image()
        t = hdr['start_time']
        int_len = hdr['int_len']
        lst = hdr['lst']
        start_freq = hdr['start_freq']/1e+6
        stop_freq = hdr['stop_freq']/1e+6
        bandwidth = hdr['bandwidth']/1e+6
        cent_az = hdr['center_az']
        cent_alt = hdr['center_alt']

        MJD = int(t)
        h = 24.0*(t - float(MJD))
        H = int(h)
        m = 60*(h - float(H))
        M = int(m)
        s = 60*(m - float(M))
        S = int(s)
           
        # Copying the array from oims to main dataset
        dset1[i,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
        dset2[i,:] = np.array([MJD,H,M,S])
        dset3[i,:] = np.array([ints,ngrid,psize,nchan,int_len,lst,start_freq,stop_freq,bandwidth,cent_az,cent_alt])
      
        
    db.close() 
         
    # closing the write out file
    outfile.close()

