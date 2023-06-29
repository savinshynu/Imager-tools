"""
Go through the oims images in a day, collect the spectral responses of 
A class sources for calibration purposes and write them into a hdf5 file

"""

import glob
import os
import numpy as np
import h5py
import time
import sys
from OrvilleImageDB import OrvilleImageDB
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
from OrvilleImageDB import OrvilleImageDB
from utils import mask_image,maxi_sub,maxi_unsub,read_oims_cont

mjd = sys.argv[1]

path3 = '/leo/savin/wide_lwasv/oims_wide/standard_sub/resp_source/'+str(mjd)+'/'
try: 
    os.mkdir(path3)
except OSError:
    print 'directry already exists'

dir = '/leo/savin/wide_lwasv/oims_wide/images/'
   
path1 = dir + str(mjd)+'/*.oims'
files1 = sorted(glob.glob(path1))

#Listing different statiion parameters

lat = 34.3491 #Station Latitude (of LWA-SV)
lon = -106.886  #Station Latitude (of LWA-SV)

# Listing the coordinates for different sources

ra_cas = 350.866417
dec_cas = 58.811778

ra_cyg = 299.868153
dec_cyg = 40.733916

ra_vir = 187.705930
dec_vir = 12.391123

ra_tau = 83.633212
dec_tau = 22.014460


# Mention the source for tracking and change based on above info for a different one
ra_cal = [ra_cyg,ra_cas,ra_vir,ra_tau]
dec_cal = [dec_cyg,dec_cas,dec_vir,dec_tau]



h5_name = str(mjd)+'cal'

max_int = 17300
nchan = 198
outfile = h5py.File(path3+h5_name+'.hdf5','a')
dset1 = outfile.create_dataset('image',(max_int,4,nchan,4),dtype='float32') #dataset with the response along 108 channels 4 and stokes
dset2 = outfile.create_dataset('time',(max_int,4),dtype='float32') # dataset with time info for each integration
dset3 = outfile.create_dataset('header',(max_int,9),dtype='float32') #dataset with header info like start,end frequency and bandwidth
dset4 = outfile.create_dataset('elev',(max_int,4),dtype='float32') #dataset with azimuth and elevation of source   
stokes = 0

capt = 0

for filename in files1:
    
    print filename
    db = OrvilleImageDB(filename,'r') #define the image data base

    #getting parameters from the input file

    ints = db.nint #number of integration
    #station =  db.header.station
    #stokes = db.header.stokes_params
    #inp_flag = db.header.flags
    file_start = db.header.start_time
    file_end = db.header.stop_time

    ngrid = db.header.ngrid #size of images (x-axis)
    psize = db.header.pixel_size #angular size of a pixel (at zenith)
    nchan = db.header.nchan # number of channels
   
    for intg in range(ints):
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
        
        
        data = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
   
        dat_im = np.mean(data[:,stokes,:,:],axis=0) # Averaged image over all channels to find the location of source
        
        for i in range(4):

            ra = ra_cal[i]
            dec = dec_cal[i]
            xp,yp,alt = trackcoord(ra, dec, MJD, H, M, S, ngrid, psize, lat, lon, cent_az, cent_alt) #tracking source location

            xa,ya = maxi_unsub(dat_im,xp,yp) # finding the peak pixels

            # Copying information to different dataset
            dset1[capt,i,:,:] = data[:,:,xa,ya] #spectral response in  4 stokes
            dset4[capt,i] = alt  # Elevation information
          

        dset2[capt,:] = [MJD,H, M, S]
        dset3[capt,:] = [start_freq, stop_freq, bandwidth, lst, cent_az, cent_alt, ngrid, psize, nchan]
        capt += 1

    db.close()    
      
        
              
# closing the write out file
outfile.close()

