"""
Plot the light curve of an event

"""
import sys
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
from windows import tukey
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
from utils import avgp, maxi_unsub

def get_light_curve(hf, ra, dec, stokes):

    dset1 = hf.get('image')
    dset2 = hf.get('time')
    dset3 = hf.get('header')
    hdr  = dset3[0,:]
    intg =  dset1.shape[0]
  
    lat1= 34.0689         # coordinates for LWA station 1
    lon1= -107.6284

    lat2 = 34.3484       #cordinates for LWA-SV sation
    lon2 = -106.8858

    size = hdr[1]
    psize = hdr[2]
    cent_az = hdr[9]
    cent_alt = hdr[10]


    light=np.zeros((intg),dtype=float)

    for j in range(intg):

        data = dset1[j,stokes,:,:]
        time = dset2[j,:]
        xpix,ypix,altpix=trackcoord(ra, dec, time[0], time[1], time[2], time[3], size, psize, lat2, lon2, cent_az, cent_alt)
        xa, ya = maxi_unsub(data,xpix,ypix)
        light[j] =  avgp(xa,ya,data)
    
        #print light[k,0], light[k,1],time[1],time[2],time[3],xpix,ypix,xa,ya,altpix

    
    #Applying tukey window filtering
    y1 = np.fft.fft(light)
    y1 = y1*(tukey(y1.shape[0],alpha = 0.015)**2)
    data_filt = np.fft.ifft(y1)

    return data_filt
    
    
    





