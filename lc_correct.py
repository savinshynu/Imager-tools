import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import glob
import h5py
#import time as t
from windows import tukey
from sliding_rfi_flagger_adv import main
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
from OrvilleImageDB import OrvilleImageDB
from utils import mjd2ut,avgp,maxi_unsub,maxi_sub
from mad import median_absolute_deviation as mad


lat1= 34.0689         # coordinates for LWA station 1
lon1= -107.6284

lat2 = 34.3484       #cordinates for LWA-SV sation
lon2 = -106.8858

def snr_in_lc(ra,dec,db,peak,stokes_par):
    
    #getting parameters from the input file
    ints = db.nint #number of integration
   
    ngrid = db.header.ngrid #size of images (x-axis)
    psize = db.header.pixel_size #angular size of a pixel (at zenith)
    nchan = db.header.nchan # number of channels
    
    peak = int(np.median(peak))
    signal = np.arange(max(0,peak-105),min(peak+105,ints))
    noise = np.where((signal < peak-5)| (signal > peak+5))[0]
    event = np.where((signal > peak-3)&(signal < peak+3))[0]
    
    
    light=np.zeros((len(signal),4),dtype=float)
    data_filt = np.zeros((len(signal),4),dtype=np.complex64)
    stokes =[0,1,2,3]

    for k,i in enumerate(signal):
        try:
            hdr, dat = db.__getitem__(i)
        except RuntimeError:
            print 'bad database at %d integration' % (i)
            continue

        t = hdr['start_time']
        #int_len = hdr['int_len']
        #lst = hdr['lst']
        #start_freq = hdr['start_freq']/1e+6
        #stop_freq = hdr['stop_freq']/1e+6
        #bandwidth = hdr['bandwidth']/1e+6
        cent_az = hdr['center_az']
        cent_alt = hdr['center_alt']
        MJD, H, M, S = mjd2ut(t) 
        data = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
       
        dat_im = np.mean(data[:,stokes[0],:,:],axis=0)
        xpix,ypix,altpix = trackcoord(ra,dec,MJD,H,M,S,ngrid,psize,lat2,lon2,cent_az,cent_alt)

        xa, ya = maxi_unsub(dat_im,xpix,ypix)

        for pol in stokes:
            gd = main(data[:,pol,xa,ya])    
            dat_im_flag = np.mean(data[gd,pol,:,:],axis=0)
            light[k,pol] =  avgp(xa,ya,dat_im_flag) 

    for pol in stokes:
        y1 = np.fft.fft(light[:,pol])
        y1 = y1*(tukey(y1.shape[0],alpha = 0.015)**2)
        data_filt[:,pol] = np.fft.ifft(y1)

    #stokes_par = 0 #for stokes I
    sig_noise_mn = np.std(data_filt[:,stokes_par].real[noise])
    sig_noise_md = mad(data_filt[:,stokes_par].real[noise])
    """
    if sig_noise_md == 0:
       sig_noise = sig_noise_mn
    else:
       sig_noise = sig_noise_md
    """
    sig_noise = sig_noise_mn #changes made for mean

    med_flag_i = np.mean(data_filt[:,0].real[noise]) #changed median to mean
    med_flag_q = np.mean(data_filt[:,1].real[noise])
    med_flag_u = np.mean(data_filt[:,2].real[noise])
    med_flag_v = np.mean(data_filt[:,3].real[noise])
    med_flag_st = np.mean(data_filt[:,stokes_par].real[noise])
   
    

    #med_st  = np.median(data_filt[:,stokes_par].real[noise])

    """
    if sig_noise == sig_noise_md:
       med_flag_st = med_st
    else:
       med_flag_st = mean_st
    """

    snr = (data_filt[:,stokes_par].real[event].max()-med_flag_st)/sig_noise

    lin_pol = (((data_filt[:,1].real[event].max()-med_flag_q)**2+(data_filt[:,2].real[event].max()-med_flag_u)**2)**0.5)/np.abs(data_filt[:,0].real[event].max()-med_flag_i)

    circ_pol = np.abs(data_filt[:,3].real[event].max()-med_flag_v)/np.abs(data_filt[:,0].real[event].max()-med_flag_i)




    #snr = data_filt[:,stokes_par].real[event].max()/sig_noise
   
    #lin_pol = (((data_filt[:,1].real[event].max())**2+(data_filt[:,2].real[event].max())**2)**0.5)/np.abs(data_filt[:,0].real[event].max())

    #circ_pol = np.abs(data_filt[:,3].real[event].max())/np.abs(data_filt[:,0].real[event].max())
 
    return snr,lin_pol,circ_pol





