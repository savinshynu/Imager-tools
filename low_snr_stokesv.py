import numpy as np
import sys
import matplotlib.pyplot as plt
import glob
import h5py
import time as t
from windows import tukey
import scipy.stats as st
from mad import median_absolute_deviation as mad
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
from utils import maxi_sub,maxi_unsub,avgp,ang_dist

vlss=np.loadtxt("/leo/savin/wide_lwasv/oims_wide/standard_sub/VLSS_10.txt")
def search_vlss(ra,dec):
    ra_in = np.ones(vlss.shape[0])*ra
    dec_in = np.ones(vlss.shape[0])*dec
    diff = ang_dist(ra_in, dec_in, vlss[:,0], vlss[:,1])
    par = np.where((vlss[:,2] > 10.0) & (diff < 3.0))[0]
    return len(par)



def filter_snr(inf,mjd_day):

    num1 = 0
    stokes = [0,1,2,3]
    list1 = np.zeros(inf.shape)
    for hour in range(0,24):
        #find the corresponding .hdf5 file for the time 
        path1 = '/leo/savin/wide_lwasv/oims_wide/standard_sub/images_hdf5/'+str(mjd_day)+'/'+str(mjd_day)+'_'+str(format(hour,'02d'))+'*.hdf5'
        files1  = sorted(glob.glob(path1))
        hour_index = np.where((inf[:,3] == hour))[0]
        if len(hour_index) == 0:
           print "no transients at utc" +str(hour)
           continue

        for filename1 in files1:
            print filename1
            hf = h5py.File(filename1,'r')
            img = hf.get('image')
            tm = hf.get('time')
            header = hf.get('header')
            head =  header[0,:]
            ints = head[0] #number of integrations
            ngrid  = int(head[1]) #size of images (x-axis)
            psize = head[2] #angular size of a pixel (at zenith)
            nchan = int(head[3])
            cent_az = head[9]
            cent_alt = head[10]
            t_comp = tm[:,1]*3600.0 + tm[:,2]*60.0 + tm[:,3]
            

            for i in hour_index:
                ra = float(inf[i,0])
                dec = float(inf[i,1])
                MJD = float(inf[i,2])
                h = float(inf[i,3])
                m = float(inf[i,4])
                s = float(inf[i,5])
       
                mint = int(m)
                sec = s+(m-mint)*60.0
                if sec >= 60.0:
                   sec = int(sec%60.0)
                   mint += 1
                   if mint >= 60:
                      mint = mint%60
                      h += 1 
                sec = int(sec)
                h  = int(h) 
                t_in = h*3600.0 + mint*60.0 + sec      
                peak = np.where((abs(t_comp-t_in) < 3))[0]
                if len(peak)==0:
                   #print"requested time not in dataset"
                   print inf[i,:]
                   continue 
                if len(peak) > 1:
                   peak = int(np.median(peak))
                    
                #signal = np.arange(max(0,peak-305),min(peak+305,ints))
                signal = np.arange(ints)
                noise = np.where((signal < peak-5)| (signal > peak+5))[0]
                event = np.where((signal > peak-3)&(signal < peak+3))[0]
    
                intg = len(signal)
                if intg == 0:
                   continue
                light=np.zeros((intg,4),dtype=float)
                data_filt = np.zeros((intg,4),dtype=np.complex64)
       
                lat1= 34.0689         # coordinates for LWA station 1
                lon1= -107.6284

                lat2 = 34.3484       #cordinates for LWA-SV sation
                lon2 = -106.8858

                for ln,j in enumerate(signal):
                    time = tm[j,:]
                    xpix,ypix,altpix=trackcoord(ra,dec,time[0],time[1],time[2],time[3],ngrid,psize,lat2,lon2,cent_az,cent_alt)
                    data_i = img[j,stokes[0],:,:]
                    xa, ya = maxi_unsub(data_i,xpix,ypix)

                    for pol in stokes:
                        data = img[j,pol,:,:]
                        light[ln,pol] =  avgp(xa,ya,data) 
            
                for pol in stokes:
                    y1 = np.fft.fft(light[:,pol])
                    y1 = y1*(tukey(y1.shape[0],alpha = 0.015)**2)
                    data_filt[:,pol] = np.fft.ifft(y1)
              
                stokes_par = 3 #stokes V
                 
                sig_noise_mn = np.std(data_filt[:,stokes_par].real[noise])
                sig_noise_md = mad(data_filt[:,stokes_par].real[noise])
                
                """
                if sig_noise_md == 0:
                   sig_noise = sig_noise_mn
                
                else:
                   sig_noise = sig_noise_md
                """
                sig_noise = sig_noise_mn #changes

                med = np.median(data_filt[:,stokes_par].real[noise])
                mean = np.mean(data_filt[:,stokes_par].real[noise])

                if sig_noise == sig_noise_md:
                   med_par = med
                else:
                   med_par = mean
                   
                snr = (data_filt[:,stokes_par].real[event].max()-med_par)/sig_noise               
                
                """
                if (snr > 5.0):
                   list1[num1,:] = inf[i,:]
                   num1 += 1
                """
                dat_noise = data_filt[:,stokes_par].real[noise]
                inp_dat_noise = dat_noise-med_par
                top = np.where((inp_dat_noise > 5*sig_noise))[0]

                if (snr > 5.0):
                   if len(top) < 2:
                      list1[num1,:] = inf[i,:]
                      num1 += 1
                   else:
                       if (search_vlss(inf[i,0],inf[i,1]) == 0) :
                          list1[num1,:] = inf[i,:]
                          num1 += 1

                  

    
            hf.close() 

    return list1[:num1,:]    


