import numpy as np
import sys, glob
import h5py
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#from scipy.stats import ks_2samp as ks
from mad import median_absolute_deviation as mad
from scipy.optimize import curve_fit

filename = sys.argv[1]

hf = h5py.File(filename,'r')

im_dat =  hf.get('image')
t_dat = hf.get('time')
head_dat = hf.get('header')
el_dat = hf.get('elev')

#print t_dat[:10,:]
time = t_dat[:,1]*60.0 + t_dat[:,2] + t_dat[:,3]/60.0

#print time[10]

trange = np.arange(0,1440,5)

cyg, cas, vir, tau = [0,1,2,3]
stokes = 0
source = cyg
freq_in = 30
el_range = np.arange(20,85)
mean_el_dat = np.zeros(el_range.shape[0])

for n,el_low in enumerate(el_range):

    ind1  = np.where((el_dat[:,source] > el_low) & (el_dat[:,source] < (el_low + 1)) & (np.int_(head_dat[:,0])==freq_in))[0]
    #print len(ind1)
    
    if len(ind1) > 0:
       med_dat_source = np.nanmedian(im_dat[ind1,source,:,stokes],axis = 0)
       sig_dev_source = np.zeros(ind1.shape)

       for i,x in enumerate(ind1):
           dev_source = im_dat[x,source,:,stokes] / med_dat_source
           sig_dev_source[i] = np.nanstd(dev_source)

       med_sig_dev_source =  np.nanmean(sig_dev_source)

       #if med_sig_dev_source < 0.05:
       #   print "good time at elavation %d" % (el_low)

       mean_el_dat[n] = med_sig_dev_source


plt.scatter(el_range,mean_el_dat) 
plt.xlabel('Elevation in degrees')
plt.ylabel('Scintillation Index (a.u.)')
plt.show()


mean_time_dat = np.zeros(trange.shape[0])
time_dat = np.zeros(trange.shape[0])
nz = 0

for m,t in enumerate(trange):

    ind2  = np.where((time > t) & (time < (t + 5)) & (np.int_(head_dat[:,0])==30) & (el_dat[:,source] > 20))[0]
    #print len(ind2)

    if len(ind2) > 0:
       med_dat_source = np.nanmedian(im_dat[ind2,source,:,stokes],axis = 0)
       sig_dev_source = np.zeros(ind2.shape)

       for k,l in enumerate(ind2):
           dev_source = im_dat[l,source,:,stokes] / med_dat_source
           sig_dev_source[k] = np.nanstd(dev_source)

       med_time =  np.mean(sig_dev_source)

       #if med_sig_dev_source < 0.05:
       #   print "good time at elavation %d" % (el_low)

       mean_time_dat[nz] = med_time
       time_dat[nz] = t+2.5
       nz += 1


plt.scatter(time_dat[:nz]/60.0,mean_time_dat[:nz])
plt.xlabel('Time in UTC')
plt.ylabel('Scintillation Index (a.u.)')
plt.show()
