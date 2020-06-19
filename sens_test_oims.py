import numpy as np
import matplotlib.pyplot as plt
#from sliding_rfi_flagger_man import main
#from lsl import astro
import glob, os
import sys
import h5py
from scipy.optimize import curve_fit
#from mad import median_absolute_deviation as mad
#from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
from OrvilleImageDB import OrvilleImageDB
#from utils import mask_image,maxi_sub,maxi_unsub,read_oims_cont,avgp,avgp_chan,avgp_chan_eq
#import matplotlib.animation as animation
#from twod_gaussian import fitgaussian
from scipy.optimize import curve_fit

filename = sys.argv[1]

hf = h5py.File(filename,'r')

im_dat =  hf.get('image')
t_dat = hf.get('time')
head_dat = hf.get('header')
el_dat = hf.get('elev')

#print t_dat[:10,:]


cyg, cas, vir, tau = [0,1,2,3]
stokes = 0
source = cyg
ind  = np.where((el_dat[:,source] > 10)&(np.int_(head_dat[:,0])==30))[0]
#print ind
ydat = np.mean(im_dat[ind,source,:,stokes],axis = 1)
xdat = el_dat[ind,source]

#print head_dat[0,0]
def func(x,a,b,c):

    return a*x**3+b*x**2+c*x


popt, pcov = curve_fit(func, xdat, ydat)

perr = np.sqrt(np.diag(pcov))

print popt
print perr

a,b,c = [ -2.60298718e-04,   2.28432424e-02,   2.15106031e+00]

range_el = np.linspace(0,90,100)

plt.scatter(xdat,ydat,c ='r')
#plt.plot(range_el,func(range_el,*popt),c ='b')
plt.plot(range_el,func(range_el,a,b,c),c ='b')
plt.xlabel("Elevation in degrees")
plt.ylabel("Power in a.u.")
plt.show()

ind2  = np.where((el_dat[:,source] > 83) & (el_dat[:,source] < 85)&(np.int_(head_dat[:,0])==30))[0]
spect_cal = np.mean(im_dat[ind2,source,:,stokes],axis=0)

plt.plot(range(len(spect_cal)), spect_cal)
plt.show()




def flux_jy(v):
    
    a,b,c = [4.695,0.085,-0.178]
    
    log_s = a + b*np.log10(v) + c*(np.log10(v)**2)
   
    return 10**log_s 

sig_noise_i = 0.0145366604753
sig_noise_v = 0.00733938874187
fr = 39.987

el = np.arange(20,91)
sens_i = (sig_noise_i*flux_jy(fr))/func(el,a,b,c)
sens_v = (sig_noise_v*flux_jy(fr))/func(el,a,b,c)

plt.subplot(1,2,1)
plt.plot(el,sens_i,c ='b', label = 'Stokes I')
plt.plot(el,sens_v,c ='r', label = 'Stokes V')
plt.xlabel('Elevation Angle (degrees)')
plt.ylabel('RMS Noise (Jy)')
plt.ylim(0,8)
plt.legend()
#plt.show()

filename2 = "/leo/savin/wide_lwasv/oims_wide/images/58939/58939_090000_30.137MHz_49.837MHz.oims"
#filename2 = "/leo/savin/wide_lwasv/oims_wide/images/58939/58939_090000_30.100MHz_49.800MHz.oims"
nchan = 198
ngrid = 128
peak = [402]
noise = range(396,402) #+range(202,208)

dset = np.zeros((len(peak),nchan,4,ngrid,ngrid))
dset_avg = np.zeros((len(noise),nchan,4,ngrid,ngrid))

db = OrvilleImageDB(filename2,'r')

print db.nint
for j,index in enumerate(peak):
    hdr, dat = db.__getitem__(index)
    dset[j,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))


ngrid = int( db.header.ngrid) #size of images (x-axis)
psize = db.header.pixel_size #angular size of a pixel (at zenith)
nchan = int(db.header.nchan) # number of channels
f1 = hdr['start_freq']/1e+6
bw = hdr['bandwidth']/1e+6
f2 = hdr['stop_freq']/1e+6 + bw
t = hdr['start_time']
cent_az = hdr['center_az']
cent_alt = hdr['center_alt']
MJD = int(t)
h = 24.0*(t - float(MJD))
H = int(h)
m = 60*(h - float(H))
M = int(m)
s = 60*(m - float(M))
S = int(s)
lat = 34.3484       #cordinates for LWA-SV sation
lon = -106.8858

f_chan = np.arange(f1,f2,bw)
flux_cal_jy = flux_jy(f_chan)
#print f_chan.size     

for j,index in enumerate(noise):
    hdr, dat = db.__getitem__(index)
    dset_avg[j,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))

dset = np.mean(dset, axis = 0)
dset_avg = np.mean(dset_avg, axis = 0)
sub_dat = dset - dset_avg
#stokes = 0
for num,stokes in enumerate([0,3]):

    #im = np.mean(sub_dat[:,stokes,:,:],axis =0)
    
    #im_chan = sub_dat[80,stokes,:,:]
    #plt.pcolormesh(np.transpose(im),cmap ='jet')
    #plt.pcolormesh(np.transpose(im_chan),cmap ='jet')
    #plt.show()



    flag = [143,144,145]
    #flag = []
    chan_range = list(set(range(nchan)) -set(flag))
    nchan1 = len(chan_range)
    nchan = nchan1
    #print nchan
    im = np.mean(sub_dat[chan_range,stokes,:,:],axis =0)
    print np.std(im[80:,80:])



    sub_dat_mod = np.zeros((nchan,ngrid,ngrid))
    for i,chan in enumerate(chan_range): #chan in range(nchan):
        #sub_dat_mod[chan,:,:] = (sub_dat[chan,stokes,:,:]/spect_cal[chan])*flux_cal_jy[chan]
        sub_dat_mod[i,:,:] = (sub_dat[chan,stokes,:,:]/spect_cal[chan])*flux_cal_jy[chan]
    

    chan_rms = np.zeros(nchan)

    for chan in range(nchan):
        chan_rms[chan] = np.std(sub_dat_mod[chan,80:,80:])

    #plt.plot(range(nchan),chan_rms)
    #plt.xlabel("Channels")
    #plt.ylabel("RMS in Jy")
    #plt.show()


    int_chan_rms = np.zeros(nchan)

    for chan in range(nchan):
    
        intdat = np.std(np.mean(sub_dat_mod[:chan+1,80:,80:],axis = 0))
        int_chan_rms[chan] = intdat 
    
    var = ['I', 'V']
    np.savetxt('noise_'+var[num]+'.txt', int_chan_rms)
    #plt.subplot(1,2,num+1)
    #plt.plot(range(nchan),int_chan_rms, label = 'Stokes ' + var[num])
    #plt.xlabel("Number of averaged channels together")
    #plt.ylabel("RMS noise in Jy")
    #plt.legend()
    #plt.ylim(0,15)

    #plt.loglog(range(nchan),int_chan_rms)
    #plt.xlabel("Number of averaged channels together in log scale")
    #plt.ylabel("RMS in Jy log scale")
    #plt.show()

ns_i = np.loadtxt('noise_I.txt')
ns_v = np.loadtxt('noise_V.txt')
plt.subplot(1,2,2)
plt.plot(range(len(ns_i)),ns_i, c ='b', label = 'Stokes I')
plt.plot(range(len(ns_v)),ns_v, c= 'r', label = 'Stokes V')
plt.xlabel("Number of averaged channels together")
plt.ylabel("RMS Noise (Jy)")
plt.legend()
plt.ylim(0,15)
plt.show()


