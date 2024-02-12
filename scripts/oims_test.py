
"""
Scripts to test the calibration of Orville imager
"""
import os
import glob
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt
from sliding_rfi_flagger_man import main
from lsl import astro
from scipy.optimize import curve_fit
from mad import median_absolute_deviation as mad
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
from OrvilleImageDB import OrvilleImageDB
from utils import mask_image,maxi_sub,maxi_unsub,read_oims_cont,avgp,avgp_chan,avgp_chan_eq
from collect_resp_data import collect_resp
import matplotlib.animation as animation
from twod_gaussian import fitgaussian

#for calibrations
"""
low = float(sys.argv[2])
high = float(sys.argv[3])

mjd = sys.argv[1]

im_dat,t_dat,head_dat,el_dat = collect_resp(mjd)
stokes_cal = 0

cyg, cas, vir, tau = [0,1,2,3]
source  = cyg
#ind =  np.where((el_dat[:,cyg]< high) & (el_dat[:,cyg] > low) &(head_dat[:,0] == 30.1375) )[0]
ind =  np.where((el_dat[:,source]< high) & (el_dat[:,source] > low) )[0]

if len(ind)==0:
   sys.exit("could not find required indices")



for x in xrange(ind.shape[0]):
    print ind[x], t_dat[ind[x],1], head_dat[ind[x],1], head_dat[ind[x],0], el_dat[ind[x],source]   

ydat = np.mean(im_dat[ind,source,:,stokes_cal],axis = 1)
xdat = el_dat[ind,source]

plt.scatter(xdat,ydat)
plt.show()
plt.plot(range(im_dat.shape[2]),np.mean(im_dat[ind,source,:,stokes_cal],axis=0))
plt.show()


indm = range(3828,3860)

xdat_ind = ind
spect_cal = np.mean(im_dat[xdat_ind,source,:,stokes_cal],axis =0)


def flux_jy(v):
    
    a,b,c = [4.695,0.085,-0.178]
    
    log_s = a + b*np.log10(v) + c*(np.log10(v)**2)
   
    return 10**log_s 

"""
   

filename2 = "/leo/savin/wide_lwasv/oims_wide/images/58850/58850_220000_30.137MHz_49.837MHz.oims"
nchan = 198
ngrid = 128
peak = [312]
#noise = range(94,100) #+range(251,256)

dset = np.zeros((len(peak),nchan,4,ngrid,ngrid))
#dset_avg = np.zeros((len(noise),nchan,4,ngrid,ngrid))
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
#flux_cal_jy = flux_jy(f_chan)
#print f_chan.size     
"""
for j,index in enumerate(noise):
    hdr, dat = db.__getitem__(index)
    dset_avg[j,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
"""
dset = np.mean(dset, axis = 0)
#dset_avg = np.mean(dset_avg, axis = 0)
sub_dat = dset #- dset_avg
stokes = 0


#plt.plot(range(nchan),sub_dat[:,stokes,86,74])
#plt.show()

im = np.mean(sub_dat[:,stokes,:,:],axis =0)
im_chan = sub_dat[80,stokes,:,:]
#plt.pcolormesh(np.transpose(im),cmap ='jet')
#plt.pcolormesh(np.transpose(im_chan),cmap ='jet')
#plt.show()

#params = fitgaussian(im_chan[54-14:54+14,88-14:88+14])

#print params
ra_cyg, dec_cyg = [299.8, 40.7]
ra_cas, dec_cas = [350.8, 58.8] 

annul_info = [ngrid, psize, MJD, H, M, S, lat, lon, cent_az, cent_alt]
resp1 = avgp_chan_eq(ra_cas,dec_cas,f_chan,annul_info,sub_dat[:,stokes,:,:])
resp2 = avgp_chan_eq(ra_cyg,dec_cyg,f_chan,annul_info,sub_dat[:,stokes,:,:])

#resp1 = avgp_chan(54,88,sub_dat[:,stokes,:,:],f_chan)
#resp2 = avgp_chan(86,74,sub_dat[:,stokes,:,:],f_chan)

#resp1 = (resp1/resp1.max())*100
#resp2 = (resp2/resp2.max())*100


ratio1 = resp1/resp2
#ratio2 = sub_dat[:,stokes,54,88]/sub_dat[:,stokes,86,74]

#plt.plot(range(len(resp1)),resp1/resp1.max())
#plt.plot(range(len(resp2)),resp2/resp2.max())
plt.plot(range(len(ratio1)),ratio1,label='Annulus subtracted and integrated over point source')
#plt.plot(range(len(ratio2)),ratio2,label="peak value from point source")
plt.xlabel('channels')
plt.ylabel('Ratio of Cas with Cyg')
#plt.ylim(0,150)
plt.legend()
plt.show()

"""
#making movie 
fig = plt.figure()
ims = []
for i in range(nchan):
    im = plt.pcolormesh(np.transpose(sub_dat[i,stokes,:,:]),cmap ='jet')
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=40, blit = True, repeat_delay=1000)
ani.save('movie_test.mp4')
plt.show()




im_chan = sub_dat[84,stokes,:,:]
print np.std(im[80:,:40]), np.std(im_chan[80:,:40])

flag = [144,145,105]
#flag = []
chan_range = list(set(range(nchan)) -set(flag))
nchan1 = len(chan_range)
nchan = nchan1
print nchan
sub_dat_mod = np.zeros((nchan,ngrid,ngrid))
for i,chan in enumerate(chan_range): #chan in range(nchan):
    #sub_dat_mod[chan,:,:] = (sub_dat[chan,stokes,:,:]/spect_cal[chan])*flux_cal_jy[chan]
    sub_dat_mod[i,:,:] = (sub_dat[chan,stokes,:,:]/spect_cal[chan])*flux_cal_jy[chan]
    

chan_rms = np.zeros(nchan)

for chan in range(nchan):
    chan_rms[chan] = np.std(sub_dat_mod[chan,80:,80:])

plt.plot(range(nchan),chan_rms)
plt.xlabel("Channels")
plt.ylabel("RMS in Jy")
plt.show()


int_chan_rms = np.zeros(nchan)

for chan in range(nchan):
    
    intdat = np.std(np.mean(sub_dat_mod[:chan+1,80:,80:],axis = 0))
    int_chan_rms[chan] = intdat 

plt.plot(range(nchan),int_chan_rms)
plt.xlabel("Number of averaged channels together")
plt.ylabel("RMS in Jy")
plt.show()

plt.loglog(range(nchan),int_chan_rms)
plt.xlabel("Number of averaged channels together in log scale")
plt.ylabel("RMS in Jy log scale")
plt.show()
"""
