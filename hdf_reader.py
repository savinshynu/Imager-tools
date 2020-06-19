import numpy as np
import sys
import matplotlib.pyplot as plt
#import pyfits
import h5py
import time
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord,pix2hrz
from utils import ang_dist_pix
from scipy.stats import scoreatpercentile as sc
import matplotlib.animation as animation
from mad import median_absolute_deviation as mad

filename = sys.argv[1]

print filename
print 'getting data from hdf5 files' 
hf=h5py.File(filename,'r')
dset1=hf.get('image')
dset2=hf.get('time')
dset3 = hf.get('header')
hdr  = dset3[0,:]
intg =  dset1.shape[0]
xSize = hdr[1]
ySize = hdr[1]

lat1= 34.0689         # coordinates for LWA station 1
lon1= -107.6284

lat2 = 34.3484       #cordinates for LWA-SV sation
lon2 = -106.8858

stokes= 0

size = hdr[1]
psize = hdr[2]
cent_az = hdr[9]
cent_alt = hdr[10]


print dset1.shape
print dset2.shape


#print dset2[100:120,:]

im = dset1[100,0,:,:] 
#data_int = np.transpose(np.mean(im[:,:,0,420:421],axis=2))
#ax.set(xlabel='x pixel',ylabel='y pixel')
#plt.pcolormesh(np.transpose(im[:,:]),cmap='jet')
#plt.imshow(np.transpose(im[:,:]),cmap='jet',origin='lower')
#plt.show()
#print np.reshape(dset2[0,:],(4,)) 

stokes = 0 
noise_dat = np.zeros(dset1.shape[0]-6)
n = 0
for x in range(dset1.shape[0]):
    if x > 5:
        sub_im  = dset1[x,stokes,:,:] - np.mean(dset1[x-6:x,stokes,:,:],axis =0) 
        noise_dat[n] = np.std(sub_im[80:,80:]) 
        n += 1

print np.median(noise_dat),np.mean(noise_dat) 
plt.plot(range(noise_dat.shape[0]),noise_dat)
plt.show()


"""
peak = range(540,564)
def on_press(event):
    if event.key.isspace():
        if anim.running:
            anim.event_source.stop()
        else:
            anim.event_source.start()
        anim.running ^= True


#making movie 
fig = plt.figure()
ax1=fig.add_subplot(1,2,1)
ax2=fig.add_subplot(1,2,2)
listims = []
fig.canvas.mpl_connect('key_press_event', on_press)

for i,index in enumerate(peak) :

    im_i = np.transpose(dset1[index,0,:,:])
    imin = sc(im_i,0.5)
    imax = sc(im_i,99.5)
    im1 = ax1.imshow(im_i,cmap ='jet',vmin = imin, vmax = imax, origin ='lower')

    im_v  = np.transpose(dset1[index,3,:,:])
    vmin = sc(im_v,0.5)
    vmax = sc(im_v,99.5)
    im2 = ax2.imshow(im_v,cmap ='jet',vmin = vmin, vmax =vmax,origin ='lower')


    listims.append([im1, im2])

anim = animation.ArtistAnimation(fig, listims, interval=500, repeat_delay=1000)
anim.save('movie_ap.mp4')
anim.running = True
plt.show()



print im[54:58,70:73]
time = dset2[50,:]
ra = 350.8
dec = 58.8
xpix,ypix,altpix=trackcoord(ra, dec, time[0], time[1], time[2], time[3], size, psize, lat2, lon2, cent_az, cent_alt)
print xpix, ypix


az1, el1 =  pix2hrz(97, 85, size, psize, cent_az, cent_alt)

az2, el2 =  pix2hrz(98.8, 81.9, size, psize, cent_az, cent_alt)

print ang_dist_pix(az1,el1,az2,el2)
"""
