"""
Create all-sky image movie of the candidate

"""


import numpy as np
import sys, os
import matplotlib.pyplot as plt
import glob
import h5py
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
from OrvilleImageDB import OrvilleImageDB
from scipy.stats import scoreatpercentile as sc
import matplotlib.animation as animation

ra, dec, mjd, h, m, s, stokes, filename = sys.argv[1:]

ra = float(ra)
dec = float(dec)
mjd = int(float(mjd))
h = int(float(h))
m = int(float(m))
s = int(float(s))
stokes = int(stokes)
#filename = sys.argv[7]
file_base = os.path.basename(filename)
file_ext  = os.path.splitext(file_base)[0]
filename2 = "/leo/savin/wide_lwasv/oims_wide/images/%d/%s.oims" %(mjd,file_ext)
print filename
print 'getting data from hdf5 files'

hf=h5py.File(filename,'r')
dset1 = hf.get('image')
dset2 = hf.get('time')
dset3 = hf.get('header')
hdr_dat  = dset3[0,:]
#intg =  dset1.shape[0]
#xSize = hdr[1]
#ySize = hdr[1]   
nchan = int(hdr_dat[3])
lat1= 34.0689         # coordinates for LWA station 1
lon1= -107.6284

lat2 = 34.3484       #cordinates for LWA-SV sation
lon2 = -106.8858

time_in = h*3600.00 + m*60.0 + s
dur  = 60.0
time_dat = dset2[:,1]*3600.00 + dset2[:,2]*60.0 + dset2[:,3]
peak = np.where((time_dat < time_in + dur) & (time_dat > time_in - dur))[0]

if len(peak) == 0:
   sys.exit(" File grabbed is not right")
      
def maxi_unsub(data,x_in,y_in):
    
    
    x_in=int(x_in)
    y_in=int(y_in)
    max=data[x_in,y_in]
    x1=x_in
    y1=y_in
    for l in range(x_in-1,x_in+1,1):
     for k in range(y_in-1,y_in+1,1):
      if data[l,k] > max:
       max=data[l,k]
       x1=l
       y1=k
 
    return x1,y1


#stokes = 0

wfall_source = np.zeros((len(peak),nchan))

db = OrvilleImageDB(filename2,'r')

for j,index in enumerate(peak):
    hdr, dat = db.__getitem__(index)
    im_dat = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))

    ngrid = int(db.header.ngrid) #size of images (x-axis)
    psize = db.header.pixel_size #angular size of a pixel (at zenith)
    nchan = int(db.header.nchan) # number of channels
    #f1 = hdr['start_freq']
    #bw = hdr['bandwidth']
    #f2 = hdr['stop_freq'] + bw
    t = hdr['start_time']
    #print hdr['lst']
    cent_az = hdr['center_az']
    cent_alt = hdr['center_alt']
    MJD = int(t)
    h = 24.0*(t - float(MJD))
    H = int(h)
    m = 60*(h - float(H))
    M = int(m)
    s = 60*(m - float(M))
    S = int(s)
    #f_chan_hz = np.arange(f1,f2,bw)
    im = dset1[index,stokes,:,:]

    xpix, ypix, altpix = trackcoord(ra, dec, MJD, H, M, S, ngrid, psize,lat2, lon2, cent_az, cent_alt)
    xa, ya = maxi_unsub(im, xpix, ypix)
    wfall_source[j,:] = im_dat[:,stokes,xa, ya]


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
    ax1.set_xlabel('Channels')
    ax1.set_ylabel('Response of Source')
    im1, = ax1.plot(range(wfall_source.shape[1]),wfall_source[i,:])

    im_int  = np.transpose(dset1[index,stokes,:,:])
    #im_dat = np.transpose(dset[i,80,stokes,:,:])
    
    vmin = sc(im_int,0)
    vmax = sc(im_int,100)
    #im = plt.pcolormesh(im_dat,cmap ='jet')
    im2 = ax2.imshow(im_int,cmap ='jet',vmin = vmin, vmax =vmax,origin ='lower')


    listims.append([im1, im2])

anim = animation.ArtistAnimation(fig, listims, interval=500, repeat_delay=1000)
#ani.save('movie_source.mp4')
anim.running = True
plt.show()

