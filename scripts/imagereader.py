"""
Plots the difference between the narrowband and broadband images

"""

import numpy as np
from OrvilleImageDB import OrvilleImageDB
import sys
import matplotlib.pyplot as plt
#from exp_functions import pix2eq,hrz2pix,eq2pix,sin_pix2sky,sin_sky2pix
#from function_ken import pix2eq,hrz2pix,eq2pix
from conversion_coord import pix2eq, eq2pix 
from scipy.stats import scoreatpercentile as sc
from mad import median_absolute_deviation as mad
from rfi_flagger_oims2hdf5 import main
import matplotlib.animation as animation

filename = sys.argv[1]

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
#print ints,ngrid,psize,nchan #,station,stokes,inp_flag,file_start,file_end



time = np.zeros((50,4))   #((ints,4))
data = np.zeros((50,nchan,4,ngrid,ngrid))    #((ints,nchan,4,ngrid,ngrid))
#stokes = 0 # right now we only do Stokes I, but we might want to do all 4 parameters, or at least I and V


lat1= 34.0689         # coordinates for LWA station 1
lon1= -107.6284

lat2 = 34.3484       #cordinates for LWA-SV sation
lon2 = -106.8858



for i in xrange(50):
    hdr, dat = db.read_image()
    t = hdr['start_time']
    int_len = hdr['int_len']
    lst = hdr['lst']
    start_freq = hdr['start_freq']/1e+6
    stop_freq = hdr['stop_freq']/1e+6 
    bandwidth = hdr['bandwidth']/1e+6
    cent_ra = hdr['center_ra']
    cent_dec = hdr['center_dec']
    cent_az = hdr['center_az']    
    cent_alt = hdr['center_alt']

    #print int_len,lst,start_freq,stop_freq,bandwidth,cent_az,cent_alt
    #print cent_ra,cent_dec,cent_az,cent_alt

    MJD = int(t)
    h = 24.0*(t - float(MJD))
    H = int(h)
    m = 60*(h - float(H))
    M = int(m)
    s = 60*(m - float(M))
    S = int(s)
    
    #xa,ya,alt_pix = eq2pix(np.array([350.8,299.8]),np.array([58.8,40.7]),np.array([MJD,MJD]),np.array([H,H]),np.array([M,M]),np.array([S,S]),ngrid,cent_az,cent_alt,psize,lat2,lon2)
    #print xa,ya,alt_pix

    #xa, ya = sin_sky2pix(np.array([350.8,299.8]), np.array([58.8,40.7]), ngrid, psize, cent_ra, cent_dec)
    #print xa,ya

    #ra, dec  = sin_pix2sky(xa, ya, ngrid, psize,cent_ra,cent_dec) 
    #print ra , dec 
    #ra,dec = pix2eq(xa,ya,ngrid,cent_az,cent_alt,psize,np.array([MJD,MJD]),np.array([H,H]),np.array([M,M]),np.array([S,S]),lat2,lon2)
    #print (ra),(dec)

    time[i,:] = np.reshape(np.array([MJD,H,M,S]),(4,))
    data[i,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
    print i, time[i,:]


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
ax1.set_title("Single 100 kHz Channel")
ax2=fig.add_subplot(1,2,2)
ax2.set_title("198 Channels combined")
listims = []
fig.canvas.mpl_connect('key_press_event', on_press)

for index in range(50) :

    im_i = np.transpose(data[index,79,0,:,:])
    imin = sc(im_i,0)
    imax = sc(im_i,100)
    im1 = ax1.imshow(im_i,cmap ='jet',vmin = imin, vmax = imax, origin ='lower')

    im_v  = np.transpose(np.mean(data[index,:,0,:,:],axis = 0))
    vmin = sc(im_v,0)
    vmax = sc(im_v,100)
    im2 = ax2.imshow(im_v,cmap ='jet',vmin = vmin, vmax =vmax,origin ='lower')


    listims.append([im1, im2])

anim = animation.ArtistAnimation(fig, listims, interval=500, repeat_delay=1000)
anim.save('movie_broadband.mp4')
anim.running = True
plt.show()





"""

#print data[:,:,5]
#for x in range(ints):

test_dat = np.mean(data[0,:,0,:,:],axis =0)

coord_az = np.zeros((128,128))
coord_alt = np.zeros((128,128))

for x in range(128):
    for y in range(128):
        az,alt = pix2eq(np.array([x,40]),np.array([y,40]),ngrid,cent_az,cent_alt,psize,np.array([MJD,MJD]),np.array([H,H]),np.array([M,M]),np.array([S,S]),lat2,lon2)
        #az,alt = sin_pix2sky(x, y, ngrid, psize, cent_ra, cent_dec)
        #print az[0]*180.0/3.14
        coord_az[x,y] = az[0]*180.0/3.14 
        coord_alt[x,y] = alt[0]*180.0/3.14 
#plt.plot(range(coord.shape[0]),coord[50,:])
#plt.show()



plt.imshow(coord_az)
plt.colorbar()
plt.title("Azimuth")
plt.legend()
plt.show()

plt.imshow(coord_alt)
plt.colorbar()
plt.title("Elevation")
plt.legend()
plt.show()

stokes = 0
im_avg = np.transpose(np.mean(data[1,:,stokes,:,:],axis =(0)))
vmin = sc(im_avg,0.5)
vmax = sc(im_avg,99.5)
plt.imshow(im_avg,vmin = vmin, vmax = vmax, cmap='jet',origin = 'lower')
plt.show()


# print time[:,x]
plt.pcolormesh(np.transpose(data[0,144,0,:,:]))
plt.show()

std_dat = np.zeros(data.shape[1])

for x in range(data.shape[0]):

    im  = data[x,:,3,:,:]
    avg = np.mean(im,axis = 0)
    print np.std(avg)

chan_im = data[1,:,stokes,:,:]
for x in range(chan_im.shape[0]):

    im  = chan_im[x,:,:]
    std_dat[x] = mad(im, axis=None) #np.std(im)

plt.plot(range(len(std_dat)),std_dat)
plt.show()

gd = main(std_dat)

print np.mean(std_dat),np.mean(std_dat[gd])

im_avg_fl = np.transpose(np.mean(data[1,gd,stokes,:,:],axis =(0)))
vmin = sc(im_avg_fl,0.5)
vmax = sc(im_avg_fl,99.5)
plt.imshow(im_avg_fl,vmin = vmin, vmax = vmax, cmap='jet', origin = 'lower')
plt.show()

#making movie
fig = plt.figure()
listims = []

for i in range(chan_im.shape[0]):

    im_dat = np.transpose(chan_im[i,:,:])
    vmin = sc(im_dat,0)
    vmax = sc(im_dat,100 )
    im = plt.imshow(im_dat,cmap ='jet',vmin = vmin, vmax =vmax,origin ='lower')


    listims.append([im])

ani = animation.ArtistAnimation(fig, listims, interval=500, repeat_delay=1000)
#ani.save('movie_vir_scint.mp4')
plt.show()

"""
