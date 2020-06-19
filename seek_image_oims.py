import numpy as np
from OrvilleImageDB import OrvilleImageDB
import sys
import matplotlib.pyplot as plt
#from exp_functions import pix2eq,hrz2pix,eq2pix,sin_pix2sky,sin_sky2pix
#from function_ken import pix2eq,hrz2pix,eq2pix
from conversion_coord import pix2eq, eq2pix 
from mad import median_absolute_deviation as mad
from utils import avgp_mean, avgp_med,avgp_test

filename = sys.argv[1]

db = OrvilleImageDB(filename,'r') #define the image data base

print db.__len__()

#getting parameters from the input file

ints = db.nint #number of integration
station =  db.header.station
#stokes = db.header.stokes_params
inp_flag = db.header.flags
file_start = db.header.start_time
file_end = db.header.stop_time  


ngrid = db.header.ngrid #size of images (x-axis)
psize = db.header.pixel_size #angular size of a pixel (at zenith)
nchan = db.header.nchan # number of channels
#print ints,ngrid,psize,nchan #,station,stokes,inp_flag,file_start,file_end



#time = np.zeros((ints,4))
#data = np.zeros((ints,nchan,4,ngrid,ngrid))
stokes = 0 # right now we only do Stokes I, but we might want to do all 4 parameters, or at least I and V


lat1= 34.0689         # coordinates for LWA station 1
lon1= -107.6284

lat2 = 34.3484       #cordinates for LWA-SV sation
lon2 = -106.8858



"""
hdr, dat = db.__getitem__(710)
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
    

#time[i,:] = np.reshape(np.array([MJD,H,M,S]),(4,))
#data[i,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
#print i, time[i,:]

print MJD,H,M,S
"""
peak = np.array(range(694,695))
noise = np.array(range(680,690)+range(700,710))
dset = np.zeros((peak.shape[0],nchan,4,ngrid,ngrid))
dset_avg = np.zeros((noise.shape[0],nchan,4,ngrid,ngrid))
    
for index in range(peak.shape[0]):
    hdr, dat = db.__getitem__(peak[index])
    dset[index,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
    #print 1.486*mad(np.mean(dset[index,:,0,:,:],axis=0),axis=None)

for index in range(noise.shape[0]):
    hdr, dat = db.__getitem__(noise[index])
    dset_avg[index,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
    #print 1.486*mad(np.mean(dset_avg[index,:,0,:,:],axis=0),axis=None)

sub_dat = np.mean(dset[:,:,stokes,:,:],axis=(0)) - np.mean(dset_avg[:,:,stokes,:,:],axis=(0))
dat_img = np.mean(sub_dat,axis=0)


plt.plot(range(sub_dat.shape[0]),sub_dat[:,97,59])
plt.show()

#print dat_img[36,56]/(1.486*mad(dat_img,axis=None)), 1.486*mad(dat_img,axis=None)

#print avgp_med(97,71,dat_img), avgp_mean(97,71,dat_img), mad(dat_img,axis = None)
print avgp_test(97,71,dat_img)

plt.pcolormesh(np.transpose(dat_img),cmap ='jet')
plt.show()


#plt.pcolormesh(np.transpose(np.mean(dat[:,0,:,:],axis =0)))
#plt.show()

#plt.pcolormesh(np.transpose(dat[100,0,:,:]))
#plt.show()



