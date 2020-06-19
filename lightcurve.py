import numpy as np
import sys
import matplotlib.pyplot as plt
import glob
import h5py
from windows import tukey
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
from utils import avgp, maxi_unsub
from scipy.signal import find_peaks, peak_widths, peak_prominences

filename = sys.argv[1]
print filename
print 'getting data from hdf5 files'

hf=h5py.File(filename,'r')
dset1 = hf.get('image')
dset2 = hf.get('time')
dset3 = hf.get('header')
hdr  = dset3[0,:]
intg =  dset1.shape[0]
xSize = hdr[1]
ySize = hdr[1]   

lat1= 34.0689         # coordinates for LWA station 1
lon1= -107.6284

lat2 = 34.3484       #cordinates for LWA-SV sation
lon2 = -106.8858



ra=float(sys.argv[2])
dec=float(sys.argv[3])
#MJD=float(sys.argv[4])
#h=float(sys.argv[5])
#m=float(sys.argv[6])
#s=float(sys.argv[7])
stokes= int(sys.argv[4])

size = hdr[1]
psize = hdr[2]
cent_az = hdr[9]
cent_alt = hdr[10]

k=0
light=np.zeros((intg,3),dtype=float)

for j in xrange(intg):

    data = dset1[j,stokes,:,:]
    time = dset2[j,:]
    xpix,ypix,altpix=trackcoord(ra, dec, time[0], time[1], time[2], time[3], size, psize, lat2, lon2, cent_az, cent_alt)
    #xpix,ypix,azpix,altpix=trackcoord_hrz(ra, dec, time[0], time[1], time[2], time[3], size, psize, lat2, lon2, cent_az, cent_alt)
    
    xa, ya = maxi_unsub(data,xpix,ypix)
   
    light[k,0] = j
    light[k,1] =  avgp(xa,ya,data)
    #light[k,1] =  avgpix(xa,ya,data)
    #light[k,1] =  data[xa,ya] 
    print light[k,0], light[k,1],time[1],time[2],time[3],xpix,ypix,xa,ya,altpix
    k=k+1


#y1 = np.fft.fft(light[:k,1])
y1 = np.fft.fft(light[:,1])
y1 = y1*(tukey(y1.shape[0],alpha = 0.015)**2)
data_filt = np.fft.ifft(y1)


# plotttingpart 

fig,ax= plt.subplots()
#ax.plot(light[:k,0],light[:k,1],color='red', marker='.', linestyle='solid')
ax.plot(range(data_filt.shape[0]),data_filt,color='blue', marker='.', linestyle='solid')
ax.set(xlabel= 'time (s)',ylabel= 'Arbitrary Power',title='Transient light curve')
ax.grid()
plt.show()


dat_noise = data_filt.real - np.mean(data_filt.real)
sig  = np.std(dat_noise)
#print sig
#print np.mean(dat_noise)

peaks,_ = find_peaks(dat_noise, height = 3.5*sig, distance = 3)
print peaks.shape
#prom = peak_prominences(dat_noise, peaks, wlen = 37 )
widths = np.zeros((len(peaks),4))
for i,pk in enumerate(peaks):
    width = peak_widths(dat_noise, np.array([pk]), rel_height= 1-((3.0*sig)/dat_noise[pk]))
    print width
    widths[i,:] = width
#width = peak_widths(dat_noise, peaks, rel_height = 0.85 )

#print widths

plt.hlines(widths[:,1],widths[:,2],widths[:,3], color="C2")
plt.plot(range(len(dat_noise)),dat_noise)
plt.scatter(peaks,dat_noise[peaks])
plt.show()
#print data_filt.real[242-2:242+2].max()/np.std(data_filt[range(242-55,242-5)+range(242+5,242+55)])
#print data_filt.real[242-2:242+2].max()/np.std(data_filt[:])



