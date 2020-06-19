import numpy as np
import sys, glob
import h5py
import matplotlib.pyplot as plt

filename = sys.argv[1]

hf = h5py.File(filename,'r')

im_dat =  hf.get('image')
t_dat = hf.get('time')
head_dat = hf.get('header')
el_dat = hf.get('elev')

print t_dat[:10,:]


cyg, cas, vir, tau = [0,1,2,3]
stokes = 0
source = cyg
ind  = np.where((el_dat[:,source] > 10))[0] # &(np.int_(head_dat[:,0])==30))[0]
print ind
ydat = np.mean(im_dat[ind,source,:,stokes],axis = 1)
xdat = el_dat[ind,source]

print head_dat[0,0]


plt.scatter(xdat,ydat)
plt.xlabel("Elevation in degrees")
plt.ylabel("Power in a.u.")
plt.show()

plt.plot(range(im_dat.shape[2]),np.mean(im_dat[ind,source,:,stokes],axis=0))
plt.show()
