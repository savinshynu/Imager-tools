"""
Filter nearby VLSS sources, low SNR scintillation events,
narrowband RFIs and airplane from the transient search list in Stokes I

"""


import time
import os
import numpy as np
import sys
import math
import glob
import matplotlib.pyplot as plt
from low_snr import filter_snr
import h5py
from narrow_rfi import flag_narrow_rfi as fnr
from utils import ang_dist
from ap_filter import filter_planes
day = sys.argv[1]

t0 = time.time() 
vlss = np.loadtxt("/leo/savin/wide_lwasv/oims_wide/standard_sub/VLSS_10.txt") 

"""
#function to search the coordinate in VLSS sources
def search_vlss(ra,dec):
    su=0
    for z in range(vlss.shape[0]):
        diff = abs(ang_dist(ra,dec,vlss[z,0],vlss[z,1]))
        if (vlss[z,2] > 50.0 and diff < 3.0): 
           su += 1
    if su > 0:
       return "true"
    else:
       return "false"
"""

def search_vlss(ra,dec):
    ra_in = np.ones(vlss.shape[0])*ra
    dec_in = np.ones(vlss.shape[0])*dec
    diff = ang_dist(ra_in, dec_in, vlss[:,0], vlss[:,1])
    #print diff
    par = np.where((vlss[:,2] > 50.0) & (diff < 3.0))[0]
    #print par 
    return len(par)

path1 = "/leo/savin/wide_lwasv/oims_wide/standard_sub/transients/"+str(int(day))+"/*.txt"
files1 = sorted(glob.glob(path1))

path2 = "/leo/savin/wide_lwasv/oims_wide/standard_sub/corr_transients/"+str(int(day))+"/"
try:
    os.mkdir(path2)
except OSError:
    print 'directry already exists'

stokes = 0

for filename1 in files1:
    
    print filename1
    list1 = np.loadtxt(filename1)
    dim = list1.shape
    if len(dim) == 1:
       list1 = np.reshape(list1,(1,dim[0]))
       
    list2 = np.zeros(list1.shape)
   
    k = 0
    for i in range(list1.shape[0]):
        if search_vlss(list1[i,0],list1[i,1]) == 0:
           list2[k,:] = list1[i,:]
           k += 1
    list2 = list2[:k,:]
    
    
    list3 = filter_snr(list2,day)
    #np.savetxt(path2+"unpol_"+os.path.basename(filename1), list3,'%10.2f')
    #np.savetxt(path2+"pol_"+os.path.basename(filename1), list4,'%10.2f')
    #list3 = np.loadtxt(path2+"unpol_"+os.path.basename(filename1))
    
     
    ext = os.path.basename(filename1)
    ext_num = int(ext.split("_")[1][:2]) 
 
    list4 = fnr(list3,ext_num,day)
    
    list5 = filter_planes(list4,day,stokes)
    
    np.savetxt(path2+os.path.basename(filename1), list5,'%10.2f')
    #np.savetxt(path2+"unpol_"+os.path.basename(filename1), list6,'%10.2f')
    #np.savetxt(path2+"pol_"+os.path.basename(filename1), list7,'%10.2f')
    

t1 = time.time()

print (t1-t0)/60.0    
