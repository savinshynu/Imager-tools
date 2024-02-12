"""
Scripts to test the spectra results from Orville
"""

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
#from sliding_rfi_flagger_man import main
#from lsl import astro
#import glob, os
#import sys
#import h5py
#import pyfits as pf
#from scipy.optimize import curve_fit
#from mad import median_absolute_deviation as mad
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord,eq2hrz,pix2hrz
from OrvilleImageDB import OrvilleImageDB
from utils import mask_image,maxi_sub,maxi_unsub,read_oims_cont,avgp,avgp_chan,ang_dist
#from collect_resp_data import collect_resp
#import matplotlib.animation as animation
from twod_gaussian import fitgaussian
from plot_antenna import response_sing_antenna
from instrumental_response import getImpedanceMisMatch,getARXResponse
from scipy.stats import scoreatpercentile as sc

def avgp_chan_eq(ra,dec,chan_info,eq_info,data,el):

    #data in the form of channels,pixels,pixels
    #chan_info has frequency values of each channels in MHz
     
    resp = np.zeros(data.shape[0]) # initilize the output result
    for chan in range(data.shape[0]): # Looking at each channel
        fwhm = 2.2*(74/(chan_info[chan]/1e+6)) # calculating fwhm and beam area at 63 degrees
        sig1 = fwhm/2.355
        sig2 = sig1/np.sin(el*np.pi/180.0)# change elevation wrt sources using
        area = np.pi*sig1*sig2
        rad_in = 1.2*fwhm # defining radius of the annulus
        rad_out = 2.2*fwhm
        
        ins = [] # array to hold pixel values inside the source area
        annul =[] # array to hold the pixel values inside the annulus
        ins_im = np.zeros((data.shape[1],data.shape[2])) # for visualizing inside source and annulus
        ann_im = np.zeros((data.shape[1],data.shape[2]))
        for x in range(eq_info.shape[0]):
            for y in range(eq_info.shape[1]):
                if ang_dist(eq_info[x,y,0],eq_info[x,y,1],ra,dec) < rad_in :  # looking with in inside radius
                   ins.append(data[chan,x,y])
                   ins_im[x,y] = data[chan,x,y]
                elif (ang_dist(eq_info[x,y,0],eq_info[x,y,1],ra,dec) < rad_out) and (ang_dist(eq_info[x,y,0],eq_info[x,y,1],ra,dec)> rad_in):
                   annul.append(data[chan,x,y]) # looking within annulus
                   ann_im[x,y] = data[chan,x,y]

        xmin = np.where((ins_im != 0))[0].min() # finding the extremum of pixel values to fit a gaussian
        xmax = np.where((ins_im != 0))[0].max()
        ymin = np.where((ins_im != 0))[1].min()
        ymax = np.where((ins_im != 0))[1].max()
        params = fitgaussian(ins_im[xmin:xmax+1,ymin:ymax+1]) # height,x,y,sigmax,sigmay from gaussian fit
        if chan ==0: #or chan ==197:
           plt.pcolormesh(ins_im,cmap='jet')
           plt.show()
           plt.pcolormesh(ann_im,cmap='jet')
           plt.show()
        #print np.median(annul)
        #plt.hist(data[chan,40:80,70:110])
        #plt.show()
        ins = np.array(ins)  # converting list to array
        annul = np.array(annul)
        resp[chan] = (params[0] - np.median(annul))*area
        #resp[chan] = ins.sum() - (np.median(annul)*len(ins))
        #print ins.max(), np.median(annul),np.median(data[chan,:,:])
    return resp

cyg, cas, vir, tau = [0,1,2,3]

par_jy = [[4.695,0.085,-0.178],[5.745,-0.770,0],[5.023,-0.856,0],[3.915,-0.299,0]]

def flux_jy_baars(v,source): # From Baars et al. paper values

    a,b,c = par_jy[source]

    log_s = a + b*np.log10(v) + c*(np.log10(v)**2)

    return 10**log_s

   

filename2 = "/leo/savin/wide_lwasv/oims_wide/images/58889/58889_190000_30.137MHz_49.837MHz.oims"
nchan = 198
ngrid = 128
peak = [552]
#noise = range(94,100) #+range(251,256)

dset = np.zeros((len(peak),nchan,4,ngrid,ngrid))
#dset_avg = np.zeros((len(noise),nchan,4,ngrid,ngrid))
db = OrvilleImageDB(filename2,'r')

print db.nint
for j,index in enumerate(peak):
    hdr, dat = db.__getitem__(index)
    dset[j,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))


ngrid = int(db.header.ngrid) #size of images (x-axis)
psize = db.header.pixel_size #angular size of a pixel (at zenith)
nchan = int(db.header.nchan) # number of channels
f1 = hdr['start_freq']
bw = hdr['bandwidth']
f2 = hdr['stop_freq'] + bw
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

f_chan_hz = np.arange(f1,f2,bw)

flux_cas_jy = flux_jy_baars(f_chan_hz/1e+6,cas)    
flux_cyg_jy = flux_jy_baars(f_chan_hz/1e+6,cyg)

"""
for j,index in enumerate(noise):
    hdr, dat = db.__getitem__(index)
    dset_avg[j,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
"""
dset = np.mean(dset, axis = 0)
#dset_avg = np.mean(dset_avg, axis = 0)
sub_dat = dset #- dset_avg
stokes = 0

#For making fits file
#hdu = pf.PrimaryHDU(np.transpose(sub_dat[:,:,:,:],axes=(0,1,3,2)))
#hdulist = pf.HDUList([hdu])
#hdu.writeto('image_test.fits')


plt.plot(range(nchan),sub_dat[:,stokes,53,88])
plt.show()

im = np.mean(sub_dat[:,stokes,:,:],axis =0)
im_chan = sub_dat[100,stokes,:,:]
#plt.pcolormesh(np.transpose(im),cmap ='jet')
plt.pcolormesh(np.transpose(im_chan),cmap ='jet')
plt.show()

ra_cyg, dec_cyg = [299.8, 40.7]
ra_cas, dec_cas = [350.8, 58.8] 

cal_dat = sub_dat[:,stokes,:,:]
cal_dat_mod = np.zeros(cal_dat.shape)
#annul_info = [ngrid, psize, MJD, H, M, S, lat, lon, cent_az, cent_alt]
eq_info = np.zeros((cal_dat.shape[1],cal_dat.shape[2],2))
hrz_info = np.zeros((cal_dat.shape[1],cal_dat.shape[2],2))
for l in range(cal_dat.shape[1]):
    for k in range(cal_dat.shape[2]):
        m_ra,m_dec = getcoord(np.array([l]),np.array([k]),ngrid, psize, np.array([MJD]), np.array([H]), np.array([M]), np.array([S]), lat, lon, cent_az, cent_alt)
        eq_info[l,k,0] = m_ra
        eq_info[l,k,1] = m_dec
        m_az, m_alt = pix2hrz(np.array([l]), np.array([k]), ngrid, psize, cent_az, cent_alt)
        hrz_info[l,k,0] = m_az
        hrz_info[l,k,1] = m_alt

#plt.imshow(hrz_info[:,:,1])
#plt.show()

el_cas = hrz_info[54,88,1]  # for cas 
el_cyg = hrz_info[86,74,1]  # for cyg

print el_cas, el_cyg

imp_resp = getImpedanceMisMatch(f_chan_hz)
imp_resp /= imp_resp.max()
arx_resp = getARXResponse(f_chan_hz, filter='split')
arx_resp /= arx_resp.max()
plt.plot(range(nchan),arx_resp)
plt.show()

ant_resp = np.zeros((nchan,ngrid,ngrid))
for chan in range(nchan):
    ant_resp_chan = response_sing_antenna(f_chan_hz[chan]/1e+6,hrz_info[:,:,0],hrz_info[:,:,1])
    ant_resp_chan /= np.nanmax(ant_resp_chan)
    cal_dat_mod_chan = cal_dat[chan,:,:]
    cal_dat_mod_chan[np.isnan(ant_resp_chan)] = 1
    ant_resp_chan[np.isnan(ant_resp_chan)] = 1
    ant_resp[chan,:,:] = ant_resp_chan
    #cal_dat_mod_chan /= imp_resp[chan]
    #cal_dat_mod_chan /= arx_resp[chan]
    #cal_dat_mod_chan /= ant_resp_chan
    cal_dat_mod[chan,:,:] = cal_dat_mod_chan

#plt.imshow(ant_resp[0,:,:])
#plt.show()

im_new = np.mean(cal_dat_mod,axis =0)
im_new_chan = cal_dat_mod[100,:,:]
vmin = sc(im_new,5)
vmax = sc(im_new,99)

plt.imshow(np.transpose(im_new),vmin = vmin,vmax =vmax)
plt.show()

#print im_new_chan[54,88], im_new_chan[86,74]

resp1 = avgp_chan_eq(ra_cas,dec_cas,f_chan_hz,eq_info,cal_dat_mod,el_cas)
resp2 = avgp_chan_eq(ra_cyg,dec_cyg,f_chan_hz,eq_info,cal_dat_mod,el_cyg)
"""
ratio_cas_jy = resp1/flux_cas_jy
ratio_cyg_jy = resp2/flux_cyg_jy
plt.plot(range(nchan),ratio_cas_jy,label = 'Ratio of Cas response wrt Jy')
plt.plot(range(nchan),ratio_cyg_jy, label = 'Ratio of Cyg response wrt Jy')
plt.title('Cyg is higher than  Cas')
plt.xlabel('Channels')
plt.ylabel('Ratio in (a.u./Jy)')
plt.legend()
plt.show()
"""


scale1 = resp1[197]
#resp1 = (resp1/scale1)
#resp2 = (resp2/scale1)
ratio1 = resp1/resp2

plt.plot(range(len(resp1)),resp1,label = 'Cas') # From the annulus method
plt.plot(range(len(resp2)),resp2,label = 'Cyg')
plt.xlabel('Channels')
plt.ylabel('Instrumental response corrected (a.u.)')
plt.title("Annulus method")
plt.legend()
plt.show()

resp3 = cal_dat_mod[:,54,88] #for cas
resp4 = cal_dat_mod[:,86,74] # for cyg   # from the maxmimum pixel value
scale2 = resp3[197]
#resp3 /= scale2
#resp4 /= scale2
ratio2 = resp3/resp4

plt.plot(range(len(resp3)),resp3,label = 'Cas') # From the peak value method
plt.plot(range(len(resp4)),resp4,label = 'Cyg')
plt.title("Peak value")
plt.xlabel('Channels')
plt.ylabel('Instrumental response corrected(a.u.)')
plt.legend()
plt.show()

ratio3 = flux_cas_jy/flux_cyg_jy #Actual value

plt.plot(range(len(ratio1)),ratio1,label='Annulus subtracted and integrated over point source')
plt.plot(range(len(ratio2)),ratio2,label="peak value from point source")
plt.plot(range(len(ratio3)),ratio3,label="Values from Baars et al paper")
plt.xlabel('channels')
plt.ylabel('Ratio of Cas with Cyg')
#plt.ylim(0,150)
plt.legend()
plt.show()

