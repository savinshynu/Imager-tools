import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord,eq2hrz,pix2hrz
from OrvilleImageDB import OrvilleImageDB
from utils import mask_image,maxi_sub,maxi_unsub,read_oims_cont,avgp,ang_dist
from twod_gaussian import fitgaussian
from plot_antenna import response_sing_antenna
from instrumental_response import getImpedanceMisMatch,getARXResponse
from scipy.stats import scoreatpercentile as sc

def avgp_chan_eq(ra,dec,chan_info,eq_info,data,el):
    #data in the form of channels,pixels,pixels
    #chan_info has frequency values of each channels in kHz
    resp = np.zeros(data.shape[0]) # initilize the output result
    for chan in range(data.shape[0]): # Looking at each channel
        fwhm = 2.2*(74/(chan_info[chan]/1e+6)) # calculating fwhm and beam area 
        sig1 = fwhm/2.355
        sig2 = sig1/np.sin(el*np.pi/180.0)# change elevation wrt sources using
        area = np.pi*sig1*sig2
        rad_in = 1.5*fwhm #1.5 defining radius of the annulus
        rad_out = 2.5*fwhm #2.5
        
        ins = [] # list to hold pixel values inside the source area
        annul =[] # list to hold the pixel values inside the annulus
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
        """
        xmin = np.where((ins_im != 0))[0].min() # finding the extremum of pixel values for anulus to fit a gaussian
        xmax = np.where((ins_im != 0))[0].max()
        ymin = np.where((ins_im != 0))[1].min()
        ymax = np.where((ins_im != 0))[1].max()
        params = fitgaussian(ins_im[xmin:xmax+1,ymin:ymax+1]) # height,x,y,sigmax,sigmay from gaussian fit
                   
        if chan ==10:  #or chan ==197:
           #print  params[0]
           plt.imshow(ins_im,cmap='jet',origin = 'lower')
           plt.show()
           plt.imshow(ann_im,cmap='jet',origin = 'lower')
           plt.show()
              
        """

        ins = np.array(ins)  # converting list to array
        annul = np.array(annul)
        #resp[chan] = (params[0] - np.median(annul)) #*area
        #resp[chan] = ins.sum() - (np.mean(annul)*len(ins))
        resp[chan] = ins.max() - np.mean(annul)
    return resp

cyg, cas, vir, tau = [0,1,2,3] #defining source values

par_jy = [[4.695,0.085,-0.178],[5.745,-0.770,0],[5.023,-0.856,0],[3.915,-0.299,0]] #a,b,c fit parameters for each source from baars etal paper

def flux_jy_baars(v,source): # From Baars et al. paper values

    a,b,c = par_jy[source]

    log_s = a + b*np.log10(v) + c*(np.log10(v)**2)

    return 10**log_s

   

filename2 = "/leo/savin/wide_lwasv/oims_wide/images/58888/58888_210000_30.137MHz_49.837MHz.oims"
nchan = 198
ngrid = 128

peak = range(330,350) #(0,20)#[552]


dset = np.zeros((len(peak),nchan,4,ngrid,ngrid))
db = OrvilleImageDB(filename2,'r')

for j,index in enumerate(peak):
    hdr, dat = db.__getitem__(index)
    dset[j,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
"""
#plt.imshow(dset[:,:,0,54,88],origin = 'lower')
dt = dset[:,:,0,54,88]
med = np.median(dset[:,:,0,54,88],axis = 0)
wfall = dt*0.0
for i in range(dset.shape[0]):
    wfall[i,:] = dt[i,:]/med
plt.pcolormesh(wfall)
plt.colorbar()
plt.show()
"""
ngrid = int(db.header.ngrid) #size of images (x-axis)
psize = db.header.pixel_size #angular size of a pixel (at zenith)
nchan = int(db.header.nchan) # number of channels
f1 = hdr['start_freq']
bw = hdr['bandwidth']
f2 = hdr['stop_freq'] + bw
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

lat = 34.3484       #cordinates for LWA-SV sation
lon = -106.8858

f_chan_hz = np.arange(f1,f2,bw)

flux_cas_jy = flux_jy_baars(f_chan_hz/1e+6,cas)    
flux_cyg_jy = flux_jy_baars(f_chan_hz/1e+6,cyg)
flux_vir_jy = flux_jy_baars(f_chan_hz/1e+6,vir)


dset = np.mean(dset, axis = 0)
sub_dat = dset #- dset_avg
stokes = 0

im = np.mean(sub_dat[:,stokes,:,:],axis =0) #averaging over channels to get a mean image to visualise
im_chan = sub_dat[1,stokes,:,:]
#plt.pcolormesh(np.transpose(im),cmap ='jet')
#plt.pcolormesh(np.transpose(im_chan),cmap ='jet')
#plt.show()

ra_cyg, dec_cyg = [299.8, 40.7] # ra dec values for cyg and cas
ra_cas, dec_cas = [350.8, 58.8] 
ra_vir, dec_vir = [187.70, 12.39]
ra_tau, dec_tau = [83.6, 22.01]

x_cas, y_cas, alt_cas = trackcoord(ra_cas, dec_cas, MJD, H, M, S, ngrid, psize, lat, lon, cent_az, cent_alt)
xa_cas, ya_cas = maxi_unsub(im, x_cas, y_cas)
print 'x,y, alt value for Cas are %f  %f %f ' % (xa_cas, ya_cas, alt_cas)

x_cyg, y_cyg, alt_cyg = trackcoord(ra_cyg, dec_cyg, MJD, H, M, S, ngrid, psize, lat, lon, cent_az, cent_alt)
xa_cyg,ya_cyg = maxi_unsub(im, x_cyg, y_cyg)
print 'x,y, alt value for Cyg are %f  %f  %f ' % (xa_cyg, ya_cyg, alt_cyg)

x_vir, y_vir, alt_vir = trackcoord(ra_vir, dec_vir, MJD, H, M, S, ngrid, psize,lat, lon, cent_az, cent_alt)
xa_vir, ya_vir = maxi_unsub(im, x_vir, y_vir)
print 'x,y, alt value for Vir are %f  %f  %f ' % (xa_vir, ya_vir, alt_vir)


x_tau, y_tau, alt_tau = trackcoord(ra_tau, dec_tau, MJD, H, M, S, ngrid, psize,lat, lon, cent_az, cent_alt)
xa_tau, ya_tau = maxi_unsub(im, x_tau, y_tau)



#plt.plot(range(nchan),sub_dat[:,0,xa_cas, ya_cas])
#plt.show()

cal_dat = sub_dat[:,stokes,:,:]
cal_dat_mod = np.zeros(cal_dat.shape)

#Taking the stokes I image and converting the corresponding pixel values to (azimuth,altitude) and (ra,dec) values
eq_info = np.zeros((cal_dat.shape[1],cal_dat.shape[2],2)) #ra,dec values
hrz_info = np.zeros((cal_dat.shape[1],cal_dat.shape[2],2)) # azimuth, altitude values
# getcoord and pix2hrz taken from conversion_coord.py
for l in range(cal_dat.shape[1]):
    for k in range(cal_dat.shape[2]):
        m_ra,m_dec = getcoord(np.array([l]),np.array([k]),ngrid, psize, np.array([MJD]), np.array([H]), np.array([M]), np.array([S]), lat, lon, cent_az, cent_alt)
        eq_info[l,k,0] = m_ra
        eq_info[l,k,1] = m_dec
        m_az, m_alt = pix2hrz(np.array([l]), np.array([k]), ngrid, psize, cent_az, cent_alt)
        hrz_info[l,k,0] = m_az
        hrz_info[l,k,1] = m_alt

#plt.imshow(hrz_info[:,:,1]) #plots the altitude distribution of each pixel
#plt.show()

el_cas = hrz_info[xa_cas,ya_cas,1]  # passing azi and elevation info for cas, this is used to estimate the beam width which changes by sin(elevation),(54,88) are the pixel values from image
el_cyg = hrz_info[xa_cyg,ya_cyg,1]  # same for cyg

el_vir = hrz_info[xa_vir,ya_vir,1]

#print el_cas, el_cyg

imp_resp = getImpedanceMisMatch(f_chan_hz) #calling impedance mismatch function
#imp_resp /= imp_resp.max()
arx_resp = getARXResponse(f_chan_hz, filter='split') # calling arx response function
#arx_resp /= arx_resp.max()
ant_resp_chan = np.zeros(nchan)

for chan in range(nchan):
    ant_resp_chan[chan] = response_sing_antenna(f_chan_hz[chan]/1e+6,hrz_info[98,81,0],hrz_info[98,81,1])
comb_resp = arx_resp*imp_resp*ant_resp_chan

#plt.plot(range(nchan),arx_resp)
plt.plot(f_chan_hz/1e+6,comb_resp/comb_resp.max(),label = 'Combined Instrumental Response')
plt.plot(f_chan_hz/1e+6,cal_dat[:,98,81]/cal_dat[:,98,81].max(),label ='Cyg')
plt.xlabel('Frequency(MHz)')
plt.ylabel('Power arb')
plt.legend()
plt.show()

ant_resp = np.zeros((nchan,ngrid,ngrid))
for chan in range(nchan):
    #ant_resp_chan = response_sing_antenna(f_chan_hz[chan]/1e+6,hrz_info[:,:,0],hrz_info[:,:,1]) # single antenna response in stokes I as a function of frequency,azimuth and elevation
    #ant_resp_chan /= np.nanmax(ant_resp_chan)
    cal_dat_mod_chan = cal_dat[chan,:,:]
    #cal_dat_mod_chan[np.isnan(ant_resp_chan)] = 1 #putting pixel value beyond horizon as 1
    #ant_resp_chan[np.isnan(ant_resp_chan)] = 1   # same for antenna reponse
    #ant_resp[chan,:,:] = ant_resp_chan
    #cal_dat_mod_chan /= imp_resp[chan]    # applying different calibrations
    #cal_dat_mod_chan /= arx_resp[chan]
    #cal_dat_mod_chan /= ant_resp_chan
    cal_dat_mod[chan,:,:] = cal_dat_mod_chan

#plt.imshow(ant_resp[0,:,:])
#plt.show()

#im_new = np.mean(cal_dat_mod,axis =0)  # Looking at instrumental response corrected mean image
#im_new_chan = cal_dat_mod[100,:,:]
#vmin = sc(im_new,5)
#vmax = sc(im_new,99)

#plt.imshow(np.transpose(im_new),vmin = vmin,vmax =vmax,origin = 'lower')
#plt.show()

#print 'peak value  Gaussian height # Cas'
#resp1 = avgp_chan_eq(ra_cas,dec_cas,f_chan_hz,eq_info,cal_dat_mod,el_cas) #Getting the response function for Cas
#print 'peak value  Gaussian height # Cyg'
resp2 = avgp_chan_eq(ra_cyg,dec_cyg,f_chan_hz,eq_info,cal_dat_mod,el_cyg) # same for cyg

resp1 = avgp_chan_eq(ra_cas,dec_cas,f_chan_hz,eq_info,cal_dat_mod,el_cas)
"""
# Finding the ratio of response in a.u. to actual Jy to see how much correction is needed
#ratio_cas_jy = resp1/flux_cas_jy
#ratio_cyg_jy = resp2/flux_cyg_jy
#plt.plot(range(nchan),ratio_cas_jy,label = 'Ratio of Cas response wrt Jy')
#plt.plot(range(nchan),ratio_cyg_jy, label = 'Ratio of Cyg response wrt Jy')
#plt.title('Cyg is higher than  Cas')
#plt.xlabel('Channels')
#plt.ylabel('Ratio in (a.u./Jy)')
#plt.legend()
#plt.show()
"""


scale1 = resp1[197]
#resp1 = (resp1/scale1) # Scaling reponse from Cyg and Cas to same level
#resp2 = (resp2/scale1)


ratio1 = resp1/resp2 #ratio of responses

plt.plot(range(len(resp1)),resp1,label = 'Cas') # From the annulus method
plt.plot(range(len(resp2)),resp2,label = 'Cyg')
plt.xlabel('Channels')
plt.ylabel('Instrumental response corrected (a.u.)')
plt.title("Annulus method")
plt.legend()
plt.show()

resp3 = cal_dat_mod[:,xa_cas,ya_cas] #for cas     (54,88) and (86,74) are the pixel locations of cyg and cas at the particular time
resp4 = cal_dat_mod[:,xa_cyg,ya_cyg] # for cyg   # from the maxmimum pixel value
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

ratio3 = flux_cas_jy/flux_cyg_jy #Actual flux value ratios for source 
ratio4 = resp1/resp4
# Now plotting the ratios of cas and cyg from different methods together

plt.plot(range(len(ratio1)),ratio1,label='Annulus subtracted and integrated over point source')
plt.plot(range(len(ratio2)),ratio2,label="peak value from point source")
plt.plot(range(len(ratio4)),ratio4,label="exp")
plt.plot(range(len(ratio3)),ratio3,label="Values from Baars et al paper")
plt.xlabel('channels')
plt.ylabel('Ratio of Cas with Cyg')
plt.legend()
plt.show()



plt.loglog((f_chan_hz/1e+6),(ratio1*flux_cyg_jy),label='Annulus subtracted and integrated over point source')
plt.loglog((f_chan_hz/1e+6),(ratio2*flux_cyg_jy),label="peak value from point source")
plt.loglog((f_chan_hz/1e+6),(ratio4*flux_cyg_jy),label="exp")
plt.loglog((f_chan_hz/1e+6),flux_cas_jy,label="Values from Baars et al paper")
plt.xlabel('channels')
plt.ylabel('Ratio of Cas with Cyg')
plt.legend()
plt.show()

