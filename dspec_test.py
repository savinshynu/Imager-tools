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
import matplotlib.animation as animation


def avgp_chan_eq(ra,dec,chan_info,eq_info,data,el):

    #data in the form of channels,pixels,pixels
    #chan_info has frequency values of each channels in kHz
     
    resp = np.zeros(data.shape[0]) # initilize the output result
    for chan in range(data.shape[0]): # Looking at each channel
        fwhm = 2.2*(74/(chan_info[chan]/1e+6)) # calculating fwhm and beam area 
        sig1 = fwhm/2.355
        sig2 = sig1/np.sin(el*np.pi/180.0)# change elevation wrt sources using
        area = np.pi*sig1*sig2
        rad_in = 1.5*fwhm # defining radius of the annulus
        rad_out = 2.5*fwhm
        
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
        #print chan
        xmin = np.where((ins_im != 0))[0].min() # finding the extremum of pixel values for anulus to fit a gaussian
        xmax = np.where((ins_im != 0))[0].max()
        ymin = np.where((ins_im != 0))[1].min()
        ymax = np.where((ins_im != 0))[1].max()
        params = fitgaussian(ins_im[xmin:xmax+1,ymin:ymax+1]) # height,x,y,sigmax,sigmay from gaussian fit
        """              
        if chan ==10: #or chan ==197:
           print  params[0]
           plt.imshow(ins_im,cmap='jet',origin = 'lower')
           plt.show()
           plt.imshow(ann_im,cmap='jet',origin = 'lower')
           plt.show()
        """
        #print ins_im.max(),params[0]
        ins = np.array(ins)  # converting list to array
        annul = np.array(annul)
        #resp[chan] = (params[0] - np.median(annul))*area
        resp[chan] = ins.sum() - (np.median(annul)*len(ins))

    return resp

cyg, cas, vir, tau = [0,1,2,3] #defining source values

par_jy = [[4.695,0.085,-0.178],[5.745,-0.770,0],[5.023,-0.856,0],[3.915,-0.299,0]] #a,b,c fit parameters for each source from baars etal paper

def flux_jy_baars(v,source): # From Baars et al. paper values

    a,b,c = par_jy[source]

    log_s = a + b*np.log10(v) + c*(np.log10(v)**2)

    return 10**log_s


filename2 = "/leo/savin/wide_lwasv/oims_wide/images/58890/58890_080000_30.137MHz_49.837MHz.oims"
nchan = 198
ngrid = 128
peak = range(0,75)#[552]
stokes = 0
lat = 34.3484       #cordinates for LWA-SV sation
lon = -106.8858
ra_cyg, dec_cyg = [299.8, 40.7] # ra dec values for cyg and cas
ra_cas, dec_cas = [350.8, 58.8]
ra_vir, dec_vir = [187.70, 12.39]
ra_tau, dec_tau = [83.6, 22.01]

dset = np.zeros((len(peak),nchan,4,ngrid,ngrid))
wfall_cas = np.zeros((len(peak),nchan))
wfall_cyg = np.zeros((len(peak),nchan))
wfall_vir = np.zeros((len(peak),nchan))
wfall_tau = np.zeros((len(peak),nchan))
wfall_cas_ann = np.zeros((len(peak),nchan))
wfall_cyg_ann = np.zeros((len(peak),nchan))

#tset = np.zeros((len(peak),4))
db = OrvilleImageDB(filename2,'r')

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
    
    f_chan_hz = np.arange(f1,f2,bw)
    im_tm = dset[j,:,stokes,:,:]
    im = np.mean(im_tm,axis =0)

    x_cas, y_cas, alt_cas = trackcoord(ra_cas, dec_cas, MJD, H, M, S, ngrid, psize,lat, lon, cent_az, cent_alt)
    xa_cas, ya_cas = maxi_unsub(im, x_cas, y_cas)
    x_cyg, y_cyg, alt_cyg = trackcoord(ra_cyg, dec_cyg, MJD, H, M, S, ngrid, psize,lat, lon, cent_az, cent_alt)
    xa_cyg,ya_cyg = maxi_unsub(im, x_cyg, y_cyg)
    #print xa_cyg, ya_cyg, alt_cas, alt_cyg
    #print H, M, S
    x_vir, y_vir, alt_vir = trackcoord(ra_vir, dec_vir, MJD, H, M, S, ngrid, psize,lat, lon, cent_az, cent_alt)
    xa_vir, ya_vir = maxi_unsub(im, x_vir, y_vir)
    x_tau, y_tau, alt_tau = trackcoord(ra_tau, dec_tau, MJD, H, M, S, ngrid, psize,lat, lon, cent_az, cent_alt)
    xa_tau, ya_tau = maxi_unsub(im, x_tau, y_tau)

    wfall_cas[j,:] = dset[j,:,stokes,xa_cas, ya_cas]
    wfall_cyg[j,:] = dset[j,:,stokes,xa_cyg, ya_cyg]
    wfall_vir[j,:] = dset[j,:,stokes,xa_vir, ya_vir]
    wfall_tau[j,:] = dset[j,:,stokes,xa_tau, ya_tau]
    """ 
    #Taking the stokes I image and converting the corresponding pixel values to (azimuth,altitude) and (ra,dec) values
    eq_info = np.zeros((im_tm.shape[1],im_tm.shape[2],2)) #ra,dec values

    #hrz_info = np.zeros((cal_dat.shape[1],cal_dat.shape[2],2)) # azimuth, altitude values
    # getcoord and pix2hrz taken from conversion_coord.py

    for l in range(im_tm.shape[1]):
        for k in range(im_tm.shape[2]):
            m_ra,m_dec = getcoord(np.array([l]),np.array([k]),ngrid, psize, np.array([MJD]), np.array([H]), np.array([M]), np.array([S]), lat, lon, cent_az, cent_alt)
            eq_info[l,k,0] = m_ra
            eq_info[l,k,1] = m_dec
            #m_az, m_alt = pix2hrz(np.array([l]), np.array([k]), ngrid, psize, cent_az, cent_alt)
            #hrz_info[l,k,0] = m_az
            #hrz_info[l,k,1] = m_alt
    
    resp1 = avgp_chan_eq(ra_cas,dec_cas,f_chan_hz,eq_info,im_tm,alt_cas) #Getting the response function for Cas
    resp2 = avgp_chan_eq(ra_cyg,dec_cyg,f_chan_hz,eq_info,im_tm,alt_cyg) # same for Cyg
    wfall_cas_ann[j,:] = resp1
    wfall_cyg_ann[j,:] = resp2
    

flux_cas_jy = flux_jy_baars(f_chan_hz/1e+6,cas)    
flux_cyg_jy = flux_jy_baars(f_chan_hz/1e+6,cyg)

"""
"""
im = np.mean(dset[:,:,stokes,:,:],axis =(0,1)) #averaging over channels to get a mean image to visualise
#im_chan = sub_dat[1,stokes,:,:]
plt.pcolormesh(np.transpose(im),cmap ='jet')
#plt.pcolormesh(np.transpose(im_chan),cmap ='jet')
plt.show()

fig = plt.figure()
for i in range(wfall_tau.shape[0]):
    plt.plot(range(wfall_tau.shape[1]),wfall_tau[i,:])
plt.title("response in successive integrations using annulus value")
plt.xlabel('Channels')
plt.ylabel('Response')
plt.show()


plt.pcolormesh(wfall_tau)
#plt.imshow(wfall_cas, origin = 'lower')
plt.colorbar()
plt.show()

lc = np.mean(wfall_tau,axis = 1)
plt.plot(range(lc.shape[0]),lc)
plt.show()
"""


#making movie 
fig = plt.figure()
ax1=fig.add_subplot(1,2,1)
ax2=fig.add_subplot(1,2,2)
listims = []

for i in range(wfall_vir.shape[0]):
    ax1.set_xlabel('Channels')
    ax1.set_ylabel('Response of Cygnus A on MJD 58890 UTC 17')
    im1, = ax1.plot(range(wfall_vir.shape[1]),wfall_vir[i,:],label = 'image %d' % (i))
    
    im_dat = np.transpose(np.mean(dset[i,:,stokes,:,:],axis = 0))
    #im_dat = np.transpose(dset[i,80,stokes,:,:])
    vmin = sc(im_dat,0)
    vmax = sc(im_dat,100)
    #im = plt.pcolormesh(im_dat,cmap ='jet')
    im2 = ax2.imshow(im_dat,cmap ='jet',vmin = vmin, vmax =vmax,origin ='lower')


    listims.append([im1, im2])
    
ani = animation.ArtistAnimation(fig, listims, interval=500, repeat_delay=1000)
ani.save('movie_vir_scint.mp4')
plt.show()


"""
#making movie 
fig = plt.figure()
ims = []
for i in range(dset.shape[0]):
    im_dat = np.transpose(np.mean(dset[i,:,stokes,:,:],axis = 0))
    #im_dat = np.transpose(dset[i,80,stokes,:,:])
    vmin = sc(im_dat,5)
    vmax = sc(im_dat,99)
    #im = plt.pcolormesh(im_dat,cmap ='jet')
    im = plt.imshow(im_dat,cmap ='jet',vmin = vmin, vmax =vmax,origin ='lower')
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=500, blit = True, repeat_delay=1000)
#ani.save('movie_test.mp4')
plt.show()

resp_cas_avg_ann = np.mean(wfall_cas_ann,axis = 0)
resp_cyg_avg_ann = np.mean(wfall_cyg_ann,axis = 0)
plt.plot(range(len(resp_cas_avg_ann)),resp_cas_avg_ann)
plt.plot(range(len(resp_cyg_avg_ann)),resp_cyg_avg_ann)
plt.show()

resp_cas_avg = np.mean(wfall_cas,axis = 0)
resp_cyg_avg = np.mean(wfall_cyg,axis = 0)


ratio1  = resp_cas_avg/resp_cyg_avg
ratio2  = resp_cas_avg_ann/resp_cyg_avg_ann
ratio3 = flux_cas_jy/flux_cyg_jy

plt.plot(range(len(ratio1)), ratio1, label = 'Peak value averaged over some integrations')
plt.plot(range(len(ratio2)), ratio2, label = 'Annulus averaged over some integrations')
plt.plot(range(len(ratio3)), ratio3, label = 'Expected from baars et al')
plt.ylabel("a.u.")
plt.xlabel("Channels")
plt.show()
"""

