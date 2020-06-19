import numpy as np
import matplotlib.pyplot as plt
from sliding_rfi_flagger_man import main
from lsl import astro
import glob
import sys
from scipy.optimize import curve_fit
from mad import median_absolute_deviation as mad
import os
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
from OrvilleImageDB import OrvilleImageDB
from utils import mask_image,maxi_sub,maxi_unsub
import h5py
from lc_correct_man import snr_in_lc
from collect_resp_data import collect_resp

mjd = int(float(sys.argv[3]))
im_dat,t_dat,head_dat,el_dat = collect_resp(mjd)
cyg, cas, vir, tau = [0,1,2,3]
stokes = 0
source = cyg

el_low = 64.4
el_up = 64.6

name = '58831'

ylow = 0
yhigh = 2500

f1 = 30.1375007629
f2 = 25.1375007629

#ind =  np.where((el_dat[:,source]< el_up) & (el_dat[:,source] > el_low)& (np.int_(head_dat[:,source])== int(f1)))[0]
ind =  np.where((el_dat[:,source]< el_up) & (el_dat[:,source] > el_low))[0]

for x in ind:
    print x,head_dat[x,0],head_dat[x,1],el_dat[x,source],t_dat[x,1],t_dat[x,2]
indm = range(13999,14010)
xdat = indm
ydat = np.mean(im_dat[xdat,source,:,stokes],axis = 1)
plt.scatter(xdat,ydat)
plt.show()

spect_cal = np.mean(im_dat[xdat,source,:,stokes],axis =0)
plt.plot(range(im_dat.shape[2]),spect_cal)
plt.show()


par_jy = [[4.695,0.085,-0.178],[5.745,-0.770,0],[5.023,-0.856,0],[3.915,-0.299,0]]

def flux_jy_baars(v): # From Baars et al. paper values

    a,b,c = par_jy[source]

    log_s = a + b*np.log10(v) + c*(np.log10(v)**2)

    return 10**log_s



def flag_narrow_rfi(capt):
       
    ra = float(sys.argv[1])
    dec = float(sys.argv[2])
    mjd = float(sys.argv[3])
    h = float(sys.argv[4])
    m = float(sys.argv[5])
    s = float(sys.argv[6])
           
    mjd =int(mjd)    
    hint = int(h)
    mint = int(m)
   
    if (m-mint) == 0:
       sint = int(s)

    elif (m-mint) != 0:
       sec = s+(m-mint)*60.0
       if sec >= 60.0:
          sec = int(sec%60.0)
          mint += 1
          if mint >= 60:
             mint = mint%60
             hint += 1
       sint = int(sec)
  
    t_in = hint*3600.0 + mint*60.0 + sint
    
    if capt not in [5,15,60]:
       sys.exit("Wrong integration length")

    if capt == 5:
       diff_tm = 3
    elif capt == 15:
       diff_tm = 8
    elif capt == 60:
       diff_tm = 30

    stokes = 0

    #find the corresponding .hdf5 file for the time and .oims file fo image
    path1 = '/leo/savin/wide_lwasv/oims_wide/standard_sub/images_hdf5/'+str(mjd)+'/'+str(mjd)+'_'+str(format(hint,'02d'))+'*.hdf5'
    files1  = sorted(glob.glob(path1))

    for filename1 in files1:
        hf=h5py.File(filename1,'r')
        tm_ind = hf.get('time')
        header = hf.get('header')
        head  =  header[0,:]
        ints = head[0] #number of integrations
        ngrid  = int(head[1]) #size of images (x-axis)
        psize = head[2] #angular size of a pixel (at zenith)
        nchan = int(head[3])
        print nchan
        cent_az = head[9]
        cent_alt = head[10]
        t_comp = tm_ind[:,1]*3600.0 + tm_ind[:,2]*60.0 + tm_ind[:,3]
        try:
           peak = np.where((abs(t_comp-t_in) < diff_tm))[0]
        except IndexError:
           continue
        if len(peak) == 0:
           continue
        peak_ind = int(np.median(peak))
        noise = np.array(range(max(0,peak_ind-15),max(0,peak_ind-5))+range(min(ints,peak_ind+5),min(ints,peak_ind+15))) 
        #noise = np.array(range(598,608)+range(630,640)) 
        print peak
        print noise
        #print peak.shape[0],noise.shape[0],nchan, ngrid
    
        file_ext1 = os.path.basename(filename1)
        file_ext1 = os.path.splitext(file_ext1)[0]
        filename2 = '/leo/savin/wide_lwasv/oims_wide/images/'+str(mjd)+'/'+file_ext1+'.oims' 
        
        dset = np.zeros((peak.shape[0],nchan,4,ngrid,ngrid))
        dset_avg = np.zeros((noise.shape[0],nchan,4,ngrid,ngrid))
        db = OrvilleImageDB(filename2,'r')
     
        for index in range(peak.shape[0]):
            
            try:
                hdr, dat = db.__getitem__(peak[index])
            except RuntimeError:
                print 'bad database at %d integration' % (i)
                continue

            dset[index,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
         
        start_freq = hdr['start_freq']/1e+6
        bandwidth = hdr['bandwidth']/1e+6
        stop_freq = hdr['stop_freq']/1e+6 + bandwidth
        freq_chan  = np.arange(start_freq,stop_freq,bandwidth)
        print start_freq, stop_freq

        for index in range(noise.shape[0]):
            
            try:
                hdr, dat = db.__getitem__(noise[index])
            except RuntimeError:
                print 'bad database at %d integration' % (i)
                continue

            dset_avg[index,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))

        if sys.argv[8]=='u':
           sub_dat = np.mean(dset[:,:,stokes,:,:],axis=(0))
        elif sys.argv[8]=='s':
           sub_dat = np.mean(dset[:,:,stokes,:,:],axis=(0)) - np.mean(dset_avg[:,:,stokes,:,:],axis=(0))

        dat_img = np.mean(sub_dat,axis=0)
        #lat = 34.3491 #Station Latitude (of LWA-SV)
        #lon = -106.886  #Station Latitude (of LWA-SV)
        #time = np.reshape(np.array([mjd,hint,mint,sint]),(4,1))
        #dat_im_mask,dat_im_mask_nohrz = mask_image(dat_img,time,psize,lat,lon,cent_az,cent_alt) 

        plt.pcolormesh(np.transpose(dat_img),cmap='jet')
        plt.show()
        time = np.reshape(np.array([mjd,hint,mint,sint]),(4,1))

        lat = 34.3491 #Station Latitude (of LWA-SV)
        lon = -106.886  #Station Latitude (of LWA-SV)

        xp,yp,alt = trackcoord(ra,dec,time[0,:],time[1,:],time[2,:],time[3,:],ngrid,psize,lat,lon,cent_az,cent_alt)

        xa,ya = maxi_sub(dat_img,xp,yp)
        print xp,yp,xa,ya
    
      
        plt.plot((np.arange(sub_dat.shape[0])),sub_dat[:,xa,ya],color='blue', marker='.', linestyle='solid')
        #plt.xlabel("Frequency (MHz)")
        #plt.ylabel("Flux_response")
        plt.show()
        gd = main(sub_dat[:,xa,ya])
        db.close()   
        hf.close()
        return freq_chan[gd], gd, sub_dat[gd,xa,ya]
               
if sys.argv[9] == 'c':

   freq_chan, spect_ind,spect_dat = flag_narrow_rfi(int(sys.argv[7]))

   flux_cal_jy = flux_jy_baars(freq_chan)
   spect_cal = spect_cal[spect_ind]

   true_spect = (spect_dat/spect_cal)*flux_cal_jy
   plt.plot(range(len(freq_chan)),spect_dat/spect_dat.max(),label='cas')
   plt.plot(range(len(freq_chan)),spect_cal/spect_cal.max(),label = 'cyg')
   plt.xlabel('Channels')
   plt.ylabel('Power a.u.')
   plt.legend()
   plt.show()
   plt.scatter(freq_chan, true_spect,color='blue', marker='.', linestyle='solid')
   plt.ylabel("Flux density (Jy)")
   plt.xlabel("Frequency (MHz)")
   plt.show()


   inc = 3
   npt = int(np.ceil(true_spect.shape[0]/float(inc)))
   result = np.zeros((npt,2),dtype =float)

   p =0
   for m in range(0,true_spect.shape[0],inc):

       result[p,0]=np.mean(freq_chan[m:m+inc])
       result[p,1]=np.mean(true_spect[m:m+inc])
       p += 1



   #np.savetxt("txt_dat/spect_"+name+".txt",result)

   # For fitting the spectrum
   def func(x,a,b):

       return (a*(x**b))

   popt, pcov = curve_fit(func,result[:,0],result[:,1],maxfev=2000)
   perr = np.sqrt(np.diag(pcov))

   print popt
   print perr

   plt.scatter(result[:,0],result[:,1],s=20,color='blue', label= 'Data')
   plt.plot(result[:,0],func(result[:,0], *popt),label='Spectral index = %5.3f $\pm$ %5.3f'%(popt[1],perr[1]),color='r',linestyle = 'solid', linewidth=2)

   #plt.plot(result[:,0],func(result[:,0], *popt),label='Spectral index = %5.3f $\pm$ %5.3f'%(popt[1]+par_jy[source][1],perr[1]),color='r',linestyle = 'solid', linewidth=2)
   plt.xlabel('Frequency (MHz)')
   plt.ylabel('Flux Density (Jy)')
   #plt.ylim(ylow,yhigh)
   #plt.xlim(22,47)
   plt.legend()
   plt.show()

   print np.mean(result[:,0]),np.mean(result[:,1])

elif sys.argv[9] == "nc":
     sys.exit()
                               
