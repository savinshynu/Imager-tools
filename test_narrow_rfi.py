import numpy as np
import matplotlib.pyplot as plt
from sliding_rfi_flagger_man import main
from lsl import astro
import glob
import sys
from mad import median_absolute_deviation as mad
import os
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
from OrvilleImageDB import OrvilleImageDB
from utils import mask_image,maxi_sub,maxi_unsub
import h5py
from lc_correct_man import snr_in_lc
from low_snr import filter_snr 
from twod_gaussian import fitgaussian
 
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

    stokes = int(sys.argv[7])

    #find the corresponding .hdf5 file for the time and .oims file fo image
    path1 = '/leo/savin/wide_lwasv/oims_wide/standard_sub/images_hdf5/'+str(mjd)+'/'+str(mjd)+'_'+str(format(hint,'02d'))+'*.hdf5'
    files1  = sorted(glob.glob(path1))

    for filename1 in files1:
        hf=h5py.File(filename1,'r')
        #dset1=hf.get('image')
        tm_ind = hf.get('time')
        header = hf.get('header')
        head  =  header[0,:]
        ints = head[0] #number of integrations
        ngrid  = int(head[1]) #size of images (x-axis)
        psize = head[2] #angular size of a pixel (at zenith)
        nchan = int(head[3])
        cent_az = head[9]
        cent_alt = head[10]
        t_comp = tm_ind[:,1]*3600.0 + tm_ind[:,2]*60.0 + tm_ind[:,3]
        try:
           #peak = np.where((tm_ind[:,0]==mjd)&(tm_ind[:,1]==hint)&(tm_ind[:,2]==mint)&(tm_ind[:,3] > sint-diff_tm) &(tm_ind[:,3] < sint+diff_tm))[0]
           peak = np.where((abs(t_comp-t_in) < diff_tm))[0]
        except IndexError:
           continue
        if len(peak) == 0:
           continue
        peak_ind = int(np.median(peak))
        #noise = np.where((tm_ind[:,0]==mjd)&(tm_ind[:,1]==hint)&((tm_ind[:,2] == mint-2)|(tm_ind[:,2] == mint+2)))[0]
        noise = np.array(range(max(0,peak_ind-15),max(0,peak_ind-5))+range(min(ints,peak_ind+5),min(ints,peak_ind+15))) 
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

        for index in range(noise.shape[0]):
            
            try:
                hdr, dat = db.__getitem__(noise[index])
            except RuntimeError:
                print 'bad database at %d integration' % (i)
                continue

            dset_avg[index,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))

        sub_dat = np.mean(dset[:,:,stokes,:,:],axis=(0)) - np.mean(dset_avg[:,:,stokes,:,:],axis=(0))
        dat_img = np.mean(sub_dat,axis=0)

        time = np.reshape(np.array([mjd,hint,mint,sint]),(4,1))

        lat = 34.3491 #Station Latitude (of LWA-SV)
        lon = -106.886  #Station Latitude (of LWA-SV)

        xp,yp,alt = trackcoord(ra,dec,time[0,:],time[1,:],time[2,:],time[3,:],ngrid,psize,lat,lon,cent_az,cent_alt)

        xai,yai = maxi_sub(dat_img,xp,yp)
        print xp,yp,alt,xai,yai
    
      
        plt.plot((np.arange(sub_dat.shape[0])),sub_dat[:,xai,yai],color='blue', marker='.', linestyle='solid')
        #plt.xlabel("Frequency (MHz)")
        #plt.ylabel("Flux_response")
        plt.show()
       

        gd = main(sub_dat[:,xai,yai])
        #print gd

        dat_img_flag = np.mean(sub_dat[gd,:,:],axis=0)
        xa,ya = maxi_unsub(dat_img_flag,xai,yai)
        
        params = fitgaussian(dat_img_flag[xa-15:xa+16,ya-15:ya+16]) # height,x,y,sigmax,sigmay from gaussian fit
        print "axial ratio is %f" % (params[3]/params[4])
        print params
        mean_flag = np.mean(dat_img_flag)
        mean_noflag = np.mean(dat_img)
        #print (dat_img[xa,ya]-med_noflag)/(mad(dat_img,axis=None)),(dat_img_flag[xa,ya]-med_flag)/(mad(dat_img_flag,axis=None)) 
        
        dat_im, dat_im_nohrz = mask_image(dat_img,time,psize,lat,lon,cent_az,cent_alt)
        dat_im_flag, dat_im_flag_nohrz = mask_image(dat_img_flag,time,psize,lat,lon,cent_az,cent_alt)


        fill_flag = np.where((dat_im_flag_nohrz != 0))
        fill = np.where((dat_im_nohrz != 0))
       
        plt.imshow(np.transpose(dat_img),cmap='jet',origin='lower')
        plt.show()   
   
        plt.imshow(np.transpose(dat_im_flag_nohrz),cmap='jet',origin= 'lower')
        plt.show()
        
        print (dat_im_nohrz[xa,ya]-mean_noflag)/np.std(dat_im_nohrz[fill]), (dat_im_flag_nohrz[xa,ya]-mean_flag)/np.std(dat_im_flag_nohrz[fill_flag]) 
        #print mad(dat_img_flag,axis=None), np.std(dat_im_flag_nohrz[fill_flag])
        
        #noise_flag = mad(dat_img_flag,axis=None)
        #im_snr = dat_img_flag[xa,ya]/noise_flag
        
        if sys.argv[9] == 'c':
           ## results from the light curve
           lc_sig_noflag,lc_sig,lc_lin_pol,lc_circ_pol = snr_in_lc(ra,dec,db,peak,stokes)
           print lc_sig_noflag,lc_sig,lc_lin_pol,lc_circ_pol
           #lc_sig,lc_lin_pol,lc_circ_pol = snr_in_lc(ra,dec,db,peak)
           #print lc_sig,lc_lin_pol,lc_circ_pol

       
        db.close()   
        hf.close()
               
flag_narrow_rfi(int(sys.argv[8]))


