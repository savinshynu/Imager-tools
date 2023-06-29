"""
Flags narrowband RFI in the candidate list in stokes I

"""


import numpy as np
import matplotlib.pyplot as plt
from sliding_rfi_flagger_adv import main
from lsl import astro
import glob
import sys
from mad import median_absolute_deviation as mad
import os
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
from OrvilleImageDB import OrvilleImageDB
from utils import mask_image,maxi_sub,maxi_unsub
import h5py
from lc_correct import snr_in_lc


def flag_narrow_rfi(inf,capt,mjd_day):
       
    
    if capt not in [5,15,60]:
       sys.exit("Wrong integration length")

    if capt == 5:
       diff_tm = 3
    elif capt == 15:
       diff_tm = 8
    elif capt == 60:
       diff_tm = 30

    stokes = 0 #For stokes I
    #print inf.shape
    list_rfi1 = np.zeros((inf.shape[0],9))
    num_rfi1 =0
    #list_rfi2 = np.zeros((inf.shape[0],9))
    #num_rfi2 =0
    for hour in range(0,24):
        #find the corresponding .hdf5 file for the time and .oims file fo image
        path1 = '/leo/savin/wide_lwasv/oims_wide/standard_sub/images_hdf5/'+str(mjd_day)+'/'+str(mjd_day)+'_'+str(format(hour,'02d'))+'*.hdf5'
        files1  = sorted(glob.glob(path1))
        hour_index = np.where((inf[:,3] == hour))[0]
        
        if len(hour_index) == 0:
           print "No transients at utc "+str(hour)
           continue

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
           
            file_ext1 = os.path.basename(filename1)
            file_ext1 = os.path.splitext(file_ext1)[0]
            filename2 = '/leo/savin/wide_lwasv/oims_wide/images/'+str(mjd_day)+'/'+file_ext1+'.oims' 
            db = OrvilleImageDB(filename2,'r')
            print filename2                  
            for ival in hour_index:
                ra = float(inf[ival,0])
                dec = float(inf[ival,1])
                mjd = float(inf[ival,2])
                h = float(inf[ival,3])
                m = float(inf[ival,4])
                s = float(inf[ival,5])

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


                try:
                   peak = np.where((abs(t_comp-t_in) < diff_tm))[0]
                except IndexError:
                   print "Error: could not find the peak image"
                   continue

                if len(peak) == 0:
                   #print "requested time not in dataset"
                   print inf[ival,:6]
                   continue
 
                peak_ind = int(np.median(peak))
                noise = np.array(range(max(0,peak_ind-15),max(0,peak_ind-5))+range(min(ints,peak_ind+5),min(ints,peak_ind+15)))
               
                if len(noise) == 0:
                   #print "requested time not in dataset"
                   print inf[ival,:6]
                   continue              
                dset = np.zeros((peak.shape[0],nchan,4,ngrid,ngrid))
                dset_avg = np.zeros((noise.shape[0],nchan,4,ngrid,ngrid))
                
     
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
                gd = main(sub_dat[:,xai,yai])
                dat_img_flag = np.mean(sub_dat[gd,:,:],axis=0)
                
                xa,ya = maxi_unsub(dat_img_flag,xai,yai)

                im_noise_md = mad(dat_img_flag,axis=None)


                dat_im_flag, dat_im_flag_nohrz = mask_image(dat_img_flag,time,psize,lat,lon,cent_az,cent_alt) # calculating noise in the standard way
                fill_flag = np.where((dat_im_flag_nohrz != 0))
                im_noise_mn = np.std(dat_im_flag_nohrz[fill_flag])
                
                """
                if im_noise_md == 0:
                   im_noise = im_noise_mn
                else:
                   im_noise = im_noise_md
                """

                im_noise = im_noise_mn # changes made

                med_noise = np.median(dat_img_flag)
                mean_noise = np.mean(dat_img_flag)

                if im_noise == im_noise_md:
                   par_noise = med_noise
                else:
                   par_noise = mean_noise

                im_snr = (dat_img_flag[xa,ya]-par_noise)/im_noise

                ## results from the light curve
                lc_snr,lc_lin_pol,lc_circ_pol = snr_in_lc(ra,dec,db,peak,stokes)
                
                """
                if (im_snr > 5.0) and (lc_snr > 5.0):
                   if  (lc_lin_pol < 0.3) & (lc_circ_pol < 0.3):
                       list_rfi1[num_rfi1,:6] = inf[ival,:]
                       list_rfi1[num_rfi1,6] = lc_snr
                       list_rfi1[num_rfi1,7] = lc_lin_pol
                       list_rfi1[num_rfi1,8] = lc_circ_pol
                       num_rfi1 += 1
                   elif (lc_lin_pol > 0.3) | (lc_circ_pol > 0.3):
                        list_rfi2[num_rfi2,:6] = inf[ival,:]
                        list_rfi2[num_rfi2,6] = lc_snr
                        list_rfi2[num_rfi2,7] = lc_lin_pol
                        list_rfi2[num_rfi2,8] = lc_circ_pol
                        num_rfi2 += 1
                """
                if (im_snr > 5.0) and (lc_snr > 5.0):
                   
                   list_rfi1[num_rfi1,:6] = inf[ival,:]
                   list_rfi1[num_rfi1,6] = lc_snr
                   list_rfi1[num_rfi1,7] = lc_lin_pol
                   list_rfi1[num_rfi1,8] = lc_circ_pol
                   num_rfi1 += 1

            db.close()
            hf.close()

    return list_rfi1[:num_rfi1,:] #,list_rfi2[:num_rfi2,:]

          




