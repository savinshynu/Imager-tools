"""
Filtering the slowly moving airplanes in the images
"""


import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import os
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord,pix2hrz
from OrvilleImageDB import OrvilleImageDB
from utils import mask_image,maxi_sub,maxi_unsub, ang_dist_pix
import h5py



def filter_planes(inf,mjd_day,stokes):
    
    #stokes = 0 #For stokes I
    list_ap = np.zeros((inf.shape[0],9))
    num_ap = 0
    list_good = np.zeros((inf.shape[0],9))
    num_good = 0

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
    
                diff_tm = 2
                try:
                   #peak = np.where((tm_ind[:,0]==mjd)&(tm_ind[:,1]==hint)&(tm_ind[:,2]==mint)&(tm_ind[:,3] > sint-diff_tm) &(tm_ind[:,3] < sint+diff_tm))[0]
                   peak = np.where((abs(t_comp-t_in) < diff_tm))[0]
                except IndexError:
                       continue
                if len(peak) == 0:
                   continue
                peak_ind = int(np.median(peak))
                #noise = np.where((tm_ind[:,0]==mjd)&(tm_ind[:,1]==hint)&((tm_ind[:,2] == mint-2)|(tm_ind[:,2] == mint+2)))[0]
                noise = np.array(range(max(0,peak_ind-24),max(0,peak_ind-12))+range(min(ints,peak_ind+12),min(ints,peak_ind+24))) 
    
                #dset = np.zeros((peak.shape[0],nchan,4,ngrid,ngrid))
                dset_avg = np.zeros((noise.shape[0],nchan,4,ngrid,ngrid))
     
                for index in range(noise.shape[0]):
            
                    try:
                       hdr, dat = db.__getitem__(noise[index])
                    except RuntimeError:
                       print 'bad database at %d integration' % (i)
                       continue

                    dset_avg[index,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
        
                pos = range(peak_ind+1,min(peak_ind+13,720),1)
                neg = range(peak_ind-1,max(peak_ind-13,0),-1)

                try:
                   hdr, dat = db.__getitem__(peak_ind)
                except RuntimeError:
                   print 'bad database at %d integration' % (i)
            

                d1 = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))

                sub_dat1 = d1[:,stokes,:,:] - np.mean(dset_avg[:,:,stokes,:,:],axis=(0))
                dat_img1 = np.mean(sub_dat1,axis=0)

                #plt.imshow(np.transpose(dat_img1),cmap='jet',origin='lower')
                #plt.show()

                time = np.reshape(np.array([mjd,hint,mint,sint]),(4,1))

                lat = 34.3491 #Station Latitude (of LWA-SV)
                lon = -106.886  #Station Latitude (of LWA-SV)

                xp,yp,alt = trackcoord(ra,dec,time[0,:],time[1,:],time[2,:],time[3,:],ngrid,psize,lat,lon,cent_az,cent_alt)

                xa,ya = maxi_sub(dat_img1,xp,yp)
                az1, alt1 = pix2hrz(xa, ya, ngrid, psize, cent_az, cent_alt) 
                #print xp,yp,xa,ya
                az_old = az1
                alt_old = alt1
                #xa_old = xa
                #ya_old =ya
                #print "forward"

                coord_pos = []
                coord_neg = []
                num_pos = 0
                num_neg = 0
                for index in pos:
                    
                    try:
                       hdr, dat = db.__getitem__(index)
                    except RuntimeError:
                       print 'bad database at %d integration' % (i)

                    d2 = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
                    sub_dat2 = d2[:,stokes,:,:] - np.mean(dset_avg[:,:,stokes,:,:],axis=(0))
                    dat_img2 = np.mean(sub_dat2,axis=0)
                    dat_im2, dat_im_nohrz2 = mask_image(dat_img2,time,psize,lat,lon,cent_az,cent_alt)
            
                    #plt.imshow(np.transpose(dat_im_nohrz2),cmap='jet',origin='lower')
                    #plt.show()

                    sig2 = np.std(dat_im_nohrz2[dat_im_nohrz2 != 0])
                    mean2 = np.mean(dat_im_nohrz2[dat_im_nohrz2 != 0])
            
                    #print sig2,mean2

                    x_list2, y_list2 = np.where(((dat_im_nohrz2 - mean2) > 5*sig2) & (dat_im_nohrz2 != 0))
            

                    az_list2, alt_list2 = pix2hrz(x_list2, y_list2,  ngrid, psize, cent_az, cent_alt)

                    az_main = np.ones(len(az_list2))*az_old
                    alt_main = np.ones(len(az_list2))*alt_old
                    diff = ang_dist_pix(az_list2,alt_list2,az_main,alt_main)
                    #diff = ((x_list2 - xa_old)**2 + (y_list2 - ya_old)**2)**0.5
                    mov = np.where((diff < 20))[0]
                    if len(mov) == 0:
                       break
                    az_mean = np.mean(az_list2[mov])
                    alt_mean = np.mean(alt_list2[mov])
                    x_mean = np.mean(x_list2[mov])
                    y_mean = np.mean(y_list2[mov])
                    #print x_mean,y_mean
            
                    if ang_dist_pix(az_mean, alt_mean, az_old, alt_old) < 1:
                       break

                    az_old = az_mean
                    alt_old  = alt_mean
                    #xa_old = x_mean
                    #ya_old = y_mean
                    num_pos += 1
                    coord_pos.append([x_mean,y_mean])
        
        
                az_old = az1
                alt_old = alt1
                #xa_old = xa
                #ya_old = ya
                #print "backward"

                for index in neg:
                    
                    try:
                       hdr, dat = db.__getitem__(index)
                    except RuntimeError:
                       print 'bad database at %d integration' % (i)

                    d3 = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
                    sub_dat3 = d3[:,stokes,:,:] - np.mean(dset_avg[:,:,stokes,:,:],axis=(0))
                    dat_img3 = np.mean(sub_dat3,axis=0)
                    dat_im3, dat_im_nohrz3 = mask_image(dat_img3,time,psize,lat,lon,cent_az,cent_alt)
            
                    #plt.imshow(np.transpose(dat_im_nohrz3),cmap='jet',origin='lower')
                    #plt.show()

                    sig3 = np.std(dat_im_nohrz3[dat_im_nohrz3 != 0])
                    mean3 = np.mean(dat_im_nohrz3[dat_im_nohrz3 != 0])
                    #print sig3, mean3
                    x_list3, y_list3 = np.where(((dat_im_nohrz3 - mean3) > 5*sig3) & (dat_im_nohrz3 != 0))
            

                    az_list3, alt_list3 = pix2hrz(x_list3, y_list3, ngrid, psize, cent_az, cent_alt)

                    az_main = np.ones(len(az_list3))*az_old
                    alt_main = np.ones(len(az_list3))*alt_old
                    diff = ang_dist_pix(az_list3,alt_list3,az_main,alt_main)
                    #diff = ((x_list3 - xa_old)**2 + (y_list3 - ya_old)**2)**0.5
                    mov = np.where((diff < 20))[0]
                    if len(mov) == 0:
                       break
                    az_mean = np.mean(az_list3[mov])
                    alt_mean = np.mean(alt_list3[mov])
                    x_mean = np.mean(x_list3[mov])
                    y_mean = np.mean(y_list3[mov])
                    #print x_mean,y_mean
                    if ang_dist_pix(az_mean, alt_mean, az_old, alt_old) < 1:
                       break
             
                    az_old = az_mean
                    alt_old  = alt_mean
                    #xa_old = x_mean
                    #ya_old = y_mean
                    num_neg += 1
                    coord_neg.append([x_mean,y_mean])


                if num_pos + num_neg < 3 :
                   list_good[num_good,:] = inf[ival,:]
                   num_good += 1
                else:
                   list_ap[num_ap,:] = inf[ival,:] 
                   num_ap += 1

            db.close()   
            hf.close()
       
    return list_good[:num_good,:] #,list_ap[:num_ap,:] 

if __name__ == "__main__":
   
   day = sys.argv[1] 
   stokes = sys.argv[2]
   if int(stokes) == 0:
      files_txt = glob.glob("/leo/savin/wide_lwasv/oims_wide/standard_sub/corr_transients/%s/*.txt" % (day))
   elif int(stokes) == 3:
      files_txt = glob.glob("/leo/savin/wide_lwasv/oims_wide/standard_sub/corr_transients_stokesv/%s/*.txt" % (day)) 
      
   for filename in files_txt:
       if os.path.getsize(filename) != 0:
          print filename
          inf = np.loadtxt(filename)
          if len(inf.shape) == 1:
             inf = np.reshape(inf,(1,9)) 
          file_ext = os.path.splitext(filename)[0] 

          first, second = filter_planes(inf,int(day),int(stokes))

          np.savetxt(file_ext+"_good.txt",first,'%10.2f')
          np.savetxt(file_ext+ "_bad.txt", second,'%10.2f')
