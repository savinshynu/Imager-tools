"""
Reads in the hdf5 file in stokes I and conducts the transient search and 
list out the candidate events as text files
"""


import sys
import glob
import os
import time as tm
import h5py
import numpy as np
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord

if sys.version_info > (2,):
    xrange = range
   
    
def findtranspim(data,time,dur,sigma,sz,size,runav,psize,lat,lon,cent_az,cent_alt):

    """
    Function to find transient candidates from the LWA images
    Returns a list of tranisent coordindates and their corresponding time as a numpy array

    Parameters:
    data - 2D image data
    time - time array
    dur - number of images to subtract
    sigma - detection threshold
    sz - Image grid size
    size - source size to mask
    runav - averaging parameter across time for finding events
    psize - physical size of each image pixel, degrees
    lat, lon - latitude and longitude of stations, in degrees
    cent_az, cent_alt - coordinates of the center of the image
    """
    
    # Masking of bright sources in the image
    tmap = np.zeros(data.shape)

    A = np.ones((data.shape[0],data.shape[1])) # mask for the horizon
    for a in xrange(A.shape[0]):
        for b in xrange(A.shape[1]):
            if (float(a)-(A.shape[0]/2.0))**2 + (float(b)-(A.shape[0]/2.0))**2 >  ((180*np.cos(30.0*(np.pi/180.0)))/(np.pi*psize))**2: # making below 30 degrees (180* cos(theta)/(np.pi * psize)
            #if ((float(a)-(A.shape[0]/2.0))/(140.0/np.pi))**2 + ((float(b)-(A.shape[0]/2.0))/(140.0/np.pi))**2 > 1.0:
               A[a,b] = 0

    #print A.shape
    badra = np.array((83.633,350.85,299.868,187.706,139.52375,252.78375,69.28,49.94,187.3,123.4,211.2,150.5,176.3,147,24.42,261.2,277.4,140.3,247.2))
    baddec = np.array((22.0145,58.815,40.734,12.391,-12.0956,4.9925,29.6,41.51,2.049,48.21,51.96,28.79,31.26,7.421,33.16,-0.9939,48.75,45.67,39.55))


    for source in xrange(badra.shape[0]):
        
        x,y,Att = trackcoord(badra[source],baddec[source],time[0,:],time[1,:],time[2,:],time[3,:],A.shape[0],psize,lat,lon,cent_az,cent_alt)

        if source == 0:
            X = np.reshape(np.round(x),(-1,1))
            Y = np.reshape(np.round(y),(-1,1))
        else:
            X = np.concatenate((X,np.reshape(np.round(x),(-1,1))),1)
            Y = np.concatenate((Y,np.reshape(np.round(y),(-1,1))),1)
   
    
    #print X[:,1]
    X = np.delete(X,np.arange(0,dur),0)
    Y = np.delete(Y,np.arange(0,dur),0)
    #print X.shape
    c = np.concatenate((np.reshape(X,(X.shape[0],1,-1)),np.reshape(Y,(Y.shape[0],1,-1))),1)
    #print c[1,:,:]
    #print ken
    A = np.tile(np.reshape(A,(A.shape[0],A.shape[1],1)),(1,1,data.shape[2]))
    #print A.shape
    
    # Averaging the data based on the runav value
    if runav == 0:
        tmap = ((data[:,:,np.arange(dur,data.shape[2])] - data[:,:,np.arange(0,data.shape[2] - dur)]) > np.tile(np.reshape(sigma*np.std(np.reshape(data[:,:,np.arange(dur,data.shape[2])] - data[:,:,np.arange(0,data.shape[2] - dur)],(-1,data.shape[2]-dur)),0) - np.mean((data[:,:,np.arange(dur,data.shape[2])] - data[:,:,np.arange(0,data.shape[2] - dur)]),(0,1)),(1,1,-1)),(data.shape[0],data.shape[1],1)))

    else:

        for i in xrange(dur):

            if i == 0:
                datasub = data[:,:,np.arange(i,data.shape[2] - (dur-i))]
            else:
                datasub =  datasub + data[:,:,np.arange(i,data.shape[2] - (dur-i))]

    
        datasub = datasub/dur
        datasub2 = np.swapaxes(data[:,:,np.arange(dur,data.shape[2])] - datasub,0,1)
        s = np.std(np.reshape(datasub2,(-1,datasub2.shape[2])),0)
        sm = np.median(s)


        for i in xrange(c.shape[0]):
            rc = np.argwhere(datasub2[:,:,i] > (50*np.std(datasub2[:,:,i] + np.mean(datasub2[:,:,i]))))
            rowtest = rc[:,0]
            coltest = rc[:,1]
            test = np.zeros((rowtest.shape[0]))

            
            for j in xrange(rowtest.shape[0]):
 
                if np.any((rowtest[j] >= (c[i,1,:]-3)) & (rowtest[j] <= (c[i,1,:]+3)) & (coltest[j] >= (c[i,0,:]-3)) & (coltest[j] <= (c[i,0,:]+3))):
                    
                    test[j] = 1
                #    print test[j]

            if np.any(test):
               # print 'yes'
                tmap[:,:,i] = 0
            elif s[i] > 2*sm:
                
                tmap[:,:,i] = 0
            else:
                for j in xrange(c.shape[2]):

                    if c[i,0,j] > 1.0:
                        
                        if (c[i,1,j]-size) < 1.0 or (c[i,0,j]-size) < 1.0 or c[i,0,j] + size > A.shape[0] or (c[i,1,j] + size) > A.shape[0]:
                            size = size - 2
 

                        maskx = np.arange(int(np.round(c[i,1,j]-size)),int(np.round(c[i,1,j] + size + 1.0)))
                        masky = np.arange(int(np.round(c[i,0,j]-size)),int(np.round(c[i,0,j] + size + 1.0)))
                        
                        for x in maskx:
                            for y in masky:
                                datasub2[x,y,i] = 0.0

                D = datasub2[:,:,i]
                

                S = np.std(D[D!=0.0])
                M = np.mean(D[D!=0.0])
                D = D*A[:,:,i]

  
                tmap[:,:,i] = D > np.tile(np.reshape(sigma*S + M,(1,1)),(data.shape[0],data.shape[1]))

   
                if np.any((A[:,:,i] == 0) & (datasub2[:,:,i]> np.tile(np.reshape(2*sigma*S + M,(1,1,-1)),(data.shape[0],data.shape[1])))):
                    tmap[:,:,i] = 0

                if np.std(data[:,:,i])==0 and np.std(data[:,:,i+dur]- data[:,:,i]) ==0:
                    tmap[:,:,i] = 0
  
                
    tt = np.argwhere(np.sum(tmap,(0,1)) > 0.0)[:,0]

 #   print tt
    tmap = tmap[:,:,tt]
    tdata = data[:,:,tt + dur]
    tmdata = data[:,:,tt]

#    
    ttime = time[:,tt+dur]
    #print ttime.shape

    # Finding sources satifying the threshold critereia
    for k in xrange(tdata.shape[2]):
        
        rc = np.argwhere(tmap[:,:,k] > 0.0)
      
        row = rc[:,0]
        col = rc[:,1]
  

        CC = (tdata[:,:,k] - tmdata[:,:,k])*A[:,:,0]

        rc = np.argwhere(CC == np.max(CC))
        row2 = rc[:,0]
        col2 = rc[:,1]


        corr = np.round(np.sqrt((row.astype('float')+1)**2.0 + (col.astype('float')+1)**2.0)/6.0)
        ucorr = np.unique(np.round(np.sqrt((row.astype('float')+1)**2.0 + (col.astype('float')+1)**2.0)/6.0))


        signifk = np.zeros((ucorr.shape[0]))
        urow = np.zeros((ucorr.shape[0]))
        ucol = np.zeros((ucorr.shape[0]))


        for rd in xrange(ucorr.shape[0]):
            signifk[rd] = np.max(CC[row[corr == ucorr[rd]],col[corr == ucorr[rd]]])/np.std(CC)
            urow[rd] = np.mean(row[corr ==ucorr[rd]])
            ucol[rd] = np.mean(col[corr ==ucorr[rd]])


        corr2 = np.round(np.sqrt((row2.astype('float') + 1.0)**2.0 + (col2.astype('float') + 1.0)**2.0)/6.0)
        ucorr2 = np.unique(np.round(np.sqrt((row2.astype('float')+1.0)**2.0 + (col2.astype('float') + 1.0)**2.0)/6.0))

        signifk2 = np.zeros((ucorr2.shape[0]))
        urow2 = np.zeros((ucorr2.shape[0]))
        ucol2 = np.zeros((ucorr2.shape[0]))

        for rd in xrange(ucorr2.shape[0]):
            signifk2[rd] = np.max(CC[row2[corr2 == ucorr2[rd]],col2[corr2 == ucorr2[rd]]])/np.std(CC)
            urow2[rd] = np.mean(row2[corr2 ==ucorr2[rd]])
            ucol2[rd] = np.mean(col2[corr2 ==ucorr2[rd]])

        # Getting the actual sky coordinates from the pixel coordinates
        if np.std(row) + np.std(col) < 20:
 
            rak,deck = getcoord(ucol,urow,sz,psize,np.tile(ttime[0,k],ucol.shape),np.tile(ttime[1,k],ucol.shape),np.tile(ttime[2,k],ucol.shape),np.tile(ttime[3,k],ucol.shape),lat,lon,cent_az,cent_alt)

            if np.any(np.isnan(rak)):
                print (urow[np.isnan(rak)])
                print (ucol[np.isnan(rak)])

        
        else:
            rak,deck = getcoord(ucol2,urow2,sz,psize,np.tile(ttime[0,k],ucol2.shape),np.tile(ttime[1,k],ucol2.shape),np.tile(ttime[2,k],ucol2.shape),np.tile(ttime[3,k],ucol2.shape),lat,lon, cent_az, cent_alt)

            signifk = signifk2
            if np.any(np.isnan(rak)):
                print (urow2[np.isnan(rak)])
                print (ucol2[np.isnan(rak)])

        if k == 0:

            ra = rak
            dec = deck
            signif = signifk

            ttimek = np.tile(np.reshape(ttime[:,k],(4,1)),(1,np.max(rak.shape)))
            #ttimek = ttime[:,k]
            #print ttimek.shape

        else:
            
            ra = np.append(ra,rak)
            dec = np.append(dec,deck)
            signif = np.append(signif,signifk)

            ttimek = np.concatenate((ttimek,np.tile(np.reshape(ttime[:,k],(4,1)),(1,np.max(rak.shape)))),1)

    try:
        ttime = ttimek

        # Galactic mask applying different threshold
        b = (180.0/np.pi)*np.arcsin(np.sin(dec*(np.pi/180.0))*np.cos(62.6*(np.pi/180.0)) - np.cos(dec*(np.pi/180.0))*np.sin(ra*(np.pi/180.0) - (282.5*(np.pi/180.0)))*np.sin(62.6*(np.pi/180.0)))
  
        l = 33.0 +(180.0/np.pi)*np.arcsin((np.cos(dec*(np.pi/180.0))*np.sin(ra*(np.pi/180.0) - (282.5*(np.pi/180.0)))*np.cos(62.6*(np.pi/180.0)) + np.sin(dec*(np.pi/180.0))*np.sin(62.6*(np.pi/180.0)))/np.cos(b*(np.pi/180.0)))


        gm = np.logical_or(np.logical_and(np.logical_and(b > -10 ,b < 10), signif < 8),np.logical_and(np.logical_and(b > -7 ,b < 7), signif < 10))

        galmask = np.logical_not(gm)
   

        ra = ra[galmask]
        dec = dec[galmask]
        ttime = ttime[:,galmask]
        
        return (ra,dec,ttime)
    
    except UnboundLocalError:

        return (np.nan,np.nan,np.nan)


