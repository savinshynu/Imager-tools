"""
A number of utility functions used in the transient search with Orville

"""


import numpy as np
from OrvilleImageDB import OrvilleImageDB
import sys
import matplotlib.pyplot as plt
from lsl import astro
from conversion_coord import pix2eq as getcoord, eq2pix as trackcoord
import math
from twod_gaussian import fitgaussian


# fo reading images from a .oims file continuously
def read_oims_cont(filename):

    db = OrvilleImageDB(filename,'r') #define the image data base

    #getting parameters from the input file

    ints = int(db.nint) #number of integration
    station =  db.header.station
    #stokes = db.header.stokes_params
    inp_flag = db.header.flags
    file_start = db.header.start_time
    file_end = db.header.stop_time
    ngrid = int( db.header.ngrid) #size of images (x-axis)
    psize = db.header.pixel_size #angular size of a pixel (at zenith)
    nchan = int(db.header.nchan) # number of channels
    #print ints,ngrid,psize,nchan #,station,stokes,inp_flag,file_start,file_end
    time = np.zeros((ints,4))
    data = np.zeros((ints,nchan,4,ngrid,ngrid))
    header = np.zeros((ints,9))
    for i in xrange(ints):
        hdr, dat = db.read_image()
        t = hdr['start_time']
        int_len = hdr['int_len']
        lst = hdr['lst']
        start_freq = hdr['start_freq']/1e+6
        stop_freq = hdr['stop_freq']/1e+6
        bandwidth = hdr['bandwidth']/1e+6
        cent_ra = hdr['center_ra']
        cent_dec = hdr['center_dec']
        cent_az = hdr['center_az']
        cent_alt = hdr['center_alt']

        #print int_len,lst,start_freq,stop_freq,bandwidth,cent_az,cent_alt
        #print cent_ra,cent_dec,cent_az,cent_alt

        MJD = int(t)
        h = 24.0*(t - float(MJD))
        H = int(h)
        m = 60*(h - float(H))
        M = int(m)
        s = 60*(m - float(M))
        S = int(s)


        time[i,:] = np.reshape(np.array([MJD,H,M,S]),(4,))
        data[i,:,:,:,:] = np.transpose(dat[:,:,:,:],axes=(0,1,3,2))
        header[i,:] = [ints, ngrid, psize, nchan, start_freq, stop_freq, bandwidth, cent_az, cent_alt]
    return time,data,header


def mask_image(im,time,psize,lat,lon,cent_az,cent_alt):

    dat_im = im
    # mask for the horizon

    A = np.ones((dat_im.shape[0],dat_im.shape[1]))
    for a in xrange(A.shape[0]):
        for b in xrange(A.shape[1]):
            if (float(a)-(A.shape[0]/2.0))**2 + (float(b)-(A.shape[0]/2.0))**2 >  ((180*np.cos(25*(np.pi/180.0)))/(np.pi*psize))**2: #(180* cos(theta)/(np.pi * psize)
                A[a,b] = 0

    sp =  astro.get_solar_equ_coords(time[0,0] + (time[1,0] + (time[2,0] + time[3,0]/60.0)/60.0)/24.0  + 2400000.5)
    Sunra = sp[0]
    Sundec = sp[1]

    badra = np.array((83.633,350.85,299.868,187.706,139.52375,252.78375,Sunra,69.28,49.94,187.3,123.4,211.2,150.5,176.3,147,24.42,261.2,277.4,140.3,247.2))
    baddec = np.array((22.0145,58.815,40.734,12.391,-12.0956,4.9925,Sundec,29.6,41.51,2.049,48.21,51.96,28.79,31.26,7.421,33.16,-0.9939,48.75,45.67,39.55))

    for source in xrange(badra.shape[0]):

        x,y,Att = trackcoord(badra[source],baddec[source],time[0,:],time[1,:],time[2,:],time[3,:],A.shape[0],psize,lat,lon,cent_az,cent_alt)
        #x = int(np.round(x))
        #y = int(np.round(y))
        if Att > 0:
           xa,ya = maxi_unsub(dat_im,x,y) 
           dat_im[xa-6:xa+7,ya-6:ya+7] = 0

    dat_im_nohrz = dat_im*A
    return dat_im,dat_im_nohrz 

#finding peak pixel in a  subtracted image
def maxi_sub(data,x_in,y_in):


    x_in=int(x_in)
    y_in=int(y_in)
    maxi=data[x_in,y_in]
    x1=x_in
    y1=y_in
    for l in range(x_in-2,x_in+3,1):
        for k in range(y_in-2,y_in+3,1):
            if data[l,k] > maxi:
               maxi=data[l,k]
               x1=l
               y1=k

    return x1,y1

#finding peak pixel in a unsubtracted image

def maxi_unsub(data,x_in,y_in):


    x_in=int(x_in)
    y_in=int(y_in)
    maxi=data[x_in,y_in]
    x1=x_in
    y1=y_in
    for l in range(x_in-1,x_in+2,1):
        for k in range(y_in-1,y_in+2,1):
            if data[l,k] > maxi:
               maxi=data[l,k]
               x1=l
               y1=k

    return x1,y1

# time from mjd day
def mjd2ut(t):

    MJD = int(t)
    h = 24.0*(t - float(MJD))
    H = int(h)
    m = 60*(h - float(H))
    M = int(m)
    s = 60*(m - float(M))
    S = int(s)

    return MJD, H, M, S

#for finding an annulus around a point in sky and subtracting the sky contribution for a single channel data

def avgp_test(xpix,ypix,data):

    xpix=int(xpix)
    ypix=int(ypix)
    sum1 = np.sum(data[xpix-4:xpix+5,ypix-4:ypix+5])
    sum2 = np.sum(data[xpix-5:xpix+6,ypix-5:ypix+6])
    resid = sum2 -sum1
    avg2 = resid/40.0
    return  sum1-(81.0*avg2)

    
    
def avgp_med(xpix,ypix,data):

    xpix=int(xpix)
    ypix=int(ypix)
    ins = data[xpix-3:xpix+4,ypix-3:ypix+4]
    ann = list(data[xpix-4:xpix+5,ypix-4]) + list(data[xpix-4:xpix+5,ypix+4]) + list(data[xpix-4,ypix-3:ypix+4]) + list(data[xpix+4,ypix-3:ypix+4])
    ann = np.array(ann)
    sum1 = np.sum(ins)
    med = np.median(ann)
    ins_len = len(np.ravel(ins))
    return  sum1-(med*ins_len)




def avgp(xpix,ypix,data):

    xpix=int(xpix)
    ypix=int(ypix)
    sum1 = np.sum(data[xpix-3:xpix+4,ypix-3:ypix+4])
    sum2 = np.sum(data[xpix-4:xpix+5,ypix-4:ypix+5])  
    resid = sum2 -sum1
    avg2 = resid/32.0
    return  sum1-(49.0*avg2)

"""


#for finding an annulus around a point in sky and subtracting the sky contribution for multi channel data
def avgp_chan(xpix,ypix,data):
    
    resp = np.zeros(data.shape[0])
    xSize,ySize = [data.shape[1],data.shape[2]]
    xpix=int(xpix)
    ypix=int(ypix)
    for chan in range(data.shape[0]):
  
        sum1 = np.sum(data[chan,xpix-3:xpix+4,ypix-3:ypix+4])
        sum2 = np.sum(data[chan,xpix-4:xpix+5,ypix-4:ypix+5])
        resid = sum2 -sum1
        avg2 = resid/32.0
        resp[chan] =  sum1-(49.0*avg2)

    return resp

#for finding an annulus around a point in sky and subtracting the sky contribution for multi channel data
def avgp_chan(xpix,ypix,data,chan_info):
    
    resp = np.zeros(data.shape[0])
    xpix=int(xpix)
    ypix=int(ypix)
    for chan in range(data.shape[0]):
        fwhm = 2.2*(74/chan_info[chan])

        rad_in = 1.5*fwhm
        xin = int(round(rad_in/1.0156))
        yin = int(round(rad_in/1.0156))
        xout = xin+1
        yout = yin+1
        sum1 = np.sum(data[chan,xpix-xin:xpix+xin+1,ypix-yin:ypix+yin+1])
        sum2 = np.sum(data[chan,xpix-xout:xpix+xout+1,ypix-yout:ypix+yout+1])
        resid = sum2 -sum1
        resid_len = ((2*xout+1)*(2*yout+1))-((2*xin+1)*(2*yin+1))
        avg2 = resid/resid_len
        in_len = (2*xin+1)*(2*yin+1)
        resp[chan] =  sum1-(in_len*avg2)

    return resp
"""

def avgp_chan_eq(ra,dec,chan_info,annul_info,data):

    size, psize, mjd, h, m, s, lat, lon, cent_az, cent_alt = annul_info
    eq_info = np.zeros((data.shape[1],data.shape[2],2))
    for l in range(data.shape[1]):
        for k in range(data.shape[2]):
            m_ra,m_dec = getcoord(np.array([l]),np.array([k]),size, psize, np.array([mjd]), np.array([h]), np.array([m]), np.array([s]), lat, lon, cent_az, cent_alt)
            eq_info[l,k,0] = m_ra
            eq_info[l,k,1] = m_dec
    resp = np.zeros(data.shape[0])
    for chan in range(data.shape[0]):
        fwhm = 2.2*(74/chan_info[chan])
        sig1 = fwhm/2.355
        sig2 = sig1/np.sin(63.5*np.pi/180.0)
        area = np.pi*sig1*sig2
        rad_in = 1.5*fwhm
        rad_out = 2.5*fwhm
        #print area,rad_in,rad_out
        ins = []
        annul =[]
        ins_im = np.zeros((data.shape[1],data.shape[2]))
        ann_im = np.zeros((data.shape[1],data.shape[2]))
        for x in range(eq_info.shape[0]):
            for y in range(eq_info.shape[1]):
                if ((eq_info[x,y,0]-ra)**2+(eq_info[x,y,1]-dec)**2 < rad_in**2):
                   ins.append(data[chan,x,y]) 
                   ins_im[x,y] = data[chan,x,y]
                elif (((eq_info[x,y,0]-ra)**2+(eq_info[x,y,1]-dec)**2 < rad_out**2) and ((eq_info[x,y,0]-ra)**2+(eq_info[x,y,1]-dec)**2 > rad_in**2)):
                   annul.append(data[chan,x,y])
                   ann_im[x,y] = data[chan,x,y]

        xmin = np.where((ins_im != 0))[0].min()
        xmax = np.where((ins_im != 0))[0].max()
        ymin = np.where((ins_im != 0))[1].min()
        ymax = np.where((ins_im != 0))[1].max()
        params = fitgaussian(ins_im[xmin:xmax+1,ymin:ymax+1])
        #plt.pcolormesh(ins_im[xmin:xmax+1,ymin:ymax+1],cmap='jet')
        #plt.show()
        #plt.pcolormesh(ann_im,cmap='jet')
        #plt.show()
        #print np.median(annul)
        #plt.hist(data[chan,40:80,70:110])
        #plt.show()
        
        ins = np.array(ins)
        annul = np.array(annul)
                
        #print len(ins),len(annul)

        resp[chan] = ins.sum() - (np.median(data[chan,:,:])*len(ins))
        #resp[chan] = params[0]*area -(np.median(annul)*area)
       
    return resp


#function to determine the angular distance between two points (ra1,dec1) and (ra2,dec2) given all of them are in degrees
def ang_dist(ra1,dec1,ra2,dec2):
    x = np.sin(dec1*(np.pi/180.0))*np.sin(dec2*(np.pi/180.0))
    y = np.cos(dec1*(np.pi/180.0))*np.cos(dec2*(np.pi/180.0))*np.cos((ra1-ra2)*(np.pi/180.0))
    z = x + y
    ang_dist = np.arccos(z)
    return ang_dist*(180.0/np.pi)


def ang_dist_pix(x1,y1,x2,y2):
    a = np.sin(y1*(np.pi/180.0))*np.sin(y2*(np.pi/180.0))
    b = np.cos(y1*(np.pi/180.0))*np.cos(y2*(np.pi/180.0))*np.cos((x1-x2)*(np.pi/180.0))
    z=a+b
    ang_dist = np.arccos(z)
    return ang_dist*(180.0/np.pi)

