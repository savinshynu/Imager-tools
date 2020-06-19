#!/usr/bin/env python

import os
import sys
import ephem
import numpy

from astropy import wcs

from lsl.astro import MJD_OFFSET, DJD_OFFSET

from OrvilleImageDB import OrvilleImageDB

from matplotlib import pyplot as plt


def separation(ra1,dec1,ra2,dec2):
    x=numpy.sin(dec1*(numpy.pi/180.0))*numpy.sin(dec2*(numpy.pi/180.0))
    y=numpy.cos(dec1*(numpy.pi/180.0))*numpy.cos(dec2*(numpy.pi/180.0))*numpy.cos((ra1-ra2)*(numpy.pi/180.0))
    z=x+y
    adist=numpy.arccos(z)
    return adist*(180.0/numpy.pi)


def filter_plane(img, db, hdr):
    t = wcs.WCS(naxis=2)
    t.wcs.crpix = [img.shape[2]//2, img.shape[3]//2]
    t.wcs.cdelt = [db.header.pixel_size, db.header.pixel_size]
    t.wcs.crval = [hdr['center_az'], hdr['center_alt']]
    t.wcs.ctype = ['RA---SIN', 'DEC--SIN']
    x, y = numpy.arange(img.shape[2]), numpy.arange(img.shape[3])
    x, y = numpy.meshgrid(x, y)
    world = t.wcs_pix2world(numpy.array([x.ravel(),y.ravel()]).T, 0)
    az, alt = (world[:,0]+90)%360, world[:,1]
    az.shape = x.shape
    alt.shape = x.shape
    
    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [img.shape[2]//2, img.shape[3]//2]
    w.wcs.cdelt = [-db.header.pixel_size, db.header.pixel_size]
    w.wcs.crval = [hdr['center_ra'], hdr['center_dec']]
    w.wcs.ctype = ['RA---SIN', 'DEC--SIN']
    epoch = hdr['start_time'] + (MJD_OFFSET - DJD_OFFSET)
    
    f = numpy.arange(img.shape[0])*hdr['bandwidth'] + hdr['start_freq']
    f.shape = (f.size,1,1)
    alt.shape = (1,)+alt.shape
    fwhm = numpy.clip(1.0/numpy.sqrt(numpy.sin(alt*numpy.pi/180)), 1.0, 1.0) * 2.2*74e6/f # geometric mean of the FWHMs
    f = f.ravel()
    
    x, y = numpy.arange(img.shape[2]), numpy.arange(img.shape[3])
    x, y = numpy.meshgrid(x, y)
    world = w.wcs_pix2world(numpy.array([x.ravel(),y.ravel()]).T, 0)
    ra, dec = world[:,0], world[:,1]
    ra.shape = x.shape
    dec.shape = x.shape
    l, b = ra*0.0, dec*0.0
    for i in xrange(ra.shape[0]):
        for j in xrange(ra.shape[1]):
            if not numpy.isfinite(ra[i,j]):
                l[i,j] = b[i,j] = numpy.nan
                continue
            eq = ephem.Equatorial(ra[i,j]*numpy.pi/180, dec[i,j]*numpy.pi/180, epoch=epoch)
            ga = ephem.Galactic(eq)
            l[i,j] = ga.lon*180/numpy.pi
            b[i,j] = ga.lat*180/numpy.pi
            
    clean = img*0.0
    for i in xrange(l.shape[0]):
        for j in xrange(l.shape[1]):
            if not numpy.isfinite(l[i,j]):
                clean[:,:,i,j] = 0.0
                continue
            dist = separation(l,b, l[i,j],b[i,j])
            for k in xrange(img.shape[0]):
                inner  = 2.50
                outer  = 3.75
                height = 1.50
                confined_ring = numpy.where((dist >= inner*fwhm[k,:,:]) \
                                            & (dist < outer*fwhm[k,:,:]) \
                                            & (numpy.abs(b-b[i,j]) <= height*fwhm[k,:,:]))
                if len(confined_ring[0]) == 0:
                    clean[k,:,i,j] = 0.0
                    continue
                confined_ring = img[k,:,confined_ring[0],confined_ring[1]]
                clean[k,:,i,j] = numpy.median(confined_ring, axis=0)
    
    #for gps in numpy.linspace(-90, 90, 181):
    #    valid = numpy.where(numpy.abs(b-gps) < 4)
    #    if len(valid[0]) < 3:
    #        continue
    #    lng = l[valid]
    #    lat = b[valid]
    #    data = img[:,:,valid[0],valid[1]]
    #    
    #    lng = lng.ravel()
    #    lat = lat.ravel()
    #    data = data.reshape(img.shape[0],img.shape[1],-1)
    #    
    #    for i in xrange(lng.size):
    #        dist = separation(lng,lat, lng[i],lat[i])
    #        valid2 = numpy.where((dist >= 10) & (dist < 15))[0]
    #        if len(valid2) < 3:
    #            continue
    #        clean[:,:,valid[0][i],valid[1][i]] = numpy.median(data[:,:,valid2], axis=2)
        
    diff = img-clean
    return clean, diff


def main(args):
    db = OrvilleImageDB(args[0], 'r')
    hdr, img = db[600]
    print(hdr, img.shape)
    
    clean, diff= filter_plane(img, db, hdr)
    
    for i in (0, 100, -1):
        fig = plt.figure()
        ax = fig.add_subplot(1, 3, 1)
        ax.imshow(img[i,0,:,:], origin='lower')
        ax = fig.add_subplot(1, 3, 2)
        ax.imshow(clean[i,0,:,:], origin='lower')
        ax = fig.add_subplot(1, 3, 3)
        ax.imshow(diff[i,0,:,:], origin='lower')
        plt.show()
    
    
    
    


if __name__ == "__main__":
    main(sys.argv[1:])
    
