"""
Different coordinate conversion scripts required for the 
transient searches with images

Author: S. S. Varghese
"""
import sys
import numpy as np

if sys.version_info > (2,):
    xrange = range


lat1= 34.0689         # coordinates for LWA station 1
lon1= -107.6284

lat2 = 34.3484       #cordinates for LWA-SV sation
lon2 = -106.8858


# Using the sin projection where you can adjust the phase center. This one converets ra,dec to azimuth,altitude which then maps to pixel using sin projection. The code needs the phase centers in azimuth and altitude (available in header info) which is fixed for the station.  

def eq2pix(ra, dec, MJD, h, m, s, size, psize, lat, lon, cent_az, cent_alt):
   """
   Conversion from equatorial coordinates to pixel coordinates
   Return pixel coordinates

   parameters:
   ra - right ascension, in degrees
   dec - declination in degrees
   MJD,H,M,S -  MJD day, hour, minute, seconds, float  
   size - grid size, float,  usually 128
   psize - pixel size,  in degrees 
   cent_az, cent_alt - phase center values in degrees 
   lat, lon - station coordinates, in degrees
   If passing an array of ra,dec also pass an array of MJD,H,M,S
   """
   RA = ra
   LAT = lat*(np.pi/180.0)
   UT = h + (m + s/60.0)/60.0
   d = MJD - 40000 - 11544 + UT/24.0
   DEC = dec*(np.pi/180)
   LST = np.mod(100.46 + 0.985647*d + lon + 15.0*UT,360) # calculation of LST and hour angle
   HA = np.mod(LST-RA,360)*(np.pi/180.0) 

   ALT = np.arcsin(np.sin(DEC)*np.sin(LAT) + np.cos(DEC)*np.cos(LAT)*np.cos(HA))   #calculation of topocentric coordinates
   az = np.arccos((np.sin(DEC) - np.sin(ALT)*np.sin(LAT))/(np.cos(ALT)*np.cos(LAT)))
   
   try:
      comp1 = np.where((np.sin(HA) >= 0.0))[0]
      comp2 = np.where((np.sin(HA) < 0.0))[0]
      AZ = HA*0.0
      AZ[comp1] = 2*np.pi - az[comp1]
      AZ[comp2] = az[comp2]
   
   except IndexError:
         if np.sin(HA) >= 0.0:
            AZ = 2*np.pi-az
         elif np.sin(HA) < 0.0:
            AZ=az

   cent_az *= np.pi/180 # Azimuth and Altitude phase center from the header information
   cent_alt *= np.pi/180
   
   # conversion to pixel coordinates  
   x = np.cos(ALT)*np.sin(AZ-(cent_az + np.pi/2.0))  # np.pi/2.0 is added to azimuth center to fix the orientation of  projection coordinates wrt to topocentric system
   y = (np.cos(cent_alt)*np.sin(ALT))-(np.sin(cent_alt)*np.cos(ALT)*np.cos(AZ-(cent_az + np.pi/2.0))) 

   x = (x*((180.0/np.pi)/psize)) + (size/2.0) -0.5  # 63,63 center pixel is used 
   y = (y*((180.0/np.pi)/psize)) + (size/2.0) -0.5
   
   Att = (ALT*180.0)/np.pi
   try:
      comp3 = np.where((Att < 0))[0]
      #print comp3
      x[comp3] = 1.0
      y[comp3] = 1.0

   except TypeError:
      if Att < 0.0:
         x = 1.0
         y = 1.0
   return (x,y,Att)


def eq2hrz(ra, dec, MJD, h, m, s, size, psize, lat, lon, cent_az, cent_alt):
   """
   Conversion from equatorial coordinates to horizontal coordinates
   Return horizontal coordinates, altitude and azimuth

   parameters:
   ra - right ascension, in degrees
   dec - declination in degrees
   MJD,H,M,S -  MJD day, hour, minute, seconds, float 
   size - grid size, float,  usually 128
   psize - pixel size,  in degrees 
   cent_az, cent_alt - phase center values in degrees 
   lat, lon - station coordinates, in degrees
   If passing an array of ra,dec also pass an array of MJD,H,M,S
   """
   RA = ra
   LAT = lat*(np.pi/180.0)
   UT = h + (m + s/60.0)/60.0
   d = MJD - 40000 - 11544 + UT/24.0
   DEC = dec*(np.pi/180)
   LST = np.mod(100.46 + 0.985647*d + lon + 15.0*UT,360) # calculation of LST and hour angle
   HA = np.mod(LST-RA,360)*(np.pi/180.0)

   ALT = np.arcsin(np.sin(DEC)*np.sin(LAT) + np.cos(DEC)*np.cos(LAT)*np.cos(HA))   #calculation of topocentric coordinates
   az = np.arccos((np.sin(DEC) - np.sin(ALT)*np.sin(LAT))/(np.cos(ALT)*np.cos(LAT)))
      
   try:
      comp1 = np.where((np.sin(HA) >= 0.0))[0]
      comp2 = np.where((np.sin(HA) < 0.0))[0]
      AZ = HA*0.0
      AZ[comp1] = 2*np.pi - az[comp1]
      AZ[comp2] = az[comp2]

   except IndexError:
         if np.sin(HA) >= 0.0:
            AZ = 2*np.pi-az
         elif np.sin(HA) < 0.0:
            AZ=az
   Att = (ALT*180.0)/np.pi
   AZ *= 180.0/np.pi 
   
   return (AZ,Att)



def hrz2pix(AZ, ALT, size, psize, cent_az, cent_alt):

   """
   Calculates pixel coordinates from horizontal coordinates
   Returns pixel coordinates

   Parameters:
   size - grid size, float,  usually 128
   psize - pixel size,  in degrees 
   cent_az, cent_alt - phase center values in degrees 
   lat, lon - station coordinates, in degrees
   """
 
   AZ = AZ*(np.pi/180.0) #Azimuth and Altitude phase center from the header information
   ALT = ALT*(np.pi/180.0)

   x = np.cos(ALT)*np.sin(AZ-(cent_az+np.pi/2.0))
   y = (np.cos(cent_alt)*np.sin(ALT))-(np.sin(cent_alt)*np.cos(ALT)*np.cos(AZ-(cent_az+np.pi/2.0)))

   
   x = (x*((180.0/np.pi)/psize)) + (size/2.0) -0.5
   y = (y*((180.0/np.pi)/psize)) + (size/2.0) -0.5

   Att = (ALT*180.0)/np.pi
   try:
      comp3 = np.where((Att < 0))[0]
      x[comp3] = 1.0
      y[comp3] = 1.0

   except TypeError:
      if Att < 0.0:
         x = 1.0
         y = 1.0

   return (x,y,Att)


def pix2eq(x, y, size, psize, MJD, h, m, s, lat, lon, cent_az, cent_alt):

   """
   Conversion from pixel coordinates to equatorial coordinates
   Return equatorial coordinates, ra, dec in degrees

   parameters:
   x, y  - pixel x, y coordinates, int
   MJD,H,M,S -  MJD day, hour, minute, seconds, float  
   size - grid size, float,  usually 128
   psize - pixel size,  in degrees 
   cent_az, cent_alt - phase center values in degrees 
   lat, lon - station coordinates, in degrees
   If passing an array of x, y also pass an array of MJD,H,M,S
   """
  
   LAT = lat*(np.pi/180.0)
   UT = h + (m + s/60.0)/60.0
   d = MJD - 40000 - 11544 + UT/24.0
   LST = np.mod(100.46 + 0.985647*d + lon + 15.0*UT,360)

   sz = (size/2.0)-0.5 
   x = ((x) - float(sz))*float(np.pi/180.0*psize)
   y = ((y) - float(sz))*float(np.pi/180.0*psize)

   cent_az *= np.pi/180
   cent_alt *= np.pi/180   # phase center values from the header info


   # The added changes for full sine projection
   rho = np.sqrt(x**2+y**2)
   c = np.arcsin(rho)

   try:
      comp = np.where( rho > 1e-10 )
      ALT = c*0.0 + cent_alt
      ALT[comp] = np.arcsin(np.cos(c)*np.sin(cent_alt)+((y*np.sin(c)*np.cos(cent_alt))/rho))[comp]

   except IndexError:
      if rho > 1e-10:
         ALT = np.arcsin(np.cos(c)*np.sin(cent_alt)+((y*np.sin(c)*np.cos(cent_alt))/rho)) 
      else:
         ALT = cent_alt


   az = cent_az + np.pi/2.0 + np.arctan2(x*np.sin(c),(rho*np.cos(c)*np.cos(cent_alt)-y*np.sin(c)*np.sin(cent_alt)))   
   az = np.mod(az,2*np.pi)

   DEC = np.arcsin(np.cos(az)*(np.cos(ALT)*np.cos(LAT)) + (np.sin(ALT)*np.sin(LAT)))
   HA = np.real(np.arccos((np.sin(ALT) - (np.sin(DEC)*np.sin(LAT)))/(np.cos(DEC)*np.cos(LAT))))
   RA = np.zeros(x.shape)

   try:
      for i in xrange(x.shape[0]):
         if x[i] < 0.0:
            RA[i] = np.mod(LST[i] + HA[i]/(np.pi/180.0),360.0)
         else:
            RA[i] = np.mod(LST[i] - HA[i]/(np.pi/180.0),360.0)
   except IndexError:
         if x < 0.0:
            RA = np.mod(LST + HA/(np.pi/180.0),360.0)
         else :
            RA = np.mod(LST - HA/(np.pi/180.0),360.0)

   dec = DEC*(180.0/np.pi)
   ra = RA

   return (ra,dec)



def pix2hrz(x, y, size, psize, cent_az, cent_alt):

   """
   Conversion from pixel coordinates to horizontal coordinates
   Return horizontal coordinates, altitude, azimuth in degrees

   parameters:
   x, y  - pixel x, y coordinates, int
   size - grid size, float,  usually 128
   psize - pixel size,  in degrees 
   cent_az, cent_alt - phase center values in degrees 
   
   """

   sz = (size/2.0)-0.5
   x = ((x) - float(sz))*float(np.pi/180.0*psize)
   y = ((y) - float(sz))*float(np.pi/180.0*psize)

   cent_az *= np.pi/180
   cent_alt *= np.pi/180   # phase center values from the header info


   # The added changes for full sine projection
   rho = np.sqrt(x**2+y**2)
   c = np.arcsin(rho)

   try:
      comp = np.where( rho > 1e-10 )
      ALT = c*0.0 + cent_alt
      ALT[comp] = np.arcsin(np.cos(c)*np.sin(cent_alt)+((y*np.sin(c)*np.cos(cent_alt))/rho))[comp]

   except IndexError:
      if rho > 1e-10:
         ALT = np.arcsin(np.cos(c)*np.sin(cent_alt)+((y*np.sin(c)*np.cos(cent_alt))/rho))
      else:
         ALT = cent_alt


   az = cent_az + np.pi/2.0 + np.arctan2(x*np.sin(c),(rho*np.cos(c)*np.cos(cent_alt)-y*np.sin(c)*np.sin(cent_alt)))
   az = np.mod(az,2*np.pi)

   return (az*180.0/np.pi,ALT*180.0/np.pi)




def eq_and_hrz(x, y, size, psize, MJD, h, m, s, lat, lon, cent_az, cent_alt):

   """
   Conversion from pixel coordinates to equatorial/horizontal coordinates
   Return equatorial and horizontal coordinates in degrees

   parameters:
   x, y  - pixel x, y coordinates, int
   MJD,H,M,S -  MJD day, hour, minute, seconds, float  
   size - grid size, float,  usually 128
   psize - pixel size,  in degrees 
   cent_az, cent_alt - phase center values in degrees 
   lat, lon - station coordinates, in degrees
   If passing an array of x, y also pass an array of MJD,H,M,S
   """
  
   LAT = lat*(np.pi/180.0)
   UT = h + (m + s/60.0)/60.0
   d = MJD - 40000 - 11544 + UT/24.0
   LST = np.mod(100.46 + 0.985647*d + lon + 15.0*UT,360)

   sz = (size/2.0)-0.5
   x = ((x) - float(sz))*float(np.pi/180.0*psize)
   y = ((y) - float(sz))*float(np.pi/180.0*psize)

   cent_az *= np.pi/180
   cent_alt *= np.pi/180   # phase center values from the header info


   # The added changes for full sine projection
   rho = np.sqrt(x**2+y**2)
   c = np.arcsin(rho)

   try:
      comp = np.where( rho > 1e-10 )
      ALT = c*0.0 + cent_alt
      ALT[comp] = np.arcsin(np.cos(c)*np.sin(cent_alt)+((y*np.sin(c)*np.cos(cent_alt))/rho))[comp]

   except IndexError:
      if rho > 1e-10:
         ALT = np.arcsin(np.cos(c)*np.sin(cent_alt)+((y*np.sin(c)*np.cos(cent_alt))/rho))
      else:
         ALT = cent_alt


   az = cent_az + np.pi/2.0 + np.arctan2(x*np.sin(c),(rho*np.cos(c)*np.cos(cent_alt)-y*np.sin(c)*np.sin(cent_alt)))
   az = np.mod(az,2*np.pi)

   DEC = np.arcsin(np.cos(az)*(np.cos(ALT)*np.cos(LAT)) + (np.sin(ALT)*np.sin(LAT)))
   HA = np.real(np.arccos((np.sin(ALT) - (np.sin(DEC)*np.sin(LAT)))/(np.cos(DEC)*np.cos(LAT))))
   RA = np.zeros(x.shape)

   try:
      for i in xrange(x.shape[0]):
         if x[i] < 0.0:
            RA[i] = np.mod(LST[i] + HA[i]/(np.pi/180.0),360.0)
         else:
            RA[i] = np.mod(LST[i] - HA[i]/(np.pi/180.0),360.0)
   except IndexError:
         if x < 0.0:
            RA = np.mod(LST + HA/(np.pi/180.0),360.0)
         else :
            RA = np.mod(LST - HA/(np.pi/180.0),360.0)

   dec = DEC*(180.0/np.pi)
   ra = RA

   return (ra,dec,az*180.0/np.pi,ALT*180.0/np.pi)

