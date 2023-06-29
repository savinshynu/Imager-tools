"""
Plots the instrumental response, taken from LSL

"""

from scipy.interpolate import interp1d
import numpy
from lsl.common import stations


def getImpedanceMisMatch(freq):

    freq4, imm4 = stations.lwasv.antennas[0].response(dB=False)
    immIntp = interp1d(freq4, imm4, kind='cubic', bounds_error=False)
    imm = immIntp(freq)
    return imm


def getARXResponse(freq, filter='full', site=stations.lwasv):
    antennas = site.getAntennas()
    f,r = antennas[0].arx.response(filter='split')
    freq2 = f
    respX2 = numpy.zeros_like(r)
    respY2 = numpy.zeros_like(r)
    for i in xrange(len(antennas)):
        if antennas[i].getStatus() != 33:
           continue
        f,r = antennas[i].arx.response(filter=filter, dB=False)

        if antennas[i].pol == 0:
           respX2 += r
        else:
           respY2 += r

    respX2 /= respX2.max()
    respY2 /= respY2.max()

    respXIntp = interp1d(freq2, respX2, kind='cubic', bounds_error=False)
    respYIntp = interp1d(freq2, respY2, kind='cubic', bounds_error=False)

    respX = respXIntp(freq)
    respY = respYIntp(freq)
    respX = numpy.where(numpy.isfinite(respX), respX, 0)
    respY = numpy.where(numpy.isfinite(respY), respY, 0)
    stokes_i = respX+respY
    return stokes_i


