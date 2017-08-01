# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 07:22:23 2014

@author: steven.hill
"""
import sys
sys.path.append('D:\\Astronomy\Python Play')
import matplotlib.pyplot as pl
import pylab
import numpy as np
import scipy
from scipy import interpolate

# SKIP PANDAS!!!

# Read and reshape spectral data files    
SaturnINT_20130612UT = scipy.fromfile(file="../../Saturn Project 2013/Spectroscopy/SaturnINTSpectrum20130612UT.txt", dtype=float, count=-1, sep='\t')    
SaturnINT_20130612UT=scipy.reshape(SaturnINT_20130612UT,[SaturnINT_20130612UT.size/2,2])
SaturnDSK_20130612UT = scipy.fromfile(file="../../Saturn Project 2013/Spectroscopy/SaturnDSKSpectrum20130612UT.txt", dtype=float, count=-1, sep='\t')    
SaturnDSK_20130612UT=scipy.reshape(SaturnDSK_20130612UT,[SaturnDSK_20130612UT.size/2,2])
SaturnRNG_20130612UT = scipy.fromfile(file="../../Saturn Project 2013/Spectroscopy/SaturnRNGSpectrum20130612UT.txt", dtype=float, count=-1, sep='\t')    
SaturnRNG_20130612UT=scipy.reshape(SaturnRNG_20130612UT,[SaturnRNG_20130612UT.size/2,2])
SaturnINT_20140629UT = scipy.fromfile(file="SaturnSpectrum20140629UT.txt", dtype=float, count=-1, sep='\t')    
SaturnINT_20140629UT=scipy.reshape(SaturnINT_20140629UT,[SaturnINT_20140629UT.size/2,2])

"""StartWV=400.
EndWV=700.
LowCutIndices=np.where((Neptune_20141028UT[:,0] < StartWV) )
Neptune_20141028UT[LowCutIndices,1]=np.nan
HighCutIndices=np.where((Neptune_20141028UT[:,0] > EndWV) )
Neptune_20141028UT[HighCutIndices,1]=np.nan

StartWV=400.
EndWV=700.
LowCutIndices=np.where((Neptune_20141029UT[:,0] < StartWV) )
Neptune_20141029UT[LowCutIndices,1]=np.nan
HighCutIndices=np.where((Neptune_20141029UT[:,0] > EndWV) )
Neptune_20141029UT[HighCutIndices,1]=np.nan

Neptune_20141028UT_Interp=interpolate.interp1d(Neptune_20141028UT[:,0],Neptune_20141028UT[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
Neptune_20141028UT_on_20141029UT=Neptune_20141028UT_Interp(Neptune_20141029UT[:,0])

Neptune_Avg_20141028and29=Neptune_20141028UT_on_20141029UT*0.25+Neptune_20141029UT[:,1]*0.75
"""

SaturnPhotometry_Labels=np.array(['NUV','BLU','GRN','CLR','RED','HAL','NIR','CH4'])
SaturnPhotometry_WaveCenters=np.array([3800.,4500.,5410.,5705.,6480.,6563.,7350.,8990.])
SaturnPhotometry_20110608UT=np.array([2211.,95670.,201177.,np.nan,138751.,np.nan,70572.,np.nan])
SaturnPhotometry_20130610UT=np.array([605.,np.nan,57777.,np.nan,np.nan,60335.,22306.,1440.]) ;

pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=375
x1=900

xtks=22
y0=1.e4
y1=1.e8

# Set x limits
pl.xlim(x0,x1)
# Set x ticks
pl.xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
pl.ylim(y0,y1)
# Set y ticks
pl.yscale('log')
pl.grid()
pl.tick_params(axis='both', which='major', labelsize=7)
pl.ylabel(r"$Counts-s^{-1}$-$m^{-2}$-$nm^{-1}$",fontsize=7)
pl.xlabel(r"$Wavelength (nm)$",fontsize=7)
pl.title("Saturn Spectrum - All Data",fontsize=9)
pl.plot(SaturnINT_20130612UT[:,0],SaturnINT_20130612UT[:,1],label='INT_20130612UT',linewidth=1)
pl.plot(SaturnDSK_20130612UT[:,0],SaturnDSK_20130612UT[:,1],label='DSK_20130612UT',linewidth=0.5)
pl.plot(SaturnRNG_20130612UT[:,0],SaturnRNG_20130612UT[:,1],label='RNG_20130612UT',linewidth=0.5)
pl.plot(SaturnINT_20140629UT[:,0],SaturnINT_20140629UT[:,1],label='INT_20140629UT',linewidth=1,color='k')

PhotometryAperture = 0.2**2.-0.07**2. #meters^2

pl.plot(SaturnPhotometry_WaveCenters/10.,0.25*10.*SaturnPhotometry_20110608UT/PhotometryAperture,
        label='20110608UT',linewidth=0,marker='+',markersize=5,markeredgewidth=1,color='b')
pl.plot(SaturnPhotometry_WaveCenters/10.,10.*SaturnPhotometry_20130610UT/PhotometryAperture,
        label='20130610UT',linewidth=0,marker='x',markersize=3.5,markeredgewidth=1,color='b')


pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})

pylab.savefig('SaturnSpectrumAllYears.png',dpi=300)