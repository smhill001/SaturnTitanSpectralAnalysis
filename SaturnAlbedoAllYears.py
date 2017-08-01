# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 07:22:23 2014

@author: steven.hill
"""
import sys
sys.path.append('g:\\Astronomy\Python Play')
import matplotlib.pyplot as pl
import pylab
import numpy as np
import scipy
from scipy import interpolate
from copy import deepcopy

# SKIP PANDAS!!!

# Read and reshape spectral data files    
SaturnINT_20130612UT = scipy.fromfile(file="../../Saturn/Spectral Data/SaturnINTAlbedo20130612UT.txt", dtype=float, count=-1, sep='\t')    
SaturnINT_20130612UT=scipy.reshape(SaturnINT_20130612UT,[SaturnINT_20130612UT.size/2,2])
SaturnDSK_20130612UT = scipy.fromfile(file="../../Saturn/Spectral Data/SaturnDSKAlbedo20130612UT.txt", dtype=float, count=-1, sep='\t')    
SaturnDSK_20130612UT=scipy.reshape(SaturnDSK_20130612UT,[SaturnDSK_20130612UT.size/2,2])
SaturnRNG_20130612UT = scipy.fromfile(file="../../Saturn/Spectral Data/SaturnRNGAlbedo20130612UT.txt", dtype=float, count=-1, sep='\t')    
SaturnRNG_20130612UT=scipy.reshape(SaturnRNG_20130612UT,[SaturnRNG_20130612UT.size/2,2])
SaturnINT_20140629UT = scipy.fromfile(file="SaturnAlbedo20140629UT.txt", dtype=float, count=-1, sep='\t')    
SaturnINT_20140629UT=scipy.reshape(SaturnINT_20140629UT,[SaturnINT_20140629UT.size/2,2])

Saturn_Karkoschka1993 = scipy.fromfile(file="../../Saturn/Spectral Data/Karkoschka/1993.tab.txt", dtype=float, count=-1, sep=" ")    
Saturn_Karkoschka1993=scipy.reshape(Saturn_Karkoschka1993,[Saturn_Karkoschka1993.size/8,8])

Saturn_KarkRef1993=np.zeros((Saturn_Karkoschka1993.size/8,2))
Saturn_KarkRef1993[:,0]=Saturn_Karkoschka1993[:,0]*10.
Saturn_KarkRef1993[:,1]=Saturn_Karkoschka1993[:,4]
print "***Saturn_KarkRef1993***", Saturn_KarkRef1993

Saturn_Avg_Smooth=SaturnINT_20140629UT[:,1]
N=5
Saturn_Avg_Smooth=scipy.convolve(Saturn_Avg_Smooth,np.ones((N,))/N)[(N-1):]
Saturn_Avg_Smooth_WV=SaturnINT_20140629UT[:,0]+7.0*(N/2-0.5)

Saturn2014Smth=deepcopy(SaturnINT_20140629UT)
Saturn2014Smth[:,0]=Saturn_Avg_Smooth_WV
Saturn2014Smth[:,1]=Saturn_Avg_Smooth

StartWV=4000.
EndWV=7500.
LowCutIndices=np.where((SaturnINT_20140629UT[:,0] < StartWV) )
SaturnINT_20140629UT[LowCutIndices,1]=np.nan
HighCutIndices=np.where((SaturnINT_20140629UT[:,0] > EndWV) )
SaturnINT_20140629UT[HighCutIndices,1]=np.nan

LowCutIndices=np.where((Saturn2014Smth[:,0] < StartWV) )
Saturn2014Smth[LowCutIndices,1]=np.nan
HighCutIndices=np.where((Saturn2014Smth[:,0] > EndWV) )
Saturn2014Smth[HighCutIndices,1]=np.nan
"""
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

pl.figure(figsize=(6.5, 1.5), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=375
x1=950

xtks=24
y0=0.0
y1=0.7

# Set x limits
pl.xlim(x0,x1)
# Set x ticks
pl.xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
pl.ylim(y0,y1)
# Set y ticks
pl.grid()
pl.tick_params(axis='both', which='major', labelsize=7)
pl.ylabel(r"$Albedo$",fontsize=7)
pl.xlabel(r"$Wavelength (nm)$",fontsize=7)
"""
pl.title("Saturn Albedo - All Data",fontsize=9)
pl.plot(SaturnINT_20130612UT[:,0]/10.,SaturnINT_20130612UT[:,1]*0.60,label='INT_20130612UT',linewidth=1)
pl.plot(SaturnDSK_20130612UT[:,0]/10.,SaturnDSK_20130612UT[:,1]*0.60,label='DSK_20130612UT',linewidth=0.5)
pl.plot(SaturnRNG_20130612UT[:,0]/10.,SaturnRNG_20130612UT[:,1]*0.60,label='RNG_20130612UT',linewidth=0.5)
pl.plot(SaturnINT_20140629UT[:,0]/10.,SaturnINT_20140629UT[:,1]*0.60,label='INT_20140629UT',linewidth=1,color='k')
pl.plot(Saturn2014Smth[:,0]/10.,Saturn2014Smth[:,1]*0.59,label='2014 Smooth',linewidth=1,color='r')

pl.plot(Saturn_KarkRef1993[:,0]/10.,Saturn_KarkRef1993[:,1],label='Karkoschka, 1993',linewidth=1,color='0.5')
"""
pl.title("Saturn",fontsize=9)
#pl.plot(SaturnINT_20130612UT[:,0]/10.,SaturnINT_20130612UT[:,1]*0.60,label='INT_20130612UT',linewidth=1)
#pl.plot(SaturnDSK_20130612UT[:,0]/10.,SaturnDSK_20130612UT[:,1]*0.60,label='DSK_20130612UT',linewidth=0.5)
#pl.plot(SaturnRNG_20130612UT[:,0]/10.,SaturnRNG_20130612UT[:,1]*0.60,label='RNG_20130612UT',linewidth=0.5)
pl.plot(SaturnINT_20140629UT[:,0]/10.,SaturnINT_20140629UT[:,1]*0.60,label='This Work (2014) 0.006m',linewidth=1)
#pl.plot(Saturn2014Smth[:,0]/10.,Saturn2014Smth[:,1]*0.59,label='2014 Smooth',linewidth=1,color='g')

pl.plot(Saturn_KarkRef1993[:,0]/10.,Saturn_KarkRef1993[:,1],label='Karkoschka, 1994',linewidth=1,color='0.5')


pl.legend(loc=0,ncol=2, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.06, bottom=0.25, right=0.98, top=0.88,
                wspace=None, hspace=None)
                
pylab.savefig('SaturnAlbedoAllYears.png',dpi=300)