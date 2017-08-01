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
from scipy import interpolate
import scipy
from copy import deepcopy
import ComputeEW

# SKIP PANDAS!!!

# Read and reshape spectral data files    
INT = scipy.fromfile(file="Data/SaturnInt-UnNorm-Cropped.dat", dtype=float, count=-1, sep='\t')    
NormResponsewithWV= scipy.fromfile(file="VegaResponse20130921UT.txt", dtype=float, count=-1, sep=" ")

INT=scipy.reshape(INT,[INT.size/2,2])
NativeDispersion=(INT[(INT.size/2.-1),0]-INT[0,0])/(INT.size/2.-1.)
NRespWV=scipy.reshape(NormResponsewithWV,[NormResponsewithWV.size/2,2])
#Load Reference Spectrum: Average G2v for albedo calculations
Ref = scipy.loadtxt("g:/Astronomy/SpectralReferenceFiles/ReferenceLibrary/g2v.dat", dtype=float, skiprows=3,usecols=(0,1))
MasterDispersion=(Ref[(Ref.size/2.-1),0]-Ref[0,0])/(Ref.size/2.-1.)

#Interpolate NIR, Response and Reference spectra onto Reference Wavelengths

INTInterp=interpolate.interp1d(INT[:,0],INT[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
INTonRef=INTInterp(Ref[:,0])

NRespInterp=interpolate.interp1d(NRespWV[:,0],NRespWV[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
NResponRef=NRespInterp(Ref[:,0])

#Create Master Observed Spectrum by merging CLR and NIR spectra

MASTER=deepcopy(Ref)
MASTER[:,1]= INTonRef[:]
#Compute EWs for telluric bands from MASTER
"""
EWFN="SaturnEW20140629UT.txt"
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca II H&K",3920.,3990.,20.,EWFN,False)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Delta",4090.,4115.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca I 'g band'",4215.,4240.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"G band 4300",4270.,4330.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Gamma",4330.,4350.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Fe I 'd band'",4365.,4405.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Beta",4815.,4890.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Mg 5170",5140.,5200.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Na D",5870.,5920.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Alpha 6563 band",6530.,6600.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca II 8498",8480.,8510.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca II 8542",8520.,8560.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca II 8662",8625.,8675.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"O2 6300 band",6260.,6360.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"O2 B band",6840.,6980.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H2O 7200a band",6910.,7090.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H2O 7200b band",7150.,7350.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"O2 A band",7560.,7700.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H2O Z band",7850.,8600.,40.,EWFN,True)

BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 6190",6130.,6210.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 7250",7100.,7400.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 8620",8480.,8660.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 8890",8700.,9200.,30.,EWFN,True)
"""

#Compute top of atmosphere spectrum

BoxcarSmoothIndices=np.where(Ref[:,0] >= 3500.)
NResponRefSmooth=deepcopy(NResponRef)
NResponRefSmooth[BoxcarSmoothIndices]=scipy.convolve(NResponRef[BoxcarSmoothIndices],np.ones((29,))/29)[(28):]
ToA=deepcopy(MASTER)
ToA=MASTER[:,1]/NResponRefSmooth

#Compute Albedo
RefSmooth=deepcopy(Ref[:,1])
RefSmooth[BoxcarSmoothIndices]=scipy.convolve(RefSmooth[BoxcarSmoothIndices],np.ones((19,))/19)[(18):]
Albedo=ToA/RefSmooth
mAlbedo = np.ma.masked_invalid(Albedo)
AlbedoNormRangeIndices=np.where((Ref[:,0] >4000.) & \
     (Ref[:,0] < 7500.))

NormAlbedo=Albedo/mAlbedo[AlbedoNormRangeIndices].max()
print NormAlbedo.max()

NormAlbedowithWV=deepcopy(Ref)
NormAlbedowithWV[:,1]=NormAlbedo

#Begin plotting 
SaturnPhotometry_Labels=np.array(['NUV','BLU','GRN','CLR','RED','HAL','NIR','CH4'])
SaturnPhotometry_WaveCenters=np.array([3800.,4500.,5410.,5705.,6480.,6563.,7350.,8990.])
SaturnPhotometry_20110608UT=np.array([2211.,95670.,201177.,np.nan,138751.,np.nan,70572.,np.nan])
SaturnPhotometry_20130610UT=np.array([605.,np.nan,57777.,np.nan,np.nan,60335.,22306.,1440.]) ;


pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=350
x1=950
xtks=13

y0=1.0e4
y1=1.0e8
#    ytks=9
ExposureCLR = 1. #seconds - already exposure normalized
Aperture = 0.2**2.-0.07**2. #meters^2
NativeDispersionNM=NativeDispersion/10.
MasterDispersionNM=MasterDispersion/10.
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
pl.title("Saturn Integrated Spectrum 20130612UT",fontsize=9)
#pl.plot(Ref[:,0],CLRonRef/150.,label='CLR',linewidth=0.5)
#pl.plot(Ref[:,0],NIRonRef/300.,label='NIRonCLR',linewidth=0.5)
pl.plot(MASTER[:,0]/10.,MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM),label='Integrated',color='k',linewidth=1)
pl.plot(MASTER[:,0]/10.,ToA/(ExposureCLR*Aperture*NativeDispersionNM),label='Top of Atm.')
pl.plot(MASTER[:,0]/10.,Ref[:,1]*1.0e6,label='Solar Ref. x 1e6')
pl.plot(MASTER[:,0]/10.,NormAlbedo*1.0e6,label='Norm. Albedo x 1e6')

pl.plot(SaturnPhotometry_WaveCenters/10.,0.25*10.*SaturnPhotometry_20110608UT/Aperture,
        label='20110608UT',linewidth=0,marker='+',markersize=5,markeredgewidth=1,color='b')
pl.plot(SaturnPhotometry_WaveCenters/10.,10.*SaturnPhotometry_20130610UT/Aperture,
        label='20130610UT',linewidth=0,marker='x',markersize=3.5,markeredgewidth=1,color='b')


pl.legend(loc=0,ncol=3, borderaxespad=0.,prop={'size':6})
pylab.savefig('SaturnINTSpectrum20130612UT.png',dpi=300)

TempMaster=MASTER
TempMaster[:,0]=MASTER[:,0]/10.
TempMaster[:,1]=MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM)
np.savetxt("SaturnINTSpectrum20130612UT.txt",MASTER,delimiter=" ",fmt="%10.3F %10.7F")
np.savetxt("SaturnINTAlbedo20130612UT.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")
