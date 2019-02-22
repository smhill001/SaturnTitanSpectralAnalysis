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
import EquivWidthUtils as EWU

# Read and reshape spectral data files    
CLR = scipy.fromfile(file="../20140629UT/Saturn-Spectrum-20140629UT-CLR-sum2m30s-Rotated-Cropped-WVCal.dat", dtype=float, count=-1, sep='\t')    
NIR = scipy.fromfile(file="../20140629UT/Saturn-Spectrum-20140629UT-NIR-sum5min0s-Rotated-Cropped-ABWvCal.dat", dtype=float, count=-1, sep='\t')    
NormResponsewithWV= scipy.fromfile(file="SpicaResponse20140629UT.txt", dtype=float, count=-1, sep=" ")
CLR=scipy.reshape(CLR,[CLR.size/2,2])
NativeDispersion=(CLR[(CLR.size/2-1),0]-CLR[0,0])/(CLR.size/2-1)
NIR=scipy.reshape(NIR,[NIR.size/2,2])
NRespWV=scipy.reshape(NormResponsewithWV,[NormResponsewithWV.size/2,2])
MasterDispersion=(NRespWV[(NRespWV.size/2-1),0]-NRespWV[0,0])/(NRespWV.size/2-1)

#Load Reference Spectrum: Average G2v for albedo calculations
RefPath="F:/Astronomy/Python Play/SPLibraries/SpectralReferenceFiles/ReferenceLibrary/"
Ref = scipy.loadtxt(RefPath+"g2v.dat", dtype=float, skiprows=3,usecols=(0,1))

#Interpolate NIR, Response and Reference spectra onto CLR Wavelengths

CLRInterp=interpolate.interp1d(CLR[:,0],CLR[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
CLRonRef=CLRInterp(Ref[:,0])

NIRInterp=interpolate.interp1d(NIR[:,0],NIR[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
NIRonRef=NIRInterp(Ref[:,0])

NRespInterp=interpolate.interp1d(NRespWV[:,0],NRespWV[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
NResponRef=NRespInterp(Ref[:,0])

#Create Master Observed Spectrum by merging CLR and NIR spectra

MergeStartWV=7200.
MergeEndWV=7500.

OverlapRegionIndices=np.where((Ref[:,0] >MergeStartWV) & \
     (Ref[:,0] < MergeEndWV))

NIRScaling2CLR= scipy.mean(CLRonRef[OverlapRegionIndices]) \
    / scipy.mean(NIRonRef[OverlapRegionIndices])

MASTER=deepcopy(Ref)
MASTER[OverlapRegionIndices,1]= \
    (CLRonRef[OverlapRegionIndices]+NIRonRef[OverlapRegionIndices]*NIRScaling2CLR)/2.

DedicatedNIRIndices=np.where(Ref[:,0] >= 7500.)
DedicatedCLRIndices=np.where(Ref[:,0] <= 7200.)
MASTER[DedicatedNIRIndices,1]= NIRonRef[DedicatedNIRIndices]*NIRScaling2CLR
MASTER[DedicatedCLRIndices,1]= CLRonRef[DedicatedCLRIndices]

#Compute EWs for telluric bands from MASTER
Bands=EWU.LinesBands_to_Measure("Saturn_ObsBands_135mm100lpm.txt")
Bands.load_records()

MASTER[:,0]=MASTER[:,0]/10.

EWFN="SaturnEW20140629UT.txt"
flag=False
for B in range(0,len(Bands.ID)):
    Temp=EWU.ComputeEW(MASTER,Bands.ID[B],Bands.WV0[B],Bands.WV1[B],Bands.WVCont[B],EWFN,flag)
    flag=True

NativeDispersionNM=NativeDispersion/10.
MasterDispersionNM=MasterDispersion/10.

#Compute top of atmosphere spectrum

BoxcarSmoothIndices=np.where(Ref[:,0] >= 8500.)
NResponRefSmooth=deepcopy(NResponRef)
NResponRefSmooth[BoxcarSmoothIndices]=scipy.convolve(NResponRef[BoxcarSmoothIndices],np.ones((29,))/29)[(28):]
ToA=deepcopy(MASTER)
ToA=MASTER[:,1]/NResponRefSmooth

#Compute Albedo
Albedo=ToA/Ref[:,1]
mAlbedo = np.ma.masked_invalid(Albedo)
AlbedoNormRangeIndices=np.where((Ref[:,0] >4000.) & \
     (Ref[:,0] < 7500.))

NormAlbedo=Albedo/mAlbedo[AlbedoNormRangeIndices].max()
print NormAlbedo.max()

NormAlbedowithWV=deepcopy(Ref)
NormAlbedowithWV[:,1]=NormAlbedo
np.savetxt("SaturnAlbedo20140629UT.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")
#Begin plotting 

pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=350
x1=950

xtks=13
y0=1.0e4
y1=1.0e8
#    ytks=9
ExposureCLR = 150. #seconds
ExposureNIR = 300. #seconds
Aperture = (0.135/22.)**2. #meters^2
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
pl.ylabel("Counts Per Second",fontsize=7)
pl.xlabel("Wavelength (A)",fontsize=7)
pl.title("Saturn Spectrum 20140629UT",fontsize=9)
pl.plot(Ref[:,0],CLRonRef/(ExposureCLR*Aperture*NativeDispersionNM),label='CLR',linewidth=0.5)
pl.plot(Ref[:,0],NIRonRef/(ExposureNIR*Aperture*NativeDispersionNM),label='NIRonCLR',linewidth=0.5)
pl.plot(MASTER[:,0],MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM),label='MASTER',color='k',linewidth=1)
pl.plot(MASTER[:,0],ToA//(ExposureCLR*Aperture*NativeDispersionNM),label='Top of Atm.')
pl.plot(MASTER[:,0],Ref[:,1]*1e6,label='Solar Ref. x 1e6')
pl.plot(MASTER[:,0],NormAlbedo*1e6,label='Norm. Albedo x 1e6')

pl.legend(loc=0,ncol=6, borderaxespad=0.,prop={'size':6})
pylab.savefig('SaturnSpectrum20140629UT.png',dpi=300)


TempMaster=MASTER
TempMaster[:,0]=MASTER[:,0]
TempMaster[:,1]=MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM)
np.savetxt("SaturnSpectrum20140629UT.txt",MASTER,delimiter=" ",fmt="%10.3F %10.7F")
np.savetxt("SaturnAlbedo20140629UT.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")