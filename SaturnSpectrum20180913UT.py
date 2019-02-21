# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 07:22:23 2014

@author: steven.hill
"""
drive="f:"
import sys
sys.path.append('f:\\Astronomy\Python Play')
sys.path.append(drive+'\\Astronomy\Python Play\SpectroPhotometry\Spectroscopy')
sys.path.append(drive+'\\Astronomy\Python Play\TechniquesLibrary')

import matplotlib.pyplot as pl
import pylab
import numpy as np
from scipy import interpolate
import scipy
from copy import deepcopy
import ComputeEW
import SpectrumFromFITS1 as SFF1
import GeneralSpecUtils as GSU
import PlotUtils as PU
import copy
import scipy.stats as ST
from PyAstronomy import pyasl

####Custom single-observation code for Saturn on 2018-09-13 UT
# Maybe a transition to newer, more modudular code?
# Read and reshape spectral data files

#This is only for a single Saturn observation (two were made on 9/13)
#I need to double check on apertures and resulting flux values
#I need to incorporate EW measurements (raw, so including solar & telluric)
#I need to write output files
#I need to have linear scales set for albedo, or maybe just do that on 
#  the all-years plot.
#I need to make sure all lines appear in the legend

VegaArray,Vega,VegaMeta=SpectrumFromFITS1("Vega","Vega","Spectra","20180913UT")
RefPath="F:/Astronomy/Python Play/SPLibraries/SpectralReferenceFiles/ReferenceLibrary/"
Ref_a0v = scipy.loadtxt(RefPath+"a0v.dat", dtype=float, skiprows=3,usecols=(0,1))
Ref_a0v[:,0]=Ref_a0v[:,0]/10.
VegaCalibrated=GSU.SpectrumMath(Vega,Ref_a0v,"Divide")
temp=PU.Draw_with_Conf_Level(VegaCalibrated,1.0,'b',"SaturnCalibrated")
VegaCalibrated[:,1]=pyasl.smooth(VegaCalibrated[:,1],3,'flat')

SaturnArray,Saturn,SaturnMeta=SpectrumFromFITS1("Saturn","Saturn","Spectra","20180913UT")
Saturn2=copy.deepcopy(Saturn)
Saturn2[:,1]=SaturnArray[:,0]
SaturnCalibrated=GSU.SpectrumMath(Saturn,VegaCalibrated,"Divide")
SaturnCalibrated2=GSU.SpectrumMath(Saturn2,VegaCalibrated,"Divide")

Ref_g2v = scipy.loadtxt(RefPath+"g2v.dat", dtype=float, skiprows=3,usecols=(0,1))
Ref_g2v[:,0]=Ref_g2v[:,0]/10.
Ref_g2v[:,1]=pyasl.smooth(Ref_g2v[:,1],3,'flat')
SaturnAlbedo=GSU.SpectrumMath(SaturnCalibrated,Ref_g2v,"Divide")
SaturnAlbedo2=GSU.SpectrumMath(SaturnCalibrated2,Ref_g2v,"Divide")

temp=PU.Draw_with_Conf_Level(SaturnCalibrated,1.0e7,'b',"Saturn Calibrated")
temp=PU.Draw_with_Conf_Level(SaturnCalibrated2,1.0e7,'b',"Saturn Calibrated2")

temp=PU.Draw_with_Conf_Level(SaturnAlbedo,1.0e7,'k',"Saturn Albedo")
temp=PU.Draw_with_Conf_Level(SaturnAlbedo2,1.0e7,'k',"Saturn Albedo2")

Albedoarray=np.zeros((SaturnAlbedo.shape[0],2),dtype=float)
Albedoarray[:,0]=SaturnAlbedo[:,1]
Albedoarray[:,1]=SaturnAlbedo2[:,1]

ZeroIndices=np.where(Albedoarray <= 0.)
Albedoarray[ZeroIndices]=np.nan
#pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")
#pl.plot(wave,signalarray[:,0])        
AvgSignal=np.nanmean(Albedoarray,axis=1)
std=np.nanstd(Albedoarray,axis=1) 
sem=ST.sem(Albedoarray,axis=1,ddof=0,nan_policy='omit')
    
MeanAlbedo=np.zeros([SaturnAlbedo.shape[0],4])
MeanAlbedo[:,0]=SaturnAlbedo[:,0]
MeanAlbedo[:,1]=AvgSignal
MeanAlbedo[:,2]=std
MeanAlbedo[:,3]=sem

Temp=PU.Draw_with_Conf_Level(MeanAlbedo,1.0e7,'r','Saturn Albedo Mean')

pl.legend(loc=0,ncol=1, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.10, bottom=0.14, right=0.98, top=0.92,
            wspace=None, hspace=None)
pl.savefig("Saturn_Mean_Albedo_20180913UT"+"_1D.png",dpi=300)

MASTER=copy.deepcopy(SaturnAlbedo2)
MASTER[:,1]=MeanAlbedo[:,1]

np.savetxt("Saturn_Mean_Albedo_20180913UT"+"_1D.txt",MASTER,delimiter=" ",fmt="%10.3F %10.7F")

#np.savetxt(path.output_path+Target+'_'+meta.DateKey+'_1D_WVCal.txt',spectrum,delimiter=" ",fmt="%10.3F %10.7F")
#Compute EWs for telluric bands from MASTER
EWFN="SaturnEW20180913UT.txt"
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca II H&K",392.,399.,2.,EWFN,False)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Delta",409.,411.5,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Ca I 'g band'",421.5,424.,2.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"G band 430",427.,433.,2.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Gamma",433.,435.,1.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Fe I 'd band'",436.5,440.5,2.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Beta",481.5,489.0,1.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Mg 517",514.,520.,2.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Na D",587.,592.,2.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Alpha 656.3 band",653.,660.,1.,EWFN,True)

BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"O2 630 band",626.,636.,2.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"O2 B band",684.,698.,2.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H2O 720a band",691.,709.,1.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H2O 720b band",715.,735.,2.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"O2 A band",756.,770.,2.,EWFN,True)

BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 543",536.,550.,2.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 619",610.,628.,2.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 725",710.,740.,3.,EWFN,True)
"""
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
np.savetxt("SaturnAlbedo20180913UT.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")
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
pl.title("Saturn Spectrum 20180913UT",fontsize=9)
pl.plot(Ref[:,0]/10.,CLRonRef/(ExposureCLR*Aperture*NativeDispersionNM),label='CLR',linewidth=0.5)
pl.plot(Ref[:,0]/10.,NIRonRef/(ExposureNIR*Aperture*NativeDispersionNM),label='NIRonCLR',linewidth=0.5)
pl.plot(MASTER[:,0]/10.,MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM),label='MASTER',color='k',linewidth=1)
pl.plot(MASTER[:,0]/10.,ToA//(ExposureCLR*Aperture*NativeDispersionNM),label='Top of Atm.')
pl.plot(MASTER[:,0]/10.,Ref[:,1]*1e6,label='Solar Ref. x 1e6')
pl.plot(MASTER[:,0]/10.,NormAlbedo*1e6,label='Norm. Albedo x 1e6')

pl.legend(loc=0,ncol=6, borderaxespad=0.,prop={'size':6})
pylab.savefig('SaturnSpectrum20180913UT.png',dpi=300)


TempMaster=MASTER
TempMaster[:,0]=MASTER[:,0]/10.
TempMaster[:,1]=MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM)
np.savetxt("SaturnSpectrum20180913UT.txt",MASTER,delimiter=" ",fmt="%10.3F %10.7F")
np.savetxt("SaturnAlbedo20180913UT.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")
"""