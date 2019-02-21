# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 07:22:23 2014

@author: steven.hill

FUNCTION: TitanSpectrum20130612UT.py

This function reads 1D spectral data from text files created in RSpect from 
FITS image finles that have been calibrated (bias, dark, but not necessarily 
flat). A 
linear wavelength calibration is applied and the resulting 1-dimensional spectrum
is re-gridded to standard 0.5nm wavelength intervals. Using reference spectra
(either instrument response functions or reference stellar spectra) the data
are converted to either Top of Atmosphere 'raw' fluxes or the relative 
(peak normalized to unity) response of the instrument. In addition, the
equivalent widths of lines or bands of interest are computed and written
to ASCII output files.

INPUT PARAMETERS:
    Target - This string identifies the target, e.g., "Jupiter", "Vega", etc.
             This parameter allows the lookup of the target type, e.g.,
             "Planet", "Star", etc. that permits construction of directory
             paths.
    DateUT - The UT date of the observation. Combined with Target, this forms
             a unique key to the observation, assuming that most unique parameters
             are invariant over a single observing night.
    FluxCalibration - This is a flag that determines whether the observation
             is to be processed to produce a system response (e.g. the target
             has a well defined spectrum) or if the calibrated spectrum is
             produced using a system spectral response function.
             
INPUT FILES:
    PlotParameters - List of plotting configuration info related only to the target,
             not the specific observation
    NNNNKeyFile.txt - List of Target-Date-Time keys
    ProcessingConfigFile.txt - Four parameters that control processing
    FileList - List of individual FITS files for processing
    ObsBands - List of spectral bands for which equiv. widths will be calculated
    ReferenceSpectrum - This may be a canonical reference, e.g., for a stellar
             spectral type for response computation, or it may be a response
             file
    FITS Images - 2D image files with necessary metadata for flux and wavelength
             calibration            
             
OUTPUTS:
    1D Spectrum - txt file
    1D Spectrum - png file
    1D Response or Albedo - txt file; depends on FluxCalibration flag - TBD 6/6/16
    Equivalent widths - txt file

NEXT STEPS as of 6/8/2016: (Need to verify which of these have been done)
1) Adapt this code to work multiple observations of Neptune in a single file list file.
2) Run through Neptune observations for the 2015-11-08UT and add keywords for processing
3) Close out raw and calibrated files for 11/08 so they can be zipped
4) Encapsulate target selection and associated file location meta data in
   text files and/or a python function
5) Bigger job: Make decisions on configuring EW parameters: Maybe as metadata in FITS files?
6) Make this so it full accommodates emission and reflected spectra, e.g., 
   for stars versus planets.
7) Make this so that the EWs can be calculated for either (or both?) albedo
   and/or total spectrum
8) Make smoothing of reference spectrum configurable (maybe in ProcessingConfig)
9) Provide options so that individual filter transmission can be assessed, e.g.,
   merge the Jovian Moons filter response code
10) Start using GIT?
11) Convert many functions to classes
12) Work on inheritance
"""

import sys
sys.path.append('f:\\Astronomy\Python Play')
import matplotlib.pyplot as pl
import pylab
import numpy as np
from scipy import interpolate
import scipy
from copy import deepcopy
import ComputeEW

datapath="f:/Astronomy/Projects/Planets/Saturn/Spectral Data/"
# Read and reshape spectral data files    
INT = scipy.fromfile(file=datapath+"1D Spectra/Titan-UnNorm-Cropped.dat", dtype=float, count=-1, sep='\t')    
INT=scipy.reshape(INT,[INT.size/2,2])
INT[:,1]=INT[:,1]/10. #Version of data read was 10x actual data values
NativeDispersion=(INT[(INT.size/2-1),0]-INT[0,0])/(INT.size/2-1)

NormResponsewithWV= scipy.fromfile(file="VegaResponse20130921UT.txt", dtype=float, count=-1, sep=" ")
NRespWV=scipy.reshape(NormResponsewithWV,[NormResponsewithWV.size/2,2])

SaturnRNG_20130612UT = scipy.fromfile(file=datapath+"SaturnTitanSpectralAnalysis/SaturnRNGSpectrum20130612UT.txt", dtype=float, count=-1, sep='\t')    
SaturnRNG_20130612UT=scipy.reshape(SaturnRNG_20130612UT,[SaturnRNG_20130612UT.size/2,2])

#Load Reference Spectrum: Average G2v for albedo calculations
Ref = scipy.loadtxt("f:/Astronomy/Python Play/SPLibraries/SpectralReferenceFiles/ReferenceLibrary/g2v.dat", dtype=float, skiprows=3,usecols=(0,1))
MasterDispersion=(Ref[(Ref.size/2-1),0]-Ref[0,0])/(Ref.size/2-1)

#Interpolate NIR, Response and Reference spectra onto Reference Wavelengths

INTInterp=interpolate.interp1d(INT[:,0],INT[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
INTonRef=INTInterp(Ref[:,0])

NRespInterp=interpolate.interp1d(NRespWV[:,0],NRespWV[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
NResponRef=NRespInterp(Ref[:,0])

RNGInterp=interpolate.interp1d(SaturnRNG_20130612UT[:,0]*10.,SaturnRNG_20130612UT[:,1],kind='linear', copy=True,
                         bounds_error=False, fill_value=0.0)  
RNGonRef=RNGInterp(Ref[:,0])


#Create Master Observed Spectrum by merging CLR and NIR spectra

MASTER=deepcopy(Ref)
MASTER[:,1]= INTonRef[:]
#Compute EWs for telluric bands from MASTER

EWFN="TitanEW20130612UT.txt"
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"G band 4300",4110.,4450.,20.,EWFN,False)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Beta",4760.,5040.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Mg 5170",5090.,5370.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"Na D",5800.,6000.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H Alpha 6563 band",6530.,6700.,10.,EWFN,True)

BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"O2 B band",6830.,6990.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H2O 7200a + CH4",6910.,7090.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H2O 7200b + CH4",7150.,7350.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"O2 A band",7540.,7760.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"H2O Z band + CH4",7850.,8600.,40.,EWFN,True)

BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 6190",6060.,6400.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 6680",6460.,6770.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 6825 + O2 B",6780.,6900.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 7050 + H2O 7200a",6920.,7120.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW.ComputeEW(MASTER,"CH4 7250 + H2O 7200b",7120.,7560.,30.,EWFN,True)


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

#Compute Rings Normalized Albedo
RNGAlbedo=INTonRef/RNGonRef
mRNGAlbedo = np.ma.masked_invalid(RNGAlbedo)
RNGAlbedoNormRangeIndices=np.where((Ref[:,0] >4000.) & \
     (Ref[:,0] < 7500.))

RNGNormAlbedo=RNGAlbedo/mRNGAlbedo[RNGAlbedoNormRangeIndices].max()
print RNGNormAlbedo.max()

RNGNormAlbedowithWV=deepcopy(Ref)
RNGNormAlbedowithWV[:,1]=RNGNormAlbedo

#Begin plotting 

TitanPhotometry_Labels=np.array(['NUV','BLU','GRN','CLR','RED','HAL','NIR','CH4'])
TitanPhotometry_WaveCenters=np.array([3800.,4500.,5410.,5705.,6480.,6563.,7350.,8990.])
TitanPhotometry_20130610UT=np.array([np.nan,np.nan,30.2,np.nan,np.nan,38.3,14.1,np.nan]) ;

pl.figure(figsize=(6.5, 2.5), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=350
x1=950
xtks=13

y0=1.0e1
y1=1.0e5
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
pl.title("Titan Spectrum 20130612UT",fontsize=9)
#pl.plot(Ref[:,0],CLRonRef/150.,label='CLR',linewidth=0.5)
#pl.plot(Ref[:,0],NIRonRef/300.,label='NIRonCLR',linewidth=0.5)
pl.plot(MASTER[:,0]/10.,MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM),label='Integrated',color='k',linewidth=1)
pl.plot(MASTER[:,0]/10.,ToA/(ExposureCLR*Aperture*NativeDispersionNM),label='Top of Atm.')
pl.plot(MASTER[:,0]/10.,Ref[:,1]*1.0e3,label='Solar Ref. x 1e3')
pl.plot(MASTER[:,0]/10.,NormAlbedo*1.0e3,label='Norm. Albedo x 1e3')
pl.plot(MASTER[:,0]/10.,RNGNormAlbedo*1.0e3,label='RNG Norm. Albedo x 1e3')




pl.plot(TitanPhotometry_WaveCenters/10.,10.*TitanPhotometry_20130610UT/Aperture,
        label='20130610UT',linewidth=0,marker='x',markersize=3.5,markeredgewidth=1,color='b')
        
pl.legend(loc=0,ncol=3, borderaxespad=0.,prop={'size':6})
pylab.savefig('TitanSpectrum20130612UT.png',dpi=300)

TempMaster=MASTER
TempMaster[:,0]=MASTER[:,0]/10.
TempMaster[:,1]=MASTER[:,1]/(ExposureCLR*Aperture*NativeDispersionNM)
np.savetxt("TitanSpectrum20130612UT.txt",MASTER,delimiter=" ",fmt="%10.3F %10.7F")
np.savetxt("TitanAlbedo20130612UT.txt",NormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")
np.savetxt("TitanRNGAlbedo20130612UT.txt",RNGNormAlbedowithWV,delimiter=" ",fmt="%10.3F %10.7F")
