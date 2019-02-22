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
import numpy as np
import scipy
from copy import deepcopy
import EquivWidthUtils as EWU
import SpectrumFromFITS1 as SFF1
import GeneralSpecUtils as GSU
import PlotUtils as PU
import ConfigFiles as CF
import scipy.stats as ST
from PyAstronomy import pyasl

####Custom single-observation code for Saturn on 2018-09-13 UT
# Maybe a transition to newer, more modudular code?
# Read and reshape spectral data files

#I need to double check on apertures and resulting flux values
#I need to have linear scales set for albedo, or maybe just do that on 
#  the all-years plot.

#######################
V=CF.Target_Parameters("f:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy/Target_Parameters.txt")
V.loadtargetparams("Vega")
VegaPath=CF.built_path(V)
VegaPath.spectra("20180913UT")

VegaArray,Vega,VegaMeta=SFF1.SpectrumFromFITS1("Vega","Vega","Spectra","20180913UT")

Ref_a0v = scipy.loadtxt(VegaPath.reference_path+V.SpecType, dtype=float, skiprows=3,usecols=(0,1))
Ref_a0v[:,0]=Ref_a0v[:,0]/10.
VegaCalibrated=GSU.SpectrumMath(Vega,Ref_a0v,"Divide")
temp=PU.Draw_with_Conf_Level(VegaCalibrated,1.0,'b',"SaturnCalibrated")
VegaCalibrated[:,1]=pyasl.smooth(VegaCalibrated[:,1],3,'flat')
#######################
S=CF.Target_Parameters("f:/Astronomy/Python Play/SpectroPhotometry/Spectroscopy/Target_Parameters.txt")
S.loadtargetparams("Saturn")
SaturnPath=CF.built_path(S)
SaturnPath.spectra("20180913UT")

SaturnArray,Saturn,SaturnMeta=SFF1.SpectrumFromFITS1("Saturn","Saturn","Spectra","20180913UT")

Saturn2=deepcopy(Saturn)
Saturn2[:,1]=SaturnArray[:,0]
SaturnCalibrated=GSU.SpectrumMath(Saturn,VegaCalibrated,"Divide")
SaturnCalibrated2=GSU.SpectrumMath(Saturn2,VegaCalibrated,"Divide")

Ref_g2v = scipy.loadtxt(SaturnPath.reference_path+S.SpecType, dtype=float, skiprows=3,usecols=(0,1))
Ref_g2v[:,0]=Ref_g2v[:,0]/10.
Ref_g2v[:,1]=pyasl.smooth(Ref_g2v[:,1],3,'flat')
SaturnAlbedo=GSU.SpectrumMath(SaturnCalibrated,Ref_g2v,"Divide")
SaturnAlbedo2=GSU.SpectrumMath(SaturnCalibrated2,Ref_g2v,"Divide")
#######################
temp=PU.Draw_with_Conf_Level(SaturnCalibrated,1.0e7,'b',"Saturn Calibrated")
temp=PU.Draw_with_Conf_Level(SaturnCalibrated2,1.0e7,'b',"Saturn Calibrated2")

temp=PU.Draw_with_Conf_Level(SaturnAlbedo,1.0e7,'k',"Saturn Albedo")
temp=PU.Draw_with_Conf_Level(SaturnAlbedo2,1.0e7,'k',"Saturn Albedo2")

Albedoarray=np.zeros((SaturnAlbedo.shape[0],2),dtype=float)
Albedoarray[:,0]=SaturnAlbedo[:,1]
Albedoarray[:,1]=SaturnAlbedo2[:,1]
#######################
ZeroIndices=np.where(Albedoarray <= 0.)
Albedoarray[ZeroIndices]=np.nan
AvgSignal=np.nanmean(Albedoarray,axis=1)
std=np.nanstd(Albedoarray,axis=1) 
sem=ST.sem(Albedoarray,axis=1,ddof=0,nan_policy='omit')
    
MeanAlbedo=np.zeros([SaturnAlbedo.shape[0],4])
MeanAlbedo[:,0]=SaturnAlbedo[:,0]
MeanAlbedo[:,1]=AvgSignal
MeanAlbedo[:,2]=std
MeanAlbedo[:,3]=sem
########################
Temp=PU.Draw_with_Conf_Level(MeanAlbedo,1.0e7,'r','Saturn Albedo Mean')

pl.legend(loc=0,ncol=1, borderaxespad=0.,prop={'size':6})
pl.subplots_adjust(left=0.10, bottom=0.14, right=0.98, top=0.92,
            wspace=None, hspace=None)
pl.savefig("Saturn_Mean_Albedo_20180913UT"+"_1D.png",dpi=300)

MASTER=deepcopy(SaturnAlbedo2)
MASTER[:,1]=MeanAlbedo[:,1]

np.savetxt("Saturn_Mean_Albedo_20180913UT"+"_1D.txt",MASTER,delimiter=" ",fmt="%10.3F %10.7F")
########################
Bands=EWU.LinesBands_to_Measure("Saturn_ObsBands_135mm100lpm.txt")
Bands.load_records(WVRange=[350.,750.])

EWFN="SaturnEW20180913UT.txt"
flag=False
for B in range(0,len(Bands.ID)):
    Temp=EWU.ComputeEW(MASTER,Bands.ID[B],Bands.WV0[B],Bands.WV1[B],Bands.WVCont[B],EWFN,flag)
    flag=True

