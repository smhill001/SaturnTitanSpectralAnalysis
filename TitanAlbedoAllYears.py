# -*- coding: utf-8 -*-
"""
Created on Tue Jul 15 07:22:23 2014

@author: steven.hill
"""
import sys
sys.path.append('f:\\Astronomy\Python Play')
import matplotlib.pyplot as pl
import pylab
import numpy as np
import scipy
from scipy import interpolate
import ComputeEW1
from copy import deepcopy
from PyAstronomy import pyasl
import datetime

# SKIP PANDAS!!!

# Read and reshape spectral data files    
Titan_20130612UT = scipy.fromfile(file="TitanAlbedo20130612UT.txt", dtype=float, count=-1, sep='\t')    
Titan_20130612UT=scipy.reshape(Titan_20130612UT,[Titan_20130612UT.size/2,2])
TitanRNG_20130612UT = scipy.fromfile(file="TitanRNGAlbedo20130612UT.txt", dtype=float, count=-1, sep='\t')    
TitanRNG_20130612UT=scipy.reshape(TitanRNG_20130612UT,[TitanRNG_20130612UT.size/2,2])
Titan_1996UT = scipy.fromfile(file="f:/Astronomy/Projects/Planets/Saturn/PROJECT/TITANS.DAT", dtype=float, count=-1, sep='\t')    
Titan_1996UT=scipy.reshape(Titan_1996UT,[Titan_1996UT.size/2,2])
Titan_1996UT[:,1]=Titan_1996UT[:,1]/Titan_1996UT[:,1].max()
Titan_1996Smooth=pyasl.smooth(Titan_1996UT[:,1],9,'flat')

Titan_1996SmoothWV=Titan_1996UT
Titan_1996SmoothWV[:,1]=Titan_1996Smooth

Titan7S_1996UT = scipy.fromfile(file="f:/Astronomy/Projects/Planets/Saturn/PROJECT/TITAN7S.DAT", dtype=float, count=-1, sep='\t')    
Titan7S_1996UT=scipy.reshape(Titan7S_1996UT,[Titan7S_1996UT.size/2,2])
Titan7S_1996UT[:,1]=Titan7S_1996UT[:,1]/Titan7S_1996UT[:,1].max()
Titan7S_1996Smooth=pyasl.smooth(Titan7S_1996UT[:,1],9,'flat') #SMOOTHING IS WRONG!
Titan7S_1996SmoothWV=Titan7S_1996UT                           #WAVELENGTH SHIFT!!!
Titan7S_1996SmoothWV[:,1]=Titan7S_1996Smooth

WVOffset=-25.0
Titan_20130612UT[:,0]=Titan_20130612UT[:,0]+WVOffset
TitanRNG_20130612UT[:,0]=TitanRNG_20130612UT[:,0]+WVOffset

Titan_Karkoschka1993 = scipy.fromfile(file="f:/Astronomy/Projects/Planets/Saturn/Spectral Data/Karkoschka/1993.tab.txt", dtype=float, count=-1, sep=" ")    
Titan_Karkoschka1993=scipy.reshape(Titan_Karkoschka1993,[Titan_Karkoschka1993.size/8,8])

Titan_KarkRef1993=np.zeros((Titan_Karkoschka1993.size/8,2))
Titan_KarkRef1993[:,0]=Titan_Karkoschka1993[:,0]*10.
Titan_KarkRef1993[:,1]=Titan_Karkoschka1993[:,7]
print "***Titan_KarkRef1993***", Titan_KarkRef1993

"""
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

EWFN="TitanKarkoschkaRefEW.txt"
Key="Titan19931101050000UT"
DateTime=datetime.datetime.strptime(Key[5:9]+"-"+Key[9:11]+"-" \
        +Key[11:13]+"T"+Key[13:15]+":"+Key[15:17]+":"+Key[17:19], \
        '%Y-%m-%dT%H:%M:%S')
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 4410",4390.,4460.,20.,EWFN,False)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 4590",4560.,4620.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 4860",4750.,4900.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 5090",5010.,5150.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 5430",5290.,5500.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 5760",5670.,5820.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 5970",5880.,6040.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 6190",6040.,6320.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 6680",6460.,6770.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 6825",6780.,6900.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 7050",6920.,7120.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 7250",7130.,7450.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4/NH3? 7900",7560.,8260.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 8420",8300.,8480.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 8620",8480.,8740.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_KarkRef1993,"Titan",DateTime,"Target","CH4 8890",8740.,9130.,30.,EWFN,True)

EWFN="Titan19961101050000UT-Albedo-EW.txt"
Key="Titan19961101050000UT"
DateTime=datetime.datetime.strptime(Key[5:9]+"-"+Key[9:11]+"-" \
        +Key[11:13]+"T"+Key[13:15]+":"+Key[15:17]+":"+Key[17:19], \
        '%Y-%m-%dT%H:%M:%S')
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_1996SmoothWV,"Titan",DateTime,"Target","CH4 5430",5390.,5500.,20.,EWFN,False)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_1996SmoothWV,"Titan",DateTime,"Target","CH4 5760",5670.,5820.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_1996SmoothWV,"Titan",DateTime,"Target","CH4 5970",5880.,6040.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan_1996SmoothWV,"Titan",DateTime,"Target","CH4 6190",6040.,6320.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(Titan7S_1996SmoothWV,"Titan",DateTime,"Target","CH4 7250",7130.,7450.,30.,EWFN,True)

EWFN="Titan20130612050000UT-Albedo-EW.txt"
Key="Titan20130612050000UT"
DateTime=datetime.datetime.strptime(Key[5:9]+"-"+Key[9:11]+"-" \
        +Key[11:13]+"T"+Key[13:15]+":"+Key[15:17]+":"+Key[17:19], \
        '%Y-%m-%dT%H:%M:%S')

BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","CH4 6190",6060.,6400.,30.,EWFN,False)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","CH4 6680",6460.,6770.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","CH4 6825",6780.,6900.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","CH4 7050",6920.,7120.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","CH4 7250",7120.,7560.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","CH4/NH3? 7900",7560.,8260.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","CH4 8420",8300.,8480.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","CH4 8620",8480.,8740.,30.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","CH4 8890",8740.,9130.,30.,EWFN,True)

BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","O2 B band",6830.,6990.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","H2O 7200a + CH4",6910.,7090.,10.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","H2O 7200b + CH4",7150.,7350.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","O2 A band",7540.,7760.,20.,EWFN,True)
BandName,BandStart,BandEnd,ContWidth,EW=ComputeEW1.ComputeEW1(TitanRNG_20130612UT,"Titan",DateTime,"Target","H2O Z band + CH4",7850.,8600.,40.,EWFN,True)

pl.figure(figsize=(8.0, 2.5), dpi=150, facecolor="white")

pl.subplot(1, 1, 1)
#Plot Layout Configuration
x0=375
x1=950
xtks=24
y0=0.0
y1=0.5
ytks=6

# Set x limits
pl.xlim(x0,x1)
# Set x ticks
pl.xticks(np.linspace(x0,x1,xtks, endpoint=True))
# Set y limits
pl.ylim(y0,y1)
# Set y ticks
pl.yticks(np.linspace(y0,y1,ytks, endpoint=True))
pl.grid()
pl.tick_params(axis='both', which='major', labelsize=7)
pl.ylabel(r"$Albedo$",fontsize=7)
pl.xlabel(r"$Wavelength (nm)$",fontsize=7)
"""
pl.title("Titan Albedo - All Data",fontsize=9)
#pl.plot(Titan_20130612UT[:,0]/10.,Titan_20130612UT[:,1]*0.33,label='Titan_20130612UT',linewidth=0.5,color='b')
pl.plot(TitanRNG_20130612UT[:,0]/10.,TitanRNG_20130612UT[:,1]*0.29,label='TitanRNG_20130612UT',linewidth=1,color='b')
pl.scatter(Titan_1996UT[:,0]/10.,Titan_1996UT[:,1]*0.29,label='Titan_1996UT',marker='.',s=0.1,color='g')
pl.scatter(Titan7S_1996UT[:,0]/10.,Titan7S_1996UT[:,1]*0.40,label='Titan_1996UT',marker='.',s=0.1,color='r')

pl.plot((Titan_1996UT[:,0]-7.)/10.,Titan_1996Smooth*0.29,label='Titan_1996UT-Smth',color='g',linewidth=1.0)
pl.plot((Titan7S_1996UT[:,0]-7.)/10.,Titan7S_1996Smooth*0.40,label='Titan_1996UT-Smth',color='r',linewidth=1.0)

pl.plot(Titan_KarkRef1993[:,0]/10.,Titan_KarkRef1993[:,1],label='Karkoschka, 1993',linewidth=1,color='0.5')
"""
pl.title("Titan",fontsize=9)
Indices=np.where((TitanRNG_20130612UT[:,0] > 8000.))
TitanRNG_20130612UT[Indices,1]=np.nan

#pl.plot(Titan_20130612UT[:,0]/10.,Titan_20130612UT[:,1]*0.33,label='Titan_20130612UT',linewidth=0.5,color='b')
pl.plot(TitanRNG_20130612UT[:,0]/10.,TitanRNG_20130612UT[:,1]*0.29,label='This Work (2013) 0.200m',linewidth=1,color='b')
pl.scatter(Titan_1996UT[:,0]/10.,Titan_1996UT[:,1]*0.29,label='Titan_1996UT',marker='.',s=0.1,color='g')
pl.scatter(Titan7S_1996UT[:,0]/10.,Titan7S_1996UT[:,1]*0.40,label='Titan_1996UT',marker='.',s=0.1,color='r')

pl.plot((Titan_1996UT[:,0]-7.)/10.,Titan_1996Smooth*0.29,label='Titan_1996UT-Smth',color='g',linewidth=1.0)
pl.plot((Titan7S_1996UT[:,0]-7.)/10.,Titan7S_1996Smooth*0.40,label='Titan_1996UT-Smth',color='r',linewidth=1.0)

pl.plot(Titan_KarkRef1993[:,0]/10.,Titan_KarkRef1993[:,1],label='Karkoschka, 1994',linewidth=1,color='0.5')

LineWVs=np.array([543.0,576.0,  #H I Balmer
                  597.0,619.0,668.0,683.0,705.0,725.0,
                  790.0,842.0,862.0,889.0,
                  
                  646.0,750.0,760.0,825.0])                    #H I Balmer
                  
                  
LineY=np.array([0.05,0.05,
                0.05,0.05,0.05,0.05,0.05,0.0,
                0.05,0.0,0.2,0.1,
                
                0.05,0.05,0.05,0.05])
                                        #H I Paschen
LineLabels=[r'$CH_4$',r'$CH_4$',
            r'$CH_4$',r'$CH_4$',r'$CH_4$',r'$CH_4$',r'$CH_4$',r'$CH_4$',
            r'$CH_4$',r'$CH_4$',r'$CH_4$',r'$CH_4$',
            
            r'$NH_3$',r'$NH_3$',r'$NH_3$',r'$NH_3$',r'$NH_3$',r'$NH_3$']

for l in range(0,LineWVs.size):                
    pl.text(LineWVs[l],LineY[l],LineWVs[l].astype('|S5')+' '+LineLabels[l],fontsize=8,
            verticalalignment='bottom',horizontalalignment='center',
            rotation='vertical')

pl.legend(loc=0,ncol=3, borderaxespad=0.,prop={'size':7})
pl.subplots_adjust(left=0.06, bottom=0.20, right=0.98, top=0.88,
                wspace=None, hspace=None)
                
pylab.savefig('TitanAlbedoAllYears.png',dpi=300)