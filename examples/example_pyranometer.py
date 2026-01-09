# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 15:31:01 2026

@author: Bkubicek
"""

#measurement equation
# E = V/S
# E = total horizontal irradiation [W/m^2]
# V = measured Voltage of pyranometer [uV]
# S = calibrated sensitivity of pyranometer [actual W/m^2 per uV]


import matplotlib.pyplot as plt

from uncertainDistribution import unDist as undi
import numpy as np

zeroOffsetBPlusMinus = 3
sensitivityPercent = 1.5/3
sensitifty = 15 #uV/(W/m^2)
dataLoggerUVoltageError = 10 #uV
nonStabilityPercent = 1
nonLinearityPercent = 0.5
spectralSelectivityPercent = 1
temperatureResponsePercent = 2
tiltResponsePercent = 0.5
maintenancePercent = 1
thermalShockAbsolut = 1 #W/m^2

samples=10000
# Distributions 
sensitivityCalibration = undi("normal", mean=1, 
                              stdDev=sensitivityPercent*0.01, maxSigma=3, samples=samples) #percent

dataLoggerVoltage = undi("rect", leftPos=-dataLoggerUVoltageError , 
                         rightPos=dataLoggerUVoltageError, samples=samples) #10uV

zeroOffsetB = undi("rect", leftPos=-zeroOffsetBPlusMinus, 
                   rightPos=zeroOffsetBPlusMinus, samples=samples) #W/m^2

nonStability = undi("rect", leftPos=1-nonStabilityPercent*0.01, 
                   rightPos=1+nonStabilityPercent*0.01, samples=samples) 

nonLinearity = undi("rect", leftPos=1-nonLinearityPercent*0.01, 
                   rightPos=1+nonLinearityPercent*0.01, samples=samples) 

spectralSelectivity = undi("rect", leftPos=1-spectralSelectivityPercent*0.01, 
                   rightPos=1+spectralSelectivityPercent*0.01, samples=samples) 

temperatureResponse = undi("rect", leftPos=1-temperatureResponsePercent*0.01, 
                   rightPos=1+temperatureResponsePercent*0.01, samples=samples) 

tiltResponse = undi("rect", leftPos=1-tiltResponsePercent*0.01, 
                   rightPos=1+tiltResponsePercent*0.01, samples=samples) 

maintenance = undi("rect", leftPos=1-maintenancePercent *0.01, 
                   rightPos=1+maintenancePercent *0.01, samples=samples) 

thermalShock = undi("rect", leftPos=-thermalShockAbsolut, 
                   rightPos=thermalShockAbsolut, samples=samples) 

DNI = 900 

sensDist = sensitivityCalibration*nonStability*nonLinearity*temperatureResponse*tiltResponse*maintenance #spectralSelectivity
idealUV = (DNI+zeroOffsetB)*sensitifty
actualUV = idealUV *sensDist + dataLoggerVoltage

wantAssurance=0.95
w= 0.5*(1-wantAssurance)
quantiles = actualUV.quantiles([w, 0.5, 1-w])
quNor = quantiles/quantiles[1]

mesUVMean = np.sum(actualUV.centers*actualUV.weights )
mesStdDev = np.sqrt(    np.sum(  ((actualUV.centers-mesUVMean)**2)*actualUV.weights ) )



print("Rel Error for %.3f confidence: -%.2f%% + %.2f%%"%(wantAssurance,100*(1-quNor[0]), 100*(quNor[2]-1)))
print("Bias: %.2f W/m^2  = %.2f%% "%(mesUVMean/sensitifty-DNI, mesUVMean/DNI/sensitifty*100-100))
seemingIrradiation = actualUV.centers/sensitifty
plt.figure()
plt.plot(seemingIrradiation,actualUV.weights, label="PDF")
plt.vlines([DNI],0,max(actualUV.weights),color="red",label="actual DNI")
plt.vlines([mesUVMean/sensitifty],0,max(actualUV.weights), color="green",label="mes DNI")
plt.vlines([quantiles[0]/sensitifty, quantiles[2]/sensitifty],0,max(actualUV.weights)*0.3, color="black",label="Confidence %.2f"%wantAssurance)
plt.vlines([(mesUVMean-i*mesStdDev)/sensitifty for i in range(-3,4)],0,max(actualUV.weights)*0.2, color="grey",label="[1-3] stdDev ")


m=mesUVMean/sensitifty
s=mesStdDev/sensitifty

xx = np.linspace(800,1000,num=200)
#0.5/np.sqrt(2*np.pi)/s
yy = max(actualUV.weights)*np.exp(-pow((xx-m)/s,2)/2)
plt.plot(xx,yy,label="real gaussian",alpha=0.5)

plt.legend()
plt.grid()
plt.xlim(800,1000)
plt.xlabel("G [W/m^2]")
plt.savefig("example_pyranometer.png")
