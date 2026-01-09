# -*- coding: utf-8 -*-
"""
Created on Fri Jan  9 13:46:42 2026

@author: Bkubicek
"""

import matplotlib.pyplot as plt
import numpy as np


import sys
sys.path.append('../') #the uncertainDistribution lib is in the folder above
from uncertainDistribution import unDist as undi


plotDir = "./plots/"

"""
we have a an electrical conductor
we put a measurement shunt resistor in series with it
and connect the plus/minus of the resistor to a voltmeter
we want to estimate the accuracy of the total measurement

The resistor has R=1 Ohm, and a relative accuracy of 0.5 %
it has a temperature coefficient tcR of 0.03%/per degree 
the accuracy of the temperature coefficient is +- 0.001 percentPoint

The temperature measurement is accurate to +-0.5°

the voltmeter has a random noise background with a sigma of 0.01 Volt
and otherwise is calibrated to 0.1 % of its output reading

The total measurement Equation is:
    I = [U+ voltageBackground] / [  R*(1+0.01*tcR*deltaT) ] 
    
lets define the distributions
"""

#R = undi("rect", leftPos=0.99, rightPos=1.01) #Ohm
R = undi("rect", leftPos=0.995, rightPos=1.005) #Ohm
tcR = undi("rect", leftPos=0.03-0.001, rightPos=0.03+0.001)
deltaT= undi("rect", leftPos=25-0.5, rightPos=25+0.5)-25
voltageBackground = undi("normal", mean=0, stdDev= 0.01, maxSigma=5)
voltageRelAcc = undi("rect", leftPos=0.999, rightPos=1.001)

"""  
lets throw all together:
"""

actualCurrent = 10 # [A]
perfectVoltage = 1*actualCurrent 
volDistr = (perfectVoltage+voltageBackground)*voltageRelAcc
Rdis = R*(1+0.01*tcR*deltaT)

I = volDistr /Rdis

"""  
lets evaluate the quantiles of the distribution for intended confidences:
Output:
  Rel Error for 0.950 confidence: -0.5487% + 0.5531%
  Rel Error for 0.990 confidence: -0.6509% + 0.6553%
  Rel Error for 0.999 confidence: -0.7530% + 0.7597%
General Bias: 0.00020 A  = 0.0020% of nominal
"""

wantConfidences=[0.95, 0.99, 0.999]
for wantAssurance in wantConfidences:
    #print("Confidence requested: %.4f"%wantAssurance)
    w= 0.5*(1-wantAssurance)
    quantiles = I.quantiles([w, 0.5, 1-w])
    quNor = quantiles/quantiles[1]
    
    mesMean, mesStdDev = I.getMeanStd()
    
    print("  Rel Error for %.3f confidence: -%.4f%% + %.4f%%"%(wantAssurance,100*(1-quNor[0]), 100*(quNor[2]-1)))

print("General Bias: %.5f A  = %.4f%% of nominal"%(mesMean-actualCurrent, mesMean/actualCurrent*100-100))


"""
and now lets plot how the distribution looks

"""
plt.figure(figsize=(11,10),dpi=200)

ymax =  max(I.weights)
plt.plot(I.centers,I.weights, label="I measured Distribution")
plt.vlines([actualCurrent],0,ymax,color="red",label="actual current")
plt.vlines([mesMean],0,ymax, color="green",label="mes most likely")
plt.vlines([quantiles[0], quantiles[2]],0,ymax*0.3, color="black",label="Confidence %.5f"%wantAssurance)
plt.vlines([(mesMean-i*mesStdDev) for i in range(-3,4)],0,ymax*0.2, color="grey",label="[1-3] stdDev ")
plt.grid()
plt.legend(bbox_to_anchor=(1, 0.5))
plt.xlabel("Current [A]")
plt.ylabel("likelyhood of measurement")
plt.tight_layout()
plt.savefig(plotDir+"example_currentMeasurement.png")
