# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 14:07:17 2026

Example for uncertainDistribution

imagine you measure n sticks of 1m with a tool that has a +- 1% accuracy

@author: Bkubicek
"""

import matplotlib.pyplot as plt
import numpy as np


import sys
sys.path.append('../') #the uncertainDistribution lib is in the folder above
from uncertainDistribution import unDist as undi


plotDir = "./plots/"

# definition of the typical Stick
plusMinusPercent = 0.1
baseLength = 1

#define the min/max limits of a single stick measurement
leftPos = baseLength*(1-0.01*plusMinusPercent)
rightPos = baseLength*(1+0.01*plusMinusPercent)

#create the Distribution as a continous likely distribution in the interval
singleStickLength = undi("rect", leftPos=leftPos, rightPos=rightPos)

# this create the triangular distribution of the total length of two sticks back2back
manual2Sticks = singleStickLength + singleStickLength

# we want to have a list of distributions for 1,2,3,4.., 10 sticks back2back

curLength = singleStickLength

totalLengths  = [curLength]
for i in range(1,10):
    #this is an addition of two distributions objects!
    curLength = singleStickLength + curLength 
    totalLengths.append( curLength)

#  now we have an array of 10 disributions


# lets see how the quantiles change:
wantQuant = np.linspace(0.01,0.99,endpoint=True,num=10)

quantTable=[]
for i,totalDis in enumerate(totalLengths):
    curMean = totalDis.quantiles(0.5)    
    qs =  totalDis.quantiles(wantQuant)    
    quantTable.append(qs/curMean)

quantTable = np.array(quantTable)
# the table will be plotted in the last image

figsize=(10,5)
if True:
    plt.figure(figsize=figsize)
    for i,totalDis in enumerate(totalLengths):
        plt.plot(totalDis.centers, totalDis.weights, label="stick count %i"%(i+1), alpha=1)
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title("Absolute Length")
    plt.xlabel("total length [m]")
    plt.ylabel("Probability")
    plt.tight_layout()
    plt.savefig(plotDir+"example_stick_absolutLength.png")

if True:
    plt.figure(figsize=figsize)
    for i,totalDis in enumerate(totalLengths):
        curMean = totalDis.quantiles(0.5)    
        plt.plot((totalDis.centers/curMean -1)*100 , totalDis.weights, label="%i sticks"%(i+1))
    plt.grid()
    plt.title("Relative Length becoming increasingly like a normal disribution")
    plt.xlabel("percent from medium length")
    plt.ylabel("Probability")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(plotDir+"example_stick_relativeLength.png")

if True:
    plt.figure(figsize=figsize)
    
    for qnr in range(len(wantQuant)):
        plt.plot(range(len(totalLengths)), (quantTable[:,qnr]-1)*100, label="quantile %.2f"%(wantQuant[qnr]))
    plt.grid()
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.title("how it develops to normal disribution: central limit theorem")
    plt.xlabel("Average Total Length")
    plt.ylabel("Percent")
    plt.tight_layout()
    plt.savefig(plotDir+"example_stick_quantiles.png")
        