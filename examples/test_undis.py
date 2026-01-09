# -*- coding: utf-8 -*-
"""
This demonstrates how to use uncertainDistribution 

@author: Bkubicek
"""

import matplotlib.pyplot as plt
import sys
sys.path.append('../') #the uncertainDistribution lib is in the folder above
from uncertainDistribution import unDist as undi

plotDir = "./plots/"


# we define three different distributions
disN = undi("normal", mean=1, stdDev= 0.1, maxSigma=5)
disR = undi("rect", leftPos=6, rightPos=7)
disT = undi("tri", leftPos=10, centerPos=11,rightPos=13) #note its asymmetric

figsize=(10,5)
if True:
    plt.figure(figsize=(20,20))
    plt.scatter(disT.centers, disT.weights, label="disT",s=0.1)
    plt.grid()
    plt.legend()
    plt.savefig(plotDir+"testUndis_tri_closeup.png")

if True:
    plt.figure(figsize=figsize)
    plt.plot(disN.centers, disN.weights, label="disN")
    plt.plot(disR.centers, disR.weights, label="disR")
    plt.plot(disT.centers, disT.weights, label="disT")
    plt.grid()
    plt.legend()
    plt.savefig(plotDir+"testUndis_basic_Distributions.png")


# now lets get the sum of two distributions
# this is the likelyhood of the sum of random picks from both distributions
# this is not an arithmetic plus but a complex algorithm 
sumExample1 = disR + disT
sumExample2 = disR + disN

if True:
    fig, (ax1, ax2) = plt.subplots(nrows=2,figsize=figsize)

    ax1.plot(disR.centers, disR.weights, label="disR", alpha=0.5)
    ax1.plot(disT.centers, disT.weights, label="disT", alpha=0.5)
    ax1.plot(sumExample1.centers, sumExample1.weights, label="Sum rect+tri")    
    
    ax1.grid()
    ax1.legend(bbox_to_anchor=(1, 0.5))
    
    
    ax2.plot(disN.centers, disN.weights, label="disN", alpha=0.5)
    ax2.plot(disR.centers, disR.weights, label="disR", alpha=0.5)
        
    ax2.plot(sumExample2.centers, sumExample2.weights, label="Sum rect+normal")
    
    ax2.grid()
    ax2.legend(bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.savefig(plotDir+"testUndis_sum.png")

# we can shift a distribution also by adding a float
sumSkalar = disR+10

if True:
    plt.figure(figsize=figsize)
    
    plt.plot(disR.centers, disR.weights, label="disR")
    plt.plot(sumSkalar.centers, sumSkalar.weights, label="disR+10")
    plt.grid()
    plt.legend()
    plt.savefig(plotDir+"testUndis_skalarAddition.png")


# we can multiply distributions
mulRT = disR*disT
mulSkalarR = disR*3

if True:
    plt.figure(figsize=figsize)    
    plt.plot(disR.centers, disR.weights, label="R")
    plt.plot(disT.centers, disT.weights, label="T")
    plt.plot(mulRT.centers, mulRT.weights, label="R * T")
    plt.plot(mulSkalarR.centers, mulSkalarR.weights, label="R * 3")
    plt.grid()
    plt.legend()
    plt.savefig(plotDir+"testUndis_mul.png")
    
    
# we can divide distributions
divRT = disR/disT
divTSkal = disT/3
invT = 1/disT

if True:
    fig, (ax1, ax2) = plt.subplots(nrows=2,figsize=(10,10))
    
    ax1.plot(disR.centers, disR.weights, label="R")
    ax1.plot(disT.centers, disT.weights, label="T")
    ax1.plot(divRT.centers, divRT.weights, label="R / T")
    ax2.plot(invT.centers, invT.weights, label="1/T")
    ax1.grid()
    ax2.grid()
    ax1.legend()
    ax2.legend()
    plt.savefig(plotDir+"testUndis_div.png")



