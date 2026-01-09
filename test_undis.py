# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 13:04:19 2026

@author: Bkubicek
"""

import matplotlib.pyplot as plt
from uncertainDistribution import unDist as undi

dis1 = undi("normal", mean=1, stdDev= 0.1, maxSigma=3)
dis2 = undi("rect", leftPos=6, rightPos=7)
dis3 = undi("tri", leftPos=10, centerPos=11,rightPos=13)

sum23 = dis2+dis3

sum2Skalar = dis2+10
sumSkalar2 = 15+dis2

mul23 = dis2*dis3
mulSkal = dis3*3

div23 = dis2/dis3
divSkal = dis3/3
divSkal2 = 1/dis3


plt.figure()
plt.plot(dis1.centers, dis1.weights, label="dis1")
plt.plot(dis2.centers, dis2.weights, label="dis2")
plt.plot(dis3.centers, dis3.weights, label="dis3")

plt.plot(sum23.centers, sum23.weights, label="sum23")
plt.plot(sum2Skalar.centers, sum2Skalar.weights, label="sum2Skalar")
plt.plot(sumSkalar2.centers, sumSkalar2.weights, label="sumSkalar2")

plt.plot(mul23.centers, mul23.weights, label="mul23")
plt.plot(mulSkal.centers, mulSkal.weights, label="mulSkal")

plt.plot(div23.centers, div23.weights, label="div23")
plt.plot(divSkal.centers, divSkal.weights, label="divSkal")
plt.plot(divSkal2.centers, divSkal2.weights, label="divSkal2")

 
plt.grid()
plt.legend()

debug= dis3.weights
