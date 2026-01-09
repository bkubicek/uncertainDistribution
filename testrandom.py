# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 09:32:44 2026

@author: Bkubicek
"""

import cppyy
import numpy as np

cppyy.include("bkrandom.h")


b= cppyy.gbl.testab(int(5))

rvec =  np.array( cppyy.gbl.getNormal(0, 1, 1000) ) 

rvec2 =  np.array( cppyy.gbl.getRect(3, 3.5, 1000) ) 

cs ,ws = cppyy.gbl.getRectW(2, 3.5, 2005)
cs=np.array(cs)
ws=np.array(ws)

cs2 ,ws2 = cppyy.gbl.getNormalW(4,0.3, 3, 2000)
cs2=np.array(cs2)
ws2=np.array(ws2)

#cs3, ws3 = cppyy.gbl.addDistributions(cs,ws, cs2,ws2, max(len(cs),len(cs2)))#
#cs3, ws3 = cppyy.gbl.multiplyDistributions(cs,ws, cs2,ws2, max(len(cs),len(cs2)))
cs3, ws3 = cppyy.gbl.divideDistributions(cs,ws, cs2,ws2, int(0.5*max(len(cs),len(cs2))))
cs3=np.array(cs3)
ws3=np.array(ws3)

#pairs=np.array(pairs)
import matplotlib.pyplot as plt
plt.figure()
#plt.hist(rvec,density=True,bins=100)
#plt.hist(rvec2,density=True,bins=100)

plt.plot(cs,ws,label="rectw")
plt.plot(cs2,ws2,label="normalw")

#plt.hist(cs3,weights=ws3,label="wsum",bins=100) #density=True
plt.plot(cs3,ws3,label="sum")
plt.legend()