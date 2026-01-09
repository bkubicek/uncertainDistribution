# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 10:05:45 2026

@author: Bkubicek
"""

import numpy as np

class Distribution:
    def __init__(self, disType:str, 
                 centralValue:float = np.nan , 
                 leftValue:float = np.nan, 
                 rightValue:float = np.nan, 
                 spread:float = 1, 
                 symmetric=True ):
        self.type = disType
        self.centralValue = centralValue
        self.spread = spread
        self.leftValue= leftValue
        self.rightValue= rightValue
        
    def sample(self, count):
        if self.type=="normal":
            return np.random.normal(self.centralValue, self.spread, count)
        elif self.type=="rect":
            if ~np.isnan(self.centralValue):
                low =self.centralValue-self.spread/2
                high=self.centralValue+self.spread/2
            else:
                low =self.leftValue
                high=self.rightValue
            return np.random.rand(count)*(high-low)+low
        #elif self.type="tri":
        else:
            print("not implemented yet: %s"%self.type)
        
            
    
    
d1 = Distribution(disType="rect",centralValue=1,spread=0.5)
d2 = Distribution(disType="rect",leftValue=5,rightValue=6)

s1 = d1.sample(10000)
s2 = d2.sample(10000)

def addSamples(s1,s2):
    return np.array( [draw1+draw2 for draw1 in s1 for draw2 in s2])

s12 = addSamples(s1, s2)
import matplotlib.pyplot as plt

plt.figure()
plt.hist(s1,label="1",density=True)
plt.hist(s2,label="2",density=True)
plt.hist(s12,label="1+2",density=True)
plt.grid()


        
            