# -*- coding: utf-8 -*-
"""
Created on Thu Jan  8 12:37:53 2026

@author: Bkubicek
"""
import numpy as np
import cppyy


cppyy.include("bkrandom.h")


class unDist:
    def __init__(self, disType="calculated", mean=np.nan, stdDev=np.nan, maxSigma=np.nan, 
                  leftPos=np.nan, centerPos=np.nan, rightPos=np.nan, samples=1001):
        self.disType = disType
        self.samples = samples
        
        self.centers = None
        self.weights = None
        
        if disType=="normal":
            if np.isnan(mean):
                raise Exception("Normal distribution needs mean")
            if np.isnan(stdDev):
                raise Exception("Normal distribution needs stdDev")
            if np.isnan(maxSigma):
                raise Exception("Normal distribution needs maxSigma, to know how much possibility to resolve")
            self.mean = mean
            self.stdDev = stdDev
            self.maxSigma = maxSigma
        elif disType=="rect":
            if np.isnan(leftPos):
                raise Exception("Rect distribution needs leftPos")
            if np.isnan(rightPos):
                raise Exception("Normal distribution needs rightPos")
            self.leftPos = leftPos
            self.rightPos = rightPos
        elif disType=="tri":
            if np.isnan(leftPos):
                raise Exception("Rect distribution needs leftPos")
            if np.isnan(centerPos):
                raise Exception("Normal distribution needs centerPos")
            if np.isnan(rightPos):
                raise Exception("Normal distribution needs rightPos")
            self.leftPos = leftPos
            self.centerPos = centerPos
            self.rightPos = rightPos
        elif disType=="calculated":
            pass
        else:
            raise Exception("Unimplemented distribution type %s"%disType)
        
        if disType!="calculated":
            self.sample()
            
    def sample(self):
        if self.disType=="normal":
            self.sampleNormal()
        elif self.disType=="rect":
            self.sampleRect()
        elif self.disType=="tri":
            self.sampleTri()
            
    def sampleNormal(self):
        centers, weights= cppyy.gbl.getNormalW(self.mean,self.stdDev, self.maxSigma, self.samples)
        self.centers = np.array(centers)
        self.weights=np.array(weights)
    def sampleRect(self):
        centers, weights= cppyy.gbl.getRectW(self.leftPos, self.rightPos, self.samples)
        self.centers = np.array(centers)
        self.weights=np.array(weights)
    def sampleTri(self):
        centers, weights= cppyy.gbl.getTriW(self.leftPos, self.centerPos, self.rightPos, self.samples)
        self.centers = np.array(centers)
        self.weights=np.array(weights)
    
    def quantiles(self,vals):
        # this would work in the curren np unstable, sadly not for everyone yet
        #return np.quantile(self.centers, vals, weights=self.weights)
        
        cum = np.zeros_like(self.weights)
        for i in range (1, len(self.weights)):
            cum[i]=cum[i-1]+self.weights[i]
        poses = np.searchsorted(cum, vals)
        
        
        
        if isinstance(vals, list) or isinstance(vals, np.ndarray):
            return np.array( [ self.centers[p] for p in poses])
        else:
            return  self.centers[int(poses)] 
        
    def __add__(self, x):
        if self.centers is None:
            raise Exception("Cannot add unsampled distributions, left one")
        if isinstance(x, float) or isinstance(x, int):
            #self.centers+=x
            resDis = unDist()
            resDis.samples = self.samples
            resDis.weights = self.weights
            resDis.centers = self.centers + x
            return resDis
        elif isinstance(x, unDist):
            if x.centers is None:
                raise Exception("Cannot add unsampled distributions, right one")
            
            newSamples = min(self.samples, x.samples) +1
            resDis = unDist()
            resDis.samples = newSamples
            resDis.centers , resDis.weights = cppyy.gbl.addDistributions(self.centers,self.weights, 
                                                  x.centers, x.weights, newSamples)
            resDis.centers = np.array(resDis.centers)
            resDis.weights = np.array(resDis.weights)
            return resDis
        else:
            raise Exception("not sure how to add type %s to an unDistributino"%type(x))
        
    def __radd__(self, x):
        return self.__add__(x)
                
    def __mul__(self, x):
        if self.centers is None:
            raise Exception("Cannot multiply unsampled distributions, left one")
        if isinstance(x, float) or isinstance(x, int):
            #self.centers+=x
            resDis = unDist()
            resDis.samples = self.samples
            resDis.weights = self.weights
            resDis.centers = self.centers * x
            return resDis
        elif isinstance(x, unDist):
            if x.centers is None:
                raise Exception("Cannot multiply unsampled distributions, right one")
            
            newSamples = min(self.samples, x.samples)+1
            resDis = unDist()
            resDis.samples = newSamples
            resDis.centers , resDis.weights = cppyy.gbl.multiplyDistributions(self.centers,self.weights, 
                                                  x.centers, x.weights, newSamples)
            resDis.centers = np.array(resDis.centers)
            resDis.weights = np.array(resDis.weights)
            return resDis
        else:
            raise Exception("not sure how to multiply type %s to an unDistributino"%type(x))
        
    def __rmul__(self, x):
        return self.__mul__(x)

    def __truediv__(self, x):
        if self.centers is None:
            raise Exception("Cannot divide unsampled distributions, left one")
        if isinstance(x, float) or isinstance(x, int):
            #self.centers+=x
            resDis = unDist()
            resDis.samples = self.samples
            resDis.weights = self.weights
            resDis.centers = self.centers / x
            return resDis
        elif isinstance(x, unDist):
            if x.centers is None:
                raise Exception("Cannot divide unsampled distributions, right one")
            
            newSamples = min(self.samples, x.samples)+1
            resDis = unDist()
            resDis.samples = newSamples
            resDis.centers , resDis.weights = cppyy.gbl.divideDistributions(self.centers,self.weights, 
                                                  x.centers, x.weights, newSamples)
            resDis.centers = np.array(resDis.centers)
            resDis.weights = np.array(resDis.weights)
            return resDis
        else:
            raise Exception("not sure how to divide type %s to an unDistributino"%type(x))
    def __rtruediv__(self, x):
        if self.centers is None:
            raise Exception("Cannot divide unsampled distributions, left one")
        if isinstance(x, float) or isinstance(x, int):
            #self.centers+=x
            resDis = unDist()
            resDis.samples = self.samples
            resDis.weights = self.weights
            resDis.centers = x/self.centers 
            return resDis
        elif isinstance(x, unDist):
            if x.centers is None:
                raise Exception("Cannot divide unsampled distributions, right one")
            
            newSamples = min(self.samples, x.samples)+1
            resDis = unDist()
            resDis.samples = newSamples
            resDis.centers , resDis.weights = cppyy.gbl.divideDistributions(x.centers, x.weights, 
                                                self.centers,self.weights, newSamples)
            resDis.centers = np.array(resDis.centers)
            resDis.weights = np.array(resDis.weights)
            return resDis
        else:
            raise Exception("not sure how to divide type %s to an unDistributino"%type(x))
      
    