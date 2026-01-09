# -*- coding: utf-8 -*-
"""
Python Library for performing calculations of distributions/error distributions

Motivation TL/DR:
    You have multiple sticks, you know how accurate they are measured.
    The accuracy is defined e.g. by a rectangular "+- 1 percent" window, 
    or a normal distribution.
    This library allows to define that accuracy types, and add multiple sticks
    to get the estimated distribution for the likely total length of the combined sticks.
    
Features:
    One can define multiple distributions:
        Normal Distribution
        Rectangular Distribution
        Triangular Distribution (possibly assymetric)
    
    one can:
        add
        multiply
        divide (if the second distribution is not zero-crossing)
    
    this allows in combination to perform estimations of the probablility 
    distribution of complex systems: 
        e.g. to proofe the measurement accuracy of an accredited measurement
            observable = bias_error + sensitivity*set_point/current/length
            where each quantity on the right hand side has individual precission intervals/types
            the bias_error might be normal distributed, while
            the other things are individually rectangular distributed
                
Technicals:
    One could do all that using monte carlo sampling and be done.
    However it would be slow.
    Instead a single distribution is defined by two equally sized arrays
    the first one, "centers", is the value of the bin center of the observable, e.g. the stick length
    the second one, "weights", is the likelyhood that the corresponding interval 
     would be hit when measuring the random stick with a perfect device.
     
    For rect and tri distributions, the centers are naturally limited
    For a normal distribution one has to choose how much of the tails are resolved
    
    when adding/multiplying/dividing two distributions, all combination of centers are looked at
    and the combined likelyhood is defining the new resulting weights.
    
    The Class uses operator overloading, so that + * / with a distributino/scalar are possible

Usage:
    One can define example distributinos:
        dis1 = undi("normal", mean=1, stdDev= 0.1, maxSigma=3)
        dis2 = undi("rect", leftPos=6, rightPos=7)
        dis3 = undi("tri", leftPos=10, centerPos=11,rightPos=13)
        
    one can look at the distribution function of three added samples of dis2:
        res = dis2+dis2+dis2
    a single dis2 is a rect function, two samples would result in a triangle disribution
    with more and more samples, the total distribution approaches a normal distribution
    see: Central Limit Theorem 
    
    to plot, one can obtain the bin centers as dis.centers, and the weights as dis.weights
    

core functionally is written in c++ for optimization
and automatically wrapped to python from "bkUncDist.hpp" using cppyy "https://cppyy.readthedocs.io/en/latest/"

License: LGPL-3.0-or-later

@author: Bkubicek
"""
import numpy as np
import cppyy
import copy

try:
    cppyy.add_include_path('./include')
except:
    cppyy.add_include_path('../include')
    
cppyy.include("bkUncDist.hpp")


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
        
        # if called with actual arguments, calculated the distribution
        if disType!="calculated":
            self.sample()
            
    def sample(self): 
        """
        create distribution, calls more specific functions

        """
        if self.disType=="normal":
            self.sampleNormal()
        elif self.disType=="rect":
            self.sampleRect()
        elif self.disType=="tri":
            self.sampleTri()
            
    def sampleNormal(self): 
        # wraps the c++ function, and converts the results to numpy 
        centers, weights= cppyy.gbl.getNormal(self.mean,self.stdDev, self.maxSigma, self.samples)
        self.centers = np.array(centers)
        self.weights=np.array(weights)
        
    def sampleRect(self):
        # wraps the c++ function, and converts the results to numpy 
        centers, weights= cppyy.gbl.getRect(self.leftPos, self.rightPos, self.samples)
        self.centers = np.array(centers)
        self.weights=np.array(weights)
        
    def sampleTri(self):
        # wraps the c++ function, and converts the results to numpy 
        centers, weights= cppyy.gbl.getTri(self.leftPos, self.centerPos, self.rightPos, self.samples)
        self.centers = np.array(centers)
        self.weights=np.array(weights)
    
    def quantiles(self,vals):
        # obtain the quantiles of the distibution
        # vals can be a single value, or a list/np.array
        # does not do anti-aliasing
        
        # this would work in the curren numpy unstable, sadly not for everyone yet
        # hence the code below is necessary
        # return np.quantile(self.centers, vals, weights=self.weights)
        
        # create cummulative distribution function
        # requiring that sum of weigths=1
        cum = np.zeros_like(self.weights)
        for i in range (1, len(self.weights)):
            cum[i]=cum[i-1]+self.weights[i]
        
        poses = np.searchsorted(cum, vals)
        
        #split based on if single of list is passed in vals
        if isinstance(vals, list) or isinstance(vals, np.ndarray):
            return np.array( [ self.centers[p] for p in poses])
        else:
            return  self.centers[int(poses)] 
    def getMeanStd(self):
        mean = np.sum(self.centers*self.weights )
        stdDev = np.sqrt(    np.sum(  ((self.centers-mean)**2)*self.weights ) )
        return mean, stdDev
        
    def __add__(self, x):
        # operater that is called if python sees code "unDist+x" 
        # will be interpred
        #  * as shift the centers, if x=float
        #  * create the combined distrbution of sample(self)+sample(x)
        
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
            resDis.centers , resDis.weights = cppyy.gbl.operateDistributionsAndResample(cppyy.gbl.OpAdd, self.centers,self.weights, 
                                                  x.centers, x.weights, newSamples)
            resDis.centers = np.array(resDis.centers)
            resDis.weights = np.array(resDis.weights)
            return resDis
        else:
            raise Exception("not sure how to add type %s to an unDistributino"%type(x))
        
    def __radd__(self, x):
        # operater that is called if python sees code "x+unDist" 
        # + is assumed commutative
        return self.__add__(x)
    
    def __sub__(self, x):
        # operater that is called if python sees code "unDist-x" 
        return self.__add__(-x)
    
    def __rsub__(self, x):
        # operater that is called if python sees code "x-unDist" 
        print("rsub")
        nd = unDist()
        nd.centers = -copy.copy(self.centers)
        nd.weights = copy.copy(self.weights)
        
        return nd.__add__(-x)
    
    
                
    def __mul__(self, x):
        # operater that is called if python sees code "unDist*x" 
        # will be interpred
        #  * as spreading of the the centers by factor x, if x=float
        #  * create the combined distrbution of sample(self)*sample(x)
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
            resDis.centers , resDis.weights = cppyy.gbl.operateDistributionsAndResample(cppyy.gbl.OpMul, self.centers,self.weights, 
                                                  x.centers, x.weights, newSamples)
            resDis.centers = np.array(resDis.centers)
            resDis.weights = np.array(resDis.weights)
            return resDis
        else:
            raise Exception("not sure how to multiply type %s to an unDistributino"%type(x))
        
    def __rmul__(self, x):
        # operater that is called if python sees code "x*unDist" 
        # * is assumed commutative
        return self.__mul__(x)

    def __truediv__(self, x):
        # operater that is called if python sees code "unDist/x" 
        # will be interpred
        #  * as divition of the the centers by factor x, if x=float
        #  * create the combined distrbution of sample(self)/sample(x)
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
            resDis.centers , resDis.weights = cppyy.gbl.operateDistributionsAndResample(cppyy.gbl.OpDiv, self.centers,self.weights, 
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
            resDis.centers , resDis.weights = cppyy.gbl.operateDistributionsAndResample(cppyy.gbl.OpDiv, 
                    x.centers, x.weights,  self.centers, self.weights, newSamples)
            resDis.centers = np.array(resDis.centers)
            resDis.weights = np.array(resDis.weights)
            return resDis
        else:
            raise Exception("not sure how to divide type %s to an unDistributino"%type(x))
      
    