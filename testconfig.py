# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 14:40:09 2026

@author: Bkubicek
"""
import math


uncerComponents={}

uncerComponents["calibration"]=\
    {"actingOn":"V", #E,V,S
    "uncertaiType":"relative", # relative,absolute
    "distribution":"normal", # normal, rectangular, triangular, uShaped
    "symmetric":True, #True False
    "correlated":True, # True False
    "valA":1.5,
    "baseUnit":"%",
    "factorUnit":1,   # mV=1e-3, muV=1e-6,...
    "analyze":True
    }

uncerComponents["datalogger"]=\
    {"actingOn":"V", #E,V,S
    "uncertaiType":"relative", # relative,absolute
    "distribution":"normal", # normal, rectangular, triangular, uShaped
    "symmetric":True, #True False
    "correlated":True, # True False
    "valA":0,
    "baseUnit":"uV",
    "factorUnit":1,   # mV=1e-3, muV=1e-6,...
    "analyze":True
    }

uncerComponents["zeroOffsetA"]=\
    {"actingOn":"V", #E,V,S
    "uncertaiType":"relative", # relative,absolute
    "distribution":"normal", # normal, rectangular, triangular, uShaped
    "symmetric":True, #True False
    "correlated":True, # True False
    "valA":7,
    "baseUnit":"W/m^2",
    "factorUnit":1,   # mV=1e-3, muV=1e-6,...
    "analyze":True
    }    

uncerComponents["zeroOffsetB"]=\
    {"actingOn":"V", #E,V,S
    "uncertaiType":"relative", # relative,absolute
    "distribution":"normal", # normal, rectangular, triangular, uShaped
    "symmetric":True, #True False
    "correlated":True, # True False
    "valA":2,
    "baseUnit":"W/m^2",
    "factorUnit":1,   # mV=1e-3, muV=1e-6,...
    "analyze":True
    }    
    
uFactors = {}
uFactors["normal"] = 1/3
uFactors["rectangular"] = 1/math.sqrt(3)
uFactors["triangular"] = 1/math.sqrt(6)
uFactors["uShaped"] = 1/math.sqrt(2)


def calculateU( component:dict):
    res = None
    compType = component["mesType"]
    res = uFactors[  MeasTypes[ compType ]  ]
    if not component["symmetry"]:
        res *= 0.5
    component["valU"]=res*component["valA"]
    
        
