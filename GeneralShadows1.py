#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 23:58:58 2017
Shadow of a rotating black hole
@author: ashcat
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

class Shadow(object):
    def __init__(self, M = 1. ,  a = 0.5, inc = np.pi/2):
        
        # Parameters of the Metric
        self.M = M          # Mass of the central body in SolarMasses
        self.a = a          # Spin Parameter (0<=a<=1)
        self.inc = inc      # Inclination angle in radians
        
        # Event Horizon
        self.rH = self.M + np.sqrt(self.M**2 - self.a**2) 
        
        # Range of the parameter r
        self.r = np.linspace(self.rH ,4*self.rH,10000000)
        
        # Constants of motion
        
        self.Delta = self.r**2 - 2*self.M*self.r + self.a**2
        # Angular Momentum
        self.xi = (self.M*(self.r**2 - self.a**2) 
                - self.r*self.Delta)/(self.a*(self.r - self.M))  
        # Carter Constant
        self.eta = ((self.r**3)*(4*self.M*self.Delta 
            - self.r*(self.r-self.M)**2))/(self.a**2 * (self.r - self.M)**2)   
        
        # Radius of the corresponding Schwarzschild BH
        self.Rcircle = np.sqrt(27 * self.M**2)
        
        #Celestial coordinates
        self.alpha = - self.xi / np.sin(self.inc)
        
        self.arg = self.eta + self.a**2 * np.cos(self.inc)**2 
        - self.xi**2 * (np.cos(self.inc)/np.sin(self.inc))**2
        
        self.beta = np.sqrt(self.arg)
        
        '''def f(r):
            arg1 =self.eta + self.a**2 * np.cos(self.inc)**2 
            - self.xi**2 * (np.cos(self.inc)/np.sin(self.inc))**2
            return arg1'''
        
        #print(f(5))

M =1.
a = [0.3, 0.6, 0.9, 0.99]

# Radius of the corresponding Schwarzschild BH
Rcircle = np.sqrt(27 * M**2)
incl = np.pi/2
print(incl)

# Plot of the orbit
fig, ax = plt.subplots()

#Axis information
ax.set_title('Shadow of a Rotating Black Hole')
ax.set_xlabel('X')
ax.set_ylabel('Y')

# Axis range and aspect (square)
axrange=1.5*Rcircle
ax.axis([-axrange,axrange,-axrange,axrange])
ax.set_aspect('equal', 'box')

# Grid ON
ax.grid(True, linestyle='-.')

# Plot of the Schwarzschild's BH shadow for comparison
t = np.linspace(0, 2*np.pi, 100000)
ax.plot(Rcircle * np.cos(t), Rcircle * np.sin(t),'k--', label='a = 0')

labels=[]
colors=[]
bh =[]
for i in range(0, len(a)):
    # Calculation of the Shadow
    bh.append(Shadow(M , a=a[i]))
    # Label of the shadow plot
    labels.append('a='+ str(bh[i].a))
    # Colors for the plots
    colors.append('C'+str(i+1)+'-')

    # Plot the shadow of the rotating BH
    ax.plot(bh[i].alpha,bh[i].beta, colors[i], label=labels[i])
    ax.plot(bh[i].alpha,-bh[i].beta, colors[i])

ax.legend(loc='upper right')
plt.show()
