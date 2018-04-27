#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 23:58:58 2017
Shadow of a rotating black hole
@author: ashcat
"""
import numpy as np
import matplotlib.pyplot as plt
#import


# Parameters of the Metric
M = 1.               # Mass of the central body in SolarMasses
a = 0.99             # Spin Parameter (0<=a<=1)
inc = np.pi / 4         # Inclination angle in radians
rH = M + np.sqrt(M**2 - a**2)
Rshadow = np.sqrt(27 * M**2)

# Range of the parameter r
r = np.linspace(rH ,4*rH,10000000)

# Constants of motion
Delta = r**2 - 2*M*r + a**2
xi = (M*(r**2 - a**2) - r*Delta)/(a*(r - M))                    # Angular Momentum 
eta = ((r**3)*(4*M*Delta - r*(r-M)**2))/(a**2 * (r - M)**2)   # Carter Constant


# Celestial Coordinates

arg = eta + a**2 * np.cos(inc)**2 - xi**2 * (np.cos(inc)/np.sin(inc))**2
alpha = - xi / np.sin(inc)
beta = np.sqrt(arg)




# Plot of the orbit
fig, ax = plt.subplots()

# Label of the shadow plot
labelshadow = 'a=' + str(a)

#Axis information
ax.set_title('Shadow of a Rotating Black Hole')
ax.set_xlabel('X')
ax.set_ylabel('Y')

# Axis range and aspect (square)
axrange=1.5*Rshadow
ax.axis([-axrange,axrange,-axrange,axrange])
ax.set_aspect('equal', 'box')

# Grid ON
ax.grid(True, linestyle='-.')


# Plot the shadow of the rotating BH
ax.plot(alpha,beta, 'b', label=labelshadow)
ax.plot(alpha,-beta, 'b')

# Plot of the Schwarzschild's BH shadow for comparison
t = np.linspace(0, 2*np.pi, 100000)
ax.plot(Rshadow* np.cos(t), Rshadow* np.sin(t),'r--', label='Schwarzschild')

ax.legend(loc='upper right')
plt.show()
