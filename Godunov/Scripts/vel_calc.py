import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.interpolate as interp
#from characteristic_plot_general import wave_speed

# Define constants

g = 9.81
f = 0.05 
A_max = 5000
L = 100 

def l(A):
    if shape == 'Rectangle':
        w = 10
        return w ** 2 + (2 * A) / w
    
    if shape == 'Wedge':
        theta = np.pi / 6
        return np.sqrt((8 * A) / (np.sin(theta)))
    
    if shape == 'Semi':
        theta = np.pi/3
        
        return np.sqrt((2 * A) / (theta - np.sin(theta))) * theta

def u_bar(A):
    if shape == 'Rectangle':
        alpha = np.arctan(0.02)
        
    if shape == 'Wedge':
        alpha = np.arctan(0.065)
        
    if shape == 'Semi':
        alpha = np.arctan(0.3)
    
    if shape == 'Parabola':
        alpha = np.arctan(0.1)
        
    return np.sqrt((g * np.sin(alpha) * A) / (f * l(A)))

def Q(A):
    return A * u_bar(A)

def A0(s):
    b = 0  
    sigma = L   
    norm = np.sqrt(2 * np.pi * (sigma ** 2))  
    AL = 3  
    return AL + (A_max/norm) * np.exp(-((s - b)**2) / (sigma**2))


def shock_velocity(A0, A1):
    return ((Q(A1) - Q(A0))/(A1 - A0))

shape = 'Rectangle'
print('rect', shock_velocity(19.2, 2.65))

shape = 'Wedge'
print('wedge', shock_velocity(19.26, 5.39))

shape = 'Semi'
print('semi', shock_velocity(19.22, 5.04))