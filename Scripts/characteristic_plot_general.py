import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv

def wave_speed(A):
    if shape == 'Rectangle':
        alpha = np.arctan(0.02)
        w = 10.0  # width of the channel (m)

        return (3/2) * np.sqrt(g * w * np.sin(alpha)/f) * ((3/2)* np.sqrt(A/(2*A + w**2)) - ((2*A + w**2) ** (-3/2)))
        
    if shape == 'Wedge':
        alpha = np.arctan(0.08)
        theta = np.pi/6
        
        return (5/4) * np.sqrt((g * np.sin(alpha)/f) * np.sqrt(np.sin(theta)/8)) * (A ** (1/4))
    
    if shape == 'Semi':
        alpha = 
        theta = 
        
        return

def A0(s):
    norm = np.sqrt(2 * np.pi * (sigma ** 2))
    return (V/norm) * np.exp(-((s - mean)**2) / (sigma**2)) + A_L


def characteristic(s0, t):
    return s0 + wave_speed(A0(s0)) * t #x(t;P) = v(c)t + P


def plot_characteristics(t, s, intersections = 'False'):
    lines = []
    plt.figure(figsize=(16, 10))
    for s0 in s:
        lines.append(characteristic(s0, t))
        plt.plot(characteristic(s0, t), t, 'b-', lw=1, alpha=0.5)

    if intersections == 'True':
        intersections = []
        for i in range(np.shape(lines[:])[0] - 1):
            for j in range(np.shape(lines[:])[1] - 1):
                if np.isclose(lines[i][j], lines[i + 1][j], rtol=1e-5):
                    intersections.append((lines[i]))
                    

        df_lines = pd.DataFrame(lines)
        df_lines.to_csv(f'Flash Floods/data/{shape}_channel.csv', index = False)

        df_intersections = pd.DataFrame(intersections)
        df_intersections.to_csv(f'Flash Floods/data/{shape}_intersections.csv', index = False)

        for i in range(len(intersections)):
            plt.plot(intersections[i], t, label = f'Intersection {i + 1}')

    plt.xlabel('Distance along river, $s$ (m)', fontsize=20)
    plt.ylabel('Time, $t$ (Hours)', fontsize=20)
    plt.title(f'Characteristic Diagram for {shape} Channel', fontsize=20)
    plt.xlim(0, 5 * sigma + mean)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.savefig(f'Flash Floods/Figures/{shape}_characteristic.pdf')
    #plt.show()

##############################################################################################################

g = 9.81
f = 0.1 

A_L = 10
V = 100

mean = 1
sigma = 2 #std dev

shape = 'Rectangle'

c0 = wave_speed(A0(0))
t_shock = 2 * sigma / c0
t_max = 3 * t_shock
print(t_max)
t = np.linspace(0, t_max, 500)
s = np.linspace(-5 * sigma, 5 * sigma + mean, 100)

plot_characteristics(t, s)