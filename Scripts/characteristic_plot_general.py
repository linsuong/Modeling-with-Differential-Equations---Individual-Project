import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv

def wave_speed(A):
    if shape == 'Rectangle':
        alpha = np.pi - np.arctan(0.02)
        w = 20.0  # width of the channel (m)

        l = w + ((2 * A)/w)
        l_prime = 2/w
        
        #return (3/2) * np.sqrt(g * w * np.sin(alpha)/f) * ((3/2)* np.sqrt(A/(2*A + w**2)) - (((A)/((2*A + w**2)) ** (3/2))))
 
    if shape == 'Wedge':
        alpha = np.arctan(0.08)
        theta = np.pi/6
        
        l = np.sqrt((8 * A)/np.sin(theta))
        l_prime = (np.sqrt(8 * A * np.sin(theta)))/(2 * np.sin(theta))
        
       #return (5/4) * np.sqrt((g * np.sin(alpha)/f) * np.sqrt(np.sin(theta)/8)) * (A ** (1/4))
    
    if shape == 'Semi':
        alpha = np.arctan(0.08)
        theta = np.pi/3
        
        l = np.sqrt((2 * A)/(theta - np.sin(theta))) * theta
        l_prime = (theta/(theta - np.sin(theta))) * np.sqrt((2 * A)/(theta - np.sin(theta)))
    
    if shape == 'Parabola':
        a = 5
        w = 10
        alpha = np.arctan(0.1)
        
        var_0 = np.sqrt(1 + (3 * A * a * w))
        var_1 = (3 * A)/ (w **2)
        var_2 = (3 * A * (w **2))/2
        
        l = var_2 * (var_1 * var_0 + np.log(1 + var_0 + var_1)) 
        l_prime =((3 * (w **2) * l)/2) + var_2 * ((3/(w ** 2)) * var_0 + (((9 * A * a)/(2 * w))/var_0) + ((var_0 + var_1)/(((3 * a * w/2)/(var_0)) + (3/(w**2)))))
    
    ratio = A/l
    
    return (1/2) * np.sqrt((g * np.sin(alpha))/f) *( (3 * np.sqrt(ratio)) - ((ratio) ** (3/2)) * l_prime)

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
    plt.ylabel('Time, $t$ (Seconds)', fontsize=20)
    plt.title(f'Characteristic Diagram for {shape} Channel', fontsize=20)
    plt.xlim(0, 5 * sigma + mean)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tick_params(axis='both', which='major', labelsize=20)
    #plt.savefig(f'Figures/{shape}_characteristic.pdf')
    plt.show()

##############################################################################################################

if __name__ == '__main__':
    g = 9.81
    f = 0.05

    A_L = 10
    V = 5000

    mean = 0
    sigma = 100 #std dev

    shape = 'Parabola'

    c0 = wave_speed(A0(0))
    t_shock = 2 * sigma / c0
    t_max = 3 * t_shock
    print(t_max)
    t = np.linspace(0, t_max, 500)
    s = np.linspace(-5 * sigma, 5 * sigma + mean, 100)

    plot_characteristics(t, s)