import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv



def wave_speed(A):
    g = 9.81
    f = 0.05

    alpha = np.arctan(0.02)
    w = 20.0  # width of the channel (m)

    l = w + ((2 * A) / w)
    l_prime = 2/w
        
    #return (3/2) * np.sqrt(g * w * np.sin(alpha)/f) * ((3/2)* np.sqrt(A/(2*A + w**2)) - (((A)/((2*A + w**2)) ** (3/2))))

    #def wave_speed(A):
    if shape == 'Wedge':
        alpha = np.arctan(0.08)
        theta = np.pi/3
        
        l = np.sqrt((8 * A)/np.sin(theta))
        l_prime = 2 / (np.sqrt(2 * A * np.sin(theta)))
        
    #return (5/4) * np.sqrt((g * np.sin(alpha)/f) * np.sqrt(np.sin(theta)/8)) * (A ** (1/4))
    
    if shape == 'Semi':
        alpha = np.arctan(0.08)
        theta = np.pi/3
        
        l = np.sqrt((2 * A)/(theta - np.sin(theta))) * theta
        l_prime = theta/(np.sqrt(2 * A * (theta - np.sin(theta))))
    
    if shape == 'Parabola':
        w = 5
        
        var_0 = (3 * A) / (2 * (w ** 2))

        # Calculate l
        l = ((2 * (w ** 3))/(3 * A)) * (var_0 * np.sqrt(1 + (var_0 ** 2)) + np.log(np.abs(np.sqrt(1 + (var_0 ** 2)) + var_0))) 
        
        # Calculate l'
        alpha = (2 * (w ** 3)) / (3 * A)
        alpha_prime = -((2 * (w ** 3)) / (3 * (A ** 2)))

        beta_prime = ((1 + 2 * (var_0 ** 2)) / np.sqrt(1 + (var_0 ** 2))) * (3 / (2 * (w ** 2)))
        gamma_prime = (9 * A) / ((4 * (w ** 4)) * np.sqrt(1 + (var_0 ** 2)))
        
        gamma = np.sqrt(1 + (var_0 ** 2)) + var_0

        l_prime = (alpha_prime * l) + alpha * (beta_prime + (gamma_prime / gamma))
        
    ratio = A/l
    
    return (1/2) * np.sqrt((g * np.sin(alpha))/f) *((3 * np.sqrt(ratio)) - (((ratio) ** (3/2)) * l_prime))

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
        df_lines.to_csv(f'data/{shape}_channel.csv', index = False)

        df_intersections = pd.DataFrame(intersections)
        df_intersections.to_csv(f'data/{shape}_intersections.csv', index = False)

        for i in range(len(intersections)):
            plt.plot(intersections[i], t, label = f'Intersection {i + 1}')

    plt.xlabel('Distance along river, $s$ (m)', fontsize=40)
    plt.ylabel('Time, $t$ (Seconds)', fontsize=40)
    plt.title(f'Characteristic Diagram for {shape} Channel', fontsize=40)
    plt.xlim(0, 5 * sigma + mean)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tick_params(axis='both', which='major', labelsize=40)
    plt.savefig(f'Figures/{shape}_characteristic.pdf')
    plt.show()

##############################################################################################################

if __name__ == '__main__':
    g = 9.81
    f = 0.05

    A_L = 3
    V = 5000

    mean = 0
    sigma = 100 #std dev

    shape = 'Wedge'

    c0 = wave_speed(A0(0))
    t_shock = 2 * sigma / c0
    t_max = 2 * t_shock
    
    t = np.linspace(0, t_max, 500)
    s = np.linspace(-5 * sigma, 5 * sigma + mean, 100)

    plot_characteristics(t, s, intersections = "False")