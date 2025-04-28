import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint

def l(A):
    if shape == 'Rectangle':
        w = 10
        return w + ((2 * A) / w)
    
    if shape == 'Wedge':
        theta = np.pi / 6
        return np.sqrt((8 * A) / (np.sin(theta)))
    
    if shape == 'Semi':
        theta = np.pi/3
        
        return np.sqrt((2 * A) / (theta - np.sin(theta))) * theta

    if shape == 'Parabola':
        w = 10
        
        var_0 = (3 * A) / (2 * (w ** 2))
        
        return ((2 * (w ** 3))/(3 * A)) * (var_0 * np.sqrt(1 + (var_0 ** 2)) + np.log(np.abs(np.sqrt(1 + (var_0 ** 2)) + var_0))) 
            
def u_bar(A):
    if shape == 'Rectangle':
        alpha = np.arctan(0.02)  # Correctly left in radians
        
    if shape == 'Wedge':
        alpha = np.arctan(0.08)  # Correctly left in radians
        
    if shape == 'Semi':
        alpha = np.arctan(0.2)
    
    if shape == 'Parabola':
        alpha = np.arctan(0.1)
        
    return np.sqrt((g * np.sin(alpha) * A) / (f * l(A)))


def int_cond(s):
    norm = np.sqrt(2 * np.pi * (sigma ** 2))
    return (V / norm) * np.exp(-((s - mean) ** 2) / (sigma ** 2)) + A_L

def Q(A):
    return A * u_bar(A)

def godunov(t, N, L):
    """
    godunov solver, plots every 1s.

    Args:
        t (_type_): duration of time for simulation
        N (_type_): number of grids
        L (_type_): total length of river
    """
    
    def conservation_law(c, Q):
        # Compute flux differences
        dc_dt = np.zeros_like(c)
        #dc_dt[0]= - Q(c[1]) + Q(c[0])
        dc_dt[0] = A_L
        
        #for i in range(1, len(c)):
        dc_dt[i] = -Q(c[i]) + Q(c[i-1])
            
        return dc_dt
    
    t = 100
    
    time = np.linspace(0, t, t)
    c0 = int_cond(1)
    sol = odeint(conservation_law, c0, time)
    
    print(sol)
    A_list = np.zeros(N)
    
    grid = np.linspace(0, L, N)
    
    for i in range(0, N):    

        A = np.zeros(N)
        c0 = int_cond(s)
        
        time = np.linspace(0, t, t)
        sol = odeint(conservation_law, c0, time)

        print(sol)
        #plt.plot(time, sol, label = 't = {t}')
        
    plt.show()



##############################################################################################################
if __name__ == "__main__":
    # Constants
    g = 9.81
    f = 0.05

    A_L = 3
    V = 5000

    mean = 0
    sigma = 10

    shape = 'Rectangle'

    N = 500
    L = 10


    #plt.close('all')
    #godunov(t, t_end, s, N, L)
    
    godunov(10, 50, 100)