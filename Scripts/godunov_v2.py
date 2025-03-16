import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp

def l(A):
    if shape == 'Rectangle':
        w = 10
        return w + (2 * A) / w
    
    if shape == 'Wedge':
        theta = np.pi / 6
        return np.sqrt((8 * A) / (np.sin(theta)))
    
    if shape == 'Semi':
        theta = np.pi / 3
        return np.sqrt((2 * A) / (theta - np.sin(theta))) * theta

def u_bar(A):
    if shape == 'Rectangle':
        alpha = np.arctan(0.02)  # Correctly left in radians
        
    if shape == 'Wedge':
        alpha = np.arctan(0.08)  # Correctly left in radians
        
    if shape == 'Semi':
        alpha = np.arctan(0.04)  # Added a value for alpha in the Semi case
    return np.sqrt((g * np.sin(alpha) * A) / (f * l(A)))

def int_cond(s):
    norm = np.sqrt(2 * np.pi * (sigma ** 2))
    return (V / norm) * np.exp(-((s - b) ** 2) / (sigma ** 2)) + A_L

def Q(A):
    return A * u_bar(A)

def cell_average(c, x, i):
    """
    Compute the cell average of c over the cell [x[i], x[i+1]].
    """
    #dx = x[i+1] - x[i]
    dx = 0.1
    return (1 / dx) * np.trapz(c[i:i+2], x[i:i+2])

def rhs(t, c):
    """
    Right-hand side of the ODE system for the cell averages.
    """
    N = len(c)
    Q_flux = np.zeros(N + 1)
    dx = x[1] - x[0]  # Assuming uniform grid
    
    # Compute fluxes at cell interfaces
    for i in range(N + 1):
        if i == 0:
            # Left boundary: use boundary condition
            Q_flux[i] = Q(c[0])
        elif i == N:
            # Right boundary: use boundary condition
            Q_flux[i] = Q(c[-1])
        else:
            # Interior interfaces: use upwind flux
            if u_bar(c[i - 1]) > 0:
                Q_flux[i] = Q(c[i - 1])
            else:
                Q_flux[i] = Q(c[i])
    
    # Compute the right-hand side of the ODE system
    dc_dt = np.zeros(N)
    for i in range(N):
        dc_dt[i] = -(Q_flux[i + 1] - Q_flux[i]) / dx
    
    return dc_dt

def godunov(t, t_end, x, N, L):
    plt.figure(figsize=(16, 10))

    dx = L / N
    
    # Initialize cell averages
    c0 = np.array([cell_average(int_cond(x), x, i) for i in range(N)])

    # Time points at which to solve the ODEs
    t_eval = np.linspace(0, t_end, int(t_end / 0.1) + 1)
    
    # Solve the ODE system using RK45
    sol = solve_ivp(rhs, [0, t_end], c0, t_eval=t_eval, method='RK45')
    
    # Plot the solution
    for i, t in enumerate(sol.t):
        if i % 100 == 0:  # Plot every 10th time step
            plt.plot(x, sol.y[:, i], label=f't={t:.1f}s')
            
    plt.xlabel('Distance along river, $s$', fontsize=20)
    plt.ylabel('Cross sectional area, $A$', fontsize=20)
    plt.xlim(0, 5 * sigma + b)
    plt.legend()
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.title('Evolution of cross sectional area, $A$ across length for different $t$', fontsize=20)
    plt.savefig(f'/Users/linusong/Repositories/Modeling-with-Differential-Equations---Individual-Project/Figures/{shape}_godunov_rk45.pdf')
    plt.show()
    
# Constants
g = 9.81
f = 0.1

A_L = 0
V = 5000

b = 200
sigma = 100

shape = 'Rectangle'

# Discretization
N = 100
L = 10
t = 0
t_end = 5000
x = np.linspace(0, 1000, N)

# Run the simulation
godunov(t, t_end, x, N, L)