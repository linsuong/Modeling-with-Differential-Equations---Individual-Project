import numpy as np
from matplotlib import pyplot as plt

L = 200
R = 500
Delta = 2*R/L

A_L = 3

shape = 'Wedge'

def l(A):
    if shape == 'Rectangle':
        w = 20
        return w + (2 * A) / w
    
    if shape == 'Wedge':
        theta = np.pi / 3
        return np.sqrt((8 * A) / (np.sin(theta)))
    
    if shape == 'Semi':
        theta = np.pi
        
        return np.sqrt((2 * A) / (theta - np.sin(theta))) * theta

    if shape == 'Parabola':
        w = 10
        
        var_0 = (3 * A) / (2 * (w ** 2))
        
        return ((2 * (w ** 3))/(3 * A)) * (var_0 * np.sqrt(1 + (var_0 ** 2)) + np.log(np.abs(np.sqrt(1 + (var_0 ** 2)) + var_0))) 
            
def u_bar(A):
    g = 9.81
    f = 0.05
    if shape == 'Rectangle':
        alpha = np.arctan(0.02)
        
    if shape == 'Wedge':
        alpha = np.arctan(0.065)
        
    if shape == 'Semi':
        alpha = np.arctan(0.3)
    
    if shape == 'Parabola':
        alpha = np.arctan(0.1)
        
    return np.sqrt((g * np.sin(alpha) * A) / (f * l(A)))

'''
def int_cond(s):
    norm = np.sqrt(2 * np.pi * (sigma ** 2))
    return (V / norm) * np.exp(-((s - mean) ** 2) / (sigma ** 2)) + A_L
'''

def Q(A):
    return A * u_bar(A)

def initial_conditions(s, A_L=3, V=5000, sigma=100, mean=0):
    """Generate initial conditions for A using a Gaussian distribution."""
    norm = np.sqrt(2 * np.pi * (sigma ** 2))
    return (V / norm) * np.exp(-((s - mean) ** 2) / (sigma ** 2)) + A_L

s = np.linspace(-R, R, L, dtype=int)
A_0 = np.array([initial_conditions(si) for si in s])

def solve_ODE(Q, A_0_array, Delta, dt, boundary_condition, t_max,):
    """
    Solve the ODE dA_i/dt + (Q(A_i) - Q(A_{i-1})) / Delta = 0 for a specified time t_max.

    Parameters:
        Q (function): Function that computes Q(A).
        A_0_array (array-like): Array of initial conditions for A at t = 0.
        Delta (float): Spatial step size.
        t_max (float): Maximum time to evolve the system.
        dt (float): Time step size for Euler's method.
        boundary_condition (function): Function that defines the boundary condition for A.

    Returns:
        A_t_array (numpy array): Approximated values of A at time t_max for each index i.
    """
    num_points = len(A_0_array)
    A_t_array = np.array(A_0_array, dtype=float)  # Initialize A_t with A_0 values
    time = 0

    while time < t_max:
        A_new = np.copy(A_t_array)  # Make a copy to update values without interference
        
        for i in range(1, num_points):
            # Finite difference scheme for the spatial derivative
            dA_dt = -(Q(A_t_array[i]) - Q(A_t_array[i-1])) / Delta
            A_new[i] = A_t_array[i] + dt * dA_dt  # Update A[i] for the current time step

        # Apply the boundary condition at i = 0 (assuming it's a fixed value or some function of time)
        A_new[0] = boundary_condition(time)

        A_t_array = np.copy(A_new)  # Update the solution for the next time step
        time += dt

    return A_t_array

def boundary_condition(t):
    """Example boundary condition function"""
    return A_L  # Example: fixed boundary condition at A_0(t) = 0


dt = 0.01

A_20_array = solve_ODE(Q, A_0, Delta, dt, boundary_condition, t_max = 20)
A_40_array = solve_ODE(Q, A_0, Delta, dt, boundary_condition, t_max = 40)
A_60_array = solve_ODE(Q, A_0, Delta, dt, boundary_condition, t_max = 60)
A_80_array = solve_ODE(Q, A_0, Delta, dt, boundary_condition, t_max = 80)
A_100_array = solve_ODE(Q, A_0, Delta, dt, boundary_condition, t_max = 100)

plt.figure(figsize=(16, 10))
plt.grid()
plt.title(f'Evolution of cross sectional area, \n $A$ across length for different $t$, {shape}', fontsize=40)

plt.xlabel('Distance along river, $s$', fontsize=40)
plt.ylabel('Cross sectional area, $A$', fontsize=40)

plt.plot(s,A_0, label = 't = 0s')
plt.plot(s,A_20_array, label = 't = 20s')
plt.plot(s,A_40_array, label = 't = 40s')
plt.plot(s,A_60_array, label = 't = 60s')
plt.plot(s,A_80_array, label = 't = 80s')
plt.plot(s,A_100_array, label = 't = 100s')
plt.legend(fontsize=40)
plt.tick_params(axis='both', which='major', labelsize=40)
plt.savefig(f'Figures/{shape}_godunov.pdf', bbox_inches="tight", pad_inches = 0.2)
#plt.show()