import numpy as np
import matplotlib.pyplot as plt

def flux(A, g=9.81, alpha=0.01, f=0.03, l=1.0):
    """
    Compute the flux function F(A) based on the given PDE.
    """
    return np.sqrt(g * np.sin(alpha) * A**3 / (f * l))

def godunov_solver(A0, dx, dt, t_max):
    """
    Solves the hyperbolic PDE using the Godunov method.
    """
    Nt = int(t_max / dt)  # Number of time steps
    Nx = len(A0)          # Number of spatial points
    A = A0.copy()
    
    for _ in range(Nt):
        A_new = A.copy()
        for i in range(1, Nx-1):
            # Compute numerical flux using Godunov's method
            F_left = flux(A[i-1])
            F_right = flux(A[i])
            
            # Update solution using conservation law
            A_new[i] -= dt/dx * (F_right - F_left)
        
        A = A_new.copy()
    
    return A

# Define simulation parameters
Nx = 100        # Number of spatial points
dx = 1.0        # Spatial step size
dt = 0.01       # Time step size
t_max = 1.0     # Maximum simulation time

# Initial condition (Gaussian bump to simulate water surge)
x = np.linspace(0, Nx*dx, Nx)
A0 = np.ones(Nx) + np.exp(-0.1*(x - 50)**2)  # Small perturbation

# Solve using Godunov's method
A_final = godunov_solver(A0, dx, dt, t_max)

# Plot results
plt.plot(x, A0, label="Initial Condition")
plt.plot(x, A_final, label="Final Condition")
plt.xlabel("Distance along river")
plt.ylabel("Cross-sectional Area A")
plt.legend()
plt.title("Flash Flood Simulation using Godunov Method")
plt.show()