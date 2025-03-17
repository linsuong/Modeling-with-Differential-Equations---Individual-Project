import numpy as np
import matplotlib.pyplot as plt
L = 200
R = 500
Delta = 2*R/L

A_L = 3

def l_rectangular(A, w=10):
    """Compute wetted perimeter for a rectangular channel."""
    return (w**2 + (2 * A)) / w

def u_bar_rectangular(A, g=9.81, f=0.05, alpha=np.arctan(0.02)):
    """Compute velocity u_bar for rectangular cross-section."""
    return np.sqrt((g * np.sin(alpha) * A) / (f * l_rectangular(A)))

def Q(A):
    """Compute discharge Q = A * u_bar."""
    return A * u_bar_rectangular(A)

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

plt.title('Godanov for Rectangle', fontsize = 20)
plt.plot(s,A_0)
plt.xlabel('Distance m')
plt.ylabel('Area $m^2$')
plt.plot(s,A_20_array, label = 't = 20s')
plt.plot(s,A_40_array, label = 't = 40s')
plt.plot(s,A_60_array, label = 't = 60s')
plt.plot(s,A_80_array, label = 't = 80s')
plt.plot(s,A_100_array, label = 't = 100s')
plt.legend()
plt.show()