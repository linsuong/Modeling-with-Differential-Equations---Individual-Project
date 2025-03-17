import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import solve_ivp

# Constants
g = 9.81
f = 0.05

A_L = 3
V = 5000

mean = 0
sigma = 100
shape = 'Wedge'

# Define l(A)
def l(A):
    if shape == 'Rectangle':
        w = 10
        return w + (2 * A) / w
    if shape == 'Wedge':
        theta = np.pi / 6
        return np.sqrt((8 * A) / np.sin(theta))
    if shape == 'Semi':
        theta = np.pi/3
        return np.sqrt((2 * A) / (theta - np.sin(theta))) * theta
    if shape == 'Parabola':
        w = 10
        var_0 = (3 * A) / (2 * w**2)
        return (2 * w**3 / (3 * A)) * (var_0 * np.sqrt(1 + var_0**2) + np.log(np.sqrt(1 + var_0**2) + var_0))

# Define u_bar(A)
def u_bar(A):
    if shape == 'Rectangle':
        alpha = np.arctan(0.02)
    if shape == 'Wedge':
        alpha = np.arctan(0.08)
    if shape == 'Semi':
        alpha = np.arctan(0.2)
    if shape == 'Parabola':
        alpha = np.arctan(0.1)
    return np.sqrt((g * np.sin(alpha) * A) / (f * l(A)))

# Define Q(A)
def Q(A):
    return A * u_bar(A)

# Initial condition for A(s, t)
def int_cond(s):
    norm = np.sqrt(2 * np.pi * sigma**2)
    return (V / norm) * np.exp(-((s - mean)**2) / sigma**2) + A_L

# Define the ODE system
def dA_dt(t, A):
    dA = np.zeros_like(A)
    dA[0] = A_L  # Boundary condition at the left end
    for i in range(1, len(A)):
        dA[i] = -Q(A[i]) - Q(A[i - 1])
    return dA

# Solver function
def solver(N):
    s = np.linspace(-5 * sigma, 5 * sigma, N)
    A0 = int_cond(s)  # Initial condition

    # Time points to evaluate the solution
    t_eval = np.linspace(0, 100, 200)

    # Solve the ODE system
    sol = solve_ivp(dA_dt, [0, t_eval[-1]], A0, t_eval=t_eval, method='RK45')

    # Plot the solution at different times
    plt.figure(figsize=(12, 8))
    for i, t in enumerate(t_eval):
        if i % 20 == 0:  # Plot every 20th time point
            plt.plot(s, sol.y[:, i], label=f't = {t:.2f}')

    plt.xlabel('Distance, s (m)')
    plt.ylabel('A(s, t)')
    plt.title(f'Evolution of A(s, t) Over Time (Shape: {shape})')
    plt.legend()
    plt.grid()
    plt.show()

# Parameters
N = 100  # Number of spatial points

# Call the solver
solver(N)