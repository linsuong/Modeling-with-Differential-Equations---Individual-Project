import numpy as np
import matplotlib.pyplot as plt

# Parameters
L = 1          # Domain length
N = 2000            # Number of cells
dx = L / N         # Cell width
cfl = 0.1          # CFL condition (must be <= 1 for stability)
total_time = 10   # Total simulation time
alpha = np.pi / 6  # Angle of inclination (30 degrees)
sin_alpha = np.sin(alpha)

# Function for v(A)
def v(A):
    # Assume l(A) is proportional to A^{1/2} (e.g., for a rectangular or triangular cross-section)
    l_A = np.sqrt(A)  # Example: l(A) = A^{1/2}
    return (3 * np.sqrt(A) / 2) * np.sqrt(sin_alpha / l_A)

# Cell centers and initial condition (Gaussian pulse)
x = np.linspace(0.5 * dx, L - 0.5 * dx, N)
A_initial = 1.0 + 0.5 * np.exp(-(x - L / 4)**2 / 0.5)  # Initial cross-sectional area
A = A_initial.copy()

# Time step based on CFL condition
dt = cfl * dx / np.max(v(A))  # Use maximum wave speed for stability

# Number of time steps
nt = int(total_time / dt)

# Godunov method (periodic boundary conditions)
for _ in range(nt):
    A_old = A.copy()
    for i in range(N):
        # Left and right cell indices (periodic boundary conditions)
        i_left = (i - 1) % N
        i_right = (i + 1) % N

        # Compute v(A) at cell centers
        v_left = v(A_old[i_left])
        v_right = v(A_old[i_right])

        # Upwind flux (Godunov method for scalar conservation laws)
        if v_left > 0:
            flux_left = v_left * A_old[i_left]
        else:
            flux_left = v_left * A_old[i]

        if v_right > 0:
            flux_right = v_right * A_old[i]
        else:
            flux_right = v_right * A_old[i_right]

        # Update A using the flux difference
        A[i] = A_old[i] - (dt / dx) * (flux_right - flux_left)

# Plot results
plt.plot(x, A_initial, 'k--', label='Initial Condition')
plt.plot(x, A, 'r-', label='Numerical Solution')
plt.xlabel('s (Distance Along River)')
plt.ylabel('A (Cross-Sectional Area)')
plt.legend()
plt.title('Godunov Method for Flash Flood PDE')
plt.show()