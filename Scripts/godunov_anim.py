import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Define simulation parameters
Nx = 100        # Number of spatial points
dx = 1.0        # Spatial step size
dt = 0.01       # Time step size
t_max = 2.0     # Maximum simulation time
frames = 100    # Number of frames in animation
time_per_frame = t_max / frames

# Define flux function
def flux(A, g=9.81, alpha=0.01, f=0.03, l=1.0):
    return np.sqrt(g * np.sin(alpha) * A**3 / (f * l))

# Godunov solver function
def godunov_solver(A, dx, dt, t_step):
    Nt = int(t_step / dt)  # Number of time steps per frame
    A_new = A.copy()
    for _ in range(Nt):
        A_temp = A_new.copy()
        for i in range(1, len(A) - 1):
            F_left = flux(A[i-1])
            F_right = flux(A[i])
            A_temp[i] -= dt/dx * (F_right - F_left)
        A_new = A_temp.copy()
    return A_new

# Initial condition (Gaussian bump to simulate water surge)
x = np.linspace(0, Nx * dx, Nx)
A = np.ones(Nx) + np.exp(-0.1 * (x - 50) ** 2)  # Small perturbation

# Set up the figure and axis
fig, ax = plt.subplots()
ax.set_xlim(0, Nx * dx)
ax.set_ylim(np.min(A), np.max(A) * 1.1)
line, = ax.plot(x, A, label="Water level")
ax.set_xlabel("Distance along river")
ax.set_ylabel("Cross-sectional Area A")
ax.set_title("Flash Flood Evolution")
ax.legend()

# Animation function
def animate(frame):
    global A
    A = godunov_solver(A, dx, dt, time_per_frame)  # Update A over a small time step
    line.set_ydata(A)  # Update plot
    return line,

# Create animation
ani = animation.FuncAnimation(fig, animate, frames=frames, interval=50, blit=True)

plt.show()
