import numpy as np
import matplotlib.pyplot as plt

g = 9.81
alpha = np.radians(5)
f = 0.02 
w = 10.0  # width of the rectangular channel (m)
A_max = 5.0  # maximum initial cross-sectional area (m^2)
L = 20.0  # half-length of the initial pulse (m)

# Wave speed function
def wave_speed(A):
    l =  w + (2 * A/w)
    return (3/2) * np.sqrt(g * np.sin(alpha) / l) * np.sqrt(A)

# Initial condition: triangular pulse
def A0(s):
    return np.where(np.abs(s) <= L, A_max * (1 - np.abs(s)/L), 0)

# Generate characteristic curves
def characteristic(s0, t):
    return s0 + wave_speed(A0(s0)) * t

# Shock formation time
t_shock = 2 * L / wave_speed(A_max)

t_max = 3 * t_shock
t = np.linspace(0, t_max, 500)
s = np.linspace(-2*L, 2*L, 1000)

plt.figure(figsize=(10, 6))
for s0 in np.linspace(-L, 2 * L, 200):
    plt.plot(characteristic(s0, t), t, 'b-', lw=1, alpha=0.5)

# Labels and title
plt.xlabel('Distance along river, $s$ (m)')
plt.ylabel('Time, $t$ (s)')
plt.title('Characteristic Diagram for Rectangular Channel')
plt.legend()
plt.grid(alpha=0.3)
plt.show()