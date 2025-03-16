from godunov import godunov
from characteristic_plot_general import plot_characteristics
from characteristic_plot_general import wave_speed
from characteristic_plot_general import A0
import numpy as np
from matplotlib import pyplot as plt

g = 9.81
f = 0.1 

A_L = 50
V = 50000

mean = 1
sigma = 2 #std dev

shape = 'Rectangle'

c0 = wave_speed(A0(0))
t_shock = 2 * sigma / c0
t_max = 3 * t_shock
print(t_max)
t = np.linspace(0, t_max, 500)
s = np.linspace(-5 * sigma, 5 * sigma + mean, 100)


shape = 'Rectangle'

# Discretization
N = 100
L = 10
t = 0
t_end = 2.5
#x = np.linspace(-sigma - mean, 5 * sigma + mean, N)

plot_characteristics(t, s)
godunov(t, t_end, s, N, L)
plt.close('all')
