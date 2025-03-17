import numpy as np
import scipy.optimize as opt
from characteristic_plot_general import wave_speed
# Constants
g = 9.81
f = 0.05 
shape = 'Rectangular'

def find_A_for_wave_speed(target_v, A_guess=10):
    def equation(A):
        return wave_speed(A) - target_v  # Solve for A where wave_speed(A) = target_v

    # Define bounds for A (must be positive)
    A_min, A_max = 0.1, 5000  # Avoid A=0 for numerical stability

    # Find all possible solutions using Brent's method
    roots = []
    try:
        root1 = opt.brentq(equation, A_min, A_max)
        roots.append(root1)
    except ValueError:
        pass  # No solution in this range

    # If multiple solutions exist, check other ranges
    A_range = np.linspace(A_min, A_max, 100)
    for i in range(len(A_range) - 1):
        try:
            root = opt.brentq(equation, A_range[i], A_range[i + 1])
            if not np.isclose(root, roots).any():  # Avoid duplicates
                roots.append(root)
        except ValueError:
            continue  # No root in this small interval

    return roots

# Example usage
target_v = float(input("Enter wave speed (m/s): "))
A_solutions = find_A_for_wave_speed(target_v)

if A_solutions:
    print(f"Possible values of A for wave speed {target_v} m/s: {A_solutions}")
else:
    print("No valid A found for the given wave speed.")
    