import numpy as np
import scipy.optimize as opt
#from characteristic_plot_general import wave_speed
# Constants
g = 9.81
f = 0.05 
shape = 'Semi'

def wave_speed(A):
    g = 9.81
    f = 0.05

    alpha = np.arctan(0.02)
    w = 20.0  # width of the channel (m)

    l = w + ((2 * A) / w)
    l_prime = 2/w
        
    #return (3/2) * np.sqrt(g * w * np.sin(alpha)/f) * ((3/2)* np.sqrt(A/(2*A + w**2)) - (((A)/((2*A + w**2)) ** (3/2))))

    #def wave_speed(A):
    if shape == 'Wedge':
        alpha = np.arctan(0.065)
        theta = np.pi/3
        
        l = np.sqrt((8 * A)/np.sin(theta))
        l_prime = 2 / (np.sqrt(2 * A * np.sin(theta)))
        
    #return (5/4) * np.sqrt((g * np.sin(alpha)/f) * np.sqrt(np.sin(theta)/8)) * (A ** (1/4))
    
    if shape == 'Semi':
        alpha = np.arctan(0.03)
        theta = np.pi
        
        l = np.sqrt((2 * A)/(theta - np.sin(theta))) * theta
        l_prime = theta/(np.sqrt(2 * A * (theta - np.sin(theta))))
    
    if shape == 'Parabola':
        alpha = np.arctan(0.1)
        w = 5
        
        var_0 = (3 * A) / (2 * (w ** 2))

        # Calculate l
        l = ((2 * (w ** 3))/(3 * A)) * (var_0 * np.sqrt(1 + (var_0 ** 2)) + np.log(np.abs(np.sqrt(1 + (var_0 ** 2)) + var_0))) 
        
        # Calculate l'
        alpha = (2 * (w ** 3)) / (3 * A)
        alpha_prime = -((2 * (w ** 3)) / (3 * (A ** 2)))

        beta_prime = ((1 + 2 * (var_0 ** 2)) / np.sqrt(1 + (var_0 ** 2))) * (3 / (2 * (w ** 2)))
        gamma_prime = (9 * A) / ((4 * (w ** 4)) * np.sqrt(1 + (var_0 ** 2)))
        
        gamma = np.sqrt(1 + (var_0 ** 2)) + var_0

        l_prime = (alpha_prime * l) + alpha * (beta_prime + (gamma_prime / gamma))
        
    ratio = A/l
    
    return (1/2) * np.sqrt((g * np.sin(alpha))/f) *((3 * np.sqrt(ratio)) - (((ratio) ** (3/2)) * l_prime))


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
    