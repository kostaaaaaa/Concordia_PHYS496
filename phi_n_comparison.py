<<<<<<< HEAD
import numpy as np
import matplotlib.pyplot as plt

def find_integer_n_with_plot(phi_range, type, step=0.001, tol=0.02):
    phi_values = []
    n_values = []
    phi_range_values = np.arange(phi_range[0], phi_range[1], step)
    
    if type == "Triangle":
        for phi in phi_range_values:
            numerator = 1 - 2 * np.cos(phi)
            denominator = 3 * (np.cos(phi) - 1) - np.sqrt(3) * np.sin(phi)
        
            if abs(denominator) < tol:  # Avoid division by zero
                continue
            
            n = numerator / denominator
            if abs(n - round(n)) < tol:  # Check if n is close to an integer
                phi_values.append(phi)
                n_values.append(int(round(n)))

    elif type == "Square":
        for phi in phi_range_values:
            numerator = -np.cos(phi)
            denominator = np.cos(phi)-np.sin(phi)-1
        
            if abs(denominator) < tol:  # Avoid division by zero
                continue
            
            n = numerator / denominator
            if abs(n - round(n)) < tol:  # Check if n is close to an integer
                phi_values.append(phi)
                n_values.append(int(round(n)))


    return phi_values, n_values

# Define the range of phi and step size
phi_range = (0, np.pi/3)
phi_range2 = (0,np.pi/2)
phi_vals, n_vals = find_integer_n_with_plot(phi_range,"Triangle", step=0.001)

phi_valss, n_valss = find_integer_n_with_plot(phi_range2, "Square", step=0.001)

# Plot the results
plt.figure(figsize=(10, 6))
plt.scatter(phi_vals, n_vals, color='blue', marker='o', label='Integer n values')
plt.scatter(np.pi/36, 6, color=(0.5,0.1,0.1,0.2)) # A point we know exists
plt.title('Integer n Values as a Function of Phi', fontsize=14)
plt.xlabel('Phi (radians)', fontsize=12)
plt.ylabel('n (Integer)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.scatter(phi_valss, n_valss, color='blue', marker='o', label='Integer n values')
plt.scatter(np.pi/36, 11, color=(0.5,0.1,0.1,0.2)) # A point we know exists
plt.title('Integer n Values as a Function of Phi', fontsize=14)
plt.xlabel('Phi (radians)', fontsize=12)
plt.ylabel('n (Integer)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()
=======
import numpy as np
import matplotlib.pyplot as plt

def find_integer_n_with_plot(phi_range, type, step=0.001, tol=0.02):
    phi_values = []
    n_values = []
    phi_range_values = np.arange(phi_range[0], phi_range[1], step)
    
    if type == "Triangle":
        for phi in phi_range_values:
            numerator = 1 - 2 * np.cos(phi)
            denominator = 3 * (np.cos(phi) - 1) - np.sqrt(3) * np.sin(phi)
        
            if abs(denominator) < tol:  # Avoid division by zero
                continue
            
            n = numerator / denominator
            if abs(n - round(n)) < tol:  # Check if n is close to an integer
                phi_values.append(phi)
                n_values.append(int(round(n)))

    elif type == "Square":
        for phi in phi_range_values:
            numerator = -np.cos(phi)
            denominator = np.cos(phi)-np.sin(phi)-1
        
            if abs(denominator) < tol:  # Avoid division by zero
                continue
            
            n = numerator / denominator
            if abs(n - round(n)) < tol:  # Check if n is close to an integer
                phi_values.append(phi)
                n_values.append(int(round(n)))


    return phi_values, n_values

# Define the range of phi and step size
phi_range = (0, np.pi/3)
phi_range2 = (0,np.pi/2)
phi_vals, n_vals = find_integer_n_with_plot(phi_range,"Triangle", step=0.001)

phi_valss, n_valss = find_integer_n_with_plot(phi_range2, "Square", step=0.001)

# Plot the results
plt.figure(figsize=(10, 6))
plt.scatter(phi_vals, n_vals, color='blue', marker='o', label='Integer n values')
plt.scatter(np.pi/36, 6, color=(0.5,0.1,0.1,0.2)) # A point we know exists
plt.title('Integer n Values as a Function of Phi', fontsize=14)
plt.xlabel('Phi (radians)', fontsize=12)
plt.ylabel('n (Integer)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()

plt.figure(figsize=(10, 6))
plt.scatter(phi_valss, n_valss, color='blue', marker='o', label='Integer n values')
plt.scatter(np.pi/36, 11, color=(0.5,0.1,0.1,0.2)) # A point we know exists
plt.title('Integer n Values as a Function of Phi', fontsize=14)
plt.xlabel('Phi (radians)', fontsize=12)
plt.ylabel('n (Integer)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.show()
>>>>>>> 696e6cb (2024-12-01)
