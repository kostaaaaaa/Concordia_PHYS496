import numpy as np
import matplotlib.pyplot as plt

n_values = np.arange(1, 31, 1)  

def compute_phi(n, s):
    numerator = 2*n**2+2*n*s
    denominator = 2*n**2+2*n*s+s**2
    return np.arccos(numerator / denominator)

s_values = range(1, 6)

plt.figure(figsize=(8, 6))

for s in s_values:
    phi_values = compute_phi(n_values, s)
    plt.plot(n_values, np.degrees(phi_values), marker='o', linestyle='-', label=f's={s}')

plt.xlabel('n')
plt.ylabel('Angle φ (degrees)')
plt.title('Angle φ vs. n for different s values (Integer n)')
plt.legend()
plt.grid(True)

plt.show()
