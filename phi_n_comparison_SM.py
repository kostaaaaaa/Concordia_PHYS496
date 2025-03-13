import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange

s_values = list(primerange(31, 61))
n_values = np.arange(1000, 1400, 10)

def compute_phi(n, s):
    numerator = 3 * n**2 + 3 * n * s + (s**2) / 2
    denominator = 3 * n**2 + 3 * n * s + s**2
    return np.arccos(numerator / denominator)

plt.figure(figsize=(8, 6))

for s in s_values:
    phi_values = compute_phi(n_values, s)
    plt.scatter([s] * len(phi_values), np.degrees(phi_values), marker='o', alpha=0.6, label=f's={s}')

plt.xlabel('s values')
plt.ylabel('Angle φ (degrees)')
plt.ylim(0, 2)  # Limit y-axis to focus on 0-2 degrees
plt.title('φ vs. s values')
plt.axhline(y=1.1, color='r', linestyle='-', label='Magic Angle (1.1°)')
plt.legend(loc='lower right')
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()
