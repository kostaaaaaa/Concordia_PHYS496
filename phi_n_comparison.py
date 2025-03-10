import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange

s_values = np.arange(1,9,1)

n_values = np.arange(-90, 71, 1)

def compute_phi(n, s):
    numerator = 3 * n**2 + 3 * n * s + (s**2) / 2
    denominator = 3 * n**2 + 3 * n * s + s**2
    return np.arccos(numerator / denominator)

plt.figure(figsize=(8, 6))

for s in s_values:
    phi_values = compute_phi(n_values, s)
    plt.plot(n_values, np.degrees(phi_values), marker='o', linestyle='-', label=f's={s}')

plt.xlabel('n')
plt.ylabel('Angle φ (degrees)')
plt.title('Angle φ vs. n for prime s values (Integer n)')
plt.axhline(y=1.1, color='r', linestyle='-')
plt.legend()
plt.grid(True)

plt.show()

