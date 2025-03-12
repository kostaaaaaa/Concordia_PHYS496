import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange

s_values = list(primerange(1, 35))
n_values = np.arange(1, 920, 10)

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
plt.ylim(0,2)
plt.title('Angle φ vs. n and s values (Integer n)')
plt.axhline(y=1.1, color='r', linestyle='-')
plt.legend()
plt.grid(True)

plt.show()

