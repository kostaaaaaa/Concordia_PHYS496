import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange


def phi(n, s):
    return np.degrees(np.arccos((3 * n**2 + 3 * n * s + (s**2) / 2) / (3 * n**2 + 3 * n * s + s**2)))


def bz_length(n, s):
    return (4 * np.pi / 3) / np.sqrt(3 * n**2 + 3 * n * s + s**2)


n_values = np.arange(1,900,1)
s_values = [1] + list(primerange(1, 20))

phi_dict = {}

for n in n_values:
    for s in s_values:
        phi_value = phi(n, s)
        if 0 <= phi_value <= 2:
            bz_len = bz_length(n, s)
            if phi_value not in phi_dict or bz_len > phi_dict[phi_value]:
                phi_dict[phi_value] = bz_len

sorted_phi = np.array(sorted(phi_dict.keys()))
max_bz_lengths = np.array([phi_dict[p] for p in sorted_phi])

plt.figure(figsize=(10, 6))
plt.scatter(sorted_phi, max_bz_lengths, s=2, color='blue', alpha=0.7)
plt.plot(sorted_phi, max_bz_lengths, linestyle='-', linewidth=0.8, color='black', alpha=0.6)
plt.axvspan(phi(30,1)-0.01, phi(30,1)+0.01, color='red', alpha=0.2, label='Magic Angle Region')
plt.xlim(0, 2)
plt.yscale('log')
plt.title("Max BZ Length vs Phi (Log Scale)", fontsize=14)
plt.xlabel("Phi (degrees)", fontsize=12)
plt.ylabel("Max BZ Length (log scale)", fontsize=12)
plt.legend()
plt.show()

