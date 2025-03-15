import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange


def phi(n, s):
    return np.degrees(np.arccos((3 * n**2 + 3 * n * s + (s**2) / 2) / (3 * n**2 + 3 * n * s + s**2)))

def bz_length(n, s):
    factor = (2 * np.pi) / (3 * (n**2) + 3 * n * s + s**2)
    b1 = factor * np.array([(2 * n + s), -s / np.sqrt(3)])
    b2 = factor * np.array([-n, (3 * n + 2 * s) / np.sqrt(3)])
    return np.linalg.norm((2 * b1 + b2) / 3)


n_values = np.arange(1, 10000, 1)
s_values = [1] + list(primerange(1, 30))

sorted_phi = []
max_bz_lengths = []
n_values_list = []
s_values_list = []

for n in n_values:
    for s in s_values:
        phi_value = phi(n, s)
        if 0 <= phi_value <= 2:  
            sorted_phi.append(phi_value)
            max_bz_lengths.append(bz_length(n, s))
            n_values_list.append(n)
            s_values_list.append(s)


sorted_phi = np.array(sorted_phi)
max_bz_lengths = np.array(max_bz_lengths)
n_values_list = np.array(n_values_list)
s_values_list = np.array(s_values_list)


sort_indices = np.argsort(sorted_phi)
sorted_phi = sorted_phi[sort_indices]
max_bz_lengths = max_bz_lengths[sort_indices]
n_values_list = n_values_list[sort_indices]
s_values_list = s_values_list[sort_indices]

colors = ['red' if s == 1 else 'blue' for s in s_values_list]

plt.figure(figsize=(10, 6))
plt.scatter(sorted_phi, max_bz_lengths, s=5, c=colors, alpha=0.7, picker=True)
plt.plot(sorted_phi, max_bz_lengths, linestyle='-', linewidth=0.8, color='black', alpha=0.6)
plt.axvspan(phi(30,1)-0.01, phi(30,1)+0.01, color='red', alpha=0.2, label='Magic Angle Region')
plt.xlim(0, 2)
plt.title("Max BZ Length vs Phi", fontsize=14)
plt.xlabel("Phi (degrees)", fontsize=12)
plt.ylabel("Max BZ Length", fontsize=12)
plt.legend()

def on_pick(event):
    ind = event.ind[0]
    print(f"n: {n_values_list[ind]}, s: {s_values_list[ind]}, Phi: {sorted_phi[ind]:.4f} degrees, Max BZ Length: {max_bz_lengths[ind]:.4f}")

plt.gcf().canvas.mpl_connect('pick_event', on_pick)
plt.show()
