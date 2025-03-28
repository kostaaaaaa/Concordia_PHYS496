
import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange

def phi(n, s):
    return np.degrees(np.arccos((3 * n**2 + 3 * n * s + (s**2) / 2) / (3 * n**2 + 3 * n * s + s**2)))

def bz_length(n, s):
    return (4 * np.pi / 3) / np.sqrt(3 * n**2 + 3 * n * s + s**2)

n_values = np.arange(1, 900, 1)
s_values = [1] + list(primerange(1, 20))

phi_dict = {}
ns_dict = {}

for n in n_values:
    for s in s_values:
        phi_value = phi(n, s)
        if 0 <= phi_value <= 2:
            bz_len = bz_length(n, s)
            if phi_value not in phi_dict or bz_len > phi_dict[phi_value]:
                phi_dict[phi_value] = bz_len
                ns_dict[phi_value] = (n, s)  

sorted_phi = np.array(sorted(phi_dict.keys()))
max_bz_lengths = np.array([phi_dict[p] for p in sorted_phi])

fig, ax = plt.subplots(figsize=(10, 6))
ax.scatter(sorted_phi, max_bz_lengths, s=2, color='blue', alpha=0.7)
ax.plot(sorted_phi, max_bz_lengths, linestyle='-', linewidth=0.8, color='black', alpha=0.6)
ax.axvspan(phi(30, 1)-0.01, phi(30, 1)+0.01, color='red', alpha=0.2, label='Magic Angle Region')
ax.set_xlim(0, 2)
ax.set_yscale('log')
ax.set_title("Max BZ Length vs Phi (Log Scale)", fontsize=14)
ax.set_xlabel("Phi (degrees)", fontsize=12)
ax.set_ylabel("Max BZ Length (log scale)", fontsize=12)
ax.legend()

def on_click(event):
    if event.xdata is None or event.ydata is None:
        return  

    closest_phi = min(sorted_phi, key=lambda x: abs(x - event.xdata))
    closest_n, closest_s = ns_dict[closest_phi]
    closest_bz_len = phi_dict[closest_phi]

    print(f"Clicked Point -> n: {closest_n}, s: {closest_s}, phi: {closest_phi:.6f}, BZ Length: {closest_bz_len:.6f}")

fig.canvas.mpl_connect('button_press_event', on_click)

plt.show()

