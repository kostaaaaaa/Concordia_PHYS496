import numpy as np
import matplotlib.pyplot as plt
from sympy import primerange

def phi(n, s):
    return np.degrees(np.arccos((3*n**2 + 3*n*s + (s**2)/2) / (3*n**2 + 3*n*s + s**2)))

n_values = np.arange(1, 100, 1)
s_values = list(primerange(1, 20))

phi_points = []
phi_1_points = []
for n in n_values:
    phi_1_points.append((phi(n,1),n,1))
    for s in s_values:
        value = phi(n,s)
        phi_points.append((value, n, s))

plt.figure(figsize=(10, 6))
plt.scatter([p[0] for p in phi_points], np.zeros(len(phi_points)), s=0.2, c='blue', alpha=0.5, picker=True)
plt.scatter([p[0] for p in phi_1_points], np.zeros(len(phi_1_points)), s=2, c='red', marker='o')
plt.yticks([])
plt.xlim(0.2, 1.8)
plt.axvspan(1.09, 1.11, color='red', alpha=0.2, label='Magic Angle Region')
plt.title("Distribution of Unique Phi Values", fontsize=14)
plt.xlabel("Phi (degrees)", fontsize=12)
plt.legend()


def on_pick(event):
    ind = event.ind[0]
    phi_value, n, s = phi_points[ind]
    print(f"Phi: {phi_value:.4f} degrees, n: {n}, s: {s}")

plt.gcf().canvas.mpl_connect('pick_event', on_pick)
plt.show()
