import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from mpl_toolkits.mplot3d import Axes3D

"""
u = 0.0797eV and u' = 0.0975eV 
Src = https://arxiv.org/pdf/1805.06819

T1(u,u') = u sigma_0+u' sigma_x

T2(u,u') = u sigma_0+u'(cos(2pi/3)sigma_x+sin(2pi/3)sigma_y)

T3(u,u') = u sigma_0+u'(cos(2pi/3)sigma_x-sin(2pi/3)sigma_y)
"""

sigma_0 = np.array([[1, 0], [0, 1]])
sigma_x = np.array([[0, 1], [1, 0]])
sigma_y = np.array([[0, -1j], [1j, 0]])

u = 0.0797  
u_p = 0.0975 

T1 = u * sigma_0 + u_p * sigma_x
T2 = u * sigma_0 + u_p * (np.cos(2 * np.pi / 3) * sigma_x + np.sin(2 * np.pi / 3) * sigma_y)
T3 = u * sigma_0 + u_p * (np.cos(2 * np.pi / 3) * sigma_x - np.sin(2 * np.pi / 3) * sigma_y)

K1 = np.array([4 * np.pi / (3 * np.sqrt(3)), 0])
K2 = -K1

num_kpoints = 25
kx_vals = np.linspace(-4 * np.pi / 3, 4 * np.pi / 3, num_kpoints)
ky_vals = np.linspace(-4 * np.pi / 3, 4 * np.pi / 3, num_kpoints)
kx_grid, ky_grid = np.meshgrid(kx_vals, ky_vals)

bands = np.zeros((num_kpoints, num_kpoints, 4))
for i, kx in enumerate(kx_vals):
    for j, ky in enumerate(ky_vals):
        k = np.array([kx, ky])
        H_k = np.block([[np.zeros((2, 2)), T1 * np.exp(1j * k[0]) + T2 * np.exp(1j * (k[0] + k[1])) + T3 * np.exp(1j * k[1])],
                         [T1.T.conj() * np.exp(-1j * k[0]) + T2.T.conj() * np.exp(-1j * (k[0] + k[1])) + T3.T.conj() * np.exp(-1j * k[1]), np.zeros((2, 2))]])
        eigvals = eigh(H_k, eigvals_only=True)
        bands[i, j, :] = eigvals

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
for band in range(4):
    ax.plot_surface(kx_grid, ky_grid, bands[:, :, band], cmap='coolwarm', alpha=0.8)

ax.set_xlabel('$k_x$')
ax.set_ylabel('$k_y$')
ax.set_zlabel('E (meV)')
ax.set_title(r'3D Bandstructure of Relaxed TBG at $\theta \approx 1.08^\circ$')
plt.show()
