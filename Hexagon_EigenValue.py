import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig

def delta(x, y):
    return np.exp(-1j * (-x)) + np.exp(-1j * (x/2 + np.sqrt(3)*y/2)) + np.exp(-1j * (x/2 - np.sqrt(3)*y/2))

def primeDelta(x, y):
    return np.exp(1j * (-x)) + np.exp(1j * (x/2 + np.sqrt(3)*y/2)) + np.exp(1j * (x/2 - np.sqrt(3)*y/2))

def H(x, y):
    t = 2
    tprime = .5
    return np.array([
        [0, -t * delta(x, y), -tprime, 0],
        [-t * primeDelta(x, y), 0, 0, -tprime],
        [-tprime, 0, 0, -t * delta(x, y)],
        [0, -tprime, -t * primeDelta(x, y), 0]
    ], dtype=complex)

def eigenvalues(x, y):
    H_matrix = H(x, y,)
    return eig(H_matrix)[0] 

x_vals = np.linspace(-4, 4, 50)
y_vals = np.linspace(-4, 4, 50)
X, Y = np.meshgrid(x_vals, y_vals)

Z1 = np.zeros_like(X, dtype=float)
Z2 = np.zeros_like(X, dtype=float)
Z3 = np.zeros_like(X, dtype=float)
Z4 = np.zeros_like(X, dtype=float)

for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        eigvals = np.real(eigenvalues(X[i, j], Y[i, j]))  
        Z1[i, j], Z2[i, j], Z3[i, j], Z4[i, j] = np.sort(eigvals)  

# Plotting
fig = plt.figure(figsize=(10, 8))
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z1, cmap='coolwarm', alpha=0.8)
ax.plot_surface(X, Y, Z2, cmap='coolwarm', alpha=0.8)
ax.plot_surface(X, Y, Z3, cmap='coolwarm', alpha=0.8)
ax.plot_surface(X, Y, Z4, cmap='coolwarm', alpha=0.8)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Energy')
ax.set_title('Eigenvalues of the Hamiltonian')

plt.show()
