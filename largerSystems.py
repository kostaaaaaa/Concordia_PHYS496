from Lattice_class import Lattice
import numpy as np
import matplotlib.pyplot as plt

a1 = np.array([1,0])
a2 = np.array([0,1])

def u1(n):
    return np.array(n*a1 + n*a2+a1)

def u2(n):
    return n*a1+n*a2+a2

u = u1(2)
v = u2(2)

dot_product = np.dot(u, v)
cross_product = np.cross(u, v)

magnitude_u = np.linalg.norm(u)
magnitude_v = np.linalg.norm(v)

cos_theta = dot_product / (magnitude_u * magnitude_v)
sin_theta = cross_product / (magnitude_u * magnitude_v)

rotation_matrix = np.array([
    [cos_theta, -sin_theta],
    [sin_theta, cos_theta]
])


###---###

lattice_square = Lattice(1, [[1,0],[0,1]], [[0,0]], 5)
lattice_square.plot_lattice()
