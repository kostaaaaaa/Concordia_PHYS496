from Lattice_class import Lattice
import numpy as np


"""Initialized Parameters"""
lattice_distance = 1.0
lattice_sites = 20
lattice_degrees = 5

vectors_triangle = [[1, 0], [0.5, np.sqrt(3)/2]]
basis_triangle = [[0, 0]]

vectors_square = [[1, 0], [0, 1]]
basis_square = [[0, 0]]

vectors_hexagon = [[1, 0], [0.5, np.sqrt(3)/2]]
basis_hexagon = [[0, 0], [0.5, np.sqrt(3)/6]]

"""Class Initialization"""
lattice_triangle = Lattice(lattice_distance, vectors_triangle, basis_triangle, lattice_sites)
lattice_square = Lattice(lattice_distance, vectors_square, basis_square, lattice_sites)
lattice_hexagon = Lattice(lattice_distance, vectors_hexagon, basis_hexagon, lattice_sites)

"""Plotting"""
"""lattice_triangle.plot_bilayer(lattice_degrees)
lattice_triangle.plot_bilayer_align(lattice_degrees)
"""
"""lattice_square.plot_bilayer(lattice_degrees)"""
"""lattice_square.plot_bilayer_align(lattice_degrees)"""
j = lattice_square.find_aligned_overlap_points(lattice_degrees)
k = lattice_square.find_overlap_points(lattice_degrees)
print(j)
print(k)
print(len(k))
"""
lattice_hexagon.plot_bilayer(lattice_degrees)
lattice_hexagon.plot_bilayer_align(lattice_degrees)
"""