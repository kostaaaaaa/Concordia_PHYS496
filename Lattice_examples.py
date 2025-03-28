from Lattice_class import Lattice
import numpy as np


"""Initialized Parameters"""
lattice_distance = 1
lattice_sites = 20
lattice_degrees = 5
save_state = True

vectors_triangle = [[1, 0], [0.5, np.sqrt(3)/2]]
basis_triangle = [[0, 0]]

vectors_square = [[1, 0], [0, 1]]
basis_square = [[0, 0]]

vectors_hexagon = [[1, 0], [0.5, np.sqrt(3)/2]]
basis_hexagon = [[0, 0], [0.5, np.sqrt(3)/6]]

vectors_kagome = [[1,0], [0.5, np.sqrt(3)/2]]
basis_kagome = [[0,0], [0.5,0],[1/4, np.sqrt(3)/4]]

vectors_lieb = [[1, 0], [0, 1]]
basis_lieb = [[0, 0], [1/2, 0], [0, 1/2]]

"""Class Initialization"""
lattice_triangle = Lattice(lattice_distance, vectors_triangle, basis_triangle, lattice_sites)
lattice_square = Lattice(lattice_distance, vectors_square, basis_square, lattice_sites)
lattice_hexagon = Lattice(lattice_distance, vectors_hexagon, basis_hexagon, lattice_sites)
lattice_kagome = Lattice(lattice_distance, vectors_kagome, basis_kagome, lattice_sites)
lattice_lieb = Lattice(lattice_distance, vectors_lieb, basis_lieb, lattice_sites)

"""Plotting"""
lattice_triangle_s1 = Lattice(lattice_distance, vectors_triangle, basis_triangle, 60)
lattice_triangle_s11= Lattice(lattice_distance, vectors_triangle, basis_triangle, 340)

lattice_triangle_s1.plot_bilayer_2D(1.084549,point_size=2,isAlign=True)
lattice_triangle_s11.plot_bilayer_2D(1.087791,point_size=0.005,isAlign=True)
