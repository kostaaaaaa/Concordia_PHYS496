from Lattice_class import Lattice
import numpy as np


"""Initialized Parameters"""
lattice_distance = 1
lattice_sites = 5
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

"""
lattice_triangle.plot_bilayer(lattice_degrees, save=save_state)
lattice_triangle.plot_bilayer_align(lattice_degrees, save=save_state)
lattice_triangle.plot_lattice_with_twist_Circles(save=save_state)

lattice_square.plot_bilayer(lattice_degrees, save=save_state)
lattice_square.plot_bilayer_align(lattice_degrees, save=save_state)
lattice_square.plot_lattice_with_twist_Circles(save=save_state)

lattice_hexagon.plot_bilayer(lattice_degrees, save=save_state)
lattice_hexagon.plot_bilayer_align(lattice_degrees, save=save_state)

"""
lattice_hexagon = Lattice(lattice_distance, vectors_hexagon, basis_hexagon, 5*lattice_sites)
lattice_hexagon.plot_bilayer(degrees=5)
lattice_hexagon.plot_bilayer_align(degrees=5)
lattice_kagome = Lattice(lattice_distance, vectors_kagome, basis_kagome, 4*lattice_sites)
lattice_kagome.plot_bilayer(degrees=5)
lattice_kagome.plot_bilayer_align(degrees=5)
lattice_lieb = Lattice(lattice_distance, vectors_lieb, basis_lieb, 6*lattice_sites)
lattice_lieb.plot_bilayer(degrees=15)
lattice_lieb.plot_bilayer_align(degrees=5)

