from Lattice_class import Lattice
import numpy as np

v = [[1,0],[0,1]]
b=[[0,0],[1/2,0],[0,1/2]]

sq = Lattice(1, v, b, 10)
sq.plot_lattice_superlattice(n=2)
