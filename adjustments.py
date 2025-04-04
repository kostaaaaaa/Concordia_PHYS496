from Lattice_class import Lattice
import numpy as np

v = [[1,0],[0,1]]
b =[[0,0]]

tr = Lattice(1, v, b, 10)
tr.plot_bz_difference(s=4,degrees=0.94701)
