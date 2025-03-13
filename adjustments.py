from Lattice_class import Lattice
import numpy as np

v = [[1,0],[1/2,np.sqrt(3)/2]]
b =[[0,0]]

tr = Lattice(1, v, b, 10)
tr.plot_bz_difference(1.1, s=47)
