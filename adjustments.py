from Lattice_class import Lattice
import numpy as np

v = [[1,0],[1/2,np.sqrt(3)/2]]
b=[[0,0]]

tr = Lattice(1, v, b, 5)
tr.plot_lattice_with_twist_Circles(save=True)
