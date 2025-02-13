from Lattice_class import Lattice
import numpy as np
import matplotlib.pyplot as plt

vectors = [[1,0],[0.5,np.sqrt(3)/2]]
basis = [[0,0],[0.5,np.sqrt(3)/6]]

hex = Lattice(1,vectors,basis,10)

hex.plot_lattice_WC_BZ_comparison(n=4)
"""
Here is the following pattern:
    n=1 points=2
    n=2 points=8
    n=3 points=18
    n=4 points=32
    n=5 points=50
    n=6 points=72

    points = 2(n^2)
"""

def atoms(s):
    return 2*(s**2)

s = np.arange(1,20,1)

"""
plt.figure()
plt.scatter(atoms(s),s)
plt.show()
"""
