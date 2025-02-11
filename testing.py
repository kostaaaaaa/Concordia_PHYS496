from Lattice_class import Lattice
import numpy as np

vectors = [[1,0],[0.5,np.sqrt(3)/2]]
basis = [[0,0],[0.5,np.sqrt(3)/6]]

hex = Lattice(1,vectors,basis,10)

hex.plot_lattice_nthWC(n=2)
"""
Here is the following pattern:
    n=1 points=2
    n=2 points=8
    n=3 points=18
    n=4 points=32
    n=5 points=50
    n=6 points=72

    up to n=4 it seems it follows a f(x)=2^{n+1} if n%2==0 or a 2^{n}+2 if else
"""
