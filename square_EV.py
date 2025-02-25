import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import eig

def squareEigenvalue(Type, bilayer=False, save=False):

    if Type == "Square":

        def delta(x, y):
            return 2*(np.cos(x)+np.cos(y))

        def H(x, y):
            t = 2
            tprime = 1
            if bilayer:
                return np.array([
                            [0, -t * delta(x, y), -tprime, 0],
                            [-t * delta(x, y), 0, 0, -tprime],
                            [-tprime, 0, 0, -t * delta(x, y)],
                            [0, -tprime, -t * delta(x, y), 0]
                            ], dtype=complex) 
            else:
               return np.array([
                            [0, -t * delta(x, y)],
                            [-t * delta(x, y), 0]
                            ], dtype=complex)  

    elif Type == "Checkered":

        def delta(x, y):
            return 2*np.cos(x/2)+2*np.cos(y/2) #???

    elif Type == "Lieb":

        def delta(x, y):
            return 2*np.cos(x/2)
        def gamma(x, y):
            return 2*np.cos(y/2)
        
        def H(x, y):
            t = 2
            tprime = 1
            if bilayer:
                return np.array([
                                [0, -t*delta(x, y), -t*gamma(x, y), -tprime, 0, 0],
                                [-t*delta(x, y), 0, 0, 0, -tprime, 0],
                                [-t*gamma(x, y), 0, 0, 0, 0, -tprime],
                                [-tprime, 0, 0, 0, -t*delta(x, y), -t*gamma(x, y)], 
                                [0, -tprime, 0, -t*delta(x, y), 0, 0],
                                [0, 0, -tprime, -t*gamma(x, y), 0, 0]
                                ], dtype=complex) 
            else:
                return np.array([
                                [0, -t*delta(x, y), -t*gamma(x, y)],
                                [-t*delta(x, y), 0, 0],
                                [-t*gamma(x, y), 0, 0]
                                ], dtype=complex) 
        
    else:
        return 0


    def eigenvalues(x, y):
        H_matrix = H(x, y)  
        return np.real(eig(H_matrix)[0])  

    x_vals = np.linspace(-6, 6, 20)
    y_vals = np.linspace(-6, 6, 20)
    X, Y = np.meshgrid(x_vals, y_vals)

    Z_list = []
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            eigvals = np.sort(eigenvalues(X[i, j], Y[i, j]))  
            if len(Z_list) < len(eigvals):
                Z_list.extend([np.zeros_like(X, dtype=float) for _ in range(len(eigvals) - len(Z_list))])
            for k in range(len(eigvals)):
                Z_list[k][i, j] = eigvals[k]

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    for Z in Z_list:
        ax.plot_surface(X, Y, Z, cmap="viridis",alpha=0.8)

    ax.set_xlabel(r'$ak_x$')
    ax.set_ylabel(r'$ak_y$')
    ax.set_zlabel(r'$\frac{Energy}{t}$')
    ax.set_title('Eigenvalues of the Hamiltonian')
    if save:
        if bilayer:
            plt.savefig(f'{Type}BilayerEigenvaluePlot.pdf')
        else:
            plt.savefig(f'{Type}EigenvaluePlot.pdf')
    plt.show()

"""
Showcase
"""

squareEigenvalue("Square", bilayer=True)
squareEigenvalue("Lieb", bilayer=True)
