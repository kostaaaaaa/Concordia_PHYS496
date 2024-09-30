import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Lattice:
    def __init__(self, vectors, basis):
        self.vectors = np.array(vectors)
        self.basis = np.array(basis)
        self.lattice_type = self.determine_lattice_type()
    
    def determine_lattice_type(self):
        angle = np.arccos(np.dot(self.vectors[0], self.vectors[1]) /(np.linalg.norm(self.vectors[0]) * np.linalg.norm(self.vectors[1])))
        angle_deg = np.degrees(angle)
        if np.isclose(angle_deg, 90):
            return "Square"
        elif np.isclose(angle_deg, 60) and len(self.basis)>1:
            return "Hexagon"
        elif np.isclose(angle_deg, 60):
            return "Triangle"
        else:
            return "Unknown Type"

    def generate_lattice_points(self, range_x=(-5, 5), range_y=(-5, 5)):
        points = []
        for i in range(range_x[0], range_x[1]):
            for j in range(range_y[0], range_y[1]):
                displacement = i * self.vectors[0] + j * self.vectors[1]
                for b in self.basis:
                    points.append(displacement + b)
        return np.array(points)
    
    def plot_lattice(self):
        points = self.generate_lattice_points()
        plt.figure(figsize=(8, 8))
        for idx, b in enumerate(self.basis):
            basis_points = points + b
            color = (0.1, 0.2, 0.5,0.5) if idx ==0 else (0.5, 0.1, 0.2,0.5)
            plt.scatter(basis_points[:,0], basis_points[:,1], color=color, s=50)
        plt.title(f'{self.lattice_type} Lattice')
        plt.axis('equal')
        plt.xticks([],[])
        plt.yticks([],[])
        plt.show()
    
    def plot_bilayer(self):
        a=0.5 #vertical distance between layers
        p1 = self.generate_lattice_points()
        p2 = self.generate_lattice_points()
        

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection="3d")

        for idx, b in enumerate(self.basis):
            basis_p1 = np.c_[p1+b, np.zeros(len(p1))]
            color = (0.1, 0.2, 0.5,0.5) if idx ==0 else (0.5, 0.1, 0.2,0.5)
            ax.scatter(basis_p1[:,0], basis_p1[:,1], basis_p1[:,2], color=color, s=50, label="Layer 1")

        for idx, b in enumerate(self.basis):
            basis_p2 = np.c_[p2+b, a*np.ones(len(p2))]
            color = (0.1, 0.5, 0.2,0.5) if idx ==0 else (0.5, 0.5, 0.2,0.5)
            ax.scatter(basis_p2[:,0]+0.5, basis_p2[:,1]+(np.sqrt(3)/6),basis_p2[:,2], color=color, s=50, label="Layer 2")
                
        ax.set_zlim(-a, a*2)
        ax.set_title("3D Bilayer Graphene")
        ax.legend()
        plt.savefig("3D_Bilayer.pdf")
        plt.show()
    
    def __str__(self):
        return "This is a "+self.lattice_type+" lattice"

vectors_triangle = [[1, 0], [0.5, np.sqrt(3)/2]]
basis_triangle = [[0, 0]]

vectors_square = [[1, 0], [0, 1]]
basis_square = [[0, 0]]

vectors_hexagon = [[1, 0], [0.5, np.sqrt(3)/2]]
basis_hexagon = [[0, 0], [0.5, np.sqrt(3)/6]]

lattice_triangle = Lattice(vectors_triangle, basis_triangle)
lattice_triangle.plot_lattice()

lattice_square = Lattice(vectors_square, basis_square)
lattice_square.plot_lattice()

lattice_hexagon = Lattice(vectors_hexagon, basis_hexagon)
lattice_hexagon.plot_lattice()
lattice_hexagon.plot_bilayer()
