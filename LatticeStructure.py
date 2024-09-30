import numpy as np
import matplotlib.pyplot as plt

class Lattice:
    def __init__(self, vectors, basis):
        self.vectors = np.array(vectors)
        self.basis = np.array(basis)
        self.lattice_type = self.determine_lattice_type()
        self.x_range = (-3,3)
        self.y_range = (-3,3)
    
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

    def generate_lattice_points(self):
        points = []
        points2 = []

        for i in range(self.x_range[0], self.x_range[1]):
            for j in range(self.y_range[0], self.y_range[1]):
                displacement = i * self.vectors[0] + j * self.vectors[1]
                points.append(displacement + self.basis[0])
                if len(self.basis) > 1:
                    points2.append(displacement + self.basis[1])

        points = np.array(points)
        if points2:
            points2 = np.array(points2)
            return points, points2
        else:
            return points

        
    
    def plot_lattice(self):
        
        if self.lattice_type=="Hexagon":
            points, points2 = self.generate_lattice_points()
        else:
            points = self.generate_lattice_points()

        plt.figure(figsize=(8, 8))
        plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50)

        if len(self.basis) > 1:
            plt.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50)
        plt.title(f'{self.lattice_type} Lattice')
        plt.axis('equal')
        plt.xticks([],[])
        plt.yticks([],[])
        plt.show()
    
    def plot_bilayer(self):
        a=0.5 #vertical distance between layers
        p11,p12 = self.generate_lattice_points()
        p21,p22 = self.generate_lattice_points()
        

        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection="3d")

        ax.scatter(p11[:, 0], p11[:, 1],0, color=(0.1, 0.2, 0.5, 0.5), s=50, label="A1")
        ax.scatter(p12[:, 0], p12[:, 1],0, color=(0.5, 0.1, 0.2, 0.5), s=50, label="B1")

        lattice_shift_x = 0.5
        lattice_shift_y = np.sqrt(3)/6

        ax.scatter(p21[:, 0]-lattice_shift_x, p21[:, 1]-lattice_shift_y,a, color=(0.1, 0.5, 0.2,0.5), s=50, label="A2")
        ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=(0.5, 0.5, 0.2,0.5), s=50, label="B2")
        
                
        ax.set_zlim(-a, a*2)
        ax.set_title("3D Bilayer Graphene")
        ax.axis('equal')
        ax.legend()
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