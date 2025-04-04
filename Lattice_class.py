import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Lattice:
    def __init__(self, lattice_distance, vectors, basis, num_sites):
        """ initialization of lattice parameters """
        self.lattice_distance = lattice_distance
        self.vectors = self.lattice_distance*np.array(vectors)
        self.basis = self.lattice_distance*np.array(basis)
        self.lattice_type = self.determine_lattice_type()
        self.reciprocal_vectors = self.get_reciprocal()
        self.nn = self.get_nn()
        self.num_sites = num_sites
        self.x_range = (-num_sites, num_sites)
        self.y_range = self.x_range

    def get_reciprocal(self):
        """ solves the reciprocal lattice with the intialized vectors """
        a1 = self.vectors[0]
        a2 = self.vectors[1]

        eq = np.array([[a1[0], a1[1], 0,     0],
                      [a2[0], a2[1], 0,     0],
                      [0,     0,     a1[0], a1[1]],
                      [0,     0,     a2[0], a2[1]]])

        ans = np.array([2 * np.pi, 0, 0, 2 * np.pi])

        solution = np.linalg.solve(eq, ans)
        b1 = solution[0:2]
        b2 = solution[2:4]

        return np.array([b1, b2])
    
    def calculate_reciprocal(self, v1, v2):
        """ calculate the reciporcal of given lattice vectors """
        eq = np.array([[v1[0], v1[1], 0,     0],
                      [v2[0], v2[1], 0,     0],
                      [0,     0,     v1[0], v1[1]],
                      [0,     0,     v2[0], v2[1]]])
        ans = np.array([2 * np.pi, 0, 0, 2 * np.pi])

        solution = np.linalg.solve(eq, ans)
        b1 = solution[0:2]
        b2 = solution[2:4]

        return np.array([b1, b2])
    
    def determine_lattice_type(self):
        """ using the angle the vectors produce the lattice type is determined """
        angle = np.arccos(np.dot(self.vectors[0], self.vectors[1]) /(np.linalg.norm(self.vectors[0]) * np.linalg.norm(self.vectors[1])))
        angle_deg = np.degrees(angle)

        if np.isclose(angle_deg, 90) and len(self.basis)==3:
            return "Lieb"
        elif np.isclose(angle_deg, 90):
            return "Square"
        elif np.isclose(angle_deg, 60) and len(self.basis)==2:
            return "Hexagon"
        elif np.isclose(angle_deg, 60) and len(self.basis)>2:
            return "Kagome"
        elif np.isclose(angle_deg, 60):
            return "Triangle"
        else:
            return "unknown type"
        
    def get_nn(self):
        """ generates the nearest neighbors of the 2d lattice strucutre """
        nearestneighbors = []
        for basis in self.basis:
            for vector in self.vectors:
                nearestneighbors.append(basis-vector)
                nearestneighbors.append(basis+vector)
        return np.array(nearestneighbors)
    
    def generate_bzedges(self, v1, v2):
        """ brillouin zone edges to plot brillouin zone within reciporcal space """
        if self.lattice_type == "Square":
            bz_edges = [
                    0.5 * (v1+v2),   
                    0.5 * (v1-v2), 
                    -0.5 * (v1+v2), 
                    -0.5 * (v1-v2) 
            ]
        else:
            bz_edges = [
                    (1/3)*(2*v1+v2),  
                    (1/3)*(v1+2*v2),  
                    (1/3)*(-v1+v2), 
                    (1/3)*(-2*v1-v2), 
                    (1/3)*(-v1-2*v2), 
                    (1/3)*(v1-v2)   
            ]
        
        return bz_edges

    def generate_lattice_points(self):
        """ generate lattice points to use for plotting, takes into account the basis """
        points = []
        points2 = []
        points3 = []

        for i in range(self.x_range[0], self.x_range[1]):
            for j in range(self.y_range[0], self.y_range[1]):
                displacement = i * self.vectors[0] + j * self.vectors[1]
                points.append(displacement + self.basis[0])
                if len(self.basis)>1:
                    points2.append(displacement + self.basis[1])
                if len(self.basis) == 3:
                    points3.append(displacement + self.basis[2])

        points = np.array(points)
        if points3:
            points3 = np.array(points3)
            points2 = np.array(points2)
            return points, points2, points3
        elif points2 and len(self.basis)==2:
            points2 = np.array(points2)
            return points, points2
        else:
            return points
    
    def generate_rotated_points(self, angle):
        """ degree rotation function for later use when discussing super lattices"""
        angle = angle *(np.pi/180)
        rot_matrix = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        points = []
        points2 = []
        points3 = []

        rot_vectors = []
        for a in self.vectors:
            rot_vectors.append(rot_matrix.dot(a))

        for i in range(self.x_range[0], self.x_range[1]):
            for j in range(self.y_range[0], self.y_range[1]):
                displacement = i * rot_vectors[0] + j * rot_vectors[1]
                points.append(displacement + self.basis[0])
                if len(self.basis) > 1:
                    points2.append(displacement + rot_matrix.dot(self.basis[1]))
                if len(self.basis) ==3:
                    points3.append(displacement + rot_matrix.dot(self.basis[2]))

        points = np.array(points)
        if points3:
            points3 = np.array(points3)
            points2 = np.array(points2)
            return points, points2, points3
        elif points2 and len(self.basis)==2:
            points2 = np.array(points2)
            return points, points2
        else:
            return points

    def generate_reciprocal_points(self):
        """ generates lattice points with reciprocal vectors"""
        points = []
        points2 = []

        for i in range(self.x_range[0], self.x_range[1]):
            for j in range(self.y_range[0], self.y_range[1]):
                displacement = i * self.reciprocal_vectors[0] + j * self.reciprocal_vectors[1]
                points.append(displacement + self.basis[0])
                if len(self.basis) > 1:
                    points2.append(displacement + self.basis[1])

        points = np.array(points)
        if points2:
            points2 = np.array(points2)
            return points, points2
        else:
            return points
        
    def generate_superlattice_reciprocal_points(self, degrees):
        """ generates reciprocal lattice points for the superlattice defined by the rotation angle """
        phi = np.radians(degrees)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        
        if self.lattice_type == "Square":
            n = np.round(-cos_phi / (cos_phi - sin_phi - 1))
            factor = 2 * np.pi / (2 * n**2 + 2 * n + 1)
            b1 = factor * np.array([(n + 1), n])
            b2 = factor * np.array([-n, (n + 1)])
        elif self.lattice_type == "Triangle":
            n = np.round((1 - 2 * cos_phi) / (3 * (cos_phi - 1) - np.sqrt(3) * sin_phi))
            factor = 2 * np.pi / (3 * n**2 + 3 * n + 1)
            b1 = factor * np.array([(2 * n + 1), -1 / np.sqrt(3)])
            b2 = factor * np.array([-n, (3 * n + 2) / np.sqrt(3)])
        else:
            raise ValueError("superlattice reciprocal not defined for this lattice type.")
        
        points = []
        for i in range(self.x_range[0], self.x_range[1]):
            for j in range(self.y_range[0], self.y_range[1]):
                displacement = i * b1 + j * b2
                points.append(displacement)
    
        return np.array(points)

    def plot_superlattice_reciprocal(self, degrees, save=False):
        """ plots the reciprocal lattice structure for the superlattice """
        points = self.generate_superlattice_reciprocal_points(degrees)

        plt.figure(figsize=(8, 8))
        plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50)
        plt.title(f'{self.lattice_type} Superlattice Reciprocal (Rotation: {degrees}°)')
        plt.axis('equal')
        if save:
            plt.savefig(f'{self.lattice_type}_superlattice_reciprocal_{degrees}.pdf')
        plt.show()

    def plot_superlattice_reciprocal_with_vectors(self, degrees, save=False):
        """ 
        plots the reciprocal lattice structure for the superlattice 
        with the reciprocal lattice vectors b1 and b2 displayed.
        """
        phi = np.radians(degrees)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        
        if self.lattice_type == "Square":
            n = np.round(-cos_phi / (cos_phi - sin_phi - 1))
            factor = 2 * np.pi / (2 * n**2 + 2 * n + 1)
            b1 = factor * np.array([(n + 1), n])
            b2 = factor * np.array([-n, (n + 1)])
        elif self.lattice_type == "Triangle":
            n = np.round((1 - 2 * cos_phi) / (3 * (cos_phi - 1) - np.sqrt(3) * sin_phi))
            factor = 2 * np.pi / (3 * n**2 + 3 * n + 1)
            b1 = factor * np.array([(2 * n + 1), -1 / np.sqrt(3)])
            b2 = factor * np.array([-n, (3 * n + 2) / np.sqrt(3)])
        else:
            raise ValueError("superlattice reciprocal not defined for this lattice type.")
        
        points = self.generate_superlattice_reciprocal_points(degrees)

        plt.figure(figsize=(8, 8))
        plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50, label="reciprocal points")

        origin = np.array([0, 0])
        plt.quiver(origin[0], origin[1], b1[0], b1[1], angles='xy', scale_units='xy', scale=1, color='r', width=0.005, label="b1")
        plt.quiver(origin[0], origin[1], b2[0], b2[1], angles='xy', scale_units='xy', scale=1, color='g', width=0.005, label="b2")

        plt.title(f'{self.lattice_type} Superlattice Reciprocal with Vectors (Rotation: {degrees}°)')
        plt.axis('equal')
        plt.legend()
        if save:
            plt.savefig(f'{self.lattice_type}_superlattice_reciprocal_vectors_{degrees}.pdf')
        plt.show()


    def plot_lattice(self, save=False, border=True):
        """ plots 2d lattice structure """
        plt.figure(figsize=(8, 8))
        if self.lattice_type=="Hexagon":
            points, points2 = self.generate_lattice_points()
            plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=25)
            plt.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=25)
        elif self.lattice_type =="Kagome":
            points, points2, points3 = self.generate_lattice_points()
            plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50)
            plt.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50)
            plt.scatter(points3[:, 0], points3[:, 1], color=(0.2, 0.5, 0.2, 0.5), s=50)
        elif self.lattice_type =="Lieb":
            points, points2, points3 = self.generate_lattice_points()
            plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50)
            plt.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50)
            plt.scatter(points3[:, 0], points3[:, 1], color=(0.2, 0.5, 0.2, 0.5), s=50)
        else:
            points = self.generate_lattice_points()
            plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50)
        
        if border:
            plt.title(f'{self.lattice_type} Lattice')
            plt.axis('equal')
            plt.xlabel('x')
            plt.ylabel('y')
        else:
            plt.axis('equal')
            plt.axis("off")


        if save:
            plt.savefig(f'{self.lattice_type}_lattice_plot.pdf')
        plt.show()

    def plot_lattice_superlattice(self, n, save=False):
        plt.figure(figsize=(8, 8))
        if self.lattice_type=="Hexagon":
            points, points2 = self.generate_lattice_points()
            plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50)
            plt.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50)
            a1 = self.vectors[0]
            a2 = self.vectors[1]           
            
            a1 = [((3*n/2)+1), np.sqrt(3)*n/2]
            a2 = [1/2, np.sqrt(3)*(n+1/2)]

            origin = np.array([0, 0])
            plt.quiver(origin[0], origin[1], a1[0], a1[1], angles='xy', scale_units='xy', scale=1, width=0.005)
            plt.quiver(origin[0], origin[1], a2[0], a2[1], angles='xy', scale_units='xy', scale=1, width=0.005)
            plt.quiver(origin[0], origin[1], a1[0], a1[1], angles='xy', scale_units='xy', scale=1, color=(0.1, 0.2, 0.5, 0.8), width=0.005)
            plt.quiver(origin[0], origin[1], a2[0], a2[1], angles='xy', scale_units='xy', scale=1, color=(0.1, 0.2, 0.5, 0.8), width=0.005)

            ws_edges = [[1/2*(n+1),np.sqrt(3)/2*(n+1/3)],
                        [-n/2, np.sqrt(3)*(n/2+1/3)],
                        [-1/2-n, 1/(2*np.sqrt(3))],
                        [-1/2*(n+1),-np.sqrt(3)/2*(n+1/3)],
                        [n/2, -np.sqrt(3)*(n/2+1/3)],
                        [1/2+n, -1/(2*np.sqrt(3))]]                      
                
            for i in range(len(ws_edges)):
                        start = ws_edges[i]
                        end = ws_edges[(i + 1) % len(ws_edges)]  
                        plt.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.1, 0.2, 0.5, 0.8), label="wigner cell" if i == 0 else "")
            

        elif self.lattice_type =="Kagome":
            points, points2, points3 = self.generate_lattice_points()
            plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50)
            plt.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50)
            plt.scatter(points3[:, 0], points3[:, 1], color=(0.2, 0.5, 0.2, 0.5), s=50)

        elif self.lattice_type == "Lieb":

            points, points2, points3 = self.generate_lattice_points()
            plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50)
            plt.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50)
            plt.scatter(points3[:, 0], points3[:, 1], color=(0.2, 0.5, 0.2, 0.5), s=50)
            a1 = self.vectors[0]
            a2 = self.vectors[1]

            A1 = [ n+1,n ]
            A2 = [ -n,n+1 ]
            ws_edges = [[1/2,n+1/2],
                        [-1/2-n,1/2],
                        [-1/2,-1/2-n],
                        [1/2+n,-1/2]]

            for i in range(len(ws_edges)):
                        start = ws_edges[i]
                        end = ws_edges[(i + 1) % len(ws_edges)]  
                        plt.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.1, 0.2, 0.5, 0.8), label="Wigner Cell" if i == 0 else "")
            
            
            origin = np.array([0, 0])
            plt.quiver(origin[0], origin[1], a1[0], a1[1], angles='xy', scale_units='xy', scale=1, width=0.005)
            plt.quiver(origin[0], origin[1], a2[0], a2[1], angles='xy', scale_units='xy', scale=1, width=0.005)
            plt.quiver(origin[0], origin[1], A1[0], A1[1], angles='xy', scale_units='xy', scale=1, color=(0.1, 0.2, 0.5, 0.8), width=0.005)
            plt.quiver(origin[0], origin[1], A2[0], A2[1], angles='xy', scale_units='xy', scale=1, color=(0.1, 0.2, 0.5, 0.8), width=0.005)            
        else:
            points = self.generate_lattice_points()
            plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=point_size0)
            a1 = self.vectors[0]
            a2 = self.vectors[1]

            a1 = [ n+1,n ]
            a2 = [ -n,n+1 ]
            ws_edges = [[1/2,n+1/2],
                        [-1/2-n,1/2],
                        [-1/2,-1/2-n],
                        [1/2+n,-1/2]]

            for i in range(len(ws_edges)):
                        start = ws_edges[i]
                        end = ws_edges[(i + 1) % len(ws_edges)]  
                        plt.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.1, 0.2, 0.5, 0.8), label="wigner cell" if i == 0 else "")
            
            
            origin = np.array([0, 0])
            plt.quiver(origin[0], origin[1], a1[0], a1[1], angles='xy', scale_units='xy', scale=1, width=0.005)
            plt.quiver(origin[0], origin[1], a2[0], a2[1], angles='xy', scale_units='xy', scale=1, width=0.005)
            plt.quiver(origin[0], origin[1], a1[0], a1[1], angles='xy', scale_units='xy', scale=1, color=(0.1, 0.2, 0.5, 0.8), width=0.005)
            plt.quiver(origin[0], origin[1], a2[0], a2[1], angles='xy', scale_units='xy', scale=1, color=(0.1, 0.2, 0.5, 0.8), width=0.005)            

        plt.title(f'{self.lattice_type} Lattice with Superlattice Formation')
        plt.axis('equal')
        if save:
            plt.savefig(f'{self.lattice_type}LatticeWSuperLattice.pdf')
        plt.show()
    

    def plot_bilayer_2D(self, degrees, point_size=2, isAlign=False, save=False, border=True):
        """ Plots 2D Bilayer of a given lattice, top-down view."""
        layer_1 = (0.1, 0.1, 0.6, 0.5)
        layer_2 = (0.6, 0.1, 0.1, 0.5)
    
        fig, ax = plt.subplots(figsize=(10, 10))
        
        if isAlign:
            shift_x = 0
            shift_y = 0
        else:
            if self.lattice_type in ["Hexagon", "Kagome", "Triangle"]:
                shift_x = self.lattice_distance * 0.5
                shift_y = self.lattice_distance * np.sqrt(3) / 4
            else:
                shift_x = 0.5*self.lattice_distance
                shift_y = 0.5*self.lattice_distance
    
        if self.lattice_type == "Hexagon":
            p11, p12 = self.generate_lattice_points()
            p21, p22 = self.generate_rotated_points(degrees)
        
            ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=point_size, label="Layer 1")
            ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=point_size)
            ax.scatter(p21[:, 0] + shift_x, p21[:, 1] + shift_y, color=layer_2, s=point_size, label="Layer 2")
            ax.scatter(p22[:, 0] + shift_x, p22[:, 1] + shift_y, color=layer_2, s=point_size)
    
        elif self.lattice_type == "Kagome":
            p11, p12, p13 = self.generate_lattice_points()
            p21, p22, p23 = self.generate_rotated_points(degrees)
        
            ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=point_size, label="Layer 1")
            ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=point_size)
            ax.scatter(p13[:, 0], p13[:, 1], color=layer_1, s=point_size)
        
            ax.scatter(p21[:, 0] - shift_x, p21[:, 1] - shift_y, color=layer_2, s=point_size, label="Layer 2")
            ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=point_size)
            ax.scatter(p23[:, 0] - shift_x, p23[:, 1] - shift_y, color=layer_2, s=point_size)
    
        elif self.lattice_type == "Triangle":
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)
        
            ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=point_size, label="Layer 1")
            ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=point_size, label="Layer 2")
    
        elif self.lattice_type == "Lieb":
            p11, p12, p13 = self.generate_lattice_points()
            p21, p22, p23 = self.generate_rotated_points(degrees)
        
            ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=point_size, label="Layer 1")
            ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=point_size)
            ax.scatter(p13[:, 0], p13[:, 1], color=layer_1, s=point_size)
        
            ax.scatter(p21[:, 0] - shift_x, p21[:, 1] - shift_y, color=layer_2, s=point_size, label="Layer 2")
            ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=point_size)
            ax.scatter(p23[:, 0] - shift_x, p23[:, 1] - shift_y, color=layer_2, s=point_size)
    
        else:
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)
        
            ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=point_size, label="a1")
            ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2,s=point_size, label="b2")
        
        if border:

            ax.set_title(f"2D Bilayer {self.lattice_type}")
            ax.axis("equal")
            ax.legend()
        else:
            ax.axis("equal")
            ax.axis("off")
    
        if save:
            if isAlign:
                plt.savefig(f"{self.lattice_type}_bilayer_{degrees}_aligned.pdf")
            else:
                plt.savefig(f"{self.lattice_type}_bilayer_{degrees}_shifted.pdf")
        plt.show()

    def plot_bilayer_2D_MM(self, degrees, point_size=2, isAlign=False, save=False, border=True):
        """ Plots 2D Bilayer of a given lattice, top-down view."""
        layer_1 = (0.1, 0.1, 0.6, 0.5)
        layer_2 = (0.6, 0.1, 0.1, 0.5)
    
        fig, ax = plt.subplots(figsize=(10, 10))
        if self.lattice_type in ["Hexagon", "Kagome", "Triangle"]:
            n = np.round(np.abs((1 - 2 * np.cos(np.radians(degrees))) / (3 * (np.cos(np.radians(degrees)) - 1) - np.sqrt(3) * np.sin(np.radians(degrees)))))
        else:
            n = np.abs(np.cos(np.radians(degrees)) / (np.cos(np.radians(degrees)) - np.sin(np.radians(degrees)) - 1))
        
        A1 = [1/2*(3*n+1),np.sqrt(3)/2*(n+1)]
        A2 = [-1/2, np.sqrt(3)/2*(2*n+1)]
        
        if isAlign:
            shift_x = 0
            shift_y = 0
        else:
            if self.lattice_type in ["Hexagon", "Kagome", "Triangle"]:
                shift_x = self.lattice_distance * 0.5
                shift_y = self.lattice_distance * np.sqrt(3) / 4
            else:
                shift_x = 0.5*self.lattice_distance
                shift_y = 0.5*self.lattice_distance
    
        if self.lattice_type == "Hexagon":
            p11, p12 = self.generate_lattice_points()
            p21, p22 = self.generate_rotated_points(degrees)
        
            ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=point_size, label="Layer 1")
            ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=point_size)
            ax.scatter(p21[:, 0] + shift_x, p21[:, 1] + shift_y, color=layer_2, s=point_size, label="Layer 2")
    
        elif self.lattice_type == "Kagome":
            p11, p12, p13 = self.generate_lattice_points()
            p21, p22, p23 = self.generate_rotated_points(degrees)
        
            ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=point_size, label="Layer 1")
            ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=point_size)
            ax.scatter(p13[:, 0], p13[:, 1], color=layer_1, s=point_size)
        
            ax.scatter(p21[:, 0] - shift_x, p21[:, 1] - shift_y, color=layer_2, s=point_size, label="Layer 2")
    
        elif self.lattice_type == "Triangle":
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)
        
            ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=point_size, label="Layer 1")
            ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=point_size, label="Layer 2")
    
        elif self.lattice_type == "Lieb":
            p11, p12, p13 = self.generate_lattice_points()
            p21, p22, p23 = self.generate_rotated_points(degrees)
        
            ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=point_size, label="Layer 1")
            ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=point_size)
            ax.scatter(p13[:, 0], p13[:, 1], color=layer_1, s=point_size)
        
            ax.scatter(p21[:, 0] - shift_x, p21[:, 1] - shift_y, color=layer_2, s=point_size, label="Layer 2")
    
        else:
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)
        
            ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=point_size, label="a1")
            ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2,s=point_size, label="b2")
        
        origin = np.array([0, 0])
        ax.quiver(origin[0], origin[1], A1[0], A1[1], angles='xy', scale_units='xy', scale=1, color=(0.1, 0.2, 0.5, 0.8), width=0.005)
        ax.quiver(origin[0], origin[1], A2[0], A2[1], angles='xy', scale_units='xy', scale=1, color=(0.1, 0.2, 0.5, 0.8), width=0.005)
        
        if border:

            ax.set_title(f"2D Bilayer {self.lattice_type}")
            ax.axis("equal")
            ax.legend()
        else:
            ax.axis("equal")
            ax.axis("off")
    
        if save:
            if isAlign:
                plt.savefig(f"{self.lattice_type}_bilayer_{degrees}_aligned.pdf")
            else:
                plt.savefig(f"{self.lattice_type}_bilayer_{degrees}_shifted.pdf")
        plt.show()

    def plot_bilayer_comparison(self, degrees, save=False):
        """ Plots 2D Bilayer of a given lattice, top-down view, showing both CW and CCW rotations."""
        layer_1 = (0.1, 0.1, 0.6, 0.5)
        layer_2 = (0.6, 0.1, 0.1, 0.5)
        
        fig, axes = plt.subplots(1, 2, figsize=(20, 10))
        rotations = [degrees, -degrees]
        titles = [f"{self.lattice_type} CCW ({degrees}°)", f"{self.lattice_type} CW ({-degrees}°)"]
        
        if self.lattice_type in ["Hexagon", "Kagome", "Triangle"]:
            n = np.abs((1 - 2 * np.cos(np.radians(degrees))) / (3 * (np.cos(np.radians(degrees)) - 1) - np.sqrt(3) * np.sin(np.radians(degrees))))
        else:
            n = np.abs(np.cos(np.radians(degrees)) / (np.cos(np.radians(degrees)) - np.sin(np.radians(degrees)) - 1))
        
        u1 = n * (self.vectors[0] + self.vectors[1]) - self.vectors[0]
        u2 = n * (self.vectors[0] + self.vectors[1]) - self.vectors[1]
        
        for ax, deg, title in zip(axes, rotations, titles):
            if self.lattice_type in ["Hexagon", "Kagome", "Triangle"]:
                shift_x = self.lattice_distance * 0.5
                shift_y = self.lattice_distance * np.sqrt(3) / 4
            else:
                shift_x = 0
                shift_y = 0.5*self.lattice_distance

            
            if self.lattice_type == "Hexagon":
                p11, p12 = self.generate_lattice_points()
                p21, p22 = self.generate_rotated_points(deg)
                
                ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=5, label="a1")
                ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=5, label="b1")
                ax.scatter(p21[:, 0] - shift_x, p21[:, 1] - shift_y, color=layer_2, s=5, label="a2")
                ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=5, label="b2")
            
            elif self.lattice_type == "Kagome":
                p11, p12, p13 = self.generate_lattice_points()
                p21, p22, p23 = self.generate_rotated_points(deg)
                
                ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=5, label="a1")
                ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=5, label="b1")
                ax.scatter(p13[:, 0], p13[:, 1], color=layer_1, s=5, label="c1")
                
                ax.scatter(p21[:, 0] - shift_x, p21[:, 1] - shift_y, color=layer_2, s=5, label="a2")
                ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=5, label="b2")
                ax.scatter(p23[:, 0] - shift_x, p23[:, 1] - shift_y, color=layer_2, s=5, label="c2")
            
            elif self.lattice_type == "Triangle":
                p11 = self.generate_lattice_points()
                p22 = self.generate_rotated_points(deg)
                
                ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=5, label="a1")
                ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=5, label="b2")
            
            elif self.lattice_type == "Lieb":
                p11, p12, p13 = self.generate_lattice_points()
                p21, p22, p23 = self.generate_rotated_points(deg)
                
                ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=5, label="Layer 1")
                ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=5)
                ax.scatter(p13[:, 0], p13[:, 1], color=layer_1, s=5)
                
                ax.scatter(p21[:, 0] - shift_x, p21[:, 1] - shift_y, color=layer_2, s=5, label="Layer 2")
                ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=5)
                ax.scatter(p23[:, 0] - shift_x, p23[:, 1] - shift_y, color=layer_2, s=5)
            
            else:
                p11 = self.generate_lattice_points()
                p22 = self.generate_rotated_points(deg)
                
                ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=5, label="a1")
                ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=10, label="b2")
            
            ax.quiver(-self.vectors[0], -self.vectors[0], u1[0], u1[1], color='r', angles='xy', scale_units='xy', scale=1, label='u1')
            ax.quiver(-self.vectors[1], -self.vectors[1], u2[0], u2[1], color='g', angles='xy', scale_units='xy', scale=1, label='u2')
            
            ax.set_title(title)
            ax.axis("equal")
            ax.legend()
        
        if save:
            plt.savefig(f"{self.lattice_type}_bilayer_CW_CCW_{degrees}.pdf")
        
        plt.show()    

    
    def plot_bilayer_align_comparison(self, degrees, save=False):
        """ Plots 2D Bilayer of a given lattice, top-down view, showing both CW and CCW rotations."""
        layer_1 = (0.1, 0.1, 0.6, 0.5)
        layer_2 = (0.6, 0.1, 0.1, 0.5)
        
        fig, axes = plt.subplots(1, 2, figsize=(20, 10))
        rotations = [degrees, -degrees]
        titles = [f"{self.lattice_type} CCW ({degrees}°)", f"{self.lattice_type} CW ({-degrees}°)"]
        
        if self.lattice_type in ["Hexagon", "Kagome", "Triangle"]:
            n = np.abs((1 - 2 * np.cos(np.radians(degrees))) / (3 * (np.cos(np.radians(degrees)) - 1) - np.sqrt(3) * np.sin(np.radians(degrees))))
        else:
            n = np.abs(np.cos(np.radians(degrees)) / (np.cos(np.radians(degrees)) - np.sin(np.radians(degrees)) - 1))
        
        u1 = n * (self.vectors[0] + self.vectors[1]) + self.vectors[0]
        u2 = n * (self.vectors[0] + self.vectors[1]) + self.vectors[1]
        
        for ax, deg, title in zip(axes, rotations, titles):
            if self.lattice_type in ["Hexagon", "Kagome", "Triangle"]:
                shift_x = 0
                shift_y = 0
            else:
                shift_x = 0
                shift_y = 0
            
            if self.lattice_type == "Hexagon":
                p11, p12 = self.generate_lattice_points()
                p21, p22 = self.generate_rotated_points(deg)
                
                ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=5, label="a1")
                ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=5, label="b1")
                ax.scatter(p21[:, 0] - shift_x, p21[:, 1] - shift_y, color=layer_2, s=5, label="a2")
                ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=5, label="b2")
            
            elif self.lattice_type == "Kagome":
                p11, p12, p13 = self.generate_lattice_points()
                p21, p22, p23 = self.generate_rotated_points(deg)
                
                ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=5, label="a1")
                ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=5, label="b1")
                ax.scatter(p13[:, 0], p13[:, 1], color=layer_1, s=5, label="c1")
                
                ax.scatter(p21[:, 0] - shift_x, p21[:, 1] - shift_y, color=layer_2, s=5, label="a2")
                ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=5, label="b2")
                ax.scatter(p23[:, 0] - shift_x, p23[:, 1] - shift_y, color=layer_2, s=5, label="c2")
            
            elif self.lattice_type == "Triangle":
                p11 = self.generate_lattice_points()
                p22 = self.generate_rotated_points(deg)
                
                ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=5)
                ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=5)
            
            elif self.lattice_type == "Lieb":
                p11, p12, p13 = self.generate_lattice_points()
                p21, p22, p23 = self.generate_rotated_points(deg)
                
                ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=5, label="Layer 1")
                ax.scatter(p12[:, 0], p12[:, 1], color=layer_1, s=5)
                ax.scatter(p13[:, 0], p13[:, 1], color=layer_1, s=5)
                
                ax.scatter(p21[:, 0] - shift_x, p21[:, 1] - shift_y, color=layer_2, s=5, label="Layer 2")
                ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=5)
                ax.scatter(p23[:, 0] - shift_x, p23[:, 1] - shift_y, color=layer_2, s=5)
            
            else:
                p11 = self.generate_lattice_points()
                p22 = self.generate_rotated_points(deg)
                
                ax.scatter(p11[:, 0], p11[:, 1], color=layer_1, s=5, label="a1")
                ax.scatter(p22[:, 0] - shift_x, p22[:, 1] - shift_y, color=layer_2, s=10, label="b2")
            
            ax.quiver(0, 0, u1[0], u1[1], color='r', angles='xy', scale_units='xy', scale=1, label='u1')
            ax.quiver(0, 0, u2[0], u2[1], color='g', angles='xy', scale_units='xy', scale=1, label='u2')
            
            ax.set_title(title)
            ax.axis("equal")
            ax.legend()
        
        if save:
            plt.savefig(f"{self.lattice_type}_bilayer_CW_CCW_{degrees}.pdf")
        
        plt.show()    

    def plot_reciprocal(self):
        """ plots reciprocal lattice strucutre """
        if self.lattice_type=="Hexagon":
            points, points2 = self.generate_reciprocal_points()
        else:
            points = self.generate_reciprocal_points()

        plt.figure(figsize=(8, 8))
        plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50)

        if len(self.basis) > 1:
            plt.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50)
        plt.title(f'{self.lattice_type} Reciprocal Lattice')
        plt.axis('equal')
        plt.show()
    
    def plot_bilayer(self, degrees, point_size=5, save=False):
        """ Plots 3D Bilayer of a given lattice, not aligned. The Second layer is shifted"""
        a = 0.5 
        layer_1 = (0.1,0.1,0.6,0.5)
        layer_2 = (0.6,0.1,0.1,0.5)
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection="3d")
        if self.lattice_type == "Hexagon":
            p11,p12 = self.generate_lattice_points()
            p21,p22 = self.generate_rotated_points(degrees)

            ax.scatter(p11[:, 0], p11[:, 1],0, color=layer_1, s=point_size, label="a1")
            ax.scatter(p12[:, 0], p12[:, 1],0, color=layer_1, s=5, label="b1")

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*np.sqrt(3)/4

            ax.scatter(p21[:, 0]-lattice_shift_x, p21[:, 1]-lattice_shift_y,a, color=layer_2, s=5, label="a2")
            ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=layer_2, s=5, label="b2")

        elif self.lattice_type == "Kagome":
            p11, p12, p13 = self.generate_lattice_points()
            p21, p22, p23 = self.generate_rotated_points(degrees)

            ax.scatter(p11[:, 0], p11[:, 1],0, color=layer_1, s=5, label="a1")
            ax.scatter(p12[:, 0], p12[:, 1],0, color=layer_1, s=5, label="b1")
            ax.scatter(p13[:, 0], p13[:, 1],0, color=layer_1, s=5, label="c1")

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*np.sqrt(3)/4

            ax.scatter(p21[:, 0]-lattice_shift_x, p21[:, 1]-lattice_shift_y,a, color=layer_2, s=5, label="a2")
            ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=layer_2, s=5, label="b2")
            ax.scatter(p23[:, 0]-lattice_shift_x, p23[:, 1]-lattice_shift_y,a, color=layer_2, s=5, label="c2")

        
        elif self.lattice_type == "Triangle":
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)
            
            ax.scatter(p11[:, 0], p11[:, 1],0, color=layer_1, s=5, label="a1")

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*np.sqrt(3)/4

            ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=layer_2, s=5, label="b2")

        elif self.lattice_type == "Lieb":
            p11, p12, p13 = self.generate_lattice_points()
            p21, p22, p23 = self.generate_rotated_points(degrees)

            ax.scatter(p11[:, 0], p11[:, 1],0, color=layer_1, s=5, label="Layer 1")
            ax.scatter(p12[:, 0], p12[:, 1],0, color=layer_1, s=5)
            ax.scatter(p13[:, 0], p13[:, 1],0, color=layer_1, s=5)

            lattice_shift_x = self.lattice_distance
            lattice_shift_y = self.lattice_distance

            ax.scatter(p21[:, 0]-lattice_shift_x, p21[:, 1]-lattice_shift_y,a, color=layer_2, s=5, label="Layer 2")
            ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=layer_2, s=5)
            ax.scatter(p23[:, 0]-lattice_shift_x, p23[:, 1]-lattice_shift_y,a, color=layer_2, s=5)

             
        else:
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)

            ax.scatter(p11[:, 0], p11[:, 1],0, color=layer_1, s=5, label="a1")

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*0.5

            ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=layer_2, s=10, label="b2")

        ax.set_zlim(-a, a*2)
        ax.set_title(f"3D Bilayer {self.lattice_type}")
        ax.axis('equal')
        ax.legend()
        ax.view_init(90,90,180)
        if(save):
            plt.savefig(f'{self.lattice_type}_bilayer_{degrees}.pdf')
        plt.show()
    
    def plot_bilayer_align(self, degrees, save=False):
        """bilayer plotting with no shift (atoms a1 on a2, b1 on b2)"""
        a = 0.5
        layer_1 = (0.1,0.1,0.6,0.5)
        layer_2 = (0.6,0.1,0.1,0.5)
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection="3d")

        if self.lattice_type == "Hexagon":
            p11,p12 = self.generate_lattice_points()
            p21,p22 = self.generate_rotated_points(degrees)
        
            ax.scatter(p11[:, 0], p11[:, 1],0, color=layer_1, s=5, label="a1")
            ax.scatter(p12[:, 0], p12[:, 1],0, color=layer_1, s=5, label="b1")

            ax.scatter(p21[:, 0], p21[:, 1],a, color=layer_2, s=5, label="a2")
            ax.scatter(p22[:, 0], p22[:, 1],a, color=layer_2, s=5, label="b2")
        
        elif self.lattice_type == "Kagome":
            p11, p12, p13 = self.generate_lattice_points()
            p21, p22, p23 = self.generate_rotated_points(degrees)

            ax.scatter(p11[:, 0], p11[:, 1],0, color=layer_1, s=5, label="a1")
            ax.scatter(p12[:, 0], p12[:, 1],0, color=layer_1, s=5, label="b1")
            ax.scatter(p13[:, 0], p13[:, 1],0, color=layer_1, s=5, label="c1")

            ax.scatter(p21[:, 0], p21[:, 1],a, color=layer_2, s=5, label="a2")
            ax.scatter(p22[:, 0], p22[:, 1],a, color=layer_2, s=5, label="b2")
            ax.scatter(p23[:, 0], p23[:, 1],a, color=layer_2, s=5, label="c2")

        elif self.lattice_type == "Lieb":

            p11, p12, p13 = self.generate_lattice_points()
            p21, p22, p23 = self.generate_rotated_points(degrees)

            ax.scatter(p11[:, 0], p11[:, 1],0, color=layer_1, s=5, label="a1")
            ax.scatter(p12[:, 0], p12[:, 1],0, color=layer_1, s=5, label="b1")
            ax.scatter(p13[:, 0], p13[:, 1],0, color=layer_1, s=5, label="c1")

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*np.sqrt(3)/4

            ax.scatter(p21[:, 0]-lattice_shift_x, p21[:, 1]-lattice_shift_y,a, color=layer_2, s=5, label="a2")
            ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=layer_2, s=5, label="b2")
            ax.scatter(p23[:, 0]-lattice_shift_x, p23[:, 1]-lattice_shift_y,a, color=layer_2, s=5, label="c2")

        else:
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)

            ax.scatter(p11[:, 0], p11[:, 1],0, color=layer_1, s=5, label="a1")

            ax.scatter(p22[:, 0], p22[:, 1],a, color=layer_2, s=5, label="a2")
            

        ax.set_zlim(-a, a*2)
        ax.set_title("3D Bilayer Graphene")
        ax.axis('equal')
        ax.legend()
        ax.view_init(90,90,180)
        if(save):
            plt.savefig(f'{self.lattice_type}_aligned_bilayer_{degrees}.pdf')
        plt.show()

    def plot_lattice_with_twist_Circles(self, save=False):
        """ 
        plots 2d lattice structure with Circles of radius n*(a_1+a_2)+a_1 
        using lattice vectors. number of Circles determined by self.num_sites.
        """
        if self.lattice_type == "Hexagon":
            points, points2 = self.generate_lattice_points()
            points3 = None
        elif self.lattice_type == "Kagome":
            points, points2, points3 = self.generate_lattice_points()
        else:
            points = self.generate_lattice_points()
            points2 = None
            points3 = None

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50, label="lattice points")
        if points2 is not None:
            ax.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50, label="basis points")
        if points3 is not None:
            ax.scatter(points3[:, 0], points3[:, 1], color=(0.2, 0.5, 0.2, 0.5), s=50, label="basis points")

        a1 = self.vectors[0]
        a2 = self.vectors[1]

        origin = np.array([0, 0])
        for n in range(0, self.num_sites):  
            radius_vector = n * (np.array(a1) + np.array(a2)) + np.array(a1)
            radius = np.linalg.norm(radius_vector)
            Circle = plt.Circle(origin, radius, color="r", fill=False, alpha=0.3, lw=1.5)
            ax.add_artist(Circle)

        max_extent = max(np.linalg.norm(n * (a1 + a2) + a1) for n in range(0, self.num_sites + 1)) + self.lattice_distance
        ax.set_xlim(-max_extent, max_extent)
        ax.set_ylim(-max_extent, max_extent)
        ax.set_aspect('equal', adjustable='datalim')  

        plt.title(f'{self.lattice_type} Lattice with Twist Circles')
        plt.legend()
        plt.axis('equal')
        if save:
            plt.savefig(f'{self.lattice_type}_twist_Circles.pdf')
        plt.show()

    def plot_lattice_with_twist_vectors(self, degrees, save=False):
        """ 
        plots 2d lattice structure with reciprocal lattice vectors b1 and b2 
        for a twist defined by the rotation angle.
        """
        phi = np.radians(degrees)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        
        if self.lattice_type == "Square":
            n = np.round(-cos_phi / (cos_phi - sin_phi - 1))
            factor = 2 * np.pi / (2 * n**2 + 2 * n + 1)
            b1 = factor * np.array([(n + 1), n])
            b2 = factor * np.array([-n, (n + 1)])
        elif self.lattice_type == "Triangle":
            n = np.round((1 - 2 * cos_phi) / (3 * (cos_phi - 1) - np.sqrt(3) * sin_phi))
            factor = 2 * np.pi / (3 * n**2 + 3 * n + 1)
            b1 = factor * np.array([(2 * n + 1), -1 / np.sqrt(3)])
            b2 = factor * np.array([-n, (3 * n + 2) / np.sqrt(3)])
        else:
            raise ValueError("twist vectors not defined for this lattice type.")

        points = self.generate_lattice_points()
        if isinstance(points, tuple):  
            points, points2 = points
        else:
            points2 = None

        plt.figure(figsize=(8, 8))
        plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50, label="lattice points")
        if points2 is not None:
            plt.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50, label="basis points")

        a1 = self.vectors[0]
        a2 = self.vectors[1]
        
        origin = np.array([0, 0])
        plt.quiver(origin[0], origin[1], a1[0], a1[1], angles='xy', scale_units='xy', scale=1, width=0.005)
        plt.quiver(origin[0], origin[1], a2[0], a2[1], angles='xy', scale_units='xy', scale=1, width=0.005)
        plt.quiver(origin[0], origin[1], b1[0], b1[1], angles='xy', scale_units='xy', scale=1, color='r', width=0.005, label="b1")
        plt.quiver(origin[0], origin[1], b2[0], b2[1], angles='xy', scale_units='xy', scale=1, color='g', width=0.005, label="b2")

        for i in range(self.num_sites):  
            radius_vector = i * (np.array(a1) + np.array(a2)) + np.array(a1)
            radius = np.linalg.norm(radius_vector)
            Circle = plt.Circle(origin, radius, color="r", fill=False, alpha=0.3, lw=1.5)
            plt.gca().add_artist(Circle)

        max_extent = max(np.linalg.norm(i * (b1 + b2)) for i in range(self.num_sites)) + self.lattice_distance
        plt.xlim(-max_extent, max_extent)
        plt.ylim(-max_extent, max_extent)
        plt.gca().set_aspect('equal', adjustable='datalim')  

        plt.title(f'{self.lattice_type} Lattice with Twist Vectors (Rotation: {degrees}°)')
        plt.legend()
        plt.axis('equal')
        if save:
            plt.savefig(f'{self.lattice_type}_twist_vectors_{degrees}.pdf')
        plt.show()

    def plot_reciprocal_difference(self, save=False):
        """ 
        animates the superlattice reciprocal lattice structure with vectors b1 and b2 
        as the twist angle spans from 0 to 90 degrees (phi: 0 to pi/2).
        """
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_title(f'{self.lattice_type} Superlattice Reciprocal Comparison')

        points_scatter = ax.scatter([], [], color=(0.1, 0.2, 0.5, 0.5), s=15, label="superlattice points")
        quiver_a1 = ax.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1, width=0.005, label="original reciprocal b1")
        quiver_a2 = ax.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1, width=0.005, label="original reciprocal b2")
        quiver_b1 = ax.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1, color='r', width=0.005, label="b1")
        quiver_b2 = ax.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1, color='g', width=0.005, label="b2")

        if self.lattice_type == "Square":
            phi_values = np.linspace(0.001, np.pi / 2, 200)
        elif self.lattice_type == "Triangle":
            phi_values = np.linspace(0.001, np.pi / 3, 200)
        b1, b2 = self.reciprocal_vectors

        def update(frame):
            phi = phi_values[frame]
            cos_phi = np.cos(phi)
            sin_phi = np.sin(phi)

            if self.lattice_type == "Square":
                n = -cos_phi / (cos_phi - sin_phi - 1)
                factor = 2 * np.pi / (2 * n**2 + 2 * n + 1)
                b1 = factor * np.array([(n + 1), n])
                b2 = factor * np.array([-n, (n + 1)])

                bz_edges = self.generate_bzedges(b1,b2)
                bz2_edges = self.generate_bzedges(b1,b2)
                
                for i in range(4):
                    start = bz_edges[i]
                    end = bz_edges[(i + 1) % 4]  
                    ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="orange", label="bz edge" if i == 0 else "")

                for i in range(4):
                    start = bz2_edges[i]
                    end = bz2_edges[(i + 1) % 4]  
                    ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="orange", label="bz2 edge" if i == 0 else "")

            elif self.lattice_type == "Triangle":
                n = (1 - 2 * cos_phi) / (3 * (cos_phi - 1) - np.sqrt(3) * sin_phi)
                factor = 2 * np.pi / (3 * n**2 + 3 * n + 1)
                b1 = factor * np.array([(2 * n + 1), -1 / np.sqrt(3)])
                b2 = factor * np.array([-n, (3 * n + 2) / np.sqrt(3)])

                bz_edges = self.generate_bzedges(b1,b2)
                bz2_edges = self.generate_bzedges(b1,b2)
                
                for i in range(len(bz_edges)):
                    start = bz_edges[i]
                    end = bz_edges[(i + 1) % len(bz_edges)]  
                    ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="orange", label="bz edge" if i == 0 else "")

                for i in range(len(bz2_edges)):
                    start = bz2_edges[i]
                    end = bz2_edges[(i + 1) % len(bz2_edges)]  
                    ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="orange", label="bz2 edge" if i == 0 else "")

                
            else:
                raise ValueError("superlattice reciprocal not defined for this lattice type.")

            points = []
            for i in range(self.x_range[0], self.x_range[1]):
                for j in range(self.y_range[0], self.y_range[1]):
                    displacement = i * b1 + j * b2 
                    points.append(displacement)
            points = np.array(points)

            points_scatter.set_offsets(points)

            quiver_a1.set_uvc(b1[0], b1[1])
            quiver_a2.set_uvc(b2[0], b2[1])
            quiver_b1.set_uvc(b1[0], b1[1])
            quiver_b2.set_uvc(b2[0], b2[1])
            
            ax.set_title(f'{self.lattice_type} Superlattice Reciprocal Comparison (Phi = {np.degrees(phi):.2f}°)')
            return points_scatter, quiver_b1, quiver_b2

        max_extent = self.lattice_distance * self.num_sites
        ax.set_xlim(-max_extent, max_extent)
        ax.set_ylim(-max_extent, max_extent)
        ax.legend(loc="lower left")

        anim = FuncAnimation(fig, update, frames=len(phi_values), blit=False)

        if save:
            anim.save(f'{self.lattice_type}_superlattice_reciprocal_difference_animation.gif', fps=20)
        else:
            plt.show()

    def plot_bz_difference(self, degrees, s=1, save=False):
        """ plot to showcase the difference in brillouin zone for the twisted system """
        b1 = self.reciprocal_vectors[0]
        b2 = self.reciprocal_vectors[1]

        a1 = self.vectors[0]
        a2 = self.vectors[1]
        angle = np.radians(degrees)
        rot_matrix = np.array([[np.cos(-angle), -np.sin(-angle)], [np.sin(-angle), np.cos(-angle)]])
        rota1 = rot_matrix @ a1
        rota2 = rot_matrix @ a2

        rotb = self.calculate_reciprocal(rota1,rota2)

        rotb1 = rotb[0]
        rotb2 = rotb[1]

        cos_phi = np.cos(angle)
        sin_phi = np.sin(angle)

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_title(f'{self.lattice_type} Brillouin Zone Pattern')

        if self.lattice_type == "Square":
            n = (-(s*cos_phi) / (cos_phi - sin_phi - 1))
            factor = 2 * np.pi / (2 * n**2 + 2 * n * s+ s**2)
            b1 = factor * np.array([n, (n+s)])
            b2 = factor * np.array([-n-s, (n)])
        elif self.lattice_type == "Triangle":
            n = s*abs((1 - 2 * cos_phi) / (3 * (cos_phi - 1) - np.sqrt(3) * sin_phi))
            factor = (2 * np.pi) / (3 * (n**2) + 3 * n * s + s**2)
            b1 = factor * np.array([(2 * n + s), -s / np.sqrt(3)])
            b2 = factor * np.array([-n, (3 * n + 2*s) / np.sqrt(3)])
        else:
            raise ValueError("twist vectors not defined for this lattice type.")
        
        grid_size = int(n+1)

        if self.lattice_type =="Square":
            
            bz_edges = self.generate_bzedges(self.reciprocal_vectors[0],self.reciprocal_vectors[1])
            bz2_edges = self.generate_bzedges(rotb1,rotb2)
            bz3_edges = self.generate_bzedges(b1,b2)
                    
            for i in range(4):
                        start = bz_edges[i]
                        end = bz_edges[(i + 1) % 4]  
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.1, 0.2, 0.5, 0.8), label="original bz" if i == 0 else "")

            for i in range(4):
                        start = bz2_edges[i]
                        end = bz2_edges[(i + 1) % 4]  
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.5, 0.1, 0.2, 0.8), label="rotated bz" if i == 0 else "")

            for i in range(grid_size-2*s, grid_size+2*s):
                for j in range(-int(grid_size/n)-4, -int(grid_size/n)+4):
                    translation = i * b1 + j * b2
                    for k in range(len(bz3_edges)):
                        start = bz3_edges[k] + translation
                        end = bz3_edges[(k + 1) % len(bz3_edges)] + translation
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="k", linewidth=0.25 , label="superlattice bz" if i == -grid_size and j == -grid_size and k == 0 else "")
        
        elif self.lattice_type=="Triangle":

            bz_edges = self.generate_bzedges(self.reciprocal_vectors[0],self.reciprocal_vectors[1])
            bz2_edges = self.generate_bzedges(rotb1,rotb2)
            bz3_edges = self.generate_bzedges(b1,b2)
            
            for i in range(len(bz_edges)):
                        start = bz_edges[i]
                        end = bz_edges[(i + 1) % len(bz_edges)]  
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="k", label="original bz" if i == 0 else "")

            for i in range(len(bz2_edges)):
                        start = bz2_edges[i]
                        end = bz2_edges[(i + 1) % len(bz2_edges)]  
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.5, 0.1, 0.2, 0.8), label="rotated bz" if i == 0 else "")

            
            for i in range(grid_size, grid_size+s):
                for j in range(grid_size, grid_size+s):
                    translation = i * b1 + j * b2
                    for k in range(len(bz3_edges)):
                        start = bz3_edges[k] + translation
                        end = bz3_edges[(k + 1) % len(bz3_edges)] + translation
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="b", linewidth=0.25 , label="superlattice bz" if i == -grid_size and j == -grid_size and k == 0 else "")
        plt.legend()
        if save:
            plt.savefig(f"{self.lattice_type}_bz_diff.pdf")
        plt.show()

    def plot_bz_difference_cw(self, degrees, s=1, save=False):
        """ plot to showcase the difference in brillouin zone for the twisted system """
        b1 = self.reciprocal_vectors[0]
        b2 = self.reciprocal_vectors[1]

        a1 = self.vectors[0]
        a2 = self.vectors[1]
        angle = np.radians(degrees)
        rot_matrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        rota1 = rot_matrix @ a1
        rota2 = rot_matrix @ a2

        rotb = self.calculate_reciprocal(rota1,rota2)

        rotb1 = rotb[0]
        rotb2 = rotb[1]

        cos_phi = np.cos(angle)
        sin_phi = np.sin(angle)

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_title(f'{self.lattice_type} Brillouin Zone Pattern')

        if self.lattice_type == "Square":
            n = (-(s*cos_phi) / (cos_phi - sin_phi - 1))
            factor = 2 * np.pi / (2 * n**2 + 2 * n * s+ s**2)
            b1_u2 = factor * np.array([(n + s), n])
            b2_u2 = factor * np.array([-n, (n + s)])
        elif self.lattice_type == "Triangle":
            n = np.round(s*((np.sqrt(3)-2*sin_phi) / (np.sqrt(3) * (cos_phi - 1) + 3 * sin_phi)))
            factor = (2 * np.pi) / (3 * (n**2) + 3 * n * s + s**2)
            b1_u2 = factor * np.array([(2 * n + s), s / np.sqrt(3)])
            b2_u2 = factor * np.array([-n-s, (3 * n + s) / np.sqrt(3)])
        else:
            raise ValueError("twist vectors not defined for this lattice type.")
        
        grid_size = int(n+2)

        if self.lattice_type =="Square":
            
            bz_edges = self.generate_bzedges(self.reciprocal_vectors[0],self.reciprocal_vectors[1])
            bz2_edges = self.generate_bzedges(rotb1,rotb2)
            bz3_edges = self.generate_bzedges(b1_u2,b2_u2)
                    
            for i in range(4):
                        start = bz_edges[i]
                        end = bz_edges[(i + 1) % 4]  
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.1, 0.2, 0.5, 0.8), label="original bz" if i == 0 else "")

            for i in range(4):
                        start = bz2_edges[i]
                        end = bz2_edges[(i + 1) % 4]  
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.5, 0.1, 0.2, 0.8), label="rotated bz" if i == 0 else "")

            for i in range(grid_size-2*s, grid_size+2*s):
                for j in range(-int(grid_size/n), -int(grid_size/n)+5):
                    translation = i * b1_u2 + j * b2_u2
                    for k in range(len(bz3_edges)):
                        start = bz3_edges[k] + translation
                        end = bz3_edges[(k + 1) % len(bz3_edges)] + translation
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="k", linewidth=0.25 , label="superlattice bz" if i == -grid_size and j == -grid_size and k == 0 else "")
        
        elif self.lattice_type=="Triangle":

            bz_edges = self.generate_bzedges(self.reciprocal_vectors[0],self.reciprocal_vectors[1])
            bz2_edges = self.generate_bzedges(rotb1,rotb2)
            bz3_edges = self.generate_bzedges(b1_u2,b2_u2)
            
            for i in range(len(bz_edges)):
                        start = bz_edges[i]
                        end = bz_edges[(i + 1) % len(bz_edges)]  
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="k", label="original bz" if i == 0 else "")

            for i in range(len(bz2_edges)):
                        start = bz2_edges[i]
                        end = bz2_edges[(i + 1) % len(bz2_edges)]  
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.5, 0.1, 0.2, 0.8), label="rotated bz" if i == 0 else "")

            
            for i in range(grid_size-s, grid_size+s):
                for j in range(grid_size-s, grid_size+s):
                    translation = i * b1_u2 + j * b2_u2
                    for k in range(len(bz3_edges)):
                        start = bz3_edges[k] + translation
                        end = bz3_edges[(k + 1) % len(bz3_edges)] + translation
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.1, 0.1, 0.6, 0.6), linewidth=0.25 , label="Superlattice BZ" if i == -grid_size and j == -grid_size and k == 0 else "")

        plt.legend()
        if save:
            plt.savefig(f"{self.lattice_type}_bz_diff.pdf")
        plt.show()
    
    def plot_bilayer_twist_animation(self, save=False):
        """ animation to showcase the twist of a bilayer """
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection="3d")
        if self.lattice_type == "Hexagon":
            return 0
        else:
            max_angle = 90
            
            p11 = self.generate_lattice_points()
            a = 1

            ax.scatter(p11[:, 0], p11[:, 1], 0, color=(0.1, 0.2, 0.5, 0.5), s=5, label="a1 (bottom layer)")
            
            top_layer_plot, = ax.plot([], [], [], 'o', color=(0.1, 0.5, 0.2, 0.5), markersize=5, label="a2 (top layer)")
            ax.axis('equal')
            ax.set_zlim(-a, a * 2)
            ax.set_title(f"{self.lattice_type} Bilayer Rotation Animation")
            ax.legend()
            ax.view_init(90,90,180)

            def update(frame):
                angle = frame  
                p22 = self.generate_rotated_points(angle)
                top_layer_plot.set_data(p22[:, 0], p22[:, 1])
                top_layer_plot.set_3d_properties(a)  
                ax.set_title(f"Rotation angle: {angle:.1f}°")
                return top_layer_plot,

            anim = FuncAnimation(fig, update, frames=np.linspace(0, max_angle, 100), blit=True)

            if save:
                anim.save(f'{self.lattice_type}_rotation_animation.mp4', fps=15)
            else:
                plt.show()

    def __str__(self):
        """ generic string output """
        return "this is a "+self.lattice_type+" lattice"
    
