import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class Lattice:
    def __init__(self, lattice_distance, vectors, basis, num_sites):
        """ Initialization of Lattice parameters """
        self.lattice_distance = lattice_distance
        self.vectors = self.lattice_distance*np.array(vectors)
        self.basis = self.lattice_distance*np.array(basis)
        self.lattice_type = self.determine_lattice_type()
        self.reciprocal_vectors = self.get_Reciprocal()
        self.NN = self.get_NN()
        self.num_sites = num_sites
        self.x_range = (-num_sites, num_sites)
        self.y_range = self.x_range

    def get_Reciprocal(self):
        """ Solves the reciprocal lattice with the intialized vectors """
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
    
    def determine_lattice_type(self):
        """ Using the angle the vectors produce the lattice type is determined """
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
        
    def get_NN(self):
        """ Generates the nearest neighbors of the 2D lattice strucutre """
        nearestNeighbors = []
        for basis in self.basis:
            for vector in self.vectors:
                nearestNeighbors.append(basis-vector)
                nearestNeighbors.append(basis+vector)
        return np.array(nearestNeighbors)

    def generate_lattice_points(self):
        """ Generate Lattice Points to use for plotting, takes into account the basis """
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
    
    def generate_rotated_points(self, angle):
        """ Degree rotation function for later use when discussing super lattices"""
        angle = angle *(np.pi/180)
        rot_matrix = np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
        points = []
        points2 = []

        rot_vectors = []
        for a in self.vectors:
            rot_vectors.append(rot_matrix.dot(a))

        for i in range(self.x_range[0], self.x_range[1]):
            for j in range(self.y_range[0], self.y_range[1]):
                displacement = i * rot_vectors[0] + j * rot_vectors[1]
                points.append(displacement + self.basis[0])
                if len(self.basis) > 1:
                    points2.append(displacement + self.basis[1])

        points = np.array(points)
        if points2:
            points2 = np.array(points2)
            return points, points2
        else:
            return points

    def generate_reciprocal_points(self):
        """ Generates lattice points with reciprocal vectors"""
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
        """ Generates reciprocal lattice points for the superlattice defined by the rotation angle """
        phi = np.radians(degrees)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        
        # Determine n based on lattice type
        if self.lattice_type == "Square":
            n = np.round(-cos_phi / (cos_phi - sin_phi - 1))
            factor = 2 * np.pi / (2 * n**2 + 2 * n + 1)
            B1 = factor * np.array([(n + 1), n])
            B2 = factor * np.array([-n, (n + 1)])
        elif self.lattice_type == "Triangle":
            n = np.round((1 - 2 * cos_phi) / (3 * (cos_phi - 1) - np.sqrt(3) * sin_phi))
            factor = 2 * np.pi / (3 * n**2 + 3 * n + 1)
            B1 = factor * np.array([(2 * n + 1), -1 / np.sqrt(3)])
            B2 = factor * np.array([-n, (3 * n + 2) / np.sqrt(3)])
        else:
            raise ValueError("Superlattice reciprocal not defined for this lattice type.")
        
        # Generate points
        points = []
        for i in range(self.x_range[0], self.x_range[1]):
            for j in range(self.y_range[0], self.y_range[1]):
                displacement = i * B1 + j * B2
                points.append(displacement)
    
        return np.array(points)

    def plot_superlattice_reciprocal(self, degrees, save=False):
        """ Plots the reciprocal lattice structure for the superlattice """
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
        Plots the reciprocal lattice structure for the superlattice 
        with the reciprocal lattice vectors B1 and B2 displayed.
        """
        phi = np.radians(degrees)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        
        # Determine n based on lattice type and calculate B1, B2
        if self.lattice_type == "Square":
            n = np.round(-cos_phi / (cos_phi - sin_phi - 1))
            factor = 2 * np.pi / (2 * n**2 + 2 * n + 1)
            B1 = factor * np.array([(n + 1), n])
            B2 = factor * np.array([-n, (n + 1)])
        elif self.lattice_type == "Triangle":
            n = np.round((1 - 2 * cos_phi) / (3 * (cos_phi - 1) - np.sqrt(3) * sin_phi))
            factor = 2 * np.pi / (3 * n**2 + 3 * n + 1)
            B1 = factor * np.array([(2 * n + 1), -1 / np.sqrt(3)])
            B2 = factor * np.array([-n, (3 * n + 2) / np.sqrt(3)])
        else:
            raise ValueError("Superlattice reciprocal not defined for this lattice type.")
        
        # Generate reciprocal points
        points = self.generate_superlattice_reciprocal_points(degrees)

        # Plot the points
        plt.figure(figsize=(8, 8))
        plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50, label="Reciprocal Points")

        # Add vectors B1 and B2 using quiver
        origin = np.array([0, 0])
        plt.quiver(origin[0], origin[1], B1[0], B1[1], angles='xy', scale_units='xy', scale=1, color='r', width=0.005, label="B1")
        plt.quiver(origin[0], origin[1], B2[0], B2[1], angles='xy', scale_units='xy', scale=1, color='g', width=0.005, label="B2")

        # Configure the plot
        plt.title(f'{self.lattice_type} Superlattice Reciprocal with Vectors (Rotation: {degrees}°)')
        plt.axis('equal')
        plt.legend()
        if save:
            plt.savefig(f'{self.lattice_type}_superlattice_reciprocal_vectors_{degrees}.pdf')
        plt.show()


    def plot_lattice(self):
        """ Plots 2D lattice structure """
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
        plt.show()

    def plot_reciprocal(self):
        """ Plots reciprocal lattice strucutre """
        if self.lattice_type=="Hexagon":
            points, points2 = self.generate_reciprocal_points()
        else:
            points = self.generate_reciprocal_points()

        plt.figure(figsize=(8, 8))
        plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50)

        if len(self.basis) > 1:
            plt.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50)
        plt.title(f'{self.lattice_type} reciprocal Lattice')
        plt.axis('equal')
        plt.show()
    
    def plot_bilayer(self, degrees, save=False):
        """ Plots 3D Bilayer Hexagon lattices (Atoms A1 and B2 overlap)"""
        a = 0.5 #vertical distance between layers
        if self.lattice_type == "Hexagon":
            p11,p12 = self.generate_lattice_points()
            p21,p22 = self.generate_rotated_points(degrees)
            

            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111, projection="3d")

            ax.scatter(p11[:, 0], p11[:, 1],0, color=(0.1, 0.2, 0.5, 0.5), s=5, label="A1")
            ax.scatter(p12[:, 0], p12[:, 1],0, color=(0.5, 0.1, 0.2, 0.5), s=5, label="B1")

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*np.sqrt(3)/4

            ax.scatter(p21[:, 0]-lattice_shift_x, p21[:, 1]-lattice_shift_y,a, color=(0.1, 0.5, 0.2,0.5), s=5, label="A2")
            ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=(0.5, 0.5, 0.2,0.5), s=5, label="B2")
            
                    
            ax.set_zlim(-a, a*2)
            ax.set_title("3D Bilayer Graphene")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,90,180)
            if(save):
                plt.savefig(f'{self.lattice_type}_bilayer_{degrees}.pdf')
            plt.show()
        elif self.lattice_type == "Triangle":
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)
            

            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111, projection="3d")

            ax.scatter(p11[:, 0], p11[:, 1],0, color=(0.1, 0.2, 0.5, 0.5), s=5, label="A1")

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*np.sqrt(3)/4

            ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=(0.5, 0.5, 0.2,0.5), s=5, label="B2")
            
                    
            ax.set_zlim(-a, a*2)
            ax.set_title("Triangular Bilayer")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,90,180)
            if(save):
                plt.savefig(f'{self.lattice_type}_bilayer_{degrees}.pdf')
            plt.show()
        else:
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)
            

            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111, projection="3d")

            ax.scatter(p11[:, 0], p11[:, 1],0, color=(0.1, 0.2, 0.5, 0.5), s=5, label="A1")

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*0.5

            ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=(0.5, 0.5, 0.2,0.5), s=10, label="B2")
            
                    
            ax.set_zlim(-a, a*2)
            ax.set_title("Square Bilayer")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,90,180)
            if(save):
                plt.savefig(f'{self.lattice_type}_bilayer_{degrees}.pdf')
            plt.show()
    
    def plot_bilayer_align(self, degrees, save=False):
        """Bilayer Plotting with no shift (Atoms A1 on A2, B1 on B2)"""
        a = 0.5
        if self.lattice_type == "Hexagon":
            p11,p12 = self.generate_lattice_points()
            p21,p22 = self.generate_rotated_points(degrees)
            

            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111, projection="3d")

            ax.scatter(p11[:, 0], p11[:, 1],0, color=(0.1, 0.2, 0.5, 0.5), s=5, label="A1")
            ax.scatter(p12[:, 0], p12[:, 1],0, color=(0.5, 0.1, 0.2, 0.5), s=5, label="B1")

            ax.scatter(p21[:, 0], p21[:, 1],a, color=(0.1, 0.5, 0.2,0.5), s=5, label="A2")
            ax.scatter(p22[:, 0], p22[:, 1],a, color=(0.5, 0.5, 0.2,0.5), s=5, label="B2")
            
                    
            ax.set_zlim(-a, a*2)
            ax.set_title("3D Bilayer Graphene")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,90,180)
            if(save):
                plt.savefig(f'{self.lattice_type}_aligned_bilayer_{degrees}.pdf')
            plt.show()
        else:
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)

            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111, projection="3d")

            ax.scatter(p11[:, 0], p11[:, 1],0, color=(0.1, 0.2, 0.5, 0.5), s=5, label="A1")

            ax.scatter(p22[:, 0], p22[:, 1],a, color=(0.1, 0.5, 0.2,0.5), s=5, label="A2")
            
                    
            ax.set_zlim(-a, a*2)
            ax.set_title(f"{self.lattice_type} Graphene Aligned")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,90,180)
            if(save):
                plt.savefig(f'{self.lattice_type}_aligned_bilayer_{degrees}.pdf')
            plt.show()

    def plot_lattice_with_twist_circles(self, save=False):
        """ 
        Plots 2D lattice structure with circles of radius n*(a_1+a_2)+a_1 
        using lattice vectors. Number of circles determined by self.num_sites.
        """
        # Determine lattice points
        if self.lattice_type == "Hexagon":
            points, points2 = self.generate_lattice_points()
        else:
            points = self.generate_lattice_points()
            points2 = None

        # Initialize plot
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50, label="Lattice Points")
        if points2 is not None:
            ax.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50, label="Basis Points")

        # Get lattice vectors
        a1 = self.vectors[0]
        a2 = self.vectors[1]

        # Draw circles for twist angles
        origin = np.array([0, 0])
        for n in range(0, self.num_sites):  # Circle range from 0 to num_sites
            radius_vector = n * (np.array(a1) + np.array(a2)) + np.array(a1)
            radius = np.linalg.norm(radius_vector)
            circle = plt.Circle(origin, radius, color="r", fill=False, alpha=0.3, lw=1.5)
            ax.add_artist(circle)

        # Set axis limits and aspect ratio
        max_extent = max(np.linalg.norm(n * (a1 + a2) + a1) for n in range(0, self.num_sites + 1)) + self.lattice_distance
        ax.set_xlim(-max_extent, max_extent)
        ax.set_ylim(-max_extent, max_extent)
        ax.set_aspect('equal', adjustable='datalim')  # Ensure equal aspect ratio

        # Finalize plot
        plt.title(f'{self.lattice_type} Lattice with Twist Circles')
        plt.legend()
        plt.axis('equal')
        if save:
            plt.savefig(f'{self.lattice_type}_twist_circles.pdf')
        plt.show()

    def plot_lattice_with_twist_vectors(self, degrees, save=False):
        """ 
        Plots 2D lattice structure with reciprocal lattice vectors B1 and B2 
        for a twist defined by the rotation angle.
        """
        phi = np.radians(degrees)
        cos_phi = np.cos(phi)
        sin_phi = np.sin(phi)
        
        # Determine n based on lattice type and calculate B1, B2
        if self.lattice_type == "Square":
            n = np.round(-cos_phi / (cos_phi - sin_phi - 1))
            factor = 2 * np.pi / (2 * n**2 + 2 * n + 1)
            B1 = factor * np.array([(n + 1), n])
            B2 = factor * np.array([-n, (n + 1)])
        elif self.lattice_type == "Triangle":
            n = np.round((1 - 2 * cos_phi) / (3 * (cos_phi - 1) - np.sqrt(3) * sin_phi))
            factor = 2 * np.pi / (3 * n**2 + 3 * n + 1)
            B1 = factor * np.array([(2 * n + 1), -1 / np.sqrt(3)])
            B2 = factor * np.array([-n, (3 * n + 2) / np.sqrt(3)])
        else:
            raise ValueError("Twist vectors not defined for this lattice type.")

        # Generate lattice points
        points = self.generate_lattice_points()
        if isinstance(points, tuple):  # Handle cases with basis points
            points, points2 = points
        else:
            points2 = None

        # Initialize plot
        plt.figure(figsize=(8, 8))
        plt.scatter(points[:, 0], points[:, 1], color=(0.1, 0.2, 0.5, 0.5), s=50, label="Lattice Points")
        if points2 is not None:
            plt.scatter(points2[:, 0], points2[:, 1], color=(0.5, 0.1, 0.2, 0.5), s=50, label="Basis Points")

        a1 = self.vectors[0]
        a2 = self.vectors[1]
        
        origin = np.array([0, 0])
        plt.quiver(origin[0], origin[1], a1[0], a1[1], angles='xy', scale_units='xy', scale=1, width=0.005)
        plt.quiver(origin[0], origin[1], a2[0], a2[1], angles='xy', scale_units='xy', scale=1, width=0.005)
        plt.quiver(origin[0], origin[1], B1[0], B1[1], angles='xy', scale_units='xy', scale=1, color='r', width=0.005, label="B1")
        plt.quiver(origin[0], origin[1], B2[0], B2[1], angles='xy', scale_units='xy', scale=1, color='g', width=0.005, label="B2")

        # Draw circles for twist visualization
        for i in range(self.num_sites):  # Circle range from 0 to num_sites
            radius_vector = i * (np.array(a1) + np.array(a2)) + np.array(a1)
            radius = np.linalg.norm(radius_vector)
            circle = plt.Circle(origin, radius, color="r", fill=False, alpha=0.3, lw=1.5)
            plt.gca().add_artist(circle)

        # Set axis limits and aspect ratio
        max_extent = max(np.linalg.norm(i * (B1 + B2)) for i in range(self.num_sites)) + self.lattice_distance
        plt.xlim(-max_extent, max_extent)
        plt.ylim(-max_extent, max_extent)
        plt.gca().set_aspect('equal', adjustable='datalim')  # Ensure equal aspect ratio

        # Finalize plot
        plt.title(f'{self.lattice_type} Lattice with Twist Vectors (Rotation: {degrees}°)')
        plt.legend()
        plt.axis('equal')
        if save:
            plt.savefig(f'{self.lattice_type}_twist_vectors_{degrees}.pdf')
        plt.show()

    def plot_reciprocal_difference(self, save=False):
        """ 
        Animates the superlattice reciprocal lattice structure with vectors B1 and B2 
        as the twist angle spans from 0 to 90 degrees (phi: 0 to pi/2).
        """
        # Prepare figure
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.set_title(f'{self.lattice_type} Superlattice Reciprocal Comparison')

        # Initialize plot elements
        points_scatter = ax.scatter([], [], color=(0.1, 0.2, 0.5, 0.5), s=15, label="Superlattice Points")
        quiver_a1 = ax.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1, width=0.005, label="Original Reciprocal b1")
        quiver_a2 = ax.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1, width=0.005, label="Original Reciprocal b2")
        quiver_b1 = ax.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1, color='r', width=0.005, label="B1")
        quiver_b2 = ax.quiver(0, 0, 0, 0, angles='xy', scale_units='xy', scale=1, color='g', width=0.005, label="B2")

        # Angle range: 0 to 90 degrees
        if self.lattice_type == "Square":
            phi_values = np.linspace(0.001, np.pi / 2, 200)
        elif self.lattice_type == "Triangle":
            phi_values = np.linspace(0.001, np.pi / 3, 200)
        b1, b2 = self.reciprocal_vectors

        def update(frame):
            phi = phi_values[frame]
            cos_phi = np.cos(phi)
            sin_phi = np.sin(phi)

            # Determine n and calculate B1, B2
            if self.lattice_type == "Square":
                n = -cos_phi / (cos_phi - sin_phi - 1)
                factor = 2 * np.pi / (2 * n**2 + 2 * n + 1)
                B1 = factor * np.array([(n + 1), n])
                B2 = factor * np.array([-n, (n + 1)])

                # Correcting the square BZ vertices calculation
                bz_edges = [
                    0.5 * (b1+b2),   
                    0.5 * (b1-b2), 
                    -0.5 * (b1+b2), 
                    -0.5 * (b1-b2) 
                ]
                bz2_edges = [
                    0.5 * (B1+B2),   
                    0.5 * (B1-B2), 
                    -0.5 * (B1+B2), 
                    -0.5 * (B1-B2) 
                ]
                
                # Plot the square BZ
                for i in range(4):
                    start = bz_edges[i]
                    end = bz_edges[(i + 1) % 4]  # Connect back to the first point to form the square
                    ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="orange", label="BZ Edge" if i == 0 else "")

                for i in range(4):
                    start = bz2_edges[i]
                    end = bz2_edges[(i + 1) % 4]  # Connect back to the first point to form the square
                    ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="orange", label="BZ2 Edge" if i == 0 else "")

            elif self.lattice_type == "Triangle":
                n = (1 - 2 * cos_phi) / (3 * (cos_phi - 1) - np.sqrt(3) * sin_phi)
                factor = 2 * np.pi / (3 * n**2 + 3 * n + 1)
                B1 = factor * np.array([(2 * n + 1), -1 / np.sqrt(3)])
                B2 = factor * np.array([-n, (3 * n + 2) / np.sqrt(3)])

                # Define the 6 hexagonal directions by rotating a1
                angles = np.linspace(0, 2 * np.pi, 7)  # 6 directions
                b_length = np.linalg.norm(b1+b2)
                bz_edges = [(np.sqrt(3)/3) * b_length * np.array([np.cos(angle), np.sin(angle)]) for angle in angles]

                bz2_edges = [
                    (2/3 * B1[0] + 1/3 * B2[0], 2/3 * B1[1] + 1/3 * B2[1]),  # Vertex 2
                    (1/3 * B1[0] + 2/3 * B2[0], 1/3 * B1[1] + 2/3 * B2[1]),  # Vertex 1
                    
                    (-1/3 * B1[0] + 1/3 * B2[0], -1/3 * B1[1] + 1/3 * B2[1]), # Vertex 3
                    (-2/3 * B1[0] - 1/3 * B2[0], -2/3 * B1[1] - 1/3 * B2[1]), # Vertex 4
                    (-1/3 * B1[0] - 2/3 * B2[0], -1/3 * B1[1] - 2/3 * B2[1]), # Vertex 5
                    (1/3 * B1[0] - 1/3 * B2[0], 1/3 * B1[1] - 1/3 * B2[1])    # Vertex 6
                ]
                
                
                # Plot the hexagon edges
                for i in range(len(bz_edges)):
                    start = bz_edges[i]
                    end = bz_edges[(i + 1) % len(bz_edges)]  # Connect back to the first point
                    ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="orange", label="BZ Edge" if i == 0 else "")

                for i in range(len(bz2_edges)):
                    start = bz2_edges[i]
                    end = bz2_edges[(i + 1) % len(bz2_edges)]  # Connect back to the first point
                    ax.plot([start[0], end[0]], [start[1], end[1]], '-', color="orange", label="BZ2 Edge" if i == 0 else "")

                
            else:
                raise ValueError("Superlattice reciprocal not defined for this lattice type.")

            # Generate superlattice reciprocal points
            points = []
            for i in range(self.x_range[0], self.x_range[1]):
                for j in range(self.y_range[0], self.y_range[1]):
                    displacement = i * B1 + j * B2 
                    points.append(displacement)
            points = np.array(points)

            # Update scatter points
            points_scatter.set_offsets(points)

            # Update quivers
            quiver_a1.set_UVC(b1[0], b1[1])
            quiver_a2.set_UVC(b2[0], b2[1])
            quiver_b1.set_UVC(B1[0], B1[1])
            quiver_b2.set_UVC(B2[0], B2[1])
            
            # Update title with current angle
            ax.set_title(f'{self.lattice_type} Superlattice Reciprocal Comparison (Phi = {np.degrees(phi):.2f}°)')
            return points_scatter, quiver_b1, quiver_b2

        # Set plot limits
        max_extent = self.lattice_distance * self.num_sites
        ax.set_xlim(-max_extent, max_extent)
        ax.set_ylim(-max_extent, max_extent)
        ax.legend(loc="lower left")

        # Create animation
        anim = FuncAnimation(fig, update, frames=len(phi_values), blit=False)

        # Save or show the animation
        if save:
            anim.save(f'{self.lattice_type}_superlattice_reciprocal_difference_animation.gif', fps=20)
        else:
            plt.show()

    def plot_BZ_difference(self, degrees, save=False):
        a1 = self.reciprocal_vectors[0]
        a2 = self.reciprocal_vectors[1]
        angle = np.radians(degrees)
        rot_matrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
        rota1 = rot_matrix.dot(a1)
        rota2 = rot_matrix.dot(a2)

        cos_phi = np.cos(angle)
        sin_phi = np.sin(angle)

        
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set_title(f'{self.lattice_type} Brillouin Zone Pattern')

        if self.lattice_type == "Square":
            n = np.round(-cos_phi / (cos_phi - sin_phi - 1))
            factor = 2 * np.pi / (2 * n**2 + 2 * n + 1)
            B1 = factor * np.array([(n + 1), n])
            B2 = factor * np.array([-n, (n + 1)])
        elif self.lattice_type == "Triangle":
            n = np.round((1 - 2 * cos_phi) / (3 * (cos_phi - 1) - np.sqrt(3) * sin_phi))
            factor = 2 * np.pi / (3 * n**2 + 3 * n + 1)
            B1 = factor * np.array([(2 * n + 1), -1 / np.sqrt(3)])
            B2 = factor * np.array([-n, (3 * n + 2) / np.sqrt(3)])
        else:
            raise ValueError("Twist vectors not defined for this lattice type.")
        
        grid_size = int(np.round(n))+1

        if self.lattice_type =="Square":
            
            bz_edges = [
                    0.5 * (a1+a2),   
                    0.5 * (a1-a2), 
                    -0.5 * (a1+a2), 
                    -0.5 * (a1-a2) 
                ]
            bz2_edges = [
                        0.5 * (rota1+rota2),   
                        0.5 * (rota1-rota2), 
                        -0.5 * (rota1+rota2), 
                        -0.5 * (rota1-rota2) 
                    ]
            bz3_edges = [
                0.5 * (B1+B2),  
                    0.5 * (B1-B2), 
                    -0.5 * (B1+B2), 
                    -0.5 * (B1-B2) 
            ]
                    
                    # Plot the square BZ
            for i in range(4):
                        start = bz_edges[i]
                        end = bz_edges[(i + 1) % 4]  # Connect back to the first point to form the square
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.1, 0.2, 0.5, 0.8), label="Original BZ" if i == 0 else "")

            for i in range(4):
                        start = bz2_edges[i]
                        end = bz2_edges[(i + 1) % 4]  # Connect back to the first point to form the square
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.5, 0.1, 0.2, 0.8), label="Rotated BZ" if i == 0 else "")

            for i in range(-grid_size, grid_size + 1):
                for j in range(-grid_size, grid_size + 1):
                    translation = i * B1 + j * B2
                    for k in range(len(bz3_edges)):
                        start = bz3_edges[k] + translation
                        end = bz3_edges[(k + 1) % len(bz3_edges)] + translation
                        ax.plot([start[0], end[0]], [start[1], end[1]], '-', color=(0.1, 0.5, 0.5, 0.5), label="Superlattice BZ" if i == -grid_size and j == -grid_size and k == 0 else "")

        ax.set_aspect('equal', 'box')
        plt.legend()
        plt.show()

    
    def plot_bilayer_twist_animation(self, save=False):
        if self.lattice_type == "Hexagon":
            return 0
        else:
            max_angle = 90
            # Set up figure and 3D plot
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111, projection="3d")

            p11 = self.generate_lattice_points()
            a = 1

            # Scatter plot for the bottom layer (fixed)
            ax.scatter(p11[:, 0], p11[:, 1], 0, color=(0.1, 0.2, 0.5, 0.5), s=5, label="A1 (Bottom Layer)")
            
            # Plot variable to store rotated points
            top_layer_plot, = ax.plot([], [], [], 'o', color=(0.1, 0.5, 0.2, 0.5), markersize=5, label="A2 (Top Layer)")

            # Set plot limits and labels
            ax.axis('equal')
            ax.set_zlim(-a, a * 2)
            ax.set_title(f"{self.lattice_type} Bilayer Rotation Animation")
            ax.legend()
            ax.view_init(90,90,180)

            # Animation function
            def update(frame):
                angle = frame  # Current rotation angle in degrees
                p22 = self.generate_rotated_points(angle)
                # Update the data for the top layer plot
                top_layer_plot.set_data(p22[:, 0], p22[:, 1])
                top_layer_plot.set_3d_properties(a)  # Set constant z for top layer
                ax.set_title(f"Rotation Angle: {angle:.1f}°")
                return top_layer_plot,

            # Create the animation
            anim = FuncAnimation(fig, update, frames=np.linspace(0, max_angle, 100), blit=True)

            # Save or show the animation
            if save:
                anim.save(f'{self.lattice_type}_rotation_animation.mp4', fps=15)
            else:
                plt.show()

    def __str__(self):
        """ Generic string output """
        return "This is a "+self.lattice_type+" lattice"
    