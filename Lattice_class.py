import numpy as np
import matplotlib.pyplot as plt
import math

class Lattice:
    def __init__(self, lattice_distance, vectors, basis, num_sites):
        """ Initialization of Lattice parameters """
        self.lattice_distance = lattice_distance
        self.vectors = self.lattice_distance*np.array(vectors)
        self.basis = self.lattice_distance*np.array(basis)
        self.lattice_type = self.determine_lattice_type()
        self.reciprocal_vectors = self.get_Reciprocal()
        self.NN = self.get_NN()
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
    
    def plot_bilayer(self, degrees, save):
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
            lattice_shift_y = self.lattice_distance*np.sqrt(3)/6

            ax.scatter(p21[:, 0]-lattice_shift_x, p21[:, 1]-lattice_shift_y,a, color=(0.1, 0.5, 0.2,0.5), s=5, label="A2")
            ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=(0.5, 0.5, 0.2,0.5), s=5, label="B2")
            
                    
            ax.set_zlim(-a, a*2)
            ax.set_title("3D Bilayer Graphene")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,0,90)
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
            lattice_shift_y = self.lattice_distance*np.sqrt(3)/6

            ax.scatter(p22[:, 0]-lattice_shift_x, p22[:, 1]-lattice_shift_y,a, color=(0.5, 0.5, 0.2,0.5), s=5, label="B2")
            
                    
            ax.set_zlim(-a, a*2)
            ax.set_title("Triangular Bilayer")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,0,90)
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
            ax.view_init(90,0,90)
            if(save):
                plt.savefig(f'{self.lattice_type}_bilayer_{degrees}.pdf')
            plt.show()
    
    def plot_bilayer_align(self, degrees, save):
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
            ax.view_init(90,0,90)
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
            ax.view_init(90,0,90)
            if(save):
                plt.savefig(f'{self.lattice_type}_aligned_bilayer_{degrees}.pdf')
            plt.show()

    def plot_bilayer_rotation_locus(self, degrees, save):
        a = 0.5 #vertical distance between layers
        if self.lattice_type == "Hexagon":
            p11,p12 = self.generate_lattice_points()
            p21,p22 = self.generate_rotated_points(degrees)
            

            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111, projection="3d")

            ax.scatter(p11[:, 0], p11[:, 1],0, color=(0.1, 0.2, 0.5, 0.5), s=5, label="A1")
            ax.scatter(p12[:, 0], p12[:, 1],0, color=(0.5, 0.1, 0.2, 0.5), s=5, label="B1")

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*np.sqrt(3)/6

            p21[:,0] = p21[:,0]-lattice_shift_x
            p21[:,1] = p21[:,1]-lattice_shift_y

            p22[:,0] = p22[:,0]-lattice_shift_x
            p22[:,1] = p22[:,1]-lattice_shift_y

            ax.scatter(p21[:, 0], p21[:, 1],a, color=(0.1, 0.5, 0.2,0.5), s=5, label="A2")
            ax.scatter(p22[:, 0], p22[:, 1],a, color=(0.5, 0.5, 0.2,0.5), s=5, label="B2")
            
            theta = np.linspace(0, np.pi/2, 20)
            for point in p22:
                radius = np.linalg.norm(point) 
                circle_x = radius * np.cos(theta)
                circle_y = radius * np.sin(theta)

                ax.plot(circle_x, circle_y, a, 'r', alpha=0.1)
                    
            ax.set_zlim(-a, a*2)
            ax.set_title("3D Bilayer Graphene with Rotations")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,0,90)
            if(save):
                plt.savefig(f'{self.lattice_type}_bilayer_{degrees}_rotation_locus.pdf')
            plt.show()
        elif self.lattice_type == "Triangle":
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)
            

            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111, projection="3d")

            ax.scatter(p11[:, 0], p11[:, 1],0, color=(0.1, 0.2, 0.5, 0.5), s=5, label="A1")

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*np.sqrt(3)/6

            p22[:,0] = p22[:,0]-lattice_shift_x
            p22[:,1] = p22[:,1]-lattice_shift_y

            ax.scatter(p22[:, 0], p22[:, 1],a, color=(0.5, 0.5, 0.2,0.5), s=5, label="B2")
            
            theta = np.linspace(0, np.pi/2, 20)
            for point in p22:
                radius = np.linalg.norm(point) 
                circle_x = radius * np.cos(theta)
                circle_y = radius * np.sin(theta)

                ax.plot(circle_x, circle_y, a, 'r', alpha=0.1)

            ax.set_zlim(-a, a*2)
            ax.set_title("Triangular Bilayer with Rotations")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,0,90)
            if(save):
                plt.savefig(f'{self.lattice_type}_bilayer_{degrees}_rotation_locus.pdf')
            plt.show()
        else:
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)
            

            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111, projection="3d")

            ax.scatter(p11[:, 0], p11[:, 1],0, color=(0.1, 0.2, 0.5, 0.5), s=5, label="A1")

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*0.5

            p22[:,0] = p22[:,0]-lattice_shift_x
            p22[:,1] = p22[:,1]-lattice_shift_y

            ax.scatter(p22[:, 0], p22[:, 1],a, color=(0.5, 0.5, 0.2,0.5), s=5, label="B2")
            
            theta = np.linspace(0, np.pi/2, 20)
            for point in p22:
                radius = np.linalg.norm(point) 
                circle_x = radius * np.cos(theta)
                circle_y = radius * np.sin(theta)

                ax.plot(circle_x, circle_y, a, 'r', alpha=0.1)
                
            ax.set_zlim(-a, a*2)
            ax.set_title("Square Bilayer with Rotations")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,0,90)
            if(save):
                plt.savefig(f'{self.lattice_type}_bilayer_{degrees}_rotation_locus.pdf')
            plt.show()

    def plot_aligned_bilayer_rotation_locus(self, degrees, save):
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

            theta = np.linspace(0, np.pi/2, 20)
            for point in p22:
                radius = np.linalg.norm(point) 
                circle_x = radius * np.cos(theta)
                circle_y = radius * np.sin(theta)

                ax.plot(circle_x, circle_y, a, 'r', alpha=0.1)
            
                    
            ax.set_zlim(-a, a*2)
            ax.set_title("3D Bilayer Graphene with rotations (Aligned)")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,0,90)
            if(save):
                plt.savefig(f'{self.lattice_type}_aligned_bilayer_{degrees}_rotation_locus.pdf')
            plt.show()
        else:
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)

            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111, projection="3d")

            ax.scatter(p11[:, 0], p11[:, 1],0, color=(0.1, 0.2, 0.5, 0.5), s=5, label="A1")
            ax.scatter(p22[:, 0], p22[:, 1],a, color=(0.1, 0.5, 0.2,0.5), s=5, label="A2")

            theta = np.linspace(0, np.pi/2, 20)
            for point in p22:
                radius = np.linalg.norm(point) 
                circle_x = radius * np.cos(theta)
                circle_y = radius * np.sin(theta)

                ax.plot(circle_x, circle_y, a, 'r', alpha=0.1)
                    
            ax.set_zlim(-a, a*2)
            ax.set_title(f"{self.lattice_type} Graphene Aligned with Rotations (Aligned)")
            ax.axis('equal')
            ax.legend()
            ax.view_init(90,0,90)
            if(save):
                plt.savefig(f'{self.lattice_type}_aligned_bilayer_{degrees}_rotation_locus.pdf')
            plt.show()

    """ To Do Later
    def find_aligned_overlap_points(self, degrees, tolerance=0.04):
        Finds overlapping points in the bilayer structure within a given tolerance.
        if self.lattice_type == "Hexagon":
            p11,p12 = self.generate_lattice_points()
            p21,p22 = self.generate_rotated_points(degrees)

            
        else:
            overlap_points = []
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)

            for point in p22:
                for point2 in p11:
                    if np.linalg.norm(point - point2)<tolerance:
                        overlap_points.append(point)
            
            return np.array(overlap_points)
    
    def find_overlap_points(self, degrees, tolerance=0.04):

        if self.lattice_type == "Hexagon":
            p11,p12 = self.generate_lattice_points()
            p21,p22 = self.generate_rotated_points(degrees)

            
        else:
            overlap_points = []
            p11 = self.generate_lattice_points()
            p22 = self.generate_rotated_points(degrees)

            lattice_shift_x = self.lattice_distance*0.5
            lattice_shift_y = self.lattice_distance*0.5

            p11[:,0]+=lattice_shift_x
            p11[:,1]+=lattice_shift_y

            for point in p22:
                for point2 in p11:
                    if math.dist(point,point2)<tolerance:
                        overlap_points.append(point)
            
            return np.array(overlap_points)
    
    """
    
    def __str__(self):
        """ Generic string output """
        return "This is a "+self.lattice_type+" lattice"
    