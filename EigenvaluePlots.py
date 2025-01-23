import numpy as np
import matplotlib.pyplot as plt

def eiganvalue(x,y):

    def euler(val_x, val_y):
        exp = 1j*(val_x*x + val_y*y)
        return (np.e**(exp))
    
    mult_term = euler(-1/2,-np.sqrt(3)/4)
    sqrt_term = np.sqrt(euler(1/2,0) + euler(3/2,0) + euler(0,np.sqrt(3)/2) + 3*euler(1, np.sqrt(3)/2) + euler(2, np.sqrt(3)/2) + euler(1/2, np.sqrt(3)) + euler(3/2, np.sqrt(3)))
    
    term = abs(mult_term * sqrt_term)

    return term

x = np.linspace(-4, 4, 100)
y = np.linspace(-4, 4, 100)
  
X, Y = np.meshgrid(x, y)
Z = eiganvalue(X, Y)
fig = plt.figure()
ax = plt.axes(projection ='3d')
ax.plot_surface(X, Y, Z, color =(0.2,0.2,0.5,1))
ax.plot_surface(X, Y, Z+2, color =(0.2,0.2,0.5,1))
ax.plot_surface(X, Y, -Z, color = (0.5,0.2,0.2,1)) 
ax.plot_surface(X, Y, -(Z+2), color =(0.5,0.2,0.2,1))
plt.title("Eigenvalues of Bilayer Graphene")
plt.show()