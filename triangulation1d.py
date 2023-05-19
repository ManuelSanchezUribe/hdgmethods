import numpy as np

# Triangulation
class mesh1d:
    def __init__(self, Coordinates, Elements):
        self.Coordinates = Coordinates
        self.Elements = Elements
        self.NE = Elements.shape[0] # number of elements
        self.NN = Coordinates.size # number of coordinates or nodes
        self.h = np.max([abs(Coordinates[j+1] - Coordinates[j]) for j in range(self.NE)]) # parameter of the triangulation

    def equispaced_mesh(xi,xf,j):
        n = 2**j;
        Coordinates = np.linspace(xi,xf,num = n+2) # nodes 
        Elements = np.vstack([range(n+1), range(1,n+2)]).T
        return mesh1d(Coordinates, Elements)