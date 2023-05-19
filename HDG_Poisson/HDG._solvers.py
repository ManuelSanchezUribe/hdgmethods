import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import spsolve
from scipy.integrate import quadrature

'''
HDG Solver class for Poisson equation in 1d
c * q   = grad(u)
-div(q) = f

'''
class HDG_solver:
    def __init__(self, mesh1d=None, poisson_problem=None, order=None):
        self.problem = poisson_problem
        self.p = order
        self.mesh = mesh1d

    def add_problem(self, poisson_problem=None):
        self.problem = poisson_problem
    
    def HDG_spaces(self):
        self.Wh = None
        self.Vh = None
        self.Mh = None
        #  we need Gauss Lobatto basis
    
    def Solver(self):

        for j in range(self.mesh.NE):
            # local element
            K = self.mesh.Coordinates[self.mesh.Elements[j,:]]

            # local system
            QK, bK = self.LocalSystem()


    def LocalSystem(self, K, cj, fj):

    
def massmatrix(self, c=None):

from scipy.special import eval_jacobi

def Vandermonde1d(npoints,r):
    V1d = np.zeros(r.size, npoints+1)

    x, w = GaussLobatto_Points_Weights(npoints)

    for j in range(npoints):
        V1d[:,j] = eval_jacobi(x, 0,0 ,j)
    return V1d