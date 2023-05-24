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
    
    def HDG_spaces(self, order):
        self.Wh = []
        self.Vh = []
        self.Mh = []

        # dimension on each element
        self.Vh['localdim'] = order+1
        self.Wh['localdim'] = order+1
        self.Mh['localdim'] = 2
        # Global dimension of the space
        self.Vh['dim'] = self.mesh.NE * self.Vh['localdim']
        self.Wh['dim'] = self.mesh.NE * self.Wh['localdim']
        self.Mh['dim'] = self.mesh.NN
        
        #  we need Gauss Lobatto basis
    
    def Solver(self, problem):

        b = np.zeros(self.Mh['dim'], dtype=np.float64)
        A = np.zeros((self.Mh['dim'], self.Mh['dim']), dtype=np.float64)
        for j in range(self.mesh.NE):
            # local element
            dof = self.mesh.Elements[j,:]
            K = self.mesh.Coordinates[dof]

            # local system
            AK, bK, _ = self.LocalSystem(K, cK, fK)

            # Global Assembling
            b[dof] += bK
            A[dof,dof] += AK 

        # Solve
        uhat = np.zeros(self.Mh['dim'])
        # Impose Dirichlet boundary condition
        Dirichletdof = [0, self.Mh['dim']-1]
        for j in Dirichletdof:
            uhat[j] = problem.uD(self.mesh.Coordinates[j])
        b -= A.dot(uhat)
        Freedof = np.setdiff1d(range(0, self.Mh['dim']), Dirichletdof)
        AFree = A[Freedof,:][:,Freedof]
        bFree = b[Freedof]
        uhat[Freedof] = np.linalg.solve(AFree, bFree)



    def LocalSystem(self, K, cj, fj):
        # Local Matricres for order = 0
        hK = np.abs(K[1]-K[0])
        # B: (uh, div(r))
        B = np.zeros((self.Vh['localdim'], self.Wh['localdim']))
        # M: (qh, r)
        M = np.zeros((self.Vh['localdim'], self.Vh['localdim']))
        Mhat = 
        M = 0.5*hK*Mhat
        # S: ( uh, v)
        S = np.zeros((self.Wh['localdim'], self.Wh['localdim']))
        # E: (uh, mu)
        E = np.zeros((self.Mh['localdim'], self.Wh['localdim']))
        # C: (qh*n, mu)
        C = np.zeros((self.Mh['localdim'], self.Vh['localdim']))
        # G: (uhat, mu)
        G = np.zeros((self.Mh['localdim'], self.Mh['localdim']))

        Q = np.vstack(( np.hstack(( M, -B.T)),  np.hstack(( B, S)) ))
        CE = np.hstack(( C, E ))
        R = np.linalg.solve(Q, np.vstack(( C.T, -E.T )))
        A =  -CE @ R - G
        Rb = np.linalg.solve(Q, np.hstack((self.Vh['localdim'], F )))
        b = -CE @ Rb
        return A, b, Q


    
def massmatrix(self, c=None):



