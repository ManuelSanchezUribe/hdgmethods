import numpy as np
import scipy.sparse as sps
from scipy.sparse.linalg import spsolve
from scipy.integrate import quadrature





def hdg1d_Poisson_hybridization(mesh, p=1, f=None, c=None):
    iA = []; jA = []; kA = [];
    iR = []; jR = []; kR = [];
    # precompute the binomial coefficients
    global Cb 
    Cb = Compute_binomial_coefficients(max(2*p,p))
    # Load vector
    b = np.zeros(mesh.NN);
    Rb = np.zeros(2*mesh.NE*(p+1))
    # General Structure of degrees of freedom: [Nodes, faces, interior]
    for j in range(mesh.NE):
        # Interval
        K = mesh.Coordinates[mesh.Elements[j,:]]
    
        # Local Matrices
        Bi = Convective1d(K, p)
        Mi = Mass1d(K, p, c)
        Si = np.zeros((p+1,p+1)); Si[0,0] = 1.0; Si[p,p] = 1.0
        Ei = np.zeros((2,p+1))  ; Ei[0,0] = 1.0; Ei[1,p] = 1.0
        Ci = np.zeros((2,p+1))  ; Ci[0,0] = -1.0; Ci[1,p] = 1.0
        Gi = np.eye(2)
        
        # Local load vector
        if f is None:
            Fi = np.zeros(p+1)
        else:
            Fi = Compute_mu(K, p, f)
        
        # condensation
        Qp = np.vstack(( np.hstack(( Mi, -Bi.T)),  np.hstack(( Bi, Si)) ))
        CEi = np.hstack(( Ci, Ei ))
        Ri = np.linalg.solve(Qp, np.vstack(( Ci.T, -Ei.T )))
        Ai =  -CEi @ Ri - Gi
        Rbi = np.linalg.solve(Qp, np.hstack((np.zeros(p+1), Fi )))
        bi = -CEi @ Rbi
        
        # Global Assembling
        dofn = mesh.Elements[j,:]
        # assembling b
        b[dofn] += bi
        # Sparse assembling of A: A[dof,dof] = AK
        rows = dofn; cols = dofn
        prod = [(x,y) for x in rows for y in cols];
        iAlocal = [x for (x,y) in prod]
        jAlocal = [y for (x,y) in prod]
        kAlocal = np.reshape(Ai, (1,Ai.size))
        iA.extend(iAlocal); jA.extend(jAlocal); kA.extend(kAlocal.flatten())
        
        # Reconstruction
        # interior degrees of freedom
        dofq = np.arange((p+1)*(j), (p+1)*(j+1)).tolist()
        dofu = np.arange(mesh.NE*(p+1)+(p+1)*(j), mesh.NE*(p+1)+(p+1)*(j+1)).tolist()
        dofi = dofq+dofu
        Rb[dofi] = Rbi
        # Sparse assembling of A: A[dof,dof] = AK
        rows = dofi; cols = dofn
        prod = [(x,y) for x in rows for y in cols];
        iRlocal = [x for (x,y) in prod]
        jRlocal = [y for (x,y) in prod]
        kRlocal = np.reshape(Ri, (1,Ri.size))
        iR.extend(iRlocal); jR.extend(jRlocal); kR.extend(kRlocal.flatten())
    # end for
    A = sps.coo_matrix(( kA,(iA, jA)), (mesh.NN, mesh.NN), dtype=np.float64)
    R = sps.coo_matrix(( kR,(iR, jR)), (2*mesh.NE*(p+1), mesh.NN), dtype=np.float64)
    #fig, axa = plt.subplots(1,1,figsize=(4,4))
    #axa.spy(A, markersize=5)
    #plt.show()
    return (A, b, R, Rb)
def hdg1d_solver(mesh, p=1, f=None, c=None, uD=None, hybridization=True):
    # Solve linear system
    if hybridization:
        (A, b, R, Rb) = hdg1d_Poisson_hybridization(mesh, p, f, c)
        uhat = np.zeros(mesh.NN);
        uhat[0] = uD(mesh.Coordinates[0])
        uhat[mesh.NN-1] = uD(mesh.Coordinates[-1])
        b = b- A.dot(uhat)
        All = range(0,mesh.NN)
        Free = np.setdiff1d(All,[0,mesh.NN-1]);
        AFree = A.tocsc()[Free, :][:, Free]
        bFree = b[Free]
        uhat[Free] = spsolve(AFree, bFree)
        
        qu = Rb-R.dot(uhat)
        qh = qu[0:mesh.NE*(p+1)]
        uh = qu[mesh.NE*(p+1):mesh.NE*2*(p+1)]
        return qh, uh, uhat