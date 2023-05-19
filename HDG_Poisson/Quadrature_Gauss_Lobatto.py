import numpy as np
from scipy.special import eval_jacobi

def GaussLobatto_Points_Weights(npoints):
    '''
    Input  : number of Gauss Lobatto points
    Output : x, Points and w, Weights
    '''
    if npoints == 3:
        x = np.array([-1.0,0,1.0], dtype=np.float64)
        w = np.array([1.0/3.0, 4.0/3.0, 1.0/3.0], dtype=np.float64)
    elif npoints == 4:
        a = np.sqrt(1.0/5.0)
        x = np.array([-1.0,-a,a,1.0], dtype=np.float64)
        w = np.array([1.0/6.0, 5.0/6.0, 5.0/6.0, 1.0/6.0], dtype=np.float64)
    elif npoints == 5:
        a = np.sqrt(3.0/7.0)
        x = np.array([-1.0,-a,0,a,1.0], dtype=np.float64)
        w = np.array([1.0/10.0, 49.0/90.0,32.0/45.0, 49.0/90.0, 1.0/10.0], dtype=np.float64)
    elif npoints == 6:
        a = np.sqrt(1.0/3.0 - 2.0*np.sqrt(7.0)/21.0)
        b = np.sqrt(1.0/3.0 + 2.0*np.sqrt(7.0)/21.0)
        x = np.array([-1.0, -b, -a, a, b, 1.0], dtype=np.float64)
        aw = (14.0+np.sqrt(7.0))/30.0
        bw = (14.0-np.sqrt(7.0))/30.0
        w = np.array([1.0/15.0, bw, aw, aw, bw, 1.0/15.0], dtype=np.float64)
    elif npoints == 7:
        a = np.sqrt(5/11.0 - 2.0*np.sqrt(5.0/3.0)/11.0)
        b = np.sqrt(5/11.0 + 2.0*np.sqrt(5.0/3.0)/11.0)
        x = np.array([-1.0, -b, -a, 0, a, b, 1.0], dtype=np.float64)
        aw = (124.0+7.0*np.sqrt(15.0))/350.0
        bw = (124.0-7.0*np.sqrt(15.0))/350.0
        w = np.array([1.0/21.0, bw, aw, 256.0/525.0, aw, bw, 1.0/21.0], dtype=np.float64)
    else:
        # if npoints < 3 or npoints > 7:
        raise ValueError(' Incorrect number of points.  3 <= npoints <= 7')
    return x, w

def GaussLobattoquad(npoints,fun):
    '''
    Gauss-Lobatto quadrature rule on [-1,1]
    This quadrature includes 
    This quadrature is exact for polynomials of order 2n-3
    Input  : npoints is the number of Gauss Lobatto quadrature points
             fun is the function to approximate
    Output : approximation to the integral of fun(x), -1 <= x <= 1.
    '''
    # obtain auss Lobatto quadrature points and weights
    x, w = GaussLobatto_Points_Weights(npoints)
    # evaluate fun at quadratue points
    fj = fun(x)
    return w.dot(fj)

def GaussLobattoquad_vals(fvals):
    '''
    Gauss-Lobatto quadrature rule on [-1,1]
    This quadrature includes 
    This quadrature is exact for polynomials of order 2n-3
    Input  : values f(xj) where xj are the Gauss Lobatto point for j=1,..., npoints 
    Output : approximation to the integral of f(x), -1 <= x <= 1, using the values f(xj)
    ''' 
    # obtain auss Lobatto quadrature points and weights
    npoints = fvals.size
    _, w = GaussLobatto_Points_Weights(npoints)
    # evaluate fun at quadratue points
    return w.dot(fvals)

def Vandermonde1d(N,r):
    '''
    Vandermonde matrix using Legendre basis
    Input  : N is polynomial degree of Legendre basis, r is point to evaluate at
    Output : Vandermode Matrix V1d = [P0(r) | P1(r) | ... | PN(r)]
    '''
    V1d = np.zeros(r.size, N+1)

    # x, w = GaussLobatto_Points_Weights(npoints)

    for j in range(N+1):
        V1d[:,j] = eval_jacobi(r, 0, 0, j)
    return V1d