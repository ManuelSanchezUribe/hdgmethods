from netgen.meshing import *
import ngsolve

# ============================
def diagonal1(N=2):
    '''
    ---------
    | \     |
    |   \   |
    |     \ |
    ---------
    '''
    ngmesh = Mesh(dim=2)
    # Nodes
    pnums = []
    for i in range(N + 1):
        for j in range(N + 1):
            pnums.append(ngmesh.Add(MeshPoint(Pnt(i / N, j / N, 0))))
    # Elements
    idx_dom = ngmesh.AddRegion("mat", dim=2)
    for j in range(N):
        for i in range(N):
            ngmesh.Add(Element2D(idx_dom, [pnums[i + j * (N + 1)],
                                           pnums[i + (j + 1) * (N + 1)],
                                           pnums[i + 1 + j * (N + 1)]] ))         
            ngmesh.Add(Element2D(idx_dom, [pnums[i + (j + 1) * (N + 1)],
                                           pnums[i + 1 + (j + 1) * (N + 1)],
                                           pnums[i + 1 + j * (N + 1)]]))
    # horizontal boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=1))
        ngmesh.Add(Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=1))
    # vertical boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[i], pnums[i + 1]], index=2))
        ngmesh.Add(Element1D([pnums[i + N * (N + 1)], pnums[i + 1 + N * (N + 1)]], index=2))
    return ngsolve.Mesh(ngmesh)

def diagonal2(N=2):
    '''
    ---------
    |     / |
    |   /   |
    | /     |
    ---------
    '''
    ngmesh = Mesh(dim=2)
    # Nodes
    pnums = []
    for i in range(N + 1):
        for j in range(N + 1):
            pnums.append(ngmesh.Add(MeshPoint(Pnt(i / N, j / N, 0))))
    # Elements
    idx_dom = ngmesh.AddRegion("mat", dim=2)
    for j in range(N):
        for i in range(N):
            ngmesh.Add(Element2D(idx_dom, [pnums[i + j * (N + 1)],
                                           pnums[i + (j + 1) * (N + 1)],
                                           pnums[i + 1 + (j + 1) * (N + 1)]] ) )           
            ngmesh.Add(Element2D(idx_dom, [pnums[i + j * (N + 1)],
                                           pnums[i + 1 + (j + 1) * (N + 1)],
                                           pnums[i + 1 + j * (N + 1)]]))
    # horizontal boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=1))
        ngmesh.Add(Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=1))
    # vertical boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[i], pnums[i + 1]], index=2))
        ngmesh.Add(Element1D([pnums[i + N * (N + 1)], pnums[i + 1 + N * (N + 1)]], index=2))
    return ngsolve.Mesh(ngmesh)

def criscross(N=2):
    '''
    ---------
    | \   / |
    |   X   |
    | /   \ |
    ---------
    '''
    ngmesh = Mesh(dim=2)
    # Nodes
    pnums = []
    for i in range(N + 1):
        for j in range(N + 1):
            pnums.append(ngmesh.Add(MeshPoint(Pnt(i / N, j / N, 0))))
    for i in range(1,N + 1):
        for j in range(1,N + 1):
            pnums.append(ngmesh.Add(MeshPoint(Pnt((i-0.5) / N, (j-0.5) / N, 0))))
    idx_dom = ngmesh.AddRegion("mat", dim=2)
    for j in range(N):
        for i in range(N):
            c = (N+1)**2+i+j*(N)
            ngmesh.Add(Element2D(idx_dom, [pnums[i + j * (N + 1)],
                                           pnums[i + (j + 1) * (N + 1)],
                                           pnums[c]] ) )           
            ngmesh.Add(Element2D(idx_dom, [pnums[i + (j + 1) * (N + 1)],
                                           pnums[i + 1 + (j + 1) * (N + 1)],
                                           pnums[c]] ) )       
            ngmesh.Add(Element2D(idx_dom, [pnums[i + 1 + (j + 1) * (N + 1)],
                                           pnums[i + 1 + j * (N + 1)],
                                           pnums[c]]))
            ngmesh.Add(Element2D(idx_dom, [pnums[i + 1 + j * (N + 1)],
                                           pnums[i + j * (N + 1)],
                                           pnums[c]]))
    # horizontal boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[N + i * (N + 1)],pnums[N + (i + 1) * (N + 1)]], index=1))
        ngmesh.Add(Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=1))
    # vertical boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[i], pnums[i + 1]], index=2))
        ngmesh.Add(Element1D([pnums[i + N * (N + 1)], pnums[i + 1 + N * (N + 1)]], index=2))
    return ngsolve.Mesh(ngmesh)

def criscross_offcenter(N=2, alpha=0.7, beta=0.25):
    '''
    ---------
    | \     |
    |   \   |
    |     X |
    ---------
    '''
    ngmesh = Mesh(dim=2)
    # Nodes
    pnums = []
    for i in range(N + 1):
        for j in range(N + 1):
            pnums.append(ngmesh.Add(MeshPoint(Pnt(i / N, j / N, 0))))
    for i in range(1,N + 1):
        for j in range(1,N + 1):
            pnums.append(ngmesh.Add(MeshPoint(Pnt((i-(1-alpha)) / N, (j-(1-beta)) / N, 0))))
    idx_dom = ngmesh.AddRegion("mat", dim=2)
    for j in range(N):
        for i in range(N):
            c = (N+1)**2+i+j*(N)
            ngmesh.Add(Element2D(idx_dom, [pnums[i + j * (N + 1)],
                                           pnums[i + (j + 1) * (N + 1)],
                                           pnums[c]] ) )           
            ngmesh.Add(Element2D(idx_dom, [pnums[i + (j + 1) * (N + 1)],
                                           pnums[i + 1 + (j + 1) * (N + 1)],
                                           pnums[c]] ) )       
            ngmesh.Add(Element2D(idx_dom, [pnums[i + 1 + (j + 1) * (N + 1)],
                                           pnums[i + 1 + j * (N + 1)],
                                           pnums[c]]))
            ngmesh.Add(Element2D(idx_dom, [pnums[i + 1 + j * (N + 1)],
                                           pnums[i + j * (N + 1)],
                                           pnums[c]]))
    # horizontal boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[N + i * (N + 1)],pnums[N + (i + 1) * (N + 1)]], index=1))
        ngmesh.Add(Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=1))
    # vertical boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[i], pnums[i + 1]], index=2))
        ngmesh.Add(Element1D([pnums[i + N * (N + 1)], pnums[i + 1 + N * (N + 1)]], index=2))
    return ngsolve.Mesh(ngmesh)

def chevron(N=2):
    ngmesh = Mesh(dim=2)
    # Nodes
    pnums = []
    for i in range(N + 1):
        for j in range(N + 1):
            pnums.append(ngmesh.Add(MeshPoint(Pnt(i / N, j / N, 0))))
    # Elements
    idx_dom = ngmesh.AddRegion("mat", dim=2)
    for j in range(N):
        for i in range(N):
            if j%2==0:
                ngmesh.Add(Element2D(idx_dom, [pnums[i + j * (N + 1)],
                                               pnums[i + (j + 1) * (N + 1)],
                                               pnums[i + 1 + (j + 1) * (N + 1)]] ) )           
                ngmesh.Add(Element2D(idx_dom, [pnums[i + j * (N + 1)],
                                               pnums[i + 1 + (j + 1) * (N + 1)],
                                               pnums[i + 1 + j * (N + 1)]]))
            else:
                ngmesh.Add(Element2D(idx_dom, [pnums[i + j * (N + 1)],
                                               pnums[i + (j + 1) * (N + 1)],
                                               pnums[i + 1 + j * (N + 1)]] ))         
                ngmesh.Add(Element2D(idx_dom, [pnums[i + (j + 1) * (N + 1)],
                                               pnums[i + 1 + (j + 1) * (N + 1)],
                                               pnums[i + 1 + j * (N + 1)]]))
    # horizontal boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=1))
        ngmesh.Add(Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=1))
    # vertical boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[i], pnums[i + 1]], index=2))
        ngmesh.Add(Element1D([pnums[i + N * (N + 1)], pnums[i + 1 + N * (N + 1)]], index=2))
    return ngsolve.Mesh(ngmesh)

def unionjack(N=2):
    ngmesh = Mesh(dim=2)
    # Nodes
    pnums = []
    for i in range(N + 1):
        for j in range(N + 1):
            pnums.append(ngmesh.Add(MeshPoint(Pnt(i / N, j / N, 0))))
    # Elements
    idx_dom = ngmesh.AddRegion("mat", dim=2)
    for j in range(N):
        for i in range(N):
            if j%2==i%2:
                ngmesh.Add(Element2D(idx_dom, [pnums[i + j * (N + 1)],
                                               pnums[i + (j + 1) * (N + 1)],
                                               pnums[i + 1 + (j + 1) * (N + 1)]] ) )           
                ngmesh.Add(Element2D(idx_dom, [pnums[i + j * (N + 1)],
                                               pnums[i + 1 + (j + 1) * (N + 1)],
                                               pnums[i + 1 + j * (N + 1)]]))
            else:
                ngmesh.Add(Element2D(idx_dom, [pnums[i + j * (N + 1)],
                                               pnums[i + (j + 1) * (N + 1)],
                                               pnums[i + 1 + j * (N + 1)]] ))         
                ngmesh.Add(Element2D(idx_dom, [pnums[i + (j + 1) * (N + 1)],
                                               pnums[i + 1 + (j + 1) * (N + 1)],
                                               pnums[i + 1 + j * (N + 1)]]))
    # horizontal boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=1))
        ngmesh.Add(Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=1))
    # vertical boundaries
    for i in range(N):
        ngmesh.Add(Element1D([pnums[i], pnums[i + 1]], index=2))
        ngmesh.Add(Element1D([pnums[i + N * (N + 1)], pnums[i + 1 + N * (N + 1)]], index=2))
    return ngsolve.Mesh(ngmesh)
