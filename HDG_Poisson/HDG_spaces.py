

# 1d 
class L2_DG_space1d:
    def __init__(self, mesh, order=0):
        self.mesh = mesh
        self.order = order
        self.local_basis = 'GaussLobatto'
