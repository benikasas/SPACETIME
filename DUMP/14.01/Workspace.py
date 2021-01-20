import numpy as np
# from gauge_latticeqcd import matrix_su2, matrix_su3, create_su3_set, lattice
# from generate import *
# import tools_v1 as tool

# beta=betas[0]
Nt=5
Nz=7
Ny=9
Nx=10
u0=1
# lattice_space
# class spacetime():
#     def __init__(self, Nt, Nx, Ny, Nz, beta, u0, U=None):
#         if None == U:
#             # Creates a spacetime grid, and at every point theres 4 
#             U = [[[[np.zeros(4) for z in range(Nz)] for y in range(Ny)] for x in range(Nx)] for t in range(Nt)]
#         # convert to numpy arrays -> significant speed up
#         self.U = np.array(U)

#         ## Unsure if needed
#         # self.beta = beta
#         # self.u0 = u0
#         # self.Nx = Nx
#         # self.Ny = Ny
#         # self.Nz = Nz
#         # self.Nt = Nt

# # def first_approx():
#     U_og=trial.U(Nt, Nx, Ny, Nz)-trial.U(Nt, Nx, Ny, Nz)
#     U_prime=trial.U(Nt, Nx, Ny, Nz)
#     approx=(U_prime-U_og)/lattice
# Trial=spacetime(Nt, Nx, Ny, Nz, beta, u0)
# print(Trial.U[1, 1, 1, 1])
partial_sum=0
jack=np.zeros((4,4))
jack[1,1]=1
for i in range(4):
    for j in range(4):
        print(jack[i][j])